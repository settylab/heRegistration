"""Stage 1: H&E ↔ Xenium DAPI registration.

Wraps VALIS (via valis_hest) directly rather than HEST's
`register_dapi_he`, because the stock function exposes only a subset of
VALIS kwargs (missing `max_processed_image_dim_px`,
`use_he_deconvolution`, `create_masks`, `align_to_reference`,
`non_rigid_registrar_cls`).

Filename trap: valis_hest's `_post_register` hardcodes the H&E slide-dict
key as `'aligned_fullres_HE'`. If the H&E basename is anything else, the
post-register step KeyErrors AFTER successful registration. This module
symlinks the input H&E to that canonical name inside the per-sample
workdir before invoking VALIS, so the operator never has to rename
source files.
"""
from __future__ import annotations

import inspect
import time
from pathlib import Path

from hexenium._internal.compat import apply_numpy_shims  # kill_jvm_safely intentionally not imported — see note in run_registration
from hexenium._internal.logging import log
from hexenium.manifest import (
    compute_params_hash,
    git_sha,
    set_default_symlink,
    utc_timestamp,
    write_manifest,
)


def _canonical_he_symlink(he_path: Path, workdir: Path) -> Path:
    """Create a symlink named aligned_fullres_HE.<ext> pointing at he_path.
    Returns the symlink path."""
    workdir.mkdir(parents=True, exist_ok=True)
    suffix = "".join(he_path.suffixes) or ".ome.tif"
    canonical = workdir / f"aligned_fullres_HE{suffix}"
    if canonical.exists() or canonical.is_symlink():
        canonical.unlink()
    canonical.symlink_to(he_path.resolve())
    return canonical


def resolve_dapi_path(
    xenium_bundle: Path,
    explicit_path: Path | None = None,
) -> Path:
    """Resolve the DAPI/morphology image used as the moving image.

    Precedence:
      1. explicit_path (from --dapi-path CLI flag or YAML dapi_path key).
      2. <xenium_bundle>/morphology_focus/morphology_focus_0000.ome.tif
         (standard Xenium output bundle layout).
    """
    if explicit_path is not None:
        if not explicit_path.exists():
            raise FileNotFoundError(
                f"DAPI path passed via --dapi-path does not exist: {explicit_path}"
            )
        return explicit_path
    dapi = xenium_bundle / "morphology_focus" / "morphology_focus_0000.ome.tif"
    if not dapi.exists():
        raise FileNotFoundError(
            f"DAPI image not found at {dapi}. "
            "Either pass --xenium-bundle pointing at the output-XETG... directory "
            "that holds morphology_focus/morphology_focus_0000.ome.tif, or "
            "pass --dapi-path pointing at the actual DAPI file (e.g. multichannel "
            "bundles use morphology_focus/ch0000_dapi.ome.tif)."
        )
    return dapi


# Backward-compat alias for internal callers.
_derive_dapi_path = resolve_dapi_path


def run_registration(
    sample_id: str,
    he_path: Path,
    xenium_bundle: Path,
    out_dir: Path,
    *,
    dapi_path: Path | None = None,
    mode: str = "rigid_only_micro",
    use_he_deconvolution: bool = False,
    check_for_reflections: bool = True,
    create_masks: bool = False,
    align_to_reference: bool = True,
    max_image_dim_px: int = 1500,
    max_processed_image_dim_px: int = 1500,
    max_non_rigid_registration_dim_px: int = 10000,
    micro_rigid_registrar_cls=None,
    micro_rigid_registrar_params: dict | None = None,
    name: str | None = None,
    symlink_to_canonical_name: bool = True,
    force_rerun: bool = False,
) -> Path:
    """Run registration. Returns the path to `_registrar.pickle`.

    `mode`:
      - rigid_only        : rigid solve, no nonrigid, no micro (fastest).
      - rigid_nonrigid    : rigid + nonrigid.
      - full_with_micro   : rigid + nonrigid + micro (stock HEST behavior).

    NOTE: `rigid_only_micro` was removed on 2026-06-30. Combining
    `non_rigid_registrar_cls=None` in the Valis() constructor (so the
    initial register() pass skips non-rigid) with a subsequent
    register_micro() call raises
    `TypeError: 'NoneType' object is not subscriptable` at
    valis_hest/registration.py:4809, because register_micro tries to
    combine `slide_obj.bk_dxdy` (the displacement field from the
    initial pass) with the new micro displacements — and bk_dxdy is
    None when the initial pass was rigid-only.

    ``out_dir`` is the layout-computed per-run register dir
    (``<output_root_he>/register/<he_job_id>``). Caller (pipeline.run)
    resolves this from the RunLayout.
    """
    apply_numpy_shims()
    from valis_hest import preprocessing, registration  # type: ignore
    from valis_hest.slide_io import BioFormatsSlideReader  # type: ignore
    from hest.SlideReaderAdapter import SlideReaderAdapter  # type: ignore
    from hest.utils import get_name_datetime  # type: ignore

    if micro_rigid_registrar_params is None:
        micro_rigid_registrar_params = {}

    if mode not in {"rigid_only", "rigid_nonrigid", "full_with_micro"}:
        raise ValueError(
            f"unknown registration mode: {mode!r}. "
            f"Valid modes: rigid_only, rigid_nonrigid, full_with_micro. "
            f"(rigid_only_micro was removed — it triggers a NoneType bug in "
            f"valis_hest's register_micro because bk_dxdy isn't populated by "
            f"a rigid-only initial pass.)"
        )
    do_nonrigid = mode in {"rigid_nonrigid", "full_with_micro"}
    do_micro = mode == "full_with_micro"

    workdir = out_dir / "_workdir"
    # container = the folder that carries the cross-run manifest.yaml.
    # This is <output_root_he>/register/ (parent of the per-he_job_id dirs).
    container = out_dir.parent

    # H&E filename rewrite — sidesteps valis_hest's "Paul fix" hardcode.
    if symlink_to_canonical_name and he_path.stem != "aligned_fullres_HE":
        he_for_valis = _canonical_he_symlink(he_path, workdir)
        log(f"[registration] symlinked H&E -> {he_for_valis} (real: {he_path})")
    else:
        he_for_valis = he_path

    dapi_path = resolve_dapi_path(xenium_bundle, explicit_path=dapi_path)

    # Resume check.
    # ``run_name`` is retained purely for provenance in the manifest —
    # the on-disk folder is <he_job_id>, so this label is diagnostic,
    # not path-carrying.
    run_name = name or f"{sample_id}_{get_name_datetime()}"
    registrar_dir = out_dir  # caller pre-encodes <he_job_id>
    he_job_id = registrar_dir.name
    sentinel = registrar_dir / "data" / "_registrar.pickle"

    # Params consumed by the registrar — used for both the manifest
    # snapshot AND the params_hash. Kept as a single dict so the hash
    # is stable across re-runs with identical params.
    reg_params = {
        "mode": mode,
        "use_he_deconvolution": use_he_deconvolution,
        "check_for_reflections": check_for_reflections,
        "create_masks": create_masks,
        "align_to_reference": align_to_reference,
        "max_image_dim_px": max_image_dim_px,
        "max_processed_image_dim_px": max_processed_image_dim_px,
        "max_non_rigid_registration_dim_px": max_non_rigid_registration_dim_px,
    }
    params_hash = compute_params_hash(reg_params)

    if sentinel.exists() and not force_rerun:
        log(f"[registration] skip — registrar already exists at {sentinel}")
        # Persist a manifest pointer if missing.
        _write_run_manifests(
            registrar_dir=registrar_dir,
            container=container,
            payload={
                "sample_id": sample_id,
                "he_job_id": he_job_id,
                "run_dir": str(registrar_dir),
                "registrar_pickle": str(sentinel),
                "skipped_resume": True,
                "params": reg_params,
                "params_hash": params_hash,
                "git_sha": git_sha(),
                "timestamp_utc": utc_timestamp(),
            },
        )
        return sentinel

    registrar_dir.mkdir(parents=True, exist_ok=True)
    container.mkdir(parents=True, exist_ok=True)

    valis_kwargs = dict(
        src_dir="",
        dst_dir=str(registrar_dir),
        reference_img_f=str(he_for_valis),
        img_list=[str(he_for_valis), str(dapi_path)],
        check_for_reflections=check_for_reflections,
        create_masks=create_masks,
        align_to_reference=align_to_reference,
        max_image_dim_px=max_image_dim_px,
        max_processed_image_dim_px=max_processed_image_dim_px,
        max_non_rigid_registration_dim_px=max_non_rigid_registration_dim_px,
        micro_rigid_registrar_cls=micro_rigid_registrar_cls,
        micro_rigid_registrar_params=micro_rigid_registrar_params,
    )
    if not do_nonrigid:
        valis_kwargs["non_rigid_registrar_cls"] = None

    # Filter to kwargs the installed Valis() actually accepts — versions drift.
    valis_sig = inspect.signature(registration.Valis)
    accepted_init = {k: v for k, v in valis_kwargs.items() if k in valis_sig.parameters}
    dropped_init = sorted(set(valis_kwargs) - set(accepted_init))
    if dropped_init:
        log(f"[registration] WARN: VALIS rejected init kwargs (version skew): {dropped_init}")

    log(f"[registration] mode={mode} (do_nonrigid={do_nonrigid}, do_micro={do_micro})")
    log(f"[registration] he_for_valis={he_for_valis}")
    log(f"[registration] dapi_path={dapi_path}")
    log(f"[registration] init_jvm() ...")
    registration.init_jvm()
    log(f"[registration] JVM started; constructing Valis(...)")

    # NOTE on JVM lifecycle: we deliberately do NOT call kill_jvm_safely()
    # after registration completes. JPype enforces ONE JVM lifecycle per
    # Python process — once kill_jvm is called, scyjava.start_jvm raises
    # `OSError: JVM cannot be restarted`. Downstream stages (warp) call
    # HEST's warp_gdf_valis which unconditionally re-runs init_jvm; that
    # OSErrors after registration kills the JVM. So no try/finally for
    # JVM cleanup here; errors still propagate as expected.
    registrar = registration.Valis(**accepted_init)
    log(f"[registration] Valis() constructed")

    register_kwargs = {
        "reader_dict": {
            str(he_for_valis): [SlideReaderAdapter],
            str(dapi_path): [BioFormatsSlideReader],
        },
    }
    if use_he_deconvolution:
        register_kwargs["brightfield_processing_cls"] = preprocessing.HEDeconvolution
        log(f"[registration] using HEDeconvolution for brightfield processing")
    else:
        log(f"[registration] brightfield processing OFF (raw RGB-grayscale path)")

    log(f"[registration] calling registrar.register() — this includes "
        f"image conversion + processing + rigid"
        f"{' + nonrigid' if do_nonrigid else ''} (typically 30-60 min)")
    t0 = time.time()
    registrar.register(**register_kwargs)
    log(f"[registration] registrar.register() done in {time.time() - t0:.1f}s")

    if do_micro:
        micro_kwargs = {
            "max_non_rigid_registration_dim_px": max_non_rigid_registration_dim_px,
            "align_to_reference": True,
            "reference_img_f": str(he_for_valis),
        }
        if use_he_deconvolution:
            micro_kwargs["brightfield_processing_cls"] = preprocessing.HEDeconvolution
        log(f"[registration] calling registrar.register_micro() at "
            f"max_non_rigid_registration_dim_px={max_non_rigid_registration_dim_px} "
            f"(typically 30 min – 2 h; this is usually the long pole)")
        t0 = time.time()
        registrar.register_micro(**micro_kwargs)
        log(f"[registration] registrar.register_micro() done in {time.time() - t0:.1f}s")

    if not sentinel.exists():
        raise RuntimeError(
            f"VALIS reported success but {sentinel} does not exist. "
            "Inspect the registrar_dir for partial output."
        )

    _write_run_manifests(
        registrar_dir=registrar_dir,
        container=container,
        payload={
            "sample_id": sample_id,
            "he_job_id": he_job_id,
            "run_name": run_name,
            "run_dir": str(registrar_dir),
            "he_input_real": str(he_path),
            "he_passed_to_valis": str(he_for_valis),
            "dapi_path": str(dapi_path),
            "registrar_pickle": str(sentinel),
            "params": reg_params,
            "params_hash": params_hash,
            "git_sha": git_sha(),
            "timestamp_utc": utc_timestamp(),
            "dropped_valis_kwargs_init": dropped_init,
        },
    )
    log(f"[registration] done -> {sentinel}")
    return sentinel


def _write_run_manifests(*, registrar_dir: Path, container: Path, payload: dict) -> None:
    """Write the per-run manifest and re-point the stage-root symlink.

    Layout:
      ``<container>/<he_job_id>/manifest.yaml`` — this run's immutable
      manifest; never overwritten by later register runs.

      ``<container>/manifest.yaml`` — relative symlink pointing at the
      per-run manifest; each successful register run rewrites this to
      itself. Warp reads through the symlink when ``--register-run-id``
      is not passed, so this pointer IS the current default.
    """
    per_run = registrar_dir / "manifest.yaml"
    write_manifest(per_run, payload)
    set_default_symlink(container / "manifest.yaml",
                        f"{registrar_dir.name}/manifest.yaml")
