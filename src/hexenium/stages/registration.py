"""Stage 1: H&E ↔ Xenium DAPI registration.

Wraps VALIS (via valis_hest) directly rather than HEST's
`register_dapi_he`, because the stock function exposes only a subset of
VALIS kwargs (missing `max_processed_image_dim_px`,
`use_he_deconvolution`, `create_masks`, `align_to_reference`,
`non_rigid_registrar_cls`). This wrapper surfaces the full set as
config parameters.

Filename trap: valis_hest's `_post_register` hardcodes the H&E slide-dict
key as `'aligned_fullres_HE'`. If the H&E basename is anything else, the
post-register step KeyErrors AFTER successful registration. This module
symlinks the input H&E to that canonical name inside the per-sample
workdir before invoking VALIS, so users never have to rename source
files.
"""
from __future__ import annotations

import inspect
import time
from pathlib import Path

from hexenium._internal.compat import apply_numpy_shims  # kill_jvm_safely intentionally not imported — see note in run_registration
from hexenium._internal.logging import log


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
    output_root: Path,
    *,
    dapi_path: Path | None = None,
    mode: str = "rigid_only_micro",
    use_he_deconvolution: bool = False,
    check_for_reflections: bool = True,
    create_masks: bool = False,
    align_to_reference: bool = True,
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
      - rigid_only        : rigid solve, no nonrigid, no micro
      - rigid_nonrigid    : rigid + nonrigid
      - full_with_micro   : rigid + nonrigid + micro (stock HEST behavior; default)

    NOTE: `rigid_only_micro` was removed. Combining
    `non_rigid_registrar_cls=None` in the Valis() constructor (so the
    initial register() pass skips non-rigid) with a subsequent
    register_micro() call raises
    `TypeError: 'NoneType' object is not subscriptable` at
    valis_hest/registration.py:4809, because register_micro tries to
    combine `slide_obj.bk_dxdy` (the displacement field from the
    initial pass) with the new micro displacements — and bk_dxdy is
    None when the initial pass was rigid-only.
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

    out_dir = output_root / sample_id / "registration"
    workdir = out_dir / "_workdir"

    # H&E filename rewrite — sidesteps valis_hest's "Paul fix" hardcode.
    if symlink_to_canonical_name and he_path.stem != "aligned_fullres_HE":
        he_for_valis = _canonical_he_symlink(he_path, workdir)
        log(f"[registration] symlinked H&E -> {he_for_valis} (real: {he_path})")
    else:
        he_for_valis = he_path

    dapi_path = resolve_dapi_path(xenium_bundle, explicit_path=dapi_path)

    # Resume check.
    run_name = name or f"{sample_id}_{get_name_datetime()}"
    out_dir.mkdir(parents=True, exist_ok=True)
    registrar_dir = out_dir / run_name
    sentinel = registrar_dir / "data" / "_registrar.pickle"
    if sentinel.exists() and not force_rerun:
        log(f"[registration] skip — registrar already exists at {sentinel}")
        # Persist a manifest pointer if missing.
        _write_manifest(out_dir / "manifest.yaml", {
            "sample_id": sample_id,
            "registrar_pickle": str(sentinel),
            "skipped_resume": True,
        })
        return sentinel

    registrar_dir.mkdir(parents=True, exist_ok=True)

    valis_kwargs = dict(
        src_dir="",
        dst_dir=str(registrar_dir),
        reference_img_f=str(he_for_valis),
        img_list=[str(he_for_valis), str(dapi_path)],
        check_for_reflections=check_for_reflections,
        create_masks=create_masks,
        align_to_reference=align_to_reference,
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
    # OSErrors after registration kills the JVM. The right pattern is:
    # init the JVM once at the start of stage 1, reuse across stages,
    # let Python's process-exit teardown shut it down. So no try/finally
    # for JVM cleanup here; errors still propagate as expected.
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

    _write_manifest(out_dir / "manifest.yaml", {
        "sample_id": sample_id,
        "run_name": run_name,
        "he_input_real": str(he_path),
        "he_passed_to_valis": str(he_for_valis),
        "dapi_path": str(dapi_path),
        "registrar_pickle": str(sentinel),
        "mode": mode,
        "use_he_deconvolution": use_he_deconvolution,
        "check_for_reflections": check_for_reflections,
        "create_masks": create_masks,
        "align_to_reference": align_to_reference,
        "max_processed_image_dim_px": max_processed_image_dim_px,
        "max_non_rigid_registration_dim_px": max_non_rigid_registration_dim_px,
        "dropped_valis_kwargs_init": dropped_init,
    })
    log(f"[registration] done -> {sentinel}")
    return sentinel


def _write_manifest(path: Path, payload: dict) -> None:
    import yaml  # PyYAML
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        yaml.safe_dump(payload, f, sort_keys=False)
