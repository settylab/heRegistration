"""Stage 2: warp Xenium objects (transcripts, cell + nucleus boundaries)
into H&E pixel space using the registrar pickle from stage 1.

Uses a locally-vendored shim of HEST's `warp_and_save_xenium_objects`
(see `hexenium._internal.hest_warp_shim`): HEST v1.2.0 — the version
pinned in `environments/heRegistration-requirements.txt` — does not
export that function upstream; it was added to HEST main in the
"hest V2" refactor (commit 1bb61159, Feb 2026). The shim composes
v1.2.0's `warp_gdf_valis` primitive to cover the cells + nuclei
targets that the hexenium default (`warp.targets: [cells, nuclei]`)
uses. The transcripts target raises a clear NotImplementedError from
the shim; enabling it requires bumping the HEST pin.

The Dask LocalCluster with JVMPlugin remains — it initialises the
BioFormats JVM on each worker for future dask-aware code paths; on
v1.2.0's non-dask `warp_gdf_valis`, the shim inits JVM in the main
process itself.
"""
from __future__ import annotations

import contextlib
from pathlib import Path

from hexenium._internal.compat import apply_numpy_shims
from hexenium._internal.logging import log
from hexenium.manifest import (
    compute_params_hash,
    git_sha,
    set_default_symlink,
    utc_timestamp,
    write_manifest,
)


VALID_TARGETS = {"cells", "nuclei", "transcripts"}


@contextlib.contextmanager
def _dask_with_jvm(
    n_workers: int,
    threads_per_worker: int,
    memory_limit: str,
    jvm_mem_gb: int,
):
    """Start a local Dask cluster whose workers each init the VALIS JVM
    via a WorkerPlugin. Tear everything down on exit."""
    import dask
    from dask.distributed import LocalCluster, Client, WorkerPlugin

    class JVMPlugin(WorkerPlugin):
        def __init__(self, mem_gb: int = 1):
            self.mem_gb = mem_gb

        def setup(self, worker):
            import jpype  # noqa: F401
            from valis_hest.registration import init_jvm  # type: ignore
            if not jpype.isJVMStarted():
                init_jvm(mem_gb=self.mem_gb)

    dask.config.set({"distributed.scheduler.worker-ttl": None})
    cluster = LocalCluster(
        n_workers=n_workers,
        threads_per_worker=threads_per_worker,
        memory_limit=memory_limit,
        dashboard_address=None,
        nanny=False,
    )
    client = Client(cluster)
    client.register_worker_plugin(JVMPlugin(mem_gb=jvm_mem_gb), name="jvm")
    try:
        log(f"[warp] dask cluster: {client}; dashboard={client.dashboard_link}")
        yield client
    finally:
        try:
            client.close()
        finally:
            cluster.close()


def run_warp(
    sample_id: str,
    registrar_pickle: Path,
    xenium_bundle: Path,
    out_dir: Path,
    *,
    dapi_path: Path | None = None,
    targets: list[str] | None = None,
    use_dask: bool = True,
    save_geojson: bool = True,
    dask_n_workers: int = 1,
    dask_threads_per_worker: int = 1,
    dask_memory_limit: str = "32GB",
    dask_jvm_mem_gb: int = 1,
    force_rerun: bool = False,
    source_register_run_id: str | None = None,
) -> Path:
    """Run warp. Returns the warp output directory.

    ``out_dir`` is the layout-computed
    ``<output_root_he>/warp/<he_job_id>/`` folder
    (caller resolves this from the RunLayout).

    ``source_register_run_id`` is the ``<he_job_id>`` of the
    ``register/<X>/`` directory whose registrar pickle was consumed.
    Recorded in the per-run warp manifest so downstream
    ``set-default-run`` can validate register↔warp lineage before
    re-pointing user-facing symlinks. When ``None`` the manifest omits
    the lineage field — expected in tests and in older no-flag flows.
    """
    apply_numpy_shims()
    # HEST v1.2.0 (installed pin) does NOT export
    # `warp_and_save_xenium_objects`; use the local shim that composes
    # v1.2.0's `warp_gdf_valis`. See module docstring above.
    from hexenium._internal.hest_warp_shim import warp_and_save_xenium_objects

    if targets is None:
        targets = ["cells", "nuclei", "transcripts"]
    bad = sorted(set(targets) - VALID_TARGETS)
    if bad:
        raise ValueError(f"unknown warp targets: {bad}; valid: {sorted(VALID_TARGETS)}")

    out_dir.mkdir(parents=True, exist_ok=True)

    # Resume guard — skip if every requested target's sentinel exists.
    target_files = {
        "cells": out_dir / "he_cell_seg.parquet",
        "nuclei": out_dir / "he_nucleus_seg.parquet",
        "transcripts": out_dir / "he_transcripts.parquet",
    }
    needed = [target_files[t] for t in targets]
    if all(p.exists() for p in needed) and not force_rerun:
        log(f"[warp] skip — all requested target files exist in {out_dir}")
        return out_dir

    # The "moving image" basename HEST passes through to VALIS's slide_dict
    # lookup. Must match the actual DAPI filename used in registration.
    # If dapi_path was explicitly set (e.g. multichannel ch0000_dapi.ome.tif),
    # use its basename; otherwise fall back to the standard Xenium layout.
    from hexenium.stages.registration import resolve_dapi_path
    dapi_resolved = resolve_dapi_path(xenium_bundle, explicit_path=dapi_path)
    dapi_filename = dapi_resolved.name

    kwargs = dict(
        path_registrar=str(registrar_pickle),
        # HEST 1.1.1 names this parameter `dapi_path` in the public signature
        # but uses it INTERNALLY as `curr_slide_name` — the slide_dict key VALIS
        # uses to look up the moving image. That key is the basename of the
        # file passed to Valis() during registration (here = the symlinked or
        # actual DAPI filename). The pipeline computed exactly that value
        # above; we just need to pass it under the public name `dapi_path`.
        dapi_path=dapi_filename,
        save_dir=str(out_dir),
        use_dask=use_dask,
        verbose=True,
        save_geojson=save_geojson,
    )

    if "cells" in targets:
        kwargs["dapi_cells"] = str(xenium_bundle / "cell_boundaries.parquet")
    if "nuclei" in targets:
        kwargs["dapi_nuclei"] = str(xenium_bundle / "nucleus_boundaries.parquet")
    if "transcripts" in targets:
        kwargs["dapi_transcripts"] = str(xenium_bundle / "transcripts.parquet")

    # Trim kwargs the installed function doesn't accept (defensive — HEST
    # versions have varied).
    import inspect
    sig = inspect.signature(warp_and_save_xenium_objects)
    if "path_registrar" not in sig.parameters and "registrar_pickle_path" in sig.parameters:
        kwargs["registrar_pickle_path"] = kwargs.pop("path_registrar")
    accepted = {k: v for k, v in kwargs.items() if k in sig.parameters}
    dropped = sorted(set(kwargs) - set(accepted))
    if dropped:
        log(f"[warp] WARN: dropped unsupported kwargs (HEST version skew): {dropped}")

    log(f"[warp] targets resolved: {targets}; outputs will land in {out_dir}")
    if use_dask:
        log(f"[warp] starting Dask LocalCluster "
            f"(n_workers={dask_n_workers}, threads={dask_threads_per_worker}, "
            f"memory_limit={dask_memory_limit}, jvm_mem_gb={dask_jvm_mem_gb})")
        with _dask_with_jvm(
            n_workers=dask_n_workers,
            threads_per_worker=dask_threads_per_worker,
            memory_limit=dask_memory_limit,
            jvm_mem_gb=dask_jvm_mem_gb,
        ):
            log(f"[warp] calling warp_and_save_xenium_objects() "
                f"(typically 20-60 min for cells+nuclei; hours if transcripts included)")
            warp_and_save_xenium_objects(**accepted)
            log(f"[warp] warp_and_save_xenium_objects() returned")
    else:
        # JVM still required even without dask.
        from valis_hest.registration import init_jvm  # type: ignore
        import jpype  # noqa: F401
        log(f"[warp] use_dask=False; initialising JVM directly")
        if not jpype.isJVMStarted():
            init_jvm(mem_gb=dask_jvm_mem_gb)
        log(f"[warp] calling warp_and_save_xenium_objects()")
        warp_and_save_xenium_objects(**accepted)
        log(f"[warp] warp_and_save_xenium_objects() returned")

    # Sanity-check outputs landed.
    missing = [str(p) for p in needed if not p.exists()]
    if missing:
        raise RuntimeError(
            f"warp_and_save_xenium_objects completed but expected outputs are missing: {missing}"
        )

    _write_run_manifests(
        out_dir=out_dir,
        sample_id=sample_id,
        targets=targets,
        use_dask=use_dask,
        save_geojson=save_geojson,
        source_register_run_id=source_register_run_id,
        registrar_pickle=registrar_pickle,
    )

    log(f"[warp] done -> {out_dir}")
    return out_dir


def _write_run_manifests(
    *,
    out_dir: Path,
    sample_id: str,
    targets: list[str],
    use_dask: bool,
    save_geojson: bool,
    source_register_run_id: str | None,
    registrar_pickle: Path,
) -> None:
    """Persist warp's per-run manifest + stage-root default symlink.

    Mirror of the register stage's writer: an immutable per-run
    manifest at ``<out_dir>/manifest.yaml`` plus a relative symlink at
    ``<container>/manifest.yaml -> <he_job_id>/manifest.yaml``. The
    lineage field ``source_register_run_id`` is what
    ``set-default-run`` reads to reject register/warp mismatches.
    """
    he_job_id = out_dir.name
    container = out_dir.parent
    warp_params = {
        "targets": sorted(targets),
        "use_dask": use_dask,
        "save_geojson": save_geojson,
    }
    payload = {
        "sample_id": sample_id,
        "he_job_id": he_job_id,
        "run_dir": str(out_dir),
        "source_register_run_id": source_register_run_id,
        "registrar_pickle": str(registrar_pickle),
        "params": warp_params,
        "params_hash": compute_params_hash(warp_params),
        "git_sha": git_sha(),
        "timestamp_utc": utc_timestamp(),
    }
    write_manifest(out_dir / "manifest.yaml", payload)
    set_default_symlink(container / "manifest.yaml",
                        f"{he_job_id}/manifest.yaml")
