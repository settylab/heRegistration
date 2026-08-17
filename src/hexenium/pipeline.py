"""Stage-orchestration loop for hexenium.

Runs the requested subset of `he_preprocess → register → warp →
celltype → viz`, each guarded by its own sentinel-existence resume check
(nuke with `force_rerun`).
"""
from __future__ import annotations

import os
import socket
import sys
import time
from pathlib import Path

import yaml

from hexenium._internal.logging import log, banner


def log_invocation_banner(argv: list[str], stages: list[str], cfg: dict) -> None:
    """Print everything anyone would want to grep out of the log if the
    job ever silently does nothing again."""
    banner("hexenium starting")
    log(f"command:      {' '.join(argv)}")
    log(f"host:         {socket.gethostname()}")
    log(f"cwd:          {os.getcwd()}")
    log(f"pid:          {os.getpid()}")
    log(f"python:       {sys.executable}")
    log(f"python ver:   {sys.version.splitlines()[0]}")
    log(f"sample_id:    {cfg.get('sample_id')}")
    log(f"he_path:      {cfg.get('he_path')}")
    log(f"xenium_bun:   {cfg.get('xenium_bundle')}")
    log(f"dapi_path:    {cfg.get('dapi_path') or '(derived from xenium_bundle)'}")
    log(f"output_root:  {cfg.get('output_root')}")
    log(f"stages:       {stages}")
    log(f"force_rerun:  {cfg.get('force_rerun', False)}")
    log(f"mode:         {cfg['registration']['mode']}")
    log(f"use_he_deconv:{cfg['registration']['use_he_deconvolution']}")
    log(f"warp.targets: {cfg['warp']['targets']}"
        f"{' +transcripts' if cfg['warp'].get('include_transcripts') else ''}")
    # Library versions — wrap each in its own try so one missing library
    # doesn't suppress the rest.
    for name in ("numpy", "yaml", "scipy", "scanpy", "anndata",
                 "geopandas", "shapely", "dask", "tifffile", "openslide",
                 "hest", "valis_hest"):
        try:
            mod = __import__(name)
            ver = getattr(mod, "__version__", "no __version__ attr")
            log(f"  {name:14s}{ver}")
        except Exception as e:
            log(f"  {name:14s}NOT IMPORTABLE ({type(e).__name__})")
    log("=" * 64)


def run(cfg: dict, stages: list[str], argv: list[str]) -> int:
    """Execute the pipeline. Assumes `cfg` has already been validated."""
    sample_id = cfg["sample_id"]
    he_path = Path(cfg["he_path"]).resolve()
    xenium_bundle = Path(cfg["xenium_bundle"]).resolve()
    dapi_path_cfg = cfg.get("dapi_path")
    dapi_path = Path(dapi_path_cfg).resolve() if dapi_path_cfg else None
    output_root = Path(cfg["output_root"]).resolve()
    force_rerun = bool(cfg.get("force_rerun", False))

    log_invocation_banner(argv, stages, cfg)

    if not he_path.exists():
        raise SystemExit(f"H&E not found: {he_path}")
    if not xenium_bundle.exists():
        raise SystemExit(f"xenium bundle not found: {xenium_bundle}")

    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / sample_id).mkdir(parents=True, exist_ok=True)

    # Snapshot resolved config to the run dir for traceability.
    snap = output_root / sample_id / "resolved_config.yaml"
    with open(snap, "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    log(f"[pipeline] resolved config snapshot -> {snap}")

    n_stages = len(stages)
    registrar_pickle = None
    warp_dir = output_root / sample_id / "warped"
    celltyped_dir = output_root / sample_id / "celltyped"
    pipeline_t0 = time.time()

    # Stage 0: VSI → OME-TIFF (no-op if input isn't a .vsi file). Runs
    # first because it substitutes he_path for the downstream stages.
    if "he_preprocess" in stages:
        idx = stages.index("he_preprocess") + 1
        banner(f"stage {idx}/{n_stages}: he_preprocess — starting")
        t0 = time.time()
        from hexenium.stages.he_preprocess import run_he_preprocess
        hp_cfg = cfg.get("he_preprocess", {}) or {}
        new_he_path = run_he_preprocess(
            sample_id=sample_id,
            he_path=he_path,
            output_root=output_root,
            tile_size=hp_cfg.get("tile_size", 1024),
            compression=hp_cfg.get("compression", "jpeg2000"),
            jpeg2000_level=hp_cfg.get("jpeg2000_level", 100),
            pyramid_levels=hp_cfg.get("pyramid_levels", 7),
            force_rerun=force_rerun,
        )
        if new_he_path != he_path:
            log(f"[pipeline] register stage will use converted OME-TIFF: {new_he_path}")
            he_path = new_he_path
        banner(f"stage {idx}/{n_stages}: he_preprocess — complete in {time.time()-t0:.1f}s")
    elif he_path.suffix.lower() == ".vsi":
        # User skipped he_preprocess but passed a .vsi — warn.
        log(f"[pipeline] WARN: --he-path is a .vsi file but stage he_preprocess "
            f"is not in --stages. Register will likely fail; add "
            f"'--stages he_preprocess register warp celltype viz' to run stage 0.")

    if "register" in stages:
        idx = stages.index("register") + 1
        banner(f"stage {idx}/{n_stages}: register — starting")
        t0 = time.time()
        from hexenium.stages.registration import run_registration
        reg = cfg["registration"]
        params = cfg["parameters"]
        he_cfg = cfg["he"]
        registrar_pickle = run_registration(
            sample_id=sample_id,
            he_path=he_path,
            xenium_bundle=xenium_bundle,
            output_root=output_root,
            dapi_path=dapi_path,
            mode=reg["mode"],
            use_he_deconvolution=reg["use_he_deconvolution"],
            check_for_reflections=reg["check_for_reflections"],
            create_masks=reg["create_masks"],
            align_to_reference=reg["align_to_reference"],
            max_processed_image_dim_px=params["max_processed_image_dim_px"],
            max_non_rigid_registration_dim_px=params["max_non_rigid_registration_dim_px"],
            micro_rigid_registrar_cls=params.get("micro_rigid_registrar_cls"),
            micro_rigid_registrar_params=params.get("micro_rigid_registrar_params") or {},
            name=params.get("name"),
            symlink_to_canonical_name=he_cfg["symlink_to_canonical_name"],
            force_rerun=force_rerun,
        )
        banner(f"stage {idx}/{n_stages}: register — complete in {time.time()-t0:.1f}s")
    else:
        # Try to discover an existing registrar from the manifest.
        manifest = output_root / sample_id / "registration" / "manifest.yaml"
        if manifest.exists():
            with open(manifest) as f:
                m = yaml.safe_load(f)
            registrar_pickle = Path(m["registrar_pickle"])
            log(f"[pipeline] using existing registrar from manifest: {registrar_pickle}")

    if "warp" in stages:
        idx = stages.index("warp") + 1
        banner(f"stage {idx}/{n_stages}: warp — starting")
        t0 = time.time()
        from hexenium.stages.warp import run_warp
        if registrar_pickle is None or not registrar_pickle.exists():
            raise SystemExit(
                "warp stage requested but no registrar pickle is available. "
                "Run the register stage first or place a manifest at "
                f"{output_root / sample_id / 'registration' / 'manifest.yaml'}."
            )
        warp_cfg = cfg["warp"]
        # Compose final targets list: start from configured list, append
        # 'transcripts' when include_transcripts is set, then dedupe while
        # preserving order.
        targets = list(warp_cfg["targets"])
        if warp_cfg.get("include_transcripts") and "transcripts" not in targets:
            targets.append("transcripts")
        seen = set()
        targets = [t for t in targets if not (t in seen or seen.add(t))]
        log(f"[warp] final targets: {targets}")
        warp_dir = run_warp(
            sample_id=sample_id,
            registrar_pickle=registrar_pickle,
            xenium_bundle=xenium_bundle,
            output_root=output_root,
            dapi_path=dapi_path,
            targets=targets,
            use_dask=warp_cfg["use_dask"],
            save_geojson=warp_cfg["save_geojson"],
            dask_n_workers=warp_cfg["dask"]["n_workers"],
            dask_threads_per_worker=warp_cfg["dask"]["threads_per_worker"],
            dask_memory_limit=warp_cfg["dask"]["memory_limit"],
            dask_jvm_mem_gb=warp_cfg["dask"]["jvm_mem_gb"],
            force_rerun=force_rerun,
        )
        banner(f"stage {idx}/{n_stages}: warp — complete in {time.time()-t0:.1f}s")

    if "celltype" in stages:
        idx = stages.index("celltype") + 1
        banner(f"stage {idx}/{n_stages}: celltype — starting")
        t0 = time.time()
        from hexenium.stages.celltyping import run_celltyping
        ct_cfg = cfg["celltype"]
        csv = ct_cfg.get("csv")
        ct_result = run_celltyping(
            sample_id=sample_id,
            warp_dir=warp_dir,
            output_root=output_root,
            celltype_csv=Path(csv) if csv else None,
            csv_id_col=ct_cfg["csv_id_col"],
            csv_group_col=ct_cfg["csv_group_col"],
            nuclei_inherit_classification=ct_cfg["nuclei_inherit_classification"],
            area_threshold_px=ct_cfg["area_threshold_px"],
            nucleus_round_ndigits=ct_cfg["nucleus_round_ndigits"],
            cell_round_ndigits=ct_cfg["cell_round_ndigits"],
            qupath_cell_round_ndigits=ct_cfg.get("qupath_cell_round_ndigits", 2),
            qupath_nucleus_round_ndigits=ct_cfg.get("qupath_nucleus_round_ndigits", 2),
            force_rerun=force_rerun,
        )
        celltyped_dir = ct_result["geojson"].parent
        banner(f"stage {idx}/{n_stages}: celltype — complete in {time.time()-t0:.1f}s")

    if "nn_celltype_mapping" in stages:
        idx = stages.index("nn_celltype_mapping") + 1
        banner(f"stage {idx}/{n_stages}: nn_celltype_mapping — starting")
        t0 = time.time()
        from hexenium.stages.nn_celltype_mapping import run_nn_celltype_mapping
        nn_cfg = cfg.get("nn_celltype_mapping", {}) or {}
        nn_in = nn_cfg.get("input", {}) or {}
        nn_alg = nn_cfg.get("nn", {}) or {}
        nn_coords = nn_cfg.get("coords", {}) or {}
        nn_out = nn_cfg.get("output", {}) or {}
        run_nn_celltype_mapping(
            sample_id=sample_id,
            xenium_bundle=xenium_bundle,
            output_root=output_root,
            proseg_source=nn_in.get("proseg_source", "h5ad/{sample_id}_purified.h5ad"),
            proseg_celltype_col=nn_in.get("proseg_celltype_col", "first_type"),
            proseg_x_col=nn_in.get("proseg_x_col", "x"),
            proseg_y_col=nn_in.get("proseg_y_col", "y"),
            proseg_id_col=nn_in.get("proseg_id_col"),
            proseg_original_id_col=nn_in.get("proseg_original_id_col", "original_cell_id"),
            xenium_cells_csv=nn_in.get("xenium_cells_csv", "cells.csv.gz"),
            xenium_id_col=nn_in.get("xenium_id_col", "cell_id"),
            xenium_x_col=nn_in.get("xenium_x_col", "x_centroid"),
            xenium_y_col=nn_in.get("xenium_y_col", "y_centroid"),
            k=nn_alg.get("k", 1),
            algorithm=nn_alg.get("algorithm", "auto"),
            metric=nn_alg.get("metric", "euclidean"),
            tiebreak=nn_alg.get("tiebreak", "min_dist"),
            distance_threshold=nn_alg.get("distance_threshold"),
            unmatched_policy=nn_alg.get("unmatched_policy", "mark_unassigned"),
            unassigned_label=nn_alg.get("unassigned_label", "Unassigned"),
            hybrid_direct_join_first=nn_alg.get("hybrid_direct_join_first", "auto"),
            hybrid_min_original_id_populated=(
                (nn_alg.get("hybrid_auto_detect") or {}).get(
                    "min_original_id_populated", 0.8
                )
            ),
            hybrid_min_xenium_overlap=(
                (nn_alg.get("hybrid_auto_detect") or {}).get(
                    "min_xenium_overlap", 0.5
                )
            ),
            proseg_scale_factor=nn_coords.get("proseg_scale_factor", 1.0),
            xenium_scale_factor=nn_coords.get("xenium_scale_factor", 1.0),
            verify_same_frame=nn_coords.get("verify_same_frame", True),
            frame_check_tolerance=nn_coords.get("frame_check_tolerance", 100.0),
            output_dir=nn_out.get("dir", "celltype_for_hexenium"),
            output_stem=nn_out.get("stem", "{sample_id}_celltype"),
            emit_inspection=nn_out.get("emit_inspection", True),
            inspection_suffix=nn_out.get("inspection_suffix", "_inspection"),
            write_compression=nn_out.get("write_compression"),
            extra_columns=nn_out.get("extra_columns"),
            force_rerun=force_rerun,
        )
        banner(f"stage {idx}/{n_stages}: nn_celltype_mapping — complete in {time.time()-t0:.1f}s")

    if "viz" in stages and cfg["viz"]["enabled"]:
        idx = stages.index("viz") + 1
        banner(f"stage {idx}/{n_stages}: viz — starting")
        t0 = time.time()
        from hexenium.stages.viz import run_viz
        v = cfg["viz"]
        run_viz(
            sample_id=sample_id,
            he_path=he_path,
            warp_dir=warp_dir,
            celltyped_dir=celltyped_dir,
            output_root=output_root,
            thumbnail_max_dim=v["thumbnail_max_dim"],
            dpi=v["dpi"],
            cell_alpha=v["cell_alpha"],
            nucleus_alpha=v["nucleus_alpha"],
            classification_palette=v.get("classification_palette"),
            palette_cmap=v.get("palette_cmap", "tab20"),
            render_boundaries=v.get("render_boundaries", "nucleus"),
            force_rerun=force_rerun,
        )
        banner(f"stage {idx}/{n_stages}: viz — complete in {time.time()-t0:.1f}s")

    banner(f"hexenium done in {time.time()-pipeline_t0:.1f}s "
           f"-> {output_root / sample_id}/")
    return 0
