"""Command-line entrypoint for hexenium.

Exposes `hexenium run …` via the `[project.scripts]` entry point in
pyproject.toml.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from hexenium import __version__
from hexenium.config import (
    DEFAULT_STAGES,
    VALID_STAGES,
    deep_update,
    load_default,
    load_yaml,
    validate,
)


def _str2bool(v):
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in {"true", "1", "yes", "y", "on"}:
        return True
    if s in {"false", "0", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError(f"expected boolean, got {v!r}")


def _str2bool_or_auto(v):
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s == "auto":
        return "auto"
    return _str2bool(v)


def _add_run_args(p: argparse.ArgumentParser) -> None:
    """Attach every `run` subcommand flag to `p`. Kept separate so the
    argparse layout mirrors run_pipeline.py's flag surface 1:1."""
    p.add_argument("--config", type=Path, default=None,
                   help="User config YAML (overrides config/default.yaml entries).")
    p.add_argument("--stages", nargs="+", choices=VALID_STAGES, default=list(DEFAULT_STAGES),
                   help="Which stages to run. Default excludes nn_celltype_mapping "
                        "(opt-in — requires an upstream proseg purified.h5ad).")
    # I/O
    p.add_argument("--sample-id", help="Sample identifier (e.g. SAMPLE1).")
    p.add_argument("--he-path", type=Path, help="Path to H&E image (any name).")
    p.add_argument("--xenium-bundle", type=Path,
                   help="Xenium output-XETG... directory with morphology_focus/, transcripts.parquet, etc.")
    p.add_argument("--dapi-path", type=Path, default=None,
                   help="Full path to the DAPI/morphology image to register against. "
                        "If omitted, derives from <xenium_bundle>/morphology_focus/morphology_focus_0000.ome.tif. "
                        "Set this for non-standard bundles, e.g. multichannel "
                        "morphology_focus/ch0000_dapi.ome.tif.")
    p.add_argument("--output-root", type=Path, help="Root output directory.")
    p.add_argument("--force-rerun", action="store_true",
                   help="Re-run all stages even if sentinel outputs exist.")
    # H&E
    p.add_argument("--symlink-to-canonical-name", type=_str2bool, default=None,
                   help="Symlink input H&E to aligned_fullres_HE.<ext> in workdir (sidesteps HEST hardcode).")
    # Registration
    p.add_argument("--mode", choices=("rigid_only", "rigid_nonrigid", "full_with_micro"),
                   default=None,
                   help="Registration mode. NOTE: rigid_only_micro was removed — "
                        "register_micro requires bk_dxdy from an initial non-rigid "
                        "pass; combining rigid-only init with micro raises "
                        "TypeError: 'NoneType' object is not subscriptable in "
                        "valis_hest/registration.py:4809.")
    p.add_argument("--use-he-deconvolution", type=_str2bool, default=None,
                   help="Apply Macenko-style HEDeconvolution to the H&E "
                        "(YAML default: true; override with false on "
                        "faint-hematoxylin samples).")
    p.add_argument("--check-for-reflections", type=_str2bool, default=None,
                   help="Try mirror/rotation orientations during rigid solve.")
    p.add_argument("--create-masks", type=_str2bool, default=None,
                   help="Let VALIS auto-mask tissue regions.")
    p.add_argument("--align-to-reference", type=_str2bool, default=None,
                   help="Treat H&E as the registration reference (default true).")
    p.add_argument("--max-processed-image-dim-px", type=int, default=None,
                   help="Max dimension for feature-detection downsample (default 1500).")
    p.add_argument("--max-non-rigid-registration-dim-px", type=int, default=None,
                   help="Max dimension for non-rigid stage (default 10000).")
    p.add_argument("--run-name", default=None,
                   help="Optional registrar run name; default <sample>_<datetime>.")
    # Warp
    p.add_argument("--warp-targets", nargs="+",
                   choices=("cells", "nuclei", "transcripts"), default=None,
                   help="Which Xenium objects to warp. Default: cells + nuclei.")
    p.add_argument("--include-transcripts", action="store_true", default=False,
                   help="Append 'transcripts' to the warp targets (opt-in; "
                        "transcripts are expensive to warp).")
    p.add_argument("--use-dask", type=_str2bool, default=None,
                   help="Use Dask LocalCluster + JVMPlugin for warp.")
    # Celltype
    p.add_argument("--celltype-csv", type=Path, default=None,
                   help="Per-sample celltype annotation CSV.")
    p.add_argument("--csv-id-col", default=None, help="Column in celltype CSV holding the Xenium cell_id.")
    p.add_argument("--csv-group-col", default=None, help="Column in celltype CSV holding the cell-type label.")
    # NN celltype mapping (proseg → Xenium)
    p.add_argument("--nn-proseg-source", default=None,
                   help="Proseg-side source (.h5ad or .csv[.gz]). Relative path is "
                        "joined onto <output_root>/<sample_id>/. "
                        "Default template: h5ad/{sample_id}_purified.h5ad.")
    p.add_argument("--nn-celltype-col", default=None,
                   help="Column in the proseg source holding the celltype label. Default: first_type.")
    p.add_argument("--nn-k", type=int, default=None,
                   help="Number of nearest neighbours. Default: 1.")
    p.add_argument("--nn-distance-threshold", type=float, default=None,
                   help="Distance beyond which no assignment is considered valid (Xenium µm). "
                        "Default: no threshold (matches notebook).")
    p.add_argument("--nn-unmatched-policy",
                   choices=("drop", "mark_unassigned", "keep"), default=None,
                   help="What to do with cells beyond nn-distance-threshold. "
                        "Ignored when threshold is unset.")
    p.add_argument("--nn-hybrid-direct-join-first", type=_str2bool_or_auto, default=None,
                   help="Hybrid direct-join policy: 'auto' (default; enable when "
                        "proseg's original_cell_id column is populated AND covers a "
                        "sufficient fraction of Xenium cells), true (force enable), "
                        "or false (force disable, pure NN).")
    # Viz
    p.add_argument("--viz-enabled", type=_str2bool, default=None)
    p.add_argument("--viz-thumbnail-max-dim", type=int, default=None)
    p.add_argument("--viz-dpi", type=int, default=None)
    p.add_argument("--viz-render-boundaries",
                   choices=("cell", "nucleus", "both"), default=None,
                   help="Which boundaries to draw on the overlay: "
                        "cell (polygons only), nucleus (outlines only), "
                        "or both. Default (from YAML): nucleus.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hexenium",
        description="H&E ↔ Xenium DAPI registration + warp + celltype + viz pipeline.",
    )
    parser.add_argument("--version", action="version", version=f"hexenium {__version__}")
    sub = parser.add_subparsers(dest="cmd", required=True)

    run_p = sub.add_parser(
        "run",
        help="Run the registration + warp + celltype + viz pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_run_args(run_p)
    return parser


def _resolve_config(args: argparse.Namespace) -> dict:
    """Load default + optional user YAML, then apply CLI overrides."""
    cfg = load_default()
    if args.config is not None:
        cfg = deep_update(cfg, load_yaml(args.config))

    # Map CLI overrides into the nested config tree.
    overrides: dict = {}
    if args.sample_id is not None:
        overrides["sample_id"] = args.sample_id
    if args.he_path is not None:
        overrides["he_path"] = str(args.he_path)
    if args.xenium_bundle is not None:
        overrides["xenium_bundle"] = str(args.xenium_bundle)
    if args.dapi_path is not None:
        overrides["dapi_path"] = str(args.dapi_path)
    if args.output_root is not None:
        overrides["output_root"] = str(args.output_root)
    if args.force_rerun:
        overrides["force_rerun"] = True
    he_over = {}
    if args.symlink_to_canonical_name is not None:
        he_over["symlink_to_canonical_name"] = args.symlink_to_canonical_name
    if he_over:
        overrides["he"] = he_over
    reg_over = {}
    if args.mode is not None:
        reg_over["mode"] = args.mode
    if args.use_he_deconvolution is not None:
        reg_over["use_he_deconvolution"] = args.use_he_deconvolution
    if args.check_for_reflections is not None:
        reg_over["check_for_reflections"] = args.check_for_reflections
    if args.create_masks is not None:
        reg_over["create_masks"] = args.create_masks
    if args.align_to_reference is not None:
        reg_over["align_to_reference"] = args.align_to_reference
    if reg_over:
        overrides["registration"] = reg_over
    param_over = {}
    if args.max_processed_image_dim_px is not None:
        param_over["max_processed_image_dim_px"] = args.max_processed_image_dim_px
    if args.max_non_rigid_registration_dim_px is not None:
        param_over["max_non_rigid_registration_dim_px"] = args.max_non_rigid_registration_dim_px
    if args.run_name is not None:
        param_over["name"] = args.run_name
    if param_over:
        overrides["parameters"] = param_over
    warp_over = {}
    if args.warp_targets is not None:
        warp_over["targets"] = list(args.warp_targets)
    if args.include_transcripts:
        warp_over["include_transcripts"] = True
    if args.use_dask is not None:
        warp_over["use_dask"] = args.use_dask
    if warp_over:
        overrides["warp"] = warp_over
    ct_over = {}
    if args.celltype_csv is not None:
        ct_over["csv"] = str(args.celltype_csv)
    if args.csv_id_col is not None:
        ct_over["csv_id_col"] = args.csv_id_col
    if args.csv_group_col is not None:
        ct_over["csv_group_col"] = args.csv_group_col
    if ct_over:
        overrides["celltype"] = ct_over
    nn_in_over: dict = {}
    nn_alg_over: dict = {}
    if args.nn_proseg_source is not None:
        nn_in_over["proseg_source"] = args.nn_proseg_source
    if args.nn_celltype_col is not None:
        nn_in_over["proseg_celltype_col"] = args.nn_celltype_col
    if args.nn_k is not None:
        nn_alg_over["k"] = args.nn_k
    if args.nn_distance_threshold is not None:
        nn_alg_over["distance_threshold"] = args.nn_distance_threshold
    if args.nn_unmatched_policy is not None:
        nn_alg_over["unmatched_policy"] = args.nn_unmatched_policy
    if args.nn_hybrid_direct_join_first is not None:
        nn_alg_over["hybrid_direct_join_first"] = args.nn_hybrid_direct_join_first
    if nn_in_over or nn_alg_over:
        nn_over: dict = {}
        if nn_in_over:
            nn_over["input"] = nn_in_over
        if nn_alg_over:
            nn_over["nn"] = nn_alg_over
        overrides["nn_celltype_mapping"] = nn_over
    viz_over = {}
    if args.viz_enabled is not None:
        viz_over["enabled"] = args.viz_enabled
    if args.viz_thumbnail_max_dim is not None:
        viz_over["thumbnail_max_dim"] = args.viz_thumbnail_max_dim
    if args.viz_dpi is not None:
        viz_over["dpi"] = args.viz_dpi
    if args.viz_render_boundaries is not None:
        viz_over["render_boundaries"] = args.viz_render_boundaries
    if viz_over:
        overrides["viz"] = viz_over

    cfg = deep_update(cfg, overrides)
    validate(cfg)
    return cfg


def main(argv: list[str] | None = None) -> int:
    import sys

    args = build_parser().parse_args(argv)
    if args.cmd == "run":
        from hexenium.pipeline import run
        cfg = _resolve_config(args)
        return run(cfg, stages=list(args.stages), argv=sys.argv)
    raise SystemExit(f"unknown command: {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
