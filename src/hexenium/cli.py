"""Command-line entrypoint for hexenium.

Exposes ``hexenium run …`` via the ``[project.scripts]`` entry point in
pyproject.toml.

Three invocation modes:

* **standalone** — ``--sample-id`` + ``--he-path`` + ``--xenium-bundle`` +
  ``--output-root`` (+ ``--proseg-purified-h5ad`` for the celltype
  stage). Outputs land at ``<output_root>/<sample_id>/{...}``.

* **integrated-by-run-id** — add ``--run-id`` on top of the standalone
  set. Derives the xenium h5ad path from the upstream layout
  (``<output-root>/<sample>/<sample>_<run-id>/spatial_adata/
  <sample>_xenium_ranger.h5ad``) and colocates H&E outputs under
  ``<output-root>/<sample>/<sample>_<run-id>/he_registration/``.

* **integrated-by-h5ad** — pass ``--xenium-h5ad`` directly. Sample
  identity read from ``.uns['sample_id']`` + ``.uns['run_id']``;
  outputs colocate under ``<xenium_run_dir>/he_registration/``.
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
from hexenium._internal.logging import log


def _str2bool(v):
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in {"true", "1", "yes", "y", "on"}:
        return True
    if s in {"false", "0", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError(f"expected boolean, got {v!r}")


def _add_run_args(p: argparse.ArgumentParser) -> None:
    """Attach every `run` subcommand flag to `p`."""
    p.add_argument("--config", type=Path, default=None,
                   help="User config YAML (overrides config/default.yaml entries).")
    p.add_argument("--stages", nargs="+", choices=VALID_STAGES, default=list(DEFAULT_STAGES),
                   help="Which stages to run.")
    # I/O — three modes: standalone (sample_id + output_root),
    # integrated-by-run-id (add --run-id), integrated-by-h5ad
    # (--xenium-h5ad).
    p.add_argument("--sample-id",
                   help="Sample identifier (e.g. SAMPLE1). "
                        "Standalone / integrated-by-run-id: required. "
                        "Integrated-by-h5ad: auto from .uns.")
    p.add_argument("--run-id", default=None,
                   help="Upstream xenium-preprocess run id (e.g. demo_v1). "
                        "When set with --sample-id + --output-root, derives the "
                        "xenium h5ad from <output-root>/<sample>/<sample>_<run-id>/"
                        "spatial_adata/<sample>_xenium_ranger.h5ad and colocates "
                        "H&E outputs under <output-root>/<sample>/<sample>_<run-id>/"
                        "he_registration/.")
    p.add_argument("--he-path", "--he-slide", dest="he_path", type=Path,
                   help="Path to H&E image (any name). Alias: --he-slide.")
    p.add_argument("--xenium-bundle", type=Path,
                   help="Xenium output-XETG... directory with morphology_focus/, "
                        "transcripts.parquet, cell_boundaries.parquet, "
                        "nucleus_boundaries.parquet.")
    p.add_argument("--xenium-h5ad", type=Path, default=None,
                   help="Xenium h5ad written by an upstream xenium-preprocess "
                        "xenium_ranger_to_anndata stage (reads .uns[sample_id]/[run_id]). "
                        "When set, activates integrated mode + colocates outputs "
                        "under <xenium_run_dir>/he_registration/. Also becomes the "
                        "celltype-stage QUERY side. Usually pass --run-id instead "
                        "so the path is derived.")
    p.add_argument("--proseg-purified-h5ad", type=Path, default=None,
                   help="Upstream proseg_purified h5ad. Its .obs carries celltypes "
                        "+ centroids that get NN-mapped onto every xenium cell. "
                        "Auto-derived from <xenium_run_dir>/spatial_adata/"
                        "<sample>_proseg_purified.h5ad when unset AND either "
                        "--xenium-h5ad or --run-id is passed (integrated modes). "
                        "Required in standalone mode when celltype is in --stages.")
    p.add_argument("--dapi-path", type=Path, default=None,
                   help="Full path to the DAPI/morphology image to register against. "
                        "If omitted, derives from <xenium_bundle>/morphology_focus/"
                        "morphology_focus_0000.ome.tif.")
    p.add_argument("--output-root", type=Path,
                   help="Root output directory (standalone / integrated-by-run-id "
                        "modes).")
    p.add_argument("--he-job-id", default=None,
                   help="Explicit run label for this invocation (overrides "
                        "$SLURM_JOB_ID). Interactive fallback: YYYYMMDDTHHMMSS "
                        "timestamp.")
    p.add_argument("--register-run-id", default=None,
                   help="When running warp against a specific prior registration, "
                        "point at register/<id>/manifest.yaml directly instead of "
                        "the register/manifest.yaml default-symlink. Use to compare "
                        "downstream results across multiple registrations of the "
                        "same sample.")
    p.add_argument("--warp-run-id", default=None,
                   help="When running celltype/viz without a preceding warp in "
                        "this invocation, point at a prior warp/<id>/ folder.")
    p.add_argument("--celltype-run-id", default=None,
                   help="When running viz without a preceding celltype in this "
                        "invocation, point viz at a prior celltyped/<id>/ folder.")
    p.add_argument("--force-rerun", action="store_true",
                   help="Re-run all stages even if sentinel outputs exist.")
    p.add_argument("--force-preprocess", action="store_true",
                   help="Force re-conversion of VSI → OME-TIFF even if a "
                        "converted OME-TIFF already exists at the canonical "
                        "next-to-source location (<vsi_dir>/<vsi_stem>.ome.tif) "
                        "or the pipeline-tree fallback. Only affects "
                        "he_preprocess; unlike --force-rerun does not re-run "
                        "register/warp/celltype/viz.")
    # H&E
    p.add_argument("--symlink-to-canonical-name", type=_str2bool, default=None,
                   help="Symlink input H&E to aligned_fullres_HE.<ext> in workdir "
                        "(sidesteps HEST hardcode).")
    # Registration
    p.add_argument("--mode", choices=("rigid_only", "rigid_nonrigid", "full_with_micro"),
                   default=None,
                   help="Registration mode. NOTE: rigid_only_micro was removed — "
                        "register_micro requires bk_dxdy from an initial non-rigid "
                        "pass; combining rigid-only init with micro raises "
                        "TypeError in valis_hest.")
    p.add_argument("--use-he-deconvolution", type=_str2bool, default=None,
                   help="Apply Macenko-style HEDeconvolution to the H&E "
                        "(YAML default: true; override with false on "
                        "faint-hematoxylin samples).")
    p.add_argument("--check-for-reflections", type=_str2bool, default=None)
    p.add_argument("--create-masks", type=_str2bool, default=None)
    p.add_argument("--align-to-reference", type=_str2bool, default=None)
    p.add_argument("--max-image-dim-px", type=int, default=None,
                   help="Max dimension for the SAVED image pyramid (VALIS "
                        "memory guard; default 1500).")
    p.add_argument("--max-processed-image-dim-px", type=int, default=None,
                   help="Max dimension for feature-detection downsample (default 1500).")
    p.add_argument("--max-non-rigid-registration-dim-px", type=int, default=None,
                   help="Max dimension for non-rigid stage (default 10000).")
    p.add_argument("--run-name", default=None,
                   help="Optional registrar run name for the manifest "
                        "(provenance only; on-disk folder is always <he_job_id>).")
    # Warp
    p.add_argument("--warp-targets", nargs="+",
                   choices=("cells", "nuclei", "transcripts"), default=None,
                   help="Which Xenium objects to warp. Default: cells + nuclei.")
    p.add_argument("--include-transcripts", action="store_true", default=False,
                   help="Append 'transcripts' to the warp targets (opt-in; "
                        "transcripts are expensive to warp).")
    p.add_argument("--use-dask", type=_str2bool, default=None,
                   help="Use Dask LocalCluster + JVMPlugin for warp.")
    # Celltype — inlined proseg → xenium NN mapping.
    p.add_argument("--celltype-col", default=None,
                   help="Column on proseg_purified.h5ad's .obs to use as celltype "
                        "label. Default 'auto' — precedence celltype > first_type > "
                        "primary_cell_type > celltype_updated.")
    p.add_argument("--id-col", default=None,
                   help="Column on xenium.h5ad's .obs holding xenium UUIDs. "
                        "Default 'auto' — shape-check first, prefers UUID-shaped "
                        "candidates. Special value '__index__' forces the index.")
    p.add_argument("--nn-k", type=int, default=None,
                   help="Neighbours per query in the proseg → xenium NN mapping.")
    # Viz
    p.add_argument("--viz-enabled", type=_str2bool, default=None)
    p.add_argument("--viz-thumbnail-max-dim", type=int, default=None)
    p.add_argument("--viz-dpi", type=int, default=None)
    p.add_argument("--viz-render-boundaries",
                   choices=("cell", "nucleus", "both"), default=None,
                   help="Which boundaries to draw on the overlay.")


def _add_set_default_run_args(p: argparse.ArgumentParser) -> None:
    """Argparse wiring for ``hexenium set-default-run``.

    Identity resolution mirrors ``run``: pass ``--sample-id`` +
    ``--output-root`` (standalone), optionally with ``--run-id`` for
    integrated-by-run-id mode. Only one ``--*-run-id`` flag is
    required; the rest keep their current pointer (partial update).
    """
    p.add_argument("--sample-id", required=True,
                   help="Sample identifier — same value used at pipeline run time.")
    p.add_argument("--run-id", default=None,
                   help="Upstream xenium-preprocess run id (integrated-by-run-id "
                        "mode). Omit for standalone-mode layouts.")
    p.add_argument("--output-root", required=True, type=Path,
                   help="Root output directory (same value passed at run time).")
    p.add_argument("--register-run-id", default=None,
                   help="Retarget output/register to register/<he_job_id>/.")
    p.add_argument("--warp-run-id", default=None,
                   help="Retarget output/warp to warp/<he_job_id>/.")
    p.add_argument("--celltyped-run-id", default=None,
                   help="Retarget output/celltyped to celltyped/<he_job_id>/.")
    p.add_argument("--viz-run-id", default=None,
                   help="Retarget output/viz to viz/<he_job_id>/.")
    p.add_argument("--force-lineage", action="store_true",
                   help="Bypass the register↔warp lineage-consistency refusal. "
                        "Use only when deliberately cross-comparing.")


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

    sdr_p = sub.add_parser(
        "set-default-run",
        help="Re-point output/<stage> symlinks after visual QA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_set_default_run_args(sdr_p)
    return parser


def _resolve_config(args: argparse.Namespace) -> dict:
    """Load default + optional user YAML, then apply CLI overrides.

    Also derives the xenium h5ad path from ``--run-id`` when
    ``--xenium-h5ad`` isn't set (integrated-by-run-id mode).
    """
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
    if args.xenium_h5ad is not None:
        overrides["xenium_h5ad"] = str(args.xenium_h5ad)
    if args.proseg_purified_h5ad is not None:
        overrides["proseg_purified_h5ad"] = str(args.proseg_purified_h5ad)
    if args.dapi_path is not None:
        overrides["dapi_path"] = str(args.dapi_path)
    if args.output_root is not None:
        overrides["output_root"] = str(args.output_root)
    if args.he_job_id is not None:
        overrides["he_job_id"] = args.he_job_id
    if args.register_run_id is not None:
        overrides["register_run_id"] = args.register_run_id
    if args.warp_run_id is not None:
        overrides["warp_run_id"] = args.warp_run_id
    if args.celltype_run_id is not None:
        overrides["celltype_run_id"] = args.celltype_run_id
    if args.force_rerun:
        overrides["force_rerun"] = True
    if args.force_preprocess:
        overrides["force_preprocess"] = True
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
    if args.max_image_dim_px is not None:
        param_over["max_image_dim_px"] = args.max_image_dim_px
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
    if args.celltype_col is not None:
        ct_over["celltype_col"] = args.celltype_col
    if args.id_col is not None:
        ct_over["id_col"] = args.id_col
    if args.nn_k is not None:
        ct_over["nn_k"] = args.nn_k
    if ct_over:
        overrides["celltype"] = ct_over
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

    # If --run-id was passed without --xenium-h5ad, derive the h5ad path
    # from the canonical upstream xenium-preprocess layout:
    #     <output_root>/<sample>/<sample>_<run_id>/spatial_adata/<sample>_xenium_ranger.h5ad
    # Only the celltype stage actually reads the h5ad's content; register,
    # warp, viz, he_preprocess don't. So the existence check is required
    # only when celltype is in --stages — otherwise register+warp can run
    # before the upstream step-1 has produced the h5ad, and identity
    # comes from the CLI args via _resolve_integrated_by_run_id.
    if args.run_id and not cfg.get("xenium_h5ad"):
        sid = cfg.get("sample_id")
        oroot = cfg.get("output_root")
        if not (sid and oroot):
            raise SystemExit(
                "--run-id requires --sample-id and --output-root "
                "(or use --xenium-h5ad directly)."
            )
        derived = (Path(oroot) / sid / f"{sid}_{args.run_id}"
                   / "spatial_adata" / f"{sid}_xenium_ranger.h5ad")
        cfg["run_id"] = args.run_id  # for resolve_layout no-h5ad branch
        if derived.exists():
            cfg["xenium_h5ad"] = str(derived)
        elif "celltype" in args.stages:
            raise SystemExit(
                f"derived xenium h5ad does not exist: {derived}\n"
                f"expected upstream xenium-preprocess step-1 output at\n"
                f"  <output-root>/<sample>/<sample>_<run-id>/spatial_adata/"
                f"<sample>_xenium_ranger.h5ad\n"
                f"the celltype stage requires this file. Run the upstream "
                f"step-1 first, or drop celltype from --stages (register+warp+viz "
                f"don't need it).\n"
                f"Check --sample-id={sid!r} / --run-id={args.run_id!r} / "
                f"--output-root={oroot!r}, or pass --xenium-h5ad explicitly."
            )
        else:
            log(f"[pipeline] xenium h5ad not present at {derived} — "
                f"proceeding without it (celltype not in --stages).")

    validate(cfg)
    return cfg


def _resolve_set_default_run_output_root_he(args: argparse.Namespace) -> Path:
    """Compute the ``he_registration/`` root from the identity flags.

    Mirrors the layout module's two integrated/standalone shapes
    (skipping the h5ad-driven mode — ``set-default-run`` is a
    lightweight symlink-toggle command, not worth another loader).
    """
    output_root = Path(args.output_root).resolve()
    sample_id = args.sample_id
    if args.run_id:
        return output_root / sample_id / f"{sample_id}_{args.run_id}" / "he_registration"
    return output_root / sample_id


def main(argv: list[str] | None = None) -> int:
    import sys

    args = build_parser().parse_args(argv)
    if args.cmd == "run":
        from hexenium.pipeline import run
        cfg = _resolve_config(args)
        return run(cfg, stages=list(args.stages), argv=sys.argv)
    if args.cmd == "set-default-run":
        from hexenium.set_default_run import format_changes, set_default_run
        output_root_he = _resolve_set_default_run_output_root_he(args)
        if not output_root_he.exists():
            raise SystemExit(
                f"set-default-run: he_registration/ tree not found at "
                f"{output_root_he}. Check --sample-id/--run-id/--output-root."
            )
        result = set_default_run(
            output_root_he=output_root_he,
            register_run_id=args.register_run_id,
            warp_run_id=args.warp_run_id,
            celltyped_run_id=args.celltyped_run_id,
            viz_run_id=args.viz_run_id,
            force_lineage=args.force_lineage,
        )
        print(format_changes(result))
        return 0
    raise SystemExit(f"unknown command: {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
