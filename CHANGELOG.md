# Changelog

All notable changes to hexenium will be documented here. Follows
[Keep a Changelog](https://keepachangelog.com/) shape.

## [0.2.0] — 2026-08-17

### Added
- New `hexenium.layout` module: `RunLayout`, `resolve_layout`,
  `resolve_he_job_id`, `compute_stages_suffix`, `STAGE_NAMES`,
  `is_xenium_uuid`. Encodes the output-tree convention (per-stage
  `<he_job_id>` sharding) and the three-mode invocation resolver.
- Three invocation modes: **standalone**, **integrated-by-run-id**
  (`--run-id`), and **integrated-by-h5ad** (`--xenium-h5ad`).
  Integrated modes colocate H&E outputs under an upstream
  xenium-preprocess run dir at `<run_dir>/he_registration/`.
- Inlined proseg → xenium NN mapping in `hexenium.stages.celltyping`.
  Auto-detects the celltype column (`celltype` > `first_type` >
  `primary_cell_type` > `celltype_updated`) and the xenium id column
  (UUID-shape-checked to reject Proseg's int64 index masquerading as
  `cell_id`). Auto-detects spatial coords via `.obsm['spatial']` and
  the common `.obs` conventions (`centroid_x/y`, `x_centroid/y_centroid`,
  `x/y`).
- OME-TIFF discovery in `hexenium.stages.he_preprocess`
  (`discover_existing_ometiff`, `_next_to_source_path`,
  `_fallback_out_path`): re-runs skip the VSI → OME-TIFF conversion when
  the converted file already exists next to the source or under the
  pipeline-tree fallback. New `--force-preprocess` flag overrides.
- New CLI flags: `--run-id`, `--xenium-h5ad`, `--proseg-purified-h5ad`,
  `--he-job-id`, `--warp-run-id`, `--celltype-run-id`,
  `--force-preprocess`, `--celltype-col`, `--id-col`, `--nn-k`,
  `--max-image-dim-px`, `--he-slide` alias for `--he-path`.
- `submit_he_registration.sh` rewrite: two-phase `sbatch --hold`
  submission with layout-integrated log routing (colocated with any
  upstream run's step logs under `<run_dir>/logs/`), plus the
  defensive `.bashrc` → micromamba/mamba/conda activation fallback
  chain from v0.1.0. The stages-suffix (`<sample>_<jobid>_<stages>`)
  records which pipeline stages ran, so re-runs at the same jobid
  don't clobber each other's logs.

### Changed
- Stage entrypoints (`run_registration`, `run_warp`, `run_celltyping`,
  `run_viz`, `run_he_preprocess`) now take a pre-computed `out_dir`
  instead of `output_root` + `sample_id`. The pipeline resolves paths
  from the `RunLayout` before calling each stage.
- Output tree switches from `he_preprocessed/ | registration/ |
  warped/ | celltyped/ | viz/` (flat, unsharded) to `converted/ |
  register/<he_job_id>/ | warp/<he_job_id>/ | celltyped/<he_job_id>/ |
  viz/<he_job_id>/` (job-id-sharded past stage 0).
- Config: mode-aware validation — `sample_id` + `output_root` are only
  required when `xenium_h5ad` isn't set (integrated-by-h5ad mode gets
  them from `.uns`).
- Default `parameters.max_image_dim_px = 1500` matches
  `max_processed_image_dim_px` so valis_hest's auto-bump warning
  doesn't fire on the stock config.
- Registration mode `full_with_micro` remains the default; the removed
  `rigid_only_micro` mode is documented in the docstring.

### Removed
- CSV-based celltype join (`--celltype-csv`, `--csv-id-col`,
  `--csv-group-col`). Celltype now sources labels directly from the
  proseg_purified h5ad via inlined NN. Legacy runs relying on the CSV
  need to re-source their labels from proseg_purified.
- Opt-in `nn_celltype_mapping` stage and module. Its role (emit a
  CSV that the CSV-based celltype path reads) is subsumed by the
  inlined NN mapping in `celltyping.py`.

### Migration notes
- Callers passing `--celltype-csv` need to switch to `--xenium-h5ad`
  + `--proseg-purified-h5ad` (or rely on integrated mode's
  auto-derivation from the upstream run dir).
- Stage entry-point signatures changed from `output_root` + `sample_id`
  to `out_dir`; direct callers of `run_registration`/`run_warp`/…
  (outside the pipeline) need to update.
- Output paths moved: `he_preprocessed/` → next-to-source or
  `converted/`; `registration/` → `register/<he_job_id>/`; `warped/` →
  `warp/<he_job_id>/`; `celltyped/` → `celltyped/<he_job_id>/`;
  `viz/` → `viz/<he_job_id>/`.

## [0.1.0] — 2026-07-09

### Added
- Initial package layout — src-layout Python package
  (`src/hexenium/`) with an installable console-script entry
  point (`hexenium`).
- Five stage modules under `hexenium.stages` (he_preprocess,
  registration, warp, celltyping, viz).
- Optional `nn_celltype_mapping` stage for propagating per-cell
  labels from an external CSV onto the warped Xenium boundaries
  via nearest-neighbour lookup.
- Config machinery split into `hexenium.config` (YAML load + deep
  merge + validation).
- CLI split into `hexenium.cli` (argparse) + `hexenium.pipeline`
  (stage-orchestration loop).
- `pyproject.toml` + `environments/heRegistration.yml` pinning
  the tested dependency set.
- `CITATION.cff`, MIT `LICENSE`, `docs/` markdown pages, and a
  suite of unit + equivalence smoke tests under `tests/`.
- `scripts/submit_he_registration.sh` — sbatch wrapper carrying a
  Slurm-context-aware log-routing fix that avoids a `tee` + `set
  -o pipefail` deadlock under sbatch.
