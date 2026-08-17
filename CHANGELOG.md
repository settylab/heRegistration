# Changelog

All notable changes to hexenium will be documented here. Follows
[Keep a Changelog](https://keepachangelog.com/) shape.

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
