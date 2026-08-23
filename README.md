# hexenium

**H&E ↔ Xenium registration pipeline.** One CLI invocation per
sample. Aligns a hematoxylin-and-eosin whole-slide image to a
matched Xenium in-situ transcriptomics run. Warps Xenium cell +
nucleus segmentations into H&E pixel space. Renders a per-slide
celltype-annotated overlay.

- **CLI**: `hexenium run …` (Python package + console script).
- **Env**: pinned conda + pip layer at
  `environments/heRegistration.yml`. Default env name
  `heRegistration`.
- **HPC**: `scripts/submit_he_registration.sh …` self-submits
  to Slurm. Handles VSI inputs via an `afterok:` conversion
  chain. Routes logs into the run tree.

Under the hood: [VALIS](https://github.com/MathOnco/valis)
(via `valis_hest`) solves the registration; HEST does the
boundary warp. Five stages
(`he_preprocess → register → warp → celltype → viz`). Each
stage is guarded by a sentinel file, so reruns pick up where
the last invocation left off.

**Contents**: [Pipeline overview](#pipeline-overview) ·
[Installation](#installation) · [Quickstart](#quickstart) ·
[Invocation modes](#invocation-modes) · [Stage details](#stage-details) ·
[Output tree](#output-tree) · [Configuration](#configuration) ·
[HPC / Slurm](#hpc-usage-slurm) ·
[Re-runs & resume](#repro--re-run-behaviour) ·
[Development & testing](#development--testing) ·
[Citation](#citation)

**Deeper docs**: [`docs/install.md`](docs/install.md)
(troubleshooting + env drift), [`docs/methods.md`](docs/methods.md)
(algorithm details), [`docs/usage.md`](docs/usage.md)
(extended examples).

## Pipeline overview

Five stages. When `--he-slide` is a `.vsi` file, the HPC driver
runs a VSI → OME-TIFF conversion job first and chains hexenium
after it with `--dependency=afterok:`. Per-stage outputs
colocate under
`<output-root>/<sample>/<sample>_<run-id>/he_registration/` —
the canonical layout that xenium-preprocess-pipeline peers with:

```
    H&E slide                    hexenium pipeline                   per-run outputs
                                                                     (<output-root>/<sample>/
    ┌────────────────┐        ┌───────────────────┐                   <sample>_<run-id>/
    │ .vsi or        │───────▶│ he_preprocess     │──┐                he_registration/)
    │ .ome.tif       │        │ VSI → OME-TIFF    │  │
    └────────────────┘        │ (no-op for TIFF)  │  │
                              └───────────────────┘  │
                                                     ▼
                                    <sample>_he.ome.tif
                                    (location depends on writability;
                                     see he_preprocess stage details)
                                                     │
                                                     │
    Xenium bundle                                    │
    ┌────────────────┐        ┌───────────────────┐  │
    │ output-XETG..  │───────▶│ register          │◀─┘
    │ + morphology_  │        │ VALIS DAPI ↔ H&E  │
    │   focus/       │        │ (rigid + non-     │
    └────────────────┘        │  rigid + micro)   │
                              └───────────────────┘
                                        │
                                        ▼
                              register/<he_job_id>/data/_registrar.pickle
                                        │
                                        │
                              ┌───────────────────┐
                              │ warp              │
                              │ Xenium cell +     │
                              │ nucleus polygons  │
                              │ → H&E pixel space │
                              └───────────────────┘
                                        │
                                        ▼
                              warp/<he_job_id>/he_cell_seg.parquet
                              warp/<he_job_id>/he_nucleus_seg.parquet
                                        │
    xenium-ranger                       │
    ┌────────────────┐        ┌───────────────────┐
    │ xenium_ranger  │───────▶│ celltype          │
    │  .h5ad         │        │ read .obs[<col>]  │
    │ .obs[<col>] =  │        │ direct on xenium  │
    │  celltype label│        │ cells             │
    └────────────────┘        └───────────────────┘
                                        │
                                        ▼
                              celltyped/<he_job_id>/
                                 <sample>_cells_analysis.geojson
                                 <sample>_nuclei_analysis.geojson
                                 <sample>_cells_qupath.geojson
                                 <sample>_nuclei_qupath.geojson
                                 <sample>_celltyped_wholeslide.parquet
                                        │
                                        ▼
                              ┌───────────────────┐
                              │ viz               │
                              │ H&E overlay PNG   │
                              │ + summary.html    │
                              └───────────────────┘
                                        │
                                        ▼
                              viz/<he_job_id>/<sample>_overlay.png

                              (only on --set-default-on-success or
                               hexenium set-default-run:)
                              output/summary.html
                              output/{register,warp,celltyped,viz}   [symlinks]
```

### H&E input handling

Hexenium accepts two H&E input types.

**OME-TIFF / TIFF** (`.ome.tif`, `.ome.tiff`, `.tif`) — used
as-is. The `he_preprocess` stage is a no-op.

**Olympus `.vsi`** — the HPC driver
(`submit_he_registration.sh`) auto-converts to OME-TIFF via
BioFormats before running hexenium. Conversion runs as a
separate Slurm job; hexenium waits on
`--dependency=afterok:<conv_jobid>`. Reruns reuse the
converted file. See [HPC usage — VSI
inputs](#vsi-inputs-automatic-bioformats-conversion-unified)
for install requirements and options.

### Choosing / promoting the default registration

Every stage writes into its own `<he_job_id>` shard. Once a
run is promoted (see below), the
`output/{register,warp,celltyped,viz}` symlinks point at
that run so `output/register/` etc. show one canonical run
for browsing.

**Auto-promote on success.** Add `--set-default-on-success`
to `hexenium run`. When the run completes cleanly, the
symlinks are repointed at its fresh `<he_job_id>`.

**Compare N registrations, then promote manually.** Run
hexenium N times without `--set-default-on-success`, review
each `register/<he_job_id>/` shard, then promote the winner:

```bash
hexenium set-default-run \
    --output-root <output_root> \
    --sample-id   <sample> \
    --run-id      <upstream_run_id> \
    --register-run-id <winning_he_job_id>
```

`hexenium set-default-run --help` lists the per-stage flags
(`--register-run-id`, `--warp-run-id`, `--celltyped-run-id`,
`--viz-run-id`) — you can retarget any subset independently.

### Celltype input

The `celltype` stage assigns a cell-type label to every
warped Xenium cell. The default (**ranger-direct**) reads
labels from `--xenium-h5ad`'s `.obs["celltype"]`; the
legacy path (`--proseg-purified-h5ad`, passed **in addition
to** `--xenium-h5ad`) does spatial NN-mapping from a proseg
h5ad instead. Missing / empty column, or neither h5ad flag
set → `unlabeled` fallback.

See the [`celltype` stage docs](#celltype--assign-cell-type-labels-to-warped-polygons)
for source-mode details, `--celltype-col` overrides,
column-resolution precedence, and the full set of knobs.

## Installation

`hexenium` targets Python 3.10+ (validated on 3.11). It
depends on VALIS, HEST, and their Java/BioFormats bridge
(`jpype1`, `openjdk=11`).

Reproducing the working env requires a **two-step install**:

1. A conda-side solve for the system libraries and
   scientific-Python core.
2. A `--no-deps` pip layer for `hest` / `valis-wsi` /
   `valis-hest` and their transitive stack.

Both steps are captured in the order they should be run in
`environments/heRegistration.yml` and
`environments/heRegistration-requirements.txt`.

`--no-deps` on the pip step is **essential** and must live
on the CLI. See
[`docs/install.md`](docs/install.md#why-no-deps-is-mandatory)
for the resolver-conflict details.

### Recommended install

```bash
# 1. Clone the repo
git clone https://github.com/settylab/heRegistration.git
cd heRegistration

# 2. Create the conda env from the pinned spec.
#    `heRegistration` is just a readable convention — the sbatch
#    wrapper activates envs by resolved absolute prefix (step 8 in
#    docs/install.md), not by name.
#    scripts/create-env.sh wraps any `micromamba` env-creation call; if
#    you set MAMBA_ROOT_PREFIX for an isolated install it also keeps
#    the package cache isolated (see docs/install.md § 2).
scripts/create-env.sh env create -n heRegistration -f environments/heRegistration.yml
# or: conda env create -n heRegistration -f environments/heRegistration.yml

# 3. Activate
micromamba activate heRegistration
# or: conda activate heRegistration

# 4. Install the pinned pip layer with `--no-deps`.
#    Ships hest v1.2.0, valis-wsi 1.1.0, valis_hest 0.0.2, and
#    their supporting stack. `--no-deps` is required — see the
#    callout below.
pip install --no-deps -r environments/heRegistration-requirements.txt

# 5. Install hexenium itself. `--no-deps` here prevents pip
#    from re-resolving pyproject.toml's deps and clobbering
#    the pinned layer from step 4.
pip install --no-deps .

# 6. Sanity-check the install.
python -c "import hest, valis_hest, dask, openslide; print('OK')"
hexenium --version
hexenium run --help

# 7. (Optional — only if you have .vsi H&E input.)
#    Download the Glencoe release zips so the pipeline can
#    auto-convert VSI to OME-TIFF via
#    submit_he_registration.sh (tested versions:
#    bioformats2raw 0.12.1 + raw2ometiff 0.9.0):
mkdir -p $HOME/opt/bftools
cd $HOME/opt/bftools
wget https://github.com/glencoesoftware/bioformats2raw/releases/download/v0.12.1/bioformats2raw-0.12.1.zip
wget https://github.com/glencoesoftware/raw2ometiff/releases/download/v0.9.0/raw2ometiff-0.9.0.zip
unzip bioformats2raw-0.12.1.zip
unzip raw2ometiff-0.9.0.zip
```

If any verification line in step 6 errors, jump to
[Troubleshooting](#troubleshooting) below.

**8. REQUIRED before any sbatch submission** — steps 1-7 only cover
interactive use. Run this in **any** shell, even a brand new one (it
does not depend on step 7's shell state — see the warning below):

```bash
micromamba env list   # find heRegistration's absolute prefix
scripts/write-env-config.sh --env-prefix /path/to/micromamba/envs/heRegistration
scripts/env-preflight.sh   # verify it resolves end-to-end
```

If you did step 7 for VSI support, pass those two paths again
explicitly — **do not reuse `$BFTOOLS_ROOT`/`$LIBBLOSC_DIR` as bare
shell variables here**: they only exist in the shell session that ran
step 7, and if that step was skipped (no `.vsi` input — the common
case) or you're in a fresh terminal, both are unset. An unset,
unquoted `$VAR` vanishes from the command line entirely rather than
expanding to an empty string, which silently shifts every argument
after it — `write-env-config.sh` then reports a confusing "does not
exist" error for a path you never intended to pass. Use the literal
paths from step 7 instead:

```bash
scripts/write-env-config.sh \
    --env-prefix    /path/to/micromamba/envs/heRegistration \
    --bftools-root  $HOME/opt/bftools \
    --libblosc-dir  $HOME/micromamba/envs/blosc/lib
scripts/env-preflight.sh
```

If you don't already have a conda env with
`conda-forge::blosc` installed (needed for
`LIBBLOSC_DIR`), [`docs/install.md`](docs/install.md) has
a minimal `blosc`-only env recipe under
`### VSI-input prerequisites`.

Step 8 is not optional if you plan to submit anything via
`scripts/submit_he_registration.sh` — see
[`docs/install.md`](docs/install.md#8-configure-slurm-env-activation-required-before-any-sbatch-submission)
for why skipping it produces a submission that *looks* successful and
fails later, silently, inside the job.

### Advanced install topics

For depth on installation edge cases, see
[`docs/install.md`](docs/install.md):

- [Why `--no-deps` is mandatory](docs/install.md#why-no-deps-is-mandatory)
  — the pandas / pyvips / scikit-image over-cap in
  `valis-wsi 1.1`'s declared metadata + why the flag must
  live on the CLI.
- [Manual install (advanced)](docs/install.md#manual-install-bypass-the-env-file)
  — reproduce the working env without the yml.
- [Troubleshooting](docs/install.md#troubleshooting) —
  the seven most common install-time failures + fixes
  (`ResolutionImpossible`, `ModuleNotFoundError`, VIPS
  warnings, `openslide-bin` sdist, etc.).
- [Env drift & recovery](docs/install.md#env-drift--recovery)
  — `--force-reinstall` gotchas around the pinned
  numpy / xarray / fastcluster / opencv-contrib layer.
- [Configure Slurm env activation](docs/install.md#8-configure-slurm-env-activation-required-before-any-sbatch-submission)
  — the required post-install step (`write-env-config.sh` +
  `env-preflight.sh`) before any `sbatch` submission works.

## Quickstart

A bare `hexenium run …` runs the full pipeline under its defaults
(`--mode full_with_micro`, `--use-he-deconvolution true`,
`--viz-render-boundaries nucleus`). The most common shape is
**integrated-by-run-id** — a single command from an upstream
`xenium-preprocess` output tree:

```bash
micromamba activate heRegistration

hexenium run \
    --sample-id       SAMPLE1 \
    --run-id          demo_v1 \
    --output-root     /data/xenium_runs \
    --he-slide        /data/SAMPLE1/HE/SAMPLE1_he.ome.tif \
    --xenium-bundle   /data/SAMPLE1/xenium/output-XETG.../
```

For a **standalone** run (no upstream xenium-preprocess tree), point at
the celltype source explicitly — pass `--xenium-h5ad` for the
ranger-direct default:

```bash
hexenium run \
    --sample-id      SAMPLE1 \
    --he-slide       /data/SAMPLE1/HE/SAMPLE1_he.ome.tif \
    --xenium-bundle  /data/SAMPLE1/xenium/output-XETG.../ \
    --dapi-path      /data/SAMPLE1/xenium/output-XETG.../morphology_focus/ch0000_dapi.ome.tif \
    --xenium-h5ad    /data/SAMPLE1/spatial_adata/SAMPLE1_xenium_ranger.h5ad \
    --output-root    /data/hexenium_runs
```

For the legacy proseg-NN mapping instead, add
`--proseg-purified-h5ad /data/SAMPLE1/spatial_adata/SAMPLE1_proseg_purified.h5ad`
to the invocation above. See [`celltype` — assign cell-type
labels to warped polygons](#celltype--assign-cell-type-labels-to-warped-polygons)
for when to reach for it.

`--dapi-path` is only needed for multichannel Xenium bundles where the
DAPI lives at `morphology_focus/ch0000_dapi.ome.tif` instead of the
standard `morphology_focus/morphology_focus_0000.ome.tif`.

Common per-run overrides:

- `--viz-render-boundaries both` — draw cell polygons and nucleus outlines
  together.
- `--use-he-deconvolution false` — faint-hematoxylin samples where
  deconvolution hurts alignment.
- `--mode rigid_only` / `--mode rigid_nonrigid` — diagnostic sweeps.
- `--stages register warp celltype viz` — restrict which stages run
  (`he_preprocess` is a no-op for OME-TIFF input; skip explicitly when the
  input is a `.vsi` you don't want re-converted).
- `--force-rerun` — remove all sentinels and redo every requested stage.
- `--force-preprocess` — force VSI → OME-TIFF re-conversion only.

## Invocation modes

Three ways to point hexenium at a sample. Pick the one that matches how
you're driving the pipeline; the on-disk output layout under the top-level
run dir is identical in all three modes.

### 1. Standalone

Drive with `--sample-id` + `--he-slide` + `--xenium-bundle` +
`--output-root`. Outputs land at `<output_root>/<sample_id>/`.

If the `celltype` stage is in `--stages`, pass `--xenium-h5ad`
so hexenium can read labels directly off its `.obs[<col>]`
(the ranger-direct default). Add `--proseg-purified-h5ad`
only if you want to opt in to the legacy proseg-NN mapping
instead. Neither flag is enforced — omitting both doesn't error,
it just labels every cell `unlabeled` (with a log line), so a
typo'd or missing path here fails silently rather than loudly.

### 2. Integrated-by-run-id

Add `--run-id <upstream_run_id>` on top of the standalone
set. Hexenium derives the xenium h5ad from the upstream
`xenium-preprocess` layout at
`<output-root>/<sample>/<sample>_<run-id>/spatial_adata/<sample>_xenium_ranger.h5ad`.
All H&E outputs colocate under
`<output-root>/<sample>/<sample>_<run-id>/he_registration/`.

The celltype stage reads its labels directly from this
ranger h5ad's `.obs["celltype"]` — no proseg lookup by
default. Pass `--proseg-purified-h5ad <path>` explicitly to
switch to the legacy proseg-NN code path.

### 3. Integrated-by-h5ad

Pass `--xenium-h5ad <path>` directly. Sample identity is
read from `.uns['sample_id']` and `.uns['run_id']` on that
h5ad. Outputs colocate under `<xenium_run_dir>/he_registration/`.
Fails LOUD if `.uns` identity is missing or disagrees with a
passed `--sample-id`.

Celltype labels still come from this h5ad's
`.obs["celltype"]`. Same ranger-direct default as mode 2;
same legacy opt-in via explicit `--proseg-purified-h5ad`.

## Stage details

The subsections below cover user-facing behavior per stage.
For algorithm details and library-level implementation
(OpenSlide / tifffile writer settings, Dask worker + JVM
lifecycle, Shapely polygon cleanup, palette assignment),
see [`docs/methods.md`](docs/methods.md).

### `he_preprocess` — VSI → OME-TIFF (no-op for OME-TIFF inputs)

If `--he-slide` points to an Olympus SlideScanner `.vsi`,
this stage converts it to a pyramidal OME-TIFF matching
10x's Xenium Explorer specification. If the input is
already an OME-TIFF, the stage returns immediately and
downstream stages consume the input as-is.

The converted OME-TIFF is written **next to the source VSI**
as `<vsi_dir>/<vsi_stem>.ome.tif` when that directory is
writable. This lets a single conversion be reused across
every pipeline invocation for that sample. When the source
dir is read-only, the pipeline falls back to
`<output_root_he>/converted/<sample_id>_he.ome.tif`. Re-runs
auto-detect either location and skip re-conversion.
`--force-preprocess` forces re-conversion of the VSI without
touching downstream stages.

### `register` — H&E ↔ Xenium DAPI alignment (VALIS)

Registers the H&E to the DAPI morphology image using VALIS via
`valis_hest`. The default `--mode full_with_micro` composes
three transforms into one registrar:

1. A rigid initial solve.
2. A non-rigid B-spline solve at
   `max_processed_image_dim_px = 1500`.
3. A micro-registration refinement at
   `max_non_rigid_registration_dim_px = 10000`.

This matches stock HEST `register_dapi_he(micro_reg=True)`.
Two diagnostic modes are available for sweeps: `rigid_only`
(fastest, rigid only) and `rigid_nonrigid` (rigid + non-rigid,
no micro).

Macenko-style H&E deconvolution runs before feature detection
by default. On faint-hematoxylin samples it can degrade
alignment. Override with `--use-he-deconvolution false` when a
run shows that failure mode.

Reflection checking is on by default so mirrored slides are
caught cheaply. To sidestep a hardcoded
`he_key='aligned_fullres_HE'` in `valis_hest`, hexenium
symlinks the user's H&E to that canonical name inside the
per-sample workdir before invoking VALIS.

**Writes:**
`register/<he_job_id>/data/_registrar.pickle` — the composed rigid ×
non-rigid × micro transforms — plus VALIS's per-stage image dumps under
sibling `rigid_registration/`, `non_rigid_registration/`,
`micro_registration/`, etc., and a cross-run
`register/manifest.yaml` symlink pointing at the
current-default registrar.

### `warp` — apply the registrar to Xenium objects

Applies the VALIS registrar to the Xenium
`cell_boundaries.parquet` and `nucleus_boundaries.parquet`
(and, optionally, `transcripts.parquet`) via HEST's
`warp_and_save_xenium_objects`.

Outputs are WKB-encoded GeoPandas parquets in H&E pixel
space. Each row's Xenium `cell_id` is preserved for
downstream joins.

**Writes:** `warp/<he_job_id>/he_cell_seg.{parquet,geojson}` and
`warp/<he_job_id>/he_nucleus_seg.{parquet,geojson}` (and
`he_transcripts.parquet` when `--include-transcripts` is set). Parquet
for programmatic downstream use, GeoJSON for viewers (QuPath,
GeoJSON.io, `napari-geojson`, …).

### `celltype` — assign cell-type labels to warped polygons

Assigns a cell-type label to every warped Xenium cell. Viz
uses these labels to colour the overlay; downstream analysis
uses them to filter by class. Two source modes are available.
The default now reads labels DIRECTLY off the query h5ad
(post `rctd-split celltype_writeback`), so no NN mapping is
needed in the normal xenium-preprocess → hexenium flow.

**Default: ranger-direct.** With just `--xenium-h5ad` (or
`--run-id` auto-deriving it), hexenium reads
`xenium_ranger.h5ad`'s `.obs[<celltype_col>]` for every cell.
Default column: `celltype` (matches
`packages/rctd-split/config/default.yaml`'s
`celltype_writeback.celltype_col`). Pass `--celltype-col <name>` to
override for a custom-built ranger h5ad. `--celltype-col auto`
falls through the historic precedence
`celltype > first_type > primary_cell_type > celltype_updated`.

No NN mapping in this mode — the label comes straight off
the query h5ad. Hexenium's overlay reads the same
`.obs["celltype"]` column that the xenium-preprocess
summary reports read, so labels match unless the column
has been rewritten between the two invocations.

**Legacy: proseg-NN.** Pass `--proseg-purified-h5ad <path>`
explicitly (in addition to `--xenium-h5ad`) to switch to the
original NN-mapping path. For every xenium cell centroid,
hexenium looks up the nearest proseg-purified cell and inherits
its celltype label. This mode is kept for two cases:

1. Standalone-mode users who don't run xenium-preprocess's
   `rctd-split celltype_writeback`.
2. Users who want to override the ranger labels with a fresh
   NN fit against a custom proseg reference.

**Missing / all-NaN column → `unlabeled`.** When the requested
column doesn't exist, or exists but is entirely NaN / empty, the
stage logs a WARN and labels every row `unlabeled`. Viz renders
`unlabeled` as grey `#888888`. Per-row NaN / empty values also
fall to `unlabeled` (only the affected rows). No crash.

Auto-detection knobs:

- **celltype column** (`--celltype-col`, default `celltype`) — set
  to `auto` for precedence walk; set to a literal column name to
  force a choice; unmatched → `unlabeled` fallback.
- **Xenium id column** (`--id-col`, default `auto`) — every candidate is
  shape-checked against the Xenium UUID regex (`^[a-z]{8}-\d+$`) before
  use, so it can't silently pick Proseg's int64 index masquerading as
  `cell_id`. Special value `__index__` reads from `.obs.index`.
- **spatial coords** (proseg-NN mode only) — checks `.obsm['spatial']`
  first, then falls back through `.obs['centroid_x'/'centroid_y']`,
  `.obs['x_centroid'/'y_centroid']`, `.obs['x'/'y']`. Override with
  `celltype.{proseg,xenium}_{x,y}_col` in the config YAML.

After the label join, polygon cleanup drops empty /
non-Polygon geometries and geometries smaller than
`area_threshold_px` (default 20 sq px). Nuclei without a
label optionally inherit their sibling cell's label
(`nuclei_inherit_classification: true`).

**Writes** under `celltyped/<he_job_id>/`:

- `<sample>_cells_analysis.geojson` — cells, raw-float precision.
- `<sample>_nuclei_analysis.geojson` — nuclei, raw-float precision.
- `<sample>_cells_qupath.geojson` — cells, rounded (default 2 decimals)
  + polygon-safety cleanup so QuPath's drag-and-drop importer accepts
  them without invalid-polygon errors.
- `<sample>_nuclei_qupath.geojson` — nuclei, ditto.
- `<sample>_celltyped_wholeslide.parquet` — merged, analysis-precision
  combined table feeding the viz stage.
- `<sample>_celltype_annotation.parquet` — per-cell label
  table (columns: `xenium_cell_id`, `group`, `nn_distance`),
  written by both the ranger-direct and legacy proseg-NN
  branches. In ranger-direct mode `nn_distance` is `NaN`
  (no NN fit occurred); in proseg-NN mode it carries the
  L2 distance to the winning proseg centroid.

### `viz` — publication-shape overlay PNG

Draws warped boundaries on a downsampled H&E thumbnail.
The default (`viz.render_boundaries = nucleus`) draws only
nucleus outlines — the cleanest read of registration
quality against H&E's hematoxylin-stained nuclei. Two
alternatives are exposed via `--viz-render-boundaries`:
`cell` draws cell polygons (semi-transparent fill) only,
and `both` draws nucleus outlines on top of the cell
polygons.

Cell-type labels drive a qualitative palette (`tab20` by
default). Override colors in `viz.classification_palette`.

**Writes:** `viz/<he_job_id>/<sample>_overlay.png` at
`viz.dpi` (default 200 DPI), with a legend of the labels
present in the data.

## Output tree

Every mode writes into a **top-level H&E-reg output dir** — call it
`<output_root_he>` — with the same shape. Integrated modes place that dir
under an upstream xenium run dir; standalone mode places it under
`<output_root>/<sample_id>/`.

```
<output_root_he>/                              ← integrated: <run_dir>/he_registration/
│                                                standalone: <output_root>/<sample_id>/
├── converted/
│   └── <sample>_he.ome.tif                    ← stage 0 fallback (only when
│                                                <vsi_dir> was read-only; the
│                                                default write is next to the VSI)
├── register/
│   ├── manifest.yaml                          ← symlink to current-default registrar
│   └── <he_job_id>/
│       ├── data/_registrar.pickle             ← VALIS composed transform (sentinel)
│       ├── rigid_registration/
│       ├── non_rigid_registration/
│       ├── micro_registration/
│       ├── deformation_fields/
│       ├── masks/
│       ├── overlaps/
│       ├── processed/
│       └── _workdir/                          ← canonical-name H&E symlink
├── warp/
│   └── <he_job_id>/
│       ├── he_cell_seg.parquet    (+ .geojson)
│       ├── he_nucleus_seg.parquet (+ .geojson)
│       └── he_transcripts.parquet             ← only with --include-transcripts
├── celltyped/
│   └── <he_job_id>/
│       ├── <sample>_cells_analysis.geojson    ← per-cell polygons, float precision
│       ├── <sample>_nuclei_analysis.geojson
│       ├── <sample>_cells_qupath.geojson      ← rounded + polygon-safe
│       ├── <sample>_nuclei_qupath.geojson
│       ├── <sample>_celltyped_wholeslide.parquet  ← merged, feeds viz
│       └── <sample>_celltype_annotation.parquet   ← label table (both modes; NaN nn_distance in ranger-direct)
└── viz/
    └── <he_job_id>/
        └── <sample>_overlay.png               ← H&E + warped boundaries
```

**`<he_job_id>` precedence** (highest wins):

1. `--he-job-id` — explicit CLI override.
2. `$SLURM_JOB_ID` — set automatically under `sbatch` /
   `srun`.
3. `YYYYMMDDTHHMMSS` — local-time timestamp
   (interactive fallback, 1-second resolution).

**Logs** — per-run logs land at:

- **Integrated modes** →
  `<run_dir>/logs/logs_heRegistration/<sample>_<he_job_id>_<stages>/`
  (colocated with any workflow-driver step logs under `<run_dir>/logs/`).
- **Standalone** →
  `<output_root_he>/logs/<sample>_<he_job_id>_<stages>/`.

The `<stages>` suffix (e.g. `all`, `register_warp`) records which pipeline
stages ran, so re-runs at the same jobid don't clobber each other's logs.
Each per-run logs dir carries the Slurm `.out`/`.err` and a
`resolved_config.yaml` snapshot with every parameter that shaped the run
(defaults + user YAML + CLI overrides, merged).

## Configuration

Every knob is defined in `src/hexenium/_defaults/default.yaml`
with a comment explaining its effect. Override any subset with
a user YAML (`hexenium run --config user.yaml …`) or
spot-override on the CLI. CLI overrides win over user YAML.
User YAML wins over the default.

### Most-tuned knobs

The knobs most runs actually touch:

| Key | Default | Meaning |
| --- | --- | --- |
| `--mode` / `registration.mode` | `full_with_micro` | `rigid_only`, `rigid_nonrigid`, or `full_with_micro`. `rigid_only_micro` is invalid — see below. |
| `--use-he-deconvolution` / `registration.use_he_deconvolution` | `true` | Macenko-style H&E stain deconvolution before registration. Override with `false` on faint-hematoxylin samples. |
| `--dapi-path` / `dapi_path` | derived from `xenium_bundle` | Pass explicitly for multichannel bundles (`morphology_focus/ch0000_dapi.ome.tif`). |
| `--proseg-purified-h5ad` / `proseg_purified_h5ad` | (unset) | **Legacy proseg-NN mode opt-in.** Explicit, **and only in addition to `--xenium-h5ad`** → celltype uses the historic NN-mapping code path. Set alone (without `--xenium-h5ad`) it does nothing — falls through to the same `unlabeled` fallback as neither flag being set. Not needed in the standard xenium-preprocess → hexenium flow (default ranger-direct). |
| `--xenium-h5ad` / `xenium_h5ad` | integrated modes: auto; standalone: `(unset)` | Query xenium h5ad. Default source for celltype labels (read directly from `.obs[<celltype_col>]`). Also the identity source in integrated-by-h5ad mode. Not enforced in standalone mode — if `celltype` is in `--stages` and neither this nor `--proseg-purified-h5ad` is set, every cell is labeled `unlabeled` (logged, not an error). |
| `--celltype-col` / `celltype.celltype_col` | `celltype` | Column on the ranger h5ad (or proseg h5ad in legacy mode) that carries per-cell labels. `auto` walks the precedence `celltype > first_type > primary_cell_type > celltype_updated`. Missing / all-NaN → `unlabeled` fallback. |
| `--id-col` / `celltype.id_col` | `auto` | Force a specific xenium-side id column. `__index__` reads from `.obs.index`. |
| `--nn-k` / `celltype.nn_k` | `1` | Neighbours per query in the proseg → xenium NN mapping. |
| `celltype.nuclei_inherit_classification` | `true` | Unlabelled nuclei inherit their sibling cell's label. |
| `celltype.area_threshold_px` | `20` | Drop polygons smaller than this (in H&E pixels²) as segmentation noise. |
| `viz.thumbnail_max_dim` | `4096` | Downsample H&E to this max dimension for the overlay. |
| `viz.dpi` | `200` | Overlay PNG DPI. |
| `--viz-render-boundaries` / `viz.render_boundaries` | `nucleus` | Which boundaries to draw: `nucleus` (outlines only), `cell` (polygons only), or `both`. |
| `--include-transcripts` / `warp.include_transcripts` | `false` | Also warp `transcripts.parquet`. Adds hours. |

### Registration algorithm knobs

The `register` stage's numeric parameters live under
`registration:` and `parameters:` in
`src/hexenium/_defaults/default.yaml`. Every one is
CLI-overridable (dash-cased flag) and YAML-overridable (dot
path).

| Key | Default | Meaning |
| --- | --- | --- |
| `--check-for-reflections` / `registration.check_for_reflections` | `true` | Have VALIS test mirrored/flipped candidates during rigid init. Cheap; catches mounted-backwards slides. |
| `--create-masks` / `registration.create_masks` | `false` | Auto-mask non-tissue regions before registration. Off by default; auto-mask can over-mask small tissue islands. |
| `--align-to-reference` / `registration.align_to_reference` | `true` | Use H&E as the reference image (rather than DAPI). |
| `--max-image-dim-px` / `parameters.max_image_dim_px` | `1500` | Cap on the SAVED image pyramid (VALIS memory guard). Constraint: `max_image_dim_px >= max_processed_image_dim_px`; if lower, valis_hest silently auto-bumps. |
| `--max-processed-image-dim-px` / `parameters.max_processed_image_dim_px` | `1500` | Long-edge pixel cap for feature detection during rigid + first-pass non-rigid. Higher = more accurate + more memory. |
| `--max-non-rigid-registration-dim-px` / `parameters.max_non_rigid_registration_dim_px` | `10000` | Long-edge pixel cap for the non-rigid B-spline solve (`rigid_nonrigid` and `full_with_micro`). |
| `parameters.micro_rigid_registrar_cls` | `null` | Override VALIS's `MicroRigidRegistrar` class. Advanced; leave `null` for the stock class. |
| `parameters.micro_rigid_registrar_params` | `{}` | Extra kwargs forwarded to the micro-rigid registrar. |
| `--run-name` / `parameters.name` | `null` (auto) | Registrar run name recorded in the manifest for provenance (the on-disk folder is always `<he_job_id>`). |

**Note on `rigid_only_micro`:** removed — `register_micro` requires
`bk_dxdy` from an initial non-rigid pass; combining rigid-only init with
micro raises `TypeError: 'NoneType' object is not subscriptable` in
`valis_hest.registration.py:register_micro`. Use `rigid_only`,
`rigid_nonrigid`, or `full_with_micro`.

See `src/hexenium/_defaults/default.yaml` for the full list.
That file covers Dask worker settings, JVM memory, palette
overrides, and per-stage rounding for QuPath vs. analysis
outputs.

## Repro / re-run behaviour

- **Sentinel-file resume.** Each stage checks for its sentinel
  output(s) before starting. If present and `--force-rerun`
  isn't set, the stage is skipped. Its outputs are threaded
  through to the next stage.
- **Per-invocation sharding.** `register/`, `warp/`,
  `celltyped/`, `viz/` all shard under `<he_job_id>/`. A new
  run under a new Slurm job id never clobbers a prior run's
  outputs. `celltype` and `viz` reading an earlier stage's
  outputs default to the same `<he_job_id>`. Point them at a
  prior run with `--warp-run-id` / `--celltype-run-id`.
- **Config snapshot.** Every run writes
  `resolved_config.yaml` into its logs dir (see [Output
  tree](#output-tree)). This is the exact merged defaults +
  user YAML + CLI overrides that shaped the run. Integrated
  modes also merge a `he_registration:` key into the xenium
  run dir's shared `resolved_config.yaml` for cross-pipeline
  provenance.
- **VSI conversion is per-sample, not per-run.** The
  converted OME-TIFF lives next to the source VSI (or the
  pipeline-tree fallback). It is reused by every invocation
  for that sample. `--force-preprocess` triggers a
  re-conversion without touching downstream sentinels.

## HPC usage (Slurm)

The bundled wrapper `scripts/submit_he_registration.sh` is
designed to be called **directly** (no `sbatch` prefix). It
self-submits under sbatch and routes its own `.out`/`.err`
under the run folder's logs dir.

Under the hood it does two-phase `sbatch --hold` → mkdir
per-job-id dir → release. This meets Slurm's
log-dir-must-exist-at-job-start requirement.

`#SBATCH` header in the wrapper (tune to your cluster):

```
--partition=campus-new  --time=2-00:00:00  --mem=196G  --cpus-per-task=12
```

Call:

```bash
./scripts/submit_he_registration.sh \
    --sample-id     SAMPLE1 \
    --run-id        demo_v1 \
    --output-root   /data/xenium_runs \
    --he-slide      /data/SAMPLE1/HE/SAMPLE1_he.ome.tif \
    --xenium-bundle /data/SAMPLE1/xenium/output-XETG.../ \
    [--stages       register warp celltype viz] \
    [--force-rerun]
```

`OUTPUT_ROOT=/data/xenium_runs` in your env is a fallback for
`--output-root`.

**Env activation is by resolved absolute prefix, never by name.** The
wrapper does **not** source `~/.bashrc` and has no
micromamba→mamba→conda fallback chain — an earlier version did, and
that `$HOME`-relative activation chain is what caused real jobs to hang
with 0-byte output (see `scripts/lib/env_config.sh`'s header comment).
It always activates the exact prefix recorded in
`scripts/env.local.conf`, written once by
`scripts/write-env-config.sh` — see [`docs/install.md` step
8](docs/install.md#8-configure-slurm-env-activation-required-before-any-sbatch-submission).
To point the wrapper at a different env, re-run `scripts/write-env-config.sh
--env-prefix /path/to/other/env` — there is no per-invocation flag for
this anymore. `--env-name` / `ENV_NAME` are still parsed (so old
invocations don't hit an "unknown flag" error) but now only emit a WARN
and change nothing.

**A misconfigured environment does not block submission — it fails
later, silently, inside the job.** `submit_he_registration.sh`'s
launcher branch (`scripts/submit_he_registration.sh:101-333` — the part
that runs when you invoke the script directly) parses flags and calls
`sbatch`/`scontrol release`; it never sources
`scripts/lib/env_config.sh`. That only happens in the **under-sbatch**
branch (`scripts/submit_he_registration.sh:388-390`), once the job
actually starts running on a compute node. So if `scripts/env.local.conf`
is missing or stale, `./scripts/submit_he_registration.sh ...` still
prints a normal `Submitted batch job NNNN` and exits `0` — the failure
only shows up afterward in that job's own `slurm-<jobid>.err`. If
you're queueing a run and stepping away, run `scripts/env-preflight.sh`
first so a bad config fails in your terminal instead.

### Where Slurm logs land

The wrapper builds the logs dir under either
`<output_root>/<sample>/<sample>_<run_id>/logs/logs_heRegistration/`
(integrated — `submit_he_registration.sh:302`) or
`<output_root>/<sample>/logs/` (standalone —
`submit_he_registration.sh:305`), with a leaf directory
`<sample>_<jobid>_<stages>` (`:307`, `:324`) — so re-runs at the same
job id don't clobber each other's logs. Slurm's own `.out`/`.err` land
inside that leaf as `slurm-<jobid>.out` / `slurm-<jobid>.err`
(`:313-314`). Concretely, for an integrated run:

```
<output_root>/<sample>/<sample>_<run_id>/logs/logs_heRegistration/<sample>_<jobid>_<stages>/slurm-<jobid>.{out,err}
```

**Integrated-by-h5ad is the one mode this wrapper doesn't specially
route.** The launcher branch has no `--xenium-h5ad` handling at all,
and still hard-requires `--output-root` + `--sample-id` to submit
(`:128-137`) even though `hexenium run --xenium-h5ad ...` itself needs
neither — identity and output location come from the h5ad's `.uns` (see
[Invocation modes](#invocation-modes) above). In practice this mode is
more often driven with a direct `sbatch`/`sbatch --wrap` call instead of
this wrapper (see "Direct `sbatch --wrap`" below); in that case Slurm's
own default applies — `slurm-<jobid>.out`/`.err` land in the submitting
shell's `cwd`, wherever that is (commonly the xenium run directory
itself, if that's what an upstream pipeline `cd`s into before
submitting). If you do drive h5ad-mode through
`submit_he_registration.sh`, you still need `--output-root`/
`--sample-id`/`--run-id` to satisfy the launcher, and the logs land per
the integrated/standalone rule above — not automatically under the
h5ad's own xenium run directory.

**Direct `sbatch --wrap`** is also fine when you want to override the
resource envelope per-run:

```bash
sbatch --partition=campus-new --time=2-00:00:00 --mem=196G --cpus-per-task=12 \
    --job-name=hexenium-SAMPLE1 \
    --wrap "micromamba activate heRegistration && hexenium run \
        --sample-id     SAMPLE1 \
        --run-id        demo_v1 \
        --output-root   /data/xenium_runs \
        --he-slide      /data/SAMPLE1/HE/SAMPLE1_he.ome.tif \
        --xenium-bundle /data/SAMPLE1/xenium/output-XETG.../"
```

Slurm's own `slurm-<jobid>.out` in `cwd` captures everything;
`PYTHONUNBUFFERED=1` (which `hexenium` sets internally) keeps stdout
line-buffered.

**Resource sizing.** Under the default `--mode
full_with_micro`, a typical whole-slide run takes 2–2.5 h
wall-clock. Breakdown:

- `register` (rigid + non-rigid): 30–60 min.
- `register_micro`: 30 min – 2 h, depending on tissue size.
- Warp + celltype + viz: another 30–60 min.

The 2-day `--time` in the wrapper is the honest ceiling for
`full_with_micro` on large samples. A `--mode rigid_only`
diagnostic sweep completes in well under an hour.

196 GB memory is sized for the non-rigid + micro solve at
`max_non_rigid_registration_dim_px = 10000` plus the Dask
warp cluster. Drop to 128 GB only if you also drop that cap.

### VSI inputs: automatic BioFormats conversion (unified)

`submit_he_registration.sh` handles VSI inputs end-to-end
without any extra flags. Pass `--he-slide /path/to/foo.vsi`
and the launcher:

1. Detects the `.vsi` extension.
2. Submits a `bioformats2raw` + `raw2ometiff` conversion
   job to Slurm as its own record.
3. Rewrites `--he-slide` internally to the converted
   OME-TIFF at
   `<run_dir>/he_registration/converted/<sample>_he.ome.tif`
   (the HPC driver always writes to this pipeline-tree
   location, whether or not the source VSI's directory is
   writable; interactive `hexenium run` on a `.vsi`
   prefers next-to-source — see the `he_preprocess`
   stage docs).
4. Submits the hexenium job with
   `--dependency=afterok:<conv_jobid>` so it runs only when the
   conversion finishes cleanly.

Both jobs are independent slurm records. A hexenium rerun
(different `--mode`, `--warp-run-id`, etc.) reuses the
converted OME-TIFF without redoing the ~15-25-min BioFormats
step. Each job's `slurm-<jobid>.out/err` lives at a stable
log path.

**Idempotency**. The converted file is reused on subsequent
invocations unless you pass `--force-preprocess` (analogous to
`--force-rerun` for the hexenium stages themselves).

**Series auto-detect**. Olympus VSI files typically carry 4
series (label / overview / 40x brightfield / macro thumbnail).
The launcher auto-picks the largest-on-disk series after
`raw2ometiff --split` — a robust proxy for pixel count that
avoids the s2-is-always-the-40x assumption. Override with
`--vsi-series <N>` if the auto-pick ever gets it wrong for
your data.

**Prerequisites** — see Step 7 of [Recommended
install](#recommended-install) above for the Glencoe-zip
setup + `BFTOOLS_ROOT` / `LIBBLOSC_DIR` env vars. Full
deep-dive (Java 11 note, minimal `blosc`-only env
recipe, release-page links) in
[`docs/install.md`](docs/install.md#vsi-input-prerequisites-only-if---he-slide-is-a-vsi-file).

**Example** — a VSI-input invocation looks identical to any
other input, just with a `.vsi` path:

```bash
./scripts/submit_he_registration.sh \
    --sample-id     SAMPLE1 \
    --run-id        demo_v1 \
    --output-root   /data/xenium_runs \
    --he-slide      /data/SAMPLE1/HE/SAMPLE1.vsi \
    --xenium-bundle /data/SAMPLE1/xenium/output-XETG.../ \
    --dapi-path     /data/SAMPLE1/xenium/output-XETG.../morphology_focus/morphology_focus_0000.ome.tif
```

Output:
```
[submit] VSI input detected — submitting conversion job:
Submitted batch job 12344 (conversion)
Submitted batch job 12345 (logs: SAMPLE1_12345_all)
  waiting on conversion job 12344 (afterok)
```

On a rerun with the same `--run-id`, the launcher notices the
canonical OME-TIFF is already present. It prints
`[submit] skipping VSI conversion`. Only the hexenium job is
submitted — no `--dependency` wait, starts immediately.

### Standalone conversion scripts (advanced)

Two auxiliary launchers cover the rare cases the unified flow
above doesn't handle cleanly:

- `scripts/submit_vsi_to_ometiff.sh` — VSI → OME-TIFF only. Use
  to prime a shared cache of converted files outside a specific
  `--run-id`, or to debug conversion parameters
  (`--vsi-series`, `--max-workers`) in isolation.
  Same install prerequisites as above.
- `scripts/submit_hexenium_from_ometiff.sh` — hexenium against
  an already-converted OME-TIFF at a non-canonical path (e.g.
  a Trident / QuPath / `bfconvert` output). Thin dispatcher to
  `submit_he_registration.sh` — forwards every hexenium flag
  untouched.

Both are called the same way, with `--help` support and
`Submitted batch job <jobid>` output routed to a `logs/`
subdirectory next to their outputs.

## Development & testing

The test suite covers:

- CLI parsing and config validation.
- Layout resolution for the three invocation modes.
- Celltyping (auto-detect, NN mapping, polygon cleanup).
- Viz `render_boundaries` variants.
- Stage imports.
- An end-to-end pipeline smoke test.

```bash
pip install -e '.[test]'
pytest tests/
```

Tests are ordinary `pytest` — no cluster or GPU required. `pytest -k
layout` and friends work for targeted runs. `conftest.py` in `tests/`
sets up shared fixtures.

## Citation

If you use hexenium in a publication, cite via the metadata in
[`CITATION.cff`](CITATION.cff) — GitHub renders that file as a
copy-paste-ready BibTeX / APA / CFF block on the repo landing page.
Methods-section template (fill in the version tag and Zenodo DOI):

> Registration of the H&E whole-slide image to the matched Xenium DAPI
> morphology channel was performed with hexenium v0.2.0 (DOI:
> <TBD>), a Python wrapper around VALIS (Gatenbee et
> al., 2023) via the `valis_hest` adapter and HEST (Jaume et al., 2024)
> that composes rigid, non-rigid, and micro-registration transforms into
> a single registrar. Xenium cell and nucleus segmentations were warped
> into H&E pixel space via HEST's `warp_and_save_xenium_objects` and
> annotated with per-cell type labels sourced from
> `xenium_ranger.h5ad`'s `.obs["celltype"]` (as written by an upstream
> `rctd-split` `celltype_writeback` step). Annotated boundaries were
> then rendered as per-slide overlays.

Please also cite VALIS and HEST directly per their upstream requests.

## License

MIT — see [`LICENSE`](LICENSE).

## Acknowledgements

- **VALIS** — Gatenbee, C. D. et al. *VALIS: Virtual Alignment of
  pathoLogy Image Series*. Nature Communications 14, 4502 (2023).
- **HEST** — Jaume, G. et al. *HEST-1k: A Dataset for Spatial
  Transcriptomics and Histology Image Analysis*. NeurIPS Datasets and
  Benchmarks (2024). Code: `mahmoodlab/HEST`.
- **Xenium In Situ platform** — 10x Genomics. The Xenium bundle layout
  (`morphology_focus/`, `cell_boundaries.parquet`,
  `nucleus_boundaries.parquet`, `transcripts.parquet`) is a 10x-Genomics
  output specification.

hexenium was adapted from an internal H&E-Xenium registration workflow.
