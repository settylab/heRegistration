# hexenium

**Single-sample H&E ↔ Xenium DAPI registration, warp, cell-type propagation, and overlay visualisation.**

> **v0.2.0** — celltype now sources labels directly from an upstream
> `proseg_purified.h5ad` via an inlined nearest-neighbour mapping (no more
> `--celltype-csv`). Three invocation modes (standalone /
> integrated-by-run-id / integrated-by-h5ad); per-stage outputs are sharded
> under `register/<he_job_id>/`, `warp/<he_job_id>/`,
> `celltyped/<he_job_id>/`, `viz/<he_job_id>/`. See
> [CHANGELOG.md](CHANGELOG.md) for the full delta and migration notes from
> v0.1.x.

At a glance:

- **Register** a hematoxylin and eosin (H&E) whole-slide image against the
  DAPI morphology channel of a matched Xenium in-situ transcriptomics run
  (VALIS via `valis_hest`), composing rigid, non-rigid, and
  micro-registration transforms into a single registrar.
- **Warp** Xenium cell + nucleus polygon boundaries from Xenium pixel space
  into H&E pixel space (HEST `warp_and_save_xenium_objects`).
- **Propagate** cell-type labels from an upstream `proseg_purified.h5ad`
  onto every Xenium cell via a spatial nearest-neighbour lookup on
  centroids; nuclei can inherit their sibling cell's label.
- **Visualise** with a publication-quality overlay: nucleus outlines on a
  downsampled H&E thumbnail, coloured by cell-type label. Add cell polygons
  with `--viz-render-boundaries both`.
- Every stage is idempotent — sentinel-file resume — and every run
  snapshots its resolved configuration to disk.

## Package name

The Python package + console script is **`hexenium`** (`pip install .`
under the tested **`heRegistration`** conda env). The env name is the
default `ENV_NAME` for the sbatch wrapper.

Table of contents: [What it does](#what-it-does) · [Installation](#installation) · [Invocation modes](#invocation-modes) · [Quickstart](#quickstart) · [Output tree](#output-tree) · [Configuration](#configuration) · [HPC (Slurm)](#hpc-usage-slurm) · [Development & testing](#development--testing) · [Citation](#citation) · [License](#license) · [Acknowledgements](#acknowledgements)

## What it does

Five stages, run in order (`he_preprocess → register → warp → celltype →
viz`). Restrict a run to a subset with `--stages`; each stage is guarded
by a sentinel file, so re-running with the same identity picks up where
the last run left off. Nuke sentinels with `--force-rerun`.

### `he_preprocess` — VSI → OME-TIFF (no-op for OME-TIFF inputs)

If `--he-path` points to an Olympus SlideScanner `.vsi`, this stage
converts it to a pyramidal OME-TIFF that matches 10x's Xenium Explorer
image-conversion specification: 1024×1024 tiles, lossless JPEG 2000 (or
ZLIB), 7-level pyramid at scale 2. `OpenSlide` reads the source;
`tifffile.TiffWriter` writes the pyramid; physical pixel size is
propagated from OpenSlide's `mpp-{x,y}` into the OME-TIFF's
`PhysicalSize{X,Y}` metadata. If the input is already an OME-TIFF, the
stage returns immediately and downstream stages consume the input as-is.

The converted OME-TIFF is written **next to the source VSI** as
`<vsi_dir>/<vsi_stem>.ome.tif` when that directory is writable, so a
single conversion is reused across every pipeline invocation for that
sample. When the source dir is read-only, the pipeline falls back to
`<output_root_he>/converted/<sample_id>_he.ome.tif`. Re-runs auto-detect
either location and skip re-conversion; `--force-preprocess` forces
re-conversion of the VSI without touching downstream stages.

### `register` — H&E ↔ Xenium DAPI alignment (VALIS)

Registers the H&E to the DAPI morphology image using VALIS via
`valis_hest`. The default `--mode full_with_micro` composes three
transforms into one registrar: a rigid initial solve, a non-rigid
B-spline solve at `max_processed_image_dim_px = 1500`, and a
micro-registration refinement at `max_non_rigid_registration_dim_px =
10000` — matching stock HEST `register_dapi_he(micro_reg=True)`. Two
diagnostic modes are available for sweeps: `rigid_only` (fastest, rigid
only) and `rigid_nonrigid` (rigid + non-rigid, no micro). Macenko-style
H&E deconvolution runs before feature detection by default; on
faint-hematoxylin samples it can degrade alignment, so override with
`--use-he-deconvolution false` when a run shows that failure mode.
Reflection checking is on by default so mirrored slides are caught
cheaply. To sidestep a hardcoded `he_key='aligned_fullres_HE'` in
`valis_hest`, hexenium symlinks the user's H&E to that canonical name
inside the per-sample workdir before invoking VALIS.

**Writes:**
`register/<he_job_id>/data/_registrar.pickle` — the composed rigid ×
non-rigid × micro transforms — plus VALIS's per-stage image dumps under
sibling `rigid_registration/`, `non_rigid_registration/`,
`micro_registration/`, etc., and a cross-run
`register/manifest.yaml` pointing at the latest registrar.

### `warp` — apply the registrar to Xenium objects

Applies the VALIS registrar to the Xenium `cell_boundaries.parquet` and
`nucleus_boundaries.parquet` (and, optionally, `transcripts.parquet`) via
HEST's `warp_and_save_xenium_objects`. Warping runs on a Dask
`LocalCluster` in which every worker initialises the BioFormats JVM
exactly once via a `WorkerPlugin`; JPype enforces a single JVM lifecycle
per Python process, so the JVM is deliberately not torn down after
registration. Outputs are WKB-encoded GeoPandas parquets in H&E pixel
space; each row's Xenium `cell_id` is preserved for downstream joins.

**Writes:** `warp/<he_job_id>/he_cell_seg.{parquet,geojson}` and
`warp/<he_job_id>/he_nucleus_seg.{parquet,geojson}` (and
`he_transcripts.parquet` when `--include-transcripts` is set). Parquet
for programmatic downstream use, GeoJSON for viewers (QuPath,
GeoJSON.io, `napari-geojson`, …).

### `celltype` — proseg → xenium NN mapping onto warped polygons

Assigns a cell-type label to every warped Xenium cell by nearest-neighbour
lookup on centroids. Source of labels is an upstream
`proseg_purified.h5ad` (its `.obs` carries the cell-type column and
per-cell centroids in Xenium µm); query is the matched xenium h5ad
(its `.obs` carries Xenium UUID `cell_id`s and centroids). No CSV, no
ID join — every Xenium cell is labelled by NN on the proseg side.

Auto-detection is deliberate:

- **celltype column** on proseg (`--celltype-col`, default `auto`) —
  precedence `celltype` > `first_type` > `primary_cell_type` >
  `celltype_updated`. Pass a literal column name to force a choice.
- **Xenium id column** (`--id-col`, default `auto`) — every candidate is
  shape-checked against the Xenium UUID regex (`^[a-z]{8}-\d+$`) before
  use, so it can't silently pick Proseg's int64 index masquerading as
  `cell_id`. Special value `__index__` reads from `.obs.index`.
- **spatial coords** — checks `.obsm['spatial']` first, then falls back
  through `.obs['centroid_x'/'centroid_y']`,
  `.obs['x_centroid'/'y_centroid']`, `.obs['x'/'y']`. Override with
  `celltype.{proseg,xenium}_{x,y}_col` in the config YAML.

After the label is joined onto the warped boundaries (via a
cumcount-augmented merge that survives Dask-emitted duplicate
`xenium_cell_id` rows), polygon cleanup runs `shapely.make_valid →
largest connected piece`, drops empty / non-Polygon /
sub-`area_threshold_px` (default 20 sq px) geometries, and sanity-checks
for finite coords, ≥3 unique vertices, and `is_valid`. Nuclei without a
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
- `<sample>_celltype_annotation.parquet` — the NN-mapping annotation
  table (xenium cell_id ↔ proseg cell_id ↔ label ↔ distance), kept as a
  debug sidecar.

### `viz` — publication-shape overlay PNG

Draws warped boundaries on a downsampled H&E thumbnail. By default
(`viz.render_boundaries = nucleus`) only nucleus outlines are drawn — the
cleanest read of registration quality against H&E's hematoxylin-stained
nuclei. Two alternatives are exposed via `--viz-render-boundaries`:
`cell` draws cell polygons (semi-transparent fill) only, and `both`
draws nucleus outlines on top of the cell polygons. Thumbnails are read
via OpenSlide with a `tifffile` fallback. Cell-type labels drive a
qualitative palette (`tab20` by default) with user-configurable
overrides in `viz.classification_palette`; any label not in the
override map is auto-assigned a colormap slot in first-appearance order
(deterministic across re-runs). The renderer is dask-friendly: polygon
batches are drawn in chunks so a whole-slide overlay does not need to
fit its geometries in memory as a single frame.

**Writes:** `viz/<he_job_id>/<sample>_overlay.png` at `viz.dpi` (default
200 DPI), with a legend of the labels present in the data.

## Installation

`hexenium` targets Python 3.10+ (validated on 3.11) and depends on
VALIS, HEST, and their Java/BioFormats bridge (`jpype1`, `openjdk=11`).
Reproducing the working env requires a **two-step install** — a
conda-side solve for the system libraries and scientific-Python core,
then a `--no-deps` pip layer for `hest` / `valis-wsi` / `valis-hest`
and their transitive stack. Both steps are captured, in the order they
should be run, in `environments/heRegistration.yml` and
`environments/heRegistration-requirements.txt`.

`--no-deps` on the pip step is **load-bearing** and must live on the
CLI (see the [why](#why---no-deps-is-mandatory) box below).

### Recommended (validated on Fred Hutch Gizmo, 2026-08)

The steps below are the exact sequence Tracy validated end-to-end on
`rhino` with `micromamba 2.6.2`. `conda` works too — swap the tool
name in step 2 (and `micromamba activate` → `conda activate` in step 3).

```bash
# 1. Clone the repo
git clone https://github.com/settylab/heRegistration.git
cd heRegistration

# 2. Create the conda env from the pinned spec. The yml intentionally
#    does NOT declare a `name:` field, so you pick the env name with
#    `-n`. `heRegistration` is the name the sbatch wrapper defaults to
#    (`ENV_NAME`), so using it keeps everything downstream working
#    without further overrides.
micromamba env create -n heRegistration -f environments/heRegistration.yml
# or with conda:
# conda env create -n heRegistration -f environments/heRegistration.yml

# 3. Activate
micromamba activate heRegistration
# or: conda activate heRegistration

# 4. Install the pinned pip layer with `--no-deps`.
#    This installs hest (git, v1.2.0), valis-wsi 1.1.0, valis_hest 0.0.2,
#    and the ~160 supporting packages (torch, transformers, ultralytics,
#    spatialdata, opencv, …) at the exact versions Tracy's working env
#    was validated against. `--no-deps` is required — see the callout
#    below. `hest @ git+https://github.com/mahmoodlab/HEST.git@v1.2.0`
#    is already the first entry in the requirements file, so no
#    separate `pip install hest` step is needed.
pip install --no-deps -r environments/heRegistration-requirements.txt

# 5. Install the hexenium package itself.
#    Non-editable is recommended for end users (see note below); use
#    `pip install -e .` if you plan to hack on hexenium's source.
pip install .

# 6. Sanity-check the install.
python -c "import hest, valis_hest, valis_hest.registration, valis_hest.slide_io, dask, openslide; print('OK')"
hexenium --version
hexenium run --help
```

If any of the verification lines errors, jump to
[Troubleshooting](#troubleshooting) below.

#### Why `--no-deps` is mandatory

`valis-wsi 1.1`'s declared `Requires-Dist` metadata caps
`pandas<2.0`, `pyvips<3.0`, and `scikit-image<0.20`. HEST has similar
declared-vs-actual divergences. Tracy's working env — the one this
package is validated against — runs the *newer* pandas 2.3.x /
pyvips 3.1.x / scikit-image 0.19.x branch. VALIS's actual code paths
never touch the pandas-1 or pyvips-2 API, so the caps are
over-defensive relative to the code that hexenium actually calls, but
pip's resolver refuses to install VALIS against pandas 2.x without
`--no-deps`. Skipping `--no-deps` yields:

```
ResolutionImpossible: valis-wsi 1.1.0 depends on pandas<2.0.0
```

`--no-deps` cannot be inlined into the yml `pip:` block (micromamba
treats each entry as a package name → `ERROR: Invalid requirement:
--no-deps`) or the top of the requirements file (pip refuses →
`no such option: --no-deps`), so it lives on the CLI in step 4.

#### Why `pip install .` (not `-e`) for end users

An editable install exposes the source tree to `sys.path`, so a stray
import via a working-directory Python — or a sibling `hexenium/`
folder in `cwd` — can shadow the installed package and silently pull
in half-updated modules. Non-editable is safer for end users. If you
are actively hacking on hexenium, the "Development & testing"
section further down covers the editable install.

### Manual install (bypass the env file)

If you can't or don't want to use the env yml, this recipe reproduces
the same working env by hand. The two-step conda solve + `--no-deps`
pip layer structure is the same.

```bash
# Conda side (system libs + scientific-Python core)
micromamba create -n heRegistration -c conda-forge -c bioconda \
    python=3.11 'numpy<2' 'pandas<3' scipy pyarrow pyyaml \
    'scanpy>=1.10' dask-geopandas shapely proj pyproj ipykernel \
    importlib_metadata 'pycparser>=2.14' \
    'libvips>=8.15' 'pyvips>=3' imagemagick openslide openjdk=11
micromamba activate heRegistration

# Pip side (VALIS + HEST + their transitive stack). See
# environments/heRegistration-requirements.txt for the exact pins
# every package should be installed at.
pip install --no-deps valis-wsi==1.1.0 valis_hest==0.0.2
pip install --no-deps "hest @ git+https://github.com/mahmoodlab/HEST.git@v1.2.0" hestcore==1.0.4
# HEST + VALIS have ~160 pip-only transitive deps (torch, transformers,
# ultralytics, spatialdata, opencv-*, anndata 0.12 override, etc.).
# The full pinned list is `environments/heRegistration-requirements.txt`
# — the easiest way to install it is just to run step 4 of the
# recommended flow above from the repo checkout.
pip install .  # or `pip install -e .` for development
```

### Troubleshooting

The env captured in the two files above is versioned deliberately to
reproduce exactly. If the install fails, the failure mode is almost
always one of the following.

- **`ResolutionImpossible: valis-wsi 1.1.0 depends on pandas<2.0.0`**
  (or similar for `pyvips` / `scikit-image`) during step 4. You forgot
  `--no-deps` — see [Why `--no-deps` is mandatory](#why---no-deps-is-mandatory)
  above. Re-run with `pip install --no-deps -r
  environments/heRegistration-requirements.txt`.
- **`ModuleNotFoundError: No module named 'valis_hest'`** (or `hest`
  / `torch` / `transformers` / `ultralytics`) when you run `hexenium`.
  Step 4 was skipped or silently failed. Re-run it and re-check with
  `pip list | grep -Ei 'valis|hest|torch|transformers'`. Expected:
  `valis_hest 0.0.2`, `valis-wsi 1.1.0`, `hest 1.1.1` (installed from
  the `v1.2.0` git tag; the package's internal version is `1.1.1`),
  `hestcore 1.0.4`, `torch 2.6.0`, `transformers 5.1.0`,
  `ultralytics 8.4.14`.
- **`VIPS-WARNING: unable to load "vips-magick.so" ... libMagickCore-7.Q16HDRI.so.10:
  cannot open shared object file`** at first `import pyvips` /
  `libvips_init`. `imagemagick` is present in the yml exactly to
  silence this — the warning should not fire on a freshly-created
  env. If it does, your conda solve is stale or pinned to an older
  yml revision: `micromamba clean -a` then re-create. The warning
  itself is benign (the pipeline never uses that loader), so
  `VIPS_WARNING=off` is a safe temporary silencer.
- **Micromamba dependency cascade** during step 2 (`libvips` /
  `gdk-pixbuf` / `librsvg` / `scanpy` "no viable options"). Almost
  always a stale local cache or a stale `~/.condarc`. Fix in order:
  1. `micromamba clean -a` (nukes cached repodata + packages).
  2. If `~/.condarc` is a stale NFS handle (`[Errno 116] Stale file
     handle`), refresh it: `rm ~/.condarc && touch ~/.condarc` (or
     restore your original file).
  3. Re-run `micromamba env create -n heRegistration -f
     environments/heRegistration.yml` on a networked node.
- **`openslide-bin` sdist build fails** at step 4 with "Install
  OpenSlide from source". This is what the `@v1.2.0` pin on HEST
  guards against — HEST HEAD depends on TRIDENT which pulls
  `openslide-bin`, and its wheels don't cover older glibc. If you're
  hitting this, verify line 1 of `environments/heRegistration-requirements.txt`
  is `hest @ git+https://github.com/mahmoodlab/HEST.git@v1.2.0` and
  not an unpinned `hest @ git+…HEST.git`.
- **`hexenium` command not found** after step 5 succeeds. The wrong
  env is active — check `which hexenium` and re-run `micromamba
  activate heRegistration`.
- **Sbatch job fails at env activation** with "environment
  `heRegistration` not found". You created the env under a different
  name in step 2. Either recreate as `heRegistration`, or pass the
  name you used to the sbatch wrapper: `./scripts/submit_he_registration.sh
  --env-name <your-env-name> …` (or export
  `ENV_NAME=<your-env-name>` in your shell).

If none of these match, `pip install --no-deps -r
environments/heRegistration-requirements.txt -v` prints per-package
progress and surfaces which entry pip is choking on — most useful
when a wheel has been yanked from PyPI or a git ref has moved.

## Invocation modes

Three ways to point hexenium at a sample. Pick the one that matches how
you're driving the pipeline; the on-disk output layout under the top-level
run dir is identical in all three modes.

### 1. Standalone

Drive with `--sample-id` + `--he-path` + `--xenium-bundle` +
`--output-root`. Outputs land at `<output_root>/<sample_id>/`. Pass
`--proseg-purified-h5ad` (and, if you have it, `--xenium-h5ad`) explicitly
if the `celltype` stage is in `--stages`.

### 2. Integrated-by-run-id

Add `--run-id <upstream_run_id>` on top of the standalone set. Derives the
xenium h5ad from the upstream `xenium-preprocess` layout at
`<output-root>/<sample>/<sample>_<run-id>/spatial_adata/<sample>_xenium_ranger.h5ad`
and colocates all H&E outputs under
`<output-root>/<sample>/<sample>_<run-id>/he_registration/`.
`proseg_purified.h5ad` is auto-derived from
`<xenium_run_dir>/spatial_adata/<sample>_proseg_purified.h5ad` unless
`--proseg-purified-h5ad` overrides.

### 3. Integrated-by-h5ad

Pass `--xenium-h5ad <path>` directly. Sample identity is read from
`.uns['sample_id']` and `.uns['run_id']` on that h5ad; outputs colocate
under `<xenium_run_dir>/he_registration/`. Fails LOUD if `.uns` identity
is missing or disagrees with a passed `--sample-id`. Same
proseg auto-derivation as mode 2.

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
the celltype sources explicitly:

```bash
hexenium run \
    --sample-id             SAMPLE1 \
    --he-slide              /data/SAMPLE1/HE/SAMPLE1_he.ome.tif \
    --xenium-bundle         /data/SAMPLE1/xenium/output-XETG.../ \
    --dapi-path             /data/SAMPLE1/xenium/output-XETG.../morphology_focus/ch0000_dapi.ome.tif \
    --xenium-h5ad           /data/SAMPLE1/spatial_adata/SAMPLE1_xenium_ranger.h5ad \
    --proseg-purified-h5ad  /data/SAMPLE1/spatial_adata/SAMPLE1_proseg_purified.h5ad \
    --output-root           /data/hexenium_runs
```

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
- `--force-rerun` — nuke all sentinels and redo every requested stage.
- `--force-preprocess` — force VSI → OME-TIFF re-conversion only.

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
│   ├── manifest.yaml                          ← cross-run provenance (latest registrar)
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
│       └── <sample>_celltype_annotation.parquet   ← NN debug sidecar
└── viz/
    └── <he_job_id>/
        └── <sample>_overlay.png               ← H&E + warped boundaries
```

**`<he_job_id>` precedence** (three tiers): `--he-job-id` > `$SLURM_JOB_ID`
> `YYYYMMDDTHHMMSS` timestamp. Under sbatch the Slurm job id is used
automatically; interactive runs fall back to a UTC timestamp so
back-to-back runs never share a folder.

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

Every knob is defined in `config/default.yaml` with a comment explaining
its effect. Override any subset with a user YAML (`hexenium run --config
user.yaml …`) or spot-override on the CLI. CLI overrides win over user
YAML which wins over the default.

### Most-tuned knobs

The knobs most runs actually touch:

| Key | Default | Meaning |
| --- | --- | --- |
| `--mode` / `registration.mode` | `full_with_micro` | `rigid_only`, `rigid_nonrigid`, or `full_with_micro`. `rigid_only_micro` is invalid — see below. |
| `--use-he-deconvolution` / `registration.use_he_deconvolution` | `true` | Macenko-style H&E stain deconvolution before registration. Override with `false` on faint-hematoxylin samples. |
| `--dapi-path` / `dapi_path` | derived from `xenium_bundle` | Pass explicitly for multichannel bundles (`morphology_focus/ch0000_dapi.ome.tif`). |
| `--proseg-purified-h5ad` / `proseg_purified_h5ad` | integrated modes: auto-derived; standalone: required for `celltype` stage | Upstream proseg h5ad whose `.obs` carries celltype labels + centroids. |
| `--xenium-h5ad` / `xenium_h5ad` | integrated modes: auto; standalone: required for `celltype` | Query xenium h5ad. Also the identity source in integrated-by-h5ad mode. |
| `--celltype-col` / `celltype.celltype_col` | `auto` | Force a specific celltype column on the proseg side. |
| `--id-col` / `celltype.id_col` | `auto` | Force a specific xenium-side id column. `__index__` reads from `.obs.index`. |
| `--nn-k` / `celltype.nn_k` | `1` | Neighbours per query in the proseg → xenium NN mapping. |
| `celltype.nuclei_inherit_classification` | `true` | Unlabelled nuclei inherit their sibling cell's label. |
| `celltype.area_threshold_px` | `20` | Drop polygons smaller than this (in H&E pixels²) as segmentation noise. |
| `viz.thumbnail_max_dim` | `4096` | Downsample H&E to this max dimension for the overlay. |
| `viz.dpi` | `200` | Overlay PNG DPI. |
| `--viz-render-boundaries` / `viz.render_boundaries` | `nucleus` | Which boundaries to draw: `nucleus` (outlines only), `cell` (polygons only), or `both`. |
| `--include-transcripts` / `warp.include_transcripts` | `false` | Also warp `transcripts.parquet`. Adds hours. |

### Registration algorithm knobs

The `register` stage's numeric parameters live under `registration:` and
`parameters:` in `config/default.yaml`. Every one is CLI-overridable
(dash-cased flag) and YAML-overridable (dot path).

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

See `config/default.yaml` for the full list (Dask worker settings, JVM
memory, palette overrides, per-stage rounding for QuPath vs. analysis
outputs).

## Repro / re-run behaviour

- **Sentinel-file resume.** Each stage checks for its sentinel output(s)
  before starting; if present and `--force-rerun` isn't set, the stage
  is skipped and its outputs are threaded through to the next stage.
- **Per-invocation sharding.** `register/`, `warp/`, `celltyped/`, `viz/`
  all shard under `<he_job_id>/`, so a new run under a new Slurm job id
  never clobbers a prior run's outputs. `celltype` and `viz` reading an
  earlier stage's outputs default to the same `<he_job_id>`; point them at
  a prior run with `--warp-run-id` / `--celltype-run-id`.
- **Config snapshot.** Every run writes `resolved_config.yaml` into its
  logs dir (see [Output tree](#output-tree)) — the exact merged defaults
  + user YAML + CLI overrides that shaped the run. Integrated modes also
  merge a `he_registration:` key into the xenium run dir's shared
  `resolved_config.yaml` for cross-pipeline provenance.
- **VSI conversion is per-sample, not per-run.** The converted OME-TIFF
  lives next to the source VSI (or the pipeline-tree fallback) and is
  reused by every invocation for that sample. `--force-preprocess`
  triggers a re-conversion without touching downstream sentinels.

## HPC usage (Slurm)

The bundled wrapper `scripts/submit_he_registration.sh` is designed to be
called **directly** (no `sbatch` prefix); it self-submits under sbatch and
routes its own `.out`/`.err` under the run folder's logs dir. Under the
hood it does two-phase `sbatch --hold` → mkdir per-job-id dir → release,
so Slurm's log-dir-must-exist-at-job-start requirement is met.

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
    [--env-name     my-custom-env] \
    [--force-rerun]
```

`OUTPUT_ROOT=/data/xenium_runs` in your env is a fallback for
`--output-root`. The conda env name defaults to `heRegistration` and can
be overridden two ways — `--env-name my-custom-env` on the CLI (wins) or
`ENV_NAME=my-custom-env` in the environment. Example:

```bash
./scripts/submit_he_registration.sh \
    --sample-id     SAMPLE1 \
    --run-id        demo_v1 \
    --output-root   /data/xenium_runs \
    --he-slide      /data/SAMPLE1/HE/SAMPLE1_he.ome.tif \
    --xenium-bundle /data/SAMPLE1/xenium/output-XETG.../ \
    --env-name      my-custom-env
```

The wrapper sources `~/.bashrc` and has a micromamba → mamba → conda
activation fallback chain so it works on any setup that has one of those
tools available.

The wrapper builds the logs dir under either
`<output_root>/<sample>/<sample>_<run_id>/logs/logs_heRegistration/`
(integrated) or `<output_root>/<sample>/logs/` (standalone), with a leaf
name `<sample>_<jobid>_<stages>` — so re-runs at the same job id don't
clobber each other's logs.

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

**Resource sizing.** Under the default `--mode full_with_micro`, a
typical whole-slide run takes 2–2.5 h wall-clock — the `register` step is
30–60 min for rigid + non-rigid, and `register_micro` adds another 30 min
– 2 h on top depending on tissue size. Warp + celltype + viz add another
30–60 min. The 2-day `--time` in the wrapper is the honest ceiling for
`full_with_micro` on large samples; a `--mode rigid_only` diagnostic
sweep completes in well under an hour. 196 GB memory is sized for the
non-rigid + micro solve at `max_non_rigid_registration_dim_px = 10000`
plus the Dask warp cluster; drop to 128 GB only if you also drop that
cap.

## Development & testing

The test suite covers CLI parsing, config validation, layout resolution
for the three invocation modes, celltyping (auto-detect, NN mapping,
polygon cleanup), viz `render_boundaries` variants, stage imports, and an
end-to-end pipeline smoke test.

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
> 10.5281/zenodo.PLACEHOLDER), a Python wrapper around VALIS (Gatenbee et
> al., 2023) via the `valis_hest` adapter and HEST (Jaume et al., 2024)
> that composes rigid, non-rigid, and micro-registration transforms into
> a single registrar. Xenium cell and nucleus segmentations were warped
> into H&E pixel space via HEST's `warp_and_save_xenium_objects`,
> assigned cell-type labels by nearest-neighbour lookup on centroids
> against an upstream `proseg_purified.h5ad`, and rendered as per-slide
> overlays.

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
