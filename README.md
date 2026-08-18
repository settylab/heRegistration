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

The Python package + console script is **`hexenium`** (`pip install -e .`
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

`hexenium` targets Python 3.10+ (validated on 3.11) and depends on VALIS,
HEST, and their Java/BioFormats bridge (`jpype1`, `openjdk=11`). The
`heRegistration` conda env in `environments/heRegistration.yml` is the
exact env this package was developed and validated against.

### Recommended: env from `environments/heRegistration.yml` + pinned pip layer

```bash
# 1. Install micromamba (skip if you already have it)
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# 2. Create the env from the pinned conda spec. The yml file does NOT
#    declare a name, so pass one with -n; pick whatever you like — the
#    wrapper below defaults to `heRegistration`, so using that keeps the
#    defaults working without further overrides.
cd /path/to/xenium-he-registration
micromamba env create -n <your-env-name> -f environments/heRegistration.yml
micromamba activate <your-env-name>

# 3. Install the pinned pip layer with `--no-deps` (torch / transformers /
#    ultralytics / hest / valis-wsi / …). `--no-deps` is load-bearing —
#    valis-wsi 1.1's declared metadata caps `pandas<2`, `pyvips<3` and
#    `scikit-image<0.20`, but its actual code paths run fine with the
#    newer versions this env has; without `--no-deps` pip's resolver
#    refuses. The flag can't be inlined in the yaml pip: block
#    (micromamba treats it as a package name) or in the requirements.txt
#    (pip itself refuses), so it lives on the CLI here.
pip install --no-deps -r environments/heRegistration-requirements.txt

# 4. Editable install of hexenium
pip install -e .

# 5. Verify
python -c "import hest, valis_hest, dask, openslide; print('OK')"
hexenium --version
hexenium run --help
```

Together, `environments/heRegistration.yml` (conda solve: system libs
+ scientific-Python core) and
`environments/heRegistration-requirements.txt` (pip layer: pinned
overrides + `hest @ git+https://github.com/mahmoodlab/HEST.git`)
reproduce the exact working env. If either verification command errors,
check that the active env is the one you created (not `base`) and that
both `pip install` steps returned successfully. If you chose a name
other than `heRegistration`, set `ENV_NAME=<your-env-name>` when using
the sbatch wrapper below.

### Manual install (if you can't use the env file)

The two `pip install --no-deps` calls below mirror what the yaml
does: `valis-wsi 1.1`'s declared metadata caps `pandas<2` and
`pyvips<3`, but its actual code runs fine with pandas 2.x / pyvips
3.x. Without `--no-deps` pip refuses to install it against a
pandas>=2 env.

```bash
micromamba create -n heRegistration -c conda-forge -c bioconda \
    python=3.11 'numpy<2' 'pandas<3' scipy pyarrow pyyaml \
    scanpy dask-geopandas shapely proj pyproj ipykernel \
    importlib_metadata 'pycparser>=2.14' \
    libvips pyvips imagemagick openslide openjdk=11
micromamba activate heRegistration

pip install --no-deps valis-wsi==1.1.0 valis_hest==0.0.2
pip install --no-deps "hest @ git+https://github.com/mahmoodlab/HEST.git@v1.2.0" hestcore==1.0.4
# hest + valis pip-only transitive deps (torch / transformers /
# ultralytics / spatialdata / opencv-*, anndata 0.12 override,
# etc.) — see the `- pip:` block of environments/heRegistration.yml
# for the full pinned set. Or just use the yml file, it's shorter.
pip install -e /path/to/xenium-he-registration
```

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
    [--force-rerun]
```

`OUTPUT_ROOT=/data/xenium_runs` in your env is a fallback for
`--output-root`. `ENV_NAME` overrides the conda env name (default:
`heRegistration`). The wrapper sources `~/.bashrc` and has a
micromamba → mamba → conda activation fallback chain so it works on any
setup that has one of those tools available.

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
