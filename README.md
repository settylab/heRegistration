# hexenium

**Single-sample H&E ↔ Xenium DAPI registration, warp, cell-type propagation, and overlay visualisation.**

> **v0.2.0 breaking changes** — celltype now sources labels directly from an
> upstream proseg-purified h5ad via inlined nearest-neighbour mapping
> (no more `--celltype-csv`); three invocation modes (standalone /
> integrated-by-run-id / integrated-by-h5ad); per-stage outputs are
> sharded under `register/<he_job_id>/`, `warp/<he_job_id>/`,
> `celltyped/<he_job_id>/`, `viz/<he_job_id>/`. See
> [CHANGELOG.md](CHANGELOG.md) for the full delta and migration notes.
> Some examples below still describe the v0.1.x flow; see the
> `--help` output and the CHANGELOG for the current interface.

At a glance:

- **Register** a hematoxylin and eosin (H&E) whole-slide image against the DAPI morphology channel of a matched Xenium in-situ transcriptomics run (VALIS via `valis_hest`), composing rigid, non-rigid, and micro-registration transforms into a single registrar.
- **Warp** Xenium cell + nucleus polygon boundaries from Xenium pixel space into H&E pixel space (HEST `warp_and_save_xenium_objects`).
- **Propagate** external per-cell type annotations onto the warped boundaries via a `how="left"` join on the Xenium UUID `cell_id`; optionally let nuclei inherit their sibling cell's label.
- **Visualise** with a publication-quality overlay: nucleus outlines on a downsampled H&E thumbnail, coloured by cell-type label. Add cell polygons with `--viz-render-boundaries both`.
- Every stage is idempotent — sentinel-file resume — and every run snapshots its resolved configuration to disk.

Table of contents: [What it does](#what-it-does) · [Installation](#installation) · [Quickstart](#quickstart) · [Example results](#example-results) · [Configuration](#configuration) · [HPC (Slurm)](#hpc-usage-slurm) · [Citation](#citation) · [License](#license) · [Acknowledgements](#acknowledgements)

## What it does

Five stages, run in order (`he_preprocess → register → warp → celltype → viz`). Restrict a run to a subset with `--stages`; each stage is guarded by a sentinel file, so re-running with the same `--sample-id` and `--output-root` picks up where the last run left off. Nuke sentinels with `--force-rerun`.

### `he_preprocess` — VSI → OME-TIFF (no-op for OME-TIFF inputs)

If `--he-path` points to an Olympus SlideScanner `.vsi`, this stage converts it to a pyramidal OME-TIFF that matches 10x's Xenium Explorer image-conversion specification: 1024×1024 tiles, lossless JPEG 2000 (or ZLIB), 7-level pyramid at scale 2. `OpenSlide` reads the source; `tifffile.TiffWriter` writes the pyramid; physical pixel size is propagated from OpenSlide's `mpp-{x,y}` into the OME-TIFF's `PhysicalSize{X,Y}` metadata. If the input is already an OME-TIFF, the stage returns immediately and downstream stages consume the input as-is.

**Writes:** `<output_root>/<sample_id>/he_preprocessed/<sample_id>_he.ome.tif` (only when the input is a `.vsi`).

### `register` — H&E ↔ Xenium DAPI alignment (VALIS)

Registers the H&E to the DAPI morphology image using VALIS via `valis_hest`. The default `--mode full_with_micro` composes three transforms into one registrar: a rigid initial solve, a non-rigid B-spline solve at `max_processed_image_dim_px = 1500`, and a micro-registration refinement at `max_non_rigid_registration_dim_px = 10000` — matching stock HEST `register_dapi_he(micro_reg=True)`. Two diagnostic modes are available for sweeps: `rigid_only` (fastest, rigid only) and `rigid_nonrigid` (rigid + non-rigid, no micro). Macenko-style H&E deconvolution runs before feature detection by default; on faint-hematoxylin samples it can degrade alignment, so override with `--use-he-deconvolution false` when a run shows that failure mode. Reflection checking is on by default so mirrored slides are caught cheaply. To sidestep a hardcoded `he_key='aligned_fullres_HE'` in `valis_hest`, hexenium symlinks the user's H&E to that canonical name inside the per-sample workdir before invoking VALIS.

**Writes:** a single VALIS registrar pickle (`registration/<run_name>/data/_registrar.pickle`) which serialises the composed rigid × non-rigid × micro transforms for later reuse, plus VALIS's per-stage image dumps under sibling `rigid_registration/`, `non_rigid_registration/`, `micro_registration/`, etc.

### `warp` — apply the registrar to Xenium objects

Applies the VALIS registrar to the Xenium `cell_boundaries.parquet` and `nucleus_boundaries.parquet` (and, optionally, `transcripts.parquet`) via HEST's `warp_and_save_xenium_objects`. Warping runs on a Dask `LocalCluster` in which every worker initialises the BioFormats JVM exactly once via a `WorkerPlugin`; JPype enforces a single JVM lifecycle per Python process, so the JVM is deliberately not torn down after registration. Outputs are WKB-encoded GeoPandas parquets in H&E pixel space; each row's Xenium `cell_id` is preserved for downstream joins.

**Writes:** `warped/he_cell_seg.{parquet,geojson}` and `warped/he_nucleus_seg.{parquet,geojson}`. Parquet for programmatic downstream use, GeoJSON for viewers (QuPath, GeoJSON.io, `napari-geojson`, …).

### `celltype` — join an annotation CSV onto the warped polygons

Left-joins an external per-cell type annotation CSV (`--celltype-csv`) onto the warped boundaries on the Xenium UUID `cell_id`. Because Dask can split a single polygon across pyramid tiles and emit duplicate `xenium_cell_id` rows, hexenium uses a `cumcount`-augmented merge (pair-by-occurrence rather than pair-by-id) to preserve pairings. Polygon cleanup runs `shapely.make_valid → largest connected piece`, drops empty / non-Polygon / sub-`area_threshold_px` (default 20 sq px) geometries, and sanity-checks for finite coords, ≥3 unique vertices, and `is_valid`. Nuclei without a label optionally inherit their sibling cell's label (`nuclei_inherit_classification: true`).

**Writes:** five files. `<sample>_cells_analysis.geojson` and `<sample>_nuclei_analysis.geojson` — raw floats for Python downstream. `<sample>_cells_qupath.geojson` and `<sample>_nuclei_qupath.geojson` — 2-decimal-rounded, minimal-property GeoJSON that QuPath's built-in import tool accepts. `<sample>_celltyped_wholeslide.parquet` — the analysis-precision combined table feeding the viz stage.

### `viz` — publication-shape overlay PNG

Draws warped boundaries on a downsampled H&E thumbnail. By default (`viz.render_boundaries = nucleus`) only nucleus outlines are drawn — the cleanest read of registration quality against H&E's hematoxylin-stained nuclei. Two alternatives are exposed via `--viz-render-boundaries`: `cell` draws cell polygons (semi-transparent fill) only, and `both` draws nucleus outlines on top of the cell polygons. Thumbnails are read via OpenSlide with a `tifffile` fallback. Cell-type labels drive a qualitative palette (`tab20` by default) with user-configurable overrides in `viz.classification_palette`; any label not in the override map is auto-assigned a colormap slot in first-appearance order (deterministic across re-runs). The renderer is dask-friendly: polygon batches are drawn in chunks so a whole-slide overlay does not need to fit its geometries in memory as a single frame.

**Writes:** `viz/<sample>_overlay.png` at `viz.dpi` (default 200 DPI), with a legend of the labels present in the data.

## Installation

`hexenium` targets Python 3.11 and depends on VALIS, HEST, and their Java/BioFormats bridge (`jpype1`, `openjdk=11`). The `heRegistration` conda environment below is the exact env this package was developed and validated against.

```bash
# 1. Install micromamba (skip if you already have it)
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# 2. Create the env
micromamba create -n heRegistration -c conda-forge -c bioconda \
    python=3.11 \
    numpy=1.26.4 \
    scipy \
    pandas \
    pyarrow=23 \
    pyyaml \
    matplotlib=3.10 \
    scikit-image=0.19 \
    scikit-learn=1.8 \
    shapely=2.1 \
    geopandas=1.1 \
    dask=2026.1 \
    dask-geopandas=0.5 \
    distributed \
    openslide-python=1.4 \
    tifffile=2026.1 \
    pyvips=2.2 \
    openjdk=11 \
    jpype1=1.6 \
    anndata=0.11 \
    scanpy=1.11

micromamba activate heRegistration

# 3. The image-registration libraries live on PyPI + GitHub
pip install valis-wsi==1.1 valis-hest==0.0.2
pip install "hest @ git+https://github.com/mahmoodlab/HEST.git"

# 4. Editable install of hexenium
pip install -e /path/to/xenium-he-registration

# 5. Verify
python -c "import hest, valis_hest, dask, openslide; print('OK')"
hexenium --help
```

Alternatively, `environment.yml` in this repo carries the full set — `micromamba env create -f environment.yml` reproduces the env in one step.

## Quickstart

A bare `hexenium run` runs the full pipeline under its defaults: `--mode full_with_micro`, `--use-he-deconvolution true`, `--viz-render-boundaries nucleus`. That yields the composed rigid × non-rigid × micro registrar, deconvolved H&E for the feature stage, and a nucleus-outline overlay. A concrete invocation (multichannel Xenium bundle, so `--dapi-path` is passed explicitly):

```bash
micromamba activate heRegistration

hexenium run \
    --sample-id       SAMPLE1 \
    --he-path         /data/SAMPLE1/HE/aligned_fullres_HE.ome.tif \
    --xenium-bundle   /data/SAMPLE1/xenium/output-XETG.../ \
    --dapi-path       /data/SAMPLE1/xenium/output-XETG.../morphology_focus/ch0000_dapi.ome.tif \
    --celltype-csv    /data/SAMPLE1/annotation/SAMPLE1_celltype.csv \
    --output-root     /data/SAMPLE1/hexenium_runs
```

For a figure that shows both cell polygons and nucleus outlines together, add `--viz-render-boundaries both`. For a faint-hematoxylin sample where deconvolution hurts alignment, flip `--use-he-deconvolution false`. For a diagnostic sweep, `--mode rigid_only` or `--mode rigid_nonrigid` restricts the transform stack.

The batch-submission form via `sbatch --wrap`:

```bash
sbatch --partition=campus-new --time=2-00:00:00 --mem=128G --cpus-per-task=10 \
    --job-name=hexenium-SAMPLE1 \
    --wrap "micromamba activate heRegistration && hexenium run \
        --sample-id       SAMPLE1 \
        --he-path         /data/SAMPLE1/HE/aligned_fullres_HE.ome.tif \
        --xenium-bundle   /data/SAMPLE1/xenium/output-XETG.../ \
        --dapi-path       /data/SAMPLE1/xenium/output-XETG.../morphology_focus/ch0000_dapi.ome.tif \
        --celltype-csv    /data/SAMPLE1/annotation/SAMPLE1_celltype.csv \
        --output-root     /data/SAMPLE1/hexenium_runs"
```

The `sbatch --wrap` form skips the tee split entirely — Slurm's own `slurm-<jobid>.out` in `cwd` captures everything, and `python -u` (which `hexenium` sets internally via `PYTHONUNBUFFERED=1`) keeps stdout line-buffered.

## Example results

A completed run at `<output_root>/<sample_id>/` writes the following tree (sizes shown for a whole-slide TMA-scale run):

```
SAMPLE1/
├── celltyped/
│   ├── SAMPLE1_cells_analysis.geojson       (~160 MB)
│   ├── SAMPLE1_cells_qupath.geojson         (~129 MB)
│   ├── SAMPLE1_celltyped_wholeslide.parquet (~ 89 MB)
│   ├── SAMPLE1_nuclei_analysis.geojson      (~148 MB)
│   └── SAMPLE1_nuclei_qupath.geojson        (~119 MB)
├── logs/
├── registration/
│   ├── manifest.yaml
│   ├── SAMPLE1_<UTC-timestamp>/
│   │   ├── data/           ← _registrar.pickle + _summary.csv
│   │   ├── deformation_fields/
│   │   ├── masks/
│   │   ├── micro_registration/
│   │   ├── non_rigid_registration/
│   │   ├── overlaps/
│   │   ├── processed/
│   │   └── rigid_registration/
│   └── _workdir/
├── resolved_config.yaml
├── viz/
│   └── SAMPLE1_overlay.png                  (~17 MB)
└── warped/
    ├── he_cell_seg/           ← Dask-partitioned parts (part.0.parquet, …)
    ├── he_cell_seg.geojson                           (~135 MB)
    ├── he_cell_seg.parquet                           (~ 46 MB)
    ├── he_nucleus_seg/        ← Dask-partitioned parts
    ├── he_nucleus_seg.geojson                        (~128 MB)
    └── he_nucleus_seg.parquet                        (~ 44 MB)
```

### `registration/`

`_workdir/` holds the (symlinked) canonical-name H&E used to work around VALIS's hardcoded `aligned_fullres_HE` key. The timestamped slot `SAMPLE1_<UTC-timestamp>/` is VALIS's output directory — each subdir is one solve stage. The one file that matters programmatically is `.../data/_registrar.pickle`: it carries the composed **rigid × non-rigid × micro** transforms, and is the sentinel every downstream stage checks. The `manifest.yaml` records the real input H&E path, the symlinked path, the DAPI path, the resolved mode, and any VALIS kwargs the installed `valis-hest` version silently dropped — useful when reproducing a run against a different `valis-hest` version.

### `warped/`

Polygon geometries of Xenium cells and nuclei aligned into H&E coordinates.

- `he_cell_seg.parquet` / `he_nucleus_seg.parquet` — WKB-encoded GeoPandas parquets. This is the format the `celltype` stage reads, and the format you want for any programmatic downstream (Python / `geopandas`, or arrow-native tools).
- `he_cell_seg.geojson` / `he_nucleus_seg.geojson` — the same polygons in GeoJSON, ready for viewers (QuPath, GeoJSON.io, `napari-geojson`).

Every polygon is tagged with its Xenium `cell_id`; the next stage joins on that. The `he_cell_seg/` and `he_nucleus_seg/` sibling directories are Dask's partitioned parquet parts written during the warp; the single `.parquet` files next to them are the consolidated view most consumers want.

### `celltyped/`

Five files. The split between `_analysis` and `_qupath` is intentional:

- `<sample>_cells_analysis.geojson`, `<sample>_nuclei_analysis.geojson` — raw-float coordinates, all properties preserved. Use these for Python downstream (`geopandas.read_file` + `.merge` back to your Xenium AnnData).
- `<sample>_cells_qupath.geojson`, `<sample>_nuclei_qupath.geojson` — 2-decimal-rounded, minimal-property GeoJSON. QuPath's built-in import tool is strict about polygon validity; these have been passed through `shapely.make_valid → largest connected piece → buffer(0)`, so QuPath loads them directly (drag-and-drop onto an opened slide) without invalid-polygon errors.
- `<sample>_celltyped_wholeslide.parquet` — the combined analysis-precision table the viz stage renders from.

The `--celltype-csv` argument is a two-column CSV whose `cell_id` column joins against the Xenium `cell_id` UUID and whose `group` column carries the label. See [Known gotchas](docs/usage.md#known-gotchas) for the CSV shape and ID-space conventions.

### `viz/`

`<sample>_overlay.png` — the overlay figure. Xenium cell polygons (semi-transparent fill, `cell_alpha=0.3`) and/or nucleus outlines (`nucleus_alpha=0.6`) drawn over a downsampled H&E thumbnail (`viz.thumbnail_max_dim=4096`), coloured per cell-type label from the annotation CSV. The palette in `viz.classification_palette` (see `config/default.yaml`) can be overridden per label; any label not in the override map is auto-assigned a `tab20` slot in first-appearance order (deterministic across re-runs).

### `resolved_config.yaml`

Snapshot of every parameter that shaped the run — CLI overrides, user-YAML overrides, and defaults, merged. Read this after any run to confirm what defaults were actually in force (e.g. `celltype.area_threshold_px`, `viz.thumbnail_max_dim`, the palette). If a downstream analysis surprises you, `resolved_config.yaml` is the first place to look.

## Configuration

Every knob is defined in `config/default.yaml` with a comment explaining its effect. Override any subset with a user YAML (`hexenium run --config user.yaml …`) or spot-override on the CLI (`--mode rigid_only`, `--use-he-deconvolution false`, …). CLI overrides win over user YAML which wins over the default.

### Most-tuned knobs

The knobs most runs actually touch:

| Key | Default | Meaning |
| --- | --- | --- |
| `--mode` / `registration.mode` | `full_with_micro` | `rigid_only`, `rigid_nonrigid`, or `full_with_micro`. `rigid_only_micro` is invalid — see `docs/usage.md`. |
| `--use-he-deconvolution` / `registration.use_he_deconvolution` | `true` | Macenko-style H&E stain deconvolution before registration. Override with `false` on faint-hematoxylin samples. |
| `--dapi-path` / `dapi_path` | derived from `xenium_bundle` | Pass explicitly for multichannel bundles (`morphology_focus/ch0000_dapi.ome.tif`). |
| `--celltype-csv` / `celltype.csv` | `null` | Two columns: `cell_id` (Xenium UUID) + `group`. Omit to skip celltyping. |
| `--csv-id-col` / `celltype.csv_id_col` | `cell_id` | Column joining to the Xenium `cell_id`. |
| `--csv-group-col` / `celltype.csv_group_col` | `group` | Column holding the cell-type label. |
| `celltype.area_threshold_px` | `20` | Drop polygons smaller than this (in H&E pixels²) as segmentation noise. |
| `celltype.nuclei_inherit_classification` | `true` | Unlabelled nuclei inherit their sibling cell's label. |
| `viz.thumbnail_max_dim` | `4096` | Downsample H&E to this max dimension for the overlay. |
| `viz.dpi` | `200` | Overlay PNG DPI. |
| `--viz-render-boundaries` / `viz.render_boundaries` | `nucleus` | Which boundaries to draw on the overlay: `nucleus` (outlines only), `cell` (polygons only), or `both`. |
| `--include-transcripts` / `warp.include_transcripts` | `false` | Also warp `transcripts.parquet`. Adds hours. |

### Registration algorithm knobs

The `register` stage's numeric parameters live under `registration:` and `parameters:` in `config/default.yaml`. Every one is CLI-overridable (dash-cased flag) and YAML-overridable (dot path).

| Key | Default | Meaning |
| --- | --- | --- |
| `--check-for-reflections` / `registration.check_for_reflections` | `true` | Have VALIS test mirrored/flipped candidates during rigid init. Cheap; catches mounted-backwards slides. |
| `--create-masks` / `registration.create_masks` | `false` | Auto-mask non-tissue regions before registration. Off by default; auto-mask can over-mask small tissue islands. |
| `--align-to-reference` / `registration.align_to_reference` | `true` | Use H&E as the reference image (rather than DAPI). |
| `--max-processed-image-dim-px` / `parameters.max_processed_image_dim_px` | `1500` | Long-edge pixel cap for the rigid + first-pass non-rigid registration. Higher = more accurate + more memory. |
| `--max-non-rigid-registration-dim-px` / `parameters.max_non_rigid_registration_dim_px` | `10000` | Long-edge pixel cap for the non-rigid B-spline solve (used by `rigid_nonrigid` and `full_with_micro`). |
| `parameters.micro_rigid_registrar_cls` | `null` | Override VALIS's `MicroRigidRegistrar` class. Advanced; leave `null` for the stock class. |
| `parameters.micro_rigid_registrar_params` | `{}` | Extra kwargs forwarded to the micro-rigid registrar. |
| `--run-name` / `parameters.name` | `null` (auto) | Registrar run-name folder inside `registration/`. Auto value is `<sample_id>_<UTC-timestamp>`. |

See `config/default.yaml` for the full list (Dask worker settings, JVM memory, palette overrides, per-stage rounding for QuPath vs. analysis outputs).

## HPC usage (Slurm)

Two supported invocation shapes:

**(a) via the wrapper** — `scripts/submit_he_registration.sh` embeds a resource-sized `#SBATCH` header (`--partition=campus-new --time=2-00:00:00 --mem=128G --cpus-per-task=10` — tune to your cluster), sources the user's `~/.bashrc` to make micromamba's shell function available in the non-interactive Slurm shell, activates `$ENV_NAME` (defaulting to `heRegistration`), and runs `hexenium run ...` with `PYTHONUNBUFFERED=1`. It splits its I/O routing based on execution context: under `sbatch` (no controlling terminal), everything after the pre-flight validation goes through a plain `exec >> "$LOG_FILE" 2>&1` — no `tee`, no process substitution — because a chatty subprocess (VALIS/Java) can burst output faster than `tee` flushes, filling `tee`'s pipe buffer and (under `set -euo pipefail`, no `SIGPIPE` escape hatch) deadlocking the whole job silently on its full allocation. Under interactive invocation (`bash submit_he_registration.sh`, or `srun`), it keeps the `stdbuf -oL tee` split so you get live console echo alongside the log file. Log lands at `<output_root>/<sample_id>/logs/<job_name>_<jobid>.log`.

```bash
OUTPUT_ROOT=/data/hexenium_runs sbatch scripts/submit_he_registration.sh \
    SAMPLE_ID HE_PATH XENIUM_BUNDLE [CELLTYPE_CSV] [--extra --flags]
```

Resource sizing: under the default `--mode full_with_micro`, a typical whole-slide run takes 2–2.5 h wall-clock — the `register` step is 30–60 min for rigid + non-rigid, and `register_micro` adds another 30 min – 2 h on top depending on tissue size. Warp + celltype + viz add another 30–60 min. The 2-day `--time` pinned in the wrapper is the honest ceiling for `full_with_micro` on large samples; a `--mode rigid_only` diagnostic sweep completes in well under an hour. 128 GB memory is sized for the non-rigid + micro solve at `max_non_rigid_registration_dim_px = 10000`; drop to 64 GB only if you also drop that cap.

**(b) direct `sbatch --wrap`** — the shape shown in [Quickstart](#quickstart). `sbatch --wrap` bypasses the wrapper entirely; you carry the `#SBATCH` flags on the `sbatch` command line and let Slurm's own `slurm-<jobid>.out` capture stdout. Prefer it when you want to override the resource envelope per-run or when the wrapper's env-activation logic doesn't fit your setup.

## Citation

If you use hexenium in a publication, cite via the metadata in [`CITATION.cff`](CITATION.cff) — GitHub renders that file as a copy-paste-ready BibTeX / APA / CFF block on the repo landing page. Methods-section template (fill in the version tag and Zenodo DOI):

> Registration of the H&E whole-slide image to the matched Xenium DAPI morphology channel was performed with hexenium v0.1.0 (DOI: 10.5281/zenodo.PLACEHOLDER), a Python wrapper around VALIS (Gatenbee et al., 2023) via the `valis_hest` adapter and HEST (Jaume et al., 2024) that composes rigid, non-rigid, and micro-registration transforms into a single registrar. Xenium cell and nucleus segmentations were warped into H&E pixel space via HEST's `warp_and_save_xenium_objects`, joined to external cell-type annotations on the Xenium `cell_id`, and rendered as per-slide overlays.

Please also cite VALIS and HEST directly per their upstream requests.

## License

MIT — see [`LICENSE`](LICENSE).

## Acknowledgements

- **VALIS** — Gatenbee, C. D. et al. *VALIS: Virtual Alignment of pathoLogy Image Series*. Nature Communications 14, 4502 (2023).
- **HEST** — Jaume, G. et al. *HEST-1k: A Dataset for Spatial Transcriptomics and Histology Image Analysis*. NeurIPS Datasets and Benchmarks (2024). Code: `mahmoodlab/HEST`.
- **Xenium In Situ platform** — 10x Genomics. The Xenium bundle layout (`morphology_focus/`, `cell_boundaries.parquet`, `nucleus_boundaries.parquet`, `transcripts.parquet`) is a 10x-Genomics output specification.

hexenium was adapted from an internal H&E-Xenium registration workflow.
