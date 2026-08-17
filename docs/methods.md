# Methods

Draft methods-section text, one paragraph per stage. Use as the
starting point for a manuscript. Cite the underlying libraries with
their peer-reviewed papers, not GitHub URLs.

## Overview

hexenium (`v0.1.0`) is a Python package that registers a hematoxylin
and eosin (H&E) whole-slide image against the DAPI morphology channel
of a matched Xenium in-situ transcriptomics run, warps the Xenium
cell + nucleus segmentations into H&E pixel space, joins them to an
external cell-type annotation table, and emits per-slide overlay
visualisations. Every stage is idempotent (sentinel-file resume);
every run snapshots its resolved parameters to disk for traceability.

## 0. VSI pre-conversion (`he_preprocess`)

Olympus SlideScanner `.vsi` inputs are converted to pyramidal
OME-TIFF using OpenSlide (`openslide.OpenSlide`) as the reader and
`tifffile.TiffWriter` as the writer. Output matches 10x's Xenium
Explorer image-conversion specification: 1024×1024 pixel tiles,
lossless JPEG 2000 (`compressionargs={level: 100}`) or ZLIB
compression, and a 7-level pyramid at scale 2. Physical pixel size is
propagated from OpenSlide's `openslide.mpp-{x,y}` properties into the
OME-TIFF's `PhysicalSize{X,Y}` metadata. Downsample interpolation is
`cv2.INTER_AREA`. OME-TIFF inputs bypass this stage.

## 1. Registration (`register`)

H&E ↔ DAPI registration is performed with VALIS via the `valis_hest`
adapter. We drive VALIS directly rather than through HEST's
`register_dapi_he` wrapper so we can expose every VALIS knob as a
config parameter. Three registration modes are supported:

- **`rigid_only`**: rigid solve only. Can match or improve over
  non-rigid variants on some sample cohorts.
- **`rigid_nonrigid`**: rigid + non-rigid B-spline solve at
  `max_processed_image_dim_px = 1500` (default).
- **`full_with_micro`** (default): rigid + non-rigid +
  micro-registration refinement at
  `max_non_rigid_registration_dim_px = 10000` (matches stock HEST
  `register_dapi_he(micro_reg=True)`).

Optional Macenko-style H&E deconvolution
(`valis_hest.preprocessing.HEDeconvolution`) is on by default; on
samples with faint hematoxylin it can degrade registration and should
be disabled with `--use-he-deconvolution false`. Reflection-checking
is on by default so mirrored slides are caught cheaply during the
rigid solve.

To sidestep a hardcoded slide-dict key in `valis_hest._post_register`
(the H&E image is expected to be named `aligned_fullres_HE`), the
pipeline symlinks the user's H&E to that canonical name inside the
per-sample workdir before invoking VALIS. This preserves the source
filename while letting VALIS's post-register housekeeping resolve
correctly.

The output is a single VALIS registrar pickle
(`_registrar.pickle`) which serialises the composed rigid + non-rigid
+ micro transforms for later reuse.

## 2. Xenium object warp (`warp`)

The VALIS registrar is applied to Xenium `cell_boundaries.parquet`,
`nucleus_boundaries.parquet`, and (optionally) `transcripts.parquet`
via HEST's `warp_and_save_xenium_objects`. Warping runs on a Dask
`LocalCluster` in which every worker initialises the BioFormats JVM
exactly once via a `WorkerPlugin`; JPype enforces a single JVM
lifecycle per Python process, so we deliberately do NOT tear down the
JVM after Stage 1 (subsequent HEST calls reuse it). Outputs are WKB-
encoded GeoPandas parquets in H&E pixel space; each row's Xenium
`cell_id` is preserved in the `__null_dask_index__` column, which the
next stage hoists to a `xenium_cell_id` column.

## 3. Cell-type annotation (`celltype`)

Warped boundaries are joined to a per-sample cell-type annotation
CSV on the Xenium UUID `cell_id`. Because Dask can split a single
polygon across pyramid tiles and emit duplicate `xenium_cell_id`
rows, we use a `cumcount`-augmented merge (pair-by-occurrence rather
than pair-by-id) to preserve pairings. Polygon cleanup runs
`shapely.make_valid → largest connected piece`, drops empty /
non-Polygon / sub-`area_threshold_px` (default 20 sq px) geometries,
optionally rounds exterior-ring coordinates, and finally sanity-
checks for finite coords, ≥3 unique vertices, and `is_valid`.

Nuclei without a label optionally inherit their sibling cell's label
(`nuclei_inherit_classification=true`). Four GeoJSONs are written per
sample — cells and nuclei, each at two precisions:

- **`*_analysis.geojson`**: raw floats, for downstream Python
  analysis.
- **`*_qupath.geojson`**: 2-decimal-rounded, minimal-property
  GeoJSON that QuPath's built-in import tool accepts.

A combined `<sample>_celltyped_wholeslide.parquet` (analysis-
precision) is also written for the visualisation stage.

## 4. Overlay visualisation (`viz`)

The visualisation stage draws warped cell polygons (semi-transparent
fill) and nucleus outlines on a downsampled H&E thumbnail. Thumbnails
are read via OpenSlide with a `tifffile` fallback. Cell classification
labels drive a qualitative palette (`tab20` by default) with
user-configurable overrides in `viz.classification_palette`; any
label not in the override map is auto-assigned a colormap slot in
first-appearance order (deterministic across re-runs). Output is a
per-sample PNG at 200 DPI, with a legend of the labels present in the
data.

## Reproducibility

All numeric libraries are pinned in `pyproject.toml` /
`environments/heRegistration.yml`.
Each run writes a `resolved_config.yaml` snapshot to
`<output_root>/<sample_id>/`, capturing every parameter — CLI
overrides, user-YAML overrides, and defaults — that shaped the run.
Registration also emits a per-run `manifest.yaml` recording the real
input H&E path, the (possibly symlinked) VALIS input path, the DAPI
path, the resolved mode, and any VALIS kwargs the installed version
did not accept.
