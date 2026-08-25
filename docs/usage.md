# Usage

Everything runs through `hexenium run …`.

## Quick reference

```bash
hexenium run --help
```

## Worked example

The whole flag surface is exercised by a run against a multichannel
Xenium bundle (where DAPI lives at `ch0000_dapi.ome.tif`, so
`--dapi-path` is required), full VALIS pipeline (`--mode
full_with_micro`), CSV annotation, and a forced re-run of every stage
after the register pass.

**Interactive (srun --pty)** — recommended when you want live output
in your terminal:

```bash
micromamba activate heRegistration

srun --time=2:00:00 --mem=128G --cpus-per-task=10 --partition=campus-new \
    --pty \
    hexenium run \
        --sample-id     SAMPLE1 \
        --he-path       /data/SAMPLE1/HE/aligned_fullres_HE.ome.tif \
        --xenium-bundle /data/SAMPLE1/xenium/output-XETG.../ \
        --dapi-path     /data/SAMPLE1/xenium/output-XETG.../morphology_focus/ch0000_dapi.ome.tif \
        --celltype-csv  /data/SAMPLE1/annotation/SAMPLE1_celltype.csv \
        --output-root   /data/SAMPLE1/hexenium_runs \
        --mode full_with_micro \
        --stages register warp celltype viz \
        --force-rerun
```

**Background (sbatch)** — same call, wrapped by the submission script:

```bash
OUTPUT_ROOT=/data/hexenium_runs sbatch scripts/submit_he_registration.sh \
    SAMPLE1 \
    /data/SAMPLE1/HE/aligned_fullres_HE.ome.tif \
    /data/SAMPLE1/xenium/output-XETG.../ \
    /data/SAMPLE1/annotation/SAMPLE1_celltype.csv \
    --dapi-path /data/SAMPLE1/xenium/output-XETG.../morphology_focus/ch0000_dapi.ome.tif \
    --mode full_with_micro \
    --stages register warp celltype viz \
    --force-rerun
```

Positional args to `submit_he_registration.sh`:
`SAMPLE_ID  HE_PATH  XENIUM_BUNDLE  [CELLTYPE_CSV]  [--extra --flags]`.
Anything after the four positionals is forwarded verbatim to
`hexenium run`.

## Inputs

| Argument | What it is |
| --- | --- |
| `--sample-id` | Any label, e.g. `SAMPLE1`. Used for output folder names. |
| `--he-path` | OME-TIFF (used as-is by register) OR Olympus `.vsi` (stage 0 converts to Xenium-Explorer-compatible OME-TIFF). |
| `--xenium-bundle` | The `output-XETG...` directory holding `morphology_focus/`, `cell_boundaries.parquet`, etc. |
| `--dapi-path` | Standard bundles: derives from `--xenium-bundle`. Multichannel bundles: `morphology_focus/ch0000_dapi.ome.tif` — pass explicitly. |
| `--celltype-csv` | Optional. Two columns: `cell_id` (Xenium UUID) + `group`. |
| `--output-root` | Where per-sample directories land. Required (no built-in default). |

## Key flags

```
--mode {rigid_only, rigid_nonrigid, full_with_micro}
                  Default: full_with_micro. Matches stock HEST
                  register_dapi_he(micro_reg=True).

--use-he-deconvolution {true,false}
                  Default: true. Macenko-style HEDeconvolution
                  preprocessing. Flip to false on samples with
                  faint hematoxylin where deconvolution hurts.

--include-transcripts
                  Default: off. Warps transcripts alongside cells +
                  nuclei. Adds hours.

--stages he_preprocess register warp celltype viz
                  Default: all five. Restrict to re-run a subset.

--force-rerun     Nuke sentinels and redo the requested stages.
```

## Config

`hexenium` looks for its default config at
`src/hexenium/_defaults/default.yaml` inside the package tree.
Override any subset with a user YAML:

```bash
hexenium run --config my_overrides.yaml --sample-id SAMPLE1 …
```

Or spot-override on the CLI (`--use-he-deconvolution true`, etc.).
CLI overrides win over user YAML which wins over the default.

## Outputs

Land under `<output_root>/<sample_id>/`:

```
he_preprocessed/<sample>_he.ome.tif             ← stage 0 (only if input was .vsi)
registration/<run_name>/data/_registrar.pickle  ← stage 1 sentinel
warped/he_cell_seg.parquet   + .geojson         ← stage 2
warped/he_nucleus_seg.parquet + .geojson
warped/he_transcripts.parquet                   ← only with --include-transcripts
celltyped/<sample>_cells_analysis.geojson       ← stage 3 (raw floats)
celltyped/<sample>_cells_qupath.geojson         ← stage 3 (QuPath-safe)
celltyped/<sample>_nuclei_analysis.geojson
celltyped/<sample>_nuclei_qupath.geojson
celltyped/<sample>_celltyped_wholeslide.parquet
viz/<sample>_overlay.png                        ← stage 4
logs/<job_name>_<jobid>.log                     ← timestamped, per-stage
resolved_config.yaml                            ← snapshot of every param used
```

`--output-root` is required (there is no built-in default). Set it on
the CLI or via the `OUTPUT_ROOT` env var used by
`scripts/submit_he_registration.sh`.

## Resume

Each stage's outputs are sentinels. Re-running with the same
`--sample-id` + `--output-root` will skip completed stages. Useful if
warp crashes — fix it, re-run, registration won't redo.

## Optional stage: nn_celltype_mapping (proseg → Xenium celltype CSV)

Generates the two-column `(cell_id, group)` celltype CSV that
`--celltype-csv` consumes, by assigning each Xenium cell the celltype
of its nearest proseg-side neighbour.

**Opt-in** — not in the default `--stages` list. Add it explicitly:

```bash
hexenium run \
    --sample-id SAMPLE1 \
    --he-path       /data/SAMPLE1/HE/aligned_fullres_HE.ome.tif \
    --xenium-bundle /data/SAMPLE1/xenium/output-XETG.../ \
    --output-root   /data/SAMPLE1/hexenium_runs \
    --stages nn_celltype_mapping \
    --nn-proseg-source /data/SAMPLE1/proseg/SAMPLE1_metadata.csv
```

Output lands at
`<output_root>/<sample_id>/celltype_for_hexenium/<sample_id>_celltype.csv`
plus an `<sample_id>_celltype_inspection.csv` sidecar with per-cell
distance + winning proseg id.

**Inputs.** The proseg source is either a `.h5ad` (canonical, from
rctd-split's `mtx_to_h5ad` stage) or a `.csv[.gz]`. Relative paths are
joined onto `<output_root>/<sample_id>/`; `{sample_id}` substitution is
supported.

**Common knobs.**

| CLI flag | Default | What it does |
| --- | --- | --- |
| `--nn-proseg-source` | `h5ad/{sample_id}_purified.h5ad` | Path to proseg .h5ad or .csv[.gz]. |
| `--nn-celltype-col` | `first_type` | Column in proseg source holding the celltype label. |
| `--nn-k` | `1` | Number of nearest neighbours (K>1 uses `nn.tiebreak` config). |
| `--nn-distance-threshold` | none | Max NN distance (Xenium µm) for a valid assignment. |
| `--nn-unmatched-policy` | `mark_unassigned` | `drop` / `mark_unassigned` / `keep` for beyond-threshold cells. |
| `--nn-hybrid-direct-join-first` | `false` | For cells with a proseg `original_cell_id`, direct-join first; NN as fallback. |

Other knobs (algorithm, tiebreak, coord scaling, output naming,
extra_columns) live in the `nn_celltype_mapping:` block of
`src/hexenium/_defaults/default.yaml` — override via
`--config user.yaml`.

**Piping into the celltype stage.** Once the CSV lands, re-run with
`--celltype-csv <path> --stages celltype viz` to consume it:

```bash
hexenium run \
    --sample-id SAMPLE1 \
    --he-path       /data/SAMPLE1/HE/aligned_fullres_HE.ome.tif \
    --xenium-bundle /data/SAMPLE1/xenium/output-XETG.../ \
    --celltype-csv  <output_root>/SAMPLE1/celltype_for_hexenium/SAMPLE1_celltype.csv \
    --stages celltype viz
```

## Known gotchas

- **`rigid_only_micro` is NOT a valid mode** — combining rigid-only
  init with `register_micro` raises `bk_dxdy is None` deep in
  valis_hest. Use `rigid_only` or `full_with_micro`.
- **Multichannel Xenium bundles** name the DAPI `ch0000_dapi.ome.tif`
  instead of `morphology_focus_0000.ome.tif`. You **must** pass
  `--dapi-path` explicitly.
- **Celltype CSV needs Xenium UUID cell_ids** (`aaaafejg-1` style),
  not integer indices. If you only have integer IDs from proseg, do
  the proseg → Xenium UUID translation outside the pipeline before
  passing the CSV, or use the `nn_celltype_mapping` stage above.
