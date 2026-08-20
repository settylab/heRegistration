# hexenium

`hexenium` takes a single-sample H&E slide and the corresponding Xenium
`output-XETG...` bundle and produces:

- a **registrar pickle** encoding the H&E ↔ Xenium DAPI transform
  (VALIS / valis_hest),
- **warped** cell + nucleus boundary parquets in H&E pixel space
  (HEST),
- **celltype-annotated GeoJSONs** (one pair for analysis, one pair for
  QuPath import),
- an **overlay PNG** for sanity-checking.

Every stage is resumable via a sentinel-file existence check; every
run snapshots its resolved config to disk for traceability.

## Pipeline shape

```
   HE slide (.ome.tif or .vsi)
                │
                ▼
   ┌─────────────────────────┐
   │ 0. he_preprocess        │  (VSI → Xenium-Explorer OME-TIFF; no-op if input is already .tif)
   └─────────────────────────┘
                │
                ▼
   ┌─────────────────────────┐
   │ 1. register             │  VALIS H&E ↔ Xenium DAPI → _registrar.pickle
   └─────────────────────────┘
                │
                ▼
   ┌─────────────────────────┐
   │ 2. warp                 │  HEST warp_and_save_xenium_objects → parquets
   └─────────────────────────┘
                │
                ▼
   ┌─────────────────────────┐
   │ 3. celltype             │  join CSV → 4× GeoJSON + combined parquet
   └─────────────────────────┘
                │
                ▼
   ┌─────────────────────────┐
   │ 4. viz                  │  matplotlib overlay PNG
   └─────────────────────────┘
```

## Reading order

New to hexenium? Read in this order:

1. [Installation](installation.md) — conda env setup, pip
   layer, common install failures.
2. [Quickstart](quickstart.md) — the canonical
   `hexenium run …` invocation.
3. [Usage](usage.md) — CLI reference, invocation modes,
   inputs, outputs.

Then, as needed:

- [Stages](stages/he_preprocess.md) — per-stage algorithm +
  parameter reference (`he_preprocess`, `register`, `warp`,
  `celltype`, `viz`).
- [Configuration](configuration.md) — full knob reference
  and override precedence.
- [HPC / Slurm](hpc.md) — the `submit_he_registration.sh`
  wrapper, unified BioFormats VSI conversion, resource
  sizing.
- [Advanced](advanced.md) — multi-registration workflows,
  `--set-default-on-success`, force-rerun semantics.
- [Troubleshooting](troubleshooting.md) — install failures,
  env drift, VSI format issues.
- [Methods](methods.md) — one-paragraph algorithm summary
  per stage, in the shape of a methods-section writeup.
- [Development](development.md) — contributing, test suite,
  code map.

!!! note "Docs are under construction"
    Most pages currently point back at the README while
    content migrates in a series of small PRs. The
    [Installation](installation.md) and
    [Methods](methods.md) pages are already complete; the
    rest will land over the next few PRs.
