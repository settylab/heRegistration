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

1. [`install.md`](install.md) — set up the `heRegistration` conda env
   (`environments/heRegistration.yml`) and `pip install -e .` the
   package.
2. [`usage.md`](usage.md) — CLI reference, config surface, worked
   examples.
3. [`methods.md`](methods.md) — one-paragraph algorithm summary per
   stage, in the shape of a methods-section writeup.
