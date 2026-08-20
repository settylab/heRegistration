# Troubleshooting

*Content coming in PR 4 (migration from README).*

Will cover:

- Install failures (`ResolutionImpossible`,
  `ModuleNotFoundError`, VIPS warnings, `openslide-bin` sdist
  builds).
- Env drift and `--force-reinstall` gotchas (numpy 1.x pin,
  xarray, fastcluster, opencv-contrib).
- VSI format issues (OpenSlide-unsupported Olympus subtypes,
  BioFormats workaround).
- HEST-shim gotchas (missing `warp_and_save_xenium_objects`,
  parquet index, `pixel_size_morph`).

In the meantime, the [Installation
guide](installation.md#troubleshooting) carries the
full postmortems and per-symptom recovery commands.
