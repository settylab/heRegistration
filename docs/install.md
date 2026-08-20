# Install

`hexenium` depends on VALIS + HEST + their Java/BioFormats bridge, so
reproducing the working env is a **two-step install**: a conda-side
solve for the system libraries and scientific-Python core, then a
`--no-deps` pip layer for `hest` / `valis-wsi` / `valis-hest` and
their transitive stack. Both steps are captured in
`environments/heRegistration.yml` and
`environments/heRegistration-requirements.txt`.

`--no-deps` on the pip step is **load-bearing** and must live on the
CLI — see [Why `--no-deps`](#why---no-deps-is-mandatory) below. This
mirrors the recommended flow in the top-level README; keep both in
sync.

## Prerequisites

- **micromamba** (or mamba/conda). Fresh install:
  ```bash
  "${SHELL}" <(curl -L micro.mamba.pm/install.sh)
  ```
- A working checkout of hexenium.

## 1. Clone the repo

```bash
git clone https://github.com/settylab/heRegistration.git
cd heRegistration
```

## 2. Create the conda env

The yml intentionally omits a `name:` field; pass one with `-n`.
`heRegistration` is the name the sbatch wrapper (`ENV_NAME`) defaults
to — use it unless you have a reason not to.

```bash
micromamba env create -n heRegistration -f environments/heRegistration.yml
# or with conda:
# conda env create -n heRegistration -f environments/heRegistration.yml
```

## 3. Activate

```bash
micromamba activate heRegistration
# or: conda activate heRegistration
```

## 4. Install the pinned pip layer with `--no-deps`

```bash
pip install --no-deps -r environments/heRegistration-requirements.txt
```

This installs `hest` (from the pinned `v1.2.0` git tag),
`valis-wsi 1.1.0`, `valis_hest 0.0.2`, and the ~160 supporting
packages (torch, transformers, ultralytics, spatialdata, opencv,
anndata 0.12 override, …) at the exact versions the working env was
validated against. `hest` is already the first entry of the
requirements file, so no separate `pip install hest` step is needed.

### Why `--no-deps` is mandatory

`valis-wsi 1.1`'s declared `Requires-Dist` metadata caps
`pandas<2.0`, `pyvips<3.0`, and `scikit-image<0.20`. HEST has similar
declared-vs-actual divergences. The working env runs the *newer*
pandas 2.3.x / pyvips 3.1.x / scikit-image 0.19.x branch; VALIS's
actual code never touches the pandas-1 or pyvips-2 API, so the caps
are over-defensive relative to the code hexenium actually calls, but
pip's resolver refuses to install VALIS against pandas 2.x without
`--no-deps`. Skipping it yields:

```
ResolutionImpossible: valis-wsi 1.1.0 depends on pandas<2.0.0
```

`--no-deps` cannot be inlined in the yml `pip:` block (micromamba
treats each entry as a package name) or the top of the requirements
file (pip refuses), so it lives on the CLI here.

## 5. Install the hexenium package

For end users:

```bash
pip install .
```

For developers (editable):

```bash
pip install -e .
```

Non-editable is recommended for end users: an editable install
exposes the source tree to `sys.path`, so a stray import via a
working-directory Python — or a sibling `hexenium/` folder in
`cwd` — can shadow the installed package and silently pull in
half-updated modules.

Both commands resolve `pyproject.toml`'s declared deps against the
env you already prepared in steps 2-4, and the pins are loose
enough to accept what's already installed (no re-downloads
expected). Notably, `pyproject.toml` declares
`opencv-contrib-python` — matching what step 2 pinned + what
`valis_hest` requires — so this step doesn't clobber the contrib
`cv2/` files with stock ones. See the OpenCV drift entry in
"Bugs & fixes" below if you're recovering an env that ended up
with `opencv-python` on top of `opencv-contrib-python`.

## 6. (Optional) HEST from an editable local clone

The requirements file already installs HEST from the `mahmoodlab/HEST`
`v1.2.0` tag. If you want an editable clone for local hacks:

```bash
git clone --branch v1.2.0 https://github.com/mahmoodlab/HEST.git ~/HEST
pip install --no-deps -e ~/HEST
```

## 7. Verify

```bash
python -c "import hest, valis_hest, valis_hest.registration, valis_hest.slide_io, dask, openslide; print('OK')"
hexenium --version
hexenium run --help
```

If any of these errors, check that the active env is the one you
created (not `base`) and see the Troubleshooting section below.

## Troubleshooting

- **`ResolutionImpossible: valis-wsi 1.1.0 depends on pandas<2.0.0`**
  (or similar for `pyvips` / `scikit-image`) during step 4. You forgot
  `--no-deps` — re-run step 4 with the flag.
- **`ModuleNotFoundError: No module named 'valis_hest'`** (or `hest`
  / `torch` / `transformers` / `ultralytics`) when you run `hexenium`.
  Step 4 was skipped or silently failed. Re-run it and re-check with
  `pip list | grep -Ei 'valis|hest|torch|transformers'`. Expected:
  `valis_hest 0.0.2`, `valis-wsi 1.1.0`, `hest 1.1.1` (installed from
  the `v1.2.0` tag — its internal version string is `1.1.1`),
  `hestcore 1.0.4`, `torch 2.6.0`, `transformers 5.1.0`,
  `ultralytics 8.4.14`.
- **`VIPS-WARNING: unable to load "vips-magick.so" ...
  libMagickCore-7.Q16HDRI.so.10: cannot open shared object file`** at
  first `import pyvips`. `imagemagick` is in the yml to prevent this;
  the warning should not fire on a freshly-created env. If it does,
  your conda solve is stale: `micromamba clean -a` then re-create.
  The warning itself is benign (the pipeline never touches that
  loader); `VIPS_WARNING=off` silences it as a temporary workaround.
- **Micromamba dependency cascade** at step 2 (`libvips` /
  `gdk-pixbuf` / `librsvg` / `scanpy` "no viable options"). Almost
  always a stale local cache or a stale `~/.condarc`. Fix in order:
  1. `micromamba clean -a`.
  2. If `~/.condarc` is a stale NFS handle (`[Errno 116] Stale file
     handle`), refresh it: `rm ~/.condarc && touch ~/.condarc` (or
     restore your original file).
  3. Re-run step 2 on a networked node.
- **`openslide-bin` sdist build fails** at step 4 with "Install
  OpenSlide from source". This is what the `@v1.2.0` pin on HEST
  guards against — HEST HEAD depends on TRIDENT which pulls
  `openslide-bin`, and its wheels don't cover older glibc. If you're
  hitting this, verify line 1 of `environments/heRegistration-requirements.txt`
  is `hest @ git+https://github.com/mahmoodlab/HEST.git@v1.2.0` and
  not an unpinned `hest @ git+…HEST.git`.
- **`hexenium` command not found** after step 5 succeeds. The wrong
  env is active — `which hexenium` and re-run `micromamba activate
  heRegistration`.
- **Sbatch job fails at env activation** with "environment
  `heRegistration` not found". You created the env under a different
  name in step 2. Either recreate as `heRegistration`, or pass the
  name to the sbatch wrapper: `./scripts/submit_he_registration.sh
  --env-name <your-env-name> …` (or export `ENV_NAME=<your-env-name>`
  in your shell).

If none of these match, `pip install --no-deps -r
environments/heRegistration-requirements.txt -v` prints per-package
progress and surfaces which entry pip is choking on — most useful
when a wheel has been yanked from PyPI or a git ref has moved.

## Env drift & recovery

These are bugs that a live env can drift into AFTER a clean install
succeeds — usually because someone ran `pip install --force-reinstall
<pkg>` (or an unpinned `pip install --upgrade <pkg>`) which
re-resolves transitive dependencies and pulls newer versions of
numpy, xarray, or their kin. The pins in `heRegistration.yml` +
`heRegistration-requirements.txt` block the initial-install form of
each drift; the recovery commands below fix an env that has already
drifted.

Each fix is diagnosed on `settylab/TracyY123-nexus#15`; comment IDs
are linked so you can retrace the debugging thread.

### The load-bearing rule

**Do NOT run `pip install --force-reinstall <pkg>` without
`--no-deps` OR a co-pinned `numpy==1.26.4`.** `--force-reinstall`
re-resolves transitive dependencies, and any package with a loose
`numpy>=1.22` bound will pull numpy 2.x — silently breaking
`fastcluster`'s compiled extension (which was built against numpy
1.x's C API) and any other C-extension package pinned against the
1.x ABI. If you must force-reinstall a package, add `--no-deps` OR
pin numpy in the same command:

```bash
# safe:
pip install --force-reinstall --no-deps "xarray==2023.10.1"
pip install --force-reinstall "xarray==2023.10.1" "numpy==1.26.4"

# UNSAFE — silently upgrades numpy to 2.x:
pip install --force-reinstall "xarray==2023.10.1"
```

The same rule applies to `pip install --upgrade` without a target
list — `pip install --upgrade` on an under-pinned env drifts numpy
just as reliably.

### Bugs & fixes

- **`AttributeError: module 'pandas.arrays' has no attribute
  'NumpyExtensionArray'`** — first `import anndata` (or any code path
  that imports `dask.array`, which optionally imports `xarray`).
  Root cause: your env's `xarray` has drifted to `2026.1.0`, which
  accesses `pd.arrays.NumpyExtensionArray` at import time — a name
  that pandas added in 2.1. On `pandas 1.5.3` the attribute is
  missing and xarray blows up at import. `dask.array` marks xarray
  as optional (`errors="ignore"`) but that swallow-branch only
  catches `ImportError`, not `AttributeError`, so a broken xarray
  crashes the whole chain even though nothing register+warp calls
  actually needs xarray. Full diagnostic:
  [comment 5335874791](https://github.com/settylab/TracyY123-nexus/issues/15#issuecomment-5335874791).

  **Fix:**
  ```bash
  pip install --force-reinstall --no-deps "xarray==2023.10.1"
  ```

  If `xarray-dataclass` / `xarray-schema` / `xarray-spatial` complain
  after the downgrade, uninstall the whole xarray-* family — nothing
  on the register+warp path calls into xarray directly:
  ```bash
  pip uninstall -y xarray xarray-dataclass xarray-schema xarray-spatial
  ```

- **`_ARRAY_API not found` / `numpy.core.multiarray failed to
  import`** — first `import fastcluster` or any downstream code path
  that pulls it (`from valis_hest import registration`,
  `import scanpy`, …). Root cause: your env's `numpy` has drifted
  to 2.x while `fastcluster`'s compiled `.so` was built against
  numpy 1.x — the classic numpy 1↔2 C-API ABI mismatch. The most
  common way an env drifts into this is a
  `pip install --force-reinstall <pkg>` where `<pkg>` has a loose
  `numpy>=1.22` bound and pip re-resolves numpy to the latest.
  Full diagnostic:
  [comment 5336001245](https://github.com/settylab/TracyY123-nexus/issues/15#issuecomment-5336001245).

  **Fix:**
  ```bash
  pip install --force-reinstall --no-deps "numpy==1.26.4"
  ```

  Verify the ABI is restored:
  ```bash
  python -c "import numpy; print(numpy.__version__)"           # expect 1.26.4
  python -c "import fastcluster; print(fastcluster.__version__)"  # expect 1.2.6
  python -c "from valis_hest import registration; print('OK')"    # expect OK
  ```

- **`AttributeError: module 'cv2.xfeatures2d' has no attribute
  'VGG_create'`** — first `from valis_hest import feature_detectors`
  (or any code path that pulls it — `from valis_hest import
  registration`, `import valis_hest`). Root cause: your env has
  more than one `opencv-*` PyPI distribution installed at the same
  time. All four (`opencv-python`, `opencv-python-headless`,
  `opencv-contrib-python`, `opencv-contrib-python-headless`)
  install into the SAME `cv2/` module directory and are meant to
  be mutually exclusive; when `pip install --no-deps` walks a
  requirements file that pins several of them, whichever variant
  lands LAST wins the shared files. If the stock `opencv-python`
  wins, its `cv2.xfeatures2d` is empty (non-free algorithms
  stripped) and `VGG_create` / `SIFT_create` / `BEBLID_create` all
  disappear. Full diagnostic:
  [comment 5349853726](https://github.com/settylab/TracyY123-nexus/issues/15#issuecomment-5349853726).

  **Fix:**
  ```bash
  pip uninstall -y opencv-python opencv-python-headless \
                   opencv-contrib-python opencv-contrib-python-headless
  pip install --no-deps "opencv-contrib-python==4.13.0.92"
  ```

  Verify:
  ```bash
  python -c "import cv2; print(cv2.__version__)"                    # expect 4.13.0
  python -c "import cv2; print(hasattr(cv2.xfeatures2d, 'VGG_create'))"  # expect True
  ```

  The current `heRegistration-requirements.txt` pins ONLY
  `opencv-contrib-python` (the other three were dropped in the
  commit that added this note) so a clean install from the
  requirements no longer triggers this bug — but an ad-hoc
  `pip install opencv-python` (or a transitive dep that specifies
  stock opencv without `--no-deps`) can still overwrite the
  contrib variant on disk.

### Verified min-blast-radius state

Tracy validated the following pin set end-to-end on Gizmo
(`heRegistration-test` env, 2026-08-18; see
[comment 5334721763](https://github.com/settylab/TracyY123-nexus/issues/15#issuecomment-5334721763)):

| Package | Version |
| --- | --- |
| `numpy` | `1.26.4` |
| `pandas` | `1.5.3` |
| `pyvips` | `2.2.3` |
| `xarray` | `2023.10.1` |
| `fastcluster` | `1.2.6` |
| `hest` | `1.1.1` (installed from the `v1.2.0` git tag) |
| `valis_hest` | `0.0.2` |

The current `heRegistration.yml` + `heRegistration-requirements.txt`
resolve to `numpy==1.26.4` + `xarray==2023.10.1` (the load-bearing
pins) and permit `pandas<3` / `pyvips>=3` / `fastcluster==1.3.0`
above them — the register+warp code paths run against both the
pandas-1/pyvips-2 set above and the pandas-2/pyvips-3 set that the
current pins permit. If you hit a bug we haven't seen and want to
narrow the surface, `pip install --force-reinstall --no-deps
"pyvips==2.2.3" "pandas==1.5.3" "fastcluster==1.2.6"` (all with
`--no-deps`) reproduces Tracy's verified state.

### Shims that are already permanent (nothing to do)

Two HEST v1.2.0 gaps landed as shims in `hexenium/_internal/`:

- `warp_and_save_xenium_objects` — HEST v1.2.0 doesn't export it,
  so `hexenium._internal.hest_warp_shim` composes v1.2.0's
  `warp_gdf_valis` primitive. Landed on commit `634158d`.
- `pixel_size_morph` — HEST's `read_gdf` needed it as a keyword
  argument; the shim plumbs it through. Landed on commit `60d97ff`.

You do not need to install anything extra for either. If you see
`ImportError: cannot import name 'warp_and_save_xenium_objects'` or
`TypeError: unsupported operand type(s) for /: 'float' and
'NoneType'` from HEST/VALIS, verify your checkout includes both
commits (`git log --oneline v0.2.0`).
