# Install

`hexenium` depends on VALIS + HEST + their Java/BioFormats bridge, so
reproducing the working env is a **two-step install**: a conda-side
solve for the system libraries and scientific-Python core, then a
`--no-deps` pip layer for `hest` / `valis-wsi` / `valis-hest` and
their transitive stack. Both steps are captured in
`environments/heRegistration.yml` and
`environments/heRegistration-requirements.txt`.

`--no-deps` on the pip step is **essential** and must live on the
CLI — see [Why `--no-deps`](#why---no-deps-is-mandatory) below. This
mirrors the recommended flow in the top-level README; keep both in
sync.

## Prerequisites

- **micromamba** (or mamba/conda). Fresh install:
  ```bash
  "${SHELL}" <(curl -L micro.mamba.pm/install.sh)
  ```
  **`micromamba activate` (step 3 below) requires the shell hook to
  already be sourced** — a fresh install of the raw `micromamba` binary
  does not wire this up by itself. If `micromamba activate
  heRegistration` fails with `critical libmamba Shell not initialized`
  / `'micromamba' is running as a subprocess and can't modify the
  parent shell`, run this once (add it to your shell rc to persist
  across sessions):
  ```bash
  eval "$(micromamba shell hook --shell bash)"   # or --shell zsh
  ```
  The install script above offers to do this for you interactively;
  if you accepted a non-interactive/scripted install, or installed
  micromamba some other way, it's easy to end up without the hook.
- A working checkout of hexenium.

## 1. Clone the repo

```bash
git clone https://github.com/settylab/heRegistration.git
cd heRegistration
```

## 2. Create the conda env

The yml intentionally omits a `name:` field; pass one with `-n`.
`heRegistration` is just a readable convention — the sbatch wrapper no
longer selects an env by name (see step 8): whatever you call it here,
you point the wrapper at its resolved prefix explicitly later.

```bash
scripts/create-env.sh env create -n heRegistration -f environments/heRegistration.yml
# or with conda:
# conda env create -n heRegistration -f environments/heRegistration.yml
```

`scripts/create-env.sh` is a thin wrapper that runs any `micromamba`
env-creation subcommand you pass it (`env create`, `create`, …). If
you don't set `MAMBA_ROOT_PREFIX`, it's a plain passthrough. If you
**do** set `MAMBA_ROOT_PREFIX` to keep this install fully isolated
(off `$HOME`, for a scratch/test install, or to keep multiple installs
from sharing state) — set these first:

```bash
export MAMBA_ROOT_PREFIX=/abs/path/to/isolated/root
export XDG_CACHE_HOME="$MAMBA_ROOT_PREFIX/xdg-cache"
export XDG_CONFIG_HOME="$MAMBA_ROOT_PREFIX/xdg-config"
```

The wrapper also exports `CONDA_PKGS_DIRS="$MAMBA_ROOT_PREFIX/pkgs"` and
asserts afterward that `~/.mamba/pkgs` was not touched. Without this,
micromamba's `pkgs_dirs` silently resolves to
`[$MAMBA_ROOT_PREFIX/pkgs, ~/.mamba/pkgs]` — an undocumented second
entry — so an "isolated" install can still quietly write package-cache
state to `~/.mamba/pkgs` with nothing erroring. A bare
`micromamba env create` still works exactly as before if you don't
need isolation. **Use `scripts/create-env.sh` (not the bare
`micromamba`/`conda` command) for every env this repo has you create**
— see also the ad-hoc `blosc` env in [VSI-input
prerequisites](#vsi-input-prerequisites-only-if---he-slide-is-a-vsi-file)
below, which uses the same wrapper.

`scripts/pip-install.sh` (step 4) covers pip's own HTTP/wheel cache,
but several importable dependencies fall back to the [XDG Base
Directory](https://specifications.freedesktop.org/basedir-spec/latest/)
spec (`$XDG_CACHE_HOME`, default `~/.cache`; `$XDG_CONFIG_HOME`,
default `~/.config`) for their OWN caches, with no repo-specific
wrapper to intercept them — e.g. `matplotlib` writes a font-list cache
(`fontlist-*.json`) and a config dir the first time it's imported,
confirmed to fire during step 7's verify import (`import hest,
valis_hest, ...`). The two `export`s above redirect this whole class
of dependency at the source instead of chasing each library that
writes to `$HOME` one at a time. Leave them unset if you didn't set
`MAMBA_ROOT_PREFIX` either — everything then falls back to the normal
`$HOME` locations.

## 3. Activate

```bash
micromamba activate heRegistration
# or: conda activate heRegistration
```

## 4. Install the pinned pip layer with `--no-deps`

```bash
scripts/pip-install.sh --no-deps -r environments/heRegistration-requirements.txt
```

`scripts/pip-install.sh` is a thin wrapper around `pip install`. Like
`scripts/create-env.sh` (step 2), it's a plain passthrough unless
`MAMBA_ROOT_PREFIX` is set — but if it IS set, plain `pip install`
still writes its HTTP + wheel cache to `~/.cache/pip` regardless
(there's no conda-side `CONDA_PKGS_DIRS` equivalent pip honors by
default), so this step alone writes 170+ files there with no way to
opt out. The wrapper overrides `PIP_CACHE_DIR` into
`$MAMBA_ROOT_PREFIX/pip-cache` and asserts afterward that
`~/.cache/pip` was not touched, the same isolation guarantee
`create-env.sh` gives the conda side. A bare `pip install --no-deps
-r ...` still works exactly as before if you don't need isolation.

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
scripts/pip-install.sh --no-deps .
```

For developers (editable):

```bash
scripts/pip-install.sh --no-deps -e .
```

`--no-deps` here is **essential**: without it, pip re-resolves
`pyproject.toml`'s declared deps against PyPI, which can
override the pinned layer from step 4. With `--no-deps`, pip
installs just the `hexenium` package + entry points and
leaves your step-4 pins alone.

Non-editable is recommended for end users: an editable install
exposes the source tree to `sys.path`, so a stray import via a
working-directory Python — or a sibling `hexenium/` folder in
`cwd` — can shadow the installed package and silently pull in
half-updated modules.

Both commands install just the `hexenium` package on top of the
env you already prepared in steps 2-4. Notably, `pyproject.toml`
declares
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
scripts/pip-install.sh --no-deps -e ~/HEST
```

## 7. Verify

```bash
python -c "import hest, valis_hest, valis_hest.registration, valis_hest.slide_io, dask, openslide; print('OK')"
hexenium --version
hexenium run --help
```

If any of these errors, check that the active env is the one you
created (not `base`) and see the Troubleshooting section below.

## 8. Configure Slurm env activation (required before any sbatch submission)

Steps 1-7 get `hexenium` running **interactively**. A separate,
one-time step is required before `scripts/submit_he_registration.sh` or
`scripts/submit_vsi_to_ometiff.sh` can run under Slurm at all.

Both sbatch wrappers activate their conda env by an **absolute,
pre-resolved prefix** — never by name, never by re-deriving anything
from `$HOME` at job time. (An earlier version sourced `~/.bashrc` and
searched `$HOME`-relative paths; that caused real jobs to hang with
0-byte output because the job's `$HOME` didn't match the submitting
shell's. See `scripts/lib/env_config.sh`'s header comment for the full
rationale.) The resolved values live in `scripts/env.local.conf` — a
gitignored, machine-local file — and every sbatch job reads it back at
job start via `scripts/lib/env_config.sh`. It does not exist yet after
a fresh clone, so this step is mandatory, not optional.

Generate it with `scripts/write-env-config.sh`:

```bash
scripts/write-env-config.sh --env-name heRegistration
```

`--env-name heRegistration` reads back the prefix that
`scripts/create-env.sh` (step 2) already recorded to
`scripts/.env-prefix-heRegistration` at env-creation time — it does
not query `micromamba env list`. **Do not use `micromamba env list` to
find this path by hand.** On any account with install history, that
list shows one row per `heRegistration`-named env ever created, across
every root ever used — only the CURRENTLY active root's rows are
name-labeled; every other root's matching row shows a blank Name
column, so a quick visual scan can easily copy a stale, foreign env's
path into this command. If you skipped `scripts/create-env.sh` and
created the env some other way (a bare `micromamba env create` /
`conda env create`), there is no receipt to read back — resolve the
prefix yourself (e.g. `conda info --envs`, cross-checking against the
root you actually just created into) and pass `--env-prefix
/abs/path/to/the/env` directly instead:

```bash
scripts/write-env-config.sh --env-prefix /home/you/micromamba/envs/heRegistration
```

`--micromamba-bin` auto-detects via `command -v micromamba` if
omitted. `--mamba-root-prefix` is **derived from `--env-prefix`** if
omitted — a micromamba env prefix is, by construction,
`<root>/envs/<name>`, so the root is inferred from the prefix you just
passed. It is deliberately **not** read from the current shell's
`$MAMBA_ROOT_PREFIX`: every `sbatch` job sources the resulting
`env.local.conf`, so a wrong value here would contaminate the real
compute job, not just this recording step. If your shell happens to
have `$MAMBA_ROOT_PREFIX` set and it disagrees with the derived value,
the script refuses to guess and exits with an error — pass
`--mamba-root-prefix` explicitly to pick one (needed only for a
nonstandard layout where `--env-prefix` isn't under `<root>/envs/`).
If you set up VSI support in step 7,
also pass `--bftools-root`/`--libblosc-dir` here — this is now the
canonical way to record those two paths (they used to be documented as
plain `export`s only; see [VSI-input
prerequisites](#vsi-input-prerequisites-only-if---he-slide-is-a-vsi-file)
below):

```bash
scripts/write-env-config.sh \
    --env-name      heRegistration \
    --bftools-root  $HOME/opt/bftools \
    --libblosc-dir  $HOME/micromamba/envs/blosc/lib
```

These are the literal paths from the [VSI-install
walkthrough](#vsi-input-prerequisites-only-if---he-slide-is-a-vsi-file)
below, not `$BFTOOLS_ROOT`/`$LIBBLOSC_DIR` read back from your shell —
don't substitute those two variable names here. They only exist in the
shell session that ran that walkthrough's own `export` lines, and are
unset in any other shell (a fresh terminal, a new SSH session, or if
you skip VSI setup entirely). An unset, unquoted `$VAR` vanishes from
the command line rather than expanding to an empty string, silently
shifting every argument after it into the wrong flag — `write-env-config.sh`
then reports a confusing "does not exist" error for a path you never
meant to pass.

Then verify it resolves end-to-end (activates the env, confirms
`hexenium` is installed, and — unless `--skip-vsi-check` — confirms
`bioformats2raw`/`raw2ometiff`/`libblosc` are reachable):

```bash
scripts/env-preflight.sh
```

**Re-run `write-env-config.sh` whenever the env moves or is recreated**
— `env.local.conf` pins the prefix that existed at the time you ran it;
nothing keeps it in sync automatically.

**Why this matters more than a normal install step:** `submit_he_registration.sh`'s
launcher branch (what runs when you invoke it directly — the part that
parses flags and calls `sbatch`) never sources `env_config.sh` at all;
that only happens later, under `sbatch`, once the job actually starts
on a compute node (`scripts/submit_he_registration.sh:388-390`). So
skipping this step does **not** make submission fail — you'll see a
normal `Submitted batch job NNNN` and the command exits 0, and the
`scripts/env.local.conf`-missing error only surfaces afterward, inside
that job's own `slurm-<jobid>.err`. If you're about to queue a run and
step away, run `scripts/env-preflight.sh` first so a bad config fails
in your terminal, not silently in a log file you have to go find. See
the README's ["Where Slurm logs
land"](../README.md#where-slurm-logs-land) for the exact log path per
invocation mode.

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
- **Interactive `micromamba activate heRegistration` fails** with
  "environment not found". You created the env under a different name
  in step 2 — activate it under that name instead, or recreate it as
  `heRegistration`.
- **Sbatch job fails** (check `slurm-<jobid>.err` under the run's logs
  dir — see step 8 for why this doesn't show up on your terminal) with
  `error: scripts/env.local.conf not found` or `MICROMAMBA_BIN ... is
  not an executable file` or similar. You haven't done step 8 yet, or
  `env.local.conf` still points at a prefix that no longer exists (env
  recreated or moved). Re-run `scripts/write-env-config.sh --env-name
  heRegistration` (or `--env-prefix <path>`) and verify with
  `scripts/env-preflight.sh`.
  **`--env-name <name>` / `ENV_NAME=<name>` passed to
  `submit_he_registration.sh` / `submit_vsi_to_ometiff.sh` do NOT fix
  this** — that's a different, deprecated flag on the sbatch launcher
  scripts (kept only so old invocations don't hit an "unknown flag"
  error; it only emits a WARN). Activation there is always by the
  absolute prefix pinned in `scripts/env.local.conf`. This is unrelated
  to `write-env-config.sh --env-name` above, which does something real:
  it resolves that absolute prefix from `scripts/create-env.sh`'s
  receipt file.

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

### The `--no-deps` rule

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
  actually needs xarray.

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
  disappear.

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

- **`ValueError: NoneType copy mode not allowed`** deep in
  VALIS's rigid-registration step, from
  `fastcluster.py:linkage`. Root cause: your env's
  `fastcluster` has drifted UP to 1.3.0 while `numpy` is
  correctly pinned at `1.26.4`. fastcluster 1.3.0's
  `linkage()` passes `copy=None` to numpy's `array()` for
  `method='single'` (the method valis_hest's
  `serial_rigid.order_Dmat` uses); numpy 2.x treats
  `copy=None` as "copy if needed", but numpy 1.x rejects it
  outright.

  **Fix:**
  ```bash
  pip install --force-reinstall --no-deps "fastcluster==1.2.6"
  ```

  Verify:
  ```bash
  python -c "import fastcluster; print(fastcluster.__version__)"  # expect 1.2.6
  python -c "from valis_hest import serial_rigid; print('OK')"    # expect OK
  ```

  The current `heRegistration-requirements.txt` pins
  `fastcluster==1.2.6` explicitly (with a comment explaining
  the numpy-1.x compat requirement), so a clean install
  from the requirements no longer triggers this — but a
  `pip install --upgrade fastcluster` or a re-freeze against
  a numpy-2.x-only wheel can drift the env back into the
  bug.

### VSI-input prerequisites (only if `--he-slide` is a `.vsi` file)

`submit_he_registration.sh` dispatches a
`bioformats2raw` + `raw2ometiff` conversion job before
running hexenium when the input is a `.vsi` file.

The two conversion tools are not included in the base
`heRegistration` environment. The tested / recommended
install path is the Glencoe binary zips; the conda-based
install did not resolve reproducibly in our tested
environment, so we don't ship it as the default.

The launcher discovers the tools through two env vars:

- `BFTOOLS_ROOT` — a directory containing the extracted
  Glencoe zips (i.e.
  `$BFTOOLS_ROOT/bioformats2raw-*/bin/` and
  `$BFTOOLS_ROOT/raw2ometiff-*/bin/` exist).
- `LIBBLOSC_DIR` — a directory containing
  `libblosc.so.1`. Any env's `lib/` that carries
  `conda-forge::blosc` works (for example a sibling env
  you already have around); the `heRegistration` env's
  imagecodecs-vendored copy is not on the discoverable
  path.

**Canonical way to set both: `scripts/write-env-config.sh
--bftools-root <dir> --libblosc-dir <dir>`** (step 8 above) — this
records the resolved absolute paths into `scripts/env.local.conf`,
which both sbatch wrappers read at job start. A plain ambient `export`
(below) still works as a fallback since `env_config.sh` never unsets a
pre-existing value, but it isn't the recommended path anymore: the
`export`-only approach is exactly the kind of `$HOME`/shell-session-
relative state this fix moved away from, even though these two specific
vars are absolute paths and so aren't unsafe the way the old
`$HOME`-relative micromamba lookup was.

Set both before running `submit_he_registration.sh` if you're using the
ambient-export fallback. The launcher (`scripts/submit_vsi_to_ometiff.sh`)
prepends `$BFTOOLS_ROOT/…/bin/` onto `PATH` and `$LIBBLOSC_DIR`
onto `LD_LIBRARY_PATH` automatically — you don't need to
manage those two variables yourself.

**Install (once per host)** — pick any writable
directory (location is your choice; the pipeline
discovers the tools via `BFTOOLS_ROOT`, not a fixed
path). Below uses `$HOME/opt/bftools/`:

```bash
# 1. Download the Glencoe release zips
#    (tested versions: bioformats2raw 0.12.1 + raw2ometiff 0.9.0):
mkdir -p $HOME/opt/bftools
cd $HOME/opt/bftools
wget https://github.com/glencoesoftware/bioformats2raw/releases/download/v0.12.1/bioformats2raw-0.12.1.zip
wget https://github.com/glencoesoftware/raw2ometiff/releases/download/v0.9.0/raw2ometiff-0.9.0.zip
unzip bioformats2raw-0.12.1.zip
unzip raw2ometiff-0.9.0.zip

# 2. Point the launcher at the extracted zips + a libblosc:
export BFTOOLS_ROOT=$HOME/opt/bftools
export LIBBLOSC_DIR=$HOME/micromamba/envs/<any-env-with-blosc>/lib
```

Release pages (newer versions also work — the launcher
globs `bioformats2raw-*/bin` and `raw2ometiff-*/bin`, so
version isn't hard-coded):

- https://github.com/glencoesoftware/bioformats2raw/releases
- https://github.com/glencoesoftware/raw2ometiff/releases

Any environment with `conda-forge::blosc` installed
provides the right `libblosc.so.1` ABI. If you don't
already have one, create a minimal env just for the
library:

```bash
scripts/create-env.sh create -n blosc -c conda-forge blosc -y
export LIBBLOSC_DIR=$HOME/micromamba/envs/blosc/lib
```

See the docstring at the top of
`scripts/submit_vsi_to_ometiff.sh` for the full env-var
contract.

Full workflow details: README ["HPC usage (Slurm) — VSI
inputs: automatic BioFormats conversion
(unified)"](../README.md#vsi-inputs-automatic-bioformats-conversion-unified).

Skip this section entirely if your H&E is already
`.ome.tif` / `.ome.tiff` / `.tif`.

### Verified min-blast-radius state

The following pin set has been validated end-to-end:

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
resolve to `numpy==1.26.4` + `xarray==2023.10.1` (the essential
pins) and permit `pandas<3` / `pyvips>=3` / `fastcluster==1.3.0`
above them — the register+warp code paths run against both the
pandas-1/pyvips-2 set above and the pandas-2/pyvips-3 set that the
current pins permit. If you hit a bug we haven't seen and want to
narrow the surface, `pip install --force-reinstall --no-deps
"pyvips==2.2.3" "pandas==1.5.3" "fastcluster==1.2.6"` (all with
`--no-deps`) reproduces the previously-verified pandas-1/pyvips-2 state.

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
