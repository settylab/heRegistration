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
