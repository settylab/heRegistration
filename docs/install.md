# Install

`hexenium` depends on VALIS + HEST, so we install the conda-forge /
bioconda stack first and then `pip install -e .` the package on top.

## Prerequisites

- **micromamba** (or mamba/conda). Fresh install:
  ```bash
  "${SHELL}" <(curl -L micro.mamba.pm/install.sh)
  ```
- A working checkout of hexenium.

## 1. Create the conda env

The `environments/heRegistration.yml` file lists the pinned
dependencies used to develop and validate this package.

```bash
cd /path/to/xenium-he-registration
micromamba env create -f environments/heRegistration.yml
micromamba activate heRegistration
```

## 2. Editable install

```bash
pip install -e .
```

This installs hexenium in `-e` (editable) mode: source edits under
`src/hexenium/` are picked up without a re-install, and the
console-script `hexenium` lands on your `PATH`.

## 3. HEST (optional editable clone)

The env file pulls HEST from GitHub via pip. If you want an editable
clone for local hacks:

```bash
git clone https://github.com/mahmoodlab/HEST.git ~/HEST
pip install -e ~/HEST
```

## 4. Verify

```bash
python -c "import hest, valis_hest, dask, openslide; print('OK')"
hexenium --version
hexenium run --help
```

If either of the last two commands errors, check that your active env
is `heRegistration` (not `base`) and that `pip install -e .` returned
successfully.
