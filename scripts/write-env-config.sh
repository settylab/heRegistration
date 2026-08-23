#!/usr/bin/env bash
# scripts/write-env-config.sh — the final step of installation: records
# RESOLVED ABSOLUTE values into the gitignored scripts/env.local.conf,
# so every sbatch job reads a fixed, already-resolved config instead of
# re-deriving `$HOME`-relative paths inside its own (possibly different)
# environment instance. See scripts/lib/env_config.sh for the consumer
# side and the full rationale.
#
# Usage:
#   scripts/write-env-config.sh \
#       --env-prefix        /abs/path/to/micromamba/envs/heRegistration \
#       [--micromamba-bin   /abs/path/to/micromamba] \
#       [--mamba-root-prefix /abs/path/to/micromamba/root] \
#       [--bftools-root     /abs/path/to/extracted/glencoe/zips] \
#       [--libblosc-dir     /abs/path/to/env/with/blosc/lib]
#
# --env-prefix must already exist (this script records, it does not
# create, the environment — run `scripts/create-env.sh env create
# -n heRegistration -f environments/heRegistration.yml` first, then
# resolve its prefix via `micromamba env list`). --micromamba-bin auto-detects via
# `command -v micromamba` if not passed. --mamba-root-prefix is DERIVED
# from --env-prefix if not passed: a micromamba env prefix is, by
# construction, `<root>/envs/<name>`, so the root is
# `${ENV_PREFIX%/envs/*}`. This is NOT --env-prefix: it's micromamba's
# own root (holds condabin/, pkgs/, etc.) — required because
# micromamba's shell hook script references `${MAMBA_ROOT_PREFIX}` with
# no fallback the moment it runs in a shell that doesn't already have
# CONDA_SHLVL set, which is every fresh Slurm job shell; see
# scripts/lib/env_config.sh for the reproduction.
#
# Deliberately NOT auto-detected from the current shell's
# $MAMBA_ROOT_PREFIX: every `sbatch` job sources the resulting
# env.local.conf, so a wrong value here contaminates the real compute
# job, not just this recording step. If the ambient shell has
# $MAMBA_ROOT_PREFIX set and it disagrees with the value derived from
# --env-prefix, this script refuses to guess and exits loud — pass
# --mamba-root-prefix explicitly to pick one (e.g. for a nonstandard
# layout where --env-prefix isn't under `<root>/envs/`).
#
# --bftools-root / --libblosc-dir are OPTIONAL — only needed if you plan
# to convert .vsi inputs via scripts/submit_vsi_to_ometiff.sh and
# bioformats2raw/raw2ometiff/libblosc aren't already installed inside
# --env-prefix (see docs/install.md "VSI-input prerequisites").

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUT="$SCRIPT_DIR/env.local.conf"

ENV_PREFIX=""
MICROMAMBA_BIN=""
MAMBA_ROOT_PREFIX_ARG=""
BFTOOLS_ROOT=""
LIBBLOSC_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --env-prefix)         ENV_PREFIX="$2"; shift 2 ;;
        --micromamba-bin)     MICROMAMBA_BIN="$2"; shift 2 ;;
        --mamba-root-prefix)  MAMBA_ROOT_PREFIX_ARG="$2"; shift 2 ;;
        --bftools-root)       BFTOOLS_ROOT="$2"; shift 2 ;;
        --libblosc-dir)       LIBBLOSC_DIR="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,34p' "$0"; exit 0 ;;
        *)
            echo "error: unknown argument: $1" >&2; exit 2 ;;
    esac
done

: "${ENV_PREFIX:?--env-prefix is required (absolute path to the micromamba env PREFIX)}"

if [[ ! -d "$ENV_PREFIX" ]]; then
    echo "error: --env-prefix does not exist: $ENV_PREFIX" >&2
    echo "       Create it first: scripts/create-env.sh env create -n heRegistration -f environments/heRegistration.yml" >&2
    echo "       then resolve its prefix with: micromamba env list" >&2
    exit 3
fi
ENV_PREFIX=$(cd "$ENV_PREFIX" && pwd)   # canonicalize to an absolute path

if [[ -z "$MICROMAMBA_BIN" ]]; then
    if ! MICROMAMBA_BIN=$(command -v micromamba); then
        echo "error: micromamba not found on PATH and --micromamba-bin not given." >&2
        exit 4
    fi
fi
if [[ ! -x "$MICROMAMBA_BIN" ]]; then
    echo "error: --micromamba-bin is not executable: $MICROMAMBA_BIN" >&2
    exit 4
fi
# Canonicalize (command -v can return a relative-looking path on some shells).
MICROMAMBA_BIN=$(cd "$(dirname "$MICROMAMBA_BIN")" && pwd)/$(basename "$MICROMAMBA_BIN")

if [[ -z "$MAMBA_ROOT_PREFIX_ARG" ]]; then
    # Derive from --env-prefix (already canonicalized above), not from
    # the ambient shell: a micromamba env prefix is, by construction,
    # <root>/envs/<name>. See the header comment for why this must not
    # fall back to reading $MAMBA_ROOT_PREFIX out of the calling shell.
    DERIVED_ROOT_PREFIX="${ENV_PREFIX%/envs/*}"
    if [[ "$DERIVED_ROOT_PREFIX" == "$ENV_PREFIX" ]]; then
        echo "error: could not derive --mamba-root-prefix from --env-prefix ($ENV_PREFIX):" >&2
        echo "       it does not look like <root>/envs/<name>." >&2
        echo "       Pass --mamba-root-prefix explicitly for this nonstandard layout." >&2
        exit 5
    fi
    if [[ -n "${MAMBA_ROOT_PREFIX:-}" && "${MAMBA_ROOT_PREFIX}" != "$DERIVED_ROOT_PREFIX" ]]; then
        echo "error: this shell's ambient \$MAMBA_ROOT_PREFIX ($MAMBA_ROOT_PREFIX) disagrees" >&2
        echo "       with the root derived from --env-prefix ($DERIVED_ROOT_PREFIX)." >&2
        echo "       Refusing to silently prefer either — pass --mamba-root-prefix" >&2
        echo "       explicitly to pick one, or unset \$MAMBA_ROOT_PREFIX if the derived" >&2
        echo "       value is correct." >&2
        exit 5
    fi
    MAMBA_ROOT_PREFIX_ARG="$DERIVED_ROOT_PREFIX"
fi
if [[ ! -d "$MAMBA_ROOT_PREFIX_ARG" ]]; then
    echo "error: --mamba-root-prefix does not exist: $MAMBA_ROOT_PREFIX_ARG" >&2
    exit 3
fi
MAMBA_ROOT_PREFIX_ARG=$(cd "$MAMBA_ROOT_PREFIX_ARG" && pwd)

if [[ -n "$BFTOOLS_ROOT" ]]; then
    if [[ ! -d "$BFTOOLS_ROOT" ]]; then
        echo "error: --bftools-root does not exist: $BFTOOLS_ROOT" >&2
        exit 3
    fi
    BFTOOLS_ROOT=$(cd "$BFTOOLS_ROOT" && pwd)
fi

if [[ -n "$LIBBLOSC_DIR" ]]; then
    if [[ ! -d "$LIBBLOSC_DIR" ]]; then
        echo "error: --libblosc-dir does not exist: $LIBBLOSC_DIR" >&2
        exit 3
    fi
    LIBBLOSC_DIR=$(cd "$LIBBLOSC_DIR" && pwd)
fi

{
    echo "# Generated by scripts/write-env-config.sh on $(date -u +%Y-%m-%dT%H:%M:%SZ)."
    echo "# Machine-local, gitignored — resolved ABSOLUTE values for THIS install."
    echo "# Do not hand-edit; re-run scripts/write-env-config.sh instead."
    echo "MICROMAMBA_BIN=$MICROMAMBA_BIN"
    echo "HEREG_ENV_PREFIX=$ENV_PREFIX"
    echo "MAMBA_ROOT_PREFIX=$MAMBA_ROOT_PREFIX_ARG"
    [[ -n "$BFTOOLS_ROOT" ]] && echo "BFTOOLS_ROOT=$BFTOOLS_ROOT"
    [[ -n "$LIBBLOSC_DIR" ]] && echo "LIBBLOSC_DIR=$LIBBLOSC_DIR"
} > "$OUT"

echo "[write-env-config] wrote $OUT:"
sed 's/^/[write-env-config]   /' "$OUT"
echo "[write-env-config] verify with: scripts/env-preflight.sh"
