#!/usr/bin/env bash
# scripts/lib/env_config.sh — load + validate scripts/env.local.conf.
#
# Meant to be SOURCED, not executed: `source ".../lib/env_config.sh"`.
# Both sbatch wrappers (submit_he_registration.sh,
# submit_vsi_to_ometiff.sh) and scripts/env-preflight.sh source this so
# there is exactly one place that knows the config's shape.
#
# Rationale for why this file exists at all: three real sbatch
# submissions of this pipeline hung before the wrapper's own script
# body logged anything (0-byte stdout, no output dir), while the SAME
# pipeline runs fine interactively outside the sandbox. The wrapper was
# built around `$HOME`: it sourced `$HOME/.bashrc`, defaulted
# `MAMBA_ROOT_PREFIX` to `$HOME/micromamba`, and searched
# `$HOME/.local/bin` / `$HOME/micromamba/bin`. Inside the sandbox
# `$HOME` is a non-persistent tmpfs, so those resolve differently
# inside the job than in the submitting shell. A file on shared
# storage, read for its literal (already-absolute) values, is immune to
# that class of bug; `sbatch --export=ALL` forwards env VAR values
# faithfully but does nothing to stop `$HOME` itself from resolving to
# a different instance once the job's own shell starts — so anything
# that re-derives a path from `$HOME` at job time is still exposed even
# with `--export=ALL`. Resolving once, at install time, and reading the
# resolved absolute values back is what actually fixes it.
#
# On success, exports:
#   MICROMAMBA_BIN     absolute path to the micromamba binary (required)
#   HEREG_ENV_PREFIX   absolute path to the micromamba env PREFIX, not
#                       a name (required)
#   BFTOOLS_ROOT        absolute path to the extracted Glencoe zips
#                       (bioformats2raw-*/bin, raw2ometiff-*/bin) — only
#                       consulted by the VSI-conversion path
#                       (submit_vsi_to_ometiff.sh) when those tools
#                       aren't already on PATH in the activated env.
#                       OPTIONAL: leave unset in env.local.conf if VSI
#                       input isn't used.
#   LIBBLOSC_DIR        absolute path to a directory holding
#                       libblosc.so* — same optionality as BFTOOLS_ROOT.
#
# Fails loud (exit 1) with a "run the install step" message if:
#   - scripts/env.local.conf is missing or unreadable
#   - MICROMAMBA_BIN or HEREG_ENV_PREFIX is missing or empty
#   - MICROMAMBA_BIN does not point at an executable file
#   - HEREG_ENV_PREFIX does not point at an existing directory
#   - BFTOOLS_ROOT or LIBBLOSC_DIR is SET in the config but empty (a
#     stale/corrupted entry) — but they may be entirely absent, since
#     they are optional.

_env_config_fail() {
    echo "error: $1" >&2
    echo "       Run scripts/write-env-config.sh to (re)generate "\
"scripts/env.local.conf — see docs/install.md." >&2
    # Unconditional exit, not `return`: `return 1` only aborts the sourcing
    # script if that script itself runs under `set -e` -- a caller contract,
    # not a property of this file. Every real call site (the two sbatch
    # wrappers, env-preflight.sh) does have `set -euo pipefail` before
    # sourcing this, but this file's whole purpose is to fail loud
    # unconditionally, so it shouldn't depend on that. This is never
    # meant to be sourced from an interactive login shell, so `exit`
    # (not `return`) is safe here.
    exit 1
}

_ENV_CONFIG_LIB_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# HEREG_ENV_LOCAL_CONF overrides the default location — used by tests to
# point at a fixture conf instead of the real scripts/env.local.conf;
# production callers should never need to set it.
ENV_LOCAL_CONF="${HEREG_ENV_LOCAL_CONF:-$_ENV_CONFIG_LIB_DIR/../env.local.conf}"

if [[ ! -f "$ENV_LOCAL_CONF" ]]; then
    _env_config_fail "scripts/env.local.conf not found (looked at $ENV_LOCAL_CONF)."
fi

if [[ ! -r "$ENV_LOCAL_CONF" ]]; then
    _env_config_fail "scripts/env.local.conf is not readable: $ENV_LOCAL_CONF"
fi

# shellcheck source=/dev/null
source "$ENV_LOCAL_CONF"

for _var in MICROMAMBA_BIN HEREG_ENV_PREFIX; do
    if [[ -z "${!_var:-}" ]]; then
        _env_config_fail "scripts/env.local.conf is missing or has an empty '$_var'."
    fi
done
unset _var

if [[ ! -x "$MICROMAMBA_BIN" ]]; then
    _env_config_fail "MICROMAMBA_BIN in scripts/env.local.conf is not an executable file: $MICROMAMBA_BIN"
fi

if [[ ! -d "$HEREG_ENV_PREFIX" ]]; then
    _env_config_fail "HEREG_ENV_PREFIX in scripts/env.local.conf does not exist: $HEREG_ENV_PREFIX"
fi

# BFTOOLS_ROOT / LIBBLOSC_DIR are OPTIONAL — only the VSI-conversion
# path needs them, and only as a fallback when bioformats2raw /
# raw2ometiff / libblosc aren't already reachable in the activated env.
# A config without either stays valid (the common non-VSI case). But if
# a key IS present, an empty value means the config is stale/corrupted
# (write-env-config.sh never writes an empty value for a key it emits),
# so that still fails loud — `-v` (not `-z`) is what tells "present but
# empty" apart from "absent".
for _var in BFTOOLS_ROOT LIBBLOSC_DIR; do
    if [[ -v "$_var" && -z "${!_var}" ]]; then
        _env_config_fail "scripts/env.local.conf sets '$_var' but leaves it empty; remove the line if VSI conversion isn't needed."
    fi
done
unset _var

if [[ -n "${BFTOOLS_ROOT:-}" && ! -d "$BFTOOLS_ROOT" ]]; then
    _env_config_fail "BFTOOLS_ROOT in scripts/env.local.conf does not exist: $BFTOOLS_ROOT"
fi

if [[ -n "${LIBBLOSC_DIR:-}" && ! -d "$LIBBLOSC_DIR" ]]; then
    _env_config_fail "LIBBLOSC_DIR in scripts/env.local.conf does not exist: $LIBBLOSC_DIR"
fi

export MICROMAMBA_BIN HEREG_ENV_PREFIX
# NOT `[[ -n ... ]] && export ...`: when sourced, this file's exit status
# is whatever its LAST executed command returned. A false `[[ ]]` as (or
# driving) the last statement makes `source lib/env_config.sh` itself
# return 1 — which silently aborts the CALLING script under `set -e`,
# even on the fully-valid "optional key absent" path, with no error
# message at all (caught by actually sourcing this file and checking
# $?, not just by reading the code). `if` avoids this: bash defines an
# `if` whose condition is false and has no matching commands executed
# as exiting 0, not the condition's own (nonzero) status.
if [[ -n "${BFTOOLS_ROOT:-}" ]]; then
    export BFTOOLS_ROOT
fi
if [[ -n "${LIBBLOSC_DIR:-}" ]]; then
    export LIBBLOSC_DIR
fi
