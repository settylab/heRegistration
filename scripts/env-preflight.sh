#!/usr/bin/env bash
# scripts/env-preflight.sh — submit-time preflight for the resolved
# environment. Run this ONCE, before dispatching submit_he_registration.sh
# / submit_vsi_to_ometiff.sh, so a broken/stale environment fails in
# seconds instead of after a queue wait — or worse, after the job's own
# script body never even logs anything (the failure mode this whole
# env-config mechanism exists to prevent; see scripts/lib/env_config.sh).
#
# Checks:
#   1. scripts/env.local.conf exists and every required key resolves to
#      a real path (delegated to scripts/lib/env_config.sh — the same
#      check both sbatch wrappers do at job start, just run here first).
#   2. The env at HEREG_ENV_PREFIX actually has `hexenium` installed —
#      the expensive check (a real `micromamba activate` + `command -v`),
#      which is why this is a separate, explicitly-invoked step rather
#      than something every sbatch job repeats redundantly.
#   3. If BFTOOLS_ROOT/LIBBLOSC_DIR are configured (VSI support), that
#      bioformats2raw + raw2ometiff actually resolve from them (or from
#      PATH inside the activated env).
#
# Usage: scripts/env-preflight.sh [--skip-vsi-check]
#   --skip-vsi-check   Only run the config-file + hexenium checks; skip
#                       the bioformats2raw/raw2ometiff resolution check
#                       (useful if you never pass a .vsi --he-slide).

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

SKIP_VSI_CHECK=0
if [[ "${1:-}" == "--skip-vsi-check" ]]; then
    SKIP_VSI_CHECK=1
fi

# shellcheck source=lib/env_config.sh
source "$SCRIPT_DIR/lib/env_config.sh"
echo "[env-preflight] scripts/env.local.conf OK:"
echo "[env-preflight]   MICROMAMBA_BIN=$MICROMAMBA_BIN"
echo "[env-preflight]   HEREG_ENV_PREFIX=$HEREG_ENV_PREFIX"
echo "[env-preflight]   BFTOOLS_ROOT=${BFTOOLS_ROOT:-<unset>}"
echo "[env-preflight]   LIBBLOSC_DIR=${LIBBLOSC_DIR:-<unset>}"

eval "$("$MICROMAMBA_BIN" shell hook --shell bash)"
if ! micromamba activate "$HEREG_ENV_PREFIX"; then
    echo "error: [env-preflight] micromamba activate failed for prefix: $HEREG_ENV_PREFIX" >&2
    exit 5
fi

if ! command -v hexenium >/dev/null 2>&1; then
    echo "error: [env-preflight] 'hexenium' is not installed in env: $HEREG_ENV_PREFIX" >&2
    echo "       Install it once: micromamba activate $HEREG_ENV_PREFIX && pip install -e ." >&2
    exit 6
fi
echo "[env-preflight] hexenium -> $(command -v hexenium)"

if [[ "$SKIP_VSI_CHECK" == "1" ]]; then
    echo "[env-preflight] --skip-vsi-check: not checking bioformats2raw/raw2ometiff resolution."
    echo "[env-preflight] all checks passed."
    exit 0
fi

if [[ -z "${BFTOOLS_ROOT:-}" && -z "${LIBBLOSC_DIR:-}" ]]; then
    echo "[env-preflight] BFTOOLS_ROOT/LIBBLOSC_DIR not configured — VSI (.vsi) input"
    echo "[env-preflight]   is only supported if bioformats2raw/raw2ometiff/libblosc"
    echo "[env-preflight]   are already on PATH in the env. Skipping the resolution check."
    echo "[env-preflight] all checks passed."
    exit 0
fi

B2R="$(command -v bioformats2raw || true)"
R2O="$(command -v raw2ometiff || true)"
if [[ -z "$B2R" && -n "${BFTOOLS_ROOT:-}" ]]; then
    # shellcheck disable=SC2010
    B2R_DIR="$(ls -d "$BFTOOLS_ROOT"/bioformats2raw-*/bin 2>/dev/null | head -1)"
    [[ -n "$B2R_DIR" && -x "$B2R_DIR/bioformats2raw" ]] && B2R="$B2R_DIR/bioformats2raw"
fi
if [[ -z "$R2O" && -n "${BFTOOLS_ROOT:-}" ]]; then
    # shellcheck disable=SC2010
    R2O_DIR="$(ls -d "$BFTOOLS_ROOT"/raw2ometiff-*/bin 2>/dev/null | head -1)"
    [[ -n "$R2O_DIR" && -x "$R2O_DIR/raw2ometiff" ]] && R2O="$R2O_DIR/raw2ometiff"
fi

if [[ -z "$B2R" ]]; then
    echo "error: [env-preflight] bioformats2raw not found on PATH or under BFTOOLS_ROOT: ${BFTOOLS_ROOT:-<unset>}" >&2
    exit 7
fi
if [[ -z "$R2O" ]]; then
    echo "error: [env-preflight] raw2ometiff not found on PATH or under BFTOOLS_ROOT: ${BFTOOLS_ROOT:-<unset>}" >&2
    exit 7
fi
echo "[env-preflight] bioformats2raw -> $B2R"
echo "[env-preflight] raw2ometiff    -> $R2O"

if [[ -n "${LIBBLOSC_DIR:-}" ]]; then
    if ! find "$LIBBLOSC_DIR" -maxdepth 1 -name 'libblosc.so*' 2>/dev/null | grep -q .; then
        echo "error: [env-preflight] LIBBLOSC_DIR contains no libblosc.so*: $LIBBLOSC_DIR" >&2
        exit 7
    fi
    echo "[env-preflight] libblosc.so* found under LIBBLOSC_DIR."
elif [[ ! -f "$CONDA_PREFIX/lib/libblosc.so" ]]; then
    echo "[env-preflight] WARN: LIBBLOSC_DIR unset and no libblosc.so in \$CONDA_PREFIX/lib." >&2
    echo "[env-preflight]   raw2ometiff may fail with 'Unable to load library blosc'." >&2
fi

echo "[env-preflight] all checks passed."
