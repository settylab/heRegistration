#!/usr/bin/env bash
# scripts/pip-install.sh — run `pip install` with an isolated
# PIP_CACHE_DIR, self-verifying that ~/.cache/pip stayed untouched.
#
# Unlike the conda/mamba side (scripts/create-env.sh, which overrides
# CONDA_PKGS_DIRS whenever MAMBA_ROOT_PREFIX is set), pip has no
# equivalent isolation by default — a plain `pip install` always
# writes its HTTP + wheel cache to `~/.cache/pip`, regardless of
# MAMBA_ROOT_PREFIX. Nothing errors when this happens: a "fully
# isolated" install (docs/install.md § 2's `MAMBA_ROOT_PREFIX`
# guidance) can still quietly leave ~170+ files behind in the user's
# home directory. This wraps `pip install` with an explicit
# `PIP_CACHE_DIR` override (when MAMBA_ROOT_PREFIX is set) and
# asserts, after the fact, that `~/.cache/pip` was not touched — same
# shape as create-env.sh's `~/.mamba/pkgs` check, because the failure
# mode is the same: silent, not loud.
#
# If MAMBA_ROOT_PREFIX is NOT set, this is a thin passthrough to `pip
# install` (no isolation was requested, so there's nothing to
# isolate or assert).
#
# Usage:
#   scripts/pip-install.sh --no-deps -r environments/heRegistration-requirements.txt
#   scripts/pip-install.sh --no-deps -e .
#   MAMBA_ROOT_PREFIX=/abs/isolated/root scripts/pip-install.sh --no-deps -e .
#
# All arguments are passed through to `pip install` verbatim.

set -euo pipefail

if ! command -v pip >/dev/null 2>&1; then
    echo "error: [pip-install] pip not found on PATH. Did you run 'micromamba activate heRegistration'?" >&2
    exit 4
fi

HOME_PIP_CACHE="$HOME/.cache/pip"

if [[ -z "${MAMBA_ROOT_PREFIX:-}" ]]; then
    echo "[pip-install] MAMBA_ROOT_PREFIX not set — using pip's default cache; nothing to isolate."
    pip install "$@"
    exit 0
fi

if [[ ! -d "$MAMBA_ROOT_PREFIX" ]]; then
    mkdir -p "$MAMBA_ROOT_PREFIX"
fi
MAMBA_ROOT_PREFIX=$(cd "$MAMBA_ROOT_PREFIX" && pwd)
export PIP_CACHE_DIR="$MAMBA_ROOT_PREFIX/pip-cache"
mkdir -p "$PIP_CACHE_DIR"

echo "[pip-install] PIP_CACHE_DIR=$PIP_CACHE_DIR"

MARKER=$(mktemp)
trap 'rm -f "$MARKER"' EXIT

pip install "$@"

if [[ -d "$HOME_PIP_CACHE" ]]; then
    # Exclude ~/.cache/pip/selfcheck/ — pip's version-check timestamp
    # file writes there unconditionally (hardcoded in
    # pip._internal.self_outdated_check, NOT controlled by
    # PIP_CACHE_DIR). It has no cache-semantics impact — it just records
    # "when did I last check whether a newer pip exists". Users who want
    # to suppress it can also set PIP_DISABLE_PIP_VERSION_CHECK=1 in
    # their environment.
    TOUCHED=$(find "$HOME_PIP_CACHE" -newer "$MARKER" ! -path '*/selfcheck*' 2>/dev/null || true)
    if [[ -n "$TOUCHED" ]]; then
        echo "error: [pip-install] $HOME_PIP_CACHE was modified during pip install" >&2
        echo "       despite PIP_CACHE_DIR=$PIP_CACHE_DIR — isolation did NOT hold." >&2
        echo "       Files touched:" >&2
        echo "$TOUCHED" | sed 's/^/       /' >&2
        exit 1
    fi
fi
echo "[pip-install] verified: $HOME_PIP_CACHE was not modified. Isolation held."
