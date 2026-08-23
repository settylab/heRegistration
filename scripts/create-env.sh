#!/usr/bin/env bash
# scripts/create-env.sh — create the micromamba env used by this
# pipeline, self-verifying that the package cache stayed isolated.
#
# micromamba's `pkgs_dirs` resolves to `[<root>/pkgs, ~/.mamba/pkgs]`
# even when `MAMBA_ROOT_PREFIX` is overridden to an isolated location
# — an undocumented second entry. Nothing errors when this fires: a
# user who sets MAMBA_ROOT_PREFIX expecting a fully self-contained
# install can still end up with package-cache state silently written
# to `~/.mamba/pkgs`. This wraps `micromamba env create` with an
# explicit `CONDA_PKGS_DIRS` override (when MAMBA_ROOT_PREFIX is set)
# and asserts, after the fact, that `~/.mamba/pkgs` was not touched —
# the failure mode is silent, so the invariant has to check itself
# rather than being trusted.
#
# If MAMBA_ROOT_PREFIX is NOT set, this is a thin passthrough to
# `micromamba env create` (the default root already IS ~/.mamba, so
# there is nothing to isolate and nothing to assert).
#
# Usage:
#   scripts/create-env.sh -n heRegistration -f environments/heRegistration.yml
#   MAMBA_ROOT_PREFIX=/abs/isolated/root scripts/create-env.sh \
#       -n heRegistration -f environments/heRegistration.yml
#
# All arguments are passed through to `micromamba env create` verbatim.

set -euo pipefail

if ! command -v micromamba >/dev/null 2>&1; then
    echo "error: [create-env] micromamba not found on PATH." >&2
    exit 4
fi

HOME_PKGS="$HOME/.mamba/pkgs"

if [[ -z "${MAMBA_ROOT_PREFIX:-}" ]]; then
    echo "[create-env] MAMBA_ROOT_PREFIX not set — using micromamba's default root; nothing to isolate."
    exec micromamba env create "$@"
fi

if [[ ! -d "$MAMBA_ROOT_PREFIX" ]]; then
    mkdir -p "$MAMBA_ROOT_PREFIX"
fi
MAMBA_ROOT_PREFIX=$(cd "$MAMBA_ROOT_PREFIX" && pwd)
export MAMBA_ROOT_PREFIX
export CONDA_PKGS_DIRS="$MAMBA_ROOT_PREFIX/pkgs"

echo "[create-env] MAMBA_ROOT_PREFIX=$MAMBA_ROOT_PREFIX"
echo "[create-env] CONDA_PKGS_DIRS=$CONDA_PKGS_DIRS"

MARKER=$(mktemp)
trap 'rm -f "$MARKER"' EXIT

micromamba env create "$@"

if [[ -d "$HOME_PKGS" ]]; then
    TOUCHED=$(find "$HOME_PKGS" -newer "$MARKER" 2>/dev/null || true)
    if [[ -n "$TOUCHED" ]]; then
        echo "error: [create-env] $HOME_PKGS was modified during env creation" >&2
        echo "       despite CONDA_PKGS_DIRS=$CONDA_PKGS_DIRS — isolation did" >&2
        echo "       NOT hold. Files touched:" >&2
        echo "$TOUCHED" | sed 's/^/       /' >&2
        exit 1
    fi
fi
echo "[create-env] verified: $HOME_PKGS was not modified. Isolation held."
