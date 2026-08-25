#!/usr/bin/env bash
# scripts/create-env.sh — run any `micromamba` env-creation subcommand
# (`create`, `env create`, …), self-verifying that the package cache
# stayed isolated.
#
# micromamba's `pkgs_dirs` resolves to `[<root>/pkgs, ~/.mamba/pkgs]`
# even when `MAMBA_ROOT_PREFIX` is overridden to an isolated location
# — an undocumented second entry. Nothing errors when this fires: a
# user who sets MAMBA_ROOT_PREFIX expecting a fully self-contained
# install can still end up with package-cache state silently written
# to `~/.mamba/pkgs`. This wraps the `micromamba` call with an
# explicit `CONDA_PKGS_DIRS` override (when MAMBA_ROOT_PREFIX is set)
# and asserts, after the fact, that `~/.mamba/pkgs` was not touched —
# the failure mode is silent, so the invariant has to check itself
# rather than being trusted.
#
# If MAMBA_ROOT_PREFIX is NOT set, this is a thin passthrough to
# `micromamba` (the default root already IS ~/.mamba, so there is
# nothing to isolate and nothing to assert).
#
# Deliberately generic (does not hardcode `create` vs `env create`):
# this repo uses both — `env create -f environments/heRegistration.yml`
# for the main env, and plain `create -n <name> -c <channel> <pkg>`
# for the ad-hoc blosc env (docs/install.md "VSI-input prerequisites")
# — and `env create` does not accept the latter's positional package
# argument.
#
# Separately, micromamba ALWAYS registers every env it creates into a
# second, global, `$HOME`-scoped registry — `~/.conda/environments.txt`
# — regardless of `MAMBA_ROOT_PREFIX` and with no flag to opt out
# (verified against `micromamba create --help`: no such flag exists).
# Unlike `~/.mamba/pkgs`, this one cannot be avoided, so this script
# does not assert it's untouched (it always will be) — instead it (a)
# says so out loud below instead of silently omitting it from the
# isolation report, and (b) resolves and records the env's own prefix
# itself (see "env prefix" below), so nothing downstream needs to
# scrape that registry to find it — the registry accumulates one row
# per env ever created on this account, across every root, so grepping
# or eyeballing it for a name match returns every prior install too,
# not just this one.
#
# Usage:
#   scripts/create-env.sh env create -n heRegistration -f environments/heRegistration.yml
#   scripts/create-env.sh create -n blosc -c conda-forge blosc -y
#   MAMBA_ROOT_PREFIX=/abs/isolated/root scripts/create-env.sh \
#       env create -n heRegistration -f environments/heRegistration.yml
#
# All arguments are passed through to `micromamba` verbatim.
#
# Side effect: if a `-n`/`--name` env name is present in the arguments,
# writes the resolved absolute env prefix to
# `scripts/.env-prefix-<name>` (gitignored) — `write-env-config.sh
# --env-name <name>` reads this back instead of asking the user to
# find the path by hand via `micromamba env list`.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Resolve the micromamba binary. Prefer $MAMBA_EXE (exported by
# micromamba's shell-hook to the exact binary that sourced the hook —
# known-working by construction) over a bare PATH lookup; a shadowed
# or wrong-arch `micromamba` earlier on PATH would otherwise fail
# with `Exec format error` on exec, well after `command -v` was happy.
if [[ -n "${MAMBA_EXE:-}" && -x "$MAMBA_EXE" ]]; then
    MAMBA_BIN="$MAMBA_EXE"
elif MAMBA_BIN=$(command -v micromamba 2>/dev/null); then
    :
else
    echo "error: [create-env] micromamba not found (\$MAMBA_EXE unset and not on PATH)." >&2
    exit 4
fi

HOME_PKGS="$HOME/.mamba/pkgs"
HOME_ENVS_REGISTRY="$HOME/.conda/environments.txt"

ENV_NAME=""
prev=""
for arg in "$@"; do
    case "$arg" in
        --name=*) ENV_NAME="${arg#--name=}" ;;
    esac
    if [[ "$prev" == "-n" || "$prev" == "--name" ]]; then
        ENV_NAME="$arg"
    fi
    prev="$arg"
done

# Resolve the env's own absolute prefix and record it — deterministic,
# so nothing downstream needs to scrape `micromamba env list` (which
# returns one row per env ever created on this account, across every
# MAMBA_ROOT_PREFIX ever used, not just this one).
record_env_prefix() {
    if [[ -z "$ENV_NAME" ]]; then
        echo "[create-env] note: no -n/--name found in arguments; skipping env-prefix receipt." >&2
        return 0
    fi
    local envs_dir
    if [[ -n "${MAMBA_ROOT_PREFIX:-}" ]]; then
        envs_dir="$MAMBA_ROOT_PREFIX/envs"
    else
        envs_dir=$("$MAMBA_BIN" info | sed -n 's/^[[:space:]]*envs directories[[:space:]]*:[[:space:]]*//p' | head -1)
    fi
    if [[ -z "$envs_dir" ]]; then
        echo "[create-env] warning: could not resolve the envs directory; skipping env-prefix receipt." >&2
        return 0
    fi
    local prefix="$envs_dir/$ENV_NAME"
    local receipt="$SCRIPT_DIR/.env-prefix-$ENV_NAME"
    echo "$prefix" > "$receipt"
    echo "[create-env] env prefix: $prefix"
    echo "[create-env] recorded to $receipt"
}

if [[ -z "${MAMBA_ROOT_PREFIX:-}" ]]; then
    echo "[create-env] MAMBA_ROOT_PREFIX not set — using micromamba's default root; nothing to isolate."
    "$MAMBA_BIN" "$@"
    record_env_prefix
    exit 0
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

"$MAMBA_BIN" "$@"

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

if [[ -f "$HOME_ENVS_REGISTRY" ]] && [[ "$HOME_ENVS_REGISTRY" -nt "$MARKER" ]]; then
    echo "[create-env] note: $HOME_ENVS_REGISTRY was appended to during env creation." >&2
    echo "             This is unavoidable — micromamba registers every env's path into" >&2
    echo "             this global, \$HOME-scoped registry regardless of" >&2
    echo "             MAMBA_ROOT_PREFIX, with no flag to opt out. Only a path string is" >&2
    echo "             written (no package data), so it does not violate this wrapper's" >&2
    echo "             package-cache isolation guarantee above — but 'micromamba env" >&2
    echo "             list' will show this env, plus every other env ever created on" >&2
    echo "             this account, going forward. Use the recorded prefix (see below)" >&2
    echo "             instead of that command." >&2
fi

record_env_prefix
