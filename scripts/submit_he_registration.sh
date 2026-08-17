#!/bin/bash -l
# ---------------------------------------------------------------------
# Slurm submission wrapper for the `hexenium` pipeline.
# `bash -l` makes this a LOGIN shell so ~/.bash_profile (and indirectly
# ~/.bashrc on most setups) gets sourced — that's what initialises
# micromamba in interactive sessions but is otherwise skipped in
# non-interactive Slurm batch jobs.
#
# Usage:
#   sbatch scripts/submit_he_registration.sh SAMPLE_ID HE_PATH XENIUM_BUNDLE [CELLTYPE_CSV] [--key value ...]
#
# Example:
#   sbatch scripts/submit_he_registration.sh SAMPLE1 \
#       /path/to/HE.ome.tif \
#       /path/to/output-XETG... \
#       /path/to/celltypes.csv
#
# Logs land at: <output_root>/<sample_id>/logs/<job_name>_<jobid>.log
#
# No #SBATCH --output / --error directives — Slurm parses those at
# submit time using your current cwd as the base, and if the cwd isn't
# writable the log file can't be created and your job fails silently.
# Instead, we redirect everything to a path under YOUR output root once
# we know it (which the pipeline owns and `mkdir -p`s on demand).
# ---------------------------------------------------------------------
#SBATCH --job-name=hexenium
#SBATCH --partition=campus-new
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00

set -euo pipefail

# ----- Positional args (parsed FIRST so we can route the log) ---------
SAMPLE_ID="${1:?usage: sbatch scripts/submit_he_registration.sh SAMPLE_ID HE_PATH XENIUM_BUNDLE [CELLTYPE_CSV] [extra --flags]}"
HE_PATH="${2:?missing HE_PATH}"
XENIUM_BUNDLE="${3:?missing XENIUM_BUNDLE}"
CELLTYPE_CSV="${4:-}"
shift $(( $# < 4 ? $# : 4 ))   # consume up to 4 positional args
EXTRA_ARGS=("$@")              # anything left = --key value overrides

# ----- Resolve OUTPUT_ROOT (mirrors hexenium's precedence) ------------
# Precedence: --output-root in EXTRA_ARGS > $OUTPUT_ROOT env > (unset -> abort).
# We intentionally do NOT ship a site-specific default; set OUTPUT_ROOT in
# your env or pass --output-root on the command line.
OUTPUT_ROOT="${OUTPUT_ROOT:-}"
for ((i=0; i<${#EXTRA_ARGS[@]}; i++)); do
    if [[ "${EXTRA_ARGS[i]}" == "--output-root" && $((i+1)) -lt ${#EXTRA_ARGS[@]} ]]; then
        OUTPUT_ROOT="${EXTRA_ARGS[i+1]}"
        break
    fi
done
if [[ -z "$OUTPUT_ROOT" ]]; then
    echo "[submit] ERROR: OUTPUT_ROOT is not set." >&2
    echo "[submit]   Set OUTPUT_ROOT=/path/to/runs in your env, or pass --output-root /path/to/runs after the positional args." >&2
    exit 2
fi

# ----- Create the user-facing log location + REDIRECT (early) --------
SAMPLE_OUT="$OUTPUT_ROOT/$SAMPLE_ID"
LOG_DIR="$SAMPLE_OUT/logs"
mkdir -p "$LOG_DIR"

JOB_TAG="${SLURM_JOB_NAME:-hexenium}_${SLURM_JOB_ID:-local-$(date +%Y%m%d-%H%M%S)}"
LOG_FILE="$LOG_DIR/${JOB_TAG}.log"

# ----- Validate paths BEFORE the tee redirect ------------------------
# Fail loud here on missing inputs so the operator sees the problem
# in stderr immediately and the Slurm job doesn't waste a time slot.
fail=0
for arg_name in HE_PATH XENIUM_BUNDLE; do
    arg_val="${!arg_name}"
    if [[ ! -e "$arg_val" ]]; then
        echo "[submit] ERROR: $arg_name does not exist: $arg_val" >&2
        fail=1
    fi
done
if [[ -n "$CELLTYPE_CSV" && ! -e "$CELLTYPE_CSV" ]]; then
    echo "[submit] ERROR: CELLTYPE_CSV does not exist: $CELLTYPE_CSV" >&2
    fail=1
fi
# --dapi-path override (if present in EXTRA_ARGS) also gets checked.
for ((i=0; i<${#EXTRA_ARGS[@]}; i++)); do
    if [[ "${EXTRA_ARGS[i]}" == "--dapi-path" && $((i+1)) -lt ${#EXTRA_ARGS[@]} ]]; then
        dp="${EXTRA_ARGS[i+1]}"
        if [[ ! -e "$dp" ]]; then
            echo "[submit] ERROR: --dapi-path does not exist: $dp" >&2
            fail=1
        fi
        break
    fi
done
if [[ $fail -ne 0 ]]; then
    echo "[submit] aborting before redirect; fix the path(s) above and resubmit." >&2
    exit 2
fi

# Redirect EVERYTHING from this point on to $LOG_FILE.
#
# TWO PATHS, chosen by execution context:
#
#   sbatch batch script  (SLURM_JOB_ID set AND SLURM_STEP_ID unset)
#       -> plain file redirect: `exec >> "$LOG_FILE" 2>&1`.
#       No tee, no process substitution. python -u keeps output
#       unbuffered, so we don't need tee's line-buffered flush.
#       Loss: no live console echo -- but sbatch has no controlling
#       terminal anyway, and slurm-<jobid>.out already captured the
#       pre-redirect lines.
#
#   interactive `bash submit.slurm.sh` or `srun bash submit.slurm.sh`
#       -> keep the tee-in-process-substitution for live console echo.
#       stdbuf -oL forces per-line flush (vs. GNU tee's default ~4 KB
#       block buffer) so tracebacks show up promptly.
#
# WHY split? Under sbatch (no controlling terminal + shared-filesystem
# latency), a chatty subprocess (VALIS/Java) can burst output faster
# than tee can flush; tee's stdin pipe buffer (64 KB kernel default)
# fills, tee's write blocks on FS latency, python's next write blocks
# on the full pipe, and with `set -euo pipefail` there is no SIGPIPE
# escape hatch -- the whole job deadlocks silently holding its full
# allocation. An interactive terminal drains fast enough to keep tee
# flowing, so the tee path is only unsafe under sbatch.
if [[ -n "${SLURM_JOB_ID:-}" && -z "${SLURM_STEP_ID:-}" ]]; then
    exec >> "$LOG_FILE" 2>&1
else
    exec > >(stdbuf -oL -eL tee -a "$LOG_FILE") 2> >(stdbuf -oL -eL tee -a "$LOG_FILE" >&2)
fi

echo "[submit] log file: $LOG_FILE"

# ----- Conda / micromamba / mamba env activation ---------------------
# Slurm batch jobs run non-interactive shells. Interactive setups put
# micromamba's init in ~/.bashrc (which defines a SHELL FUNCTION
# `micromamba`, not a binary) — non-interactive shells don't source
# .bashrc, so the function is missing and `micromamba activate` fails.
#
# Fix in three layers:
#   1) source ~/.bashrc explicitly (most direct).
#   2) ensure MAMBA_ROOT_PREFIX is set + the micromamba binary is on PATH.
#   3) source the micromamba hook script to (re)define the shell function.
#
# ENV_NAME defaults to heRegistration; override via
# `ENV_NAME=other_env sbatch scripts/submit_he_registration.sh ...`.
ENV_NAME="${ENV_NAME:-heRegistration}"
activated=0

# 1) Replay the interactive shell init. set +e so a noisy .bashrc
# doesn't kill us under `set -euo pipefail`.
if [[ -f "$HOME/.bashrc" ]]; then
    echo "[submit] sourcing $HOME/.bashrc"
    set +e; set +u
    # shellcheck disable=SC1091
    source "$HOME/.bashrc"
    set -e; set -u
fi

# 2) Make sure MAMBA_ROOT_PREFIX is set and the binary is reachable.
: "${MAMBA_ROOT_PREFIX:=$HOME/micromamba}"
export MAMBA_ROOT_PREFIX
if ! command -v micromamba >/dev/null 2>&1; then
    for bindir in \
        "$MAMBA_ROOT_PREFIX/bin" \
        "$HOME/.local/bin" \
        "$HOME/micromamba/bin" \
        "/app/software/micromamba/bin"; do
        if [[ -x "$bindir/micromamba" ]]; then
            echo "[submit] adding micromamba bin dir to PATH: $bindir"
            export PATH="$bindir:$PATH"
            break
        fi
    done
fi

# 3) Ensure the shell function exists (not just the binary).
if [[ "$(type -t micromamba 2>/dev/null)" != "function" ]]; then
    for hook in \
        "$MAMBA_ROOT_PREFIX/etc/profile.d/micromamba.sh" \
        "$HOME/micromamba/etc/profile.d/micromamba.sh" \
        "$HOME/.local/share/mamba/etc/profile.d/micromamba.sh" \
        "/app/software/micromamba/etc/profile.d/micromamba.sh"; do
        if [[ -f "$hook" ]]; then
            echo "[submit] sourcing micromamba hook: $hook"
            # shellcheck disable=SC1090
            source "$hook"
            break
        fi
    done
fi

# Now try to activate.
if command -v micromamba >/dev/null 2>&1; then
    if micromamba activate "$ENV_NAME" 2>/dev/null; then
        activated=1
        echo "[submit] activated env via micromamba: $ENV_NAME"
    else
        echo "[submit] WARN: micromamba activate '$ENV_NAME' failed; trying mamba/conda"
    fi
fi

# 4) mamba fallback.
if [[ $activated -eq 0 ]] && command -v mamba >/dev/null 2>&1; then
    eval "$(mamba shell hook -s bash 2>/dev/null)" || true
    if mamba activate "$ENV_NAME" 2>/dev/null; then
        activated=1
        echo "[submit] activated env via mamba: $ENV_NAME"
    fi
fi

# 5) conda fallback.
if [[ $activated -eq 0 ]] && command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" 2>/dev/null || true
    if conda activate "$ENV_NAME" 2>/dev/null; then
        activated=1
        echo "[submit] activated env via conda: $ENV_NAME"
    fi
fi

if [[ $activated -eq 0 ]]; then
    echo "[submit] ERROR: could not activate env '$ENV_NAME'." >&2
    echo "[submit]   Discovered tools:" >&2
    echo "[submit]     micromamba: $(command -v micromamba 2>/dev/null || echo 'NOT FOUND')" >&2
    echo "[submit]     mamba:      $(command -v mamba      2>/dev/null || echo 'NOT FOUND')" >&2
    echo "[submit]     conda:      $(command -v conda      2>/dev/null || echo 'NOT FOUND')" >&2
    echo "[submit]     MAMBA_ROOT_PREFIX=$MAMBA_ROOT_PREFIX" >&2
    echo "[submit]   Available envs (best effort):" >&2
    micromamba env list 2>/dev/null || true
    echo "[submit]   Override env name: ENV_NAME=your_env sbatch scripts/submit_he_registration.sh ..." >&2
    exit 1
fi
# Verify python is the env's python, not the system one.
echo "[submit] which python: $(command -v python)"
echo "[submit] which hexenium: $(command -v hexenium 2>/dev/null || echo 'NOT FOUND — did you `pip install -e .` inside the env?')"
echo "[submit] CONDA_PREFIX: ${CONDA_PREFIX:-unset}"

# ----- Slurm sanity ---------------------------------------------------
echo "===== SLURM INFO ====="
echo "JOBID:        ${SLURM_JOB_ID:-NA}"
echo "NODELIST:     ${SLURM_NODELIST:-NA}"
echo "HOST:         $(hostname)"
echo "PWD:          $(pwd)"
echo "SUBMIT_DIR:   ${SLURM_SUBMIT_DIR:-NA}"
echo "DATE:         $(date)"
echo "PYTHON:       $(which python)"
echo "OUTPUT_ROOT:  $OUTPUT_ROOT"
echo "SAMPLE_OUT:   $SAMPLE_OUT"
echo "LOG_FILE:     $LOG_FILE"
python --version || true
echo "======================"

# ----- Run ------------------------------------------------------------
ARGS=(
    run
    --sample-id "$SAMPLE_ID"
    --he-path "$HE_PATH"
    --xenium-bundle "$XENIUM_BUNDLE"
)
if [[ -n "$CELLTYPE_CSV" ]]; then
    ARGS+=( --celltype-csv "$CELLTYPE_CSV" )
fi
ARGS+=( "${EXTRA_ARGS[@]}" )

echo "[submit] running:"
echo "    hexenium ${ARGS[*]}"
# PYTHONUNBUFFERED=1 + python -u together force unbuffered stdout/stderr.
# Without this, Python's prints sit in OS pipe buffers and are LOST when
# Slurm SIGKILLs at time-limit. Belt + braces because some libraries
# reopen FDs.
export PYTHONUNBUFFERED=1
hexenium "${ARGS[@]}"

echo "[submit] DONE at $(date)"
