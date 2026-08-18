#!/bin/bash -l
# ---------------------------------------------------------------------
# Slurm submission wrapper for the `hexenium` pipeline.
#
# `bash -l` makes this a LOGIN shell so ~/.bash_profile (and indirectly
# ~/.bashrc on most setups) gets sourced — that's what initialises
# micromamba in interactive sessions but is otherwise skipped in
# non-interactive Slurm batch jobs.
#
# Usage — call directly (no `sbatch` prefix); the script self-submits
# to slurm and routes its own .out/.err under the run folder:
#
#   ./scripts/submit_he_registration.sh \
#       --sample-id     SAMPLE1 \
#       --run-id        demo_v1 \
#       --output-root   /path/to/runs \
#       --he-slide      /path/to/aligned_fullres_HE.ome.tif \
#       --xenium-bundle /path/to/output-XETG... \
#       [--stages       register warp celltype viz] \
#       [--force-rerun]
#
# Legacy `sbatch scripts/submit_he_registration.sh ...` also works; slurm's
# .out/.err then land in the submit directory instead of the integrated
# run folder.
#
# Environment variables:
#   OUTPUT_ROOT — fallback for --output-root when not passed as a flag.
#   ENV_NAME    — conda/micromamba env name (default: heRegistration).
# ---------------------------------------------------------------------
#SBATCH --job-name=hexenium
#SBATCH --partition=campus-new
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=196G
#SBATCH --time=2-00:00:00

set -euo pipefail

# Canonical stage list — MUST stay in sync with hexenium.layout.STAGE_NAMES.
DEFAULT_STAGES=(he_preprocess register warp celltype viz)

# Compute the stages-suffix from the CLI args. Mirrors
# hexenium.layout.compute_stages_suffix: empty (argparse default) → "all";
# full canonical set (any order) → "all"; otherwise underscore-joined.
_compute_stages_suffix() {
    local -a stages=("$@")
    if [[ ${#stages[@]} -eq 0 ]]; then
        echo "all"
        return
    fi
    local sorted_actual
    sorted_actual=$(printf '%s\n' "${stages[@]}" | sort | tr '\n' ' ')
    local sorted_default
    sorted_default=$(printf '%s\n' "${DEFAULT_STAGES[@]}" | sort | tr '\n' ' ')
    if [[ "$sorted_actual" == "$sorted_default" ]]; then
        echo "all"
    else
        (IFS=_; echo "${stages[*]}")
    fi
}

# ---------------------------------------------------------------------
# Launcher branch: when called directly (SLURM_JOB_ID unset), parse the
# routing-relevant flags, self-submit under sbatch with --output/--error
# routed to the integrated run folder, then exit. Under sbatch the
# script re-enters with SLURM_JOB_ID set and skips this block.
# ---------------------------------------------------------------------
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    SAMPLE_ID=""; RUN_ID=""; OUTPUT_ROOT_ARG=""
    STAGES=()
    args=("$@")
    for ((i=0; i<${#args[@]}; i++)); do
        case "${args[i]}" in
            --sample-id)   SAMPLE_ID="${args[i+1]:-}" ;;
            --run-id)      RUN_ID="${args[i+1]:-}" ;;
            --output-root) OUTPUT_ROOT_ARG="${args[i+1]:-}" ;;
            --stages)
                # nargs="+" — consume values until next flag or end.
                j=$((i+1))
                while [[ $j -lt ${#args[@]} && "${args[j]}" != --* ]]; do
                    STAGES+=("${args[j]}")
                    j=$((j+1))
                done
                ;;
        esac
    done
    OUTPUT_ROOT="${OUTPUT_ROOT_ARG:-${OUTPUT_ROOT:-}}"
    STAGES_SUFFIX=$(_compute_stages_suffix "${STAGES[@]}")

    if [[ -z "$OUTPUT_ROOT" ]]; then
        echo "[submit] ERROR: OUTPUT_ROOT is not set." >&2
        echo "[submit]   Set OUTPUT_ROOT=/path/to/runs in your env, or pass" >&2
        echo "[submit]   --output-root /path/to/runs on the command line." >&2
        exit 2
    fi
    if [[ -z "$SAMPLE_ID" ]]; then
        echo "[submit] ERROR: --sample-id is required." >&2
        exit 2
    fi

    # Two-phase submission: sbatch will silently drop .out/.err if the
    # --output parent dir doesn't exist at job start, but the job id is
    # only known AFTER sbatch. So submit --hold, mkdir the per-job dir,
    # then release.
    if [[ -n "$RUN_ID" ]]; then
        # Integrated: colocate under the run dir with any upstream logs.
        LOG_DIR_BASE="$OUTPUT_ROOT/$SAMPLE_ID/${SAMPLE_ID}_${RUN_ID}/logs/logs_heRegistration"
    else
        # Standalone: no shared run dir.
        LOG_DIR_BASE="$OUTPUT_ROOT/$SAMPLE_ID/logs"
    fi
    LOG_LEAF_TMPL="${SAMPLE_ID}_%j_${STAGES_SUFFIX}"
    mkdir -p "$LOG_DIR_BASE"
    JOBID=$(sbatch --parsable --hold \
        --output="$LOG_DIR_BASE/${LOG_LEAF_TMPL}/slurm-%j.out" \
        --error="$LOG_DIR_BASE/${LOG_LEAF_TMPL}/slurm-%j.err" \
        "$0" "$@")
    if [[ -z "$JOBID" ]]; then
        echo "[submit] sbatch failed to return a job id" >&2
        exit 1
    fi
    mkdir -p "$LOG_DIR_BASE/${SAMPLE_ID}_${JOBID}_${STAGES_SUFFIX}"
    scontrol release "$JOBID"
    echo "Submitted batch job $JOBID (logs: ${SAMPLE_ID}_${JOBID}_${STAGES_SUFFIX})"
    exit 0
fi

# ---------------------------------------------------------------------
# Under-sbatch branch: activate the env and run hexenium.
# ---------------------------------------------------------------------

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
# `ENV_NAME=other_env ./scripts/submit_he_registration.sh ...`.
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
    echo "[submit]   Override env name: ENV_NAME=your_env ./scripts/submit_he_registration.sh ..." >&2
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
python --version || true
echo "======================"

# ----- Run ------------------------------------------------------------
# PYTHONUNBUFFERED=1 keeps stdout/stderr flushed so nothing is lost when
# Slurm SIGKILLs at time-limit.
export PYTHONUNBUFFERED=1
echo "[submit] running: hexenium run $*"
hexenium run "$@"

echo "[submit] DONE at $(date)"
