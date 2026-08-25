#!/usr/bin/env bash
# ---------------------------------------------------------------------
# Slurm submission wrapper for the `hexenium` pipeline.
#
# NOT a login shell (`bash -l`): a login shell implicitly sources
# ~/.bash_profile (and, on most setups, indirectly ~/.bashrc) via
# `$HOME` at process startup — before this script's own body runs at
# all. That's exactly the class of bug scripts/lib/env_config.sh exists
# to eliminate (see its rationale comment), so it can't be the
# mechanism that initialises micromamba here. Env activation instead
# uses the resolved absolute values in scripts/env.local.conf — see the
# "Env activation" section below.
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
# Auto-VSI (unified Path B) — if --he-slide is a .vsi file AND
# he_preprocess is in --stages (default), the launcher automatically
# submits a bioformats2raw + raw2ometiff conversion job FIRST,
# rewrites --he-slide in the hexenium args to the converted OME-TIFF,
# and submits hexenium with --dependency=afterok:<conv_jobid>. On
# successful reruns the converted file is reused unless
# --force-preprocess is set. The conversion and hexenium jobs remain
# independent slurm records so hexenium re-runs don't need to redo
# the conversion. See:
#   scripts/submit_vsi_to_ometiff.sh (standalone conversion, advanced)
#   scripts/submit_hexenium_from_ometiff.sh (standalone OME-TIFF entry
#     point, advanced)
#
# Legacy `sbatch scripts/submit_he_registration.sh ...` also works; slurm's
# .out/.err then land in the submit directory instead of the integrated
# run folder.
#
# Environment variables:
#   OUTPUT_ROOT — fallback for --output-root when not passed as a flag.
#
# Env activation: always by the absolute prefix pinned in
# scripts/env.local.conf (scripts/write-env-config.sh) — never by name,
# never re-derived from $HOME. See scripts/lib/env_config.sh for why.
# --env-name / $ENV_NAME are still parsed (so old invocations don't hit
# an "unknown flag" error) but only WARN if given; they no longer
# select which env is activated. Re-run write-env-config.sh
# --env-prefix <path> to point at a different env.
#
# VSI-specific flags (only honoured when --he-slide is a .vsi file):
#   --force-preprocess  — re-run the VSI → OME-TIFF conversion even if
#                         a converted file already exists in the
#                         canonical location.
#   --vsi-series <N>    — force a specific bioformats2raw series index.
#                         Default: auto-pick the largest-on-disk
#                         per-series file (a robust proxy for the
#                         full-resolution 40x brightfield image).
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
    HE_SLIDE=""; ENV_NAME_CLI=""; VSI_SERIES=""; FORCE_PREPROCESS=0
    STAGES=()
    args=("$@")
    for ((i=0; i<${#args[@]}; i++)); do
        case "${args[i]}" in
            --sample-id)        SAMPLE_ID="${args[i+1]:-}" ;;
            --run-id)           RUN_ID="${args[i+1]:-}" ;;
            --output-root)      OUTPUT_ROOT_ARG="${args[i+1]:-}" ;;
            --he-slide|--he-path) HE_SLIDE="${args[i+1]:-}" ;;
            --env-name)         ENV_NAME_CLI="${args[i+1]:-}" ;;
            --vsi-series)       VSI_SERIES="${args[i+1]:-}" ;;
            --force-preprocess) FORCE_PREPROCESS=1 ;;
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

    # -----------------------------------------------------------------
    # Auto-VSI (unified Path B).
    #
    # If the H&E input is a .vsi file, submit a conversion job first
    # and rewrite --he-slide to the converted OME-TIFF's canonical
    # location BEFORE submitting the hexenium job. The hexenium job
    # then runs with --dependency=afterok:<conv_jobid>, so
    # (a) hexenium doesn't start until conversion finishes cleanly, and
    # (b) hexenium reruns don't touch the conversion job's state.
    #
    # Idempotent: if the converted file already exists in the canonical
    # location (successful prior conversion), skip conversion entirely
    # and just submit hexenium. --force-preprocess opts back into a
    # full re-conversion.
    #
    # OME-TIFF / .tif inputs bypass this block entirely (no conversion
    # needed; the below sbatch call runs hexenium as-is).
    # -----------------------------------------------------------------
    CONV_JOBID=""
    _he_lc="${HE_SLIDE,,}"
    if [[ "$_he_lc" == *.vsi ]]; then
        # Refuse if user opted out of he_preprocess with a .vsi input —
        # register/warp would crash at OpenSlide open on the raw .vsi.
        _stages_for_check=("${STAGES[@]:-${DEFAULT_STAGES[@]}}")
        _has_he_preprocess=0
        for s in "${_stages_for_check[@]}"; do
            [[ "$s" == "he_preprocess" ]] && _has_he_preprocess=1
        done
        if [[ $_has_he_preprocess -eq 0 ]]; then
            echo "[submit] ERROR: --he-slide is a .vsi file but he_preprocess" >&2
            echo "[submit]        was dropped from --stages. Registration would" >&2
            echo "[submit]        crash on the raw VSI. Either:" >&2
            echo "[submit]          - keep he_preprocess in --stages, OR" >&2
            echo "[submit]          - pre-convert to OME-TIFF and pass that path." >&2
            exit 2
        fi

        if [[ ! -f "$HE_SLIDE" ]]; then
            echo "[submit] ERROR: --he-slide does not exist: $HE_SLIDE" >&2
            exit 2
        fi

        # Canonical location for the converted OME-TIFF — matches
        # hexenium.layout.RunLayout.converted_dir + the pipeline-tree
        # fallback used by hexenium.stages.he_preprocess.discover_existing_ometiff.
        if [[ -n "$RUN_ID" ]]; then
            _run_dir="$OUTPUT_ROOT/$SAMPLE_ID/${SAMPLE_ID}_${RUN_ID}"
        else
            _run_dir="$OUTPUT_ROOT/$SAMPLE_ID"
        fi
        CONV_DIR="$_run_dir/he_registration/converted"
        CONV_OMETIFF="$CONV_DIR/${SAMPLE_ID}_he.ome.tif"
        mkdir -p "$CONV_DIR"

        # Rewrite --he-slide in the argv we pass to the hexenium sbatch
        # so the pipeline sees the OME-TIFF, not the raw .vsi. Under
        # sbatch, hexenium.stages.he_preprocess.is_vsi_input() returns
        # False for .ome.tif and the stage becomes a no-op.
        _rewritten_args=()
        _j=0
        _rewrote=0
        while [[ $_j -lt ${#args[@]} ]]; do
            case "${args[$_j]}" in
                --he-slide|--he-path)
                    _rewritten_args+=("${args[$_j]}" "$CONV_OMETIFF")
                    _j=$((_j+2))
                    _rewrote=1
                    ;;
                --vsi-series)
                    # Consume + drop (conversion-only flag; not for
                    # hexenium's argparse).
                    _j=$((_j+2))
                    ;;
                --force-preprocess)
                    # Consume + drop (this launcher's flag).
                    _j=$((_j+1))
                    ;;
                *)
                    _rewritten_args+=("${args[$_j]}")
                    _j=$((_j+1))
                    ;;
            esac
        done
        if [[ $_rewrote -eq 0 ]]; then
            # Defensive — should never trigger since $HE_SLIDE was
            # extracted from these same args.
            echo "[submit] internal error: --he-slide vanished during arg rewrite" >&2
            exit 1
        fi

        # Do we need to convert, or does a prior successful conversion
        # cover us?
        _need_convert=1
        if [[ $FORCE_PREPROCESS -eq 0 ]] && [[ -f "$CONV_OMETIFF" ]] &&
           [[ -s "$CONV_OMETIFF" ]]; then
            # -f: exists; -s: non-empty. Also accept symlinks to
            # non-empty files (submit_vsi_to_ometiff.sh writes the
            # canonical name as a symlink to the picked per-series file).
            _target="$(readlink -f "$CONV_OMETIFF" 2>/dev/null || echo "$CONV_OMETIFF")"
            if [[ -s "$_target" ]]; then
                echo "[submit] skipping VSI conversion — canonical OME-TIFF exists:"
                echo "[submit]   $CONV_OMETIFF"
                echo "[submit]   → $_target ($(du -h "$_target" | cut -f1))"
                echo "[submit]   Pass --force-preprocess to override."
                _need_convert=0
            fi
        fi

        if [[ $_need_convert -eq 1 ]]; then
            # Log dir for the conversion job. Sits alongside the
            # hexenium logs under logs_heRegistration/ so a triage
            # sweep of the run dir gathers everything in one place.
            if [[ -n "$RUN_ID" ]]; then
                _log_dir_base="$OUTPUT_ROOT/$SAMPLE_ID/${SAMPLE_ID}_${RUN_ID}/logs/logs_heRegistration"
            else
                _log_dir_base="$OUTPUT_ROOT/$SAMPLE_ID/logs"
            fi
            mkdir -p "$_log_dir_base"

            SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
            CONV_LAUNCHER="$SCRIPT_DIR/submit_vsi_to_ometiff.sh"
            if [[ ! -x "$CONV_LAUNCHER" ]]; then
                echo "[submit] ERROR: sibling launcher not found: $CONV_LAUNCHER" >&2
                exit 1
            fi

            # `submit_vsi_to_ometiff.sh` self-submits and prints
            # "Submitted batch job <id>" — capture the id from that.
            _conv_cmd=("$CONV_LAUNCHER"
                       --vsi "$HE_SLIDE"
                       --out-dir "$CONV_DIR"
                       --sample-id "$SAMPLE_ID")
            [[ -n "$ENV_NAME_CLI" ]] && _conv_cmd+=(--env-name "$ENV_NAME_CLI")
            [[ -n "$VSI_SERIES" ]] && _conv_cmd+=(--series "$VSI_SERIES")

            echo "[submit] VSI input detected — submitting conversion job:"
            printf '[submit]   %s\n' "${_conv_cmd[*]}"
            _conv_out="$("${_conv_cmd[@]}")"
            echo "$_conv_out"
            CONV_JOBID="$(printf '%s\n' "$_conv_out" \
                | grep -oE 'Submitted batch job [0-9]+' \
                | awk '{print $NF}' | head -1)"
            if [[ -z "$CONV_JOBID" ]]; then
                echo "[submit] ERROR: could not extract conversion job id" >&2
                echo "[submit]        from submit_vsi_to_ometiff.sh output." >&2
                exit 1
            fi
            echo "[submit] conversion job id: $CONV_JOBID"
            echo "[submit] hexenium job will wait for it via --dependency=afterok"
        fi

        # Hand off the rewritten argv (with --he-slide swapped for the
        # converted path, --vsi-series / --force-preprocess dropped) to
        # the sbatch submission below.
        set -- "${_rewritten_args[@]}"
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
    # Resolve this script's own scripts/ dir here (interactive shell —
    # BASH_SOURCE[0] is reliable), and forward it via --export so the
    # sbatch body (see line ~388 below) can source lib/env_config.sh
    # from the real repo tree, not from Slurm's per-job tmp copy of the
    # script at /var/tmp/slurmd/job<id>/slurm_script/ (which does NOT
    # include the surrounding scripts/lib/ dir).
    HEREG_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    export HEREG_SCRIPT_DIR
    # Assemble the sbatch call. If conversion was queued, wire the
    # afterok dependency in so hexenium starts only after the
    # conversion job SUCCEEDS (afterok, not afterany).
    _sbatch_args=(--parsable --hold
                  --output="$LOG_DIR_BASE/${LOG_LEAF_TMPL}/slurm-%j.out"
                  --error="$LOG_DIR_BASE/${LOG_LEAF_TMPL}/slurm-%j.err"
                  --export=ALL,HEREG_SCRIPT_DIR="$HEREG_SCRIPT_DIR")
    if [[ -n "$CONV_JOBID" ]]; then
        _sbatch_args+=(--dependency="afterok:$CONV_JOBID"
                       --kill-on-invalid-dep=yes)
    fi
    JOBID=$(sbatch "${_sbatch_args[@]}" "$0" "$@")
    if [[ -z "$JOBID" ]]; then
        echo "[submit] sbatch failed to return a job id" >&2
        exit 1
    fi
    mkdir -p "$LOG_DIR_BASE/${SAMPLE_ID}_${JOBID}_${STAGES_SUFFIX}"
    scontrol release "$JOBID"
    if [[ -n "$CONV_JOBID" ]]; then
        echo "Submitted batch job $JOBID (logs: ${SAMPLE_ID}_${JOBID}_${STAGES_SUFFIX})"
        echo "  waiting on conversion job $CONV_JOBID (afterok)"
    else
        echo "Submitted batch job $JOBID (logs: ${SAMPLE_ID}_${JOBID}_${STAGES_SUFFIX})"
    fi
    exit 0
fi

# ---------------------------------------------------------------------
# Under-sbatch branch: activate the env and run hexenium.
# ---------------------------------------------------------------------

# ----- Env activation ---------------------------------------------------
# Resolved from scripts/env.local.conf (scripts/write-env-config.sh),
# NEVER re-derived from `$HOME` here: three real sbatch submissions of
# this pipeline hung before the wrapper's own script body logged
# anything (0-byte stdout, no output dir), while the SAME pipeline runs
# fine interactively outside the sandbox — because inside the sandbox
# `$HOME` is a non-persistent tmpfs, so a `$HOME`-relative micromamba
# root / binary search resolves differently inside the job than in the
# submitting shell. `sbatch --export=ALL` does NOT fix this: it forwards
# env VAR values faithfully but does nothing to stop `$HOME` itself from
# resolving to a different instance once the job's own shell starts.
# See scripts/lib/env_config.sh for the full rationale. Activation is by
# absolute PREFIX, never by name, for the same reason.
#
# --env-name / $ENV_NAME previously selected a *named* env, activated
# via a $HOME/MAMBA_ROOT_PREFIX-relative lookup — exactly the
# sandbox-unsafe pattern this fix removes. Name-based activation is no
# longer supported: this wrapper always activates the single prefix
# pinned in scripts/env.local.conf. The flag is still parsed here (so
# old invocations don't hit an "unknown flag" error from `hexenium run
# "$@"` below) but only WARNS if given.
_env_name_cli=""
_filtered_args=()
_args=("$@")
_i=0
while [[ $_i -lt ${#_args[@]} ]]; do
    case "${_args[_i]}" in
        --env-name)
            _env_name_cli="${_args[$((_i+1))]:-}"
            _i=$((_i+2))
            ;;
        --env-name=*)
            _env_name_cli="${_args[_i]#--env-name=}"
            _i=$((_i+1))
            ;;
        *)
            _filtered_args+=("${_args[_i]}")
            _i=$((_i+1))
            ;;
    esac
done
set -- "${_filtered_args[@]+"${_filtered_args[@]}"}"
if [[ -n "$_env_name_cli" || -n "${ENV_NAME:-}" ]]; then
    echo "[submit] WARN: --env-name/\$ENV_NAME ('${_env_name_cli:-${ENV_NAME:-}}') is ignored." >&2
    echo "[submit]   Activation is now always by the absolute prefix pinned in" >&2
    echo "[submit]   scripts/env.local.conf. Re-run scripts/write-env-config.sh" >&2
    echo "[submit]   --env-prefix <path> to point at a different env." >&2
fi

# Under Slurm the script body runs from a per-job tmp copy at
# /var/tmp/slurmd/job<id>/slurm_script; BASH_SOURCE[0] there points at
# that copy, so `dirname` gives a path with no sibling lib/ dir. Prefer
# HEREG_SCRIPT_DIR (forwarded via --export by the pre-sbatch section)
# when it's set AND names a scripts/ dir that actually has our
# lib/env_config.sh; fall back to BASH_SOURCE for the case where this
# body runs interactively (e.g. `bash scripts/submit_he_registration.sh
# --dry-run` outside of a Slurm dispatch).
if [[ -n "${HEREG_SCRIPT_DIR:-}" && -f "$HEREG_SCRIPT_DIR/lib/env_config.sh" ]]; then
    SCRIPT_DIR="$HEREG_SCRIPT_DIR"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
# shellcheck source=lib/env_config.sh
source "$SCRIPT_DIR/lib/env_config.sh"

eval "$("$MICROMAMBA_BIN" shell hook --shell bash)"
if ! micromamba activate "$HEREG_ENV_PREFIX"; then
    echo "[submit] ERROR: micromamba activate failed for prefix: $HEREG_ENV_PREFIX" >&2
    echo "[submit]   Re-run scripts/write-env-config.sh to verify/regenerate env.local.conf." >&2
    exit 1
fi

# Verify python is the env's python, not the system one.
echo "[submit] env prefix: $HEREG_ENV_PREFIX"
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
