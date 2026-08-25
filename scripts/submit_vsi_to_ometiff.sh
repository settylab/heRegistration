#!/usr/bin/env bash
# ---------------------------------------------------------------------
# NOT a login shell (`bash -l`): a login shell implicitly sources
# ~/.bash_profile / ~/.bashrc via `$HOME` at process startup, before
# this script's own body runs — see scripts/lib/env_config.sh for why
# that class of bug is exactly what this pipeline can no longer risk.
# Env activation uses the resolved absolute values in
# scripts/env.local.conf instead — see the "Env activation" section
# below.
#
# Slurm submission wrapper: convert an Olympus VSI (or any BioFormats-
# readable WSI format) to pyramidal OME-TIFF using bioformats2raw +
# raw2ometiff. Independent of the hexenium pipeline — the produced
# OME-TIFF is then handed to `submit_hexenium_from_ometiff.sh` (or
# passed directly as `--he-slide` to `submit_he_registration.sh`).
#
# When to use this instead of hexenium's built-in he_preprocess:
#   - Your VSI file isn't understood by OpenSlide's Olympus reader
#     (`openslide.OpenSlide(...)` raises OpenSlideUnsupportedFormatError
#     or the pipeline crashes with vendor=None). BioFormats has much
#     wider VSI coverage than OpenSlide.
#   - You want to pre-convert once and reuse the OME-TIFF across many
#     hexenium invocations / samples / registration modes.
#
# Usage — call directly (no `sbatch` prefix); the script self-submits
# to slurm and routes its .out/.err under `<out-dir>/logs/`:
#
#   ./scripts/submit_vsi_to_ometiff.sh \
#       --vsi         /path/to/SAMPLE.vsi \
#       --out-dir     /path/to/output_dir \
#       --sample-id   SAMPLE \
#       [--env-name   heRegistration] \
#       [--series     <N>]           # keep only zarr series N; default:
#                                    # write all with --split and let the
#                                    # user pick <sample>_he_s2.ome.tiff
#                                    # (Olympus 40x lives in s2 usually)
#       [--max-workers 4]
#
# Env activation: always by the absolute prefix pinned in
# scripts/env.local.conf (scripts/write-env-config.sh) — never by name,
# never re-derived from $HOME. See scripts/lib/env_config.sh for why.
# --env-name / $ENV_NAME are still parsed but only WARN if given; they
# no longer select which env is activated.
#
# BFTOOLS_ROOT / LIBBLOSC_DIR are normally set via scripts/env.local.conf
# (scripts/write-env-config.sh --bftools-root/--libblosc-dir), which is
# sourced before this section and so takes precedence; exporting them
# directly still works as a one-off override when env.local.conf leaves
# them unset:
#   BFTOOLS_ROOT   — directory containing extracted Glencoe zips
#                    (bioformats2raw-*/bin/bioformats2raw +
#                    raw2ometiff-*/bin/raw2ometiff). Only consulted
#                    when the tools aren't on PATH in the env.
#   LIBBLOSC_DIR   — directory containing libblosc.so if the env
#                    doesn't ship one (raw2ometiff needs it via JNA).
#
# Prerequisites (install ONCE per env):
#   micromamba install -p <env-prefix> -c conda-forge \
#       bioformats2raw raw2ometiff c-blosc -y
#
# If conda-forge is unreachable, drop the Glencoe zips into a scratch
# dir and export BFTOOLS_ROOT + LIBBLOSC_DIR:
#   TOOL=/path/to/scratch/bftools
#   mkdir -p $TOOL && cd $TOOL
#   curl -L -O https://github.com/glencoesoftware/bioformats2raw/releases/download/v0.12.1/bioformats2raw-0.12.1.zip
#   curl -L -O https://github.com/glencoesoftware/raw2ometiff/releases/download/v0.9.0/raw2ometiff-0.9.0.zip
#   unzip -q bioformats2raw-0.12.1.zip && unzip -q raw2ometiff-0.9.0.zip
#   export BFTOOLS_ROOT=$TOOL
#   export LIBBLOSC_DIR=/path/to/env/with/libblosc/lib
# ---------------------------------------------------------------------
#SBATCH --job-name=vsi2ometiff
#SBATCH --partition=campus-new
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=6:00:00

set -euo pipefail

# ---------------------------------------------------------------------
# Launcher branch: self-submit under sbatch.
# ---------------------------------------------------------------------
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    VSI=""; OUT_DIR=""; SAMPLE_ID=""; ENV_NAME_CLI=""
    args=("$@")
    for ((i=0; i<${#args[@]}; i++)); do
        case "${args[i]}" in
            --vsi)         VSI="${args[i+1]:-}" ;;
            --out-dir)     OUT_DIR="${args[i+1]:-}" ;;
            --sample-id)   SAMPLE_ID="${args[i+1]:-}" ;;
            --env-name)    ENV_NAME_CLI="${args[i+1]:-}" ;;
        esac
    done

    if [[ -z "$VSI" || -z "$OUT_DIR" || -z "$SAMPLE_ID" ]]; then
        echo "ERROR: --vsi, --out-dir, and --sample-id are all required." >&2
        echo "Usage:" >&2
        echo "  ./scripts/submit_vsi_to_ometiff.sh \\" >&2
        echo "      --vsi <path> --out-dir <path> --sample-id <name> \\" >&2
        echo "      [--env-name <env>] [--series <N>] [--max-workers <N>]" >&2
        exit 2
    fi

    if [[ ! -f "$VSI" ]]; then
        echo "ERROR: --vsi does not exist: $VSI" >&2
        exit 2
    fi

    # Pre-create the output tree so sbatch's --output/--error can land.
    LOG_DIR="$OUT_DIR/logs"
    mkdir -p "$LOG_DIR"

    LEAF="${SAMPLE_ID}_vsi2ometiff_%j"
    JOBID=$(sbatch --parsable --hold \
        --output="$LOG_DIR/${LEAF}.out" \
        --error="$LOG_DIR/${LEAF}.err" \
        "$0" "$@")
    if [[ -z "$JOBID" ]]; then
        echo "ERROR: sbatch failed to return a job id" >&2
        exit 1
    fi
    scontrol release "$JOBID"
    echo "Submitted batch job $JOBID"
    echo "  logs: $LOG_DIR/${SAMPLE_ID}_vsi2ometiff_${JOBID}.{out,err}"
    exit 0
fi

# ---------------------------------------------------------------------
# Under-sbatch branch: activate env + convert.
# ---------------------------------------------------------------------

# Parse the CLI args again (self-submitted with the same argv).
VSI=""; OUT_DIR=""; SAMPLE_ID=""; ENV_NAME_CLI=""
SERIES=""; MAX_WORKERS="4"
args=("$@")
for ((i=0; i<${#args[@]}; i++)); do
    case "${args[i]}" in
        --vsi)         VSI="${args[i+1]:-}" ;;
        --out-dir)     OUT_DIR="${args[i+1]:-}" ;;
        --sample-id)   SAMPLE_ID="${args[i+1]:-}" ;;
        --env-name)    ENV_NAME_CLI="${args[i+1]:-}" ;;
        --series)      SERIES="${args[i+1]:-}" ;;
        --max-workers) MAX_WORKERS="${args[i+1]:-4}" ;;
    esac
done

if [[ -n "$ENV_NAME_CLI" || -n "${ENV_NAME:-}" ]]; then
    echo "[vsi2ometiff] WARN: --env-name/\$ENV_NAME ('${ENV_NAME_CLI:-${ENV_NAME:-}}') is ignored." >&2
    echo "[vsi2ometiff]   Activation is now always by the absolute prefix pinned in" >&2
    echo "[vsi2ometiff]   scripts/env.local.conf. Re-run scripts/write-env-config.sh" >&2
    echo "[vsi2ometiff]   --env-prefix <path> to point at a different env." >&2
fi

# ---- Env activation ---------------------------------------------------
# Resolved from scripts/env.local.conf, NEVER re-derived from `$HOME`
# here. See scripts/lib/env_config.sh for the full rationale and
# scripts/submit_he_registration.sh for the identical pattern.
# Activation is by absolute PREFIX, never by name.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/env_config.sh
source "$SCRIPT_DIR/lib/env_config.sh"

eval "$("$MICROMAMBA_BIN" shell hook --shell bash)"
if ! micromamba activate "$HEREG_ENV_PREFIX"; then
    echo "ERROR: micromamba activate failed for prefix: $HEREG_ENV_PREFIX" >&2
    echo "  Re-run scripts/write-env-config.sh to verify/regenerate env.local.conf." >&2
    exit 1
fi

echo "[vsi2ometiff] env activated: $HEREG_ENV_PREFIX (CONDA_PREFIX=$CONDA_PREFIX)"
echo "[vsi2ometiff] python:        $(command -v python)"

# ---- Locate the bftools binaries + libblosc -------------------------
# Prefer the env's PATH (canonical: mamba install ...). Fall back to
# BFTOOLS_ROOT-provided zips if the env is missing them.
B2R="$(command -v bioformats2raw || true)"
R2O="$(command -v raw2ometiff || true)"
if [[ -z "$B2R" || -z "$R2O" ]]; then
    if [[ -n "${BFTOOLS_ROOT:-}" ]]; then
        # shellcheck disable=SC2010
        B2R_DIR="$(ls -d "$BFTOOLS_ROOT"/bioformats2raw-*/bin 2>/dev/null | head -1)"
        R2O_DIR="$(ls -d "$BFTOOLS_ROOT"/raw2ometiff-*/bin 2>/dev/null | head -1)"
        if [[ -n "$B2R_DIR" && -n "$R2O_DIR" ]]; then
            export PATH="$B2R_DIR:$R2O_DIR:$PATH"
            B2R="$(command -v bioformats2raw)"
            R2O="$(command -v raw2ometiff)"
        fi
    fi
fi
if [[ -z "$B2R" || -z "$R2O" ]]; then
    echo "ERROR: bioformats2raw / raw2ometiff not on PATH." >&2
    echo "  Install once with:" >&2
    echo "    micromamba install -p $HEREG_ENV_PREFIX -c conda-forge \\" >&2
    echo "        bioformats2raw raw2ometiff c-blosc -y" >&2
    echo "  Or set BFTOOLS_ROOT via scripts/write-env-config.sh --bftools-root." >&2
    exit 1
fi
echo "[vsi2ometiff] bioformats2raw: $B2R ($($B2R --version 2>&1 | grep '^Version' | head -1))"
echo "[vsi2ometiff] raw2ometiff:    $R2O ($($R2O --version 2>&1 | grep '^Version' | head -1))"

# libblosc: raw2ometiff loads it via JNA. Env-installed c-blosc puts it
# under $CONDA_PREFIX/lib; user-provided override via LIBBLOSC_DIR.
if [[ -n "${LIBBLOSC_DIR:-}" ]]; then
    export LD_LIBRARY_PATH="$LIBBLOSC_DIR:${LD_LIBRARY_PATH:-}"
    echo "[vsi2ometiff] LIBBLOSC_DIR override: $LIBBLOSC_DIR"
elif [[ -f "$CONDA_PREFIX/lib/libblosc.so" ]]; then
    export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
    echo "[vsi2ometiff] libblosc.so: $CONDA_PREFIX/lib/libblosc.so"
else
    echo "[vsi2ometiff] WARN: no libblosc.so found in env or LIBBLOSC_DIR." >&2
    echo "[vsi2ometiff]  raw2ometiff may fail with 'Unable to load library blosc'." >&2
fi

echo "[vsi2ometiff] slurm:  jobid=${SLURM_JOB_ID} host=$(hostname)"
echo "[vsi2ometiff] input:  $VSI ($(stat -c%s "$VSI") bytes)"
echo "[vsi2ometiff] output: $OUT_DIR"
echo "[vsi2ometiff] sample: $SAMPLE_ID"

# ---- Convert --------------------------------------------------------
ZARR="$OUT_DIR/${SAMPLE_ID}.zarr"
OMETIFF_STEM="$OUT_DIR/${SAMPLE_ID}_he.ome.tif"

if [[ -e "$ZARR" ]]; then
    echo "[vsi2ometiff] existing zarr at $ZARR — removing (idempotent rerun)"
    rm -rf "$ZARR"
fi

echo "[vsi2ometiff] running: bioformats2raw --max-workers $MAX_WORKERS $VSI $ZARR"
t0=$(date +%s)
if [[ -n "$SERIES" ]]; then
    "$B2R" --max-workers "$MAX_WORKERS" --series "$SERIES" "$VSI" "$ZARR"
else
    "$B2R" --max-workers "$MAX_WORKERS" "$VSI" "$ZARR"
fi
echo "[vsi2ometiff] bioformats2raw done in $(( $(date +%s) - t0 ))s; zarr size = $(du -sh "$ZARR" | cut -f1)"

echo "[vsi2ometiff] running: raw2ometiff --rgb --split --max_workers $MAX_WORKERS $ZARR $OMETIFF_STEM"
t0=$(date +%s)
"$R2O" --rgb --split --max_workers "$MAX_WORKERS" "$ZARR" "$OMETIFF_STEM"
echo "[vsi2ometiff] raw2ometiff done in $(( $(date +%s) - t0 ))s"

echo "[vsi2ometiff] per-series OME-TIFF files:"
ls -lh "${OMETIFF_STEM%.ome.tif}"_s*.ome.tiff 2>/dev/null || true

# ---- Canonical output: symlink <sample>_he.ome.tif → full-res series ---
# The pipeline (submit_he_registration.sh's auto-VSI path AND
# hexenium.stages.he_preprocess.discover_existing_ometiff) expects the
# converted file at a stable name — the actual full-res image lives in
# one of the s<N>.ome.tiff files (usually s2 for Olympus 40x, but not
# guaranteed). Auto-detect by picking the largest-on-disk per-series
# file (a robust proxy for pixel count). --series <N> override wins.
if [[ -n "$SERIES" ]]; then
    PICK="${OMETIFF_STEM%.ome.tif}_s${SERIES}.ome.tiff"
    if [[ ! -f "$PICK" ]]; then
        echo "[vsi2ometiff] ERROR: --series $SERIES asked for $PICK" >&2
        echo "[vsi2ometiff]        but that file wasn't produced." >&2
        exit 1
    fi
    echo "[vsi2ometiff] --series $SERIES override → $PICK"
else
    # `ls -S` sorts by size descending.
    PICK=$(ls -S -1 "${OMETIFF_STEM%.ome.tif}"_s*.ome.tiff 2>/dev/null | head -1)
    if [[ -z "$PICK" ]]; then
        echo "[vsi2ometiff] ERROR: no per-series files produced." >&2
        exit 1
    fi
    echo "[vsi2ometiff] auto-detected full-res series (largest on disk):"
    echo "[vsi2ometiff]   $PICK ($(du -h "$PICK" | cut -f1))"
fi
# Relative symlink so the canonical path stays valid on move.
ln -sfT "$(basename "$PICK")" "$OMETIFF_STEM"
echo "[vsi2ometiff] canonical link: $OMETIFF_STEM -> $(readlink "$OMETIFF_STEM")"

# All series with the picked one starred — recorded in the log so a
# later triage can grep it without re-inspecting the OME-TIFF.
echo ""
echo "[vsi2ometiff] series metadata (* = picked):"
python -c "
import tifffile, glob, os
canonical = os.path.realpath('$OMETIFF_STEM')
for p in sorted(glob.glob('${OMETIFF_STEM%.ome.tif}_s*.ome.tiff')):
    with tifffile.TiffFile(p) as t:
        s = t.series[0]
    mark = '  *' if os.path.realpath(p) == canonical else '   '
    print(f'{mark} {os.path.basename(p)}: {s.shape} name={s.name!r}')
" 2>&1 || echo "[vsi2ometiff]   (tifffile unavailable; skipping metadata dump)"

echo "[vsi2ometiff] DONE at $(date)"
