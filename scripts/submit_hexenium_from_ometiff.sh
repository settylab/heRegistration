#!/bin/bash -l
# ---------------------------------------------------------------------
# Slurm submission wrapper: run hexenium (register + warp + celltype +
# viz) on a pre-converted OME-TIFF — the "Path B" second step.
# Independent of `submit_vsi_to_ometiff.sh`; there's no shared state
# beyond the OME-TIFF file itself, so you can rerun this stage as many
# times as you like (different modes, different celltype sources, etc.)
# without re-running the conversion.
#
# When to use this instead of `submit_he_registration.sh` with
# `--stages he_preprocess register warp celltype viz`:
#   - Your H&E is already OME-TIFF (either you converted a VSI via
#     `submit_vsi_to_ometiff.sh` first, or you have a
#     Trident/QuPath/bfconvert output).
#   - You want a script whose args make it OBVIOUS the OME-TIFF
#     step was completed and hexenium starts at `register`.
#
# Under the hood this is a thin dispatcher to
# `submit_he_registration.sh` with `--stages` forced to
# `register warp celltype viz`. If you need a different stage subset
# (e.g. just register+warp), pass `--stages` explicitly — it wins.
#
# Usage — call directly (no `sbatch` prefix); this script hands off
# to `submit_he_registration.sh` which self-submits under sbatch:
#
#   ./scripts/submit_hexenium_from_ometiff.sh \
#       --he-ometiff  /path/to/converted/SAMPLE_he_s2.ome.tiff \
#       --sample-id   SAMPLE \
#       --run-id      demo_v1 \
#       --output-root /path/to/runs \
#       --xenium-bundle /path/to/output-XETG... \
#       --dapi-path   /path/to/output-XETG.../morphology_focus/morphology_focus_0000.ome.tif \
#       [--proseg-purified-h5ad /path/to/proseg_purified.h5ad] \
#       [--env-name   heRegistration] \
#       [--stages     register warp celltype viz] \
#       [--set-default-on-success] \
#       [<any-other-hexenium-flag>...]
#
# The `--he-ometiff` flag is the ONE thing this script names
# differently from the underlying hexenium CLI; it just gets
# forwarded as `--he-slide`. Every OTHER flag passes through
# unchanged to `hexenium run`, so anything the CLI accepts
# (`--mode`, `--warp-run-id`, `--celltype-col`, `--viz-dpi`, …) works
# here too.
#
# Environment variables (fallbacks; CLI flags win):
#   ENV_NAME     — conda/micromamba env with hexenium installed.
#                  Default: heRegistration.
#   OUTPUT_ROOT  — fallback for --output-root.
# ---------------------------------------------------------------------

set -euo pipefail

# ---------------------------------------------------------------------
# Parse --he-ometiff (special) + rewrite it to --he-slide. Everything
# else passes through to submit_he_registration.sh untouched.
# ---------------------------------------------------------------------
_ometiff=""
_stages_set=0
_forwarded_args=()
_args=("$@")
_i=0
while [[ $_i -lt ${#_args[@]} ]]; do
    case "${_args[_i]}" in
        --he-ometiff)
            _ometiff="${_args[$((_i+1))]:-}"
            _i=$((_i+2))
            ;;
        --he-ometiff=*)
            _ometiff="${_args[_i]#--he-ometiff=}"
            _i=$((_i+1))
            ;;
        --stages)
            # Consume all nargs="+" values so we can detect "user gave
            # their own --stages". Then re-forward them unchanged.
            _forwarded_args+=("--stages")
            _i=$((_i+1))
            while [[ $_i -lt ${#_args[@]} && "${_args[$_i]}" != --* ]]; do
                _forwarded_args+=("${_args[$_i]}")
                _i=$((_i+1))
            done
            _stages_set=1
            ;;
        *)
            _forwarded_args+=("${_args[_i]}")
            _i=$((_i+1))
            ;;
    esac
done

if [[ -z "$_ometiff" ]]; then
    echo "ERROR: --he-ometiff is required." >&2
    echo "Usage:" >&2
    echo "  ./scripts/submit_hexenium_from_ometiff.sh \\" >&2
    echo "      --he-ometiff <path.ome.tif> --sample-id <name> \\" >&2
    echo "      --run-id <id> --output-root <path> \\" >&2
    echo "      --xenium-bundle <path> [--dapi-path <path>] \\" >&2
    echo "      [--proseg-purified-h5ad <path>] [--env-name <env>] \\" >&2
    echo "      [--stages register warp celltype viz]" >&2
    exit 2
fi

if [[ ! -f "$_ometiff" ]]; then
    echo "ERROR: --he-ometiff does not exist: $_ometiff" >&2
    exit 2
fi

# Warn if the file extension isn't one hexenium's is_vsi_input() will
# short-circuit through. `.ome.tif` / `.ome.tiff` / `.tif` all pass;
# `.vsi` would cause the he_preprocess stage to try to CONVERT the
# file. Since the whole point of this script is "input is already
# OME-TIFF", refuse.
case "${_ometiff,,}" in
    *.vsi)
        echo "ERROR: --he-ometiff looks like a .vsi file: $_ometiff" >&2
        echo "  This script is for pre-converted OME-TIFFs. Convert first" >&2
        echo "  via ./scripts/submit_vsi_to_ometiff.sh, then pass the" >&2
        echo "  resulting *.ome.tiff (usually the s2 series for Olympus" >&2
        echo "  40x scans) as --he-ometiff." >&2
        exit 2
        ;;
esac

# Default --stages when user didn't set one: skip he_preprocess
# (input is already OME-TIFF; no conversion needed).
if [[ $_stages_set -eq 0 ]]; then
    _forwarded_args+=("--stages" "register" "warp" "celltype" "viz")
fi

# Locate the sibling launcher script (same dir as this one).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LAUNCHER="$SCRIPT_DIR/submit_he_registration.sh"
if [[ ! -x "$LAUNCHER" ]]; then
    echo "ERROR: sibling launcher not found or not executable: $LAUNCHER" >&2
    exit 1
fi

echo "[hexenium-from-ometiff] dispatching to $LAUNCHER"
echo "[hexenium-from-ometiff] --he-slide will be:  $_ometiff"
echo "[hexenium-from-ometiff] forwarded args:      ${_forwarded_args[*]}"
exec "$LAUNCHER" --he-slide "$_ometiff" "${_forwarded_args[@]}"
