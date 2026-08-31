"""Pre-registration OME-metadata validation for H&E + DAPI inputs.

Runs after paths are resolved but before VALIS starts, so a
metadata problem fails fast (or logs loudly, when strict mode is
off) rather than eating a 30-60 min VALIS run only to produce a
misaligned result.

The load-bearing check is **"did we actually read PhysicalSize
from the source OME, or did we fall back to HEST's hard-coded
0.25 default?"** — that is the failure class the OME-aware adapter
fix targets from a different angle, and the one we want to catch
one layer earlier here. Range and ratio checks are additional
sanity guards for obviously-suspicious metadata (e.g., a 100×
scale ratio would flag a µm-vs-cm unit mixup) — they are NOT the
primary defense against the H&E-0.25-hardcoded class of bug, which
is a small ratio (~1.2× typical) that this validator would only
detect via the "fallback used" signal.
"""
from __future__ import annotations

from pathlib import Path
from typing import Literal, NamedTuple, Optional

from hexenium._internal.he_slide_reader import _read_ome_physical_size


# Bounds on plausible optical-microscope pixel sizes (µm/px).
# 0.05 µm/px is roughly 40× immersion at the Nyquist limit for
# visible light; 1.0 µm/px is well past what any Xenium or standard
# H&E scanner produces. Values outside this range are almost
# certainly a units mixup (metres, cm) or corrupted metadata.
MIN_PLAUSIBLE_UM_PER_PX = 0.05
MAX_PLAUSIBLE_UM_PER_PX = 1.0

# Beyond this scale ratio the two images almost certainly disagree
# about their units; VALIS's feature matcher can compensate for
# ~2-4× scale differences, so 10× leaves an honest safety margin.
MAX_PLAUSIBLE_SCALE_RATIO = 10.0

# The value HEST's stock SlideReaderAdapter hard-codes when it can't
# read PhysicalSize. Not a valid Xenium H&E scale in practice — see
# he_slide_reader.py for the full story.
HEST_FALLBACK_UM_PER_PX = 0.25

# Tolerance for spotting "the HEST fallback was used" when we see
# PhysicalSize == 0.25 come out of the OME reader for a file that
# also has NO OME PhysicalSize. Real 0.25 µm/px scanners exist —
# the trigger is "OME metadata was missing AND the reader ended up
# at 0.25," not "the file reports 0.25."
_HEST_FALLBACK_EQ_TOL = 1e-9


SourceTag = Literal["ome_metadata", "hest_fallback"]


class MetadataValidationError(ValueError):
    """Raised by :func:`validate_registration_inputs` in strict mode when
    OME PhysicalSize is missing / unusable / clearly wrong.
    """


class PhysicalSize(NamedTuple):
    px: float          # PhysicalSizeX in µm/px
    py: float          # PhysicalSizeY in µm/px
    source: SourceTag  # where the value came from


def read_physical_size(path: Path | str) -> PhysicalSize:
    """Return the source's PhysicalSize with a source tag.

    ``ome_metadata`` when the source file carries usable
    PhysicalSize; ``hest_fallback`` when it does not (and HEST's
    hard-coded default would apply). Never raises for a missing
    source — the raising happens in :func:`check_physical_size`.
    """
    phys = _read_ome_physical_size(path)
    if phys is None:
        return PhysicalSize(
            HEST_FALLBACK_UM_PER_PX, HEST_FALLBACK_UM_PER_PX, "hest_fallback",
        )
    px, py = phys
    return PhysicalSize(float(px), float(py), "ome_metadata")


def check_physical_size(phys: PhysicalSize, label: str) -> None:
    """Fail-fast on a missing / out-of-range PhysicalSize.

    Called from :func:`validate_registration_inputs`; kept public so
    callers can also invoke it directly on files that don't reach the
    registration path (e.g., a preflight CLI verb).
    """
    if phys.source == "hest_fallback":
        raise MetadataValidationError(
            f"{label}: OME PhysicalSizeX/Y could not be read; HEST's stock "
            f"SlideReaderAdapter would fall back to {HEST_FALLBACK_UM_PER_PX} µm/px, "
            "which is almost certainly wrong for a Xenium H&E scan (VS-series "
            "40× native is 0.1369 µm/px). Fix: regenerate the source with a "
            "tool that preserves OME PhysicalSize (hexenium.he_preprocess, A1's "
            "dapi_crop, or bfconvert with --overwrite-ome-meta). Override: "
            "pass --strict-metadata-check false to force registration to run "
            "with the fallback scale (accuracy will be degraded)."
        )
    for axis, v in (("PhysicalSizeX", phys.px), ("PhysicalSizeY", phys.py)):
        if not (MIN_PLAUSIBLE_UM_PER_PX <= v <= MAX_PLAUSIBLE_UM_PER_PX):
            raise MetadataValidationError(
                f"{label}: {axis}={v} µm/px is outside the plausible range "
                f"[{MIN_PLAUSIBLE_UM_PER_PX}, {MAX_PLAUSIBLE_UM_PER_PX}] "
                "µm/px. This usually indicates a units mixup (metres, cm) "
                "or corrupted OME metadata."
            )


def check_ratio(he_phys: PhysicalSize, dapi_phys: PhysicalSize) -> None:
    """Fail on an obviously-suspicious H&E ↔ DAPI scale ratio.

    Catches unit-error footguns (µm-vs-cm, µm-vs-mm) — NOT the small
    (~1.2×) ratio produced by the HEST-0.25-hardcoded bug that the
    OME-aware adapter targets. Use :func:`check_physical_size`'s
    ``hest_fallback`` branch to catch the latter.
    """
    scales = (he_phys.px, he_phys.py, dapi_phys.px, dapi_phys.py)
    ratio = max(scales) / min(scales)
    if ratio > MAX_PLAUSIBLE_SCALE_RATIO:
        raise MetadataValidationError(
            f"H&E ↔ DAPI scale ratio {ratio:.2f}× exceeds the plausible "
            f"maximum {MAX_PLAUSIBLE_SCALE_RATIO}×. H&E "
            f"({he_phys.px}, {he_phys.py}) µm/px, DAPI "
            f"({dapi_phys.px}, {dapi_phys.py}) µm/px. This usually "
            "indicates one file's PhysicalSize is off by an order of "
            "magnitude (a units mixup — µm vs mm vs cm)."
        )


def validate_registration_inputs(
    he_path: Path | str,
    dapi_path: Path | str,
    *,
    strict: bool = True,
    log=None,
) -> tuple[PhysicalSize, PhysicalSize]:
    """Log resolved paths + PhysicalSizes, then run all checks.

    Returns ``(he_phys, dapi_phys)`` unconditionally (also when
    ``strict=False`` and a check would have failed — the failure is
    logged as a WARNING and the caller proceeds at their own risk).

    ``log`` is any callable taking a single string (defaults to
    :func:`print`). :mod:`hexenium.stages.registration` passes its
    stage-local ``log`` here so the emitted lines land in the slurm
    job's `.out` stream alongside the rest of the stage's output.
    """
    if log is None:
        log = print

    he_phys = read_physical_size(he_path)
    dapi_phys = read_physical_size(dapi_path)

    log(f"[pre-registration] H&E : path={he_path}")
    _log_physical_size(log, "H&E", he_phys)
    log(f"[pre-registration] DAPI: path={dapi_path}")
    _log_physical_size(log, "DAPI", dapi_phys)

    if (
        he_phys.source == "ome_metadata"
        and dapi_phys.source == "ome_metadata"
    ):
        ratio = max(he_phys.px, dapi_phys.px) / min(he_phys.px, dapi_phys.px)
        log(
            f"[pre-registration] scale ratio (H&E / DAPI) = "
            f"{he_phys.px / dapi_phys.px:.3f} (both from OME metadata; "
            f"max/min = {ratio:.2f})"
        )

    checks = (
        ("H&E", he_phys, lambda: check_physical_size(he_phys, "H&E")),
        ("DAPI", dapi_phys, lambda: check_physical_size(dapi_phys, "DAPI")),
        ("ratio", None, lambda: check_ratio(he_phys, dapi_phys)),
    )
    for _label, _phys, run_check in checks:
        try:
            run_check()
        except MetadataValidationError as e:
            if strict:
                raise
            log(f"[pre-registration] WARNING (strict off): {e}")

    return he_phys, dapi_phys


def _log_physical_size(log, label: str, phys: PhysicalSize) -> None:
    if phys.source == "ome_metadata":
        log(
            f"[pre-registration] {label}: PhysicalSizeX={phys.px:.6f} µm/px, "
            f"PhysicalSizeY={phys.py:.6f} µm/px (source: OME metadata)"
        )
        return
    log(
        f"[pre-registration] WARNING: {label} OME metadata missing "
        "PhysicalSizeX/Y; the OME-aware reader will fall back to HEST's "
        f"hard-coded [{HEST_FALLBACK_UM_PER_PX}, {HEST_FALLBACK_UM_PER_PX}, "
        "PIXEL_UNIT] default (source: HEST fallback). This is almost "
        "certainly wrong for a Xenium H&E scan (VS-series 40× native is "
        "0.1369 µm/px). Verify the source file, or regenerate the crop "
        "with a tool that preserves OME PhysicalSize."
    )


__all__ = [
    "MetadataValidationError",
    "PhysicalSize",
    "check_physical_size",
    "check_ratio",
    "read_physical_size",
    "validate_registration_inputs",
]
