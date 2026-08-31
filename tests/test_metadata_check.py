"""Pre-registration OME-metadata validator tests.

Verifies the two failure modes we most want to catch before a
30-60 min VALIS run:

* HEST-fallback detection — the source has no readable OME
  PhysicalSize, so the reader would fall back to HEST's hard-coded
  0.25 µm/px default. This is the failure class the OME-aware
  adapter fix targets from the reader side; the validator is the
  cheap pre-flight version of the same check.
* Obviously-suspicious scale ratios (a units mixup) — a separate
  sanity guard, NOT the primary defense against the 0.25-hardcoded
  bug (that ratio is small).

Also covers the fail-fast vs. warn-only knob and the successful
"both metadata are fine" path.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import tifffile

from hexenium._internal.metadata_check import (
    HEST_FALLBACK_UM_PER_PX,
    MAX_PLAUSIBLE_SCALE_RATIO,
    MAX_PLAUSIBLE_UM_PER_PX,
    MIN_PLAUSIBLE_UM_PER_PX,
    MetadataValidationError,
    check_physical_size,
    check_ratio,
    read_physical_size,
    validate_registration_inputs,
)


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

def _write_ome_tif(path: Path, phys_x: float, phys_y: float,
                   rgb: bool = True) -> None:
    if rgb:
        img = np.zeros((32, 32, 3), dtype=np.uint8)
        photometric = "rgb"
        axes = "YXS"
    else:
        img = np.zeros((32, 32), dtype=np.uint16)
        photometric = "minisblack"
        axes = "YX"
    tifffile.imwrite(
        path, img, photometric=photometric,
        metadata={
            "axes": axes,
            "PhysicalSizeX": phys_x, "PhysicalSizeXUnit": "µm",
            "PhysicalSizeY": phys_y, "PhysicalSizeYUnit": "µm",
        }, ome=True,
    )


def _write_plain_tif(path: Path) -> None:
    tifffile.imwrite(
        path, np.zeros((16, 16, 3), dtype=np.uint8), photometric="rgb",
    )


class _LogCapture:
    def __init__(self):
        self.lines: list[str] = []

    def __call__(self, s: str) -> None:
        self.lines.append(s)

    @property
    def text(self) -> str:
        return "\n".join(self.lines)


# ---------------------------------------------------------------------------
# read_physical_size — pure helper, no policy
# ---------------------------------------------------------------------------

def test_read_physical_size_returns_ome_when_present(tmp_path: Path):
    p = tmp_path / "he.ome.tif"
    _write_ome_tif(p, 0.1369, 0.1369)
    phys = read_physical_size(p)
    assert phys.source == "ome_metadata"
    assert phys.px == pytest.approx(0.1369)
    assert phys.py == pytest.approx(0.1369)


def test_read_physical_size_returns_hest_fallback_when_missing(tmp_path: Path):
    p = tmp_path / "flat.tif"
    _write_plain_tif(p)
    phys = read_physical_size(p)
    assert phys.source == "hest_fallback"
    assert phys.px == HEST_FALLBACK_UM_PER_PX
    assert phys.py == HEST_FALLBACK_UM_PER_PX


# ---------------------------------------------------------------------------
# check_physical_size
# ---------------------------------------------------------------------------

def test_check_physical_size_passes_for_typical_values(tmp_path: Path):
    p = tmp_path / "he.ome.tif"
    _write_ome_tif(p, 0.1369, 0.1369)
    check_physical_size(read_physical_size(p), "H&E")  # no raise


def test_check_physical_size_fails_when_hest_fallback_used(tmp_path: Path):
    p = tmp_path / "flat.tif"
    _write_plain_tif(p)
    with pytest.raises(MetadataValidationError, match="OME PhysicalSizeX/Y could not be read"):
        check_physical_size(read_physical_size(p), "H&E")


def test_check_physical_size_fails_for_out_of_range_scale(tmp_path: Path):
    # 10 µm/px — well past MAX_PLAUSIBLE_UM_PER_PX.
    p = tmp_path / "he.ome.tif"
    _write_ome_tif(p, MAX_PLAUSIBLE_UM_PER_PX * 10, MAX_PLAUSIBLE_UM_PER_PX * 10)
    with pytest.raises(MetadataValidationError, match="outside the plausible range"):
        check_physical_size(read_physical_size(p), "H&E")


def test_check_physical_size_fails_for_tiny_scale(tmp_path: Path):
    p = tmp_path / "he.ome.tif"
    _write_ome_tif(p, MIN_PLAUSIBLE_UM_PER_PX / 10, MIN_PLAUSIBLE_UM_PER_PX / 10)
    with pytest.raises(MetadataValidationError, match="outside the plausible range"):
        check_physical_size(read_physical_size(p), "H&E")


# ---------------------------------------------------------------------------
# check_ratio
# ---------------------------------------------------------------------------

def test_check_ratio_passes_on_realistic_TMA1A_case(tmp_path: Path):
    """H&E 0.1369, DAPI 0.2125 → 1.55× ratio, well under the 10× guard."""
    he = tmp_path / "he.ome.tif"; _write_ome_tif(he, 0.1369, 0.1369)
    da = tmp_path / "da.ome.tif"; _write_ome_tif(da, 0.2125, 0.2125, rgb=False)
    check_ratio(read_physical_size(he), read_physical_size(da))  # no raise


def test_check_ratio_flags_100x_units_mixup(tmp_path: Path):
    """H&E 0.13 µm/px vs DAPI mislabeled 13 (units off by ×100) → flagged."""
    he = tmp_path / "he.ome.tif"; _write_ome_tif(he, 0.13, 0.13)
    da = tmp_path / "da.ome.tif"; _write_ome_tif(da, 13.0, 13.0, rgb=False)
    with pytest.raises(MetadataValidationError, match="scale ratio.*exceeds"):
        check_ratio(read_physical_size(he), read_physical_size(da))


def test_check_ratio_does_not_flag_TMA1A_hest_fallback_ratio(tmp_path: Path):
    """The 0.25-vs-0.1369 case that motivated all this work has a ~1.2×
    ratio — well under the 10× threshold. The ratio check is NOT the
    defense against the HEST-hardcoded bug; the hest_fallback branch of
    check_physical_size is.
    """
    he_phys = read_physical_size(_ome_pair(tmp_path, "he_wrong.ome.tif",
                                          HEST_FALLBACK_UM_PER_PX)[0])
    da_phys = read_physical_size(_ome_pair(tmp_path, "da.ome.tif", 0.2125)[0])
    check_ratio(he_phys, da_phys)  # no raise — 1.18× ratio is not suspicious


def _ome_pair(tmp_path: Path, name: str, s: float) -> tuple[Path, float]:
    p = tmp_path / name
    _write_ome_tif(p, s, s)
    return p, s


# ---------------------------------------------------------------------------
# validate_registration_inputs — orchestration + strict knob
# ---------------------------------------------------------------------------

def test_validate_logs_paths_and_scales_on_success(tmp_path: Path):
    he = tmp_path / "he.ome.tif"; _write_ome_tif(he, 0.1369, 0.1369)
    da = tmp_path / "da.ome.tif"; _write_ome_tif(da, 0.2125, 0.2125, rgb=False)
    logs = _LogCapture()

    he_phys, da_phys = validate_registration_inputs(he, da, log=logs)

    assert he_phys.source == "ome_metadata"
    assert da_phys.source == "ome_metadata"
    text = logs.text
    assert "[pre-registration] H&E : path=" in text
    assert "[pre-registration] DAPI: path=" in text
    assert "PhysicalSizeX=0.136900" in text
    assert "PhysicalSizeX=0.212500" in text
    assert "(source: OME metadata)" in text
    assert "scale ratio" in text


def test_validate_strict_raises_on_hest_fallback(tmp_path: Path):
    he = tmp_path / "he_bad.tif"; _write_plain_tif(he)  # no OME
    da = tmp_path / "da.ome.tif"; _write_ome_tif(da, 0.2125, 0.2125, rgb=False)
    logs = _LogCapture()

    with pytest.raises(MetadataValidationError, match="HEST"):
        validate_registration_inputs(he, da, strict=True, log=logs)

    assert "WARNING: H&E OME metadata missing" in logs.text
    assert "source: HEST fallback" in logs.text


def test_validate_nonstrict_warns_but_returns(tmp_path: Path):
    """Escape hatch: strict=False downgrades the failures to WARNs."""
    he = tmp_path / "he_bad.tif"; _write_plain_tif(he)
    da = tmp_path / "da.ome.tif"; _write_ome_tif(da, 0.2125, 0.2125, rgb=False)
    logs = _LogCapture()

    he_phys, da_phys = validate_registration_inputs(
        he, da, strict=False, log=logs,
    )
    assert he_phys.source == "hest_fallback"
    assert da_phys.source == "ome_metadata"
    assert "WARNING (strict off):" in logs.text


def test_validate_strict_raises_on_extreme_ratio(tmp_path: Path):
    """Extreme ratio inside per-file range → ratio check fires.

    Both values individually pass ``check_physical_size`` (0.08 and
    0.9 are inside [0.05, 1.0]), so the failure comes from
    ``check_ratio``: 0.9 / 0.08 = 11.25 > 10.
    """
    he = tmp_path / "he.ome.tif"; _write_ome_tif(he, 0.08, 0.08)
    da = tmp_path / "da.ome.tif"; _write_ome_tif(da, 0.9, 0.9, rgb=False)
    logs = _LogCapture()
    with pytest.raises(MetadataValidationError, match="scale ratio"):
        validate_registration_inputs(he, da, strict=True, log=logs)


def test_validate_ratio_boundary_at_10x_does_not_fail(tmp_path: Path):
    """MAX_PLAUSIBLE_SCALE_RATIO = 10 — a 10× ratio should pass; only
    ratios strictly greater than the threshold trip the check.
    """
    he = tmp_path / "he.ome.tif"; _write_ome_tif(he, 0.05, 0.05)
    da = tmp_path / "da.ome.tif"; _write_ome_tif(da, 0.5, 0.5, rgb=False)
    logs = _LogCapture()
    assert MAX_PLAUSIBLE_SCALE_RATIO == 10.0
    validate_registration_inputs(he, da, strict=True, log=logs)  # no raise
