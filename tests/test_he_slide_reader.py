"""OMEAwareSlideReaderAdapter tests.

Verifies the adapter picks up ``PhysicalSizeX/Y`` from the source
OME-TIFF instead of HEST's hard-coded ``[0.25, 0.25, PIXEL_UNIT]``
default.

The full :class:`OMEAwareSlideReaderAdapter` class needs a real WSI
backend (HEST's ``wsi_factory``), so these tests exercise the pure
helpers ``_read_ome_physical_size`` + ``_apply_ome_physical_size``,
which contain all the OME-parsing logic. A ``SimpleNamespace`` stand-in
represents the ``MetaData`` object HEST would hand us.
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import tifffile
from valis.slide_io import MICRON_UNIT, PIXEL_UNIT

from hexenium._internal.he_slide_reader import (
    _apply_ome_physical_size,
    _read_ome_physical_size,
)


def _write_rgb_ome_tiff(path: Path, phys_x: float, phys_y: float) -> None:
    img = np.zeros((32, 32, 3), dtype=np.uint8)
    tifffile.imwrite(
        path, img, photometric="rgb",
        metadata={
            "axes": "YXS",
            "PhysicalSizeX": phys_x,
            "PhysicalSizeXUnit": "µm",
            "PhysicalSizeY": phys_y,
            "PhysicalSizeYUnit": "µm",
        },
        ome=True,
    )


def _fake_meta():
    return SimpleNamespace(
        pixel_physical_size_xyu=[0.25, 0.25, PIXEL_UNIT],
    )


# ---------------------------------------------------------------------------
# _read_ome_physical_size
# ---------------------------------------------------------------------------

def test_read_ome_physical_size_returns_ome_values(tmp_path: Path):
    p = tmp_path / "he.ome.tif"
    _write_rgb_ome_tiff(p, 0.1369, 0.1369)
    phys = _read_ome_physical_size(p)
    assert phys is not None
    px, py = phys
    assert abs(px - 0.1369) < 1e-9
    assert abs(py - 0.1369) < 1e-9


def test_read_ome_physical_size_returns_none_on_plain_tiff(tmp_path: Path):
    p = tmp_path / "flat.tif"
    tifffile.imwrite(p, np.zeros((16, 16, 3), dtype=np.uint8), photometric="rgb")
    assert _read_ome_physical_size(p) is None


def test_read_ome_physical_size_returns_none_on_nonexistent(tmp_path: Path):
    assert _read_ome_physical_size(tmp_path / "nope.ome.tif") is None


# ---------------------------------------------------------------------------
# _apply_ome_physical_size — the adapter's create_metadata delegate
# ---------------------------------------------------------------------------

def test_apply_ome_physical_size_overwrites_hest_default(tmp_path: Path):
    p = tmp_path / "he.ome.tif"
    _write_rgb_ome_tiff(p, 0.1369, 0.1369)
    meta = _fake_meta()

    ok = _apply_ome_physical_size(meta, str(p))

    assert ok is True
    assert abs(meta.pixel_physical_size_xyu[0] - 0.1369) < 1e-9
    assert abs(meta.pixel_physical_size_xyu[1] - 0.1369) < 1e-9
    assert meta.pixel_physical_size_xyu[2] == MICRON_UNIT
    assert meta.pixel_physical_size_xyu[2] != PIXEL_UNIT


def test_apply_ome_physical_size_falls_back_on_missing_ome(
    tmp_path: Path, caplog: pytest.LogCaptureFixture,
):
    """No PhysicalSize → warn loudly, leave HEST default in place."""
    p = tmp_path / "flat.tif"
    tifffile.imwrite(p, np.zeros((16, 16, 3), dtype=np.uint8), photometric="rgb")
    meta = _fake_meta()

    with caplog.at_level("WARNING"):
        ok = _apply_ome_physical_size(meta, str(p))

    assert ok is False
    assert meta.pixel_physical_size_xyu == [0.25, 0.25, PIXEL_UNIT]
    assert "could not read OME" in caplog.text


def test_apply_ome_physical_size_matches_vs_40x_scan(tmp_path: Path):
    """Regression test against the exact Olympus VS-series 40× values.

    ``PhysicalSizeX=0.1369049680068557``, ``PhysicalSizeY=0.136906127921726``
    are what the TMA1A_bladder H&E scans carry. Round-trip must
    preserve them to VALIS.
    """
    p = tmp_path / "he_vs40x.ome.tif"
    _write_rgb_ome_tiff(p, 0.1369049680068557, 0.136906127921726)
    meta = _fake_meta()

    ok = _apply_ome_physical_size(meta, str(p))

    assert ok is True
    assert meta.pixel_physical_size_xyu[0] == pytest.approx(0.1369049680068557)
    assert meta.pixel_physical_size_xyu[1] == pytest.approx(0.136906127921726)
    assert meta.pixel_physical_size_xyu[2] == MICRON_UNIT
