"""Bounded-memory DAPI crop tests.

Xenium DAPI OME-TIFFs are (4, ~52k, ~40k) uint16 — ~8 GB fully
materialised. A1's :func:`hexenium.tools.roi_subset._crop_dapi`
must decode only the tiles overlapping the ROI bbox, never the
whole page. These tests build a tiled synthetic tiff much larger
than the ROI and assert that (a) the crop is correct and (b) the
peak Python allocation during the crop is a small multiple of the
ROI, not of the full page.
"""
from __future__ import annotations

import gc
import tracemalloc
from pathlib import Path

import numpy as np
import pytest
import tifffile
from shapely.geometry import box

from hexenium.tools.roi_subset import _crop_dapi


SLIDE_PX = 4000                # 4000×4000 uint16 page = 32 MB nominal
TILE_PX = 512
PIXEL_SIZE = 1.0               # µm/px — keeps arithmetic 1:1
ROI_PX = 100                   # 100×100 px crop
ROI_ORIGIN_PX = 100


@pytest.fixture(scope="module")
def big_tiled_dapi(tmp_path_factory) -> Path:
    path = tmp_path_factory.mktemp("bmem") / "big_dapi.ome.tif"
    y, x = np.indices((SLIDE_PX, SLIDE_PX), dtype=np.uint32)
    img = (y * SLIDE_PX + x).astype(np.uint16)
    tifffile.imwrite(
        path, img, tile=(TILE_PX, TILE_PX), photometric="minisblack",
    )
    return path


def _roi_polygon():
    return box(
        ROI_ORIGIN_PX * PIXEL_SIZE,
        ROI_ORIGIN_PX * PIXEL_SIZE,
        (ROI_ORIGIN_PX + ROI_PX) * PIXEL_SIZE,
        (ROI_ORIGIN_PX + ROI_PX) * PIXEL_SIZE,
    )


def test_crop_dapi_correctness_on_tiled_page(big_tiled_dapi: Path):
    crop, bounds = _crop_dapi(big_tiled_dapi, _roi_polygon(), PIXEL_SIZE)
    assert crop.shape == (ROI_PX, ROI_PX)
    assert crop.dtype == np.uint16
    assert bounds == (
        ROI_ORIGIN_PX, ROI_ORIGIN_PX,
        ROI_ORIGIN_PX + ROI_PX, ROI_ORIGIN_PX + ROI_PX,
    )
    full = tifffile.imread(big_tiled_dapi)
    expected = full[
        ROI_ORIGIN_PX : ROI_ORIGIN_PX + ROI_PX,
        ROI_ORIGIN_PX : ROI_ORIGIN_PX + ROI_PX,
    ]
    assert np.array_equal(crop, expected)


def test_crop_dapi_returns_clamped_bounds_when_roi_overshoots(
    tmp_path: Path,
):
    """The returned bbox is the ACTUAL clamped slice, not the pre-clamped
    ROI request. When the ROI polygon extends past the DAPI's right/bottom
    edge, the returned ``(x0, y0, x1, y1)`` must match what was sliced.
    """
    slide = 40
    path = tmp_path / "small_dapi.ome.tif"
    y, x = np.indices((slide, slide), dtype=np.uint32)
    tifffile.imwrite(
        path, (y * slide + x).astype(np.uint16), photometric="minisblack",
    )
    # ROI extends 20 px past the slide on the right and bottom.
    poly = box(10.0, 10.0, 60.0, 60.0)  # at pixel_size 1.0
    crop, bounds = _crop_dapi(path, poly, 1.0)
    assert crop.shape == (30, 30)              # clamped: 40 - 10 = 30
    assert bounds == (10, 10, slide, slide)    # x1, y1 == clamped, not 60


def test_crop_dapi_peak_memory_is_bounded(big_tiled_dapi: Path):
    """Peak alloc during _crop_dapi is a small multiple of the ROI,
    not of the full page — otherwise we regressed to ``asarray()``.
    """
    full_page_bytes = SLIDE_PX * SLIDE_PX * 2  # uint16
    bound_bytes = full_page_bytes // 4          # 8 MB for a 32 MB page

    gc.collect()
    tracemalloc.start()
    tracemalloc.reset_peak()
    crop, _ = _crop_dapi(big_tiled_dapi, _roi_polygon(), PIXEL_SIZE)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    assert crop.nbytes == ROI_PX * ROI_PX * 2
    assert peak < bound_bytes, (
        f"peak={peak/1e6:.1f} MB exceeds bound {bound_bytes/1e6:.1f} MB; "
        f"full page is {full_page_bytes/1e6:.1f} MB — did _crop_dapi "
        "regress to asarray()?"
    )
