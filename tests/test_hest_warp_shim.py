"""Tests for the vendored HEST warp shim.

Covers the pixel_size_morph resolution + reader-plumbing that works
around HEST v1.2.0's `warp_gdf_valis` bug (calls `read_gdf(shapes)`
with no reader_kwargs, so `XeniumParquetCellReader.pixel_size_morph`
stays None and the reader crashes on the first divide).

Skips if hest / geopandas aren't installed — matches the pattern in
`test_stage_imports.py`.
"""
from __future__ import annotations

import json
from unittest import mock

import pytest


hest_seg_readers = pytest.importorskip("hest.io.seg_readers")


def test_resolve_pixel_size_morph_reads_experiment_xenium(tmp_path):
    """When `experiment.xenium` lives next to the parquet, its
    `pixel_size` field wins over the fallback default."""
    from hexenium._internal.hest_warp_shim import _resolve_pixel_size_morph

    bundle = tmp_path
    (bundle / "experiment.xenium").write_text(
        json.dumps({"pixel_size": 0.2125, "z_step_size": 3.0})
    )
    parquet = bundle / "cell_boundaries.parquet"
    parquet.write_bytes(b"")  # doesn't need to be a real parquet; only
                              # its .parent is consulted here.

    assert _resolve_pixel_size_morph(str(parquet)) == 0.2125


def test_resolve_pixel_size_morph_uses_metadata_value_when_nondefault(tmp_path):
    """If a future Xenium slide reports a different `pixel_size`, we
    honour it rather than silently clamping to 0.2125."""
    from hexenium._internal.hest_warp_shim import _resolve_pixel_size_morph

    bundle = tmp_path
    (bundle / "experiment.xenium").write_text(
        json.dumps({"pixel_size": 0.5})
    )
    parquet = bundle / "nucleus_boundaries.parquet"
    parquet.write_bytes(b"")

    assert _resolve_pixel_size_morph(str(parquet)) == 0.5


def test_resolve_pixel_size_morph_falls_back_when_metadata_missing(tmp_path):
    """No `experiment.xenium` next to the parquet → fall back to the
    canonical Xenium morphology default (0.2125)."""
    from hexenium._internal.hest_warp_shim import (
        _XENIUM_PIXEL_SIZE_MORPH_DEFAULT,
        _resolve_pixel_size_morph,
    )

    parquet = tmp_path / "cell_boundaries.parquet"
    parquet.write_bytes(b"")

    assert _resolve_pixel_size_morph(str(parquet)) == _XENIUM_PIXEL_SIZE_MORPH_DEFAULT


def test_resolve_pixel_size_morph_falls_back_on_bad_json(tmp_path):
    """Corrupt `experiment.xenium` (or one missing the `pixel_size`
    field) must not crash the shim — fall back gracefully."""
    from hexenium._internal.hest_warp_shim import (
        _XENIUM_PIXEL_SIZE_MORPH_DEFAULT,
        _resolve_pixel_size_morph,
    )

    bundle = tmp_path
    (bundle / "experiment.xenium").write_text("not valid json {")
    parquet = bundle / "cell_boundaries.parquet"
    parquet.write_bytes(b"")

    assert _resolve_pixel_size_morph(str(parquet)) == _XENIUM_PIXEL_SIZE_MORPH_DEFAULT

    (bundle / "experiment.xenium").write_text(json.dumps({"z_step_size": 3.0}))
    assert _resolve_pixel_size_morph(str(parquet)) == _XENIUM_PIXEL_SIZE_MORPH_DEFAULT


def test_load_xenium_boundaries_forwards_pixel_size_to_reader(tmp_path):
    """Regression guard against the original v1.2.0 bug: the shim's
    loader MUST pass `pixel_size_morph` into `read_gdf`'s reader
    kwargs. Otherwise HEST v1.2.0's
    `XeniumParquetCellReader.pixel_size_morph` stays at its
    constructor default of `None` and the next `df['vertex_x'] / None`
    crashes."""
    from hexenium._internal.hest_warp_shim import _load_xenium_boundaries

    bundle = tmp_path
    (bundle / "experiment.xenium").write_text(json.dumps({"pixel_size": 0.2125}))
    parquet = bundle / "cell_boundaries.parquet"
    parquet.write_bytes(b"")

    with mock.patch("hest.io.seg_readers.read_gdf") as mocked:
        mocked.return_value = "sentinel-gdf"
        result = _load_xenium_boundaries(str(parquet))

    assert result == "sentinel-gdf"
    mocked.assert_called_once()
    args, kwargs = mocked.call_args
    # First positional: the parquet path.
    assert args[0] == str(parquet)
    # reader_kwargs must carry pixel_size_morph — the whole point.
    assert kwargs.get("reader_kwargs") == {"pixel_size_morph": 0.2125}
