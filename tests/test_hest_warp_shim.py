"""Tests for the vendored HEST warp shim.

Covers the pixel_size_morph resolution + reader-plumbing that works
around HEST v1.2.0's `warp_gdf_valis` bug (calls `read_gdf(shapes)`
with no reader_kwargs, so `XeniumParquetCellReader.pixel_size_morph`
stays None and the reader crashes on the first divide), and the
`__null_dask_index__` marker the shim now writes so
`celltyping._read_warped_gdf` can decode the parquet.

Skips if hest / geopandas aren't installed — matches the pattern in
`test_stage_imports.py`.
"""
from __future__ import annotations

import json
from unittest import mock

import pytest


hest_seg_readers = pytest.importorskip("hest.io.seg_readers")
gpd = pytest.importorskip("geopandas")
pd = pytest.importorskip("pandas")
_shapely = pytest.importorskip("shapely")


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


def _fake_warped_gdf() -> gpd.GeoDataFrame:
    """Build a small GeoDataFrame in the exact shape v1.2.0's
    ``warp_gdf_valis`` returns: cell IDs live in an UNNAMED index,
    only a ``geometry`` column."""
    from shapely.geometry import Polygon

    polys = [
        Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
        Polygon([(2, 2), (3, 2), (3, 3), (2, 3)]),
        Polygon([(5, 5), (6, 5), (6, 6), (5, 6)]),
    ]
    gdf = gpd.GeoDataFrame(geometry=polys)
    gdf.index = ["aaaaegdo-1", "aaaagjph-1", "aaaakpdi-1"]
    assert gdf.index.name is None  # matches v1.2.0 warp_gdf_valis output
    return gdf


def test_shim_parquet_has_null_dask_index_marker(tmp_path):
    """Regression guard: the shim MUST write the warp parquet with
    the anonymous cell-id index named ``__null_dask_index__``. Under
    HEST main's dask path, ``to_parquet`` writes that sentinel
    automatically; under our v1.2.0 non-dask shim, we set it
    explicitly so ``celltyping._read_warped_gdf``'s downstream
    ``rename(columns={"__null_dask_index__": "xenium_cell_id"})``
    hits its target."""
    from hexenium._internal.hest_warp_shim import (
        _write_warped_parquet_and_geojson,
    )

    _write_warped_parquet_and_geojson(
        _fake_warped_gdf(),
        save_dir=str(tmp_path),
        stem="he_cell_seg",
        save_parquet=True,
        save_geojson=False,
    )

    parquet = tmp_path / "he_cell_seg.parquet"
    assert parquet.exists()
    read = pd.read_parquet(parquet)
    assert read.index.name == "__null_dask_index__", (
        f"expected the anonymous cell-id index to be named "
        f"'__null_dask_index__' so celltyping._read_warped_gdf's "
        f"rename target lands; got {read.index.name!r}"
    )
    # After reset_index (what _read_warped_gdf does), the marker
    # column MUST be present so its subsequent rename hits.
    assert "__null_dask_index__" in read.reset_index().columns


def test_celltyping_reader_accepts_shim_parquet(tmp_path):
    """End-to-end regression against the KeyError that surfaced on
    Tracy's `65099653` celltype run: read a shim-written parquet
    with the exact reader `celltyping.run_celltyping` uses. The
    reader must succeed and expose a `xenium_cell_id` column with
    the original cell IDs."""
    from hexenium._internal.hest_warp_shim import (
        _write_warped_parquet_and_geojson,
    )

    _write_warped_parquet_and_geojson(
        _fake_warped_gdf(),
        save_dir=str(tmp_path),
        stem="he_cell_seg",
        save_parquet=True,
        save_geojson=False,
    )

    try:
        from hexenium.stages.celltyping import _read_warped_gdf
    except ModuleNotFoundError as e:  # celltyping pulls heavier deps
        pytest.skip(f"celltyping module unavailable: {e.name}")

    gdf = _read_warped_gdf(tmp_path / "he_cell_seg.parquet")
    assert "xenium_cell_id" in gdf.columns
    assert list(gdf["xenium_cell_id"]) == [
        "aaaaegdo-1",
        "aaaagjph-1",
        "aaaakpdi-1",
    ]
    # Sanity: geometry column reconstructed via from_wkb.
    assert gdf.geometry.iloc[0].is_valid
