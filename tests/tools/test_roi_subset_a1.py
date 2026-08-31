"""A1 tests: coord-frame, multi-ROI + overlap, crop origin, sparse ROI.

The tests rely on the synthetic fixtures in :mod:`tests.tools.conftest`;
they exercise the actual :func:`hexenium.tools.roi_subset.run_a1`
against small on-disk bundles + ROI files.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tifffile
import yaml

from hexenium._internal.roi_config import RoiSpec, SubsetConfig, Outputs
from hexenium._internal.coord_frames import (
    convert_polygon,
    read_pixel_size_morph,
    um_to_dapi_px,
)
from hexenium._internal.roi_readers import (
    read_qupath_csv,
    read_geojson,
)
from hexenium.tools.roi_subset import run_a1, _slug


# ---------------------------------------------------------------------------
# Coord-frame + reader unit tests
# ---------------------------------------------------------------------------

def test_read_pixel_size_falls_back_when_missing(tmp_path: Path):
    assert read_pixel_size_morph(tmp_path / "nope.xenium") == 0.2125


def test_qupath_csv_reader_parses_first_selection(roi_csv_um: Path):
    poly = read_qupath_csv(roi_csv_um)
    # Left selection default: bbox (0,0)-(10,10) in microns
    assert poly.bounds == (0.0, 0.0, 10.0, 10.0)


def test_qupath_csv_reader_picks_by_selection_name(roi_csv_um: Path):
    poly = read_qupath_csv(roi_csv_um, selection="Center")
    assert poly.bounds == (5.0, 5.0, 15.0, 15.0)


def test_qupath_csv_reader_raises_on_unknown_selection(roi_csv_um: Path):
    with pytest.raises(ValueError, match="Nope"):
        read_qupath_csv(roi_csv_um, selection="Nope")


def test_geojson_reader_parses_polygon(roi_geojson_dapi_px: Path):
    poly = read_geojson(roi_geojson_dapi_px)
    assert poly.bounds == (0.0, 0.0, 20.0, 20.0)


def test_um_to_dapi_px_scale():
    from shapely.geometry import box
    poly_um = box(0, 0, 10, 10)
    poly_px = um_to_dapi_px(poly_um, pixel_size_morph=0.5)
    assert poly_px.bounds == (0.0, 0.0, 20.0, 20.0)


def test_convert_polygon_rejects_he_px():
    from shapely.geometry import box
    with pytest.raises(ValueError, match="he_px"):
        convert_polygon(box(0, 0, 1, 1), "he_px", "microns", 0.2125)


# ---------------------------------------------------------------------------
# Config validation
# ---------------------------------------------------------------------------

def test_config_requires_explicit_frame(tmp_path: Path):
    with pytest.raises(ValueError, match="frame must be"):
        SubsetConfig(
            mode="a1",
            xenium_bundle=tmp_path,
            output_dir=tmp_path / "out",
            rois=[RoiSpec(name="X", path=tmp_path / "roi.csv",
                          frame="whatever", source_tool="qupath_csv")],
        ).validate()


def test_config_a1_rejects_he_px_frame(tmp_path: Path):
    with pytest.raises(ValueError, match="he_px"):
        SubsetConfig(
            mode="a1",
            xenium_bundle=tmp_path,
            output_dir=tmp_path / "out",
            rois=[RoiSpec(name="X", path=tmp_path / "r.geojson",
                          frame="he_px", source_tool="qupath_geojson")],
        ).validate()


# ---------------------------------------------------------------------------
# A1 end-to-end
# ---------------------------------------------------------------------------

def _cfg_for(bundle: Path, tmp_path: Path, roi_specs: list[RoiSpec],
             geometry_frame: str = "global",
             outputs: Outputs | None = None) -> SubsetConfig:
    return SubsetConfig(
        mode="a1",
        xenium_bundle=bundle,
        output_dir=tmp_path / "subset_run",
        rois=roi_specs,
        outputs=outputs or Outputs(geometry_frame=geometry_frame),
    )


def test_a1_multi_roi_with_overlap(
    synthetic_bundle: Path, roi_csv_um: Path, tmp_path: Path,
):
    """Both selections in the same CSV → two output dirs.

    A cell inside both selections is present in BOTH per-ROI outputs
    (no auto-dedup).
    """
    specs = [
        RoiSpec(name="Left",   path=roi_csv_um, frame="microns",
                source_tool="qupath_csv", selection="Left"),
        RoiSpec(name="Center", path=roi_csv_um, frame="microns",
                source_tool="qupath_csv", selection="Center"),
    ]
    dirs = run_a1(_cfg_for(synthetic_bundle, tmp_path, specs))

    assert [d.name for d in dirs] == ["left", "center"]

    left_cells  = pd.read_parquet(dirs[0] / "cells.parquet")
    center_cells = pd.read_parquet(dirs[1] / "cells.parquet")

    # Left bbox (0,0)-(10,10) covers the four cells with centroids
    # < 10 in both axes: (2.5,2.5), (7.5,2.5), (2.5,7.5), (7.5,7.5).
    assert set(left_cells["cell_id"]) == {
        "cell_0_0", "cell_0_1", "cell_1_0", "cell_1_1",
    }
    # Center bbox (5,5)-(15,15) covers (7.5,7.5), (12.5,7.5), (7.5,12.5), (12.5,12.5).
    assert set(center_cells["cell_id"]) == {
        "cell_1_1", "cell_1_2", "cell_2_1", "cell_2_2",
    }
    # cell_1_1 appears in BOTH — no dedup.
    assert "cell_1_1" in set(left_cells["cell_id"])
    assert "cell_1_1" in set(center_cells["cell_id"])


def test_a1_boundary_parquets_subset_by_cell_id(
    synthetic_bundle: Path, roi_csv_um: Path, tmp_path: Path,
):
    specs = [RoiSpec(name="Left", path=roi_csv_um, frame="microns",
                     source_tool="qupath_csv", selection="Left")]
    [roi_dir] = run_a1(_cfg_for(synthetic_bundle, tmp_path, specs))

    cb = pd.read_parquet(roi_dir / "cell_boundaries.parquet")
    nb = pd.read_parquet(roi_dir / "nucleus_boundaries.parquet")

    # 4 cells × 8 vertices each — schema preserved.
    assert list(cb.columns) == ["cell_id", "vertex_x", "vertex_y", "label_id"]
    assert list(nb.columns) == ["cell_id", "vertex_x", "vertex_y", "label_id"]
    assert cb["cell_id"].nunique() == 4
    assert nb["cell_id"].nunique() == 4
    assert len(cb) == 32
    assert len(nb) == 32


def test_a1_experiment_xenium_copied_verbatim(
    synthetic_bundle: Path, roi_csv_um: Path, tmp_path: Path,
):
    specs = [RoiSpec(name="Left", path=roi_csv_um, frame="microns",
                     source_tool="qupath_csv", selection="Left")]
    [roi_dir] = run_a1(_cfg_for(synthetic_bundle, tmp_path, specs))

    src = (synthetic_bundle / "experiment.xenium").read_text()
    dst = (roi_dir / "experiment.xenium").read_text()
    assert src == dst


def test_a1_dapi_crop_origin_and_shape(
    synthetic_bundle: Path, roi_csv_um: Path, tmp_path: Path,
):
    """Left roi bbox (0,0)-(10,10) µm at pixel_size 0.5 = (0,0)-(20,20) px.

    Crop shape should be (20, 20) uint16 and pixel (0, 0) of the crop
    should equal pixel (0, 0) of the source (gradient image, so
    ``value == y * W + x`` and origin px value == 0).
    """
    specs = [RoiSpec(name="Left", path=roi_csv_um, frame="microns",
                     source_tool="qupath_csv", selection="Left")]
    [roi_dir] = run_a1(_cfg_for(synthetic_bundle, tmp_path, specs))

    crop = tifffile.imread(roi_dir / "dapi_subset.ome.tif")
    assert crop.shape == (20, 20)
    assert crop.dtype == np.uint16
    assert crop[0, 0] == 0


def test_a1_manifest_records_crop_origin_and_frame(
    synthetic_bundle: Path, roi_csv_um: Path, tmp_path: Path,
):
    specs = [RoiSpec(name="Left", path=roi_csv_um, frame="microns",
                     source_tool="qupath_csv", selection="Left")]
    [roi_dir] = run_a1(_cfg_for(synthetic_bundle, tmp_path, specs))

    manifest = yaml.safe_load((roi_dir / "manifest.yaml").read_text())
    assert manifest["roi_name"] == "Left"
    assert manifest["geometry_frame"] == "global"
    assert manifest["pixel_size_morph_um_per_px"] == 0.5
    assert manifest["dapi_crop"]["origin_dapi_px"] == [0, 0]
    assert manifest["dapi_crop"]["origin_dapi_um"] == [0.0, 0.0]
    assert manifest["dapi_crop"]["size_dapi_px"] == [20, 20]
    assert manifest["roi_bbox_microns"] == [0.0, 0.0, 10.0, 10.0]
    assert manifest["n_cells_kept"] == 4
    assert manifest["n_cells_total"] == 16


def test_a1_geometry_frame_local_shifts_roi_geojson(
    synthetic_bundle: Path, roi_csv_um: Path, tmp_path: Path,
):
    """When geometry_frame='local', the emitted roi.geojson is shifted so
    (0, 0) matches the crop origin. For the Left roi (crop origin
    (0, 0) µm), this is a no-op — pick a ROI with a non-zero crop
    origin to make the shift observable.
    """
    specs = [RoiSpec(name="Center", path=roi_csv_um, frame="microns",
                     source_tool="qupath_csv", selection="Center")]
    outputs = Outputs(geometry_frame="local")
    [roi_dir] = run_a1(
        _cfg_for(synthetic_bundle, tmp_path, specs, outputs=outputs),
    )

    import geopandas as gpd
    gdf = gpd.read_file(roi_dir / "roi.geojson")
    # Center in global microns: (5,5)-(15,15). Crop origin px = (5/0.5, 5/0.5) = (10, 10),
    # crop origin µm = (5, 5). Local frame: bbox becomes (0,0)-(10,10).
    assert gdf.geometry.iloc[0].bounds == (0.0, 0.0, 10.0, 10.0)


def test_a1_geojson_dapi_px_roi_matches_csv_microns_roi(
    synthetic_bundle: Path, roi_csv_um: Path, roi_geojson_dapi_px: Path,
    tmp_path: Path,
):
    """A QuPath GeoJSON (DAPI px) and a QuPath CSV (µm) representing the
    same physical region produce identical cell subsets.
    """
    csv_dirs = run_a1(_cfg_for(
        synthetic_bundle, tmp_path / "csv",
        [RoiSpec(name="Left", path=roi_csv_um, frame="microns",
                 source_tool="qupath_csv", selection="Left")],
    ))
    gj_dirs = run_a1(_cfg_for(
        synthetic_bundle, tmp_path / "gj",
        [RoiSpec(name="Left", path=roi_geojson_dapi_px, frame="dapi_px",
                 source_tool="qupath_geojson")],
    ))
    csv_cells = pd.read_parquet(csv_dirs[0] / "cells.parquet")
    gj_cells  = pd.read_parquet(gj_dirs[0]  / "cells.parquet")
    assert set(csv_cells["cell_id"]) == set(gj_cells["cell_id"])


def test_a1_slug_normalizes_names():
    assert _slug("Selection 1") == "selection_1"
    assert _slug("Left") == "left"
    assert _slug("") == "roi"
    assert _slug("A/B c") == "a_b_c"
