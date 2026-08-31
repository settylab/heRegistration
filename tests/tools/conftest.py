"""Shared fixtures for tools/roi_subset tests.

Builds a tiny synthetic Xenium bundle on disk under a tmp path:
``experiment.xenium`` with a known ``pixel_size``, a 40×40 uint16
"DAPI" OME-TIFF (with an image big enough to catch off-by-one crop
bugs), ``cells.parquet`` with centroids on a grid, and long-format
``cell_boundaries.parquet`` / ``nucleus_boundaries.parquet``.

Fixture ROI files:

* ``roi.csv`` — QuPath format (microns) with two overlapping selections
  named ``Left`` and ``Center``.
* ``roi.geojson`` — QuPath-style GeoJSON of the same regions in DAPI px.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tifffile


PIXEL_SIZE = 0.5  # µm/px — tiny but easy arithmetic in tests
SLIDE_H_PX = 40
SLIDE_W_PX = 40
SLIDE_H_UM = SLIDE_H_PX * PIXEL_SIZE  # 20 µm
SLIDE_W_UM = SLIDE_W_PX * PIXEL_SIZE  # 20 µm


def _write_experiment_xenium(path: Path) -> None:
    with path.open("w") as f:
        json.dump({"pixel_size": PIXEL_SIZE, "num_cells": 16}, f)


def _write_dapi(path: Path) -> None:
    # Gradient image so a crop is visually distinct from the whole.
    y, x = np.indices((SLIDE_H_PX, SLIDE_W_PX), dtype=np.uint16)
    img = (y * SLIDE_W_PX + x).astype(np.uint16)
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(path, img, photometric="minisblack")


def _synthetic_cells() -> pd.DataFrame:
    """16 cells on a 4×4 grid inside the 20×20 µm slide."""
    rows = []
    for iy in range(4):
        for ix in range(4):
            cx = 2.5 + 5.0 * ix  # 2.5, 7.5, 12.5, 17.5
            cy = 2.5 + 5.0 * iy
            rows.append({
                "cell_id": f"cell_{iy}_{ix}",
                "x_centroid": cx,
                "y_centroid": cy,
                "transcript_counts": 100 + iy * 10 + ix,
                "cell_area": 4.0,
                "nucleus_area": 2.0,
            })
    return pd.DataFrame(rows)


def _boundary_long(cells: pd.DataFrame, radius: float) -> pd.DataFrame:
    """Emit 8-vertex square polygons around each cell centroid."""
    rows = []
    for _, c in cells.iterrows():
        for dx, dy in [(-1, -1), (0, -1), (1, -1), (1, 0),
                       (1, 1), (0, 1), (-1, 1), (-1, 0)]:
            rows.append({
                "cell_id": c["cell_id"],
                "vertex_x": c["x_centroid"] + dx * radius,
                "vertex_y": c["y_centroid"] + dy * radius,
                "label_id": 1,
            })
    return pd.DataFrame(rows)


@pytest.fixture
def synthetic_bundle(tmp_path: Path) -> Path:
    bundle = tmp_path / "TMA_synth"
    (bundle / "morphology_focus").mkdir(parents=True)
    _write_experiment_xenium(bundle / "experiment.xenium")
    _write_dapi(bundle / "morphology_focus" / "ch0000_dapi.ome.tif")

    cells = _synthetic_cells()
    cells.to_parquet(bundle / "cells.parquet", index=False)
    _boundary_long(cells, radius=1.5).to_parquet(
        bundle / "cell_boundaries.parquet", index=False)
    _boundary_long(cells, radius=1.0).to_parquet(
        bundle / "nucleus_boundaries.parquet", index=False)
    return bundle


# ---------------------------------------------------------------------------
# Synthetic AnnData fixtures for A1's h5ad path
# ---------------------------------------------------------------------------

def _make_synthetic_adata(cells: pd.DataFrame, obs_indexed_by_cell_id: bool = True):
    """Build a small AnnData mirroring cells.parquet's 16-cell grid.

    Preserves the same cell_ids so the A1 subset can cross-check the
    h5ad-derived selection against cells.parquet's selection.

    Populates enough surface (obs, var, obsm, uns, layers, .raw) that
    the preservation tests are meaningful.
    """
    import anndata as ad

    n_cells = len(cells)
    n_genes = 5
    rng = np.random.default_rng(42)
    X = rng.integers(0, 100, size=(n_cells, n_genes)).astype(np.float32)
    layer_norm = X / (X.sum(axis=1, keepdims=True) + 1e-9)

    obs = pd.DataFrame({
        "cell_id": cells["cell_id"].to_numpy(),
        "x_centroid": cells["x_centroid"].to_numpy(),
        "y_centroid": cells["y_centroid"].to_numpy(),
        "transcript_counts": cells["transcript_counts"].to_numpy(),
    })
    if obs_indexed_by_cell_id:
        obs.index = obs["cell_id"].astype(str)

    var = pd.DataFrame({
        "gene_symbol": [f"GENE_{i}" for i in range(n_genes)],
        "highly_variable": [True] * n_genes,
    }, index=[f"ENSG_{i:04d}" for i in range(n_genes)])

    obsm = {"spatial": obs[["x_centroid", "y_centroid"]].to_numpy(np.float64)}
    uns = {
        "sample_id": "TMA_synth",
        "run_id": "run_synth_0",
        "provenance": {"tool": "conftest-synth", "version": "0"},
    }

    adata = ad.AnnData(X=X, obs=obs, var=var, obsm=obsm, uns=uns,
                       layers={"normalized": layer_norm})
    adata.raw = adata.copy()
    return adata


@pytest.fixture
def synthetic_h5ad(tmp_path: Path, synthetic_bundle: Path) -> Path:
    """A minimal xenium-ranger-style h5ad on the same 16-cell grid."""
    cells = _synthetic_cells()
    adata = _make_synthetic_adata(cells)
    path = tmp_path / "TMA_synth_xenium_ranger.h5ad"
    adata.write_h5ad(str(path))
    return path


@pytest.fixture
def synthetic_h5ad_no_obsm(tmp_path: Path) -> Path:
    """Same 16-cell grid but coords via .obs[x_centroid/y_centroid] only.

    Exercises the `_resolve_spatial_coords` fallback branch when
    ``.obsm['spatial']`` is absent (proseg_raw / early xenium h5ads).
    """
    import anndata as ad
    cells = _synthetic_cells()
    adata = _make_synthetic_adata(cells)
    del adata.obsm["spatial"]
    path = tmp_path / "TMA_synth_no_obsm.h5ad"
    adata.write_h5ad(str(path))
    return path


@pytest.fixture
def roi_csv_um(tmp_path: Path) -> Path:
    """QuPath CSV with two overlapping polygons — Left and Center."""
    path = tmp_path / "rois.csv"
    lines = [
        "#Selection names: Left, Center\n",
        "#Areas (um^2): 100.0, 100.0\n",
        "Selection,X,Y,Class,Color\n",
        # Left: x=[0,10], y=[0,10] — covers cells at (2.5, 2.5), (7.5, 2.5), (2.5, 7.5), (7.5, 7.5)
        "Left,0.0,0.0,Unclassified,#c2c2c2\n",
        "Left,10.0,0.0,Unclassified,#c2c2c2\n",
        "Left,10.0,10.0,Unclassified,#c2c2c2\n",
        "Left,0.0,10.0,Unclassified,#c2c2c2\n",
        # Center: x=[5,15], y=[5,15] — overlaps Left at (7.5, 7.5)
        "Center,5.0,5.0,Unclassified,#c2c2c2\n",
        "Center,15.0,5.0,Unclassified,#c2c2c2\n",
        "Center,15.0,15.0,Unclassified,#c2c2c2\n",
        "Center,5.0,15.0,Unclassified,#c2c2c2\n",
    ]
    path.write_text("".join(lines))
    return path


@pytest.fixture
def roi_geojson_dapi_px(tmp_path: Path) -> Path:
    """QuPath-style GeoJSON of a rectangle in DAPI-pixel space.

    Same physical region as ``roi_csv_um``'s Left selection, but
    expressed in pixel coordinates (µm / 0.5 = ×2).
    """
    path = tmp_path / "roi_left.geojson"
    payload = {
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [[
                    [0.0, 0.0], [20.0, 0.0], [20.0, 20.0], [0.0, 20.0], [0.0, 0.0],
                ]],
            },
            "properties": {"name": "Left", "objectType": "annotation"},
        }],
    }
    path.write_text(json.dumps(payload))
    return path
