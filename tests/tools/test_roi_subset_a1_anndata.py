"""A1 AnnData ROI subset tests.

Covers the per-ROI subset path: reads the optional
``xenium_ranger_h5ad`` / ``proseg_purified_h5ad`` / ``proseg_raw_h5ad``
inputs in backed mode, materialises only the ROI-kept rows via
``to_memory()``, writes them to ``<roi_dir>/<stem>.h5ad``, and
preserves obs/var/obsm/uns/layers/.raw. Coord frame stays global.
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import yaml

from hexenium._internal.roi_config import Outputs, RoiSpec, SubsetConfig
from hexenium.tools.roi_subset import run_a1


def _cfg(
    bundle: Path,
    tmp_path: Path,
    specs: list[RoiSpec],
    *,
    xenium_ranger_h5ad: Path | None = None,
    proseg_purified_h5ad: Path | None = None,
    proseg_raw_h5ad: Path | None = None,
    outputs: Outputs | None = None,
) -> SubsetConfig:
    return SubsetConfig(
        mode="a1",
        xenium_bundle=bundle,
        output_dir=tmp_path / "subset_run",
        rois=specs,
        outputs=outputs or Outputs(),
        xenium_ranger_h5ad=xenium_ranger_h5ad,
        proseg_purified_h5ad=proseg_purified_h5ad,
        proseg_raw_h5ad=proseg_raw_h5ad,
    )


def test_h5ad_subset_agrees_with_cells_parquet(
    synthetic_bundle: Path, synthetic_h5ad: Path, roi_csv_um: Path, tmp_path: Path,
):
    """The h5ad subset for the xenium_ranger source keeps the same cells
    that the cells.parquet subset keeps — same ROI, same predicate.
    """
    specs = [RoiSpec(name="Left", path=roi_csv_um, frame="microns",
                     source_tool="qupath_csv", selection="Left")]
    [roi_dir] = run_a1(_cfg(
        synthetic_bundle, tmp_path, specs,
        xenium_ranger_h5ad=synthetic_h5ad,
    ))

    parquet_ids = set(
        pd.read_parquet(roi_dir / "cells.parquet")["cell_id"].astype(str)
    )
    subset = ad.read_h5ad(str(roi_dir / "xenium_ranger.h5ad"))
    h5ad_ids = set(subset.obs["cell_id"].astype(str))
    assert h5ad_ids == parquet_ids
    assert subset.n_obs == 4  # Left ROI covers 4 cells on the 4×4 grid


def test_h5ad_subset_preserves_obs_var_obsm_uns_layers_raw(
    synthetic_bundle: Path, synthetic_h5ad: Path, roi_csv_um: Path, tmp_path: Path,
):
    specs = [RoiSpec(name="Left", path=roi_csv_um, frame="microns",
                     source_tool="qupath_csv", selection="Left")]
    [roi_dir] = run_a1(_cfg(
        synthetic_bundle, tmp_path, specs,
        xenium_ranger_h5ad=synthetic_h5ad,
    ))

    src = ad.read_h5ad(str(synthetic_h5ad))
    sub = ad.read_h5ad(str(roi_dir / "xenium_ranger.h5ad"))

    # obs — same columns, cell-count matches ROI membership
    assert list(sub.obs.columns) == list(src.obs.columns)
    assert sub.n_obs == 4
    # var — untouched
    assert list(sub.var.columns) == list(src.var.columns)
    assert sub.var.index.equals(src.var.index)
    assert sub.n_vars == src.n_vars
    # obsm — same key set, first-axis subset
    assert set(sub.obsm.keys()) == set(src.obsm.keys())
    assert sub.obsm["spatial"].shape == (sub.n_obs, 2)
    # uns — pass-through
    assert sub.uns["sample_id"] == src.uns["sample_id"]
    assert sub.uns["run_id"] == src.uns["run_id"]
    assert sub.uns["provenance"]["tool"] == "conftest-synth"
    # layers — sliced consistently
    assert set(sub.layers.keys()) == {"normalized"}
    assert sub.layers["normalized"].shape == (sub.n_obs, sub.n_vars)
    # .raw — preserved AND row-sliced symmetrically with .X. Backed
    # mode's biggest known compatibility caveat; the .to_memory() call
    # inside _subset_anndata must propagate the row mask through
    # adata.raw. If it doesn't, sub.raw.X.shape[0] would still equal
    # src.n_obs (the whole population).
    assert sub.raw is not None
    assert sub.raw.n_vars == src.raw.n_vars
    assert sub.raw.n_obs == sub.n_obs
    assert sub.raw.X.shape == (sub.n_obs, sub.raw.n_vars)
    assert list(sub.raw.obs_names) == list(sub.obs_names)


def test_h5ad_subset_raw_values_match_source_rows(
    synthetic_bundle: Path, synthetic_h5ad: Path, roi_csv_um: Path, tmp_path: Path,
):
    """Backed-mode ``.raw`` slice materialises the exact source rows.

    Shape checks in the preservation test catch structural
    regressions; this test catches a subtler case where the backed
    slice materialised the wrong rows (e.g., an off-by-one on the
    HDF5-backed index) — values must line up cell-for-cell.
    """
    specs = [RoiSpec(name="Left", path=roi_csv_um, frame="microns",
                     source_tool="qupath_csv", selection="Left")]
    [roi_dir] = run_a1(_cfg(
        synthetic_bundle, tmp_path, specs,
        xenium_ranger_h5ad=synthetic_h5ad,
    ))

    src = ad.read_h5ad(str(synthetic_h5ad))
    sub = ad.read_h5ad(str(roi_dir / "xenium_ranger.h5ad"))
    src_kept = src[src.obs_names.isin(list(sub.obs_names)), :].copy()

    np.testing.assert_array_equal(np.asarray(sub.X), np.asarray(src_kept.X))
    np.testing.assert_array_equal(
        np.asarray(sub.raw.X), np.asarray(src_kept.raw.X),
    )


def test_h5ad_subset_keeps_global_coords(
    synthetic_bundle: Path, synthetic_h5ad: Path, roi_csv_um: Path, tmp_path: Path,
):
    """Default anndata_geometry_frame='global' → no shift. The four
    Left-cell centroids in the subset match the source verbatim.
    """
    specs = [RoiSpec(name="Left", path=roi_csv_um, frame="microns",
                     source_tool="qupath_csv", selection="Left")]
    [roi_dir] = run_a1(_cfg(
        synthetic_bundle, tmp_path, specs,
        xenium_ranger_h5ad=synthetic_h5ad,
    ))

    src = ad.read_h5ad(str(synthetic_h5ad))
    sub = ad.read_h5ad(str(roi_dir / "xenium_ranger.h5ad"))
    kept_ids = set(sub.obs["cell_id"].astype(str))
    src_kept = src[src.obs["cell_id"].astype(str).isin(kept_ids), :].copy()

    # Coord columns and obsm both unchanged.
    np.testing.assert_array_equal(
        sub.obs[["x_centroid", "y_centroid"]].to_numpy(),
        src_kept.obs[["x_centroid", "y_centroid"]].to_numpy(),
    )
    np.testing.assert_array_equal(sub.obsm["spatial"], src_kept.obsm["spatial"])


def test_h5ad_subset_falls_back_to_obs_when_obsm_missing(
    synthetic_bundle: Path, synthetic_h5ad_no_obsm: Path,
    roi_csv_um: Path, tmp_path: Path,
):
    """The reused celltyping resolver handles h5ads with no obsm['spatial']
    by falling back to .obs[x_centroid/y_centroid]. A1 must inherit that
    behavior so proseg-style inputs work.
    """
    specs = [RoiSpec(name="Left", path=roi_csv_um, frame="microns",
                     source_tool="qupath_csv", selection="Left")]
    [roi_dir] = run_a1(_cfg(
        synthetic_bundle, tmp_path, specs,
        proseg_raw_h5ad=synthetic_h5ad_no_obsm,
    ))

    sub = ad.read_h5ad(str(roi_dir / "proseg_raw.h5ad"))
    assert sub.n_obs == 4


def test_missing_optional_h5ad_inputs_do_not_break_a1(
    synthetic_bundle: Path, roi_csv_um: Path, tmp_path: Path,
):
    """A1 must remain runnable when the three h5ad fields stay unset."""
    specs = [RoiSpec(name="Left", path=roi_csv_um, frame="microns",
                     source_tool="qupath_csv", selection="Left")]
    [roi_dir] = run_a1(_cfg(synthetic_bundle, tmp_path, specs))  # no h5ads

    # Parquet outputs still emit; no h5ad files appear.
    assert (roi_dir / "cells.parquet").exists()
    assert not (roi_dir / "xenium_ranger.h5ad").exists()
    assert not (roi_dir / "proseg_purified.h5ad").exists()
    assert not (roi_dir / "proseg_raw.h5ad").exists()


def test_manifest_records_anndata_subset_counts(
    synthetic_bundle: Path, synthetic_h5ad: Path, roi_csv_um: Path, tmp_path: Path,
):
    specs = [RoiSpec(name="Left", path=roi_csv_um, frame="microns",
                     source_tool="qupath_csv", selection="Left")]
    [roi_dir] = run_a1(_cfg(
        synthetic_bundle, tmp_path, specs,
        xenium_ranger_h5ad=synthetic_h5ad,
    ))
    manifest = yaml.safe_load((roi_dir / "manifest.yaml").read_text())
    entries = manifest["anndata_subsets"]
    assert "xenium_ranger" in entries
    assert entries["xenium_ranger"]["n_kept"] == 4
    assert entries["xenium_ranger"]["n_total"] == 16
    assert entries["xenium_ranger"]["output"].endswith(
        "/xenium_ranger.h5ad"
    )


def test_anndata_geometry_frame_local_raises_notimplemented(
    synthetic_bundle: Path, synthetic_h5ad: Path, roi_csv_um: Path, tmp_path: Path,
):
    """local coord-shift for AnnData is deferred; explicit not-yet-implemented."""
    specs = [RoiSpec(name="Left", path=roi_csv_um, frame="microns",
                     source_tool="qupath_csv", selection="Left")]
    outputs = Outputs(anndata_geometry_frame="local")
    with pytest.raises(NotImplementedError, match="anndata_geometry_frame"):
        run_a1(_cfg(
            synthetic_bundle, tmp_path, specs,
            xenium_ranger_h5ad=synthetic_h5ad,
            outputs=outputs,
        ))


