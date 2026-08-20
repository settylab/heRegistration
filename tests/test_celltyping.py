"""Unit tests for the inlined celltype NN mapping (proseg_purified →
xenium), auto-detect helpers with the shape-checked id-col precedence
(guards against Proseg's int64 index).

Run:
    pytest tests/test_celltyping.py -v
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import anndata as ad
from shapely.geometry import Polygon

from hexenium.stages.celltyping import (
    UNLABELED,
    _CELLTYPE_COL_CANDIDATES,
    _SPATIAL_COORD_OBS_PAIRS,
    _XENIUM_ID_COL_CANDIDATES,
    _celltype_from_nn,
    _celltype_from_xenium_ranger_direct,
    _get_xenium_ids,
    _resolve_celltype_col,
    _resolve_spatial_coords,
    _resolve_xenium_id_col,
    run_celltyping,
)


# ---------------------------------------------------------------------
# _resolve_celltype_col — auto-detect precedence.
# ---------------------------------------------------------------------
class TestResolveCelltypeCol:
    def test_prefers_celltype_over_first_type(self):
        cols = ["celltype", "first_type", "primary_cell_type"]
        assert _resolve_celltype_col(cols, "auto") == "celltype"

    def test_falls_back_to_first_type(self):
        # Today's proseg_purified.h5ad has first_type + primary_cell_type
        # but no `celltype` (that lands after the roadmap's writeback stage).
        cols = ["first_type", "primary_cell_type"]
        assert _resolve_celltype_col(cols, "auto") == "first_type"

    def test_falls_back_to_primary_cell_type(self):
        cols = ["primary_cell_type", "spot_class"]
        assert _resolve_celltype_col(cols, "auto") == "primary_cell_type"

    def test_falls_back_to_celltype_updated(self):
        cols = ["celltype_updated", "spot_class"]
        assert _resolve_celltype_col(cols, "auto") == "celltype_updated"

    def test_none_treated_as_auto(self):
        cols = ["first_type"]
        assert _resolve_celltype_col(cols, None) == "first_type"

    def test_explicit_col_returned_verbatim(self):
        cols = ["custom_label", "first_type"]
        assert _resolve_celltype_col(cols, "custom_label") == "custom_label"

    def test_explicit_col_missing_fails_loud(self):
        cols = ["first_type"]
        with pytest.raises(KeyError, match="missing celltype column 'nope'"):
            _resolve_celltype_col(cols, "nope")

    def test_no_candidates_fails_loud(self):
        cols = ["x", "y"]
        with pytest.raises(KeyError, match="could not auto-detect"):
            _resolve_celltype_col(cols, "auto")

    def test_candidate_order_is_authoritative(self):
        # Even if `celltype_updated` comes first alphabetically, the
        # precedence order in _CELLTYPE_COL_CANDIDATES wins.
        cols = ["celltype_updated", "celltype", "first_type"]
        assert _resolve_celltype_col(cols, "auto") == "celltype"
        assert _CELLTYPE_COL_CANDIDATES[0] == "celltype"


# ---------------------------------------------------------------------
# _resolve_xenium_id_col — THE shape-check trap from the design-pass
# skeptic. Must NEVER pick Proseg's int64 "cell_id" over xenium UUIDs.
# ---------------------------------------------------------------------
class TestResolveXeniumIdCol:
    @staticmethod
    def _obs(cols: dict, index=None) -> pd.DataFrame:
        df = pd.DataFrame(cols)
        if index is not None:
            df.index = index
        return df

    def test_uuid_shaped_cell_id_accepted(self):
        # Post-writeback xenium_ranger.h5ad has UUID-shaped cell_id.
        obs = self._obs({"cell_id": ["aaaaafep-1", "lkkgeihi-1", "knjdobdk-1"]})
        assert _resolve_xenium_id_col(obs, "auto") == "cell_id"

    def test_proseg_int64_cell_id_rejected_falls_through(self):
        # THE bug the design-pass skeptic caught: on proseg_purified.h5ad,
        # cell_id is int64 [0, 1, 2, ...]. Auto-detect must skip it.
        obs = self._obs({
            "cell_id": [0, 1, 2, 3, 5],                              # Proseg int64 index
            "xenium_cell_id_nn": ["aaaaafep-1", "lkkgeihi-1",
                                  "knjdobdk-1", "degchhml-1", "goidknaf-1"],
        })
        # Must pick xenium_cell_id_nn (UUID-shaped) not cell_id (int64).
        assert _resolve_xenium_id_col(obs, "auto") == "xenium_cell_id_nn"

    def test_int_cast_to_string_still_rejected(self):
        # If someone .astype(str)-ed a Proseg cell_id before landing it
        # in the h5ad, the values look like "0", "1", "2" — still fails
        # the UUID regex.
        obs = self._obs({
            "cell_id": ["0", "1", "2", "3", "5"],
            "xenium_cell_id_nn": ["aaaaafep-1"] * 5,
        })
        assert _resolve_xenium_id_col(obs, "auto") == "xenium_cell_id_nn"

    def test_precedence_prefers_xenium_cell_id(self):
        # Both xenium_cell_id AND cell_id UUID-shaped: precedence wins.
        obs = self._obs({
            "xenium_cell_id": ["aaaaafep-1", "lkkgeihi-1", "knjdobdk-1"],
            "cell_id":        ["mmmmmfep-1", "abcdefgh-1", "ijklmnop-1"],
        })
        assert _resolve_xenium_id_col(obs, "auto") == "xenium_cell_id"

    def test_falls_back_to_index_when_columns_absent(self):
        # xenium_ranger.h5ad might carry UUIDs on the index (not a column).
        obs = self._obs(
            {"x_centroid": [1.0, 2.0, 3.0]},
            index=["aaaaafep-1", "lkkgeihi-1", "knjdobdk-1"],
        )
        assert _resolve_xenium_id_col(obs, "auto") == "__index__"

    def test_explicit_col_bypasses_shape_check(self):
        # If the user KNOWS they want cell_id, they can force it.
        obs = self._obs({"cell_id": [0, 1, 2, 3, 5]})
        assert _resolve_xenium_id_col(obs, "cell_id") == "cell_id"

    def test_explicit_missing_col_fails_loud(self):
        obs = self._obs({"cell_id": ["aaaaafep-1"]})
        with pytest.raises(KeyError, match="missing id column"):
            _resolve_xenium_id_col(obs, "nonexistent")

    def test_explicit_index_marker(self):
        obs = self._obs({"x": [1.0]}, index=["aaaaafep-1"])
        assert _resolve_xenium_id_col(obs, "__index__") == "__index__"

    def test_no_uuid_shaped_source_fails_loud(self):
        # Nothing looks like a xenium UUID → fail with actionable message.
        obs = self._obs(
            {"x_centroid": [1.0, 2.0]},
            index=["0", "1"],
        )
        with pytest.raises(KeyError, match="Pass --id-col explicitly"):
            _resolve_xenium_id_col(obs, "auto")


# ---------------------------------------------------------------------
# _get_xenium_ids — id extraction after resolution.
# ---------------------------------------------------------------------
class TestGetXeniumIds:
    def test_named_column(self):
        obs = pd.DataFrame({"xenium_cell_id": ["aaaaafep-1", "lkkgeihi-1"]})
        ids = _get_xenium_ids(obs, "xenium_cell_id")
        assert ids.tolist() == ["aaaaafep-1", "lkkgeihi-1"]

    def test_index_sentinel(self):
        obs = pd.DataFrame(
            {"x": [1.0, 2.0]},
            index=["aaaaafep-1", "lkkgeihi-1"],
        )
        ids = _get_xenium_ids(obs, "__index__")
        assert ids.tolist() == ["aaaaafep-1", "lkkgeihi-1"]

    def test_whitespace_stripped(self):
        obs = pd.DataFrame({"cell_id": ["  aaaaafep-1  ", "lkkgeihi-1"]})
        ids = _get_xenium_ids(obs, "cell_id")
        assert ids.tolist() == ["aaaaafep-1", "lkkgeihi-1"]


# ---------------------------------------------------------------------
# _resolve_spatial_coords — obsm['spatial'] → centroid_x/y → x_centroid/y_centroid
# → x/y fallback chain. Tracy's second #15 failure (slurm-64031325.err):
# an h5ad with coords in .obsm['spatial'] and NONE of the obs pairs.
# ---------------------------------------------------------------------
class TestResolveSpatialCoords:
    @staticmethod
    def _make_adata(*, obs_cols=None, obsm=None, n=3):
        obs = pd.DataFrame(obs_cols or {}, index=[str(i) for i in range(n)])
        # AnnData requires at least one obs column to infer n; pad if empty.
        if obs.shape[1] == 0:
            obs = pd.DataFrame({"_dummy": [0] * n}, index=obs.index)
        adata = ad.AnnData(
            X=np.zeros((n, 2), dtype=np.float32),
            obs=obs,
        )
        for k, v in (obsm or {}).items():
            adata.obsm[k] = np.asarray(v, dtype=float)
        return adata

    def test_obsm_spatial_only(self):
        # Xenium/scanpy convention — the case Tracy hit in slurm-64031325.
        adata = self._make_adata(
            obs_cols={"celltype": ["A", "B", "C"]},
            obsm={"spatial": [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]},
            n=3,
        )
        xy = _resolve_spatial_coords(adata)
        assert xy.shape == (3, 2)
        assert xy[0].tolist() == [1.0, 2.0]
        assert xy[2].tolist() == [5.0, 6.0]

    def test_obs_centroid_x_y_pair(self):
        # proseg_to_anndata convention.
        adata = self._make_adata(
            obs_cols={"centroid_x": [10.0, 20.0], "centroid_y": [11.0, 21.0]},
            n=2,
        )
        xy = _resolve_spatial_coords(adata)
        assert xy.tolist() == [[10.0, 11.0], [20.0, 21.0]]

    def test_obs_x_centroid_y_centroid_pair(self):
        # xenium_ranger convention (the previous hard-coded xenium side).
        adata = self._make_adata(
            obs_cols={"x_centroid": [7.0, 8.0], "y_centroid": [70.0, 80.0]},
            n=2,
        )
        xy = _resolve_spatial_coords(adata)
        assert xy.tolist() == [[7.0, 70.0], [8.0, 80.0]]

    def test_obs_x_y_generic_pair(self):
        # Generic fallback — the ORIGINAL default before the first #15 fix.
        adata = self._make_adata(
            obs_cols={"x": [0.5, 1.5], "y": [0.6, 1.6]},
            n=2,
        )
        xy = _resolve_spatial_coords(adata)
        assert xy.tolist() == [[0.5, 0.6], [1.5, 1.6]]

    def test_none_present_raises_value_error(self):
        adata = self._make_adata(
            obs_cols={"transcript_counts": [1, 2, 3], "cell_area": [10.0, 20.0, 30.0]},
            n=3,
        )
        with pytest.raises(ValueError, match="could not resolve spatial coords"):
            _resolve_spatial_coords(adata, side_name="xenium.h5ad")

    def test_obsm_spatial_beats_obs_pairs(self):
        # Precedence: if BOTH .obsm['spatial'] and an obs pair are present,
        # .obsm wins (matches the fallback-chain order).
        adata = self._make_adata(
            obs_cols={"centroid_x": [99.0], "centroid_y": [99.0]},
            obsm={"spatial": [[1.0, 2.0]]},
            n=1,
        )
        xy = _resolve_spatial_coords(adata)
        assert xy.tolist() == [[1.0, 2.0]]

    def test_obs_pair_precedence_order(self):
        # centroid_x/y beats x_centroid/y_centroid beats x/y when all present.
        adata = self._make_adata(
            obs_cols={
                "centroid_x": [1.0], "centroid_y": [2.0],
                "x_centroid": [10.0], "y_centroid": [20.0],
                "x": [100.0], "y": [200.0],
            },
            n=1,
        )
        xy = _resolve_spatial_coords(adata)
        assert xy.tolist() == [[1.0, 2.0]]
        # Confirm _SPATIAL_COORD_OBS_PAIRS enshrines the same order.
        assert _SPATIAL_COORD_OBS_PAIRS[0] == ("centroid_x", "centroid_y")

    def test_explicit_override_wins(self):
        # Explicit x_col + y_col hints skip the fallback chain entirely.
        adata = self._make_adata(
            obs_cols={
                "custom_x": [7.7], "custom_y": [8.8],
                "centroid_x": [1.0], "centroid_y": [2.0],   # decoy
            },
            obsm={"spatial": [[99.0, 99.0]]},                # decoy
            n=1,
        )
        xy = _resolve_spatial_coords(adata, "custom_x", "custom_y")
        assert xy.tolist() == [[7.7, 8.8]]

    def test_explicit_override_missing_col_fails_loud(self):
        adata = self._make_adata(
            obs_cols={"centroid_x": [1.0], "centroid_y": [2.0]},
            n=1,
        )
        with pytest.raises(KeyError, match="missing centroid column"):
            _resolve_spatial_coords(adata, "nope_x", "nope_y")

    def test_partial_obs_pair_not_matched(self):
        # Only ONE side of a pair present → keep walking the chain.
        adata = self._make_adata(
            obs_cols={"centroid_x": [1.0], "y": [2.0], "x": [3.0]},
            n=1,
        )
        # centroid_x present but centroid_y absent → skip.
        # x_centroid/y_centroid both absent → skip.
        # x/y both present → use them.
        xy = _resolve_spatial_coords(adata)
        assert xy.tolist() == [[3.0, 2.0]]

    def test_obsm_spatial_wrong_shape_raises(self):
        adata = self._make_adata(
            obs_cols={"celltype": ["A"]},
            obsm={"spatial": [[1.0]]},   # only one column
            n=1,
        )
        with pytest.raises(ValueError, match="expected \\(n, >=2\\)"):
            _resolve_spatial_coords(adata)


# ---------------------------------------------------------------------
# _celltype_from_nn — end-to-end NN mapping.
#
# Verifies the CORRECTED direction (Tracy's #15 amendment):
#   fit NN on proseg_purified centroids;
#   query with xenium centroids;
#   each xenium cell inherits its nearest proseg cell's celltype.
# ---------------------------------------------------------------------
def _write_proseg(tmp_path: Path, centroids, labels, celltype_col="first_type"):
    obs = pd.DataFrame({
        "x": [c[0] for c in centroids],
        "y": [c[1] for c in centroids],
        celltype_col: labels,
    })
    obs.index = [str(i) for i in range(len(centroids))]  # Proseg int64 index cast to str
    adata = ad.AnnData(
        X=np.zeros((len(centroids), 2), dtype=np.float32),
        obs=obs,
    )
    path = tmp_path / "proseg_purified.h5ad"
    adata.write_h5ad(path)
    return path


def _write_xenium(tmp_path: Path, centroids, uuids):
    obs = pd.DataFrame({
        "cell_id": uuids,
        "x_centroid": [c[0] for c in centroids],
        "y_centroid": [c[1] for c in centroids],
    })
    obs.index = uuids
    adata = ad.AnnData(
        X=np.zeros((len(centroids), 2), dtype=np.float32),
        obs=obs,
    )
    path = tmp_path / "xenium_ranger.h5ad"
    adata.write_h5ad(path)
    return path


class TestCelltypeFromNn:
    def test_nn_direction_proseg_to_xenium(self, tmp_path):
        # Two proseg cells: Tumor at (0, 0), Fibroblast at (100, 100).
        # Three xenium cells: (5, 5) → near Tumor; (95, 95) → near Fibro;
        # (50, 50) → equidistant → tie broken by NN sklearn ordering.
        proseg = _write_proseg(
            tmp_path,
            centroids=[(0.0, 0.0), (100.0, 100.0)],
            labels=["Tumor", "Fibroblast"],
        )
        xenium = _write_xenium(
            tmp_path,
            centroids=[(5.0, 5.0), (95.0, 95.0), (50.0, 50.0)],
            uuids=["aaaaafep-1", "lkkgeihi-1", "knjdobdk-1"],
        )
        annot = _celltype_from_nn(
            proseg_purified_h5ad=proseg,
            xenium_h5ad=xenium,
            celltype_col="auto",
            id_col="auto",
            k=1,
        )
        assert list(annot.columns) == ["xenium_cell_id", "group", "nn_distance"]
        assert len(annot) == 3
        # First xenium cell → nearest proseg is Tumor.
        assert annot.loc[annot["xenium_cell_id"] == "aaaaafep-1", "group"].iloc[0] == "Tumor"
        # Second xenium cell → nearest is Fibroblast.
        assert annot.loc[annot["xenium_cell_id"] == "lkkgeihi-1", "group"].iloc[0] == "Fibroblast"

    def test_every_xenium_cell_labelled_no_unassigned(self, tmp_path):
        # Design constraint from Tracy: NO mark_unassigned policy.
        # Every xenium cell — even ones very far from proseg cells —
        # gets a label.
        proseg = _write_proseg(
            tmp_path,
            centroids=[(0.0, 0.0)],
            labels=["Tumor"],
        )
        xenium = _write_xenium(
            tmp_path,
            centroids=[(999_999.0, 999_999.0)],   # extremely far
            uuids=["aaaaafep-1"],
        )
        annot = _celltype_from_nn(
            proseg_purified_h5ad=proseg,
            xenium_h5ad=xenium,
            celltype_col="auto",
            id_col="auto",
        )
        assert len(annot) == 1
        # Even the very far cell gets Tumor (nearest available label).
        # No "Unassigned" or "Unclassified" — those come later at the
        # Unclassified-fillna stage on cells missing from the annot.
        assert annot["group"].iloc[0] == "Tumor"

    def test_id_col_auto_detects_uuid_shaped(self, tmp_path):
        # Cross-check that the shape-check integration works end-to-end.
        # If auto-detect picked Proseg's int64-derived cell_id, this
        # would return "0" strings that won't join to any warped parquet.
        proseg = _write_proseg(
            tmp_path,
            centroids=[(10.0, 10.0), (20.0, 20.0)],
            labels=["A", "B"],
        )
        xenium = _write_xenium(
            tmp_path,
            centroids=[(11.0, 11.0), (21.0, 21.0)],
            uuids=["aaaaafep-1", "lkkgeihi-1"],
        )
        annot = _celltype_from_nn(
            proseg_purified_h5ad=proseg,
            xenium_h5ad=xenium,
            celltype_col="auto",
            id_col="auto",
        )
        assert set(annot["xenium_cell_id"]) == {"aaaaafep-1", "lkkgeihi-1"}
        # Every returned id looks like a xenium UUID (regex).
        import re
        for i in annot["xenium_cell_id"]:
            assert re.match(r"^[a-z]{8}-\d+$", i), f"non-UUID id leaked into output: {i!r}"

    def test_auto_detects_celltype_col_first_type(self, tmp_path):
        # Today's proseg_purified.h5ad carries `first_type`, not `celltype`.
        # Auto-detect must fall through the precedence chain.
        proseg = _write_proseg(
            tmp_path,
            centroids=[(0.0, 0.0)],
            labels=["Tumor"],
            celltype_col="first_type",   # no `celltype` column
        )
        xenium = _write_xenium(
            tmp_path,
            centroids=[(1.0, 1.0)],
            uuids=["aaaaafep-1"],
        )
        annot = _celltype_from_nn(
            proseg_purified_h5ad=proseg,
            xenium_h5ad=xenium,
            celltype_col="auto",
            id_col="auto",
        )
        assert annot["group"].iloc[0] == "Tumor"

    def test_xenium_coords_from_obsm_spatial(self, tmp_path):
        # Regression for Tracy's slurm-64031325 failure: xenium.h5ad has
        # spatial coords in .obsm['spatial'] and NONE of the obs pairs.
        # Old code hard-required .obs['x_centroid'] and crashed.
        proseg_path = tmp_path / "proseg_purified.h5ad"
        proseg_obs = pd.DataFrame(
            {"centroid_x": [0.0, 100.0], "centroid_y": [0.0, 100.0],
             "first_type": ["Tumor", "Fibroblast"]},
            index=["0", "1"],
        )
        ad.AnnData(
            X=np.zeros((2, 2), dtype=np.float32), obs=proseg_obs
        ).write_h5ad(proseg_path)

        xenium_path = tmp_path / "xenium.h5ad"
        xenium_obs = pd.DataFrame(
            # None of the historical obs pairs — coords must come from obsm.
            {"celltype": ["a", "b", "c"]},
            index=["aaaaafep-1", "lkkgeihi-1", "knjdobdk-1"],
        )
        xenium_adata = ad.AnnData(
            X=np.zeros((3, 2), dtype=np.float32), obs=xenium_obs,
        )
        xenium_adata.obsm["spatial"] = np.array(
            [[5.0, 5.0], [95.0, 95.0], [50.0, 50.0]], dtype=float,
        )
        xenium_adata.write_h5ad(xenium_path)

        annot = _celltype_from_nn(
            proseg_purified_h5ad=proseg_path,
            xenium_h5ad=xenium_path,
            celltype_col="auto",
            id_col="auto",
        )
        assert len(annot) == 3
        assert annot.loc[annot["xenium_cell_id"] == "aaaaafep-1", "group"].iloc[0] == "Tumor"
        assert annot.loc[annot["xenium_cell_id"] == "lkkgeihi-1", "group"].iloc[0] == "Fibroblast"


def _make_warped_parquet(path: Path, xenium_ids, polygons):
    """Emit a warp-stage-shaped parquet — mirror of the helper in
    ``test_pipeline_e2e.py``. Kept local so this module has no
    cross-file test-fixture import.
    """
    from shapely import to_wkb
    df = pd.DataFrame({
        "geometry": [to_wkb(p) for p in polygons],
    })
    df.index = xenium_ids
    df.index.name = "__null_dask_index__"
    df.to_parquet(path)


# ---------------------------------------------------------------------
# Ranger-direct celltype source (Tracy `5352008910` Ask 1).
#
# These tests cover the NEW default path introduced on commit
# ``26e78b0`` and the routing fix from skeptic F1 (commit that lands
# on top): with just ``--xenium-h5ad``, hexenium reads
# ``.obs[<celltype_col>]`` directly instead of NN-mapping from
# proseg_purified. Adds direct test coverage that skeptic F2 flagged
# as absent from ``test_celltyping.py``.
# ---------------------------------------------------------------------
class TestRangerDirect:
    def _ranger_h5ad(
        self,
        tmp_path: Path,
        celltype_values,
        celltype_col: str = "celltype",
        uuids=("aaaaafep-1", "lkkgeihi-1", "knjdobdk-1"),
    ) -> Path:
        """Build a xenium_ranger-shaped h5ad with a celltype column."""
        n = len(uuids)
        obs = pd.DataFrame(
            {celltype_col: list(celltype_values)},
            index=list(uuids),
        )
        adata = ad.AnnData(
            X=np.zeros((n, 2), dtype=np.float32), obs=obs,
        )
        adata.obsm["spatial"] = np.array(
            [[10.0 * i, 10.0 * i] for i in range(n)], dtype=float,
        )
        path = tmp_path / "xenium_ranger.h5ad"
        adata.write_h5ad(path)
        return path

    def test_ranger_direct_reads_obs_celltype_column(self, tmp_path: Path):
        # Default source path: three cells labelled by their own
        # .obs["celltype"] value. No NN — labels come straight off
        # the query h5ad.
        h5ad_path = self._ranger_h5ad(
            tmp_path, celltype_values=["Tumor", "Fibroblast", "Bcell"],
        )
        annot = _celltype_from_xenium_ranger_direct(
            xenium_h5ad=h5ad_path,
            celltype_col="celltype",
            id_col="auto",
        )
        assert len(annot) == 3
        assert set(annot["group"]) == {"Tumor", "Fibroblast", "Bcell"}
        # nn_distance column is present + all-NaN — no NN happened.
        assert "nn_distance" in annot.columns
        assert annot["nn_distance"].isna().all()
        # ids preserved (whitespace-stripped by the helper).
        assert set(annot["xenium_cell_id"]) == {
            "aaaaafep-1", "lkkgeihi-1", "knjdobdk-1",
        }

    def test_ranger_direct_all_nan_column_falls_to_unlabeled(
        self, tmp_path: Path,
    ):
        # Ask 2c: an existing but entirely-NaN celltype column also
        # falls back to UNLABELED for every row (no crash).
        h5ad_path = self._ranger_h5ad(
            tmp_path, celltype_values=[np.nan, np.nan, np.nan],
        )
        annot = _celltype_from_xenium_ranger_direct(
            xenium_h5ad=h5ad_path,
            celltype_col="celltype",
            id_col="auto",
        )
        assert len(annot) == 3
        assert (annot["group"] == UNLABELED).all()

    def test_ranger_direct_mixed_nan_falls_to_unlabeled_per_row(
        self, tmp_path: Path,
    ):
        # Companion to Ask 2c: per-row NaN falls to UNLABELED while
        # non-null rows keep their labels.
        h5ad_path = self._ranger_h5ad(
            tmp_path, celltype_values=["Tumor", np.nan, "Bcell"],
        )
        annot = _celltype_from_xenium_ranger_direct(
            xenium_h5ad=h5ad_path,
            celltype_col="celltype",
            id_col="auto",
        )
        groups = dict(zip(annot["xenium_cell_id"], annot["group"]))
        assert groups["aaaaafep-1"] == "Tumor"
        assert groups["lkkgeihi-1"] == UNLABELED
        assert groups["knjdobdk-1"] == "Bcell"

    def test_ranger_direct_missing_column_falls_to_unlabeled(
        self, tmp_path: Path,
    ):
        # Ask 2 headline case: the requested column doesn't exist on
        # the h5ad → WARN + all rows UNLABELED. Was KeyError.
        h5ad_path = self._ranger_h5ad(
            tmp_path, celltype_values=["Tumor", "Fibroblast", "Bcell"],
            celltype_col="not_the_celltype_col",  # column named differently
        )
        annot = _celltype_from_xenium_ranger_direct(
            xenium_h5ad=h5ad_path,
            celltype_col="celltype",       # asking for a column that doesn't exist
            id_col="auto",
        )
        assert len(annot) == 3
        assert (annot["group"] == UNLABELED).all()

    def test_run_celltyping_routes_to_ranger_direct_when_only_xenium_h5ad(
        self, tmp_path: Path,
    ):
        # The routing test that would have caught skeptic F1.
        # Given `xenium_h5ad` but NOT `proseg_purified_h5ad`, run_celltyping
        # must go through _celltype_from_xenium_ranger_direct (no NN).
        # We monkeypatch both entry points to spy on which fires.
        from hexenium.stages import celltyping as _ct

        h5ad_path = self._ranger_h5ad(
            tmp_path, celltype_values=["Tumor", "Fibroblast", "Bcell"],
        )
        # Minimal warp parquet so run_celltyping doesn't refuse on
        # missing warp inputs. Reuse the shape from the test helper
        # at the top of this module (WKB geometry + __null_dask_index__).
        warp_dir = tmp_path / "warp"; warp_dir.mkdir()
        _make_warped_parquet(
            warp_dir / "he_cell_seg.parquet",
            xenium_ids=["aaaaafep-1", "lkkgeihi-1", "knjdobdk-1"],
            # Polygons must be area > `area_threshold_px` (default 20 px²)
            # or `_clean_boundary_gdf` drops them. 10x10 squares (area
            # 100) match real warp output scale + skeptic F-A's verified
            # E2E fixture.
            polygons=[
                Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]),
                Polygon([(20, 20), (30, 20), (30, 30), (20, 30)]),
                Polygon([(40, 40), (50, 40), (50, 50), (40, 50)]),
            ],
        )
        out_dir = tmp_path / "celltyped"; out_dir.mkdir()

        calls = {"nn": 0, "ranger": 0}
        real_ranger = _ct._celltype_from_xenium_ranger_direct
        real_nn = _ct._celltype_from_nn

        def spy_ranger(**kw):
            calls["ranger"] += 1
            return real_ranger(**kw)
        def spy_nn(**kw):
            calls["nn"] += 1
            return real_nn(**kw)

        import pytest as _pt
        monkey = _pt.MonkeyPatch()
        try:
            monkey.setattr(_ct, "_celltype_from_xenium_ranger_direct", spy_ranger)
            monkey.setattr(_ct, "_celltype_from_nn", spy_nn)
            run_celltyping(
                sample_id="SAMPLE1",
                warp_dir=warp_dir,
                out_dir=out_dir,
                xenium_h5ad=h5ad_path,
                proseg_purified_h5ad=None,   # THE routing input
                celltype_col="celltype",
                id_col="auto",
            )
        finally:
            monkey.undo()

        assert calls["ranger"] == 1, (
            "expected _celltype_from_xenium_ranger_direct to fire (new default)"
        )
        assert calls["nn"] == 0, (
            "legacy _celltype_from_nn should NOT fire when "
            "proseg_purified_h5ad is None (F1 routing fix)"
        )

    def test_run_celltyping_routes_to_legacy_nn_when_proseg_purified_set(
        self, tmp_path: Path,
    ):
        # Companion to the above: passing --proseg-purified-h5ad
        # explicitly opts back into the legacy proseg-NN path.
        from hexenium.stages import celltyping as _ct

        h5ad_path = self._ranger_h5ad(
            tmp_path, celltype_values=["Tumor", "Fibroblast", "Bcell"],
        )
        # Minimal proseg h5ad with celltype + centroids so the legacy
        # NN path can fit + run without erroring.
        proseg_path = tmp_path / "proseg_purified.h5ad"
        proseg_obs = pd.DataFrame(
            {"celltype": ["Tumor", "Fibroblast", "Bcell"],
             "centroid_x": [0.0, 10.0, 20.0],
             "centroid_y": [0.0, 10.0, 20.0]},
            index=["0", "1", "2"],
        )
        ad.AnnData(
            X=np.zeros((3, 2), dtype=np.float32), obs=proseg_obs,
        ).write_h5ad(proseg_path)

        warp_dir = tmp_path / "warp"; warp_dir.mkdir()
        _make_warped_parquet(
            warp_dir / "he_cell_seg.parquet",
            xenium_ids=["aaaaafep-1", "lkkgeihi-1", "knjdobdk-1"],
            # Polygons must be area > `area_threshold_px` (default 20 px²)
            # or `_clean_boundary_gdf` drops them. 10x10 squares (area
            # 100) match real warp output scale + skeptic F-A's verified
            # E2E fixture.
            polygons=[
                Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]),
                Polygon([(20, 20), (30, 20), (30, 30), (20, 30)]),
                Polygon([(40, 40), (50, 40), (50, 50), (40, 50)]),
            ],
        )
        out_dir = tmp_path / "celltyped"; out_dir.mkdir()

        calls = {"nn": 0, "ranger": 0}
        real_ranger = _ct._celltype_from_xenium_ranger_direct
        real_nn = _ct._celltype_from_nn

        def spy_ranger(**kw):
            calls["ranger"] += 1
            return real_ranger(**kw)
        def spy_nn(**kw):
            calls["nn"] += 1
            return real_nn(**kw)

        import pytest as _pt
        monkey = _pt.MonkeyPatch()
        try:
            monkey.setattr(_ct, "_celltype_from_xenium_ranger_direct", spy_ranger)
            monkey.setattr(_ct, "_celltype_from_nn", spy_nn)
            run_celltyping(
                sample_id="SAMPLE1",
                warp_dir=warp_dir,
                out_dir=out_dir,
                xenium_h5ad=h5ad_path,
                proseg_purified_h5ad=proseg_path,   # explicit legacy opt-in
                celltype_col="auto",
                id_col="auto",
            )
        finally:
            monkey.undo()

        assert calls["nn"] == 1, (
            "legacy _celltype_from_nn should fire when "
            "proseg_purified_h5ad is set explicitly"
        )
        assert calls["ranger"] == 0, (
            "ranger-direct path should NOT fire when the legacy proseg-NN "
            "was explicitly opted into via --proseg-purified-h5ad"
        )
