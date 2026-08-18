"""End-to-end smoke test for the celltype+viz stages against a small
synthetic warped-parquet fixture. Verifies the reorg layout writes
files under the expected paths + that the celltype refactor's NN
mapping stitches into the existing cumcount-merge flow correctly.

We stub out the registrar+warp stages (they require VALIS + a real
Xenium bundle + tens of minutes). But the celltype+viz surface is
fixture-friendly since they consume plain parquets.

Run:
    /home/ryang/micromamba/envs/heRegistration/bin/python -m pytest \
        test_pipeline_e2e.py -v
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Polygon
from shapely import to_wkb

from hexenium.layout import resolve_layout, resolve_he_job_id


def _make_warped_parquet(path: Path, xenium_ids, polygons):
    """Emit a warp-stage-shaped parquet: WKB geometry, __null_dask_index__ = xenium_cell_id.

    Matches the dask.dataframe.to_parquet output that `_read_warped_gdf`
    reverses in celltyping.py.
    """
    df = pd.DataFrame({
        "geometry": [to_wkb(p) for p in polygons],
    })
    df.index = xenium_ids
    df.index.name = "__null_dask_index__"
    df.to_parquet(path)


def _make_proseg(path: Path, centroids, labels):
    obs = pd.DataFrame({
        "x": [c[0] for c in centroids],
        "y": [c[1] for c in centroids],
        "first_type": labels,
    })
    obs.index = [str(i) for i in range(len(centroids))]
    adata = ad.AnnData(X=np.zeros((len(centroids), 2), dtype=np.float32), obs=obs)
    adata.write_h5ad(path)
    return path


def _make_xenium(path: Path, centroids, uuids):
    obs = pd.DataFrame({
        "cell_id": uuids,
        "x_centroid": [c[0] for c in centroids],
        "y_centroid": [c[1] for c in centroids],
    })
    obs.index = uuids
    adata = ad.AnnData(X=np.zeros((len(centroids), 2), dtype=np.float32), obs=obs)
    adata.write_h5ad(path)
    return path


def _make_integrated_xenium(path: Path, centroids, uuids, sample_id, run_id):
    obs = pd.DataFrame({
        "cell_id": uuids,
        "x_centroid": [c[0] for c in centroids],
        "y_centroid": [c[1] for c in centroids],
    })
    obs.index = uuids
    adata = ad.AnnData(X=np.zeros((len(centroids), 2), dtype=np.float32), obs=obs)
    adata.uns["sample_id"] = sample_id
    adata.uns["run_id"] = run_id
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
    return path


class TestCelltypeStageStandalone:
    def test_writes_outputs_under_he_job_id_folder(self, tmp_path):
        # Two "cells" (polygons) with xenium UUIDs.
        uuids = ["aaaaafep-1", "lkkgeihi-1"]
        polys = [
            Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]),
            Polygon([(50, 50), (60, 50), (60, 60), (50, 60)]),
        ]

        # Set up the layout as if we ran through the reorg.
        layout = resolve_layout(
            sample_id="SAMPLE1",
            output_root=tmp_path,
            xenium_h5ad=None,
            he_job_id="myjob",
        )
        layout.ensure_dirs()

        warp_dir = layout.warp_dir()
        warp_dir.mkdir(parents=True, exist_ok=True)
        _make_warped_parquet(warp_dir / "he_cell_seg.parquet", uuids, polys)
        _make_warped_parquet(warp_dir / "he_nucleus_seg.parquet", uuids, polys)

        proseg = _make_proseg(
            tmp_path / "proseg.h5ad",
            centroids=[(5.0, 5.0), (55.0, 55.0)],
            labels=["Tumor", "Fibroblast"],
        )
        xenium = _make_xenium(
            tmp_path / "xenium.h5ad",
            centroids=[(5.0, 5.0), (55.0, 55.0)],
            uuids=uuids,
        )

        from hexenium.stages.celltyping import run_celltyping
        result = run_celltyping(
            sample_id="SAMPLE1",
            warp_dir=warp_dir,
            out_dir=layout.celltyped_dir(),
            xenium_h5ad=xenium,
            proseg_purified_h5ad=proseg,
        )

        # Every output landed under celltyped/<he_job_id>/.
        ct_dir = layout.celltyped_dir()
        assert ct_dir.exists()
        for expected in (
            f"SAMPLE1_cells_analysis.geojson",
            f"SAMPLE1_cells_qupath.geojson",
            f"SAMPLE1_nuclei_analysis.geojson",
            f"SAMPLE1_nuclei_qupath.geojson",
            f"SAMPLE1_celltyped_wholeslide.parquet",
            f"SAMPLE1_celltype_annotation.parquet",
        ):
            assert (ct_dir / expected).exists(), f"missing: {ct_dir / expected}"
        # And nothing landed under the OLD pre-flatten wrapper path
        # (dropped 2026-08-12 — was <root>/registration/celltyped/<id>/).
        assert not (tmp_path / "SAMPLE1" / "registration").exists()
        # Nor directly at <root>/celltyped/ without the he_job_id shard.
        assert not (tmp_path / "SAMPLE1" / "celltyped" / "SAMPLE1_cells_analysis.geojson").exists()

    def test_annotation_parquet_carries_nn_labels(self, tmp_path):
        uuids = ["aaaaafep-1", "lkkgeihi-1"]
        polys = [
            Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]),
            Polygon([(50, 50), (60, 50), (60, 60), (50, 60)]),
        ]
        layout = resolve_layout(
            sample_id="SAMPLE1",
            output_root=tmp_path,
            xenium_h5ad=None,
            he_job_id="myjob",
        )
        layout.ensure_dirs()
        warp_dir = layout.warp_dir()
        warp_dir.mkdir(parents=True, exist_ok=True)
        _make_warped_parquet(warp_dir / "he_cell_seg.parquet", uuids, polys)

        proseg = _make_proseg(
            tmp_path / "proseg.h5ad",
            centroids=[(5.0, 5.0), (55.0, 55.0)],
            labels=["Tumor", "Fibroblast"],
        )
        xenium = _make_xenium(
            tmp_path / "xenium.h5ad",
            centroids=[(5.0, 5.0), (55.0, 55.0)],
            uuids=uuids,
        )
        from hexenium.stages.celltyping import run_celltyping
        run_celltyping(
            sample_id="SAMPLE1",
            warp_dir=warp_dir,
            out_dir=layout.celltyped_dir(),
            xenium_h5ad=xenium,
            proseg_purified_h5ad=proseg,
        )
        annot = pd.read_parquet(layout.celltyped_dir() / "SAMPLE1_celltype_annotation.parquet")
        assert set(annot["xenium_cell_id"]) == set(uuids)
        # Cell at (5,5) → nearest proseg is (5,5)="Tumor"
        assert annot.loc[annot["xenium_cell_id"] == "aaaaafep-1", "group"].iloc[0] == "Tumor"


class TestIntegratedMode:
    def test_layout_colocates_under_xenium_run_dir(self, tmp_path):
        # Build the upstream canonical layout.
        sample_id, run_id = "SAMPLE1", "run42"
        xenium_run_dir = tmp_path / sample_id / f"{sample_id}_{run_id}"
        xenium_h5ad = _make_integrated_xenium(
            xenium_run_dir / "spatial_adata" / f"{sample_id}_xenium_ranger.h5ad",
            centroids=[(5.0, 5.0), (55.0, 55.0)],
            uuids=["aaaaafep-1", "lkkgeihi-1"],
            sample_id=sample_id, run_id=run_id,
        )
        # Prior xenium-preprocess resolved_config.yaml, to prove
        # sibling-preserving merge works.
        prior_cfg = xenium_run_dir / "resolved_config.yaml"
        prior_cfg.write_text("driver:\n  version: 1.2\nstep4:\n  n_iters: 5\n")

        layout = resolve_layout(
            sample_id=None,
            output_root=None,
            xenium_h5ad=xenium_h5ad,
            he_job_id="he-42",
        )
        assert layout.integrated is True
        assert layout.output_root_he == (xenium_run_dir / "he_registration").resolve()

        # Simulate the pipeline's merge step.
        from hexenium.pipeline import _merge_run_dir_config
        _merge_run_dir_config(layout, {
            "sample_id": sample_id,
            "he_job_id": "he-42",
        })
        # Merged file still carries the pre-existing keys.
        import yaml
        with open(prior_cfg) as f:
            merged = yaml.safe_load(f)
        assert "driver" in merged
        assert merged["driver"]["version"] == 1.2
        assert "step4" in merged
        assert merged["step4"]["n_iters"] == 5
        # And gained the he_registration: key.
        assert "he_registration" in merged
        assert merged["he_registration"]["sample_id"] == sample_id
        assert merged["he_registration"]["he_job_id"] == "he-42"

    def test_integrated_mode_writes_logs_under_run_level_logs(self, tmp_path):
        sample_id, run_id = "SAMPLE1", "run42"
        xenium_run_dir = tmp_path / sample_id / f"{sample_id}_{run_id}"
        xenium_h5ad = _make_integrated_xenium(
            xenium_run_dir / "spatial_adata" / f"{sample_id}_xenium_ranger.h5ad",
            centroids=[(5.0, 5.0)],
            uuids=["aaaaafep-1"],
            sample_id=sample_id, run_id=run_id,
        )
        layout = resolve_layout(
            sample_id=None,
            output_root=None,
            xenium_h5ad=xenium_h5ad,
            he_job_id="he-42",
        )
        layout.ensure_dirs()
        # Logs share the RUN-LEVEL logs/ dir with step1/3/4, in a
        # logs_heRegistration/<sid>_<jobid>/ subfolder
        # ().
        assert layout.logs_dir == (xenium_run_dir / "logs"
                                   / "logs_heRegistration"
                                   / "SAMPLE1_he-42").resolve()
        assert layout.logs_dir.exists()
        # A second concurrent he-run would land in a DIFFERENT subfolder.
        layout2 = resolve_layout(
            sample_id=None,
            output_root=None,
            xenium_h5ad=xenium_h5ad,
            he_job_id="he-99",
        )
        assert layout2.logs_dir != layout.logs_dir


class TestProsegPurifiedAutoDerive:
    """Auto-derive of --proseg-purified-h5ad from the xenium run dir
    (). Mirrors the xenium h5ad
    auto-derive from --run-id already in place."""

    def _run_id_layout(self, tmp_path, sample_id="SAMPLE1", run_id="demo_v1"):
        return resolve_layout(
            sample_id=sample_id,
            output_root=tmp_path,
            xenium_h5ad=None,
            he_job_id="he-1",
            run_id=run_id,
        )

    def test_derives_from_run_dir_when_file_exists(self, tmp_path):
        from hexenium.pipeline import _maybe_derive_proseg_purified
        layout = self._run_id_layout(tmp_path)
        target = (layout.xenium_run_dir / "spatial_adata"
                  / f"{layout.sample_id}_proseg_purified.h5ad")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("")  # existence-only; not read here
        cfg: dict = {}
        _maybe_derive_proseg_purified(["celltype"], cfg, layout)
        assert cfg["proseg_purified_h5ad"] == str(target)

    def test_fails_loud_when_derived_missing(self, tmp_path):
        from hexenium.pipeline import _maybe_derive_proseg_purified
        layout = self._run_id_layout(tmp_path)
        cfg: dict = {}
        with pytest.raises(SystemExit, match="proseg_purified h5ad not found"):
            _maybe_derive_proseg_purified(["celltype"], cfg, layout)

    def test_noop_when_celltype_not_in_stages(self, tmp_path):
        from hexenium.pipeline import _maybe_derive_proseg_purified
        layout = self._run_id_layout(tmp_path)
        cfg: dict = {}
        # No file, no celltype stage → do nothing, no raise.
        _maybe_derive_proseg_purified(["register", "warp"], cfg, layout)
        assert "proseg_purified_h5ad" not in cfg

    def test_explicit_override_wins(self, tmp_path):
        from hexenium.pipeline import _maybe_derive_proseg_purified
        layout = self._run_id_layout(tmp_path)
        cfg: dict = {"proseg_purified_h5ad": "/explicit/path.h5ad"}
        # No file at derived path, but explicit wins → no raise, no change.
        _maybe_derive_proseg_purified(["celltype"], cfg, layout)
        assert cfg["proseg_purified_h5ad"] == "/explicit/path.h5ad"

    def test_noop_in_standalone_mode_no_run_dir(self, tmp_path):
        from hexenium.pipeline import _maybe_derive_proseg_purified
        # Standalone: no --run-id → xenium_run_dir is None.
        layout = resolve_layout(
            sample_id="SAMPLE1",
            output_root=tmp_path,
            xenium_h5ad=None,
            he_job_id="he-1",
        )
        assert layout.xenium_run_dir is None
        cfg: dict = {}
        # No derive possible; must not raise (falls through to whatever
        # the celltype stage does when proseg is None — today that's
        # skip the NN mapping; not this function's concern).
        _maybe_derive_proseg_purified(["celltype"], cfg, layout)
        assert "proseg_purified_h5ad" not in cfg

    def test_integrated_by_h5ad_derives_via_run_dir(self, tmp_path):
        from hexenium.pipeline import _maybe_derive_proseg_purified
        sample_id, run_id = "SAMPLE1", "run42"
        run_dir = tmp_path / sample_id / f"{sample_id}_{run_id}"
        xenium_h5ad = _make_integrated_xenium(
            run_dir / "spatial_adata" / f"{sample_id}_xenium_ranger.h5ad",
            centroids=[(5.0, 5.0)],
            uuids=["aaaaafep-1"],
            sample_id=sample_id, run_id=run_id,
        )
        target = (run_dir / "spatial_adata"
                  / f"{sample_id}_proseg_purified.h5ad")
        target.write_text("")
        layout = resolve_layout(
            sample_id=None,
            output_root=None,
            xenium_h5ad=xenium_h5ad,
            he_job_id="he-42",
        )
        cfg: dict = {}
        _maybe_derive_proseg_purified(["celltype", "viz"], cfg, layout)
        assert cfg["proseg_purified_h5ad"] == str(target.resolve())
