"""Unit tests for hexenium.layout — the RunLayout module that computes
all output paths + the he_job_id fallback chain.

Run:
    pytest tests/test_layout.py -v
"""
from __future__ import annotations

from pathlib import Path

import pytest

from hexenium.layout import (
    RunLayout,
    STAGE_NAMES,
    compute_stages_suffix,
    is_xenium_uuid,
    resolve_he_job_id,
    resolve_layout,
)


# ---------------------------------------------------------------------
# is_xenium_uuid — shape-check guard used by celltyping.py's id-col
# auto-detect to reject Proseg's int64 "cell_id" masquerading as a UUID.
# ---------------------------------------------------------------------
class TestIsXeniumUuid:
    def test_matches_real_uuid(self):
        for s in ("aaaaafep-1", "lkkgeihi-1", "knjdobdk-1", "degchhml-1"):
            assert is_xenium_uuid(s), f"expected {s!r} to match"

    def test_rejects_proseg_int64_index(self):
        # This is the SILENT-FAILURE mode the design-pass skeptic flagged:
        # <S>_step4_purified.h5ad's `cell_id` column is `[0, 1, 2, 3, 5, ...]`
        # cast to int64. If we don't shape-check, the auto-detect picks
        # it and every join misses.
        for v in (0, 1, "0", "42", "1234"):
            assert not is_xenium_uuid(v), f"expected {v!r} to be rejected"

    def test_rejects_uppercase(self):
        # UUIDs are always lowercase in xenium output.
        assert not is_xenium_uuid("AAAAAFEP-1")

    def test_rejects_wrong_length(self):
        assert not is_xenium_uuid("aaaa-1")       # too short
        assert not is_xenium_uuid("aaaaaaaaa-1")  # too long

    def test_rejects_none_and_nan(self):
        import math
        assert not is_xenium_uuid(None)
        # NaN gets str()-ified to "nan" which shouldn't match.
        assert not is_xenium_uuid(math.nan)


# ---------------------------------------------------------------------
# compute_stages_suffix — folder-name suffix from --stages.
# Baked into logs_dir so re-runs at the
# same jobid don't clobber and operators can tell what a run covered.
# ---------------------------------------------------------------------
class TestComputeStagesSuffix:
    def test_all_stages_collapses_to_all(self):
        assert compute_stages_suffix(list(STAGE_NAMES)) == "all"

    def test_all_stages_reordered_still_all(self):
        # Order-insensitive: if every canonical stage is present, the
        # suffix collapses to "all" no matter what order they came in.
        reordered = list(reversed(STAGE_NAMES))
        assert compute_stages_suffix(reordered) == "all"

    def test_subset_preserves_cli_order(self):
        assert compute_stages_suffix(["celltype", "viz"]) == "celltype_viz"
        assert compute_stages_suffix(["viz", "celltype"]) == "viz_celltype"

    def test_single_stage(self):
        assert compute_stages_suffix(["register"]) == "register"

    def test_empty_input_is_none_sentinel(self):
        # Defensive — argparse default ensures this branch isn't hit in
        # normal use, but if someone builds a RunLayout by hand with an
        # empty stages list they get a clearly-labeled folder rather
        # than a silent duplicate of the pre-follow-up shape.
        assert compute_stages_suffix([]) == "none"
        assert compute_stages_suffix(None) == "none"


# ---------------------------------------------------------------------
# resolve_he_job_id — three-tier precedence:
#   --he-job-id > $SLURM_JOB_ID > YYYYMMDDTHHMMSS
# ---------------------------------------------------------------------
class TestResolveHeJobId:
    def test_explicit_wins_over_slurm(self, monkeypatch):
        monkeypatch.setenv("SLURM_JOB_ID", "12345")
        assert resolve_he_job_id("manual-run-1") == "manual-run-1"

    def test_explicit_wins_over_timestamp(self, monkeypatch):
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        assert resolve_he_job_id("manual-run-1") == "manual-run-1"

    def test_slurm_env_var_falls_through(self, monkeypatch):
        monkeypatch.setenv("SLURM_JOB_ID", "61438029")
        assert resolve_he_job_id() == "61438029"
        assert resolve_he_job_id(None) == "61438029"

    def test_empty_slurm_var_treated_as_unset(self, monkeypatch):
        # SLURM_JOB_ID="" (empty string) should NOT be honoured.
        # Otherwise back-to-back interactive runs would collide on "".
        monkeypatch.setenv("SLURM_JOB_ID", "")
        result = resolve_he_job_id()
        # Should have fallen through to timestamp — 15 chars YYYYMMDDTHHMMSS.
        assert len(result) == 15
        assert "T" in result

    def test_timestamp_fallback_shape(self, monkeypatch):
        monkeypatch.delenv("SLURM_JOB_ID", raising=False)
        result = resolve_he_job_id()
        assert len(result) == 15                    # YYYYMMDDTHHMMSS
        assert result[8] == "T"
        # First 8 chars are digits (date).
        assert result[:8].isdigit()
        assert result[9:].isdigit()

    def test_whitespace_stripped(self, monkeypatch):
        # Explicit override with padding.
        assert resolve_he_job_id("  61438029  ") == "61438029"


# ---------------------------------------------------------------------
# resolve_layout — standalone mode (sample_id + output_root).
# ---------------------------------------------------------------------
class TestResolveLayoutStandalone:
    def test_basic(self, tmp_path):
        layout = resolve_layout(
            sample_id="SAMPLE1",
            output_root=tmp_path,
            xenium_h5ad=None,
            he_job_id="job-999",
        )
        assert layout.sample_id == "SAMPLE1"
        assert layout.he_job_id == "job-999"
        assert layout.integrated is False
        assert layout.xenium_h5ad is None
        assert layout.output_root_he == tmp_path.resolve() / "SAMPLE1"

    def test_stage_dirs_include_he_job_id(self, tmp_path):
        layout = resolve_layout(
            sample_id="SAMPLE1",
            output_root=tmp_path,
            xenium_h5ad=None,
            he_job_id="61438029",
        )
        root = tmp_path.resolve() / "SAMPLE1"
        # Per-run sharded stages (flattened 2026-08-12 — no outer wrapper).
        assert layout.registration_dir == root / "register" / "61438029"
        assert layout.warp_dir() == root / "warp" / "61438029"
        assert layout.celltyped_dir() == root / "celltyped" / "61438029"
        assert layout.viz_dir == root / "viz" / "61438029"
        # Non-sharded top-level dirs.
        assert layout.converted_dir == root / "converted"
        assert layout.test_samples_dir == root / "test_samples"
        # Logs are also sharded per-run (logs shard per-run).
        # No stages passed → pre-follow-up shape (no suffix).
        assert layout.logs_dir == root / "logs" / "SAMPLE1_61438029"

    def test_standalone_logs_dir_with_stages_suffix(self, tmp_path):
        # --stages baked into logs
        # folder name so re-runs at the same jobid don't clobber and
        # operators can tell at a glance what a run covered.
        layout = resolve_layout(
            sample_id="SAMPLE1",
            output_root=tmp_path,
            xenium_h5ad=None,
            he_job_id="61438029",
            stages=["celltype", "viz"],
        )
        root = tmp_path.resolve() / "SAMPLE1"
        assert layout.stages_suffix == "celltype_viz"
        assert layout.logs_dir == root / "logs" / "SAMPLE1_61438029_celltype_viz"

    def test_standalone_logs_dir_all_stages_collapses_to_all(self, tmp_path):
        layout = resolve_layout(
            sample_id="SAMPLE1",
            output_root=tmp_path,
            xenium_h5ad=None,
            he_job_id="61438029",
            stages=list(STAGE_NAMES),
        )
        root = tmp_path.resolve() / "SAMPLE1"
        assert layout.stages_suffix == "all"
        assert layout.logs_dir == root / "logs" / "SAMPLE1_61438029_all"

    def test_run_id_overrides_for_viz_only(self, tmp_path):
        # --warp-run-id / --celltype-run-id let viz-only reruns
        # target a PRIOR he_job_id's outputs.
        layout = resolve_layout(
            sample_id="SAMPLE1",
            output_root=tmp_path,
            xenium_h5ad=None,
            he_job_id="new-job",
        )
        root = tmp_path.resolve() / "SAMPLE1"
        assert layout.warp_dir("old-warp") == root / "warp" / "old-warp"
        assert layout.celltyped_dir("old-ct") == root / "celltyped" / "old-ct"

    def test_missing_sample_id_and_no_h5ad(self, tmp_path):
        with pytest.raises(SystemExit, match="either --xenium-h5ad"):
            resolve_layout(
                sample_id=None,
                output_root=tmp_path,
                xenium_h5ad=None,
                he_job_id="j",
            )

    def test_ensure_dirs_creates_top_level(self, tmp_path):
        layout = resolve_layout(
            sample_id="SAMPLE1",
            output_root=tmp_path,
            xenium_h5ad=None,
            he_job_id="j",
        )
        layout.ensure_dirs()
        assert layout.output_root_he.exists()
        assert layout.logs_dir.exists()


# ---------------------------------------------------------------------
# resolve_layout — integrated mode: reads .uns from a real h5ad,
# derives run_dir from the h5ad's on-disk location.
# ---------------------------------------------------------------------
def _write_test_h5ad(
    path: Path,
    *,
    sample_id: str | None,
    run_id: str | None,
    n_cells: int = 5,
) -> Path:
    """Build a tiny h5ad with the .uns metadata contract."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    obs = pd.DataFrame({
        "cell_id": [f"aaaaafep-{i}" for i in range(n_cells)],
        "x_centroid": np.linspace(100.0, 200.0, n_cells),
        "y_centroid": np.linspace(100.0, 200.0, n_cells),
    })
    obs.index = obs["cell_id"]
    adata = ad.AnnData(
        X=np.zeros((n_cells, 3), dtype=np.float32),
        obs=obs,
    )
    if sample_id is not None:
        adata.uns["sample_id"] = sample_id
    if run_id is not None:
        adata.uns["run_id"] = run_id
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
    return path


class TestResolveLayoutIntegrated:
    def test_reads_uns_and_derives_run_dir(self, tmp_path):
        # Build the canonical upstream layout:
        # <root>/<sample>/<sample>_<run_id>/spatial_adata/<sample>_xenium_ranger.h5ad
        run_dir = tmp_path / "SAMPLE1" / "SAMPLE1_run42"
        h5ad = run_dir / "spatial_adata" / "SAMPLE1_xenium_ranger.h5ad"
        _write_test_h5ad(h5ad, sample_id="SAMPLE1", run_id="run42")

        layout = resolve_layout(
            sample_id=None,
            output_root=None,
            xenium_h5ad=h5ad,
            he_job_id="he-1",
        )
        assert layout.sample_id == "SAMPLE1"
        assert layout.he_job_id == "he-1"
        assert layout.integrated is True
        assert layout.xenium_run_id == "run42"
        assert layout.xenium_run_dir == run_dir.resolve()
        assert layout.output_root_he == run_dir.resolve() / "he_registration"

    def test_integrated_layout_paths_all_under_he_registration(self, tmp_path):
        run_dir = tmp_path / "SAMPLE1" / "SAMPLE1_run42"
        h5ad = run_dir / "spatial_adata" / "SAMPLE1_xenium_ranger.h5ad"
        _write_test_h5ad(h5ad, sample_id="SAMPLE1", run_id="run42")

        layout = resolve_layout(
            sample_id=None,
            output_root=None,
            xenium_h5ad=h5ad,
            he_job_id="he-1",
        )
        base = run_dir.resolve() / "he_registration"
        # Logs live under the RUN-LEVEL logs/ folder (peer of step1/3/4),
        # NOT under he_registration/ ().
        # No stages passed → pre-follow-up shape (no suffix).
        assert layout.logs_dir == (run_dir.resolve() / "logs"
                                   / "logs_heRegistration" / "SAMPLE1_he-1")
        assert layout.converted_dir == base / "converted"
        assert layout.registration_dir == base / "register" / "he-1"
        assert layout.warp_dir() == base / "warp" / "he-1"
        assert layout.celltyped_dir() == base / "celltyped" / "he-1"
        assert layout.viz_dir == base / "viz" / "he-1"

    def test_integrated_logs_dir_with_stages_suffix(self, tmp_path):
        run_dir = tmp_path / "SAMPLE1" / "SAMPLE1_run42"
        h5ad = run_dir / "spatial_adata" / "SAMPLE1_xenium_ranger.h5ad"
        _write_test_h5ad(h5ad, sample_id="SAMPLE1", run_id="run42")
        layout = resolve_layout(
            sample_id=None,
            output_root=None,
            xenium_h5ad=h5ad,
            he_job_id="he-1",
            stages=["register", "warp"],
        )
        assert layout.stages_suffix == "register_warp"
        assert layout.logs_dir == (run_dir.resolve() / "logs"
                                   / "logs_heRegistration"
                                   / "SAMPLE1_he-1_register_warp")

    def test_missing_uns_sample_id_fails_loud(self, tmp_path):
        run_dir = tmp_path / "SAMPLE1" / "SAMPLE1_run42"
        h5ad = run_dir / "spatial_adata" / "SAMPLE1_xenium_ranger.h5ad"
        _write_test_h5ad(h5ad, sample_id=None, run_id="run42")
        with pytest.raises(SystemExit, match="missing .uns"):
            resolve_layout(
                sample_id=None,
                output_root=None,
                xenium_h5ad=h5ad,
                he_job_id="he-1",
            )

    def test_missing_uns_run_id_fails_loud(self, tmp_path):
        run_dir = tmp_path / "SAMPLE1" / "SAMPLE1_run42"
        h5ad = run_dir / "spatial_adata" / "SAMPLE1_xenium_ranger.h5ad"
        _write_test_h5ad(h5ad, sample_id="SAMPLE1", run_id=None)
        with pytest.raises(SystemExit, match="missing .uns"):
            resolve_layout(
                sample_id=None,
                output_root=None,
                xenium_h5ad=h5ad,
                he_job_id="he-1",
            )

    def test_h5ad_outside_spatial_adata_fails_loud(self, tmp_path):
        # Wrong parent name — must be exactly "spatial_adata".
        run_dir = tmp_path / "SAMPLE1" / "SAMPLE1_run42"
        h5ad = run_dir / "wrong_parent" / "SAMPLE1_xenium_ranger.h5ad"
        _write_test_h5ad(h5ad, sample_id="SAMPLE1", run_id="run42")
        with pytest.raises(SystemExit, match="spatial_adata"):
            resolve_layout(
                sample_id=None,
                output_root=None,
                xenium_h5ad=h5ad,
                he_job_id="he-1",
            )

    def test_run_dir_name_mismatch_fails_loud(self, tmp_path):
        # .uns says run42 but folder is named differently → fail.
        run_dir = tmp_path / "SAMPLE1" / "renamed_by_hand"
        h5ad = run_dir / "spatial_adata" / "SAMPLE1_xenium_ranger.h5ad"
        _write_test_h5ad(h5ad, sample_id="SAMPLE1", run_id="run42")
        with pytest.raises(SystemExit, match="run dir name mismatch"):
            resolve_layout(
                sample_id=None,
                output_root=None,
                xenium_h5ad=h5ad,
                he_job_id="he-1",
            )

    def test_explicit_sample_id_must_agree_with_uns(self, tmp_path):
        run_dir = tmp_path / "SAMPLE1" / "SAMPLE1_run42"
        h5ad = run_dir / "spatial_adata" / "SAMPLE1_xenium_ranger.h5ad"
        _write_test_h5ad(h5ad, sample_id="SAMPLE1", run_id="run42")
        with pytest.raises(SystemExit, match="disagrees with"):
            resolve_layout(
                sample_id="SAMPLE99",   # disagrees with .uns
                output_root=None,
                xenium_h5ad=h5ad,
                he_job_id="he-1",
            )

    def test_explicit_sample_id_matching_uns_is_accepted(self, tmp_path):
        run_dir = tmp_path / "SAMPLE1" / "SAMPLE1_run42"
        h5ad = run_dir / "spatial_adata" / "SAMPLE1_xenium_ranger.h5ad"
        _write_test_h5ad(h5ad, sample_id="SAMPLE1", run_id="run42")
        layout = resolve_layout(
            sample_id="SAMPLE1",   # agrees
            output_root=None,
            xenium_h5ad=h5ad,
            he_job_id="he-1",
        )
        assert layout.sample_id == "SAMPLE1"
        assert layout.integrated is True


# ---------------------------------------------------------------------
# resolve_layout — integrated-by-run-id mode: CLI args carry identity,
# no h5ad read. Lets register+warp run before xenium-preprocess step-1
# has produced the h5ad.
# ---------------------------------------------------------------------
class TestResolveLayoutIntegratedByRunId:
    def test_derives_run_dir_from_cli_args(self, tmp_path):
        layout = resolve_layout(
            sample_id="SAMPLE2",
            output_root=tmp_path,
            xenium_h5ad=None,
            he_job_id="he-1",
            run_id="demo_v1",
        )
        run_dir = tmp_path.resolve() / "SAMPLE2" / "SAMPLE2_demo_v1"
        assert layout.integrated is True
        assert layout.sample_id == "SAMPLE2"
        assert layout.xenium_run_id == "demo_v1"
        assert layout.xenium_run_dir == run_dir
        assert layout.xenium_h5ad is None
        assert layout.output_root_he == run_dir / "he_registration"

    def test_layout_paths_all_under_he_registration(self, tmp_path):
        layout = resolve_layout(
            sample_id="SAMPLE2",
            output_root=tmp_path,
            xenium_h5ad=None,
            he_job_id="he-1",
            run_id="demo_v1",
        )
        run_dir = tmp_path.resolve() / "SAMPLE2" / "SAMPLE2_demo_v1"
        base = run_dir / "he_registration"
        # Logs live under the RUN-LEVEL logs/ folder (peer of step1/3/4),
        # NOT under he_registration/ ().
        # No stages passed → pre-follow-up shape (no suffix).
        assert layout.logs_dir == (run_dir / "logs" / "logs_heRegistration"
                                   / "SAMPLE2_he-1")
        assert layout.registration_dir == base / "register" / "he-1"
        assert layout.warp_dir() == base / "warp" / "he-1"
        assert layout.celltyped_dir() == base / "celltyped" / "he-1"
        assert layout.viz_dir == base / "viz" / "he-1"

    def test_integrated_by_run_id_logs_dir_with_stages_suffix(self, tmp_path):
        layout = resolve_layout(
            sample_id="SAMPLE2",
            output_root=tmp_path,
            xenium_h5ad=None,
            he_job_id="64027728",
            run_id="demo_v1",
            stages=["celltype", "viz"],
        )
        run_dir = tmp_path.resolve() / "SAMPLE2" / "SAMPLE2_demo_v1"
        # This is the exact shape matches the target layout
        # SAMPLE2_64027728_celltype_viz.
        assert layout.logs_dir == (run_dir / "logs" / "logs_heRegistration"
                                   / "SAMPLE2_64027728_celltype_viz")

    def test_no_uns_read_h5ad_absent_is_fine(self, tmp_path):
        # Explicit invariant: the branch must NOT touch a filesystem path
        # for the h5ad. tmp_path is empty; layout resolves anyway.
        layout = resolve_layout(
            sample_id="SAMPLE2",
            output_root=tmp_path,
            xenium_h5ad=None,
            he_job_id="he-1",
            run_id="demo_v1",
        )
        assert not (tmp_path / "SAMPLE2" / "SAMPLE2_demo_v1"
                    / "spatial_adata").exists()
        assert layout.output_root_he.name == "he_registration"

    def test_h5ad_wins_over_run_id_when_both_present(self, tmp_path):
        # If a real h5ad is passed, we take the .uns route even if
        # run_id was also given — the .uns is the source of truth.
        run_dir = tmp_path / "SAMPLE1" / "SAMPLE1_run42"
        h5ad = run_dir / "spatial_adata" / "SAMPLE1_xenium_ranger.h5ad"
        _write_test_h5ad(h5ad, sample_id="SAMPLE1", run_id="run42")
        layout = resolve_layout(
            sample_id="SAMPLE1",
            output_root=tmp_path,
            xenium_h5ad=h5ad,
            he_job_id="he-1",
            run_id="demo_v1",   # would-be conflict, ignored
        )
        assert layout.xenium_run_id == "run42"
        assert layout.xenium_run_dir == run_dir.resolve()
