"""Unit tests for hexenium.set_default_run — the CLI subcommand that
retargets ``output/<stage>`` symlinks after visual QA.

Coverage matrix (per Tracy's spec in
``settylab/TracyY123-nexus#15``):

* Happy path — partial + full updates land correctly.
* Lineage validation — register↔warp mismatches refuse.
* ``--force-lineage`` — bypass path emits a warning but does the
  update.
* Idempotency — re-running with the same targets is a no-op.
* Preservation — non-default ``<he_job_id>`` folders stay on disk.
* Refusal — target dir missing → no symlinks touched.

Run:
    pytest tests/test_set_default_run.py -v
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from hexenium.manifest import output_symlink_run_id, write_manifest
from hexenium.set_default_run import set_default_run


def _seed_stage(root: Path, stage: str, run_id: str, *,
                extra: dict | None = None) -> Path:
    """Build ``<root>/<stage>/<run_id>/`` with a plausible manifest.

    Tests use this to lay down as many parallel runs as needed under
    one stage, then flip the ``output/`` pointer with set_default_run.
    ``extra`` is merged into the manifest payload so tests can set
    ``source_register_run_id`` on warp runs.
    """
    d = root / stage / run_id
    d.mkdir(parents=True, exist_ok=True)
    payload = {
        "sample_id": "S1",
        "he_job_id": run_id,
        "params_hash": f"hash-{stage}-{run_id}",
    }
    if extra:
        payload.update(extra)
    write_manifest(d / "manifest.yaml", payload)
    return d


@pytest.fixture
def he_root(tmp_path: Path) -> Path:
    """A minimal ``he_registration/`` tree with two registrations and
    one warp that consumed the first.

    Layout:
        <he_root>/
          register/reg_A/manifest.yaml
          register/reg_B/manifest.yaml
          warp/warp_A/manifest.yaml   (source_register_run_id: reg_A)
          celltyped/ct_A/manifest.yaml
          viz/viz_A/manifest.yaml
    """
    root = tmp_path / "he_registration"
    _seed_stage(root, "register", "reg_A")
    _seed_stage(root, "register", "reg_B")
    _seed_stage(root, "warp", "warp_A",
                extra={"source_register_run_id": "reg_A"})
    _seed_stage(root, "celltyped", "ct_A")
    _seed_stage(root, "viz", "viz_A")
    return root


# ---------------------------------------------------------------------
# Happy path — partial + full updates.
# ---------------------------------------------------------------------
class TestHappyPath:
    def test_creates_output_dir_on_first_invocation(self, he_root: Path):
        # A fresh he_registration tree has no output/ folder yet; the
        # first set-default-run should create it and populate it.
        assert not (he_root / "output").exists()
        set_default_run(
            output_root_he=he_root,
            register_run_id="reg_A",
            warp_run_id="warp_A",
            celltyped_run_id="ct_A",
            viz_run_id="viz_A",
        )
        assert (he_root / "output").is_dir()
        for stage, expected in [("register", "reg_A"), ("warp", "warp_A"),
                                ("celltyped", "ct_A"), ("viz", "viz_A")]:
            link = he_root / "output" / stage
            assert link.is_symlink()
            assert output_symlink_run_id(link) == expected

    def test_partial_update_leaves_others_intact(self, he_root: Path):
        set_default_run(
            output_root_he=he_root,
            register_run_id="reg_A",
            warp_run_id="warp_A",
        )
        set_default_run(
            output_root_he=he_root,
            celltyped_run_id="ct_A",
        )
        # Register/warp untouched; celltype now added; viz still absent.
        assert output_symlink_run_id(he_root / "output" / "register") == "reg_A"
        assert output_symlink_run_id(he_root / "output" / "warp") == "warp_A"
        assert output_symlink_run_id(he_root / "output" / "celltyped") == "ct_A"
        assert not (he_root / "output" / "viz").exists()

    def test_relative_symlink_target(self, he_root: Path):
        set_default_run(output_root_he=he_root, register_run_id="reg_A")
        raw = os.readlink(he_root / "output" / "register")
        assert raw == "../register/reg_A"  # portable when tree is moved.


# ---------------------------------------------------------------------
# Lineage validation — register ↔ warp only, per Tracy's spec (C).
# ---------------------------------------------------------------------
class TestLineageValidation:
    def test_refuses_register_mismatched_with_current_warp(self, he_root: Path):
        # Seed the current pointer at warp_A (which consumed reg_A).
        set_default_run(
            output_root_he=he_root,
            register_run_id="reg_A",
            warp_run_id="warp_A",
        )
        # Now try to point register at reg_B — mismatches warp_A's source.
        with pytest.raises(SystemExit) as exc:
            set_default_run(output_root_he=he_root, register_run_id="reg_B")
        assert "lineage mismatch" in str(exc.value)
        # Nothing changed.
        assert output_symlink_run_id(he_root / "output" / "register") == "reg_A"

    def test_refuses_warp_mismatched_with_current_register(self, he_root: Path):
        # Seed pointers.
        set_default_run(
            output_root_he=he_root,
            register_run_id="reg_A",
            warp_run_id="warp_A",
        )
        # Introduce a second warp whose source is reg_B, then try to
        # point warp at it while output/register is still reg_A.
        _seed_stage(he_root, "warp", "warp_B",
                    extra={"source_register_run_id": "reg_B"})
        with pytest.raises(SystemExit) as exc:
            set_default_run(output_root_he=he_root, warp_run_id="warp_B")
        assert "lineage mismatch" in str(exc.value)
        assert output_symlink_run_id(he_root / "output" / "warp") == "warp_A"

    def test_both_at_once_consistent_is_allowed(self, he_root: Path):
        # Passing register + warp together that DO agree — no refusal.
        _seed_stage(he_root, "warp", "warp_B",
                    extra={"source_register_run_id": "reg_B"})
        set_default_run(
            output_root_he=he_root,
            register_run_id="reg_B",
            warp_run_id="warp_B",
        )
        assert output_symlink_run_id(he_root / "output" / "register") == "reg_B"
        assert output_symlink_run_id(he_root / "output" / "warp") == "warp_B"

    def test_force_lineage_bypasses_refusal(self, he_root: Path):
        set_default_run(
            output_root_he=he_root,
            register_run_id="reg_A",
            warp_run_id="warp_A",
        )
        # Deliberately point register at reg_B despite warp_A recording
        # reg_A — with --force-lineage this must succeed.
        set_default_run(
            output_root_he=he_root,
            register_run_id="reg_B",
            force_lineage=True,
        )
        assert output_symlink_run_id(he_root / "output" / "register") == "reg_B"

    def test_no_lineage_field_skips_check(self, he_root: Path):
        # Older warp outputs may lack source_register_run_id — the
        # check should skip silently rather than refuse.
        _seed_stage(he_root, "warp", "warp_legacy", extra={})
        set_default_run(
            output_root_he=he_root,
            register_run_id="reg_A",
            warp_run_id="warp_legacy",
        )
        assert output_symlink_run_id(he_root / "output" / "warp") == "warp_legacy"


# ---------------------------------------------------------------------
# Idempotency + preservation.
# ---------------------------------------------------------------------
class TestIdempotencyAndPreservation:
    def test_repeat_invocation_is_noop(self, he_root: Path):
        set_default_run(
            output_root_he=he_root,
            register_run_id="reg_A",
            warp_run_id="warp_A",
        )
        result = set_default_run(
            output_root_he=he_root,
            register_run_id="reg_A",
            warp_run_id="warp_A",
        )
        for change in result.changes:
            assert "unchanged" in change

    def test_non_default_run_dirs_preserved(self, he_root: Path):
        # Before: both reg_A and reg_B exist on disk. After we point at
        # reg_A, reg_B must still be there — this is the D contract.
        set_default_run(output_root_he=he_root, register_run_id="reg_A")
        assert (he_root / "register" / "reg_A").is_dir()
        assert (he_root / "register" / "reg_B").is_dir()
        assert (he_root / "register" / "reg_B" / "manifest.yaml").exists()

    def test_re_pointing_preserves_previous_target(self, he_root: Path):
        set_default_run(output_root_he=he_root, register_run_id="reg_A")
        set_default_run(output_root_he=he_root, register_run_id="reg_B")
        assert (he_root / "register" / "reg_A").is_dir()
        assert (he_root / "register" / "reg_B").is_dir()
        assert output_symlink_run_id(he_root / "output" / "register") == "reg_B"


# ---------------------------------------------------------------------
# Refusal on missing target — atomicity of the "no half-updates" rule.
# ---------------------------------------------------------------------
class TestMissingTargets:
    def test_refuses_when_any_target_missing(self, he_root: Path):
        with pytest.raises(SystemExit) as exc:
            set_default_run(
                output_root_he=he_root,
                register_run_id="reg_A",
                warp_run_id="warp_DOES_NOT_EXIST",
            )
        assert "target directories not found" in str(exc.value)
        # Neither symlink was touched.
        assert not (he_root / "output" / "register").exists()
        assert not (he_root / "output" / "warp").exists()

    def test_refuses_when_no_flags_given(self, he_root: Path):
        with pytest.raises(SystemExit) as exc:
            set_default_run(output_root_he=he_root)
        assert "at least one" in str(exc.value)
