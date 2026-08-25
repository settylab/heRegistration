"""Unit tests for the ``--set-default-on-success`` promotion hook.

Coverage matrix from Tracy's spec in
``settylab/TracyY123-nexus#15`` comment ``5336523381``:

* **Full 5-stage run** — every promotable stage (register, warp,
  celltype, viz) that ran gets its ``output/<stage>`` symlink pointed
  at the current invocation's ``<he_job_id>``.
* **Partial celltype+viz rerun with older register/warp lineage** —
  only the stages that ran are promoted; existing ``output/register``
  + ``output/warp`` are left alone. Lineage of the older register/warp
  must match the warp that celltype/viz consumed.
* **Lineage-validation failure** — celltype/viz consumed a warp
  different from the current ``output/warp``. Promote refuses;
  ``output/`` unchanged.
* **Stage failure** — never exercised in-pipeline (stage errors
  propagate out of ``run()`` before the hook fires), but the
  behavioural guarantee is asserted by construction: the hook is the
  LAST statement in ``run()``. Documented in
  ``test_hook_is_last_in_pipeline_run``.
* **Opt-out** — no flag = ``output/`` never touched, regardless of
  what ran.

The tests exercise the module-private helpers
``_resolve_promotion_ids`` + ``_validate_promotion_lineage`` +
``_promote_completed_run`` directly, plus a full-round-trip against
the on-disk tree via ``_promote_completed_run`` calling into the
real ``set_default_run``. VALIS/HEST are not touched.

Run:
    pytest tests/test_set_default_on_success.py -v
"""
from __future__ import annotations

from pathlib import Path

import pytest

from hexenium.layout import RunLayout
from hexenium.manifest import output_symlink_run_id, write_manifest
from hexenium.pipeline import (
    _promote_completed_run,
    _resolve_promotion_ids,
    _validate_promotion_lineage,
    _PROMOTE_STAGE_DIR,
)


def _seed_stage(root: Path, stage: str, run_id: str, *,
                extra: dict | None = None) -> Path:
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
    """Minimal ``he_registration/`` with two registrations and one warp
    that consumed the first (``source_register_run_id: reg_A``)."""
    root = tmp_path / "he_registration"
    _seed_stage(root, "register", "reg_A")
    _seed_stage(root, "register", "reg_B")
    _seed_stage(root, "warp", "warp_A",
                extra={"source_register_run_id": "reg_A"})
    _seed_stage(root, "celltyped", "ct_A")
    _seed_stage(root, "viz", "viz_A")
    return root


def _layout(he_root: Path, he_job_id: str) -> RunLayout:
    return RunLayout(
        sample_id="S1", he_job_id=he_job_id, output_root_he=he_root,
        integrated=False, stages_suffix="",
    )


# ---------------------------------------------------------------------
# _resolve_promotion_ids — only stages that ran + are promotable.
# ---------------------------------------------------------------------
class TestResolvePromotionIds:
    def test_maps_ran_stages_to_fresh_he_job_id(self, he_root: Path):
        layout = _layout(he_root, "job_new")
        ids = _resolve_promotion_ids(
            layout=layout, stages_ran=["register", "warp", "celltype", "viz"],
        )
        assert ids == {
            "register": "job_new", "warp": "job_new",
            "celltype": "job_new", "viz": "job_new",
        }

    def test_partial_stage_list(self, he_root: Path):
        layout = _layout(he_root, "job_new")
        ids = _resolve_promotion_ids(
            layout=layout, stages_ran=["celltype", "viz"],
        )
        assert ids == {"celltype": "job_new", "viz": "job_new"}

    def test_he_preprocess_is_not_promotable(self, he_root: Path):
        # he_preprocess writes to converted/ (un-sharded), so it must
        # never appear in the promote dict.
        layout = _layout(he_root, "job_new")
        ids = _resolve_promotion_ids(
            layout=layout, stages_ran=["he_preprocess", "register"],
        )
        assert ids == {"register": "job_new"}


# ---------------------------------------------------------------------
# Full 5-stage promotion — happy path.
# ---------------------------------------------------------------------
class TestFullRunPromotion:
    def test_all_symlinks_updated_to_fresh(self, he_root: Path):
        _seed_stage(he_root, "warp", "job_new",
                    extra={"source_register_run_id": "job_new"})
        _seed_stage(he_root, "register", "job_new")
        _seed_stage(he_root, "celltyped", "job_new")
        _seed_stage(he_root, "viz", "job_new")

        layout = _layout(he_root, "job_new")
        _promote_completed_run(
            layout=layout, cfg={},
            stages_ran=["register", "warp", "celltype", "viz"],
            source_register_run_id="job_new",
        )
        for stage_dirname in ("register", "warp", "celltyped", "viz"):
            link = he_root / "output" / stage_dirname
            assert link.is_symlink()
            assert output_symlink_run_id(link) == "job_new"


# ---------------------------------------------------------------------
# Partial celltype+viz rerun — the load-bearing lineage case Tracy
# called out explicitly ("a celltype+viz rerun may use an older
# register/warp run").
# ---------------------------------------------------------------------
class TestPartialCelltypeVizRerun:
    def test_leaves_register_warp_alone_when_lineage_matches(self, he_root: Path):
        # Seed the current output pointer at register=reg_A + warp=warp_A
        # (consistent lineage) via a prior full set_default_run.
        from hexenium.set_default_run import set_default_run
        set_default_run(
            output_root_he=he_root,
            register_run_id="reg_A",
            warp_run_id="warp_A",
        )
        # Now a celltype+viz rerun lands fresh outputs under job_new,
        # consuming the same warp_A.
        _seed_stage(he_root, "celltyped", "job_new")
        _seed_stage(he_root, "viz", "job_new")
        layout = _layout(he_root, "job_new")
        _promote_completed_run(
            layout=layout,
            cfg={"warp_run_id": "warp_A"},
            stages_ran=["celltype", "viz"],
            source_register_run_id=None,   # register wasn't in this run.
        )
        assert output_symlink_run_id(he_root / "output" / "register") == "reg_A"
        assert output_symlink_run_id(he_root / "output" / "warp") == "warp_A"
        assert output_symlink_run_id(he_root / "output" / "celltyped") == "job_new"
        assert output_symlink_run_id(he_root / "output" / "viz") == "job_new"

    def test_refuses_when_output_warp_not_seeded(self, he_root: Path):
        # No prior set-default-run — output/ doesn't exist. A celltype+viz
        # rerun cannot know which warp the CURRENT default is, so the
        # promote refuses with an actionable message.
        _seed_stage(he_root, "celltyped", "job_new")
        _seed_stage(he_root, "viz", "job_new")
        layout = _layout(he_root, "job_new")
        with pytest.raises(SystemExit) as exc:
            _promote_completed_run(
                layout=layout,
                cfg={"warp_run_id": "warp_A"},
                stages_ran=["celltype", "viz"],
                source_register_run_id=None,
            )
        assert "output/warp does not exist" in str(exc.value)
        assert not (he_root / "output").exists()

    def test_refuses_when_downstream_used_different_warp(self, he_root: Path):
        # output/warp is warp_A, but this celltype+viz rerun consumed
        # warp_B (which also happens to exist on disk). Promoting
        # would leave output/ lineage-orphaned.
        from hexenium.set_default_run import set_default_run
        set_default_run(
            output_root_he=he_root,
            register_run_id="reg_A",
            warp_run_id="warp_A",
        )
        _seed_stage(he_root, "warp", "warp_B",
                    extra={"source_register_run_id": "reg_B"})
        _seed_stage(he_root, "celltyped", "job_new")
        _seed_stage(he_root, "viz", "job_new")
        layout = _layout(he_root, "job_new")
        with pytest.raises(SystemExit) as exc:
            _promote_completed_run(
                layout=layout,
                cfg={"warp_run_id": "warp_B"},
                stages_ran=["celltype", "viz"],
                source_register_run_id=None,
            )
        assert "lineage mismatch" in str(exc.value)
        # output/ untouched — refusal is atomic.
        assert output_symlink_run_id(he_root / "output" / "register") == "reg_A"
        assert output_symlink_run_id(he_root / "output" / "warp") == "warp_A"
        assert not (he_root / "output" / "celltyped").exists()
        assert not (he_root / "output" / "viz").exists()


# ---------------------------------------------------------------------
# --stages warp only — the second lineage-extension case (warp
# without register).
# ---------------------------------------------------------------------
class TestWarpOnlyPromotion:
    def test_refuses_when_output_register_not_seeded(self, he_root: Path):
        # No prior set-default-run + warp-only invocation: the warp
        # consumed register/reg_A but there's no output/register yet,
        # so lineage cannot be validated. Refuse.
        _seed_stage(he_root, "warp", "job_new",
                    extra={"source_register_run_id": "reg_A"})
        layout = _layout(he_root, "job_new")
        with pytest.raises(SystemExit) as exc:
            _promote_completed_run(
                layout=layout, cfg={},
                stages_ran=["warp"],
                source_register_run_id="reg_A",
            )
        assert "output/register does not exist" in str(exc.value)
        assert not (he_root / "output").exists()

    def test_refuses_when_output_register_mismatches_warp_source(self, he_root: Path):
        from hexenium.set_default_run import set_default_run
        set_default_run(
            output_root_he=he_root,
            register_run_id="reg_B",  # NOTE: not the source of warp we're about to promote
        )
        _seed_stage(he_root, "warp", "job_new",
                    extra={"source_register_run_id": "reg_A"})
        layout = _layout(he_root, "job_new")
        with pytest.raises(SystemExit) as exc:
            _promote_completed_run(
                layout=layout, cfg={},
                stages_ran=["warp"],
                source_register_run_id="reg_A",
            )
        assert "lineage mismatch" in str(exc.value)
        # output/register untouched.
        assert output_symlink_run_id(he_root / "output" / "register") == "reg_B"
        assert not (he_root / "output" / "warp").exists()


# ---------------------------------------------------------------------
# The register↔warp check inherited from set_default_run — verifies
# we still catch what the existing validation catches (defence in depth).
# ---------------------------------------------------------------------
class TestInheritedRegisterWarpLineage:
    def test_refuses_full_run_with_inconsistent_source(self, tmp_path: Path):
        # Fresh he_root with a warp that CLAIMS a different source than
        # the register being promoted. set_default_run's own lineage
        # check should reject this.
        root = tmp_path / "he_registration"
        _seed_stage(root, "register", "job_new")
        _seed_stage(root, "warp", "job_new",
                    extra={"source_register_run_id": "reg_OTHER"})
        _seed_stage(root, "celltyped", "job_new")
        _seed_stage(root, "viz", "job_new")
        layout = _layout(root, "job_new")
        with pytest.raises(SystemExit) as exc:
            _promote_completed_run(
                layout=layout, cfg={},
                stages_ran=["register", "warp", "celltype", "viz"],
                source_register_run_id="job_new",
            )
        # Message comes from set_default_run.
        assert "lineage mismatch" in str(exc.value)
        assert not (root / "output").exists()


# ---------------------------------------------------------------------
# _validate_promotion_lineage — direct test of the pipeline-side
# extension helper (surfaces the errors before set_default_run runs).
# ---------------------------------------------------------------------
class TestValidatePromotionLineageDirect:
    def test_ok_when_only_register_promoted(self, he_root: Path):
        layout = _layout(he_root, "job_new")
        _seed_stage(he_root, "register", "job_new")
        _validate_promotion_lineage(
            layout=layout, cfg={},
            stages_ran=["register"], promote_ids={"register": "job_new"},
            source_register_run_id=None,
        )  # no raise expected

    def test_ok_when_downstream_matches_current_output_warp(self, he_root: Path):
        from hexenium.set_default_run import set_default_run
        set_default_run(
            output_root_he=he_root,
            register_run_id="reg_A",
            warp_run_id="warp_A",
        )
        layout = _layout(he_root, "job_new")
        _validate_promotion_lineage(
            layout=layout,
            cfg={"warp_run_id": "warp_A"},
            stages_ran=["celltype", "viz"],
            promote_ids={"celltype": "job_new", "viz": "job_new"},
            source_register_run_id=None,
        )  # no raise expected


# ---------------------------------------------------------------------
# Opt-out — no flag = output/ never touched (integration-shape test
# that exercises the flag gate, without actually running any stage).
# ---------------------------------------------------------------------
class TestOptOut:
    """Opt-out semantics — no flag = ``output/`` never touched even
    across a full ``pipeline.run()`` invocation.

    F2 skeptic finding on ``settylab/TracyY123-nexus#15``: the earlier
    version of this test asserted trivial facts and never called
    ``pipeline.run()`` — so a regression that moved / mis-gated the
    promote hook could pass CI silently. These tests invoke the real
    ``pipeline.run`` with an empty ``--stages`` list (which exercises
    the tail — layout + banner + config snapshot + the promote gate
    — without needing VALIS/HEST). The gate is verified by spying on
    ``_promote_completed_run``.
    """

    def _minimal_pipeline_inputs(self, tmp_path: Path) -> tuple:
        """Shared fixture body — the smallest set of on-disk inputs
        that lets ``pipeline.run(cfg, stages=[])`` run to completion
        (registers no stage, but exercises the tail-of-run promote
        gate we're testing)."""
        from hexenium.config import deep_update, load_default
        he_path = tmp_path / "he.ome.tif"
        he_path.touch()
        bundle = tmp_path / "output-XETG"
        bundle.mkdir()
        output_root = tmp_path / "runs"
        output_root.mkdir()
        cfg = deep_update(load_default(), {
            "sample_id": "S1",
            "output_root": str(output_root),
            "he_path": str(he_path),
            "xenium_bundle": str(bundle),
        })
        return cfg, output_root

    def test_flag_off_pipeline_run_never_calls_promote(
        self, tmp_path: Path, monkeypatch,
    ):
        from hexenium import pipeline
        called: list[dict] = []
        monkeypatch.setattr(
            pipeline, "_promote_completed_run",
            lambda **kw: called.append(kw),
        )

        cfg, output_root = self._minimal_pipeline_inputs(tmp_path)
        # Flag explicitly UNSET.
        assert "set_default_on_success" not in cfg

        rc = pipeline.run(cfg, stages=[], argv=["hexenium", "run"])
        assert rc == 0
        assert called == []  # promote hook never fired
        # No output/ subdir created either.
        assert not (output_root / "S1" / "output").exists()

    def test_flag_off_explicit_false_pipeline_run_never_calls_promote(
        self, tmp_path: Path, monkeypatch,
    ):
        # Same as above but with the key EXPLICITLY set to False —
        # covers the edge case where CLI wiring changes to always
        # write the key.
        from hexenium import pipeline
        called: list[dict] = []
        monkeypatch.setattr(
            pipeline, "_promote_completed_run",
            lambda **kw: called.append(kw),
        )

        cfg, output_root = self._minimal_pipeline_inputs(tmp_path)
        cfg["set_default_on_success"] = False

        pipeline.run(cfg, stages=[], argv=["hexenium", "run"])
        assert called == []
        assert not (output_root / "S1" / "output").exists()

    def test_flag_on_pipeline_run_calls_promote(
        self, tmp_path: Path, monkeypatch,
    ):
        # Parity check: with the flag ON, ``pipeline.run`` reaches the
        # promote hook (regardless of what it does inside — we spy on
        # entry, not on ``set_default_run`` behavior which has its
        # own suite).
        from hexenium import pipeline
        called: list[dict] = []
        monkeypatch.setattr(
            pipeline, "_promote_completed_run",
            lambda **kw: called.append(kw),
        )

        cfg, output_root = self._minimal_pipeline_inputs(tmp_path)
        cfg["set_default_on_success"] = True

        pipeline.run(cfg, stages=[], argv=["hexenium", "run"])
        assert len(called) == 1
        # Kwargs shape as the real hook expects.
        assert set(called[0].keys()) == {
            "layout", "cfg", "stages_ran", "source_register_run_id",
        }
        assert called[0]["stages_ran"] == []
        assert called[0]["source_register_run_id"] is None

    def test_promote_helper_still_writes_output_when_invoked_directly(
        self, he_root: Path,
    ):
        # Kept from the previous version — verifies the helper itself
        # works when its inputs are valid. Complements the gate tests
        # above (gate + downstream, both covered).
        _seed_stage(he_root, "register", "job_new")
        _seed_stage(he_root, "warp", "job_new",
                    extra={"source_register_run_id": "job_new"})
        layout = _layout(he_root, "job_new")
        _promote_completed_run(
            layout=layout, cfg={"set_default_on_success": True},
            stages_ran=["register", "warp"],
            source_register_run_id="job_new",
        )
        assert (he_root / "output" / "register").is_symlink()
        assert (he_root / "output" / "warp").is_symlink()


# ---------------------------------------------------------------------
# Structural / regression guards — the flag's contract lives partly
# in code structure (the hook is the LAST thing in run()) and partly
# in data structure (_PROMOTE_STAGE_DIR shape).
# ---------------------------------------------------------------------
class TestStructuralGuards:
    def test_promote_stage_dir_shape_stable(self):
        # If a stage is added/renamed, this test forces a conscious
        # update of the promote surface (set_default_run kwargs would
        # otherwise silently drop).
        assert _PROMOTE_STAGE_DIR == {
            "register": "register",
            "warp": "warp",
            "celltype": "celltyped",
            "viz": "viz",
        }

    def test_hook_is_last_in_pipeline_run(self):
        # Structural guarantee: the ``--set-default-on-success`` hook
        # is the LAST statement in ``pipeline.run()`` — a stage error
        # raises before it, so ``output/`` stays put on stage failure.
        # This is not exhaustively enforceable, but grepping for the
        # hook in the tail of ``run()`` catches the common refactor
        # accident of moving it earlier.
        import inspect
        from hexenium import pipeline
        src = inspect.getsource(pipeline.run)
        hook_line = src.rfind("_promote_completed_run(")
        done_line = src.rfind('hexenium done in')
        assert hook_line > 0, "promote hook missing from pipeline.run()"
        assert hook_line < done_line, (
            "promote hook must run before the final banner (both live "
            "in the tail of run()); code path moved unexpectedly."
        )
