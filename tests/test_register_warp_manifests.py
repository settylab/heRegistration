"""Unit tests for the per-stage manifest writers.

These exercise the module-private ``_write_run_manifests`` helpers in
``registration.py`` and ``warp.py``, and the pipeline's
``_resolve_existing_registrar``. Full stage execution is out of scope
(both stages call VALIS/HEST) — this file locks in the manifest
schema + symlink invariants that set-default-run depends on.

Run:
    pytest tests/test_register_warp_manifests.py -v
"""
from __future__ import annotations

from pathlib import Path

import pytest

from hexenium.manifest import (
    output_symlink_run_id,
    read_manifest,
    symlink_target_run_id,
)


# ---------------------------------------------------------------------
# register stage: registration._write_run_manifests
# ---------------------------------------------------------------------
class TestRegisterManifest:
    def test_writes_per_run_and_symlink(self, tmp_path: Path):
        from hexenium.stages.registration import _write_run_manifests

        register_root = tmp_path / "register"
        registrar_dir = register_root / "reg_A"
        registrar_dir.mkdir(parents=True)
        _write_run_manifests(
            registrar_dir=registrar_dir,
            container=register_root,
            payload={"sample_id": "S1", "he_job_id": "reg_A",
                     "params_hash": "abc123"},
        )
        # Per-run file exists AND has the expected content.
        per_run = registrar_dir / "manifest.yaml"
        assert per_run.exists()
        assert read_manifest(per_run)["he_job_id"] == "reg_A"
        # Root symlink points at the per-run.
        root_link = register_root / "manifest.yaml"
        assert root_link.is_symlink()
        assert symlink_target_run_id(root_link) == "reg_A"
        # Reading through the symlink returns the same manifest.
        assert read_manifest(root_link)["he_job_id"] == "reg_A"

    def test_second_run_repoints_root_but_preserves_first(self, tmp_path: Path):
        # The critical invariant: reg_B's write must NOT destroy reg_A's
        # manifest — that's what makes lineage checks reproducible.
        from hexenium.stages.registration import _write_run_manifests

        register_root = tmp_path / "register"
        for run in ("reg_A", "reg_B"):
            d = register_root / run
            d.mkdir(parents=True)
            _write_run_manifests(
                registrar_dir=d,
                container=register_root,
                payload={"sample_id": "S1", "he_job_id": run,
                         "params_hash": f"hash-{run}"},
            )
        assert (register_root / "reg_A" / "manifest.yaml").exists()
        assert (register_root / "reg_B" / "manifest.yaml").exists()
        # Root symlink now points at reg_B (last wins).
        assert symlink_target_run_id(register_root / "manifest.yaml") == "reg_B"


# ---------------------------------------------------------------------
# warp stage: warp._write_run_manifests
# ---------------------------------------------------------------------
class TestWarpManifest:
    def test_records_source_register_run_id(self, tmp_path: Path):
        from hexenium.stages.warp import _write_run_manifests

        warp_root = tmp_path / "warp"
        out_dir = warp_root / "warp_A"
        out_dir.mkdir(parents=True)
        _write_run_manifests(
            out_dir=out_dir,
            sample_id="S1",
            targets=["cells", "nuclei"],
            use_dask=True,
            save_geojson=True,
            source_register_run_id="reg_A",
            registrar_pickle=Path("/tmp/_registrar.pickle"),
        )
        m = read_manifest(out_dir / "manifest.yaml")
        assert m["source_register_run_id"] == "reg_A"
        assert m["he_job_id"] == "warp_A"
        assert sorted(m["params"]["targets"]) == ["cells", "nuclei"]
        # Root symlink present.
        assert symlink_target_run_id(warp_root / "manifest.yaml") == "warp_A"

    def test_source_register_run_id_none_is_serialisable(self, tmp_path: Path):
        # Path B without --register-run-id (older/legacy) — the field
        # must land as YAML null, not blow up serialisation.
        from hexenium.stages.warp import _write_run_manifests

        out_dir = tmp_path / "warp" / "warp_legacy"
        out_dir.mkdir(parents=True)
        _write_run_manifests(
            out_dir=out_dir,
            sample_id="S1",
            targets=["cells"],
            use_dask=False,
            save_geojson=False,
            source_register_run_id=None,
            registrar_pickle=Path("/tmp/_registrar.pickle"),
        )
        m = read_manifest(out_dir / "manifest.yaml")
        assert m["source_register_run_id"] is None


# ---------------------------------------------------------------------
# pipeline._resolve_existing_registrar
# ---------------------------------------------------------------------
class TestResolveExistingRegistrar:
    def _seed_registration(self, register_root: Path, run: str,
                           pickle_path: Path) -> None:
        d = register_root / run
        d.mkdir(parents=True)
        from hexenium.manifest import set_default_symlink, write_manifest
        write_manifest(d / "manifest.yaml", {
            "sample_id": "S1",
            "he_job_id": run,
            "registrar_pickle": str(pickle_path),
        })
        set_default_symlink(register_root / "manifest.yaml",
                            f"{run}/manifest.yaml")

    def test_default_symlink_path(self, tmp_path: Path):
        from hexenium.pipeline import _resolve_existing_registrar

        register_root = tmp_path / "register"
        self._seed_registration(register_root, "reg_A",
                                pickle_path=tmp_path / "pA.pkl")
        pickle, source = _resolve_existing_registrar(
            register_root=register_root, explicit_register_run_id=None,
        )
        assert pickle == tmp_path / "pA.pkl"
        assert source == "reg_A"

    def test_explicit_flag_bypasses_default(self, tmp_path: Path):
        from hexenium.pipeline import _resolve_existing_registrar

        register_root = tmp_path / "register"
        # Default points at reg_A, but we explicitly ask for reg_B.
        self._seed_registration(register_root, "reg_A",
                                pickle_path=tmp_path / "pA.pkl")
        self._seed_registration(register_root, "reg_B",
                                pickle_path=tmp_path / "pB.pkl")
        # (reg_B write repoints the default; ensure we still bypass it.)
        pickle, source = _resolve_existing_registrar(
            register_root=register_root, explicit_register_run_id="reg_A",
        )
        assert pickle == tmp_path / "pA.pkl"
        assert source == "reg_A"

    def test_explicit_flag_fails_loud_on_typo(self, tmp_path: Path):
        # If the operator mistypes --register-run-id, we MUST NOT fall
        # back to the default — that would poison downstream lineage.
        from hexenium.pipeline import _resolve_existing_registrar

        register_root = tmp_path / "register"
        self._seed_registration(register_root, "reg_A",
                                pickle_path=tmp_path / "pA.pkl")
        with pytest.raises(SystemExit) as exc:
            _resolve_existing_registrar(
                register_root=register_root,
                explicit_register_run_id="reg_TYPO",
            )
        assert "reg_TYPO" in str(exc.value)

    def test_empty_register_returns_none(self, tmp_path: Path):
        # No prior registration and no explicit flag — return (None,
        # None) so the caller's downstream guard raises the actionable
        # "run register first" error.
        from hexenium.pipeline import _resolve_existing_registrar

        register_root = tmp_path / "register"
        register_root.mkdir()
        pickle, source = _resolve_existing_registrar(
            register_root=register_root, explicit_register_run_id=None,
        )
        assert pickle is None
        assert source is None


# Unused shape helper used by TestSetDefaultRun — kept here so its
# import doesn't fail when running this file alone.
_ = output_symlink_run_id
