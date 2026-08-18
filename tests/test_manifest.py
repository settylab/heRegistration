"""Unit tests for hexenium.manifest — the shared per-run manifest,
params_hash, atomic-symlink, and lineage helpers used by the register
and warp stages plus the set-default-run subcommand.

Run:
    pytest tests/test_manifest.py -v
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from hexenium.manifest import (
    compute_params_hash,
    output_symlink_run_id,
    read_manifest,
    read_symlink_target,
    set_default_symlink,
    symlink_target_run_id,
    utc_timestamp,
    write_manifest,
)


# ---------------------------------------------------------------------
# compute_params_hash — must be stable across key order + Path types.
# Reproducibility hinges on the invariant "same effective params →
# same digest," so semantically-equivalent invocations MUST agree.
# ---------------------------------------------------------------------
class TestComputeParamsHash:
    def test_deterministic_across_calls(self):
        params = {"mode": "rigid_only", "max_image_dim_px": 1500}
        assert compute_params_hash(params) == compute_params_hash(params)

    def test_key_order_insensitive(self):
        # This is the surface where a naive `json.dumps` (no sort_keys)
        # would silently produce two hashes for the same params.
        a = {"mode": "rigid_only", "max_image_dim_px": 1500}
        b = {"max_image_dim_px": 1500, "mode": "rigid_only"}
        assert compute_params_hash(a) == compute_params_hash(b)

    def test_different_values_produce_different_hashes(self):
        a = {"mode": "rigid_only"}
        b = {"mode": "rigid_nonrigid"}
        assert compute_params_hash(a) != compute_params_hash(b)

    def test_path_values_serialise_via_str_default(self):
        # Paths mustn't crash the hash, and their string form should be
        # what's hashed — so `Path("/a") == "/a"` at the digest level.
        via_path = compute_params_hash({"input": Path("/a/b")})
        via_str = compute_params_hash({"input": "/a/b"})
        assert via_path == via_str

    def test_hex_length_is_stable(self):
        # If the truncation length ever changes, downstream greps break;
        # this test locks in the current contract.
        h = compute_params_hash({"foo": 1})
        assert len(h) == 16
        assert all(c in "0123456789abcdef" for c in h)


# ---------------------------------------------------------------------
# write_manifest / read_manifest — atomic write + symlink transparency.
# ---------------------------------------------------------------------
class TestManifestIO:
    def test_write_then_read_roundtrip(self, tmp_path):
        path = tmp_path / "manifest.yaml"
        payload = {"sample_id": "S1", "params_hash": "abc123",
                   "params": {"mode": "rigid_only"}}
        write_manifest(path, payload)
        assert read_manifest(path) == payload

    def test_read_follows_symlink(self, tmp_path):
        real = tmp_path / "run_A" / "manifest.yaml"
        write_manifest(real, {"sample_id": "S1", "he_job_id": "A"})
        link = tmp_path / "manifest.yaml"
        os.symlink("run_A/manifest.yaml", link)
        assert read_manifest(link)["he_job_id"] == "A"

    def test_read_missing_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            read_manifest(tmp_path / "nope.yaml")

    def test_atomic_write_no_dangling_tmp(self, tmp_path):
        # The atomic-write path uses <name>.tmp.<pid>; a successful
        # write should not leave that tmp file behind.
        path = tmp_path / "manifest.yaml"
        write_manifest(path, {"foo": "bar"})
        stray = list(tmp_path.glob("manifest.yaml.tmp.*"))
        assert stray == []


# ---------------------------------------------------------------------
# set_default_symlink — the primitive that both register/warp writers
# and set-default-run share. Must be safe to run repeatedly.
# ---------------------------------------------------------------------
class TestSetDefaultSymlink:
    def test_creates_new_symlink(self, tmp_path):
        link = tmp_path / "manifest.yaml"
        set_default_symlink(link, "run_A/manifest.yaml")
        assert link.is_symlink()
        assert read_symlink_target(link) == "run_A/manifest.yaml"

    def test_repoints_existing_symlink(self, tmp_path):
        link = tmp_path / "manifest.yaml"
        set_default_symlink(link, "run_A/manifest.yaml")
        set_default_symlink(link, "run_B/manifest.yaml")
        assert read_symlink_target(link) == "run_B/manifest.yaml"

    def test_idempotent(self, tmp_path):
        link = tmp_path / "manifest.yaml"
        set_default_symlink(link, "run_A/manifest.yaml")
        set_default_symlink(link, "run_A/manifest.yaml")
        assert read_symlink_target(link) == "run_A/manifest.yaml"

    def test_refuses_to_clobber_real_file(self, tmp_path):
        # Guard against silently destroying operator-managed content —
        # e.g. someone hand-copied a real manifest.yaml here.
        real = tmp_path / "manifest.yaml"
        real.write_text("hand-written content\n")
        with pytest.raises(FileExistsError):
            set_default_symlink(real, "run_A/manifest.yaml")
        assert real.read_text() == "hand-written content\n"


# ---------------------------------------------------------------------
# symlink_target_run_id / output_symlink_run_id — the two schemas.
# ---------------------------------------------------------------------
class TestSymlinkRunIdHelpers:
    def test_stage_root_scheme(self, tmp_path):
        link = tmp_path / "manifest.yaml"
        os.symlink("run_A/manifest.yaml", link)
        assert symlink_target_run_id(link) == "run_A"

    def test_output_scheme(self, tmp_path):
        # output/<stage> -> ../<stage>/<he_job_id>
        (tmp_path / "output").mkdir()
        link = tmp_path / "output" / "register"
        os.symlink("../register/run_A", link)
        assert output_symlink_run_id(link) == "run_A"

    def test_helpers_return_none_when_no_link(self, tmp_path):
        assert symlink_target_run_id(tmp_path / "nope") is None
        assert output_symlink_run_id(tmp_path / "nope") is None


# ---------------------------------------------------------------------
# Timestamp shape — sortable + timezone-explicit.
# ---------------------------------------------------------------------
class TestUtcTimestamp:
    def test_iso_format_with_z_suffix(self):
        ts = utc_timestamp()
        # 2026-08-18T12:34:56Z
        assert len(ts) == 20
        assert ts[4] == "-"
        assert ts[10] == "T"
        assert ts.endswith("Z")
