"""Regression tests for hexenium.stages.registration._canonical_he_symlink.

Anchors the fix for `#15` — filenames with extra dots before the
image extension (e.g. ``MH9_0088754_SR25-1407.40x_BF_01.ome.tif``)
used to produce ``aligned_fullres_HE.40x_BF_01.ome.tif``, which broke
valis_hest's hardcoded ``slide_dict['aligned_fullres_HE']`` lookup at
``valis_hest/registration.py:4644``.

Run:
    pytest tests/test_canonical_he_symlink.py -v
"""
from __future__ import annotations

from pathlib import Path

import pytest

from hexenium.stages.registration import (
    _canonical_he_suffix,
    _canonical_he_symlink,
)


# ---------------------------------------------------------------------
# _canonical_he_suffix — pure suffix picker, no filesystem
# ---------------------------------------------------------------------
class TestCanonicalHeSuffix:
    def test_plain_ome_tif(self):
        assert _canonical_he_suffix("sample.ome.tif") == ".ome.tif"

    def test_multi_dot_ome_tif_mh9(self):
        # The concrete case from `#15` slurm-5541605.
        assert (
            _canonical_he_suffix("MH9_0088754_SR25-1407.40x_BF_01.ome.tif")
            == ".ome.tif"
        )

    def test_many_middle_dots(self):
        assert (
            _canonical_he_suffix("sample_with.dots.and.more.dots.ome.tif")
            == ".ome.tif"
        )

    def test_bigtiff_middle_dot(self):
        # TMA1A bladder case: `BottomSubset_offaligned_core.bigtiff.ome.tif`.
        # `.bigtiff` is NOT a known suffix, so the picker must walk past
        # it and land on `.ome.tif`.
        assert (
            _canonical_he_suffix("BottomSubset_offaligned_core.bigtiff.ome.tif")
            == ".ome.tif"
        )

    def test_ome_tiff_double_f(self):
        assert _canonical_he_suffix("sample.ome.tiff") == ".ome.tiff"

    def test_plain_tif(self):
        assert _canonical_he_suffix("sample.tif") == ".tif"

    def test_plain_tiff(self):
        assert _canonical_he_suffix("sample.tiff") == ".tiff"

    def test_qptiff(self):
        assert _canonical_he_suffix("scan.qptiff") == ".qptiff"

    def test_svs(self):
        assert _canonical_he_suffix("slide.svs") == ".svs"

    def test_uppercase_ome_tif_normalises(self):
        # Case-insensitive match; the canonical suffix is always the
        # lowercase form so the emitted symlink is stable regardless
        # of source casing.
        assert _canonical_he_suffix("SAMPLE.OME.TIF") == ".ome.tif"

    def test_mixed_case_tif(self):
        assert _canonical_he_suffix("Sample.Tif") == ".tif"

    def test_unknown_suffix_raises(self):
        with pytest.raises(ValueError, match="no recognised image suffix"):
            _canonical_he_suffix("data.h5ad")

    def test_no_suffix_raises(self):
        with pytest.raises(ValueError, match="no recognised image suffix"):
            _canonical_he_suffix("bare_name")

    def test_ome_tif_wins_over_tif(self):
        # Explicit ordering guard: `.ome.tif` must be matched before
        # `.tif`, otherwise plain `.tif` would win and leave `.ome` in
        # the canonical stem.
        assert _canonical_he_suffix("x.ome.tif") == ".ome.tif"


# ---------------------------------------------------------------------
# _canonical_he_symlink — filesystem integration
# ---------------------------------------------------------------------
class TestCanonicalHeSymlink:
    def test_multi_dot_produces_canonical_name(self, tmp_path: Path):
        # Reproduces Tracy's MH9 case end-to-end.
        source = tmp_path / "src" / "MH9_0088754_SR25-1407.40x_BF_01.ome.tif"
        source.parent.mkdir()
        source.write_bytes(b"fake tiff")
        workdir = tmp_path / "workdir"

        link = _canonical_he_symlink(source, workdir)

        assert link.name == "aligned_fullres_HE.ome.tif"
        assert link.is_symlink()
        assert link.resolve() == source.resolve()

    def test_plain_ome_tif_regression(self, tmp_path: Path):
        # Pre-existing single-dot case must still resolve identically.
        source = tmp_path / "src" / "plain_sample.ome.tif"
        source.parent.mkdir()
        source.write_bytes(b"fake tiff")
        workdir = tmp_path / "workdir"

        link = _canonical_he_symlink(source, workdir)

        assert link.name == "aligned_fullres_HE.ome.tif"
        assert link.resolve() == source.resolve()

    def test_plain_tif_preserves_single_extension(self, tmp_path: Path):
        source = tmp_path / "src" / "sample.tif"
        source.parent.mkdir()
        source.write_bytes(b"fake tiff")
        workdir = tmp_path / "workdir"

        link = _canonical_he_symlink(source, workdir)

        assert link.name == "aligned_fullres_HE.tif"

    def test_idempotent_second_call_replaces_link(self, tmp_path: Path):
        source = tmp_path / "src" / "sample.ome.tif"
        source.parent.mkdir()
        source.write_bytes(b"fake tiff")
        workdir = tmp_path / "workdir"

        first = _canonical_he_symlink(source, workdir)
        second = _canonical_he_symlink(source, workdir)

        assert first == second
        assert second.is_symlink()
        assert second.resolve() == source.resolve()

    def test_replaces_stale_link_pointing_elsewhere(self, tmp_path: Path):
        # If the workdir already holds a stale canonical link pointing
        # at an old H&E, a fresh call must repoint it at the new source.
        old = tmp_path / "src" / "old.ome.tif"
        new = tmp_path / "src" / "new.ome.tif"
        old.parent.mkdir()
        old.write_bytes(b"old tiff")
        new.write_bytes(b"new tiff")
        workdir = tmp_path / "workdir"

        _canonical_he_symlink(old, workdir)
        link = _canonical_he_symlink(new, workdir)

        assert link.resolve() == new.resolve()

    def test_unknown_extension_raises(self, tmp_path: Path):
        source = tmp_path / "src" / "data.h5ad"
        source.parent.mkdir()
        source.write_bytes(b"not an image")
        workdir = tmp_path / "workdir"

        with pytest.raises(ValueError, match="no recognised image suffix"):
            _canonical_he_symlink(source, workdir)
        # The workdir may be created by the mkdir() call before the
        # suffix check — but the canonical symlink must NOT exist.
        assert not (workdir / "aligned_fullres_HE.h5ad").exists()
