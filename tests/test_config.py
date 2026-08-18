"""Smoke test: default YAML loads + carries the required top-level keys."""
from __future__ import annotations


def test_default_yaml_loads():
    from hexenium.config import load_default

    cfg = load_default()
    assert isinstance(cfg, dict)
    assert cfg  # non-empty


def test_default_yaml_has_expected_top_level_sections():
    from hexenium.config import load_default

    cfg = load_default()
    for section in (
        "he_preprocess", "he", "registration", "parameters",
        "warp", "celltype", "viz",
    ):
        assert section in cfg, f"missing section: {section}"


def test_default_yaml_has_required_stubs():
    """Required top-level keys are present in default.yaml but null —
    user must set them via CLI/user YAML (except in integrated-by-h5ad
    mode where identity comes from .uns)."""
    from hexenium.config import REQUIRED_KEYS, load_default

    cfg = load_default()
    for key in REQUIRED_KEYS:
        assert key in cfg, f"missing required-key stub: {key}"


def test_default_yaml_has_integrated_mode_stubs():
    """Integrated-mode fields (xenium_h5ad, proseg_purified_h5ad, he_job_id,
    warp_run_id, celltype_run_id) are stubbed to null so users can enable
    integrated mode by setting them on the CLI."""
    from hexenium.config import load_default

    cfg = load_default()
    for key in ("xenium_h5ad", "proseg_purified_h5ad", "he_job_id",
                "warp_run_id", "celltype_run_id"):
        assert key in cfg, f"missing integrated-mode stub: {key}"
        assert cfg[key] is None, f"{key} should default to null"


def test_validate_rejects_missing_required_standalone():
    """Validation raises when a required standalone key is null."""
    import pytest
    from hexenium.config import validate

    with pytest.raises(SystemExit):
        validate({})


def test_validate_accepts_integrated_mode_without_output_root():
    """In integrated-by-h5ad mode, sample_id + output_root are optional
    (identity read from .uns). Only he_path + xenium_bundle are required."""
    from hexenium.config import validate

    # Not raising is the assertion here.
    validate({
        "xenium_h5ad": "/path/to.h5ad",
        "he_path": "/path/to/he.ome.tif",
        "xenium_bundle": "/path/to/bundle",
    })


def test_deep_update_merges_recursively():
    from hexenium.config import deep_update

    base = {"a": 1, "nested": {"x": 10, "y": 20}}
    override = {"nested": {"y": 99, "z": 100}, "b": 2}
    got = deep_update(base, override)
    assert got == {"a": 1, "b": 2, "nested": {"x": 10, "y": 99, "z": 100}}
