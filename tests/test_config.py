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
        "warp", "celltype", "nn_celltype_mapping", "viz",
    ):
        assert section in cfg, f"missing section: {section}"


def test_default_yaml_has_required_stubs():
    """The four required-at-runtime keys are present in default.yaml but
    null — user must set them via CLI/user YAML."""
    from hexenium.config import REQUIRED_KEYS, load_default

    cfg = load_default()
    for key in REQUIRED_KEYS:
        assert key in cfg, f"missing required-key stub: {key}"


def test_validate_rejects_missing_required():
    """Validation raises when a required key is null."""
    import pytest
    from hexenium.config import validate

    with pytest.raises(SystemExit):
        validate({})


def test_deep_update_merges_recursively():
    from hexenium.config import deep_update

    base = {"a": 1, "nested": {"x": 10, "y": 20}}
    override = {"nested": {"y": 99, "z": 100}, "b": 2}
    got = deep_update(base, override)
    assert got == {"a": 1, "b": 2, "nested": {"x": 10, "y": 99, "z": 100}}
