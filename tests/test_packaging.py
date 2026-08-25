"""Packaging regression tests.

Cross-user finding on ``settylab/TracyY123-nexus#15`` comment
``5360185231``: `config/default.yaml` didn't ship in a non-editable
`pip install .` because (a) it lived at the repo-root
``config/default.yaml`` (outside the package tree), and (b)
``pyproject.toml``'s ``["../../config/default.yaml"]`` package-data
pattern is invalid setuptools syntax. Editable installs "worked"
purely because `Path(__file__).parents[2]` in the source tree
happened to walk up to the repo root; wheel installs broke.

Fix: yaml moved to ``src/hexenium/_defaults/default.yaml`` and the
loader now uses ``importlib.resources.files("hexenium._defaults")``.

These tests lock the fix in — they do NOT need a fresh non-editable
install to run (they exercise the code path directly), but the
verification section of the follow-up report also runs `pip install
--no-deps .` in a fresh venv to belt-and-braces the packaging.
"""
from __future__ import annotations

import pytest


def test_default_yaml_ships_as_package_resource():
    """``hexenium._defaults`` package exposes ``default.yaml`` via
    ``importlib.resources``. This is the API that survives both
    editable and non-editable installs — if this test passes on the
    source layout AND in a wheel, both install modes work."""
    from importlib.resources import files
    resource = files("hexenium._defaults").joinpath("default.yaml")
    assert resource.is_file(), (
        f"expected default.yaml to ship inside hexenium._defaults; "
        f"got resource={resource!r} that isn't a file. Likely "
        f"pyproject.toml [tool.setuptools.package-data] regressed."
    )
    text = resource.read_text()
    assert text.strip(), "default.yaml is empty; packaging shipped an empty file"


def test_load_default_returns_populated_dict():
    """``load_default()`` reads the yaml + returns the known top-level
    schema. Guards against the historic path-resolution bug (walked
    out of the package tree in wheel installs) coming back."""
    from hexenium.config import load_default
    cfg = load_default()
    assert isinstance(cfg, dict) and cfg, (
        "load_default() returned empty; the yaml either didn't ship "
        "or the importlib.resources anchor is wrong."
    )
    # Sample of top-level keys the pipeline depends on — presence is
    # enough (values may drift; the shape is what we lock).
    for key in ("registration", "warp", "celltype", "viz"):
        assert key in cfg, (
            f"default.yaml is missing top-level {key!r} — did we ship "
            f"the wrong file, or overwrite it with a stub?"
        )


def test_default_config_path_module_attr_is_stable():
    """``DEFAULT_CONFIG_PATH`` module attribute survives — some
    external callers / log lines reference it. Locks in the
    contract that it resolves to a real file after the packaging
    refactor."""
    from hexenium.config import DEFAULT_CONFIG_PATH
    from pathlib import Path
    assert isinstance(DEFAULT_CONFIG_PATH, Path)
    assert DEFAULT_CONFIG_PATH.exists(), (
        f"DEFAULT_CONFIG_PATH={DEFAULT_CONFIG_PATH} does not resolve to "
        f"a real file. Wheel install likely missing the resource."
    )
    # Path ends at the shipped yaml, wherever the anchor resolves.
    assert DEFAULT_CONFIG_PATH.name == "default.yaml"


def test_config_module_never_references_repo_root_walk():
    """Historic regression guard: the old loader used
    ``Path(__file__).resolve().parents[2] / "config" / "default.yaml"``
    to reach the repo-root config. That walk breaks in wheel installs
    (``site-packages/hexenium/../..`` isn't the repo root). Guard
    against the fragile pattern coming back in ACTUAL CODE — parse
    the module's AST and inspect subscript expressions, so the
    docstring's narrative mention of ``parents[2]`` (as documentation
    of the historic bug) doesn't false-positive the check.
    """
    import ast
    import inspect
    from hexenium import config as _cfg
    tree = ast.parse(inspect.getsource(_cfg))
    for node in ast.walk(tree):
        # Match ``.parents[<INT>]`` subscript accesses in code.
        if not isinstance(node, ast.Subscript):
            continue
        val = node.value
        if not isinstance(val, ast.Attribute) or val.attr != "parents":
            continue
        # Some Python versions wrap the index in ast.Index; unwrap.
        idx_node = node.slice
        if hasattr(ast, "Index") and isinstance(idx_node, ast.Index):
            idx_node = idx_node.value
        if isinstance(idx_node, ast.Constant) and idx_node.value == 2:
            raise AssertionError(
                "hexenium.config regressed: contains a "
                "`.parents[2]` code expression. That walk lands "
                "outside site-packages in a wheel install — use "
                "importlib.resources.files(...) instead."
            )
    # Positive guard: importlib.resources.files IS in use.
    src = inspect.getsource(_cfg)
    assert "files(" in src and "_defaults" in src, (
        "hexenium.config no longer uses importlib.resources / the "
        "_defaults subpackage anchor — the packaging fix has been "
        "reverted."
    )
