"""YAML config loading, merging, and validation for hexenium.

The canonical default config ships inside the package at
``src/hexenium/_defaults/default.yaml`` and is loaded via
``importlib.resources`` — so both editable AND non-editable
(``pip install .`` from a clone or a wheel) installs resolve it
correctly. Cross-user finding on
``settylab/TracyY123-nexus#15`` comment ``5360185231``: the prior
layout kept the yaml at the repo-root ``config/`` and referenced it
via ``Path(__file__).resolve().parents[2]``, which walks OUT of the
package tree in a wheel install (``site-packages/hexenium/`` has no
sibling ``config/`` directory). The `pyproject.toml`'s
``../../config/default.yaml`` package-data pattern is not valid
setuptools syntax either, so the wheel simply omitted the file.

``load_default()`` returns it as a dict; ``load_yaml(path)`` returns
a user override YAML; ``deep_update(base, override)`` does a
recursive merge. ``validate(cfg)`` raises with a readable list of
missing required top-level keys.

Required-key semantics track the three-mode CLI:

* standalone AND integrated-by-run-id — need ``sample_id`` +
  ``output_root`` (plus ``he_path`` + ``xenium_bundle`` always).
* integrated-by-h5ad — ``xenium_h5ad`` is present, so ``sample_id`` +
  ``output_root`` are optional; identity comes from ``.uns``.

``validate()`` inspects the config to decide which mode's required-set
applies.
"""
from __future__ import annotations

from importlib.resources import files as _pkg_files
from pathlib import Path

import yaml

# Package layout:
#   src/hexenium/config.py                <- this file
#   src/hexenium/_defaults/__init__.py    <- makes it a subpackage so
#                                            importlib.resources can find it
#   src/hexenium/_defaults/default.yaml   <- the shipped default config
#
# ``importlib.resources.files("hexenium._defaults")`` returns a
# ``Traversable`` — a file-like anchor that works for editable
# installs (source layout) AND non-editable installs (zipfile / wheel)
# without touching sys.path or __file__.
_DEFAULTS_PKG = "hexenium._defaults"
_DEFAULT_YAML_NAME = "default.yaml"

#: Config keys required at ALL invocation modes.
CORE_REQUIRED_KEYS = ("he_path", "xenium_bundle")
#: Additional keys required in standalone / integrated-by-run-id modes
#: (integrated-by-h5ad derives them from .uns).
STANDALONE_REQUIRED_KEYS = ("sample_id", "output_root")
#: Retained for tests / external callers that want the pre-refactor
#: "everything required" set (standalone view).
REQUIRED_KEYS = CORE_REQUIRED_KEYS + STANDALONE_REQUIRED_KEYS

#: Canonical stage names, in execution order. Kept aligned with
#: :data:`hexenium.layout.STAGE_NAMES`.
VALID_STAGES = (
    "he_preprocess", "register", "warp", "celltype", "viz",
)
DEFAULT_STAGES = VALID_STAGES


def deep_update(base: dict, override: dict) -> dict:
    """Recursive dict merge — override wins."""
    out = dict(base)
    for k, v in override.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = deep_update(out[k], v)
        else:
            out[k] = v
    return out


def load_yaml(path: Path) -> dict:
    """Read a YAML file, return `{}` if the file is empty."""
    with open(path) as f:
        return yaml.safe_load(f) or {}


def default_config_path() -> Path:
    """Return the on-disk path to the shipped default config.

    Uses ``importlib.resources`` so the path resolves correctly for
    BOTH editable installs (points at the source-tree location) and
    non-editable wheel installs (points at the copy shipped inside
    the wheel). Callers should treat the result as read-only.
    """
    # ``Traversable`` supports ``.__fspath__`` on real-filesystem
    # backends (the common case). For unusual backends (e.g. zipfile),
    # cast via ``str()`` since ``Path`` accepts strings.
    resource = _pkg_files(_DEFAULTS_PKG).joinpath(_DEFAULT_YAML_NAME)
    return Path(str(resource))


#: Retained as a module attribute for callers that peek at the path
#: (tests, log messages, ``--help`` output). Wraps the ``importlib``
#: resolution so the value is stable at import time.
DEFAULT_CONFIG_PATH = default_config_path()


def load_default() -> dict:
    """Load the package-shipped default config via ``importlib.resources``.

    Works in both editable AND non-editable installs — see the module
    docstring for the fix history.
    """
    resource = _pkg_files(_DEFAULTS_PKG).joinpath(_DEFAULT_YAML_NAME)
    # ``Traversable.read_text`` is the portable read API across
    # filesystem and zipfile backends. Keep it here rather than
    # dispatching through ``load_yaml`` so a zipped install works too.
    return yaml.safe_load(resource.read_text()) or {}


def validate(cfg: dict) -> None:
    """Raise SystemExit if a required key is missing/null.

    Mode-aware: when ``xenium_h5ad`` is set (integrated-by-h5ad),
    ``sample_id`` + ``output_root`` are optional (identity comes from
    the h5ad's ``.uns``). Otherwise they're required.
    """
    missing = [k for k in CORE_REQUIRED_KEYS if not cfg.get(k)]
    if not cfg.get("xenium_h5ad"):
        missing += [k for k in STANDALONE_REQUIRED_KEYS if not cfg.get(k)]
    if missing:
        raise SystemExit(f"missing required config keys: {missing}")
