"""YAML config loading, merging, and validation for hexenium.

The canonical default config lives at `<repo-root>/config/default.yaml`.
`load_default()` returns it as a dict; `load_yaml(path)` returns a user
override YAML; `deep_update(base, override)` does a recursive merge.
`validate(cfg)` raises with a readable list of missing required top-level
keys.

Required-key semantics track the three-mode CLI:

* standalone AND integrated-by-run-id — need ``sample_id`` +
  ``output_root`` (plus ``he_path`` + ``xenium_bundle`` always).
* integrated-by-h5ad — ``xenium_h5ad`` is present, so ``sample_id`` +
  ``output_root`` are optional; identity comes from ``.uns``.

``validate()`` inspects the config to decide which mode's required-set
applies.
"""
from __future__ import annotations

from pathlib import Path

import yaml

# Package layout:
#   src/hexenium/config.py         <- this file
#   src/hexenium/
#   config/default.yaml            <- default config (repo top-level)
#
# From this file, walk up: parents[0]=hexenium, [1]=src, [2]=repo root.
_REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = _REPO_ROOT / "config" / "default.yaml"

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


def load_default() -> dict:
    """Load the package-shipped default config."""
    return load_yaml(DEFAULT_CONFIG_PATH)


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
