"""YAML config loading, merging, and validation for hexenium.

The canonical default config lives at `<package>/config/default.yaml`
(top-level of the source tree). `load_default()` returns it as a dict;
`load_user(path)` returns a user override YAML; `deep_update(base, override)`
does a recursive merge. `validate(cfg)` raises with a readable list of
missing required top-level keys.
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

REQUIRED_KEYS = ("sample_id", "he_path", "xenium_bundle", "output_root")

VALID_STAGES = (
    "he_preprocess", "register", "warp", "celltype", "viz",
    "nn_celltype_mapping",
)
# `nn_celltype_mapping` is opt-in: it needs a proseg-side purified.h5ad
# and produces the hexenium-input CSV — not every run has that upstream
# artifact yet. Pass it explicitly via --stages when needed.
DEFAULT_STAGES = ("he_preprocess", "register", "warp", "celltype", "viz")


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
    """Raise SystemExit if any required top-level key is missing/null."""
    missing = [k for k in REQUIRED_KEYS if not cfg.get(k)]
    if missing:
        raise SystemExit(f"missing required config keys: {missing}")
