"""Cross-module compatibility shims and helpers.

Apply numpy compatibility shims BEFORE importing HEST/VALIS.
"""
from __future__ import annotations

import numpy as np


def apply_numpy_shims() -> None:
    """Restore numpy attributes that older valis-hest deps still expect."""
    if not hasattr(np, "bool8"):
        np.bool8 = np.bool_
    if not hasattr(np, "float_"):
        np.float_ = np.float64


def kill_jvm_safely() -> None:
    """Best-effort JVM shutdown. Both valis_hest.registration and valis.slide_io
    expose kill_jvm; we try both because installed versions vary."""
    for mod in ("valis_hest.registration", "valis.slide_io"):
        try:
            module = __import__(mod, fromlist=["kill_jvm"])
            if hasattr(module, "kill_jvm"):
                module.kill_jvm()
                return
        except Exception:
            continue
