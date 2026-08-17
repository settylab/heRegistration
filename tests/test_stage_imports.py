"""Smoke test: each stage module imports without ImportError.

Each stage's TOP-LEVEL imports declare its dependency contract. VALIS
/ HEST / cv2 / dask / openslide are imported LAZILY inside the
`run_<stage>` bodies, so this smoke test needs only the numeric-stack
subset (numpy for _internal.compat, plus the celltyping stack for
that module's module-level pandas/geopandas/shapely imports).

Missing science stack? The test is skipped with a clear reason — do
NOT block CI on the full heRegistration env being present.
"""
from __future__ import annotations

import importlib

import pytest


STAGES = (
    "hexenium.stages.he_preprocess",
    "hexenium.stages.registration",
    "hexenium.stages.warp",
    "hexenium.stages.celltyping",
    "hexenium.stages.viz",
    "hexenium.stages.nn_celltype_mapping",
)


@pytest.mark.parametrize("modname", STAGES)
def test_stage_module_imports(modname: str):
    try:
        importlib.import_module(modname)
    except ModuleNotFoundError as e:
        # Missing heRegistration env dep — skip, don't fail.
        pytest.skip(f"missing dependency for {modname}: {e.name}")


def test_stages_export_run_functions():
    """Each stage exposes a `run_<stage>` entry point that pipeline.py
    calls into. This test skips a stage if its module-level imports
    can't be satisfied in the current env."""
    entrypoints = {
        "hexenium.stages.he_preprocess": "run_he_preprocess",
        "hexenium.stages.registration": "run_registration",
        "hexenium.stages.warp": "run_warp",
        "hexenium.stages.celltyping": "run_celltyping",
        "hexenium.stages.viz": "run_viz",
        "hexenium.stages.nn_celltype_mapping": "run_nn_celltype_mapping",
    }
    for modname, fname in entrypoints.items():
        try:
            mod = importlib.import_module(modname)
        except ModuleNotFoundError as e:
            pytest.skip(f"missing dependency for {modname}: {e.name}")
        assert hasattr(mod, fname), f"{modname} is missing entry point {fname}"
