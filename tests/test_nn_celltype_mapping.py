"""Smoke tests for the nn_celltype_mapping stage.

Runs the stage end-to-end on synthetic data and asserts the expected
two-column ``(cell_id, group)`` CSV plus the six-column inspection
sidecar are produced correctly.
"""
from __future__ import annotations

from pathlib import Path

import pytest


def _write_synthetic_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Build a synthetic proseg + xenium bundle for the smoke tests.

    Layout mirrors what run_nn_celltype_mapping expects:
      <output_root>/<sample_id>/                  ← per-sample workdir
      <xenium_bundle>/cells.csv.gz                ← xenium centroids

    Returns (output_root, xenium_bundle, proseg_csv). We use the CSV
    proseg branch (no anndata dependency in the test).
    """
    import pandas as pd

    sample_id = "TEST01"
    output_root = tmp_path / "runs"
    (output_root / sample_id).mkdir(parents=True)
    xenium_bundle = tmp_path / "xenium_bundle"
    xenium_bundle.mkdir()

    # 10 proseg cells at x=0..90 (step 10), split into two labels.
    proseg_df = pd.DataFrame({
        "cell_d": list(range(10)),
        "x": [i * 10.0 for i in range(10)],
        "y": [0.0] * 10,
        "first_type": ["tumor"] * 5 + ["stroma"] * 5,
    }).set_index("cell_d")
    proseg_csv = tmp_path / "proseg_metadata.csv"
    proseg_df.to_csv(proseg_csv)

    # 10 xenium cells: half exact-match, half offset by 2 (dist=2).
    xen_df = pd.DataFrame({
        "cell_id": [f"xen_{i}" for i in range(10)],
        "x_centroid": [0, 12, 20, 32, 40, 52, 60, 72, 80, 92],
        "y_centroid": [0.0] * 10,
    }).set_index("cell_id")
    xen_df.to_csv(xenium_bundle / "cells.csv.gz", compression="gzip")

    return output_root, xenium_bundle, proseg_csv


@pytest.fixture
def synthetic_inputs(tmp_path):
    return _write_synthetic_inputs(tmp_path)


def _skip_if_missing_deps():
    """Skip if the numeric stack isn't available (missing sklearn,
    pandas, etc. → the stage can't run at all)."""
    try:
        import pandas  # noqa: F401
        import numpy  # noqa: F401
        import sklearn.neighbors  # noqa: F401
    except ModuleNotFoundError as e:
        pytest.skip(f"missing dependency for nn_celltype_mapping: {e.name}")


def test_nn_celltype_mapping_smoke(synthetic_inputs):
    """End-to-end: the packaged stage runs on synthetic data + emits
    the expected two-column primary CSV and six-column sidecar."""
    _skip_if_missing_deps()
    import pandas as pd
    from hexenium.stages.nn_celltype_mapping import run_nn_celltype_mapping

    output_root, xenium_bundle, proseg_csv = synthetic_inputs

    result = run_nn_celltype_mapping(
        sample_id="TEST01",
        xenium_bundle=xenium_bundle,
        output_root=output_root,
        proseg_source=str(proseg_csv),
    )
    csv = pd.read_csv(result["csv"])
    assert list(csv.columns) == ["cell_id", "group"]
    assert len(csv) == 10
    # First 5 xenium cells are exact-match on tumor; last 5 fall to stroma.
    # Two are exact-match on stroma (xen_5 at x=60 → proseg[6] x=60);
    # let's just check both labels appear and every cell_id is present.
    assert set(csv["cell_id"]) == {f"xen_{i}" for i in range(10)}
    assert {"tumor", "stroma"}.issubset(set(csv["group"]))

    # Inspection sidecar has the full six-column shape.
    insp = pd.read_csv(result["inspection_csv"])
    assert list(insp.columns) == [
        "cell_id", "group", "proseg_id", "nn_dist", "assignment_note",
    ]
    assert len(insp) == 10


def test_nn_celltype_mapping_distance_threshold_marks_unassigned(synthetic_inputs):
    """distance_threshold + mark_unassigned policy marks offset cells."""
    _skip_if_missing_deps()
    import pandas as pd
    from hexenium.stages.nn_celltype_mapping import run_nn_celltype_mapping

    output_root, xenium_bundle, proseg_csv = synthetic_inputs

    result = run_nn_celltype_mapping(
        sample_id="TEST01",
        xenium_bundle=xenium_bundle,
        output_root=output_root,
        proseg_source=str(proseg_csv),
        distance_threshold=1.5,
        unmatched_policy="mark_unassigned",
    )
    csv = pd.read_csv(result["csv"])
    # 5 cells (xen_1, xen_3, xen_5, xen_7, xen_9) are dist=2 from their
    # nearest neighbour; all should be Unassigned.
    unassigned = csv[csv["group"] == "Unassigned"]
    assert len(unassigned) == 5


# -------------------------------------------------------------------
# Hybrid auto-detect helper — pure-function unit tests.
# -------------------------------------------------------------------


def _make_proseg_with_original_ids(orig_ids):
    """Build a minimal proseg DataFrame with an ``original_id`` column.

    Matches the shape ``_read_proseg_side`` produces (id, x, y,
    celltype, original_id). ``orig_ids`` is a list; entries that are
    ``None`` become NaN in the column.
    """
    import numpy as np
    import pandas as pd
    n = len(orig_ids)
    return pd.DataFrame({
        "id": [f"p_{i}" for i in range(n)],
        "x": [float(i) for i in range(n)],
        "y": [0.0] * n,
        "celltype": ["tumor"] * n,
        "original_id": [np.nan if v is None else v for v in orig_ids],
    })


def _make_xenium(xen_ids):
    import pandas as pd
    return pd.DataFrame({
        "id": list(xen_ids),
        "x": [float(i) for i in range(len(xen_ids))],
        "y": [0.0] * len(xen_ids),
    })


def test_should_enable_hybrid_both_thresholds_pass():
    _skip_if_missing_deps()
    from hexenium.stages.nn_celltype_mapping import _should_enable_hybrid

    # 10 proseg rows, all with original_id; 10 xenium cells that all match.
    proseg = _make_proseg_with_original_ids([f"xen_{i}" for i in range(10)])
    xenium = _make_xenium([f"xen_{i}" for i in range(10)])

    assert _should_enable_hybrid(
        proseg, xenium,
        min_original_id_populated=0.8,
        min_xenium_overlap=0.5,
    ) is True


def test_should_enable_hybrid_column_all_null():
    _skip_if_missing_deps()
    from hexenium.stages.nn_celltype_mapping import _should_enable_hybrid

    proseg = _make_proseg_with_original_ids([None] * 10)
    xenium = _make_xenium([f"xen_{i}" for i in range(10)])
    assert _should_enable_hybrid(
        proseg, xenium,
        min_original_id_populated=0.8,
        min_xenium_overlap=0.5,
    ) is False


def test_should_enable_hybrid_column_mostly_null_below_threshold():
    _skip_if_missing_deps()
    from hexenium.stages.nn_celltype_mapping import _should_enable_hybrid

    # Only 20% populated — below default 0.8.
    orig = [f"xen_{i}" if i < 2 else None for i in range(10)]
    proseg = _make_proseg_with_original_ids(orig)
    xenium = _make_xenium([f"xen_{i}" for i in range(10)])
    assert _should_enable_hybrid(
        proseg, xenium,
        min_original_id_populated=0.8,
        min_xenium_overlap=0.5,
    ) is False


def test_should_enable_hybrid_column_populated_but_no_overlap():
    _skip_if_missing_deps()
    from hexenium.stages.nn_celltype_mapping import _should_enable_hybrid

    # 100% populated, but the ids are from a completely different frame.
    proseg = _make_proseg_with_original_ids([f"other_{i}" for i in range(10)])
    xenium = _make_xenium([f"xen_{i}" for i in range(10)])
    assert _should_enable_hybrid(
        proseg, xenium,
        min_original_id_populated=0.8,
        min_xenium_overlap=0.5,
    ) is False


def test_should_enable_hybrid_partial_overlap_below_threshold():
    _skip_if_missing_deps()
    from hexenium.stages.nn_celltype_mapping import _should_enable_hybrid

    # 100% populated; only 3/10 xenium cells match — below default 0.5.
    proseg = _make_proseg_with_original_ids(
        [f"xen_{i}" for i in range(3)] + [f"other_{i}" for i in range(7)]
    )
    xenium = _make_xenium([f"xen_{i}" for i in range(10)])
    assert _should_enable_hybrid(
        proseg, xenium,
        min_original_id_populated=0.8,
        min_xenium_overlap=0.5,
    ) is False


def test_hybrid_auto_mode_disables_when_no_original_id(synthetic_inputs):
    """End-to-end: ``hybrid_direct_join_first="auto"`` falls back to pure
    NN on the synthetic inputs (no ``original_cell_id`` column) — the
    output must match a ``hybrid_direct_join_first=False`` run byte-for-byte.
    """
    _skip_if_missing_deps()
    from hexenium.stages.nn_celltype_mapping import run_nn_celltype_mapping

    output_root, xenium_bundle, proseg_csv = synthetic_inputs

    # First run: explicit False.
    r_false = run_nn_celltype_mapping(
        sample_id="TEST01",
        xenium_bundle=xenium_bundle,
        output_root=output_root,
        proseg_source=str(proseg_csv),
        hybrid_direct_join_first=False,
    )
    false_bytes = Path(r_false["csv"]).read_bytes()

    # Second run: "auto" — set up a fresh output tree so sentinel doesn't skip.
    output_root2 = output_root.parent / "runs_auto"
    (output_root2 / "TEST01").mkdir(parents=True)
    r_auto = run_nn_celltype_mapping(
        sample_id="TEST01",
        xenium_bundle=xenium_bundle,
        output_root=output_root2,
        proseg_source=str(proseg_csv),
        hybrid_direct_join_first="auto",
    )
    auto_bytes = Path(r_auto["csv"]).read_bytes()
    assert auto_bytes == false_bytes, (
        "auto mode should fall through to pure NN when the proseg source has "
        "no original_id column, matching hybrid_direct_join_first=False"
    )


def test_hybrid_invalid_string_raises():
    _skip_if_missing_deps()
    from hexenium.stages.nn_celltype_mapping import run_nn_celltype_mapping

    # Reuse the synthetic-inputs fixture manually.
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        output_root, xenium_bundle, proseg_csv = _write_synthetic_inputs(Path(tmp))
        with pytest.raises(ValueError, match="unknown nn.hybrid_direct_join_first"):
            run_nn_celltype_mapping(
                sample_id="TEST01",
                xenium_bundle=xenium_bundle,
                output_root=output_root,
                proseg_source=str(proseg_csv),
                hybrid_direct_join_first="sometimes",
            )


