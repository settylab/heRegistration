"""Tests for the legacy-``Unclassified`` → ``unlabeled`` viz-time fold.

Locked in by Tracy on ``settylab/TracyY123-nexus#15`` comment
``5360206242`` (B+fold): legacy celltyped parquets that still carry
``"Unclassified"`` on disk render at the same ``#888888`` as fresh
``"unlabeled"`` runs, with no separate palette entry. New runs
continue to write ``"unlabeled"``; only VIZ-side rendering treats
the legacy label as a synonym.
"""
from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest


matplotlib_pyplot = pytest.importorskip("matplotlib.pyplot")


# ---------- _canonicalize_classification ----------------------------


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("Unclassified", "unlabeled"),     # the fold
        ("unlabeled", "unlabeled"),         # canonical stays canonical
        ("T", "T"),                         # regular label passthrough
        ("Fibroblast", "Fibroblast"),       # regular label passthrough
        (None, "unlabeled"),                # None → UNLABELED
        (float("nan"), "unlabeled"),        # NaN (pandas' None-in-object)
                                            # → UNLABELED
    ],
)
def test_canonicalize_classification(raw, expected):
    from hexenium.stages.viz import _canonicalize_classification

    assert _canonicalize_classification(raw) == expected


# ---------- _build_full_palette folds legacy Unclassified -----------


def test_build_palette_legacy_unclassified_folds_to_unlabeled():
    """A raw label list containing 'Unclassified' produces a palette
    keyed by 'unlabeled' (not 'Unclassified'), and the value is the
    canonical #888888. The legacy [180,180,180] shade is gone."""
    from hexenium.stages.viz import (
        UNLABELED,
        UNLABELED_RGB,
        _build_full_palette,
    )

    palette = _build_full_palette(labels=["Unclassified"])

    assert UNLABELED in palette
    assert palette[UNLABELED] == list(UNLABELED_RGB)
    # No separate legacy entry.
    assert "Unclassified" not in palette
    # And the value is NOT the old [180,180,180] shade.
    assert palette[UNLABELED] != [180, 180, 180]


def test_build_palette_mixed_unclassified_and_unlabeled_collapse():
    """A dataset with BOTH legacy 'Unclassified' AND fresh 'unlabeled'
    rows collapses to ONE palette entry (canonical), both rendering
    identically."""
    from hexenium.stages.viz import (
        UNLABELED,
        UNLABELED_RGB,
        _build_full_palette,
    )

    palette = _build_full_palette(labels=["Unclassified", "unlabeled"])

    assert list(palette.keys()) == [UNLABELED]
    assert palette[UNLABELED] == list(UNLABELED_RGB)


def test_build_palette_yaml_override_wins_over_fold():
    """A YAML override for 'unlabeled' still wins after the fold —
    legacy 'Unclassified' rows pick up the OVERRIDDEN color, not the
    default #888888. Users can still pin the sentinel's colour."""
    from hexenium.stages.viz import UNLABELED, _build_full_palette

    palette = _build_full_palette(
        labels=["Unclassified"],
        user_overrides={UNLABELED: [10, 20, 30]},
    )
    assert palette[UNLABELED] == [10, 20, 30]


def test_build_palette_yaml_unclassified_override_also_folds():
    """A YAML override keyed on the LEGACY 'Unclassified' label also
    lands on the canonical entry — the fold treats the two names as
    synonyms at every stage, including user overrides."""
    from hexenium.stages.viz import UNLABELED, _build_full_palette

    # User's YAML still says "Unclassified: [10, 20, 30]" — the fold
    # should route that override to the canonical UNLABELED slot.
    palette = _build_full_palette(
        labels=["Unclassified"],
        user_overrides={"Unclassified": [10, 20, 30]},
    )
    # The palette key is canonical.
    assert UNLABELED in palette
    # The value was routed from the user's YAML override.
    assert palette[UNLABELED] == [10, 20, 30]
    assert "Unclassified" not in palette


# ---------- End-to-end through run_viz ------------------------------


def _legend_from_run_viz(labels, tmp_path):
    """Drive run_viz with a mocked H&E loader + parquet + matplotlib
    plot; capture the (label, color) legend the palette produces."""
    import geopandas as gpd
    from shapely.geometry import Polygon
    from hexenium.stages import viz as viz_mod

    poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    gdf = gpd.GeoDataFrame(
        {"boundary_type": ["cell"] * len(labels),
         "classification": labels},
        geometry=[poly] * len(labels),
    )
    ax = MagicMock()
    legend = []

    def _fake_plot(x, y, *_a, color=None, label=None, **_kw):
        if label is not None:
            legend.append((label, tuple(color)))
        return [SimpleNamespace()]

    ax.plot.side_effect = _fake_plot
    fig = MagicMock()

    celltyped_dir = tmp_path / "celltyped"
    celltyped_dir.mkdir()
    (celltyped_dir / "S_celltyped_wholeslide.parquet").touch()
    out_dir = tmp_path / "viz_out"

    with (
        patch.object(viz_mod, "_load_he_thumbnail",
                     return_value=(np.zeros((4, 4, 3), np.uint8), 1.0)),
        patch("geopandas.read_parquet", return_value=gdf),
        patch("matplotlib.pyplot.subplots", return_value=(fig, ax)),
        patch("matplotlib.collections.PolyCollection",
              side_effect=lambda *_a, **_kw: SimpleNamespace()),
    ):
        viz_mod.run_viz(
            sample_id="S",
            he_path=tmp_path / "he.tif",
            warp_dir=tmp_path / "warped",
            celltyped_dir=celltyped_dir,
            out_dir=out_dir,
            render_boundaries="cell",
        )
    return dict(legend)


def test_run_viz_legacy_unclassified_renders_as_888888(tmp_path):
    """End-to-end regression: a parquet with 'Unclassified' rows
    renders those cells at #888888 (canonical unlabeled color) —
    NOT the historic #b4b4b4 shade. Confirms the fold reaches all
    the way from disk read to matplotlib color arg."""
    legend = _legend_from_run_viz(["Unclassified"], tmp_path)
    # Legend key: canonicalised label.
    assert "unlabeled" in legend
    assert "Unclassified" not in legend
    # Legend color: #888888 in matplotlib's 0-1 float space.
    assert legend["unlabeled"] == pytest.approx(
        (0x88 / 255, 0x88 / 255, 0x88 / 255), rel=1e-3,
    )
    # Explicitly NOT the historic shade.
    assert legend["unlabeled"] != pytest.approx(
        (180 / 255, 180 / 255, 180 / 255), rel=1e-3,
    )


def test_run_viz_mixed_dataset_unifies_greys(tmp_path):
    """A parquet with BOTH legacy 'Unclassified' AND fresh 'unlabeled'
    rows renders them under a single legend entry with the same
    color — no visual distinction between old and new."""
    legend = _legend_from_run_viz(
        ["Unclassified", "unlabeled"], tmp_path,
    )
    # Both fold into one legend entry.
    assert list(legend.keys()) == ["unlabeled"]
    assert legend["unlabeled"] == pytest.approx(
        (0x88 / 255, 0x88 / 255, 0x88 / 255), rel=1e-3,
    )


def test_run_viz_none_classification_falls_into_unlabeled(tmp_path):
    """A row with a null classification (None) also renders as
    UNLABELED — matches the fallback default at the render site."""
    import geopandas as gpd
    from shapely.geometry import Polygon
    from hexenium.stages import viz as viz_mod

    poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    # Include a plain label so the palette has more than one entry
    # (checks the None → unlabeled fold picks up the right color).
    gdf = gpd.GeoDataFrame(
        {"boundary_type": ["cell", "cell"],
         "classification": [None, "T"]},
        geometry=[poly, poly],
    )
    ax = MagicMock()
    legend = []

    def _fake_plot(x, y, *_a, color=None, label=None, **_kw):
        if label is not None:
            legend.append((label, tuple(color)))
        return [SimpleNamespace()]

    ax.plot.side_effect = _fake_plot
    fig = MagicMock()

    celltyped_dir = tmp_path / "celltyped"
    celltyped_dir.mkdir()
    (celltyped_dir / "S_celltyped_wholeslide.parquet").touch()
    with (
        patch.object(viz_mod, "_load_he_thumbnail",
                     return_value=(np.zeros((4, 4, 3), np.uint8), 1.0)),
        patch("geopandas.read_parquet", return_value=gdf),
        patch("matplotlib.pyplot.subplots", return_value=(fig, ax)),
        patch("matplotlib.collections.PolyCollection",
              side_effect=lambda *_a, **_kw: SimpleNamespace()),
    ):
        viz_mod.run_viz(
            sample_id="S",
            he_path=tmp_path / "he.tif",
            warp_dir=tmp_path / "warped",
            celltyped_dir=celltyped_dir,
            out_dir=tmp_path / "viz_out",
            render_boundaries="cell",
        )

    legend_d = dict(legend)
    assert "unlabeled" in legend_d
    assert legend_d["unlabeled"] == pytest.approx(
        (0x88 / 255, 0x88 / 255, 0x88 / 255), rel=1e-3,
    )
