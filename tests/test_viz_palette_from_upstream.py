"""Tests for hexenium viz palette inheritance from the upstream
xenium-preprocess color map (`<xenium_run_dir>/summary/
<sample_id>_color_map.json`).

Precedence enforced:
  1. YAML ``viz.classification_palette`` override.
  2. Upstream ``first_type`` map from the xenium-preprocess summary.
  3. Cmap fallback.

Every failure mode of the upstream read must fall through to the
cmap, not crash.
"""
from __future__ import annotations

import json
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest


matplotlib_pyplot = pytest.importorskip("matplotlib.pyplot")


# ---------- _hex_to_rgb ----------------------------------------------


@pytest.mark.parametrize(
    "hex_str,expected",
    [
        ("#FF1493", [255, 20, 147]),   # canonical form
        ("FF1493", [255, 20, 147]),    # no leading '#'
        ("#a9a9a9", [169, 169, 169]),  # lowercase
        ("  #808000  ", [128, 128, 0]),  # whitespace tolerated
    ],
)
def test_hex_to_rgb_valid(hex_str, expected):
    from hexenium.stages.viz import _hex_to_rgb

    assert _hex_to_rgb(hex_str) == expected


@pytest.mark.parametrize(
    "bad",
    ["", "not hex", "#12345", "#GGGGGG", None, 12345, "#12345678"],
)
def test_hex_to_rgb_bad_input_returns_none(bad):
    from hexenium.stages.viz import _hex_to_rgb

    assert _hex_to_rgb(bad) is None


# ---------- _load_upstream_color_map ---------------------------------


def _write_color_map_json(dir_path, sample_id, first_type=None,
                          purification_status=None):
    dir_path.mkdir(parents=True, exist_ok=True)
    payload = {}
    if first_type is not None:
        payload["first_type"] = first_type
    if purification_status is not None:
        payload["purification_status"] = purification_status
    path = dir_path / f"{sample_id}_color_map.json"
    path.write_text(json.dumps(payload))
    return path


def test_load_upstream_color_map_reads_first_type(tmp_path):
    from hexenium.stages.viz import _load_upstream_color_map

    _write_color_map_json(
        tmp_path / "summary",
        "S",
        first_type={"Fibroblast": "#FF1493", "Myeloid": "#E69F00"},
        purification_status={"singlet": "#117733"},
    )
    got = _load_upstream_color_map(tmp_path, "S")
    assert got == {"Fibroblast": [255, 20, 147], "Myeloid": [230, 159, 0]}


def test_load_upstream_color_map_none_xenium_run_dir_returns_empty():
    from hexenium.stages.viz import _load_upstream_color_map

    assert _load_upstream_color_map(None, "S") == {}


def test_load_upstream_color_map_missing_file_returns_empty(tmp_path):
    from hexenium.stages.viz import _load_upstream_color_map

    # No `summary/` under xenium_run_dir at all.
    assert _load_upstream_color_map(tmp_path, "S") == {}


def test_load_upstream_color_map_malformed_json_returns_empty(tmp_path):
    from hexenium.stages.viz import _load_upstream_color_map

    (tmp_path / "summary").mkdir()
    (tmp_path / "summary" / "S_color_map.json").write_text("not valid json {")
    assert _load_upstream_color_map(tmp_path, "S") == {}


def test_load_upstream_color_map_missing_first_type_returns_empty(tmp_path):
    from hexenium.stages.viz import _load_upstream_color_map

    _write_color_map_json(
        tmp_path / "summary", "S",
        purification_status={"singlet": "#117733"},
        # note: no first_type
    )
    assert _load_upstream_color_map(tmp_path, "S") == {}


def test_load_upstream_color_map_drops_malformed_hex_entries(tmp_path):
    """Corrupt individual hex codes must not tank the whole read; the
    good entries survive."""
    from hexenium.stages.viz import _load_upstream_color_map

    _write_color_map_json(
        tmp_path / "summary", "S",
        first_type={
            "GoodCell": "#a9a9a9",
            "BadHexTooShort": "#12345",
            "NotEvenAString": ["not", "a", "string"],
            "AlsoGood": "#808000",
        },
    )
    got = _load_upstream_color_map(tmp_path, "S")
    assert got == {"GoodCell": [169, 169, 169], "AlsoGood": [128, 128, 0]}


# ---------- _build_full_palette precedence --------------------------


def test_build_palette_yaml_override_takes_priority():
    """Manual YAML `viz.classification_palette` wins over the upstream
    JSON per Tracy's spec on comment 5349046505."""
    from hexenium.stages.viz import _build_full_palette

    palette = _build_full_palette(
        labels=["Fibroblast"],
        user_overrides={"Fibroblast": [10, 20, 30]},
        upstream_palette={"Fibroblast": [255, 20, 147]},
    )
    assert palette["Fibroblast"] == [10, 20, 30]


def test_build_palette_upstream_used_when_no_yaml_override():
    """No YAML override → upstream JSON supplies the color."""
    from hexenium.stages.viz import _build_full_palette

    palette = _build_full_palette(
        labels=["Fibroblast"],
        user_overrides=None,
        upstream_palette={"Fibroblast": [255, 20, 147]},
    )
    assert palette["Fibroblast"] == [255, 20, 147]


def test_build_palette_cmap_fallback_when_no_yaml_and_no_upstream():
    """No YAML, no upstream → cmap fallback assigns a color from
    `tab20` (specifically the first cmap slot for the first label)."""
    from hexenium.stages.viz import _build_full_palette

    palette = _build_full_palette(labels=["Fibroblast"])
    color = palette["Fibroblast"]
    assert isinstance(color, list) and len(color) == 3
    assert all(0 <= c <= 255 for c in color)


def test_build_palette_mixed_sources_per_label():
    """A label in YAML wins; a label only in upstream inherits from
    upstream; a label in neither draws from the cmap. Missing labels
    across sources never crash."""
    from hexenium.stages.viz import _build_full_palette

    palette = _build_full_palette(
        labels=["Yaml", "Upstream", "Neither"],
        user_overrides={"Yaml": [1, 2, 3]},
        upstream_palette={"Upstream": [255, 20, 147]},
    )
    assert palette["Yaml"] == [1, 2, 3]
    assert palette["Upstream"] == [255, 20, 147]
    # Neither in YAML nor upstream — draws from cmap.
    assert palette["Neither"] not in (
        [1, 2, 3], [255, 20, 147]
    )
    assert isinstance(palette["Neither"], list) and len(palette["Neither"]) == 3


def test_build_palette_unclassified_folds_to_unlabeled_default():
    """Legacy ``"Unclassified"`` folds into the ``UNLABELED`` canonical
    slot (Tracy `5360206242`, B+fold). The historic
    ``[180, 180, 180]`` shade is gone; the returned palette carries
    ``unlabeled → #888888`` and no separate ``Unclassified`` key.
    See ``tests/test_viz_unclassified_fold.py`` for the full fold
    matrix."""
    from hexenium.stages.viz import (
        UNLABELED,
        UNLABELED_RGB,
        _build_full_palette,
    )

    palette = _build_full_palette(
        labels=["Unclassified"],
        user_overrides=None,
        upstream_palette=None,
    )
    assert UNLABELED in palette
    assert palette[UNLABELED] == list(UNLABELED_RGB)
    assert "Unclassified" not in palette


def test_build_palette_yaml_unclassified_override_routes_to_unlabeled():
    """A YAML override keyed on the legacy ``"Unclassified"`` label
    now routes to the ``UNLABELED`` canonical slot after the B+fold —
    same precedence rule (YAML wins), same visual outcome. See
    ``tests/test_viz_unclassified_fold.py`` for parity coverage."""
    from hexenium.stages.viz import UNLABELED, _build_full_palette

    palette = _build_full_palette(
        labels=["Unclassified"],
        user_overrides={"Unclassified": [77, 77, 77]},
        upstream_palette={"Unclassified": [0, 0, 0]},
    )
    assert palette[UNLABELED] == [77, 77, 77]
    assert "Unclassified" not in palette


# ---------- End-to-end through run_viz ------------------------------


def _run_viz_with_labels(tmp_path, labels, palette_kwargs):
    """Call `run_viz` with a mocked H&E loader + parquet + matplotlib
    calls, and capture the (RGB, label) pairs handed to
    `ax.plot` for the legend. That's the cleanest read of the final
    palette without re-running matplotlib rendering."""
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
    legend_entries = []

    def _fake_plot(x, y, *_a, color=None, label=None, **_kw):
        # ax.plot([], [], "s", color=rgb, label=label)
        if label is not None:
            legend_entries.append((label, tuple(color)))
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
            render_boundaries="cell",  # ensures ax.plot legend is populated
            **palette_kwargs,
        )
    return dict(legend_entries)


def test_run_viz_end_to_end_upstream_palette_flows_through(tmp_path):
    """Integration: run_viz reads the summary JSON from
    ``xenium_run_dir/summary/`` and colors labels accordingly."""
    _write_color_map_json(
        tmp_path / "xenium_run" / "summary", "S",
        first_type={"Fibroblast": "#FF1493"},
    )
    legend = _run_viz_with_labels(
        tmp_path=tmp_path,
        labels=["Fibroblast"],
        palette_kwargs={"xenium_run_dir": tmp_path / "xenium_run"},
    )
    assert legend["Fibroblast"] == pytest.approx(
        (255 / 255, 20 / 255, 147 / 255), rel=1e-3,
    )


def test_run_viz_end_to_end_yaml_beats_upstream(tmp_path):
    """Integration: an explicit YAML `classification_palette` entry
    trumps the upstream JSON for that label."""
    _write_color_map_json(
        tmp_path / "xenium_run" / "summary", "S",
        first_type={"Fibroblast": "#FF1493"},
    )
    legend = _run_viz_with_labels(
        tmp_path=tmp_path,
        labels=["Fibroblast"],
        palette_kwargs={
            "xenium_run_dir": tmp_path / "xenium_run",
            "classification_palette": {"Fibroblast": [0, 0, 0]},
        },
    )
    assert legend["Fibroblast"] == pytest.approx((0.0, 0.0, 0.0))


def test_run_viz_end_to_end_missing_summary_falls_back_to_cmap(tmp_path):
    """Integration: no summary JSON → cmap fallback, no crash."""
    # No summary/ dir created under xenium_run
    legend = _run_viz_with_labels(
        tmp_path=tmp_path,
        labels=["NovelCell"],
        palette_kwargs={"xenium_run_dir": tmp_path / "xenium_run"},
    )
    assert "NovelCell" in legend
    # cmap-drawn color is a 3-tuple of floats in [0, 1].
    color = legend["NovelCell"]
    assert len(color) == 3 and all(0.0 <= c <= 1.0 for c in color)


def test_run_viz_end_to_end_no_xenium_run_dir_falls_back(tmp_path):
    """Integration: xenium_run_dir=None (standalone mode) → cmap
    fallback, no crash."""
    legend = _run_viz_with_labels(
        tmp_path=tmp_path,
        labels=["Fibroblast"],
        palette_kwargs={"xenium_run_dir": None},
    )
    assert "Fibroblast" in legend
