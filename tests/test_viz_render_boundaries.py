"""Thin viz test: `render_boundaries` gates which boundary_type branch runs.

Mocks the H&E loader and geopandas read; records the facecolors of each
`PolyCollection` build to tell cell (filled RGBA) vs nucleus (`"none"`)
apart. Confirms `nucleus` (2026-07-09 default) skips cell polygons.
"""
from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest


def _kinds_drawn(render_boundaries, tmp_path):
    import geopandas as gpd
    from shapely.geometry import Polygon
    from hexenium.stages import viz as viz_mod

    poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    gdf = gpd.GeoDataFrame(
        {"boundary_type": ["cell", "nucleus"], "classification": ["T", "T"]},
        geometry=[poly, poly],
    )
    seen = []

    def _fake_pc(verts, **kw):
        seen.append("nucleus" if kw.get("facecolors") == "none" else "cell")
        return SimpleNamespace()

    celltyped_dir = tmp_path / "celltyped"
    celltyped_dir.mkdir()
    (celltyped_dir / "S_celltyped_wholeslide.parquet").touch()
    out_dir = tmp_path / "viz_out"
    with (
        patch.object(viz_mod, "_load_he_thumbnail",
                     return_value=(np.zeros((4, 4, 3), np.uint8), 1.0)),
        patch("geopandas.read_parquet", return_value=gdf),
        patch("matplotlib.pyplot.subplots", return_value=(MagicMock(), MagicMock())),
        patch("matplotlib.collections.PolyCollection", side_effect=_fake_pc),
    ):
        viz_mod.run_viz(
            sample_id="S",
            he_path=tmp_path / "he.tif",
            warp_dir=tmp_path / "warped",
            celltyped_dir=celltyped_dir,
            out_dir=out_dir,
            render_boundaries=render_boundaries,
        )
    return set(seen)


@pytest.mark.parametrize(
    "render_boundaries,expected",
    [("nucleus", {"nucleus"}), ("cell", {"cell"}), ("both", {"cell", "nucleus"})],
)
def test_render_boundaries_gates_draw_loop(render_boundaries, expected, tmp_path):
    assert _kinds_drawn(render_boundaries, tmp_path) == expected


def test_render_boundaries_rejects_invalid(tmp_path):
    from hexenium.stages import viz as viz_mod
    with pytest.raises(ValueError):
        viz_mod.run_viz(
            sample_id="S", he_path=tmp_path / "h",
            warp_dir=tmp_path / "w", celltyped_dir=tmp_path / "c",
            out_dir=tmp_path / "o", render_boundaries="everything",
        )
