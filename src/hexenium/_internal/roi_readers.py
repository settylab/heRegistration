"""ROI polygon readers by source-tool convention.

Each reader returns a Shapely :class:`~shapely.geometry.Polygon` in the
frame declared by the caller. The reader itself is frame-agnostic — it
just decodes the file's vertex list. The caller (typically
:mod:`hexenium.tools.roi_subset`) is responsible for asserting that the
declared frame in the config matches what the source tool actually
exports.

**Known conventions** (empirically verified, see the ``#29`` design
thread):

* QuPath CSV — three ``#``-prefixed header lines, then
  ``Selection,X,Y,Class,Color`` rows. Coordinates in **microns**.
* QuPath GeoJSON — vertex coordinates in **DAPI pixels**
  (0.2125 µm/px × the CSV values, vertex-by-vertex verified).
* Xenium Explorer GeoJSON — coordinate frame **not yet empirically
  verified**; the 10x ``subset2zarr`` tutorial's use of
  ``target_coordinate_system="global"`` implies **microns**. Do NOT
  hard-code a default until an actual XE export is inspected.

The reader signatures are uniform so the dispatcher (:func:`read_roi`)
can route by ``source_tool`` in the config.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterator, Literal, Union

import pandas as pd
from shapely.geometry import MultiPolygon, Polygon
from shapely.geometry.base import BaseGeometry


SourceTool = Literal["qupath_csv", "qupath_geojson", "xenium_explorer"]


def _polygon_from_coords(coords: list) -> Polygon:
    """Build a Shapely Polygon from a list of ``(x, y)`` vertices."""
    if len(coords) < 3:
        raise ValueError(f"polygon needs >= 3 vertices, got {len(coords)}")
    return Polygon(coords)


def read_qupath_csv(
    path: Union[str, Path],
    selection: str | None = None,
) -> Polygon:
    """Read a QuPath-style annotation CSV and return one Shapely polygon.

    The CSV has three ``#``-prefixed comment lines followed by
    ``Selection,X,Y,Class,Color`` rows; each ``Selection`` group is one
    polygon's ordered vertices. When ``selection`` is ``None``, the
    first selection in the file wins; when set, that name is picked.

    Coordinates are in **microns**.
    """
    df = pd.read_csv(path, comment="#")
    required = {"Selection", "X", "Y"}
    if not required <= set(df.columns):
        raise ValueError(
            f"QuPath CSV {path} missing required columns "
            f"{sorted(required - set(df.columns))}"
        )
    if selection is None:
        selection = str(df["Selection"].iloc[0])
    block = df[df["Selection"].astype(str) == str(selection)]
    if len(block) == 0:
        raise ValueError(
            f"QuPath CSV {path} has no Selection={selection!r} "
            f"(available: {sorted(df['Selection'].astype(str).unique())})"
        )
    verts = list(zip(block["X"].tolist(), block["Y"].tolist()))
    return _polygon_from_coords(verts)


def _iter_geojson_features(path: Union[str, Path]) -> Iterator[dict]:
    """Yield features from a GeoJSON Feature or FeatureCollection."""
    with Path(path).open() as f:
        data = json.load(f)
    if data.get("type") == "FeatureCollection":
        yield from data.get("features", [])
    elif data.get("type") == "Feature":
        yield data
    else:
        raise ValueError(
            f"GeoJSON {path}: expected FeatureCollection or Feature, "
            f"got {data.get('type')!r}"
        )


def _feature_name(feature: dict, default: str) -> str:
    """Extract the annotation name from a GeoJSON feature."""
    props = feature.get("properties") or {}
    name = props.get("name") or props.get("classification", {}).get("name")
    return str(name) if name else default


def _feature_to_polygon(feature: dict) -> BaseGeometry:
    """Materialize a GeoJSON feature's geometry as a Shapely object."""
    geom = feature.get("geometry") or {}
    gtype = geom.get("type")
    coords = geom.get("coordinates")
    if gtype == "Polygon":
        return _polygon_from_coords([tuple(v) for v in coords[0]])
    if gtype == "MultiPolygon":
        parts = [_polygon_from_coords([tuple(v) for v in poly[0]]) for poly in coords]
        return MultiPolygon(parts)
    raise ValueError(f"Unsupported GeoJSON geometry type: {gtype!r}")


def read_geojson(
    path: Union[str, Path],
    selection: str | None = None,
) -> BaseGeometry:
    """Read a GeoJSON polygon (or MultiPolygon).

    Works for both QuPath and Xenium Explorer exports; the FRAME is
    caller-declared (via config), not detected here. When ``selection``
    is ``None``, the first feature wins; otherwise features are matched
    by ``properties.name``.

    QuPath GeoJSON: vertices in DAPI pixels. XE GeoJSON: vertices likely
    in microns (unverified — see module docstring).
    """
    for feature in _iter_geojson_features(path):
        name = _feature_name(feature, default="unnamed")
        if selection is None or name == selection:
            return _feature_to_polygon(feature)
    raise ValueError(
        f"GeoJSON {path} has no feature named {selection!r}"
    )


def iter_geojson_features(path: Union[str, Path]) -> Iterator[tuple[str, BaseGeometry]]:
    """Yield ``(name, geometry)`` for every feature in a GeoJSON.

    Used by callers that want to expand a multi-feature file into
    multiple ROIs without needing to know each name upfront.
    """
    for i, feature in enumerate(_iter_geojson_features(path)):
        name = _feature_name(feature, default=f"feature_{i}")
        yield name, _feature_to_polygon(feature)


def iter_qupath_csv_selections(
    path: Union[str, Path],
) -> Iterator[tuple[str, BaseGeometry]]:
    """Yield ``(selection_name, polygon)`` for every Selection in a QuPath CSV."""
    df = pd.read_csv(path, comment="#")
    if not {"Selection", "X", "Y"} <= set(df.columns):
        raise ValueError(f"QuPath CSV {path} missing required columns")
    for name, block in df.groupby("Selection", sort=False):
        verts = list(zip(block["X"].tolist(), block["Y"].tolist()))
        if len(verts) >= 3:
            yield str(name), _polygon_from_coords(verts)


READERS = {
    "qupath_csv": read_qupath_csv,
    "qupath_geojson": read_geojson,
    "xenium_explorer": read_geojson,
}


def read_roi(
    path: Union[str, Path],
    source_tool: SourceTool,
    selection: str | None = None,
) -> BaseGeometry:
    """Dispatch to the reader for ``source_tool`` and return the polygon."""
    reader = READERS.get(source_tool)
    if reader is None:
        raise ValueError(
            f"Unknown source_tool={source_tool!r}; "
            f"expected one of {sorted(READERS)}"
        )
    return reader(path, selection=selection)
