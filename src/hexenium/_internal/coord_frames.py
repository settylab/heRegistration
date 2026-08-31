"""Scalar coordinate-frame conversions for ROI subsetting.

The Xenium pipeline touches four coordinate frames:

* ``microns``   — physical unit; the frame of ``cell_boundaries.parquet``,
  ``nucleus_boundaries.parquet``, ``cells.parquet`` (``x_centroid``,
  ``y_centroid``), and Xenium Explorer GeoJSON exports.
* ``dapi_px``   — pixel index on the DAPI morphology OME-TIFF. Scalar
  factor ``pixel_size_morph`` µm/px between this and microns; default
  ``0.2125`` for Xenium chemistry v2.
* ``he_px``     — pixel index on the H&E OME-TIFF. Non-linearly related
  to the other two via the VALIS registrar. **Not handled here** — A1
  operates entirely in the DAPI / microns half of the diagram; A2 will
  add the H&E half via ``_internal.hest_warp_shim``.
* ``he_um``     — mostly conceptual; H&E `openslide.mpp-{x,y}` gives the
  µm/px factor if a caller ever needs to bridge microns↔H&E-px linearly.
  Not used in Phase 2a.

Only the µm ↔ DAPI-px direction is exposed here; A2 imports the H&E
side from ``hest_warp_shim`` (which owns the registrar-loading code).

Every function is pure and stateless — pass ``pixel_size_morph`` in.
Callers resolve it once from ``experiment.xenium`` via
``read_pixel_size_morph`` and thread it through.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Literal, Union

from shapely.affinity import scale as shapely_scale
from shapely.geometry.base import BaseGeometry


DEFAULT_PIXEL_SIZE_MORPH_UM_PER_PX = 0.2125

Frame = Literal["microns", "dapi_px", "he_px"]
"""Coordinate frames a ROI polygon may be expressed in.

A1 accepts ``microns`` and ``dapi_px``. ``he_px`` is reserved for B1 /
A2 / B2 and will raise from :func:`convert_polygon` in A1 code paths.
"""


def read_pixel_size_morph(experiment_xenium_path: Union[str, Path]) -> float:
    """Read ``pixel_size`` (µm/px) from an ``experiment.xenium`` file.

    Falls back to :data:`DEFAULT_PIXEL_SIZE_MORPH_UM_PER_PX` when the
    file is absent, unreadable, or the ``pixel_size`` field is missing
    or non-positive. The fallback matches Xenium chemistry v2's canonical
    value and is what ``spatialdata_io`` and ``hest`` use themselves.
    """
    path = Path(experiment_xenium_path)
    if not path.exists():
        return DEFAULT_PIXEL_SIZE_MORPH_UM_PER_PX
    try:
        with path.open() as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError):
        return DEFAULT_PIXEL_SIZE_MORPH_UM_PER_PX
    value = data.get("pixel_size")
    if not isinstance(value, (int, float)) or value <= 0:
        return DEFAULT_PIXEL_SIZE_MORPH_UM_PER_PX
    return float(value)


def um_to_dapi_px(polygon: BaseGeometry, pixel_size_morph: float) -> BaseGeometry:
    """Scale a Shapely geometry from microns to DAPI pixel space."""
    factor = 1.0 / pixel_size_morph
    return shapely_scale(polygon, xfact=factor, yfact=factor, origin=(0, 0))


def dapi_px_to_um(polygon: BaseGeometry, pixel_size_morph: float) -> BaseGeometry:
    """Scale a Shapely geometry from DAPI pixel space to microns."""
    return shapely_scale(polygon, xfact=pixel_size_morph, yfact=pixel_size_morph,
                         origin=(0, 0))


def convert_polygon(
    polygon: BaseGeometry,
    src_frame: Frame,
    dst_frame: Frame,
    pixel_size_morph: float,
) -> BaseGeometry:
    """Convert a polygon between the two frames A1 touches.

    ``he_px`` is intentionally rejected here — A1 has no registrar; the
    conversion belongs on the A2/B2 code path once ``hest_warp_shim``
    grows an inverse warp. Callers that need ``he_px`` in A1 should raise
    a configuration error, not silently pretend microns.
    """
    if src_frame == dst_frame:
        return polygon
    if src_frame == "he_px" or dst_frame == "he_px":
        raise ValueError(
            f"A1 does not convert to/from 'he_px' (got {src_frame!r} -> {dst_frame!r}). "
            "H&E-space conversion needs the VALIS registrar; use A2 / B2."
        )
    if src_frame == "microns" and dst_frame == "dapi_px":
        return um_to_dapi_px(polygon, pixel_size_morph)
    if src_frame == "dapi_px" and dst_frame == "microns":
        return dapi_px_to_um(polygon, pixel_size_morph)
    raise ValueError(f"Unhandled conversion: {src_frame!r} -> {dst_frame!r}")
