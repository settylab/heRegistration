"""OME-metadata-aware wrapper around HEST's ``SlideReaderAdapter``.

Upstream ``hest.SlideReaderAdapter.create_metadata`` hard-codes the
H&E slide's ``pixel_physical_size_xyu`` to ``[0.25, 0.25, PIXEL_UNIT]``
regardless of what the source OME-TIFF's metadata says. For any scan
that isn't natively 0.25 µm/px — including the Olympus VS-series 40× BF
we get from the H&E scanner (~0.1369 µm/px) — VALIS then aligns against
a wrong H&E scale, inflating the H&E's implied physical size by
``0.25 / PhysicalSize`` (~1.83× for the VS 40× case) and breaking the
DAPI ↔ H&E alignment.

This wrapper subclasses HEST's adapter, keeps its WSI/pyramid handling
untouched, and only overrides ``create_metadata`` to inject the true
PhysicalSizeX/Y from the source file's OME-XML.

Fallback semantics: if PhysicalSize cannot be read (non-OME TIFF, no
``Image[0].Pixels`` block, unit that isn't microns, or a parse
failure), we log a loud WARNING and leave HEST's default in place.
This preserves backward-compatibility for any pipeline invocation
that was working before this file existed — the failure mode is at
worst identical to today's behavior, plus an actionable log line the
operator can grep for in slurm ``.err`` files.
"""
from __future__ import annotations

import logging
from typing import Optional

import tifffile
from hest.SlideReaderAdapter import SlideReaderAdapter
from valis.slide_io import MICRON_UNIT

try:  # pragma: no cover — env-dep
    from ome_types import from_xml as _ome_from_xml
except ImportError:  # pragma: no cover — env-dep
    _ome_from_xml = None


logger = logging.getLogger(__name__)


_MICRON_UNIT_STRINGS = frozenset({"µm", "micrometer", "um", "micron"})


class OMEAwareSlideReaderAdapter(SlideReaderAdapter):
    """SlideReaderAdapter that reads OME PhysicalSize from the input.

    Every other behavior is inherited from
    :class:`hest.SlideReaderAdapter` unchanged.
    """

    def create_metadata(self):
        meta = super().create_metadata()
        _apply_ome_physical_size(meta, self.src_f)
        return meta


def _apply_ome_physical_size(meta, src_f) -> bool:
    """Mutate ``meta.pixel_physical_size_xyu`` in place from OME PhysicalSize.

    Returns ``True`` if the OME PhysicalSize was applied, ``False`` if
    we fell back to whatever the passed-in ``meta`` already had. In
    the fallback branch, logs a WARNING pointing at the file.
    """
    phys = _read_ome_physical_size(src_f)
    if phys is None:
        logger.warning(
            "[OMEAwareSlideReaderAdapter] could not read OME "
            "PhysicalSizeX/Y from %s; falling back to HEST's "
            "hard-coded [0.25, 0.25, PIXEL_UNIT] default. "
            "Registration accuracy will be degraded if the scan is "
            "not natively 0.25 µm/px.", src_f,
        )
        return False
    phys_x, phys_y = phys
    meta.pixel_physical_size_xyu = [phys_x, phys_y, MICRON_UNIT]
    logger.info(
        "[OMEAwareSlideReaderAdapter] using OME PhysicalSize "
        "%.6f × %.6f µm/px for %s", phys_x, phys_y, src_f,
    )
    return True


def _read_ome_physical_size(src_f) -> Optional[tuple[float, float]]:
    """Return ``(PhysicalSizeX, PhysicalSizeY)`` in µm from an OME-TIFF, or None.

    Uses :func:`ome_types.from_xml` to parse the OME-XML block that
    ``tifffile.TiffFile.ome_metadata`` extracts. Returns ``None`` for
    any of: non-OME input, no ``Image[0].Pixels`` block, missing
    PhysicalSizeX/Y, non-micron unit, or a parse error.
    """
    if _ome_from_xml is None:
        return None
    try:
        with tifffile.TiffFile(str(src_f)) as tf:
            xml = tf.ome_metadata
        if not xml:
            return None
        ome = _ome_from_xml(xml)
        if not ome.images:
            return None
        pixels = ome.images[0].pixels
        px = pixels.physical_size_x
        py = pixels.physical_size_y
        if px is None or py is None:
            return None
        for u in (
            getattr(pixels, "physical_size_x_unit", None),
            getattr(pixels, "physical_size_y_unit", None),
        ):
            if u is None:
                continue
            u_str = getattr(u, "value", str(u)).lower()
            if u_str not in _MICRON_UNIT_STRINGS:
                return None
        return float(px), float(py)
    except Exception as e:
        logger.warning(
            "[OMEAwareSlideReaderAdapter] OME PhysicalSize read failed "
            "for %s: %s", src_f, e,
        )
        return None


__all__ = ["OMEAwareSlideReaderAdapter"]
