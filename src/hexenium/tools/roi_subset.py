"""2×2 ROI subsetting entry points.

The 2×2 catalogue (see the ``#29`` design thread):

+----+--------------------------+----------------------------------+----------------------+
|    | ROI drawn on             | Applied to                       | Transform            |
+====+==========================+==================================+======================+
| A1 | original DAPI (µm or px) | original Xenium bundle           | identity (scalar)    |
+----+--------------------------+----------------------------------+----------------------+
| A2 | original DAPI (µm or px) | H&E-registered outputs           | forward VALIS        |
+----+--------------------------+----------------------------------+----------------------+
| B1 | H&E (px)                 | warped H&E-space parquets        | identity             |
+----+--------------------------+----------------------------------+----------------------+
| B2 | H&E (px)                 | original Xenium bundle           | inverse VALIS        |
+----+--------------------------+----------------------------------+----------------------+

This module owns **A1** end-to-end. B1 and A2 live alongside once their
diffs land; B2 is Phase 2b.

Each run consumes a :class:`SubsetConfig` and materialises one
directory per ROI under ``config.output_dir/rois/<slug(name)>/`` with
a ``manifest.yaml`` that records crop origins, source files, and the
geometry frame the emitted GeoDataFrames sit in.
"""
from __future__ import annotations

import json
import logging
import re
import shutil
from dataclasses import asdict
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import tifffile
import yaml
import zarr
from shapely.geometry import Polygon
from shapely.geometry.base import BaseGeometry

from hexenium._internal.coord_frames import (
    DEFAULT_PIXEL_SIZE_MORPH_UM_PER_PX,
    convert_polygon,
    read_pixel_size_morph,
)
from hexenium._internal.roi_config import RoiSpec, SubsetConfig
from hexenium._internal.roi_readers import read_roi


logger = logging.getLogger(__name__)


_SLUG_RE = re.compile(r"[^a-z0-9]+")


def _slug(name: str) -> str:
    """Slug a ROI name for filesystem use. Empty result -> ``roi``."""
    s = _SLUG_RE.sub("_", name.strip().lower()).strip("_")
    return s or "roi"


# ---------------------------------------------------------------------------
# ROI polygon loading + frame conversion
# ---------------------------------------------------------------------------

def _load_roi_polygon_um(
    spec: RoiSpec,
    pixel_size_morph: float,
) -> BaseGeometry:
    """Read a ROI polygon and return it in **microns**.

    A1 operates in microns end-to-end because the boundary parquets and
    the ``cells.parquet`` centroids are natively in microns. DAPI-pixel
    conversion happens only at DAPI-crop time.
    """
    polygon = read_roi(spec.path, source_tool=spec.source_tool, selection=spec.selection)
    return convert_polygon(polygon, spec.frame, "microns", pixel_size_morph)


# ---------------------------------------------------------------------------
# Cell / boundary subsetting (all in microns)
# ---------------------------------------------------------------------------

def _load_cells(xenium_bundle: Path) -> pd.DataFrame:
    return pd.read_parquet(xenium_bundle / "cells.parquet")


def _centroid_mask_xy(
    xs, ys, polygon_um: BaseGeometry,
) -> np.ndarray:
    """Point-in-polygon predicate over arbitrary x/y iterables (microns)."""
    return np.asarray(gpd.points_from_xy(xs, ys).within(polygon_um))


def _centroid_mask(cells: pd.DataFrame, polygon_um: BaseGeometry) -> np.ndarray:
    return _centroid_mask_xy(cells["x_centroid"], cells["y_centroid"], polygon_um)


# ---------------------------------------------------------------------------
# AnnData subsetting (Xenium/Proseg h5ads — all in microns)
# ---------------------------------------------------------------------------

# Canonical mapping from a config field name to the on-disk stem used
# for A1's per-ROI + annotated outputs. Ordered — the annotation step
# iterates in this order for deterministic multi-ROI comma joins.
ANNDATA_INPUTS: tuple[tuple[str, str], ...] = (
    ("xenium_ranger_h5ad", "xenium_ranger"),
    ("proseg_purified_h5ad", "proseg_purified"),
    ("proseg_raw_h5ad", "proseg_raw"),
)


def _iter_configured_anndata_inputs(config: SubsetConfig):
    """Yield ``(stem, path)`` for h5ad inputs the caller populated."""
    for field_name, stem in ANNDATA_INPUTS:
        path = getattr(config, field_name)
        if path is not None:
            yield stem, Path(path)


def _read_anndata(src_path: Path):
    """Lazy import wrapper — keeps anndata out of module-load cost."""
    import anndata as ad
    return ad.read_h5ad(str(src_path))


def _resolve_centroids_um(adata, side_name: str) -> np.ndarray:
    """Reuse the celltyping-stage resolver so A1 sees the same coords
    downstream stages do. Returns (n, 2) [x, y] in microns.
    """
    # Local import: hexenium.stages.celltyping pulls in scanpy at import
    # time; keep it out of A1's module load path.
    from hexenium.stages.celltyping import _resolve_spatial_coords
    return _resolve_spatial_coords(adata, side_name=side_name)


def _subset_anndata(
    src_path: Path,
    polygon_um: BaseGeometry,
    out_path: Path,
    anndata_geometry_frame: str,
) -> tuple[int, int]:
    """Subset an h5ad by ROI-polygon centroid membership.

    Opens the source in ``backed='r'`` so ``.X``, ``.layers[*]``, and
    ``.raw`` are HDF5 references rather than materialised arrays.
    Slicing on the backed view is lazy; ``to_memory()`` materialises
    only the kept rows — peak RAM tracks the subset size, not the
    full-file size. On the largest real proseg_raw observed in
    ``metx_liver_met/anndata/LiverMetx_samples/`` (~33 GB in memory
    on a full read), a ~4% ROI subset peaks at roughly 1-2 GB.

    Preserves ``obs``, ``var``, ``obsm``, ``uns``, ``layers``, and
    ``.raw`` via AnnData's own backed slice + ``to_memory()`` +
    ``write_h5ad`` chain. Coord frame is ``global`` today — no shift.
    Returns ``(n_kept, n_total)``.
    """
    if anndata_geometry_frame != "global":
        raise NotImplementedError(
            f"anndata_geometry_frame={anndata_geometry_frame!r} not implemented; "
            "only 'global' is supported today. Use the parquet/GeoJSON "
            "``geometry_frame`` for polygon-frame shifting; AnnData coord "
            "recentering is a future extension."
        )
    import anndata as ad

    adata = ad.read_h5ad(str(src_path), backed="r")
    try:
        coords = _resolve_centroids_um(adata, side_name=src_path.name)
        mask = _centroid_mask_xy(coords[:, 0], coords[:, 1], polygon_um)
        # to_memory() materialises the slice out of the HDF5 backing;
        # afterwards `sub` no longer references the source file.
        sub = adata[mask, :].to_memory()
    finally:
        f = getattr(adata, "file", None)
        if f is not None:
            f.close()
    sub.write_h5ad(str(out_path))
    return int(mask.sum()), int(len(mask))


def _annotate_source_h5ad(
    src_path: Path,
    rois_polygons_ordered: list[tuple[str, BaseGeometry]],
    out_path: Path,
) -> tuple[int, int]:
    """Write ``<out_path>`` = full uncut h5ad + ``obs['roi_annotation']``.

    Never mutates the source file. Reads it once, iterates ROIs in the
    order the caller supplies, computes per-cell membership against
    each ROI polygon, and joins names with commas for cells inside
    multiple ROIs. Cells outside every ROI get the literal ``outside``.

    Returns ``(n_annotated_inside_any_roi, n_total)``.
    """
    adata = _read_anndata(src_path)
    coords = _resolve_centroids_um(adata, side_name=src_path.name)
    n = int(len(coords))
    xs, ys = coords[:, 0], coords[:, 1]

    per_cell: list[list[str]] = [[] for _ in range(n)]
    for name, poly in rois_polygons_ordered:
        mask = _centroid_mask_xy(xs, ys, poly)
        for i in np.where(mask)[0]:
            per_cell[i].append(name)

    annotation = np.array(
        [",".join(m) if m else "outside" for m in per_cell],
        dtype=object,
    )
    adata.obs["roi_annotation"] = pd.Categorical(annotation)
    adata.write_h5ad(str(out_path))
    return int(sum(1 for m in per_cell if m)), n


def _subset_boundary_parquet(
    bundle: Path,
    stem: str,
    keep_ids: set[str],
    out_dir: Path,
) -> tuple[Path, int, int]:
    src = bundle / f"{stem}.parquet"
    df = pd.read_parquet(src)
    sub = df[df["cell_id"].astype(str).isin(keep_ids)]
    dst = out_dir / f"{stem}.parquet"
    sub.to_parquet(dst, index=False)
    return dst, len(sub), int(sub["cell_id"].nunique())


# ---------------------------------------------------------------------------
# DAPI crop
# ---------------------------------------------------------------------------

def _dapi_pixel_bounds(polygon_um: BaseGeometry, pixel_size_morph: float,
                      ) -> tuple[int, int, int, int]:
    """Return integer (x0, y0, x1, y1) px, clamped to non-negative x0/y0."""
    minx_um, miny_um, maxx_um, maxy_um = polygon_um.bounds
    x0 = max(int(minx_um / pixel_size_morph), 0)
    y0 = max(int(miny_um / pixel_size_morph), 0)
    x1 = int(maxx_um / pixel_size_morph)
    y1 = int(maxy_um / pixel_size_morph)
    return x0, y0, x1, y1


def _crop_dapi(
    dapi_path: Path,
    polygon_um: BaseGeometry,
    pixel_size_morph: float,
) -> tuple[np.ndarray, tuple[int, int, int, int]]:
    """Read only the tiles covering the ROI bbox — bounded memory.

    ``page.asarray()`` materialises the entire full-res DAPI page
    before slicing — ~8 GB in RAM for a Xenium (4, ~52k, ~40k) slide.
    ``page.aszarr()`` exposes the page as a zarr store backed by the
    on-disk tile decoders, so ``z[y0:y1, x0:x1]`` decodes only the
    tiles that intersect the bbox. Peak memory is
    O(roi_tiles × tile_size), independent of full-slide size.

    Xenium v4 splits morphology into per-channel OME-TIFFs
    (``morphology_focus/ch0000_dapi.ome.tif`` etc.); page 0 is 2-D.
    Legacy single-file ``morphology.ome.tif`` bundles four channels
    in one 3-D page — we take channel 0 to match the v4 output.
    """
    x0, y0, x1, y1 = _dapi_pixel_bounds(polygon_um, pixel_size_morph)
    with tifffile.TiffFile(dapi_path) as tf:
        page = tf.pages[0]
        store = page.aszarr()
        try:
            z = zarr.open(store, mode="r")
            shape = z.shape
            if len(shape) == 2:
                H, W = shape
                y1c, x1c = min(y1, H), min(x1, W)
                crop = np.asarray(z[y0:y1c, x0:x1c])
            elif len(shape) == 3:
                _, H, W = shape
                y1c, x1c = min(y1, H), min(x1, W)
                crop = np.asarray(z[0, y0:y1c, x0:x1c])
            else:
                raise ValueError(
                    f"unexpected DAPI page ndim={len(shape)} shape={shape}"
                )
        finally:
            store.close()
    return crop, (x0, y0, x1c, y1c)


def _resolve_dapi_path(bundle: Path) -> Path:
    """Locate the DAPI OME-TIFF within an XOA bundle.

    v4 layout: ``morphology_focus/ch0000_dapi.ome.tif``. Older bundles
    have a single-file ``morphology.ome.tif`` — that's the fallback.
    """
    v4 = bundle / "morphology_focus" / "ch0000_dapi.ome.tif"
    if v4.exists():
        return v4
    legacy = bundle / "morphology.ome.tif"
    if legacy.exists():
        return legacy
    raise FileNotFoundError(f"No DAPI OME-TIFF found under {bundle}")


# ---------------------------------------------------------------------------
# Manifest + geometry-frame shift
# ---------------------------------------------------------------------------

def _shift_to_local(polygon: BaseGeometry, origin: tuple[float, float]) -> BaseGeometry:
    from shapely.affinity import translate

    return translate(polygon, xoff=-origin[0], yoff=-origin[1])


def _write_manifest(
    manifest_path: Path,
    roi_name: str,
    polygon_um: BaseGeometry,
    dapi_crop_px_bounds: tuple[int, int, int, int] | None,
    pixel_size_morph: float,
    n_cells_kept: int,
    n_cells_total: int,
    sources: dict[str, Any],
    geometry_frame: str,
    anndata_subsets: dict[str, dict[str, Any]] | None = None,
) -> None:
    minx_um, miny_um, maxx_um, maxy_um = polygon_um.bounds
    manifest = {
        "roi_name": roi_name,
        "geometry_frame": geometry_frame,
        "pixel_size_morph_um_per_px": pixel_size_morph,
        "roi_bbox_microns": [minx_um, miny_um, maxx_um, maxy_um],
        "roi_area_microns_sq": float(polygon_um.area),
        "n_cells_kept": n_cells_kept,
        "n_cells_total": n_cells_total,
        "sources": {k: str(v) for k, v in sources.items()},
    }
    if dapi_crop_px_bounds is not None:
        x0, y0, x1, y1 = dapi_crop_px_bounds
        manifest["dapi_crop"] = {
            "origin_dapi_px": [x0, y0],
            "origin_dapi_um": [x0 * pixel_size_morph, y0 * pixel_size_morph],
            "size_dapi_px": [x1 - x0, y1 - y0],
        }
    if anndata_subsets:
        manifest["anndata_subsets"] = anndata_subsets
    with manifest_path.open("w") as f:
        yaml.safe_dump(manifest, f, sort_keys=False)


# ---------------------------------------------------------------------------
# A1 entry point
# ---------------------------------------------------------------------------

def run_a1(config: SubsetConfig) -> list[Path]:
    """Execute A1 — lightweight registration-input subset.

    For every ROI in the config, emits a directory containing the
    quartet a subsequent registration run needs: cropped DAPI OME-TIFF,
    subsetted cell + nucleus boundary parquets, and the source bundle's
    ``experiment.xenium`` copied verbatim. Also writes an optional
    subsetted ``cells.parquet`` and a per-ROI ``manifest.yaml``.

    Returns the list of per-ROI directories written.
    """
    if config.mode != "a1":
        raise ValueError(f"run_a1() called with mode={config.mode!r}")

    bundle = config.xenium_bundle.resolve()
    if not bundle.exists():
        raise FileNotFoundError(f"xenium_bundle not found: {bundle}")

    experiment_xenium = bundle / "experiment.xenium"
    pixel_size_morph = (
        config.pixel_size_morph
        if config.pixel_size_morph is not None
        else read_pixel_size_morph(experiment_xenium)
    )
    logger.info("[a1] pixel_size_morph = %s µm/px", pixel_size_morph)

    dapi_path = _resolve_dapi_path(bundle) if config.outputs.dapi_crop else None

    logger.info("[a1] loading cells.parquet")
    cells = _load_cells(bundle)

    out_root = config.output_dir.resolve()
    rois_dir = out_root / "rois"
    rois_dir.mkdir(parents=True, exist_ok=True)

    written_dirs: list[Path] = []
    seen_slugs: dict[str, int] = {}

    for spec in config.rois:
        slug = _slug(spec.name)
        if slug in seen_slugs:
            seen_slugs[slug] += 1
            slug = f"{slug}_{seen_slugs[slug]}"
        else:
            seen_slugs[slug] = 0

        roi_dir = rois_dir / slug
        roi_dir.mkdir(parents=True, exist_ok=True)

        polygon_um = _load_roi_polygon_um(spec, pixel_size_morph)
        logger.info(
            "[a1] roi=%r bbox_um=%s",
            spec.name, tuple(round(v, 1) for v in polygon_um.bounds),
        )

        mask = _centroid_mask(cells, polygon_um)
        kept_cells = cells.loc[mask]
        keep_ids = set(kept_cells["cell_id"].astype(str))
        logger.info(
            "[a1] roi=%r cells: %d of %d kept",
            spec.name, len(kept_cells), len(cells),
        )

        if config.outputs.cells_parquet:
            kept_cells.to_parquet(roi_dir / "cells.parquet", index=False)

        for stem in ("cell_boundaries", "nucleus_boundaries"):
            enabled = getattr(config.outputs, stem)
            if enabled:
                _subset_boundary_parquet(bundle, stem, keep_ids, roi_dir)

        if config.outputs.experiment_xenium and experiment_xenium.exists():
            shutil.copyfile(experiment_xenium, roi_dir / "experiment.xenium")

        dapi_crop_bounds: tuple[int, int, int, int] | None = None
        if config.outputs.dapi_crop and dapi_path is not None:
            crop, dapi_crop_bounds = _crop_dapi(dapi_path, polygon_um, pixel_size_morph)
            tifffile.imwrite(
                roi_dir / "dapi_subset.ome.tif",
                crop,
                photometric="minisblack",
                metadata={
                    "axes": "YX",
                    "PhysicalSizeX": pixel_size_morph,
                    "PhysicalSizeY": pixel_size_morph,
                    "PhysicalSizeXUnit": "µm",
                    "PhysicalSizeYUnit": "µm",
                },
                ome=True,
            )
            if config.outputs.dapi_polygon_mask:
                _write_polygon_mask(
                    roi_dir / "dapi_subset_mask.tif",
                    polygon_um,
                    pixel_size_morph,
                    dapi_crop_bounds,
                )

        # ROI polygon: emit in the caller's chosen geometry frame.
        # global = source-slide microns; local = shifted so (0, 0) is the crop origin.
        polygon_for_geojson = polygon_um
        if config.outputs.geometry_frame == "local" and dapi_crop_bounds is not None:
            x0, y0 = dapi_crop_bounds[0], dapi_crop_bounds[1]
            polygon_for_geojson = _shift_to_local(
                polygon_um, (x0 * pixel_size_morph, y0 * pixel_size_morph)
            )
        gpd.GeoDataFrame(
            {"roi_name": [spec.name]},
            geometry=[polygon_for_geojson],
        ).to_file(roi_dir / "roi.geojson", driver="GeoJSON")

        # Per-ROI h5ad subsets — only the inputs the caller declared.
        anndata_subsets: dict[str, dict[str, Any]] = {}
        for stem, src in _iter_configured_anndata_inputs(config):
            out = roi_dir / f"{stem}.h5ad"
            n_kept, n_total = _subset_anndata(
                src, polygon_um, out,
                anndata_geometry_frame=config.outputs.anndata_geometry_frame,
            )
            # QC: the h5ad subset for the xenium_ranger source should
            # match the cells.parquet subset exactly (same cells,
            # same ROI polygon). Log a mismatch; don't hard-fail
            # (proseg h5ads use different segmentations and are
            # expected to differ from cells.parquet).
            if stem == "xenium_ranger" and n_kept != len(kept_cells):
                logger.warning(
                    "[a1] roi=%r xenium_ranger.h5ad subset kept %d cells but "
                    "cells.parquet subset kept %d — mismatch suggests a "
                    "coord-frame drift between the h5ad and cells.parquet",
                    spec.name, n_kept, len(kept_cells),
                )
            anndata_subsets[stem] = {
                "source": str(src),
                "output": str(out),
                "n_kept": n_kept,
                "n_total": n_total,
            }
            logger.info(
                "[a1] roi=%r %s.h5ad cells: %d of %d kept",
                spec.name, stem, n_kept, n_total,
            )

        _write_manifest(
            roi_dir / "manifest.yaml",
            roi_name=spec.name,
            polygon_um=polygon_um,
            dapi_crop_px_bounds=dapi_crop_bounds,
            pixel_size_morph=pixel_size_morph,
            n_cells_kept=len(kept_cells),
            n_cells_total=len(cells),
            anndata_subsets=anndata_subsets,
            sources={
                "xenium_bundle": bundle,
                "roi_path": spec.path,
                "roi_source_tool": spec.source_tool,
                "roi_frame_declared": spec.frame,
            },
            geometry_frame=config.outputs.geometry_frame,
        )

        written_dirs.append(roi_dir)
        logger.info("[a1] roi=%r wrote %s", spec.name, roi_dir)

    if config.outputs.annotate_source_h5ad:
        _write_annotated_source_copies(config, out_root, pixel_size_morph)

    _write_run_summary(out_root, config, written_dirs)
    return written_dirs


def _write_annotated_source_copies(
    config: SubsetConfig,
    out_root: Path,
    pixel_size_morph: float,
) -> None:
    """One annotated copy per configured h5ad, holding all cells.

    Iterates ROIs in config order so overlaps stringify deterministically
    (`obs['roi_annotation'] == "Left,Center"`, never `"Center,Left"`).
    Source files are never modified — the copy lives at
    ``<output_dir>/<stem>_roi_annotated.h5ad``.
    """
    inputs = list(_iter_configured_anndata_inputs(config))
    if not inputs:
        return
    rois_polygons = [
        (spec.name, _load_roi_polygon_um(spec, pixel_size_morph))
        for spec in config.rois
    ]
    for stem, src in inputs:
        out = out_root / f"{stem}_roi_annotated.h5ad"
        n_inside, n_total = _annotate_source_h5ad(src, rois_polygons, out)
        logger.info(
            "[a1] annotated %s → %s (%d of %d cells inside at least one ROI)",
            src.name, out, n_inside, n_total,
        )


def _write_run_summary(
    out_root: Path,
    config: SubsetConfig,
    written_dirs: list[Path],
) -> None:
    summary = {
        "mode": config.mode,
        "xenium_bundle": str(config.xenium_bundle),
        "n_rois": len(written_dirs),
        "rois": [d.name for d in written_dirs],
        "outputs": asdict(config.outputs),
    }
    with (out_root / "resolved_config.yaml").open("w") as f:
        yaml.safe_dump(summary, f, sort_keys=False)


def _write_polygon_mask(
    dst: Path,
    polygon_um: BaseGeometry,
    pixel_size_morph: float,
    crop_bounds_px: tuple[int, int, int, int],
) -> None:
    """Rasterise the polygon into a uint8 mask matching the DAPI crop."""
    try:
        from rasterio.features import rasterize
        from rasterio.transform import Affine
    except ImportError as e:
        raise ImportError(
            "outputs.dapi_polygon_mask requires 'rasterio'; install it or "
            "disable the mask output."
        ) from e

    x0, y0, x1, y1 = crop_bounds_px
    height, width = y1 - y0, x1 - x0
    transform = Affine(
        pixel_size_morph, 0, x0 * pixel_size_morph,
        0, pixel_size_morph, y0 * pixel_size_morph,
    )
    mask = rasterize(
        [(polygon_um, 1)],
        out_shape=(height, width),
        transform=transform,
        fill=0,
        dtype=np.uint8,
    )
    tifffile.imwrite(dst, mask, photometric="minisblack")


__all__ = ["run_a1"]
