"""Vendored shim for HEST's ``warp_and_save_xenium_objects``.

HEST v1.2.0 (pinned via ``hest @ v1.2.0``, commit ``b87046f``) does NOT
export ``warp_and_save_xenium_objects`` — that function was added to
HEST main in the "hest V2" refactor (commit ``1bb61159``, Feb 2026),
which also pulls in a wider dependency shift (trident + hestcore
deprecation, spatialdata becoming optional, new dask/transcript
readers). Rather than force an env rebuild past a stable pin, this
shim composes v1.2.0's ``warp_gdf_valis`` primitive into the
higher-level entrypoint hexenium's warp stage expects.

Also works around a second v1.2.0 upstream gap: its ``warp_gdf_valis``
calls ``read_gdf(shapes)`` with NO reader kwargs, so the underlying
``XeniumParquetCellReader.pixel_size_morph`` stays at its constructor
default of ``None`` and the reader crashes on the first divide
(``TypeError: unsupported operand type(s) for /: 'float' and
'NoneType'`` at ``hest/io/seg_readers.py:129``). HEST main fixed
this by threading ``XENIUM_PIXEL_SIZE_MORPH = 0.2125`` through the
reader kwargs; v1.2.0 didn't. We pre-load the parquet ourselves with
the correct pixel size (read from the sample's ``experiment.xenium``
metadata when available, else the canonical Xenium morphology
default) and pass the resulting ``GeoDataFrame`` to
``warp_gdf_valis``, which routes through its
``elif isinstance(shapes, gpd.GeoDataFrame)`` branch and skips the
broken str-loading path entirely.

Third v1.2.0 gap: HEST main's ``warp_and_save_xenium_objects`` uses
dask under the hood, and its ``to_parquet`` writes the anonymous
cell-id index under Dask's ``__null_dask_index__`` sentinel column.
The hexenium celltype stage
(``hexenium.stages.celltyping._read_warped_gdf``) is written against
that shape: ``df.reset_index().rename(columns={"__null_dask_index__":
"xenium_cell_id"})``. Under v1.2.0's non-dask ``warp_gdf_valis``, our
plain-pandas ``to_parquet`` produces an unnamed index → the rename
misses → ``KeyError: 'xenium_cell_id'`` downstream. Fix: name the
warped GDF's index ``__null_dask_index__`` before writing, so the
resulting parquet carries the same shape the reader expects.

Scope:
    - Fully implements the cells / nuclei paths — matches the
      hexenium default (``warp.targets: [cells, nuclei]``).
    - Raises ``NotImplementedError`` when ``dapi_transcripts`` is
      set: HEST v1.2.0 has no ``XeniumTranscriptsReader`` / dask
      transcript path, and reimplementing that here would fork half
      of HEST V2. To enable transcripts, bump the HEST pin in
      ``environments/heRegistration-requirements.txt`` to a commit
      that includes the upstream function (>= ``1bb61159``).
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Optional

import geopandas as gpd

from hest.registration import warp_gdf_valis

from hexenium._internal.logging import log


# Xenium morphology sensor pitch (µm/pixel). Canonical fallback when
# a bundle's ``experiment.xenium`` metadata isn't reachable from the
# input parquet's directory. HEST main hard-codes the same constant
# in its ``warp_gdf_valis`` (``registration.py:166``).
_XENIUM_PIXEL_SIZE_MORPH_DEFAULT = 0.2125


def _ensure_jvm(mem_gb: int = 1) -> None:
    """Start the BioFormats JVM on the current process if it isn't up.

    v1.2.0's ``warp_gdf_valis`` runs single-threaded on the main
    process (no dask path in this tag) and loads slides via
    ``valis.registration``, which needs the JVM before any
    ``Slide.warp_xy_from_to`` call. The outer hexenium
    ``_dask_with_jvm`` context inits JVM on dask workers only.
    """
    import jpype
    if jpype.isJVMStarted():
        return
    try:
        from valis.registration import init_jvm
    except Exception:
        from valis_hest.registration import init_jvm
    init_jvm(mem_gb=mem_gb)


def _resolve_pixel_size_morph(parquet_path: str) -> float:
    """Resolve the Xenium morphology pixel size (µm/pixel) for a
    cell/nucleus boundaries parquet.

    Reads ``experiment.xenium`` next to the parquet (both live at the
    Xenium bundle root by 10x's output layout) and returns its
    ``pixel_size`` field. Falls back to the canonical Xenium
    morphology default (0.2125 µm/pixel) if the metadata file is
    absent, unreadable, or missing the field. Every Xenium slide
    shipped to date reports this value; the read-from-metadata path
    is future-proofing for slides that ever report differently.
    """
    parent = Path(parquet_path).parent
    meta = parent / "experiment.xenium"
    if not meta.exists():
        return _XENIUM_PIXEL_SIZE_MORPH_DEFAULT
    try:
        with open(meta) as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError):
        return _XENIUM_PIXEL_SIZE_MORPH_DEFAULT
    value = data.get("pixel_size")
    if not isinstance(value, (int, float)) or value <= 0:
        return _XENIUM_PIXEL_SIZE_MORPH_DEFAULT
    return float(value)


def _load_xenium_boundaries(parquet_path: str) -> gpd.GeoDataFrame:
    """Load a Xenium cell/nucleus boundaries parquet as a
    ``GeoDataFrame``, plumbing the resolved ``pixel_size_morph``
    into the reader.

    Bypasses HEST v1.2.0's ``warp_gdf_valis`` file-loading branch
    (which drops ``reader_kwargs``); the caller passes the returned
    ``GeoDataFrame`` into ``warp_gdf_valis`` instead.
    """
    from hest.io.seg_readers import read_gdf

    pixel_size_morph = _resolve_pixel_size_morph(parquet_path)
    return read_gdf(
        parquet_path,
        reader_kwargs={"pixel_size_morph": pixel_size_morph},
    )


def _write_warped_parquet_and_geojson(
    warped: gpd.GeoDataFrame,
    *,
    save_dir: str,
    stem: str,
    save_parquet: bool,
    save_geojson: bool,
) -> None:
    """Write the shim's warped ``GeoDataFrame`` to disk with the
    ``__null_dask_index__`` marker on the anonymous cell-id index.

    HEST main's ``warp_and_save_xenium_objects`` runs under Dask, and
    dask's ``to_parquet`` writes anonymous indices under the sentinel
    ``__null_dask_index__``. Downstream ``celltyping._read_warped_gdf``
    reads the parquet with ``df.reset_index().rename(columns={
    "__null_dask_index__": "xenium_cell_id"})`` — so we need to write
    that exact column name for the rename to hit.

    Plain-pandas ``to_parquet`` produces an unnamed index (→ column
    ``index`` after reset), which slips past the rename and crashes
    at ``gdf["xenium_cell_id"]`` in celltyping. Naming the index
    matches dask's shape without pulling dask into the write path.
    """
    # Copy so we don't mutate the caller's frame; index-name assignment
    # is cheap but still a mutation.
    df = warped.copy()
    df.index.name = "__null_dask_index__"
    if save_parquet:
        df.to_parquet(os.path.join(save_dir, f"{stem}.parquet"))
    if save_geojson:
        # v1.2.0's `hest.io.seg_readers.write_geojson` requires a
        # `category_key` column we don't carry; write via geopandas.
        df.to_file(
            os.path.join(save_dir, f"{stem}.geojson"),
            driver="GeoJSON",
        )


def warp_and_save_xenium_objects(
    path_registrar: str,
    dapi_path: str,
    save_dir: str,
    dapi_cells: Optional[str] = None,
    dapi_transcripts: Optional[str] = None,
    dapi_nuclei: Optional[str] = None,
    use_dask: bool = True,  # accepted for signature parity; ignored on v1.2.0
    verbose: bool = True,
    save_parquet: bool = True,
    save_geojson: bool = True,
) -> None:
    """Warp Xenium cells / nuclei to H&E and save. Signature-compatible
    with HEST main's ``warp_and_save_xenium_objects``.

    ``use_dask`` is accepted for signature parity but is a no-op on
    HEST v1.2.0's ``warp_gdf_valis`` (single-threaded main-process
    warp). The dask cluster hexenium spins up outside this call is
    still useful for its JVMPlugin on future dask-aware paths, but
    the actual warp here runs on the caller.
    """
    if not os.path.exists(save_dir):
        raise ValueError(f"Save directory '{save_dir}' doesn't exist.")

    if dapi_transcripts:
        raise NotImplementedError(
            "Warping xenium transcripts requires HEST >= 'hest V2' "
            "(commit 1bb61159, Feb 2026) for its "
            "`XeniumTranscriptsReader` + dask transcript path. This "
            "env pins hest @ v1.2.0 (commit b87046f), which does not "
            "expose the transcript reader classes. To enable "
            "transcripts, bump the HEST pin in "
            "environments/heRegistration-requirements.txt and rebuild "
            "the env; the hexenium wiring is otherwise ready."
        )

    _ensure_jvm()

    if dapi_cells is not None:
        if verbose:
            log("[hexenium hest-warp shim] warping cells DAPI -> H&E")
        cells_gdf = _load_xenium_boundaries(dapi_cells)
        warped_cells: gpd.GeoDataFrame = warp_gdf_valis(
            cells_gdf,
            path_registrar=path_registrar,
            curr_slide_name=dapi_path,
        )
        _write_warped_parquet_and_geojson(
            warped_cells,
            save_dir=save_dir,
            stem="he_cell_seg",
            save_parquet=save_parquet,
            save_geojson=save_geojson,
        )

    if dapi_nuclei is not None:
        if verbose:
            log("[hexenium hest-warp shim] warping nuclei DAPI -> H&E")
        nuclei_gdf = _load_xenium_boundaries(dapi_nuclei)
        warped_nuclei: gpd.GeoDataFrame = warp_gdf_valis(
            nuclei_gdf,
            path_registrar=path_registrar,
            curr_slide_name=dapi_path,
        )
        _write_warped_parquet_and_geojson(
            warped_nuclei,
            save_dir=save_dir,
            stem="he_nucleus_seg",
            save_parquet=save_parquet,
            save_geojson=save_geojson,
        )
