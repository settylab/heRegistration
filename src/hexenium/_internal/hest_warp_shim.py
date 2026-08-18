"""Vendored shim for HEST's ``warp_and_save_xenium_objects``.

HEST v1.2.0 (pinned via ``hest @ v1.2.0``, commit ``b87046f``) does NOT
export ``warp_and_save_xenium_objects`` — that function was added to
HEST main in the "hest V2" refactor (commit ``1bb61159``, Feb 2026),
which also pulls in a wider dependency shift (trident + hestcore
deprecation, spatialdata becoming optional, new dask/transcript
readers). Rather than force an env rebuild past a stable pin, this
shim composes v1.2.0's ``warp_gdf_valis`` primitive into the
higher-level entrypoint hexenium's warp stage expects.

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

import os
from typing import Optional

import geopandas as gpd

from hest.registration import warp_gdf_valis

from hexenium._internal.logging import log


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
        warped_cells: gpd.GeoDataFrame = warp_gdf_valis(
            dapi_cells,
            path_registrar=path_registrar,
            curr_slide_name=dapi_path,
        )
        if save_parquet:
            warped_cells.to_parquet(
                os.path.join(save_dir, "he_cell_seg.parquet")
            )
        if save_geojson:
            # v1.2.0's `hest.io.seg_readers.write_geojson` requires a
            # `category_key` column we don't carry; write via geopandas.
            warped_cells.to_file(
                os.path.join(save_dir, "he_cell_seg.geojson"),
                driver="GeoJSON",
            )

    if dapi_nuclei is not None:
        if verbose:
            log("[hexenium hest-warp shim] warping nuclei DAPI -> H&E")
        warped_nuclei: gpd.GeoDataFrame = warp_gdf_valis(
            dapi_nuclei,
            path_registrar=path_registrar,
            curr_slide_name=dapi_path,
        )
        if save_parquet:
            warped_nuclei.to_parquet(
                os.path.join(save_dir, "he_nucleus_seg.parquet")
            )
        if save_geojson:
            warped_nuclei.to_file(
                os.path.join(save_dir, "he_nucleus_seg.geojson"),
                driver="GeoJSON",
            )
