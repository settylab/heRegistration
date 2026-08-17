"""Stage 0: convert an Olympus VSI slide into an Xenium-Explorer-
compatible pyramidal OME-TIFF.

Follows 10x's Xenium Explorer image-conversion tutorial spec:
  - Tile size: 1024 px
  - Compression: JPEG 2000 (lossless, level=100) OR ZLIB
  - Pyramid scale: 2 (each level downsampled 2×)
  - Format: pyramidal OME-TIFF

Reader:  openslide (already in env; opens Olympus VSI natively).
Writer:  tifffile (matches 10x's tutorial writer path, level=100 for
         lossless JPEG 2000 — Xenium Explorer rejects lossy variants).

Stage 0 gets skipped (no-op) if `--he-path` isn't a `.vsi` file OR
if the output OME-TIFF already exists. Downstream stages register /
warp / celltype / viz then run against the converted OME-TIFF.
"""
from __future__ import annotations

import time
from pathlib import Path

from hexenium._internal.compat import apply_numpy_shims
from hexenium._internal.logging import log


VSI_EXTENSIONS = {".vsi"}


def is_vsi_input(he_path: Path) -> bool:
    """Return True iff we should route he_path through stage 0."""
    return he_path.suffix.lower() in VSI_EXTENSIONS


def _derive_ometiff_out_path(vsi_path: Path, out_root: Path, sample_id: str) -> Path:
    """Where the converted OME-TIFF lands under the pipeline output tree."""
    out_dir = out_root / sample_id / "he_preprocessed"
    return out_dir / f"{sample_id}_he.ome.tif"


def convert_vsi_to_ometiff(
    vsi_path: Path,
    out_path: Path,
    *,
    tile_size: int = 1024,
    compression: str = "jpeg2000",
    jpeg2000_level: int = 100,
    pyramid_levels: int = 7,
    force_rerun: bool = False,
) -> Path:
    """Convert a VSI file to a pyramidal OME-TIFF matching Xenium Explorer's
    spec. Returns the output path.

    Parameters
    ----------
    vsi_path
        Path to the Olympus VSI file.
    out_path
        Where the OME-TIFF should be written.
    tile_size
        OME-TIFF tile edge (px). 1024 per 10x's Explorer spec.
    compression
        `jpeg2000` (with `jpeg2000_level=100` for lossless) or `zlib`.
        Explorer's requirement is "ZLIB or lossless JPEG 2000".
    jpeg2000_level
        Only used when `compression='jpeg2000'`. 100 = lossless.
        Anything <100 is treated as lossy — Xenium Explorer rejects it.
    pyramid_levels
        Number of downsampled subresolutions to write (each half the
        edge of the prior). 7 matches Xenium Explorer's expectations.
    force_rerun
        If False (default) and out_path already exists, skip and
        return it — matches the resume-guard pattern used by the
        other stages.
    """
    apply_numpy_shims()

    if out_path.exists() and not force_rerun:
        log(f"[he_preprocess] skip — output already exists at {out_path}")
        return out_path

    import numpy as np
    import openslide
    import tifffile as tf
    import cv2

    if not vsi_path.exists():
        raise FileNotFoundError(f"VSI input not found: {vsi_path}")

    slide = openslide.OpenSlide(str(vsi_path))
    W, H = slide.dimensions
    vendor = slide.properties.get("openslide.vendor", "unknown")
    log(f"[he_preprocess] opened {vsi_path.name}: {W}×{H} px, "
        f"{slide.level_count} native pyramid levels, vendor={vendor}")

    # Physical pixel size from openslide's mpp-{x,y} properties (µm/px).
    # Xenium H&E is typically 0.25 µm/px (40x scan on 4k CMOS).
    try:
        mpp_x = float(slide.properties["openslide.mpp-x"])
        mpp_y = float(slide.properties["openslide.mpp-y"])
    except (KeyError, ValueError, TypeError):
        log("[he_preprocess] WARN: openslide.mpp-{x,y} missing/invalid; "
            "defaulting to 0.25 µm/px")
        mpp_x = mpp_y = 0.25
    px_per_cm_x = 1e4 / mpp_x
    px_per_cm_y = 1e4 / mpp_y
    log(f"[he_preprocess] pixel size: {mpp_x:.4f}×{mpp_y:.4f} µm/px "
        f"(→ {px_per_cm_x:.1f}×{px_per_cm_y:.1f} px/cm)")

    # Read level 0 into memory. On 128 GB request this is fine for
    # typical H&E (5-30 GB decoded RGB). If future VSI files are much
    # larger we can switch to tile-streaming — see NOTE below.
    n_gb = W * H * 3 / (1024**3)
    log(f"[he_preprocess] reading level 0 ({n_gb:.1f} GB decoded RGB)")
    t0 = time.time()
    img = np.asarray(slide.read_region((0, 0), 0, (W, H)))
    if img.shape[-1] == 4:
        img = img[..., :3]  # drop alpha
    log(f"[he_preprocess] read done in {time.time()-t0:.1f}s; shape={img.shape}")
    slide.close()

    # NOTE on memory: openslide.read_region returns a PIL RGBA which
    # np.asarray copies. Peak RSS is ~2× the decoded image size during
    # this step. For images >40 GB decoded, refactor to iterate over
    # 1024-tile grid: openslide.read_region((x, y), 0, (tile, tile)) in
    # a generator fed to tifffile.TiffWriter.write(data=<gen>).

    # Build tifffile options that match 10x's Explorer spec.
    options = dict(
        photometric="rgb",
        tile=(tile_size, tile_size),
        compression=compression,
        resolutionunit="CENTIMETER",
    )
    if compression == "jpeg2000":
        options["compressionargs"] = {"level": jpeg2000_level}
        if jpeg2000_level < 100:
            log(f"[he_preprocess] WARN: jpeg2000_level={jpeg2000_level} is LOSSY; "
                f"Xenium Explorer requires lossless (level=100 or compression='zlib')")

    metadata = {
        "PhysicalSizeX": mpp_x,
        "PhysicalSizeXUnit": "µm",
        "PhysicalSizeY": mpp_y,
        "PhysicalSizeYUnit": "µm",
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    log(f"[he_preprocess] writing OME-TIFF: {out_path.name} "
        f"(tile={tile_size}, compression={compression}"
        f"{f', level={jpeg2000_level}' if compression == 'jpeg2000' else ''}, "
        f"pyramid_levels={pyramid_levels})")

    with tf.TiffWriter(str(out_path), bigtiff=True) as writer:
        # Level 0 (main IFD) with SubIFDs for the pyramid.
        t0 = time.time()
        log(f"[he_preprocess]   level 0 ({img.shape[1]}×{img.shape[0]})")
        writer.write(
            img,
            subifds=pyramid_levels,
            resolution=(px_per_cm_x, px_per_cm_y),
            metadata=metadata,
            **options,
        )
        log(f"[he_preprocess]   level 0 done in {time.time()-t0:.1f}s")

        # Downsampled pyramid — each level is half the edge of the last.
        scale = 1.0
        for i in range(pyramid_levels):
            scale *= 0.5
            new_w = max(1, int(img.shape[1] * 0.5))
            new_h = max(1, int(img.shape[0] * 0.5))
            t0 = time.time()
            # cv2.INTER_AREA is the standard downsample interpolation.
            img = cv2.resize(img, (new_w, new_h), interpolation=cv2.INTER_AREA)
            log(f"[he_preprocess]   level {i+1} ({img.shape[1]}×{img.shape[0]}) "
                f"resize took {time.time()-t0:.1f}s")
            t0 = time.time()
            writer.write(
                img,
                subfiletype=1,
                resolution=(scale * px_per_cm_x, scale * px_per_cm_y),
                **options,
            )
            log(f"[he_preprocess]   level {i+1} write done in {time.time()-t0:.1f}s")

    file_gb = out_path.stat().st_size / (1024**3)
    log(f"[he_preprocess] done -> {out_path} ({file_gb:.2f} GB on disk)")
    return out_path


def run_he_preprocess(
    sample_id: str,
    he_path: Path,
    output_root: Path,
    *,
    tile_size: int = 1024,
    compression: str = "jpeg2000",
    jpeg2000_level: int = 100,
    pyramid_levels: int = 7,
    force_rerun: bool = False,
) -> Path:
    """Pipeline entrypoint. If he_path is a VSI, convert to OME-TIFF and
    return the converted path. If he_path is already an OME-TIFF, return
    it unchanged (no-op)."""
    if not is_vsi_input(he_path):
        log(f"[he_preprocess] input is not a VSI ({he_path.suffix}); skipping stage 0 — "
            f"downstream stages will use {he_path} as-is")
        return he_path

    out_path = _derive_ometiff_out_path(he_path, output_root, sample_id)
    return convert_vsi_to_ometiff(
        he_path, out_path,
        tile_size=tile_size,
        compression=compression,
        jpeg2000_level=jpeg2000_level,
        pyramid_levels=pyramid_levels,
        force_rerun=force_rerun,
    )
