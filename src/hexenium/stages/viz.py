"""Stage 4: visualise warped cell + nucleus boundaries on the H&E.

Produces a per-sample PNG: downsampled H&E in the background, cell
polygons coloured by classification, nucleus polygons drawn as a darker
outline. Uses matplotlib + openslide for the H&E thumbnail.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from hexenium._internal.logging import log


def _load_he_thumbnail(he_path: Path, max_dim: int) -> tuple[np.ndarray, float]:
    """Return (RGB ndarray, downsample_factor). downsample_factor maps
    full-res H&E pixels -> thumbnail pixels (multiply polygon coords by it)."""
    # Try OpenSlide first (handles big OME-TIFF / SVS), fall back to tifffile.
    try:
        import openslide
        slide = openslide.OpenSlide(str(he_path))
        full_w, full_h = slide.dimensions
        downsample = max(full_w, full_h) / max_dim
        if downsample < 1:
            downsample = 1.0
        level = slide.get_best_level_for_downsample(downsample)
        lvl_w, lvl_h = slide.level_dimensions[level]
        # Read whole level then resize to target.
        img = slide.read_region((0, 0), level, (lvl_w, lvl_h)).convert("RGB")
        target_w = int(full_w / downsample)
        target_h = int(full_h / downsample)
        img = img.resize((target_w, target_h))
        return np.asarray(img), 1.0 / downsample
    except Exception as exc:
        log(f"[viz] OpenSlide failed ({exc!r}), trying tifffile.")
        import tifffile
        with tifffile.TiffFile(str(he_path)) as tif:
            series = tif.series[0]
            # Pick the smallest pyramid level >= max_dim.
            best = series.levels[-1]  # smallest
            for lvl in series.levels:
                if max(lvl.shape[-2:]) >= max_dim:
                    best = lvl
                    break
            arr = best.asarray()
        # Normalize to HxWx3.
        if arr.ndim == 2:
            arr = np.stack([arr] * 3, axis=-1)
        elif arr.shape[0] in (3, 4) and arr.ndim == 3:
            arr = np.moveaxis(arr, 0, -1)
        if arr.shape[-1] == 4:
            arr = arr[..., :3]
        full_max = max(arr.shape[:2])
        # We don't actually know full-resolution here without OME metadata,
        # so assume polygons are in this thumbnail's pixel space already.
        return arr, full_max / max_dim if full_max > max_dim else 1.0


def _build_full_palette(
    labels,
    user_overrides: dict | None = None,
    cmap_name: str = "tab20",
) -> dict:
    """Return a dict mapping every label to an [R,G,B] triple in 0-255.

    Labels get colors in this precedence order:
      1. User-provided override (from YAML `viz.classification_palette`).
      2. "Unclassified" defaults to gray [180,180,180] if unspecified.
      3. Any remaining label draws the next unused color from a
         qualitative matplotlib colormap (`tab20` by default → 20 hues,
         wraps around for more). Assignment order = first-appearance in
         `labels`, so runs are deterministic.
    """
    import matplotlib.pyplot as plt

    user_overrides = dict(user_overrides or {})
    # Sensible default for "Unclassified" — only used if user hasn't set one.
    user_overrides.setdefault("Unclassified", [180, 180, 180])

    cmap = plt.get_cmap(cmap_name)
    n_cmap = cmap.N if hasattr(cmap, "N") else 20

    out: dict = {}
    cmap_i = 0
    for lbl in labels:
        if lbl is None:
            continue
        lbl = str(lbl)
        if lbl in out:
            continue
        if lbl in user_overrides:
            out[lbl] = list(user_overrides[lbl])
        else:
            rgba = cmap(cmap_i % n_cmap)
            out[lbl] = [int(round(rgba[0] * 255)),
                        int(round(rgba[1] * 255)),
                        int(round(rgba[2] * 255))]
            cmap_i += 1
    return out


def run_viz(
    sample_id: str,
    he_path: Path,
    warp_dir: Path,
    celltyped_dir: Path,
    out_dir: Path,
    *,
    thumbnail_max_dim: int = 4096,
    dpi: int = 200,
    cell_alpha: float = 0.3,
    nucleus_alpha: float = 0.6,
    classification_palette: dict | None = None,
    palette_cmap: str = "tab20",
    render_boundaries: str = "nucleus",
    force_rerun: bool = False,
) -> Path:
    """Draw warped boundaries on the H&E overlay.

    ``out_dir`` is the layout-computed
    ``<output_root_he>/viz/<he_job_id>/`` folder.
    ``celltyped_dir`` may point at THIS run's celltype output OR (via
    ``--celltype-run-id``) at an earlier run's — the caller resolves
    that choice from the RunLayout before invoking.
    """
    if render_boundaries not in ("cell", "nucleus", "both"):
        raise ValueError(
            f"render_boundaries must be one of cell|nucleus|both, got {render_boundaries!r}"
        )
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    import geopandas as gpd

    out_dir.mkdir(parents=True, exist_ok=True)
    sentinel = out_dir / f"{sample_id}_overlay.png"
    if sentinel.exists() and not force_rerun:
        log(f"[viz] skip — {sentinel} already exists")
        return sentinel

    parquet_path = celltyped_dir / f"{sample_id}_celltyped_wholeslide.parquet"
    if not parquet_path.exists():
        raise FileNotFoundError(
            f"Celltyped parquet not found at {parquet_path}; run celltype stage first."
        )
    gdf = gpd.read_parquet(parquet_path)

    # Read every actual classification label from the data first (in
    # first-appearance order), then build a palette that covers all of
    # them — using YAML overrides where present, auto-assigning from
    # `palette_cmap` where not.
    labels_in_data = gdf["classification"].dropna().unique().tolist()
    full_palette = _build_full_palette(
        labels_in_data,
        user_overrides=classification_palette,
        cmap_name=palette_cmap,
    )
    auto_labels = [lbl for lbl in labels_in_data
                   if lbl not in (classification_palette or {})]
    if auto_labels:
        log(f"[viz] auto-assigned {len(auto_labels)} celltype(s) from "
            f"colormap '{palette_cmap}': {auto_labels}")
    else:
        log(f"[viz] all {len(labels_in_data)} celltype(s) matched user palette")
    log(f"[viz] final palette: {full_palette}")
    classification_palette = full_palette

    log(f"[viz] loading H&E thumbnail (max dim {thumbnail_max_dim} px)…")
    img, scale = _load_he_thumbnail(he_path, thumbnail_max_dim)
    h, w = img.shape[:2]
    log(f"[viz] thumbnail {w}×{h}; polygon scale factor = {scale:.5f}")

    fig, ax = plt.subplots(figsize=(w / dpi, h / dpi), dpi=dpi)
    ax.imshow(img, origin="upper")

    def _scaled_xy(poly):
        try:
            return [(x * scale, y * scale) for x, y in poly.exterior.coords]
        except Exception:
            return []

    # Draw cell polygons first, nucleus outlines on top. `render_boundaries`
    # gates each branch.
    active = {"cell": ("cell", "both"), "nucleus": ("nucleus", "both")}
    seen_labels = set()
    for boundary_type, alpha in [("cell", cell_alpha), ("nucleus", nucleus_alpha)]:
        if render_boundaries not in active[boundary_type]:
            continue
        sub = gdf[gdf["boundary_type"] == boundary_type]
        if sub.empty:
            continue
        verts_by_class: dict[str, list] = {}
        for _, row in sub.iterrows():
            label = str(row.get("classification", "Unclassified"))
            verts_by_class.setdefault(label, []).append(_scaled_xy(row.geometry))
        for label, verts in verts_by_class.items():
            # classification_palette is now guaranteed to cover every
            # label present in the data (see _build_full_palette above).
            color = classification_palette.get(label) or classification_palette.get("Unclassified", [180, 180, 180])
            rgb = tuple(c / 255 for c in color)
            pc = PolyCollection(
                verts,
                facecolors=(*rgb, alpha) if boundary_type == "cell" else "none",
                edgecolors=(*rgb, min(1.0, alpha + 0.3)) if boundary_type == "nucleus" else "none",
                linewidths=0.2 if boundary_type == "nucleus" else 0,
            )
            ax.add_collection(pc)
            if label not in seen_labels:
                ax.plot([], [], "s", color=rgb, label=label)
                seen_labels.add(label)

    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.set_axis_off()
    _title_kind = {
        "cell": "cell polygons",
        "nucleus": "nucleus outlines",
        "both": "cell + nucleus boundaries",
    }[render_boundaries]
    ax.set_title(f"{sample_id} — warped {_title_kind} on H&E")
    if seen_labels:
        ax.legend(loc="upper right", fontsize=6, framealpha=0.8, ncol=1)
    fig.tight_layout(pad=0)
    fig.savefig(sentinel, bbox_inches="tight", pad_inches=0)
    plt.close(fig)

    log(f"[viz] wrote {sentinel}")
    return sentinel
