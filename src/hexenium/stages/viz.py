"""Stage 4: visualise warped cell + nucleus boundaries on the H&E.

Produces a per-sample PNG: downsampled H&E in the background, cell
polygons coloured by classification, nucleus polygons drawn as a darker
outline. Uses matplotlib + openslide for the H&E thumbnail.

Palette source (precedence, high→low):

  1. ``viz.classification_palette`` YAML override — user's explicit
     per-label colors always win. Manual overrides let the caller pin
     specific hues regardless of the upstream palette.
  2. Upstream xenium-preprocess color map — for the same sample/run
     pair, ``rctd_split`` writes
     ``<xenium_run_dir>/summary/<sample_id>_color_map.json`` with a
     ``first_type`` dict of ``label -> #rrggbb`` used to colour the
     summary UMAP. When available, hexenium reads it here so the H&E
     overlay uses the SAME celltype colors as the UMAP that already
     lives beside the run (Tracy's ask on
     settylab/TracyY123-nexus#15 comments 5348966457 / 5349046505).
  3. ``palette_cmap`` fallback (default ``tab20``) — labels not
     covered by (1) or (2) draw the next unused color from the cmap.

Missing / malformed upstream JSON is a graceful fall-through to the
cmap; the shim never crashes the viz stage on a palette-source read.
"""
from __future__ import annotations

import json
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


def _hex_to_rgb(hex_str: str) -> list[int] | None:
    """Decode ``#rrggbb`` / ``rrggbb`` to ``[R, G, B]`` (0-255).

    Returns ``None`` on any malformed input rather than raising —
    the caller silently drops bad entries and falls through to the
    cmap fallback, so a single corrupt row in the upstream JSON
    doesn't crash the viz stage.
    """
    if not isinstance(hex_str, str):
        return None
    s = hex_str.strip().lstrip("#")
    if len(s) != 6:
        return None
    try:
        return [int(s[0:2], 16), int(s[2:4], 16), int(s[4:6], 16)]
    except ValueError:
        return None


def _load_upstream_color_map(
    xenium_run_dir: Path | None,
    sample_id: str,
) -> dict[str, list[int]]:
    """Read the celltype palette that ``rctd_split.qc_report`` persists
    beside the xenium-preprocess UMAP.

    Path: ``<xenium_run_dir>/summary/<sample_id>_color_map.json``
    (write site: ``rctd_split.stages.qc_report`` around line 1409).
    The file's ``first_type`` key holds a ``{label: '#rrggbb'}`` dict.

    Returns a ``{label: [R,G,B]}`` dict when the file is readable and
    non-empty, else an empty dict. Every failure mode (no
    ``xenium_run_dir``, missing file, unreadable file, bad JSON,
    missing ``first_type`` key, malformed hex entries) logs and
    returns ``{}`` so ``_build_full_palette`` cleanly falls through
    to the cmap fallback.
    """
    if xenium_run_dir is None:
        return {}
    path = Path(xenium_run_dir) / "summary" / f"{sample_id}_color_map.json"
    if not path.exists():
        log(f"[viz] palette: no upstream color map at {path} — "
            f"falling back to cmap.")
        return {}
    try:
        with open(path) as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError) as exc:
        log(f"[viz] palette: could not read upstream color map "
            f"{path} ({exc!r}) — falling back to cmap.")
        return {}
    raw = data.get("first_type") if isinstance(data, dict) else None
    if not isinstance(raw, dict) or not raw:
        log(f"[viz] palette: upstream color map at {path} has no "
            f"non-empty 'first_type' block — falling back to cmap.")
        return {}
    out: dict[str, list[int]] = {}
    for label, hex_str in raw.items():
        rgb = _hex_to_rgb(hex_str)
        if rgb is not None:
            out[str(label)] = rgb
    log(f"[viz] palette: loaded {len(out)} celltype color(s) from "
        f"upstream {path}.")
    return out


# RGB triple for the "unlabeled" sentinel. Hex #888888 = [136,136,136].
# Locked in by Tracy on ``settylab/TracyY123-nexus#15`` comment
# ``5352008910`` — keep in sync with ``celltyping.UNLABELED`` so the
# viz palette entry matches whatever the celltype stage writes.
UNLABELED = "unlabeled"
UNLABELED_RGB = [0x88, 0x88, 0x88]

# Historic label from the pre-``26e78b0`` codebase. Legacy celltyped
# parquets still on disk carry ``"Unclassified"`` rows. We do NOT
# rewrite those files — instead we fold the legacy label into the
# canonical ``UNLABELED`` at viz time so old and new runs render the
# same shade of grey. Tracy on ``settylab/TracyY123-nexus#15`` comment
# ``5360206242``.
_LEGACY_UNCLASSIFIED = "Unclassified"


def _canonicalize_classification(label) -> str:
    """Return the canonical viz label for a ``.classification`` value.

    Folds the historic ``"Unclassified"`` sentinel into ``UNLABELED``
    so all "no celltype" cells resolve to the same palette entry
    (``#888888``) regardless of when the parquet was written.
    ``None`` and NaN (pandas' float representation of a missing entry
    in an object column) also fold to ``UNLABELED``.
    """
    if label is None:
        return UNLABELED
    # NaN is the only float where x != x — dodge pandas' None→NaN
    # coercion in mixed-type columns without importing pandas here.
    try:
        if label != label:
            return UNLABELED
    except Exception:
        pass
    s = str(label)
    if s == _LEGACY_UNCLASSIFIED:
        return UNLABELED
    return s


def _build_full_palette(
    labels,
    user_overrides: dict | None = None,
    upstream_palette: dict | None = None,
    cmap_name: str = "tab20",
) -> dict:
    """Return a dict mapping every label to an [R,G,B] triple in 0-255.

    Labels get colors in this precedence order:
      1. User-provided override (from YAML `viz.classification_palette`).
      2. ``"unlabeled"`` defaults to ``#888888`` (medium grey) if
         unspecified by (1) — Tracy `5352008910`.
      3. Upstream cross-pipeline palette (`upstream_palette` — the
         xenium-preprocess color map, keyed by celltype label). Fills
         labels not covered by (1) so the H&E overlay matches the UMAP
         for the same sample/run.
      4. Any remaining label draws the next unused color from a
         qualitative matplotlib colormap (`tab20` by default → 20 hues,
         wraps around for more). Assignment order = first-appearance in
         `labels`, so runs are deterministic.

    Labels are canonicalized via :func:`_canonicalize_classification`
    before palette lookup — the historic ``"Unclassified"`` sentinel
    folds into ``UNLABELED`` so old celltyped parquets render at the
    same ``#888888`` as fresh ``unlabeled`` runs (B+fold decision,
    Tracy `5360206242`).
    """
    import matplotlib.pyplot as plt

    # Canonicalize the overrides dict too so a YAML entry keyed on the
    # legacy "Unclassified" name still lands on the UNLABELED slot.
    user_overrides = {
        _canonicalize_classification(k): v
        for k, v in (user_overrides or {}).items()
    }
    user_overrides.setdefault(UNLABELED, list(UNLABELED_RGB))
    upstream_palette = {
        _canonicalize_classification(k): v
        for k, v in (upstream_palette or {}).items()
    }

    cmap = plt.get_cmap(cmap_name)
    n_cmap = cmap.N if hasattr(cmap, "N") else 20

    out: dict = {}
    cmap_i = 0
    for lbl in labels:
        lbl = _canonicalize_classification(lbl)
        if lbl in out:
            continue
        if lbl in user_overrides:
            out[lbl] = list(user_overrides[lbl])
        elif lbl in upstream_palette:
            out[lbl] = list(upstream_palette[lbl])
        else:
            rgba = cmap(cmap_i % n_cmap)
            out[lbl] = [int(round(rgba[0] * 255)),
                        int(round(rgba[1] * 255)),
                        int(round(rgba[2] * 255))]
            cmap_i += 1
    return out


def _synthesise_unlabeled_from_warp(
    *, warp_dir: Path, sample_id: str, render_boundaries: str,
):
    """Build the wholeslide GeoDataFrame the celltype stage WOULD have
    written, but with every row labelled :data:`UNLABELED`.

    Called from :func:`run_viz` when the celltype stage did NOT run in
    this invocation (Tracy `5352008910` Ask 4). Reads whichever warp
    parquets exist under ``warp_dir`` and reconstructs geometry from
    WKB the same way ``celltyping._read_warped_gdf`` does — but
    without importing the celltype stage (which would pull sklearn
    etc. for a viz-only run).

    ``render_boundaries`` selects which of ``he_cell_seg.parquet`` /
    ``he_nucleus_seg.parquet`` to load. Errors LOUD if the requested
    boundary type has no warp parquet on disk — a viz-only invocation
    that asks for boundaries the warp stage never wrote is a config
    mistake worth surfacing early.
    """
    import geopandas as gpd
    import pandas as pd
    from shapely import from_wkb

    from hexenium.stages.celltyping import UNLABELED as _CT_UNLABELED

    # Match constants across modules — if they ever diverge, prefer
    # the celltyping side since that's where the label is written to
    # disk in the normal flow.
    assert _CT_UNLABELED == UNLABELED, (
        "viz.UNLABELED and celltyping.UNLABELED disagree; sentinel drift "
        "would produce a broken palette. Keep them in sync."
    )

    cell_pq = warp_dir / "he_cell_seg.parquet"
    nuc_pq = warp_dir / "he_nucleus_seg.parquet"

    want_cells = render_boundaries in ("cell", "both")
    want_nuclei = render_boundaries in ("nucleus", "both")

    frames: list[gpd.GeoDataFrame] = []
    if want_cells:
        if not cell_pq.exists():
            raise FileNotFoundError(
                f"viz --render-boundaries={render_boundaries!r} requires "
                f"a cell warp parquet; {cell_pq} is missing. Include "
                f"`warp` in --stages or re-run warp with cells in --warp-targets."
            )
        df = pd.read_parquet(cell_pq).reset_index().rename(
            columns={"__null_dask_index__": "cell_id"},
        )
        df["boundary_type"] = "cell"
        df["classification"] = UNLABELED
        df["cell_id"] = df["cell_id"].astype(str)
        frames.append(gpd.GeoDataFrame(
            df[["cell_id", "classification", "boundary_type"]],
            geometry=from_wkb(df["geometry"].values), crs=None,
        ))
    if want_nuclei:
        if not nuc_pq.exists():
            raise FileNotFoundError(
                f"viz --render-boundaries={render_boundaries!r} requires "
                f"a nucleus warp parquet; {nuc_pq} is missing. Include "
                f"`warp` in --stages or re-run warp with nuclei in --warp-targets."
            )
        df = pd.read_parquet(nuc_pq).reset_index().rename(
            columns={"__null_dask_index__": "cell_id"},
        )
        df["boundary_type"] = "nucleus"
        df["classification"] = UNLABELED
        df["cell_id"] = df["cell_id"].astype(str)
        frames.append(gpd.GeoDataFrame(
            df[["cell_id", "classification", "boundary_type"]],
            geometry=from_wkb(df["geometry"].values), crs=None,
        ))
    combined = gpd.GeoDataFrame(
        pd.concat(frames, ignore_index=True),
        geometry="geometry", crs=None,
    )
    log(f"[viz] synthesised {len(combined)} boundaries from warp "
        f"(render_boundaries={render_boundaries!r}, all labelled "
        f"{UNLABELED!r})")
    return combined


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
    xenium_run_dir: Path | None = None,
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
    if parquet_path.exists():
        gdf = gpd.read_parquet(parquet_path)
    else:
        # Fallback path per Tracy `5352008910` Ask 4: viz can run
        # without a preceding celltype stage as long as the warp
        # parquets are on disk. Synthesise a wholeslide GeoDataFrame
        # from the warp outputs; every polygon gets
        # ``classification = UNLABELED`` so the palette maps to grey.
        # ``render_boundaries`` is respected (Tracy's 4a: default
        # ``nucleus`` unless the operator explicitly overrode).
        log(f"[viz] celltyped parquet absent at {parquet_path} — "
            f"falling back to warp parquets directly (all polygons "
            f"labelled {UNLABELED!r}).")
        gdf = _synthesise_unlabeled_from_warp(
            warp_dir=warp_dir, sample_id=sample_id,
            render_boundaries=render_boundaries,
        )

    # Read every actual classification label from the data first (in
    # first-appearance order), then build a palette that covers all of
    # them. Precedence: YAML overrides → upstream xenium-preprocess
    # color map (matches the summary UMAP for the same sample/run) →
    # `palette_cmap` fallback for anything still uncovered.
    # Canonicalize BEFORE the palette build so a legacy "Unclassified"
    # row folds into UNLABELED and hits the same #888888 entry as a
    # fresh "unlabeled" row (B+fold, Tracy `5360206242`).
    raw_labels_in_data = gdf["classification"].dropna().unique().tolist()
    labels_in_data = [_canonicalize_classification(lbl)
                      for lbl in raw_labels_in_data]
    # Preserve first-appearance order without duplicates.
    seen = set()
    labels_in_data = [lbl for lbl in labels_in_data
                      if not (lbl in seen or seen.add(lbl))]
    upstream_palette = _load_upstream_color_map(xenium_run_dir, sample_id)
    full_palette = _build_full_palette(
        labels_in_data,
        user_overrides=classification_palette,
        upstream_palette=upstream_palette,
        cmap_name=palette_cmap,
    )
    yaml_overrides = classification_palette or {}
    yaml_matched = [lbl for lbl in labels_in_data if lbl in yaml_overrides]
    upstream_matched = [lbl for lbl in labels_in_data
                        if lbl not in yaml_overrides and lbl in upstream_palette]
    auto_labels = [lbl for lbl in labels_in_data
                   if lbl not in yaml_overrides and lbl not in upstream_palette
                   and lbl != UNLABELED]
    if yaml_matched:
        log(f"[viz] palette: {len(yaml_matched)} label(s) from YAML "
            f"classification_palette: {yaml_matched}")
    if upstream_matched:
        log(f"[viz] palette: {len(upstream_matched)} label(s) inherited "
            f"from upstream xenium color map: {upstream_matched}")
    if UNLABELED in labels_in_data and UNLABELED not in yaml_overrides:
        legacy = _LEGACY_UNCLASSIFIED in raw_labels_in_data
        log(f"[viz] palette: 'unlabeled' → #888888 (default sentinel)"
            + (f" — folded {_LEGACY_UNCLASSIFIED!r} legacy rows in as well"
               if legacy else ""))
    if auto_labels:
        log(f"[viz] palette: {len(auto_labels)} label(s) auto-assigned "
            f"from cmap '{palette_cmap}': {auto_labels}")
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
            # Canonicalize BEFORE palette lookup so a legacy
            # "Unclassified" row hits the same UNLABELED palette entry
            # as a fresh "unlabeled" row — no separate rendering path.
            label = _canonicalize_classification(row.get("classification"))
            verts_by_class.setdefault(label, []).append(_scaled_xy(row.geometry))
        for label, verts in verts_by_class.items():
            # classification_palette is now guaranteed to cover every
            # label present in the data (see _build_full_palette above,
            # which canonicalises labels the same way).
            color = classification_palette.get(label) or classification_palette.get(UNLABELED, list(UNLABELED_RGB))
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
