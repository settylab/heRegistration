"""Stage 3: assign cell-type labels to warped boundaries + write
four GeoJSON output files.

Sources of celltype label. The default (ranger-direct) reads labels
directly off ``xenium_ranger.h5ad``'s ``.obs[<celltype_col>]`` — no NN
involved. This module's inlined NN mapping is the LEGACY path: it only
runs when a ``proseg_purified.h5ad`` is passed *in addition to* the
ranger h5ad (see :func:`run_celltyping`'s source-mode routing). Every
Xenium cell is labelled by nearest-neighbour lookup on the proseg-side
centroids only in that legacy path.

  * ``proseg_purified.h5ad`` — SOURCE (legacy path only). Its ``.obs``
    carries centroids (``x``/``y``, Xenium µm frame) and a celltype
    column (``celltype``, ``first_type``, ``primary_cell_type``, or
    ``celltype_updated`` — auto-detected by precedence).
  * ``xenium_ranger.h5ad`` — QUERY. Its ``.obs`` carries the Xenium
    per-cell UUIDs (``cell_id``) and centroids for the ORIGINAL Xenium
    segmentation. In the legacy path, every Xenium cell gets a celltype
    label via spatial NN on the proseg centroids (no ID join); in the
    default ranger-direct path, labels come straight from this h5ad's
    own celltype column instead.

Algorithm:

  1. Read warped cells + warped nuclei parquets from the warp stage;
     reconstruct geometry from WKB.
  2. Read the two h5ads. NN-fit on proseg spatial coords, query with
     xenium spatial coords. Every xenium cell gets its nearest
     proseg cell's label. NO ``mark_unassigned`` policy — every
     Xenium cell is labelled.
  3. Auto-detect the celltype column (``celltype`` > ``first_type`` >
     ``primary_cell_type`` > ``celltype_updated``) and the Xenium-side
     id column (UUID-shape-first, to reject Proseg's int64-index masquerading
     as ``cell_id``).
  4. Join annotation to warped cells+nuclei on ``xenium_cell_id`` via
     the cumcount-merge trick (survives duplicate rows from dask).
  5. Emit FOUR GeoJSON outputs on two axes (target × precision):
       - ``<sample>_cells_analysis.geojson``   cells,   raw floats
       - ``<sample>_nuclei_analysis.geojson``  nuclei,  raw floats
       - ``<sample>_cells_qupath.geojson``     cells,   rounded + polygon-safe
       - ``<sample>_nuclei_qupath.geojson``    nuclei,  rounded + polygon-safe
     Plus a combined ``<sample>_celltyped_wholeslide.parquet`` for viz.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely import from_wkb, make_valid
from shapely.geometry import (
    GeometryCollection,
    MultiPolygon,
    Polygon,
)

from hexenium.layout import is_xenium_uuid
from hexenium._internal.logging import log


# Column-name candidates for auto-detection. Order = precedence.
_CELLTYPE_COL_CANDIDATES = (
    "celltype",            # roadmap's celltype_writeback output
    "first_type",          # RCTD top call (proseg_purified today)
    "primary_cell_type",   # purified-only alias of first_type
    "celltype_updated",    # allow-listed for future manual curation
)

# Xenium-side id-column candidates. Order deliberately prefers
# xenium-UUID-shaped column names — ``cell_id`` on Proseg's h5ad is an
# int64 index, so preferring ``cell_id`` first would silently produce
# all-Unclassified. Every candidate is ALSO shape-checked against
# ``is_xenium_uuid`` before use; a candidate whose values don't look
# like xenium UUIDs is rejected and the next candidate tried.
_XENIUM_ID_COL_CANDIDATES = (
    "xenium_cell_id",
    "xenium_cell_id_nn",
    "cell_id",
    "original_cell_id",
)


# Spatial-coord auto-detect: (x_col, y_col) pairs, checked against
# .obs in order after .obsm['spatial'] is checked first. Covers the
# three input conventions we've observed:
#   * scanpy / xenium AnnData in-place: coords in .obsm['spatial']
#   * proseg_to_anndata:  .obs['centroid_x'] / .obs['centroid_y']
#   * xenium_ranger:      .obs['x_centroid'] / .obs['y_centroid']
#   * generic fallback:   .obs['x'] / .obs['y']
_SPATIAL_COORD_OBS_PAIRS = (
    ("centroid_x", "centroid_y"),
    ("x_centroid", "y_centroid"),
    ("x", "y"),
)


#: The sentinel value hexenium uses in-memory (and downstream in
#: viz palettes) when a source h5ad has no usable celltype column.
#: Value + grey hex fixed by Tracy on ``settylab/TracyY123-nexus#15``
#: comment ``5352008910`` — do not rename without a matching palette
#: entry in ``viz._build_full_palette`` + doc update.
UNLABELED = "unlabeled"


def _resolve_celltype_col(adata_obs_cols, spec: str | None) -> str | None:
    """Pick the celltype column on the source h5ad.

    ``spec`` = ``None`` / ``"auto"`` triggers precedence-based
    auto-detect from ``_CELLTYPE_COL_CANDIDATES``. Anything else is
    treated as a literal column name and returned verbatim after an
    existence check.

    Returns the chosen column name OR ``None`` when nothing matched.
    Callers that need "give me labels" behavior should treat ``None``
    as "fall back to :data:`UNLABELED` for every row". Prior behavior
    was to raise ``KeyError`` here; the fallback semantic was requested
    by Tracy (`5351815301`) to avoid crashing celltype on samples that
    lack the writeback column.
    """
    cols = list(adata_obs_cols)
    if spec and spec != "auto":
        if spec not in cols:
            log(f"[celltype] WARN: source h5ad is missing requested "
                f"celltype column {spec!r}. Available .obs columns "
                f"(first 20): {cols[:20]}. Falling back to "
                f"{UNLABELED!r} for every row.")
            return None
        return spec
    for c in _CELLTYPE_COL_CANDIDATES:
        if c in cols:
            log(f"[celltype] auto-detected celltype column: {c!r}")
            return c
    log(f"[celltype] WARN: no celltype column found on source h5ad. "
        f"Tried (in order): {_CELLTYPE_COL_CANDIDATES}. "
        f"Available .obs columns (first 20): {cols[:20]}. "
        f"Falling back to {UNLABELED!r} for every row.")
    return None


def _read_labels_or_unlabeled(
    obs: pd.DataFrame, celltype_col: str | None, n_expected: int,
) -> np.ndarray:
    """Return an ``n_expected``-length object array of celltype labels.

    Three cases produce the ``UNLABELED`` sentinel across all rows:

    1. ``celltype_col is None`` — auto-detect / explicit-lookup failed;
       the resolver already logged WARN.
    2. Column exists but every value is null / empty (``pd.isna`` or
       stripped-empty string). Tracy's 2c: "an entirely-NaN/empty
       existing column also falls back to ``unlabeled``".
    3. Column exists but a specific row is null — that row falls to
       ``UNLABELED`` while non-null rows keep their labels.
    """
    if celltype_col is None:
        return np.full(n_expected, UNLABELED, dtype=object)
    raw = obs[celltype_col].astype(object).to_numpy()
    if len(raw) != n_expected:
        # Defensive — misalignment would corrupt the downstream join.
        raise ValueError(
            f"celltype column {celltype_col!r} has {len(raw)} values but "
            f"expected {n_expected} rows on the h5ad. Refusing to guess."
        )
    # Per-row NaN → sentinel. Stripped-empty strings also count as
    # missing so downstream palettes / GeoJSON consumers don't render
    # a blank legend entry.
    mask_missing = np.array(
        [(v is None) or (isinstance(v, float) and np.isnan(v))
         or (isinstance(v, str) and not v.strip())
         for v in raw],
        dtype=bool,
    )
    if mask_missing.all():
        log(f"[celltype] WARN: source h5ad column {celltype_col!r} is "
            f"entirely NaN / empty across {len(raw)} rows. Falling "
            f"back to {UNLABELED!r} for every row.")
        return np.full(n_expected, UNLABELED, dtype=object)
    if mask_missing.any():
        log(f"[celltype] {mask_missing.sum()}/{len(raw)} rows on column "
            f"{celltype_col!r} were NaN/empty; falling back to "
            f"{UNLABELED!r} for those rows.")
        raw = raw.copy()
        raw[mask_missing] = UNLABELED
    return raw


def _resolve_xenium_id_col(obs: pd.DataFrame, spec: str | None) -> str:
    """Pick the xenium-side id column, shape-checked against the xenium UUID regex.

    Returns the special sentinel ``"__index__"`` when the h5ad's
    ``.obs.index`` is UUID-shaped and no candidate column matches
    (fallback for xenium_ranger.h5ad where the UUID lives on the index).
    """
    cols = list(obs.columns)
    # Explicit spec: no shape-check (user knows what they're doing) but
    # still verify existence so a typo fails loud.
    if spec and spec != "auto":
        if spec == "__index__":
            return spec
        if spec not in cols:
            raise KeyError(
                f"xenium.h5ad is missing id column {spec!r}. "
                f"Available .obs columns: {cols}"
            )
        return spec

    # Auto-detect: try each candidate, shape-check the first ~10 values.
    def _probe(sample) -> bool:
        vals = [v for v in list(sample)[:10] if v is not None and pd.notna(v)]
        if not vals:
            return False
        # Require the majority (all in the sample) to look like xenium UUIDs.
        return all(is_xenium_uuid(v) for v in vals)

    for c in _XENIUM_ID_COL_CANDIDATES:
        if c not in cols:
            continue
        if _probe(obs[c].astype(str).values):
            log(f"[celltype] auto-detected xenium id column: {c!r} (UUID-shape check passed)")
            return c
        log(f"[celltype] xenium id column candidate {c!r} present but values do NOT "
            f"look like xenium UUIDs (first: {obs[c].astype(str).head(3).tolist()!r}); "
            "trying next candidate")

    # Last resort: the index itself.
    idx_head = obs.index.astype(str).tolist()[:10]
    if _probe(idx_head):
        log(f"[celltype] auto-detected xenium id from .obs.index (UUID-shape check passed)")
        return "__index__"

    raise KeyError(
        "could not auto-detect a xenium UUID-shaped id column in "
        f"xenium.h5ad. Tried columns {list(_XENIUM_ID_COL_CANDIDATES)!r} and "
        f".obs.index. First few index values: {idx_head!r}. "
        "Pass --id-col explicitly (or --id-col __index__ to force the index)."
    )


def _get_xenium_ids(obs: pd.DataFrame, id_col: str) -> np.ndarray:
    """Extract xenium UUID strings from ``obs`` given a resolved id-col."""
    if id_col == "__index__":
        return obs.index.astype(str).str.strip().to_numpy()
    return obs[id_col].astype(str).str.strip().to_numpy()


def _resolve_spatial_coords(
    adata,
    x_col: str | None = None,
    y_col: str | None = None,
    *,
    side_name: str = "h5ad",
) -> np.ndarray:
    """Return (n, 2) float xy centroids with an .obsm['spatial'] → .obs fallback chain.

    Fallback order (skipped entirely when both ``x_col`` and ``y_col``
    are given — explicit overrides win):

      1. ``adata.obsm['spatial']`` (first two cols) — the scanpy /
         Xenium-in-place convention.
      2. ``adata.obs[x_col] / .obs[y_col]`` for each pair in
         ``_SPATIAL_COORD_OBS_PAIRS`` (``centroid_x/y``, ``x_centroid/y_centroid``,
         ``x/y``).

    Raises ``ValueError`` with an actionable message listing what was
    tried and what's actually available if nothing matches.
    """
    obs = adata.obs
    if x_col and y_col:
        for c in (x_col, y_col):
            if c not in obs.columns:
                raise KeyError(
                    f"{side_name} is missing centroid column {c!r} "
                    f"(explicit override). Available .obs columns: "
                    f"{list(obs.columns)[:20]}. Available .obsm keys: "
                    f"{list(adata.obsm.keys())}. Unset the override to "
                    "fall back to auto-detect."
                )
        return obs[[x_col, y_col]].astype(float).to_numpy()

    if "spatial" in adata.obsm:
        arr = np.asarray(adata.obsm["spatial"])
        if arr.ndim != 2 or arr.shape[1] < 2:
            raise ValueError(
                f"{side_name} .obsm['spatial'] has shape {arr.shape}; "
                "expected (n, >=2). Cannot use as spatial coords."
            )
        log(f"[celltype] {side_name}: coords from .obsm['spatial'] (n={len(arr)})")
        return arr[:, :2].astype(float)

    for xc, yc in _SPATIAL_COORD_OBS_PAIRS:
        if xc in obs.columns and yc in obs.columns:
            log(f"[celltype] {side_name}: coords from .obs[{xc!r}] / .obs[{yc!r}]")
            return obs[[xc, yc]].astype(float).to_numpy()

    raise ValueError(
        f"could not resolve spatial coords for {side_name}. Tried "
        f".obsm['spatial'] and .obs pairs {list(_SPATIAL_COORD_OBS_PAIRS)!r}. "
        f"Available .obs cols (first 20): {list(obs.columns)[:20]}. "
        f"Available .obsm keys: {list(adata.obsm.keys())}. "
        "Pass explicit x_col / y_col to override."
    )


def _celltype_from_nn(
    *,
    proseg_purified_h5ad: Path,
    xenium_h5ad: Path,
    celltype_col: str | None,
    id_col: str | None,
    proseg_x_col: str | None = None,
    proseg_y_col: str | None = None,
    xenium_x_col: str | None = None,
    xenium_y_col: str | None = None,
    k: int = 1,
    algorithm: str = "auto",
    metric: str = "euclidean",
) -> pd.DataFrame:
    """NN-map celltype labels from proseg_purified → xenium.

    Returns a DataFrame with columns ``(xenium_cell_id, group,
    nn_distance)``. Every xenium cell that has a finite centroid is
    labelled (nearest proseg_purified cell's celltype). No unmatched
    policy — every xenium cell is labelled.
    """
    import anndata as ad
    from sklearn.neighbors import NearestNeighbors

    log(f"[celltype] reading proseg_purified: {proseg_purified_h5ad}")
    ap = ad.read_h5ad(proseg_purified_h5ad)
    obs_p = ap.obs
    ct_col = _resolve_celltype_col(obs_p.columns, celltype_col)
    xy_p = _resolve_spatial_coords(
        ap, proseg_x_col, proseg_y_col, side_name="proseg_purified.h5ad"
    )
    labels_p = _read_labels_or_unlabeled(obs_p, ct_col, len(xy_p))
    ok_p = np.isfinite(xy_p).all(axis=1)
    if not ok_p.all():
        log(f"[celltype] proseg: dropped {(~ok_p).sum()} rows with NaN centroid "
            f"(kept {ok_p.sum()}/{len(xy_p)})")
    xy_p = xy_p[ok_p]
    labels_p = labels_p[ok_p]
    if len(xy_p) == 0:
        raise ValueError(
            f"proseg_purified.h5ad {proseg_purified_h5ad} has 0 rows with finite "
            "centroids; cannot fit NN."
        )
    log(f"[celltype] proseg centroids: {len(xy_p)} rows; "
        f"celltype col {ct_col!r} → {pd.Series(labels_p).nunique()} unique labels")

    log(f"[celltype] reading xenium h5ad: {xenium_h5ad}")
    ax = ad.read_h5ad(xenium_h5ad)
    obs_x = ax.obs
    xid_col = _resolve_xenium_id_col(obs_x, id_col)
    ids_x = _get_xenium_ids(obs_x, xid_col)
    xy_x = _resolve_spatial_coords(
        ax, xenium_x_col, xenium_y_col, side_name="xenium.h5ad"
    )
    ok_x = np.isfinite(xy_x).all(axis=1)
    if not ok_x.all():
        log(f"[celltype] xenium: dropped {(~ok_x).sum()} rows with NaN centroid "
            f"(kept {ok_x.sum()}/{len(xy_x)})")
    xy_x = xy_x[ok_x]
    ids_x = ids_x[ok_x]
    log(f"[celltype] xenium centroids: {len(xy_x)} rows; id col {xid_col!r}")

    # Cheap frame-mismatch guard — points in genuinely different frames
    # (µm vs px, or px-in-full-res vs px-in-thumbnail) show wildly
    # different bounding boxes. Warn loudly rather than fail: false-
    # positives are cheap (log line) and false-negatives are silent
    # coord-frame bugs.
    def _bbox(arr: np.ndarray) -> tuple:
        return (
            float(arr[:, 0].min()), float(arr[:, 0].max()),
            float(arr[:, 1].min()), float(arr[:, 1].max()),
        )
    bp, bx = _bbox(xy_p), _bbox(xy_x)
    span_p = max(bp[1] - bp[0], bp[3] - bp[2])
    span_x = max(bx[1] - bx[0], bx[3] - bx[2])
    ratio = max(span_p, span_x) / max(min(span_p, span_x), 1e-9)
    log(f"[celltype] centroid bboxes: proseg={bp}  xenium={bx}  span_ratio={ratio:.2f}")
    if ratio > 10:
        log(f"[celltype] WARN: proseg/xenium centroid spans disagree by "
            f"{ratio:.1f}× — coord frames may not match. Cross-check the "
            "resolved centroid source on each side (see the two 'coords from …' "
            "log lines above).")

    log(f"[celltype] fitting NN (k={k}, algorithm={algorithm}, metric={metric}) on proseg")
    nbrs = NearestNeighbors(n_neighbors=k, algorithm=algorithm, metric=metric)
    nbrs.fit(xy_p)
    dists, idx = nbrs.kneighbors(xy_x)
    if k == 1:
        winning_idx = idx[:, 0]
        winning_dist = dists[:, 0]
    else:
        # Majority vote across the k, break ties by min-dist.
        from collections import Counter
        winning_idx = np.empty(len(idx), dtype=idx.dtype)
        winning_dist = np.empty(len(idx), dtype=dists.dtype)
        for i in range(len(idx)):
            row_labels = labels_p[idx[i]]
            counts = Counter(row_labels.tolist())
            top = max(counts.values())
            for j in range(k):
                if counts[row_labels[j]] == top:
                    winning_idx[i] = idx[i, j]
                    winning_dist[i] = dists[i, j]
                    break
    labels_out = labels_p[winning_idx]
    log(f"[celltype] NN done: {len(labels_out)} xenium cells labelled; "
        f"median distance {np.median(winning_dist):.2f}, "
        f"p95 {np.percentile(winning_dist, 95):.2f}")

    annot = pd.DataFrame({
        "xenium_cell_id": pd.Series(ids_x).astype(str).str.strip(),
        "group": pd.Series(labels_out).astype(object),
        "nn_distance": winning_dist.astype(float),
    })
    return annot


def _celltype_from_xenium_ranger_direct(
    *,
    xenium_h5ad: Path,
    celltype_col: str | None,
    id_col: str | None,
) -> pd.DataFrame:
    """Read celltype labels DIRECTLY from ``xenium_ranger.h5ad``'s ``.obs``.

    Post-xenium-preprocess rctd-split's ``celltype_writeback`` stage,
    ``xenium_ranger.h5ad`` carries per-cell celltype labels at
    ``.obs["celltype"]`` (default column from
    ``packages/rctd-split/config/default.yaml:172``). Reading them
    directly saves a whole NN fit + guarantees hexenium sees the
    same labels the upstream summary reports use — no drift risk.

    Returns a DataFrame with columns ``(xenium_cell_id, group,
    nn_distance)``. ``nn_distance`` is populated as ``NaN`` because
    this path involves no NN mapping — the label came straight off
    the query h5ad. Missing / all-NaN column falls to ``UNLABELED``
    per Tracy's spec (`5352008910`).
    """
    import anndata as ad

    log(f"[celltype] reading xenium_ranger (ranger-direct): {xenium_h5ad}")
    ax = ad.read_h5ad(xenium_h5ad)
    obs_x = ax.obs
    xid_col = _resolve_xenium_id_col(obs_x, id_col)
    ids_x = _get_xenium_ids(obs_x, xid_col)
    ct_col = _resolve_celltype_col(obs_x.columns, celltype_col)
    labels_x = _read_labels_or_unlabeled(obs_x, ct_col, len(ids_x))

    log(f"[celltype] ranger-direct: {len(ids_x)} cells labelled from "
        f".obs[{ct_col!r}] "
        f"→ {pd.Series(labels_x).nunique()} unique labels")

    return pd.DataFrame({
        "xenium_cell_id": pd.Series(ids_x).astype(str).str.strip(),
        "group": pd.Series(labels_x).astype(object),
        # No NN in this path; keep the column shape stable for
        # downstream callers that expect it.
        "nn_distance": np.full(len(ids_x), np.nan, dtype=float),
    })


def _read_warped_gdf(parquet_path: Path) -> gpd.GeoDataFrame:
    """Read a warped boundary parquet, hoist the dask index to a column,
    convert WKB → shapely geometry."""
    df = pd.read_parquet(parquet_path)
    df = df.reset_index().rename(columns={"__null_dask_index__": "xenium_cell_id"})
    gdf = gpd.GeoDataFrame(
        df,
        geometry=from_wkb(df["geometry"].values),
        crs=None,
    )
    gdf["xenium_cell_id"] = gdf["xenium_cell_id"].astype(str).str.strip()
    return gdf


def _cumcount_merge(gdf: gpd.GeoDataFrame, annot: pd.DataFrame) -> gpd.GeoDataFrame:
    """Pair-by-occurrence merge that survives duplicate xenium_cell_id rows
    coming out of the dask warp."""
    g = gdf.copy()
    a = annot.copy()
    g["_k"] = g.groupby("xenium_cell_id").cumcount()
    a["_k"] = a.groupby("xenium_cell_id").cumcount()
    out = (
        g[["xenium_cell_id", "_k"]]
        .merge(a, on=["xenium_cell_id", "_k"], how="left")
        .drop(columns="_k")
    )
    out["geometry"] = g["geometry"].values
    out = gpd.GeoDataFrame(out, geometry="geometry", crs=None)
    return out


def _fix_to_polygon(g):
    """make_valid → return a single Polygon (largest piece if multi)."""
    if g is None:
        return None
    try:
        if g.is_empty:
            return None
    except Exception:
        return None
    try:
        g2 = make_valid(g)
    except Exception:
        return None
    if g2 is None or g2.is_empty:
        return None
    if isinstance(g2, Polygon):
        return g2
    if isinstance(g2, MultiPolygon):
        polys = [p for p in g2.geoms if p is not None and not p.is_empty]
        return max(polys, key=lambda p: p.area) if polys else None
    if isinstance(g2, GeometryCollection):
        polys = []
        for x in g2.geoms:
            if isinstance(x, Polygon) and not x.is_empty:
                polys.append(x)
            elif isinstance(x, MultiPolygon):
                polys.extend(p for p in x.geoms if p is not None and not p.is_empty)
        return max(polys, key=lambda p: p.area) if polys else None
    return None


# -------------------------------------------------------------------
# Cleanup helpers. This ordering (make_valid → largest piece → drop
# empties/tiny → optional round → final validity check) has been
# validated end-to-end against QuPath's polygon importer.
# -------------------------------------------------------------------


def _round_polygon(poly: Polygon, ndigits: int) -> Polygon | None:
    """Round exterior ring coords to ndigits. Nothing else — no dedup,
    no close-ring, no buffer(0), no orient (all of which either
    over-processed or wrong-directioned in prior iterations)."""
    try:
        ext = [(round(x, ndigits), round(y, ndigits)) for x, y in poly.exterior.coords]
        return Polygon(ext)
    except Exception:
        return None


def _finite_coords(poly) -> bool:
    """True iff every exterior-ring coordinate is finite (no NaN/Inf)."""
    try:
        arr = np.asarray(poly.exterior.coords)
        return np.isfinite(arr).all()
    except Exception:
        return False


def _enough_vertices(poly) -> bool:
    """True iff exterior ring has >=4 total coords (3 corners + closure)
    AND >=3 unique coords."""
    try:
        coords = list(poly.exterior.coords)
        uniq = set(coords)
        return len(coords) >= 4 and len(uniq) >= 3
    except Exception:
        return False


def _refilter_gdf(
    subset: gpd.GeoDataFrame, mask, *, geometry_col: str = "geometry",
) -> gpd.GeoDataFrame:
    """``subset[mask].copy()`` re-typed as a GeoDataFrame.

    Guards against a geopandas quirk (surfaced by
    ``settylab/msetty-nexus#35`` comment ``5356740940`` — cross-user
    finding on real data): certain boolean-mask filters degrade an
    EMPTY result set from ``GeoDataFrame`` to plain ``DataFrame``,
    which then blows up the very next ``.geometry.apply()`` call
    with ``AttributeError: 'DataFrame' object has no attribute
    'geometry'``. Realistic triggers include a warp output whose
    polygons all fall below ``area_threshold_px``, small-sample
    runs, and corrupted / partial warp outputs.

    **Short-circuit on empty subset** (skeptic-suspect follow-up on
    ``5359956566``): the initial fix rewrapped TYPE but not SCHEMA.
    When ``subset`` is already empty, ``subset.geometry.apply(fn)``
    returns a geometry-dtype Series (pandas can't infer bool return
    dtype from a zero-row input), and using that as a boolean
    indexer collapses the frame to zero columns. The subsequent
    ``gpd.GeoDataFrame(zero_col_df, geometry="geometry")`` then
    raises ``ValueError: Unknown column geometry``. Returning
    ``subset`` unchanged when it's already empty avoids the whole
    apply-then-index-then-wrap chain and preserves the invariant
    (still a GeoDataFrame, still has ``geometry`` column, still
    zero rows).

    Non-empty subset: rewrap-if-needed via ``gpd.GeoDataFrame(...)``
    is cheap (constructor is a no-op on the already-typed input).
    """
    if len(subset) == 0:
        # Already empty and (by invariant, since this function is
        # the only writer) already a GeoDataFrame with a geometry
        # column. Skip the mask filter entirely — evaluating a
        # boolean mask against zero rows is meaningless anyway.
        return subset
    filtered = subset[mask].copy()
    if not isinstance(filtered, gpd.GeoDataFrame):
        filtered = gpd.GeoDataFrame(
            filtered, geometry=geometry_col, crs=None,
        )
    return filtered


def _clean_boundary_gdf(
    gdf: gpd.GeoDataFrame,
    id_col: str,
    class_col: str | None,
    boundary_type: str,
    area_threshold_px: float,
    round_ndigits: int | None,
    roi=None,
) -> gpd.GeoDataFrame:
    """Cleanup + optional rounding of a warped-boundary GeoDataFrame.

    Pipeline order (matters):
      1. Standardise `cell_id` (astype str).
      2. Standardise `classification` (NaN → "Unclassified"; astype str).
      3. Optional ROI filter (spatial index + intersects).
      4. First `fix_to_polygon` pass (make_valid → largest piece if multi).
      5. Optional ROI clip (intersection).
      6. Second `fix_to_polygon` pass (fixes anything intersection broke).
      7. Drop na / empty / invalid; keep only Polygon type.
      8. Drop `area <= area_threshold_px` (default 20 sq px).
      9. If `round_ndigits is not None`, round then re-drop
         na/empty/invalid.
     10. Final: `finite_coords + enough_vertices + is_valid`.
     11. Keep only `{cell_id, classification, boundary_type, geometry}`.

    `roi` is optional; when None the ROI-branch is dormant (whole-slide
    use).
    """
    subset = gdf.copy()
    subset = gpd.GeoDataFrame(subset, geometry="geometry", crs=None)

    # 1. Standardise id.
    subset["cell_id"] = subset[id_col].astype(str)

    # 2. Standardise classification. Historic behavior mapped NaN →
    # "Unclassified"; Tracy `5352008910` unified this to the single
    # ``UNLABELED`` sentinel so the viz palette only needs one grey
    # entry regardless of how the row ended up unlabeled (missing
    # column, all-NaN column, orphan cell id from a warp/h5ad
    # mismatch, or the no-h5ad path that synthesises UNLABELED here).
    if class_col is not None and class_col in subset.columns:
        subset["classification"] = subset[class_col]
    else:
        subset["classification"] = np.nan
    subset["classification"] = subset["classification"].astype("object")
    subset["classification"] = subset["classification"].where(
        ~pd.isna(subset["classification"]), UNLABELED
    )
    # Empty strings can slip through when a real h5ad has "" values
    # (Xenium ranger's celltype col sometimes carries empty strings
    # for cells that fell out of the reference-mapped RCTD grid).
    subset["classification"] = subset["classification"].astype(str).where(
        subset["classification"].astype(str).str.strip().astype(bool), UNLABELED
    )

    # 3. Optional ROI filter.
    if roi is not None:
        try:
            cand_idx = subset.sindex.query(roi, predicate="intersects")
            subset = subset.iloc[cand_idx].copy()
            subset = _refilter_gdf(subset, subset.intersects(roi))
        except Exception:
            subset = _refilter_gdf(subset, subset.intersects(roi))

    # 4. First fix_to_polygon.
    subset["geometry"] = subset["geometry"].apply(_fix_to_polygon)
    subset = _refilter_gdf(subset, subset.geometry.notna())

    # 5. Optional ROI clip.
    if roi is not None:
        subset["geometry"] = subset.geometry.intersection(roi)
        subset = _refilter_gdf(subset, subset.geometry.notna())
        subset = _refilter_gdf(subset, ~subset.geometry.is_empty)

    # 6. Re-fix after clip.
    subset["geometry"] = subset["geometry"].apply(_fix_to_polygon)
    subset = _refilter_gdf(subset, subset.geometry.notna())
    subset = _refilter_gdf(subset, ~subset.geometry.is_empty)
    subset = _refilter_gdf(subset, subset.is_valid)
    subset = _refilter_gdf(subset, subset.geometry.geom_type.isin(["Polygon"]))

    # 7. Drop small junk (config-driven, default 20 sq px).
    n_pre_area = len(subset)
    subset = _refilter_gdf(subset, subset.geometry.area > area_threshold_px)
    if len(subset) != n_pre_area:
        log(f"[celltype]   {boundary_type} area filter (>{area_threshold_px}) "
            f"dropped {n_pre_area - len(subset)}")

    # 8. Optional rounding.
    if round_ndigits is not None:
        n_pre_round = len(subset)
        subset["geometry"] = subset.geometry.apply(
            lambda g: _round_polygon(g, ndigits=round_ndigits))
        subset = _refilter_gdf(subset, subset.geometry.notna())
        subset = _refilter_gdf(subset, ~subset.geometry.is_empty)
        subset = _refilter_gdf(subset, subset.is_valid)
        if len(subset) != n_pre_round:
            log(f"[celltype]   {boundary_type} rounding@{round_ndigits} "
                f"dropped {n_pre_round - len(subset)}")

    # 9. Final sanity: finite coords + enough vertices + valid.
    # Guard the .geometry.apply() sites with _refilter_gdf so an empty
    # subset (all polygons dropped by earlier filters) doesn't crash
    # here with AttributeError — see cross-user bug on
    # settylab/msetty-nexus#35 comment 5356740940.
    n_pre_final = len(subset)
    subset = _refilter_gdf(subset, subset.geometry.apply(_finite_coords))
    subset = _refilter_gdf(subset, subset.geometry.apply(_enough_vertices))
    subset = _refilter_gdf(subset, subset.is_valid)
    if len(subset) != n_pre_final:
        log(f"[celltype]   {boundary_type} final sanity dropped {n_pre_final - len(subset)}")

    subset["boundary_type"] = boundary_type
    keep_cols = [c for c in ["cell_id", "classification", "boundary_type", "geometry"]
                 if c in subset.columns]
    return gpd.GeoDataFrame(subset[keep_cols].copy(), geometry="geometry", crs=None)


def _write_geojson(
    gdf: gpd.GeoDataFrame,
    path: Path,
    *,
    schema: str = "analysis",  # kept for API compat; both branches produce identical output now
    palette: dict | None = None,  # unused; retained for signature stability
) -> int:
    """Write a GeoDataFrame as a GeoJSON FeatureCollection.

    The format targeted here is a MINIMAL GeoJSON that QuPath's
    built-in Import Objects tool accepts:
      - `type: "FeatureCollection"` at root
      - Each feature: `{type: "Feature", properties: {...}, geometry: {...}}`
      - Properties: `{cell_id, classification (bare string), boundary_type}`
      - NO top-level `id`, NO `objectType`, NO `isLocked`, NO nested
        classification object with color.

    Uses geopandas' `.to_file(driver='GeoJSON')` (fiona/pyogrio under
    the hood). Returns feature count.
    """
    # Keep only the three properties that end up in QuPath / downstream:
    write_gdf = gdf[["cell_id", "classification", "boundary_type", "geometry"]].copy()
    # Silence GeoJSON-driver's "engine will infer CRS" warning by
    # setting a placeholder — QuPath ignores CRS anyway.
    try:
        write_gdf = write_gdf.set_crs(None, allow_override=True)
    except Exception:
        pass
    # Overwrite if exists (matches force_rerun semantics).
    if path.exists():
        path.unlink()
    write_gdf.to_file(str(path), driver="GeoJSON")
    return len(write_gdf)


def run_celltyping(
    sample_id: str,
    warp_dir: Path,
    out_dir: Path,
    *,
    xenium_h5ad: Path | None = None,
    proseg_purified_h5ad: Path | None = None,
    celltype_col: str | None = "auto",
    id_col: str | None = "auto",
    proseg_x_col: str | None = None,
    proseg_y_col: str | None = None,
    xenium_x_col: str | None = None,
    xenium_y_col: str | None = None,
    nn_k: int = 1,
    nn_algorithm: str = "auto",
    nn_metric: str = "euclidean",
    nuclei_inherit_classification: bool = True,
    area_threshold_px: float = 20.0,
    nucleus_round_ndigits: int | None = 2,
    cell_round_ndigits: int | None = None,
    qupath_cell_round_ndigits: int = 2,
    qupath_nucleus_round_ndigits: int = 2,
    force_rerun: bool = False,
) -> dict:
    """Run cell-typing + write four GeoJSONs + a combined parquet.

    ``xenium_h5ad`` (query) and ``proseg_purified_h5ad`` (source) drive
    the NN mapping. Both are optional at this level, and at the
    pipeline level too — ``pipeline.py`` does not require or validate
    either before calling this function (its former auto-derive/require
    helper was removed; see the comment at its celltype call site). If
    neither is set, cells go out with ``classification=UNLABELED``.

    ``out_dir`` is the layout-computed
    ``<output_root_he>/celltyped/<he_job_id>/`` folder (caller resolves
    this from the RunLayout).

    Parameters:
      cell_round_ndigits, nucleus_round_ndigits:
        Rounding applied to the **analysis** GeoJSONs. Default `None`
        for cells + `2` for nuclei.
      qupath_cell_round_ndigits, qupath_nucleus_round_ndigits:
        Rounding applied to the **QuPath** GeoJSONs. Default `2` for
        both — coarse enough that QuPath's polygon renderer accepts
        them, plus we run the polygon-safety cleanup.

    Returns:
      Dict of output paths + per-file counts.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # Sentinel = first of the four geojsons.
    sentinel = out_dir / f"{sample_id}_cells_analysis.geojson"
    if sentinel.exists() and not force_rerun:
        log(f"[celltype] skip — outputs already exist at {out_dir}")
        return {"geojson": sentinel, "out_dir": out_dir, "skipped_resume": True}

    cell_parquet = warp_dir / "he_cell_seg.parquet"
    nuc_parquet = warp_dir / "he_nucleus_seg.parquet"
    has_cells = cell_parquet.exists()
    has_nuclei = nuc_parquet.exists()
    if not (has_cells or has_nuclei):
        raise FileNotFoundError(
            f"Neither {cell_parquet} nor {nuc_parquet} exists; nothing to celltype."
        )

    # Source-mode routing (Tracy `5352008910`):
    #
    # * If ``--proseg-purified-h5ad`` was passed explicitly (present
    #   in ``proseg_purified_h5ad``) AND ``--xenium-h5ad`` is available,
    #   use the LEGACY proseg → xenium NN-mapping path. Kept for
    #   backwards compat with users who don't run xenium-preprocess's
    #   rctd-split celltype_writeback (or who want to override the
    #   ranger labels with a fresh NN fit).
    # * If only ``--xenium-h5ad`` is set (the typical post-rctd-split
    #   flow), read celltype labels DIRECTLY from ``.obs[<celltype_col>]``.
    #   No NN, no proseg. This is the new default.
    # * If NEITHER h5ad is set, we can't label anything — every row
    #   goes out as ``UNLABELED``. Preserves the previous "cells go
    #   out unclassified" ergonomics without crashing.
    if xenium_h5ad is not None and proseg_purified_h5ad is not None:
        annot = _celltype_from_nn(
            proseg_purified_h5ad=Path(proseg_purified_h5ad),
            xenium_h5ad=Path(xenium_h5ad),
            celltype_col=celltype_col,
            id_col=id_col,
            proseg_x_col=proseg_x_col,
            proseg_y_col=proseg_y_col,
            xenium_x_col=xenium_x_col,
            xenium_y_col=xenium_y_col,
            k=nn_k,
            algorithm=nn_algorithm,
            metric=nn_metric,
        )
        # Persist annotation for downstream inspection (also lets
        # subsequent viz-only runs skip the NN step by re-reading this).
        annot_path = out_dir / f"{sample_id}_celltype_annotation.parquet"
        annot.to_parquet(annot_path, index=False)
        log(f"[celltype] wrote annotation parquet ({len(annot)} rows) -> {annot_path.name}")
        # Downstream _cumcount_merge only needs (xenium_cell_id, group).
        annot = annot[["xenium_cell_id", "group"]]
    elif xenium_h5ad is not None:
        annot = _celltype_from_xenium_ranger_direct(
            xenium_h5ad=Path(xenium_h5ad),
            celltype_col=celltype_col,
            id_col=id_col,
        )
        annot_path = out_dir / f"{sample_id}_celltype_annotation.parquet"
        annot.to_parquet(annot_path, index=False)
        log(f"[celltype] wrote annotation parquet ({len(annot)} rows) -> {annot_path.name}")
        annot = annot[["xenium_cell_id", "group"]]
    else:
        annot = None
        log(f"[celltype] no --xenium-h5ad; classification will be "
            f"{UNLABELED!r} everywhere")

    # ------------------------------------------------------------------
    # Load + annotate + dedupe (once per boundary type, before precision
    # branching). Precision-specific cleanup happens in _clean_boundary_gdf.
    # ------------------------------------------------------------------
    cells_annotated = None
    nuclei_annotated = None

    if has_cells:
        cells_annotated = _read_warped_gdf(cell_parquet)
        log(f"[celltype] cells: {len(cells_annotated)} polygons, "
            f"{cells_annotated['xenium_cell_id'].nunique()} unique ids "
            f"({cells_annotated['xenium_cell_id'].duplicated().sum()} duplicate rows)")
        if annot is None:
            # No h5ad at all — synthesise "unlabeled" for every row so
            # the downstream viz stage still gets a `classification`
            # column (Tracy's Ask 2 semantics on no-source-h5ad).
            cells_annotated = cells_annotated.assign(group=UNLABELED)
        else:
            cells_annotated = _cumcount_merge(cells_annotated, annot)
        cells_annotated = cells_annotated.drop_duplicates(
            subset="xenium_cell_id", keep="first")

    if has_nuclei:
        nuclei_annotated = _read_warped_gdf(nuc_parquet)
        log(f"[celltype] nuclei: {len(nuclei_annotated)} polygons, "
            f"{nuclei_annotated['xenium_cell_id'].nunique()} unique ids "
            f"({nuclei_annotated['xenium_cell_id'].duplicated().sum()} duplicate rows)")
        if annot is None:
            nuclei_annotated = nuclei_annotated.assign(group=UNLABELED)
        else:
            nuclei_annotated = _cumcount_merge(nuclei_annotated, annot)
        nuclei_annotated = nuclei_annotated.drop_duplicates(
            subset="xenium_cell_id", keep="first")

    # `group` is now guaranteed to be present on both annotated frames
    # (either from the resolver / synthetic UNLABELED here, OR from a
    # real h5ad via _celltype_from_*). Downstream renderers can rely
    # on `class_col="group"` unconditionally — no more "unclassified
    # missing-column" branches to reason about.
    class_col = "group"

    # ------------------------------------------------------------------
    # Emit the four files. Each pass through _clean_boundary_gdf is
    # independent so raw analysis and QuPath variants don't interfere.
    # ------------------------------------------------------------------
    outputs: dict = {"out_dir": out_dir}

    if has_cells:
        cells_analysis = _clean_boundary_gdf(
            cells_annotated, id_col="xenium_cell_id", class_col=class_col,
            boundary_type="cell", area_threshold_px=area_threshold_px,
            round_ndigits=cell_round_ndigits,
        )
        cells_qupath = _clean_boundary_gdf(
            cells_annotated, id_col="xenium_cell_id", class_col=class_col,
            boundary_type="cell", area_threshold_px=area_threshold_px,
            round_ndigits=qupath_cell_round_ndigits,
        )
        p_ca = out_dir / f"{sample_id}_cells_analysis.geojson"
        p_cq = out_dir / f"{sample_id}_cells_qupath.geojson"
        outputs["cells_analysis_n"] = _write_geojson(cells_analysis, p_ca, schema="analysis")
        outputs["cells_qupath_n"] = _write_geojson(cells_qupath, p_cq, schema="qupath")
        outputs["cells_analysis"] = p_ca
        outputs["cells_qupath"] = p_cq
        log(f"[celltype] wrote {outputs['cells_analysis_n']} cell features "
            f"(no round) -> {p_ca.name}")
        log(f"[celltype] wrote {outputs['cells_qupath_n']} cell features "
            f"(round@{qupath_cell_round_ndigits}, qupath-safe) -> {p_cq.name}")

    if has_nuclei:
        # For nuclei, optionally inherit cell classification for
        # Unclassified rows before splitting into analysis + qupath.
        nuclei_to_write = nuclei_annotated
        if nuclei_inherit_classification and has_cells:
            # Build the cell classification lookup from the (unrounded)
            # annotated cells table.
            if class_col is not None:
                cell_lookup = cells_annotated[["xenium_cell_id", "group"]].copy()
                cell_lookup = cell_lookup.drop_duplicates(subset="xenium_cell_id")
                nuclei_to_write = nuclei_to_write.merge(
                    cell_lookup, on="xenium_cell_id", how="left",
                    suffixes=("", "_cell"),
                )
                nuclei_to_write["group"] = nuclei_to_write["group"].fillna(
                    nuclei_to_write["group_cell"])
                nuclei_to_write = nuclei_to_write.drop(columns=["group_cell"])

        nuclei_analysis = _clean_boundary_gdf(
            nuclei_to_write, id_col="xenium_cell_id", class_col=class_col,
            boundary_type="nucleus", area_threshold_px=area_threshold_px,
            round_ndigits=nucleus_round_ndigits,
        )
        nuclei_qupath = _clean_boundary_gdf(
            nuclei_to_write, id_col="xenium_cell_id", class_col=class_col,
            boundary_type="nucleus", area_threshold_px=area_threshold_px,
            round_ndigits=qupath_nucleus_round_ndigits,
        )
        p_na = out_dir / f"{sample_id}_nuclei_analysis.geojson"
        p_nq = out_dir / f"{sample_id}_nuclei_qupath.geojson"
        outputs["nuclei_analysis_n"] = _write_geojson(nuclei_analysis, p_na, schema="analysis")
        outputs["nuclei_qupath_n"] = _write_geojson(nuclei_qupath, p_nq, schema="qupath")
        outputs["nuclei_analysis"] = p_na
        outputs["nuclei_qupath"] = p_nq
        log(f"[celltype] wrote {outputs['nuclei_analysis_n']} nucleus features "
            f"(round@{nucleus_round_ndigits}) -> {p_na.name}")
        log(f"[celltype] wrote {outputs['nuclei_qupath_n']} nucleus features "
            f"(round@{qupath_nucleus_round_ndigits}, qupath-safe) -> {p_nq.name}")

    # ------------------------------------------------------------------
    # Combined parquet — kept for the viz stage which reads a single
    # parquet with boundary_type column to slice. Uses ANALYSIS variants
    # (raw floats) so viz gets the highest-precision geometry.
    # ------------------------------------------------------------------
    pieces = []
    if has_cells:
        pieces.append(cells_analysis)
    if has_nuclei:
        pieces.append(nuclei_analysis)
    combined = pd.concat(pieces, ignore_index=True)
    combined = gpd.GeoDataFrame(combined, geometry="geometry", crs=None)
    parquet_out = out_dir / f"{sample_id}_celltyped_wholeslide.parquet"
    combined.to_parquet(parquet_out, index=False)
    outputs["parquet"] = parquet_out
    log(f"[celltype] combined parquet (analysis-precision) -> {parquet_out.name}")

    # For backwards-compat with any downstream reader still expecting the
    # old sentinel name, keep pointing at the first file.
    outputs["geojson"] = outputs.get("cells_analysis") or outputs.get("nuclei_analysis")
    return outputs
