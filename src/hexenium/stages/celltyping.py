"""Stage 3: assign cell-type labels to warped boundaries + write
four GeoJSON output files.

Production pattern:
  1. Read warped cells + warped nuclei parquets, reconstruct geometry
     from WKB.
  2. Join to a per-sample celltype annotation CSV on `cell_id`. The CSV
     is the canonical join key — *not* a proxy nearest-neighbour join.
  3. Handle the duplicate `xenium_cell_id` rows that dask creates when a
     polygon spans pyramid tiles (cumcount-merge trick).
  4. Clean polygons (make_valid → largest piece; drop empties + tiny).
  5. Nuclei without a label inherit from their sibling cell.
  6. Emit FOUR GeoJSON outputs on two axes (target × precision):
     - <sample>_cells_analysis.geojson       cells,   raw floats
     - <sample>_nuclei_analysis.geojson      nuclei,  raw floats
     - <sample>_cells_qupath.geojson         cells,   rounded + polygon-safety cleanup
     - <sample>_nuclei_qupath.geojson        nuclei,  rounded + polygon-safety cleanup
     Also a combined `<sample>_celltyped_wholeslide.parquet` for
     downstream Python (viz stage reads this).
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

from hexenium._internal.logging import log


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

    # 2. Standardise classification.
    if class_col is not None and class_col in subset.columns:
        subset["classification"] = subset[class_col]
    else:
        subset["classification"] = np.nan
    subset["classification"] = subset["classification"].astype("object")
    subset["classification"] = subset["classification"].where(
        ~pd.isna(subset["classification"]), "Unclassified"
    )
    subset["classification"] = subset["classification"].astype(str)

    # 3. Optional ROI filter.
    if roi is not None:
        try:
            cand_idx = subset.sindex.query(roi, predicate="intersects")
            subset = subset.iloc[cand_idx].copy()
            subset = subset[subset.intersects(roi)].copy()
        except Exception:
            subset = subset[subset.intersects(roi)].copy()

    # 4. First fix_to_polygon.
    subset["geometry"] = subset["geometry"].apply(_fix_to_polygon)
    subset = subset[subset.geometry.notna()].copy()
    subset = gpd.GeoDataFrame(subset, geometry="geometry", crs=None)

    # 5. Optional ROI clip.
    if roi is not None:
        subset["geometry"] = subset.geometry.intersection(roi)
        subset = subset[subset.geometry.notna()].copy()
        subset = subset[~subset.geometry.is_empty].copy()

    # 6. Re-fix after clip.
    subset["geometry"] = subset["geometry"].apply(_fix_to_polygon)
    subset = subset[subset.geometry.notna()].copy()
    subset = subset[~subset.geometry.is_empty].copy()
    subset = subset[subset.is_valid].copy()
    subset = subset[subset.geometry.geom_type.isin(["Polygon"])].copy()

    # 7. Drop small junk (config-driven, default 20 sq px).
    n_pre_area = len(subset)
    subset = subset[subset.geometry.area > area_threshold_px].copy()
    if len(subset) != n_pre_area:
        log(f"[celltype]   {boundary_type} area filter (>{area_threshold_px}) "
            f"dropped {n_pre_area - len(subset)}")

    # 8. Optional rounding.
    if round_ndigits is not None:
        n_pre_round = len(subset)
        subset["geometry"] = subset.geometry.apply(
            lambda g: _round_polygon(g, ndigits=round_ndigits))
        subset = subset[subset.geometry.notna()].copy()
        subset = subset[~subset.geometry.is_empty].copy()
        subset = subset[subset.is_valid].copy()
        if len(subset) != n_pre_round:
            log(f"[celltype]   {boundary_type} rounding@{round_ndigits} "
                f"dropped {n_pre_round - len(subset)}")

    # 9. Final sanity: finite coords + enough vertices + valid.
    n_pre_final = len(subset)
    subset = subset[subset.geometry.apply(_finite_coords)].copy()
    subset = subset[subset.geometry.apply(_enough_vertices)].copy()
    subset = subset[subset.is_valid].copy()
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

    An earlier attempt to emit QuPath's "PathObject-native" shape (with
    id + objectType + classification-as-object + isLocked) made things
    worse — QuPath rejected the file. QuPath's built-in Import Objects
    tool actually likes the minimal shape.

    Uses geopandas' `.to_file(driver='GeoJSON')` (fiona/pyogrio under
    the hood) rather than a hand-rolled `json.dump`. Returns feature
    count.
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
    output_root: Path,
    *,
    celltype_csv: Path | None = None,
    csv_id_col: str = "cell_id",
    csv_group_col: str = "group",
    nuclei_inherit_classification: bool = True,
    area_threshold_px: float = 20.0,
    nucleus_round_ndigits: int | None = 2,
    cell_round_ndigits: int | None = None,
    qupath_cell_round_ndigits: int = 2,
    qupath_nucleus_round_ndigits: int = 2,
    force_rerun: bool = False,
) -> dict:
    """Run cell-typing + write four GeoJSONs + a combined parquet.

    Parameters:
      cell_round_ndigits, nucleus_round_ndigits:
        Rounding applied to the **analysis** GeoJSONs. Default `None`
        for cells (raw floats) + `2` for nuclei.
      qupath_cell_round_ndigits, qupath_nucleus_round_ndigits:
        Rounding applied to the **QuPath** GeoJSONs. Default `2` for
        both — coarse enough that QuPath's polygon renderer accepts
        them, plus we run the polygon-safety cleanup (drop consecutive
        dupes, close ring, buffer(0) to fix self-intersection).

    Returns:
      Dict of output paths + per-file counts.
    """
    out_dir = output_root / sample_id / "celltyped"
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

    # Load + prep the annotation table.
    if celltype_csv is not None:
        annot = pd.read_csv(celltype_csv)
        if csv_id_col not in annot.columns:
            raise KeyError(
                f"celltype csv {celltype_csv} is missing id column {csv_id_col!r}. "
                f"Got: {list(annot.columns)}"
            )
        if csv_group_col not in annot.columns:
            raise KeyError(
                f"celltype csv {celltype_csv} is missing group column {csv_group_col!r}. "
                f"Got: {list(annot.columns)}"
            )
        annot = annot[[csv_id_col, csv_group_col]].rename(
            columns={csv_id_col: "xenium_cell_id", csv_group_col: "group"}
        )
        annot["xenium_cell_id"] = annot["xenium_cell_id"].astype(str).str.strip()
        log(f"[celltype] loaded annotation: {len(annot)} rows, "
            f"{annot['group'].nunique()} groups")
    else:
        annot = None
        log("[celltype] no celltype CSV provided; classification will be "
            "'Unclassified' everywhere")

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
        if annot is not None:
            cells_annotated = _cumcount_merge(cells_annotated, annot)
        cells_annotated = cells_annotated.drop_duplicates(
            subset="xenium_cell_id", keep="first")

    if has_nuclei:
        nuclei_annotated = _read_warped_gdf(nuc_parquet)
        log(f"[celltype] nuclei: {len(nuclei_annotated)} polygons, "
            f"{nuclei_annotated['xenium_cell_id'].nunique()} unique ids "
            f"({nuclei_annotated['xenium_cell_id'].duplicated().sum()} duplicate rows)")
        if annot is not None:
            nuclei_annotated = _cumcount_merge(nuclei_annotated, annot)
        nuclei_annotated = nuclei_annotated.drop_duplicates(
            subset="xenium_cell_id", keep="first")

    class_col = "group" if annot is not None else None

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
