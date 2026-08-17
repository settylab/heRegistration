"""Stage: map proseg-derived cell types to Xenium cells via nearest-neighbor.

Config-driven; every sample-specific literal (paths, column names, K,
distance policy) is a config knob.

Algorithm:

  1. Load proseg-side centroids + celltype (from purified.h5ad `.obs`
     by default; a CSV also works — auto-detected from extension).
  2. Load Xenium-side centroids from ``<xenium_bundle>/cells.csv.gz``
     (index = xenium ``cell_id``; columns ``x_centroid``, ``y_centroid``).
  3. Fit ``sklearn.neighbors.NearestNeighbors`` on the proseg centroids.
  4. Query with the Xenium centroids; propagate proseg's celltype
     label back onto each Xenium cell.
  5. Emit ``<sample_id>_celltype.csv`` with columns ``(cell_id, group)``
     — the exact shape hexenium's ``--celltype-csv`` consumes.

Optional knobs beyond the notebook default:

  * ``nn.k`` + ``nn.tiebreak`` — K>1 with majority-vote tie-break.
  * ``nn.distance_threshold`` + ``nn.unmatched_policy`` — cap the max
    NN distance; drop, keep, or mark ``Unassigned`` cells beyond it.
  * ``nn.hybrid_direct_join_first`` — for Xenium cells that already
    have a matching ``original_cell_id`` on the proseg side, use the
    direct join; fall back to NN only for the rest. Accepts ``True``
    / ``False`` for manual control, or ``"auto"`` (default) to
    auto-detect viability from the loaded data via
    :func:`_should_enable_hybrid` — enables the hybrid path only when
    proseg's ``original_cell_id`` column is populated AND has enough
    overlap with the Xenium ids (thresholds under
    ``nn.hybrid_auto_detect`` in the config).
  * ``coords.*_scale_factor`` + ``coords.verify_same_frame`` — cheap
    frame-mismatch guardrails.
  * ``output.emit_inspection`` — write a sidecar CSV with distance +
    proseg id per Xenium cell (default on; ~free at TMA-scale).

Sentinel-resume matches the other stages: on success, write
``_done.sentinel`` under the output dir. ``force_rerun`` re-runs.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from hexenium._internal.logging import log


# -------------------------------------------------------------------
# Input readers — each returns a normalised (index, x, y[, celltype,
# original_id]) DataFrame so the NN core doesn't have to branch on the
# input file type.
# -------------------------------------------------------------------


def _read_proseg_side(
    source_path: Path,
    *,
    celltype_col: str,
    x_col: str,
    y_col: str,
    id_col: str | None,
    original_id_col: str | None,
) -> pd.DataFrame:
    """Read proseg-side inputs from either an .h5ad or a .csv/.csv.gz.

    Returns a DataFrame with columns
      ``id``: proseg cell identifier (from ``id_col`` if given, else the
        underlying index — cast to str for merge stability).
      ``x``, ``y``: proseg centroid, floats.
      ``celltype``: label from ``celltype_col``.
      ``original_id``: if ``original_id_col`` is given and present,
        the direct-join Xenium UUID that proseg preserved; else NaN.
    """
    suffix = source_path.suffix.lower()
    if suffix == ".h5ad":
        # Lazy import — anndata pulls in h5py which is only needed on this branch.
        import anndata as ad
        adata = ad.read_h5ad(source_path)
        obs = adata.obs
    elif suffix in (".csv", ".gz", ".tsv"):
        sep = "\t" if suffix == ".tsv" else ","
        obs = pd.read_csv(source_path, sep=sep, index_col=0)
    else:
        raise ValueError(
            f"unsupported proseg source_h5ad extension {suffix!r} for {source_path} "
            f"(expected .h5ad, .csv, .csv.gz, or .tsv)"
        )

    for name, col in (("celltype_col", celltype_col), ("x_col", x_col), ("y_col", y_col)):
        if col not in obs.columns:
            raise KeyError(
                f"proseg source {source_path} is missing {name}={col!r}. "
                f"Available columns: {list(obs.columns)}"
            )

    if id_col is not None:
        if id_col not in obs.columns:
            raise KeyError(
                f"proseg source {source_path} is missing id_col={id_col!r}. "
                f"Available columns: {list(obs.columns)}"
            )
        proseg_id = obs[id_col].astype(str).to_numpy()
    else:
        proseg_id = obs.index.astype(str).to_numpy()

    out = pd.DataFrame({
        "id": proseg_id,
        "x": obs[x_col].astype(float).to_numpy(),
        "y": obs[y_col].astype(float).to_numpy(),
        "celltype": obs[celltype_col].astype(object).to_numpy(),
    })
    if original_id_col is not None and original_id_col in obs.columns:
        out["original_id"] = obs[original_id_col].astype(object).to_numpy()
    else:
        out["original_id"] = np.nan

    # Drop rows with NaN x/y — NN sees them as fatal and it's a real
    # data-quality signal worth logging.
    n_pre = len(out)
    out = out.dropna(subset=["x", "y"]).reset_index(drop=True)
    if len(out) != n_pre:
        log(f"[nn_celltype] proseg: dropped {n_pre - len(out)} rows with NaN x/y "
            f"(kept {len(out)})")
    return out


def _read_xenium_side(
    cells_csv_path: Path,
    *,
    id_col: str,
    x_col: str,
    y_col: str,
) -> pd.DataFrame:
    """Read Xenium ``cells.csv.gz`` and normalise to (id, x, y).

    Matches the notebook: ``pd.read_csv(cells.csv.gz, index_col=0)``.
    ``id_col`` names the CSV's index (typically ``cell_id``); if the
    index name doesn't match, we still fall back to the CSV's index
    (matches the notebook's implicit assumption).
    """
    df = pd.read_csv(cells_csv_path, index_col=0)

    for name, col in (("x_col", x_col), ("y_col", y_col)):
        if col not in df.columns:
            raise KeyError(
                f"xenium cells csv {cells_csv_path} is missing {name}={col!r}. "
                f"Available columns: {list(df.columns)}"
            )

    if df.index.name is not None and df.index.name != id_col:
        log(f"[nn_celltype] xenium cells csv index is named "
            f"{df.index.name!r}, config asks for {id_col!r}; using the "
            "CSV's own index (matches notebook's index_col=0 pattern).")

    out = pd.DataFrame({
        "id": df.index.astype(str).to_numpy(),
        "x": df[x_col].astype(float).to_numpy(),
        "y": df[y_col].astype(float).to_numpy(),
    })
    n_pre = len(out)
    out = out.dropna(subset=["x", "y"]).reset_index(drop=True)
    if len(out) != n_pre:
        log(f"[nn_celltype] xenium: dropped {n_pre - len(out)} rows with NaN x/y "
            f"(kept {len(out)})")
    return out


# -------------------------------------------------------------------
# NN core.
# -------------------------------------------------------------------


def _tiebreak_labels(
    labels_kn: np.ndarray,  # (n_query, k), object dtype
    dists_kn: np.ndarray,   # (n_query, k), float
    *,
    tiebreak: str,
) -> np.ndarray:
    """Reduce K neighbours per query to one label per query.

    ``min_dist``: take the first (already sorted by distance ascending
                  in sklearn output).
    ``majority``: mode across the K, break majority ties by nearest.
    """
    if tiebreak == "min_dist":
        return labels_kn[:, 0]
    if tiebreak == "majority":
        out = np.empty(len(labels_kn), dtype=object)
        for i in range(len(labels_kn)):
            counts = Counter(labels_kn[i].tolist())
            top_count = max(counts.values())
            # Preserve neighbour-order (== distance-order) among the
            # majority-tied labels; the first one wins.
            for j in range(labels_kn.shape[1]):
                lab = labels_kn[i, j]
                if counts[lab] == top_count:
                    out[i] = lab
                    break
        return out
    raise ValueError(
        f"unknown nn.tiebreak={tiebreak!r} (expected 'min_dist' or 'majority')"
    )


def _run_nn(
    proseg: pd.DataFrame,
    xenium: pd.DataFrame,
    *,
    k: int,
    algorithm: str,
    metric: str,
    tiebreak: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fit on proseg centroids, query with xenium centroids.

    Returns ``(label_per_xenium, min_dist_per_xenium, proseg_idx_per_xenium)``.
    ``proseg_idx_per_xenium`` is the index into ``proseg`` of the
    tie-break-winning neighbour, useful for the inspection CSV.
    """
    # Lazy import — sklearn is only needed on this branch.
    from sklearn.neighbors import NearestNeighbors

    nbrs = NearestNeighbors(n_neighbors=k, algorithm=algorithm, metric=metric)
    nbrs.fit(proseg[["x", "y"]].to_numpy())

    dists_kn, idx_kn = nbrs.kneighbors(xenium[["x", "y"]].to_numpy())

    proseg_celltype = proseg["celltype"].to_numpy(dtype=object)
    labels_kn = proseg_celltype[idx_kn]

    # For k=1, tiebreak is a no-op; still call so the same code path
    # is exercised regardless of k.
    labels = _tiebreak_labels(labels_kn, dists_kn, tiebreak=tiebreak)

    # For the inspection sidecar: which proseg index actually won.
    if tiebreak == "min_dist" or k == 1:
        winning_idx = idx_kn[:, 0]
        winning_dist = dists_kn[:, 0]
    else:
        winning_idx = np.empty(len(labels), dtype=idx_kn.dtype)
        winning_dist = np.empty(len(labels), dtype=dists_kn.dtype)
        for i in range(len(labels)):
            for j in range(k):
                if labels_kn[i, j] == labels[i]:
                    winning_idx[i] = idx_kn[i, j]
                    winning_dist[i] = dists_kn[i, j]
                    break

    return labels, winning_dist, winning_idx


# -------------------------------------------------------------------
# Hybrid auto-detect: decide whether direct-join is viable per sample.
# -------------------------------------------------------------------


def _should_enable_hybrid(
    proseg: pd.DataFrame,
    xenium: pd.DataFrame,
    *,
    min_original_id_populated: float,
    min_xenium_overlap: float,
) -> bool:
    """Auto-detect whether the direct-join hybrid path is viable.

    Called only when ``hybrid_direct_join_first == "auto"``. Returns
    True iff both:
      (a) ``proseg["original_id"]`` is non-null on at least
          ``min_original_id_populated`` (fraction) of proseg rows.
      (b) The proseg ``original_id`` set covers at least
          ``min_xenium_overlap`` (fraction) of the Xenium cell ids
          (catches the "column present but frames don't overlap" case).

    Logs the decision so it's auditable in run logs.
    """
    populated_mask = proseg["original_id"].notna()
    populated_frac = float(populated_mask.mean()) if len(proseg) else 0.0
    if populated_frac < min_original_id_populated:
        log(
            f"[nn_celltype] hybrid auto-detect: original_id populated="
            f"{populated_frac*100:.1f}%, below "
            f"{min_original_id_populated*100:.0f}% threshold "
            "-> DISABLING hybrid, using pure NN"
        )
        return False

    orig_ids = set(proseg.loc[populated_mask, "original_id"].astype(str).to_numpy())
    xen_hits = xenium["id"].astype(str).isin(orig_ids)
    overlap_frac = float(xen_hits.mean()) if len(xenium) else 0.0
    if overlap_frac < min_xenium_overlap:
        log(
            f"[nn_celltype] hybrid auto-detect: original_id populated="
            f"{populated_frac*100:.1f}%, xenium overlap="
            f"{overlap_frac*100:.1f}% (below "
            f"{min_xenium_overlap*100:.0f}% threshold) "
            "-> DISABLING hybrid, using pure NN"
        )
        return False

    log(
        f"[nn_celltype] hybrid auto-detect: original_id populated="
        f"{populated_frac*100:.1f}%, xenium overlap="
        f"{overlap_frac*100:.1f}%; both above threshold "
        "-> ENABLING hybrid direct-join-first pass"
    )
    return True


# -------------------------------------------------------------------
# Public entrypoint.
# -------------------------------------------------------------------


def _resolve_input_path(template: str | Path, sample_id: str, output_root: Path) -> Path:
    """Resolve a config-provided path.

    Absolute path → used verbatim. Otherwise treated as a template
    (with ``{sample_id}``/``{output_root}`` substitutions) rooted at
    ``<output_root>/<sample_id>/``.
    """
    s = str(template).format(sample_id=sample_id, output_root=str(output_root))
    p = Path(s)
    if not p.is_absolute():
        p = output_root / sample_id / p
    return p


def run_nn_celltype_mapping(
    sample_id: str,
    xenium_bundle: Path,
    output_root: Path,
    *,
    # Proseg side.
    proseg_source: str | Path,
    proseg_celltype_col: str = "first_type",
    proseg_x_col: str = "x",
    proseg_y_col: str = "y",
    proseg_id_col: str | None = None,
    proseg_original_id_col: str | None = "original_cell_id",
    # Xenium side.
    xenium_cells_csv: str = "cells.csv.gz",
    xenium_id_col: str = "cell_id",
    xenium_x_col: str = "x_centroid",
    xenium_y_col: str = "y_centroid",
    # NN algorithm.
    k: int = 1,
    algorithm: str = "auto",
    metric: str = "euclidean",
    tiebreak: str = "min_dist",
    distance_threshold: float | None = None,
    unmatched_policy: str = "mark_unassigned",
    unassigned_label: str = "Unassigned",
    hybrid_direct_join_first: bool | str = "auto",
    hybrid_min_original_id_populated: float = 0.8,
    hybrid_min_xenium_overlap: float = 0.5,
    # Coordinate frame.
    proseg_scale_factor: float = 1.0,
    xenium_scale_factor: float = 1.0,
    verify_same_frame: bool = True,
    frame_check_tolerance: float = 100.0,
    # Output.
    output_dir: str = "celltype_for_hexenium",
    output_stem: str = "{sample_id}_celltype",
    emit_inspection: bool = True,
    inspection_suffix: str = "_inspection",
    write_compression: str | None = None,
    extra_columns: list[str] | None = None,
    # Runtime.
    force_rerun: bool = False,
) -> dict:
    """Map proseg-derived celltypes to Xenium cells via NN.

    Returns a dict:
      ``csv``: path to the two-column hexenium-input CSV.
      ``inspection_csv``: path to the sidecar CSV (if emit_inspection).
      ``sentinel``: sentinel-file path.
      ``n_xenium``: total Xenium cells in.
      ``n_proseg``: total proseg cells in.
      ``n_matched``: rows written to the primary CSV
                     (excludes ``drop`` policy skips).
      ``n_unassigned``: rows given ``unassigned_label`` under
                        ``mark_unassigned`` policy.
      ``n_direct_join``: cells resolved by the direct-join hybrid pass
                         (0 if ``hybrid_direct_join_first`` is off).
    """
    out_dir = output_root / sample_id / output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = output_stem.format(sample_id=sample_id)
    ext = ".csv.gz" if write_compression == "gzip" else ".csv"
    primary_csv = out_dir / f"{stem}{ext}"
    inspection_csv = out_dir / f"{stem}{inspection_suffix}{ext}" if emit_inspection else None
    sentinel = out_dir / "_done.sentinel"

    if sentinel.exists() and not force_rerun:
        log(f"[nn_celltype] skip — sentinel present at {sentinel}")
        return {
            "csv": primary_csv, "inspection_csv": inspection_csv,
            "sentinel": sentinel, "skipped_resume": True,
        }

    # --- Load inputs ------------------------------------------------
    proseg_path = _resolve_input_path(proseg_source, sample_id, output_root)
    if not proseg_path.exists():
        raise FileNotFoundError(f"proseg source not found: {proseg_path}")

    xenium_cells_path = Path(xenium_cells_csv)
    if not xenium_cells_path.is_absolute():
        xenium_cells_path = xenium_bundle / xenium_cells_path
    if not xenium_cells_path.exists():
        raise FileNotFoundError(f"xenium cells csv not found: {xenium_cells_path}")

    log(f"[nn_celltype] proseg source: {proseg_path}")
    log(f"[nn_celltype] xenium cells: {xenium_cells_path}")

    proseg = _read_proseg_side(
        proseg_path,
        celltype_col=proseg_celltype_col,
        x_col=proseg_x_col,
        y_col=proseg_y_col,
        id_col=proseg_id_col,
        original_id_col=proseg_original_id_col,
    )
    xenium = _read_xenium_side(
        xenium_cells_path,
        id_col=xenium_id_col,
        x_col=xenium_x_col,
        y_col=xenium_y_col,
    )

    log(f"[nn_celltype] loaded: n_proseg={len(proseg)} n_xenium={len(xenium)} "
        f"celltype_col={proseg_celltype_col!r}")

    if len(proseg) == 0:
        raise ValueError(f"proseg source {proseg_path} yielded 0 cells after cleanup")
    if len(xenium) == 0:
        raise ValueError(f"xenium cells csv {xenium_cells_path} yielded 0 cells")

    # --- Coordinate rescale + frame check ---------------------------
    if proseg_scale_factor != 1.0:
        proseg[["x", "y"]] *= proseg_scale_factor
        log(f"[nn_celltype] proseg coords scaled by {proseg_scale_factor}")
    if xenium_scale_factor != 1.0:
        xenium[["x", "y"]] *= xenium_scale_factor
        log(f"[nn_celltype] xenium coords scaled by {xenium_scale_factor}")

    if verify_same_frame:
        px_lo, px_hi = float(proseg["x"].min()), float(proseg["x"].max())
        py_lo, py_hi = float(proseg["y"].min()), float(proseg["y"].max())
        xx_lo, xx_hi = float(xenium["x"].min()), float(xenium["x"].max())
        xy_lo, xy_hi = float(xenium["y"].min()), float(xenium["y"].max())
        log(f"[nn_celltype] proseg x=[{px_lo:.1f},{px_hi:.1f}] y=[{py_lo:.1f},{py_hi:.1f}]")
        log(f"[nn_celltype] xenium x=[{xx_lo:.1f},{xx_hi:.1f}] y=[{xy_lo:.1f},{xy_hi:.1f}]")
        gap_x = max(abs(px_lo - xx_lo), abs(px_hi - xx_hi))
        gap_y = max(abs(py_lo - xy_lo), abs(py_hi - xy_hi))
        if gap_x > frame_check_tolerance or gap_y > frame_check_tolerance:
            log(f"[nn_celltype] WARN: proseg / xenium bounding boxes differ by "
                f"gap_x={gap_x:.1f} gap_y={gap_y:.1f} > tol={frame_check_tolerance} — "
                "coordinate-frame mismatch possible. If NN distances look absurd, "
                "check proseg_scale_factor / xenium_scale_factor.")

    # --- Hybrid direct-join first pass (opt-in) ---------------------
    # Resolve "auto" -> True/False before entering the hybrid branch.
    # True/False bypass the detector and behave exactly as before.
    if isinstance(hybrid_direct_join_first, str):
        if hybrid_direct_join_first.lower() == "auto":
            hybrid_direct_join_first = _should_enable_hybrid(
                proseg, xenium,
                min_original_id_populated=hybrid_min_original_id_populated,
                min_xenium_overlap=hybrid_min_xenium_overlap,
            )
        else:
            raise ValueError(
                f"unknown nn.hybrid_direct_join_first="
                f"{hybrid_direct_join_first!r} (expected True, False, or 'auto')"
            )

    n_direct_join = 0
    direct_labels = pd.Series(index=xenium["id"].values, dtype=object)
    direct_proseg_id = pd.Series(index=xenium["id"].values, dtype=object)
    if hybrid_direct_join_first and proseg["original_id"].notna().any():
        # Build a lookup: original_id (xenium uuid) → (celltype, proseg id)
        p_with_orig = proseg[proseg["original_id"].notna()].copy()
        p_with_orig["original_id_str"] = p_with_orig["original_id"].astype(str)
        # Keep first proseg row per original_id (should already be unique).
        p_lookup = p_with_orig.drop_duplicates(subset="original_id_str", keep="first")
        p_lookup = p_lookup.set_index("original_id_str")
        # Intersect with xenium ids.
        hit_mask = xenium["id"].isin(p_lookup.index)
        n_direct_join = int(hit_mask.sum())
        if n_direct_join:
            hit_ids = xenium.loc[hit_mask, "id"].to_numpy()
            hit_labels = p_lookup.loc[hit_ids, "celltype"].to_numpy(dtype=object)
            hit_proseg = p_lookup.loc[hit_ids, "id"].to_numpy(dtype=object)
            direct_labels.loc[hit_ids] = hit_labels
            direct_proseg_id.loc[hit_ids] = hit_proseg
        log(f"[nn_celltype] hybrid direct-join first pass: matched {n_direct_join} "
            f"of {len(xenium)} xenium cells via original_id "
            f"({100.0*n_direct_join/len(xenium):.1f}%)")

    # --- NN for the rest (or all, if hybrid off) --------------------
    if hybrid_direct_join_first and n_direct_join > 0:
        need_nn = xenium[direct_labels.reindex(xenium["id"]).isna().to_numpy()].reset_index(drop=True)
    else:
        need_nn = xenium

    if len(need_nn) > 0:
        if k > len(proseg):
            raise ValueError(
                f"nn.k={k} exceeds proseg cell count {len(proseg)}"
            )
        nn_labels, nn_dist, nn_idx = _run_nn(
            proseg, need_nn,
            k=k, algorithm=algorithm, metric=metric, tiebreak=tiebreak,
        )
        d = nn_dist
        log(f"[nn_celltype] nn distances (n={len(d)}): "
            f"min={d.min():.3f} median={np.median(d):.3f} "
            f"p95={np.percentile(d, 95):.3f} max={d.max():.3f}")
        nn_result = pd.DataFrame({
            "id": need_nn["id"].to_numpy(),
            "celltype": nn_labels,
            "nn_dist": nn_dist,
            "proseg_idx": nn_idx,
        })
    else:
        nn_result = pd.DataFrame(
            {"id": [], "celltype": [], "nn_dist": [], "proseg_idx": []}
        )

    # Add proseg_id from the winning proseg_idx.
    if len(nn_result):
        proseg_ids_arr = proseg["id"].to_numpy()
        nn_result["proseg_id"] = proseg_ids_arr[nn_result["proseg_idx"].to_numpy()]
    else:
        nn_result["proseg_id"] = pd.Series(dtype=object)

    # --- Distance threshold policy ----------------------------------
    n_unassigned = 0
    assignment_note = pd.Series("", index=nn_result.index, dtype=object)
    keep_mask_nn = pd.Series(True, index=nn_result.index)

    if distance_threshold is not None and len(nn_result):
        far_mask = nn_result["nn_dist"].to_numpy() > float(distance_threshold)
        n_far = int(far_mask.sum())
        if n_far:
            if unmatched_policy == "drop":
                keep_mask_nn = pd.Series(~far_mask, index=nn_result.index)
                log(f"[nn_celltype] distance_threshold={distance_threshold}: "
                    f"dropped {n_far} rows (unmatched_policy=drop)")
            elif unmatched_policy == "mark_unassigned":
                nn_result.loc[far_mask, "celltype"] = unassigned_label
                assignment_note.loc[far_mask] = "beyond_threshold_unassigned"
                n_unassigned = n_far
                log(f"[nn_celltype] distance_threshold={distance_threshold}: "
                    f"marked {n_far} rows as {unassigned_label!r} (unmatched_policy=mark_unassigned)")
            elif unmatched_policy == "keep":
                assignment_note.loc[far_mask] = "beyond_threshold_kept"
                log(f"[nn_celltype] distance_threshold={distance_threshold}: "
                    f"{n_far} rows are beyond threshold, kept anyway (unmatched_policy=keep)")
            else:
                raise ValueError(
                    f"unknown nn.unmatched_policy={unmatched_policy!r} "
                    "(expected 'drop', 'mark_unassigned', or 'keep')"
                )

    # --- Assemble final tables --------------------------------------
    # Primary output: (cell_id, group) — the hexenium input shape.
    rows = []
    if hybrid_direct_join_first:
        # Direct-join cells first.
        hit_ids = direct_labels.dropna().index.to_numpy()
        for xid in hit_ids:
            rows.append({
                "cell_id": xid,
                "group": direct_labels.loc[xid],
                "proseg_id": direct_proseg_id.loc[xid],
                "nn_dist": 0.0,
                "assignment_note": "direct_join",
                "keep": True,
            })
    for _, r in nn_result.iterrows():
        rows.append({
            "cell_id": r["id"],
            "group": r["celltype"],
            "proseg_id": r["proseg_id"],
            "nn_dist": float(r["nn_dist"]),
            "assignment_note": assignment_note.loc[r.name] or ("nn" if not hybrid_direct_join_first else "nn_fallback"),
            "keep": bool(keep_mask_nn.loc[r.name]),
        })
    combined = pd.DataFrame(rows)

    if len(combined) == 0:
        raise RuntimeError(
            "nn_celltype_mapping produced 0 rows — check inputs (proseg/xenium)."
        )

    # Deduplicate on cell_id — a Xenium cell id shouldn't appear
    # twice; if it does (e.g. a duplicated cells.csv.gz), keep first.
    n_before_dedup = len(combined)
    combined = combined.drop_duplicates(subset="cell_id", keep="first").reset_index(drop=True)
    if len(combined) != n_before_dedup:
        log(f"[nn_celltype] dropped {n_before_dedup - len(combined)} duplicate "
            f"cell_id rows (kept first)")

    kept = combined[combined["keep"]].reset_index(drop=True)
    n_matched = len(kept)

    # Primary CSV (only cell_id + group + any explicitly requested extras).
    primary_cols = ["cell_id", "group"]
    extra_columns = extra_columns or []
    for c in extra_columns:
        if c not in {"nn_dist", "proseg_id", "assignment_note"}:
            log(f"[nn_celltype] WARN: extra column {c!r} not recognised, skipping")
            continue
        if c not in primary_cols:
            primary_cols.append(c)
    kept_out = kept[primary_cols]
    kept_out.to_csv(
        primary_csv, index=False,
        compression="gzip" if write_compression == "gzip" else None,
    )
    log(f"[nn_celltype] wrote {len(kept_out)} rows -> {primary_csv.name} "
        f"(cols: {list(kept_out.columns)})")

    # Inspection sidecar (always writes the full six-column shape when on).
    if emit_inspection:
        insp_cols = ["cell_id", "group", "proseg_id", "nn_dist", "assignment_note"]
        combined[insp_cols].to_csv(
            inspection_csv, index=False,
            compression="gzip" if write_compression == "gzip" else None,
        )
        log(f"[nn_celltype] wrote inspection sidecar -> {inspection_csv.name}")

    # --- Sentinel + return ------------------------------------------
    sentinel.write_text("ok\n")
    log(f"[nn_celltype] sentinel -> {sentinel.name}")

    label_counts = kept["group"].value_counts()
    log(f"[nn_celltype] label distribution (top 10):")
    for lab, n in label_counts.head(10).items():
        log(f"[nn_celltype]   {lab}: {n}")

    return {
        "csv": primary_csv,
        "inspection_csv": inspection_csv,
        "sentinel": sentinel,
        "n_xenium": len(xenium),
        "n_proseg": len(proseg),
        "n_matched": n_matched,
        "n_unassigned": n_unassigned,
        "n_direct_join": n_direct_join,
    }
