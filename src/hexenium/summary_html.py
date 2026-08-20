"""HTML summary for the promoted default run.

Renders ``<sample>/he_registration/output/summary.html`` after every
successful promotion (both ``hexenium run --set-default-on-success``
and ``hexenium set-default-run``). Reads what ``output/<stage>``
currently points at, pulls each stage's per-run ``manifest.yaml``,
detects two flavors of warning (missing-file, param-drift), and
emits a single **self-contained** HTML file — every ``<img>`` on
the page is an inline ``data:image/png;base64,…`` URI, downscaled
via :func:`_image_to_data_uri` to keep the whole document under a
sane size cap. That means the HTML can be moved, emailed, or
attached to a report without losing the figures, which is the
promise Tracy asked for on
``settylab/TracyY123-nexus#15`` comment ``5346844112``.

Design references (all on ``settylab/TracyY123-nexus#15``):

* Comment ``5334834202`` — Tracy's original ``output/`` +
  ``summary.html`` design proposal.
* Comment ``5334874062`` — the interpretation reply with
  clarifying questions.
* Comment ``5334969034`` — Tracy's answers (both warnings kept
  with distinct labels; overlap-quality metric deferred).
* Comment ``5337905668`` — green-light for the standard-content
  renderer.
* Comment ``5338129631`` — green-light for A3 params sections +
  C2 dynamic overlap-figure discovery.
* Comment ``5346844112`` — the reversal of the earlier
  linked-thumbnails choice; every figure now embeds inline.

**Atomicity**: the renderer writes via ``tmp + os.replace`` so a
mid-write crash never leaves a half-serialised HTML behind. The
symlinks under ``output/`` are the source of truth; render
failures LOG loudly and re-raise but do NOT touch the symlinks
(that would violate the earlier atomicity contract from
``a92cd8e``).
"""
from __future__ import annotations

import base64
import html
import io
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from hexenium._internal.logging import log
from hexenium.manifest import (
    compute_params_hash,
    output_symlink_run_id,
    read_manifest,
)


#: Max side length in px for embedded overlap / viz thumbnails.
#: VALIS overlap PNGs at full-res are 10-30MB each; embedding four of
#: them raw would produce 100+MB HTML that crashes browsers. 800px on
#: the longer edge is a QA-legible size at ~200KB per PNG post-encode
#: (rough — depends on content complexity). Tracy asked for self-contained
#: HTML on settylab/TracyY123-nexus#15 comment 5346844112; this cap is
#: what makes that promise practical.
_EMBED_MAX_DIM_PX = 800


#: Stage → dirname mapping; matches ``pipeline._PROMOTE_STAGE_DIR`` but
#: kept independent so this module has no pipeline import (avoids
#: import cycles and keeps the renderer usable from set-default-run
#: without pulling the pipeline entry point).
STAGE_DIRS: tuple[tuple[str, str], ...] = (
    ("register",  "register"),
    ("warp",      "warp"),
    ("celltype",  "celltyped"),
    ("viz",       "viz"),
)


#: Per-stage artifacts that MUST exist under ``<stage>/<he_job_id>/`` for
#: the run to be considered whole. Missing → ``missing-file`` warning.
#: Kept as a set so ordering doesn't affect log stability.
_EXPECTED_ARTIFACTS: dict[str, tuple[str, ...]] = {
    "register":  ("data/_registrar.pickle",),
    "warp":      ("he_cell_seg.parquet", "he_nucleus_seg.parquet"),
    "celltype":  (
        "{sample}_cells_analysis.geojson",
        "{sample}_cells_qupath.geojson",
        "{sample}_nuclei_analysis.geojson",
        "{sample}_nuclei_qupath.geojson",
        "{sample}_celltyped_wholeslide.parquet",
    ),
    "viz":       ("{sample}_overlay.png",),
}


#: Register "highlight" params — the ones a scanner-of-summaries cares
#: about at a glance. Everything ELSE in ``manifest["params"]`` still
#: lands in the collapsible full-dump ``<details>`` block. Priority
#: order = display order.
_REGISTER_HIGHLIGHT_KEYS: tuple[str, ...] = (
    "mode",
    "use_he_deconvolution",
    "check_for_reflections",
    "align_to_reference",
)

#: Warp highlights — the load-bearing ones are the target list + dask
#: switch. save_geojson toggles a bulky per-target sidecar that some
#: downstream QuPath workflows depend on, so it earns a highlight slot.
_WARP_HIGHLIGHT_KEYS: tuple[str, ...] = (
    "targets",
    "use_dask",
    "save_geojson",
)

#: VALIS overlap-diagnostic filenames in REFINEMENT order — used to
#: give the discovered overlaps a stable, semantically-meaningful
#: display order. Filenames NOT in this list fall to alphabetical
#: order after these (dynamic discovery is the source of truth for
#: which images render, per Tracy's constraint in
#: ``settylab/TracyY123-nexus#15`` comment ``5338129631``).
_OVERLAP_PRIORITY: tuple[str, ...] = (
    "_original_overlap.png",
    "_rigid_overlap.png",
    "_non_rigid_overlap.png",
    "_micro_reg.png",
)

#: Human-readable labels for the known VALIS overlap filenames. Files
#: not in this map render with their bare filename as the label.
_OVERLAP_LABELS: dict[str, str] = {
    "_original_overlap.png":  "Original (before registration)",
    "_rigid_overlap.png":     "After rigid solve",
    "_non_rigid_overlap.png": "After non-rigid solve",
    "_micro_reg.png":         "After micro solve",
}


@dataclass
class _StageInfo:
    """Everything the renderer needs about one promoted stage.

    ``he_job_id`` = ``None`` when ``output/<stage>`` is missing (a
    partial promotion — celltype/viz may be un-promoted while
    register/warp are set). ``manifest`` = ``None`` when the linked
    ``<stage>/<he_job_id>/manifest.yaml`` doesn't exist (celltype
    and viz stages don't currently write one — that's expected).
    """
    stage: str
    stage_dir_name: str            # "celltyped" for stage "celltype", etc.
    he_job_id: str | None
    stage_dir: Path | None
    manifest: dict[str, Any] | None
    missing_artifacts: list[str] = field(default_factory=list)
    param_drift: bool = False


def render_summary_html(
    *,
    output_root_he: Path,
    sample_id: str,
    run_id: str | None = None,
) -> Path:
    """Render + atomically write ``output/summary.html``.

    ``output_root_he`` is the ``he_registration/`` root (either
    ``<output_root>/<sample>/`` standalone, or
    ``<xenium_run_dir>/he_registration/`` integrated). Reads the
    live ``output/<stage>`` symlinks to determine what's promoted;
    per-run manifests are read from wherever each symlink resolves
    to.

    Returns the written path (``<output_root_he>/output/summary.html``).
    Raises the underlying ``OSError`` if the write itself fails; the
    caller is responsible for logging LOUDLY and NOT rolling back
    the symlinks that were already committed atomically by
    ``set_default_run``.
    """
    output_dir = output_root_he / "output"
    stages = [
        _collect_stage_info(
            stage=stage, dir_name=dir_name,
            output_root_he=output_root_he, output_dir=output_dir,
            sample_id=sample_id,
        )
        for stage, dir_name in STAGE_DIRS
    ]
    html_body = _render(sample_id=sample_id, run_id=run_id, stages=stages)
    dest = output_dir / "summary.html"
    _atomic_write_text(dest, html_body)
    log(f"[summary-html] wrote {dest}")
    return dest


# ---------------------------------------------------------------------
# Per-stage collectors
# ---------------------------------------------------------------------
def _collect_stage_info(
    *,
    stage: str,
    dir_name: str,
    output_root_he: Path,
    output_dir: Path,
    sample_id: str,
) -> _StageInfo:
    """Build a ``_StageInfo`` for one stage from what's on disk NOW."""
    info = _StageInfo(stage=stage, stage_dir_name=dir_name,
                      he_job_id=None, stage_dir=None, manifest=None)

    link = output_dir / dir_name
    he_job_id = output_symlink_run_id(link)
    info.he_job_id = he_job_id
    if he_job_id is None:
        return info

    stage_dir = output_root_he / dir_name / he_job_id
    info.stage_dir = stage_dir

    manifest_path = stage_dir / "manifest.yaml"
    if manifest_path.exists():
        try:
            info.manifest = read_manifest(manifest_path)
        except Exception as exc:  # noqa: BLE001 — best-effort read
            log(f"[summary-html] WARN: failed to read {manifest_path}: {exc!r}")

    info.missing_artifacts = _check_missing_artifacts(
        stage=stage, stage_dir=stage_dir, sample_id=sample_id,
    )
    info.param_drift = _check_param_drift(manifest=info.manifest)
    return info


def _check_missing_artifacts(
    *, stage: str, stage_dir: Path, sample_id: str,
) -> list[str]:
    """Return the list of expected artifacts missing under ``stage_dir``.

    Absence of the stage-dir itself → all artifacts marked missing.
    ``{sample}`` placeholder is substituted per stage from the
    caller's ``sample_id``.
    """
    if not stage_dir.is_dir():
        return list(_EXPECTED_ARTIFACTS.get(stage, ()))
    missing: list[str] = []
    for tmpl in _EXPECTED_ARTIFACTS.get(stage, ()):
        rel = tmpl.format(sample=sample_id)
        if not (stage_dir / rel).exists():
            missing.append(rel)
    return missing


def _check_param_drift(*, manifest: dict[str, Any] | None) -> bool:
    """Self-consistency check on the stored ``params_hash``.

    Recomputes ``compute_params_hash(manifest["params"])`` and
    compares to the stored ``manifest["params_hash"]``. A mismatch
    means the manifest was manually edited (or the hashing algorithm
    changed under it) — flag it. Missing either field → no drift
    (can't check).
    """
    if manifest is None:
        return False
    stored = manifest.get("params_hash")
    params = manifest.get("params")
    if stored is None or params is None:
        return False
    return compute_params_hash(params) != stored


# ---------------------------------------------------------------------
# HTML rendering — stdlib only, no Jinja2.
# ---------------------------------------------------------------------
def _render(*, sample_id: str, run_id: str | None,
            stages: list[_StageInfo]) -> str:
    """Compose the final HTML document from per-stage info."""
    ts = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    title = html.escape(f"hexenium — {sample_id}"
                        + (f" / {run_id}" if run_id else ""))
    pointer_line = _format_output_pointer(stages)
    per_stage_rows = "\n".join(_render_stage_row(s) for s in stages)
    lineage_html = _render_lineage(stages)
    thumb_html = _render_thumbnails(sample_id=sample_id, stages=stages)
    stage_by_name = {s.stage: s for s in stages}
    register_params_html = _render_stage_params(
        stage_by_name.get("register"),
        heading="Registration parameters",
        highlight_keys=_REGISTER_HIGHLIGHT_KEYS,
        resolution_keys=(
            "max_image_dim_px",
            "max_processed_image_dim_px",
            "max_non_rigid_registration_dim_px",
        ),
    )
    warp_params_html = _render_stage_params(
        stage_by_name.get("warp"),
        heading="Warp parameters",
        highlight_keys=_WARP_HIGHLIGHT_KEYS,
        resolution_keys=(),
    )
    overlap_html = _render_overlap_figures(stage_by_name.get("register"))

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI",
                       sans-serif; margin: 2em; max-width: 1100px; }}
  h1, h2 {{ margin-bottom: 0.2em; }}
  .meta {{ color: #666; font-size: 0.9em; margin: 0.2em 0 1em 0; }}
  table {{ border-collapse: collapse; margin: 0.5em 0 1.5em 0;
           font-size: 0.92em; }}
  th, td {{ border: 1px solid #ccc; padding: 0.4em 0.7em;
            text-align: left; vertical-align: top; }}
  th {{ background: #f0f0f0; }}
  code {{ background: #f6f8fa; padding: 0.05em 0.35em; border-radius: 3px;
           font-size: 0.9em; }}
  .warn {{ color: #b71c1c; font-weight: bold; }}
  .warn-drift {{ color: #b71c1c; font-weight: bold;
                 background: #ffebee; padding: 0.1em 0.4em; }}
  .warn-missing {{ color: #b71c1c; font-weight: bold;
                   background: #fff3e0; padding: 0.1em 0.4em; }}
  .ok {{ color: #2e7d32; }}
  .none {{ color: #999; font-style: italic; }}
  .thumb {{ max-width: 640px; max-height: 480px; border: 1px solid #ccc;
            display: block; margin: 0.3em 0 1em 0; }}
  ul {{ margin: 0.3em 0 1em 1.5em; padding: 0; }}
  .chips {{ margin: 0.3em 0 0.5em 0; line-height: 1.9; }}
  .chip {{ display: inline-block; background: #eef4fb; border: 1px solid #cfd8dc;
           border-radius: 3px; padding: 0.15em 0.5em; margin-right: 0.4em;
           font-size: 0.9em; }}
  .params-dump {{ font-size: 0.88em; margin: 0.3em 0 0.8em 0; }}
  .params-dump th {{ text-align: right; background: #fafafa;
                     min-width: 12em; }}
  details {{ margin: 0.3em 0 1em 0; }}
  details > summary {{ cursor: pointer; color: #1565c0;
                       font-size: 0.9em; padding: 0.2em 0; }}
  .overlap-grid {{ display: grid; grid-template-columns: repeat(auto-fill,
                     minmax(320px, 1fr)); gap: 1em; margin: 0.3em 0 1em 0; }}
  .overlap-tile {{ margin: 0; }}
  .overlap-tile .thumb.overlap {{ max-width: 100%; max-height: 320px;
                                  margin-bottom: 0.2em; }}
  .overlap-tile .meta {{ margin: 0; }}
</style>
</head>
<body>
<h1>{title}</h1>
<p class="meta">Rendered {html.escape(ts)}</p>
<p class="meta">{pointer_line}</p>

<h2>Per-stage summary</h2>
<table>
  <thead><tr>
    <th>Stage</th><th>he_job_id</th><th>params_hash</th>
    <th>Timestamp (UTC)</th><th>git_sha</th><th>Flags</th>
  </tr></thead>
  <tbody>
{per_stage_rows}
  </tbody>
</table>

<h2>Lineage</h2>
{lineage_html}

{register_params_html}

{warp_params_html}

{overlap_html}

<h2>Overlay</h2>
{thumb_html}
</body>
</html>
"""


def _format_output_pointer(stages: list[_StageInfo]) -> str:
    """One-line summary of what ``output/`` currently points at."""
    parts = []
    for s in stages:
        if s.he_job_id is None:
            parts.append(f"{html.escape(s.stage_dir_name)}: <span class='none'>(unset)</span>")
        else:
            parts.append(
                f"{html.escape(s.stage_dir_name)}: "
                f"<code>{html.escape(s.he_job_id)}</code>"
            )
    return "output/ &rarr; " + " &middot; ".join(parts)


def _render_stage_row(s: _StageInfo) -> str:
    """One table row per stage, warnings folded into the Flags column."""
    if s.he_job_id is None:
        return (
            "    <tr>"
            f"<td>{html.escape(s.stage_dir_name)}</td>"
            "<td colspan='5' class='none'>output symlink not set</td>"
            "</tr>"
        )

    m = s.manifest or {}
    params_hash = m.get("params_hash") or ""
    timestamp = m.get("timestamp_utc") or ""
    git_sha = m.get("git_sha") or ""
    flags: list[str] = []
    for missing in s.missing_artifacts:
        flags.append(
            "<span class='warn-missing' title='expected artifact absent under "
            "the linked stage folder'>missing-file: "
            f"{html.escape(missing)}</span>"
        )
    if s.param_drift:
        flags.append(
            "<span class='warn-drift' title='stored params_hash disagrees with "
            "compute_params_hash(manifest.params)'>param-drift</span>"
        )
    if not flags and not m:
        flags.append("<span class='none'>no manifest</span>")
    elif not flags:
        flags.append("<span class='ok'>ok</span>")

    return (
        "    <tr>"
        f"<td>{html.escape(s.stage_dir_name)}</td>"
        f"<td><code>{html.escape(s.he_job_id)}</code></td>"
        f"<td><code>{html.escape(str(params_hash))}</code></td>"
        f"<td>{html.escape(str(timestamp))}</td>"
        f"<td><code>{html.escape(str(git_sha)[:12])}</code></td>"
        f"<td>{' '.join(flags)}</td>"
        "</tr>"
    )


def _render_lineage(stages: list[_StageInfo]) -> str:
    """Render the register↔warp lineage line, plus a mismatch flag.

    Reads the warp stage's ``source_register_run_id`` from its
    manifest and compares against what ``output/register`` currently
    points at. Reports OK or mismatch inline.
    """
    stage_by_name = {s.stage: s for s in stages}
    register = stage_by_name.get("register")
    warp = stage_by_name.get("warp")
    if register is None or warp is None or warp.manifest is None:
        return "<p class='none'>register↔warp lineage unavailable (missing manifest).</p>"

    source = warp.manifest.get("source_register_run_id")
    r_id = register.he_job_id
    items: list[str] = []
    items.append(
        f"warp/{html.escape(warp.he_job_id or '?')}/manifest.yaml recorded "
        f"<code>source_register_run_id</code>: "
        f"<code>{html.escape(str(source) if source else '(none)')}</code>"
    )
    items.append(
        f"output/register &rarr; register/<code>{html.escape(str(r_id) if r_id else '(unset)')}</code>"
    )
    if source and r_id and source == r_id:
        items.append("<span class='ok'>lineage OK</span>")
    elif source and r_id:
        items.append(
            "<span class='warn'>lineage MISMATCH — "
            "warp's source_register_run_id does not match output/register</span>"
        )
    return "<ul>\n" + "\n".join(f"  <li>{it}</li>" for it in items) + "\n</ul>"


def _render_stage_params(
    info: _StageInfo | None,
    *,
    heading: str,
    highlight_keys: tuple[str, ...],
    resolution_keys: tuple[str, ...],
) -> str:
    """Render one stage's params as ``<highlights row>`` + ``<details>``
    with a full key/value dump.

    A3 shape (per Tracy's spec on
    ``settylab/TracyY123-nexus#15`` comment ``5338129631``):
    a compact highlights line for the load-bearing keys the operator
    scans at a glance, plus a collapsible ``<details>`` with every
    remaining key from ``manifest["params"]`` so nothing is silently
    hidden.

    ``resolution_keys`` (when non-empty) fold into a single
    ``resolution: <k1>/<k2>/<k3>`` one-liner in the highlights,
    because "1500/1500/10000" reads better than three separate rows.
    """
    if info is None or info.he_job_id is None or info.manifest is None:
        return (
            f"<h2>{html.escape(heading)}</h2>\n"
            "<p class='none'>(no manifest — stage may not have been "
            "promoted or was written by an older pipeline version)</p>"
        )
    params = info.manifest.get("params") or {}
    if not params:
        return (
            f"<h2>{html.escape(heading)}</h2>\n"
            "<p class='none'>(manifest present but no params block)</p>"
        )

    highlight_html = _render_highlight_line(
        params, highlight_keys, resolution_keys,
    )
    full_dump_html = _render_full_params_dump(params)

    return (
        f"<h2>{html.escape(heading)}</h2>\n"
        f"{highlight_html}\n"
        "<details>\n"
        "  <summary>Full manifest params</summary>\n"
        f"{full_dump_html}\n"
        "</details>"
    )


def _render_highlight_line(
    params: dict[str, Any],
    highlight_keys: tuple[str, ...],
    resolution_keys: tuple[str, ...],
) -> str:
    """One-line dense summary of the highlight-tier params.

    Missing keys are silently omitted so a run that DIDN'T set a
    given knob doesn't fabricate a fake value. Resolution keys are
    folded into a single ``resolution`` chip because their numeric
    combination is what an operator eyeballs.
    """
    chips: list[str] = []
    for k in highlight_keys:
        if k not in params:
            continue
        v = params[k]
        chips.append(
            f"<span class='chip'><code>{html.escape(k)}</code> = "
            f"<code>{html.escape(_stringify_param_value(v))}</code></span>"
        )
    if resolution_keys:
        present = [k for k in resolution_keys if k in params]
        if present:
            joined = "/".join(str(params[k]) for k in present)
            chips.append(
                "<span class='chip'><code>resolution</code> = "
                f"<code>{html.escape(joined)}</code> "
                "<span class='meta'>("
                + ", ".join(html.escape(k) for k in present) +
                ")</span></span>"
            )
    if not chips:
        return "<p class='none'>(no highlight params matched this stage)</p>"
    return "<p class='chips'>" + " ".join(chips) + "</p>"


def _render_full_params_dump(params: dict[str, Any]) -> str:
    """Two-column table with EVERY key/value from ``params``.

    Sorted by key so the dump is stable across re-renders and
    diffable. Lists / dicts get repr'd via ``_stringify_param_value``.
    """
    rows: list[str] = []
    for k in sorted(params.keys()):
        v = params[k]
        rows.append(
            "    <tr>"
            f"<th><code>{html.escape(k)}</code></th>"
            f"<td><code>{html.escape(_stringify_param_value(v))}</code></td>"
            "</tr>"
        )
    return (
        "  <table class='params-dump'>\n"
        + "\n".join(rows)
        + "\n  </table>"
    )


def _stringify_param_value(v: Any) -> str:
    """Compact string form of a param value for a table cell.

    ``bool``/``int``/``float``/``str`` render as-is; lists get a
    compact ``[a, b, c]`` form; anything else falls to ``repr``.
    """
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float, str)):
        return str(v)
    if isinstance(v, list):
        return "[" + ", ".join(_stringify_param_value(x) for x in v) + "]"
    return repr(v)


def _render_overlap_figures(register: _StageInfo | None) -> str:
    """Discover + link ``register/<he_job_id>/overlaps/*.png`` dynamically.

    Per Tracy's C2 spec: NEVER hardcode the list of overlap files —
    walk the directory at render time and include whichever PNGs
    exist. Files in :data:`_OVERLAP_PRIORITY` render in that
    refinement order first; anything else falls to alphabetical
    order after. Empty state renders an explicit
    "no overlap diagnostics" line rather than a broken/blank block.
    """
    if register is None or register.he_job_id is None or register.stage_dir is None:
        return (
            "<h2>Registration overlap figures</h2>\n"
            "<p class='none'>(no promoted register run &mdash; "
            "no overlaps to display)</p>"
        )
    overlaps_dir = register.stage_dir / "overlaps"
    if not overlaps_dir.is_dir():
        return (
            "<h2>Registration overlap figures</h2>\n"
            f"<p class='none'>No overlap diagnostics on disk under "
            f"<code>{html.escape(str(overlaps_dir.name))}/</code>.</p>"
        )
    pngs = sorted(overlaps_dir.glob("*.png"), key=_overlap_sort_key)
    if not pngs:
        return (
            "<h2>Registration overlap figures</h2>\n"
            "<p class='none'>No overlap diagnostics on disk "
            "(register may have failed before writing them).</p>"
        )

    tiles: list[str] = []
    for png in pngs:
        label = _OVERLAP_LABELS.get(png.name, png.name)
        data_uri = _image_to_data_uri(png)
        if not data_uri:
            # Graceful degradation: PIL failed on this specific PNG.
            # Emit a placeholder tile instead of a broken <img> so the
            # rest of the section still reads cleanly.
            tiles.append(
                "<div class='overlap-tile'>\n"
                "  <p class='none'>(could not embed "
                f"<code>{html.escape(png.name)}</code> &mdash; see log)</p>\n"
                "</div>"
            )
            continue
        tiles.append(
            "<div class='overlap-tile'>\n"
            f"  <img class='thumb overlap' src='{data_uri}' "
            f"alt='{html.escape(label)}'>\n"
            f"  <p class='meta'>{html.escape(label)} &mdash; "
            f"<code>{html.escape(png.name)}</code></p>\n"
            "</div>"
        )
    return (
        "<h2>Registration overlap figures</h2>\n"
        "<div class='overlap-grid'>\n"
        + "\n".join(tiles) +
        "\n</div>"
    )


def _overlap_sort_key(png_path: Path) -> tuple[int, str]:
    """Sort key: (priority index, filename). Priority-listed files
    first in the documented refinement order; unknowns alphabetical
    after (all get ``len(_OVERLAP_PRIORITY)`` as index so ``sorted``
    breaks ties by filename)."""
    name = png_path.name
    try:
        return (_OVERLAP_PRIORITY.index(name), name)
    except ValueError:
        return (len(_OVERLAP_PRIORITY), name)


def _render_thumbnails(*, sample_id: str, stages: list[_StageInfo]) -> str:
    """Render the viz overlay as the single top-level thumbnail.

    Register / warp / celltyped don't currently produce a canonical
    single image, so we only thumbnail viz. If viz isn't promoted
    (or its overlay is missing), report so; the ``missing-file``
    flag on the viz row will name it.
    """
    stage_by_name = {s.stage: s for s in stages}
    v = stage_by_name.get("viz")
    if v is None or v.he_job_id is None or v.stage_dir is None:
        return (
            "<p class='none'>viz stage was not run in this invocation.</p>"
        )
    overlay_path = v.stage_dir / f"{sample_id}_overlay.png"
    if not overlay_path.exists():
        return (
            f"<p class='none'>viz overlay expected at "
            f"<code>{html.escape(overlay_path.name)}</code> but the file is "
            "missing on disk &mdash; see the missing-file flag on the "
            "per-stage summary above.</p>"
        )
    data_uri = _image_to_data_uri(overlay_path)
    if not data_uri:
        return (
            "<p class='none'>viz overlay present on disk but could not be "
            "embedded &mdash; see log for details.</p>"
        )
    return f"<img class='thumb' src='{data_uri}' alt='overlay'>"


# ---------------------------------------------------------------------
# Atomic write
# ---------------------------------------------------------------------
def _image_to_data_uri(path: Path, *, max_dim: int = _EMBED_MAX_DIM_PX) -> str:
    """Downscale + re-encode ``path`` as a ``data:image/png;base64,…`` URI.

    Tracy's self-contained-HTML ask (comment ``5346844112``) means we
    embed every figure inline. The downscale keeps output size sane —
    a full-res VALIS overlap PNG can be 20+ MB; at 800px longer edge
    it's typically <300 KB after PNG re-encode.

    Returns ``""`` on any read/decode failure — callers should detect
    the empty result and either skip the image or emit a placeholder.
    Logs the pre/post byte sizes so a future compression regression is
    visible in the pipeline log without a rerun. All exceptions from
    PIL are caught and surfaced as WARN + empty string, so a corrupted
    or truncated PNG can't take down the whole HTML render — the rest
    of the summary still lands.
    """
    try:
        from PIL import Image  # imported lazily so tests that don't
                               # touch images don't need PIL loaded.
    except ImportError as exc:
        log(f"[summary-html] WARN: PIL not available; cannot embed "
            f"{path}: {exc!r}")
        return ""
    try:
        raw_size = path.stat().st_size
        with Image.open(path) as img:
            img.load()
            w, h = img.size
            if max(w, h) > max_dim:
                scale = max_dim / float(max(w, h))
                # Round to int, keep aspect ratio.
                new_size = (max(1, int(w * scale)), max(1, int(h * scale)))
                img = img.resize(new_size, Image.LANCZOS)
            # Force-drop alpha for PNG optim; VALIS overlaps are RGB
            # composites but be defensive against RGBA/P modes.
            if img.mode not in ("RGB", "L"):
                img = img.convert("RGB")
            buf = io.BytesIO()
            img.save(buf, format="PNG", optimize=True)
            encoded = base64.b64encode(buf.getvalue()).decode("ascii")
        log(f"[summary-html] embedded {path.name}: "
            f"{raw_size / 1024:.0f} KiB on disk -> "
            f"{len(buf.getvalue()) / 1024:.0f} KiB embedded "
            f"({img.size[0]}x{img.size[1]}px after downscale)")
        return f"data:image/png;base64,{encoded}"
    except Exception as exc:  # noqa: BLE001 — one image failing must
        # not sink the whole HTML render; degrade gracefully.
        log(f"[summary-html] WARN: could not embed {path}: {exc!r}")
        return ""


def _atomic_write_text(path: Path, text: str) -> None:
    """Write ``text`` to ``path`` via tmp + ``os.replace``.

    Mirrors ``manifest.write_manifest``'s atomicity guarantee: a
    crash mid-write never leaves a half-serialised HTML that would
    be misread by a downstream viewer / grep.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with open(tmp, "w") as f:
        f.write(text)
    os.replace(tmp, path)
