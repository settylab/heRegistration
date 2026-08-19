"""HTML summary for the promoted default run.

Renders ``<sample>/he_registration/output/summary.html`` after every
successful promotion (both ``hexenium run --set-default-on-success``
and ``hexenium set-default-run``). Reads what ``output/<stage>``
currently points at, pulls each stage's per-run ``manifest.yaml``,
detects two flavors of warning (missing-file, param-drift), and
emits a single self-contained HTML file. Thumbnails LINK to
``<stage>/<he_job_id>/`` outputs (relative refs) rather than
embedding — moving the HTML out of the tree breaks the thumbnails,
which is the desired signal that the tree owns the artifacts.

Design references (all on ``settylab/TracyY123-nexus#15``):

* Comment ``5334834202`` — Tracy's original ``output/`` +
  ``summary.html`` design proposal.
* Comment ``5334874062`` — the interpretation reply with
  clarifying questions.
* Comment ``5334969034`` — Tracy's answers (both warnings kept
  with distinct labels; overlap-quality metric deferred; thumbnails).
* Comment ``5337905668`` — the green-light for this iteration.

**Atomicity**: the renderer writes via ``tmp + os.replace`` so a
mid-write crash never leaves a half-serialised HTML behind. The
symlinks under ``output/`` are the source of truth; render
failures LOG loudly and re-raise but do NOT touch the symlinks
(that would violate the earlier atomicity contract from
``a92cd8e``).
"""
from __future__ import annotations

import html
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


def _render_thumbnails(*, sample_id: str, stages: list[_StageInfo]) -> str:
    """Render the viz overlay as the single top-level thumbnail.

    Register / warp / celltyped don't currently produce a canonical
    single image, so we only thumbnail viz. If viz isn't promoted
    (or its overlay is missing), report so; the ``missing-file``
    flag on the viz row will name it.
    """
    stage_by_name = {s.stage: s for s in stages}
    v = stage_by_name.get("viz")
    if v is None or v.he_job_id is None:
        return "<p class='none'>viz not promoted &mdash; no overlay to show.</p>"
    overlay_rel = f"../viz/{v.he_job_id}/{sample_id}_overlay.png"
    return (
        f"<a href='{html.escape(overlay_rel)}'>"
        f"<img class='thumb' src='{html.escape(overlay_rel)}' alt='overlay'>"
        f"</a>"
    )


# ---------------------------------------------------------------------
# Atomic write
# ---------------------------------------------------------------------
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
