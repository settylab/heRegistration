"""Per-run manifests and lineage helpers for hexenium stages.

Each stage writes an immutable per-run manifest to
``<stage>/<he_job_id>/manifest.yaml``. A top-level
``<stage>/manifest.yaml`` is maintained as a **relative symlink** to
the current-default run's manifest. This gives us three things:

* **Immutability** — a per-run manifest is never overwritten after its
  stage completes, so downstream reproducibility (params_hash lookup,
  lineage checks) can always resolve what a specific ``<he_job_id>``
  ran with, even after the "default" pointer moves.
* **A single default pointer per stage** — resolving ``manifest.yaml``
  at the stage root goes through the symlink to whatever the current
  default is. Warp-only invocations without ``--register-run-id``
  follow this pointer.
* **Cross-stage lineage** — the warp manifest records
  ``source_register_run_id`` (which registration it consumed), so a
  later ``set-default-run`` can validate consistency between the
  register and warp choices before re-pointing user-facing symlinks.

Manifest fields (both stages share ``sample_id``, ``he_job_id``,
``params_hash``, ``git_sha``, ``timestamp_utc``; each adds stage-
specific keys). ``params_hash`` is a stable digest over the resolved
stage-parameter dict — key order and value types are canonicalised
before hashing so semantically-equivalent invocations agree.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import time
from pathlib import Path
from typing import Any

import yaml


PARAMS_HASH_LEN = 16


def compute_params_hash(params: dict[str, Any]) -> str:
    """Deterministic hex digest over ``params``.

    Canonicalises via ``json.dumps(..., sort_keys=True, default=str)``
    so semantically-equivalent invocations produce equal hashes:
    ``Path`` values, ints, floats, and nested dicts all round-trip
    through the ``default=str`` conversion. Returns the first
    :data:`PARAMS_HASH_LEN` hex chars of SHA-256 — full 256 bits is
    unnecessary for a provenance tag, and 16 chars greps cleanly.
    """
    canonical = json.dumps(params, sort_keys=True, default=str, ensure_ascii=False)
    digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    return digest[:PARAMS_HASH_LEN]


def git_sha(source_path: Path | None = None) -> str | None:
    """Return HEAD SHA of the git repo containing ``source_path``.

    Defaults to the repo containing this module — so a checkout of
    ``xenium-he-registration`` under the ``hexenium.` `` package resolves
    correctly regardless of the caller's cwd. Returns ``None`` if not
    inside a git repo or if ``git`` isn't on PATH — provenance is
    best-effort, never a hard failure.
    """
    if source_path is None:
        source_path = Path(__file__).resolve().parent
    try:
        result = subprocess.run(
            ["git", "-C", str(source_path), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None
    if result.returncode != 0:
        return None
    sha = result.stdout.strip()
    return sha or None


def utc_timestamp() -> str:
    """ISO-8601 UTC timestamp — stable, sortable, greppable."""
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def write_manifest(path: Path, payload: dict[str, Any]) -> None:
    """Atomically write ``payload`` as YAML to ``path``.

    Uses a tmp-file + ``os.replace`` so a crash mid-write never leaves
    the manifest half-serialised (which would poison downstream
    ``read_manifest`` calls with a YAML parse error).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    with open(tmp, "w") as f:
        yaml.safe_dump(payload, f, sort_keys=False)
    os.replace(tmp, path)


def read_manifest(path: Path) -> dict[str, Any]:
    """Read a YAML manifest, following symlinks transparently.

    Raises ``FileNotFoundError`` with a message that names the resolved
    target — critical when the caller passed the stage-root symlink and
    the target's per-run dir was deleted.
    """
    p = Path(path)
    if not p.exists():
        resolved = p.resolve() if p.is_symlink() else p
        raise FileNotFoundError(
            f"manifest not found: {p} (resolves to {resolved})"
        )
    with open(p) as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"manifest at {p} did not parse as a YAML mapping")
    return data


def set_default_symlink(link_path: Path, target_relative: str) -> None:
    """Atomically point ``link_path`` at ``target_relative`` (relative).

    Uses the tmp-symlink + ``os.replace`` idiom so an interrupted
    ``set-default-run`` never leaves ``link_path`` missing.
    ``target_relative`` is stored verbatim so pointers stay portable
    when the stage tree is moved or symlinked into place.

    Refuses to clobber a non-symlink at ``link_path`` — that would be a
    real file or directory the operator manages themselves; silently
    replacing it with a symlink would destroy content.
    """
    link_path.parent.mkdir(parents=True, exist_ok=True)
    if link_path.exists() and not link_path.is_symlink():
        raise FileExistsError(
            f"cannot re-point {link_path}: exists and is not a symlink. "
            f"Investigate and remove it manually before re-running."
        )
    tmp = link_path.with_name(link_path.name + f".tmp.{os.getpid()}")
    if tmp.exists() or tmp.is_symlink():
        tmp.unlink()
    os.symlink(target_relative, tmp)
    os.replace(tmp, link_path)


def read_symlink_target(link_path: Path) -> str | None:
    """Return the raw target of a symlink, or ``None`` if not a symlink."""
    if not link_path.is_symlink():
        return None
    return os.readlink(link_path)


def symlink_target_run_id(link_path: Path) -> str | None:
    """Extract the ``<he_job_id>`` component from a **stage-root**
    manifest symlink, i.e. one whose target is ``<he_job_id>/manifest.yaml``.

    Returns the leading path component (after stripping any leading
    ``..`` segments) — for that scheme, the leading component IS the
    run id. Returns ``None`` when the link is missing or malformed so
    callers can treat "no default set" as a first-class state.

    Note: does NOT work on ``output/<stage>`` symlinks (whose target is
    ``../<stage>/<he_job_id>``); use :func:`output_symlink_run_id` for
    those.
    """
    raw = read_symlink_target(link_path)
    if raw is None:
        return None
    parts = [p for p in Path(raw).parts if p not in ("..", ".", "/")]
    if not parts:
        return None
    return parts[0]


def output_symlink_run_id(link_path: Path) -> str | None:
    """Extract the ``<he_job_id>`` from an ``output/<stage>`` symlink.

    Scheme: ``output/<stage> -> ../<stage>/<he_job_id>`` — so the run
    id is the LAST path component of the target. Returns ``None`` when
    the link doesn't exist or isn't a symlink. Robust to dangling links
    (does not resolve() through the filesystem — pure textual parse).
    """
    raw = read_symlink_target(link_path)
    if raw is None:
        return None
    return Path(raw).name or None
