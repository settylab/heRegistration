"""``hexenium set-default-run`` — retarget the browsable ``output/``
symlinks after Tracy visually verifies a run.

Semantics (green-lit in ``settylab/TracyY123-nexus#15``):

* **Partial updates** — every ``--*-run-id`` flag is independent; only
  the stages the operator names are re-pointed. Others keep pointing
  wherever they were.
* **Lineage validation (register ↔ warp)** — before re-pointing
  ``output/register`` or ``output/warp``, the effective (post-update)
  register-id and warp-id are cross-checked against the warp's
  recorded ``source_register_run_id``. Mismatches refuse with an
  actionable error unless ``--force-lineage`` is passed.
* **Non-destructive** — only symlinks under ``output/`` are touched.
  Nothing under ``register/``, ``warp/``, ``celltyped/``, ``viz/`` is
  moved or deleted; every historical run stays on disk.
* **Idempotent** — re-running with the same target is a no-op with a
  clear "unchanged" message.

Symlink scheme: ``output/<stage> -> ../<stage>/<he_job_id>`` (relative,
so the ``he_registration/`` tree stays movable).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from hexenium._internal.logging import log
from hexenium.manifest import (
    output_symlink_run_id,
    read_manifest,
    set_default_symlink,
)


STAGES: tuple[str, ...] = ("register", "warp", "celltyped", "viz")


@dataclass(frozen=True)
class SetDefaultResult:
    """Human-readable summary of what set-default-run did.

    ``changes`` is one line per stage the operator asked to update
    (including no-op unchanged lines, so an idempotent second run
    still logs something). ``output_dir`` is the ``output/`` folder
    that was re-populated.
    """
    output_dir: Path
    changes: list[str]


def set_default_run(
    *,
    output_root_he: Path,
    register_run_id: str | None = None,
    warp_run_id: str | None = None,
    celltyped_run_id: str | None = None,
    viz_run_id: str | None = None,
    force_lineage: bool = False,
) -> SetDefaultResult:
    """Re-point ``output/<stage>`` symlinks and return a change summary.

    ``output_root_he`` is the ``he_registration/`` directory (either
    ``<output_root>/<sample>/`` in standalone mode, or
    ``<xenium_run_dir>/he_registration/`` in integrated mode).
    """
    requested = {
        "register": register_run_id,
        "warp": warp_run_id,
        "celltyped": celltyped_run_id,
        "viz": viz_run_id,
    }
    provided = {s: rid for s, rid in requested.items() if rid}
    if not provided:
        raise SystemExit(
            "set-default-run: at least one of --register-run-id, "
            "--warp-run-id, --celltyped-run-id, --viz-run-id is required."
        )

    output_dir = output_root_he / "output"

    _validate_targets_exist(output_root_he, provided)

    _validate_register_warp_lineage(
        output_root_he=output_root_he,
        output_dir=output_dir,
        provided=provided,
        force_lineage=force_lineage,
    )

    changes = _apply_updates(output_dir=output_dir, provided=provided)
    return SetDefaultResult(output_dir=output_dir, changes=changes)


# ---------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------
def _validate_targets_exist(
    output_root_he: Path, provided: dict[str, str]
) -> None:
    """Refuse the whole invocation if any named ``<he_job_id>`` is missing.

    Half-relinking then failing on the last stage would leave the
    ``output/`` view inconsistent — better to refuse upfront and
    require the operator to fix all missing dirs before we touch
    anything.
    """
    missing = []
    for stage, run_id in provided.items():
        target = output_root_he / stage / run_id
        if not target.is_dir():
            missing.append(f"{stage}/{run_id} (looked for {target})")
    if missing:
        raise SystemExit(
            "set-default-run: target directories not found:\n  "
            + "\n  ".join(missing)
            + "\nNo symlinks were touched. Fix the run ids (or run the "
            "missing stages first) and re-run."
        )


def _effective_ids(
    *, output_dir: Path, provided: dict[str, str]
) -> dict[str, str | None]:
    """The ``<he_job_id>`` each output symlink WILL point at after this run.

    Uses the requested value where given; otherwise the CURRENT symlink
    target's leading ``<he_job_id>`` (``None`` if there is no current
    symlink — first-time invocation).
    """
    current = {
        stage: output_symlink_run_id(output_dir / stage) for stage in STAGES
    }
    return {stage: provided.get(stage, current[stage]) for stage in STAGES}


def _validate_register_warp_lineage(
    *,
    output_root_he: Path,
    output_dir: Path,
    provided: dict[str, str],
    force_lineage: bool,
) -> None:
    """Refuse register↔warp mismatches (unless ``--force-lineage``).

    The check reads ``warp/<W>/manifest.yaml``'s recorded
    ``source_register_run_id`` and compares it to the effective
    register id. Runs when EITHER register or warp is being changed
    (an unchanged output still has to stay consistent with the updated
    counterpart).

    Skips silently when the warp manifest doesn't have
    ``source_register_run_id`` — older warp outputs written before
    this iteration have no lineage to check.
    """
    if "register" not in provided and "warp" not in provided:
        return
    effective = _effective_ids(output_dir=output_dir, provided=provided)
    r_id = effective["register"]
    w_id = effective["warp"]
    if r_id is None or w_id is None:
        return  # nothing to compare against yet.

    warp_manifest_path = output_root_he / "warp" / w_id / "manifest.yaml"
    if not warp_manifest_path.exists():
        return
    m = read_manifest(warp_manifest_path)
    recorded = m.get("source_register_run_id")
    if recorded is None:
        return
    if recorded == r_id:
        return

    msg = (
        f"lineage mismatch: warp/{w_id}/manifest.yaml records "
        f"source_register_run_id={recorded!r}, but the effective "
        f"output/register would point at {r_id!r}.\n"
        f"This means output/warp={w_id!r} was NOT generated against "
        f"register={r_id!r}. Either:\n"
        f"  - pick a warp that consumed register {r_id!r}, or\n"
        f"  - pick a register that warp {w_id!r} actually consumed "
        f"(recorded: {recorded!r}), or\n"
        f"  - pass --force-lineage if you deliberately want inconsistent "
        f"pointers (e.g. cross-comparing overlays)."
    )
    if not force_lineage:
        raise SystemExit(msg)
    log("[set-default-run] WARN: --force-lineage bypass: " + msg)


def _apply_updates(*, output_dir: Path, provided: dict[str, str]) -> list[str]:
    """Re-point each provided stage's ``output/<stage>`` symlink.

    Iteration is stable across ``STAGES`` (not the ``provided`` dict)
    so the log order is predictable and diffable. Each update is a
    tmp-symlink + ``os.replace`` — atomic per-stage; a KeyboardInterrupt
    between stages leaves each finished stage consistent.
    """
    changes: list[str] = []
    for stage in STAGES:
        if stage not in provided:
            continue
        target_run_id = provided[stage]
        link = output_dir / stage
        current = output_symlink_run_id(link)
        target_relative = f"../{stage}/{target_run_id}"
        if current == target_run_id:
            msg = f"{stage}: unchanged ({target_run_id})"
        else:
            set_default_symlink(link, target_relative)
            msg = f"{stage}: {current!r} -> {target_run_id!r}"
        changes.append(msg)
        log(f"[set-default-run] {msg}")
    return changes


def format_changes(result: SetDefaultResult) -> str:
    """Human-readable one-block summary — used by the CLI epilogue."""
    lines = [f"output dir: {result.output_dir}"]
    lines.extend(f"  {c}" for c in result.changes)
    return "\n".join(lines)


def stage_flag_names() -> Iterable[tuple[str, str]]:
    """Yield ``(cli_flag, attr_name)`` pairs — one per stage.

    Kept as a helper so the CLI argparse wiring and this module's
    keyword arguments stay in sync. Attr names match the
    ``set_default_run`` kwargs so ``**vars(args)`` filtering works.
    """
    for stage in STAGES:
        yield f"--{stage}-run-id", f"{stage}_run_id"
