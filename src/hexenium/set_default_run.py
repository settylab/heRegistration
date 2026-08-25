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

import os
from dataclasses import dataclass

from hexenium._internal.logging import log
from hexenium.manifest import (
    output_symlink_run_id,
    read_manifest,
    read_symlink_target,
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


@dataclass
class _PendingChange:
    """A single ``output/<stage>`` symlink about to be re-pointed.

    Captured up-front in Phase 0 so both the pre-write guard AND the
    Phase-2 rollback path have everything they need without another
    filesystem probe (which would introduce a TOCTOU window).
    """
    stage: str
    link_path: Path
    tmp_path: Path
    target_relative: str
    current_run_id: str | None      # ``output_symlink_run_id`` — for the log msg
    prior_symlink_target: str | None  # raw ``readlink`` — for rollback


def _apply_updates(*, output_dir: Path, provided: dict[str, str]) -> list[str]:
    """Re-point each provided stage's ``output/<stage>`` symlink.

    Batch-atomic in three phases (skeptic finding F1 on
    ``settylab/TracyY123-nexus#15``):

    * **Phase 0 — plan + guard.** Build the list of changes in
      ``STAGES`` order. For each stage that needs a new target, verify
      the destination is either missing or already a symlink (a
      real file/dir there would trigger the pre-existing
      ``set_default_symlink`` ``FileExistsError``, but only after
      earlier stages had already been re-pointed — so we surface
      the refusal BEFORE any write). Snapshot each stage's current
      ``readlink`` target so rollback can restore it byte-identical.
    * **Phase 1 — prep.** Create ALL tmp symlinks
      (``output/<stage>.tmp.<pid>``). No live ``output/<stage>`` is
      touched yet. If any tmp-create fails, unlink the ones we made
      and raise — ``output/`` unchanged.
    * **Phase 2 — commit.** ``os.replace`` each tmp → live path in
      order. On any exception mid-loop, walk back the completed
      renames and restore each to its prior target (or unlink if
      the stage had no prior symlink). Best-effort rollback: a
      rollback-of-rollback failure logs and continues so we don't
      lose the original exception context.

    ``KeyboardInterrupt`` between Phase-2 renames is caught by the
    same rollback branch; the operator sees ``output/`` untouched
    or in its pre-invocation state.
    """
    plan, unchanged = _build_plan(output_dir=output_dir, provided=provided)

    if not plan:
        changes = _format_change_log(plan=plan, unchanged=unchanged, provided=provided)
        for m in changes:
            log(f"[set-default-run] {m}")
        return changes

    output_dir.mkdir(parents=True, exist_ok=True)

    # Phase 1 — create ALL tmp symlinks; on any failure remove partial tmps.
    created_tmps: list[Path] = []
    try:
        for change in plan:
            _create_tmp_symlink(change.tmp_path, change.target_relative)
            created_tmps.append(change.tmp_path)
    except BaseException:
        for tmp in created_tmps:
            _best_effort_unlink(tmp)
        raise

    # Phase 2 — commit each tmp → live path; on any failure roll back.
    committed: list[_PendingChange] = []
    try:
        for change in plan:
            os.replace(change.tmp_path, change.link_path)
            committed.append(change)
    except BaseException:
        _rollback_committed(committed)
        # Clean up any tmps that never got renamed (e.g. we failed on
        # the middle one and later ones are still on disk).
        renamed_tmps = {c.tmp_path for c in committed}
        for tmp in created_tmps:
            if tmp not in renamed_tmps:
                _best_effort_unlink(tmp)
        raise

    changes = _format_change_log(plan=plan, unchanged=unchanged, provided=provided)
    for m in changes:
        log(f"[set-default-run] {m}")
    return changes


def _build_plan(
    *, output_dir: Path, provided: dict[str, str],
) -> tuple[list[_PendingChange], list[tuple[str, str]]]:
    """Turn ``provided`` into a list of ``_PendingChange`` (needs a write)
    plus a list of ``(stage, run_id)`` that are already at target.

    Also enforces the "non-symlink obstacle" guard: refuses UP-FRONT
    if any target ``output/<stage>`` exists but isn't a symlink. This
    is the F1 skeptic's empirically-most-likely trigger — an operator
    manually planted a real dir under ``output/`` to inspect a
    specific run. Refusing here (before any write) preserves batch
    atomicity that the old per-write ``FileExistsError`` in
    ``set_default_symlink`` couldn't.
    """
    to_change: list[_PendingChange] = []
    unchanged: list[tuple[str, str]] = []
    for stage in STAGES:
        if stage not in provided:
            continue
        target_run_id = provided[stage]
        link_path = output_dir / stage
        current = output_symlink_run_id(link_path)
        if current == target_run_id:
            unchanged.append((stage, target_run_id))
            continue
        if link_path.exists() and not link_path.is_symlink():
            raise FileExistsError(
                f"cannot re-point {link_path}: exists and is not a symlink. "
                f"Refusing to touch any output/ symlink for this batch — "
                f"investigate and remove {link_path} manually, then re-run."
            )
        prior = read_symlink_target(link_path) if link_path.is_symlink() else None
        to_change.append(_PendingChange(
            stage=stage,
            link_path=link_path,
            tmp_path=link_path.with_name(link_path.name + f".tmp.{os.getpid()}"),
            target_relative=f"../{stage}/{target_run_id}",
            current_run_id=current,
            prior_symlink_target=prior,
        ))
    return to_change, unchanged


def _create_tmp_symlink(tmp_path: Path, target_relative: str) -> None:
    """Create the Phase-1 tmp symlink (unlinking any leftover first)."""
    if tmp_path.exists() or tmp_path.is_symlink():
        tmp_path.unlink()
    os.symlink(target_relative, tmp_path)


def _rollback_committed(committed: list[_PendingChange]) -> None:
    """Restore each committed change to its pre-invocation state.

    Best-effort: if a single rollback fails (disk-full mid-rollback,
    permissions changed after we started), we log and continue so
    the original Phase-2 exception surfaces with as much of
    ``output/`` restored as possible.
    """
    for change in committed:
        try:
            if change.prior_symlink_target is None:
                change.link_path.unlink()
            else:
                _atomic_restore_symlink(
                    change.link_path, change.prior_symlink_target,
                )
        except Exception as exc:  # noqa: BLE001 — best-effort teardown
            log(f"[set-default-run] WARN: rollback of {change.stage!r} "
                f"({change.link_path}) failed: {exc!r}")


def _atomic_restore_symlink(link_path: Path, prior_target: str) -> None:
    """Point ``link_path`` back at ``prior_target`` via tmp + replace."""
    rb_tmp = link_path.with_name(link_path.name + f".rollback.{os.getpid()}")
    if rb_tmp.exists() or rb_tmp.is_symlink():
        rb_tmp.unlink()
    os.symlink(prior_target, rb_tmp)
    os.replace(rb_tmp, link_path)


def _best_effort_unlink(path: Path) -> None:
    """``unlink`` that swallows ``FileNotFoundError`` — for cleanup paths."""
    try:
        path.unlink()
    except FileNotFoundError:
        pass


def _format_change_log(
    *,
    plan: list[_PendingChange],
    unchanged: list[tuple[str, str]],
    provided: dict[str, str],
) -> list[str]:
    """Emit change-log lines in canonical ``STAGES`` order.

    Merges the "unchanged" and "changed" sets back into a single
    stable-ordered list so the log reads the same regardless of
    which stages were passed. Preserves the log surface the
    existing tests assert on (``'unchanged' in change``).
    """
    unchanged_map = dict(unchanged)
    plan_map = {c.stage: c for c in plan}
    lines: list[str] = []
    for stage in STAGES:
        if stage in unchanged_map:
            lines.append(f"{stage}: unchanged ({unchanged_map[stage]})")
        elif stage in plan_map:
            c = plan_map[stage]
            lines.append(f"{stage}: {c.current_run_id!r} -> {provided[stage]!r}")
    return lines


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
