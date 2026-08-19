"""Stage-orchestration loop for hexenium.

Runs the requested subset of ``he_preprocess → register → warp →
celltype → viz``, each guarded by its own sentinel-existence resume check
(nuke with ``force_rerun``).

Three invocation modes (resolved via ``hexenium.layout.resolve_layout``):

* **standalone** — ``sample_id`` + ``output_root``.
* **integrated-by-run-id** — add ``run_id``; outputs colocate under
  ``<output_root>/<sample>/<sample>_<run_id>/he_registration/``.
* **integrated-by-h5ad** — pass ``xenium_h5ad``; identity read from
  ``.uns``; outputs colocate under ``<xenium_run_dir>/he_registration/``.

Path assembly is delegated to a :class:`RunLayout` — every stage
receives a pre-computed ``out_dir`` so the writer side never has to
know which mode produced it.
"""
from __future__ import annotations

import copy
import os
import socket
import sys
import time
from pathlib import Path

import yaml

from hexenium._internal.logging import log, banner
from hexenium.layout import (
    RunLayout,
    resolve_he_job_id,
    resolve_layout,
)
from hexenium.manifest import output_symlink_run_id


def _snapshotable(cfg: dict) -> dict:
    """Return a copy of ``cfg`` safe for YAML serialisation."""
    def _fix(v):
        if isinstance(v, Path):
            return str(v)
        if isinstance(v, dict):
            return {k: _fix(x) for k, x in v.items()}
        if isinstance(v, list):
            return [_fix(x) for x in v]
        return v
    return _fix(copy.deepcopy(cfg))


def _log_invocation_banner(argv: list[str], stages: list[str], cfg: dict,
                            layout: RunLayout) -> None:
    """Print everything anyone would want to grep out of the log."""
    banner("hexenium starting")
    log(f"command:      {' '.join(argv)}")
    log(f"host:         {socket.gethostname()}")
    log(f"cwd:          {os.getcwd()}")
    log(f"pid:          {os.getpid()}")
    log(f"python:       {sys.executable}")
    log(f"python ver:   {sys.version.splitlines()[0]}")
    log(f"sample_id:    {layout.sample_id}")
    log(f"he_job_id:    {layout.he_job_id}")
    log(f"mode:         {'INTEGRATED' if layout.integrated else 'STANDALONE'}")
    if layout.integrated:
        log(f"xenium_run_id: {layout.xenium_run_id}")
        log(f"xenium_h5ad:   {layout.xenium_h5ad}")
        log(f"xenium_run_dir: {layout.xenium_run_dir}")
    log(f"output_root_he: {layout.output_root_he}")
    log(f"he_path:      {cfg.get('he_path')}")
    log(f"xenium_bun:   {cfg.get('xenium_bundle')}")
    log(f"dapi_path:    {cfg.get('dapi_path') or '(derived from xenium_bundle)'}")
    log(f"stages:       {stages}")
    log(f"force_rerun:  {cfg.get('force_rerun', False)}")
    log(f"reg mode:     {cfg['registration']['mode']}")
    log(f"use_he_deconv:{cfg['registration']['use_he_deconvolution']}")
    log(f"warp.targets: {cfg['warp']['targets']}"
        f"{' +transcripts' if cfg['warp'].get('include_transcripts') else ''}")
    # Library versions — wrap each in its own try so one missing library
    # doesn't suppress the rest.
    for name in ("numpy", "yaml", "scipy", "scanpy", "anndata",
                 "geopandas", "shapely", "dask", "tifffile", "openslide",
                 "hest", "valis_hest"):
        try:
            mod = __import__(name)
            ver = getattr(mod, "__version__", "no __version__ attr")
            log(f"  {name:14s}{ver}")
        except Exception as e:
            log(f"  {name:14s}NOT IMPORTABLE ({type(e).__name__})")
    log("=" * 64)


def _merge_run_dir_config(layout: RunLayout, cfg: dict) -> None:
    """Write the ``he_registration:`` key into the xenium run-dir's
    ``resolved_config.yaml`` (integrated mode only).

    Preserves any pre-existing keys (``driver``, ``step1``, ``step4``, …)
    written by an upstream xenium-preprocess pipeline. Creates the file
    if the upstream side hasn't written it yet.
    """
    if not layout.integrated or layout.xenium_run_dir is None:
        return
    path = layout.xenium_run_dir / "resolved_config.yaml"
    existing = {}
    if path.exists():
        with open(path) as f:
            existing = yaml.safe_load(f) or {}
    existing["he_registration"] = {
        **layout.as_dict(),
        # A trimmed copy of the H&E-side config. Full snapshot still
        # lands in <layout.logs_dir>/resolved_config.yaml.
        "cfg": _snapshotable(cfg),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        yaml.safe_dump(existing, f, sort_keys=False)
    log(f"[pipeline] merged he_registration: key into {path}")


#: Stage name (from `--stages`) → dirname under `he_registration/`.
#: Only the four sharded stages participate in promotion; ``he_preprocess``
#: writes to the un-sharded ``converted/`` dir and has no ``output/`` link.
_PROMOTE_STAGE_DIR = {
    "register": "register",
    "warp": "warp",
    "celltype": "celltyped",
    "viz": "viz",
}


def _promote_completed_run(
    *,
    layout: RunLayout,
    cfg: dict,
    stages_ran: list[str],
    source_register_run_id: str | None,
) -> None:
    """Post-success promotion hook for ``--set-default-on-success``.

    Fires ONLY after all requested stages complete successfully. Called
    from the tail of :func:`run` — a stage failure raises before this
    runs, so any exception in the pipeline body leaves ``output/`` un-
    touched (the flag is opt-in and the promote is the LAST thing).

    Lineage resolution (per Tracy's ask in
    ``settylab/TracyY123-nexus#15`` comment ``5336523381``):

    * For each stage that ran this invocation → the fresh
      ``layout.he_job_id`` (the ``<he_job_id>`` this run's outputs
      landed under).
    * Un-run stages are NOT promoted — their ``output/<stage>``
      symlink keeps whatever it pointed at before, if anything.
      "Promote the completed run" means what actually completed.

    Lineage validation extension (see
    :func:`_validate_promotion_lineage`):

    * ``set_default_run`` already enforces the register↔warp check
      (celltyped/viz not in scope there because celltype/viz stages
      write no per-run manifest and thus record no ``source_warp_run_id``).
    * The extension here catches the "``--stages celltype viz`` with
      ``--warp-run-id X`` while ``output/warp`` points at ``Y != X``"
      case that the register↔warp check misses — a downstream-only
      re-promote that would leave ``output/`` inconsistent.

    The whole helper is atomic-by-refusal: it validates everything
    upfront and calls :func:`set_default_run` (which itself validates
    ALL targets exist before touching any symlink). If any check fails
    it raises ``SystemExit``; ``output/`` is unchanged.
    """
    from hexenium.set_default_run import format_changes, set_default_run

    promote_ids = _resolve_promotion_ids(
        layout=layout, stages_ran=stages_ran,
    )
    if not promote_ids:
        log("[promote] --set-default-on-success: no promotable stage in "
            f"--stages={stages_ran!r} (only he_preprocess ran? skipping).")
        return

    _validate_promotion_lineage(
        layout=layout,
        cfg=cfg,
        stages_ran=stages_ran,
        promote_ids=promote_ids,
        source_register_run_id=source_register_run_id,
    )

    kwargs = {
        "register_run_id":  promote_ids.get("register"),
        "warp_run_id":      promote_ids.get("warp"),
        "celltyped_run_id": promote_ids.get("celltype"),
        "viz_run_id":       promote_ids.get("viz"),
    }
    log(f"[promote] --set-default-on-success: promoting "
        f"{ {k: v for k, v in kwargs.items() if v} }")
    result = set_default_run(output_root_he=layout.output_root_he, **kwargs)
    log("[promote] output/ updated:\n" + format_changes(result))

    # Regenerate summary.html — the promotion above is atomic and
    # already committed; if the renderer fails, log LOUDLY but do NOT
    # roll back the symlinks (they're the source of truth). Tracy's
    # green-light: settylab/TracyY123-nexus#15 comment 5337905668.
    try:
        from hexenium.summary_html import render_summary_html
        render_summary_html(
            output_root_he=layout.output_root_he,
            sample_id=layout.sample_id,
            run_id=layout.xenium_run_id,
        )
    except Exception as exc:  # noqa: BLE001 — post-commit, log-and-continue
        log(f"[promote] ERROR: summary.html render failed after successful "
            f"symlink promotion: {exc!r}. The output/ symlinks are the "
            f"source of truth and remain unchanged; re-run "
            f"`hexenium set-default-run` to regenerate summary.html.")


def _resolve_promotion_ids(
    *, layout: RunLayout, stages_ran: list[str],
) -> dict[str, str]:
    """Which promotable stages ran and under which ``<he_job_id>``.

    Every stage in :data:`_PROMOTE_STAGE_DIR` that appears in
    ``stages_ran`` maps to the current invocation's fresh
    ``layout.he_job_id`` — that's the ``<he_job_id>`` its outputs just
    landed under. Un-run stages are omitted from the dict entirely,
    which is how ``set_default_run`` learns to leave their
    ``output/<stage>`` symlink alone.
    """
    return {
        stage: layout.he_job_id
        for stage in stages_ran
        if stage in _PROMOTE_STAGE_DIR
    }


def _validate_promotion_lineage(
    *,
    layout: RunLayout,
    cfg: dict,
    stages_ran: list[str],
    promote_ids: dict[str, str],
    source_register_run_id: str | None,
) -> None:
    """Cross-stage lineage checks that ``set_default_run`` alone can't do.

    ``set_default_run`` reads warp's per-run ``source_register_run_id``
    to catch register↔warp mismatches, but celltype/viz stages don't
    write per-run manifests, so a downstream-only promote where the
    warp INPUT (the stage that celltype/viz consumed) disagrees with
    ``output/warp`` would slip through. This function reads the
    pipeline's OWN resolved view (``cfg['warp_run_id']`` /
    ``layout.he_job_id``) to close that gap.

    Rules enforced:

    1. If ``celltype`` or ``viz`` is being promoted AND ``warp`` is
       NOT, ``output/warp`` MUST already exist and point at whichever
       warp celltype/viz just consumed. Otherwise the promoted
       downstream would be lineage-orphaned.
    2. If ``register`` is NOT being promoted but ``warp`` IS
       (rare: ``--stages warp`` alone), ``output/register`` MUST
       already exist and match the warp's
       ``source_register_run_id``. ``set_default_run``'s existing
       register↔warp check catches this too, but we duplicate the
       preflight here so we fail BEFORE handing off to set-default-run.
    """
    output_dir = layout.output_root_he / "output"

    downstream_ran = {"celltype", "viz"} & set(stages_ran)
    if downstream_ran and "warp" not in promote_ids:
        warp_consumed = cfg.get("warp_run_id") or layout.he_job_id
        current_output_warp = output_symlink_run_id(output_dir / "warp")
        if current_output_warp is None:
            raise SystemExit(
                "--set-default-on-success: cannot promote "
                f"{sorted(downstream_ran)!r} without also promoting warp — "
                f"output/warp does not exist yet.\n"
                f"Either include `warp` in --stages, or run "
                f"`hexenium set-default-run --warp-run-id {warp_consumed}` "
                f"first (which will also seed output/register)."
            )
        if current_output_warp != warp_consumed:
            raise SystemExit(
                "--set-default-on-success: lineage mismatch.\n"
                f"  This run's {sorted(downstream_ran)!r} consumed "
                f"warp/{warp_consumed}/\n"
                f"  but output/warp currently points at "
                f"warp/{current_output_warp}/.\n"
                f"Promoting would leave output/ inconsistent. Fix:\n"
                f"  * re-run with `--warp-run-id {current_output_warp}` "
                f"(match the current default warp), OR\n"
                f"  * `hexenium set-default-run --warp-run-id {warp_consumed}` "
                f"before this run (make warp={warp_consumed!r} the default), OR\n"
                f"  * include `warp` in --stages so a fresh warp is promoted "
                f"alongside downstream."
            )

    if "warp" in promote_ids and "register" not in promote_ids:
        current_output_register = output_symlink_run_id(output_dir / "register")
        if source_register_run_id is None:
            # Should be unreachable — warp always resolves a source. Guard
            # so a future refactor can't silently drop the field.
            raise SystemExit(
                "--set-default-on-success: internal error — warp promoted "
                "but source_register_run_id is unset. Refusing to promote."
            )
        if current_output_register is None:
            raise SystemExit(
                "--set-default-on-success: cannot promote warp without "
                "also promoting register — output/register does not exist "
                "yet.\n"
                f"Either include `register` in --stages, or run "
                f"`hexenium set-default-run --register-run-id "
                f"{source_register_run_id}` first."
            )
        if current_output_register != source_register_run_id:
            raise SystemExit(
                "--set-default-on-success: lineage mismatch.\n"
                f"  This run's warp consumed register/{source_register_run_id}/\n"
                f"  but output/register currently points at "
                f"register/{current_output_register}/.\n"
                f"Include `register` in --stages, or set-default-run the "
                f"correct register first."
            )


def _resolve_existing_registrar(
    *,
    register_root: Path,
    explicit_register_run_id: str | None,
) -> tuple[Path | None, str | None]:
    """Locate the registrar pickle to feed into warp when register is
    NOT in this invocation's ``--stages``.

    Precedence:

    1. ``explicit_register_run_id`` (from ``--register-run-id``) —
       read ``<register_root>/<X>/manifest.yaml`` directly. Fail LOUD
       if that manifest is missing rather than silently falling back to
       the default symlink; the operator asked for a specific pick and
       hiding a typo behind the fallback would poison downstream
       lineage checks.
    2. ``<register_root>/manifest.yaml`` — the default-symlink pointer
       the register stage maintains. Reads through the symlink to the
       per-run manifest. Returned ``source_register_run_id`` is
       resolved from the manifest's ``he_job_id`` field (falling back
       to the symlink target's leading path component for older
       manifests written before this field existed).

    Returns ``(registrar_pickle, source_register_run_id)``. Both are
    ``None`` when no register output exists yet — the warp stage's
    subsequent ``registrar_pickle is None`` guard raises the actionable
    error message.
    """
    from hexenium.manifest import read_manifest, symlink_target_run_id

    if explicit_register_run_id:
        per_run = register_root / explicit_register_run_id / "manifest.yaml"
        if not per_run.exists():
            raise SystemExit(
                f"--register-run-id {explicit_register_run_id!r}: "
                f"manifest not found at {per_run}. Either fix the run id "
                f"or drop the flag to fall back to the default (the last "
                f"successful register writes {register_root}/manifest.yaml)."
            )
        m = read_manifest(per_run)
        log(f"[pipeline] --register-run-id {explicit_register_run_id!r} -> {per_run}")
        return Path(m["registrar_pickle"]), explicit_register_run_id

    default_link = register_root / "manifest.yaml"
    if not default_link.exists():
        return None, None
    m = read_manifest(default_link)
    source_id = m.get("he_job_id") or symlink_target_run_id(default_link)
    log(f"[pipeline] using default registrar (he_job_id={source_id!r}) "
        f"from {default_link}")
    return Path(m["registrar_pickle"]), source_id


def _maybe_derive_proseg_purified(stages, cfg: dict, layout: RunLayout) -> None:
    """Auto-derive ``proseg_purified_h5ad`` from the xenium run dir when
    the user didn't set it explicitly.

    Only fires when celltype is in ``stages`` AND we're in an integrated
    mode (either --xenium-h5ad or --run-id gives us a xenium_run_dir).
    Standalone mode has no run_dir; the user must pass the path there.
    Explicit ``--proseg-purified-h5ad`` always wins (short-circuited before
    this runs).

    Mutates ``cfg`` in place. Raises ``SystemExit`` with an actionable
    message if celltype requires proseg but the derived path is missing.
    """
    if "celltype" not in stages:
        return
    if cfg.get("proseg_purified_h5ad"):
        return
    if layout.xenium_run_dir is None:
        return
    derived = (layout.xenium_run_dir / "spatial_adata"
               / f"{layout.sample_id}_proseg_purified.h5ad")
    if derived.exists():
        cfg["proseg_purified_h5ad"] = str(derived)
        log(f"[pipeline] auto-detected proseg_purified h5ad -> {derived}")
        return
    raise SystemExit(
        f"proseg_purified h5ad not found at expected path:\n"
        f"  {derived}\n"
        f"the celltype stage requires this file. Either run an upstream "
        f"proseg-purified export first, or pass --proseg-purified-h5ad "
        f"<path> to override.\n"
        f"Identity: sample_id={layout.sample_id!r} "
        f"xenium_run_dir={layout.xenium_run_dir}."
    )


def run(cfg: dict, stages: list[str], argv: list[str]) -> int:
    """Execute the pipeline. Assumes `cfg` has already been validated."""
    he_job_id = resolve_he_job_id(cfg.get("he_job_id"))
    layout = resolve_layout(
        sample_id=cfg.get("sample_id"),
        output_root=Path(cfg["output_root"]) if cfg.get("output_root") else None,
        xenium_h5ad=Path(cfg["xenium_h5ad"]) if cfg.get("xenium_h5ad") else None,
        he_job_id=he_job_id,
        run_id=cfg.get("run_id"),
        stages=stages,
    )
    # Back-populate cfg with the resolved sample_id (integrated-by-h5ad
    # mode reads it from .uns so the config may still have sample_id=null).
    cfg["sample_id"] = layout.sample_id
    cfg["he_job_id"] = layout.he_job_id

    _maybe_derive_proseg_purified(stages, cfg, layout)

    layout.ensure_dirs()

    he_path = Path(cfg["he_path"]).resolve()
    xenium_bundle = Path(cfg["xenium_bundle"]).resolve()
    dapi_path_cfg = cfg.get("dapi_path")
    dapi_path = Path(dapi_path_cfg).resolve() if dapi_path_cfg else None
    force_rerun = bool(cfg.get("force_rerun", False))

    _log_invocation_banner(argv, stages, cfg, layout)

    if not he_path.exists():
        raise SystemExit(f"H&E not found: {he_path}")
    if not xenium_bundle.exists():
        raise SystemExit(f"xenium bundle not found: {xenium_bundle}")

    # Snapshot resolved config to the per-<he_job_id> logs dir.
    snap = layout.logs_dir / "resolved_config.yaml"
    with open(snap, "w") as f:
        yaml.safe_dump(_snapshotable(cfg), f, sort_keys=False)
    log(f"[pipeline] resolved config snapshot -> {snap}")

    # In integrated mode, ALSO merge a he_registration: key into the
    # xenium run-dir's shared resolved_config.yaml.
    _merge_run_dir_config(layout, cfg)

    n_stages = len(stages)
    registrar_pickle = None
    # Tracks which register/<he_job_id>/ the current invocation is
    # consuming. Populated by same-invocation register OR by the
    # manifest-resolution branch below (--register-run-id explicit /
    # default-symlink implicit). Passed into warp's per-run manifest as
    # `source_register_run_id` so set-default-run can validate lineage.
    source_register_run_id: str | None = None
    # Default warp/celltype dirs for the current run — may be overridden
    # by --warp-run-id / --celltype-run-id below to point at prior runs.
    warp_dir = layout.warp_dir(cfg.get("warp_run_id"))
    celltyped_dir = layout.celltyped_dir(cfg.get("celltype_run_id"))
    pipeline_t0 = time.time()

    force_preprocess = bool(cfg.get("force_preprocess", False))

    # Stage 0: VSI → OME-TIFF (no-op if input isn't a .vsi file).
    if "he_preprocess" in stages:
        idx = stages.index("he_preprocess") + 1
        banner(f"stage {idx}/{n_stages}: he_preprocess — starting")
        t0 = time.time()
        from hexenium.stages.he_preprocess import run_he_preprocess
        hp_cfg = cfg.get("he_preprocess", {}) or {}
        new_he_path = run_he_preprocess(
            sample_id=layout.sample_id,
            he_path=he_path,
            out_dir=layout.converted_dir,
            tile_size=hp_cfg.get("tile_size", 1024),
            compression=hp_cfg.get("compression", "jpeg2000"),
            jpeg2000_level=hp_cfg.get("jpeg2000_level", 100),
            pyramid_levels=hp_cfg.get("pyramid_levels", 7),
            force_rerun=force_rerun or force_preprocess,
        )
        if new_he_path != he_path:
            log(f"[pipeline] register stage will use converted OME-TIFF: {new_he_path}")
            he_path = new_he_path
        banner(f"stage {idx}/{n_stages}: he_preprocess — complete in {time.time()-t0:.1f}s")
    elif he_path.suffix.lower() == ".vsi":
        # Auto-detect an existing per-sample OME-TIFF from a prior run so
        # `--stages register warp ...` (i.e. skipping he_preprocess) works
        # without threading the converted path back through the CLI.
        from hexenium.stages.he_preprocess import discover_existing_ometiff
        existing = discover_existing_ometiff(
            he_path, layout.converted_dir, layout.sample_id,
        )
        if existing is not None:
            log(f"[pipeline] --he-path is a .vsi and he_preprocess is not in "
                f"--stages, but a converted OME-TIFF exists at {existing} — "
                f"using it as the H&E input for downstream stages")
            he_path = existing
        else:
            log(f"[pipeline] WARN: --he-path is a .vsi file but stage he_preprocess "
                f"is not in --stages and no converted OME-TIFF was found next to "
                f"the source or under {layout.converted_dir}. Register will likely fail.")

    if "register" in stages:
        idx = stages.index("register") + 1
        banner(f"stage {idx}/{n_stages}: register — starting")
        t0 = time.time()
        from hexenium.stages.registration import run_registration
        reg = cfg["registration"]
        params = cfg["parameters"]
        he_cfg = cfg["he"]
        registrar_pickle = run_registration(
            sample_id=layout.sample_id,
            he_path=he_path,
            xenium_bundle=xenium_bundle,
            out_dir=layout.registration_dir,
            dapi_path=dapi_path,
            mode=reg["mode"],
            use_he_deconvolution=reg["use_he_deconvolution"],
            check_for_reflections=reg["check_for_reflections"],
            create_masks=reg["create_masks"],
            align_to_reference=reg["align_to_reference"],
            max_image_dim_px=params["max_image_dim_px"],
            max_processed_image_dim_px=params["max_processed_image_dim_px"],
            max_non_rigid_registration_dim_px=params["max_non_rigid_registration_dim_px"],
            micro_rigid_registrar_cls=params.get("micro_rigid_registrar_cls"),
            micro_rigid_registrar_params=params.get("micro_rigid_registrar_params") or {},
            name=params.get("name"),
            symlink_to_canonical_name=he_cfg["symlink_to_canonical_name"],
            force_rerun=force_rerun,
        )
        # Same-invocation register+warp: source register IS this run's
        # <he_job_id>. Recorded in the warp manifest for lineage.
        source_register_run_id = layout.he_job_id
        banner(f"stage {idx}/{n_stages}: register — complete in {time.time()-t0:.1f}s")
    else:
        # Warp-without-register: resolve which register/<X>/ to consume.
        registrar_pickle, source_register_run_id = _resolve_existing_registrar(
            register_root=layout.registration_dir.parent,
            explicit_register_run_id=cfg.get("register_run_id"),
        )

    if "warp" in stages:
        idx = stages.index("warp") + 1
        banner(f"stage {idx}/{n_stages}: warp — starting")
        t0 = time.time()
        from hexenium.stages.warp import run_warp
        if registrar_pickle is None or not registrar_pickle.exists():
            raise SystemExit(
                "warp stage requested but no registrar pickle is available. "
                "Run the register stage first or place a manifest at "
                f"{layout.registration_dir.parent / 'manifest.yaml'}."
            )
        warp_cfg = cfg["warp"]
        # Compose final targets list: start from configured list, append
        # 'transcripts' when include_transcripts is set, then dedupe while
        # preserving order.
        targets = list(warp_cfg["targets"])
        if warp_cfg.get("include_transcripts") and "transcripts" not in targets:
            targets.append("transcripts")
        seen = set()
        targets = [t for t in targets if not (t in seen or seen.add(t))]
        log(f"[warp] final targets: {targets}")
        warp_dir = run_warp(
            sample_id=layout.sample_id,
            registrar_pickle=registrar_pickle,
            xenium_bundle=xenium_bundle,
            out_dir=layout.warp_dir(),
            dapi_path=dapi_path,
            targets=targets,
            source_register_run_id=source_register_run_id,
            use_dask=warp_cfg["use_dask"],
            save_geojson=warp_cfg["save_geojson"],
            dask_n_workers=warp_cfg["dask"]["n_workers"],
            dask_threads_per_worker=warp_cfg["dask"]["threads_per_worker"],
            dask_memory_limit=warp_cfg["dask"]["memory_limit"],
            dask_jvm_mem_gb=warp_cfg["dask"]["jvm_mem_gb"],
            force_rerun=force_rerun,
        )
        banner(f"stage {idx}/{n_stages}: warp — complete in {time.time()-t0:.1f}s")

    if "celltype" in stages:
        idx = stages.index("celltype") + 1
        banner(f"stage {idx}/{n_stages}: celltype — starting")
        t0 = time.time()
        from hexenium.stages.celltyping import run_celltyping
        ct_cfg = cfg["celltype"]
        # Resolve h5ad paths: celltype.<key> > top-level <key> > None.
        xh5 = ct_cfg.get("xenium_h5ad") or cfg.get("xenium_h5ad")
        ph5 = ct_cfg.get("proseg_purified_h5ad") or cfg.get("proseg_purified_h5ad")
        ct_result = run_celltyping(
            sample_id=layout.sample_id,
            warp_dir=warp_dir,
            out_dir=layout.celltyped_dir(),
            xenium_h5ad=Path(xh5) if xh5 else None,
            proseg_purified_h5ad=Path(ph5) if ph5 else None,
            celltype_col=ct_cfg.get("celltype_col", "auto"),
            id_col=ct_cfg.get("id_col", "auto"),
            proseg_x_col=ct_cfg.get("proseg_x_col"),
            proseg_y_col=ct_cfg.get("proseg_y_col"),
            xenium_x_col=ct_cfg.get("xenium_x_col"),
            xenium_y_col=ct_cfg.get("xenium_y_col"),
            nn_k=ct_cfg.get("nn_k", 1),
            nn_algorithm=ct_cfg.get("nn_algorithm", "auto"),
            nn_metric=ct_cfg.get("nn_metric", "euclidean"),
            nuclei_inherit_classification=ct_cfg["nuclei_inherit_classification"],
            area_threshold_px=ct_cfg["area_threshold_px"],
            nucleus_round_ndigits=ct_cfg["nucleus_round_ndigits"],
            cell_round_ndigits=ct_cfg["cell_round_ndigits"],
            qupath_cell_round_ndigits=ct_cfg.get("qupath_cell_round_ndigits", 2),
            qupath_nucleus_round_ndigits=ct_cfg.get("qupath_nucleus_round_ndigits", 2),
            force_rerun=force_rerun,
        )
        celltyped_dir = ct_result["out_dir"]
        banner(f"stage {idx}/{n_stages}: celltype — complete in {time.time()-t0:.1f}s")

    if "viz" in stages and cfg["viz"]["enabled"]:
        idx = stages.index("viz") + 1
        banner(f"stage {idx}/{n_stages}: viz — starting")
        t0 = time.time()
        from hexenium.stages.viz import run_viz
        v = cfg["viz"]
        run_viz(
            sample_id=layout.sample_id,
            he_path=he_path,
            warp_dir=warp_dir,
            celltyped_dir=celltyped_dir,
            out_dir=layout.viz_dir,
            thumbnail_max_dim=v["thumbnail_max_dim"],
            dpi=v["dpi"],
            cell_alpha=v["cell_alpha"],
            nucleus_alpha=v["nucleus_alpha"],
            classification_palette=v.get("classification_palette"),
            palette_cmap=v.get("palette_cmap", "tab20"),
            render_boundaries=v.get("render_boundaries", "nucleus"),
            force_rerun=force_rerun,
        )
        banner(f"stage {idx}/{n_stages}: viz — complete in {time.time()-t0:.1f}s")

    # Post-success promotion (opt-in). Fires ONLY after all stages
    # above completed without raising — a stage-level SystemExit /
    # RuntimeError propagates out of this function before this point,
    # so the promote never sees a partial-success run and `output/`
    # stays put on failure. Any SystemExit inside the promote itself
    # (lineage mismatch, missing target) also propagates out unchanged
    # — set_default_run validates ALL targets up-front so we're
    # atomic-by-refusal.
    if cfg.get("set_default_on_success"):
        _promote_completed_run(
            layout=layout,
            cfg=cfg,
            stages_ran=list(stages),
            source_register_run_id=source_register_run_id,
        )

    banner(f"hexenium done in {time.time()-pipeline_t0:.1f}s "
           f"-> {layout.output_root_he}/")
    return 0
