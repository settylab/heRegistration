"""Output layout for the hexenium pipeline.

Three invocation modes:

- **standalone** — driven by ``--sample-id`` + ``--output-root``. Outputs
  land at ``<output_root>/<sample_id>/{logs,converted,register,...}``.

- **integrated-by-h5ad** — driven by ``--xenium-h5ad`` pointing at the
  xenium h5ad written by an upstream xenium-preprocess pipeline. Sample
  identity comes from ``.uns['sample_id']`` + ``.uns['run_id']`` on that
  h5ad; the H&E-reg run dir is derived from the h5ad's on-disk location
  (``h5ad.parent.parent``) and outputs are colocated under
  ``<run_dir>/he_registration/{logs,converted,...}``.

- **integrated-by-run-id** — driven by ``--sample-id`` + ``--run-id`` +
  ``--output-root``. Same colocated layout as integrated-by-h5ad, but
  without reading the h5ad (identity comes from the CLI args). Lets
  stages that don't need the h5ad content (register, warp) run before
  the upstream step-1 has produced it.

In BOTH modes the internal layout under the top-level output dir is
identical:

    <output_root_he>/
      converted/                — VSI → OME-TIFF conversion (stage 0)
      register/<he_job_id>/     — VALIS registrar + workdir (stage 1)
      warp/<he_job_id>/         — warped cell/nucleus/transcript parquets (stage 2)
      celltyped/<he_job_id>/    — GeoJSONs + combined parquet (stage 3)
      viz/<he_job_id>/          — overlay PNGs (stage 4)
      test_samples/             — operator-managed manual samples

Per-run **logs** land at a different place per mode:

* **integrated** → ``<run_dir>/logs/logs_heRegistration/<sample_id>_<he_job_id>_<stages>/``
  (colocated with any workflow-driver step logs under ``<run_dir>/logs/``).
* **standalone** → ``<output_root_he>/logs/<sample_id>_<he_job_id>_<stages>/``
  (no shared ``<run_dir>`` in standalone mode).

``<stages>`` is the underscore-joined list of stages that ran (from
``--stages``), or ``all`` when every stage in :data:`STAGE_NAMES` runs.
The suffix is determined ONCE at submission time and reused by both the
sbatch launcher (``submit_he_registration.sh``) and the pipeline
(:func:`compute_stages_suffix`), so the folder is stable across the
whole run.

Each per-run logs dir holds the slurm ``.out``/``.err`` and the
``resolved_config.yaml`` snapshot for that ``<he_job_id>``.

``<he_job_id>`` precedence (three tiers): ``--he-job-id`` > ``$SLURM_JOB_ID``
> ``YYYYMMDDTHHMMSS`` timestamp. The interactive-runs fallback is
essential so re-runs from an operator shell (no Slurm) don't collide.
"""
from __future__ import annotations

import os
import re
import time
from dataclasses import dataclass
from pathlib import Path


_XENIUM_UUID_RE = re.compile(r"^[a-z]{8}-\d+$")

#: Canonical stage names, in execution order. Kept here (not in
#: ``pipeline.py``) so the sbatch launcher can compute the SAME
#: stages-suffix as the pipeline without importing pipeline code.
STAGE_NAMES: tuple[str, ...] = (
    "he_preprocess", "register", "warp", "celltype", "viz",
)


def compute_stages_suffix(stages) -> str:
    """Folder-name suffix from a list of stages that ran.

    ``all`` when every canonical stage is present (regardless of CLI
    order); otherwise underscore-joined in CLI order. Empty input →
    ``none`` (a defensive sentinel — the pipeline's argparse default
    ensures this branch is not hit in normal use).
    """
    s = list(stages or [])
    if not s:
        return "none"
    if set(s) == set(STAGE_NAMES) and len(s) == len(STAGE_NAMES):
        return "all"
    return "_".join(s)


def resolve_he_job_id(explicit: str | None = None) -> str:
    """Pick the H&E-reg run id (``<he_job_id>``) at invocation time.

    Precedence:
      1. ``explicit`` — the ``--he-job-id`` CLI flag.
      2. ``$SLURM_JOB_ID`` — set automatically by every ``sbatch`` /
         ``srun`` submission on a Slurm cluster.
      3. ``time.strftime("%Y%m%dT%H%M%S")`` — timestamp fallback for
         interactive invocations, so two back-to-back runs never share
         a folder.
    """
    if explicit:
        return str(explicit).strip()
    slurm = os.environ.get("SLURM_JOB_ID", "").strip()
    if slurm:
        return slurm
    return time.strftime("%Y%m%dT%H%M%S")


@dataclass(frozen=True)
class RunLayout:
    """All output paths for one hexenium pipeline invocation."""

    sample_id: str
    he_job_id: str
    output_root_he: Path
    integrated: bool = False
    xenium_run_id: str | None = None
    xenium_h5ad: Path | None = None
    xenium_run_dir: Path | None = None
    #: Stages-suffix baked into ``logs_dir``. Empty string preserves a
    #: shape-agnostic layout for any test/caller that doesn't yet supply
    #: stages.
    stages_suffix: str = ""

    # ---- top-level dirs ------------------------------------------------
    @property
    def _logs_leaf(self) -> str:
        base = f"{self.sample_id}_{self.he_job_id}"
        return f"{base}_{self.stages_suffix}" if self.stages_suffix else base

    @property
    def logs_dir(self) -> Path:
        # Integrated: colocate under <run_dir>/logs/, in a
        # logs_heRegistration/ subfolder that disambiguates from any
        # workflow-driver's slurm-*-step*.out files. Standalone has no
        # shared <run_dir>.
        if self.integrated and self.xenium_run_dir is not None:
            return (self.xenium_run_dir / "logs" / "logs_heRegistration"
                    / self._logs_leaf)
        return self.output_root_he / "logs" / self._logs_leaf

    @property
    def converted_dir(self) -> Path:
        return self.output_root_he / "converted"

    @property
    def test_samples_dir(self) -> Path:
        return self.output_root_he / "test_samples"

    # ---- per-stage dirs (default = this run's <he_job_id>) ------------
    @property
    def registration_dir(self) -> Path:
        return self.output_root_he / "register" / self.he_job_id

    def warp_dir(self, run_id: str | None = None) -> Path:
        """Path to the warp stage's output.

        Defaults to this invocation's ``<he_job_id>``. Pass ``run_id``
        (via ``--warp-run-id``) to point at a prior warp when running
        celltype/viz stages standalone against an earlier warp.
        """
        return self.output_root_he / "warp" / (run_id or self.he_job_id)

    def celltyped_dir(self, run_id: str | None = None) -> Path:
        """Path to the celltype stage's output.

        Defaults to this invocation's ``<he_job_id>``. Pass ``run_id``
        (via ``--celltype-run-id``) to point viz at a prior celltype run.
        """
        return self.output_root_he / "celltyped" / (run_id or self.he_job_id)

    @property
    def viz_dir(self) -> Path:
        return self.output_root_he / "viz" / self.he_job_id

    def ensure_dirs(self) -> None:
        """Create the top-level directories eagerly.

        Individual stages still create their own sub-run-id folders when
        they actually run, but pre-creating ``logs/<he_job_id>/`` is
        important so any stage can drop a resolved-config snapshot there
        before its first output arrives.
        """
        self.output_root_he.mkdir(parents=True, exist_ok=True)
        self.logs_dir.mkdir(parents=True, exist_ok=True)

    def as_dict(self) -> dict:
        """Serialise for the resolved-config snapshot's ``he_registration:`` key."""
        d = {
            "sample_id": self.sample_id,
            "he_job_id": self.he_job_id,
            "output_root_he": str(self.output_root_he),
            "integrated": self.integrated,
            "stages_suffix": self.stages_suffix,
        }
        if self.integrated:
            d["xenium_run_id"] = self.xenium_run_id
            d["xenium_h5ad"] = str(self.xenium_h5ad) if self.xenium_h5ad else None
            d["xenium_run_dir"] = str(self.xenium_run_dir) if self.xenium_run_dir else None
        return d


def resolve_layout(
    *,
    sample_id: str | None,
    output_root: Path | None,
    xenium_h5ad: Path | None,
    he_job_id: str,
    run_id: str | None = None,
    stages: list[str] | tuple[str, ...] | None = None,
) -> RunLayout:
    """Pick the right layout for the invocation.

    Three modes:

    * ``xenium_h5ad`` given → integrated mode: read
      ``.uns['sample_id']`` + ``.uns['run_id']`` and colocate outputs
      under ``<xenium_run_dir>/he_registration/``.
    * ``xenium_h5ad`` absent AND ``sample_id + output_root + run_id`` all
      given → integrated-by-run-id: construct the same colocation without
      touching the h5ad. Lets stages that don't need the h5ad content
      (register, warp) run before the upstream step-1 has produced it —
      identity comes from the CLI args, no ``.uns`` validation.
    * Otherwise → standalone mode driven by ``sample_id + output_root``,
      outputs at ``<output_root>/<sample_id>/``.

    ``stages``: when passed, becomes the logs-folder suffix (see
    :func:`compute_stages_suffix`). When absent, ``logs_dir`` falls back
    to the pre-suffix shape without any suffix — retained for backward
    compatibility with test callers.
    """
    stages_suffix = compute_stages_suffix(stages) if stages is not None else ""
    if xenium_h5ad is not None:
        return _resolve_integrated(
            sample_id=sample_id,
            xenium_h5ad=Path(xenium_h5ad),
            he_job_id=he_job_id,
            stages_suffix=stages_suffix,
        )
    if run_id and sample_id and output_root is not None:
        return _resolve_integrated_by_run_id(
            sample_id=sample_id,
            run_id=run_id,
            output_root=Path(output_root),
            he_job_id=he_job_id,
            stages_suffix=stages_suffix,
        )
    if not sample_id or output_root is None:
        raise SystemExit(
            "one of the two invocation modes must be complete: pass "
            "either --xenium-h5ad (integrated mode; sample identity is "
            "read from .uns metadata) OR --sample-id + --output-root "
            "(standalone mode)."
        )
    output_root_he = Path(output_root).resolve() / sample_id
    return RunLayout(
        sample_id=sample_id,
        he_job_id=he_job_id,
        output_root_he=output_root_he,
        stages_suffix=stages_suffix,
    )


def _resolve_integrated_by_run_id(
    *,
    sample_id: str,
    run_id: str,
    output_root: Path,
    he_job_id: str,
    stages_suffix: str = "",
) -> RunLayout:
    """Integrated mode without an h5ad: identity comes from CLI args.

    Layout matches ``_resolve_integrated`` — outputs land under
    ``<output_root>/<sample_id>/<sample_id>_<run_id>/he_registration/`` —
    but no ``.uns`` read, so this succeeds before the upstream step-1 has
    written the h5ad. The celltype stage will still need the h5ad at
    run-time; other stages don't.
    """
    run_dir = (Path(output_root).resolve() / sample_id
               / f"{sample_id}_{run_id}")
    return RunLayout(
        sample_id=sample_id,
        he_job_id=he_job_id,
        output_root_he=run_dir / "he_registration",
        integrated=True,
        xenium_run_id=run_id,
        xenium_h5ad=None,
        xenium_run_dir=run_dir,
        stages_suffix=stages_suffix,
    )


def _resolve_integrated(
    *,
    sample_id: str | None,
    xenium_h5ad: Path,
    he_job_id: str,
    stages_suffix: str = "",
) -> RunLayout:
    """Derive the integrated-mode layout from a xenium-preprocess h5ad.

    Reads ``adata.uns['sample_id']`` and ``adata.uns['run_id']`` (both
    are expected to be written by the upstream xenium-preprocess pipeline
    on its xenium_ranger_to_anndata stage). Fails LOUD if either is
    missing rather than falling back — silent identity drift between the
    two pipelines is exactly the class of bug this design exists to
    prevent.
    """
    h5ad = Path(xenium_h5ad).resolve()
    if not h5ad.exists():
        raise SystemExit(f"--xenium-h5ad path does not exist: {h5ad}")

    uns_sid, uns_rid = _read_uns_identity(h5ad)
    if not uns_sid or not uns_rid:
        raise SystemExit(
            f"--xenium-h5ad {h5ad} is missing .uns['sample_id']={uns_sid!r} "
            f"and/or .uns['run_id']={uns_rid!r}. The upstream "
            "xenium_ranger_to_anndata stage must write these before "
            "hexenium can integrate with it."
        )
    if sample_id is not None and sample_id != uns_sid:
        raise SystemExit(
            f"--sample-id {sample_id!r} disagrees with --xenium-h5ad's "
            f".uns['sample_id']={uns_sid!r}. In integrated mode, either "
            "omit --sample-id (auto-detect) or make them agree."
        )

    if h5ad.parent.name != "spatial_adata":
        raise SystemExit(
            f"--xenium-h5ad is expected to live at "
            f"<run_dir>/spatial_adata/<file>.h5ad; got parent="
            f"{h5ad.parent.name!r} for {h5ad}. Was the h5ad moved out "
            "of the upstream run dir?"
        )
    run_dir = h5ad.parent.parent
    expected_run_name = f"{uns_sid}_{uns_rid}"
    if run_dir.name != expected_run_name:
        raise SystemExit(
            f"xenium run dir name mismatch: expected {expected_run_name!r} "
            f"(from .uns[sample_id]={uns_sid!r} + .uns[run_id]={uns_rid!r}), "
            f"got {run_dir.name!r} at {run_dir}. Was the run dir manually "
            "renamed? Fix the folder name or rewrite the .uns fields."
        )

    output_root_he = run_dir / "he_registration"
    return RunLayout(
        sample_id=uns_sid,
        he_job_id=he_job_id,
        output_root_he=output_root_he,
        integrated=True,
        xenium_run_id=uns_rid,
        xenium_h5ad=h5ad,
        xenium_run_dir=run_dir,
        stages_suffix=stages_suffix,
    )


def _read_uns_identity(h5ad: Path) -> tuple[str | None, str | None]:
    """Return ``(sample_id, run_id)`` from ``.uns``.

    Uses backed='r' + immediate close to avoid loading X into memory
    — the h5ad can be tens of GB.
    """
    import anndata as ad
    adata = ad.read_h5ad(h5ad, backed="r")
    try:
        uns = adata.uns
        sid = uns.get("sample_id")
        rid = uns.get("run_id")
        # anndata sometimes wraps small scalars in numpy arrays.
        sid = str(sid) if sid is not None and str(sid).strip() else None
        rid = str(rid) if rid is not None and str(rid).strip() else None
        return sid, rid
    finally:
        try:
            adata.file.close()
        except Exception:
            pass


def is_xenium_uuid(value: object) -> bool:
    """Cheap shape-check: does ``value`` look like a xenium cell UUID?

    Xenium's per-cell UUIDs have the form ``aaaaafep-1`` — 8 lowercase
    letters, a dash, one digit. This is stable across every bundle
    since Xenium Ranger 1.x. Used by the celltype stage's id-column
    auto-detection to reject ambiguous alternatives (e.g. Proseg's
    ``cell_id`` column, which is an int64 index that would silently
    coerce to string ``"0"``, ``"1"``, … and cause every join row to
    miss).
    """
    if value is None:
        return False
    return bool(_XENIUM_UUID_RE.match(str(value)))
