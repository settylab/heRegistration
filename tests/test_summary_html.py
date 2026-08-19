"""Unit tests for ``hexenium.summary_html`` — the ``output/summary.html``
renderer wired into both promotion paths.

Coverage from Tracy's spec on ``settylab/TracyY123-nexus#15``
comment ``5337905668``:

* **Standard content** — per-stage table with he_job_id + params_hash
  + timestamp + git_sha + flags; lineage line for register↔warp.
* **Linked thumbnails** — the viz overlay is referenced via a
  relative path, not embedded.
* **Missing-file warning** — plant a symlink to a stage dir missing
  an expected artifact; assert the flag appears.
* **Param-drift warning** — mutate a manifest's ``params_hash`` so
  it disagrees with ``compute_params_hash(params)``; assert the flag
  appears.
* **No overlap-quality section** — the deferred column stays absent
  (regression guard).
* **Atomic write** — write twice, no stray ``.tmp.*`` residue.
* **Both promotion paths** — the two callers exercise the renderer.

Run:
    pytest tests/test_summary_html.py -v
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from hexenium.manifest import (
    compute_params_hash,
    set_default_symlink,
    write_manifest,
)
from hexenium.summary_html import (
    _EXPECTED_ARTIFACTS,
    render_summary_html,
)


def _seed_stage(root: Path, stage_dir_name: str, run_id: str, *,
                manifest_extra: dict | None = None,
                artifact_files: list[str] | None = None) -> Path:
    """Build a ``<root>/<stage>/<run_id>/`` folder with an optional
    manifest and optional artifact files.

    ``manifest_extra`` is merged into a baseline manifest with a
    self-consistent ``params_hash``; pass ``None`` to skip the
    manifest entirely (celltyped/viz stages don't write one today).
    """
    d = root / stage_dir_name / run_id
    d.mkdir(parents=True, exist_ok=True)
    if manifest_extra is not None:
        params = manifest_extra.pop("params", {"foo": "bar"})
        payload = {
            "sample_id": "S1",
            "he_job_id": run_id,
            "params": params,
            "params_hash": compute_params_hash(params),
            "git_sha": "0123456789abcdef",
            "timestamp_utc": "2026-08-18T18:00:00Z",
        }
        payload.update(manifest_extra)
        write_manifest(d / "manifest.yaml", payload)
    for rel in artifact_files or []:
        art = d / rel
        art.parent.mkdir(parents=True, exist_ok=True)
        art.touch()
    return d


def _seed_output_symlink(output_dir: Path, stage_dir_name: str,
                          run_id: str) -> None:
    """Point ``output/<stage>`` at ``../<stage>/<run_id>/``."""
    output_dir.mkdir(parents=True, exist_ok=True)
    link = output_dir / stage_dir_name
    set_default_symlink(link, f"../{stage_dir_name}/{run_id}")


@pytest.fixture
def he_root(tmp_path: Path) -> Path:
    """Full happy-path tree: all four stages seeded, all promoted, all
    artifacts present, manifests self-consistent."""
    root = tmp_path / "he_registration"
    output_dir = root / "output"
    output_dir.mkdir(parents=True)

    _seed_stage(
        root, "register", "reg_A",
        manifest_extra={
            "params": {"mode": "rigid_only", "max_image_dim_px": 1500},
        },
        artifact_files=["data/_registrar.pickle"],
    )
    _seed_stage(
        root, "warp", "warp_A",
        manifest_extra={
            "params": {"targets": ["cells", "nuclei"], "use_dask": True,
                       "save_geojson": True},
            "source_register_run_id": "reg_A",
        },
        artifact_files=["he_cell_seg.parquet", "he_nucleus_seg.parquet"],
    )
    _seed_stage(
        root, "celltyped", "ct_A",
        artifact_files=[
            "S1_cells_analysis.geojson",
            "S1_cells_qupath.geojson",
            "S1_nuclei_analysis.geojson",
            "S1_nuclei_qupath.geojson",
            "S1_celltyped_wholeslide.parquet",
        ],
    )
    _seed_stage(
        root, "viz", "viz_A",
        artifact_files=["S1_overlay.png"],
    )

    for stage_dir_name, rid in [
        ("register", "reg_A"), ("warp", "warp_A"),
        ("celltyped", "ct_A"), ("viz", "viz_A"),
    ]:
        _seed_output_symlink(output_dir, stage_dir_name, rid)

    return root


# ---------------------------------------------------------------------
# 1. Happy path — four stages, all sections present.
# ---------------------------------------------------------------------
class TestHappyPath:
    def test_renders_all_four_stages_with_manifests(self, he_root: Path):
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        assert dest == he_root / "output" / "summary.html"
        content = dest.read_text()

        # Title carries sample + run.
        assert "S1" in content
        assert "demo_v1" in content

        # Stable-order pointer summary line.
        assert "output/" in content
        for rid in ("reg_A", "warp_A", "ct_A", "viz_A"):
            assert rid in content

        # Per-stage table shows params_hash + timestamp.
        assert "params_hash" in content
        assert "2026-08-18T18:00:00Z" in content

        # Lineage line agrees.
        assert "lineage OK" in content
        assert "MISMATCH" not in content

    def test_linked_thumbnail_relative_path(self, he_root: Path):
        # Tracy chose LINKED thumbnails (not self-contained base64)
        # so the HTML must reference the overlay via a relative
        # ``../viz/<he_job_id>/<sample>_overlay.png`` path.
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        assert "src='../viz/viz_A/S1_overlay.png'" in content
        assert "href='../viz/viz_A/S1_overlay.png'" in content
        # No base64-embedded image.
        assert "base64," not in content


# ---------------------------------------------------------------------
# 2. Missing-file warning.
# ---------------------------------------------------------------------
class TestMissingFileWarning:
    def test_missing_registrar_pickle_flags_stage(self, he_root: Path):
        # Delete the expected artifact under the promoted register
        # dir; renderer must surface a missing-file flag on that row
        # without touching the other three stages.
        (he_root / "register" / "reg_A" / "data" / "_registrar.pickle").unlink()
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        assert "missing-file" in content
        assert "data/_registrar.pickle" in content

    def test_missing_flag_is_labeled_distinctly_from_drift(
        self, he_root: Path,
    ):
        # Tracy's spec (#15 comment 5334969034 answer 2): both warnings
        # kept with distinct labels — CSS classes differ so they can
        # be styled / grepped independently. Match against the SPAN
        # opening tag (not the class name alone) since the CSS style
        # block also mentions both class names.
        (he_root / "warp" / "warp_A" / "he_cell_seg.parquet").unlink()
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        assert "<span class='warn-missing'" in content
        assert "<span class='warn-drift'" not in content

    def test_no_warning_when_all_artifacts_present(self, he_root: Path):
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        assert "missing-file" not in content
        # register and warp both have manifests → "ok"; celltyped +
        # viz have no manifest → "no manifest".
        assert ">ok<" in content


# ---------------------------------------------------------------------
# 3. Param-drift warning — self-consistency check on stored hash.
# ---------------------------------------------------------------------
class TestParamDriftWarning:
    def test_mutated_hash_triggers_param_drift(self, he_root: Path):
        # Write a plausible-looking hash that doesn't match the params
        # — simulates a manually-edited manifest.
        import yaml as _yaml
        manifest_path = he_root / "register" / "reg_A" / "manifest.yaml"
        m = _yaml.safe_load(manifest_path.read_text())
        m["params_hash"] = "deadbeefdeadbeef"
        manifest_path.write_text(_yaml.safe_dump(m, sort_keys=False))
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        assert "<span class='warn-drift'" in content
        assert "param-drift</span>" in content

    def test_no_drift_when_hash_matches_params(self, he_root: Path):
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        # Match the SPAN opening tag (the class name itself lives in
        # the CSS block).
        assert "<span class='warn-drift'" not in content
        assert "param-drift</span>" not in content

    def test_missing_hash_field_does_not_flag_drift(self, he_root: Path):
        # If the manifest never had params_hash (older manifests),
        # skip the check silently rather than mis-flag.
        import yaml as _yaml
        manifest_path = he_root / "register" / "reg_A" / "manifest.yaml"
        m = _yaml.safe_load(manifest_path.read_text())
        m.pop("params_hash", None)
        manifest_path.write_text(_yaml.safe_dump(m, sort_keys=False))
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        assert "param-drift</span>" not in content


# ---------------------------------------------------------------------
# 7. Registration + warp params sections (Tracy's A3 + B2 on
# ``settylab/TracyY123-nexus#15`` comment ``5338129631``). Highlights
# line + collapsible full dump; no schema change.
# ---------------------------------------------------------------------
class TestParamsSections:
    def test_registration_params_section_renders_from_manifest(
        self, he_root: Path,
    ):
        # The fixture's register manifest has ``mode: rigid_only``.
        # Overwrite it with the load-bearing set so all A3 highlight
        # keys are exercised.
        import yaml as _yaml
        manifest_path = he_root / "register" / "reg_A" / "manifest.yaml"
        params = {
            "mode": "full_with_micro",
            "use_he_deconvolution": True,
            "check_for_reflections": True,
            "create_masks": False,
            "align_to_reference": True,
            "max_image_dim_px": 1500,
            "max_processed_image_dim_px": 1500,
            "max_non_rigid_registration_dim_px": 10000,
        }
        m = _yaml.safe_load(manifest_path.read_text())
        m["params"] = params
        m["params_hash"] = compute_params_hash(params)
        manifest_path.write_text(_yaml.safe_dump(m, sort_keys=False))
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()

        assert "<h2>Registration parameters</h2>" in content
        # A3 highlights line — the four load-bearing keys.
        for k in ("mode", "use_he_deconvolution",
                  "check_for_reflections", "align_to_reference"):
            assert k in content, f"expected {k} in the register highlights"
        # Values render — the specific ones on the fixture.
        assert "full_with_micro" in content
        # Resolution chip folds the three dim keys into one row.
        assert "1500/1500/10000" in content
        # Collapsible <details> block contains the full dump.
        assert "<details>" in content
        assert "Full manifest params" in content
        # `create_masks: false` NOT in highlights but IS in the full
        # dump. `content.split("<details>")` yields multiple details
        # blocks (register + warp); locate the register block by
        # heading position, then slice from there.
        reg_heading = content.index("<h2>Registration parameters</h2>")
        reg_details_open = content.index("<details>", reg_heading)
        reg_details_close = content.index("</details>", reg_details_open)
        register_details_block = content[reg_details_open:reg_details_close]
        assert "create_masks" in register_details_block
        assert ">false<" in register_details_block

    def test_warp_params_section_renders_from_manifest(
        self, he_root: Path,
    ):
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        assert "<h2>Warp parameters</h2>" in content
        # B2 highlights.
        for k in ("targets", "use_dask", "save_geojson"):
            assert k in content, f"expected {k} in the warp highlights"
        # Targets list renders as a compact form.
        assert "[cells, nuclei]" in content or "['cells', 'nuclei']" in content

    def test_no_params_section_when_manifest_absent(self, tmp_path: Path):
        # celltype/viz stages don't write per-run manifests today —
        # B2 defers them. Register-side manifest missing = the section
        # shows a placeholder instead of a broken layout.
        root = tmp_path / "he_registration"
        (root / "output").mkdir(parents=True)
        # Only celltyped seeded, no register (no manifest).
        _seed_stage(
            root, "celltyped", "ct_A",
            artifact_files=[
                "S1_cells_analysis.geojson",
                "S1_cells_qupath.geojson",
                "S1_nuclei_analysis.geojson",
                "S1_nuclei_qupath.geojson",
                "S1_celltyped_wholeslide.parquet",
            ],
        )
        _seed_output_symlink(root / "output", "celltyped", "ct_A")
        dest = render_summary_html(
            output_root_he=root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        assert "<h2>Registration parameters</h2>" in content
        assert "no manifest" in content  # placeholder text

    def test_lists_render_as_compact_form(self, he_root: Path):
        # Regression guard for the list-value stringifier — a Python
        # `repr([...])` would render `['cells', 'nuclei']` with quotes.
        # We want the compact `[cells, nuclei]` form (see
        # ``_stringify_param_value``).
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        # The compact form is in the highlights chip.
        assert "[cells, nuclei]" in content


# ---------------------------------------------------------------------
# 8. Dynamic overlap-figure discovery (Tracy's C2 constraint on
# ``settylab/TracyY123-nexus#15`` comment ``5338129631``). Discover
# via ``glob`` at render time; never assume a specific set of files.
# ---------------------------------------------------------------------
class TestOverlapFigures:
    def test_overlap_images_render_from_actual_disk_layout(
        self, he_root: Path,
    ):
        # Plant the four canonical VALIS overlap PNGs under the
        # promoted register run; assert every one gets an <img>.
        overlaps = he_root / "register" / "reg_A" / "overlaps"
        overlaps.mkdir(parents=True)
        for name in ("_original_overlap.png", "_rigid_overlap.png",
                     "_non_rigid_overlap.png", "_micro_reg.png"):
            (overlaps / name).touch()
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()

        assert "<h2>Registration overlap figures</h2>" in content
        # Each file gets a linked thumbnail via a relative path.
        for name in ("_original_overlap.png", "_rigid_overlap.png",
                     "_non_rigid_overlap.png", "_micro_reg.png"):
            assert f"../register/reg_A/overlaps/{name}" in content
        # Labels come from _OVERLAP_LABELS.
        assert "Original (before registration)" in content
        assert "After rigid solve" in content
        assert "After non-rigid solve" in content
        assert "After micro solve" in content

    def test_only_available_images_shown(self, he_root: Path):
        # rigid_only mode produces only two overlays; C2 says the
        # HTML must show ONLY those two, not fabricate the missing pair.
        overlaps = he_root / "register" / "reg_A" / "overlaps"
        overlaps.mkdir(parents=True)
        (overlaps / "_original_overlap.png").touch()
        (overlaps / "_rigid_overlap.png").touch()
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        assert "_original_overlap.png" in content
        assert "_rigid_overlap.png" in content
        # These MUST be absent — no fake references.
        assert "_non_rigid_overlap.png" not in content
        assert "_micro_reg.png" not in content
        assert "After non-rigid solve" not in content
        assert "After micro solve" not in content

    def test_priority_order_stable_across_filesystem_orderings(
        self, he_root: Path,
    ):
        # Filesystem readdir order is unspecified; priority sort must
        # produce a deterministic display order (rigid before nonrigid,
        # etc.). Plant in reverse alphabetical to force the point.
        overlaps = he_root / "register" / "reg_A" / "overlaps"
        overlaps.mkdir(parents=True)
        for name in ("_rigid_overlap.png", "_original_overlap.png",
                     "_non_rigid_overlap.png", "_micro_reg.png"):
            (overlaps / name).touch()
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        positions = {
            n: content.index(n)
            for n in ("_original_overlap.png", "_rigid_overlap.png",
                      "_non_rigid_overlap.png", "_micro_reg.png")
        }
        # Refinement order: original → rigid → non_rigid → micro.
        assert positions["_original_overlap.png"] < positions["_rigid_overlap.png"]
        assert positions["_rigid_overlap.png"] < positions["_non_rigid_overlap.png"]
        assert positions["_non_rigid_overlap.png"] < positions["_micro_reg.png"]

    def test_unknown_filename_falls_to_alphabetical_after_known(
        self, he_root: Path,
    ):
        # A future VALIS release might add a new overlap type; C2
        # requires we still surface it (dynamic discovery is source
        # of truth), just after the known ones in alphabetical order.
        overlaps = he_root / "register" / "reg_A" / "overlaps"
        overlaps.mkdir(parents=True)
        (overlaps / "_original_overlap.png").touch()
        (overlaps / "_rigid_overlap.png").touch()
        (overlaps / "aaa_unknown_extra.png").touch()  # sorts BEFORE known if pure alpha
        (overlaps / "zzz_unknown_extra.png").touch()
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        # Known-first: originals + rigid render before the unknowns.
        p_orig = content.index("_original_overlap.png")
        p_rigid = content.index("_rigid_overlap.png")
        p_aaa = content.index("aaa_unknown_extra.png")
        p_zzz = content.index("zzz_unknown_extra.png")
        assert p_orig < p_aaa
        assert p_rigid < p_aaa  # known priority beats alpha ordering
        assert p_aaa < p_zzz    # alphabetical among unknowns
        # Unknown filename renders bare (no fake label).
        assert "aaa_unknown_extra.png" in content

    def test_overlap_section_dynamic_when_no_images(self, he_root: Path):
        # register promoted, but the overlaps/ dir is empty (or
        # doesn't exist). Section must render the "no overlap
        # diagnostics" line rather than an empty/broken block.
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        assert "<h2>Registration overlap figures</h2>" in content
        assert "No overlap diagnostics on disk" in content
        # No dangling <img> tag from a phantom overlay.
        assert "src='../register/reg_A/overlaps/" not in content

    def test_no_promoted_register_shows_placeholder(self, tmp_path: Path):
        # Fresh he_root with no promoted register at all; the section
        # must not blow up on the missing symlink.
        root = tmp_path / "he_registration"
        (root / "output").mkdir(parents=True)
        _seed_stage(
            root, "viz", "viz_A",
            artifact_files=["S1_overlay.png"],
        )
        _seed_output_symlink(root / "output", "viz", "viz_A")
        dest = render_summary_html(
            output_root_he=root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        assert "<h2>Registration overlap figures</h2>" in content
        assert "no promoted register run" in content.lower()


# ---------------------------------------------------------------------
# 4. No overlap-quality section — regression guard for the deferred
# metric staying entirely absent.
# ---------------------------------------------------------------------
class TestOverlapQualityDeferred:
    def test_no_overlap_quality_section_present(self, he_root: Path):
        # Tracy's spec: "no placeholder needed until the metric is
        # defined." The HTML must not contain the section header or
        # a stub column — every mention of that string comes back to
        # THIS test's search token.
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text().lower()
        assert "overlap-quality" not in content
        assert "overlap quality" not in content
        assert "iou" not in content
        assert "ssim" not in content

    def test_column_count_matches_documented_fields(self, he_root: Path):
        # Per-stage table headers are exactly: Stage, he_job_id,
        # params_hash, Timestamp, git_sha, Flags. Six columns — no
        # overlap-quality slot. Ties the deferred spec to a numeric
        # column count that changes noisily if someone adds a stub.
        dest = render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        content = dest.read_text()
        header_block = content.split("<thead>")[1].split("</thead>")[0]
        assert header_block.count("<th>") == 6


# ---------------------------------------------------------------------
# 5. Atomicity + determinism.
# ---------------------------------------------------------------------
class TestAtomicWriteAndDeterminism:
    def test_output_path_deterministic_and_atomic(self, he_root: Path):
        # First write produces summary.html; second write cleanly
        # replaces it with no stray .tmp.<pid> residue.
        render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        render_summary_html(
            output_root_he=he_root, sample_id="S1", run_id="demo_v1",
        )
        assert (he_root / "output" / "summary.html").exists()
        stray = sorted(
            p.name for p in (he_root / "output").iterdir()
            if p.name.startswith("summary.html.tmp.")
        )
        assert stray == []

    def test_write_survives_crash_between_stages(
        self, he_root: Path, monkeypatch,
    ):
        # Simulate an exception AFTER the tmp is written but before
        # os.replace. The live summary.html must be either the prior
        # state or absent — never a half-written new file.
        from hexenium import summary_html as _sh
        original_replace = os.replace
        def crashing_replace(src, dst):
            if str(dst).endswith("/summary.html"):
                # Simulate a rename fault (disk-full / permission).
                raise OSError("simulated post-tmp fault")
            return original_replace(src, dst)
        monkeypatch.setattr(_sh.os, "replace", crashing_replace)
        with pytest.raises(OSError, match="simulated post-tmp"):
            _sh.render_summary_html(
                output_root_he=he_root, sample_id="S1", run_id="demo_v1",
            )
        # summary.html was never present, so shouldn't be now.
        assert not (he_root / "output" / "summary.html").exists()


# ---------------------------------------------------------------------
# 6. Both promotion paths — the renderer is reached from
# --set-default-on-success AND from hexenium set-default-run.
# ---------------------------------------------------------------------
class TestBothPromotionPaths:
    def test_promote_hook_calls_renderer_on_success(self, tmp_path: Path):
        # The promote hook writes to register/<layout.he_job_id>/ and
        # warp/<layout.he_job_id>/ — both stages share the fresh id.
        # Seed the on-disk layout matching what a real invocation
        # would have produced.
        he_root = tmp_path / "he_registration"
        (he_root / "output").mkdir(parents=True)
        _seed_stage(
            he_root, "register", "job_new",
            manifest_extra={"params": {"mode": "rigid_only"}},
            artifact_files=["data/_registrar.pickle"],
        )
        _seed_stage(
            he_root, "warp", "job_new",
            manifest_extra={
                "params": {"targets": ["cells"], "use_dask": True,
                           "save_geojson": True},
                "source_register_run_id": "job_new",
            },
            artifact_files=["he_cell_seg.parquet", "he_nucleus_seg.parquet"],
        )
        from hexenium import pipeline
        from hexenium.layout import RunLayout
        layout = RunLayout(
            sample_id="S1", he_job_id="job_new", output_root_he=he_root,
            integrated=False, stages_suffix="",
        )
        # Short-circuit path: no promotable stages_ran → no side effects.
        pipeline._promote_completed_run(
            layout=layout, cfg={"set_default_on_success": True},
            stages_ran=[], source_register_run_id=None,
        )
        assert not (he_root / "output" / "summary.html").exists()

        # Real promotion: register + warp both under the same job_new.
        pipeline._promote_completed_run(
            layout=layout, cfg={"set_default_on_success": True},
            stages_ran=["register", "warp"],
            source_register_run_id="job_new",
        )
        assert (he_root / "output" / "summary.html").exists()

    def test_cli_set_default_run_writes_summary_html(
        self, tmp_path: Path, capsys,
    ):
        # Build a minimal integrated-by-run-id layout under an
        # <output_root>/<sample>/<sample>_<run>/he_registration/ tree
        # so `hexenium set-default-run --sample-id ... --run-id ...`
        # resolves the right directory.
        output_root = tmp_path / "runs"
        he_root = (output_root / "S1" / "S1_demo_v1" / "he_registration")
        (he_root / "output").mkdir(parents=True)
        _seed_stage(
            he_root, "register", "reg_A",
            manifest_extra={"params": {"mode": "rigid_only"}},
            artifact_files=["data/_registrar.pickle"],
        )
        _seed_stage(
            he_root, "warp", "warp_A",
            manifest_extra={
                "params": {"targets": ["cells"], "use_dask": True,
                           "save_geojson": True},
                "source_register_run_id": "reg_A",
            },
            artifact_files=["he_cell_seg.parquet", "he_nucleus_seg.parquet"],
        )
        from hexenium.cli import main
        rc = main([
            "set-default-run",
            "--sample-id", "S1",
            "--run-id", "demo_v1",
            "--output-root", str(output_root),
            "--register-run-id", "reg_A",
            "--warp-run-id", "warp_A",
        ])
        assert rc == 0
        assert (he_root / "output" / "summary.html").exists()

    def test_promote_hook_survives_summary_render_failure(
        self, tmp_path: Path, monkeypatch,
    ):
        # If the renderer explodes AFTER a successful symlink
        # promotion, the hook must log-and-continue rather than
        # rolling back the symlinks (they're the source of truth).
        he_root = tmp_path / "he_registration"
        (he_root / "output").mkdir(parents=True)
        _seed_stage(
            he_root, "register", "job_new",
            manifest_extra={"params": {"mode": "rigid_only"}},
            artifact_files=["data/_registrar.pickle"],
        )
        _seed_stage(
            he_root, "warp", "job_new",
            manifest_extra={
                "params": {"targets": ["cells"], "use_dask": True,
                           "save_geojson": True},
                "source_register_run_id": "job_new",
            },
            artifact_files=["he_cell_seg.parquet", "he_nucleus_seg.parquet"],
        )
        from hexenium import pipeline
        from hexenium.layout import RunLayout

        def exploding(**kw):
            raise RuntimeError("simulated summary crash")
        monkeypatch.setattr(
            "hexenium.summary_html.render_summary_html", exploding,
        )
        layout = RunLayout(
            sample_id="S1", he_job_id="job_new", output_root_he=he_root,
            integrated=False, stages_suffix="",
        )
        pipeline._promote_completed_run(
            layout=layout, cfg={"set_default_on_success": True},
            stages_ran=["register", "warp"],
            source_register_run_id="job_new",
        )
        # Symlinks committed atomically even though summary render blew up.
        assert (he_root / "output" / "register").is_symlink()
        assert (he_root / "output" / "warp").is_symlink()
        # summary.html was NOT written (the exploding mock intercepted).
        assert not (he_root / "output" / "summary.html").exists()
