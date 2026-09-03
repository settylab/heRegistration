"""Smoke test: the argparse tree assembles and `--help` runs."""
from __future__ import annotations


def test_top_level_help_runs():
    """`hexenium --help` — invoked in-process — exits 0 without exploding."""
    import pytest
    from hexenium.cli import build_parser

    parser = build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["--help"])
    assert exc.value.code == 0


def test_run_subcommand_help_runs():
    """`hexenium run --help` also exits 0 (dodges an argparse subparser
    regression that we've had before)."""
    import pytest
    from hexenium.cli import build_parser

    parser = build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["run", "--help"])
    assert exc.value.code == 0


def test_run_subcommand_parses_minimal_args():
    """The full CLI surface accepts a minimal invocation without raising
    during arg-parsing (validation happens later, in _resolve_config)."""
    from hexenium.cli import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "run",
        "--sample-id", "SAMPLE1",
        "--he-path", "/tmp/HE.ome.tif",
        "--xenium-bundle", "/tmp/output-XETG...",
    ])
    assert args.cmd == "run"
    assert args.sample_id == "SAMPLE1"


def test_run_id_persisted_in_cfg_when_h5ad_also_present(tmp_path):
    """``--run-id`` must land in ``cfg["run_id"]`` even when
    ``--xenium-h5ad`` is also supplied — otherwise ``pipeline.py``
    forwards ``run_id=None`` to :func:`hexenium.layout.resolve_layout`,
    the trio-precedence branch is skipped, and outputs fall back to the
    h5ad-derived location. This is the exact regression fixed by
    hoisting the persistence assignment out of the h5ad-derivation
    block in :func:`hexenium.cli._resolve_config`.
    """
    from hexenium.cli import _resolve_config, build_parser

    # A path that just needs to LOOK like a real h5ad to argparse
    # + _resolve_config (which doesn't check existence — layout does).
    h5ad = tmp_path / "xenium_ranger.h5ad"
    h5ad.write_bytes(b"")

    parser = build_parser()
    args = parser.parse_args([
        "run",
        "--sample-id", "SAMPLE1",
        "--output-root", str(tmp_path),
        "--run-id", "downstream_v1",
        "--xenium-h5ad", str(h5ad),
        "--he-path", str(tmp_path / "HE.ome.tif"),
        "--xenium-bundle", str(tmp_path),
    ])
    cfg = _resolve_config(args)
    assert cfg["run_id"] == "downstream_v1", (
        "regression: --run-id was dropped from cfg because "
        "--xenium-h5ad was also supplied. pipeline.py would then "
        "call resolve_layout with run_id=None, silently bypassing "
        "the trio-precedence branch of the output-root fix."
    )
    # And the h5ad passes through unchanged as an input reference.
    assert cfg["xenium_h5ad"] == str(h5ad)


def test_run_id_still_persisted_without_h5ad(tmp_path):
    """The existing "integrated-by-run-id" path (no h5ad supplied)
    keeps working — ``cfg["run_id"]`` is populated and the h5ad path
    is derived from the canonical upstream layout when it exists on
    disk (or omitted when it doesn't and celltype is out of scope).
    """
    from hexenium.cli import _resolve_config, build_parser

    parser = build_parser()
    args = parser.parse_args([
        "run",
        "--sample-id", "SAMPLE1",
        "--output-root", str(tmp_path),
        "--run-id", "run42",
        "--he-path", str(tmp_path / "HE.ome.tif"),
        "--xenium-bundle", str(tmp_path),
        "--stages", "register", "warp",   # skip celltype so derived-missing is fine
    ])
    cfg = _resolve_config(args)
    assert cfg["run_id"] == "run42"


def test_run_subcommand_accepts_register_run_id():
    """``--register-run-id`` is the flag Tracy uses to compare multiple
    registrations. Regressions here would silently drop back to the
    default-symlink registration and confound downstream diffs."""
    from hexenium.cli import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "run",
        "--sample-id", "SAMPLE1",
        "--he-path", "/tmp/HE.ome.tif",
        "--xenium-bundle", "/tmp/output-XETG...",
        "--register-run-id", "reg_abc123",
    ])
    assert args.register_run_id == "reg_abc123"


def test_set_default_run_subcommand_help_runs():
    """``hexenium set-default-run --help`` — a smoke on the second
    subparser Tracy uses to pick a run after visual QA."""
    import pytest
    from hexenium.cli import build_parser

    parser = build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["set-default-run", "--help"])
    assert exc.value.code == 0


def test_set_default_run_parses_minimal_args():
    from hexenium.cli import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "set-default-run",
        "--sample-id", "SAMPLE1",
        "--output-root", "/tmp/out",
        "--register-run-id", "reg_A",
    ])
    assert args.cmd == "set-default-run"
    assert args.sample_id == "SAMPLE1"
    assert args.register_run_id == "reg_A"
    assert args.warp_run_id is None
    assert args.force_lineage is False


def test_version_flag_reports_package_version():
    import pytest
    from hexenium import __version__
    from hexenium.cli import build_parser

    parser = build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["--version"])
    # argparse --version exits with code 0.
    assert exc.value.code == 0
    assert __version__  # non-empty
