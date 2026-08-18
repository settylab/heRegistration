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
