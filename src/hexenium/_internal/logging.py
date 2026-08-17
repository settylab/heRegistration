"""Minimal logging helper for hexenium.

`log(msg)` prints a timestamped, flushed line to stdout. Use it for
every per-stage progress line so the log file shows where time goes
and (critically) so output never sits in stdio buffers when Slurm
SIGKILLs a job at the time limit.

Belt + braces vs `python -u` / PYTHONUNBUFFERED=1: each `log()` call
also explicitly flushes. If a downstream library reopens stdout, this
still works.
"""
from __future__ import annotations

import sys
import time


def log(msg: str) -> None:
    """Print msg with an ISO-8601 timestamp prefix, flushed."""
    ts = time.strftime("%Y-%m-%dT%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)
    try:
        sys.stdout.flush()
    except Exception:
        pass


def banner(title: str) -> None:
    """Heavy visual separator + log line. Use for stage boundaries."""
    log("=" * 64)
    log(title)
    log("=" * 64)
