#!/usr/bin/env python3
"""Thin CLI for :func:`hexenium.tools.roi_subset.run_a1` (and future b1/a2).

Reads a YAML config that declares mode, xenium_bundle, output_dir, and
one or more ROIs with explicit ``frame`` and ``source_tool``. See
``src/hexenium/tools/roi_subset.py`` for the entry-point docstring and
``src/hexenium/_internal/roi_config.py`` for the schema.

Example
-------

.. code-block:: bash

   python scripts/subset_by_roi.py --config configs/tma111a_a1.yaml

The script does not open a Slurm job or spin up a JVM — A1 is a pure
CPU/IO pass. A2 will grow those requirements in a follow-up.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True,
                        help="Path to the YAML subset config.")
    parser.add_argument("--log-level", default="INFO",
                        choices=("DEBUG", "INFO", "WARNING", "ERROR"))
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    # Import here so a bad config error surfaces before the heavy imports.
    from hexenium._internal.roi_config import SubsetConfig
    from hexenium.tools.roi_subset import run_a1

    config = SubsetConfig.from_yaml(args.config)
    if config.mode == "a1":
        dirs = run_a1(config)
    else:
        raise SystemExit(
            f"mode={config.mode!r} not yet implemented; only A1 is available "
            "in this branch. b1 / a2 land in follow-up commits."
        )
    print(f"[subset_by_roi] wrote {len(dirs)} ROI dir(s) under {config.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
