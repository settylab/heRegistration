"""User-facing tools that consume hexenium pipeline outputs.

Currently exposes :mod:`~hexenium.tools.roi_subset` — the 2×2 workflow
for subsetting Xenium data by user-drawn ROI polygons. A1 (lightweight
registration-input subset) and B1 (H&E-space sjoin) land in Phase 2a;
A2 (forward VALIS warp) follows in the same phase; B2 (inverse warp)
lands in Phase 2b.
"""
from __future__ import annotations
