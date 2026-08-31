"""Config dataclasses for the ROI-subset pipeline.

The config is YAML-first and mode-scoped: ``mode: a1`` for the
lightweight registration-input subset. ``b1`` and ``a2`` share the
same top-level schema with different required fields (populated as
each mode's implementation lands). ``b2`` deliberately absent from
this phase.

Frame declarations are **explicit** — every ROI carries ``frame:
microns | dapi_px | he_px``. There is no auto-detection and no default,
per the corrections addressed on ``#29``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Union

import yaml


Mode = Literal["a1", "b1", "a2"]
GeometryFrame = Literal["global", "local"]


@dataclass
class RoiSpec:
    """One user-declared ROI."""
    name: str
    path: Path
    frame: str            # microns | dapi_px | he_px (validated at run time)
    source_tool: str      # qupath_csv | qupath_geojson | xenium_explorer
    selection: str | None = None
    """When the file carries multiple features, pick this one by name.
    ``None`` = first feature wins."""


@dataclass
class Outputs:
    dapi_crop: bool = True
    dapi_polygon_mask: bool = False
    he_crop: bool = True             # only honored by A2
    he_polygon_mask: bool = False    # only honored by A2
    cells_parquet: bool = True
    cell_boundaries: bool = True
    nucleus_boundaries: bool = True
    experiment_xenium: bool = True   # A1: verbatim copy
    geometry_frame: GeometryFrame = "global"


@dataclass
class WarpKnobs:
    """Only meaningful for A2 (forward warp)."""
    densify_max_segment_length_dapi_px: float = 10.0
    """~2.125 µm at the default pixel_size_morph. Matches the deformation-
    field resolution used at register time
    (``max_non_rigid_registration_dim_px=10000`` downsampled ~4×)."""
    make_valid: bool = True
    keep_largest_component: bool = True


@dataclass
class SubsetConfig:
    mode: Mode
    xenium_bundle: Path
    output_dir: Path
    rois: list[RoiSpec] = field(default_factory=list)
    outputs: Outputs = field(default_factory=Outputs)
    warp: WarpKnobs = field(default_factory=WarpKnobs)
    he_bundle: Path | None = None            # required for A2 (H&E crop)
    registrar_pickle: Path | None = None     # required for A2 (forward warp)
    pixel_size_morph: float | None = None
    """Override for pixel_size_morph. When None, read from
    ``<xenium_bundle>/experiment.xenium``; falls back to 0.2125."""

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> "SubsetConfig":
        with Path(path).open() as f:
            raw = yaml.safe_load(f)
        return cls.from_dict(raw)

    @classmethod
    def from_dict(cls, raw: dict) -> "SubsetConfig":
        rois = [
            RoiSpec(
                name=str(r["name"]),
                path=Path(r["path"]),
                frame=str(r["frame"]),
                source_tool=str(r["source_tool"]),
                selection=r.get("selection"),
            )
            for r in raw.get("rois", [])
        ]
        outputs_raw = raw.get("outputs", {}) or {}
        warp_raw = raw.get("warp", {}) or {}
        cfg = cls(
            mode=str(raw["mode"]),                       # type: ignore[arg-type]
            xenium_bundle=Path(raw["xenium_bundle"]),
            output_dir=Path(raw["output_dir"]),
            rois=rois,
            outputs=Outputs(**outputs_raw),
            warp=WarpKnobs(**warp_raw),
            he_bundle=Path(raw["he_bundle"]) if raw.get("he_bundle") else None,
            registrar_pickle=(
                Path(raw["registrar_pickle"]) if raw.get("registrar_pickle") else None
            ),
            pixel_size_morph=raw.get("pixel_size_morph"),
        )
        cfg.validate()
        return cfg

    def validate(self) -> None:
        if self.mode not in ("a1", "b1", "a2"):
            raise ValueError(f"mode must be one of a1|b1|a2, got {self.mode!r}")
        if len(self.rois) == 0:
            raise ValueError("at least one ROI is required")
        for r in self.rois:
            if r.frame not in ("microns", "dapi_px", "he_px"):
                raise ValueError(
                    f"ROI {r.name!r}: frame must be microns|dapi_px|he_px, "
                    f"got {r.frame!r} (no default — declare explicitly)"
                )
            if r.source_tool not in ("qupath_csv", "qupath_geojson", "xenium_explorer"):
                raise ValueError(
                    f"ROI {r.name!r}: unknown source_tool={r.source_tool!r}"
                )
        if self.mode == "a1":
            for r in self.rois:
                if r.frame == "he_px":
                    raise ValueError(
                        f"A1: ROI {r.name!r} has frame='he_px', which requires the "
                        "VALIS registrar. Use A2/B2, or convert to microns first."
                    )
        if self.mode == "a2":
            if self.registrar_pickle is None:
                raise ValueError("A2 requires 'registrar_pickle' in config")
            if self.outputs.he_crop and self.he_bundle is None:
                raise ValueError(
                    "A2 with outputs.he_crop=true requires 'he_bundle' in config"
                )
