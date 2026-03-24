from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class GasCase:
    temperature_k: float
    pressure_torr: float
    mole_fraction: float
    path_length_cm: float

    def __post_init__(self) -> None:
        if self.temperature_k <= 0.0:
            raise ValueError("temperature_k must be positive")
        if self.pressure_torr < 0.0:
            raise ValueError("pressure_torr must be non-negative")
        if not 0.0 <= self.mole_fraction <= 1.0:
            raise ValueError("mole_fraction must be between 0 and 1")
        if self.path_length_cm < 0.0:
            raise ValueError("path_length_cm must be non-negative")

    @property
    def pressure_atm(self) -> float:
        return self.pressure_torr / 760.0


@dataclass(frozen=True)
class SpectralWindow:
    wn_min: float
    wn_max: float
    wn_step: float

    def __post_init__(self) -> None:
        if self.wn_max <= self.wn_min:
            raise ValueError("wn_max must be greater than wn_min")
        if self.wn_step <= 0.0:
            raise ValueError("wn_step must be positive")


@dataclass(frozen=True)
class DatasetPaths:
    root_dir: Path
    hitran_db_dir: Path
    exomol_db_dir: Path
    artifacts_dir: Path


@dataclass(frozen=True)
class BandSelection:
    species: str
    mode_label: str
    lower_band: str | None = None
    upper_band: str | None = None
    require_same_other_modes: bool = True
    require_unit_step: bool = False


@dataclass
class SpectrumResult:
    wavenumber: np.ndarray
    values: np.ndarray
    quantity: str
    metadata: dict[str, Any] = field(default_factory=dict)
    csv_path: Path | None = None
    html_path: Path | None = None


@dataclass
class BandExportResult:
    rows: list[dict[str, Any]]
    output_dir: Path | None = None
    manifest_path: Path | None = None
    count: int = 0
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class SummaryResult:
    rows: list[dict[str, Any]]
    csv_path: Path | None = None
    html_path: Path | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class ComparisonResult:
    band_label: str | None = None
    matched_band_label: str | None = None
    metrics: dict[str, Any] = field(default_factory=dict)
    report_path: Path | None = None
    figure_paths: list[Path] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
