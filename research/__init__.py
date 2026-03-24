from __future__ import annotations

from importlib import import_module

from .models import (
    BandExportResult,
    BandSelection,
    ComparisonResult,
    DatasetPaths,
    GasCase,
    SpectralWindow,
    SpectrumResult,
    SummaryResult,
)


_LAZY_MODULES = {"absorbance", "bands", "compare", "exomol", "hitran", "io", "models", "spectra"}


def __getattr__(name: str):
    if name in _LAZY_MODULES:
        module = import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "bands",
    "compare",
    "absorbance",
    "exomol",
    "hitran",
    "io",
    "models",
    "spectra",
    "BandExportResult",
    "BandSelection",
    "ComparisonResult",
    "DatasetPaths",
    "GasCase",
    "SpectralWindow",
    "SpectrumResult",
    "SummaryResult",
]
