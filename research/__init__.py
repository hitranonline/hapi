from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING

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

if TYPE_CHECKING:
    from . import (
        absorbance,
        bands,
        ch4_treanor,
        combined,
        compare,
        compare_cases,
        exomol,
        hitran,
        io,
        models,
        nonlte,
        spectra,
        treanor,
    )


_LAZY_MODULES = {
    "absorbance",
    "bands",
    "combined",
    "ch4_treanor",
    "compare",
    "compare_cases",
    "exomol",
    "hitran",
    "io",
    "models",
    "nonlte",
    "spectra",
    "treanor",
}


def __getattr__(name: str):
    if name in _LAZY_MODULES:
        module = import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "bands",
    "combined",
    "ch4_treanor",
    "compare",
    "compare_cases",
    "absorbance",
    "exomol",
    "hitran",
    "io",
    "models",
    "nonlte",
    "spectra",
    "treanor",
    "BandExportResult",
    "BandSelection",
    "ComparisonResult",
    "DatasetPaths",
    "GasCase",
    "SpectralWindow",
    "SpectrumResult",
    "SummaryResult",
]
