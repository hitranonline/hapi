"""ExoMol non-LTE intensity scaling for CH₄ ν₃.

This module applies an effective single-mode CH₄ ν₃ Treanor population
approximation to LTE ExoMol line intensities by ratio scaling of each grouped
transition set. It is not a full polyatomic non-LTE model.
"""

from __future__ import annotations

from pathlib import Path
from typing import TypedDict, cast

from .ch4_treanor import CH4_NU3, CH4Nu3Constants, nonlte_intensity_scale_factor
from .exomol import collect_nu3_transitions_by_jpair

TransitionKey = tuple[int, int, int, int]


class TransitionGroup(TypedDict):
    wavenumber: list[float]
    intensity: list[float]


TransitionGroups = dict[TransitionKey, TransitionGroup]
RunMetadata = dict[str, object]


def collect_nu3_transitions_nonlte(
    *,
    data_dir: Path,
    temperature_k: float,
    vibrational_temperature_k: float,
    mol: CH4Nu3Constants = CH4_NU3,
    wn_min: float = 2500.0,
    wn_max: float = 3500.0,
    wing: float = 0.5,
    intensity_threshold: float = 0.0,
) -> tuple[TransitionGroups, RunMetadata]:
    """Collect CH₄ ν₃ ExoMol transitions with Treanor-scaled intensities."""

    grouped_lte_raw, run_metadata = collect_nu3_transitions_by_jpair(
        data_dir=data_dir,
        temperature_k=temperature_k,
        wn_min=wn_min,
        wn_max=wn_max,
        wing=wing,
        intensity_threshold=intensity_threshold,
    )
    grouped_lte = cast(TransitionGroups, grouped_lte_raw)

    grouped_nonlte: TransitionGroups = {}
    for key, group in grouped_lte.items():
        lower_n3, _upper_n3, _lower_j, _upper_j = key
        scale = float(
            nonlte_intensity_scale_factor(
                lower_n3,
                temperature_k,
                vibrational_temperature_k,
                mol,
            )
        )
        grouped_nonlte[key] = {
            "wavenumber": [float(value) for value in group["wavenumber"]],
            "intensity": [float(value) * scale for value in group["intensity"]],
        }

    return grouped_nonlte, run_metadata
