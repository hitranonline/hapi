from __future__ import annotations

import math

import numpy as np

from .models import GasCase, SpectralWindow, SpectrumResult
from .spectra import BAR_PER_ATM, build_grid, cross_section_to_absorbance, doppler_hwhm_cm, voigt_profile_cm


T_REF_K = 296.0


def render_cross_section_from_lines(
    grid: np.ndarray,
    *,
    line_centers: np.ndarray,
    line_intensities: np.ndarray,
    mass_da: float,
    temperature_k: float,
    pressure_torr: float,
    gamma0: float,
    n_exponent: float,
    line_cutoff: float,
) -> np.ndarray:
    centers = np.asarray(line_centers, dtype=np.float64)
    intensities = np.asarray(line_intensities, dtype=np.float64)
    if centers.shape != intensities.shape:
        raise ValueError("line_centers and line_intensities must have the same shape")
    if line_cutoff <= 0.0:
        raise ValueError("line_cutoff must be positive")
    if temperature_k <= 0.0:
        raise ValueError("temperature_k must be positive")
    if pressure_torr < 0.0:
        raise ValueError("pressure_torr must be non-negative")
    if mass_da <= 0.0:
        raise ValueError("mass_da must be positive")

    pressure_bar = pressure_torr / 760.0 * BAR_PER_ATM
    lorentz_hwhm = gamma0 * pressure_bar * (T_REF_K / temperature_k) ** n_exponent
    doppler_hwhm = doppler_hwhm_cm(centers, temperature_k=temperature_k, mass_da=mass_da)

    cross_section = np.zeros_like(grid, dtype=np.float64)
    for center, strength, alpha_d in zip(centers, intensities, doppler_hwhm):
        local_cutoff = max(line_cutoff, 25.0 * max(alpha_d, lorentz_hwhm))
        left = np.searchsorted(grid, center - local_cutoff, side="left")
        right = np.searchsorted(grid, center + local_cutoff, side="right")
        if left == right:
            continue
        profile = voigt_profile_cm(
            grid=grid[left:right],
            center_cm=float(center),
            doppler_hwhm_cm_value=float(alpha_d),
            lorentz_hwhm_cm_value=float(lorentz_hwhm),
        )
        cross_section[left:right] += float(strength) * profile
    return cross_section


def render_absorbance_on_grid(
    grid: np.ndarray,
    *,
    case: GasCase,
    line_centers: np.ndarray,
    line_intensities: np.ndarray,
    mass_da: float,
    gamma0: float,
    n_exponent: float,
    line_cutoff: float,
) -> np.ndarray:
    cross_section = render_cross_section_from_lines(
        grid,
        line_centers=line_centers,
        line_intensities=line_intensities,
        mass_da=mass_da,
        temperature_k=case.temperature_k,
        pressure_torr=case.pressure_torr,
        gamma0=gamma0,
        n_exponent=n_exponent,
        line_cutoff=line_cutoff,
    )
    return cross_section_to_absorbance(cross_section, case)


def render_absorbance_from_lines(
    *,
    case: GasCase,
    window: SpectralWindow,
    line_centers: np.ndarray,
    line_intensities: np.ndarray,
    mass_da: float,
    gamma0: float,
    n_exponent: float,
    line_cutoff: float,
) -> SpectrumResult:
    grid = build_grid(window)
    absorbance = render_absorbance_on_grid(
        grid,
        case=case,
        line_centers=line_centers,
        line_intensities=line_intensities,
        mass_da=mass_da,
        gamma0=gamma0,
        n_exponent=n_exponent,
        line_cutoff=line_cutoff,
    )
    return SpectrumResult(
        wavenumber=grid,
        values=absorbance,
        quantity="absorbance",
        metadata={
            "line_count": int(np.asarray(line_centers).size),
            "mass_da": mass_da,
            "gamma0": gamma0,
            "n_exponent": n_exponent,
            "line_cutoff": line_cutoff,
        },
    )
