from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from scipy.special import wofz

from .io import ensure_directory
from .models import GasCase, SpectralWindow


AMU_TO_KG = 1.66053906660e-27
BAR_PER_ATM = 1.01325
BOLTZMANN_J_K = 1.380649e-23
LIGHT_SPEED_M_S = 2.99792458e8
SQRT_2PI = math.sqrt(2.0 * math.pi)
SQRT_2LN2 = math.sqrt(2.0 * math.log(2.0))


def build_grid(window: SpectralWindow) -> np.ndarray:
    return np.arange(window.wn_min, window.wn_max + 0.5 * window.wn_step, window.wn_step, dtype=np.float64)


def to_absorbance(transmittance: np.ndarray) -> np.ndarray:
    clipped = np.clip(transmittance, 1.0e-300, None)
    return -np.log(clipped)


def number_density_cm3(pressure_torr: float, temperature_k: float) -> float:
    pressure_pa = pressure_torr * 133.32236842105263
    density_m3 = pressure_pa / (BOLTZMANN_J_K * temperature_k)
    return density_m3 / 1.0e6


def cross_section_to_absorbance(cross_section: np.ndarray, case: GasCase) -> np.ndarray:
    absorber_density_cm3 = number_density_cm3(case.pressure_torr, case.temperature_k) * case.mole_fraction
    return cross_section * absorber_density_cm3 * case.path_length_cm


def doppler_hwhm_cm(wavenumbers: np.ndarray, temperature_k: float, mass_da: float) -> np.ndarray:
    mass_kg = mass_da * AMU_TO_KG
    factor = math.sqrt(2.0 * BOLTZMANN_J_K * temperature_k * math.log(2.0) / (mass_kg * LIGHT_SPEED_M_S**2))
    return wavenumbers * factor


def voigt_profile_cm(
    grid: np.ndarray,
    center_cm: float,
    doppler_hwhm_cm_value: float,
    lorentz_hwhm_cm_value: float,
) -> np.ndarray:
    sigma = max(doppler_hwhm_cm_value / SQRT_2LN2, 1.0e-12)
    z = ((grid - center_cm) + 1j * lorentz_hwhm_cm_value) / (sigma * math.sqrt(2.0))
    return np.real(wofz(z)) / (sigma * SQRT_2PI)


def save_spectrum_csv(path: Path, wavenumber: np.ndarray, values: np.ndarray, y_label: str) -> Path:
    ensure_directory(path.parent)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["wavenumber_cm-1", y_label])
        for x_value, y_value in zip(wavenumber, values):
            writer.writerow([x_value, y_value])
    return path


def save_spectrum_html(
    path: Path,
    wavenumber: np.ndarray,
    values: np.ndarray,
    *,
    title: str,
    trace_name: str,
    y_label: str,
) -> Path:
    ensure_directory(path.parent)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=wavenumber, y=values, mode="lines", name=trace_name))
    fig.update_layout(
        title=title,
        xaxis_title="Wavenumber (cm^-1)",
        yaxis_title=y_label,
        template="plotly_white",
    )
    fig.update_xaxes(range=[float(wavenumber[0]), float(wavenumber[-1])])
    fig.write_html(path, include_plotlyjs="cdn")
    return path
