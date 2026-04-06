"""Non-equilibrium vibrational distributions.

Treanor distribution: Fridman (2008), Plasma Chemistry, §3.1.8, Eq. (3-37).
Original: Treanor, Rich & Rehm (1968), J. Chem. Phys. 48(4), 1798-1807.

The Treanor distribution describes vibrational populations when V-V exchange
dominates over V-T relaxation. Two temperatures govern the distribution:
T_v (vibrational) controls the harmonic part, T_0 (translational) controls
the anharmonic correction.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# hc/k_B conversion factor: wavenumber (cm^-1) -> temperature (K)
HC_OVER_KB = 1.4388  # cm·K


@dataclass(frozen=True)
class DiatomicConstants:
    """Spectroscopic constants for an anharmonic diatomic oscillator."""

    omega_e: float  # harmonic frequency, cm^-1
    omega_e_x_e: float  # anharmonicity constant, cm^-1
    label: str = ""

    @property
    def x_e(self) -> float:
        return self.omega_e_x_e / self.omega_e

    @property
    def theta_v(self) -> float:
        """Characteristic vibrational temperature [K]."""
        return HC_OVER_KB * self.omega_e


# Pre-built constants for common molecules (NIST values)
CO = DiatomicConstants(omega_e=2169.81, omega_e_x_e=13.29, label="CO")
N2 = DiatomicConstants(omega_e=2358.57, omega_e_x_e=14.32, label="N₂")


def energy_v(v: np.ndarray | float, mol: DiatomicConstants) -> np.ndarray:
    """Vibrational energy E_v in temperature units [K].

    Uses standard anharmonic oscillator: E_v = hc [omega_e v - omega_e x_e v(v+1)]
    """
    v = np.asarray(v, dtype=np.float64)
    return HC_OVER_KB * (mol.omega_e * v - mol.omega_e_x_e * v * (v + 1))


def ln_treanor(
    v: np.ndarray | float,
    Tv: float,
    T0: float,
    mol: DiatomicConstants,
) -> np.ndarray:
    """Natural log of the Treanor distribution — Fridman Eq. (3-37).

    Uses the compact v^2 form:  ln f = -theta_v * v / Tv + x_e * theta_v * v^2 / T0

    Equivalent to the v(v-1) derivation form; the difference is absorbed
    into the linear coefficient.
    """
    v = np.asarray(v, dtype=np.float64)
    theta = mol.theta_v
    xe = mol.x_e
    return -v * theta / Tv + xe * theta * v**2 / T0


def treanor(
    v: np.ndarray | float,
    Tv: float,
    T0: float,
    mol: DiatomicConstants,
) -> np.ndarray:
    """Treanor distribution N_v / N_0 — Fridman Eq. (3-37)."""
    return np.exp(ln_treanor(v, Tv, T0, mol))


def ln_boltzmann(
    v: np.ndarray | float,
    T: float,
    mol: DiatomicConstants,
) -> np.ndarray:
    """Natural log of Boltzmann distribution with anharmonic energy levels."""
    return -energy_v(v, mol) / T


def boltzmann(
    v: np.ndarray | float,
    T: float,
    mol: DiatomicConstants,
) -> np.ndarray:
    """Boltzmann distribution N_v / N_0 with anharmonic energy levels."""
    return np.exp(ln_boltzmann(v, T, mol))


def treanor_minimum(Tv: float, T0: float, mol: DiatomicConstants) -> float:
    """Treanor minimum quantum number — Fridman Eq. (3-38).

    Returns continuous value; the physical minimum is at the nearest integer.
    """
    return T0 / (2 * mol.x_e * Tv)
