"""Non-equilibrium vibrational distributions for CH₄ ν₃.

Effective single-mode approximation for the CH₄ ν₃ ladder (0,0,v₃,0),
with the other vibrational modes held at v = 0. This is not a full polyatomic
non-LTE model. Treanor distribution follows Fridman (2008), Plasma Chemistry,
§3.1.8, Eq. (3-37).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

SECOND_RADIATION_CONSTANT_CM_K = 1.438776877


@dataclass(frozen=True)
class CH4Nu3Constants:
    """Spectroscopic constants for the CH₄ ν₃ anharmonic oscillator."""

    omega_3: float
    x_33: float
    label: str = ""

    @property
    def x_e(self) -> float:
        return self.x_33 / self.omega_3

    @property
    def theta_v(self) -> float:
        return SECOND_RADIATION_CONSTANT_CM_K * self.omega_3


# ω₃ ≈ 3019 cm⁻¹, ω₃x₃₃ ≈ 62.2 cm⁻¹ (Herzberg 1945; HITRAN CH₄ ν₃ mode)
CH4_NU3 = CH4Nu3Constants(omega_3=3019.0, x_33=62.2, label="CH₄ ν₃")


def energy_v3(v: np.ndarray | float, mol: CH4Nu3Constants) -> np.ndarray:
    """CH₄ ν₃ energy E_v in temperature units [K]."""

    v = np.asarray(v, dtype=np.float64)
    return SECOND_RADIATION_CONSTANT_CM_K * (mol.omega_3 * v - mol.x_33 * v * (v + 1))


def ln_treanor_nu3(
    v: np.ndarray | float,
    Tv: float,
    T0: float,
    mol: CH4Nu3Constants,
) -> np.ndarray:
    """Natural log of the CH₄ ν₃ Treanor distribution."""

    v = np.asarray(v, dtype=np.float64)
    return -v * mol.theta_v / Tv + mol.x_e * mol.theta_v * v**2 / T0


def treanor_nu3(
    v: np.ndarray | float,
    Tv: float,
    T0: float,
    mol: CH4Nu3Constants,
) -> np.ndarray:
    """CH₄ ν₃ Treanor population ratio N_v / N_0."""

    return np.exp(ln_treanor_nu3(v, Tv, T0, mol))


def ln_boltzmann_nu3(
    v: np.ndarray | float,
    T: float,
    mol: CH4Nu3Constants,
) -> np.ndarray:
    """Natural log of the Boltzmann distribution for CH₄ ν₃.

    Uses compact v² energy form consistent with Treanor normalization; ensures
    nonlte_intensity_scale_factor = 1.0 when Tv = T0.
    """

    v = np.asarray(v, dtype=np.float64)
    return -v * mol.theta_v / T + mol.x_e * mol.theta_v * v**2 / T


def boltzmann_nu3(
    v: np.ndarray | float,
    T: float,
    mol: CH4Nu3Constants,
) -> np.ndarray:
    """CH₄ ν₃ Boltzmann population ratio N_v / N_0."""

    return np.exp(ln_boltzmann_nu3(v, T, mol))


def treanor_minimum_nu3(Tv: float, T0: float, mol: CH4Nu3Constants) -> float:
    """Continuous Treanor minimum v* for CH₄ ν₃."""

    return T0 / (2 * mol.x_e * Tv)


def nonlte_intensity_scale_factor(
    v3_lower: np.ndarray | float,
    temperature_k: float,
    vibrational_temperature_k: float,
    mol: CH4Nu3Constants,
) -> np.ndarray:
    """Non-LTE intensity scale factor for CH₄ ν₃ transitions.

    Uses compact v² Boltzmann form consistent with Treanor normalization;
    ensures equilibrium identity when vibrational_temperature_k equals
    temperature_k.
    """

    v3_lower = np.asarray(v3_lower, dtype=np.float64)
    return treanor_nu3(
        v3_lower, vibrational_temperature_k, temperature_k, mol
    ) / boltzmann_nu3(v3_lower, temperature_k, mol)
