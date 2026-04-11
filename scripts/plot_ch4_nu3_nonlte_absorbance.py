"""Plot CH4 nu3 LTE vs non-LTE absorbance spectra side by side.

Saves a PNG comparison figure to ch4_nu3_treanor/.

Configuration constants are at the top of the file; no CLI arguments needed.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _bootstrap import ensure_repo_root

ROOT_DIR = ensure_repo_root()

from research.exomol import collect_nu3_transitions_by_jpair, render_cross_section
from research.ch4_treanor import CH4_NU3, nonlte_intensity_scale_factor
from research.spectra import build_grid, cross_section_to_absorbance
from research.models import GasCase, SpectralWindow

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
T0 = 600.0  # translational / gas temperature [K]
TV = 3000.0  # vibrational temperature for non-LTE Treanor scaling [K]
DATA_DIR = Path("exomol_db/CH4/12C-1H4/MM")
PRESSURE_TORR = 3.0
MOLE_FRACTION = 0.008
PATH_CM = 100.0
WN_MIN = 2500.0
WN_MAX = 3500.0
WN_STEP = 0.1

OUT_DIR = ROOT_DIR / "ch4_nu3_treanor"
OUT_FILE = OUT_DIR / "ch4_nu3_nonlte_absorbance.png"

# ---------------------------------------------------------------------------
# Step A — Collect LTE lines once; derive non-LTE via Treanor ratio scaling
#          (single ExoMol bz2 parse pass instead of two)
# ---------------------------------------------------------------------------
print("Collecting LTE transitions (single ExoMol parse pass)...")
grouped_lte, metadata = collect_nu3_transitions_by_jpair(
    data_dir=DATA_DIR,
    temperature_k=T0,
    wn_min=WN_MIN,
    wn_max=WN_MAX,
)

grouped_nonlte = {}
for key, group in grouped_lte.items():
    lower_n3 = key[0]
    scale = float(nonlte_intensity_scale_factor(lower_n3, T0, TV, CH4_NU3))
    grouped_nonlte[key] = {
        "wavenumber": list(group["wavenumber"]),
        "intensity": [float(v) * scale for v in group["intensity"]],
    }

# ---------------------------------------------------------------------------
# Step B — Flatten grouped lines into arrays
# ---------------------------------------------------------------------------
lte_centers = np.array([wn for g in grouped_lte.values() for wn in g["wavenumber"]])
lte_intensities = np.array([i for g in grouped_lte.values() for i in g["intensity"]])
nonlte_centers = np.array(
    [wn for g in grouped_nonlte.values() for wn in g["wavenumber"]]
)
nonlte_intensities = np.array(
    [i for g in grouped_nonlte.values() for i in g["intensity"]]
)

print(
    f"Lines: {len(lte_centers)} (positions identical; only intensities differ by Treanor scale)"
)

# ---------------------------------------------------------------------------
# Step C — Build grid and render cross-sections
# ---------------------------------------------------------------------------
window = SpectralWindow(wn_min=WN_MIN, wn_max=WN_MAX, wn_step=WN_STEP)
grid = build_grid(window)

mass_da = float(metadata["mass_da"])
gamma0 = float(metadata["gamma0"])
n_exponent = float(metadata["n_exponent"])

print("Rendering LTE cross-section...")
cs_lte = render_cross_section(
    grid,
    line_centers=lte_centers,
    line_intensities=lte_intensities,
    mass_da=mass_da,
    temperature_k=T0,
    pressure_torr=PRESSURE_TORR,
    gamma0=gamma0,
    n_exponent=n_exponent,
    line_cutoff=25.0,
)

print("Rendering non-LTE cross-section...")
cs_nonlte = render_cross_section(
    grid,
    line_centers=nonlte_centers,
    line_intensities=nonlte_intensities,
    mass_da=mass_da,
    temperature_k=T0,
    pressure_torr=PRESSURE_TORR,
    gamma0=gamma0,
    n_exponent=n_exponent,
    line_cutoff=25.0,
)

# ---------------------------------------------------------------------------
# Step D — Convert to absorbance
# ---------------------------------------------------------------------------
case = GasCase(
    temperature_k=T0,
    pressure_torr=PRESSURE_TORR,
    mole_fraction=MOLE_FRACTION,
    path_length_cm=PATH_CM,
)

abs_lte = cross_section_to_absorbance(cs_lte, case)
abs_nonlte = cross_section_to_absorbance(cs_nonlte, case)

print(f"LTE max absorbance: {abs_lte.max():.6e}")
print(f"Non-LTE max absorbance: {abs_nonlte.max():.6e}")

# ---------------------------------------------------------------------------
# Step E — Plot and save
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 5))
ax.plot(grid, abs_lte, color="steelblue", linewidth=0.5, label=f"LTE  T={T0:.0f} K")
ax.plot(
    grid,
    abs_nonlte,
    color="firebrick",
    linewidth=0.5,
    alpha=0.8,
    label=f"Non-LTE  Tv={TV:.0f} K, T₀={T0:.0f} K",
)
ax.set_xlabel("Wavenumber (cm⁻¹)")
ax.set_ylabel("Absorbance")
ax.set_title(
    f"CH₄ ν₃ LTE vs Non-LTE Absorbance"
    f"  (P={case.pressure_torr} Torr, L={case.path_length_cm} cm, x={case.mole_fraction})"
)
ax.legend()
fig.tight_layout()

OUT_DIR.mkdir(exist_ok=True)
fig.savefig(OUT_FILE, dpi=120)
print(f"Saved: {OUT_FILE}")
