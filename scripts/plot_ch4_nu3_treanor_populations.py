import matplotlib.pyplot as plt
import numpy as np

from _bootstrap import ensure_repo_root

ROOT_DIR = ensure_repo_root()

from research.ch4_treanor import (
    CH4_NU3,
    ln_boltzmann_nu3,
    ln_treanor_nu3,
    treanor_minimum_nu3,
)

OUTPUT_PATH = ROOT_DIR / "ch4_nu3_treanor" / "ch4_nu3_treanor_populations.png"

T0 = 600.0
TV_VALUES = [1000, 2000, 3000, 5000]
V_MAX = 10
v = np.arange(0, V_MAX + 1, dtype=float)

colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

OUTPUT_PATH.parent.mkdir(exist_ok=True)

fig, ax = plt.subplots(figsize=(9, 6))

for Tv, color in zip(TV_VALUES, colors):
    ln_tr = ln_treanor_nu3(v, Tv, T0, CH4_NU3)
    ln_bz = ln_boltzmann_nu3(v, Tv, CH4_NU3)
    v_star = treanor_minimum_nu3(Tv, T0, CH4_NU3)

    ax.plot(v, ln_tr, "-", color=color, lw=1.8, label=f"Treanor $T_v$={Tv} K")
    ax.plot(
        v, ln_bz, "--", color=color, lw=1.2, alpha=0.7, label=f"Boltzmann $T$={Tv} K"
    )

    if 1 < v_star < V_MAX:
        v_star_int = int(round(v_star))
        ax.plot(
            v_star_int,
            ln_treanor_nu3(v_star_int, Tv, T0, CH4_NU3),
            "o",
            color=color,
            ms=7,
            zorder=5,
        )

ln_bz_T0 = ln_boltzmann_nu3(v, T0, CH4_NU3)
ax.plot(v, ln_bz_T0, "k:", lw=1.5, alpha=0.5, label=f"Boltzmann $T$={T0:.0f} K (ref)")

ax.set_xlabel(r"Vibrational quantum number $v_3$", fontsize=12)
ax.set_ylabel(r"$\ln(N_v\, /\, N_0)$", fontsize=12)
ax.set_title(f"CH₄ ν₃ Treanor vs Boltzmann — T₀ = {T0:.0f} K", fontsize=13)

ax.set_xlim(0, V_MAX)
ax.legend(fontsize=8, ncol=2, loc="lower left", framealpha=0.9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_PATH, dpi=200)
print(f"Saved: {OUTPUT_PATH}")
for Tv in TV_VALUES:
    v_star = treanor_minimum_nu3(Tv, T0, CH4_NU3)
    print(f"  Tv={Tv:5d} K  v*={v_star:.1f}")
