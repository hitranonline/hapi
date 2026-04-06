import matplotlib.pyplot as plt
import numpy as np

from _bootstrap import ensure_repo_root

ROOT_DIR = ensure_repo_root()

from research.treanor import CO, ln_boltzmann, ln_treanor, treanor_minimum

OUTPUT_PATH = ROOT_DIR / "treanor_distribution" / "treanor_fig3_3_vib_number.png"

T0 = 300  # translational temperature [K]
Tv_values = [1000, 2000, 3000, 5000, 8000]  # vibrational temperatures [K]
v_max = 40
v = np.arange(0, v_max + 1, dtype=float)

colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

fig, ax = plt.subplots(figsize=(9, 6))

for Tv, color in zip(Tv_values, colors):
    ln_tr = ln_treanor(v, Tv, T0, CO)
    ln_bz = ln_boltzmann(v, Tv, CO)
    v_star = treanor_minimum(Tv, T0, CO)

    ax.plot(v, ln_tr, "-", color=color, lw=1.8, label=f"Treanor $T_v$={Tv} K")
    ax.plot(
        v, ln_bz, "--", color=color, lw=1.2, alpha=0.7, label=f"Boltzmann $T$={Tv} K"
    )

    if 1 < v_star < v_max:
        v_star_int = int(round(v_star))
        ax.plot(
            v_star_int,
            ln_treanor(v_star_int, Tv, T0, CO),
            "o",
            color=color,
            ms=6,
            zorder=5,
        )

ln_bz_T0 = ln_boltzmann(v, T0, CO)
ax.plot(v, ln_bz_T0, "k:", lw=1.5, alpha=0.5, label=f"Boltzmann $T$={T0} K (ref)")

ax.set_xlabel("Vibrational quantum number $v$", fontsize=12)
ax.set_ylabel(r"$\ln(N_v\, /\, N_0)$", fontsize=12)
ax.set_title(f"Treanor vs Boltzmann — CO, $T_0$ = {T0} K", fontsize=13)

ax.set_xlim(0, v_max)
ax.set_ylim(-60, 30)
ax.legend(fontsize=8, ncol=2, loc="lower left", framealpha=0.9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_PATH, dpi=200)
print(f"Saved: {OUTPUT_PATH}")
for Tv in Tv_values:
    v_star = treanor_minimum(Tv, T0, CO)
    print(f"  Tv={Tv:5d} K  v*={v_star:.1f}")
