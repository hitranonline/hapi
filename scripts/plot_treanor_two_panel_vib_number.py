import matplotlib.pyplot as plt
import numpy as np

from _bootstrap import ensure_repo_root

ROOT_DIR = ensure_repo_root()

from research.treanor import CO, boltzmann, treanor, treanor_minimum

OUTPUT_PATH = ROOT_DIR / "treanor_distribution" / "treanor_two_panel_vib_number.png"

v = np.arange(0, 45, dtype=float)
T0 = 300
Tv = 5000

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

axes[0].semilogy(
    v,
    treanor(v, Tv, T0, CO),
    "b-o",
    ms=4,
    label=f"Treanor ($T_v$={Tv} K, $T_0$={T0} K)",
)
axes[0].semilogy(v, boltzmann(v, Tv, CO), "r--s", ms=4, label=f"Boltzmann ($T$={Tv} K)")
axes[0].semilogy(
    v, boltzmann(v, T0, CO), "k:^", ms=4, alpha=0.5, label=f"Boltzmann ($T$={T0} K)"
)

v_star = treanor_minimum(Tv, T0, CO)
axes[0].axvline(
    v_star,
    color="green",
    ls="--",
    alpha=0.7,
    label=f"$v^*_{{\\mathrm{{Tr}}}}$ = {v_star:.1f}",
)

axes[0].set_xlabel("Vibrational quantum number $v$")
axes[0].set_ylabel("$N_v / N_0$")
axes[0].set_title("Treanor vs Boltzmann (CO)")
axes[0].legend(fontsize=9)
axes[0].set_ylim(1e-30, 10)
axes[0].grid(True, alpha=0.3)

for Tv_val in [1000, 2000, 3000, 5000, 8000]:
    axes[1].semilogy(
        v, treanor(v, Tv_val, T0, CO), "-o", ms=3, label=f"$T_v$ = {Tv_val} K"
    )

axes[1].set_xlabel("Vibrational quantum number $v$")
axes[1].set_ylabel("$N_v / N_0$")
axes[1].set_title(f"Treanor Distribution ($T_0$ = {T0} K, CO)")
axes[1].legend(fontsize=9)
axes[1].set_ylim(1e-30, 10)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_PATH, dpi=150)
print(f"Saved: {OUTPUT_PATH}")
plt.show()
