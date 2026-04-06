import matplotlib.pyplot as plt
import numpy as np

from _bootstrap import ensure_repo_root

ROOT_DIR = ensure_repo_root()
OUTPUT_PATH = ROOT_DIR / "treanor_distribution" / "treanor_fig3_3_ev.png"

omega_e = 2169.81  # cm^-1
omega_exe = 13.29  # cm^-1
hc_k = 1.4388  # hc/k_B in cm·K

theta_v = hc_k * omega_e
x_e = omega_exe / omega_e


def energy_v(v):
    """Vibrational energy E_v in temperature units [K]."""
    return hc_k * (omega_e * v - omega_exe * v * (v + 1))


def ln_treanor(v, Tv, T0):
    """ln f(E_v, T_v, T_0) — Fridman Eq. (3-37), v^2 form."""
    return -v * theta_v / Tv + x_e * theta_v * v**2 / T0


def ln_boltzmann(v, Tv):
    """ln of Boltzmann distribution at temperature Tv."""
    return -energy_v(v) / Tv


Tv = 5000
T0 = 300
v_max = 40
v = np.arange(0, v_max + 1, dtype=float)
Ev = energy_v(v)

ln_tr = ln_treanor(v, Tv, T0)
ln_bz = ln_boltzmann(v, Tv)

v_star = T0 / (2 * x_e * Tv)

fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(Ev, ln_tr, "k-", lw=1.8, label="Treanor Distribution")
ax.plot(Ev, ln_bz, "k--", lw=1.5, label="Boltzmann Distribution")

ax.set_xlabel(r"$E_v$", fontsize=12)
ax.set_ylabel(r"$\ln\, f(E_v,\, T_v,\, T_0)$", fontsize=12)
ax.set_title("Treanor vs Boltzmann — CO", fontsize=13)

v_plot_max = min(v_max, int(2.5 * v_star) + 5)
Ev_cutoff = energy_v(v_plot_max)
ax.set_xlim(0, Ev_cutoff)

y_min = min(ln_tr[: v_plot_max + 1].min(), ln_bz[: v_plot_max + 1].min())
ax.set_ylim(y_min * 1.15, 0.5)

ax.legend(fontsize=8, loc="lower left", framealpha=0.9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_PATH, dpi=150)
print(f"Saved: {OUTPUT_PATH}")
print(f"  Tv={Tv} K, T0={T0} K, v*={v_star:.1f}, v_plot_max={v_plot_max}")
print(f"  E_v range: 0 to {Ev_cutoff:.0f} K")
plt.show()
