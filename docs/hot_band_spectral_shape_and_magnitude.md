# Hot Band Spectral Shape and Magnitude

## Why the fundamental is comb-shaped but hot bands are not

### The comb shape of $\nu_3$: $0 \to 1$

The fundamental band has a textbook rovibrational structure with three distinct branches:

| Branch | $\Delta J$ | Spectral position | Shape |
| --- | --- | --- | --- |
| **P branch** | $-1$ | Below band center, marching downward | Evenly spaced peaks |
| **Q branch** | $0$ | Piled up near band center (~3019 cm$^{-1}$) | Sharp central spike |
| **R branch** | $+1$ | Above band center, marching upward | Evenly spaced peaks — the "comb teeth" |

The R branch is the most prominent comb. Its peaks are spaced ~9.4 cm$^{-1}$ apart with a spread of only 1.06 cm$^{-1}$ — nearly perfectly regular. The Q branch clusters within ~11 cm$^{-1}$, forming a single sharp spike at the band center.

This regularity arises because:

1. **The upper state ($v_3 = 1$, mode 0010) has simple symmetry splitting.** Under CH$_4$'s tetrahedral ($T_d$) point group, $v_3 = 1$ transforms as $F_2$ — a single irreducible representation. Each $J$ manifold splits into a modest number of sub-levels.
2. **HITRAN provides experimentally measured line positions.** The strongest R-branch peaks (J $\geq$ 3) are HITRAN-sourced, with positions accurate to ~0.001 cm$^{-1}$. Lines within each J-pair stack coherently at precise positions.
3. **Relatively few lines per J-pair.** The R branch averages 87 lines per J-pair, so the absorption concentrates into narrow spectral features.

### Why the $1 \to 2$ hot band loses the comb shape

The $1 \to 2$ hot band (0010 $\to$ 0020) looks qualitatively different — filled-in rather than comb-shaped. This is physically correct, not an artifact.

**Three factors compound:**

#### 1. The $v_3 = 2$ upper state is more degenerate

The overtone state 0020 splits into multiple vibrational sub-levels under $T_d$ symmetry:

$$v_3 = 1\ (F_2):\ \text{one sub-level}$$
$$v_3 = 2\ (A_1 + E + F_2):\ \text{three sub-levels}$$

Each sub-level has slightly different rotational constants, so transitions from a given lower $J$ to different sub-levels land at different wavenumbers. This multiplies the number of spectral features per J-pair.

#### 2. ExoMol resolves all sub-level transitions

HITRAN curates dominant transitions; ExoMol's ab-initio calculation includes every symmetry-allowed transition. The result:

| | $0 \to 1$ (fundamental) | $1 \to 2$ (hot band) |
| --- | --- | --- |
| Total lines | 42,699 | 282,540 |
| Lines per J-pair (avg) | 378 | 2,974 |
| R-branch lines per J-pair (avg) | 87 | 135 |
| P-branch lines per J-pair (avg) | 448 | 3,728 |
| Q-branch lines per J-pair (avg) | 473 | 4,034 |
| Source | HITRAN + ExoMol | **ExoMol only** |

The $1 \to 2$ band has **6.6× more lines** overall. These additional lines fill the gaps between what would otherwise be comb teeth.

#### 3. Ab-initio line positions scatter

ExoMol positions are computed from a potential energy surface. While globally accurate, individual line positions can scatter by ~0.01–0.1 cm$^{-1}$ relative to true positions. With thousands of lines per J-pair, these small positional uncertainties spread the absorption over a wider wavenumber range instead of stacking into sharp peaks.

### Branch-by-branch comparison at $T = 2500$ K

**R branch ($\Delta J = +1$):**

| | $0 \to 1$ | $1 \to 2$ |
| --- | --- | --- |
| Peak spacing | 9.25–9.75 cm$^{-1}$ (spread 1.1) | 6.96–13.84 cm$^{-1}$ (spread 6.9) |
| Regularity | Extremely regular | **6× more irregular** |
| Strongest source | HITRAN (experimental) | ExoMol (ab-initio) |

The $0 \to 1$ R-branch marches outward in near-perfect ~9.4 cm$^{-1}$ steps. The $1 \to 2$ R-branch has similar average spacing (~10.4 cm$^{-1}$) but 6× more scatter, so the teeth are uneven.

**Q branch ($\Delta J = 0$):**

| | $0 \to 1$ | $1 \to 2$ |
| --- | --- | --- |
| Peak position range | 3007–3019 cm$^{-1}$ (11 cm$^{-1}$) | **2975–3024 cm$^{-1}$ (49 cm$^{-1}$)** |
| Spacing spread | 1.9 cm$^{-1}$ | **80.8 cm$^{-1}$** |
| Structure | Single sharp spike | **Splits into two sub-branches** |

The $0 \to 1$ Q-branch lines pile up near 3019 cm$^{-1}$ into a sharp central spike. The $1 \to 2$ Q-branch splits into two groups — odd-$J$ peaks near ~2980 cm$^{-1}$ and even-$J$ peaks near ~3020 cm$^{-1}$ — scattered across a ~40 cm$^{-1}$ range. This is the clearest fingerprint of the $v_3 = 2$ degeneracy splitting.

**P branch ($\Delta J = -1$):**

| | $0 \to 1$ | $1 \to 2$ |
| --- | --- | --- |
| Peak spacing | $-9$ to $-11.5$ cm$^{-1}$ (spread 2.5) | **Chaotic** (spread 107) |
| Pattern | Smooth downward march | Peak positions jump by tens of cm$^{-1}$ |

The $0 \to 1$ P-branch marches smoothly downward from the band center. The $1 \to 2$ P-branch jumps erratically — e.g. J 3$\to$2 peaks at 2991, J 4$\to$3 jumps to 2911, J 5$\to$4 at 2932.

### The non-comb shape is physically correct

The filled-in appearance of the $1 \to 2$ hot band is real physics:

- The $v_3 = 2$ overtone state has more vibrational sub-levels than $v_3 = 1$
- ExoMol resolves all symmetry-allowed transitions between these sub-levels
- The resulting dense line forest is what the real spectrum would look like in a laboratory measurement at high temperature

The comb shape of the $0 \to 1$ fundamental is the special case — it appears clean because the $v_3 = 1$ upper state has minimal splitting and HITRAN's curated experimental lines reinforce the regular structure.

## Magnitude differences between progressions

### The Boltzmann constraint

Absorbance is proportional to the population in the lower vibrational state. At any given temperature, higher progressions are exponentially weaker:

$$A(v \to v+1) \propto \exp\!\left(-\frac{v \cdot 4343}{T}\right) \cdot x \cdot P \cdot L$$

where $x$ is mole fraction, $P$ is pressure, and $L$ is path length (via number density $n = Px / k_B T$).

At $T = 2500$ K with $x = 1.0$, $P = 3$ Torr, $L = 1000$ cm:

| Progression | Peak absorbance | Ratio to fundamental | Boltzmann factor |
| --- | --- | --- | --- |
| $\nu_3$: $0 \to 1$ | $5.71 \times 10^{-1}$ | 1 | 1.0 |
| $\nu_3$: $1 \to 2$ | $8.26 \times 10^{-2}$ | 0.14 | 0.18 |
| $\nu_3$: $2 \to 3$ | $1.16 \times 10^{-2}$ | 0.020 | 0.031 |
| $\nu_3$: $3 \to 4$ | $1.00 \times 10^{-3}$ | 0.0018 | 0.0054 |

The peak absorbance ratios are smaller than the Boltzmann factors because the hot band intensity is spread across more lines over a wider wavenumber range (the sub-level splitting discussed above). The total integrated absorbance would track the Boltzmann ratio more closely, but the peak is diluted by the denser line forest.

### No single condition is ideal for all progressions

The four progressions span ~4 orders of magnitude in peak absorbance. An "ideal" absorbance for well-resolved plots is roughly $0.01$–$0.1$. Since absorbance scales linearly with $x \cdot P \cdot L$, any multiplicative boost that brings the $3 \to 4$ band into range simultaneously pushes the $0 \to 1$ band into saturation ($A \gg 1$).

| Target | Conditions needed | $A(0 \to 1)$ | Problem |
| --- | --- | --- | --- |
| $0 \to 1$ ideal | $x = 0.008$, $L = 100$ cm, $T = 1000$ K | $5.93 \times 10^{-2}$ | Hot bands invisible |
| $1 \to 2$ ideal | $x = 1.0$, $L = 100$ cm, $T = 1500$ K | $\sim 2.3 \times 10^{-2}$ | Fundamental at 0.39 (OK) |
| $2 \to 3$ ideal | $x = 1.0$, $L = 1000$ cm, $T = 2500$ K | $\sim 5.7 \times 10^{-1}$ | Fundamental approaching saturation |
| $3 \to 4$ ideal | $x = 1.0$, $L = 10000$ cm, $T = 2500$ K | $\sim 5.7$ | **Fundamental saturated** |

In practice, experimentalists measure different progressions in separate runs with conditions optimized for each band of interest.

### Realistic mole fractions

The mole fraction $x = 0.008$ (0.8% CH$_4$) used in the reference runs is physically reasonable:

| Scenario | Typical CH$_4$ mole fraction |
| --- | --- |
| Earth's atmosphere | $\sim 2 \times 10^{-6}$ (1.9 ppm) |
| Combustion exhaust | 0.001–0.01 (0.1–1%) |
| Laboratory gas cell (diluted in N$_2$) | 0.001–0.05 |
| **Reference runs ($x = 0.008$)** | **0.8% — realistic combustion/lab value** |
| Pure CH$_4$ ($x = 1.0$) | Pure gas cell with no diluent |

At $x = 0.008$, the fundamental band absorbance falls in the ideal measurement range ($A \sim 0.01$–$1$) across a wide temperature range. The weak hot band signal at this mole fraction is a real experimental challenge — it is why hot band spectroscopy of polyatomic molecules typically requires either long path lengths (multi-pass cells), high concentrations, or high-sensitivity detection techniques.

### Absorbance and optical thickness

When absorbance exceeds $\sim 1$, the sample becomes optically thick — most photons at that wavenumber are already absorbed, so additional absorbers contribute diminishing returns. The Beer-Lambert law ($A = \sigma n L$) remains mathematically valid, but in practice:

- Detector noise dominates the transmitted signal
- Small baseline errors become amplified
- Spectral features appear to "flatten" at the top

The $0 \to 1$ band at $x = 1.0$, $L = 1000$ cm, $T = 2500$ K has $A \approx 0.57$ — still in the linear regime but approaching the practical limit. This is acceptable for visualizing the hot bands, but not a condition one would use for quantitative fundamental-band analysis.

## Observed results

All figures are from `plot_combined_exomol_i1_absorbance_progressions.py` using ExoMol MM I1 and HITRAN through the HAPI temp-table rendering path. Wavenumber window: 2500–3500 cm$^{-1}$, step 0.001 cm$^{-1}$.

### Reference conditions: $T = 1000$ K, $x = 0.008$, $P = 3$ Torr, $L = 100$ cm

| Progression | Peak absorbance | Shape |
| --- | --- | --- |
| $\nu_3$: $0 \to 1$ | $5.93 \times 10^{-2}$ | Clean comb (HITRAN R-branch dominates) |
| $\nu_3$: $1 \to 2$ | $6.45 \times 10^{-4}$ | Weak, filled-in (ExoMol only) |
| $\nu_3$: $2 \to 3$ | $1.22 \times 10^{-5}$ | Barely detectable |
| $\nu_3$: $3 \to 4$ | $8.55 \times 10^{-8}$ | Noise-level |

Figures: `artifacts/combined_exomol_i1_absorbance_T1000K_P3Torr_x0p008_L100cm_step0p001/`

### Boosted conditions: $T = 2500$ K, $x = 1.0$, $P = 3$ Torr, $L = 1000$ cm

| Progression | Peak absorbance | Shape |
| --- | --- | --- |
| $\nu_3$: $0 \to 1$ | $5.71 \times 10^{-1}$ | Strong comb, approaching optical thickness |
| $\nu_3$: $1 \to 2$ | $8.26 \times 10^{-2}$ | Ideal range, filled-in PQR structure |
| $\nu_3$: $2 \to 3$ | $1.16 \times 10^{-2}$ | Good, visible J-pair features |
| $\nu_3$: $3 \to 4$ | $1.00 \times 10^{-3}$ | Detectable, individual J-pairs resolved |

Figures: `artifacts/combined_exomol_i1_absorbance_T2500K_P3Torr_x1p0_L1000cm_step0p001/`

![nu3 0→1 at 2500 K, x=1.0, L=1000 cm](../artifacts/combined_exomol_i1_absorbance_T2500K_P3Torr_x1p0_L1000cm_step0p001/nu3_0_to_1_absorbance.png)

![nu3 1→2 at 2500 K, x=1.0, L=1000 cm](../artifacts/combined_exomol_i1_absorbance_T2500K_P3Torr_x1p0_L1000cm_step0p001/nu3_1_to_2_absorbance.png)

![nu3 2→3 at 2500 K, x=1.0, L=1000 cm](../artifacts/combined_exomol_i1_absorbance_T2500K_P3Torr_x1p0_L1000cm_step0p001/nu3_2_to_3_absorbance.png)

![nu3 3→4 at 2500 K, x=1.0, L=1000 cm](../artifacts/combined_exomol_i1_absorbance_T2500K_P3Torr_x1p0_L1000cm_step0p001/nu3_3_to_4_absorbance.png)
