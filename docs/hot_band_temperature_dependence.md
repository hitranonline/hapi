# Why Hot Bands Vanish at Low Temperature

## The short answer

Hot bands originate from vibrationally excited lower states. At low temperature, almost no molecules occupy those states, so the absorption is negligible.

## Boltzmann population

The fraction of molecules in a vibrational state with energy $E$ is proportional to:

$$f \propto \exp\!\left(-\frac{E}{k_B T}\right)$$

where $k_B$ is Boltzmann's constant and $T$ is temperature in Kelvin.

For $\text{CH}_4$ $\nu_3$, one quantum of vibration corresponds to $\nu_3 = 3019\ \text{cm}^{-1}$. Converting to a characteristic vibrational temperature using $c_2 = hc/k_B = 1.4388\ \text{cm} \cdot \text{K}$:

$$\theta_{\text{vib}} = c_2 \cdot \nu_3 = 1.4388 \times 3019 \approx 4343\ \text{K}$$

So the population ratio between $v_3 = n$ and $v_3 = 0$ is:

$$\frac{N(v_3 = n)}{N(v_3 = 0)} = \exp\!\left(-\frac{n \cdot 4343}{T}\right)$$

At 1000 K, $kT/hc \approx 695\ \text{cm}^{-1}$, which is only 23% of $\nu_3$. The mode is still largely "frozen out" — most molecules remain in $v_3 = 0$.

| $T$ (K) | $v_3 = 1$ / $v_3 = 0$ | $v_3 = 2$ / $v_3 = 0$ | Meaning |
| --- | --- | --- | --- |
| 296 | $4.2 \times 10^{-7}$ | $1.8 \times 10^{-13}$ | Essentially zero — less than 1 in a million in $v_3 = 1$ |
| 600 | $7.2 \times 10^{-4}$ | $5.1 \times 10^{-7}$ | About 1 in 1400 — small but measurable |
| 1000 | $1.3 \times 10^{-2}$ | $1.7 \times 10^{-4}$ | Only 1.3% in $v_3 = 1$; explains ~100× drop from fundamental |
| 1500 | $5.5 \times 10^{-2}$ | $3.1 \times 10^{-3}$ | About 1 in 18 — hot bands becoming significant |
| 2000 | $1.1 \times 10^{-1}$ | $1.3 \times 10^{-2}$ | 11% in $v_3 = 1$ — 1→2 hot band clearly visible |
| 2500 | $1.8 \times 10^{-1}$ | $3.1 \times 10^{-2}$ | About 1 in 6 — hot bands are strong |

To get the $1 \to 2$ hot band absorbance within an order of magnitude of the fundamental requires $T \gtrsim 2000\ \text{K}$, where 11% of molecules occupy $v_3 = 1$. For the $2 \to 3$ band to reach comparable relative strength would require temperatures beyond the TIPS 2500 K limit.

## Full line intensity temperature scaling

The Boltzmann factor alone does not determine the observed line strength. HAPI scales each line intensity from the reference temperature $T_{\text{ref}} = 296\ \text{K}$ to the target temperature $T$ using:

$$S(T) = S(T_{\text{ref}}) \cdot \frac{Q(T_{\text{ref}})}{Q(T)} \cdot \frac{\exp(-c_2 E'' / T)}{\exp(-c_2 E'' / T_{\text{ref}})} \cdot \frac{1 - \exp(-c_2 \nu_0 / T)}{1 - \exp(-c_2 \nu_0 / T_{\text{ref}})}$$

where:

- $S(T_{\text{ref}})$ is the reference line intensity (from the database, at 296 K)
- $Q(T)$ is the total internal partition function at temperature $T$ (TIPS2025)
- $E''$ is the lower-state energy in $\text{cm}^{-1}$
- $\nu_0$ is the line center in $\text{cm}^{-1}$
- $c_2 = hc/k_B \approx 1.4388\ \text{cm} \cdot \text{K}$

Each factor has a distinct physical role:

| Factor | Physical role | Effect on hot bands |
| --- | --- | --- |
| $Q(T_{\text{ref}}) / Q(T)$ | Partition function ratio | Always $\leq 1$ for $T > T_{\text{ref}}$. Dilutes all transitions at high $T$ because population spreads across more states. |
| $\exp(-c_2 E'' / T) / \exp(-c_2 E'' / T_{\text{ref}})$ | Lower-state population ratio | This is the Boltzmann factor. For hot bands, $E''$ is large ($\geq 3000\ \text{cm}^{-1}$), so this factor grows enormously with $T$. |
| $[1 - \exp(-c_2 \nu_0 / T)] / [1 - \exp(-c_2 \nu_0 / T_{\text{ref}})]$ | Stimulated emission correction | Close to 1 for $\nu_3$ transitions at all relevant temperatures. Minor effect. |

### Why the fundamental weakens at high T

For the fundamental band ($0 \to 1$), $E'' \approx 0$ (ground state), so the Boltzmann ratio is $\approx 1$ at all temperatures. The only temperature-dependent factor that matters is $Q(T_{\text{ref}}) / Q(T)$, which decreases as $T$ rises. This is why the fundamental peak absorbance drops from $1.61$ at 300 K to $2.99 \times 10^{-4}$ at 2500 K.

### Why hot bands peak and then decline

For the $1 \to 2$ hot band, two competing effects determine $S(T)$:

1. The Boltzmann factor $\exp(-c_2 E'' / T)$ **grows** with $T$ (more molecules in $v_3 = 1$)
2. The partition function ratio $Q(T_{\text{ref}}) / Q(T)$ **shrinks** with $T$ (population spreads across all states)

At moderate temperatures (600–1000 K), effect 1 dominates and the hot band strengthens. At very high temperatures ($\geq 2000$ K), effect 2 overtakes and the hot band weakens. This explains why $\nu_3$: $1 \to 2$ peaks around 1000 K ($6.45 \times 10^{-4}$) and then drops at 2500 K ($6.64 \times 10^{-5}$).

For the higher hot bands ($2 \to 3$, $3 \to 4$), the Boltzmann factor has a steeper exponential rise (larger $E''$), so their peak shifts to even higher temperatures. The $3 \to 4$ band is still growing at 2500 K.

## What this means for each progression

- **$\nu_3$: $0 \to 1$ (fundamental)**: Lower state is the ground state ($v_3 = 0$). Fully populated at all temperatures. Weakens at high $T$ due to partition function dilution.
- **$\nu_3$: $1 \to 2$ (first hot band)**: Lower state is $v_3 = 1$. Requires $\exp(-4343/T)$ population. Negligible at 296 K ($4.2 \times 10^{-7}$), visible at 600 K, peaks near 1000 K — but even there, only 1.3% of molecules are in $v_3 = 1$, producing absorbance ~100× weaker than the fundamental.
- **$\nu_3$: $2 \to 3$ (second hot band)**: Lower state is $v_3 = 2$. Requires $\exp(-8686/T)$ population. The square of the first hot band factor. At 1000 K this is $1.7 \times 10^{-4}$, yielding absorbance ~5000× weaker than the fundamental.
- **$\nu_3$: $3 \to 4$ (third hot band)**: Lower state is $v_3 = 3$. Requires $\exp(-13029/T)$. At 296 K this is effectively zero. Even at 1000 K it is $\sim 2.2 \times 10^{-6}$. It only becomes detectable above $\sim 1500$ K.

## Data source

The hot band lines used in these results come from **ExoMol MM I1** — not HITRAN. HITRAN's $\text{CH}_4$ database does not include $\nu_3$ hot band transitions ($1 \to 2$, $2 \to 3$, $3 \to 4$). Only the fundamental ($0 \to 1$) has HITRAN coverage.

The ExoMol MM (methane) line list includes all vibrational bands up to high excitation, which is why the hot band progressions appear exclusively as `exomol_only` in the J-pair CSV (`hitran_line_count = 0` for all $\nu_3 \geq 1$ transitions).

Both sources are rendered through HAPI's Voigt profile engine using the temperature scaling formula above. The ExoMol intensity threshold is set to $0.0$ (no filtering) to ensure weak hot-band reference-temperature intensities ($\sim 10^{-26}\ \text{cm/molecule}$) are not discarded before the temperature scaling can amplify them.

## Observed results from this workflow

Running `plot_combined_exomol_i1_absorbance_progressions.py` at four temperatures confirms the physics. All cases use $P = 3\ \text{Torr}$, $x = 0.008$, $L = 100\ \text{cm}$, $\Delta\nu = 0.001\ \text{cm}^{-1}$.

| Progression | 300 K | 600 K | 1000 K | 2500 K |
| --- | --- | --- | --- | --- |
| $\nu_3$: $0 \to 1$ | $1.61$ | $2.92 \times 10^{-1}$ | $5.93 \times 10^{-2}$ | $2.99 \times 10^{-4}$ |
| $\nu_3$: $1 \to 2$ | $9.97 \times 10^{-7}$ | $1.97 \times 10^{-4}$ | $6.45 \times 10^{-4}$ | $6.64 \times 10^{-5}$ |
| $\nu_3$: $2 \to 3$ | $7.17 \times 10^{-13}$ | $2.39 \times 10^{-7}$ | $1.22 \times 10^{-5}$ | $9.35 \times 10^{-6}$ |
| $\nu_3$: $3 \to 4$ | $0.0$ | $1.02 \times 10^{-10}$ | $8.55 \times 10^{-8}$ | $8.03 \times 10^{-7}$ |

The corresponding figures are in the artifact directories:

- `artifacts/combined_exomol_i1_absorbance_T300K_P3Torr_x0p008_L100cm_step0p001/`
- `artifacts/combined_exomol_i1_absorbance_T600K_P3Torr_x0p008_L100cm_step0p001/`
- `artifacts/combined_exomol_i1_absorbance_T1000K_P3Torr_x0p008_L100cm_step0p001/`
- `artifacts/combined_exomol_i1_absorbance_T2500K_P3Torr_x0p008_L100cm_step0p001/`

Each directory contains one PNG and one HTML per progression (`nu3_0_to_1`, `nu3_1_to_2`, `nu3_2_to_3`, `nu3_3_to_4`) plus a J-pair CSV, a summary CSV, and a report.

## Why the CSV shows `peak_source = "none"`

When HAPI computes the temperature-scaled line intensities and all of them produce absorbance below floating-point significance, the resulting spectrum is identically zero across the entire grid. The peak-finding logic then returns:

- `peak_wavenumber = 2500.0` (the grid start, from `np.argmax` on an all-zero array)
- `peak_absorbance = 0.0`
- `peak_source = "none"`

This is not a bug. It correctly reflects that no absorption feature exists at that temperature.

## Boosting hot band absorbance to match the fundamental

Because absorbance scales linearly with mole fraction ($x$), pressure ($P$), and path length ($L$) — via $A = \sigma \cdot n \cdot L$ where $n = Px/(k_B T)$ — the Boltzmann penalty on hot bands can be compensated by adjusting gas conditions. Temperature also plays a role but has competing effects (more population in excited states vs. partition function dilution and lower number density).

### Which knobs preserve spectral shape

| Parameter | Effect on absorbance | Effect on line shape |
| --- | --- | --- |
| Path length $L$ | Linear scaling | None — cleanest knob |
| Mole fraction $x$ | Linear scaling | None |
| Pressure $P$ | Linear scaling (via $n$) | Broadens Lorentz wings — smears J-pair structure |
| Temperature $T$ | Complex (Boltzmann × partition function × $1/T$) | Changes Doppler width; shifts J-pair intensity envelope |

To preserve the resolved J-pair structure visible in the fundamental band, prefer boosting $x$ and $L$ over $P$.

### Recommended conditions for ideal hot band plots

At the reference conditions ($P = 3\ \text{Torr}$, $x = 0.008$, $L = 100\ \text{cm}$), the $1 \to 2$ hot band is ~100× weaker than the fundamental at 1000 K. The product $x \cdot L$ must increase by that factor to compensate.

The best verified condition uses **$T = 2500\ \text{K}$, $x = 1.0$ (pure $\text{CH}_4$), $P = 3\ \text{Torr}$, $L = 1000\ \text{cm}$**. This combines three effects:

1. High temperature pushes 18% of molecules into $v_3 = 1$ and 3.1% into $v_3 = 2$
2. Pure methane gives 125× boost over $x = 0.008$
3. The 10 m path gives 10× boost over 100 cm

| Progression | Peak absorbance | Y-axis scale | Quality |
| --- | --- | --- | --- |
| $\nu_3$: $0 \to 1$ | $5.71 \times 10^{-1}$ | $10^{-1}$ | Strong, well-resolved J-pairs |
| $\nu_3$: $1 \to 2$ | $8.26 \times 10^{-2}$ | $10^{-2}$ | Ideal — clear PQR branch structure |
| $\nu_3$: $2 \to 3$ | $1.16 \times 10^{-2}$ | $10^{-2}$ | Good — visible J-pair structure |
| $\nu_3$: $3 \to 4$ | $1.00 \times 10^{-3}$ | $10^{-3}$ | Detectable, individual J-pairs resolved |

The figures are in `artifacts/combined_exomol_i1_absorbance_T2500K_P3Torr_x1p0_L1000cm_step0p001/`.

![nu3 0→1 (fundamental)](../artifacts/combined_exomol_i1_absorbance_T2500K_P3Torr_x1p0_L1000cm_step0p001/nu3_0_to_1_absorbance.png)

![nu3 1→2 (first hot band)](../artifacts/combined_exomol_i1_absorbance_T2500K_P3Torr_x1p0_L1000cm_step0p001/nu3_1_to_2_absorbance.png)

![nu3 2→3 (second hot band)](../artifacts/combined_exomol_i1_absorbance_T2500K_P3Torr_x1p0_L1000cm_step0p001/nu3_2_to_3_absorbance.png)

![nu3 3→4 (third hot band)](../artifacts/combined_exomol_i1_absorbance_T2500K_P3Torr_x1p0_L1000cm_step0p001/nu3_3_to_4_absorbance.png)

### Why not just raise temperature?

Temperature alone cannot bring the hot band up to fundamental-band levels. The $1 \to 2$ hot band peak absorbance at reference conditions ($x = 0.008$, $L = 100\ \text{cm}$) peaks near 1000 K at $6.45 \times 10^{-4}$ and then *declines* at higher $T$ as the partition function dilution overtakes the Boltzmann gain. Even at the optimum, the hot band is ~100× below the fundamental. Raising $T$ helps but must be combined with larger $x$ or $L$ to reach the $10^{-2}$ absorbance range.

### Alternative condition sets for the $1 \to 2$ hot band

All target peak absorbance $\sim 0.01$–$0.06$ (the "ideal" fundamental range):

| $T$ (K) | $x$ | $P$ (Torr) | $L$ (cm) | Estimated $A(1 \to 2)$ | Trade-off |
| --- | --- | --- | --- | --- | --- |
| 1500 | 1.0 | 3 | 100 | $\sim 2.3 \times 10^{-2}$ | Both fundamental ($0.39$) and hot band in a useful range |
| 1000 | 0.1 | 3 | 1000 | $\sim 8.1 \times 10^{-2}$ | Moderate mix, 10 m path |
| 2500 | 1.0 | 3 | 1000 | $8.3 \times 10^{-2}$ | Best for all four progressions simultaneously |
| 600 | 1.0 | 3 | 100 | $\sim 2.5 \times 10^{-2}$ | Fundamental saturated ($A \approx 37$) |

## Practical implication

When running multi-temperature comparisons:

- At 296 K, only the fundamental band ($0 \to 1$) produces meaningful absorption
- Hot bands ($1 \to 2$, $2 \to 3$, $3 \to 4$) require elevated temperatures to appear
- Each hot band has a temperature at which it peaks, determined by the competition between Boltzmann population growth and partition function dilution
- **At 1000 K**, the $1 \to 2$ hot band is visible but still ~100× weaker than the fundamental because only 1.3% of molecules populate $v_3 = 1$. The $2 \to 3$ band is ~5000× weaker. This is real physics — the $\nu_3$ mode characteristic temperature ($\theta_{\text{vib}} = 4343\ \text{K}$) is far above 1000 K, so the mode remains largely frozen out
- **At 2000–2500 K**, hot bands become strong enough for quantitative comparison with the fundamental. This is the regime where $kT/hc$ approaches $\nu_3$ and Boltzmann population of excited states becomes significant
- The HAPI TIPS2025 partition function supports temperatures up to 2500 K, which limits how high the hot bands can be pushed
