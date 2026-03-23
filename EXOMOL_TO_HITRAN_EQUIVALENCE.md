# ExoMol to HITRAN: Core-Logic Equivalence

This note explains why the original script at
[xnx/ExoMol_to_HITRAN.py](https://raw.githubusercontent.com/xnx/ExoMol_to_HITRAN/master/ExoMol_to_HITRAN.py)
and the newer CH4 MM converter in this repo perform the same core
ExoMol-to-line-list calculation, even though they do not produce the same
output format or metadata.

## Short Answer

They are equivalent in the core spectroscopy pipeline:

1. Read ExoMol `.states`, `.trans`, and `.pf`.
2. Compute transition wavenumber as upper-state energy minus lower-state energy.
3. Compute LTE line intensity from Einstein-A, upper-state degeneracy, lower-state energy, and partition function.
4. Filter by wavenumber and intensity.
5. Write surviving transitions to an output database.

The newer script adds more metadata, more validation, and a different output
format, but it does not change the underlying line-position or line-intensity
physics.

## One-to-One Mapping

### 1. Input data

Original script:

- Reads `.states` in `read_states`
- Reads `.trans` / `.trans.bz2` in `read_trans`
- Reads `.pf` in `read_Q` / `get_Q`

New script:

- Reads `.states` / `.states.bz2` in `load_state_arrays`
- Reads `.trans` / `.trans.bz2` via `discover_transition_files` and `iter_text_lines`
- Reads `.pf` in `load_partition_function`

Conclusion: both scripts use the same ExoMol source files.

### 2. Transition wavenumber

Original script:

```python
Ep, gp = states[stateIDp-1]
Epp, gpp = states[stateIDpp-1]
nu = Ep - Epp
```

New script:

```python
wavenumber = float(energies[upper_id] - energies[lower_id])
```

These are the same calculation:

\[
\nu = E_{\mathrm{upper}} - E_{\mathrm{lower}}
\]

The indexing differs only because the original script stores states in a
0-based NumPy table, while the newer script stores arrays sized as
`nstates + 1` and uses ExoMol state IDs directly.

### 3. Partition function at temperature T

Original script:

- Reads tabulated `Q(T)` from the `.pf` file.
- Uses linear interpolation in `get_Q`.

New script:

- Reads tabulated `Q(T)` from the `.pf` file.
- Uses `np.interp` in `interpolate_partition_function`.

These are the same mathematical operation.

### 4. LTE line intensity

Original script:

```python
def calc_S(nu0, Epp, gp, A, Q, T):
    c2oT = c2 / T
    S = A * gp * np.exp(-c2oT * Epp) * (1 - np.exp(-c2oT * nu0))\
             / pic8 / nu0**2 / Q
    return S
```

New script:

```python
def lte_line_intensity_cm_per_molecule(
    wavenumber,
    a_coefficient,
    g_upper,
    lower_energy_cm,
    temperature_k,
    partition_function,
):
    boltzmann = math.exp(-SECOND_RADIATION_CONSTANT_CM_K * lower_energy_cm / temperature_k)
    stimulated = 1.0 - math.exp(-SECOND_RADIATION_CONSTANT_CM_K * wavenumber / temperature_k)
    return (
        g_upper
        * a_coefficient
        * boltzmann
        * stimulated
        / (8.0 * math.pi * LIGHT_SPEED_CM_S * wavenumber * wavenumber * partition_function)
    )
```

Term-by-term mapping:

- `nu0` -> `wavenumber`
- `A` -> `a_coefficient`
- `gp` -> `g_upper`
- `Epp` -> `lower_energy_cm`
- `Q` -> `partition_function`
- `T` -> `temperature_k`

Both implement:

$$
S(\nu, T) =
\frac{
g_{\mathrm{upper}} \, A \,
\exp\!\left(-\frac{c_2 E_{\mathrm{lower}}}{T}\right)
\left(1 - \exp\!\left(-\frac{c_2 \nu}{T}\right)\right)
}{
8 \pi c \, \nu^2 \, Q(T)
}
$$

The constant handling is also equivalent:

- Original script computes `c2 = h * c * 100 / kB`
- New script hardcodes `SECOND_RADIATION_CONSTANT_CM_K = 1.438776877`
- Original script uses `c * 100` to convert m/s to cm/s
- New script hardcodes `LIGHT_SPEED_CM_S = 2.99792458e10`

So the newer function is not merely similar in spirit; it is the same LTE
intensity formula written in a cleaner form.

### 5. Wavenumber and intensity filtering

Original script:

```python
if nu < numin:
    continue
if nu > numax:
    ...
S = calc_S(...)
S *= abundance
if S < Smin:
    continue
```

New script:

```python
if args.wn_min is not None and wavenumber < args.wn_min:
    continue
if args.wn_max is not None and wavenumber > args.wn_max:
    ...
sw = lte_line_intensity_cm_per_molecule(...)
sw *= args.abundance
if sw < args.intensity_threshold:
    continue
```

This is the same filtering pipeline:

1. keep only lines in the requested wavenumber window
2. compute LTE intensity
3. scale by abundance
4. discard lines below the intensity threshold

## What Is Not the Same

The scripts are not identical overall.

The newer script adds:

- `.def` parsing
- CH4-specific quantum-label extraction
- symmetry mapping to HITRAN-style labels
- HAPI `.data` and `.header` output instead of a single `.par`
- default broadening parameters such as `gamma_air`, `gamma_self`, `n_air`
- stricter argument validation
- optional transition-file selection and progress reporting

These are output-schema and metadata improvements. They do not alter the core
calculation of line positions or LTE intensities.

## Bottom Line

If the question is:

`Do both scripts perform the same core ExoMol-to-HITRAN line calculation?`

the answer is:

`Yes.`

If the question is:

`Are the two scripts otherwise the same program?`

the answer is:

`No.`

The newer script is a more specialized and more complete implementation built
around the same physics pipeline.
