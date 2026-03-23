# ExoMol Database and Spectrum Modeling Notes

This note combines the ExoMol database overview and the ExoMol spectrum-modeling explanation for this repo.

It focuses on the local CH4 MM workflow implemented in [scripts/plot_exomol_ch4_mm_absorbance.py](../scripts/plot_exomol_ch4_mm_absorbance.py), with related usage in [scripts/scan_exomol_ch4_nu3_band_types.py](../scripts/scan_exomol_ch4_nu3_band_types.py) and the ExoMol-to-HITRAN-style exporter [scripts/build_exomol_ch4_mm_hitran_db.py](../scripts/build_exomol_ch4_mm_hitran_db.py).

Note: this note predates the `scripts/` reorganization, so some inline file links further down may still reflect the old root-level locations.

## Big Picture

ExoMol is not arranged like a HITRAN-ready line list where each row already gives:

- line center `nu`
- reference line intensity `sw`
- pressure-broadening fields

Instead, ExoMol splits the information across separate files:

- `.states`: one row per quantum state
- `.trans`: one row per transition between two states
- `.pf`: partition function table
- `.def`: dataset metadata

The main idea is:

- each state is stored once in `.states`
- transitions refer to states by ID in `.trans`
- transition-specific radiative information is stored in `.trans`
- temperature-dependent line quantities are derived later

So ExoMol is state-plus-transition oriented, not line-list oriented in the same way as HITRAN.

## Why ExoMol Uses Separate `states` and `trans`

The same state can appear in many transitions.

If every transition row repeated the full upper-state and lower-state information, the same energies and quantum labels would be duplicated many times. For a very large line list, that would waste a lot of space.

So ExoMol uses a normalized structure:

- `.states` stores each state once
- `.trans` stores links between states plus the Einstein-A value

This means each state energy only needs to be declared once.

## What Each ExoMol File Means

### `.states`

The `.states` file is the state table.

Each row defines a state such as:

- `state_id`
- energy
- degeneracy
- quantum labels
- symmetry labels

In this repo, state energies and degeneracies are loaded by `load_state_arrays(...)` in [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py#L248) and more detailed state-label parsing is used in [extract_exomol_ch4_mm_pure_nu3_band_texts.py](/f:/GitHub/hapi/extract_exomol_ch4_mm_pure_nu3_band_texts.py).

### `.trans`

The `.trans` file is the transition table.

For the workflow used here, each row is treated as:

- `upper_id`
- `lower_id`
- `A`

where `A` is the Einstein-A coefficient.

The `.trans` row does not directly store:

- line-center wavenumber
- temperature-dependent line intensity

Those are derived later by combining `.trans` with `.states` and `.pf`.

In this repo, transition rows are streamed in:

- `collect_relevant_transitions(...)` in [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py#L350)
- the main scan loop in [scan_exomol_ch4_nu3_band_types.py](/f:/GitHub/hapi/scan_exomol_ch4_nu3_band_types.py)
- the main conversion loop in [build_exomol_ch4_mm_hitran_db.py](/f:/GitHub/hapi/build_exomol_ch4_mm_hitran_db.py)

### `.pf`

The `.pf` file contains the partition function `Q(T)` as a function of temperature.

In this repo it is read by:

- `load_partition_function(...)`
- `interpolate_partition_function(...)`

in both [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py#L232) and [build_exomol_ch4_mm_hitran_db.py](/f:/GitHub/hapi/build_exomol_ch4_mm_hitran_db.py).

### `.def`

The `.def` file contains dataset metadata such as:

- isotopologue mass
- number of states
- default Lorentz broadening width `gamma0`
- temperature exponent `n_exponent`
- quantum-label schema for the `.states` file

In this repo it is read by:

- `parse_def_file(...)` in [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py#L206)
- `parse_exomol_def(...)` in [extract_exomol_ch4_mm_pure_nu3_band_texts.py](/f:/GitHub/hapi/extract_exomol_ch4_mm_pure_nu3_band_texts.py#L172)
- `parse_exomol_def(...)` in [build_exomol_ch4_mm_hitran_db.py](/f:/GitHub/hapi/build_exomol_ch4_mm_hitran_db.py)

## How ExoMol Defines a Spectral Line

For a transition, ExoMol gives you:

- which upper state
- which lower state
- the Einstein-A value

From the two state IDs, you can look up the state energies and compute the line center.

### Line Center

The line-center wavenumber is:

```text
nu0 = E_upper - E_lower
```

If the state energies are in `cm^-1`, then `nu0` is also in `cm^-1`.

In this repo, that happens in:

- `collect_relevant_transitions(...)` in [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py#L350)
- the scan loop in [scan_exomol_ch4_nu3_band_types.py](/f:/GitHub/hapi/scan_exomol_ch4_nu3_band_types.py)
- the output-row builder in [build_exomol_ch4_mm_hitran_db.py](/f:/GitHub/hapi/build_exomol_ch4_mm_hitran_db.py)

This is why ExoMol does not need to store the line center directly in `.trans`: it is already implied by the two state energies.

### Why `.trans` Stores `A` Instead of the Line Center

The line center can be reconstructed from `.states`.

The Einstein-A value cannot.

So ExoMol stores the transition-specific quantity that is not redundant:

- `A`, the spontaneous-emission probability for that exact upper-to-lower transition

Short version:

- `Delta E` tells you where the line is
- `A` tells you how strongly that transition radiates

They belong to the same transition, but they are not the same physical quantity.

## How Line Intensity Is Obtained

ExoMol does not usually give a ready-to-use HITRAN-style line intensity in each transition row.

Instead, the line intensity is derived from:

- Einstein-A
- upper-state degeneracy
- lower-state energy
- line-center wavenumber
- gas temperature
- partition function `Q(T)`

In this repo, that conversion is done by:

- `lte_line_intensity_cm_per_molecule(...)` in [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py#L271)
- `lte_line_intensity_cm_per_molecule(...)` in [build_exomol_ch4_mm_hitran_db.py](/f:/GitHub/hapi/build_exomol_ch4_mm_hitran_db.py)

The output of that function is:

- integrated line intensity `S(T)`
- units: `cm/molecule`

This is not absorbance yet. It is the per-line strength before broadening is applied.

## Spectrum Modeling Pipeline in This Repo

The CH4 ExoMol MM modeling workflow in [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py) is:

1. Read `.def`, `.states`, `.pf`, and the overlapping `.trans` files.
2. Reconstruct each line center from state energies.
3. Convert Einstein-A into LTE line intensity `S(T)`.
4. Keep only transitions in the requested spectral window and above the intensity threshold.
5. Apply a Voigt profile to each kept line.
6. Sum all broadened lines into a cross section.
7. Convert cross section to absorbance for the chosen gas condition.

### Step 1: Stream the compressed ExoMol files

In this repo, `.states.bz2` and `.trans.bz2` are not converted into permanent `.txt` files on disk.

They are decompressed in memory and parsed line by line by:

- `iter_bz2_text_lines(...)` in [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py#L133)
- equivalent streaming logic in [build_exomol_ch4_mm_hitran_db.py](/f:/GitHub/hapi/build_exomol_ch4_mm_hitran_db.py)

### Step 2: Keep only relevant transitions

The script does not use every transition in the full database for a windowed spectrum.

`collect_relevant_transitions(...)` in [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py#L350):

- scans only the `.trans.bz2` chunks overlapping the requested spectral window
- reconstructs the wavenumber of each transition
- keeps only transitions in the requested range plus a line-wing margin
- discards lines weaker than the chosen intensity threshold

The kept outputs from that step are arrays of:

- line centers
- line intensities

### Step 3: Convert line intensity into a broadened line shape

Once a line has:

- center `nu0`
- line intensity `S(T)`
- broadening parameters

the script applies a profile.

In this repo:

- `doppler_hwhm_cm(...)` computes Doppler broadening in [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py#L327)
- `voigt_profile_cm(...)` defines the Voigt line shape in [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py#L333)
- `render_cross_section(...)` applies the profile and accumulates the result on the grid in [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py#L414)

Conceptually, one line contributes:

```text
sigma_i(nu) = S_i(T) * profile_i(nu - nu0_i)
```

where `profile_i` is normalized and centered at that line's `nu0_i`.

### Step 4: Sum all lines into cross section

The total cross section is the sum of all broadened lines:

```text
sigma(nu) = sum_i sigma_i(nu)
```

In the code, `render_cross_section(...)` loops over kept transitions and adds each broadened contribution onto the spectral grid.

The spectral grid is built in `main()` with:

```text
grid = arange(wn_min, wn_max, wn_step)
```

See [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py#L551).

### Step 5: Convert cross section to absorbance

After the cross section is built, the script converts it to absorbance using Beer-Lambert scaling:

```text
absorbance = cross_section * absorber_number_density * path_length
```

In this repo:

- total number density comes from `number_density_cm3(...)` in [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py#L344)
- absorber number density is total density times CH4 mole fraction
- absorbance is computed in `main()` in [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py#L564)

## Where the Case Settings Are Set

In [plot_exomol_ch4_mm_absorbance.py](/f:/GitHub/hapi/plot_exomol_ch4_mm_absorbance.py), the main case settings are command-line arguments in `parse_args(...)`.

Important ones are:

- `--wn-min`
- `--wn-max`
- `--wn-step`
- `--temperature-k`
- `--pressure-torr`
- `--mole-fraction`
- `--path-length-cm`
- `--intensity-threshold`
- `--line-cutoff`
- optional overrides for `--gamma0` and `--n-exponent`

These settings control:

- which transitions are scanned
- how strong the lines are at the chosen temperature
- how wide the lines become after broadening
- how the final cross section is scaled into absorbance

If concentration is known in ppm, convert it to mole fraction before passing it in:

```text
1000 ppm = 1000e-6 = 0.001
```

## How This Compares With HITRAN

For HITRAN, each row is already much closer to a spectroscopy-ready line record.

A HITRAN row directly gives:

- line-center wavenumber `nu`
- reference line intensity `sw` at 296 K
- lower-state energy
- broadening parameters

So HITRAN is more transition-centered and spectrum-oriented.

ExoMol is more state-plus-transition oriented:

- `.states` defines the states
- `.trans` links them and stores `A`
- the line center and line intensity are reconstructed later

With ExoMol, the repo has to derive important quantities first:

- line center from state energies
- line intensity from Einstein `A` and `Q(T)`
- pressure-broadened line shape from metadata plus case conditions

So the ExoMol workflow is more of a reconstruction-and-modeling pipeline than a direct read-and-plot workflow.

## What Is the Same Between HITRAN and ExoMol

Once you already have:

- line center
- line intensity
- broadening parameters

the downstream spectral physics is basically the same:

1. choose a line profile
2. evaluate it on a spectral grid
3. sum the contributions
4. convert to cross section, transmittance, or absorbance

So the main difference between HITRAN and ExoMol is earlier in the pipeline:

- HITRAN starts closer to a ready-made line list
- ExoMol starts from states, transitions, and Einstein-A values

## Practical Summary

The shortest possible mental model:

- HITRAN: each row already looks like a spectral line
- ExoMol: reconstruct each spectral line from state IDs plus Einstein-A

For ExoMol in this repo:

- `.states` gives the energies and degeneracies
- `.trans` gives `upper_id`, `lower_id`, and `A`
- line center comes from `E_upper - E_lower`
- line intensity comes from `A` plus temperature-dependent physics
- Voigt broadening turns each line into a smooth function on the grid
- summed cross section is converted to absorbance using number density and path length

For the local CH4 MM modeling script, the short workflow answer is:

ExoMol gives state energies, degeneracies, Einstein `A`, and `Q(T)`. The repo reconstructs each line center, computes temperature-dependent line intensity, broadens each line with a Voigt profile, sums the contributions into cross section, and then converts cross section to absorbance.
