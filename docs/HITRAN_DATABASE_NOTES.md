# HITRAN Database Notes

This note explains how the HITRAN database is organized and how it is used in this repo, especially in [scripts/plot_vibrational_mode_progressions.py](../scripts/plot_vibrational_mode_progressions.py), [scripts/extract_ch4_nu3_band_texts.py](../scripts/extract_ch4_nu3_band_texts.py), and the local HAPI implementation in [hapi/hapi.py](../hapi/hapi.py).

Note: this note predates the `scripts/` reorganization, so some inline file links further down may still reflect the old root-level locations.

## Big Picture

HITRAN is arranged much closer to a spectroscopy-ready line list than ExoMol.

Each HITRAN row already contains transition-level quantities such as:

- line-center wavenumber `nu`
- reference line intensity `sw`
- lower-state energy
- broadening parameters
- quantum labels

So a HITRAN row already looks like one spectral line.

That is the main structural difference from ExoMol:

- HITRAN stores line-oriented spectroscopic parameters directly
- ExoMol stores states and transitions separately and reconstructs line quantities later

## What a HITRAN Row Contains

In this repo, the table layout is described by files such as:

- [CH4_M6_I1.header](/f:/GitHub/hapi/hitran_db/CH4_M6_I1.header)
- [CO_M5_I1.header](/f:/GitHub/hapi/hitran_db/CO_M5_I1.header)

Important fields include:

- `nu`: transition wavenumber, the line center
- `sw`: line intensity at the HITRAN reference temperature
- `elower`: lower-state energy
- `gamma_air`: air-broadened Lorentz half-width
- `gamma_self`: self-broadened Lorentz half-width
- `n_air`: temperature exponent for air broadening
- `delta_air`: pressure shift
- `global_upper_quanta`, `global_lower_quanta`: vibrational/electronic labels
- `local_upper_quanta`, `local_lower_quanta`: rotational or local quantum labels

So for a HITRAN row:

- `x` is directly available as `nu`
- the line-strength input is directly available as `sw`

But the final plotted spectrum is still not stored in the row. You still need line-shape modeling and summation over all rows.

## Why HITRAN Feels More Direct Than ExoMol

HITRAN is already transition-centered.

That means a row already gives you:

- where the line is
- how strong the line is at the reference temperature
- what broadening and shift parameters to use

This is why HAPI can take HITRAN rows and go directly into spectral calculations such as:

- absorption coefficient
- transmittance
- absorption spectrum

without first reconstructing line centers from a separate state table.

## What `sw` Means

The HITRAN line intensity `sw` is not temperature-free.

It is the line intensity at the reference temperature, usually `296 K`.

In this repo, the header description in [CH4_M6_I1.header](/f:/GitHub/hapi/hitran_db/CH4_M6_I1.header) states that `sw` is the line intensity at `T = 296 K`.

So:

- HITRAN does not store intensity for every possible temperature
- it stores a standard reference intensity
- HAPI rescales that intensity when the target temperature is different from `296 K`

## How Temperature Affects a HITRAN Spectrum

If the target temperature is not `296 K`, the line intensity must be adjusted.

That scaling depends on:

- the stored reference intensity `sw`
- lower-state energy
- partition function ratio
- stimulated-emission correction

So HITRAN still needs partition-function information when calculating spectra at temperatures other than `296 K`.

Pressure mainly affects:

- line broadening
- line shift

Pressure does not usually redefine the intrinsic line intensity itself.

## How HITRAN Is Used in This Repo

The main HITRAN workflow in this repo is in [plot_vibrational_mode_progressions.py](/f:/GitHub/hapi/plot_vibrational_mode_progressions.py).

That script:

1. opens a local HITRAN table with HAPI
2. reads the global quantum labels
3. groups rows by exact vibrational transition labels
4. builds temporary HAPI tables for selected band groups
5. asks HAPI to compute an absorption coefficient with a Voigt profile
6. converts that result to transmittance
7. converts transmittance to absorbance

## How the Repo Detects Bands From HITRAN Rows

Band classification in this repo is based on the quantum-label fields, not on the line-center value alone.

For example:

- [extract_ch4_nu3_band_texts.py](/f:/GitHub/hapi/extract_ch4_nu3_band_texts.py) slices `global_lower_quanta` and `global_upper_quanta` from the fixed-width HITRAN text file and builds band labels like:
  `lower -> upper`
- [plot_vibrational_mode_progressions.py](/f:/GitHub/hapi/plot_vibrational_mode_progressions.py) parses those quanta labels to determine mode progression categories such as:
  `fundamental_band`
  `first_overtone_band`
  `hot_fundamental_band`

So in this repo:

- `nu` tells you where a line is on the spectral axis
- the quantum-label columns tell you which vibrational band the line belongs to

## How HAPI Builds the Spectrum

In [plot_vibrational_mode_progressions.py](/f:/GitHub/hapi/plot_vibrational_mode_progressions.py), the main spectrum calculation is:

- `hapi.absorptionCoefficient_Voigt(...)`
- `hapi.transmittanceSpectrum(...)`

inside `render_absorbance_table(...)`.

Conceptually, the process is:

1. take each HITRAN row
2. start from its line center `nu`
3. start from its line intensity `sw`
4. rescale the intensity for the target temperature
5. apply a line profile, here Voigt
6. sum all line contributions onto the spectral grid
7. return the absorption coefficient

This is standard line-by-line spectroscopy.

## Absorption Coefficient, Transmittance, and Absorbance

There is an important distinction between these quantities.

### Absorption coefficient

HAPI computes the absorption coefficient using functions such as:

- `absorptionCoefficient_Voigt(...)`

in [hapi.py](/f:/GitHub/hapi/hapi/hapi.py).

### Transmittance

HAPI can convert the absorption coefficient to transmittance using:

- `transmittanceSpectrum(...)`

which computes:

```text
T = exp(-kL)
```

This function exists in [hapi.py](/f:/GitHub/hapi/hapi/hapi.py).

### Absorption spectrum

The local HAPI copy also has:

- `absorptionSpectrum(...)`

which computes:

```text
1 - exp(-kL)
```

This is absorption fraction, not the same as spectroscopic absorbance.

### Absorbance

In this repo, absorbance is computed manually after transmittance:

```text
A = -ln(T)
```

That is what [plot_vibrational_mode_progressions.py](/f:/GitHub/hapi/plot_vibrational_mode_progressions.py) does in `render_absorbance_table(...)`.

So the accurate statement is:

- HAPI gives you absorption coefficient and transmittance
- this repo computes absorbance from transmittance

## Why HITRAN Rows Are Usually Sorted by Wavenumber

Because HITRAN is a line list, it is naturally organized along the spectral axis.

In practice, the rows are typically ordered by increasing `nu`, meaning:

- smaller line-center wavenumbers first
- larger line-center wavenumbers later

So it is reasonable to think of a HITRAN file as moving from lower spectral position to higher spectral position.

## What `extract_ch4_nu3_band_texts.py` Is Doing

[extract_ch4_nu3_band_texts.py](/f:/GitHub/hapi/extract_ch4_nu3_band_texts.py) does not calculate spectra.

It does something simpler:

- reads the fixed-width HITRAN text file
- extracts `global_lower_quanta` and `global_upper_quanta`
- builds exact band labels
- copies matching rows into separate text files
- checks the exported row count against the expected `line_count`

In that script:

- `band_label` is used for row selection
- `line_count` is used only for verification

There, `line_count` is effectively the number of transition rows written for that band.

## Standard HITRAN Workflow

The standard process for using HITRAN data is:

1. choose the molecule and spectral range
2. read the HITRAN line list
3. set the case conditions:
   temperature, pressure, composition, path length
4. rescale line intensity from the reference temperature if needed
5. apply a line profile such as Voigt
6. sum all broadened lines onto the grid
7. convert the result to absorption coefficient, transmittance, or absorbance

This is the same downstream spectral physics used for ExoMol once line center and line intensity are known.

## How HITRAN and ExoMol Differ

The key difference is where the workflow starts.

For HITRAN:

- line center is already stored as `nu`
- reference line intensity is already stored as `sw`

For ExoMol:

- line center is reconstructed from state energies
- line intensity is reconstructed from Einstein-A and thermodynamic quantities

After that point, the downstream line-shape and grid-summation logic is basically the same for both databases.

## Practical Summary

If you want the shortest possible mental model:

- HITRAN gives you a line list that already looks like spectroscopy input
- HAPI takes those lines and calculates the spectrum under your chosen conditions
- in this repo, absorbance is obtained by converting HAPI transmittance with `-ln(T)`
