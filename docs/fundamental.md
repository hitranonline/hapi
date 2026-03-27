# Fundamentals of Line Intensity, Broadening, and Absorbance

This note explains the basic spectroscopy-modeling ideas behind the local
ExoMol, HITRAN, and combined workflows in this repo.

The main goal is to separate two things that are easy to mix up:

- quantities that belong to individual spectral lines
- quantities that belong to the final modeled spectrum

The core pipeline is:

```text
line list -> broaden each line -> sum on a grid -> convert to absorbance
```

This distinction matters when reading the current CSV outputs and when
comparing ExoMol and HITRAN curves in the combined workflow.

## What a Spectral Line Row Gives You

A spectral line row describes one transition.

Typical line-level quantities are:

- line center `nu`
- line intensity such as `sw`
- lower-state energy `elower`
- broadening parameters such as `gamma_air`, `gamma_self`, `n_air`, `delta_air`

For ExoMol-style reconstructed rows, the same idea still applies even if some
fields were derived rather than stored directly in the original source files.

The important point is:

- line intensity is a property of one transition
- it is not the plotted peak height of the final spectrum

So if one row has a strong intensity, that does not mean the final absorbance
curve will have a peak equal to that intensity value.

In plain language, these parameters mean:

- `nu`
  - where the line is centered on the wavenumber axis
- `sw` or line intensity
  - how strong the transition is overall
  - this is more like the total strength or area of the line, not the final
    peak height of the plotted curve
- `elower`
  - the energy of the lower state
  - this matters when intensity is rescaled with temperature
- `gamma_air` or `gamma0`
  - how wide the collision-broadened line is at reference conditions
- `n_air` or `n_exponent`
  - how that collision broadening changes with temperature
- `delta_air`
  - how much the line center shifts under pressure

## From Line Intensity to Absorbance

The modeled spectrum is not built by plotting raw line intensities directly.

An ideal mathematical line before broadening is often thought of as a spike at
one exact line center:

$$
\sigma_i(\nu) = S_i \, \delta(\nu - \nu_i)
$$

where:

- $S_i$ is the line intensity
- $\nu_i$ is the line center
- $\delta(\nu - \nu_i)$ is an ideal delta-function spike

That ideal spike is not what we plot in the real modeling workflow.

Instead, each line is broadened into a line-shape function around its center,
and then all of those broadened contributions are added together on a
wavenumber grid.

A compact way to write that is:

$$
\sigma(\nu) = \sum_i S_i \, \phi_i(\nu)
$$

where:

- $S_i$ is the intensity of line $i$
- $\phi_i(\nu)$ is the broadened profile of line $i$ at wavenumber $\nu$
- $\sigma(\nu)$ is the modeled cross section at that grid point

The key idea is:

- the line intensity $S_i$ gives the total strength
- the profile $\phi_i(\nu)$ spreads that strength across nearby wavenumbers

So one row gives more than just "a value at one point". It gives:

- a center
- a strength
- enough information to turn that strength into a distribution around the center

For Gaussian, Lorentzian, and Voigt profiles, the support is mathematically not
strictly limited to a finite interval. In theory, the tails extend forever.
But in practical code, the profile is evaluated only on a finite local window
because far-away tails become negligible. In this repo, the ExoMol renderer
uses a finite `line_cutoff` window for that reason.

After the cross section is built, it is converted to absorbance with
Beer-Lambert scaling:

$$
A(\nu) = \sigma(\nu)\, n\, L
$$

where:

- $A(\nu)$ is absorbance
- $n$ is absorber number density
- $L$ is path length

So absorbance at one exact wavenumber is an accumulated curve value. It is not
automatically a direct copy of one line row.

## What the X-Y Plot Means

The x-axis is usually wavenumber.

The y-axis depends on which stage of the workflow you are plotting.

### Line-list view

If you are looking at raw line rows, then you may plot:

- x = line center `nu`
- y = line intensity such as `sw`

That is still a line-list representation, not yet the final modeled spectrum.

### Spectrum after broadening

After broadening and summation, the natural curve is:

- x = wavenumber
- y = cross section $\sigma(\nu)$

Then, after Beer-Lambert scaling, the final absorbance curve is:

- x = wavenumber
- y = absorbance $A(\nu)$

So the workflow is not:

$$
(\text{wavenumber}, \text{cross section}) \rightarrow (\text{wavenumber}, \text{intensity})
$$

It is the other way around:

$$
\text{line intensity} \rightarrow \text{profile-broadened cross section} \rightarrow \text{absorbance}
$$

## Why Two Widths Are Needed

Real gas-phase lines are broadened by more than one physical effect.

### Doppler Width

Doppler broadening comes from thermal motion of molecules.

- hotter gas means faster molecular motion
- motion shifts the apparent line position
- the resulting profile is Gaussian-like

This width depends mainly on:

- line center
- temperature
- molecular mass

### Lorentz Width

Lorentz broadening comes mainly from collisions and finite-state lifetime.

- more pressure usually means more collisions
- collisions broaden the line wings
- the resulting profile is Lorentzian-like

This width depends mainly on:

- pressure
- broadening coefficients such as `gamma_air` or a default `gamma0`
- temperature exponent such as `n_air` or a default `n_exponent`

### Voigt Profile

In practice, both effects matter, so the final line shape is usually modeled as
a Voigt profile:

- Gaussian part from Doppler broadening
- Lorentzian part from collisional broadening

That is why the code needs both width types.

## Why a Peak Is Not Necessarily One Line

One of the most common misunderstandings is:

> one absorbance peak must come from one line

That is not generally true.

Because each line is broadened, one line contributes not only at its exact
center but also to nearby grid points. If multiple lines are close together,
their broadened tails overlap.

So at a peak wavenumber:

- one line may dominate the peak
- or several nearby lines may overlap and create a blended peak

This matters for the combined J-pair CSV fields:

- `peak_exomol_absorbance`
- `peak_hitran_absorbance`
- `peak_absorbance`
- `peak_exomol_line_intensity`
- `peak_hitran_line_intensity`
- `total_exomol_line_intensity`
- `total_hitran_line_intensity`

The absorbance fields are curve values evaluated at a peak wavenumber. They are
not guaranteed to be one-row quantities.

So if a J pair has more than one ExoMol line or more than one HITRAN line, its
peak absorbance may be the result of multiple contributing lines.

By contrast, the `*_line_intensity` fields are line-list summaries:

- `peak_*_line_intensity`
  - strongest single line in that source's J pair
- `total_*_line_intensity`
  - sum of all grouped line intensities in that source's J pair

So the absorbance fields and the line-intensity fields are not interchangeable.

## How This Maps to the Repo

## ExoMol Path

For the current sorted ExoMol absorbance workflow, grouped rows provide:

- `wavenumber`
- `intensity`

The ExoMol absorbance path then:

1. takes those line centers and line intensities
2. computes Doppler and Lorentz widths
3. builds Voigt profiles
4. sums the broadened contributions on the grid
5. converts cross section to absorbance

In the current repo implementation:

- the custom ExoMol renderer uses the passed line intensities directly
- one shared ExoMol default broadening pair is applied to the whole set:
  - `gamma0`
  - `n_exponent`

That means the ExoMol path is using one default Lorentz broadening model for
the selected ExoMol lines in that rendering call.

## HITRAN Path

For the HITRAN path, HAPI works from line-level HITRAN fields such as:

- `sw`
- `gamma_air`
- `gamma_self`
- `n_air`
- `delta_air`
- `elower`

HAPI then:

1. starts from the HITRAN line list
2. rescales intensity from the reference temperature when needed
3. applies line-level broadening and line shift
4. builds the spectral coefficient
5. converts the result to transmittance
6. converts transmittance to absorbance

So HITRAN is closer to a spectroscopy-ready line list, while ExoMol in this
repo is still partially reconstructed and then modeled through a custom path.

## Combined Workflow

The current combined workflow does not yet use one identical modeling path for
both sources.

Instead:

- ExoMol contribution uses the custom ExoMol renderer
- HITRAN contribution uses the HAPI Voigt path

That means even if ExoMol and HITRAN rows have similar line centers and similar
line intensities, the final absorbance curves can still differ because the
downstream modeling inputs are not identical.

Examples of differences include:

- ExoMol uses shared default broadening values such as `gamma0` and `n_exponent`
- HITRAN uses per-line broadening and shift values such as `gamma_air`,
  `gamma_self`, `n_air`, and `delta_air`
- the HITRAN path includes HAPI temperature-dependent intensity scaling
- the ExoMol path in the current workflow uses the provided ExoMol intensities
  directly for that rendering pass

So equal or similar line intensities do not guarantee equal absorbance curves.

## How to Read the Combined J-Pair CSV

The combined J-pair CSV is a curve-summary table, not a raw line-list table.

Important fields include:

- `peak_exomol_absorbance`
- `peak_hitran_absorbance`
- `peak_absorbance`
- `peak_exomol_line_intensity`
- `peak_hitran_line_intensity`
- `total_exomol_line_intensity`
- `total_hitran_line_intensity`
- `total_intensity_ratio_exomol_to_hitran`

These should be read as:

- `peak_exomol_absorbance`
  - ExoMol absorbance value at the combined peak wavenumber
- `peak_hitran_absorbance`
  - HITRAN absorbance value at the combined peak wavenumber
- `peak_absorbance`
  - the final combined peak value for that J pair
- `peak_exomol_line_intensity`
  - strongest single ExoMol line intensity inside that J pair
- `peak_hitran_line_intensity`
  - strongest single HITRAN line intensity inside that J pair
- `total_exomol_line_intensity`
  - sum of all grouped ExoMol line intensities inside that J pair
- `total_hitran_line_intensity`
  - sum of all grouped HITRAN line intensities inside that J pair
- `total_intensity_ratio_exomol_to_hitran`
  - ratio of overall grouped line strength, not ratio of absorbance peaks

Only the `peak_*_absorbance` and `peak_absorbance` fields are evaluated on the
modeled curve. They are not copied directly from one line row.

The `*_line_intensity` and `total_intensity_ratio_exomol_to_hitran` fields are
line-list summary values, not curve values.

So when reading the CSV:

- treat `sw` or ExoMol line intensity as line-list input
- treat `peak_*_absorbance` as modeled spectrum output
- treat `peak_*_line_intensity` and `total_*_line_intensity` as grouped
  line-strength summaries

Those are different stages of the workflow.

## Common Misunderstandings

### “One absorbance peak must come from one line.”

Not necessarily. Peaks can be blended from several broadened lines.

### “Same intensity means same absorbance.”

Not necessarily. The final absorbance also depends on broadening, temperature,
pressure, line overlap, mole fraction, path length, and the modeling path.

### “A symmetry-label mismatch proves the transition is different.”

Not necessarily. A label mismatch may reflect mapping differences,
representation differences, or a genuine assignment issue. It is not by itself
enough to prove that two modeled rows are different physical transitions.

## Short Version

The simplest mental model is:

- a line list tells you where transitions are and how strong they are
- modeling turns each transition into a broadened curve
- the final spectrum is the sum of those curves
- absorbance peaks belong to the modeled spectrum, not automatically to one
  individual line row

That is why a peak value in the combined workflow should be interpreted as a
curve result, not as a direct line-intensity value.
