# Four nu3 Absorbance Scripts

This note is about these four scripts:

- `scripts/plot_exomol_sorted_nu3_absorbance_progressions.py`
- `scripts/plot_hitran_band_text_absorbance_progressions.py`
- `scripts/plot_combined_pure_nu3_absorbance_progressions.py`
- `scripts/plot_combined_exomol_i1_absorbance_progressions.py`

This version is written as direct answers to the main questions in review:

- what is already grouped before plotting
- what "J-pair group" means
- why a label can appear even when the peak is not obvious
- why PNG labels and HTML labels behave differently
- which code controls colors
- how the combined script now merges shared J pairs
- what the combined CSV now records for each merged J pair
- where the case parameters are controlled

## 1. Real inputs

For the current plotting workflow, the real ExoMol input is already preprocessed.

### ExoMol side

The ExoMol and combined scripts read:

- `exomol_ch4_mm_pure_nu3_band_texts_hitran_style_2500_3500_sorted`

That is the practical line-list input folder.

The code does not rescan the full ExoMol transition database during plotting.

It only reads one raw ExoMol metadata file:

- `exomol_db/CH4/12C-1H4/MM/12C-1H4__MM.def`

That file is used to get:

- molecular mass
- default Lorentz width `gamma0`
- temperature exponent `n_exponent`

Those numbers are needed for ExoMol absorbance modeling.

### HITRAN side

The HITRAN and combined scripts read line rows from:

- `ch4_nu3_progressions/band_line_texts`

But they still need:

- `hitran_db/CH4_M6_I1.header`
- `hitran_db/CH4_M6_I1.data`

Reason:

- the `band_line_texts` folder gives the filtered line rows
- HAPI still needs the original schema/bootstrap table files
- the code copies those files into a runtime folder, bootstraps the table schema, and then rebuilds one temporary table per J pair from the text rows

So your understanding is correct:

- the band text folder is the main line input
- `hitran_db` is still needed for metadata/runtime setup

## 2. What is already grouped before plotting

This is the first important correction.

At the band level, the grouping is already finished before these three scripts run.

In practical terms:

- one `.txt` file already represents one band export
- the plotting scripts do not redo band discovery

What the plotting scripts still do is:

1. gather the `.txt` files that belong to the same `nu3 X->Y` progression
2. read the rows inside those files
3. separate those rows into J-pair buckets
4. separate those J-pair curves into the three `delta J` classes

So the plot-level grouping is:

- figure level: one `nu3` progression
- curve level: one J pair
- panel level: one `delta J` class

For the review comment "different `.txt` means different band", yes, that is true.

The extra grouping done by the plotting code is not another band grouping. It is only:

- J-pair grouping inside the already prepared band/progression inputs
- panel split by `delta J = -1, 0, +1`

## 3. What a "J-pair group" means

A J-pair group means:

- all rows in the current progression whose lower and upper rotational quantum numbers are the same

Example:

- every row with `lower J = 2` and `upper J = 3` belongs to the J-pair group `(2, 3)`
- that group becomes one plotted curve labeled `J 2->3`
- since `3 - 2 = +1`, that curve belongs to the `delta J = +1` panel

Another example:

- every row with `lower J = 2` and `upper J = 2` belongs to the J-pair group `(2, 2)`
- that becomes one plotted curve labeled `J 2->2`
- since `2 - 2 = 0`, it goes into the `delta J = 0` panel

So "J-pair group" does not mean a new band.

It means:

- a subset of rows inside the current plotted progression
- grouped only by the pair `(lower_j, upper_j)`

## 4. Which module computes absorbance

You asked "how to keep, how to compute, which module". That is the right question.

### ExoMol-only script

Wrapper:

- `scripts/plot_exomol_sorted_nu3_absorbance_progressions.py`

Main function:

- `research.exomol.plot_sorted_nu3_absorbance_progressions`

Absorbance computation:

- `research.absorbance.render_absorbance_on_grid`

That function uses:

- `research.absorbance.render_cross_section_from_lines`
- `research.spectra.cross_section_to_absorbance`

This is a direct line-by-line absorbance modeling path from:

- line centers
- line intensities
- broadening parameters
- gas conditions

### HITRAN-only script

Wrapper:

- `scripts/plot_hitran_band_text_absorbance_progressions.py`

Main function:

- `research.hitran.plot_band_text_absorbance_progressions`

Absorbance computation:

- HAPI Voigt path

For each J pair, the code:

1. rebuilds a temporary HAPI table from the raw rows
2. calls `hapi.absorptionCoefficient_Voigt(...)`
3. calls `hapi.transmittanceSpectrum(...)`
4. converts transmittance to absorbance with `research.spectra.to_absorbance`

### Combined script

Wrapper:

- `scripts/plot_combined_pure_nu3_absorbance_progressions.py`

Main function:

- `research.combined.plot_combined_pure_nu3_absorbance_progressions`

Absorbance computation:

- ExoMol contribution uses the ExoMol modeling path
- HITRAN contribution uses the HAPI Voigt path

So "combined (pure nu3)" now means one merged J-pair curve.

It means:

- if the same J pair exists in both datasets, model both contributions
- take the pointwise maximum of those two contributions into one uniform J-pair curve
- if the J pair exists in only one dataset, use that one contribution alone

This is closer to what you described in review:

- the key is not to choose one source
- the key is to make one uniform J pair so the peak is less likely to disappear when one dataset misses part of the set

### Combined ExoMol MM I1 script

Wrapper:

- `scripts/plot_combined_exomol_i1_absorbance_progressions.py`

Main function:

- `research.combined.plot_combined_exomol_i1_absorbance_progressions`

ExoMol input:

- `exomol_ch4_mm_i1_pure_nu3_band_texts_hitran_style`
- header: `hitran_exomolCH4_db/CH4_EXOMOL_MM_I1.header`

This is a different ExoMol source than the combined pure nu3 script. The MM I1 source uses fixed-width HITRAN-style band texts extracted from ExoMol with isotopologue I1.

Absorbance computation:

- both ExoMol MM I1 and HITRAN contributions use the same HAPI Voigt path
- for each J pair, the code rebuilds a temporary HAPI table from raw rows, then calls `hapi.absorptionCoefficient_Voigt`
- this is different from the combined pure nu3 script, where ExoMol uses `render_absorbance_on_grid` (direct Voigt sum)

Key difference: separate intensity thresholds.

The ExoMol MM I1 source has much weaker reference-temperature line intensities for hot bands (nu3 1->2, 2->3, 3->4), often around 1e-26 to 1e-25 cm/molecule. HAPI's `IntensityThreshold` filter operates on reference-temperature intensities. If a single threshold is used, hot-band lines get filtered out entirely, producing zero absorbance even at high temperatures where those bands would be physically significant.

To solve this, the function takes two separate thresholds:

- `intensity_threshold` (default `0.0`) — for ExoMol MM I1 temp-table rendering
- `hitran_intensity_threshold` (default `1e-23`) — for HITRAN temp-table rendering

CLI arguments:

- `--intensity-threshold` controls ExoMol (default `0.0`)
- `--hitran-intensity-threshold` controls HITRAN (default `1e-23`)

The merge rule is the same as the combined pure nu3 script: pointwise maximum when both sources exist.

PNG subplot colors:

- all three delta J panels use the same branch color as each other:
  - delta J = -1: blue
  - delta J = 0: orange
  - delta J = +1: green
- this matches the overlay convention so traces are visually consistent across panels

## 5. Why a label can appear even if there is no obvious peak

This is the main behavior problem.

The current code does not ask:

- "is there a clearly visible peak?"
- "is the peak prominent enough?"
- "does this J pair create a distinguishable feature in the figure?"

Instead, the current code asks only:

1. does this J pair have rows
2. was a curve built for this J pair
3. was this curve selected for labeling

That is different.

### What label presence really means today

If the figure shows label `J 2->3`, it means:

- there was a J-pair bucket `(2, 3)`
- the code computed a curve for that bucket
- that curve was selected by the label logic

It does not mean:

- `J 2->3` makes a strong isolated visible peak

### Why this happens

For every selected trace, the code finds a label anchor by taking:

- the maximum point of the curve with `np.argmax(...)`

That always returns some index for a non-empty curve.

So every non-empty modeled curve has a "peak location" in the code, even if:

- the curve is extremely weak
- the curve is buried under stronger traces
- the curve is broad and visually unclear
- the y-axis scale hides it

### Forced labels still make this stronger

The ExoMol and combined scripts have built-in forced labels:

- `J 2->3`
- `J 3->4`

Those are defined by:

- `DEFAULT_FORCED_ABSORBANCE_J_PAIRS = ((2, 3), (3, 4))`

So if those traces exist, the label logic can add them even when they are not among the strongest traces.

That is why your observation is still possible:

- the `.txt` rows contain `J 2->3`
- the merged J-pair curve exists
- but the visible peak is still weak or hard to see
- the label can still appear

So the first question to answer before any fix is:

- should a label mean "this J pair exists in the data"
- or should a label mean "this J pair creates a visible peak under the current plotting conditions"

Right now the code uses the first meaning, not the second.

## 6. PNG label behavior versus HTML label behavior

You pointed out an important difference:

- PNG labels are on the right side, which you think is correct
- HTML labels are not behaving the same way

I agree this is a separate problem from the label-selection logic.

### Intended behavior

The code computes label positions in `_label_positions(...)`.

That function pushes label text toward the right side by using a fixed `anchor_x` inside the plotting window.

So the intended label style is:

- peak point near the actual peak
- text placed on the right side
- a connector from peak to label text

### PNG behavior

The PNG writers use Matplotlib text and connector lines.

That is why the PNG version can look correct.

### HTML behavior

The HTML writers use Plotly annotations.

In practice, your observation is:

- the HTML text is not staying on the right side the same way as the PNG
- it can appear on or near the peak itself

So there are really two different issues:

1. label selection:
   why is this J pair labeled at all
2. HTML placement:
   why is the text not displayed on the right like the PNG

Those should be debugged separately.

For the documentation, the important conclusion is:

- the current HTML output is not matching the intended right-side label style seen in the PNG output

## 7. Which code controls the colors

You asked about the color transition.

The responsible function is:

- `_color_for_index(...)`

It appears in:

- `research.exomol`
- `research.hitran`

What it does:

- it chooses one color per J-pair trace from the Matplotlib `turbo` colormap
- the color depends on the trace index among all traces in that figure

So the colors do not change continuously with wavelength along one single line.

Instead:

- each J-pair curve gets one fixed color
- the collection of curves shows a sweep through the colormap

If you like the way the colors change across traces, `_color_for_index(...)` is the code responsible for it.

## 8. How the combined script now handles one shared J pair

This is the first problem we solved.

The combined script no longer chooses one winning source for a shared J pair.

Now it does this:

1. check whether ExoMol has the J pair
2. check whether HITRAN has the J pair
3. model the ExoMol contribution if present
4. model the HITRAN contribution if present
5. take the pointwise maximum of the two absorbance contributions into one uniform J-pair curve

### Example 1

If both datasets contain `J 2->3`:

- ExoMol contributes to that J pair
- HITRAN contributes to that J pair
- the plotted `J 2->3` curve is the pointwise maximum of both modeled contributions

### Example 2

If only ExoMol contains `J 2->3`:

- the plotted `J 2->3` curve is ExoMol only

### Example 3

If only HITRAN contains `J 2->3`:

- the plotted `J 2->3` curve is HITRAN only

So the current combined logic is now:

- merged when both exist
- single-source only when just one source exists

## 9. What the combined CSV now records

For each J pair in the combined workflow, the CSV now records:

- `source_mode`
- `available_sources`
- `peak_exomol_line_intensity`
- `peak_hitran_line_intensity`
- `total_exomol_line_intensity`
- `total_hitran_line_intensity`
- `total_intensity_ratio_exomol_to_hitran`
- `exomol_line_count`
- `hitran_line_count`
- `line_count`
- `exomol_source_file_count`
- `hitran_source_file_count`
- `peak_wavenumber_cm-1`
- `peak_source`
- `peak_exomol_absorbance`
- `peak_hitran_absorbance`

This is the direct answer to your review requirement:

- how many lines came from ExoMol
- how many lines came from HITRAN
- how strong the strongest ExoMol line is in that J pair
- how strong the strongest HITRAN line is in that J pair
- how the total ExoMol line strength compares with the total HITRAN line strength
- which source dominates at the plotted peak

Important detail:

- `peak_source` means the source contributing more absorbance at the combined peak wavenumber
- `peak_source` can be `exomol`, `hitran`, `both`, or `none`
- it is not a raw database row ID

The new intensity fields are line-list summaries, not curve summaries:

- `peak_exomol_line_intensity`
  - the strongest single ExoMol line intensity inside that J pair
- `peak_hitran_line_intensity`
  - the strongest single HITRAN line intensity inside that J pair
- `total_exomol_line_intensity`
  - the sum of all grouped ExoMol line intensities inside that J pair
- `total_hitran_line_intensity`
  - the sum of all grouped HITRAN line intensities inside that J pair
- `total_intensity_ratio_exomol_to_hitran`
  - `total_exomol_line_intensity / total_hitran_line_intensity`

Behavior rules:

- if a source is missing for that J pair, its `peak_*_line_intensity` and `total_*_line_intensity` values are `0.0`
- if HITRAN is missing or its total intensity is not positive, `total_intensity_ratio_exomol_to_hitran` is written as `nan`

This gives two different comparison layers in the same CSV:

- absorbance comparison:
  - `peak_exomol_absorbance`
  - `peak_hitran_absorbance`
  - `peak_absorbance`
- line-list intensity comparison:
  - `peak_exomol_line_intensity`
  - `peak_hitran_line_intensity`
  - `total_exomol_line_intensity`
  - `total_hitran_line_intensity`
  - `total_intensity_ratio_exomol_to_hitran`

So the practical reading rule is:

- use `peak_*_absorbance` when you want to compare modeled spectra
- use `peak_*_line_intensity` and `total_*_line_intensity` when you want to compare the underlying grouped line strengths

### Important comparison rule

When you compare HITRAN-only output against combined output, do not compare:

- HITRAN-only `peak_absorbance`
with
- combined `peak_absorbance`

unless both runs use the same simulation case and the same `wn_step`.

Reason:

- in the HITRAN-only workflow, `peak_absorbance` means the peak of the HITRAN curve
- in the combined workflow, `peak_absorbance` means the peak of the merged curve `max(exomol, hitran)`

So for source-to-source comparison, the more correct comparison is:

- HITRAN-only `peak_absorbance`
with
- combined `peak_hitran_absorbance`

## 10. What we learned from the direct comparison experiments

We ran direct comparison experiments for:

- `nu3 0->1`
- `J 2->3`
- `T = 600 K`
- `P = 3 Torr`
- `x = 0.008`
- `L = 100 cm`

The output folders were named to include the main parameters:

- `artifacts/hitran_band_text_absorbance_cmp_T600K_P3Torr_x0p008_L100cm_step0p01_nu3_0to1`
- `artifacts/hitran_band_text_absorbance_cmp_T600K_P3Torr_x0p008_L100cm_step0p1_nu3_0to1`
- `artifacts/combined_pure_nu3_absorbance_cmp_T600K_P3Torr_x0p008_L100cm_step0p01_minI0_hitranI1e23`
- `artifacts/combined_pure_nu3_absorbance_cmp_T600K_P3Torr_x0p008_L100cm_step0p1_minI0_hitranI1e23`

### Experiment result for `J 2->3`

At `wn_step = 0.01`:

- HITRAN-only `peak_absorbance = 8.196447e-02`
- combined `peak_hitran_absorbance = 8.196445e-02`
- combined `peak_exomol_absorbance = 2.681899e-01`
- combined `peak_absorbance = 2.681899e-01`

At `wn_step = 0.1`:

- HITRAN-only `peak_absorbance = 8.879971e-07`
- combined `peak_hitran_absorbance = 8.879971e-07`
- combined `peak_exomol_absorbance = 2.715447e-04`
- combined `peak_absorbance = 2.715447e-04`

### What this proves

First:

- the HITRAN contribution inside the combined workflow is consistent with the HITRAN-only workflow when the case and `wn_step` are matched

Second:

- the large mismatch we saw earlier was not because the combined code was using the wrong HITRAN rows
- it was mainly because the comparison used different `wn_step` values and different CSV columns

Third:

- `wn_step = 0.1 cm^-1` is too coarse for this narrow `J 2->3` feature
- the peak collapses by many orders of magnitude compared with `wn_step = 0.01 cm^-1`

So for narrow-line comparison, the practical rule is:

- use the same `wn_step`
- prefer `0.01 cm^-1` or finer
- compare HITRAN-only `peak_absorbance` against combined `peak_hitran_absorbance`

## 11. Why two datasets are still important

Your comment here is important.

If some ExoMol J-pair curves do not show a clear visible peak, that is exactly why the second dataset can matter.

The current combined code now uses the two datasets in a more useful way:

- not as two separate competing plotted curves
- but as two contributions to the same plotted J pair when both exist

So the practical benefit now is:

- if one dataset is missing important lines for a J pair, the other dataset can still contribute to the same final curve
- the CSV still tells us how much each source contributed

This is still not a side-by-side comparison plot, but it is now much closer to the "uniform J pair" idea from your comment.

## 12. Where case parameters are controlled

You asked where the user can control the case.

The answer depends on which script.

### ExoMol-only script

User can control:

- `--wn-min`
- `--wn-max`
- `--wn-step`
- `--pressure-torr`
- `--mole-fraction`
- `--path-length-cm`
- `--line-cutoff`
- `--min-line-intensity`
- label-count settings

But temperature is not exposed there.

Reason:

- the ExoMol line strengths in this workflow are already exported at `296 K`
- the core function hard-codes `temperature_k = 296.0`

### HITRAN-only script

User can control:

- `--wn-min`
- `--wn-max`
- `--wn-step`
- `--temperature-k`
- `--pressure-torr`
- `--mole-fraction`
- `--path-length-cm`
- `--intensity-threshold`
- label-count settings

So HITRAN already exposes the full case more directly.

### Combined pure nu3 script

User can control:

- `--wn-min`
- `--wn-max`
- `--wn-step`
- `--pressure-torr`
- `--mole-fraction`
- `--path-length-cm`
- `--line-cutoff`
- `--min-line-intensity`
- `--hitran-intensity-threshold`
- `--sources`
- label-count settings

The combined script now also exposes:

- `--temperature-k`

So the combined workflow now takes a full simulation case from CLI, like the HITRAN absorbance workflow.

Important detail:

- the same CLI case is applied to both ExoMol and HITRAN inside the combined run
- this includes temperature, pressure, mole fraction, and path length
- the existing ExoMol exported text files still carry `T296.0K` in their filenames, so changing `--temperature-k` changes the modeled absorbance case without changing those source file names

### Combined ExoMol MM I1 script

User can control:

- `--wn-min`
- `--wn-max`
- `--wn-step`
- `--temperature-k`
- `--pressure-torr`
- `--mole-fraction`
- `--path-length-cm`
- `--intensity-threshold` (ExoMol, default `0.0`)
- `--hitran-intensity-threshold` (HITRAN, default `1e-23`)
- `--sources`
- label-count settings

The two intensity thresholds are independent:

- `--intensity-threshold` only affects ExoMol MM I1 lines going through HAPI
- `--hitran-intensity-threshold` only affects HITRAN lines going through HAPI
- setting ExoMol to `0.0` ensures weak hot-band lines are not filtered out at reference temperature

## 13. The most important current discussion point

The first issue to settle is this one:

- why can a label appear when there is no visible peak

The code-level answer is now clear:

- because label selection is based on J-pair existence plus ranking plus forced labels
- not on visible-peak detection

So if we want the figures to behave more like human interpretation, the next design question is:

- should the label logic include a visibility test or prominence threshold

That is the real decision point before fixing the code.

## 14. Practical code map

If we continue the discussion, these are the most important code areas to read next.

### J-pair grouping

- `research.exomol._collect_sorted_progression_groups`
- `research.hitran._parse_band_text_groups`

### Absorbance modeling

- `research.absorbance.render_absorbance_on_grid`
- `research.hitran.plot_band_text_absorbance_progressions`
- `research.combined.plot_combined_pure_nu3_absorbance_progressions`

### Label selection

- `research.exomol._branch_label_candidates`
- `research.hitran._branch_label_candidates`
- `DEFAULT_FORCED_ABSORBANCE_J_PAIRS`

### Label placement

- `research.exomol._label_positions`
- `research.hitran._label_positions`
- the HTML writer functions in `research.exomol`, `research.hitran`, and `research.combined`

### Color assignment

- `_color_for_index(...)`

### Combined source merge logic

- merged-source J-pair logic inside `research.combined.plot_combined_pure_nu3_absorbance_progressions`
- merged-source J-pair logic inside `research.combined.plot_combined_exomol_i1_absorbance_progressions`

If needed, the next pass should be line-by-line on:

- why `J 2->3` is labeled
- why HTML labels are not staying on the right side
