# Absorbance Module Validation

- Validation target: `nu3 0->1` from the sorted ExoMol folder
- Comparison: new `research.absorbance` module vs existing `research.exomol.render_cross_section` path
- Note: the sorted line intensities are exported at `296 K`, so the validation case keeps `T = 296 K`
- Window: `2500` to `3500 cm^-1` with `step = 0.25 cm^-1`
- Broadening cutoff: `0.5 cm^-1`
- Minimum line intensity kept: `0.000e+00 cm/molecule`

- Metrics CSV: [validation_metrics.csv](validation_metrics.csv)
- Validation HTML: [validation_nu3_0_to_1.html](validation_nu3_0_to_1.html)

## Metrics

- progression_label: `nu3 0->1`
- line_count: `55558`
- grid_point_count: `4001`
- max_abs_difference: `0.000000000000e+00`
- mean_abs_difference: `0.000000000000e+00`
- max_relative_difference: `0.000000000000e+00`
- mean_relative_difference: `0.000000000000e+00`
- legacy_peak_absorbance: `7.494702898863e-01`
- module_peak_absorbance: `7.494702898863e-01`
- legacy_peak_wavenumber_cm-1: `3016.500000`
- module_peak_wavenumber_cm-1: `3016.500000`
- peak_wavenumber_delta_cm-1: `0.000000e+00`

![Validation overlay](validation_nu3_0_to_1.png)
