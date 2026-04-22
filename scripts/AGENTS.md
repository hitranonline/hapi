<!-- Parent: ../AGENTS.md -->
<!-- Generated: 2026-04-03 | Updated: 2026-04-22 -->

# scripts/ — CLI Entry Points

## OVERVIEW

Standalone Python scripts run from repo root. Thin wrappers around `research/` functions. Treat existing scripts as archived reference implementations.

## HOW TO RUN

```bash
# Always from repo root
python scripts/<script_name>.py
python scripts/<script_name>.py --help   # Most support --help
```

## BOOTSTRAP PATTERN

Every script that needs `hapi` or `research` imports starts with:
```python
from _bootstrap import ensure_repo_root
ensure_repo_root()
```
This inserts repo root into `sys.path` so packages are importable without installation.

## SCRIPT CATEGORIES

| Category | Scripts | Notes |
|----------|---------|-------|
| Download | `download_exomol_ch4_mm.py`, `download_hitemp_ch4.py`, `download_hitran_ch4.py` | Fetch remote data into local DB dirs |
| Build tables | `build_exomol_ch4_mm_hitran_db.py`, `build_band_line_export_table.py`, `build_exomol_band_line_export_table.py`, `build_nu3_absorbance_slides.py` | Convert/export derived tables and presentation slides |
| Extract bands | `extract_ch4_nu3_band_texts.py`, `extract_exomol_ch4_mm_pure_nu3_band_texts.py` | Parse line data into per-band text files |
| Sort / merge | `sort_exomol_hitran_style_band_texts_by_intensity.py`, `sort_hitemp_band_texts_by_intensity.py`, `merge_exomol_pure_nu3_band_texts_hitran_style.py` | Reorder/combine band exports |
| Select | `select_co_hot_band.py` | Subset a HITRAN-style band file to hot-band rows only |
| Plot progressions | `plot_vibrational_mode_progressions.py`, `plot_hitran_band_text_absorbance_progressions.py`, `plot_exomol_direct_nu3_absorbance_progressions.py`, `plot_exomol_sorted_nu3_absorbance_progressions.py`, `plot_exomol_sorted_nu3_intensity_progressions.py` | Per-source PNG / HTML / CSV progression outputs |
| Combined workflows | `plot_combined_pure_nu3_absorbance_progressions.py`, `plot_combined_exomol_i1_absorbance_progressions.py` | Merge ExoMol + HITRAN absorbance |
| Compare | `compare_hitran_exomol_nearest_lines.py`, `plot_comparison_two_cases.py`, `plot_cases_one_species_per_plot.py`, `plot_abs_five_species_one_case.py`, `plot_overlay_all_species.py`, `plot_spectrum_2500_3500.py` | Source, case, and multi-species comparisons |
| ExoMol rendering | `plot_exomol_ch4_mm_absorbance.py`, `render_band_line_text_curves.py`, `plot_methane_vibrational_bands.py` | Render absorbance from ExoMol / band-text line lists |
| Non-LTE / Treanor | `plot_ch4_nu3_nonlte_absorbance.py`, `plot_ch4_nu3_treanor_populations.py`, `validate_ch4_nu3_treanor_equilibrium.py`, `plot_treanor_fig3_3_ev.py`, `plot_treanor_fig3_3_vib_number.py`, `plot_treanor_two_panel_vib_number.py` | Two-temperature vibrational distribution plots + LTE equivalence validation |
| Reports / scans | `report_and_plot_band_intensity_compare.py`, `scan_exomol_ch4_nu3_band_types.py` | Analysis reports and band-type scans |

## NEW SCRIPT RULES

A new script should only:
1. Parse CLI arguments (argparse)
2. Create config objects (`GasCase`, `SpectralWindow`, `DatasetPaths`)
3. Call one or two `research/` package functions
4. Print summary or write outputs

Never reimplement workflow logic that belongs in `research/`.

## GOTCHAS

- Older scripts use **module-level constants** near top-of-file instead of CLI args; edit those to change defaults
- `download_hitemp_ch4.py` requires HITRAN credentials (env vars `HITRAN_EMAIL`/`HITRAN_PASSWORD`)
- Scripts that print "warning: line cutoff reduced" silently change Voigt rendering width; inspect results carefully
- `_bootstrap.py` is a shared helper, not a runnable script
- Treanor plotting scripts (`plot_treanor_fig3_3_*.py`, `plot_treanor_two_panel_vib_number.py`) hard-code `T0`, `Tv`, and `v_max` at top-of-file — edit constants to change the figure parameters
- `validate_ch4_nu3_treanor_equilibrium.py` is an executable regression check: it asserts that `vibrational_temperature_k == temperature_k` reproduces LTE intensities within `RELATIVE_ERROR_LIMIT = 1e-10`. Run it after any change to `research/nonlte.py` or `research/ch4_treanor.py`.
- `plot_ch4_nu3_nonlte_absorbance.py` outputs to `ch4_nu3_treanor/` (not `artifacts/`) — that directory lives at the repo root and holds just the two non-LTE PNG figures
