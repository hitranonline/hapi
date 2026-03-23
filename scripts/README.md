# Scripts

All standalone utilities now live in this directory. Run them from the repository root, for example:

```powershell
python scripts/plot_vibrational_mode_progressions.py
```

## Download and ingest

- `download_exomol_ch4_mm.py`: download ExoMol CH4 MM files
- `download_hitemp_ch4.py`: download HITEMP CH4 data
- `fetch_hydrocarbon_tables.py`: fetch hydrocarbon HITRAN tables into the local HAPI database

## Build derived tables and exports

- `build_exomol_ch4_mm_hitran_db.py`: convert ExoMol CH4 MM data into a HITRAN-style table
- `build_exomol_band_line_export_table.py`: build ExoMol per-band absorbance exports and summaries
- `build_band_line_export_table.py`: build summary tables from previously generated band exports

## Extract, merge, and sort band text files

- `extract_ch4_nu3_band_texts.py`
- `extract_exomol_ch4_mm_pure_nu3_band_texts.py`
- `merge_exomol_pure_nu3_band_texts_hitran_style.py`
- `sort_exomol_hitran_style_band_texts_by_intensity.py`
- `sort_hitemp_band_texts_by_intensity.py`

## Plotting and spectrum generation

- `plot_vibrational_mode_progressions.py`
- `plot_exomol_ch4_mm_absorbance.py`
- `plot_methane_vibrational_bands.py`
- `plot_abs_five_species_one_case.py`
- `plot_cases_one_species_per_plot.py`
- `plot_overlay_all_species.py`
- `plot_spectrum_2500_3500.py`
- `render_band_line_text_curves.py`

## Reporting and analysis helpers

- `report_and_plot_band_intensity_compare.py`
- `scan_exomol_ch4_nu3_band_types.py`
- `select_co_hot_band.py`

## Notes

- The scripts are still a mixed workspace rather than a polished CLI package.
- The first cleanup pass groups them in one place and fixes repo-root path handling.
- Some scripts expose `--help`; older plotting scripts still use constants near the top of the file instead of a full CLI.
- A later pass can split these further by workflow once the interfaces settle down.
