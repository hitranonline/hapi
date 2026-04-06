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
| Build tables | `build_exomol_ch4_mm_hitran_db.py`, `build_*_export_table.py` | Convert/export derived tables |
| Extract bands | `extract_ch4_nu3_band_texts.py`, `extract_exomol_ch4_mm_pure_nu3_band_texts.py` | Parse line data into per-band text files |
| Sort/merge | `sort_*_by_intensity.py`, `merge_exomol_pure_nu3_band_texts_hitran_style.py` | Reorder/combine band exports |
| Plot progressions | `plot_vibrational_mode_progressions.py`, `plot_hitran_band_text_absorbance_progressions.py`, `plot_exomol_*` | Generate PNG/HTML/CSV progression outputs |
| Combined workflows | `plot_combined_*_absorbance_progressions.py` | Merge ExoMol + HITRAN absorbance |
| Compare | `compare_hitran_exomol_nearest_lines.py`, `plot_comparison_two_cases.py` | Source and case comparisons |
| Reports | `report_and_plot_band_intensity_compare.py`, `scan_exomol_ch4_nu3_band_types.py` | Analysis reports |

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
