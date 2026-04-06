# research/ — Workflow Package (Active Development Target)

## OVERVIEW

Repo-specific Python package built on top of `hapi`. Functions and small dataclasses, not large orchestration classes. This is where all new research code should go.

## DESIGN RULES (from `docs/framework.md`)

- Functions and small dataclasses as default abstraction
- Pure functions where practical; file I/O explicit and optional
- Return Python objects; caller decides whether to save
- No module-level mutable state or hidden globals
- Use `pathlib.Path`, dataclasses, built-in collections, standard exceptions
- Descriptive domain names: `extract_bands`, `render_absorbance`, `build_summary`
- Avoid: `manager`, `engine`, `processor`, `pipeline` unless technically accurate

## MODULE MAP

| Module | Role | Key Functions |
|--------|------|---------------|
| `models.py` | Shared dataclasses | `GasCase`, `SpectralWindow`, `DatasetPaths`, `BandSelection`, `SpectrumResult`, `BandExportResult`, `SummaryResult`, `ComparisonResult` |
| `io.py` | Path/file helpers | `default_paths()`, `ensure_directory()`, `iter_bz2_text_lines()`, CSV/HTML/MD writers |
| `spectra.py` | Grid + profiles | `build_grid()`, `doppler_hwhm_cm()`, `voigt_profile_cm()`, `number_density()`, `save_spectrum_csv/html()` |
| `absorbance.py` | Cross-section rendering | `render_cross_section_from_lines()`, `render_absorbance_on_grid()`, `render_absorbance_from_lines()` |
| `exomol.py` | ExoMol workflows | `parse_def_file()`, `load_partition_function()`, `collect_relevant_transitions()`, `render_absorbance()`, `scan_band_types()`, `plot_sorted_nu3_*()` |
| `hitran.py` | HITRAN workflows | `fetch_tables()`, `load_table()`, `extract_bands()`, `render_absorbance()`, `plot_band_text_absorbance_progressions()` |
| `combined.py` | Merge workflows | `plot_combined_pure_nu3_absorbance_progressions()`, `plot_combined_exomol_i1_absorbance_progressions()` |
| `bands.py` | Band utilities | `clean_quanta_label()`, `safe_label_fragment()`, `build_summary()` |
| `compare.py` | Source comparison | `nearest_lines()`, report writer |
| `compare_cases.py` | Case comparison | `load_raw_case()`, `plot_dual_case_comparison()` |

## TYPICAL USAGE

```python
from research.models import DatasetPaths, GasCase, SpectralWindow
from research.exomol import render_absorbance

paths = DatasetPaths(root_dir=Path("."), ...)
case = GasCase(temperature_k=600.0, pressure_torr=3.0, mole_fraction=0.008, path_length_cm=100.0)
window = SpectralWindow(wn_min=3000.0, wn_max=3010.0, wn_step=0.01)
result = render_absorbance(paths=paths, case=case, window=window)
```

## TWO RENDERING PATHS

- **ExoMol path**: `research/spectra.py` + `research/absorbance.py` (pure-Python Voigt via `scipy.special.wofz`)
- **HITRAN path**: `research/hitran.py` delegates to `hapi.absorptionCoefficient_Voigt()` + `hapi.transmittanceSpectrum()`
- **Combined**: `research/combined.py` merges both per J-pair (pointwise maximum where both sources exist)

## ANTI-PATTERNS

- Do not put low-level spectroscopy code here; that belongs in `hapi/`
- Do not add large orchestration classes
- Do not encode every default into function names (bad: `build_exomol_ch4_mm_hitran_db`)
- Do not let exploratory scripts dictate package structure

## GOTCHAS

- `__init__.py` lazy-loads submodules via `__getattr__` — import errors surface on first access, not at package import
- HITRAN workflows create/drop temporary HAPI tables (`_runtime_db/`); interrupted runs may leave remnants
- `combined.py` saves raw `.npz` + sidecar CSV metadata for cross-case comparison
- ExoMol `collect_relevant_transitions()` warns if line cutoff is reduced due to missing transition chunks
