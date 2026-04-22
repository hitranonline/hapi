<!-- Parent: ../AGENTS.md -->
<!-- Generated: 2026-04-03 | Updated: 2026-04-22 -->

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

### Spectroscopy / LTE core

| Module | Role | Key Functions |
|--------|------|---------------|
| `models.py` | Shared dataclasses | `GasCase`, `SpectralWindow`, `DatasetPaths`, `BandSelection`, `SpectrumResult`, `BandExportResult`, `SummaryResult`, `ComparisonResult` |
| `io.py` | Path/file helpers | `default_paths()`, `ensure_directory()`, `iter_bz2_text_lines()`, CSV/HTML/MD writers |
| `spectra.py` | Grid + profiles | `build_grid()`, `doppler_hwhm_cm()`, `voigt_profile_cm()`, `number_density()`, `save_spectrum_csv/html()`, `cross_section_to_absorbance()` |
| `absorbance.py` | Cross-section rendering | `render_cross_section_from_lines()`, `render_absorbance_on_grid()`, `render_absorbance_from_lines()` |
| `exomol.py` | ExoMol workflows | `parse_def_file()`, `load_partition_function()`, `collect_relevant_transitions()`, `collect_nu3_transitions_by_jpair()`, `render_absorbance()`, `render_cross_section()`, `scan_band_types()`, `plot_sorted_nu3_*()` |
| `hitran.py` | HITRAN workflows | `fetch_tables()`, `load_table()`, `extract_bands()`, `render_absorbance()`, `plot_band_text_absorbance_progressions()` |
| `combined.py` | Merge workflows | `plot_combined_pure_nu3_absorbance_progressions()`, `plot_combined_exomol_i1_absorbance_progressions()` |
| `bands.py` | Band utilities | `clean_quanta_label()`, `safe_label_fragment()`, `build_summary()` |
| `compare.py` | Source comparison | `nearest_lines()`, report writer |
| `compare_cases.py` | Case comparison | `load_raw_case()`, `plot_dual_case_comparison()` |

### Non-equilibrium / two-temperature vibrational distributions

| Module | Role | Key Functions |
|--------|------|---------------|
| `treanor.py` | Generic diatomic Treanor distribution (Fridman 2008, §3.1.8, Eq. 3-37) | `DiatomicConstants`, `CO`, `N2`, `O2`, `ln_treanor()`, `treanor()`, `ln_boltzmann()`, `boltzmann()`, `treanor_minimum()` |
| `ch4_treanor.py` | Effective single-mode CH4 nu3 Treanor (0,0,v3,0 ladder) | `CH4Nu3Constants`, `CH4_NU3`, `energy_v3()`, `ln_treanor_nu3()`, `treanor_nu3()`, `ln_boltzmann_nu3()`, `nonlte_intensity_scale_factor()` |
| `nonlte.py` | Apply Treanor non-LTE scaling to ExoMol line intensities | `collect_nu3_transitions_nonlte()`, `TransitionGroup` / `TransitionGroups` typed dicts |

## TYPICAL USAGE

### LTE absorbance

```python
from research.models import DatasetPaths, GasCase, SpectralWindow
from research.exomol import render_absorbance

paths = DatasetPaths(root_dir=Path("."), ...)
case = GasCase(temperature_k=600.0, pressure_torr=3.0, mole_fraction=0.008, path_length_cm=100.0)
window = SpectralWindow(wn_min=3000.0, wn_max=3010.0, wn_step=0.01)
result = render_absorbance(paths=paths, case=case, window=window)
```

### Non-LTE (two-temperature) absorbance

```python
from research.nonlte import collect_nu3_transitions_nonlte
from research.ch4_treanor import CH4_NU3

groups, meta = collect_nu3_transitions_nonlte(
    data_dir=Path("exomol_db/CH4/12C-1H4/MM"),
    temperature_k=600.0,             # T_0 (translational)
    vibrational_temperature_k=3000.0,  # T_v
    mol=CH4_NU3,
    wn_min=2500.0, wn_max=3500.0,
)
# groups: dict keyed by (v3_lower, v3_upper, J_lower, J_upper) → list of (nu, intensity)
```

## THREE RENDERING PATHS

- **ExoMol LTE path**: `research/spectra.py` + `research/absorbance.py` (pure-Python Voigt via `scipy.special.wofz`)
- **HITRAN LTE path**: `research/hitran.py` delegates to `hapi.absorptionCoefficient_Voigt()` + `hapi.transmittanceSpectrum()`
- **Combined LTE**: `research/combined.py` merges both per J-pair (pointwise maximum where both sources exist)
- **ExoMol non-LTE**: `research/nonlte.py` scales each grouped transition's LTE intensity by the Treanor/LTE population ratio from `ch4_treanor.nonlte_intensity_scale_factor()`

## ANTI-PATTERNS

- Do not put low-level spectroscopy code here; that belongs in `hapi/`
- Do not add large orchestration classes
- Do not encode every default into function names (bad: `build_exomol_ch4_mm_hitran_db`)
- Do not let exploratory scripts dictate package structure
- Do not treat `ch4_treanor.py` as a general polyatomic non-LTE model — it is an **effective single-mode approximation** for the (0,0,v3,0) ladder only. Other vibrational modes are held at v=0.

## GOTCHAS

- `__init__.py` lazy-loads submodules via `__getattr__` — import errors surface on first access, not at package import
- HITRAN workflows create/drop temporary HAPI tables (`_runtime_db/`); interrupted runs may leave remnants
- `combined.py` saves raw `.npz` + sidecar CSV metadata for cross-case comparison
- ExoMol `collect_relevant_transitions()` warns if line cutoff is reduced due to missing transition chunks
- The Treanor distribution is **unnormalized**; `treanor()` returns population ratios N_v / N_0, not absolute populations
- `treanor_minimum()` returns the (possibly non-integer) location of the distribution's minimum in v-space — useful for diagnosing the anharmonic turnover, not for selecting physical states
- Setting `vibrational_temperature_k == temperature_k` in `nonlte.collect_nu3_transitions_nonlte` must reproduce LTE within ~1e-10 relative error — this identity is tested by `scripts/validate_ch4_nu3_treanor_equilibrium.py`

<!-- MANUAL: -->
