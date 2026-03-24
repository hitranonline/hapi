# Framework Plan: Python-Standard Research Package

This document defines the target structure for rebuilding the current workflow code as a normal Python package.

It is a design spec, not an implementation guide for modifying the existing scripts in place.

## Summary

The repository currently mixes three things:

- the low-level `hapi/` package
- workflow scripts in `scripts/`
- generated artifacts and local datasets

The future direction should be:

- keep `hapi/` as the low-level spectroscopy dependency
- treat `scripts/` as archived reference implementations
- build a new separate package for repo-specific research workflows
- write new scripts only as thin wrappers around that package

The goal is a package that feels like normal Python:

- import modules and functions
- pass clear arguments or small dataclasses
- get Python objects back
- optionally write files afterwards

This should be readable by humans first. The structure should help someone learn normal Python package habits while still being practical for research work.

## Core Decisions

- `hapi/` remains the underlying spectroscopy engine.
- New workflow code lives outside `hapi/`.
- Old scripts are archived reference, not the target architecture.
- Useful workflows should be rebuilt as functions before any new wrapper scripts are added.
- Naming and module structure should follow ordinary Python conventions.
- Defaults are acceptable when they reduce noise and remain scientifically reasonable.
- Simplicity and clarity matter more than full abstraction or enterprise-style layering.

## Design Rules

The new package should follow these rules.

### General rules

- Use modules and functions as the default abstraction.
- Use small dataclasses only when they make calling code clearer.
- Prefer pure functions or mostly-pure functions where practical.
- Return Python objects first; saving files should be explicit and optional.
- Avoid module-level mutable workflow settings.
- Avoid hidden global state.
- Avoid large orchestration classes unless a real need appears later.
- Prefer `pathlib.Path`, dataclasses, built-in collections, and normal exceptions.

### Naming rules

- Use short domain names that humans already recognize: `hitran`, `exomol`, `bands`, `spectra`, `compare`.
- Prefer names like `extract_bands`, `render_absorbance`, `build_summary`, `scan_band_types`.
- Avoid vague names like `manager`, `engine`, `processor`, or `pipeline` unless they describe a real technical concept.
- Avoid script-style public names such as `build_exomol_ch4_mm_hitran_db` for the new package API.
- Do not encode every default into the function name.

### Usage rules

Normal usage should look like:

1. import a function or module
2. create a small config object or pass clear keyword arguments
3. call the function
4. inspect the returned object
5. optionally save outputs

## Target Package Layout

Use one separate top-level package:

```text
research/
    __init__.py
    io.py
    models.py
    hitran.py
    exomol.py
    bands.py
    spectra.py
    compare.py
```

This layout is intentionally flat. If the package grows later, modules can be split into subpackages, but that is not the default starting point.

### Module responsibilities

#### `research.io`

Responsibilities:

- dataset discovery
- path normalization
- shared file reading and writing helpers
- output-directory helpers
- manifest/CSV/HTML save helpers

Keep this module practical. It should not contain domain physics.

#### `research.models`

Responsibilities:

- shared dataclasses
- result containers
- small configuration objects

This module provides a common vocabulary across the package.

#### `research.hitran`

Responsibilities:

- load or fetch local HITRAN/HAPI tables
- extract HITRAN band subsets
- render HITRAN-based spectra and progression outputs
- wrap HAPI-specific workflow logic in cleaner functions

#### `research.exomol`

Responsibilities:

- convert ExoMol data into local workflow objects
- build HITRAN-style tables from ExoMol
- extract ExoMol bands
- render ExoMol absorbance
- scan ExoMol band types

#### `research.bands`

Responsibilities:

- band grouping
- band labeling
- merge operations
- sort operations
- shared band-summary construction helpers

#### `research.spectra`

Responsibilities:

- build spectral grids
- shared coefficient/transmittance/absorbance helpers
- common render-output helpers
- functions used by both HITRAN and ExoMol workflows

#### `research.compare`

Responsibilities:

- ExoMol-vs-HITRAN band comparison
- nearest-line matching
- report-ready comparison summaries

## Shared Models

The package should define these shared types in `research.models`.

### `GasCase`

Purpose:

- describe thermodynamic and composition conditions for a calculation

Suggested fields:

- `temperature_k`
- `pressure_torr`
- `mole_fraction`
- `path_length_cm`

### `SpectralWindow`

Purpose:

- describe the spectral range and sampling

Suggested fields:

- `wn_min`
- `wn_max`
- `wn_step`

### `DatasetPaths`

Purpose:

- carry resolved locations for local data sources and outputs

Suggested fields:

- `root_dir`
- `hitran_db_dir`
- `exomol_db_dir`
- `artifacts_dir`

Optional workflow-specific fields are acceptable if added later.

### `BandSelection`

Purpose:

- describe how bands are selected or grouped

Suggested fields:

- `species`
- `mode_label`
- `lower_band`
- `upper_band`
- `require_same_other_modes`
- `require_unit_step`

### `SpectrumResult`

Purpose:

- represent an in-memory spectral result

Suggested fields:

- `wavenumber`
- `values`
- `quantity`
- `metadata`

### `BandExportResult`

Purpose:

- represent the result of extracting, merging, or sorting band files

Suggested fields:

- `rows`
- `output_dir`
- `manifest_path`
- `count`

### `SummaryResult`

Purpose:

- represent a generated summary table before or after optional file export

Suggested fields:

- `rows`
- `csv_path`
- `html_path`

### `ComparisonResult`

Purpose:

- represent a comparison between ExoMol and HITRAN results

Suggested fields:

- `band_label`
- `matched_band_label`
- `metrics`
- `report_path`
- `figure_paths`

These dataclasses should stay small. If a type becomes too large, it is a sign the workflow should be split further.

## Public API Surface

The public API should stay compact and descriptive.

### `research.exomol`

Public function family:

- `build_hitran_table(...)`
- `extract_bands(...)`
- `render_absorbance(...)`
- `scan_band_types(...)`

Expected behavior:

- inputs are paths, a `GasCase`, a `SpectralWindow`, or a `BandSelection`
- outputs are result objects, arrays, rows, or paths wrapped in dataclasses
- file writing is optional and explicit

### `research.hitran`

Public function family:

- `fetch_tables(...)`
- `load_table(...)`
- `extract_bands(...)`
- `render_absorbance(...)`
- `render_progressions(...)`

Expected behavior:

- hide direct HAPI workflow details behind clearer functions
- keep HAPI table names and local storage details available but not dominant

### `research.bands`

Public function family:

- `merge_bands(...)`
- `sort_bands(...)`
- `build_summary(...)`
- `format_band_label(...)`

Expected behavior:

- provide source-agnostic helpers where possible
- keep band grouping logic in one place instead of reimplementing it in multiple scripts

### `research.spectra`

Public function family:

- `build_grid(...)`
- `to_absorbance(...)`
- `save_csv(...)`
- `save_html(...)`

Expected behavior:

- contain shared spectrum utilities only
- not absorb domain-specific responsibilities from `research.hitran` or `research.exomol`

### `research.compare`

Public function family:

- `compare_bands(...)`
- `nearest_lines(...)`
- `build_report(...)`

Expected behavior:

- return structured comparison data
- keep plotting and report generation usable without forcing file output

## Legacy Script Migration Map

The existing `scripts/` directory should be treated as archived reference code. It remains useful for:

- understanding current behavior
- checking edge-case handling
- validating rebuilt workflows

It should not define the future public interface.

### Ingest and reference scripts

Archived scripts:

- `download_exomol_ch4_mm.py`
- `download_hitemp_ch4.py`
- `download_hitran_ch4.py`
- `fetch_hydrocarbon_tables.py`

Target module:

- primarily `research.io`
- HITRAN-specific fetch helpers in `research.hitran`
- ExoMol-specific dataset helpers in `research.exomol`

### ExoMol rebuild candidates

Archived scripts:

- `build_exomol_ch4_mm_hitran_db.py`
- `extract_exomol_ch4_mm_pure_nu3_band_texts.py`
- `build_exomol_band_line_export_table.py`
- `plot_exomol_ch4_mm_absorbance.py`
- `scan_exomol_ch4_nu3_band_types.py`

Target modules:

- `research.exomol`
- shared band logic in `research.bands`
- shared spectrum logic in `research.spectra`

### HITRAN rebuild candidates

Archived scripts:

- `extract_ch4_nu3_band_texts.py`
- `render_band_line_text_curves.py`
- `plot_vibrational_mode_progressions.py`
- `plot_methane_vibrational_bands.py`
- `select_co_hot_band.py`

Target modules:

- `research.hitran`
- shared band logic in `research.bands`
- shared spectrum logic in `research.spectra`

### Comparison and reporting rebuild candidates

Archived scripts:

- `report_and_plot_band_intensity_compare.py`
- `compare_hitran_exomol_nearest_lines.py`
- `build_band_line_export_table.py`

Target modules:

- `research.compare`
- shared band-summary helpers in `research.bands`
- shared file-output helpers in `research.io`

### Lower-priority exploratory plotting scripts

Archived scripts:

- `plot_abs_five_species_one_case.py`
- `plot_cases_one_species_per_plot.py`
- `plot_overlay_all_species.py`
- `plot_spectrum_2500_3500.py`

Target modules:

- only rebuild pieces that prove reusable
- do not treat every exploratory script as required public API

## Rebuild Priority

Rebuild in this order.

### 1. Core reusable workflows first

Highest priority:

- ExoMol absorbance workflow
- ExoMol-to-HITRAN table build workflow
- HITRAN band extraction workflow
- HITRAN progression and curve rendering workflow
- shared band labeling and summary logic

These are the most reusable parts and already show script-to-script dependency patterns.

### 2. Comparison and reporting second

Next priority:

- ExoMol-vs-HITRAN comparison
- nearest-line matching
- report-ready summary building

These should sit on top of stable extraction and spectrum functions.

### 3. Exploratory plotting last

Lowest priority:

- one-off exploratory plotting scripts
- multi-species convenience plots
- script-specific visualizations that do not define core reusable workflows

These can remain archived reference longer without harming the package design.

## Wrapper Script Rule

If a new script is added later, it should be thin.

Allowed responsibilities for a new script:

1. parse CLI arguments
2. create config objects
3. call one or two package functions
4. print a short summary or write requested outputs

A new script should not reimplement workflow logic that belongs in the package.

## Example Usage

These are examples of the intended Python style.

### Example 1: render ExoMol absorbance

```python
from pathlib import Path

from research.models import DatasetPaths, GasCase, SpectralWindow
from research.exomol import render_absorbance

paths = DatasetPaths(
    root_dir=Path("."),
    exomol_db_dir=Path("exomol_db"),
    hitran_db_dir=Path("hitran_db"),
    artifacts_dir=Path("artifacts"),
)

case = GasCase(
    temperature_k=600.0,
    pressure_torr=3.0,
    mole_fraction=0.008,
    path_length_cm=100.0,
)

window = SpectralWindow(
    wn_min=3000.0,
    wn_max=3010.0,
    wn_step=0.01,
)

result = render_absorbance(paths=paths, case=case, window=window)
print(result.values[:5])
```

### Example 2: extract HITRAN bands and build a summary

```python
from research.models import BandSelection
from research.hitran import extract_bands
from research.bands import build_summary

selection = BandSelection(
    species="CH4",
    mode_label="nu3",
    require_same_other_modes=True,
    require_unit_step=False,
)

exports = extract_bands(table_name="CH4_M6_I1", selection=selection)
summary = build_summary(exports.rows)
print(summary.rows[0])
```

### Example 3: compare ExoMol and HITRAN outputs

```python
from research.compare import compare_bands

comparison = compare_bands(
    exomol_band="0_0_0_0_1A1_to_0_0_1_0_1A2",
    hitran_table="CH4_M6_I1",
)

print(comparison.metrics)
```

These examples are intentionally simple. They show normal Python package usage rather than chained shell scripts.

## Acceptance Checklist

This framework direction is correct only if all of the following remain true:

- old scripts are clearly treated as archived reference
- the new workflow library is separate from `hapi/`
- the package layout is simple and concrete
- the public API is function-oriented and Python-standard
- shared models provide a common vocabulary without becoming heavy abstractions
- file output remains optional and explicit
- exploratory scripts are not allowed to dictate the package structure
- a human reader can understand how the package should be organized and how to use it

## Non-Goals

This document does not require:

- refactoring the current scripts in place
- moving generated artifact directories immediately
- forcing all exploratory analysis into the first package version
- building an enterprise-style framework

The package should start small, clear, and useful.
