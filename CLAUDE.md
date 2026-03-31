# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Purpose

This is a spectroscopy research workspace for computing molecular absorption spectra using HITRAN and ExoMol line-list databases. The primary molecules of interest are CH4 (methane) and CO, focused on the nu3 vibrational mode (~3000 cm⁻¹).

It serves dual roles:
- An importable Python spectroscopy library (`hapi/`)
- A research workspace with local databases, scripts, and generated outputs

## Running Scripts

All scripts are run from the repository root:

```bash
python scripts/<script_name>.py
```

Some scripts accept `--help`. Older plotting scripts use module-level constants near the top of the file for configuration instead of CLI arguments.

## Package Installation

```bash
pip install -e .
```

This installs both the `hapi` and `research` packages in editable mode.

## No Automated Tests

There is no test suite. Correctness is validated by inspecting output artifacts in `artifacts/` and comparing against known reference spectra.

## Architecture

The codebase has three layers, each with a specific role:

### Layer 1: `hapi/hapi.py` — Low-level spectroscopy engine

A single large file (~54k lines) that is the upstream HITRAN API. It provides:
- `fetch()` — download line data from the HITRAN web server
- `absorptionCoefficient_Voigt()` / `_HT()` / `_Lorentz()` etc. — compute absorption cross-sections using various line profiles
- `transmittanceSpectrum()`, `absorptionSpectrum()` — full spectrum generation
- `partitionSum()` — temperature-dependent partition functions (TIPS 2017/2021/2025)
- SQL-like table operations (`createTable`, `select`, etc.) for the local line database

Do not modify `hapi/hapi.py` for research work; it tracks the upstream hitranonline/hapi repo.

### Layer 2: `research/` — Workflow package (the active development target)

This is the repo-specific Python package built on top of `hapi`. The design rules are:
- Functions and small dataclasses, not large orchestration classes
- Pure functions where practical; file I/O is explicit and optional
- Return Python objects; the caller decides whether to save them

**Module responsibilities:**

| Module | Role |
|---|---|
| `research/models.py` | Shared dataclasses: `GasCase`, `SpectralWindow`, `DatasetPaths`, `BandSelection`, `SpectrumResult`, `BandExportResult`, `SummaryResult` |
| `research/io.py` | Path handling, file discovery, CSV/HTML/Markdown writing helpers |
| `research/spectra.py` | Grid construction, Voigt profiles, physical constants, cross-section → absorbance conversion |
| `research/absorbance.py` | Render cross-sections from individual line lists (Doppler + Lorentz → Voigt sum on grid) |
| `research/hitran.py` | HITRAN table loading, band text extraction, J-pair analysis, progression rendering |
| `research/exomol.py` | ExoMol dataset loading, band text parsing, nu3 progression rendering |
| `research/bands.py` | Source-agnostic band grouping, labeling, sorting, and summary construction |
| `research/compare.py` | ExoMol vs HITRAN band comparison and nearest-line matching |
| `research/combined.py` | Combined ExoMol + HITRAN absorbance workflow (pointwise maximum where both sources exist) |

**Typical usage pattern:**
```python
from research.models import DatasetPaths, GasCase, SpectralWindow
from research.exomol import render_absorbance

paths = DatasetPaths(root_dir=Path("."), ...)
case = GasCase(temperature_k=600.0, pressure_torr=3.0, mole_fraction=0.008, path_length_cm=100.0)
window = SpectralWindow(wn_min=3000.0, wn_max=3010.0, wn_step=0.01)
result = render_absorbance(paths=paths, case=case, window=window)
```

### Layer 3: `scripts/` — Standalone CLI entry points

Scripts are thin wrappers around `research/` functions. Treat existing scripts as archived reference implementations. New scripts should only: parse args, create config objects, call one or two package functions, write requested outputs.

## Data Layout

```
hitran_db/          # Local HITRAN line tables (.header + .data pairs)
exomol_db/          # ExoMol source files (.def metadata + line data)
artifacts/          # Generated outputs (HTML plots, CSV exports, reports)
ch4_nu3_progressions/   # CH4 nu3 vibrational progressions
co_nu1_progressions/    # CO nu1 vibrational progressions
docs/               # Design specs, physics background, workflow notes
```

The `artifacts/` subdirectory names encode the run conditions (e.g. `_T600K_P3Torr_x0p008_L100cm`). Each artifact directory typically contains `.html` plots, `.csv` J-pair tables, `report.md`, and a `_runtime_db/` with the HAPI table snapshot used.

## Key Physics Concepts

- **Cross-section**: σ(ν) = Σᵢ Sᵢ · φᵢ(ν), where Sᵢ is line intensity and φᵢ is the broadened line profile
- **Absorbance**: A = -ln(transmittance) = σ · n · L, where n is number density and L is path length
- **Voigt profile**: convolution of Doppler (temperature-dependent) and Lorentz (pressure-dependent) broadening
- Line intensity is temperature-dependent; the HITRAN reference is 296 K; ExoMol is often different
- Wavenumber units throughout (cm⁻¹); pressure in Torr; temperature in K; path length in cm

## Python Environment

The project uses a venv at `~/hapi-venv` (on the home drive, not the T7):

```bash
~/hapi-venv/bin/python3 scripts/<script_name>.py
```

Install with `--no-build-isolation` since `hapi/hapi.py` imports numpy at setup time:

```bash
~/hapi-venv/bin/pip install --no-build-isolation -e .
```

Dependencies: numpy, matplotlib, plotly, scipy, pandas, setuptools.

## Key Constraints

- **HAPI temperature range**: TIPS2025 partition functions only cover 1–2500 K. Runs above 2500 K will fail.
- **T7 external drive**: The repo lives on `/Volumes/T7` (exFAT). macOS creates `._*` resource fork files on non-APFS volumes — these are gitignored and hidden in VS Code.
- **Artifact output naming**: Output folder names must encode run conditions (e.g. `combined_exomol_i1_absorbance_T600K_P3Torr_x0p008_L100cm_step0p001`).
- **Intensity thresholds differ by source**: In `plot_combined_exomol_i1_absorbance_progressions`, ExoMol uses `--intensity-threshold` (default 0.0) and HITRAN uses `--hitran-intensity-threshold` (default 1e-23). ExoMol must default to 0.0 because hot-band reference-temperature intensities (~1e-26 cm/mol) are far below HITRAN's typical threshold.
- **gitignored directories**: `artifacts/`, `hitran_db/`, `ch4_nu3_progressions/`, `co_nu1_progressions/`, `exomol_db/`, `.vscode/` are all gitignored. If files in these directories show up in git status, they need `git rm -r --cached <dir>` to untrack.

## Combined Workflow Scripts

There are two combined ExoMol+HITRAN absorbance scripts with different rendering paths:

| Script | ExoMol source | ExoMol rendering | Line filter |
|---|---|---|---|
| `plot_combined_pure_nu3_absorbance_progressions.py` | sorted pure-nu3 band texts | `render_absorbance_on_grid` (direct Voigt) | `--min-line-intensity` (default 0.0) |
| `plot_combined_exomol_i1_absorbance_progressions.py` | MM I1 HITRAN-style band texts | HAPI temp-table Voigt | `--intensity-threshold` (default 0.0) |

Both use the HAPI temp-table path for HITRAN lines with `--hitran-intensity-threshold` (default 1e-23).

## PNG Plot Convention

- PNG figures have 4 subplots: dJ=-1, dJ=0, dJ=+1, and an "All ΔJ overlaid" panel
- All subplots use consistent branch colors: blue (dJ=-1), orange (dJ=0), green (dJ=+1)
- The overlay panel shows all J-pair curves colored by branch with a legend

## Git Remotes

- `origin`: fork at Ezeki3lRchy/hapi
- `upstream`: official hitranonline/hapi

Changes to `hapi/hapi.py` should track upstream. All new research code goes in `research/` and `scripts/`.
