<!-- Generated: 2026-04-03 | Updated: 2026-04-22 -->

# PROJECT KNOWLEDGE BASE

**Branch:** codex/reorganize-repo

## OVERVIEW

Python spectroscopy research workspace: an importable HAPI library (`hapi/`), a research workflow package (`research/`), CLI scripts (`scripts/`), and local HITRAN/ExoMol databases. Primary molecules: CH4 and CO, focused on nu3 vibrational mode (~3000 cm-1). Recent work has extended the workspace beyond pure LTE to cover two-temperature vibrational non-equilibrium (Treanor distribution).

## STRUCTURE

```
./
├── hapi/                   # Upstream HITRAN API engine (DO NOT MODIFY for research)
│   └── AGENTS.md
├── research/               # Active dev target: workflow package built on hapi
│   └── AGENTS.md
├── scripts/                # CLI entry points (thin wrappers around research/)
│   └── AGENTS.md
├── docs/                   # Design specs, physics notes, workflow docs
│   └── AGENTS.md
├── classify/               # Legacy analysis (Python + MATLAB)
│   └── AGENTS.md
├── hitran_db/              # [gitignored] Local HITRAN .header/.data tables (INPUT)
├── exomol_db/              # [gitignored] ExoMol .def/.pf/.states/.trans files (INPUT)
├── artifacts/              # [gitignored] Generated outputs (HTML/CSV/PNG/NPZ)
├── ch4_nu3_progressions/   # [gitignored] Generated CH4 band exports
├── co_nu1_progressions/    # [gitignored] Generated CO band exports
├── band_intensity_comparisons/  # [gitignored] ExoMol vs HITRAN compare reports
├── ch4_nu3_treanor/        # Generated non-LTE figures (2 PNGs)
├── treanor_distribution/   # Generated Treanor textbook-recreation figures (3 PNGs)
├── methane_vibrational_bands/  # Generated: reference CSV + HTML for CH4 bands
├── exomol_ch4_mm_scans/    # Generated: band-type scan CSVs
├── species_plots_2500_3500_multiTP/  # Generated multi-T/P overlay HTMLs
├── presentation/           # Slide deck (.pptx)
├── setup.py                # Package: "hitran-api", includes hapi + research
└── CLAUDE.md               # Authoritative AI coding guidance (read first)
```

## HIERARCHICAL AGENTS.md COVERAGE

Every code-bearing directory has its own `AGENTS.md`. Each child file begins with `<!-- Parent: ../AGENTS.md -->` so the hierarchy is navigable. Pure-output directories (`artifacts/`, `*_db/`, `*_progressions/`, `treanor_distribution/`, `ch4_nu3_treanor/`, etc.) intentionally do not have `AGENTS.md` — their purpose is recorded here in the root map.

## WHERE TO LOOK

| Task | Location | Notes |
|------|----------|-------|
| Understand repo purpose & AI rules | `CLAUDE.md` | Read before any work |
| Design rules for research/ | `docs/framework.md` | Authoritative API style guide |
| Code review standards | `docs/code_review_checklist.md` | AI-generated code validation |
| Repo layout reference | `docs/REPOSITORY_STRUCTURE.md` | What goes where |
| Physics background | `docs/HITRAN_DATABASE_NOTES.md`, `docs/EXOMOL_DATABASE_NOTES.md`, `docs/fundamental.md` | Domain context |
| Hot-band physics | `docs/hot_band_spectral_shape_and_magnitude.md`, `docs/hot_band_temperature_dependence.md` | Why hot bands behave as they do |
| Non-LTE / two-temperature theory | `docs/TREANOR_DISTRIBUTION_NOTES.md` | Treanor (1968), Fridman (2008) |
| Run a workflow | `scripts/README.md`, `scripts/AGENTS.md` | Script inventory and usage |
| Shared dataclasses | `research/models.py` | GasCase, SpectralWindow, DatasetPaths, etc. |
| Path/IO helpers | `research/io.py` | default_paths(), iter_bz2_text_lines() |
| HITRAN rendering | `research/hitran.py` | Calls hapi for Voigt/transmittance |
| ExoMol rendering | `research/exomol.py` | Pure-Python Voigt from line lists |
| Combined workflows | `research/combined.py` | Merge ExoMol + HITRAN per J-pair |
| Treanor distribution (diatomic) | `research/treanor.py` | CO / N2 / O2 constants and population ratios |
| CH4 nu3 non-LTE | `research/ch4_treanor.py`, `research/nonlte.py` | Effective single-mode (0,0,v3,0) ladder |

## THREE-LAYER ARCHITECTURE

```
Layer 3: scripts/        Thin CLI wrappers (argparse -> config -> call -> save)
Layer 2: research/       Workflow functions, dataclasses, domain logic
Layer 1: hapi/hapi.py    Low-level spectroscopy engine (upstream, read-only)
```

- `scripts/` calls `research/` functions; never reimplements workflow logic
- `research/` calls `hapi` for HITRAN computations; implements ExoMol path independently
- `hapi/` tracks upstream `hitranonline/hapi` and must not be modified locally

## CONVENTIONS

- **Run from repo root**: `python scripts/<name>.py`
- **Install**: `pip install --no-build-isolation -e .` (hapi imports numpy at setup time)
- **Venv**: `~/hapi-venv` (on home drive, not T7)
- **No test suite**: correctness validated by inspecting output artifacts
- **No CI/CD**: no GitHub Actions, no automated checks
- **No linter/formatter config**: standards enforced via docs and review
- **Artifact naming**: folder names encode run conditions (e.g. `_T600K_P3Torr_x0p008_L100cm_step0p001`)
- **Units**: wavenumber cm-1, pressure Torr, temperature K, path length cm throughout
- **PNG plots**: 4 subplots (dJ=-1, dJ=0, dJ=+1, overlay); branch colors: blue/orange/green

## ANTI-PATTERNS (THIS PROJECT)

- **NEVER modify `hapi/hapi.py`** for research work. It tracks upstream. New code goes in `research/` or `scripts/`.
- **NEVER suppress types** with `as any`, `@ts-ignore` (N/A for Python, but: no `# type: ignore` without justification)
- **NEVER commit secrets** (HITRAN credentials go in env vars: `HITRAN_EMAIL`/`HITRAN_PASSWORD`)
- **NEVER use `shell=True`** in subprocess calls
- **NEVER hardcode absolute paths** (use `research.io.default_paths()` or argparse)
- **NEVER commit gitignored data** (artifacts/, hitran_db/, exomol_db/). If they show in `git status`: `git rm -r --cached <dir>`
- **NEVER run above 2500 K**: TIPS2025 partition functions only cover 1-2500 K
- **NEVER treat `research/ch4_treanor.py` as a full polyatomic non-LTE model** — it is an effective single-mode approximation for the CH4 (0,0,v3,0) ladder only

## EXTERNAL DEPENDENCIES

numpy, matplotlib, plotly, scipy, pandas, pillow, setuptools. No requirements.txt — install manually or from setup.py.

## COMMANDS

```bash
# Install (editable)
pip install --no-build-isolation -e .

# Run a script
python scripts/plot_vibrational_mode_progressions.py
python scripts/plot_combined_exomol_i1_absorbance_progressions.py --help

# Validate Treanor → LTE identity (regression check for nonlte.py)
python scripts/validate_ch4_nu3_treanor_equilibrium.py

# Package (legacy)
python setup.py sdist && python setup.py bdist_wheel --universal
```

## GIT

- `origin`: fork at Ezeki3lRchy/hapi
- `upstream`: official hitranonline/hapi
- Changes to `hapi/hapi.py` should track upstream via merge/rebase

## KEY PHYSICS

- Cross-section: sigma(nu) = Sum_i S_i * phi_i(nu) (intensity * broadened profile)
- Absorbance: A = sigma * n * L (number density * path length)
- Voigt profile: convolution of Doppler (temperature) and Lorentz (pressure) broadening
- Line intensity is temperature-dependent; HITRAN reference is 296 K
- ExoMol intensity thresholds default to 0.0 (hot-band reference intensities ~1e-26 cm/mol)
- Treanor distribution (non-LTE): ln(N_v / N_0) = −v·θ_v/T_v + x_e·θ_v·v²/T_0 (Fridman 2008, Eq. 3-37)

## NOTES

- T7 external drive (exFAT): macOS creates `._*` resource fork files; these are gitignored
- `setup.py` imports `hapi.hapi.HAPI_VERSION` at install time (triggers module-level prints)
- `hapi/hapi.py` prints version/license text on import; suppress stdout if needed
- `research/__init__.py` uses lazy `__getattr__` loading; import errors surface on first submodule access
- Older scripts use module-level constants instead of CLI args; check top-of-file for config
- `classify/` contains legacy MATLAB + Python analysis; not part of the active workflow (see `classify/AGENTS.md`)

<!-- MANUAL: Custom project notes can be added below this line. The MANUAL marker tells the deepinit skill to preserve them on regeneration. -->
