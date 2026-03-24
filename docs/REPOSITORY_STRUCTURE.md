# Repository Structure

This document explains what lives where after the first cleanup pass on the `codex/reorganize-repo` branch.

## Core source

- `hapi/`
  The importable HAPI package code.
- `research/`
  The new importable research workflow package built on top of HAPI.
- `scripts/`
  Standalone Python entry points for downloading data, extracting bands, building derived tables, plotting, and generating reports.
- `setup.py`
  Packaging entry point for the `hapi` package.

## Documentation

- `README.md`
  High-level repository guide.
- `docs/framework.md`
  Design spec for the future research workflow package.
- `docs/HITRAN_DATABASE_NOTES.md`
  Detailed notes about the local HITRAN/HAPI workflow.
- `docs/EXOMOL_DATABASE_NOTES.md`
  Detailed notes about the local ExoMol workflow.
- `CITATION.md`
  Citation text kept at the repo root for visibility.

## Local data and reference tables

- `hitran_db/`
  Local HITRAN-style `.header` and `.data` tables.
- `exomol_db/`
  Local ExoMol datasets.
- `plot_trans_all_species.ipynb`
  Legacy notebook still at the repo root.
- `[0]CH4,HITRAN_ ... .csv`
  Reference CSV kept at the repo root from earlier analysis work.

## Generated and legacy artifact directories

These directories are mostly outputs, intermediates, or workflow-specific exports. They remain in place for now to avoid a risky bulk move during the first pass.

- `ch4_nu3_progressions/`
- `co_nu1_progressions/`
- `species_plots_2500_3500_multiTP/`
- `methane_vibrational_bands/`
- `band_intensity_comparisons/`
- `exomol_ch4_mm_band_exports/`
- `exomol_ch4_mm_scans/`
- `exomol_ch4_mm_pure_nu3_band_texts/`
- `exomol_ch4_mm_pure_nu3_band_texts_hitran_style_2500_3500_sorted/`
- `exomol_test_band_texts/`
- `exomol_i1_smoke_band_texts/`

## Conventions for future cleanup

- New low-level spectroscopy code should go into `hapi/`.
- New repo-specific workflow code should go into `research/`.
- New one-shot or CLI-style tooling should go into `scripts/`.
- New documentation should go into `docs/`.
- New generated outputs should stay in workflow-specific artifact folders and be ignored in `.gitignore` when they are reproducible.

## Known remaining mess

- Several generated artifact directories still live at the repo root.
- Some script names are still workflow-specific rather than domain-grouped.
- The notebook and a few reference files still sit at the repo root.

Those are intentionally left for a later pass so this branch improves clarity without mixing structural cleanup with broad data migration.
