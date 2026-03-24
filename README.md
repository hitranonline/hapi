# hapi Workspace

This repository currently serves two roles:

- the importable `hapi/` Python package
- a spectroscopy workspace with local databases, analysis scripts, and generated outputs

The recent cleanup starts separating those concerns so the repo is easier to read.

## Current layout

- `hapi/`: the library/package code
- `research/`: the new repo-specific research workflow package
- `scripts/`: standalone workflow utilities and analysis entry points
- `docs/`: repository notes and structure documentation
- `hitran_db/`: local HITRAN-style tables used by HAPI and the analysis scripts
- `exomol_db/`: local ExoMol source files
- `*_progressions/`, `band_intensity_comparisons/`, `exomol_*`: generated workflow outputs and intermediate artifacts

## Where to start

- Read [docs/REPOSITORY_STRUCTURE.md](docs/REPOSITORY_STRUCTURE.md) for the top-level map.
- Read [docs/framework.md](docs/framework.md) for the planned Python package direction.
- Read [scripts/README.md](scripts/README.md) for the script inventory and workflow entry points.
- Read [docs/HITRAN_DATABASE_NOTES.md](docs/HITRAN_DATABASE_NOTES.md) or [docs/EXOMOL_DATABASE_NOTES.md](docs/EXOMOL_DATABASE_NOTES.md) for the spectroscopy-specific background.

## Working convention going forward

- Put low-level spectroscopy code in `hapi/`.
- Put repo-specific research workflows in `research/`.
- Put runnable utilities in `scripts/`.
- Put explanatory Markdown in `docs/`.
- Treat generated outputs as artifacts, not as primary source code.

This branch does not try to relocate every historical artifact directory yet. It establishes a clearer source/docs layout first, then documents the remaining legacy directories so the next cleanup pass can be more deliberate.
