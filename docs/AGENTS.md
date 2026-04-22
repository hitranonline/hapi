<!-- Parent: ../AGENTS.md -->
<!-- Generated: 2026-04-22 | Updated: 2026-04-22 -->

# docs/ — Design Specs, Physics Notes, Workflow Guides

## Purpose

Human-readable markdown documentation for the research workspace. Three kinds of content live here: (1) authoritative design specs (`framework.md`, `code_review_checklist.md`) that define how this repo is built, (2) physics/domain notes explaining the spectroscopy and distribution models, and (3) workflow narratives explaining how specific scripts work.

Read `docs/framework.md` before adding code to `research/`. Read `docs/code_review_checklist.md` before reviewing AI-generated code.

## Key Files

| File | Kind | Description |
|------|------|-------------|
| `framework.md` | Design spec | Python-standard package design rules for `research/` (functions + dataclasses, no orchestration classes) |
| `code_review_checklist.md` | Design spec | Scientific-computing code review checklist — AI validation, correctness, security, reproducibility |
| `REPOSITORY_STRUCTURE.md` | Layout reference | What lives where after the `codex/reorganize-repo` cleanup (core source, docs, data, artifacts) |
| `HITRAN_DATABASE_NOTES.md` | Physics | How HITRAN line tables are structured and consumed by `hapi` / `research/hitran.py` |
| `EXOMOL_DATABASE_NOTES.md` | Physics | ExoMol `.states` + `.trans` + `.pf` + `.def` decomposition; how `research/exomol.py` reconstructs line quantities |
| `NU3_PROGRESSIONS_WORKFLOW.md` | Workflow narrative | The four nu3 absorbance scripts (exomol-sorted, hitran-band-text, combined-pure-nu3, combined-exomol-i1) — inputs, grouping, labels |
| `fundamental.md` | Physics | Line intensity vs. final spectrum distinction; broadening → sum → absorbance pipeline |
| `TREANOR_DISTRIBUTION_NOTES.md` | Physics | Two-temperature Treanor distribution for vibrational non-LTE; maps to `research/treanor.py` |
| `hot_band_spectral_shape_and_magnitude.md` | Physics | Why nu3 fundamental is comb-shaped but hot bands fill in (symmetry, line density, curation) |
| `hot_band_temperature_dependence.md` | Physics | Boltzmann population math for hot bands; TIPS 2500 K ceiling explained |
| `explicit_j_6_6_tol_1e4.md` | Report | Frozen run-log of a specific HITRAN vs. ExoMol nearest-line comparison (J=(6,6), tol=1e-4) |

## Subdirectories

None.

## For AI Agents

### Working In This Directory

- These are **prose docs**, not code. No tests, no linting. Treat them as authoritative references.
- When `research/` or `scripts/` behavior changes, check here for stale claims (file paths, function names, default values) and update them.
- Inline file links in physics notes (`HITRAN_DATABASE_NOTES.md`, `EXOMOL_DATABASE_NOTES.md`) predate the `scripts/` reorganization — some point to old root-level script locations and have not been fully migrated.
- Math is written in LaTeX inside `$...$` / `$$...$$`. Do not rewrite unless you are correcting a physical error.
- When adding a new physics note, link to it from root `AGENTS.md` → "WHERE TO LOOK" and from the relevant code module's docstring.

### Testing Requirements

None. Prose-only.

### Common Patterns

- Each note starts with a one-paragraph framing, then develops mathematical or structural detail, then links to specific code files that implement the ideas.
- Tables are preferred over bullet lists for parameter sweeps and branch summaries.
- Uppercase section headers in code-oriented docs (`# OVERVIEW`, `# STRUCTURE`, `# CONVENTIONS`) match the style of the AGENTS.md files and let agents scan both sets uniformly.

## Dependencies

### Internal

- Most physics notes link back to specific modules in `research/` and scripts in `scripts/`. Keep those links current.
- `framework.md` is the canonical design spec referenced by `research/AGENTS.md`.

### External

- None — pure markdown.

<!-- MANUAL: Custom notes about docs workflow can be added below this line -->
