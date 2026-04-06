## Repository Structure Plan

This note captures a possible cleanup direction for the repo. It is intentionally a discussion draft, not a migration commit plan.

## Current Problems

1. Source code, runnable scripts, and generated outputs are mixed at the repo root.
2. Some workflow folders contain both code and artifacts, which makes ownership unclear.
3. Script naming is inconsistent across older exploratory work and newer reusable workflows.
4. Generated results are organized by historical task rather than by a stable top-level convention.

## Design Goals

1. Keep reusable code separate from one-shot runnable entry points.
2. Keep generated outputs separate from source code.
3. Preserve working historical workflows while making new work easier to locate.
4. Favor incremental migration over a risky one-shot reorganization.

## Recommended Direction

### 1. Keep the core code zones strict

- `hapi/` for upstream/core package code
- `research/` for reusable research functions and small dataclasses
- `scripts/` for runnable entry points only
- `docs/` for repository guidance and scientific notes

### 2. Standardize generated outputs

Use workflow-specific artifact folders for outputs. Long term, consider centralizing under a top-level `artifacts/` directory.

Example future direction:

```text
artifacts/
  treanor/
    figures/
    tables/
    runs/
  ch4_nu3/
  co_nu1/
  comparisons/
```

Short term, it is acceptable to keep existing workflow output folders if new code stops mixing scripts into them.

### 3. Make workflow folders either code-only or artifact-only

Avoid folders that contain both scripts and PNG/CSV outputs.

For example, `treanor_distribution/` should ideally be treated as an artifact folder, while runnable scripts live in `scripts/`.

### 4. Group script names by workflow and variant

Continue using lowercase, underscore-separated script names with a `plot_` prefix for plotting entry points.

Pattern:

```text
plot_<workflow>_<figure_or_variant>.py
```

Examples:

- `plot_treanor_fig3_3_ev.py`
- `plot_treanor_fig3_3_vib_number.py`
- `plot_treanor_two_panel_vib_number.py`

## Incremental Migration Strategy

### Phase 1: Stop adding new mixed-purpose folders

- Put new runnable code in `scripts/`
- Put new reusable logic in `research/`
- Keep new generated outputs in dedicated workflow artifact folders

### Phase 2: Clean one workflow at a time

- Choose one workflow such as Treanor
- Move runnable scripts into `scripts/`
- Leave outputs in the workflow artifact folder
- Update documentation and output paths

### Phase 3: Revisit historical root-level artifact folders

- Migrate only when there is a clear benefit
- Avoid breaking older workflows without a replacement path

## Open Questions

1. Should the repo eventually adopt a single top-level `artifacts/` directory?
2. Should `scripts/` remain flat, or be grouped by workflow in subdirectories?
3. Which historical output folders are still actively used and should remain where they are?
4. Which workflows deserve promotion from exploratory script form into reusable `research/` APIs?

## Recommended Next Discussion Topics

1. Decide whether `artifacts/` should be introduced now or later.
2. Decide whether `scripts/` should stay flat or become grouped by workflow.
3. Identify the first 2-3 historical artifact folders that are safe to normalize.
4. Define a naming policy for figures, tables, and run directories.
