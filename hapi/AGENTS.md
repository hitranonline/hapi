<!-- Parent: ../AGENTS.md -->
<!-- Generated: 2026-04-03 | Updated: 2026-04-22 -->

# hapi/ — Upstream HITRAN API Engine

## OVERVIEW

Single monolithic module (`hapi.py`, ~54k lines) tracking the upstream `hitranonline/hapi` repository. Provides the low-level spectroscopy engine: line data fetching, table management, absorption coefficient computation, and partition functions.

## DO NOT MODIFY THIS PACKAGE FOR RESEARCH WORK

All research code belongs in `research/`. If upstream changes are needed, coordinate with `hitranonline/hapi` or fork explicitly.

## STRUCTURE

```
hapi/
├── __init__.py    # `from .hapi import *` (wildcard re-export, no __all__)
└── hapi.py        # Everything: constants, I/O, queries, profiles, ISO tables
```

## KEY API FUNCTIONS

| Function | Purpose |
|----------|---------|
| `databaseBegin(db)` | Initialize local DB directory, load cache |
| `databaseCommit()` | Save in-memory tables to disk |
| `queryHITRAN(TableName, iso_ids, numin, numax)` | Fetch from remote HITRAN server |
| `storage2cache(TableName)` / `cache2storage(TableName)` | File <-> memory conversion |
| `absorptionCoefficient_Voigt()` | Compute absorption cross-sections |
| `transmittanceSpectrum()` / `absorptionSpectrum()` | Full spectrum generation |
| `partitionSum()` | Temperature-dependent partition functions (TIPS2025) |
| `select(TableName, Conditions=...)` | SQL-like table filtering |
| `createTable()` / `dropTable()` | Table lifecycle |
| `getColumn()` / `addColumn()` | Column access and computed columns |

## GOTCHAS

- **Import side effects**: prints version/license text on `import hapi`
- **Global mutable state**: `LOCAL_TABLE_CACHE` dict is the in-memory DB; `VARIABLES` dict holds config
- **Legacy naming**: tutorial text says `fetch()`, `db_begin()`, `tableList()` — actual names differ (see table above)
- **Condition language**: list/tuple AST, e.g. `('AND', ('>=', 'nu', 3000), ('<=', 'nu', 3100))`
- **XXX hack in `setRowObject`**: temporary insertion that increments header row count — upstream technical debt
- **`storage2cache` warning**: prints "WARNING: reading block of maximum N lines" for large tables
- **Format quirks**: handles Fortran D-exponent notation (e.g. "2.700-164" -> 2.700E-164)
- **TIPS2025 ceiling**: partition functions valid only for 1-2500 K

## CONVENTIONS

- All column data stored as masked numpy arrays (`np.ma.array`)
- Fixed-width 160-char HITRAN format (`HITRAN_DEFAULT_HEADER`)
- `CaseInsensitiveDict` used for all header/parameter dictionaries
- Procedural style throughout (free functions, no classes for domain logic)
