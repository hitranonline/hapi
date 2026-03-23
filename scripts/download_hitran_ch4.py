import re
from pathlib import Path
import pandas as pd

from _bootstrap import ensure_repo_root

ensure_repo_root()

import hapi  # classic HAPI module

# -------------------------
# USER SETTINGS
# -------------------------
# For this classic hapi package, fetch() does not consume an API key directly.
# Keep this variable only for your own reference if needed.
API_KEY = "86677861-553f-4190-b391-950f9a40c94c"
ROOT_DIR = Path(__file__).resolve().parents[1]
DB_DIR = str(ROOT_DIR / "hitran_db")
NU_MIN = 0.0      # cm^-1 Lower wavenumber bound (cm^-1).
NU_MAX = 10000.0  # cm^-1 Upper wavenumber bound (cm^-1).
I = 1             # main isotopologue

# -------------------------
# 1) Init local HAPI database
# -------------------------
hapi.db_begin(DB_DIR) # Initializes/opens local HAPI database directory.

# -------------------------
# 2) Load HITRAN molecule metadata (id + formula)
# -------------------------
url = "https://hitran.org/docs/molec-meta/" # Source page containing molecule metadata table.
mol_df = pd.read_html(url)[0] # Reads first HTML table from that page into a dataframe.
mol_df.columns = [str(c).strip().lower().replace(" ", "_") for c in mol_df.columns] 
# Makes column names lowercase and replaces spaces with _ to reduce schema issues.

# Standardize key column names
if "molecule_id" not in mol_df.columns and "id" in mol_df.columns:
    mol_df = mol_df.rename(columns={"id": "molecule_id"})

if "formula" not in mol_df.columns:
    for c in mol_df.columns:
        if "formula" in c:
            mol_df = mol_df.rename(columns={c: "formula"})
            break

# Hard fail early if table structure changed
missing = [c for c in ("molecule_id", "formula") if c not in mol_df.columns]
if missing:
    raise RuntimeError(
        f"Missing required columns from {url}: {missing}. Found: {list(mol_df.columns)}"
    )

# -------------------------
# 3) Filter formulas that contain ONLY C and H
# -------------------------
# Accept: CH4, C2H2, C2H6, C10H8, CH3CH3
# Reject: CO2, CH3OH, HCN, C2H4O
pat = re.compile(r"^C\d*H\d*$")
hc = mol_df[mol_df["formula"].astype(str).str.match(pat, na=False)].copy()

print("Hydrocarbon-only formulas found:", len(hc))
print(hc[["molecule_id", "formula"]].head(20).to_string(index=False))

# -------------------------
# 4) Download line lists
# -------------------------
for _, row in hc.iterrows():
    M = int(row["molecule_id"])
    formula = str(row["formula"])
    table_name = f"{formula}_M{M}_I{I}"

    print(f"Fetching {table_name}: M={M}, I={I}, nu=[{NU_MIN},{NU_MAX}] cm^-1")
    try:
        hapi.fetch(table_name, M, I, NU_MIN, NU_MAX)
    except Exception as exc:
        print(f"FAILED {table_name}: {exc}")

print("Done.")
