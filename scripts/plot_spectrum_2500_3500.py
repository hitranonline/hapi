import glob
import os
import re
from pathlib import Path

from _bootstrap import ensure_repo_root

ensure_repo_root()

import hapi

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise RuntimeError("matplotlib is required. Install with: pip install matplotlib") from exc

# -------------------------
# USER SETTINGS
# -------------------------
ROOT_DIR = Path(__file__).resolve().parents[1]
DB_DIR = str(ROOT_DIR / "hitran_db")
WN_MIN = 2500.0
WN_MAX = 3500.0
WN_STEP = 0.01
TEMPERATURE_K = 296.0
PRESSURE_ATM = 1.0
OUTPUT_PNG = str(ROOT_DIR / "spectrum_2500_3500.png")

# -------------------------
# 1) Open local HAPI database
# -------------------------
hapi.db_begin(DB_DIR)

# -------------------------
# 2) Discover downloaded hydrocarbon table names
# -------------------------
header_files = glob.glob(os.path.join(DB_DIR, "*.header"))

# Match names like CH4_M6_I1.header, C2H2_M26_I1.header, CH3_M57_I1.header
name_pat = re.compile(r"^(?P<formula>C\d*H\d*)_M\d+_I1\.header$", re.IGNORECASE)

source_tables = []
for path in header_files:
    fname = os.path.basename(path)
    if name_pat.match(fname):
        source_tables.append(os.path.splitext(fname)[0])

source_tables = sorted(set(source_tables))

if not source_tables:
    raise RuntimeError(
        f"No matching hydrocarbon tables found in '{DB_DIR}'. "
        "Expected files like CH4_M6_I1.header"
    )

print(f"Using {len(source_tables)} source tables:")
for t in source_tables:
    print(f"  - {t}")

# -------------------------
# 3) Calculate absorption coefficient (Voigt)
# -------------------------
nu, coef = hapi.absorptionCoefficient_Voigt(
    SourceTables=source_tables,
    WavenumberRange=[WN_MIN, WN_MAX],
    WavenumberStep=WN_STEP,
    Environment={"T": TEMPERATURE_K, "p": PRESSURE_ATM},
    HITRAN_units=False,
)

# -------------------------
# 4) Plot and save
# -------------------------
plt.figure(figsize=(10, 5))
plt.plot(nu, coef, lw=1.0)
plt.xlim(WN_MIN, WN_MAX)
plt.xlabel("Wavenumber (cm^-1)")
plt.ylabel("Absorption coefficient")
plt.title("Hydrocarbon Spectrum (Voigt) 2500-3500 cm^-1")
plt.grid(alpha=0.25)
plt.tight_layout()
plt.savefig(OUTPUT_PNG, dpi=160)

print(f"Saved plot: {OUTPUT_PNG}")
