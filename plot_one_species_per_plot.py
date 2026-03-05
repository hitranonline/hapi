import glob
import os
import re

import hapi

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise RuntimeError("matplotlib is required. Install with: pip install matplotlib") from exc

# -------------------------
# USER SETTINGS
# -------------------------
DB_DIR = "hitran_db"
WN_MIN = 2500.0
WN_MAX = 3500.0
WN_STEP = 0.05  # faster default for per-species plotting
TEMPERATURE_K = 296.0
PRESSURE_ATM = 1.0
OUTPUT_DIR = "species_plots_2500_3500"

# -------------------------
# 1) Open local HAPI database
# -------------------------
hapi.db_begin(DB_DIR)

# -------------------------
# 2) Discover downloaded hydrocarbon table names
# -------------------------
header_files = glob.glob(os.path.join(DB_DIR, "*.header"))
name_pat = re.compile(r"^(C\d*H\d*)_M\d+_I1\.header$", re.IGNORECASE)

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

os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Found {len(source_tables)} species tables")

# -------------------------
# 3) One plot per species
# -------------------------
for table in source_tables:
    print(f"Computing: {table}")
    try:
        nu, coef = hapi.absorptionCoefficient_Voigt(
            SourceTables=[table],
            WavenumberRange=[WN_MIN, WN_MAX],
            WavenumberStep=WN_STEP,
            Environment={"T": TEMPERATURE_K, "p": PRESSURE_ATM},
            HITRAN_units=False,
        )

        out_png = os.path.join(OUTPUT_DIR, f"{table}_{int(WN_MIN)}_{int(WN_MAX)}.png")
        plt.figure(figsize=(10, 5))
        plt.plot(nu, coef, lw=1.0)
        plt.xlim(WN_MIN, WN_MAX)
        plt.xlabel("Wavenumber (cm^-1)")
        plt.ylabel("Absorption coefficient")
        plt.title(f"{table} (Voigt) {int(WN_MIN)}-{int(WN_MAX)} cm^-1")
        plt.grid(alpha=0.25)
        plt.tight_layout()
        plt.savefig(out_png, dpi=300)
        plt.close()

        print(f"Saved: {out_png}")
    except Exception as exc:
        print(f"FAILED {table}: {exc}")

print("Done.")
