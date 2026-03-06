import glob
import os
import re

import hapi
import plotly.graph_objects as go

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
WN_STEP = 0.01
PATH_LENGTH_CM = 100.0 # 100 cm
OUTPUT_DIR = "species_plots_2500_3500_multiTP"

# Define the (T, p) cases you want on the same plot
CASES = [
    {"T": 296.0, "p": 1.0,  "label": "T=296K, p=1atm"},
    {"T": 350.0, "p": 1.0,  "label": "T=350K, p=1atm"},
    {"T": 400.0, "p": 1.0, "label": "T=400K, p=1atm"},
]

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
# 3) One plot per species, multiple (T,p) curves on same axes
# -------------------------
for table in source_tables:
    print(f"Computing: {table}")

    # plt.figure(figsize=(10, 5))
    fig = go.Figure()
    any_curve_plotted = False

    for case in CASES:
        T_K = case["T"]
        P_ATM = case["p"]
        label = case["label"]

        try:
            nu, coef = hapi.absorptionCoefficient_Voigt(
                SourceTables=[table],
                WavenumberRange=[WN_MIN, WN_MAX],
                WavenumberStep=WN_STEP,
                Environment={"T": T_K, "p": P_ATM},
                HITRAN_units=False,
            )

            # coef must be in cm^-1 for transmittanceSpectrum
            _, trans = hapi.transmittanceSpectrum(
                nu,
                coef,
                Environment={"l": PATH_LENGTH_CM}  # cm
            )
            absorptance = 1.0 - trans

            # plt.plot(nu, trans, lw=1.0, label=label)
            fig.add_trace(
                go.Scatter(
                    x=nu, 
                    y=absorptance, 
                    mode='lines', 
                    name=label
                )
            )
            any_curve_plotted = True

        except Exception as exc:
            print(f"  FAILED {table} at {label}: {exc}")

    if not any_curve_plotted:
        plt.close()
        continue

    # out_png = os.path.join(OUTPUT_DIR, f"{table}_{int(WN_MIN)}_{int(WN_MAX)}.png")
    # plt.xlim(WN_MIN, WN_MAX)
    # plt.xlabel("Wavenumber (cm^-1)")
    # plt.ylabel(f"Absorptance (path length = {PATH_LENGTH_CM} cm)")
    # plt.title(f"{table} (Voigt) {int(WN_MIN)}-{int(WN_MAX)} cm^-1")
    # plt.grid(alpha=0.25)
    # plt.legend(loc="best", frameon=True)
    # plt.tight_layout()
    # plt.savefig(out_png, dpi=1000)
    # plt.close()

    # print(f"Saved: {out_png}")
    fig.update_layout(
    title=f"{table} (Voigt) {int(WN_MIN)}-{int(WN_MAX)} cm^-1",
    xaxis_title="Wavenumber (cm^-1)",
    yaxis_title=f"Absorptance (path length = {PATH_LENGTH_CM} cm)",
    template="plotly_white",
    )
    fig.update_xaxes(range=[WN_MIN, WN_MAX])

    out_html = os.path.join(OUTPUT_DIR, f"{table}_{int(WN_MIN)}_{int(WN_MAX)}.html")
    fig.write_html(out_html, include_plotlyjs="cdn")
    print(f"Saved: {out_html}")

print("Done.")