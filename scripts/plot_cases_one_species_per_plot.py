import glob
import os
import re
from pathlib import Path
import numpy as np

from _bootstrap import ensure_repo_root

ensure_repo_root()

import hapi
import plotly.graph_objects as go

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

PATH_LENGTH_CM = 100.0 # 100 cm

OUTPUT_DIR = str(ROOT_DIR / "species_plots_2500_3500_test")
SELECT_SPECIES = ["CH3", "CH4"]  # e.g. ["CH3"] or ["CH3_M57_I1", "C2H2"]
CSV_OVERLAY_FILES = [
    str(ROOT_DIR / "[0]CH4,HITRAN_ X = 0.0001, T = 300 K, P = 0.13157894736842105 atm, L = 100 cm .csv")
]

X = 1e-4 # 100 ppm, x = 1e-4

P_ATM_TOTAL = 100 / 760 # 100 Torr in atm



# Define the (T, p) cases you want on the same plot
# CASES = [
#     {"T": 296.0, "p": 1.0,  "label": "T=296K, p=1atm"},
#     {"T": 350.0, "p": 1.0,  "label": "T=350K, p=1atm"},
#     {"T": 400.0, "p": 1.0, "label": "T=400K, p=1atm"},
# ]
CASES = [
    {"T": 300.0, "p": P_ATM_TOTAL,  "label": f"T=300K, p={P_ATM_TOTAL}atm"},
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

if SELECT_SPECIES:
    selected = {s.strip().upper() for s in SELECT_SPECIES if s.strip()}
    source_tables = [
        table for table in source_tables
        if table.upper() in selected or table.split("_")[0].upper() in selected
    ]
    if not source_tables:
        raise RuntimeError(
            f"No tables matched SELECT_SPECIES={SELECT_SPECIES}. "
            "Use formula (e.g. CH3) or full table name (e.g. CH3_M57_I1)."
        )

os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Found {len(source_tables)} species tables")

# Build species -> CSV overlays map from names like "[0]CH4,...csv"
overlay_curves_by_species = {}
overlay_species_pat = re.compile(r"\[\d+\](C\d*H\d+)", re.IGNORECASE)
for csv_path in CSV_OVERLAY_FILES:
    if not os.path.isfile(csv_path):
        print(f"CSV overlay not found, skipping: {csv_path}")
        continue

    base_name = os.path.basename(csv_path)
    species_match = overlay_species_pat.search(base_name)
    if not species_match:
        print(f"Could not infer species from CSV name, skipping: {csv_path}")
        continue

    species = species_match.group(1).upper()
    try:
        arr = np.loadtxt(csv_path, delimiter=",", skiprows=1)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        if arr.shape[1] < 2:
            raise ValueError("CSV must have at least two columns (nu, y)")

        overlay_curves_by_species.setdefault(species, []).append(
            {
                "x": arr[:, 0],
                "y": arr[:, 1],
                "label": f"CSV: {base_name}",
            }
        )
        print(f"Loaded CSV overlay for {species}: {base_name}")
    except Exception as exc:
        print(f"Failed to load CSV overlay '{csv_path}': {exc}")

# -------------------------
# 3) One plot per species, multiple (T,p) curves on same axes
# -------------------------
for table in source_tables:
    print(f"Computing: {table}")

    # plt.figure(figsize=(10, 5))
    fig = go.Figure()
    any_curve_plotted = False
    m = re.search(r"_M(\d+)_I(\d+)", table)
    mol_id = int(m.group(1))
    iso_id = int(m.group(2))    
    species_formula = table.split("_")[0].upper()

    for case in CASES:
        T_K_I = case["T"]
        P_ATM_I = case["p"]
        label = case["label"]

        try:
            nu, coef = hapi.absorptionCoefficient_Voigt(
                Components=[(mol_id, iso_id, X)],
                SourceTables=[table],
                WavenumberRange=[WN_MIN, WN_MAX],
                WavenumberStep=WN_STEP,
                Environment={"T": T_K_I, "p": P_ATM_I},
                Diluent={"self": X, "air": 1 - X},
                HITRAN_units=False,
            )

            # coef must be in cm^-1 for transmittanceSpectrum
            _, trans = hapi.transmittanceSpectrum(
                nu,
                coef,
                Environment={"l": PATH_LENGTH_CM}  # cm
            )
            # absorptance = 1.0 - trans
            # absorbance = - np.log(trans)
            absorbance = - np.log(trans)
            # plt.plot(nu, trans, lw=1.0, label=label)
            fig.add_trace(
                go.Scatter(
                    x=nu, 
                    y=absorbance, 
                    mode='lines', 
                    name=label
                )
            )
            fig.update_layout(showlegend=True)
            any_curve_plotted = True

        except Exception as exc:
            print(f"  FAILED {table} at {label}: {exc}")

    if not any_curve_plotted:
        plt.close()
        continue

    # Add optional external CSV overlays for this species (if any)
    for overlay in overlay_curves_by_species.get(species_formula, []):
        fig.add_trace(
            go.Scatter(
                x=overlay["x"],
                y=overlay["y"],
                mode="lines",
                name=overlay["label"],
                line={"dash": "dash"},
            )
        )
        any_curve_plotted = True

    # print(f"Saved: {out_png}")
    fig.update_layout(
    title=f"{table} (Voigt) {int(WN_MIN)}-{int(WN_MAX)} cm^-1",
    xaxis_title="Wavenumber [cm^-1]",
    yaxis_title=f"Absorbance [-] (path length = {PATH_LENGTH_CM} cm)",
    template="plotly_white",
    )
    fig.update_xaxes(range=[WN_MIN, WN_MAX])

    out_html = os.path.join(OUTPUT_DIR, f"{table}_{int(WN_MIN)}_{int(WN_MAX)}.html")
    fig.write_html(out_html, include_plotlyjs="cdn")
    print(f"Saved: {out_html}")

print("Done.")
