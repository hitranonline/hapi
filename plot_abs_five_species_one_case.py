import os
import re

import hapi
import plotly.graph_objects as go

# -------------------------
# USER SETTINGS
# -------------------------
DB_DIR = "hitran_db"
OUTPUT_DIR = "5_species"

WN_MIN = 2500.0
WN_MAX = 3500.0
WN_STEP = 0.001

T_K = 600.0
P_TORR = 4.0
P_ATM = P_TORR / 760.0
PATH_LENGTH_CM = 100.0

SPECIES = [
    {"name": "Methane", "formula": "CH4", "table": "CH4_M6_I1", "x": 0.10},
    {"name": "Methyl", "formula": "CH3", "table": "CH3_M57_I1", "x": 1000e-6},
    {"name": "Ethane", "formula": "C2H6", "table": "C2H6_M27_I1", "x": 0.05},
    {"name": "Acetylene", "formula": "C2H2", "table": "C2H2_M26_I1", "x": 0.01},
    {"name": "Ethylene", "formula": "C2H4", "table": "C2H4_M38_I1", "x": 0.01},
]


def main() -> None:
    hapi.db_begin(DB_DIR)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    fig = go.Figure()
    any_curve_plotted = False

    for sp in SPECIES:
        table = sp["table"]
        
        m = re.search(r"_M(\d+)_I(\d+)", table)
        mol_id = int(m.group(1)) # type: ignore
        iso_id = int(m.group(2)) # type: ignore

        x_i = sp["x"]  # mole fraction
        # p_i_atm = x_i * P_ATM  # partial pressure

        header_path = os.path.join(DB_DIR, f"{table}.header")
        if not os.path.exists(header_path):
            print(f"Missing table header: {header_path}")
            continue

        try:
            nu, coef = hapi.absorptionCoefficient_Voigt(
                Components=[(mol_id, iso_id, x_i)],
                SourceTables=[table],
                WavenumberRange=[WN_MIN, WN_MAX],
                WavenumberStep=WN_STEP,
                Environment={"T": T_K, "p": P_ATM},
                Diluent={"self": x_i, "He": 1 - x_i},
                HITRAN_units=False,
            )

            _, trans = hapi.transmittanceSpectrum(
                nu,
                coef,
                Environment={"l": PATH_LENGTH_CM},
            )
            absorptance = 1.0 - trans

            label = f"{sp['name']} ({sp['formula']}, x={x_i:g})"
            fig.add_trace(
                go.Scatter(
                    x=nu,
                    y=absorptance,
                    mode="lines",
                    name=label,
                )
            )
            any_curve_plotted = True

        except Exception as exc:
            print(f"FAILED {table}: {exc}")

    if not any_curve_plotted:
        raise RuntimeError("No curves were plotted.")

    fig.update_layout(
        title=(
            f"5-Species Absorptance (Voigt) "
            f"{int(WN_MIN)}-{int(WN_MAX)} cm^-1 | T={T_K:g} K, P={P_TORR:g} Torr"
        ),
        xaxis_title="Wavenumber (cm^-1)",
        yaxis_title=f"Absorptance (path length = {PATH_LENGTH_CM:g} cm)",
        template="plotly_white",
    )
    fig.update_xaxes(range=[WN_MIN, WN_MAX])

    out_html = os.path.join(
        OUTPUT_DIR,
        (
            f"five_species_abs_{int(WN_MIN)}_{int(WN_MAX)}"
            f"_T{int(T_K)}K_P{P_TORR:g}Torr_L{int(PATH_LENGTH_CM)}cm.html"
        ),
    )
    fig.write_html(out_html, include_plotlyjs="cdn")
    print(f"Saved: {out_html}")


if __name__ == "__main__":
    main()
