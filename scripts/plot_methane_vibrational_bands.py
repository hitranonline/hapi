import csv
import os
import re
from collections import defaultdict
from pathlib import Path

from _bootstrap import ensure_repo_root

ensure_repo_root()

import hapi
import numpy as np
import plotly.graph_objects as go

# -------------------------
# USER SETTINGS
# -------------------------
ROOT_DIR = Path(__file__).resolve().parents[1]
DB_DIR = str(ROOT_DIR / "hitran_db")
OUTPUT_DIR = str(ROOT_DIR / "methane_vibrational_bands")

SOURCE_TABLE = "CH4_M6_I1"

WN_MIN = 2500.0
WN_MAX = 3500.0
WN_STEP = 0.001

T_K = 600.0
P_TORR = 3.0
P_ATM = P_TORR / 760.0

MOLE_FRACTION = 0.008
PATH_LENGTH_CM = 100.0
LINE_INTENSITY_THRESHOLD = 1.0e-23

# Use None to plot every detected band.
MAX_BANDS_TO_PLOT = 12


def clean_quanta_label(value: str) -> str:
    text = value.strip()
    return text if text else "000"


def safe_label_fragment(value: str) -> str:
    text = clean_quanta_label(value)
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_") or "000"


def main() -> None:
    hapi.db_begin(DB_DIR)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    header_path = os.path.join(DB_DIR, f"{SOURCE_TABLE}.header")
    if not os.path.exists(header_path):
        raise FileNotFoundError(f"Missing table header: {header_path}")

    m = re.search(r"_M(\d+)_I(\d+)", SOURCE_TABLE)
    if m is None:
        raise ValueError(f"Could not parse molecule/isotopologue IDs from {SOURCE_TABLE}")

    mol_id = int(m.group(1))
    iso_id = int(m.group(2))

    global_lower = hapi.getColumn(SOURCE_TABLE, "global_lower_quanta")
    global_upper = hapi.getColumn(SOURCE_TABLE, "global_upper_quanta")
    line_strengths = hapi.getColumn(SOURCE_TABLE, "sw")

    band_stats = defaultdict(lambda: {"sum_sw": 0.0, "count": 0})
    for lo_raw, up_raw, sw in zip(global_lower, global_upper, line_strengths):
        key = (lo_raw, up_raw)
        band_stats[key]["sum_sw"] += float(sw)
        band_stats[key]["count"] += 1

    sorted_bands = sorted(
        band_stats.items(),
        key=lambda item: item[1]["sum_sw"],
        reverse=True,
    )

    if MAX_BANDS_TO_PLOT is not None:
        sorted_bands = sorted_bands[:MAX_BANDS_TO_PLOT]

    fig = go.Figure()
    any_curve_plotted = False
    csv_wavenumbers = None
    csv_columns = []

    for band_index, ((lo_raw, up_raw), stats) in enumerate(sorted_bands, start=1):
        lo_label = clean_quanta_label(lo_raw)
        up_label = clean_quanta_label(up_raw)
        band_label = f"{lo_label} -> {up_label}"
        temp_table = (
            f"__CH4_BAND_{band_index:03d}_"
            f"{safe_label_fragment(lo_raw)}_TO_{safe_label_fragment(up_raw)}__"
        )

        try:
            hapi.select(
                SOURCE_TABLE,
                DestinationTableName=temp_table,
                Conditions=(
                    "and",
                    ("=", "global_lower_quanta", lo_raw),
                    ("=", "global_upper_quanta", up_raw),
                ),
                Output=False,
            )

            if hapi.length(temp_table) == 0:
                continue

            nu, coef = hapi.absorptionCoefficient_Voigt(
                Components=[(mol_id, iso_id, MOLE_FRACTION)],
                SourceTables=[temp_table],
                WavenumberRange=[WN_MIN, WN_MAX],
                WavenumberStep=WN_STEP,
                Environment={"T": T_K, "p": P_ATM},
                Diluent={"self": MOLE_FRACTION, "He": 1 - MOLE_FRACTION},
                IntensityThreshold=LINE_INTENSITY_THRESHOLD,
                HITRAN_units=False,
            )

            _, trans = hapi.transmittanceSpectrum(
                nu,
                coef,
                Environment={"l": PATH_LENGTH_CM},
            )
            absorbance = -np.log(trans)

            if csv_wavenumbers is None:
                csv_wavenumbers = nu
            elif len(nu) != len(csv_wavenumbers) or any(a != b for a, b in zip(nu, csv_wavenumbers)):
                raise RuntimeError(
                    f"Wavenumber grid mismatch for band {band_label}; CSV export requires a shared grid."
                )

            column_name = f"band_{band_index:02d}_{safe_label_fragment(lo_raw)}_to_{safe_label_fragment(up_raw)}"
            csv_columns.append((column_name, absorbance))

            fig.add_trace(
                go.Scatter(
                    x=nu,
                    y=absorbance,
                    mode="lines",
                    name=f"{band_label} (N={stats['count']}, sum S={stats['sum_sw']:.2e})",
                )
            )
            any_curve_plotted = True
        except Exception as exc:
            print(f"FAILED {band_label}: {exc}")
        finally:
            hapi.dropTable(temp_table)

    if not any_curve_plotted:
        raise RuntimeError("No methane vibrational bands were plotted.")

    fig.update_layout(
        title=(
            f"CH4 Vibrational-Band Absorbance (Voigt) "
            f"{int(WN_MIN)}-{int(WN_MAX)} cm^-1 | T={T_K:g} K, P={P_TORR:g} Torr"
            f" | x={MOLE_FRACTION:g}, L={PATH_LENGTH_CM:g} cm, S>{LINE_INTENSITY_THRESHOLD:.0e}"
        ),
        xaxis_title="Wavenumber (cm^-1)",
        yaxis_title="Absorbance",
        template="plotly_white",
    )
    fig.update_xaxes(range=[WN_MIN, WN_MAX])

    out_stem = (
        f"CH4_vibrational_bands_abs_{int(WN_MIN)}_{int(WN_MAX)}"
        f"_T{int(T_K)}K_P{P_TORR:g}Torr_L{int(PATH_LENGTH_CM)}cm"
    )
    out_html = os.path.join(OUTPUT_DIR, f"{out_stem}.html")
    out_csv = os.path.join(OUTPUT_DIR, f"{out_stem}.csv")

    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["wavenumber_cm-1"] + [name for name, _ in csv_columns])
        for i, wn in enumerate(csv_wavenumbers):
            writer.writerow([wn] + [values[i] for _, values in csv_columns])

    fig.write_html(out_html, include_plotlyjs="cdn")
    print(f"Saved: {out_csv}")
    print(f"Saved: {out_html}")


if __name__ == "__main__":
    main()
