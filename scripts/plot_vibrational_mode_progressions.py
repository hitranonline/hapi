"""
How to use this script
----------------------
1. Choose `SPECIES_KEY` such as `"CH4"` or `"CO"`.
2. Normally you do not need to set a table name manually. The default table for that species is used automatically.
3. Only use `SOURCE_TABLE_OVERRIDE` if you want a non-default table.
4. Adjust the spectral range and gas conditions below.
5. Run:
       python scripts/plot_vibrational_mode_progressions.py

Outputs
-------
- Exact-band outputs: one CSV/HTML per exact lower/upper vibrational-state pair
- Category outputs: merged CSV/HTML per category such as `fundamental_band`

How classification works
------------------------
- The script reads `global_lower_quanta` and `global_upper_quanta`
- It extracts the vibrational quantum numbers for the chosen species
- It tracks only the selected mode:
  - CH4 uses `nu3` -> mode index 2 in `(v1, v2, v3, v4)`
  - CO uses `nu1` -> mode index 0 in `(v,)`
- Category rules:
  - `0->1` => `fundamental_band`
  - `0->2` => `first_overtone_band`
  - `0->3` => `second_overtone_band`
  - `1->2`, `2->3`, ... => `hot_fundamental_band`
  - `1->3`, `2->4`, ... => `hot_first_overtone_band`

Important interpretation
------------------------
- `REQUIRE_SAME_OTHER_MODES = True` means all other vibrational modes must stay unchanged
- `REQUIRE_UNIT_STEP = False` means overtone-type jumps like `0->2` are allowed
"""

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

SPECIES_KEY = "CO"
SOURCE_TABLE_OVERRIDE = ""

WN_MIN = 2500.0
WN_MAX = 3500.0
WN_STEP = 0.001

T_K = 600.0
P_TORR = 3.0
P_ATM = P_TORR / 760.0

MOLE_FRACTION = 0.008
PATH_LENGTH_CM = 100.0
LINE_INTENSITY_THRESHOLD = 1.0e-23

REQUIRE_SAME_OTHER_MODES = True
REQUIRE_UNIT_STEP = False
MIN_LINES_PER_GROUP = 1

GENERATE_EXACT_BAND_OUTPUTS = True
GENERATE_CATEGORY_OUTPUTS = True


def clean_quanta_label(value: str) -> str:
    text = value.strip()
    return text if text else "000"


def safe_label_fragment(value: str) -> str:
    text = clean_quanta_label(value)
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_") or "000"


def parse_ch4_global_quanta(value: str) -> tuple[int, ...]:
    parts = value.split()
    if len(parts) != 5:
        raise ValueError(f"Unexpected CH4 global quanta label: {value!r}")
    return tuple(int(parts[i]) for i in range(4))


def extract_vibrational_quantum(value: object) -> int | None:
    text = str(value).strip()
    if not text:
        return None
    if re.fullmatch(r"[+-]?\d+", text):
        return int(text)

    labeled_patterns = (
        r"(?i)\bv(?:p|pp)?\s*(?:=|:)\s*([+-]?\d+)\b",
        r"(?i)\bv['\"]{0,2}\s*(?:=|:)\s*([+-]?\d+)\b",
    )
    for pattern in labeled_patterns:
        match = re.search(pattern, text)
        if match is not None:
            return int(match.group(1))

    numbers = re.findall(r"[+-]?\d+", text)
    if len(numbers) == 1:
        return int(numbers[0])
    return None


def parse_co_global_quanta(value: str) -> tuple[int, ...]:
    quantum = extract_vibrational_quantum(value)
    if quantum is None:
        raise ValueError(f"Unexpected CO global quanta label: {value!r}")
    return (quantum,)


SPECIES_CONFIGS = {
    "CH4": {
        "display_name": "CH4",
        "default_table": "CH4_M6_I1",
        "mode_label": "nu3",
        "mode_index": 2,
        "parse_quanta": parse_ch4_global_quanta,
    },
    "CO": {
        "display_name": "CO",
        "default_table": "CO_M5_I1",
        "mode_label": "nu1",
        "mode_index": 0,
        "parse_quanta": parse_co_global_quanta,
    },
}

if SPECIES_KEY not in SPECIES_CONFIGS:
    raise ValueError(f"Unsupported SPECIES_KEY: {SPECIES_KEY}")

SPECIES_CONFIG = SPECIES_CONFIGS[SPECIES_KEY]
DISPLAY_NAME = SPECIES_CONFIG["display_name"]
MODE_LABEL = SPECIES_CONFIG["mode_label"]
MODE_INDEX = SPECIES_CONFIG["mode_index"]
PARSE_QUANTA = SPECIES_CONFIG["parse_quanta"]

SOURCE_TABLE = SOURCE_TABLE_OVERRIDE.strip() or SPECIES_CONFIG["default_table"]

OUTPUT_DIR = str(ROOT_DIR / f"{DISPLAY_NAME.lower()}_{MODE_LABEL}_progressions")
CATEGORY_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "categories")


def parse_global_quanta(value: str) -> tuple[int, ...]:
    return PARSE_QUANTA(value)


def mode_pair(lower_raw: str, upper_raw: str) -> tuple[int, int] | None:
    lower_state = parse_global_quanta(lower_raw)
    upper_state = parse_global_quanta(upper_raw)

    if MODE_INDEX >= len(lower_state) or MODE_INDEX >= len(upper_state):
        raise ValueError(
            f"Mode index {MODE_INDEX} is out of range for labels {lower_raw!r} -> {upper_raw!r}"
        )

    lower_mode_q = lower_state[MODE_INDEX]
    upper_mode_q = upper_state[MODE_INDEX]

    if upper_mode_q <= lower_mode_q:
        return None

    if REQUIRE_SAME_OTHER_MODES:
        lower_other = tuple(value for index, value in enumerate(lower_state) if index != MODE_INDEX)
        upper_other = tuple(value for index, value in enumerate(upper_state) if index != MODE_INDEX)
        if lower_other != upper_other:
            return None

    if REQUIRE_UNIT_STEP and upper_mode_q != lower_mode_q + 1:
        return None

    return lower_mode_q, upper_mode_q


def exact_band_key(lower_raw: str, upper_raw: str) -> tuple[str, str] | None:
    pair = mode_pair(lower_raw, upper_raw)
    if pair is None:
        return None
    return lower_raw, upper_raw


def progression_category(lower_q: int, upper_q: int) -> str:
    delta_v = upper_q - lower_q
    if delta_v <= 0:
        return f"{MODE_LABEL}_transition"

    if delta_v == 1:
        base_label = "fundamental"
    else:
        overtone_order = delta_v - 1
        overtone_names = {
            1: "first",
            2: "second",
            3: "third",
            4: "fourth",
            5: "fifth",
        }
        order_label = overtone_names.get(overtone_order, f"{overtone_order}th")
        base_label = f"{order_label}_overtone"

    if lower_q > 0:
        return f"hot_{base_label}_band"
    return f"{base_label}_band"


def category_output_stem(category: str, lower_q: int, upper_q: int) -> str:
    return (
        f"{DISPLAY_NAME}_{MODE_LABEL}_{category}_{lower_q}_to_{upper_q}"
        f"_abs_{int(WN_MIN)}_{int(WN_MAX)}"
        f"_T{int(T_K)}K_P{P_TORR:g}Torr_L{int(PATH_LENGTH_CM)}cm"
    )


def progression_output_stem(lower_q: int, upper_q: int, lower_raw: str, upper_raw: str) -> str:
    return (
        f"{DISPLAY_NAME}_{MODE_LABEL}_{lower_q}_to_{upper_q}"
        f"_{safe_label_fragment(lower_raw)}_to_{safe_label_fragment(upper_raw)}"
        f"_abs_{int(WN_MIN)}_{int(WN_MAX)}"
        f"_T{int(T_K)}K_P{P_TORR:g}Torr_L{int(PATH_LENGTH_CM)}cm"
    )


def write_summary_csv(summary_rows: list[dict[str, object]]) -> str:
    summary_paths = [
        os.path.join(OUTPUT_DIR, f"{DISPLAY_NAME}_{MODE_LABEL}_progression_summary.csv"),
        os.path.join(OUTPUT_DIR, f"{DISPLAY_NAME}_{MODE_LABEL}_progression_summary_new.csv"),
        os.path.join(OUTPUT_DIR, f"{DISPLAY_NAME}_{MODE_LABEL}_progression_summary_latest.csv"),
    ]

    for summary_path in summary_paths:
        try:
            with open(summary_path, "w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "band_label",
                        "mode_pair",
                        "mode_pair_label",
                        "category",
                        "category_label",
                        "line_count",
                        "csv_path",
                        "html_path",
                    ],
                )
                writer.writeheader()
                writer.writerows(summary_rows)
            return summary_path
        except PermissionError:
            continue

    raise PermissionError(
        f"Could not write the {DISPLAY_NAME} {MODE_LABEL} progression summary CSV because all target files are locked."
    )


def write_category_summary_csv(summary_rows: list[dict[str, object]]) -> str:
    summary_paths = [
        os.path.join(CATEGORY_OUTPUT_DIR, f"{DISPLAY_NAME}_{MODE_LABEL}_category_summary.csv"),
        os.path.join(CATEGORY_OUTPUT_DIR, f"{DISPLAY_NAME}_{MODE_LABEL}_category_summary_new.csv"),
        os.path.join(CATEGORY_OUTPUT_DIR, f"{DISPLAY_NAME}_{MODE_LABEL}_category_summary_latest.csv"),
    ]

    for summary_path in summary_paths:
        try:
            with open(summary_path, "w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "category",
                        "mode_pair",
                        "mode_pair_label",
                        "category_label",
                        "line_count",
                        "band_count",
                        "csv_path",
                        "html_path",
                    ],
                )
                writer.writeheader()
                writer.writerows(summary_rows)
            return summary_path
        except PermissionError:
            continue

    raise PermissionError(
        f"Could not write the {DISPLAY_NAME} {MODE_LABEL} category summary CSV because all target files are locked."
    )


def collect_grouped_bands(
    global_lower: list[str],
    global_upper: list[str],
    line_strengths: list[float],
) -> dict[tuple[str, str], dict[str, float]]:
    grouped_bands = defaultdict(lambda: {"sum_sw": 0.0, "count": 0})

    for lower_raw, upper_raw, sw in zip(global_lower, global_upper, line_strengths):
        try:
            key = exact_band_key(lower_raw, upper_raw)
        except ValueError:
            continue

        if key is None:
            continue

        grouped_bands[key]["sum_sw"] += float(sw)
        grouped_bands[key]["count"] += 1

    return grouped_bands


def detect_progressions(grouped_bands: dict[tuple[str, str], dict[str, float]]) -> list[tuple[str, str]]:
    if not grouped_bands:
        raise RuntimeError(f"No exact {DISPLAY_NAME} {MODE_LABEL} transition groups matched the current filters.")
    return sorted(grouped_bands)


def collect_category_groups(
    grouped_bands: dict[tuple[str, str], dict[str, float]],
) -> dict[tuple[str, int, int], list[tuple[str, str, dict[str, float]]]]:
    category_groups = defaultdict(list)

    for (lower_raw, upper_raw), band_stats in grouped_bands.items():
        pair = mode_pair(lower_raw, upper_raw)
        if pair is None:
            continue
        lower_q, upper_q = pair
        category = progression_category(lower_q, upper_q)
        category_groups[(category, lower_q, upper_q)].append((lower_raw, upper_raw, band_stats))

    return category_groups


def render_absorbance_table(
    mol_id: int,
    iso_id: int,
    temp_table: str,
) -> tuple[np.ndarray, np.ndarray]:
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
    return nu, -np.log(trans)


def save_curve(out_csv: str, out_html: str, x_values: np.ndarray, y_values: np.ndarray, trace_name: str, title: str) -> None:
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x_values, y=y_values, mode="lines", name=trace_name))
    fig.update_layout(
        title=title,
        xaxis_title="Wavenumber (cm^-1)",
        yaxis_title="Absorbance",
        template="plotly_white",
    )
    fig.update_xaxes(range=[WN_MIN, WN_MAX])
    fig.write_html(out_html, include_plotlyjs="cdn")

    with open(out_csv, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["wavenumber_cm-1", "absorbance"])
        for wavenumber, value in zip(x_values, y_values):
            writer.writerow([wavenumber, value])


def main() -> None:
    hapi.db_begin(DB_DIR)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(CATEGORY_OUTPUT_DIR, exist_ok=True)

    header_path = os.path.join(DB_DIR, f"{SOURCE_TABLE}.header")
    if not os.path.exists(header_path):
        raise FileNotFoundError(f"Missing table header: {header_path}")

    match = re.search(r"_M(\d+)_I(\d+)", SOURCE_TABLE)
    if match is None:
        raise ValueError(f"Could not parse molecule/isotopologue IDs from {SOURCE_TABLE}")

    mol_id = int(match.group(1))
    iso_id = int(match.group(2))

    global_lower = hapi.getColumn(SOURCE_TABLE, "global_lower_quanta")
    global_upper = hapi.getColumn(SOURCE_TABLE, "global_upper_quanta")
    line_strengths = hapi.getColumn(SOURCE_TABLE, "sw")

    grouped_bands = collect_grouped_bands(global_lower, global_upper, line_strengths)
    progressions = detect_progressions(grouped_bands)
    category_groups = collect_category_groups(grouped_bands)

    exact_summary_rows = []

    if GENERATE_EXACT_BAND_OUTPUTS:
        for lower_raw, upper_raw in progressions:
            pair = mode_pair(lower_raw, upper_raw)
            if pair is None:
                continue
            lower_q, upper_q = pair
            total_lines = int(grouped_bands[(lower_raw, upper_raw)]["count"])
            if total_lines < MIN_LINES_PER_GROUP:
                continue

            category = progression_category(lower_q, upper_q)
            band_label = f"{clean_quanta_label(lower_raw)} -> {clean_quanta_label(upper_raw)}"
            pair_label = f"{MODE_LABEL} {lower_q}->{upper_q}"
            temp_table = (
                f"__{DISPLAY_NAME}_{MODE_LABEL}_{lower_q}_TO_{upper_q}_"
                f"{safe_label_fragment(lower_raw)}_TO_{safe_label_fragment(upper_raw)}__"
            )
            out_stem = progression_output_stem(lower_q, upper_q, lower_raw, upper_raw)
            out_csv = os.path.join(OUTPUT_DIR, f"{out_stem}.csv")
            out_html = os.path.join(OUTPUT_DIR, f"{out_stem}.html")

            try:
                hapi.select(
                    SOURCE_TABLE,
                    DestinationTableName=temp_table,
                    Conditions=(
                        "AND",
                        ("=", "global_lower_quanta", lower_raw),
                        ("=", "global_upper_quanta", upper_raw),
                    ),
                    Output=False,
                )

                if hapi.length(temp_table) == 0:
                    continue

                nu, absorbance = render_absorbance_table(mol_id, iso_id, temp_table)
                save_curve(
                    out_csv=out_csv,
                    out_html=out_html,
                    x_values=nu,
                    y_values=absorbance,
                    trace_name=f"{band_label} (N={total_lines})",
                    title=(
                        f"{DISPLAY_NAME} {band_label} | {pair_label} Absorbance (Voigt) "
                        f"{int(WN_MIN)}-{int(WN_MAX)} cm^-1 | T={T_K:g} K, P={P_TORR:g} Torr"
                        f" | x={MOLE_FRACTION:g}, L={PATH_LENGTH_CM:g} cm, S>{LINE_INTENSITY_THRESHOLD:.0e}"
                    ),
                )

                exact_summary_rows.append(
                    {
                        "band_label": band_label,
                        "mode_pair": pair_label,
                        "mode_pair_label": f"{MODE_LABEL} {lower_q}->{upper_q}",
                        "category": category,
                        "category_label": f"{MODE_LABEL} {lower_q}->{upper_q} | {category}",
                        "line_count": total_lines,
                        "csv_path": out_csv,
                        "html_path": out_html,
                    }
                )
                print(f"Saved: {out_csv}")
                print(f"Saved: {out_html}")
            finally:
                hapi.dropTable(temp_table)

        if exact_summary_rows:
            summary_path = write_summary_csv(exact_summary_rows)
            print(f"Saved: {summary_path}")

    if GENERATE_CATEGORY_OUTPUTS:
        category_summary_rows = []
        for (category, lower_q, upper_q), grouped_category_bands in sorted(category_groups.items()):
            temp_table = f"__{DISPLAY_NAME}_{MODE_LABEL}_{category}_{lower_q}_TO_{upper_q}__"
            out_stem = category_output_stem(category, lower_q, upper_q)
            out_csv = os.path.join(CATEGORY_OUTPUT_DIR, f"{out_stem}.csv")
            out_html = os.path.join(CATEGORY_OUTPUT_DIR, f"{out_stem}.html")
            total_lines = sum(int(band_stats["count"]) for _, _, band_stats in grouped_category_bands)
            pair_label = f"{MODE_LABEL} {lower_q}->{upper_q}"

            conditions = ["OR"]
            for lower_raw, upper_raw, _ in grouped_category_bands:
                conditions.append(
                    (
                        "AND",
                        ("=", "global_lower_quanta", lower_raw),
                        ("=", "global_upper_quanta", upper_raw),
                    )
                )

            try:
                hapi.select(
                    SOURCE_TABLE,
                    DestinationTableName=temp_table,
                    Conditions=tuple(conditions),
                    Output=False,
                )

                nu, absorbance = render_absorbance_table(mol_id, iso_id, temp_table)
                save_curve(
                    out_csv=out_csv,
                    out_html=out_html,
                    x_values=nu,
                    y_values=absorbance,
                    trace_name=(
                        f"{category} {MODE_LABEL} {lower_q}->{upper_q} "
                        f"(N={total_lines}, bands={len(grouped_category_bands)})"
                    ),
                    title=(
                        f"{DISPLAY_NAME} {MODE_LABEL} {category} | delta_v={upper_q-lower_q} ({lower_q}->{upper_q}) Absorbance (Voigt) "
                        f"{int(WN_MIN)}-{int(WN_MAX)} cm^-1 | T={T_K:g} K, P={P_TORR:g} Torr"
                        f" | x={MOLE_FRACTION:g}, L={PATH_LENGTH_CM:g} cm, S>{LINE_INTENSITY_THRESHOLD:.0e}"
                    ),
                )

                category_summary_rows.append(
                    {
                        "category": category,
                        "mode_pair": pair_label,
                        "mode_pair_label": pair_label,
                        "category_label": f"{pair_label} | {category}",
                        "line_count": total_lines,
                        "band_count": len(grouped_category_bands),
                        "csv_path": out_csv,
                        "html_path": out_html,
                    }
                )
                print(f"Saved: {out_csv}")
                print(f"Saved: {out_html}")
            finally:
                hapi.dropTable(temp_table)

        if category_summary_rows:
            category_summary_path = write_category_summary_csv(category_summary_rows)
            print(f"Saved: {category_summary_path}")


if __name__ == "__main__":
    main()
