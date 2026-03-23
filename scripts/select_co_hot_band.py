import csv
import os
import re
from pathlib import Path
from typing import Any

from _bootstrap import ensure_repo_root

ensure_repo_root()

import hapi

# -------------------------
# USER SETTINGS
# -------------------------
ROOT_DIR = Path(__file__).resolve().parents[1]
DB_DIR = str(ROOT_DIR / "hitran_db")
OUTPUT_DIR = str(ROOT_DIR / "co_hot_bands")

MOLECULE_ID = 5
ISOTOPOLOGUE_ID = 1
MOLECULE_FORMULA = "CO"

SOURCE_TABLE = f"{MOLECULE_FORMULA}_M{MOLECULE_ID}_I{ISOTOPOLOGUE_ID}"
DESTINATION_TABLE = "HotBands"

NU_MIN = 2000.0
NU_MAX = 4500.0

LOWER_V = 2
UPPER_V = 3

FORCE_FETCH = False


def nested_condition(operator: str, conditions: list[tuple[Any, ...]]) -> tuple[Any, ...]:
    if not conditions:
        raise ValueError("At least one condition is required.")

    combined = conditions[0]
    for condition in conditions[1:]:
        combined = (operator, combined, condition)
    return combined


def extract_vibrational_quantum(value: Any) -> int | None:
    if isinstance(value, bool):
        return None

    if isinstance(value, int):
        return value

    if isinstance(value, float):
        rounded = int(round(value))
        if abs(value - rounded) < 1.0e-9:
            return rounded
        return None

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


def unique_in_order(values: list[Any]) -> list[Any]:
    unique_values: list[Any] = []
    seen = set()
    for value in values:
        key = repr(value)
        if key in seen:
            continue
        seen.add(key)
        unique_values.append(value)
    return unique_values


def ensure_source_table() -> None:
    hapi.db_begin(DB_DIR)

    header_path = os.path.join(DB_DIR, f"{SOURCE_TABLE}.header")
    data_path = os.path.join(DB_DIR, f"{SOURCE_TABLE}.data")
    table_exists = os.path.exists(header_path) and os.path.exists(data_path)

    if table_exists and not FORCE_FETCH:
        print(f"Using existing table: {SOURCE_TABLE}")
        return

    print(
        f"Fetching {SOURCE_TABLE}: molecule={MOLECULE_ID}, isotopologue={ISOTOPOLOGUE_ID}, "
        f"nu=[{NU_MIN}, {NU_MAX}] cm^-1"
    )
    hapi.fetch(SOURCE_TABLE, MOLECULE_ID, ISOTOPOLOGUE_ID, NU_MIN, NU_MAX)


def get_column_map(table_name: str) -> dict[str, str]:
    header = hapi.getTableHeader(table_name)
    return {column_name.lower(): column_name for column_name in header["order"]}


def export_table_to_csv(table_name: str, csv_path: str) -> None:
    header = hapi.getTableHeader(table_name)
    column_names = list(header["order"])
    column_data = [hapi.getColumn(table_name, column_name) for column_name in column_names]
    row_count = hapi.length(table_name)

    with open(csv_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(column_names)
        for row_index in range(row_count):
            writer.writerow([column[row_index] for column in column_data])


def select_with_direct_v_columns(upper_column: str, lower_column: str) -> int:
    hapi.select(
        SOURCE_TABLE,
        DestinationTableName=DESTINATION_TABLE,
        Conditions=(
            "AND",
            ("=", upper_column, UPPER_V),
            ("=", lower_column, LOWER_V),
        ),
        Output=False,
    )
    return hapi.length(DESTINATION_TABLE)


def select_with_quantum_labels(upper_column: str, lower_column: str) -> int:
    upper_values = hapi.getColumn(SOURCE_TABLE, upper_column)
    lower_values = hapi.getColumn(SOURCE_TABLE, lower_column)

    matched_upper_raw: list[Any] = []
    matched_lower_raw: list[Any] = []

    for upper_raw, lower_raw in zip(upper_values, lower_values):
        upper_v = extract_vibrational_quantum(upper_raw)
        lower_v = extract_vibrational_quantum(lower_raw)
        if upper_v == UPPER_V and lower_v == LOWER_V:
            matched_upper_raw.append(upper_raw)
            matched_lower_raw.append(lower_raw)

    if not matched_upper_raw:
        return 0

    upper_conditions = [("=", upper_column, value) for value in unique_in_order(matched_upper_raw)]
    lower_conditions = [("=", lower_column, value) for value in unique_in_order(matched_lower_raw)]

    condition = (
        "AND",
        nested_condition("OR", upper_conditions),
        nested_condition("OR", lower_conditions),
    )

    hapi.select(
        SOURCE_TABLE,
        DestinationTableName=DESTINATION_TABLE,
        Conditions=condition,
        Output=False,
    )
    return hapi.length(DESTINATION_TABLE)


def main() -> None:
    ensure_source_table()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    column_map = get_column_map(SOURCE_TABLE)
    available_columns = list(hapi.getTableHeader(SOURCE_TABLE)["order"])
    print("Available columns:")
    print(", ".join(available_columns))

    row_count = 0
    if "vp" in column_map and "vpp" in column_map:
        print(f"Selecting with direct vibrational columns: {column_map['vpp']} -> {column_map['vp']}")
        row_count = select_with_direct_v_columns(
            upper_column=column_map["vp"],
            lower_column=column_map["vpp"],
        )
    elif "global_upper_quanta" in column_map and "global_lower_quanta" in column_map:
        print(
            "Selecting by parsing standard HITRAN quantum-label columns: "
            f"{column_map['global_lower_quanta']} -> {column_map['global_upper_quanta']}"
        )
        row_count = select_with_quantum_labels(
            upper_column=column_map["global_upper_quanta"],
            lower_column=column_map["global_lower_quanta"],
        )
    else:
        raise RuntimeError(
            "Could not find a usable vibrational-quantum column pair. "
            f"Columns found: {available_columns}"
        )

    if row_count == 0:
        print(f"No {DESTINATION_TABLE} rows matched v''={LOWER_V} and v'={UPPER_V}.")
        return

    txt_path = os.path.join(OUTPUT_DIR, f"{DESTINATION_TABLE}.txt")
    csv_path = os.path.join(OUTPUT_DIR, f"{DESTINATION_TABLE}.csv")

    hapi.outputTable(DESTINATION_TABLE, File=txt_path)
    export_table_to_csv(DESTINATION_TABLE, csv_path)

    print(f"Created table: {DESTINATION_TABLE} ({row_count} rows)")
    print(f"Saved: {txt_path}")
    print(f"Saved: {csv_path}")


if __name__ == "__main__":
    main()
