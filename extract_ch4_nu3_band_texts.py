"""Extract exact CH4 nu3 band lines into per-band text files.

Inputs:
- ``hitran_db/CH4_M6_I1.header``: JSON metadata that defines the fixed-width
  field positions for the HITRAN text records.
- ``hitran_db/CH4_M6_I1.data``: source fixed-width line list to scan.
- ``ch4_nu3_progressions/CH4_nu3_progression_summary.csv``: summary table that
  lists the exact target ``band_label`` values and expected ``line_count`` for
  the 9 CH4 nu3 bands to export. The extractor uses ``band_label`` to decide
  which lines to copy; ``line_count`` is used only afterward as a verification
  check that each exported text file contains the full expected number of rows.

Outputs:
- ``ch4_nu3_progressions/band_line_texts/*.txt``: one plain-text file per
  target band, each containing the original matching records copied verbatim
  from ``CH4_M6_I1.data``.
- Console status lines reporting the expected and actual line count for each
  exported band file.

The script raises ``RuntimeError`` if the summary CSV does not contain exactly
9 target bands or if any exported file's line count does not match the expected
count from the CSV.
"""

import csv
import json
import os
import re
from collections import Counter


DB_DIR = "hitran_db"
SOURCE_TABLE = "CH4_M6_I1"
SUMMARY_CSV = os.path.join("ch4_nu3_progressions", "CH4_nu3_progression_summary.csv")
OUTPUT_DIR = os.path.join("ch4_nu3_progressions", "band_line_texts")


def clean_quanta_label(value: str) -> str:
    text = value.strip()
    return text if text else "000"


def safe_label_fragment(value: str) -> str:
    text = clean_quanta_label(value)
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_") or "000"


def parse_fixed_width(format_string: str) -> int:
    match = re.fullmatch(r"%(\d+)(?:\.\d+)?[A-Za-z]", format_string)
    if match is None:
        raise ValueError(f"Unsupported fixed-width format: {format_string}")
    return int(match.group(1))


def load_header_metadata() -> tuple[dict[str, int], dict[str, int]]:
    header_path = os.path.join(DB_DIR, f"{SOURCE_TABLE}.header")
    with open(header_path, "r", encoding="utf-8") as handle:
        header = json.load(handle)

    positions = header["position"]
    widths = {name: parse_fixed_width(fmt) for name, fmt in header["format"].items()}
    return positions, widths


def load_target_bands() -> tuple[list[dict[str, object]], dict[str, int]]:
    rows: list[dict[str, object]] = []
    expected_counts: dict[str, int] = {}

    with open(SUMMARY_CSV, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            band_label = str(row["band_label"]).strip()
            line_count = int(row["line_count"])
            rows.append(row)
            expected_counts[band_label] = line_count

    if len(expected_counts) != 9:
        raise RuntimeError(f"Expected 9 exact CH4 nu3 band labels, found {len(expected_counts)}")

    return rows, expected_counts


def extract_band_label(line: str, positions: dict[str, int], widths: dict[str, int]) -> str:
    upper_start = positions["global_upper_quanta"]
    upper_end = upper_start + widths["global_upper_quanta"]
    lower_start = positions["global_lower_quanta"]
    lower_end = lower_start + widths["global_lower_quanta"]

    upper = clean_quanta_label(line[upper_start:upper_end])
    lower = clean_quanta_label(line[lower_start:lower_end])
    return f"{lower} -> {upper}"


def output_path_for_band(band_label: str) -> str:
    lower_label, upper_label = [part.strip() for part in band_label.split("->", maxsplit=1)]
    filename = (
        f"{SOURCE_TABLE}_{safe_label_fragment(lower_label)}"
        f"_to_{safe_label_fragment(upper_label)}.txt"
    )
    return os.path.join(OUTPUT_DIR, filename)


def main() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    _, expected_counts = load_target_bands()
    positions, widths = load_header_metadata()

    matched_counts: Counter[str] = Counter()
    output_handles: dict[str, object] = {}

    try:
        for band_label in expected_counts:
            output_handles[band_label] = open(
                output_path_for_band(band_label),
                "w",
                encoding="utf-8",
                newline="",
            )

        data_path = os.path.join(DB_DIR, f"{SOURCE_TABLE}.data")
        with open(data_path, "r", encoding="utf-8", newline="") as handle:
            for line in handle:
                band_label = extract_band_label(line, positions, widths)
                if band_label not in output_handles:
                    continue
                output_handles[band_label].write(line)
                matched_counts[band_label] += 1
    finally:
        for handle in output_handles.values():
            handle.close()

    mismatches = []
    for band_label, expected_count in expected_counts.items():
        actual_count = matched_counts[band_label]
        out_path = output_path_for_band(band_label)
        # band_label controls extraction; line_count only validates completeness.
        status = "OK" if actual_count == expected_count else "MISMATCH"
        print(f"{status}: {band_label} | expected={expected_count} actual={actual_count} | {out_path}")
        if actual_count != expected_count:
            mismatches.append((band_label, expected_count, actual_count))

    if mismatches:
        raise RuntimeError(f"Line-count mismatches found for {len(mismatches)} band files.")


if __name__ == "__main__":
    main()
