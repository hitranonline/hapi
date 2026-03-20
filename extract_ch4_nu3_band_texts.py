"""Extract exact CH4 nu3 HITEMP band lines into per-band text files.

Inputs:
- ``hitemp_db/CH4/06_HITEMP2020.par_2500-3500.header``: JSON metadata that
  defines the fixed-width field positions for the extracted HITEMP text records.
- ``hitemp_db/CH4/06_HITEMP2020.par_2500-3500.par``: source fixed-width HITEMP
  line list to scan.
- ``ch4_nu3_progressions/CH4_nu3_progression_summary.csv``: summary table that
  lists the exact target ``band_label`` values to export. For the HITEMP
  workflow, the CSV's ``line_count`` values are treated as historical HITRAN
  reference counts only and are not enforced.

Outputs:
- ``ch4_nu3_progressions/hitemp_band_line_texts/*.txt``: one plain-text file
  per target band, each containing the original matching records copied
  verbatim from the HITEMP `.par` file.
- Console status lines reporting the extracted line count for each exported
  band file.

The script raises ``RuntimeError`` if the summary CSV does not contain exactly
9 target bands.
"""

import csv
import json
import os
import re
from collections import Counter
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parent
SOURCE_HEADER = ROOT_DIR / "hitemp_db" / "CH4" / "06_HITEMP2020.par_2500-3500.header"
SOURCE_DATA = ROOT_DIR / "hitemp_db" / "CH4" / "06_HITEMP2020.par_2500-3500.par"
SOURCE_STEM = SOURCE_DATA.stem
SUMMARY_CSV = ROOT_DIR / "ch4_nu3_progressions" / "CH4_nu3_progression_summary.csv"
OUTPUT_DIR = ROOT_DIR / "ch4_nu3_progressions" / "hitemp_band_line_texts"


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
    with SOURCE_HEADER.open("r", encoding="utf-8") as handle:
        header = json.load(handle)

    positions = header["position"]
    widths = {name: parse_fixed_width(fmt) for name, fmt in header["format"].items()}
    return positions, widths


def load_target_bands() -> tuple[list[dict[str, object]], dict[str, int]]:
    rows: list[dict[str, object]] = []
    expected_counts: dict[str, int] = {}

    with SUMMARY_CSV.open("r", encoding="utf-8-sig", newline="") as handle:
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
        f"{SOURCE_STEM}_{safe_label_fragment(lower_label)}"
        f"_to_{safe_label_fragment(upper_label)}.txt"
    )
    return str(OUTPUT_DIR / filename)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

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

        with SOURCE_DATA.open("r", encoding="utf-8", newline="") as handle:
            for line in handle:
                band_label = extract_band_label(line, positions, widths)
                if band_label not in output_handles:
                    continue
                output_handles[band_label].write(line)
                matched_counts[band_label] += 1
    finally:
        for handle in output_handles.values():
            handle.close()

    for band_label, expected_count in expected_counts.items():
        actual_count = matched_counts[band_label]
        out_path = output_path_for_band(band_label)
        # For HITEMP, the CSV's line_count is a historical HITRAN reference only.
        print(
            f"EXTRACTED: {band_label} | "
            f"hitemp_actual={actual_count} | hitran_reference={expected_count} | {out_path}"
        )


if __name__ == "__main__":
    main()
