"""
Sort each extracted HITEMP CH4 band-text file by line intensity.

Inputs:
- `ch4_nu3_progressions/hitemp_band_line_texts/*.txt`: source per-band HITEMP
  text files containing fixed-width HITRAN/HITEMP rows.
- `hitemp_db/CH4/06_HITEMP2020.par_2500-3500.header`: header metadata used to
  locate the fixed-width `sw` and `nu` columns.

Outputs:
- `ch4_nu3_progressions/hitemp_band_line_texts_sorted_by_intensity/*.txt`:
  one output file per source file, with rows sorted by descending `sw`.
- `ch4_nu3_progressions/hitemp_band_line_texts_sorted_by_intensity/strongest_lines_summary.csv`:
  one row per source file with the strongest line information.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT_DIR / "ch4_nu3_progressions" / "hitemp_band_line_texts"
OUTPUT_DIR = ROOT_DIR / "ch4_nu3_progressions" / "hitemp_band_line_texts_sorted_by_intensity"
HEADER_PATH = ROOT_DIR / "hitemp_db" / "CH4" / "06_HITEMP2020.par_2500-3500.header"
SUMMARY_PATH = OUTPUT_DIR / "strongest_lines_summary.csv"


def parse_fixed_width(format_string: str) -> int:
    match = re.fullmatch(r"%(\d+)(?:\.\d+)?[A-Za-z]", format_string)
    if match is None:
        raise ValueError(f"Unsupported fixed-width format: {format_string}")
    return int(match.group(1))


def load_positions_and_widths() -> tuple[dict[str, int], dict[str, int]]:
    with HEADER_PATH.open("r", encoding="utf-8") as handle:
        header = json.load(handle)
    positions = header["position"]
    widths = {name: parse_fixed_width(fmt) for name, fmt in header["format"].items()}
    return positions, widths


def parse_float_field(line: str, start: int, width: int, field_name: str) -> float:
    raw = line[start : start + width]
    try:
        return float(raw)
    except ValueError as exc:
        raise RuntimeError(f"Could not parse {field_name} from line: {line!r}") from exc


def sort_one_file(
    source_path: Path,
    destination_path: Path,
    sw_start: int,
    sw_width: int,
    nu_start: int,
    nu_width: int,
) -> dict[str, object]:
    rows: list[tuple[float, float, str]] = []
    with source_path.open("r", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            sw = parse_float_field(line, sw_start, sw_width, "sw")
            nu = parse_float_field(line, nu_start, nu_width, "nu")
            rows.append((sw, nu, line))

    rows.sort(key=lambda item: item[0], reverse=True)

    destination_path.parent.mkdir(parents=True, exist_ok=True)
    with destination_path.open("w", encoding="utf-8", newline="\n") as handle:
        for _sw, _nu, line in rows:
            handle.write(line)
            handle.write("\n")

    if not rows:
        return {
            "source_file": source_path.name,
            "sorted_file": destination_path.name,
            "line_count": 0,
            "strongest_sw": "",
            "strongest_nu": "",
            "strongest_line": "",
        }

    strongest_sw, strongest_nu, strongest_line = rows[0]
    return {
        "source_file": source_path.name,
        "sorted_file": destination_path.name,
        "line_count": len(rows),
        "strongest_sw": f"{strongest_sw:.6e}",
        "strongest_nu": f"{strongest_nu:.6f}",
        "strongest_line": strongest_line,
    }


def main() -> None:
    if not INPUT_DIR.exists():
        raise FileNotFoundError(f"Input directory does not exist: {INPUT_DIR}")
    if not HEADER_PATH.exists():
        raise FileNotFoundError(f"Header file does not exist: {HEADER_PATH}")

    positions, widths = load_positions_and_widths()
    sw_start = positions["sw"]
    sw_width = widths["sw"]
    nu_start = positions["nu"]
    nu_width = widths["nu"]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    summaries: list[dict[str, object]] = []
    for source_path in sorted(INPUT_DIR.glob("*.txt")):
        destination_path = OUTPUT_DIR / source_path.name
        summary = sort_one_file(
            source_path=source_path,
            destination_path=destination_path,
            sw_start=sw_start,
            sw_width=sw_width,
            nu_start=nu_start,
            nu_width=nu_width,
        )
        summaries.append(summary)
        print(
            f"SORTED: {source_path.name} | "
            f"lines={summary['line_count']} | strongest_sw={summary['strongest_sw']} | "
            f"strongest_nu={summary['strongest_nu']}"
        )

    with SUMMARY_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "source_file",
                "sorted_file",
                "line_count",
                "strongest_sw",
                "strongest_nu",
                "strongest_line",
            ],
        )
        writer.writeheader()
        writer.writerows(summaries)

    print(f"WROTE: {SUMMARY_PATH}")


if __name__ == "__main__":
    main()
