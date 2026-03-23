"""
Merge fine-grained ExoMol pure-nu3 band text files into HITRAN-style groups.

Input
-----
- `exomol_ch4_mm_pure_nu3_band_texts/exomol_pure_nu3_band_text_summary.csv`
  produced by `extract_exomol_ch4_mm_pure_nu3_band_texts.py`

Behavior
--------
- Groups exported ExoMol band text files by `band_label_hitran_style`
- Ignores all fine-grained MM labels other than `n1`, `n2`, `n3`, `n4`, and
  `Gtot` when collapsing to the coarse HITRAN-style grouping
- Concatenates data rows from each matching source file into one merged file

Outputs
-------
- `exomol_ch4_mm_pure_nu3_band_texts_hitran_style/*.txt`
- `exomol_ch4_mm_pure_nu3_band_texts_hitran_style/exomol_pure_nu3_band_text_summary_hitran_style.csv`
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]


DEFAULT_INPUT_DIR = ROOT_DIR / "exomol_ch4_mm_pure_nu3_band_texts"
DEFAULT_OUTPUT_DIR = ROOT_DIR / "exomol_ch4_mm_pure_nu3_band_texts_hitran_style"
DEFAULT_SUMMARY_CSV = DEFAULT_INPUT_DIR / "exomol_pure_nu3_band_text_summary.csv"
OUTPUT_SUMMARY_NAME = "exomol_pure_nu3_band_text_summary_hitran_style.csv"
HITRAN_STYLE_MAP = {
    "A1": "1A1",
    "A2": "1A2",
    "E": "1E",
    "F1": "1F1",
    "F2": "1F2",
    "T1": "1F1",
    "T2": "1F2",
}
DATA_COLUMNS = [
    "upper_id",
    "lower_id",
    "upper_energy_cm-1",
    "lower_energy_cm-1",
    "wavenumber_cm-1",
    "einstein_A_s-1",
    "line_intensity_cm_per_molecule",
    "local_upper_quanta",
    "local_lower_quanta",
]
DATA_HEADER = "\t".join(DATA_COLUMNS)
STATE_LABEL_PATTERN = re.compile(
    r"^\s*"
    r"(?P<n1>\d+)\s+(?P<n2>\d+)\s+(?P<n3>\d+)\s+(?P<n4>\d+)\s+"
    r"Gtot=(?P<gtot>\S+)"
    r"(?:\s+.*)?\s*$"
)
COARSE_STATE_LABEL_PATTERN = re.compile(
    r"^\s*"
    r"(?P<n1>\d+)\s+(?P<n2>\d+)\s+(?P<n3>\d+)\s+(?P<n4>\d+)\s+"
    r"(?P<gtot>\S+)\s*$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge fine-grained ExoMol pure-nu3 band text files into HITRAN-style band groups.",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=DEFAULT_SUMMARY_CSV,
        help="Manifest CSV produced by extract_exomol_ch4_mm_pure_nu3_band_texts.py.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for merged HITRAN-style text files.",
    )
    return parser.parse_args()


def safe_label_fragment(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_") or "000"


def output_path_for_band(output_dir: Path, band_label_hitran_style: str, reference_temperature_k: str) -> Path:
    lower_label, upper_label = [part.strip() for part in band_label_hitran_style.split("->", maxsplit=1)]
    filename = (
        "EXOMOL_CH4_MM_pure_nu3_hitran_style_"
        f"{safe_label_fragment(lower_label)}_to_{safe_label_fragment(upper_label)}"
        f"_T{reference_temperature_k}K.txt"
    )
    return output_dir / filename


def iter_data_lines(txt_path: Path) -> tuple[list[str], str]:
    preamble: list[str] = []
    found_header = False
    with txt_path.open("r", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("#"):
                preamble.append(line)
                continue

            if not found_header:
                if line != DATA_HEADER:
                    raise RuntimeError(f"Unexpected data header in {txt_path}: {line}")
                found_header = True
                continue
            yield preamble, line

    if not found_header:
        raise RuntimeError(f"Missing data header in {txt_path}")


def load_rows(summary_csv: Path) -> list[dict[str, str]]:
    with summary_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def coarse_labels_from_state_label(state_label: str) -> tuple[str, str]:
    match = STATE_LABEL_PATTERN.fullmatch(state_label)
    if match is None:
        match = COARSE_STATE_LABEL_PATTERN.fullmatch(state_label)
    if match is None:
        raise ValueError(f"Unexpected ExoMol state label: {state_label}")

    exomol_label = (
        f"{match.group('n1')} {match.group('n2')} {match.group('n3')} "
        f"{match.group('n4')} {match.group('gtot')}"
    )
    hitran_label = (
        f"{match.group('n1')} {match.group('n2')} {match.group('n3')} "
        f"{match.group('n4')} {HITRAN_STYLE_MAP.get(match.group('gtot'), match.group('gtot'))}"
    )
    return exomol_label, hitran_label


def coarse_labels_from_band_label(band_label: str) -> tuple[str, str]:
    lower_label, upper_label = [part.strip() for part in band_label.split("->", maxsplit=1)]
    lower_exomol, lower_hitran = coarse_labels_from_state_label(lower_label)
    upper_exomol, upper_hitran = coarse_labels_from_state_label(upper_label)
    return f"{lower_exomol} -> {upper_exomol}", f"{lower_hitran} -> {upper_hitran}"


def validate_group(rows: list[dict[str, str]]) -> tuple[str, str, str]:
    reference_temperatures = {row["reference_temperature_k"].strip() for row in rows}
    if len(reference_temperatures) != 1:
        raise RuntimeError(f"Expected one reference temperature per merged band, found {sorted(reference_temperatures)}")

    wn_mins = {row["wn_min_cm-1"].strip() for row in rows}
    wn_maxs = {row["wn_max_cm-1"].strip() for row in rows}
    if len(wn_mins) != 1 or len(wn_maxs) != 1:
        raise RuntimeError("Expected one wavenumber window per merged band group")

    return reference_temperatures.pop(), wn_mins.pop(), wn_maxs.pop()


def write_merged_band(
    output_path: Path,
    band_label_exomol: str,
    band_label_hitran_style: str,
    rows: list[dict[str, str]],
    *,
    reference_temperature_k: str,
    wn_min: str,
    wn_max: str,
) -> int:
    total_lines = 0
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("# ExoMol MM pure-nu3 HITRAN-style merged band export\n")
        handle.write("# source=merge_exomol_pure_nu3_band_texts_hitran_style.py\n")
        handle.write(f"# reference_temperature_k={reference_temperature_k}\n")
        handle.write(f"# band_label_exomol={band_label_exomol}\n")
        handle.write(f"# hitran_style_band_label={band_label_hitran_style}\n")
        if wn_min and wn_max:
            handle.write(f"# wavenumber_window_cm-1={wn_min} to {wn_max}\n")
        handle.write(f"# merged_source_file_count={len(rows)}\n")
        handle.write(DATA_HEADER + "\n")

        for row in rows:
            txt_path = Path(row["txt_path"])
            for _, data_line in iter_data_lines(txt_path):
                handle.write(data_line + "\n")
                total_lines += 1

    return total_lines


def write_summary_csv(output_dir: Path, rows: list[dict[str, object]]) -> Path:
    summary_path = output_dir / OUTPUT_SUMMARY_NAME
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "band_label_hitran_style",
                "band_label_exomol",
                "mode_pair",
                "line_count",
                "reference_temperature_k",
                "wn_min_cm-1",
                "wn_max_cm-1",
                "merged_source_file_count",
                "txt_path",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return summary_path


def main() -> None:
    args = parse_args()
    if not args.summary_csv.exists():
        raise FileNotFoundError(f"Missing summary CSV: {args.summary_csv}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    source_rows = load_rows(args.summary_csv)
    if not source_rows:
        raise RuntimeError(f"No rows found in {args.summary_csv}")

    grouped_rows: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in source_rows:
        band_label_exomol, band_label_hitran_style = coarse_labels_from_band_label(row["band_label"].strip())
        grouped_rows[(band_label_exomol, band_label_hitran_style)].append(row)

    summary_rows: list[dict[str, object]] = []
    for band_label_exomol, band_label_hitran_style in sorted(grouped_rows):
        rows = grouped_rows[(band_label_exomol, band_label_hitran_style)]
        reference_temperature_k, wn_min, wn_max = validate_group(rows)
        out_path = output_path_for_band(args.output_dir, band_label_hitran_style, reference_temperature_k)
        merged_line_count = write_merged_band(
            out_path,
            band_label_exomol,
            band_label_hitran_style,
            rows,
            reference_temperature_k=reference_temperature_k,
            wn_min=wn_min,
            wn_max=wn_max,
        )

        source_line_count = sum(int(row["line_count"]) for row in rows)
        if merged_line_count != source_line_count:
            raise RuntimeError(
                f"Merged line-count mismatch for {band_label_hitran_style}: "
                f"expected {source_line_count}, wrote {merged_line_count}"
            )

        summary_rows.append(
            {
                "band_label_hitran_style": band_label_hitran_style,
                "band_label_exomol": band_label_exomol,
                "mode_pair": rows[0]["mode_pair"].strip(),
                "line_count": merged_line_count,
                "reference_temperature_k": reference_temperature_k,
                "wn_min_cm-1": wn_min,
                "wn_max_cm-1": wn_max,
                "merged_source_file_count": len(rows),
                "txt_path": str(out_path),
            }
        )

    summary_path = write_summary_csv(args.output_dir, summary_rows)
    print(f"merged HITRAN-style bands: {len(summary_rows):,}")
    print(f"saved {summary_path}")


if __name__ == "__main__":
    main()
