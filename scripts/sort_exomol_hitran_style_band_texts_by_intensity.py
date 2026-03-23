"""
Sort merged ExoMol HITRAN-style band text files by line intensity.

The script reads all `.txt` files from a source directory, preserves the
comment preamble and tab-delimited header row, and sorts the data rows by
`line_intensity_cm_per_molecule` in descending order.
"""

from __future__ import annotations

import argparse
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[1]


DEFAULT_INPUT_DIR = ROOT_DIR / "exomol_ch4_mm_pure_nu3_band_texts_hitran_style"
DEFAULT_OUTPUT_DIR = ROOT_DIR / "exomol_ch4_mm_pure_nu3_band_texts_hitran_style_sorted"
INTENSITY_COLUMN = "line_intensity_cm_per_molecule"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sort HITRAN-style ExoMol band text files by line intensity.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Directory containing HITRAN-style ExoMol band text files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for the sorted text files.",
    )
    return parser.parse_args()


def sort_one_file(input_path: Path, output_path: Path) -> int:
    preamble: list[str] = []
    header_line: str | None = None
    data_rows: list[tuple[float, str]] = []

    with input_path.open("r", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("#"):
                preamble.append(line)
                continue
            if header_line is None:
                header_line = line
                columns = header_line.split("\t")
                try:
                    intensity_index = columns.index(INTENSITY_COLUMN)
                except ValueError as exc:
                    raise RuntimeError(f"Missing {INTENSITY_COLUMN!r} in header of {input_path}") from exc
                continue

            parts = line.split("\t")
            if len(parts) <= intensity_index:
                raise RuntimeError(f"Malformed data row in {input_path}: {line}")
            data_rows.append((float(parts[intensity_index]), line))

    if header_line is None:
        raise RuntimeError(f"No header row found in {input_path}")

    data_rows.sort(key=lambda item: item[0], reverse=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        for line in preamble:
            handle.write(line + "\n")
        handle.write(header_line + "\n")
        for _, line in data_rows:
            handle.write(line + "\n")

    return len(data_rows)


def main() -> None:
    args = parse_args()
    if not args.input_dir.exists():
        raise FileNotFoundError(f"Missing input directory: {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    txt_files = sorted(args.input_dir.glob("*.txt"))
    if not txt_files:
        raise RuntimeError(f"No .txt files found in {args.input_dir}")

    total_rows = 0
    for txt_path in txt_files:
        out_path = args.output_dir / txt_path.name
        row_count = sort_one_file(txt_path, out_path)
        total_rows += row_count
        print(f"sorted {txt_path.name} ({row_count:,} rows)")

    print(f"sorted files: {len(txt_files):,}")
    print(f"sorted rows: {total_rows:,}")
    print(f"saved {args.output_dir}")


if __name__ == "__main__":
    main()
