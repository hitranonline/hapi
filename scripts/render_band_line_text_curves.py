"""
Build absorbance CSV/HTML outputs from raw band-line text files.

This script leaves the existing plotting scripts untouched. It reads the
fixed-width rows stored in `ch4_nu3_progressions/band_line_texts`, rebuilds a
temporary in-memory HAPI table for each text file, and renders one absorbance
curve per file.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from _bootstrap import ensure_repo_root

ensure_repo_root()

import hapi
import plot_vibrational_mode_progressions as progression_plot


ROOT_DIR = Path(__file__).resolve().parents[1]
DB_DIR = ROOT_DIR / "hitran_db"
SOURCE_TABLE = "CH4_M6_I1"
TXT_DIR = ROOT_DIR / "ch4_nu3_progressions" / "band_line_texts"
OUTPUT_DIR = ROOT_DIR / "ch4_nu3_progressions" / "band_line_text_exports"
EXPECTED_FILE_COUNT = 9

WN_MIN = 2500.0
WN_MAX = 3500.0
WN_STEP = 0.001

T_K = 600.0
P_TORR = 3.0
P_ATM = P_TORR / 760.0

MOLE_FRACTION = 0.008
PATH_LENGTH_CM = 100.0
LINE_INTENSITY_THRESHOLD = 1.0e-23


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render one absorbance CSV/HTML pair per band-line text file.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Render only the first N text files after sorting by name.",
    )
    return parser.parse_args()


def parse_source_ids(table_name: str) -> tuple[int, int]:
    match = re.fullmatch(r".*_M(\d+)_I(\d+)", table_name)
    if match is None:
        raise ValueError(f"Could not parse molecule/isotopologue IDs from {table_name}")
    return int(match.group(1)), int(match.group(2))


def band_label_from_path(txt_path: Path) -> str:
    stem = txt_path.stem
    prefix = f"{SOURCE_TABLE}_"
    if stem.startswith(prefix):
        stem = stem[len(prefix) :]
    return stem.replace("_to_", " -> ").replace("_", " ")


def temp_table_name(txt_path: Path) -> str:
    return "__band_text_" + re.sub(r"[^A-Za-z0-9]+", "_", txt_path.stem).strip("_") + "__"


def bootstrap_source_schema() -> None:
    hapi.db_begin(str(DB_DIR))
    # Load only the header/default schema; the full source data are not needed here.
    hapi.storage2cache(SOURCE_TABLE, nlines=0)


def configure_progression_helpers() -> None:
    progression_plot.WN_MIN = WN_MIN
    progression_plot.WN_MAX = WN_MAX
    progression_plot.WN_STEP = WN_STEP
    progression_plot.T_K = T_K
    progression_plot.P_TORR = P_TORR
    progression_plot.P_ATM = P_ATM
    progression_plot.MOLE_FRACTION = MOLE_FRACTION
    progression_plot.PATH_LENGTH_CM = PATH_LENGTH_CM
    progression_plot.LINE_INTENSITY_THRESHOLD = LINE_INTENSITY_THRESHOLD


def build_temp_table_from_text(txt_path: Path, temp_table: str) -> int:
    hapi.dropTable(temp_table)
    hapi.createTable(temp_table, hapi.getDefaultRowObject(SOURCE_TABLE))

    line_count = 0
    with txt_path.open("r", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line.strip():
                continue
            row_object = hapi.getRowObjectFromString(line, SOURCE_TABLE)
            hapi.addRowObject(row_object, temp_table)
            line_count += 1

    hapi.LOCAL_TABLE_CACHE[temp_table]["header"]["number_of_rows"] = line_count
    if line_count == 0:
        raise RuntimeError(f"No usable HITRAN rows found in {txt_path}")
    return line_count


def main() -> None:
    args = parse_args()

    txt_files = sorted(TXT_DIR.glob("*.txt"))
    if len(txt_files) != EXPECTED_FILE_COUNT:
        raise RuntimeError(
            f"Expected {EXPECTED_FILE_COUNT} text files in {TXT_DIR}, found {len(txt_files)}"
        )
    if args.limit is not None:
        txt_files = txt_files[: args.limit]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    bootstrap_source_schema()
    configure_progression_helpers()
    mol_id, iso_id = parse_source_ids(SOURCE_TABLE)

    for txt_path in txt_files:
        band_label = band_label_from_path(txt_path)
        temp_table = temp_table_name(txt_path)
        out_csv = OUTPUT_DIR / f"{txt_path.stem}.csv"
        out_html = OUTPUT_DIR / f"{txt_path.stem}.html"

        try:
            line_count = build_temp_table_from_text(txt_path, temp_table)
            nu, absorbance = progression_plot.render_absorbance_table(mol_id, iso_id, temp_table)
            progression_plot.save_curve(
                out_csv=out_csv,
                out_html=out_html,
                x_values=nu,
                y_values=absorbance,
                trace_name=f"{band_label} (N={line_count})",
                title=(
                    f"CH4 {band_label} Absorbance (Voigt) "
                    f"{int(WN_MIN)}-{int(WN_MAX)} cm^-1 | T={T_K:g} K, P={P_TORR:g} Torr"
                    f" | x={MOLE_FRACTION:g}, L={PATH_LENGTH_CM:g} cm, S>{LINE_INTENSITY_THRESHOLD:.0e}"
                ),
            )
            print(f"Saved: {out_csv}")
            print(f"Saved: {out_html}")
        finally:
            hapi.dropTable(temp_table)


if __name__ == "__main__":
    main()
