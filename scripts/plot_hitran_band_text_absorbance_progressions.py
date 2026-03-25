"""
Render absorbance progressions from HITRAN band-line text files.

This is a thin wrapper around `research.hitran.plot_band_text_absorbance_progressions`.
It groups the fixed-width `band_line_texts` rows by vibrational progression and
full J pair `(lower J, upper J)`, then renders one absorbance figure per
progression plus one J-pair CSV per figure.

Each progression figure is split into 3 stacked branch panels:
- `delta J = -1`
- `delta J = 0`
- `delta J = +1`

Inputs
------
- Band-line text files in `ch4_nu3_progressions/band_line_texts`
- HITRAN header/schema in `hitran_db/CH4_M6_I1.header`
- Case settings from CLI arguments:
  `--wn-min`, `--wn-max`, `--wn-step`,
  `--temperature-k`, `--pressure-torr`, `--mole-fraction`, `--path-length-cm`
- Optional progression filter:
  `--lower-mode` and `--upper-mode`

Outputs
-------
- One PNG per progression:
  `nu3_X_to_Y_absorbance.png`
- One HTML figure per progression:
  `nu3_X_to_Y_absorbance.html`
- One J-pair CSV per progression:
  `nu3_X_to_Y_jpairs.csv`
- One summary CSV for the whole run:
  `progression_summary.csv`
- One markdown report for the whole run:
  `report.md`

Default output folder
---------------------
- `artifacts/hitran_band_text_absorbance`
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from _bootstrap import ensure_repo_root


ROOT_DIR = ensure_repo_root()

from research.hitran import band_line_text_dir, plot_band_text_absorbance_progressions


J_PAIR_PATTERN = re.compile(r"^\s*(\d+)\s*->\s*(\d+)\s*$")


def parse_mode(text: str) -> tuple[int, int, int, int]:
    parts = text.split()
    if len(parts) != 4:
        raise argparse.ArgumentTypeError("mode must contain four integers, for example: '0 0 0 0'")
    try:
        return tuple(int(part) for part in parts)  # type: ignore[return-value]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("mode must contain only integers") from exc


def parse_j_pair(text: str) -> tuple[int, int]:
    match = J_PAIR_PATTERN.fullmatch(text)
    if match is None:
        raise argparse.ArgumentTypeError("J pair must look like '2->3'")
    return int(match.group(1)), int(match.group(2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot HITRAN band-line-text absorbance progressions with one line per J pair.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=band_line_text_dir(ROOT_DIR),
        help="Folder containing the HITRAN band-line text exports.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT_DIR / "artifacts" / "hitran_band_text_absorbance",
        help="Directory for the generated figures, J-pair CSVs, and report.",
    )
    parser.add_argument("--wn-min", type=float, default=2500.0, help="Minimum wavenumber in cm^-1.")
    parser.add_argument("--wn-max", type=float, default=3500.0, help="Maximum wavenumber in cm^-1.")
    parser.add_argument("--wn-step", type=float, default=0.01, help="Wavenumber step in cm^-1.")
    parser.add_argument("--temperature-k", type=float, default=600.0, help="Gas temperature in K.")
    parser.add_argument("--pressure-torr", type=float, default=3.0, help="Total pressure in Torr.")
    parser.add_argument("--mole-fraction", type=float, default=0.008, help="CH4 mole fraction.")
    parser.add_argument("--path-length-cm", type=float, default=100.0, help="Optical path length in cm.")
    parser.add_argument(
        "--lower-mode",
        type=parse_mode,
        default=None,
        help="Optional lower vibrational mode filter, for example '0 0 0 0'.",
    )
    parser.add_argument(
        "--upper-mode",
        type=parse_mode,
        default=None,
        help="Optional upper vibrational mode filter, for example '0 0 1 0'.",
    )
    parser.add_argument(
        "--intensity-threshold",
        type=float,
        default=1.0e-23,
        help="HAPI intensity threshold used while rendering each J-pair table.",
    )
    parser.add_argument(
        "--label-top-n-per-delta-j",
        dest="label_top_n_per_delta_j",
        type=int,
        default=8,
        help="Number of strongest J-pair curves to label inside each ΔJ panel.",
    )
    parser.add_argument(
        "--label-top-n",
        dest="label_top_n_per_delta_j",
        type=int,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--forced-j-pairs",
        type=parse_j_pair,
        nargs="*",
        default=None,
        help="Extra J pairs to force onto the figures, for example: 2->3 3->4 6->7.",
    )
    parser.add_argument(
        "--html-max-points",
        type=int,
        default=5000,
        help="Maximum points per J-pair trace in the HTML figures.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = plot_band_text_absorbance_progressions(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        wn_min=args.wn_min,
        wn_max=args.wn_max,
        wn_step=args.wn_step,
        temperature_k=args.temperature_k,
        pressure_torr=args.pressure_torr,
        mole_fraction=args.mole_fraction,
        path_length_cm=args.path_length_cm,
        intensity_threshold=args.intensity_threshold,
        label_top_n_per_delta_j=args.label_top_n_per_delta_j,
        html_max_points=args.html_max_points,
        lower_mode_filter=args.lower_mode,
        upper_mode_filter=args.upper_mode,
        forced_j_pairs=None if args.forced_j_pairs is None else tuple(args.forced_j_pairs),
    )
    print(f"progressions written: {result.count}")
    print(f"report: {result.manifest_path}")
    print(f"summary csv: {result.metadata['summary_csv_path']}")
    for row in result.rows:
        print(
            f"{row['progression_label']}: "
            f"files={row['source_file_count']} rows={row['row_count']} "
            f"jpairs={row['jpair_count']} grid={row['grid_point_count']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
