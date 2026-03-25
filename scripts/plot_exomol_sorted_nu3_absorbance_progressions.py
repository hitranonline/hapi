"""
Render one absorbance figure per ExoMol sorted pure-nu3 progression.

This is a thin wrapper around `research.exomol.plot_sorted_nu3_absorbance_progressions`.
It also runs the module validation step for `nu3 0->1` before building the
per-progression absorbance figures and final report.

Each progression figure is split into 3 stacked branch panels:
- `delta J = -1`
- `delta J = 0`
- `delta J = +1`

Notes
-----
- The sorted ExoMol exports in this workflow carry `T296.0K` line intensities.
- The absorbance workflow therefore keeps the line-strength temperature fixed at
  296 K and applies the requested pressure, mole fraction, path length, and
  broadening/grid settings on top of those exported line strengths.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from _bootstrap import ensure_repo_root


ROOT_DIR = ensure_repo_root()

from research.exomol import plot_sorted_nu3_absorbance_progressions, sorted_nu3_band_dir


J_PAIR_PATTERN = re.compile(r"^\s*(\d+)\s*->\s*(\d+)\s*$")


def parse_j_pair(text: str) -> tuple[int, int]:
    match = J_PAIR_PATTERN.fullmatch(text)
    if match is None:
        raise argparse.ArgumentTypeError("J pair must look like '2->3'")
    return int(match.group(1)), int(match.group(2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot ExoMol sorted pure-nu3 absorbance progressions with one line per J pair.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=sorted_nu3_band_dir(ROOT_DIR),
        help="Folder containing the sorted ExoMol pure-nu3 HITRAN-style band texts.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT_DIR / "artifacts" / "exomol_sorted_nu3_absorbance",
        help="Directory for the generated validation outputs, figures, CSV mappings, and report.",
    )
    parser.add_argument("--wn-min", type=float, default=2500.0, help="Minimum wavenumber in cm^-1.")
    parser.add_argument("--wn-max", type=float, default=3500.0, help="Maximum wavenumber in cm^-1.")
    parser.add_argument("--wn-step", type=float, default=0.25, help="Wavenumber step in cm^-1.")
    parser.add_argument("--pressure-torr", type=float, default=3.0, help="Total pressure in Torr.")
    parser.add_argument("--mole-fraction", type=float, default=0.008, help="CH4 mole fraction.")
    parser.add_argument("--path-length-cm", type=float, default=100.0, help="Optical path length in cm.")
    parser.add_argument(
        "--line-cutoff",
        type=float,
        default=0.5,
        help="Voigt half-width cutoff in cm^-1 used when rendering each exported line list.",
    )
    parser.add_argument(
        "--min-line-intensity",
        type=float,
        default=0.0,
        help="Discard exported lines weaker than this value in cm/molecule before rendering.",
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
    result = plot_sorted_nu3_absorbance_progressions(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        wn_min=args.wn_min,
        wn_max=args.wn_max,
        wn_step=args.wn_step,
        pressure_torr=args.pressure_torr,
        mole_fraction=args.mole_fraction,
        path_length_cm=args.path_length_cm,
        line_cutoff=args.line_cutoff,
        min_line_intensity=args.min_line_intensity,
        label_top_n_per_delta_j=args.label_top_n_per_delta_j,
        html_max_points=args.html_max_points,
        forced_j_pairs=None if args.forced_j_pairs is None else tuple(args.forced_j_pairs),
    )

    print(f"progressions written: {result.metadata['progression_count']}")
    print(f"summary csv: {result.csv_path}")
    print(f"report: {result.metadata['report_path']}")
    print(f"validation report: {result.metadata['validation_report_path']}")
    for row in result.rows:
        print(
            f"{row['progression_label']}: "
            f"files={row['source_file_count']} lines={row['row_count']} "
            f"jpairs={row['jpair_count']} grid={row['grid_point_count']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
