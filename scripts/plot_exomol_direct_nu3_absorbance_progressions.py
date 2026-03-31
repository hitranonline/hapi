"""
Render nu3 absorbance progressions directly from the ExoMol MM database.

Unlike the sorted-band-text workflow (which uses pre-exported 296 K line
intensities), this script reads the raw ExoMol .trans.bz2 and .states.bz2
files and computes line intensities at the requested temperature using the
LTE formula.  Broadening uses the default gamma0 and n_exponent from the
ExoMol .def file.

The ExoMol partition function covers up to 5000 K, so this workflow can
run at temperatures beyond HAPI's TIPS2025 limit of 2500 K.

Each progression figure has 4 stacked panels:
- delta J = -1
- delta J = 0
- delta J = +1
- All delta J overlaid
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from _bootstrap import ensure_repo_root


ROOT_DIR = ensure_repo_root()

from research.exomol import dataset_dir, plot_exomol_direct_nu3_absorbance_progressions


J_PAIR_PATTERN = re.compile(r"^\s*(\d+)\s*->\s*(\d+)\s*$")


def parse_j_pair(text: str) -> tuple[int, int]:
    match = J_PAIR_PATTERN.fullmatch(text)
    if match is None:
        raise argparse.ArgumentTypeError("J pair must look like '2->3'")
    return int(match.group(1)), int(match.group(2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot nu3 absorbance progressions directly from the ExoMol MM database at arbitrary temperature.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=dataset_dir(),
        help="ExoMol MM dataset directory containing .def, .states.bz2, .pf, and .trans.bz2 files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for generated figures, CSV mappings, and report. Default includes temperature in the name.",
    )
    parser.add_argument("--wn-min", type=float, default=2500.0, help="Minimum wavenumber in cm^-1.")
    parser.add_argument("--wn-max", type=float, default=3500.0, help="Maximum wavenumber in cm^-1.")
    parser.add_argument("--wn-step", type=float, default=0.001, help="Wavenumber step in cm^-1.")
    parser.add_argument("--temperature-k", type=float, default=600.0, help="Gas temperature in K.")
    parser.add_argument("--pressure-torr", type=float, default=3.0, help="Total pressure in Torr.")
    parser.add_argument("--mole-fraction", type=float, default=0.008, help="CH4 mole fraction.")
    parser.add_argument("--path-length-cm", type=float, default=100.0, help="Optical path length in cm.")
    parser.add_argument("--line-cutoff", type=float, default=0.5, help="Voigt half-width cutoff in cm^-1.")
    parser.add_argument("--intensity-threshold", type=float, default=0.0, help="Discard lines weaker than this (cm/molecule).")
    parser.add_argument(
        "--label-top-n-per-delta-j",
        dest="label_top_n_per_delta_j",
        type=int,
        default=8,
        help="Number of strongest J-pair curves to label inside each delta-J panel.",
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
    output_dir = args.output_dir
    if output_dir is None:
        t = args.temperature_k
        p = args.pressure_torr
        x = args.mole_fraction
        l = args.path_length_cm
        step = args.wn_step
        slug = f"exomol_direct_nu3_absorbance_T{t:g}K_P{p:g}Torr_x{str(x).replace('.', 'p')}_L{l:g}cm_step{str(step).replace('.', 'p')}"
        output_dir = ROOT_DIR / "artifacts" / slug

    result = plot_exomol_direct_nu3_absorbance_progressions(
        data_dir=args.data_dir,
        output_dir=output_dir,
        wn_min=args.wn_min,
        wn_max=args.wn_max,
        wn_step=args.wn_step,
        temperature_k=args.temperature_k,
        pressure_torr=args.pressure_torr,
        mole_fraction=args.mole_fraction,
        path_length_cm=args.path_length_cm,
        line_cutoff=args.line_cutoff,
        intensity_threshold=args.intensity_threshold,
        label_top_n_per_delta_j=args.label_top_n_per_delta_j,
        html_max_points=args.html_max_points,
        forced_j_pairs=None if args.forced_j_pairs is None else tuple(args.forced_j_pairs),
    )
    print(f"progressions written: {len(result.rows)}")
    print(f"summary csv: {result.csv_path}")
    print(f"report: {result.metadata['report_path']}")
    for row in result.rows:
        print(
            f"{row['progression_label']}: "
            f"jpairs={row['jpair_count']} grid={row['grid_point_count']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
