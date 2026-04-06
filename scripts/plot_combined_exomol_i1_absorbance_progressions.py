"""
Render one combined absorbance figure per pure-nu3 progression using:

- ExoMol MM I1 fixed-width band texts
- HITRAN fixed-width band texts

Both sources are rendered through the same HAPI temp-table absorbance flow.
The current combined script remains unchanged; this wrapper targets the
MM I1 ExoMol source specifically.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from _bootstrap import ensure_repo_root


ROOT_DIR = ensure_repo_root()

from research.combined import SOURCE_CHOICES, plot_combined_exomol_i1_absorbance_progressions
from research.hitran import band_line_text_dir


J_PAIR_PATTERN = re.compile(r"^\s*(\d+)\s*<-\s*(\d+)\s*$")


def parse_j_pair(text: str) -> tuple[int, int]:
    match = J_PAIR_PATTERN.fullmatch(text)
    if match is None:
        raise argparse.ArgumentTypeError("J pair must look like '2<-3'")
    return int(match.group(1)), int(match.group(2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot combined pure-nu3 absorbance progressions using ExoMol MM I1 and HITRAN through HAPI.",
    )
    parser.add_argument(
        "--exomol-input-dir",
        type=Path,
        default=ROOT_DIR / "exomol_ch4_mm_i1_pure_nu3_band_texts_hitran_style",
        help="Folder containing the ExoMol MM I1 fixed-width pure-nu3 band texts.",
    )
    parser.add_argument(
        "--exomol-header-path",
        type=Path,
        default=ROOT_DIR / "hitran_exomolCH4_db" / "CH4_EXOMOL_MM_I1.header",
        help="Header schema for the ExoMol MM I1 fixed-width rows.",
    )
    parser.add_argument(
        "--hitran-input-dir",
        type=Path,
        default=band_line_text_dir(ROOT_DIR),
        help="Folder containing the HITRAN band-line text exports.",
    )
    parser.add_argument(
        "--hitran-header-path",
        type=Path,
        default=ROOT_DIR / "hitran_db" / "CH4_M6_I1.header",
        help="Header schema for the HITRAN fixed-width rows.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT_DIR / "artifacts" / "combined_exomol_i1_absorbance_progressions",
        help="Directory for the generated figures, CSV mappings, and report.",
    )
    parser.add_argument(
        "--sources",
        nargs="*",
        default=list(SOURCE_CHOICES),
        choices=SOURCE_CHOICES,
        help="Subset of sources to load. Defaults to both exomol and hitran.",
    )
    parser.add_argument("--wn-min", type=float, default=2500.0, help="Minimum wavenumber in cm^-1.")
    parser.add_argument("--wn-max", type=float, default=3500.0, help="Maximum wavenumber in cm^-1.")
    parser.add_argument("--wn-step", type=float, default=0.001, help="Wavenumber step in cm^-1.")
    parser.add_argument("--temperature-k", type=float, default=600.0, help="Gas temperature in K.")
    parser.add_argument("--pressure-torr", type=float, default=3.0, help="Total pressure in Torr.")
    parser.add_argument("--mole-fraction", type=float, default=0.008, help="CH4 mole fraction.")
    parser.add_argument("--path-length-cm", type=float, default=100.0, help="Optical path length in cm.")
    parser.add_argument(
        "--intensity-threshold",
        type=float,
        default=0.0,
        help="HAPI intensity threshold used for ExoMol MM I1 temp-table rendering.",
    )
    parser.add_argument(
        "--hitran-intensity-threshold",
        type=float,
        default=1.0e-23,
        help="HAPI intensity threshold used for HITRAN temp-table rendering.",
    )
    parser.add_argument(
        "--label-top-n-per-delta-j",
        dest="label_top_n_per_delta_j",
        type=int,
        default=8,
        help="Number of strongest J-pair curves to label inside each delta-J panel.",
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
    parser.add_argument(
        "--save-raw-data",
        action="store_true",
        default=False,
        help="Save raw wavenumber + absorbance arrays per progression as .npz + CSV sidecar.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = plot_combined_exomol_i1_absorbance_progressions(
        exomol_input_dir=args.exomol_input_dir,
        exomol_header_path=args.exomol_header_path,
        hitran_input_dir=args.hitran_input_dir,
        hitran_header_path=args.hitran_header_path,
        output_dir=args.output_dir,
        sources=tuple(args.sources),
        wn_min=args.wn_min,
        wn_max=args.wn_max,
        wn_step=args.wn_step,
        temperature_k=args.temperature_k,
        pressure_torr=args.pressure_torr,
        mole_fraction=args.mole_fraction,
        path_length_cm=args.path_length_cm,
        intensity_threshold=args.intensity_threshold,
        hitran_intensity_threshold=args.hitran_intensity_threshold,
        label_top_n_per_delta_j=args.label_top_n_per_delta_j,
        html_max_points=args.html_max_points,
        forced_j_pairs=None if args.forced_j_pairs is None else tuple(args.forced_j_pairs),
        save_raw_data=args.save_raw_data,
    )
    print(f"progressions written: {len(result.rows)}")
    print(f"summary csv: {result.csv_path}")
    print(f"report: {result.metadata['report_path']}")
    for row in result.rows:
        print(
            f"{row['progression_label']}: "
            f"shared={row['shared_jpair_count']} "
            f"hitran_only={row['hitran_only_jpair_count']} "
            f"exomol_only={row['exomol_only_jpair_count']} "
            f"jpairs={row['jpair_count']} grid={row['grid_point_count']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
