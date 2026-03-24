"""
Render one intensity figure per ExoMol sorted pure-nu3 progression.

This is a thin wrapper around `research.exomol.plot_sorted_nu3_intensity_progressions`.
It scans the sorted ExoMol folder, merges the symmetry-resolved files for each
nu3 progression, groups rows by full J pair `(lower J, upper J)`, and writes:

- one PNG per progression
- one HTML figure per progression
- one J-pair mapping CSV per progression
- one report.md summarizing all generated figures
"""

from __future__ import annotations

import argparse
from pathlib import Path

from _bootstrap import ensure_repo_root


ROOT_DIR = ensure_repo_root()

from research.exomol import plot_sorted_nu3_intensity_progressions, sorted_nu3_band_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot ExoMol sorted pure-nu3 intensity progressions with one line per J pair.",
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
        default=ROOT_DIR / "artifacts" / "exomol_sorted_nu3_intensity",
        help="Directory for the generated figures, CSV mappings, and report.",
    )
    parser.add_argument("--wn-min", type=float, default=2500.0, help="Minimum wavenumber in cm^-1.")
    parser.add_argument("--wn-max", type=float, default=3500.0, help="Maximum wavenumber in cm^-1.")
    parser.add_argument(
        "--label-top-n",
        type=int,
        default=8,
        help="Number of strongest J-pair curves to label directly on each figure.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = plot_sorted_nu3_intensity_progressions(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        wn_min=args.wn_min,
        wn_max=args.wn_max,
        label_top_n=args.label_top_n,
    )

    report_path = Path(result.metadata["report_path"])
    print(f"progressions written: {result.metadata['progression_count']}")
    print(f"summary csv: {result.csv_path}")
    print(f"report: {report_path}")
    for row in result.rows:
        print(
            f"{row['progression_label']}: "
            f"files={row['source_file_count']} rows={row['row_count']} jpairs={row['jpair_count']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
