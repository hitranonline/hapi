"""Run the combined ExoMol I1 absorbance workflow at two temperatures
and produce dual-case comparison figures.

Usage::

    python scripts/plot_comparison_two_cases.py \\
        --temperature-a 300 --temperature-b 2500 --mole-fraction 0.1
"""

from __future__ import annotations

import argparse
from pathlib import Path

from _bootstrap import ensure_repo_root

ROOT_DIR = ensure_repo_root()

from research.combined import (
    SOURCE_CHOICES,
    plot_combined_exomol_i1_absorbance_progressions,
    PURE_NU3_PROGRESSIONS,
)
from research.compare_cases import CaseData, load_raw_case, plot_dual_case_comparison
from research.hitran import band_line_text_dir


def _progression_slug(lower_mode: tuple, upper_mode: tuple) -> str:
    return f"nu3_{lower_mode[2]}_to_{upper_mode[2]}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare ExoMol I1 absorbance at two temperatures side by side.",
    )
    parser.add_argument("--temperature-a", type=float, default=300.0, help="Temperature for case A (K).")
    parser.add_argument("--temperature-b", type=float, default=2500.0, help="Temperature for case B (K).")
    parser.add_argument("--mole-fraction", type=float, default=0.1, help="CH4 mole fraction.")
    parser.add_argument("--pressure-torr", type=float, default=3.0, help="Total pressure (Torr).")
    parser.add_argument("--path-length-cm", type=float, default=100.0, help="Optical path length (cm).")
    parser.add_argument("--wn-min", type=float, default=2500.0)
    parser.add_argument("--wn-max", type=float, default=3500.0)
    parser.add_argument("--wn-step", type=float, default=0.001)
    parser.add_argument("--color-a", type=str, default="#1f77b4", help="Hex color for case A.")
    parser.add_argument("--color-b", type=str, default="#ff7f0e", help="Hex color for case B.")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--exomol-input-dir", type=Path,
        default=ROOT_DIR / "exomol_ch4_mm_i1_pure_nu3_band_texts_hitran_style",
    )
    parser.add_argument(
        "--exomol-header-path", type=Path,
        default=ROOT_DIR / "hitran_exomolCH4_db" / "CH4_EXOMOL_MM_I1.header",
    )
    parser.add_argument(
        "--hitran-input-dir", type=Path,
        default=band_line_text_dir(ROOT_DIR),
    )
    parser.add_argument(
        "--hitran-header-path", type=Path,
        default=ROOT_DIR / "hitran_db" / "CH4_M6_I1.header",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    ta, tb = args.temperature_a, args.temperature_b
    output_dir = args.output_dir or (
        ROOT_DIR / "artifacts" / f"comparison_T{ta:.0f}K_vs_T{tb:.0f}K"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    common_kwargs = dict(
        exomol_input_dir=args.exomol_input_dir,
        exomol_header_path=args.exomol_header_path,
        hitran_input_dir=args.hitran_input_dir,
        hitran_header_path=args.hitran_header_path,
        sources=tuple(SOURCE_CHOICES),
        wn_min=args.wn_min,
        wn_max=args.wn_max,
        wn_step=args.wn_step,
        pressure_torr=args.pressure_torr,
        mole_fraction=args.mole_fraction,
        path_length_cm=args.path_length_cm,
        intensity_threshold=0.0,
        hitran_intensity_threshold=1e-23,
        save_raw_data=True,
    )

    case_dirs = {}
    for temp, label in [(ta, "a"), (tb, "b")]:
        case_dir = output_dir / f"case_{label}_T{temp:.0f}K"
        print(f"\n{'='*60}")
        print(f"Running case {label}: T = {temp} K")
        print(f"{'='*60}")
        plot_combined_exomol_i1_absorbance_progressions(
            **common_kwargs,
            temperature_k=temp,
            output_dir=case_dir,
        )
        case_dirs[label] = case_dir

    print(f"\n{'='*60}")
    print("Generating comparison figures")
    print(f"{'='*60}")

    png_paths = []
    for lower_mode, upper_mode in PURE_NU3_PROGRESSIONS:
        slug = _progression_slug(lower_mode, upper_mode)
        progression_label = f"nu3 {upper_mode[2]}<-{lower_mode[2]}"

        npz_a = case_dirs["a"] / f"{slug}_raw.npz"
        csv_a = case_dirs["a"] / f"{slug}_raw_meta.csv"
        npz_b = case_dirs["b"] / f"{slug}_raw.npz"
        csv_b = case_dirs["b"] / f"{slug}_raw_meta.csv"

        if not npz_a.exists() or not npz_b.exists():
            print(f"  Skipping {progression_label}: raw data not found")
            continue

        case_a = load_raw_case(
            npz_a, csv_a,
            label=f"T = {ta:.0f} K",
            color=args.color_a,
            temperature_k=ta,
        )
        case_b = load_raw_case(
            npz_b, csv_b,
            label=f"T = {tb:.0f} K",
            color=args.color_b,
            temperature_k=tb,
        )

        png = plot_dual_case_comparison(
            case_a, case_b,
            progression_label=progression_label,
            output_dir=output_dir,
            wn_min=args.wn_min,
            wn_max=args.wn_max,
        )
        png_paths.append(png)
        print(f"  {progression_label}: {png.name}")

    print(f"\nDone. {len(png_paths)} comparison figures in {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
