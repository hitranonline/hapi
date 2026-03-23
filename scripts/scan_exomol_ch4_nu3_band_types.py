"""
Scan CH4 ExoMol MM transition files and count distinct nu3 transition types.

The script reads the ExoMol MM `.def`, `.states.bz2`, and overlapping
`.trans.bz2` files, derives transition wavenumbers from upper/lower state
energies, and classifies transitions by lower/upper vibrational labels.

It does not calculate line intensities, cross sections, or absorbance.

Outputs
-------
- A category CSV grouped only by nu3 quantum-number changes, such as `0->1`.
- An exact-pair CSV grouped by full ExoMol lower/upper vibrational+symmetry
  labels, plus HITRAN-style symmetry labels.
- A console summary with total transitions scanned, transitions kept, and the
  top exact pairs.

Default filter
--------------
By default, only upward pure-nu3 transitions are counted:
- `upper n3 > lower n3`
- `n1 = n2 = n4 = 0` for both lower and upper states

Use `--allow-other-modes` to keep upward `n3` transitions even when `n1`, `n2`,
or `n4` are non-zero.
"""

from __future__ import annotations

import argparse
import bz2
import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT_DIR / "exomol_db" / "CH4" / "12C-1H4" / "MM"
OUTPUT_DIR = ROOT_DIR / "exomol_ch4_mm_scans"
DATASET_STEM = "12C-1H4__MM"

EXOMOL_SYMMETRIES = ("A1", "A2", "E", "F1", "F2", "T1", "T2")
HITRAN_STYLE_MAP = {
    "A1": "1A1",
    "A2": "1A2",
    "E": "1E",
    "F1": "1F1",
    "F2": "1F2",
    "T1": "1F1",
    "T2": "1F2",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count distinct CH4 ExoMol MM nu3 band types in a wavenumber range.",
    )
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR, help="Directory containing ExoMol MM files.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Directory for summary CSV outputs.")
    parser.add_argument("--wn-min", type=float, default=2500.0, help="Minimum wavenumber in cm^-1.")
    parser.add_argument("--wn-max", type=float, default=3500.0, help="Maximum wavenumber in cm^-1.")
    parser.add_argument(
        "--allow-other-modes",
        action="store_true",
        help="Disable the pure-nu3 filter and allow non-zero n1, n2, or n4 states.",
    )
    parser.add_argument(
        "--print-top-exact",
        type=int,
        default=20,
        help="Print the top N exact nu3 band pairs by transition count.",
    )
    return parser.parse_args()


def iter_bz2_text_lines(path: Path, chunk_size: int = 8 * 1024 * 1024):
    decompressor = bz2.BZ2Decompressor()
    pending = b""

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break

            data = decompressor.decompress(chunk)
            if not data:
                continue

            pending += data
            lines = pending.split(b"\n")
            pending = lines.pop()
            for line in lines:
                yield line.decode("utf-8").rstrip("\r")

        if pending:
            yield pending.decode("utf-8").rstrip("\r")


def transition_filename(start_cm: int) -> str:
    return f"{DATASET_STEM}__{start_cm:05d}-{start_cm + 100:05d}.trans.bz2"


def available_transition_starts(data_dir: Path) -> list[int]:
    pattern = re.compile(rf"^{re.escape(DATASET_STEM)}__(\d{{5}})-(\d{{5}})\.trans\.bz2$")
    starts: list[int] = []
    for path in data_dir.glob(f"{DATASET_STEM}__*.trans.bz2"):
        match = pattern.match(path.name)
        if match is None:
            continue
        start_cm = int(match.group(1))
        end_cm = int(match.group(2))
        if end_cm == start_cm + 100:
            starts.append(start_cm)
    if not starts:
        raise FileNotFoundError(f"No transition files found in {data_dir}")
    return sorted(starts)


def overlapping_transition_files(data_dir: Path, wn_min: float, wn_max: float) -> list[Path]:
    available_starts = available_transition_starts(data_dir)
    available_min = available_starts[0]
    available_max = available_starts[-1] + 100

    start_chunk = int(math.floor(max(0.0, wn_min) / 100.0) * 100)
    end_chunk = int(math.ceil(wn_max / 100.0) * 100)

    clipped_start = max(start_chunk, available_min)
    clipped_end = min(end_chunk, available_max)
    if clipped_end <= clipped_start:
        raise FileNotFoundError(
            f"No available transition files overlap {wn_min:g}-{wn_max:g} cm^-1 in {data_dir}"
        )

    files: list[Path] = []
    for start_cm in range(clipped_start, clipped_end, 100):
        path = data_dir / transition_filename(start_cm)
        if not path.exists():
            raise FileNotFoundError(f"Missing transition file: {path}")
        files.append(path)
    return files


def parse_def_file(def_path: Path) -> dict[str, int]:
    metadata: dict[str, int] = {}
    with def_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            value_text, _, comment = raw_line.partition("#")
            value_text = value_text.strip()
            comment = comment.strip()
            if not value_text:
                continue
            if comment.startswith("No. of states in .states file"):
                metadata["nstates"] = int(value_text.split()[0])
    if "nstates" not in metadata:
        raise RuntimeError(f"Could not find state count in {def_path}")
    return metadata


def load_state_arrays(states_path: Path, nstates: int):
    energies = np.empty(nstates + 1, dtype=np.float64)
    n1 = np.empty(nstates + 1, dtype=np.int16)
    n2 = np.empty(nstates + 1, dtype=np.int16)
    n3 = np.empty(nstates + 1, dtype=np.int16)
    n4 = np.empty(nstates + 1, dtype=np.int16)
    symmetry = np.empty(nstates + 1, dtype=np.int8)

    energies.fill(np.nan)
    n1.fill(-1)
    n2.fill(-1)
    n3.fill(-1)
    n4.fill(-1)
    symmetry.fill(0)

    symmetry_to_code = {label: index + 1 for index, label in enumerate(EXOMOL_SYMMETRIES)}

    for line_number, raw_line in enumerate(iter_bz2_text_lines(states_path), start=1):
        parts = raw_line.split()
        if len(parts) < 19:
            continue

        state_id = int(parts[0])
        energies[state_id] = float(parts[1])
        symmetry[state_id] = symmetry_to_code.get(parts[6], 0)
        n1[state_id] = int(parts[9])
        n2[state_id] = int(parts[10])
        n3[state_id] = int(parts[12])
        n4[state_id] = int(parts[15])

        if line_number % 1_000_000 == 0:
            print(f"loaded {line_number:,} states")

    if np.isnan(energies[1:]).any():
        raise RuntimeError(f"State energies were not fully populated from {states_path}")

    return energies, n1, n2, n3, n4, symmetry


def symmetry_label(code: int) -> str:
    if code <= 0 or code > len(EXOMOL_SYMMETRIES):
        return "?"
    return EXOMOL_SYMMETRIES[code - 1]


def exomol_band_label(n1_value: int, n2_value: int, n3_value: int, n4_value: int, sym_code: int) -> str:
    return f"{n1_value} {n2_value} {n3_value} {n4_value} {symmetry_label(sym_code)}"


def hitran_style_band_label(n1_value: int, n2_value: int, n3_value: int, n4_value: int, sym_code: int) -> str:
    raw_symmetry = symmetry_label(sym_code)
    return f"{n1_value} {n2_value} {n3_value} {n4_value} {HITRAN_STYLE_MAP.get(raw_symmetry, raw_symmetry)}"


def output_stem(wn_min: float, wn_max: float, require_zero_other_modes: bool) -> str:
    mode_text = "pure_nu3_only" if require_zero_other_modes else "mixed_other_modes"
    return f"CH4_MM_nu3_band_types_{int(wn_min)}_{int(wn_max)}_{mode_text}"


def save_category_csv(path: Path, category_counts: Counter[tuple[int, int]], exact_to_category: dict[tuple[str, str], tuple[int, int]]) -> None:
    exact_count_per_category: defaultdict[tuple[int, int], int] = defaultdict(int)
    for category in exact_to_category.values():
        exact_count_per_category[category] += 1

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "category",
                "nu3_lower",
                "nu3_upper",
                "transition_count",
                "distinct_exact_pair_count",
            ]
        )
        for (lower_q, upper_q), count in sorted(category_counts.items()):
            writer.writerow(
                [
                    f"nu3 {lower_q}->{upper_q}",
                    lower_q,
                    upper_q,
                    count,
                    exact_count_per_category[(lower_q, upper_q)],
                ]
            )


def save_exact_csv(path: Path, exact_counts: Counter[tuple[str, str]], exact_to_hitran: dict[tuple[str, str], tuple[str, str]], exact_to_category: dict[tuple[str, str], tuple[int, int]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "category",
                "nu3_lower",
                "nu3_upper",
                "lower_label_exomol",
                "upper_label_exomol",
                "lower_label_hitran_style",
                "upper_label_hitran_style",
                "transition_count",
            ]
        )

        for (lower_label, upper_label), count in exact_counts.most_common():
            lower_hitran, upper_hitran = exact_to_hitran[(lower_label, upper_label)]
            lower_q, upper_q = exact_to_category[(lower_label, upper_label)]
            writer.writerow(
                [
                    f"nu3 {lower_q}->{upper_q}",
                    lower_q,
                    upper_q,
                    lower_label,
                    upper_label,
                    lower_hitran,
                    upper_hitran,
                    count,
                ]
            )


def validate_args(args: argparse.Namespace) -> None:
    if args.wn_max <= args.wn_min:
        raise ValueError("--wn-max must be greater than --wn-min")


def main() -> None:
    args = parse_args()
    validate_args(args)
    require_zero_other_modes = not args.allow_other_modes

    def_path = args.data_dir / f"{DATASET_STEM}.def"
    states_path = args.data_dir / f"{DATASET_STEM}.states.bz2"

    metadata = parse_def_file(def_path)
    print("loading states")
    energies, n1, n2, n3, n4, symmetry = load_state_arrays(states_path, int(metadata["nstates"]))

    transition_files = overlapping_transition_files(args.data_dir, args.wn_min, args.wn_max)
    print(f"scanning {len(transition_files)} transition files")

    category_counts: Counter[tuple[int, int]] = Counter()
    exact_counts: Counter[tuple[str, str]] = Counter()
    exact_to_hitran: dict[tuple[str, str], tuple[str, str]] = {}
    exact_to_category: dict[tuple[str, str], tuple[int, int]] = {}

    total_in_window = 0
    total_kept = 0

    for path in transition_files:
        print(f"scan {path.name}")
        for raw_line in iter_bz2_text_lines(path):
            parts = raw_line.split()
            if len(parts) != 3:
                continue

            upper_id = int(parts[0])
            lower_id = int(parts[1])
            wavenumber = energies[upper_id] - energies[lower_id]
            if wavenumber < args.wn_min or wavenumber > args.wn_max:
                continue
            total_in_window += 1

            lower_n3 = int(n3[lower_id])
            upper_n3 = int(n3[upper_id])
            if upper_n3 <= lower_n3:
                continue

            if require_zero_other_modes:
                if int(n1[lower_id]) != 0 or int(n1[upper_id]) != 0:
                    continue
                if int(n2[lower_id]) != 0 or int(n2[upper_id]) != 0:
                    continue
                if int(n4[lower_id]) != 0 or int(n4[upper_id]) != 0:
                    continue

            lower_label = exomol_band_label(int(n1[lower_id]), int(n2[lower_id]), lower_n3, int(n4[lower_id]), int(symmetry[lower_id]))
            upper_label = exomol_band_label(int(n1[upper_id]), int(n2[upper_id]), upper_n3, int(n4[upper_id]), int(symmetry[upper_id]))
            lower_hitran = hitran_style_band_label(int(n1[lower_id]), int(n2[lower_id]), lower_n3, int(n4[lower_id]), int(symmetry[lower_id]))
            upper_hitran = hitran_style_band_label(int(n1[upper_id]), int(n2[upper_id]), upper_n3, int(n4[upper_id]), int(symmetry[upper_id]))

            category = (lower_n3, upper_n3)
            exact_key = (lower_label, upper_label)

            category_counts[category] += 1
            exact_counts[exact_key] += 1
            exact_to_hitran[exact_key] = (lower_hitran, upper_hitran)
            exact_to_category[exact_key] = category
            total_kept += 1

            if total_in_window % 1_000_000 == 0:
                print(f"  in-window {total_in_window:,}, kept {total_kept:,}, exact pairs {len(exact_counts):,}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_stem(args.wn_min, args.wn_max, require_zero_other_modes)
    category_csv = args.output_dir / f"{stem}_categories.csv"
    exact_csv = args.output_dir / f"{stem}_exact_pairs.csv"

    save_category_csv(category_csv, category_counts, exact_to_category)
    save_exact_csv(exact_csv, exact_counts, exact_to_hitran, exact_to_category)

    print()
    print(f"window transitions found: {total_in_window:,}")
    print(f"nu3 transitions kept: {total_kept:,}")
    print(f"distinct exact nu3 band types: {len(exact_counts):,}")
    print()
    print("category summary:")
    for (lower_q, upper_q), count in sorted(category_counts.items()):
        print(f"  nu3 {lower_q}->{upper_q}: {count:,} transitions")

    if args.print_top_exact > 0:
        print()
        print("top exact nu3 band types:")
        for rank, ((lower_label, upper_label), count) in enumerate(exact_counts.most_common(args.print_top_exact), start=1):
            print(f"  {rank:2d}. {lower_label} -> {upper_label}: {count:,}")

    print()
    print(f"saved {category_csv}")
    print(f"saved {exact_csv}")


if __name__ == "__main__":
    main()
