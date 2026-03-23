"""Extract CH4 nu3 band rows from a HITRAN-style fixed-width database.

This script now works from any local HAPI/HITRAN-style `.header` + `.data`
pair that exposes:

- `nu`
- `global_lower_quanta`
- `global_upper_quanta`

It supports two workflows:

1. Summary-driven export
   Read an external CSV manifest containing exact `band_label` values and copy
   matching fixed-width rows into one `.txt` file per target band. This
   preserves the legacy HITEMP workflow used in this repo.

2. Auto-discovered pure-nu3 export
   Scan a HITRAN-style table, detect upward CH4 pure-`nu3` transitions directly
   from `global_*_quanta`, and write one `.txt` file per discovered band.
   This is intended for the local `CH4_EXOMOL_MM_I1` table.

All exported `.txt` files contain only the original fixed-width rows so they
remain compatible with the local HAPI text-import helpers.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, OrderedDict
from pathlib import Path
from typing import TextIO


ROOT_DIR = Path(__file__).resolve().parents[1]
MAX_OPEN_OUTPUTS = 128
FIELD_FORMAT_RE = re.compile(r"%(\d+)(?:\.\d+)?[A-Za-z]")

PRESETS = {
    "hitemp-ch4-nu3": {
        "source_header": ROOT_DIR / "hitemp_db" / "CH4" / "06_HITEMP2020.par_2500-3500.header",
        "source_data": ROOT_DIR / "hitemp_db" / "CH4" / "06_HITEMP2020.par_2500-3500.par",
        "summary_csv": ROOT_DIR / "ch4_nu3_progressions" / "CH4_nu3_progression_summary.csv",
        "output_dir": ROOT_DIR / "ch4_nu3_progressions" / "hitemp_band_line_texts",
        "selection_mode": "summary",
        "expected_band_count": 9,
        "wn_min": None,
        "wn_max": None,
        "assume_ordered": False,
    },
    "exomol-mm-i1-pure-nu3": {
        "source_header": ROOT_DIR / "hitran_db" / "CH4_EXOMOL_MM_I1.header",
        "source_data": ROOT_DIR / "hitran_db" / "CH4_EXOMOL_MM_I1.data",
        "summary_csv": None,
        "output_dir": ROOT_DIR / "exomol_ch4_mm_i1_pure_nu3_band_texts_hitran_style",
        "selection_mode": "pure-nu3",
        "expected_band_count": None,
        "wn_min": 2500.0,
        "wn_max": 3500.0,
        "assume_ordered": True,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract CH4 nu3 band rows from a HITRAN-style fixed-width table.",
    )
    parser.add_argument(
        "--preset",
        choices=tuple(PRESETS),
        default="exomol-mm-i1-pure-nu3",
        help="Repo-local configuration preset.",
    )
    parser.add_argument("--source-header", type=Path, default=None, help="Path to the HAPI/HITRAN-style header JSON.")
    parser.add_argument("--source-data", type=Path, default=None, help="Path to the fixed-width source data file.")
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=None,
        help="CSV containing exact `band_label` values. Required for --selection-mode=summary.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for exported per-band `.txt` files and the summary CSV.",
    )
    parser.add_argument(
        "--selection-mode",
        choices=("summary", "pure-nu3"),
        default=None,
        help="Select exact bands from a summary CSV or auto-discover pure-nu3 CH4 bands.",
    )
    parser.add_argument(
        "--expected-band-count",
        type=int,
        default=None,
        help="Optional validation for the number of unique summary-csv target bands.",
    )
    parser.add_argument("--wn-min", type=float, default=None, help="Optional minimum wavenumber in cm^-1.")
    parser.add_argument("--wn-max", type=float, default=None, help="Optional maximum wavenumber in cm^-1.")
    parser.add_argument(
        "--assume-ordered",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Assume the source file is ordered by `nu` and break after passing --wn-max.",
    )
    parser.add_argument(
        "--break-margin-cm",
        type=float,
        default=10.0,
        help="With --assume-ordered, stop after `nu > wn_max + break_margin_cm`.",
    )
    parser.add_argument(
        "--require-zero-other-modes",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="For pure-nu3 mode, require both states to satisfy n1=n2=n4=0.",
    )
    parser.add_argument(
        "--require-unit-step",
        action="store_true",
        default=False,
        help="For pure-nu3 mode, require unit-step progression in nu3, e.g. 0->1 or 1->2.",
    )
    parser.add_argument(
        "--max-written-lines",
        type=int,
        default=None,
        help="Optional hard stop after writing this many matched rows.",
    )
    parser.add_argument(
        "--stop-after-bands",
        type=int,
        default=None,
        help="Optional early stop after this many distinct exported bands have at least one row.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1_000_000,
        help="Print progress every N scanned input rows.",
    )
    return parser.parse_args()


def apply_preset(args: argparse.Namespace) -> argparse.Namespace:
    preset = PRESETS[args.preset]
    for key, value in preset.items():
        if getattr(args, key) is None:
            setattr(args, key, value)
    return args


def validate_args(args: argparse.Namespace) -> None:
    if args.source_header is None or args.source_data is None or args.output_dir is None:
        raise ValueError("source header, source data, and output directory must be resolved before running")
    if args.selection_mode == "summary" and args.summary_csv is None:
        raise ValueError("--summary-csv is required when --selection-mode=summary")
    if (args.wn_min is None) != (args.wn_max is None):
        raise ValueError("--wn-min and --wn-max must be provided together")
    if args.wn_min is not None and args.wn_max is not None and args.wn_max <= args.wn_min:
        raise ValueError("--wn-max must be greater than --wn-min")
    if args.break_margin_cm < 0.0:
        raise ValueError("--break-margin-cm must be non-negative")
    if args.max_written_lines is not None and args.max_written_lines <= 0:
        raise ValueError("--max-written-lines must be positive")
    if args.stop_after_bands is not None and args.stop_after_bands <= 0:
        raise ValueError("--stop-after-bands must be positive")
    if args.progress_every <= 0:
        raise ValueError("--progress-every must be positive")


def clean_quanta_label(value: str) -> str:
    text = value.strip()
    return text if text else "000"


def safe_label_fragment(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", clean_quanta_label(value)).strip("_") or "000"


def parse_fixed_width(format_string: str) -> int:
    match = FIELD_FORMAT_RE.fullmatch(format_string)
    if match is None:
        raise ValueError(f"Unsupported fixed-width format: {format_string}")
    return int(match.group(1))


def load_header_metadata(header_path: Path) -> tuple[dict[str, int], dict[str, int]]:
    with header_path.open("r", encoding="utf-8") as handle:
        header = json.load(handle)

    positions = header["position"]
    widths = {name: parse_fixed_width(fmt) for name, fmt in header["format"].items()}
    required_fields = {"nu", "global_lower_quanta", "global_upper_quanta"}
    missing = [name for name in required_fields if name not in positions or name not in widths]
    if missing:
        raise RuntimeError(f"Header {header_path} is missing required fields: {sorted(missing)}")
    return positions, widths


def load_target_bands(summary_csv: Path, expected_band_count: int | None) -> tuple[list[str], dict[str, int]]:
    target_labels: list[str] = []
    expected_counts: dict[str, int] = {}

    with summary_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if "band_label" not in reader.fieldnames:
            raise RuntimeError(f"Summary CSV {summary_csv} must contain a `band_label` column")
        has_line_count = reader.fieldnames is not None and "line_count" in reader.fieldnames
        for row in reader:
            band_label = str(row["band_label"]).strip()
            if not band_label:
                continue
            if band_label not in expected_counts:
                target_labels.append(band_label)
            expected_counts[band_label] = int(row["line_count"]) if has_line_count and row["line_count"].strip() else -1

    if expected_band_count is not None and len(target_labels) != expected_band_count:
        raise RuntimeError(
            f"Expected {expected_band_count} unique target bands in {summary_csv}, found {len(target_labels)}"
        )
    return target_labels, expected_counts


def extract_field(line: str, field_name: str, positions: dict[str, int], widths: dict[str, int]) -> str:
    start = positions[field_name]
    end = start + widths[field_name]
    return line[start:end]


def extract_band_label(line: str, positions: dict[str, int], widths: dict[str, int]) -> str:
    lower = clean_quanta_label(extract_field(line, "global_lower_quanta", positions, widths))
    upper = clean_quanta_label(extract_field(line, "global_upper_quanta", positions, widths))
    return f"{lower} -> {upper}"


def extract_wavenumber(line: str, positions: dict[str, int], widths: dict[str, int]) -> float:
    return float(extract_field(line, "nu", positions, widths).strip())


def parse_ch4_global_quanta(value: str) -> tuple[int, int, int, int, str]:
    parts = clean_quanta_label(value).split()
    if len(parts) != 5:
        raise ValueError(f"Unexpected CH4 global quanta label: {value!r}")
    return int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]), parts[4]


def pure_nu3_mode_pair(band_label: str, *, require_zero_other_modes: bool, require_unit_step: bool) -> tuple[int, int] | None:
    lower_raw, upper_raw = [part.strip() for part in band_label.split("->", maxsplit=1)]
    lower_n1, lower_n2, lower_n3, lower_n4, _ = parse_ch4_global_quanta(lower_raw)
    upper_n1, upper_n2, upper_n3, upper_n4, _ = parse_ch4_global_quanta(upper_raw)

    if upper_n3 <= lower_n3:
        return None
    if require_zero_other_modes:
        if any(value != 0 for value in (lower_n1, lower_n2, lower_n4, upper_n1, upper_n2, upper_n4)):
            return None
    if require_unit_step and upper_n3 != lower_n3 + 1:
        return None
    return lower_n3, upper_n3


def output_path_for_band(output_dir: Path, source_stem: str, band_label: str) -> Path:
    lower_label, upper_label = [part.strip() for part in band_label.split("->", maxsplit=1)]
    filename = f"{source_stem}_{safe_label_fragment(lower_label)}_to_{safe_label_fragment(upper_label)}.txt"
    return output_dir / filename


def output_summary_path(output_dir: Path, source_stem: str) -> Path:
    return output_dir / f"{source_stem}_band_text_summary.csv"


def mode_pair_label_from_band_label(band_label: str) -> str:
    lower_raw, upper_raw = [part.strip() for part in band_label.split("->", maxsplit=1)]
    lower_state = parse_ch4_global_quanta(lower_raw)
    upper_state = parse_ch4_global_quanta(upper_raw)
    return f"nu3 {lower_state[2]}->{upper_state[2]}"


def open_band_handle(
    band_label: str,
    *,
    band_paths: dict[str, Path],
    output_handles: OrderedDict[str, TextIO],
) -> TextIO:
    handle = output_handles.get(band_label)
    if handle is not None:
        output_handles.move_to_end(band_label)
        return handle

    if len(output_handles) >= MAX_OPEN_OUTPUTS:
        _, old_handle = output_handles.popitem(last=False)
        old_handle.close()

    handle = band_paths[band_label].open("a", encoding="utf-8", newline="")
    output_handles[band_label] = handle
    return handle


def write_summary_csv(
    path: Path,
    *,
    matched_counts: Counter[str],
    band_paths: dict[str, Path],
    expected_counts: dict[str, int],
) -> None:
    summary_rows = sorted(
        matched_counts.keys(),
        key=lambda label: (
            pure_nu3_mode_pair(label, require_zero_other_modes=False, require_unit_step=False) or (999, 999),
            label,
        ),
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "band_label",
                "mode_pair",
                "line_count",
                "reference_line_count",
                "txt_path",
            ],
        )
        writer.writeheader()
        for band_label in summary_rows:
            writer.writerow(
                {
                    "band_label": band_label,
                    "mode_pair": mode_pair_label_from_band_label(band_label),
                    "line_count": matched_counts[band_label],
                    "reference_line_count": "" if expected_counts[band_label] < 0 else expected_counts[band_label],
                    "txt_path": str(band_paths[band_label]),
                }
            )


def main() -> None:
    args = apply_preset(parse_args())
    validate_args(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    positions, widths = load_header_metadata(args.source_header)
    source_stem = args.source_data.stem

    target_bands: list[str]
    expected_counts: dict[str, int]
    if args.selection_mode == "summary":
        target_bands, expected_counts = load_target_bands(args.summary_csv, args.expected_band_count)
    else:
        target_bands = []
        expected_counts = {}

    band_paths: dict[str, Path] = {}
    matched_counts: Counter[str] = Counter()
    output_handles: OrderedDict[str, TextIO] = OrderedDict()
    exported_band_count = 0
    stats = {
        "scanned": 0,
        "in_window": 0,
        "band_matches": 0,
        "written": 0,
        "bad_quanta": 0,
    }

    if args.selection_mode == "summary":
        for band_label in target_bands:
            band_paths[band_label] = output_path_for_band(args.output_dir, source_stem, band_label)
            band_paths[band_label].write_text("", encoding="utf-8")

    try:
        with args.source_data.open("r", encoding="utf-8", newline="") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\r\n")
                if not line:
                    continue

                stats["scanned"] += 1
                if stats["scanned"] % args.progress_every == 0:
                    print(
                        f"scanned {stats['scanned']:,} rows, "
                        f"in-window {stats['in_window']:,}, "
                        f"matched {stats['band_matches']:,}, "
                        f"written {stats['written']:,}, "
                        f"bands {exported_band_count:,}"
                    )

                wavenumber = extract_wavenumber(line, positions, widths)
                if args.wn_min is not None and wavenumber < args.wn_min:
                    continue
                if args.wn_max is not None and wavenumber > args.wn_max:
                    if args.assume_ordered and wavenumber > args.wn_max + args.break_margin_cm:
                        print(
                            f"breaking early at nu={wavenumber:.6f} cm^-1 because "
                            f"--assume-ordered is enabled and the scan passed the selected window"
                        )
                        break
                    continue
                stats["in_window"] += 1

                band_label = extract_band_label(line, positions, widths)
                if args.selection_mode == "summary":
                    if band_label not in band_paths:
                        continue
                else:
                    try:
                        if pure_nu3_mode_pair(
                            band_label,
                            require_zero_other_modes=args.require_zero_other_modes,
                            require_unit_step=args.require_unit_step,
                        ) is None:
                            continue
                    except ValueError:
                        stats["bad_quanta"] += 1
                        continue
                    if band_label not in band_paths:
                        band_paths[band_label] = output_path_for_band(args.output_dir, source_stem, band_label)
                        expected_counts[band_label] = -1
                        target_bands.append(band_label)

                stats["band_matches"] += 1
                handle_out = open_band_handle(band_label, band_paths=band_paths, output_handles=output_handles)
                handle_out.write(line + "\n")
                if matched_counts[band_label] == 0:
                    exported_band_count += 1
                matched_counts[band_label] += 1
                stats["written"] += 1

                if args.max_written_lines is not None and stats["written"] >= args.max_written_lines:
                    print(f"stopped after writing {stats['written']:,} rows because --max-written-lines was reached")
                    break
                if args.stop_after_bands is not None and exported_band_count >= args.stop_after_bands:
                    print(f"stopped after exporting {exported_band_count:,} bands because --stop-after-bands was reached")
                    break
    finally:
        for handle in output_handles.values():
            handle.close()

    if args.selection_mode == "summary":
        for band_label in target_bands:
            if band_label not in matched_counts:
                matched_counts[band_label] = 0

    if not band_paths:
        raise RuntimeError("No output files were created for the current configuration.")
    if exported_band_count == 0:
        raise RuntimeError("No matching CH4 nu3 rows were found for the current filters.")

    summary_path = output_summary_path(args.output_dir, source_stem)
    write_summary_csv(
        summary_path,
        matched_counts=matched_counts,
        band_paths=band_paths,
        expected_counts=expected_counts,
    )

    print()
    print(f"source header: {args.source_header}")
    print(f"source data: {args.source_data}")
    print(f"selection mode: {args.selection_mode}")
    if args.wn_min is not None and args.wn_max is not None:
        print(f"wavenumber window: {args.wn_min:g} to {args.wn_max:g} cm^-1")
    print(f"scanned rows: {stats['scanned']:,}")
    print(f"in-window rows: {stats['in_window']:,}")
    print(f"matching band rows: {stats['band_matches']:,}")
    print(f"written rows: {stats['written']:,}")
    print(f"exported bands: {exported_band_count:,}")
    if stats["bad_quanta"] > 0:
        print(f"rows skipped due to unparseable CH4 global quanta: {stats['bad_quanta']:,}")
    if args.selection_mode == "summary":
        for band_label in target_bands:
            expected_count = expected_counts[band_label]
            if expected_count >= 0:
                print(
                    f"EXTRACTED: {band_label} | actual={matched_counts[band_label]} | "
                    f"reference={expected_count} | {band_paths[band_label]}"
                )
    print(f"saved {summary_path}")


if __name__ == "__main__":
    main()
