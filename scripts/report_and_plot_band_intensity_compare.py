"""Generate an ExoMol-vs-HITRAN intensity comparison report and HTML plot.

Usage
-----
Run the script with one ExoMol band text file:

    python scripts/report_and_plot_band_intensity_compare.py ^
        --exomol-txt exomol_ch4_mm_pure_nu3_band_texts_hitran_style_2500_3500_sorted\\
EXOMOL_CH4_MM_pure_nu3_hitran_style_0_0_0_0_1F2_to_0_0_1_0_1F1_T296.0K.txt

What it does
------------
- Reads the ExoMol band file and extracts `wavenumber_cm-1` and
  `line_intensity_cm_per_molecule`.
- Searches `ch4_nu3_progressions/band_line_texts` for the exact matching
  HITRAN band label.
- If no exact HITRAN label exists, falls back to the nearest band with the
  same lower/upper vibrational quanta `n1 n2 n3 n4`.
- Writes both:
  - a markdown report describing the chosen comparison
  - an HTML Plotly overlay with two intensity curves

Outputs
-------
By default outputs are written to `band_intensity_comparisons/`.

The script now writes two HTML plots:
- one for the selected band comparison
- one residual plot for the selected band comparison
- one aggregate overlay combining all available `nu3 0->1` bands from both
  ExoMol and HITRAN inputs
- one aggregate residual plot for that `nu3 0->1` comparison
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_HITRAN_DIR = ROOT_DIR / "ch4_nu3_progressions" / "band_line_texts"
DEFAULT_HITRAN_HEADER = ROOT_DIR / "hitran_db" / "CH4_M6_I1.header"
DEFAULT_OUTPUT_DIR = ROOT_DIR / "band_intensity_comparisons"

EXOMOL_PREFIX = "EXOMOL_CH4_MM_pure_nu3_hitran_style_"
HITRAN_PREFIX = "CH4_M6_I1_"
LABEL_RE = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s*$")
FORMAT_RE = re.compile(r"%(\d+)(?:\.\d+)?[A-Za-z]")


@dataclass(frozen=True)
class BandLabel:
    lower: str
    upper: str
    lower_quanta: tuple[int, int, int, int]
    upper_quanta: tuple[int, int, int, int]
    lower_symmetry: str
    upper_symmetry: str


@dataclass(frozen=True)
class CandidateBand:
    path: Path
    label: BandLabel


@dataclass(frozen=True)
class SpectrumStats:
    line_count: int
    strongest_wavenumber: float
    strongest_intensity: float


@dataclass(frozen=True)
class AggregateStats:
    band_count: int
    spectrum: SpectrumStats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a band-comparison report and intensity overlay plot for ExoMol vs HITRAN CH4 exports.",
    )
    parser.add_argument("--exomol-txt", type=Path, required=True, help="Path to one ExoMol band text file.")
    parser.add_argument(
        "--hitran-dir",
        type=Path,
        default=DEFAULT_HITRAN_DIR,
        help="Directory containing HITRAN per-band text files.",
    )
    parser.add_argument(
        "--hitran-header",
        type=Path,
        default=DEFAULT_HITRAN_HEADER,
        help="Path to the CH4 HITRAN fixed-width header JSON.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for the generated report and HTML plot.",
    )
    return parser.parse_args()


def parse_fixed_width(format_string: str) -> int:
    match = FORMAT_RE.fullmatch(format_string)
    if match is None:
        raise ValueError(f"Unsupported fixed-width format: {format_string}")
    return int(match.group(1))


def parse_state_label(value: str) -> tuple[tuple[int, int, int, int], str]:
    match = LABEL_RE.fullmatch(value)
    if match is None:
        raise ValueError(f"Unexpected band-state label: {value!r}")
    quanta = tuple(int(match.group(i)) for i in range(1, 5))
    symmetry = match.group(5)
    return quanta, symmetry


def build_band_label(lower: str, upper: str) -> BandLabel:
    lower_quanta, lower_symmetry = parse_state_label(lower)
    upper_quanta, upper_symmetry = parse_state_label(upper)
    return BandLabel(
        lower=lower,
        upper=upper,
        lower_quanta=lower_quanta,
        upper_quanta=upper_quanta,
        lower_symmetry=lower_symmetry,
        upper_symmetry=upper_symmetry,
    )


def safe_label_fragment(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_") or "000"


def find_exomol_band_label(path: Path) -> BandLabel:
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith("# hitran_style_band_label="):
                payload = line.split("=", maxsplit=1)[1].strip()
                lower, upper = [part.strip() for part in payload.split("->", maxsplit=1)]
                return build_band_label(lower, upper)

    stem = path.stem
    if not stem.startswith(EXOMOL_PREFIX):
        raise ValueError(f"Could not infer ExoMol band label from file name: {path.name}")
    payload = stem[len(EXOMOL_PREFIX) :]
    payload = re.sub(r"_T\d+(?:\.\d+)?K$", "", payload)
    lower, upper = [part.replace("_", " ").strip() for part in payload.split("_to_", maxsplit=1)]
    return build_band_label(lower, upper)


def hitran_band_label_from_path(path: Path) -> BandLabel:
    stem = path.stem
    if not stem.startswith(HITRAN_PREFIX):
        raise ValueError(f"Unexpected HITRAN band file name: {path.name}")
    payload = stem[len(HITRAN_PREFIX) :]
    lower, upper = [part.replace("_", " ").strip() for part in payload.split("_to_", maxsplit=1)]
    return build_band_label(lower, upper)


def load_hitran_candidates(hitran_dir: Path) -> list[CandidateBand]:
    candidates: list[CandidateBand] = []
    for path in sorted(hitran_dir.glob("*.txt")):
        candidates.append(CandidateBand(path=path, label=hitran_band_label_from_path(path)))
    if not candidates:
        raise RuntimeError(f"No HITRAN band text files found in {hitran_dir}")
    return candidates


def load_exomol_candidates(exomol_dir: Path) -> list[CandidateBand]:
    candidates: list[CandidateBand] = []
    for path in sorted(exomol_dir.glob("*.txt")):
        candidates.append(CandidateBand(path=path, label=find_exomol_band_label(path)))
    if not candidates:
        raise RuntimeError(f"No ExoMol band text files found in {exomol_dir}")
    return candidates


def choose_hitran_candidate(exomol_label: BandLabel, candidates: list[CandidateBand]) -> tuple[CandidateBand, str]:
    for candidate in candidates:
        if candidate.label == exomol_label:
            return candidate, "exact_match"

    same_quanta = [
        candidate
        for candidate in candidates
        if candidate.label.lower_quanta == exomol_label.lower_quanta
        and candidate.label.upper_quanta == exomol_label.upper_quanta
    ]
    if not same_quanta:
        raise RuntimeError(
            "No HITRAN candidate shares the same lower/upper vibrational quanta as the ExoMol band label "
            f"{exomol_label.lower} -> {exomol_label.upper}"
        )

    def fallback_key(candidate: CandidateBand) -> tuple[int, int, str]:
        symmetry_matches = 0
        if candidate.label.lower_symmetry == exomol_label.lower_symmetry:
            symmetry_matches += 1
        if candidate.label.upper_symmetry == exomol_label.upper_symmetry:
            symmetry_matches += 1
        lower_swapped = int(candidate.label.lower_symmetry == exomol_label.upper_symmetry)
        upper_swapped = int(candidate.label.upper_symmetry == exomol_label.lower_symmetry)
        return (-symmetry_matches, -(lower_swapped + upper_swapped), candidate.path.name)

    return sorted(same_quanta, key=fallback_key)[0], "nearest_quanta_fallback"


def read_exomol_spectrum(path: Path) -> tuple[list[float], list[float]]:
    x_values: list[float] = []
    y_values: list[float] = []
    header_line: str | None = None

    with path.open("r", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if header_line is None:
                header_line = line
                columns = header_line.split("\t")
                try:
                    x_index = columns.index("wavenumber_cm-1")
                    y_index = columns.index("line_intensity_cm_per_molecule")
                except ValueError as exc:
                    raise RuntimeError(f"Unexpected ExoMol header in {path}") from exc
                continue

            parts = line.split("\t")
            x_values.append(float(parts[x_index]))
            y_values.append(float(parts[y_index]))

    if header_line is None or not x_values:
        raise RuntimeError(f"No ExoMol data rows found in {path}")
    return sort_spectrum(x_values, y_values)


def load_hitran_header(path: Path) -> tuple[dict[str, int], dict[str, int]]:
    with path.open("r", encoding="utf-8") as handle:
        header = json.load(handle)
    positions = header["position"]
    widths = {name: parse_fixed_width(fmt) for name, fmt in header["format"].items()}
    return positions, widths


def read_hitran_spectrum(path: Path, positions: dict[str, int], widths: dict[str, int]) -> tuple[list[float], list[float]]:
    x_values: list[float] = []
    y_values: list[float] = []

    nu_start = positions["nu"]
    nu_end = nu_start + widths["nu"]
    sw_start = positions["sw"]
    sw_end = sw_start + widths["sw"]

    with path.open("r", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line.strip():
                continue
            x_values.append(float(line[nu_start:nu_end]))
            y_values.append(float(line[sw_start:sw_end]))

    if not x_values:
        raise RuntimeError(f"No HITRAN data rows found in {path}")
    return sort_spectrum(x_values, y_values)


def sort_spectrum(x_values: list[float], y_values: list[float]) -> tuple[list[float], list[float]]:
    pairs = sorted(zip(x_values, y_values), key=lambda item: item[0])
    sorted_x = [item[0] for item in pairs]
    sorted_y = [item[1] for item in pairs]
    return sorted_x, sorted_y


def is_nu3_zero_to_one(label: BandLabel) -> bool:
    return label.lower_quanta[2] == 0 and label.upper_quanta[2] == 1


def aggregate_candidate_spectra(
    candidates: list[CandidateBand],
    reader,
    *args,
) -> tuple[list[float], list[float], AggregateStats]:
    x_values: list[float] = []
    y_values: list[float] = []
    used_bands = 0

    for candidate in candidates:
        if not is_nu3_zero_to_one(candidate.label):
            continue
        band_x, band_y = reader(candidate.path, *args)
        x_values.extend(band_x)
        y_values.extend(band_y)
        used_bands += 1

    if used_bands == 0:
        raise RuntimeError("No `nu3 0->1` bands found for aggregate comparison")

    sorted_x, sorted_y = sort_spectrum(x_values, y_values)
    return sorted_x, sorted_y, AggregateStats(
        band_count=used_bands,
        spectrum=spectrum_stats(sorted_x, sorted_y),
    )


def spectrum_stats(x_values: list[float], y_values: list[float]) -> SpectrumStats:
    strongest_index = max(range(len(y_values)), key=y_values.__getitem__)
    return SpectrumStats(
        line_count=len(y_values),
        strongest_wavenumber=x_values[strongest_index],
        strongest_intensity=y_values[strongest_index],
    )


def write_report(
    path: Path,
    *,
    requested_label: BandLabel,
    chosen_hitran_label: BandLabel,
    exomol_path: Path,
    hitran_path: Path,
    status: str,
    exomol_stats: SpectrumStats,
    hitran_stats: SpectrumStats,
    exomol_aggregate_stats: AggregateStats,
    hitran_aggregate_stats: AggregateStats,
) -> None:
    requested_band = f"{requested_label.lower} -> {requested_label.upper}"
    chosen_hitran_band = f"{chosen_hitran_label.lower} -> {chosen_hitran_label.upper}"
    exact_match = status == "exact_match"
    interpretation = (
        "Exact symmetry-resolved comparison. This is the strongest available check for the selected band."
        if exact_match
        else "Nearest same-quanta fallback. This supports band-level intensity validation, but it is not a symmetry-resolved one-to-one check."
    )

    strongest_wavenumber_delta = abs(exomol_stats.strongest_wavenumber - hitran_stats.strongest_wavenumber)
    if hitran_stats.strongest_intensity > 0.0:
        strongest_intensity_ratio = exomol_stats.strongest_intensity / hitran_stats.strongest_intensity
    else:
        strongest_intensity_ratio = float("inf")

    aggregate_strongest_wavenumber_delta = abs(
        exomol_aggregate_stats.spectrum.strongest_wavenumber - hitran_aggregate_stats.spectrum.strongest_wavenumber
    )
    if hitran_aggregate_stats.spectrum.strongest_intensity > 0.0:
        aggregate_strongest_intensity_ratio = (
            exomol_aggregate_stats.spectrum.strongest_intensity / hitran_aggregate_stats.spectrum.strongest_intensity
        )
    else:
        aggregate_strongest_intensity_ratio = float("inf")

    report = f"""# ExoMol vs HITRAN Intensity Comparison Report

## Executive Summary

- Requested ExoMol band: `{requested_band}`
- Chosen HITRAN band: `{chosen_hitran_band}`
- Match status: `{status}`
- Interpretation: {interpretation}

## Inputs

- ExoMol source: `{exomol_path}`
- HITRAN source: `{hitran_path}`

## Selected Band Comparison

### ExoMol

- Line count: {exomol_stats.line_count}
- Strongest line wavenumber: {exomol_stats.strongest_wavenumber:.6f} cm^-1
- Strongest line intensity: {exomol_stats.strongest_intensity:.6e}

### HITRAN

- Line count: {hitran_stats.line_count}
- Strongest line wavenumber: {hitran_stats.strongest_wavenumber:.6f} cm^-1
- Strongest line intensity: {hitran_stats.strongest_intensity:.6e}

### Band-Level Comparison Notes

- Strongest-line wavenumber difference: {strongest_wavenumber_delta:.6f} cm^-1
- Strongest-line intensity ratio ExoMol/HITRAN: {strongest_intensity_ratio:.6f}
- Plot file: `{path.with_name(path.stem.replace("_comparison_report", "_intensity_overlay") + ".html")}`
- Residual plot file: `{path.with_name(path.stem.replace("_comparison_report", "_intensity_residual") + ".html")}`

## Aggregate `nu3 0->1` Comparison

### ExoMol Aggregate

- Band count: {exomol_aggregate_stats.band_count}
- Total line count: {exomol_aggregate_stats.spectrum.line_count}
- Strongest line wavenumber: {exomol_aggregate_stats.spectrum.strongest_wavenumber:.6f} cm^-1
- Strongest line intensity: {exomol_aggregate_stats.spectrum.strongest_intensity:.6e}

### HITRAN Aggregate

- Band count: {hitran_aggregate_stats.band_count}
- Total line count: {hitran_aggregate_stats.spectrum.line_count}
- Strongest line wavenumber: {hitran_aggregate_stats.spectrum.strongest_wavenumber:.6f} cm^-1
- Strongest line intensity: {hitran_aggregate_stats.spectrum.strongest_intensity:.6e}

### Aggregate Comparison Notes

- Strongest-line wavenumber difference: {aggregate_strongest_wavenumber_delta:.6f} cm^-1
- Strongest-line intensity ratio ExoMol/HITRAN: {aggregate_strongest_intensity_ratio:.6f}
- Plot file: `{path.with_name("nu3_0_to_1_intensity_overlay.html")}`
- Residual plot file: `{path.with_name("nu3_0_to_1_intensity_residual.html")}`

## Assessment

{"- The selected comparison is an exact label match and can be used as a direct band-level validity check." if exact_match else "- The selected comparison uses the nearest same-quanta HITRAN band, so the result should be treated as a same-band sanity check rather than an exact label validation."}
- The aggregate `nu3 0->1` figure is the better high-level check of whether the ExoMol intensity distribution is broadly consistent with HITRAN in this spectral region.
"""
    path.write_text(report, encoding="utf-8")


def write_plot(
    path: Path,
    *,
    exomol_label: BandLabel,
    hitran_label: BandLabel,
    status: str,
    exomol_x: list[float],
    exomol_y: list[float],
    hitran_x: list[float],
    hitran_y: list[float],
) -> None:
    x_min = min(min(exomol_x), min(hitran_x))
    x_max = max(max(exomol_x), max(hitran_x))
    title = (
        "Band Intensity Comparison"
        "<br>Requested: "
        f"{exomol_label.lower} -> {exomol_label.upper}"
        "<br>Chosen HITRAN: "
        f"{hitran_label.lower} -> {hitran_label.upper}"
        "<br>Status: "
        f"{status}"
    )

    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Band Intensity Comparison</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    body {{
      margin: 0;
      font-family: Segoe UI, Arial, sans-serif;
      background: #ffffff;
      color: #111111;
    }}
    #plot {{
      width: 100%;
      height: 720px;
    }}
  </style>
</head>
<body>
  <div id="plot"></div>
  <script>
    const data = [
      {{
        x: {json.dumps(exomol_x)},
        y: {json.dumps(exomol_y)},
        mode: "lines",
        name: {json.dumps(f"ExoMol: {exomol_label.lower} -> {exomol_label.upper}")},
        line: {{ color: "#1f77b4", width: 1.5 }}
      }},
      {{
        x: {json.dumps(hitran_x)},
        y: {json.dumps(hitran_y)},
        mode: "lines",
        name: {json.dumps(f"HITRAN: {hitran_label.lower} -> {hitran_label.upper}")},
        line: {{ color: "#d62728", width: 1.5 }}
      }}
    ];

    const layout = {{
      title: {json.dumps(title)},
      template: "plotly_white",
      xaxis: {{
        title: "Wavenumber (cm^-1)",
        range: [{x_min}, {x_max}]
      }},
      yaxis: {{
        title: "Intensity"
      }},
      margin: {{
        l: 80,
        r: 30,
        t: 100,
        b: 70
      }}
    }};

    Plotly.newPlot("plot", data, layout, {{ responsive: true }});
  </script>
</body>
</html>
"""
    path.write_text(html_text, encoding="utf-8")


def residual_series(
    exomol_x: list[float],
    exomol_y: list[float],
    hitran_x: list[float],
    hitran_y: list[float],
) -> tuple[list[float], list[float]]:
    grid = np.array(sorted(set(exomol_x) | set(hitran_x)), dtype=float)
    exomol_interp = np.interp(grid, np.array(exomol_x, dtype=float), np.array(exomol_y, dtype=float), left=0.0, right=0.0)
    hitran_interp = np.interp(grid, np.array(hitran_x, dtype=float), np.array(hitran_y, dtype=float), left=0.0, right=0.0)
    residual = exomol_interp - hitran_interp
    return grid.tolist(), residual.tolist()


def write_residual_plot(
    path: Path,
    *,
    exomol_label: BandLabel,
    hitran_label: BandLabel,
    status: str,
    exomol_x: list[float],
    exomol_y: list[float],
    hitran_x: list[float],
    hitran_y: list[float],
) -> None:
    residual_x, residual_y = residual_series(exomol_x, exomol_y, hitran_x, hitran_y)
    x_min = min(residual_x)
    x_max = max(residual_x)
    title = (
        "Band Intensity Residual"
        "<br>Residual = ExoMol - HITRAN"
        "<br>Requested: "
        f"{exomol_label.lower} -> {exomol_label.upper}"
        "<br>Chosen HITRAN: "
        f"{hitran_label.lower} -> {hitran_label.upper}"
        "<br>Status: "
        f"{status}"
    )

    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Band Intensity Residual</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    body {{
      margin: 0;
      font-family: Segoe UI, Arial, sans-serif;
      background: #ffffff;
      color: #111111;
    }}
    #plot {{
      width: 100%;
      height: 720px;
    }}
  </style>
</head>
<body>
  <div id="plot"></div>
  <script>
    const data = [
      {{
        x: {json.dumps(residual_x)},
        y: {json.dumps(residual_y)},
        mode: "lines",
        name: "ExoMol - HITRAN",
        line: {{ color: "#2ca02c", width: 1.5 }}
      }},
      {{
        x: [{x_min}, {x_max}],
        y: [0, 0],
        mode: "lines",
        name: "Zero",
        line: {{ color: "#444444", width: 1, dash: "dash" }}
      }}
    ];

    const layout = {{
      title: {json.dumps(title)},
      template: "plotly_white",
      xaxis: {{
        title: "Wavenumber (cm^-1)",
        range: [{x_min}, {x_max}]
      }},
      yaxis: {{
        title: "Residual Intensity"
      }},
      margin: {{
        l: 80,
        r: 30,
        t: 110,
        b: 70
      }}
    }};

    Plotly.newPlot("plot", data, layout, {{ responsive: true }});
  </script>
</body>
</html>
"""
    path.write_text(html_text, encoding="utf-8")


def output_stem(exomol_label: BandLabel) -> str:
    label_text = f"{exomol_label.lower}_to_{exomol_label.upper}"
    return safe_label_fragment(label_text)


def main() -> None:
    args = parse_args()
    if not args.exomol_txt.exists():
        raise FileNotFoundError(f"Missing ExoMol text file: {args.exomol_txt}")
    if not args.hitran_dir.exists():
        raise FileNotFoundError(f"Missing HITRAN band-text directory: {args.hitran_dir}")
    if not args.hitran_header.exists():
        raise FileNotFoundError(f"Missing HITRAN header JSON: {args.hitran_header}")

    exomol_label = find_exomol_band_label(args.exomol_txt)
    exomol_candidates = load_exomol_candidates(args.exomol_txt.parent)
    hitran_candidates = load_hitran_candidates(args.hitran_dir)
    chosen_hitran, status = choose_hitran_candidate(exomol_label, hitran_candidates)

    exomol_x, exomol_y = read_exomol_spectrum(args.exomol_txt)
    positions, widths = load_hitran_header(args.hitran_header)
    hitran_x, hitran_y = read_hitran_spectrum(chosen_hitran.path, positions, widths)
    exomol_aggregate_x, exomol_aggregate_y, exomol_aggregate_stats = aggregate_candidate_spectra(
        exomol_candidates,
        read_exomol_spectrum,
    )
    hitran_aggregate_x, hitran_aggregate_y, hitran_aggregate_stats = aggregate_candidate_spectra(
        hitran_candidates,
        read_hitran_spectrum,
        positions,
        widths,
    )

    exomol_stats = spectrum_stats(exomol_x, exomol_y)
    hitran_stats = spectrum_stats(hitran_x, hitran_y)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    stem = output_stem(exomol_label)
    report_path = args.out_dir / f"{stem}_comparison_report.md"
    plot_path = args.out_dir / f"{stem}_intensity_overlay.html"
    residual_plot_path = args.out_dir / f"{stem}_intensity_residual.html"
    aggregate_plot_path = args.out_dir / "nu3_0_to_1_intensity_overlay.html"
    aggregate_residual_plot_path = args.out_dir / "nu3_0_to_1_intensity_residual.html"

    write_report(
        report_path,
        requested_label=exomol_label,
        chosen_hitran_label=chosen_hitran.label,
        exomol_path=args.exomol_txt.resolve(),
        hitran_path=chosen_hitran.path.resolve(),
        status=status,
        exomol_stats=exomol_stats,
        hitran_stats=hitran_stats,
        exomol_aggregate_stats=exomol_aggregate_stats,
        hitran_aggregate_stats=hitran_aggregate_stats,
    )
    write_plot(
        plot_path,
        exomol_label=exomol_label,
        hitran_label=chosen_hitran.label,
        status=status,
        exomol_x=exomol_x,
        exomol_y=exomol_y,
        hitran_x=hitran_x,
        hitran_y=hitran_y,
    )
    write_residual_plot(
        residual_plot_path,
        exomol_label=exomol_label,
        hitran_label=chosen_hitran.label,
        status=status,
        exomol_x=exomol_x,
        exomol_y=exomol_y,
        hitran_x=hitran_x,
        hitran_y=hitran_y,
    )
    write_plot(
        aggregate_plot_path,
        exomol_label=build_band_label("0 0 0 0 agg", "0 0 1 0 agg"),
        hitran_label=build_band_label("0 0 0 0 agg", "0 0 1 0 agg"),
        status="aggregate_nu3_0_to_1",
        exomol_x=exomol_aggregate_x,
        exomol_y=exomol_aggregate_y,
        hitran_x=hitran_aggregate_x,
        hitran_y=hitran_aggregate_y,
    )
    write_residual_plot(
        aggregate_residual_plot_path,
        exomol_label=build_band_label("0 0 0 0 agg", "0 0 1 0 agg"),
        hitran_label=build_band_label("0 0 0 0 agg", "0 0 1 0 agg"),
        status="aggregate_nu3_0_to_1",
        exomol_x=exomol_aggregate_x,
        exomol_y=exomol_aggregate_y,
        hitran_x=hitran_aggregate_x,
        hitran_y=hitran_aggregate_y,
    )

    print(f"status: {status}")
    print(f"report: {report_path}")
    print(f"plot: {plot_path}")
    print(f"residual_plot: {residual_plot_path}")
    print(f"aggregate_plot: {aggregate_plot_path}")
    print(f"aggregate_residual_plot: {aggregate_residual_plot_path}")


if __name__ == "__main__":
    main()
