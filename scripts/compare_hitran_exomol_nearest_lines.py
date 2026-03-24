# Compare CH4 `nu3` lines from HITRAN and two ExoMol-derived exports.
#
# Implementation summary
# ----------------------
# - Sources:
#   - `exomol_ch4_mm_pure_nu3_band_texts_hitran_style_2500_3500_sorted`
#   - `exomol_ch4_mm_i1_pure_nu3_band_texts_hitran_style`
#   - `ch4_nu3_progressions/band_line_texts`
# - Default comparison window: `3000-3100 cm^-1`
# - Default vibrational mode: `0 0 0 0 -> 0 0 1 0`
# - Symmetry labels are ignored when grouping by vibrational mode.
# - Local labels are normalized to the first integer only, so rows are grouped
#   by `(lower J, upper J)`.
# - Because alpha/symmetry omission can merge multiple distinct rows into one J
#   pair, rows are not treated as one-to-one matches by J alone. Instead, each
#   ExoMol subset is matched to HITRAN by nearest wavenumber within a
#   tolerance.
# - Matching is done separately for:
#   - ExoMol sorted vs HITRAN
#   - ExoMol I1 vs HITRAN
# - Residuals are defined as:
#   - `ExoMol sorted intensity - HITRAN intensity`
#   - `ExoMol I1 intensity - HITRAN intensity`
# - Outputs:
#   - one markdown summary
#   - one HTML figure with a top intensity panel and a bottom residual panel

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt


ROOT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_EXOMOL_SORTED_DIR = ROOT_DIR / "exomol_ch4_mm_pure_nu3_band_texts_hitran_style_2500_3500_sorted"
DEFAULT_EXOMOL_I1_DIR = ROOT_DIR / "exomol_ch4_mm_i1_pure_nu3_band_texts_hitran_style"
DEFAULT_HITRAN_DIR = ROOT_DIR / "ch4_nu3_progressions" / "band_line_texts"
DEFAULT_HITRAN_HEADER = ROOT_DIR / "hitran_db" / "CH4_M6_I1.header"
DEFAULT_EXOMOL_I1_HEADER = ROOT_DIR / "hitran_db" / "CH4_EXOMOL_MM_I1.header"
DEFAULT_OUTPUT_DIR = ROOT_DIR / "band_intensity_comparisons"

EXOMOL_SORTED_PREFIX = "EXOMOL_CH4_MM_pure_nu3_hitran_style_"
EXOMOL_I1_PREFIX = "CH4_EXOMOL_MM_I1_"
HITRAN_PREFIX = "CH4_M6_I1_"

FIELD_FORMAT_RE = re.compile(r"%(\d+)(?:\.\d+)?[A-Za-z]")
FIRST_INTEGER_RE = re.compile(r"(\d+)")


@dataclass(frozen=True)
class LineRecord:
    dataset: str
    source_path: Path
    wavenumber: float
    intensity: float
    lower_mode: tuple[int, int, int, int]
    upper_mode: tuple[int, int, int, int]
    lower_j: int
    upper_j: int
    local_upper_label: str
    local_lower_label: str


@dataclass(frozen=True)
class MatchPair:
    exomol: LineRecord
    hitran: LineRecord
    delta_wavenumber: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare HITRAN and ExoMol CH4 lines by mode, J pair, and nearest wavenumber.",
    )
    parser.add_argument("--exomol-sorted-dir", type=Path, default=DEFAULT_EXOMOL_SORTED_DIR)
    parser.add_argument("--exomol-i1-dir", type=Path, default=DEFAULT_EXOMOL_I1_DIR)
    parser.add_argument("--hitran-dir", type=Path, default=DEFAULT_HITRAN_DIR)
    parser.add_argument("--hitran-header", type=Path, default=DEFAULT_HITRAN_HEADER)
    parser.add_argument("--lower-mode", type=str, default="0 0 0 0")
    parser.add_argument("--upper-mode", type=str, default="0 0 1 0")
    parser.add_argument("--wn-min", type=float, default=3000.0)
    parser.add_argument("--wn-max", type=float, default=3100.0)
    parser.add_argument("--lower-j", type=int, default=None)
    parser.add_argument("--upper-j", type=int, default=None)
    parser.add_argument("--match-tol", type=float, default=0.001)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--out-name", type=str, default=None)
    return parser.parse_args()


def parse_mode_text(value: str) -> tuple[int, int, int, int]:
    parts = value.split()
    if len(parts) != 4:
        raise ValueError(f"Expected four integers for vibrational mode, got: {value!r}")
    try:
        return tuple(int(part) for part in parts)  # type: ignore[return-value]
    except ValueError as exc:
        raise ValueError(f"Expected four integers for vibrational mode, got: {value!r}") from exc


def validate_args(args: argparse.Namespace) -> None:
    required_paths = [
        args.exomol_sorted_dir,
        args.exomol_i1_dir,
        args.hitran_dir,
        args.hitran_header,
        DEFAULT_EXOMOL_I1_HEADER,
    ]
    for path in required_paths:
        if not path.exists():
            raise FileNotFoundError(f"Missing required path: {path}")
    if args.wn_max <= args.wn_min:
        raise ValueError("--wn-max must be greater than --wn-min")
    if args.match_tol <= 0.0:
        raise ValueError("--match-tol must be positive")
    if (args.lower_j is None) != (args.upper_j is None):
        raise ValueError("--lower-j and --upper-j must be provided together")


def parse_fixed_width(format_string: str) -> int:
    match = FIELD_FORMAT_RE.fullmatch(format_string)
    if match is None:
        raise ValueError(f"Unsupported fixed-width format: {format_string}")
    return int(match.group(1))


def load_header_metadata(path: Path) -> tuple[dict[str, int], dict[str, int]]:
    with path.open("r", encoding="utf-8") as handle:
        header = json.load(handle)
    positions = header["position"]
    widths = {name: parse_fixed_width(fmt) for name, fmt in header["format"].items()}
    return positions, widths


def slice_field(line: str, positions: dict[str, int], widths: dict[str, int], field: str) -> str:
    start = positions[field]
    end = start + widths[field]
    return line[start:end].strip()


def parse_local_j(label: str) -> int:
    match = FIRST_INTEGER_RE.search(label.strip())
    if match is None:
        raise ValueError(f"Could not parse J from local label: {label!r}")
    return int(match.group(1))


def safe_label_fragment(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_") or "comparison"


def mode_pair_from_name(path: Path, prefix: str, *, strip_temperature: bool = False) -> tuple[tuple[int, int, int, int], tuple[int, int, int, int]]:
    stem = path.stem
    if not stem.startswith(prefix):
        raise ValueError(f"Unexpected file name for prefix {prefix!r}: {path.name}")
    payload = stem[len(prefix) :]
    if strip_temperature:
        payload = re.sub(r"_T\d+(?:\.\d+)?K$", "", payload)
    lower_raw, upper_raw = payload.split("_to_", maxsplit=1)
    lower_mode = tuple(int(part) for part in lower_raw.split("_")[:4])
    upper_mode = tuple(int(part) for part in upper_raw.split("_")[:4])
    return lower_mode, upper_mode


def read_exomol_sorted_rows(
    input_dir: Path,
    *,
    target_lower_mode: tuple[int, int, int, int],
    target_upper_mode: tuple[int, int, int, int],
    wn_min: float,
    wn_max: float,
) -> dict[tuple[int, int], list[LineRecord]]:
    grouped: dict[tuple[int, int], list[LineRecord]] = {}
    txt_files = sorted(input_dir.glob("*.txt"))
    if not txt_files:
        raise RuntimeError(f"No ExoMol sorted text files found in {input_dir}")

    for path in txt_files:
        lower_mode, upper_mode = mode_pair_from_name(path, EXOMOL_SORTED_PREFIX, strip_temperature=True)
        if lower_mode != target_lower_mode or upper_mode != target_upper_mode:
            continue

        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader((line for line in handle if line.strip() and not line.startswith("#")), delimiter="\t")
            for row in reader:
                wavenumber = float(row["wavenumber_cm-1"])
                if wavenumber < wn_min or wavenumber > wn_max:
                    continue
                intensity = float(row["line_intensity_cm_per_molecule"])
                local_upper = row["local_upper_quanta"].strip()
                local_lower = row["local_lower_quanta"].strip()
                lower_j = parse_local_j(local_lower)
                upper_j = parse_local_j(local_upper)
                key = (lower_j, upper_j)
                grouped.setdefault(key, []).append(
                    LineRecord(
                        dataset="ExoMol Sorted",
                        source_path=path,
                        wavenumber=wavenumber,
                        intensity=intensity,
                        lower_mode=lower_mode,
                        upper_mode=upper_mode,
                        lower_j=lower_j,
                        upper_j=upper_j,
                        local_upper_label=local_upper,
                        local_lower_label=local_lower,
                    )
                )

    return sort_grouped_rows(grouped)


def read_fixed_width_rows(
    input_dir: Path,
    *,
    header_path: Path,
    dataset_name: str,
    file_prefix: str,
    target_lower_mode: tuple[int, int, int, int],
    target_upper_mode: tuple[int, int, int, int],
    wn_min: float,
    wn_max: float,
) -> dict[tuple[int, int], list[LineRecord]]:
    grouped: dict[tuple[int, int], list[LineRecord]] = {}
    positions, widths = load_header_metadata(header_path)
    txt_files = sorted(input_dir.glob("*.txt"))
    if not txt_files:
        raise RuntimeError(f"No fixed-width text files found in {input_dir}")

    for path in txt_files:
        lower_mode, upper_mode = mode_pair_from_name(path, file_prefix)
        if lower_mode != target_lower_mode or upper_mode != target_upper_mode:
            continue

        with path.open("r", encoding="utf-8", newline="") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\r\n")
                if not line.strip():
                    continue
                wavenumber = float(slice_field(line, positions, widths, "nu"))
                if wavenumber < wn_min or wavenumber > wn_max:
                    continue
                intensity = float(slice_field(line, positions, widths, "sw"))
                local_upper = slice_field(line, positions, widths, "local_upper_quanta")
                local_lower = slice_field(line, positions, widths, "local_lower_quanta")
                lower_j = parse_local_j(local_lower)
                upper_j = parse_local_j(local_upper)
                key = (lower_j, upper_j)
                grouped.setdefault(key, []).append(
                    LineRecord(
                        dataset=dataset_name,
                        source_path=path,
                        wavenumber=wavenumber,
                        intensity=intensity,
                        lower_mode=lower_mode,
                        upper_mode=upper_mode,
                        lower_j=lower_j,
                        upper_j=upper_j,
                        local_upper_label=local_upper,
                        local_lower_label=local_lower,
                    )
                )

    return sort_grouped_rows(grouped)


def sort_grouped_rows(grouped: dict[tuple[int, int], list[LineRecord]]) -> dict[tuple[int, int], list[LineRecord]]:
    for rows in grouped.values():
        rows.sort(key=lambda row: row.wavenumber)
    return grouped


def select_j_pair(
    sorted_rows: dict[tuple[int, int], list[LineRecord]],
    i1_rows: dict[tuple[int, int], list[LineRecord]],
    hitran_rows: dict[tuple[int, int], list[LineRecord]],
    *,
    lower_j: int | None,
    upper_j: int | None,
) -> tuple[int, int]:
    common_pairs = set(sorted_rows) & set(i1_rows) & set(hitran_rows)
    if not common_pairs:
        raise RuntimeError("No common J pair exists across all three datasets in the requested mode/window.")

    if lower_j is not None and upper_j is not None:
        requested = (lower_j, upper_j)
        if requested not in common_pairs:
            raise RuntimeError(
                f"Requested J pair {requested} is not available across all three datasets in the requested mode/window."
            )
        return requested

    def pair_score(pair: tuple[int, int]) -> tuple[float, float, int, int]:
        totals = [
            sum(row.intensity for row in sorted_rows[pair]),
            sum(row.intensity for row in i1_rows[pair]),
            sum(row.intensity for row in hitran_rows[pair]),
        ]
        return (
            min(totals),
            sum(totals) / len(totals),
            -pair[0],
            -pair[1],
        )

    return max(common_pairs, key=pair_score)


def nearest_observed_delta(exomol_rows: list[LineRecord], hitran_rows: list[LineRecord]) -> float:
    return min(abs(exomol.wavenumber - hitran.wavenumber) for exomol in exomol_rows for hitran in hitran_rows)


def match_rows(exomol_rows: list[LineRecord], hitran_rows: list[LineRecord], tolerance: float) -> list[MatchPair]:
    candidates: list[tuple[float, int, int]] = []
    for ex_index, ex_row in enumerate(exomol_rows):
        for hit_index, hit_row in enumerate(hitran_rows):
            delta = abs(ex_row.wavenumber - hit_row.wavenumber)
            if delta <= tolerance:
                candidates.append((delta, ex_index, hit_index))

    candidates.sort(key=lambda item: (item[0], exomol_rows[item[1]].wavenumber, hitran_rows[item[2]].wavenumber))
    used_exomol: set[int] = set()
    used_hitran: set[int] = set()
    matches: list[MatchPair] = []
    for delta, ex_index, hit_index in candidates:
        if ex_index in used_exomol or hit_index in used_hitran:
            continue
        used_exomol.add(ex_index)
        used_hitran.add(hit_index)
        matches.append(
            MatchPair(
                exomol=exomol_rows[ex_index],
                hitran=hitran_rows[hit_index],
                delta_wavenumber=delta,
            )
        )

    matches.sort(key=lambda pair: pair.hitran.wavenumber)
    return matches


def relative_difference(exomol_intensity: float, hitran_intensity: float) -> float:
    if hitran_intensity == 0.0:
        return math.inf if exomol_intensity != 0.0 else 0.0
    return abs(exomol_intensity - hitran_intensity) / abs(hitran_intensity)


def summarize_matches(matches: list[MatchPair]) -> dict[str, float]:
    deltas = [match.delta_wavenumber for match in matches]
    relative_diffs = [relative_difference(match.exomol.intensity, match.hitran.intensity) for match in matches]
    finite_relative_diffs = [value for value in relative_diffs if math.isfinite(value)]
    return {
        "count": float(len(matches)),
        "max_delta_wavenumber": max(deltas),
        "mean_delta_wavenumber": sum(deltas) / len(deltas),
        "max_relative_intensity_difference": max(relative_diffs),
        "mean_relative_intensity_difference": (
            sum(finite_relative_diffs) / len(finite_relative_diffs) if finite_relative_diffs else math.inf
        ),
    }


def output_stem(
    lower_mode: tuple[int, int, int, int],
    upper_mode: tuple[int, int, int, int],
    j_pair: tuple[int, int],
    out_name: str | None,
) -> str:
    if out_name:
        return safe_label_fragment(out_name)
    lower_text = "_".join(str(value) for value in lower_mode)
    upper_text = "_".join(str(value) for value in upper_mode)
    return f"compare_{lower_text}_to_{upper_text}_J{j_pair[0]}_{j_pair[1]}"


def hitran_union_points(sorted_matches: list[MatchPair], i1_matches: list[MatchPair]) -> tuple[list[float], list[float]]:
    union: dict[tuple[float, float], None] = {}
    for match in sorted_matches + i1_matches:
        union[(match.hitran.wavenumber, match.hitran.intensity)] = None
    ordered = sorted(union)
    return [item[0] for item in ordered], [item[1] for item in ordered]


def reference_percent_bound(residual_bound: float, reference_intensity: float) -> float:
    if reference_intensity <= 0.0:
        return 100.0
    return 100.0 * residual_bound / reference_intensity


def build_plot_html(
    *,
    lower_mode: tuple[int, int, int, int],
    upper_mode: tuple[int, int, int, int],
    j_pair: tuple[int, int],
    sorted_matches: list[MatchPair],
    i1_matches: list[MatchPair],
) -> str:
    sorted_x = [match.hitran.wavenumber for match in sorted_matches]
    sorted_y = [match.exomol.intensity for match in sorted_matches]
    sorted_residual = [match.exomol.intensity - match.hitran.intensity for match in sorted_matches]

    i1_x = [match.hitran.wavenumber for match in i1_matches]
    i1_y = [match.exomol.intensity for match in i1_matches]
    i1_residual = [match.exomol.intensity - match.hitran.intensity for match in i1_matches]

    hitran_x, hitran_y = hitran_union_points(sorted_matches, i1_matches)
    x_min = min(hitran_x)
    x_max = max(hitran_x)
    residual_bound = max(abs(value) for value in sorted_residual + i1_residual)
    reference_intensity = max(abs(value) for value in hitran_y)
    percent_bound = reference_percent_bound(residual_bound * 1.05, reference_intensity)

    title = (
        "HITRAN vs ExoMol CH4 Line Comparison"
        f"<br>Mode: {' '.join(map(str, lower_mode))} -> {' '.join(map(str, upper_mode))}"
        f"<br>J pair: ({j_pair[0]}, {j_pair[1]})"
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>HITRAN vs ExoMol CH4 Comparison</title>
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
      height: 900px;
    }}
  </style>
</head>
<body>
  <div id="plot"></div>
  <script>
    const data = [
      {{
        x: {json.dumps(hitran_x)},
        y: {json.dumps(hitran_y)},
        mode: "markers",
        name: "HITRAN",
        marker: {{ size: 5, symbol: "x", color: "#111111" }}
      }},
      {{
        x: {json.dumps(i1_x)},
        y: {json.dumps(i1_y)},
        mode: "markers",
        name: "ExoMol I1",
        marker: {{ size: 7, symbol: "square", color: "#d62728" }}
      }},
      {{
        x: {json.dumps(sorted_x)},
        y: {json.dumps(sorted_y)},
        mode: "markers",
        name: "ExoMol Sorted",
        marker: {{ size: 11, symbol: "circle-open", color: "#1f77b4", line: {{ width: 2 }} }}
      }},
      {{
        x: {json.dumps(sorted_x)},
        y: {json.dumps(sorted_residual)},
        mode: "markers",
        name: "Residual: Sorted - HITRAN",
        marker: {{ size: 6, symbol: "circle", color: "#2ca02c" }},
        xaxis: "x2",
        yaxis: "y2"
      }},
      {{
        x: {json.dumps(i1_x)},
        y: {json.dumps(i1_residual)},
        mode: "markers",
        name: "Residual: I1 - HITRAN",
        marker: {{ size: 7, symbol: "square", color: "#ff7f0e" }},
        xaxis: "x2",
        yaxis: "y2"
      }},
      {{
        x: [{x_min}, {x_max}],
        y: [0, 0],
        mode: "lines",
        name: "Zero (Residual)",
        line: {{ color: "#666666", width: 1.0, dash: "dash" }},
        xaxis: "x2",
        yaxis: "y2"
      }}
    ];

    const layout = {{
      title: {json.dumps(title)},
      template: "plotly_white",
      height: 900,
      margin: {{ l: 80, r: 90, t: 110, b: 70 }},
      legend: {{ orientation: "h", x: 0, y: 1.10 }},
      xaxis: {{
        domain: [0, 1],
        anchor: "y",
        range: [{x_min}, {x_max}],
        title: "HITRAN Wavenumber (cm^-1)"
      }},
      yaxis: {{
        domain: [0.46, 1],
        title: "Intensity"
      }},
      xaxis2: {{
        domain: [0, 1],
        anchor: "y2",
        range: [{x_min}, {x_max}],
        title: "HITRAN Wavenumber (cm^-1)"
      }},
      yaxis2: {{
        domain: [0, 0.28],
        title: "Residual Intensity",
        range: [{-residual_bound * 1.05}, {residual_bound * 1.05}]
      }},
      yaxis4: {{
        overlaying: "y2",
        anchor: "x2",
        side: "right",
        title: "Residual / peak HITRAN intensity (%)",
        range: [{-percent_bound}, {percent_bound}]
      }},
      annotations: [
        {{
          text: "Matched Intensities",
          xref: "paper",
          yref: "paper",
          x: 0.5,
          y: 1.02,
          showarrow: false,
          font: {{ size: 14 }}
        }},
        {{
          text: "Residuals with shared percent right axis",
          xref: "paper",
          yref: "paper",
          x: 0.5,
          y: 0.33,
          showarrow: false,
          font: {{ size: 14 }}
        }}
      ]
    }};

    Plotly.newPlot("plot", data, layout, {{ responsive: true }});
  </script>
</body>
</html>
"""


def write_plot_png(
    path: Path,
    *,
    lower_mode: tuple[int, int, int, int],
    upper_mode: tuple[int, int, int, int],
    j_pair: tuple[int, int],
    sorted_matches: list[MatchPair],
    i1_matches: list[MatchPair],
) -> None:
    sorted_x = [match.hitran.wavenumber for match in sorted_matches]
    sorted_y = [match.exomol.intensity for match in sorted_matches]
    sorted_residual = [match.exomol.intensity - match.hitran.intensity for match in sorted_matches]

    i1_x = [match.hitran.wavenumber for match in i1_matches]
    i1_y = [match.exomol.intensity for match in i1_matches]
    i1_residual = [match.exomol.intensity - match.hitran.intensity for match in i1_matches]

    hitran_x, hitran_y = hitran_union_points(sorted_matches, i1_matches)
    residual_bound = max(abs(value) for value in sorted_residual + i1_residual)
    reference_intensity = max(abs(value) for value in hitran_y)
    percent_bound = reference_percent_bound(residual_bound * 1.05, reference_intensity)

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        figsize=(12, 9),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.15]},
    )

    ax_top.scatter(hitran_x, hitran_y, color="#111111", marker="x", s=30, label="HITRAN", zorder=2)
    ax_top.scatter(i1_x, i1_y, color="#d62728", marker="s", s=36, label="ExoMol I1", zorder=3)
    ax_top.scatter(sorted_x, sorted_y, facecolors="none", edgecolors="#1f77b4", marker="o", s=90, linewidths=1.8, label="ExoMol Sorted", zorder=4)
    ax_top.set_ylabel("Intensity")
    ax_top.legend(loc="best")
    ax_top.grid(True, alpha=0.25)
    ax_top.set_title(
        "HITRAN vs ExoMol CH4 Line Comparison\n"
        f"Mode: {' '.join(map(str, lower_mode))} -> {' '.join(map(str, upper_mode))} | "
        f"J pair: ({j_pair[0]}, {j_pair[1]})"
    )

    ax_bottom.scatter(sorted_x, sorted_residual, color="#2ca02c", marker="o", s=36, label="Residual: Sorted - HITRAN")
    ax_bottom.scatter(i1_x, i1_residual, color="#ff7f0e", marker="s", s=36, label="Residual: I1 - HITRAN")
    ax_bottom.axhline(0.0, color="#666666", linewidth=1.0, linestyle="--")
    ax_bottom.set_xlabel("HITRAN Wavenumber (cm^-1)")
    ax_bottom.set_ylabel("Residual Intensity")
    ax_bottom.set_ylim(-residual_bound * 1.05, residual_bound * 1.05)
    ax_bottom.grid(True, alpha=0.25)
    ax_bottom.legend(loc="upper left")

    ax_bottom_right = ax_bottom.twinx()
    ax_bottom_right.set_ylabel("Residual / peak HITRAN intensity (%)")
    ax_bottom_right.set_ylim(-percent_bound, percent_bound)

    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def build_summary_text(
    *,
    lower_mode: tuple[int, int, int, int],
    upper_mode: tuple[int, int, int, int],
    j_pair: tuple[int, int],
    wn_min: float,
    wn_max: float,
    match_tol: float,
    sorted_candidates: list[LineRecord],
    i1_candidates: list[LineRecord],
    hitran_candidates: list[LineRecord],
    sorted_matches: list[MatchPair],
    i1_matches: list[MatchPair],
    plot_path: Path,
    figure_path: Path,
) -> str:
    sorted_stats = summarize_matches(sorted_matches)
    i1_stats = summarize_matches(i1_matches)
    shared_hitran_points = len(set(match.hitran.wavenumber for match in sorted_matches) & set(match.hitran.wavenumber for match in i1_matches))
    figure_ref = figure_path.name
    return f"""# HITRAN vs ExoMol CH4 Nearest-Line Comparison

## Selection

- Mode: `{' '.join(map(str, lower_mode))} -> {' '.join(map(str, upper_mode))}`
- J pair: `({j_pair[0]}, {j_pair[1]})`
- Wavenumber window: `{wn_min:.1f} to {wn_max:.1f} cm^-1`
- Nearest-match tolerance: `{match_tol:.6f} cm^-1`

## Candidate Rows Before Matching

- ExoMol sorted: `{len(sorted_candidates)}`
- ExoMol I1: `{len(i1_candidates)}`
- HITRAN: `{len(hitran_candidates)}`

## Accepted Matches

### ExoMol Sorted vs HITRAN

- Accepted matches: `{int(sorted_stats["count"])}`
- Max `|dnu|`: `{sorted_stats["max_delta_wavenumber"]:.6f} cm^-1`
- Mean `|dnu|`: `{sorted_stats["mean_delta_wavenumber"]:.6f} cm^-1`
- Max relative intensity difference: `{sorted_stats["max_relative_intensity_difference"]:.6e}`
- Mean relative intensity difference: `{sorted_stats["mean_relative_intensity_difference"]:.6e}`

### ExoMol I1 vs HITRAN

- Accepted matches: `{int(i1_stats["count"])}`
- Max `|dnu|`: `{i1_stats["max_delta_wavenumber"]:.6f} cm^-1`
- Mean `|dnu|`: `{i1_stats["mean_delta_wavenumber"]:.6f} cm^-1`
- Max relative intensity difference: `{i1_stats["max_relative_intensity_difference"]:.6e}`
- Mean relative intensity difference: `{i1_stats["mean_relative_intensity_difference"]:.6e}`

## Plot

- Shared matched HITRAN wavenumbers between both comparisons: `{shared_hitran_points}`
- The right-side residual axis is a relabeled scale of the same residual points using `100 * residual / peak(|HITRAN intensity|)`.
- The top panel uses scatter markers, not connected lines, for the matched line comparison.
- Figure image: `{figure_path}`
- HTML plot: `{plot_path}`

## Figure

![HITRAN vs ExoMol comparison]({figure_ref})
"""


def main() -> None:
    args = parse_args()
    validate_args(args)

    lower_mode = parse_mode_text(args.lower_mode)
    upper_mode = parse_mode_text(args.upper_mode)

    sorted_grouped = read_exomol_sorted_rows(
        args.exomol_sorted_dir,
        target_lower_mode=lower_mode,
        target_upper_mode=upper_mode,
        wn_min=args.wn_min,
        wn_max=args.wn_max,
    )
    i1_grouped = read_fixed_width_rows(
        args.exomol_i1_dir,
        header_path=DEFAULT_EXOMOL_I1_HEADER,
        dataset_name="ExoMol I1",
        file_prefix=EXOMOL_I1_PREFIX,
        target_lower_mode=lower_mode,
        target_upper_mode=upper_mode,
        wn_min=args.wn_min,
        wn_max=args.wn_max,
    )
    hitran_grouped = read_fixed_width_rows(
        args.hitran_dir,
        header_path=args.hitran_header,
        dataset_name="HITRAN",
        file_prefix=HITRAN_PREFIX,
        target_lower_mode=lower_mode,
        target_upper_mode=upper_mode,
        wn_min=args.wn_min,
        wn_max=args.wn_max,
    )

    j_pair = select_j_pair(
        sorted_grouped,
        i1_grouped,
        hitran_grouped,
        lower_j=args.lower_j,
        upper_j=args.upper_j,
    )
    sorted_candidates = sorted_grouped[j_pair]
    i1_candidates = i1_grouped[j_pair]
    hitran_candidates = hitran_grouped[j_pair]

    sorted_matches = match_rows(sorted_candidates, hitran_candidates, args.match_tol)
    if not sorted_matches:
        observed = nearest_observed_delta(sorted_candidates, hitran_candidates)
        raise RuntimeError(
            f"No ExoMol sorted vs HITRAN matches survived |dnu| <= {args.match_tol:.6f} cm^-1 for J pair {j_pair}. "
            f"Nearest observed |dnu| was {observed:.6f} cm^-1."
        )

    i1_matches = match_rows(i1_candidates, hitran_candidates, args.match_tol)
    if not i1_matches:
        observed = nearest_observed_delta(i1_candidates, hitran_candidates)
        raise RuntimeError(
            f"No ExoMol I1 vs HITRAN matches survived |dnu| <= {args.match_tol:.6f} cm^-1 for J pair {j_pair}. "
            f"Nearest observed |dnu| was {observed:.6f} cm^-1."
        )

    stem = output_stem(lower_mode, upper_mode, j_pair, args.out_name)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    plot_path = args.out_dir / f"{stem}.html"
    figure_path = args.out_dir / f"{stem}.png"
    summary_path = args.out_dir / f"{stem}.md"

    plot_path.write_text(
        build_plot_html(
            lower_mode=lower_mode,
            upper_mode=upper_mode,
            j_pair=j_pair,
            sorted_matches=sorted_matches,
            i1_matches=i1_matches,
        ),
        encoding="utf-8",
    )
    write_plot_png(
        figure_path,
        lower_mode=lower_mode,
        upper_mode=upper_mode,
        j_pair=j_pair,
        sorted_matches=sorted_matches,
        i1_matches=i1_matches,
    )
    summary_path.write_text(
        build_summary_text(
            lower_mode=lower_mode,
            upper_mode=upper_mode,
            j_pair=j_pair,
            wn_min=args.wn_min,
            wn_max=args.wn_max,
            match_tol=args.match_tol,
            sorted_candidates=sorted_candidates,
            i1_candidates=i1_candidates,
            hitran_candidates=hitran_candidates,
            sorted_matches=sorted_matches,
            i1_matches=i1_matches,
            plot_path=plot_path.resolve(),
            figure_path=figure_path.resolve(),
        ),
        encoding="utf-8",
    )

    print(f"mode: {' '.join(map(str, lower_mode))} -> {' '.join(map(str, upper_mode))}")
    print(f"j_pair: {j_pair}")
    print(f"window_cm-1: {args.wn_min:.1f} to {args.wn_max:.1f}")
    print(f"sorted_candidates: {len(sorted_candidates)}")
    print(f"i1_candidates: {len(i1_candidates)}")
    print(f"hitran_candidates: {len(hitran_candidates)}")
    print(f"sorted_matches: {len(sorted_matches)}")
    print(f"i1_matches: {len(i1_matches)}")
    print(f"summary: {summary_path}")
    print(f"figure: {figure_path}")
    print(f"plot: {plot_path}")


if __name__ == "__main__":
    main()
