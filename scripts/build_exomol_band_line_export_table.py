"""
Build per-band ExoMol absorbance exports and a summary table directly from
raw ExoMol `.states/.trans/.pf/.def` files.

This script is focused on CH4 `nu3` progressions in the MM dataset. It does
not require any intermediate HITRAN-style text exports or manifest CSVs.

Outputs
-------
- One CSV/HTML absorbance curve per detected band
- One summary CSV
- One summary HTML table

Notes
-----
- The script assumes the updated ExoMol `.states` structure described by the
  ExoMol 2020 format update: the first columns are `i, E, gtot, J`, followed
  by optional uncertainty / lifetime / Lande-g columns, then the extra quantum
  labels listed in the `.def` file.
- For "pure nu3" filtering, the default is to keep only states with
  `n1 = n2 = n4 = 0` and upward transitions in `n3`.
"""

from __future__ import annotations

import argparse
import csv
import html
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import plotly.graph_objects as go

from plot_exomol_ch4_mm_absorbance import (
    DATASET_STEM,
    DATA_DIR,
    ROOT_DIR,
    interpolate_partition_function,
    iter_bz2_text_lines,
    load_partition_function,
    lte_line_intensity_cm_per_molecule,
    number_density_cm3,
    overlapping_transition_files,
    render_cross_section,
)


OUTPUT_DIR = ROOT_DIR / "exomol_ch4_mm_band_exports"
FLAT_CURVE_THRESHOLD = 1.0e-12
MODE_LABEL = "nu3"


@dataclass(frozen=True)
class BandSignature:
    n1: int
    n2: int
    n3: int
    n4: int
    l3: int
    m3: int
    gvib: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build CH4 ExoMol per-band absorbance curves and a summary table from raw MM files.",
    )
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR, help="Directory containing the ExoMol MM files.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Directory for CSV/HTML outputs.")
    parser.add_argument("--wn-min", type=float, default=3000.0, help="Minimum wavenumber in cm^-1.")
    parser.add_argument("--wn-max", type=float, default=3100.0, help="Maximum wavenumber in cm^-1.")
    parser.add_argument("--wn-step", type=float, default=0.01, help="Wavenumber step in cm^-1.")
    parser.add_argument("--temperature-k", type=float, default=600.0, help="Gas temperature in K.")
    parser.add_argument("--pressure-torr", type=float, default=3.0, help="Total pressure in Torr.")
    parser.add_argument("--mole-fraction", type=float, default=0.008, help="CH4 mole fraction.")
    parser.add_argument("--path-length-cm", type=float, default=100.0, help="Optical path length in cm.")
    parser.add_argument(
        "--intensity-threshold",
        type=float,
        default=1.0e-23,
        help="Discard transitions weaker than this LTE intensity in cm/molecule.",
    )
    parser.add_argument(
        "--line-cutoff",
        type=float,
        default=25.0,
        help="Half-width cutoff in cm^-1 used when adding each Voigt line to the grid.",
    )
    parser.add_argument(
        "--gamma0",
        type=float,
        default=None,
        help="Override the default Lorentz HWHM at 1 bar from the .def file (cm^-1/bar).",
    )
    parser.add_argument(
        "--n-exponent",
        type=float,
        default=None,
        help="Override the default temperature exponent from the .def file.",
    )
    parser.add_argument(
        "--allow-nonzero-other-modes",
        action="store_false",
        dest="require_zero_other_modes",
        help="Allow bands where n1, n2, or n4 are nonzero.",
    )
    parser.add_argument(
        "--allow-other-mode-changes",
        action="store_false",
        dest="require_same_other_modes",
        help="Allow n1, n2, or n4 to change across the transition.",
    )
    parser.add_argument(
        "--require-unit-step",
        action="store_true",
        default=False,
        help="Require unit-step progression in nu3, e.g. 0->1, 1->2.",
    )
    parser.set_defaults(
        require_zero_other_modes=True,
        require_same_other_modes=True,
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.wn_max <= args.wn_min:
        raise ValueError("--wn-max must be greater than --wn-min")
    if args.wn_step <= 0.0:
        raise ValueError("--wn-step must be positive")
    if args.temperature_k <= 0.0:
        raise ValueError("--temperature-k must be positive")
    if args.pressure_torr < 0.0:
        raise ValueError("--pressure-torr must be non-negative")
    if not (0.0 <= args.mole_fraction <= 1.0):
        raise ValueError("--mole-fraction must be between 0 and 1")
    if args.path_length_cm < 0.0:
        raise ValueError("--path-length-cm must be non-negative")
    if args.line_cutoff <= 0.0:
        raise ValueError("--line-cutoff must be positive")


def parse_exomol_def(def_path: Path) -> dict[str, object]:
    metadata: dict[str, object] = {
        "quantum_labels": {},
        "auxiliary_titles": {},
    }

    quantum_label_pattern = re.compile(r"^Quantum label (\d+)$")
    auxiliary_title_pattern = re.compile(r"^Auxiliary title (\d+)$")

    with def_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            value_text, _, comment = raw_line.partition("#")
            value_text = value_text.strip()
            comment = comment.strip()
            if not value_text:
                continue

            if comment.startswith("Isotopologue mass"):
                metadata["mass_da"] = float(value_text.split()[0])
            elif comment.startswith("No. of states in .states file"):
                metadata["nstates"] = int(value_text.split()[0])
            elif comment.startswith("Default value of Lorentzian half-width"):
                metadata["gamma0"] = float(value_text.split()[0])
            elif comment.startswith("Default value of temperature exponent"):
                metadata["n_exponent"] = float(value_text.split()[0])
            elif comment == "Uncertainty availability (1=yes, 0=no)":
                metadata["has_uncertainty"] = bool(int(value_text.split()[0]))
            elif comment == "Lifetime availability (1=yes, 0=no)":
                metadata["has_lifetime"] = bool(int(value_text.split()[0]))
            elif comment == "Lande g-factor availability (1=yes, 0=no)":
                metadata["has_lande_g"] = bool(int(value_text.split()[0]))
            else:
                quantum_match = quantum_label_pattern.match(comment)
                if quantum_match is not None:
                    metadata["quantum_labels"][int(quantum_match.group(1))] = value_text
                    continue
                auxiliary_match = auxiliary_title_pattern.match(comment)
                if auxiliary_match is not None:
                    metadata["auxiliary_titles"][int(auxiliary_match.group(1))] = value_text

    required_keys = {"mass_da", "nstates", "gamma0", "n_exponent", "has_uncertainty", "has_lifetime", "has_lande_g"}
    missing = required_keys - metadata.keys()
    if missing:
        raise RuntimeError(f"Missing required metadata in {def_path}: {sorted(missing)}")

    quantum_labels = metadata["quantum_labels"]
    if not quantum_labels:
        raise RuntimeError(f"No quantum labels found in {def_path}")

    ordered_quantum_labels = [quantum_labels[index] for index in sorted(quantum_labels)]
    ordered_auxiliary_titles = [metadata["auxiliary_titles"][index] for index in sorted(metadata["auxiliary_titles"])]
    metadata["ordered_quantum_labels"] = ordered_quantum_labels
    metadata["ordered_auxiliary_titles"] = ordered_auxiliary_titles
    return metadata


def states_column_indices(def_metadata: dict[str, object]) -> dict[str, int]:
    base_columns = ["i", "E", "gtot", "J"]
    if bool(def_metadata["has_uncertainty"]):
        base_columns.append("unc")
    if bool(def_metadata["has_lifetime"]):
        base_columns.append("lifetime")
    if bool(def_metadata["has_lande_g"]):
        base_columns.append("lande_g")

    all_columns = base_columns + list(def_metadata["ordered_quantum_labels"]) + list(def_metadata["ordered_auxiliary_titles"])
    return {name: index for index, name in enumerate(all_columns)}


def progression_category(lower_q: int, upper_q: int) -> str:
    delta_v = upper_q - lower_q
    if delta_v <= 0:
        raise ValueError(f"Expected upward {MODE_LABEL} transition, got {lower_q}->{upper_q}")
    if delta_v == 1:
        base_label = "fundamental"
    else:
        overtone_order = delta_v - 1
        overtone_names = {
            1: "first",
            2: "second",
            3: "third",
            4: "fourth",
            5: "fifth",
        }
        order_label = overtone_names.get(overtone_order, f"{overtone_order}th")
        base_label = f"{order_label}_overtone"
    if lower_q > 0:
        return f"hot_{base_label}_band"
    return f"{base_label}_band"


def build_signature(parts: list[str], column_indices: dict[str, int]) -> BandSignature:
    return BandSignature(
        n1=int(parts[column_indices["n1"]]),
        n2=int(parts[column_indices["n2"]]),
        n3=int(parts[column_indices["n3"]]),
        n4=int(parts[column_indices["n4"]]),
        l3=int(parts[column_indices["L3"]]),
        m3=int(parts[column_indices["M3"]]),
        gvib=parts[column_indices["Gvib"]],
    )


def signature_label(signature: BandSignature) -> str:
    return (
        f"{signature.n1} {signature.n2} {signature.n3} {signature.n4} "
        f"{signature.gvib} L3={signature.l3} M3={signature.m3}"
    )


def signature_stem(signature: BandSignature) -> str:
    return (
        f"{signature.n1}_{signature.n2}_{signature.n3}_{signature.n4}_"
        f"{signature.gvib}_L3_{signature.l3}_M3_{signature.m3}"
    )


def load_state_arrays(
    states_path: Path,
    nstates: int,
    column_indices: dict[str, int],
    *,
    require_zero_other_modes: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[BandSignature | None]]:
    energies = np.empty(nstates + 1, dtype=np.float64)
    g_totals = np.empty(nstates + 1, dtype=np.float64)
    signature_ids = np.zeros(nstates + 1, dtype=np.int32)
    energies.fill(np.nan)
    g_totals.fill(np.nan)

    signature_lookup: list[BandSignature | None] = [None]
    signature_to_id: dict[BandSignature, int] = {}

    required_labels = {"n1", "n2", "n3", "n4", "L3", "M3", "Gvib"}
    missing_labels = required_labels - column_indices.keys()
    if missing_labels:
        raise RuntimeError(f"Missing required state columns: {sorted(missing_labels)}")

    for line_number, raw_line in enumerate(iter_bz2_text_lines(states_path), start=1):
        parts = raw_line.split()
        if len(parts) <= max(column_indices.values()):
            continue

        state_id = int(parts[column_indices["i"]])
        energies[state_id] = float(parts[column_indices["E"]])
        g_totals[state_id] = float(parts[column_indices["gtot"]])

        signature = build_signature(parts, column_indices)
        if require_zero_other_modes and (signature.n1 != 0 or signature.n2 != 0 or signature.n4 != 0):
            continue

        signature_id = signature_to_id.get(signature)
        if signature_id is None:
            signature_id = len(signature_lookup)
            signature_lookup.append(signature)
            signature_to_id[signature] = signature_id
        signature_ids[state_id] = signature_id

        if line_number % 1_000_000 == 0:
            print(f"loaded {line_number:,} states")

    if np.isnan(energies[1:]).any() or np.isnan(g_totals[1:]).any():
        raise RuntimeError(f"State table {states_path} did not fill all expected state IDs")

    return energies, g_totals, signature_ids, signature_lookup


def is_allowed_progression(
    lower_signature: BandSignature,
    upper_signature: BandSignature,
    *,
    require_same_other_modes: bool,
    require_unit_step: bool,
) -> bool:
    if upper_signature.n3 <= lower_signature.n3:
        return False

    if require_same_other_modes:
        if (lower_signature.n1, lower_signature.n2, lower_signature.n4) != (
            upper_signature.n1,
            upper_signature.n2,
            upper_signature.n4,
        ):
            return False

    if require_unit_step and upper_signature.n3 != lower_signature.n3 + 1:
        return False

    return True


def collect_band_groups(
    transition_files: list[Path],
    energies: np.ndarray,
    g_totals: np.ndarray,
    signature_ids: np.ndarray,
    signature_lookup: list[BandSignature | None],
    *,
    wn_min: float,
    wn_max: float,
    wing: float,
    temperature_k: float,
    partition_function: float,
    intensity_threshold: float,
    require_same_other_modes: bool,
    require_unit_step: bool,
) -> tuple[dict[tuple[int, int], dict[str, object]], dict[str, int]]:
    groups: dict[tuple[int, int], dict[str, object]] = defaultdict(
        lambda: {
            "centers": [],
            "strengths": [],
            "line_count": 0,
        }
    )
    stats = {
        "files": len(transition_files),
        "parsed": 0,
        "matching_states": 0,
        "in_window": 0,
        "kept": 0,
    }

    search_min = wn_min - wing
    search_max = wn_max + wing

    for path in transition_files:
        print(f"scan {path.name}")
        for raw_line in iter_bz2_text_lines(path):
            parts = raw_line.split()
            if len(parts) != 3:
                continue

            upper_id = int(parts[0])
            lower_id = int(parts[1])
            a_value = float(parts[2])
            stats["parsed"] += 1

            upper_signature_id = int(signature_ids[upper_id])
            lower_signature_id = int(signature_ids[lower_id])
            if upper_signature_id == 0 or lower_signature_id == 0:
                continue

            upper_signature = signature_lookup[upper_signature_id]
            lower_signature = signature_lookup[lower_signature_id]
            if upper_signature is None or lower_signature is None:
                continue
            if not is_allowed_progression(
                lower_signature,
                upper_signature,
                require_same_other_modes=require_same_other_modes,
                require_unit_step=require_unit_step,
            ):
                continue
            stats["matching_states"] += 1

            wavenumber = energies[upper_id] - energies[lower_id]
            if wavenumber < search_min or wavenumber > search_max:
                continue
            stats["in_window"] += 1

            intensity = lte_line_intensity_cm_per_molecule(
                wavenumber=wavenumber,
                a_coefficient=a_value,
                g_upper=g_totals[upper_id],
                lower_energy_cm=energies[lower_id],
                temperature_k=temperature_k,
                partition_function=partition_function,
            )
            if intensity < intensity_threshold:
                continue

            key = (lower_signature_id, upper_signature_id)
            groups[key]["centers"].append(wavenumber)
            groups[key]["strengths"].append(intensity)
            groups[key]["line_count"] += 1
            stats["kept"] += 1

            if stats["parsed"] % 1_000_000 == 0:
                print(
                    f"  parsed {stats['parsed']:,} transitions, "
                    f"matching {stats['matching_states']:,}, "
                    f"window {stats['in_window']:,}, kept {stats['kept']:,}"
                )

    return groups, stats


def wavelength_um(wavenumber_cm_1: float) -> float:
    return 10000.0 / wavenumber_cm_1


def curve_peak_metrics(grid: np.ndarray, absorbance: np.ndarray) -> tuple[float, float]:
    peak_index = int(np.argmax(absorbance))
    return float(grid[peak_index]), float(absorbance[peak_index])


def weighted_line_center(line_centers: np.ndarray, line_strengths: np.ndarray) -> float:
    return float(np.sum(line_centers * line_strengths) / np.sum(line_strengths))


def center_choice(
    curve_peak_wavenumber: float,
    peak_absorbance: float,
    weighted_center: float,
) -> tuple[float, str, str]:
    if abs(peak_absorbance) > FLAT_CURVE_THRESHOLD:
        return curve_peak_wavenumber, "curve_peak", "Center from absorbance maximum in the rendered window."
    return (
        weighted_center,
        "sw_weighted_line_center",
        "Rendered curve is flat in the current window; center falls back to the intensity-weighted line position.",
    )


def save_curve(
    out_csv: Path,
    out_html: Path,
    x_values: np.ndarray,
    y_values: np.ndarray,
    trace_name: str,
    title: str,
) -> None:
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x_values, y=y_values, mode="lines", name=trace_name))
    fig.update_layout(
        title=title,
        xaxis_title="Wavenumber (cm^-1)",
        yaxis_title="Absorbance",
        template="plotly_white",
    )
    fig.update_xaxes(range=[float(x_values[0]), float(x_values[-1])])
    fig.write_html(out_html, include_plotlyjs="cdn")

    with out_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["wavenumber_cm-1", "absorbance"])
        for wavenumber, value in zip(x_values, y_values):
            writer.writerow([wavenumber, value])


def summary_html_row(row: dict[str, object]) -> str:
    return f"""
      <tr>
        <td class="num">{row["index"]}</td>
        <td>{html.escape(str(row["band_label"]))}<div class="muted">{html.escape(str(row["mode_pair"]))} | {html.escape(str(row["category"]))}</div></td>
        <td class="metric">{float(row["display_wavelength_um"]):.4f} <span class="unit">um</span></td>
        <td class="metric">{float(row["display_wavenumber_cm-1"]):.2f} <span class="unit">cm^-1</span></td>
        <td class="metric">{float(row["peak_absorbance"]):.4e}</td>
        <td>{html.escape(str(row["notes"]))}<div class="muted">Lines: {row["line_count"]} | Method: {html.escape(str(row["center_method"]))}</div></td>
        <td><a href="{html.escape(str(row["curve_html"]))}">curve html</a><br><a href="{html.escape(str(row["curve_csv"]))}">curve csv</a></td>
      </tr>"""


def write_summary_csv(rows: list[dict[str, object]], output_dir: Path) -> Path:
    summary_path = output_dir / "exomol_band_export_summary_table.csv"
    fieldnames = [
        "index",
        "band_label",
        "mode_pair",
        "category",
        "line_count",
        "display_wavelength_um",
        "display_wavenumber_cm-1",
        "center_method",
        "curve_peak_wavelength_um",
        "curve_peak_wavenumber_cm-1",
        "peak_absorbance",
        "strongest_line_wavelength_um",
        "strongest_line_wavenumber_cm-1",
        "sw_weighted_line_center_um",
        "sw_weighted_line_center_cm-1",
        "curve_csv",
        "curve_html",
        "notes",
    ]
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fieldnames})
    return summary_path


def write_summary_html(rows: list[dict[str, object]], output_dir: Path, args: argparse.Namespace) -> Path:
    html_path = output_dir / "exomol_band_export_summary_table.html"
    body_rows = "\n".join(summary_html_row(row) for row in rows)
    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>CH4 ExoMol {MODE_LABEL} Band Export Summary</title>
  <style>
    :root {{
      --bg: #f4f4f1;
      --panel: #fcfcfa;
      --ink: #2c2c2c;
      --muted: #787878;
      --grid: #d6d6cf;
      --accent: #1f4f46;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: Georgia, "Times New Roman", serif;
      color: var(--ink);
      background: radial-gradient(circle at top left, #ffffff 0, #efefe9 42%, #e8e6df 100%);
    }}
    .page {{
      max-width: 1380px;
      margin: 0 auto;
      padding: 28px 30px 40px;
    }}
    h1 {{
      margin: 0 0 10px;
      font-size: 34px;
      font-weight: 600;
      letter-spacing: -0.02em;
    }}
    .lead {{
      max-width: 980px;
      margin: 0 0 22px;
      color: var(--muted);
      font-size: 17px;
      line-height: 1.45;
    }}
    .card {{
      overflow-x: auto;
      border: 1px solid rgba(44, 44, 44, 0.08);
      border-radius: 18px;
      background: rgba(252, 252, 250, 0.86);
      box-shadow: 0 24px 60px rgba(0, 0, 0, 0.08);
      backdrop-filter: blur(8px);
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      min-width: 1120px;
    }}
    thead th {{
      padding: 18px 16px;
      text-align: left;
      font-size: 15px;
      font-weight: 600;
      border-bottom: 1px solid var(--grid);
      background: rgba(250, 250, 247, 0.96);
      position: sticky;
      top: 0;
      z-index: 1;
    }}
    tbody td {{
      padding: 20px 16px;
      vertical-align: top;
      border-bottom: 1px solid var(--grid);
      font-size: 15px;
      line-height: 1.35;
    }}
    tbody tr:hover {{
      background: rgba(31, 79, 70, 0.045);
    }}
    .num {{
      width: 58px;
      color: #9d9d96;
      font-size: 18px;
    }}
    .metric {{
      white-space: nowrap;
      font-weight: 600;
      font-size: 18px;
    }}
    .unit {{
      font-weight: 400;
      font-style: italic;
    }}
    .muted {{
      margin-top: 6px;
      color: var(--muted);
      font-size: 13px;
    }}
    a {{
      color: var(--accent);
      text-decoration: none;
    }}
    a:hover {{
      text-decoration: underline;
    }}
    .footer {{
      margin-top: 16px;
      color: var(--muted);
      font-size: 13px;
    }}
  </style>
</head>
<body>
  <div class="page">
    <h1>CH4 ExoMol {MODE_LABEL} Band Export Summary</h1>
    <p class="lead">
      This table is built directly from the raw ExoMol MM `.states/.trans/.pf/.def` files.
      Each row corresponds to one vibrational-signature band export, and the displayed
      center uses the absorbance-curve maximum when the rendered window contains signal.
    </p>
    <div class="card">
      <table>
        <thead>
          <tr>
            <th>#</th>
            <th>Band Transition</th>
            <th>Wavelength (um)</th>
            <th>Wavenumber (cm^-1)</th>
            <th>Peak Absorbance</th>
            <th>Notes</th>
            <th>Files</th>
          </tr>
        </thead>
        <tbody>
{body_rows}
        </tbody>
      </table>
    </div>
    <div class="footer">
      Source window: {args.wn_min:.0f}-{args.wn_max:.0f} cm^-1 |
      T={args.temperature_k:g} K |
      P={args.pressure_torr:g} Torr |
      x={args.mole_fraction:g} |
      L={args.path_length_cm:g} cm
    </div>
  </div>
</body>
</html>
"""
    html_path.write_text(html_text, encoding="utf-8")
    return html_path


def main() -> None:
    args = parse_args()
    validate_args(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    def_path = args.data_dir / f"{DATASET_STEM}.def"
    pf_path = args.data_dir / f"{DATASET_STEM}.pf"
    states_path = args.data_dir / f"{DATASET_STEM}.states.bz2"

    def_metadata = parse_exomol_def(def_path)
    column_indices = states_column_indices(def_metadata)

    gamma0 = args.gamma0 if args.gamma0 is not None else float(def_metadata["gamma0"])
    n_exponent = args.n_exponent if args.n_exponent is not None else float(def_metadata["n_exponent"])

    pf_temperatures, pf_values = load_partition_function(pf_path)
    partition_function = interpolate_partition_function(pf_temperatures, pf_values, args.temperature_k)

    print("loading states")
    energies, g_totals, signature_ids, signature_lookup = load_state_arrays(
        states_path=states_path,
        nstates=int(def_metadata["nstates"]),
        column_indices=column_indices,
        require_zero_other_modes=args.require_zero_other_modes,
    )

    transition_files, effective_wing = overlapping_transition_files(
        data_dir=args.data_dir,
        wn_min=args.wn_min,
        wn_max=args.wn_max,
        wing=args.line_cutoff,
    )
    if effective_wing < args.line_cutoff:
        print(
            f"warning: requested line cutoff {args.line_cutoff:g} cm^-1 reduced to "
            f"{effective_wing:g} cm^-1 because neighboring transition chunks are not available"
        )

    groups, stats = collect_band_groups(
        transition_files=transition_files,
        energies=energies,
        g_totals=g_totals,
        signature_ids=signature_ids,
        signature_lookup=signature_lookup,
        wn_min=args.wn_min,
        wn_max=args.wn_max,
        wing=effective_wing,
        temperature_k=args.temperature_k,
        partition_function=partition_function,
        intensity_threshold=args.intensity_threshold,
        require_same_other_modes=args.require_same_other_modes,
        require_unit_step=args.require_unit_step,
    )
    if not groups:
        raise RuntimeError("No ExoMol bands survived the current filters.")

    grid = np.arange(args.wn_min, args.wn_max + 0.5 * args.wn_step, args.wn_step, dtype=np.float64)
    absorber_density_cm3 = number_density_cm3(args.pressure_torr, args.temperature_k) * args.mole_fraction

    rows: list[dict[str, object]] = []
    for band_index, key in enumerate(sorted(groups), start=1):
        lower_signature = signature_lookup[key[0]]
        upper_signature = signature_lookup[key[1]]
        if lower_signature is None or upper_signature is None:
            continue

        lower_q = lower_signature.n3
        upper_q = upper_signature.n3
        category = progression_category(lower_q, upper_q)

        line_centers = np.asarray(groups[key]["centers"], dtype=np.float64)
        line_strengths = np.asarray(groups[key]["strengths"], dtype=np.float64)
        cross_section = render_cross_section(
            grid=grid,
            line_centers=line_centers,
            line_intensities=line_strengths,
            mass_da=float(def_metadata["mass_da"]),
            temperature_k=args.temperature_k,
            pressure_torr=args.pressure_torr,
            gamma0=gamma0,
            n_exponent=n_exponent,
            line_cutoff=args.line_cutoff,
        )
        absorbance = cross_section * absorber_density_cm3 * args.path_length_cm

        band_label = f"{signature_label(lower_signature)} -> {signature_label(upper_signature)}"
        pair_label = f"{MODE_LABEL} {lower_q}->{upper_q}"
        stem = (
            f"EXOMOL_CH4_{MODE_LABEL}_{lower_q}_to_{upper_q}_"
            f"{signature_stem(lower_signature)}_to_{signature_stem(upper_signature)}_"
            f"abs_{int(args.wn_min)}_{int(args.wn_max)}_"
            f"T{args.temperature_k:g}K_P{args.pressure_torr:g}Torr_"
            f"x{args.mole_fraction:g}_L{args.path_length_cm:g}cm"
        )
        out_csv = args.output_dir / f"{stem}.csv"
        out_html = args.output_dir / f"{stem}.html"
        save_curve(
            out_csv=out_csv,
            out_html=out_html,
            x_values=grid,
            y_values=absorbance,
            trace_name=f"{pair_label} (N={groups[key]['line_count']})",
            title=(
                f"CH4 ExoMol {band_label} | {pair_label} Absorbance "
                f"{int(args.wn_min)}-{int(args.wn_max)} cm^-1 | "
                f"T={args.temperature_k:g} K, P={args.pressure_torr:g} Torr, "
                f"x={args.mole_fraction:g}, L={args.path_length_cm:g} cm"
            ),
        )

        curve_peak_wavenumber, peak_absorbance = curve_peak_metrics(grid, absorbance)
        strongest_idx = int(np.argmax(line_strengths))
        strongest_line_wavenumber = float(line_centers[strongest_idx])
        sw_center = weighted_line_center(line_centers, line_strengths)
        display_wavenumber, center_method, notes = center_choice(
            curve_peak_wavenumber=curve_peak_wavenumber,
            peak_absorbance=peak_absorbance,
            weighted_center=sw_center,
        )

        rows.append(
            {
                "index": band_index,
                "band_label": band_label,
                "mode_pair": pair_label,
                "category": category,
                "line_count": int(groups[key]["line_count"]),
                "display_wavelength_um": wavelength_um(display_wavenumber),
                "display_wavenumber_cm-1": display_wavenumber,
                "center_method": center_method,
                "curve_peak_wavelength_um": wavelength_um(curve_peak_wavenumber),
                "curve_peak_wavenumber_cm-1": curve_peak_wavenumber,
                "peak_absorbance": peak_absorbance,
                "strongest_line_wavelength_um": wavelength_um(strongest_line_wavenumber),
                "strongest_line_wavenumber_cm-1": strongest_line_wavenumber,
                "sw_weighted_line_center_um": wavelength_um(sw_center),
                "sw_weighted_line_center_cm-1": sw_center,
                "curve_csv": out_csv.name,
                "curve_html": out_html.name,
                "notes": notes,
                "sort_delta_v": upper_q - lower_q,
                "sort_lower_q": lower_q,
            }
        )
        print(f"saved {out_csv}")
        print(f"saved {out_html}")

    rows.sort(
        key=lambda row: (
            int(row["sort_delta_v"]),
            int(row["sort_lower_q"]),
            -float(row["display_wavenumber_cm-1"]),
            str(row["band_label"]),
        )
    )
    for index, row in enumerate(rows, start=1):
        row["index"] = index

    summary_csv_path = write_summary_csv(rows, args.output_dir)
    summary_html_path = write_summary_html(rows, args.output_dir, args)

    print(
        f"parsed {stats['parsed']:,} transitions across {stats['files']} file(s); "
        f"matching states {stats['matching_states']:,}; window {stats['in_window']:,}; kept {stats['kept']:,}"
    )
    print(f"partition function Q({args.temperature_k:g} K) = {partition_function:.6g}")
    print(f"saved {summary_csv_path}")
    print(f"saved {summary_html_path}")


if __name__ == "__main__":
    main()
