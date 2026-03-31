"""
Build a band-summary table from previously generated per-band curve exports.

Inputs
------
This script supports two input modes:

1. `--input-kind=hitran`
   Reads the existing CH4 nu3 HITRAN workflow outputs:
   - per-band line-text files in `ch4_nu3_progressions/band_line_texts`
   - per-band curve CSV/HTML files produced by `render_band_line_text_curves.py`

2. `--input-kind=exomol`
   Reads a manifest CSV with one exported band per row. This is intended for
   ExoMol-style workflows and other manifest-driven band exports.

Manifest columns for `--input-kind=exomol`
------------------------------------------
Required:
- `band_label`
- `mode_pair` or `mode_pair_label`   example: `nu3 0->1`
- `line_count`
- `csv_path` or `curve_csv`

Optional:
- `html_path` or `curve_html`
- `category`
- `strongest_line_wavenumber_cm-1`
- `sw_weighted_line_center_cm-1`
- `notes`

Outputs
-------
The script writes a summary table for the selected input source:
- summary CSV in the output directory
- summary HTML in the output directory

Each output row represents one band and includes metadata such as band label,
mode pair, line count, category, linked curve paths, and peak/center metrics.
"""

from __future__ import annotations

import argparse
import csv
import html
import os
import re
from pathlib import Path

from _bootstrap import ensure_repo_root

ensure_repo_root()

import hapi
import render_band_line_text_curves as rendered_curves


ROOT_DIR = Path(__file__).resolve().parents[1]

HITRAN_TXT_DIR = rendered_curves.TXT_DIR
HITRAN_CURVE_DIR = rendered_curves.OUTPUT_DIR

DEFAULT_INPUT_KIND = "hitran"
DEFAULT_EXOMOL_MANIFEST = ROOT_DIR / "exomol_ch4_mm_plots" / "exomol_band_summary.csv"

FLAT_CURVE_THRESHOLD = 1.0e-12
HITRAN_MODE_LABEL = "nu3"
HITRAN_MODE_INDEX = 2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a band summary table from HITRAN text exports or an ExoMol manifest CSV.",
    )
    parser.add_argument(
        "--input-kind",
        choices=("hitran", "exomol"),
        default=DEFAULT_INPUT_KIND,
        help="Input source for the band summary table.",
    )
    parser.add_argument(
        "--manifest-csv",
        type=Path,
        default=None,
        help="Manifest CSV for ExoMol-style inputs. Required when --input-kind=exomol.",
    )
    parser.add_argument(
        "--curve-dir",
        type=Path,
        default=None,
        help="Override the directory used to resolve relative curve paths from the manifest.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for the generated summary CSV/HTML. Defaults to the active curve directory.",
    )
    return parser.parse_args()


def parse_ch4_global_quanta(value: str) -> tuple[int, ...]:
    parts = value.split()
    if len(parts) != 5:
        raise ValueError(f"Unexpected CH4 global quanta label: {value!r}")
    return tuple(int(parts[i]) for i in range(4))


def hitran_mode_pair(lower_raw: str, upper_raw: str) -> tuple[int, int]:
    lower_state = parse_ch4_global_quanta(lower_raw)
    upper_state = parse_ch4_global_quanta(upper_raw)
    lower_mode_q = lower_state[HITRAN_MODE_INDEX]
    upper_mode_q = upper_state[HITRAN_MODE_INDEX]
    if upper_mode_q <= lower_mode_q:
        raise ValueError(f"Expected upward nu3 transition, got {lower_raw!r} -> {upper_raw!r}")
    return lower_mode_q, upper_mode_q


def parse_mode_pair_label(mode_pair_label: str) -> tuple[str, int, int]:
    match = re.fullmatch(r"\s*([A-Za-z0-9_]+)\s+(\d+)\s*->\s*(\d+)\s*", mode_pair_label)
    if match is None:
        raise ValueError(f"Unexpected mode-pair label: {mode_pair_label!r}")
    return match.group(1), int(match.group(2)), int(match.group(3))


def progression_category(mode_label: str, lower_q: int, upper_q: int) -> str:
    delta_v = upper_q - lower_q
    if delta_v <= 0:
        raise ValueError(f"Expected upward {mode_label} transition, got {lower_q}->{upper_q}")
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


def curve_metrics(curve_csv_path: Path) -> tuple[float, float]:
    peak_wavenumber = 0.0
    peak_absorbance = float("-inf")

    with curve_csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            wavenumber = float(row["wavenumber_cm-1"])
            absorbance = float(row["absorbance"])
            if absorbance > peak_absorbance:
                peak_wavenumber = wavenumber
                peak_absorbance = absorbance

    if peak_absorbance == float("-inf"):
        raise RuntimeError(f"No data rows found in {curve_csv_path}")
    return peak_wavenumber, peak_absorbance


def hitran_line_metrics(txt_path: Path) -> tuple[int, float, float]:
    temp_table = rendered_curves.temp_table_name(txt_path)
    try:
        line_count = rendered_curves.build_temp_table_from_text(txt_path, temp_table)
        nus = [float(value) for value in hapi.getColumn(temp_table, "nu")]
        strengths = [float(value) for value in hapi.getColumn(temp_table, "sw")]
    finally:
        hapi.dropTable(temp_table)

    strongest_index = max(range(len(strengths)), key=strengths.__getitem__)
    strongest_line_wavenumber = nus[strongest_index]
    weighted_line_center = sum(nu * sw for nu, sw in zip(nus, strengths)) / sum(strengths)
    return line_count, strongest_line_wavenumber, weighted_line_center


def center_choice(
    curve_peak_wavenumber: float,
    peak_absorbance: float,
    weighted_line_center: float,
) -> tuple[float, str, str]:
    if abs(peak_absorbance) > FLAT_CURVE_THRESHOLD:
        return curve_peak_wavenumber, "curve_peak", "Center from absorbance maximum in the rendered window."
    return (
        weighted_line_center,
        "sw_weighted_line_center",
        "Rendered curve is flat in the current window; center falls back to the SW-weighted line position.",
    )


def wavelength_um(wavenumber_cm_1: float) -> float:
    return 10000.0 / wavenumber_cm_1


def parse_band_labels(txt_path: Path) -> tuple[str, str]:
    stem = txt_path.stem
    prefix = f"{rendered_curves.SOURCE_TABLE}_"
    if not stem.startswith(prefix):
        raise ValueError(f"Unexpected file name: {txt_path.name}")
    payload = stem[len(prefix) :]
    lower_token, upper_token = payload.split("_to_", maxsplit=1)
    return lower_token.replace("_", " "), upper_token.replace("_", " ")


def manifest_path_value(row: dict[str, str], keys: tuple[str, ...], base_dir: Path) -> Path | None:
    for key in keys:
        value = row.get(key, "").strip()
        if not value:
            continue
        path = Path(value)
        if not path.is_absolute():
            path = base_dir / path
        return path
    return None


def manifest_float_value(row: dict[str, str], keys: tuple[str, ...]) -> float | None:
    for key in keys:
        value = row.get(key, "").strip()
        if not value:
            continue
        return float(value)
    return None


def relative_path_text(path: Path | None, output_dir: Path) -> str:
    if path is None:
        return ""
    return Path(os.path.relpath(path, start=output_dir)).as_posix()


def make_row(
    *,
    band_label: str,
    mode_label: str,
    lower_q: int,
    upper_q: int,
    category: str,
    line_count: int,
    curve_csv_path: Path,
    curve_html_path: Path | None,
    strongest_line_wavenumber: float,
    weighted_line_center: float,
    band_text_path: Path | None = None,
    notes: str | None = None,
) -> dict[str, object]:
    curve_peak_wavenumber, peak_absorbance = curve_metrics(curve_csv_path)
    display_wavenumber, center_method, fallback_notes = center_choice(
        curve_peak_wavenumber=curve_peak_wavenumber,
        peak_absorbance=peak_absorbance,
        weighted_line_center=weighted_line_center,
    )

    return {
        "band_label": band_label,
        "mode_pair": f"{mode_label} {lower_q}->{upper_q}",
        "category": category,
        "line_count": line_count,
        "display_wavenumber_cm-1": display_wavenumber,
        "display_wavelength_um": wavelength_um(display_wavenumber),
        "center_method": center_method,
        "curve_peak_wavenumber_cm-1": curve_peak_wavenumber,
        "curve_peak_wavelength_um": wavelength_um(curve_peak_wavenumber),
        "peak_absorbance": peak_absorbance,
        "strongest_line_wavenumber_cm-1": strongest_line_wavenumber,
        "strongest_line_wavelength_um": wavelength_um(strongest_line_wavenumber),
        "sw_weighted_line_center_cm-1": weighted_line_center,
        "sw_weighted_line_center_um": wavelength_um(weighted_line_center),
        "curve_csv_path": curve_csv_path,
        "curve_html_path": curve_html_path,
        "band_text_path": band_text_path,
        "notes": notes.strip() if notes else fallback_notes,
        "sort_delta_v": upper_q - lower_q,
        "sort_lower_q": lower_q,
    }


def build_hitran_rows() -> list[dict[str, object]]:
    rendered_curves.bootstrap_source_schema()

    rows: list[dict[str, object]] = []
    for txt_path in sorted(HITRAN_TXT_DIR.glob("*.txt")):
        curve_csv_path = HITRAN_CURVE_DIR / f"{txt_path.stem}.csv"
        curve_html_path = HITRAN_CURVE_DIR / f"{txt_path.stem}.html"
        if not curve_csv_path.exists():
            raise FileNotFoundError(f"Missing exported curve CSV: {curve_csv_path}")
        if not curve_html_path.exists():
            raise FileNotFoundError(f"Missing exported curve HTML: {curve_html_path}")

        lower_raw, upper_raw = parse_band_labels(txt_path)
        lower_q, upper_q = hitran_mode_pair(lower_raw, upper_raw)
        category = progression_category(HITRAN_MODE_LABEL, lower_q, upper_q)
        line_count, strongest_line_wavenumber, weighted_line_center = hitran_line_metrics(txt_path)

        rows.append(
            make_row(
                band_label=f"{lower_raw} -> {upper_raw}",
                mode_label=HITRAN_MODE_LABEL,
                lower_q=lower_q,
                upper_q=upper_q,
                category=category,
                line_count=line_count,
                curve_csv_path=curve_csv_path,
                curve_html_path=curve_html_path,
                strongest_line_wavenumber=strongest_line_wavenumber,
                weighted_line_center=weighted_line_center,
                band_text_path=txt_path,
            )
        )

    return rows


def build_manifest_rows(manifest_path: Path, curve_dir: Path) -> list[dict[str, object]]:
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing manifest CSV: {manifest_path}")

    rows: list[dict[str, object]] = []
    with manifest_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for manifest_row in reader:
            band_label = manifest_row.get("band_label", "").strip()
            if not band_label:
                raise ValueError(f"Manifest row is missing band_label: {manifest_row}")

            mode_pair_label = (
                manifest_row.get("mode_pair", "").strip()
                or manifest_row.get("mode_pair_label", "").strip()
            )
            if not mode_pair_label:
                raise ValueError(f"Manifest row is missing mode_pair/mode_pair_label: {manifest_row}")
            mode_label, lower_q, upper_q = parse_mode_pair_label(mode_pair_label)

            category = manifest_row.get("category", "").strip()
            if not category:
                category = progression_category(mode_label, lower_q, upper_q)

            line_count_text = manifest_row.get("line_count", "").strip()
            if not line_count_text:
                raise ValueError(f"Manifest row is missing line_count: {manifest_row}")
            line_count = int(line_count_text)

            curve_csv_path = manifest_path_value(manifest_row, ("csv_path", "curve_csv"), curve_dir)
            if curve_csv_path is None:
                raise ValueError(f"Manifest row is missing csv_path/curve_csv: {manifest_row}")
            if not curve_csv_path.exists():
                raise FileNotFoundError(f"Missing curve CSV referenced by manifest: {curve_csv_path}")

            curve_html_path = manifest_path_value(manifest_row, ("html_path", "curve_html"), curve_dir)
            if curve_html_path is not None and not curve_html_path.exists():
                raise FileNotFoundError(f"Missing curve HTML referenced by manifest: {curve_html_path}")

            band_text_path = manifest_path_value(manifest_row, ("band_text_path", "band_text"), curve_dir)
            if band_text_path is not None and not band_text_path.exists():
                band_text_path = None

            strongest_line_wavenumber = manifest_float_value(
                manifest_row,
                ("strongest_line_wavenumber_cm-1", "strongest_line_wavenumber"),
            )
            weighted_line_center = manifest_float_value(
                manifest_row,
                ("sw_weighted_line_center_cm-1", "sw_weighted_line_center"),
            )

            curve_peak_wavenumber, _ = curve_metrics(curve_csv_path)
            if strongest_line_wavenumber is None:
                strongest_line_wavenumber = curve_peak_wavenumber
            if weighted_line_center is None:
                weighted_line_center = strongest_line_wavenumber

            rows.append(
                make_row(
                    band_label=band_label,
                    mode_label=mode_label,
                    lower_q=lower_q,
                    upper_q=upper_q,
                    category=category,
                    line_count=line_count,
                    curve_csv_path=curve_csv_path,
                    curve_html_path=curve_html_path,
                    strongest_line_wavenumber=strongest_line_wavenumber,
                    weighted_line_center=weighted_line_center,
                    band_text_path=band_text_path,
                    notes=manifest_row.get("notes", ""),
                )
            )

    return rows


def sort_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
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
    return rows


def write_csv(rows: list[dict[str, object]], output_dir: Path) -> Path:
    output_path = output_dir / "band_line_export_summary_table.csv"
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
        "band_text",
        "notes",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "index": row["index"],
                    "band_label": row["band_label"],
                    "mode_pair": row["mode_pair"],
                    "category": row["category"],
                    "line_count": row["line_count"],
                    "display_wavelength_um": row["display_wavelength_um"],
                    "display_wavenumber_cm-1": row["display_wavenumber_cm-1"],
                    "center_method": row["center_method"],
                    "curve_peak_wavelength_um": row["curve_peak_wavelength_um"],
                    "curve_peak_wavenumber_cm-1": row["curve_peak_wavenumber_cm-1"],
                    "peak_absorbance": row["peak_absorbance"],
                    "strongest_line_wavelength_um": row["strongest_line_wavelength_um"],
                    "strongest_line_wavenumber_cm-1": row["strongest_line_wavenumber_cm-1"],
                    "sw_weighted_line_center_um": row["sw_weighted_line_center_um"],
                    "sw_weighted_line_center_cm-1": row["sw_weighted_line_center_cm-1"],
                    "curve_csv": relative_path_text(row["curve_csv_path"], output_dir),
                    "curve_html": relative_path_text(row["curve_html_path"], output_dir),
                    "band_text": relative_path_text(row["band_text_path"], output_dir),
                    "notes": row["notes"],
                }
            )
    return output_path


def html_row(row: dict[str, object], output_dir: Path) -> str:
    curve_csv = html.escape(relative_path_text(row["curve_csv_path"], output_dir))
    curve_html = html.escape(relative_path_text(row["curve_html_path"], output_dir))
    band_text = html.escape(relative_path_text(row["band_text_path"], output_dir))
    note_text = html.escape(str(row["notes"]))

    links: list[str] = []
    if curve_html:
        links.append(f'<a href="{curve_html}">curve html</a>')
    if curve_csv:
        links.append(f'<a href="{curve_csv}">curve csv</a>')
    if band_text:
        links.append(f'<a href="{band_text}">band text</a>')
    files_html = "<br>".join(links) if links else '<span class="muted">No file links</span>'

    return f"""
      <tr>
        <td class="num">{row["index"]}</td>
        <td>{html.escape(str(row["band_label"]))}<div class="muted">{html.escape(str(row["mode_pair"]))} | {html.escape(str(row["category"]))}</div></td>
        <td class="metric">{float(row["display_wavelength_um"]):.4f} <span class="unit">um</span></td>
        <td class="metric">{float(row["display_wavenumber_cm-1"]):.2f} <span class="unit">cm^-1</span></td>
        <td class="metric">{float(row["peak_absorbance"]):.4e}</td>
        <td>{note_text}<div class="muted">Lines: {row["line_count"]} | Method: {html.escape(str(row["center_method"]))}</div></td>
        <td>{files_html}</td>
      </tr>"""


def write_html(
    rows: list[dict[str, object]],
    output_dir: Path,
    *,
    title: str,
    lead: str,
    footer: str,
) -> Path:
    output_path = output_dir / "band_line_export_summary_table.html"
    body_rows = "\n".join(html_row(row, output_dir) for row in rows)
    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(title)}</title>
  <style>
    :root {{
      --bg: #f4f4f1;
      --panel: #fcfcfa;
      --ink: #2c2c2c;
      --muted: #787878;
      --grid: #d6d6cf;
      --accent: #1f4f46;
    }}
    * {{
      box-sizing: border-box;
    }}
    body {{
      margin: 0;
      font-family: Georgia, "Times New Roman", serif;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, #ffffff 0, #efefe9 42%, #e8e6df 100%);
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
    <h1>{html.escape(title)}</h1>
    <p class="lead">{html.escape(lead)}</p>
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
    <div class="footer">{html.escape(footer)}</div>
  </div>
</body>
</html>
"""
    output_path.write_text(html_text, encoding="utf-8")
    return output_path


def build_rows(args: argparse.Namespace) -> tuple[list[dict[str, object]], Path, str, str, str]:
    if args.input_kind == "hitran":
        rows = build_hitran_rows()
        curve_dir = HITRAN_CURVE_DIR
        lead = (
            "This table is derived from the generated files in "
            "`ch4_nu3_progressions/band_line_text_exports`. The displayed center uses the "
            "absorbance-curve maximum when the rendered window contains signal; otherwise it "
            "falls back to the SW-weighted line-list center from the corresponding band text file."
        )
        footer = (
            f"Source window: {rendered_curves.WN_MIN:.0f}-{rendered_curves.WN_MAX:.0f} cm^-1 | "
            f"T={rendered_curves.T_K:g} K | "
            f"P={rendered_curves.P_TORR:g} Torr | "
            f"x={rendered_curves.MOLE_FRACTION:g} | "
            f"L={rendered_curves.PATH_LENGTH_CM:g} cm"
        )
        return sort_rows(rows), curve_dir, "CH4 nu3 Band Line Export Summary", lead, footer

    manifest_path = args.manifest_csv or DEFAULT_EXOMOL_MANIFEST
    curve_dir = args.curve_dir or manifest_path.parent
    rows = build_manifest_rows(manifest_path, curve_dir)
    lead = (
        f"This table is derived from the manifest `{manifest_path.name}`. "
        "Curve centers use the absorbance maximum when the rendered window contains signal; "
        "otherwise they fall back to the SW-weighted or supplied line-list center."
    )
    footer = f"Manifest source: {manifest_path}"
    return sort_rows(rows), curve_dir, "Band Line Export Summary", lead, footer


def main() -> None:
    args = parse_args()
    rows, default_output_dir, title, lead, footer = build_rows(args)
    output_dir = args.output_dir or default_output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    out_csv = write_csv(rows, output_dir)
    out_html = write_html(rows, output_dir, title=title, lead=lead, footer=footer)
    print(f"Saved: {out_csv}")
    print(f"Saved: {out_html}")


if __name__ == "__main__":
    main()
