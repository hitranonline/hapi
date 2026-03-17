"""
Build a summary table from the generated band-line export curves.

Outputs:
- CSV summary with derived peak/center metrics
- HTML table styled for quick browsing
"""

from __future__ import annotations

import csv
import html
from pathlib import Path

import hapi
import render_band_line_text_curves as rendered_curves


ROOT_DIR = Path(__file__).resolve().parent
TXT_DIR = rendered_curves.TXT_DIR
CURVE_DIR = rendered_curves.OUTPUT_DIR
OUT_CSV = CURVE_DIR / "band_line_export_summary_table.csv"
OUT_HTML = CURVE_DIR / "band_line_export_summary_table.html"

FLAT_CURVE_THRESHOLD = 1.0e-12
MODE_LABEL = "nu3"
MODE_INDEX = 2


def parse_band_labels(txt_path: Path) -> tuple[str, str]:
    stem = txt_path.stem
    prefix = f"{rendered_curves.SOURCE_TABLE}_"
    if not stem.startswith(prefix):
        raise ValueError(f"Unexpected file name: {txt_path.name}")
    payload = stem[len(prefix) :]
    lower_token, upper_token = payload.split("_to_", maxsplit=1)
    return lower_token.replace("_", " "), upper_token.replace("_", " ")


def parse_ch4_global_quanta(value: str) -> tuple[int, ...]:
    parts = value.split()
    if len(parts) != 5:
        raise ValueError(f"Unexpected CH4 global quanta label: {value!r}")
    return tuple(int(parts[i]) for i in range(4))


def mode_pair(lower_raw: str, upper_raw: str) -> tuple[int, int]:
    lower_state = parse_ch4_global_quanta(lower_raw)
    upper_state = parse_ch4_global_quanta(upper_raw)
    lower_mode_q = lower_state[MODE_INDEX]
    upper_mode_q = upper_state[MODE_INDEX]
    if upper_mode_q <= lower_mode_q:
        raise ValueError(f"Expected upward nu3 transition, got {lower_raw!r} -> {upper_raw!r}")
    return lower_mode_q, upper_mode_q


def progression_category(lower_q: int, upper_q: int) -> str:
    delta_v = upper_q - lower_q
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


def line_metrics(txt_path: Path) -> tuple[int, float, float]:
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
        "Rendered curve is flat in the current window; center falls back to SW-weighted line position.",
    )


def wavelength_um(wavenumber_cm_1: float) -> float:
    return 10000.0 / wavenumber_cm_1


def build_rows() -> list[dict[str, object]]:
    rendered_curves.bootstrap_source_schema()

    rows: list[dict[str, object]] = []
    for txt_path in sorted(TXT_DIR.glob("*.txt")):
        curve_csv_path = CURVE_DIR / f"{txt_path.stem}.csv"
        curve_html_path = CURVE_DIR / f"{txt_path.stem}.html"
        if not curve_csv_path.exists():
            raise FileNotFoundError(f"Missing exported curve CSV: {curve_csv_path}")
        if not curve_html_path.exists():
            raise FileNotFoundError(f"Missing exported curve HTML: {curve_html_path}")

        lower_raw, upper_raw = parse_band_labels(txt_path)
        lower_q, upper_q = mode_pair(lower_raw, upper_raw)
        category = progression_category(lower_q, upper_q)
        pair_label = f"{MODE_LABEL} {lower_q}->{upper_q}"

        curve_peak_wavenumber, peak_absorbance = curve_metrics(curve_csv_path)
        line_count, strongest_line_wavenumber, weighted_line_center = line_metrics(txt_path)
        display_wavenumber, center_method, notes = center_choice(
            curve_peak_wavenumber=curve_peak_wavenumber,
            peak_absorbance=peak_absorbance,
            weighted_line_center=weighted_line_center,
        )

        rows.append(
            {
                "band_label": f"{lower_raw} -> {upper_raw}",
                "mode_pair": pair_label,
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
                "curve_csv": curve_csv_path.name,
                "curve_html": curve_html_path.name,
                "band_text": txt_path.name,
                "notes": notes,
                "sort_delta_v": upper_q - lower_q,
                "sort_lower_q": lower_q,
            }
        )

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


def write_csv(rows: list[dict[str, object]]) -> None:
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
    with OUT_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row[name] for name in fieldnames})


def html_row(row: dict[str, object]) -> str:
    curve_csv = html.escape(str(row["curve_csv"]))
    curve_html = html.escape(str(row["curve_html"]))
    band_text = html.escape(str(row["band_text"]))
    note_text = html.escape(str(row["notes"]))
    return f"""
      <tr>
        <td class="num">{row["index"]}</td>
        <td>{html.escape(str(row["band_label"]))}<div class="muted">{html.escape(str(row["mode_pair"]))} | {html.escape(str(row["category"]))}</div></td>
        <td class="metric">{float(row["display_wavelength_um"]):.4f} <span class="unit">um</span></td>
        <td class="metric">{float(row["display_wavenumber_cm-1"]):.2f} <span class="unit">cm^-1</span></td>
        <td class="metric">{float(row["peak_absorbance"]):.4e}</td>
        <td>{note_text}<div class="muted">Lines: {row["line_count"]} | Method: {html.escape(str(row["center_method"]))}</div></td>
        <td><a href="{curve_html}">curve html</a><br><a href="{curve_csv}">curve csv</a><br><a href="../band_line_texts/{band_text}">band txt</a></td>
      </tr>"""


def write_html(rows: list[dict[str, object]]) -> None:
    body_rows = "\n".join(html_row(row) for row in rows)
    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>CH4 nu3 Band Line Export Summary</title>
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
    <h1>CH4 nu3 Band Line Export Summary</h1>
    <p class="lead">
      This table is derived from the generated files in <code>ch4_nu3_progressions/band_line_text_exports</code>.
      The displayed center uses the absorbance-curve maximum when the rendered window contains signal; otherwise it
      falls back to the SW-weighted line-list center from the corresponding band text file.
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
      Source window: {rendered_curves.WN_MIN:.0f}-{rendered_curves.WN_MAX:.0f} cm^-1 |
      T={rendered_curves.T_K:g} K |
      P={rendered_curves.P_TORR:g} Torr |
      x={rendered_curves.MOLE_FRACTION:g} |
      L={rendered_curves.PATH_LENGTH_CM:g} cm
    </div>
  </div>
</body>
</html>
"""
    OUT_HTML.write_text(html_text, encoding="utf-8")


def main() -> None:
    CURVE_DIR.mkdir(parents=True, exist_ok=True)
    rows = build_rows()
    write_csv(rows)
    write_html(rows)
    print(f"Saved: {OUT_CSV}")
    print(f"Saved: {OUT_HTML}")


if __name__ == "__main__":
    main()
