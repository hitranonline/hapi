from __future__ import annotations

import csv
import json
import re
from collections import Counter, OrderedDict
from pathlib import Path
from typing import TextIO

import hapi
import numpy as np

from .bands import clean_quanta_label, format_band_label, safe_label_fragment
from .io import default_paths, ensure_directory, write_rows_csv
from .models import BandExportResult, BandSelection, GasCase, SpectrumResult, SpectralWindow
from .spectra import save_spectrum_csv, save_spectrum_html, to_absorbance


FIELD_FORMAT_RE = re.compile(r"%(\d+)(?:\.\d+)?[A-Za-z]")
TABLE_ID_RE = re.compile(r".*_M(\d+)_I(\d+)")
MAX_OPEN_OUTPUTS = 128


def fetch_tables(
    requests: list[tuple[str, int, int, float, float]],
    *,
    db_dir: Path | None = None,
) -> list[str]:
    paths = default_paths()
    resolved_db_dir = db_dir or paths.hitran_db_dir
    hapi.db_begin(str(resolved_db_dir))
    fetched: list[str] = []
    for table_name, molecule_id, isotopologue_id, nu_min, nu_max in requests:
        hapi.fetch(table_name, molecule_id, isotopologue_id, nu_min, nu_max)
        fetched.append(table_name)
    return fetched


def load_table(table_name: str, *, db_dir: Path | None = None) -> dict[str, object]:
    paths = default_paths()
    resolved_db_dir = db_dir or paths.hitran_db_dir
    header_path = resolved_db_dir / f"{table_name}.header"
    data_path = resolved_db_dir / f"{table_name}.data"
    if not header_path.exists():
        raise FileNotFoundError(f"Missing HITRAN header: {header_path}")
    if not data_path.exists():
        raise FileNotFoundError(f"Missing HITRAN data: {data_path}")

    with header_path.open("r", encoding="utf-8") as handle:
        header = json.load(handle)
    return {
        "table_name": table_name,
        "header_path": header_path,
        "data_path": data_path,
        "header": header,
    }


def render_absorbance(
    table_name: str,
    *,
    case: GasCase,
    window: SpectralWindow,
    db_dir: Path | None = None,
    lower_quanta: str | None = None,
    upper_quanta: str | None = None,
    intensity_threshold: float = 1.0e-23,
    diluent: dict[str, float] | None = None,
    output_dir: Path | None = None,
    output_stem: str | None = None,
) -> SpectrumResult:
    paths = default_paths()
    resolved_db_dir = db_dir or paths.hitran_db_dir
    header_path = resolved_db_dir / f"{table_name}.header"
    if not header_path.exists():
        raise FileNotFoundError(f"Missing HITRAN header: {header_path}")

    match = TABLE_ID_RE.fullmatch(table_name)
    if match is None:
        raise ValueError(f"Could not parse molecule/isotopologue IDs from {table_name}")
    molecule_id = int(match.group(1))
    isotopologue_id = int(match.group(2))

    hapi.db_begin(str(resolved_db_dir))
    source_table = table_name
    temp_table: str | None = None

    if lower_quanta is not None and upper_quanta is not None:
        temp_table = "__research_" + re.sub(r"[^A-Za-z0-9]+", "_", f"{table_name}_{lower_quanta}_{upper_quanta}").strip("_")
        hapi.select(
            table_name,
            DestinationTableName=temp_table,
            Conditions=(
                "AND",
                ("=", "global_lower_quanta", lower_quanta),
                ("=", "global_upper_quanta", upper_quanta),
            ),
            Output=False,
        )
        if hapi.length(temp_table) == 0:
            hapi.dropTable(temp_table)
            raise RuntimeError("No rows matched the requested quantum-label filter")
        source_table = temp_table

    if diluent is None:
        diluent = {"self": case.mole_fraction, "He": 1.0 - case.mole_fraction}

    try:
        wavenumber, coefficient = hapi.absorptionCoefficient_Voigt(
            Components=[(molecule_id, isotopologue_id, case.mole_fraction)],
            SourceTables=[source_table],
            WavenumberRange=[window.wn_min, window.wn_max],
            WavenumberStep=window.wn_step,
            Environment={"T": case.temperature_k, "p": case.pressure_atm},
            Diluent=diluent,
            IntensityThreshold=intensity_threshold,
            HITRAN_units=False,
        )
        _, transmittance = hapi.transmittanceSpectrum(
            wavenumber,
            coefficient,
            Environment={"l": case.path_length_cm},
        )
        absorbance = to_absorbance(np.asarray(transmittance))
        result = SpectrumResult(
            wavenumber=np.asarray(wavenumber),
            values=absorbance,
            quantity="absorbance",
            metadata={
                "table_name": table_name,
                "lower_quanta": lower_quanta,
                "upper_quanta": upper_quanta,
                "temperature_k": case.temperature_k,
                "pressure_torr": case.pressure_torr,
                "mole_fraction": case.mole_fraction,
                "path_length_cm": case.path_length_cm,
            },
        )

        if output_dir is not None:
            ensure_directory(output_dir)
            stem = output_stem or table_name
            result.csv_path = save_spectrum_csv(output_dir / f"{stem}.csv", result.wavenumber, result.values, "absorbance")
            result.html_path = save_spectrum_html(
                output_dir / f"{stem}.html",
                result.wavenumber,
                result.values,
                title=f"{table_name} Absorbance",
                trace_name=table_name,
                y_label="Absorbance",
            )
        return result
    finally:
        if temp_table is not None:
            hapi.dropTable(temp_table)


def render_progressions(*args, **kwargs):
    raise NotImplementedError(
        "render_progressions is deferred in the first framework pass. "
        "Use extract_bands plus render_absorbance as the stable path for now."
    )


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
    return format_band_label(lower, upper)


def extract_wavenumber(line: str, positions: dict[str, int], widths: dict[str, int]) -> float:
    return float(extract_field(line, "nu", positions, widths).strip())


def parse_ch4_global_quanta(value: str) -> tuple[int, int, int, int, str]:
    parts = clean_quanta_label(value).split()
    if len(parts) != 5:
        raise ValueError(f"Unexpected CH4 global quanta label: {value!r}")
    return int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]), parts[4]


def pure_nu3_mode_pair(
    band_label: str,
    *,
    require_zero_other_modes: bool,
    require_unit_step: bool,
) -> tuple[int, int] | None:
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


def extract_bands(
    *,
    table_name: str | None = None,
    selection: BandSelection | None = None,
    db_dir: Path | None = None,
    source_header: Path | None = None,
    source_data: Path | None = None,
    summary_csv: Path | None = None,
    output_dir: Path | None = None,
    expected_band_count: int | None = None,
    wn_min: float | None = None,
    wn_max: float | None = None,
    assume_ordered: bool = False,
    break_margin_cm: float = 10.0,
    require_zero_other_modes: bool | None = None,
    require_unit_step: bool | None = None,
    max_written_lines: int | None = None,
    stop_after_bands: int | None = None,
    progress_every: int = 1_000_000,
) -> BandExportResult:
    paths = default_paths()
    resolved_db_dir = db_dir or paths.hitran_db_dir
    if table_name is not None:
        source_header = source_header or resolved_db_dir / f"{table_name}.header"
        source_data = source_data or resolved_db_dir / f"{table_name}.data"
    if source_header is None or source_data is None:
        raise ValueError("source_header/source_data or table_name must be provided")
    if not source_header.exists():
        raise FileNotFoundError(f"Missing source header: {source_header}")
    if not source_data.exists():
        raise FileNotFoundError(f"Missing source data: {source_data}")
    if (wn_min is None) != (wn_max is None):
        raise ValueError("wn_min and wn_max must be provided together")
    if wn_min is not None and wn_max is not None and wn_max <= wn_min:
        raise ValueError("wn_max must be greater than wn_min")
    if progress_every <= 0:
        raise ValueError("progress_every must be positive")

    if selection is not None:
        if selection.species != "CH4" or selection.mode_label != "nu3":
            raise NotImplementedError("The first framework pass supports CH4 nu3 extraction only")
        if require_zero_other_modes is None:
            require_zero_other_modes = selection.require_same_other_modes
        if require_unit_step is None:
            require_unit_step = selection.require_unit_step
    if require_zero_other_modes is None:
        require_zero_other_modes = True
    if require_unit_step is None:
        require_unit_step = False

    selection_mode = "summary" if summary_csv is not None else "pure-nu3"
    output_dir = output_dir or (paths.artifacts_dir / "bands" / source_data.stem)
    ensure_directory(output_dir)

    positions, widths = load_header_metadata(source_header)
    source_stem = source_data.stem
    target_bands: list[str]
    expected_counts: dict[str, int]
    if selection_mode == "summary":
        target_bands, expected_counts = load_target_bands(summary_csv, expected_band_count)
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

    if selection_mode == "summary":
        for band_label in target_bands:
            band_paths[band_label] = output_path_for_band(output_dir, source_stem, band_label)
            band_paths[band_label].write_text("", encoding="utf-8")

    try:
        with source_data.open("r", encoding="utf-8", newline="") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\r\n")
                if not line:
                    continue

                stats["scanned"] += 1
                if stats["scanned"] % progress_every == 0:
                    print(
                        f"scanned {stats['scanned']:,} rows, in-window {stats['in_window']:,}, "
                        f"matched {stats['band_matches']:,}, written {stats['written']:,}"
                    )

                wavenumber = extract_wavenumber(line, positions, widths)
                if wn_min is not None and wavenumber < wn_min:
                    continue
                if wn_max is not None and wavenumber > wn_max:
                    if assume_ordered and wavenumber > wn_max + break_margin_cm:
                        break
                    continue
                stats["in_window"] += 1

                band_label = extract_band_label(line, positions, widths)
                if selection_mode == "summary":
                    if band_label not in band_paths:
                        continue
                else:
                    try:
                        if pure_nu3_mode_pair(
                            band_label,
                            require_zero_other_modes=require_zero_other_modes,
                            require_unit_step=require_unit_step,
                        ) is None:
                            continue
                    except ValueError:
                        stats["bad_quanta"] += 1
                        continue
                    if band_label not in band_paths:
                        band_paths[band_label] = output_path_for_band(output_dir, source_stem, band_label)
                        expected_counts[band_label] = -1
                        target_bands.append(band_label)

                stats["band_matches"] += 1
                handle_out = open_band_handle(band_label, band_paths=band_paths, output_handles=output_handles)
                handle_out.write(line + "\n")
                if matched_counts[band_label] == 0:
                    exported_band_count += 1
                matched_counts[band_label] += 1
                stats["written"] += 1

                if max_written_lines is not None and stats["written"] >= max_written_lines:
                    break
                if stop_after_bands is not None and exported_band_count >= stop_after_bands:
                    break
    finally:
        for handle in output_handles.values():
            handle.close()

    if selection_mode == "summary":
        for band_label in target_bands:
            if band_label not in matched_counts:
                matched_counts[band_label] = 0

    if not band_paths:
        raise RuntimeError("No output files were created for the current configuration")
    if exported_band_count == 0:
        raise RuntimeError("No matching CH4 nu3 rows were found for the current filters")

    summary_rows = []
    for band_label in sorted(
        matched_counts.keys(),
        key=lambda label: (
            pure_nu3_mode_pair(label, require_zero_other_modes=False, require_unit_step=False) or (999, 999),
            label,
        ),
    ):
        summary_rows.append(
            {
                "band_label": band_label,
                "mode_pair": mode_pair_label_from_band_label(band_label),
                "line_count": matched_counts[band_label],
                "reference_line_count": "" if expected_counts[band_label] < 0 else expected_counts[band_label],
                "txt_path": str(band_paths[band_label]),
            }
        )

    manifest_path = write_rows_csv(output_summary_path(output_dir, source_stem), summary_rows)
    return BandExportResult(
        rows=summary_rows,
        output_dir=output_dir,
        manifest_path=manifest_path,
        count=exported_band_count,
        metadata=stats,
    )
