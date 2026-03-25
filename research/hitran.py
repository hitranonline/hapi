from __future__ import annotations

import csv
import json
import math
import re
import shutil
from collections import Counter, OrderedDict, defaultdict
from pathlib import Path
from typing import TextIO

import hapi
import numpy as np
import plotly.graph_objects as go
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
from plotly.subplots import make_subplots

from .bands import clean_quanta_label, format_band_label, safe_label_fragment
from .io import default_paths, ensure_directory, write_markdown, write_rows_csv
from .models import BandExportResult, BandSelection, GasCase, SpectrumResult, SpectralWindow
from .spectra import compute_panel_y_limits, save_spectrum_csv, save_spectrum_html, to_absorbance


FIELD_FORMAT_RE = re.compile(r"%(\d+)(?:\.\d+)?[A-Za-z]")
TABLE_ID_RE = re.compile(r".*_M(\d+)_I(\d+)")
MAX_OPEN_OUTPUTS = 128
HITRAN_BAND_TEXT_FILENAME_RE = re.compile(
    r"^CH4_M6_I1_"
    r"(?P<lower_n1>\d+)_(?P<lower_n2>\d+)_(?P<lower_n3>\d+)_(?P<lower_n4>\d+)_(?P<lower_sym>[^_]+)"
    r"_to_"
    r"(?P<upper_n1>\d+)_(?P<upper_n2>\d+)_(?P<upper_n3>\d+)_(?P<upper_n4>\d+)_(?P<upper_sym>[^.]+)"
    r"\.txt$"
)
J_VALUE_RE = re.compile(r"(\d+)")
REQUIRED_HITRAN_LABELS = ("J 2->3", "J 3->4")
ABSORBANCE_DELTA_J_VALUES = (-1, 0, 1)
DEFAULT_FORCED_HITRAN_J_PAIRS = ((2, 3), (3, 4))


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

    runtime_db_dir = ensure_directory((paths.artifacts_dir / "_runtime_hitran_db" / table_name).resolve())
    prepared_runtime_db = _prepare_runtime_db(
        source_db_dir=resolved_db_dir,
        runtime_db_dir=runtime_db_dir,
        source_table=table_name,
    )
    hapi.db_begin(str(prepared_runtime_db))
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


def band_line_text_dir(root_dir: Path | None = None) -> Path:
    paths = default_paths(root_dir=root_dir)
    return paths.root_dir / "ch4_nu3_progressions" / "band_line_texts"


def _parse_band_text_filename(path: Path) -> dict[str, object]:
    match = HITRAN_BAND_TEXT_FILENAME_RE.match(path.name)
    if match is None:
        raise ValueError(f"Unrecognized HITRAN band text filename: {path.name}")
    groups = match.groupdict()
    lower_mode = tuple(int(groups[f"lower_n{index}"]) for index in range(1, 5))
    upper_mode = tuple(int(groups[f"upper_n{index}"]) for index in range(1, 5))
    return {
        "path": path,
        "lower_mode": lower_mode,
        "upper_mode": upper_mode,
        "lower_sym": str(groups["lower_sym"]),
        "upper_sym": str(groups["upper_sym"]),
    }


def _mode_label(mode: tuple[int, int, int, int]) -> str:
    return " ".join(str(value) for value in mode)


def _progression_label(lower_mode: tuple[int, int, int, int], upper_mode: tuple[int, int, int, int]) -> str:
    return f"nu3 {lower_mode[2]}->{upper_mode[2]} | {_mode_label(lower_mode)} -> {_mode_label(upper_mode)}"


def _progression_slug(lower_mode: tuple[int, int, int, int], upper_mode: tuple[int, int, int, int]) -> str:
    return (
        f"nu3_{lower_mode[2]}_to_{upper_mode[2]}__"
        f"{lower_mode[0]}_{lower_mode[1]}_{lower_mode[2]}_{lower_mode[3]}_to_"
        f"{upper_mode[0]}_{upper_mode[1]}_{upper_mode[2]}_{upper_mode[3]}"
    )


def _extract_j_value(local_label: str) -> int:
    match = J_VALUE_RE.search(local_label.strip())
    if match is None:
        raise ValueError(f"Could not extract J from local label: {local_label!r}")
    return int(match.group(1))


def _format_jpair_label(lower_j: int, upper_j: int) -> str:
    return f"J {lower_j}->{upper_j}"


def _delta_j_value(lower_j: int, upper_j: int) -> int:
    return int(upper_j) - int(lower_j)


def _delta_j_branch_label(delta_j: int) -> str:
    return f"dJ_{delta_j:+d}"


def _delta_j_panel_title(delta_j: int, trace_count: int) -> str:
    return f"delta J = {delta_j:+d} ({trace_count} J pairs)"


def _merge_forced_j_pairs(
    base_pairs: tuple[tuple[int, int], ...],
    extra_pairs: tuple[tuple[int, int], ...] | None,
) -> tuple[tuple[int, int], ...]:
    merged: list[tuple[int, int]] = list(base_pairs)
    for pair in extra_pairs or ():
        normalized = (int(pair[0]), int(pair[1]))
        if normalized not in merged:
            merged.append(normalized)
    return tuple(merged)


def _rank_branch_traces(
    traces: list[dict[str, object]],
    *,
    peak_key: str,
    total_key: str,
) -> list[dict[str, object]]:
    return sorted(
        traces,
        key=lambda trace: (
            float(trace[peak_key]),
            float(trace[total_key]),
            -int(trace["lower_j"]),
            -int(trace["upper_j"]),
        ),
        reverse=True,
    )


def _branch_label_candidates(
    traces: list[dict[str, object]],
    label_top_n_per_delta_j: int,
    *,
    peak_key: str,
    total_key: str,
    forced_j_pairs: tuple[tuple[int, int], ...],
) -> tuple[dict[int, list[dict[str, object]]], dict[int, set[str]], dict[tuple[int, int], int]]:
    candidates_by_delta_j: dict[int, list[dict[str, object]]] = {delta_j: [] for delta_j in ABSORBANCE_DELTA_J_VALUES}
    labeled_names_by_delta_j: dict[int, set[str]] = {delta_j: set() for delta_j in ABSORBANCE_DELTA_J_VALUES}
    rank_lookup: dict[tuple[int, int], int] = {}

    for delta_j in sorted({int(trace["delta_j"]) for trace in traces}):
        branch_traces = [trace for trace in traces if int(trace["delta_j"]) == delta_j]
        ranked = _rank_branch_traces(branch_traces, peak_key=peak_key, total_key=total_key)
        for rank, trace in enumerate(ranked, start=1):
            rank_lookup[(int(trace["lower_j"]), int(trace["upper_j"]))] = rank
        if delta_j not in ABSORBANCE_DELTA_J_VALUES:
            continue

        selected: list[dict[str, object]] = ranked[: max(0, label_top_n_per_delta_j)]
        selected_names = {str(trace["jpair_label"]) for trace in selected}
        for forced_lower_j, forced_upper_j in forced_j_pairs:
            if _delta_j_value(forced_lower_j, forced_upper_j) != delta_j:
                continue
            forced_label = _format_jpair_label(forced_lower_j, forced_upper_j)
            if forced_label in selected_names:
                continue
            for trace in ranked:
                if int(trace["lower_j"]) == forced_lower_j and int(trace["upper_j"]) == forced_upper_j:
                    selected.append(trace)
                    selected_names.add(forced_label)
                    break

        branch_candidates: list[dict[str, object]] = []
        for trace in selected:
            x_values = np.asarray(trace["wavenumber"], dtype=float)
            y_values = np.asarray(trace["absorbance"], dtype=float)
            peak_index = int(np.argmax(y_values))
            branch_candidates.append(
                {
                    "trace": trace,
                    "peak_x": float(x_values[peak_index]),
                    "peak_y": float(y_values[peak_index]),
                    "text": str(trace["jpair_label"]),
                    "color": str(trace["color"]),
                    "delta_j": delta_j,
                }
            )
        candidates_by_delta_j[delta_j] = branch_candidates
        labeled_names_by_delta_j[delta_j] = {str(trace["jpair_label"]) for trace in selected}

    return candidates_by_delta_j, labeled_names_by_delta_j, rank_lookup


def _summarize_delta_j_counts(traces: list[dict[str, object]]) -> dict[int, int]:
    counts = {delta_j: 0 for delta_j in ABSORBANCE_DELTA_J_VALUES}
    for trace in traces:
        delta_j = int(trace["delta_j"])
        if delta_j in counts and bool(trace["plotted_in_figure"]):
            counts[delta_j] += 1
    return counts


def _format_branch_label_summary(labeled_by_delta_j: dict[int, list[dict[str, object]]]) -> str:
    parts: list[str] = []
    for delta_j in ABSORBANCE_DELTA_J_VALUES:
        labels = ", ".join(item["text"] for item in labeled_by_delta_j.get(delta_j, []))
        parts.append(f"{_delta_j_branch_label(delta_j)}: {labels or 'none'}")
    return "; ".join(parts)


def _color_for_index(index: int, total: int) -> str:
    if total <= 1:
        return "#1f77b4"
    cmap = plt.get_cmap("turbo")
    return mcolors.to_hex(cmap(index / max(1, total - 1)), keep_alpha=False)


def _decimate_series(x_values: np.ndarray, y_values: np.ndarray, max_points: int) -> tuple[np.ndarray, np.ndarray]:
    if max_points <= 0 or len(x_values) <= max_points:
        return x_values, y_values
    indices = np.linspace(0, len(x_values) - 1, max_points, dtype=np.int64)
    return x_values[indices], y_values[indices]


def _label_candidates(
    traces: list[dict[str, object]],
    label_top_n: int,
    *,
    required_labels: tuple[str, ...] = (),
) -> list[dict[str, object]]:
    if label_top_n <= 0 and not required_labels:
        return []
    ranked = sorted(
        traces,
        key=lambda trace: (
            float(trace["peak_absorbance"]),
            float(trace["integrated_absorbance"]),
            -int(trace["lower_j"]),
            -int(trace["upper_j"]),
        ),
        reverse=True,
    )
    selected: list[dict[str, object]] = ranked[: max(0, label_top_n)]
    if required_labels:
        selected_names = {str(trace["jpair_label"]) for trace in selected}
        for label in required_labels:
            if label in selected_names:
                continue
            for trace in ranked:
                if str(trace["jpair_label"]) == label:
                    selected.append(trace)
                    selected_names.add(label)
                    break

    candidates: list[dict[str, object]] = []
    for trace in selected:
        x_values = np.asarray(trace["wavenumber"], dtype=float)
        y_values = np.asarray(trace["absorbance"], dtype=float)
        peak_index = int(np.argmax(y_values))
        candidates.append(
            {
                "trace": trace,
                "peak_x": float(x_values[peak_index]),
                "peak_y": float(y_values[peak_index]),
                "text": str(trace["jpair_label"]),
                "color": str(trace["color"]),
            }
        )
    return candidates


def _label_positions(candidates: list[dict[str, object]], wn_min: float, wn_max: float) -> list[dict[str, object]]:
    if not candidates:
        return []
    ordered = sorted(candidates, key=lambda item: item["peak_y"])
    peak_values = [float(item["peak_y"]) for item in ordered]
    ymin = min(peak_values)
    ymax = max(peak_values)
    if math.isclose(ymin, ymax):
        ymax = ymin * 1.05 if ymin > 0.0 else 1.0
    lower_bound = ymin + 0.05 * (ymax - ymin)
    upper_bound = ymax + 0.08 * (ymax - ymin)
    if len(ordered) == 1:
        target_ys = [peak_values[0]]
    else:
        target_ys = np.linspace(lower_bound, upper_bound, len(ordered))

    anchor_x = wn_min + 0.86 * (wn_max - wn_min)
    positioned: list[dict[str, object]] = []
    for item, target_y in zip(ordered, target_ys):
        positioned.append({**item, "label_x": anchor_x, "label_y": float(target_y)})
    return positioned


def _save_progression_png(
    path: Path,
    *,
    progression_label: str,
    lower_mode: tuple[int, int, int, int],
    upper_mode: tuple[int, int, int, int],
    traces: list[dict[str, object]],
    labeled_traces_by_delta_j: dict[int, list[dict[str, object]]],
    wn_min: float,
    wn_max: float,
) -> Path:
    ensure_directory(path.parent)
    figure, axes = plt.subplots(3, 1, figsize=(16, 12), sharex=True, constrained_layout=True)
    figure.suptitle(f"HITRAN {progression_label} absorbance by J pair, split by delta J\n{_mode_label(lower_mode)} -> {_mode_label(upper_mode)}")
    for axis, delta_j in zip(axes, ABSORBANCE_DELTA_J_VALUES):
        branch_traces = [trace for trace in traces if int(trace["delta_j"]) == delta_j and bool(trace["plotted_in_figure"])]
        for trace in branch_traces:
            axis.plot(trace["wavenumber"], trace["absorbance"], color=trace["color"], linewidth=1.0, alpha=0.95)
        label_items = labeled_traces_by_delta_j.get(delta_j, [])
        for item in label_items:
            axis.plot([item["peak_x"], item["label_x"]], [item["peak_y"], item["label_y"]], color=item["color"], linewidth=0.8, alpha=0.8)
            axis.text(
                item["label_x"],
                item["label_y"],
                item["text"],
                color=item["color"],
                fontsize=8,
                ha="left",
                va="center",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.0},
            )
        axis.set_xlim(wn_min, wn_max)
        y_min, y_max = compute_panel_y_limits(branch_traces, y_key="absorbance", label_items=label_items)
        axis.set_ylim(y_min, y_max)
        axis.set_ylabel("Absorbance")
        axis.set_title(_delta_j_panel_title(delta_j, len(branch_traces)))
        axis.grid(True, alpha=0.25, linewidth=0.5)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        if not branch_traces:
            axis.text(
                0.5,
                0.5,
                "No J pairs in this delta J class",
                transform=axis.transAxes,
                ha="center",
                va="center",
                fontsize=10,
                color="#666666",
            )
    axes[-1].set_xlabel("Wavenumber (cm^-1)")
    figure.savefig(path, dpi=200)
    plt.close(figure)
    return path


def _save_progression_html(
    path: Path,
    *,
    progression_label: str,
    lower_mode: tuple[int, int, int, int],
    upper_mode: tuple[int, int, int, int],
    traces: list[dict[str, object]],
    labeled_traces_by_delta_j: dict[int, list[dict[str, object]]],
    wn_min: float,
    wn_max: float,
    html_max_points: int,
) -> Path:
    ensure_directory(path.parent)
    subplot_titles = [
        _delta_j_panel_title(
            delta_j,
            len([trace for trace in traces if int(trace["delta_j"]) == delta_j and bool(trace["plotted_in_figure"])]),
        )
        for delta_j in ABSORBANCE_DELTA_J_VALUES
    ]
    figure = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.04, subplot_titles=subplot_titles)
    labeled_keys = {
        item["text"]
        for delta_j in ABSORBANCE_DELTA_J_VALUES
        for item in labeled_traces_by_delta_j.get(delta_j, [])
    }
    for row_index, delta_j in enumerate(ABSORBANCE_DELTA_J_VALUES, start=1):
        branch_traces = [trace for trace in traces if int(trace["delta_j"]) == delta_j and bool(trace["plotted_in_figure"])]
        for trace in branch_traces:
            x_plot, y_plot = _decimate_series(
                np.asarray(trace["wavenumber"], dtype=float),
                np.asarray(trace["absorbance"], dtype=float),
                html_max_points,
            )
            figure.add_trace(
                go.Scattergl(
                    x=x_plot,
                    y=y_plot,
                    mode="lines",
                    name=str(trace["jpair_label"]),
                    meta=[str(trace["jpair_label"]), _delta_j_branch_label(delta_j), int(trace["line_count"]), int(trace["source_file_count"])],
                    legendgroup=str(trace["jpair_label"]),
                    showlegend=str(trace["jpair_label"]) in labeled_keys,
                    line={"color": trace["color"], "width": 1.3},
                    hovertemplate=(
                        "J pair: %{meta[0]}<br>"
                        "Branch: %{meta[1]}<br>"
                        "Wavenumber: %{x:.4f} cm^-1<br>"
                        "Absorbance: %{y:.4e}<br>"
                        "Line count: %{meta[2]}<br>"
                        "Source files: %{meta[3]}<extra></extra>"
                    ),
                ),
                row=row_index,
                col=1,
            )
        axis_suffix = "" if row_index == 1 else str(row_index)
        xref = f"x{axis_suffix}"
        yref = f"y{axis_suffix}"
        label_items = labeled_traces_by_delta_j.get(delta_j, [])
        for item in label_items:
            figure.add_annotation(
                x=item["label_x"],
                y=item["label_y"],
                xref=xref,
                yref=yref,
                text=item["text"],
                showarrow=True,
                arrowhead=0,
                arrowcolor=item["color"],
                arrowwidth=1.0,
                axref=xref,
                ayref=yref,
                ax=item["peak_x"],
                ay=item["peak_y"],
                font={"color": item["color"], "size": 11},
                bgcolor="rgba(255,255,255,0.7)",
            )
        y_min, y_max = compute_panel_y_limits(branch_traces, y_key="absorbance", label_items=label_items)
        figure.update_yaxes(
            title_text="Absorbance",
            tickformat=".2e",
            range=[y_min, y_max],
            row=row_index,
            col=1,
        )
        if not branch_traces:
            figure.add_annotation(
                x=0.5,
                y=0.5,
                xref=f"{xref} domain",
                yref=f"{yref} domain",
                text="No J pairs in this delta J class",
                showarrow=False,
                font={"color": "#666666", "size": 11},
            )
    figure.update_layout(
        title=f"HITRAN {progression_label} absorbance by J pair, split by delta J<br><sup>{_mode_label(lower_mode)} -> {_mode_label(upper_mode)}</sup>",
        template="plotly_white",
        xaxis_title="Wavenumber (cm^-1)",
        hovermode="x unified",
        width=1400,
        height=1100,
        legend_title_text="Labeled J pairs",
    )
    figure.update_xaxes(range=[wn_min, wn_max], row=3, col=1)
    figure.write_html(path, include_plotlyjs="cdn")
    return path


def _prepare_runtime_db(*, source_db_dir: Path, runtime_db_dir: Path, source_table: str) -> Path:
    ensure_directory(runtime_db_dir)
    source_header_path = source_db_dir / f"{source_table}.header"
    source_data_path = source_db_dir / f"{source_table}.data"
    if not source_header_path.exists():
        raise FileNotFoundError(f"Missing HITRAN header: {source_header_path}")
    if not source_data_path.exists():
        raise FileNotFoundError(f"Missing HITRAN data: {source_data_path}")
    shutil.copy2(source_header_path, runtime_db_dir / source_header_path.name)
    shutil.copy2(source_data_path, runtime_db_dir / source_data_path.name)
    return runtime_db_dir


def _bootstrap_source_schema(*, source_db_dir: Path, runtime_db_dir: Path, source_table: str) -> Path:
    prepared_runtime_db = _prepare_runtime_db(
        source_db_dir=source_db_dir,
        runtime_db_dir=runtime_db_dir,
        source_table=source_table,
    )
    hapi.db_begin(str(prepared_runtime_db))
    hapi.storage2cache(source_table, nlines=0)
    return prepared_runtime_db


def _temp_table_name(source_table: str, progression_slug: str, lower_j: int, upper_j: int) -> str:
    return f"__research_{source_table}_{progression_slug}_J{lower_j}_{upper_j}__"


def _parse_band_text_groups(
    *,
    input_dir: Path,
    header_path: Path,
    wn_min: float,
    wn_max: float,
    lower_mode_filter: tuple[int, int, int, int] | None = None,
    upper_mode_filter: tuple[int, int, int, int] | None = None,
) -> list[dict[str, object]]:
    positions, widths = load_header_metadata(header_path)
    grouped_files: dict[tuple[tuple[int, int, int, int], tuple[int, int, int, int]], list[dict[str, object]]] = defaultdict(list)
    for path in sorted(input_dir.glob("*.txt")):
        parsed = _parse_band_text_filename(path)
        grouped_files[(parsed["lower_mode"], parsed["upper_mode"])].append(parsed)

    progression_groups: list[dict[str, object]] = []
    for (lower_mode, upper_mode), file_infos in sorted(grouped_files.items(), key=lambda item: (item[0][0][2], item[0][1][2], item[0][0], item[0][1])):
        if lower_mode_filter is not None and lower_mode != lower_mode_filter:
            continue
        if upper_mode_filter is not None and upper_mode != upper_mode_filter:
            continue
        grouped_rows: dict[tuple[int, int], dict[str, object]] = defaultdict(lambda: {"raw_lines": [], "source_files": set()})
        row_count = 0
        for file_info in file_infos:
            path = Path(file_info["path"])
            with path.open("r", encoding="utf-8") as handle:
                for raw_line in handle:
                    line = raw_line.rstrip("\r\n")
                    if not line.strip():
                        continue
                    wavenumber = extract_wavenumber(line, positions, widths)
                    if wavenumber < wn_min or wavenumber > wn_max:
                        continue
                    lower_j = _extract_j_value(extract_field(line, "local_lower_quanta", positions, widths))
                    upper_j = _extract_j_value(extract_field(line, "local_upper_quanta", positions, widths))
                    group = grouped_rows[(lower_j, upper_j)]
                    group["raw_lines"].append(line)
                    group["source_files"].add(path.name)
                    row_count += 1
        progression_groups.append(
            {
                "lower_mode": lower_mode,
                "upper_mode": upper_mode,
                "progression_label": _progression_label(lower_mode, upper_mode),
                "progression_slug": _progression_slug(lower_mode, upper_mode),
                "file_infos": file_infos,
                "grouped_rows": grouped_rows,
                "row_count": row_count,
            }
        )
    return progression_groups


def _build_temp_table_from_lines(temp_table: str, source_table: str, raw_lines: list[str]) -> int:
    hapi.dropTable(temp_table)
    hapi.createTable(temp_table, hapi.getDefaultRowObject(source_table))
    line_count = 0
    for line in raw_lines:
        row_object = hapi.getRowObjectFromString(line, source_table)
        hapi.addRowObject(row_object, temp_table)
        line_count += 1
    hapi.LOCAL_TABLE_CACHE[temp_table]["header"]["number_of_rows"] = line_count
    return line_count


def plot_band_text_absorbance_progressions(
    *,
    input_dir: Path | None = None,
    output_dir: Path | None = None,
    db_dir: Path | None = None,
    source_table: str = "CH4_M6_I1",
    header_path: Path | None = None,
    wn_min: float = 2500.0,
    wn_max: float = 3500.0,
    wn_step: float = 0.01,
    temperature_k: float = 600.0,
    pressure_torr: float = 3.0,
    mole_fraction: float = 0.008,
    path_length_cm: float = 100.0,
    intensity_threshold: float = 1.0e-23,
    label_top_n_per_delta_j: int = 8,
    html_max_points: int = 5000,
    lower_mode_filter: tuple[int, int, int, int] | None = None,
    upper_mode_filter: tuple[int, int, int, int] | None = None,
    forced_j_pairs: tuple[tuple[int, int], ...] | None = None,
) -> BandExportResult:
    if wn_max <= wn_min:
        raise ValueError("wn_max must be greater than wn_min")
    if wn_step <= 0.0:
        raise ValueError("wn_step must be positive")
    if label_top_n_per_delta_j < 0:
        raise ValueError("label_top_n_per_delta_j must be non-negative")

    paths = default_paths()
    resolved_input_dir = (input_dir or band_line_text_dir()).resolve()
    resolved_db_dir = (db_dir or paths.hitran_db_dir).resolve()
    resolved_header_path = (header_path or (resolved_db_dir / f"{source_table}.header")).resolve()
    resolved_output_dir = ensure_directory((output_dir or (paths.artifacts_dir / "hitran_band_text_absorbance")).resolve())
    if not resolved_input_dir.exists():
        raise FileNotFoundError(f"HITRAN band-line text folder not found: {resolved_input_dir}")
    if not resolved_header_path.exists():
        raise FileNotFoundError(f"HITRAN header not found: {resolved_header_path}")
    effective_forced_j_pairs = _merge_forced_j_pairs(DEFAULT_FORCED_HITRAN_J_PAIRS, forced_j_pairs)

    progression_groups = _parse_band_text_groups(
        input_dir=resolved_input_dir,
        header_path=resolved_header_path,
        wn_min=wn_min,
        wn_max=wn_max,
        lower_mode_filter=lower_mode_filter,
        upper_mode_filter=upper_mode_filter,
    )
    prepared_runtime_db = _bootstrap_source_schema(
        source_db_dir=resolved_db_dir,
        runtime_db_dir=(resolved_output_dir / "_runtime_db").resolve(),
        source_table=source_table,
    )
    match = TABLE_ID_RE.fullmatch(source_table)
    if match is None:
        raise ValueError(f"Could not parse molecule/isotopologue IDs from {source_table}")
    molecule_id = int(match.group(1))
    isotopologue_id = int(match.group(2))
    case = GasCase(temperature_k=temperature_k, pressure_torr=pressure_torr, mole_fraction=mole_fraction, path_length_cm=path_length_cm)
    window = SpectralWindow(wn_min=wn_min, wn_max=wn_max, wn_step=wn_step)
    diluent = {"self": case.mole_fraction, "He": 1.0 - case.mole_fraction}

    summary_rows: list[dict[str, object]] = []
    for progression in progression_groups:
        traces: list[dict[str, object]] = []
        for index, ((lower_j, upper_j), payload) in enumerate(sorted(progression["grouped_rows"].items())):
            temp_table = _temp_table_name(source_table, str(progression["progression_slug"]), lower_j, upper_j)
            try:
                line_count = _build_temp_table_from_lines(temp_table, source_table, list(payload["raw_lines"]))
                wavenumber, coefficient = hapi.absorptionCoefficient_Voigt(
                    Components=[(molecule_id, isotopologue_id, case.mole_fraction)],
                    SourceTables=[temp_table],
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
            finally:
                hapi.dropTable(temp_table)

            x_values = np.asarray(wavenumber, dtype=np.float64)
            traces.append(
                {
                    "trace_index": index + 1,
                    "lower_j": lower_j,
                    "upper_j": upper_j,
                    "delta_j": _delta_j_value(lower_j, upper_j),
                    "branch_label": _delta_j_branch_label(_delta_j_value(lower_j, upper_j)),
                    "jpair_label": _format_jpair_label(lower_j, upper_j),
                    "wavenumber": x_values,
                    "absorbance": absorbance,
                    "line_count": line_count,
                    "grid_point_count": int(x_values.size),
                    "peak_absorbance": float(np.max(absorbance)) if len(absorbance) else 0.0,
                    "integrated_absorbance": float(np.trapezoid(absorbance, x_values)) if len(absorbance) else 0.0,
                    "source_file_count": len(payload["source_files"]),
                    "plotted_in_figure": _delta_j_value(lower_j, upper_j) in ABSORBANCE_DELTA_J_VALUES,
                }
            )

        if not traces:
            continue

        for index, trace in enumerate(traces):
            trace["color"] = _color_for_index(index, len(traces))

        (
            labeled_candidates_by_delta_j,
            labeled_names_by_delta_j,
            branch_rank_lookup,
        ) = _branch_label_candidates(
            traces,
            label_top_n_per_delta_j,
            peak_key="peak_absorbance",
            total_key="integrated_absorbance",
            forced_j_pairs=effective_forced_j_pairs,
        )
        labeled_traces_by_delta_j = {
            delta_j: _label_positions(labeled_candidates_by_delta_j.get(delta_j, []), wn_min, wn_max)
            for delta_j in ABSORBANCE_DELTA_J_VALUES
        }
        labeled_names = {
            label
            for branch_names in labeled_names_by_delta_j.values()
            for label in branch_names
        }

        ranked_for_mapping = sorted(
            traces,
            key=lambda trace: (
                float(trace["peak_absorbance"]),
                float(trace["integrated_absorbance"]),
                -int(trace["lower_j"]),
                -int(trace["upper_j"]),
            ),
            reverse=True,
        )
        rank_lookup = {
            (int(trace["lower_j"]), int(trace["upper_j"])): rank
            for rank, trace in enumerate(ranked_for_mapping, start=1)
        }
        mapping_rows: list[dict[str, object]] = []
        for trace in traces:
            trace_key = (int(trace["lower_j"]), int(trace["upper_j"]))
            mapping_rows.append(
                {
                    "trace_index": trace["trace_index"],
                    "strength_rank": rank_lookup[(int(trace["lower_j"]), int(trace["upper_j"]))],
                    "delta_j": trace["delta_j"],
                    "branch_label": trace["branch_label"],
                    "delta_j_strength_rank": branch_rank_lookup.get(trace_key, 0),
                    "lower_j": trace["lower_j"],
                    "upper_j": trace["upper_j"],
                    "jpair_label": trace["jpair_label"],
                    "line_count": trace["line_count"],
                    "grid_point_count": trace["grid_point_count"],
                    "peak_absorbance": f"{float(trace['peak_absorbance']):.12e}",
                    "integrated_absorbance": f"{float(trace['integrated_absorbance']):.12e}",
                    "color_hex": trace["color"],
                    "source_file_count": trace["source_file_count"],
                    "plotted_in_figure": "yes" if bool(trace["plotted_in_figure"]) else "no",
                    "labeled_on_figure": "yes" if trace["jpair_label"] in labeled_names else "no",
                }
            )

        branch_counts = _summarize_delta_j_counts(traces)
        skipped_traces = [trace for trace in traces if not bool(trace["plotted_in_figure"])]
        labeled_summary = _format_branch_label_summary(labeled_traces_by_delta_j)
        png_path = _save_progression_png(
            resolved_output_dir / f"{progression['progression_slug']}_absorbance.png",
            progression_label=str(progression["progression_label"]),
            lower_mode=progression["lower_mode"],
            upper_mode=progression["upper_mode"],
            traces=traces,
            labeled_traces_by_delta_j=labeled_traces_by_delta_j,
            wn_min=wn_min,
            wn_max=wn_max,
        )
        html_path = _save_progression_html(
            resolved_output_dir / f"{progression['progression_slug']}_absorbance.html",
            progression_label=str(progression["progression_label"]),
            lower_mode=progression["lower_mode"],
            upper_mode=progression["upper_mode"],
            traces=traces,
            labeled_traces_by_delta_j=labeled_traces_by_delta_j,
            wn_min=wn_min,
            wn_max=wn_max,
            html_max_points=html_max_points,
        )
        mapping_path = write_rows_csv(
            resolved_output_dir / f"{progression['progression_slug']}_jpairs.csv",
            mapping_rows,
            fieldnames=[
                "trace_index",
                "strength_rank",
                "delta_j",
                "branch_label",
                "delta_j_strength_rank",
                "lower_j",
                "upper_j",
                "jpair_label",
                "line_count",
                "grid_point_count",
                "peak_absorbance",
                "integrated_absorbance",
                "color_hex",
                "source_file_count",
                "plotted_in_figure",
                "labeled_on_figure",
            ],
        )

        summary_rows.append(
            {
                "progression_label": progression["progression_label"],
                "lower_mode": _mode_label(progression["lower_mode"]),
                "upper_mode": _mode_label(progression["upper_mode"]),
                "source_file_count": len(progression["file_infos"]),
                "row_count": progression["row_count"],
                "jpair_count": len(traces),
                "grid_point_count": traces[0]["grid_point_count"],
                "delta_j_minus_1_count": branch_counts[-1],
                "delta_j_0_count": branch_counts[0],
                "delta_j_plus_1_count": branch_counts[1],
                "skipped_jpair_count": len(skipped_traces),
                "labeled_j_pairs": labeled_summary,
                "png_file": png_path.name,
                "html_file": html_path.name,
                "mapping_csv": mapping_path.name,
            }
        )

    summary_csv_path = write_rows_csv(
        resolved_output_dir / "progression_summary.csv",
        summary_rows,
        fieldnames=[
            "progression_label",
            "lower_mode",
            "upper_mode",
            "source_file_count",
            "row_count",
            "jpair_count",
            "grid_point_count",
            "delta_j_minus_1_count",
            "delta_j_0_count",
            "delta_j_plus_1_count",
            "skipped_jpair_count",
            "labeled_j_pairs",
            "png_file",
            "html_file",
            "mapping_csv",
        ],
    )

    report_lines = [
        "# HITRAN Band Text Absorbance Progressions",
        "",
        f"- Input folder: `{resolved_input_dir}`",
        f"- Source table schema: `{source_table}`",
        f"- Isolated runtime DB: `{prepared_runtime_db}`",
        f"- Wavenumber window: `{wn_min:g}` to `{wn_max:g} cm^-1` with `step = {wn_step:g} cm^-1`",
        "- Y axis: `absorbance`",
        "- Curve grouping: full J pair `(lower J, upper J)`, split into `delta J = -1, 0, +1` panels",
        f"- Temperature: `{temperature_k:g} K`",
        f"- Pressure: `{pressure_torr:g} Torr`",
        f"- Mole fraction: `{mole_fraction:g}`",
        f"- Path length: `{path_length_cm:g} cm`",
        f"- HAPI intensity threshold: `{intensity_threshold:.3e}`",
        f"- On-figure labels: strongest `{label_top_n_per_delta_j}` J pairs per `delta J` panel, plus forced labels for `{', '.join(_format_jpair_label(lower_j, upper_j) for lower_j, upper_j in effective_forced_j_pairs)}` when present",
        f"- Summary CSV: [{summary_csv_path.name}]({summary_csv_path.name})",
        "",
    ]
    for row in summary_rows:
        report_lines.extend(
            [
                f"## {row['progression_label']}",
                "",
                f"- Modes: `{row['lower_mode']} -> {row['upper_mode']}`",
                f"- Files merged: `{row['source_file_count']}`",
                f"- HITRAN rows used: `{row['row_count']}`",
                f"- J-pair curves: `{row['jpair_count']}`",
                f"- Grid points per curve: `{row['grid_point_count']}`",
                f"- Plotted branch counts: `delta J=-1: {row['delta_j_minus_1_count']}`, `delta J=0: {row['delta_j_0_count']}`, `delta J=+1: {row['delta_j_plus_1_count']}`",
                f"- Skipped J pairs outside plotted branches: `{row['skipped_jpair_count']}`",
                f"- Labeled J pairs by branch: `{row['labeled_j_pairs'] or 'none'}`",
                f"- Outputs: [PNG]({row['png_file']}), [HTML]({row['html_file']}), [J-pair CSV]({row['mapping_csv']})",
                "",
                f"![{row['progression_label']}]({row['png_file']})",
                "",
            ]
        )
    report_path = write_markdown(resolved_output_dir / "report.md", "\n".join(report_lines))
    return BandExportResult(
        rows=summary_rows,
        output_dir=resolved_output_dir,
        manifest_path=report_path,
        count=len(summary_rows),
        metadata={
            "input_dir": resolved_input_dir,
            "db_dir": resolved_db_dir,
            "header_path": resolved_header_path,
            "summary_csv_path": summary_csv_path,
            "report_path": report_path,
            "case": case,
            "window": window,
        },
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
