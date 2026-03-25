"""
Combined pure-nu3 absorbance workflow for ExoMol sorted exports and HITRAN band texts.

Workflow structure:
- collect pure-nu3 progression groups from ExoMol, when enabled
- collect pure-nu3 progression groups from HITRAN, when enabled
- build the J-pair union for each pure-nu3 progression
- choose one source per J pair using the ordered `source_priority`
- render the chosen J pair with that source's native absorbance path
- write one PNG, one HTML, one J-pair CSV, one summary CSV, and one report
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import hapi
import numpy as np
import plotly.graph_objects as go
from matplotlib import pyplot as plt
from plotly.subplots import make_subplots

from .exomol import (
    ABSORBANCE_DELTA_J_VALUES,
    DATASET_STEM,
    DEFAULT_FORCED_ABSORBANCE_J_PAIRS,
    _branch_label_candidates,
    _collect_sorted_progression_groups,
    _color_for_index,
    _delta_j_value,
    _delta_j_panel_title,
    _format_branch_label_summary,
    _format_jpair_label,
    _label_positions,
    _merge_forced_j_pairs,
    _mode_label,
    _summarize_delta_j_counts,
    dataset_dir,
    parse_def_file,
    render_absorbance_on_grid,
    sorted_nu3_band_dir,
)
from .hitran import (
    TABLE_ID_RE,
    _bootstrap_source_schema,
    _build_temp_table_from_lines,
    _parse_band_text_groups,
    _temp_table_name,
    band_line_text_dir,
)
from .io import default_paths, ensure_directory, write_markdown, write_rows_csv
from .models import GasCase, SummaryResult, SpectralWindow
from .spectra import build_grid, compute_panel_y_limits, to_absorbance


SOURCE_CHOICES = ("exomol", "hitran")
DEFAULT_SOURCE_PRIORITY = ("hitran", "exomol")
PURE_NU3_PROGRESSIONS = tuple(
    ((0, 0, lower_nu3, 0), (0, 0, lower_nu3 + 1, 0))
    for lower_nu3 in range(4)
)


def _combined_progression_label(lower_mode: tuple[int, int, int, int], upper_mode: tuple[int, int, int, int]) -> str:
    return f"nu3 {lower_mode[2]}->{upper_mode[2]}"


def _combined_progression_slug(lower_mode: tuple[int, int, int, int], upper_mode: tuple[int, int, int, int]) -> str:
    return f"nu3_{lower_mode[2]}_to_{upper_mode[2]}"


def _source_label(source_name: str) -> str:
    return "ExoMol" if source_name == "exomol" else "HITRAN"


def _normalize_source_priority(
    normalized_sources: tuple[str, ...],
    source_priority: tuple[str, ...],
) -> tuple[str, ...]:
    loaded_sources = set(normalized_sources)
    normalized_priority: list[str] = []
    for source_name in source_priority:
        if source_name in loaded_sources and source_name not in normalized_priority:
            normalized_priority.append(source_name)
    for source_name in normalized_sources:
        if source_name not in normalized_priority:
            normalized_priority.append(source_name)
    return tuple(normalized_priority)


def _collect_exomol_pure_nu3_groups(
    input_dir: Path,
    *,
    wn_min: float,
    wn_max: float,
    min_line_intensity: float,
) -> dict[tuple[tuple[int, int, int, int], tuple[int, int, int, int]], dict[str, object]]:
    progression_groups = _collect_sorted_progression_groups(
        input_dir,
        wn_min=wn_min,
        wn_max=wn_max,
        min_line_intensity=min_line_intensity,
    )
    return {
        (tuple(group["lower_mode"]), tuple(group["upper_mode"])): group
        for group in progression_groups
        if (tuple(group["lower_mode"]), tuple(group["upper_mode"])) in PURE_NU3_PROGRESSIONS
    }


def _collect_hitran_pure_nu3_groups(
    input_dir: Path,
    *,
    header_path: Path,
    wn_min: float,
    wn_max: float,
) -> dict[tuple[tuple[int, int, int, int], tuple[int, int, int, int]], dict[str, object]]:
    progression_groups = _parse_band_text_groups(
        input_dir=input_dir,
        header_path=header_path,
        wn_min=wn_min,
        wn_max=wn_max,
    )
    return {
        (tuple(group["lower_mode"]), tuple(group["upper_mode"])): group
        for group in progression_groups
        if (tuple(group["lower_mode"]), tuple(group["upper_mode"])) in PURE_NU3_PROGRESSIONS
    }


def _save_combined_progression_png(
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
    figure.suptitle(
        f"{progression_label} combined absorbance by J pair, split by delta J\n"
        f"{_mode_label(lower_mode)} -> {_mode_label(upper_mode)}"
    )
    for axis, delta_j in zip(axes, ABSORBANCE_DELTA_J_VALUES):
        branch_traces = [trace for trace in traces if int(trace["delta_j"]) == delta_j and bool(trace["plotted_in_figure"])]
        label_items = labeled_traces_by_delta_j.get(delta_j, [])
        for trace in branch_traces:
            axis.plot(trace["wavenumber"], trace["absorbance"], color=trace["color"], linewidth=1.0, alpha=0.95)
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


def _save_combined_progression_html(
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
            x_values = np.asarray(trace["wavenumber"], dtype=float)
            y_values = np.asarray(trace["absorbance"], dtype=float)
            if html_max_points > 0 and x_values.size > html_max_points:
                indices = np.linspace(0, x_values.size - 1, html_max_points, dtype=np.int64)
                x_values = x_values[indices]
                y_values = y_values[indices]
            figure.add_trace(
                go.Scattergl(
                    x=x_values,
                    y=y_values,
                    mode="lines",
                    name=str(trace["jpair_label"]),
                    meta=[
                        str(trace["jpair_label"]),
                        str(trace["selected_source"]),
                        str(trace["available_sources_text"]),
                        str(trace["fallback_used_text"]),
                        int(trace["line_count"]),
                        int(trace["source_file_count"]),
                    ],
                    legendgroup=str(trace["jpair_label"]),
                    showlegend=str(trace["jpair_label"]) in labeled_keys,
                    line={"color": trace["color"], "width": 1.3},
                    hovertemplate=(
                        "J pair: %{meta[0]}<br>"
                        "Selected source: %{meta[1]}<br>"
                        "Available sources: %{meta[2]}<br>"
                        "Fallback used: %{meta[3]}<br>"
                        "Wavenumber: %{x:.4f} cm^-1<br>"
                        "Absorbance: %{y:.4e}<br>"
                        "Line count: %{meta[4]}<br>"
                        "Source files: %{meta[5]}<extra></extra>"
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
        title=(
            f"{progression_label} combined absorbance by J pair"
            f"<br><sup>{_mode_label(lower_mode)} -> {_mode_label(upper_mode)}</sup>"
        ),
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


def plot_combined_pure_nu3_absorbance_progressions(
    *,
    exomol_input_dir: Path | None = None,
    hitran_input_dir: Path | None = None,
    output_dir: Path | None = None,
    hitran_db_dir: Path | None = None,
    hitran_header_path: Path | None = None,
    source_table: str = "CH4_M6_I1",
    sources: tuple[str, ...] = SOURCE_CHOICES,
    source_priority: tuple[str, ...] = DEFAULT_SOURCE_PRIORITY,
    wn_min: float = 2500.0,
    wn_max: float = 3500.0,
    wn_step: float = 0.25,
    pressure_torr: float = 3.0,
    mole_fraction: float = 0.008,
    path_length_cm: float = 100.0,
    line_cutoff: float = 0.5,
    min_line_intensity: float = 0.0,
    hitran_intensity_threshold: float = 1.0e-23,
    label_top_n_per_delta_j: int = 8,
    html_max_points: int = 5000,
    forced_j_pairs: tuple[tuple[int, int], ...] | None = None,
) -> SummaryResult:
    if wn_max <= wn_min:
        raise ValueError("wn_max must be greater than wn_min")
    if wn_step <= 0.0:
        raise ValueError("wn_step must be positive")
    if line_cutoff <= 0.0:
        raise ValueError("line_cutoff must be positive")
    if label_top_n_per_delta_j < 0:
        raise ValueError("label_top_n_per_delta_j must be non-negative")

    normalized_sources = tuple(source for source in SOURCE_CHOICES if source in set(sources))
    if not normalized_sources:
        raise ValueError("At least one source must be selected from exomol/hitran")
    normalized_source_priority = _normalize_source_priority(normalized_sources, source_priority)

    paths = default_paths()
    resolved_exomol_input_dir = (exomol_input_dir or sorted_nu3_band_dir()).resolve()
    resolved_hitran_input_dir = (hitran_input_dir or band_line_text_dir()).resolve()
    resolved_output_dir = ensure_directory((output_dir or (paths.artifacts_dir / "combined_pure_nu3_absorbance")).resolve())
    resolved_hitran_db_dir = (hitran_db_dir or paths.hitran_db_dir).resolve()
    resolved_hitran_header_path = (hitran_header_path or (resolved_hitran_db_dir / f"{source_table}.header")).resolve()
    effective_forced_j_pairs = _merge_forced_j_pairs(DEFAULT_FORCED_ABSORBANCE_J_PAIRS, forced_j_pairs)

    exomol_groups: dict[tuple[tuple[int, int, int, int], tuple[int, int, int, int]], dict[str, object]] = {}
    hitran_groups: dict[tuple[tuple[int, int, int, int], tuple[int, int, int, int]], dict[str, object]] = {}

    if "exomol" in normalized_sources:
        if not resolved_exomol_input_dir.exists():
            raise FileNotFoundError(f"Sorted ExoMol folder not found: {resolved_exomol_input_dir}")
        exomol_groups = _collect_exomol_pure_nu3_groups(
            resolved_exomol_input_dir,
            wn_min=wn_min,
            wn_max=wn_max,
            min_line_intensity=min_line_intensity,
        )

    if "hitran" in normalized_sources:
        if not resolved_hitran_input_dir.exists():
            raise FileNotFoundError(f"HITRAN band-line text folder not found: {resolved_hitran_input_dir}")
        if not resolved_hitran_header_path.exists():
            raise FileNotFoundError(f"HITRAN header not found: {resolved_hitran_header_path}")
        hitran_groups = _collect_hitran_pure_nu3_groups(
            resolved_hitran_input_dir,
            header_path=resolved_hitran_header_path,
            wn_min=wn_min,
            wn_max=wn_max,
        )

    metadata = parse_def_file(dataset_dir() / f"{DATASET_STEM}.def")
    case = GasCase(
        temperature_k=296.0,
        pressure_torr=pressure_torr,
        mole_fraction=mole_fraction,
        path_length_cm=path_length_cm,
    )
    window = SpectralWindow(wn_min=wn_min, wn_max=wn_max, wn_step=wn_step)
    grid = build_grid(window)

    prepared_runtime_db: Path | None = None
    molecule_id: int | None = None
    isotopologue_id: int | None = None
    diluent = {"self": case.mole_fraction, "He": 1.0 - case.mole_fraction}
    if "hitran" in normalized_sources:
        prepared_runtime_db = _bootstrap_source_schema(
            source_db_dir=resolved_hitran_db_dir,
            runtime_db_dir=(resolved_output_dir / "_runtime_db").resolve(),
            source_table=source_table,
        )
        match = TABLE_ID_RE.fullmatch(source_table)
        if match is None:
            raise ValueError(f"Could not parse molecule/isotopologue IDs from {source_table}")
        molecule_id = int(match.group(1))
        isotopologue_id = int(match.group(2))

    summary_rows: list[dict[str, object]] = []
    for lower_mode, upper_mode in PURE_NU3_PROGRESSIONS:
        progression_key = (lower_mode, upper_mode)
        progression_label = _combined_progression_label(lower_mode, upper_mode)
        progression_slug = _combined_progression_slug(lower_mode, upper_mode)
        exomol_payload = exomol_groups.get(progression_key)
        hitran_payload = hitran_groups.get(progression_key)

        available_jpairs = set()
        if exomol_payload is not None:
            available_jpairs.update(exomol_payload["grouped_rows"].keys())
        if hitran_payload is not None:
            available_jpairs.update(hitran_payload["grouped_rows"].keys())
        if not available_jpairs:
            continue

        traces: list[dict[str, object]] = []
        selected_counts = {source_name: 0 for source_name in SOURCE_CHOICES}
        for index, (lower_j, upper_j) in enumerate(sorted(available_jpairs)):
            available_payloads: dict[str, dict[str, object]] = {}
            selected_source: str | None = None
            selected_payload: dict[str, object] | None = None
            render_mode: str | None = None

            if exomol_payload is not None and (lower_j, upper_j) in exomol_payload["grouped_rows"]:
                available_payloads["exomol"] = exomol_payload["grouped_rows"][(lower_j, upper_j)]
            if hitran_payload is not None and (lower_j, upper_j) in hitran_payload["grouped_rows"]:
                available_payloads["hitran"] = hitran_payload["grouped_rows"][(lower_j, upper_j)]
            available_sources = [
                source_name for source_name in normalized_source_priority if source_name in available_payloads
            ]
            if available_sources:
                selected_source = available_sources[0]
                selected_payload = available_payloads[selected_source]
                render_mode = selected_source

            if selected_source is None or selected_payload is None or render_mode is None:
                continue

            if selected_source == "exomol":
                selected_counts["exomol"] += 1
                line_centers = np.asarray(selected_payload["wavenumber"], dtype=np.float64)
                line_intensities = np.asarray(selected_payload["intensity"], dtype=np.float64)
                absorbance = render_absorbance_on_grid(
                    grid,
                    case=case,
                    line_centers=line_centers,
                    line_intensities=line_intensities,
                    mass_da=float(metadata["mass_da"]),
                    gamma0=float(metadata["gamma0"]),
                    n_exponent=float(metadata["n_exponent"]),
                    line_cutoff=line_cutoff,
                )
                line_count = int(line_centers.size)
                source_file_count = len(selected_payload["source_files"])
                x_values = grid
            else:
                selected_counts["hitran"] += 1
                if molecule_id is None or isotopologue_id is None:
                    raise RuntimeError("HITRAN runtime database was not initialized")
                temp_table = _temp_table_name(source_table, progression_slug, lower_j, upper_j)
                try:
                    line_count = _build_temp_table_from_lines(temp_table, source_table, list(selected_payload["raw_lines"]))
                    wavenumber, coefficient = hapi.absorptionCoefficient_Voigt(
                        Components=[(molecule_id, isotopologue_id, case.mole_fraction)],
                        SourceTables=[temp_table],
                        WavenumberRange=[window.wn_min, window.wn_max],
                        WavenumberStep=window.wn_step,
                        Environment={"T": case.temperature_k, "p": case.pressure_atm},
                        Diluent=diluent,
                        IntensityThreshold=hitran_intensity_threshold,
                        HITRAN_units=False,
                    )
                    _, transmittance = hapi.transmittanceSpectrum(
                        wavenumber,
                        coefficient,
                        Environment={"l": case.path_length_cm},
                    )
                    absorbance = to_absorbance(np.asarray(transmittance))
                    x_values = np.asarray(wavenumber, dtype=np.float64)
                finally:
                    hapi.dropTable(temp_table)
                source_file_count = len(selected_payload["source_files"])

            traces.append(
                {
                    "trace_index": index + 1,
                    "lower_j": lower_j,
                    "upper_j": upper_j,
                    "delta_j": _delta_j_value(lower_j, upper_j),
                    "branch_label": f"dJ_{_delta_j_value(lower_j, upper_j):+d}",
                    "jpair_label": _format_jpair_label(lower_j, upper_j),
                    "wavenumber": x_values,
                    "absorbance": absorbance,
                    "grid_point_count": int(np.asarray(x_values).size),
                    "line_count": line_count,
                    "peak_absorbance": float(np.max(absorbance)) if len(absorbance) else 0.0,
                    "integrated_absorbance": float(np.trapezoid(absorbance, x_values)) if len(absorbance) else 0.0,
                    "source_file_count": source_file_count,
                    "selected_source": selected_source,
                    "selected_source_label": _source_label(selected_source),
                    "available_sources": tuple(available_sources),
                    "available_sources_text": ", ".join(_source_label(source_name) for source_name in available_sources),
                    "fallback_used": bool(normalized_source_priority) and selected_source != normalized_source_priority[0],
                    "fallback_used_text": (
                        "yes"
                        if bool(normalized_source_priority) and selected_source != normalized_source_priority[0]
                        else "no"
                    ),
                    "render_mode": render_mode,
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
            y_key="absorbance",
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
                    "strength_rank": rank_lookup[trace_key],
                    "delta_j": trace["delta_j"],
                    "branch_label": trace["branch_label"],
                    "delta_j_strength_rank": branch_rank_lookup.get(trace_key, 0),
                    "lower_j": trace["lower_j"],
                    "upper_j": trace["upper_j"],
                    "jpair_label": trace["jpair_label"],
                    "selected_source": trace["selected_source"],
                    "available_sources": ", ".join(trace["available_sources"]),
                    "fallback_used": "yes" if trace["fallback_used"] else "no",
                    "render_mode": trace["render_mode"],
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
        png_path = _save_combined_progression_png(
            resolved_output_dir / f"{progression_slug}_absorbance.png",
            progression_label=progression_label,
            lower_mode=lower_mode,
            upper_mode=upper_mode,
            traces=traces,
            labeled_traces_by_delta_j=labeled_traces_by_delta_j,
            wn_min=wn_min,
            wn_max=wn_max,
        )
        html_path = _save_combined_progression_html(
            resolved_output_dir / f"{progression_slug}_absorbance.html",
            progression_label=progression_label,
            lower_mode=lower_mode,
            upper_mode=upper_mode,
            traces=traces,
            labeled_traces_by_delta_j=labeled_traces_by_delta_j,
            wn_min=wn_min,
            wn_max=wn_max,
            html_max_points=html_max_points,
        )
        mapping_path = write_rows_csv(
            resolved_output_dir / f"{progression_slug}_jpairs.csv",
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
                "selected_source",
                "available_sources",
                "fallback_used",
                "render_mode",
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
                "progression_label": progression_label,
                "lower_mode": _mode_label(lower_mode),
                "upper_mode": _mode_label(upper_mode),
                "hitran_selected_jpair_count": selected_counts["hitran"],
                "exomol_selected_jpair_count": selected_counts["exomol"],
                "jpair_count": len(traces),
                "grid_point_count": int(max(int(trace["grid_point_count"]) for trace in traces)),
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
            "hitran_selected_jpair_count",
            "exomol_selected_jpair_count",
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
        "# Combined Pure nu3 Absorbance Progressions",
        "",
        f"- ExoMol input folder: `{resolved_exomol_input_dir}`" if "exomol" in normalized_sources else "- ExoMol input folder: not used",
        f"- HITRAN input folder: `{resolved_hitran_input_dir}`" if "hitran" in normalized_sources else "- HITRAN input folder: not used",
        f"- HITRAN source table schema: `{source_table}`" if "hitran" in normalized_sources else "- HITRAN source table schema: not used",
        f"- Wavenumber window: `{wn_min:g}` to `{wn_max:g} cm^-1` with `step = {wn_step:g} cm^-1`",
        "- Y axis: `absorbance`",
        "- Progression family: pure `nu3` only with `n1 = n2 = n4 = 0`",
        "- J-pair source priority: "
        + " -> ".join(f"`{_source_label(source_name)}`" for source_name in normalized_source_priority),
        f"- Pressure: `{pressure_torr:g} Torr`",
        f"- Mole fraction: `{mole_fraction:g}`",
        f"- Path length: `{path_length_cm:g} cm`",
        f"- ExoMol broadening cutoff: `{line_cutoff:g} cm^-1`",
        f"- ExoMol minimum line intensity kept: `{min_line_intensity:.3e} cm/molecule`",
        f"- HITRAN intensity threshold: `{hitran_intensity_threshold:.3e}`",
        f"- On-figure labels: strongest `{label_top_n_per_delta_j}` J pairs per `delta J` panel, plus forced labels for `{', '.join(_format_jpair_label(lower_j, upper_j) for lower_j, upper_j in effective_forced_j_pairs)}` when present",
        f"- HTML traces are decimated to at most `{html_max_points}` points per J pair for responsiveness",
        f"- Summary CSV: [{summary_csv_path.name}]({summary_csv_path.name})",
        "",
    ]
    for row in summary_rows:
        report_lines.extend(
            [
                f"## {row['progression_label']}",
                "",
                f"- Modes: `{row['lower_mode']} -> {row['upper_mode']}`",
                f"- Selected J-pair curves: `{row['jpair_count']}`",
                f"- Selected from HITRAN: `{row['hitran_selected_jpair_count']}`",
                f"- Selected from ExoMol: `{row['exomol_selected_jpair_count']}`",
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
    return SummaryResult(
        rows=summary_rows,
        csv_path=summary_csv_path,
        metadata={
            "output_dir": resolved_output_dir,
            "report_path": report_path,
            "sources": normalized_sources,
            "source_priority": normalized_source_priority,
            "hitran_runtime_db": prepared_runtime_db,
        },
    )
