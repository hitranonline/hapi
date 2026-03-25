from __future__ import annotations

import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
from plotly.subplots import make_subplots

from .absorbance import render_absorbance_on_grid
from .io import default_paths, ensure_directory, iter_bz2_text_lines, write_markdown, write_rows_csv
from .models import BandSelection, GasCase, SpectrumResult, SpectralWindow, SummaryResult
from .spectra import (
    BAR_PER_ATM,
    build_grid,
    compute_panel_y_limits,
    cross_section_to_absorbance,
    doppler_hwhm_cm,
    save_spectrum_csv,
    save_spectrum_html,
    voigt_profile_cm,
)


DATASET_STEM = "12C-1H4__MM"
SECOND_RADIATION_CONSTANT_CM_K = 1.438776877
LIGHT_SPEED_CM_S = 2.99792458e10
T_REF_K = 296.0
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
SORTED_NU3_FOLDER_NAME = "exomol_ch4_mm_pure_nu3_band_texts_hitran_style_2500_3500_sorted"
SORTED_NU3_FILENAME_PATTERN = re.compile(
    r"^EXOMOL_CH4_MM_pure_nu3_hitran_style_"
    r"(?P<lower_n1>\d+)_(?P<lower_n2>\d+)_(?P<lower_n3>\d+)_(?P<lower_n4>\d+)_(?P<lower_sym>[^_]+)"
    r"_to_"
    r"(?P<upper_n1>\d+)_(?P<upper_n2>\d+)_(?P<upper_n3>\d+)_(?P<upper_n4>\d+)_(?P<upper_sym>[^_]+)"
    r"_T(?P<temperature>[0-9.]+)K\.txt$"
)
J_VALUE_PATTERN = re.compile(r"(\d+)")
REQUIRED_ABSORBANCE_LABELS = ("J 2->3", "J 3->4")
ABSORBANCE_DELTA_J_VALUES = (-1, 0, 1)
DEFAULT_FORCED_ABSORBANCE_J_PAIRS = ((2, 3), (3, 4))


def dataset_dir(data_dir: Path | None = None) -> Path:
    paths = default_paths()
    return data_dir or (paths.exomol_db_dir / "CH4" / "12C-1H4" / "MM")


def transition_filename(start_cm: int, *, dataset_stem: str = DATASET_STEM) -> str:
    return f"{dataset_stem}__{start_cm:05d}-{start_cm + 100:05d}.trans.bz2"


def available_transition_starts(data_dir: Path, *, dataset_stem: str = DATASET_STEM) -> list[int]:
    starts: list[int] = []
    pattern = re.compile(rf"^{re.escape(dataset_stem)}__(\d{{5}})-(\d{{5}})\.trans\.bz2$")
    for path in data_dir.glob(f"{dataset_stem}__*.trans.bz2"):
        match = pattern.match(path.name)
        if match is None:
            continue
        start_cm = int(match.group(1))
        end_cm = int(match.group(2))
        if end_cm == start_cm + 100:
            starts.append(start_cm)
    if not starts:
        raise FileNotFoundError(f"No ExoMol transition files found in {data_dir}")
    return sorted(starts)


def overlapping_transition_files(
    data_dir: Path,
    *,
    wn_min: float,
    wn_max: float,
    wing: float = 0.0,
    dataset_stem: str = DATASET_STEM,
) -> tuple[list[Path], float]:
    available_starts = available_transition_starts(data_dir, dataset_stem=dataset_stem)
    available_min = available_starts[0]
    available_max = available_starts[-1] + 100

    start_chunk = int(math.floor(max(0.0, wn_min - wing) / 100.0) * 100)
    end_chunk = int(math.ceil((wn_max + wing) / 100.0) * 100)
    clipped_start = max(start_chunk, available_min)
    clipped_end = min(end_chunk, available_max)
    if clipped_end <= clipped_start:
        raise FileNotFoundError(
            f"No available transition files overlap {wn_min:g}-{wn_max:g} cm^-1 in {data_dir}"
        )

    effective_wing = wing
    if clipped_start > start_chunk:
        effective_wing = min(effective_wing, max(0.0, wn_min - clipped_start))
    if clipped_end < end_chunk:
        effective_wing = min(effective_wing, max(0.0, clipped_end - wn_max))

    files: list[Path] = []
    for start_cm in range(clipped_start, clipped_end, 100):
        path = data_dir / transition_filename(start_cm, dataset_stem=dataset_stem)
        if not path.exists():
            raise FileNotFoundError(f"Missing transition file: {path}")
        files.append(path)
    return files, effective_wing


def parse_def_file(def_path: Path) -> dict[str, float]:
    metadata: dict[str, float] = {}
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
    required = {"mass_da", "nstates", "gamma0", "n_exponent"}
    missing = required - metadata.keys()
    if missing:
        raise RuntimeError(f"Missing required metadata in {def_path}: {sorted(missing)}")
    return metadata


def load_partition_function(pf_path: Path) -> tuple[np.ndarray, np.ndarray]:
    pf_table = np.loadtxt(pf_path)
    if pf_table.ndim != 2 or pf_table.shape[1] < 2:
        raise RuntimeError(f"Unexpected partition-function format in {pf_path}")
    return pf_table[:, 0], pf_table[:, 1]


def interpolate_partition_function(temperatures: np.ndarray, partition_values: np.ndarray, temperature_k: float) -> float:
    if temperature_k < temperatures[0] or temperature_k > temperatures[-1]:
        raise ValueError(
            f"Temperature {temperature_k:g} K is outside the partition-function range "
            f"{temperatures[0]:g}-{temperatures[-1]:g} K"
        )
    return float(np.interp(temperature_k, temperatures, partition_values))


def load_state_arrays(states_path: Path, nstates: int) -> tuple[np.ndarray, np.ndarray]:
    energies = np.empty(nstates + 1, dtype=np.float64)
    g_totals = np.empty(nstates + 1, dtype=np.float64)
    energies.fill(np.nan)
    g_totals.fill(np.nan)

    for raw_line in iter_bz2_text_lines(states_path):
        parts = raw_line.split(None, 4)
        if len(parts) < 3:
            continue
        state_id = int(parts[0])
        energies[state_id] = float(parts[1])
        g_totals[state_id] = float(parts[2])

    if np.isnan(energies[1:]).any() or np.isnan(g_totals[1:]).any():
        raise RuntimeError(f"State table {states_path} did not fill all expected state IDs")
    return energies, g_totals


def load_scan_state_arrays(states_path: Path, nstates: int):
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

    for raw_line in iter_bz2_text_lines(states_path):
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

    if np.isnan(energies[1:]).any():
        raise RuntimeError(f"State energies were not fully populated from {states_path}")
    return energies, n1, n2, n3, n4, symmetry


def lte_line_intensity_cm_per_molecule(
    wavenumber: float,
    a_coefficient: float,
    g_upper: float,
    lower_energy_cm: float,
    temperature_k: float,
    partition_function: float,
) -> float:
    if wavenumber <= 0.0:
        return 0.0
    boltzmann = math.exp(-SECOND_RADIATION_CONSTANT_CM_K * lower_energy_cm / temperature_k)
    stimulated = 1.0 - math.exp(-SECOND_RADIATION_CONSTANT_CM_K * wavenumber / temperature_k)
    return (
        g_upper
        * a_coefficient
        * boltzmann
        * stimulated
        / (8.0 * math.pi * LIGHT_SPEED_CM_S * wavenumber * wavenumber * partition_function)
    )


def collect_relevant_transitions(
    transition_files: list[Path],
    *,
    energies: np.ndarray,
    g_totals: np.ndarray,
    wn_min: float,
    wn_max: float,
    wing: float,
    temperature_k: float,
    partition_function: float,
    intensity_threshold: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    kept_wavenumbers: list[float] = []
    kept_intensities: list[float] = []
    stats = {"files": len(transition_files), "parsed": 0, "in_window": 0, "kept": 0}

    search_min = wn_min - wing
    search_max = wn_max + wing
    for path in transition_files:
        for raw_line in iter_bz2_text_lines(path):
            parts = raw_line.split()
            if len(parts) != 3:
                continue
            upper_id = int(parts[0])
            lower_id = int(parts[1])
            a_value = float(parts[2])
            stats["parsed"] += 1

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

            kept_wavenumbers.append(wavenumber)
            kept_intensities.append(intensity)
            stats["kept"] += 1

    return np.asarray(kept_wavenumbers), np.asarray(kept_intensities), stats


def render_cross_section(
    grid: np.ndarray,
    *,
    line_centers: np.ndarray,
    line_intensities: np.ndarray,
    mass_da: float,
    temperature_k: float,
    pressure_torr: float,
    gamma0: float,
    n_exponent: float,
    line_cutoff: float,
) -> np.ndarray:
    pressure_bar = pressure_torr / 760.0 * BAR_PER_ATM
    lorentz_hwhm = gamma0 * pressure_bar * (T_REF_K / temperature_k) ** n_exponent
    doppler_hwhm = doppler_hwhm_cm(line_centers, temperature_k=temperature_k, mass_da=mass_da)

    cross_section = np.zeros_like(grid)
    for center, strength, alpha_d in zip(line_centers, line_intensities, doppler_hwhm):
        local_cutoff = max(line_cutoff, 25.0 * max(alpha_d, lorentz_hwhm))
        left = np.searchsorted(grid, center - local_cutoff, side="left")
        right = np.searchsorted(grid, center + local_cutoff, side="right")
        if left == right:
            continue
        profile = voigt_profile_cm(
            grid=grid[left:right],
            center_cm=center,
            doppler_hwhm_cm_value=alpha_d,
            lorentz_hwhm_cm_value=lorentz_hwhm,
        )
        cross_section[left:right] += strength * profile
    return cross_section


def render_absorbance(
    *,
    case: GasCase,
    window: SpectralWindow,
    data_dir: Path | None = None,
    intensity_threshold: float = 1.0e-23,
    line_cutoff: float = 25.0,
    gamma0: float | None = None,
    n_exponent: float | None = None,
    output_dir: Path | None = None,
    output_stem: str | None = None,
) -> SpectrumResult:
    resolved_data_dir = dataset_dir(data_dir)
    def_path = resolved_data_dir / f"{DATASET_STEM}.def"
    pf_path = resolved_data_dir / f"{DATASET_STEM}.pf"
    states_path = resolved_data_dir / f"{DATASET_STEM}.states.bz2"

    metadata = parse_def_file(def_path)
    gamma0_value = gamma0 if gamma0 is not None else metadata["gamma0"]
    n_exponent_value = n_exponent if n_exponent is not None else metadata["n_exponent"]

    pf_temperatures, pf_values = load_partition_function(pf_path)
    partition_function = interpolate_partition_function(pf_temperatures, pf_values, case.temperature_k)
    energies, g_totals = load_state_arrays(states_path, int(metadata["nstates"]))
    transition_files, effective_wing = overlapping_transition_files(
        resolved_data_dir,
        wn_min=window.wn_min,
        wn_max=window.wn_max,
        wing=line_cutoff,
    )
    line_centers, line_intensities, stats = collect_relevant_transitions(
        transition_files,
        energies=energies,
        g_totals=g_totals,
        wn_min=window.wn_min,
        wn_max=window.wn_max,
        wing=effective_wing,
        temperature_k=case.temperature_k,
        partition_function=partition_function,
        intensity_threshold=intensity_threshold,
    )
    if len(line_centers) == 0:
        raise RuntimeError("No transitions survived the current spectral/intensity filters")

    grid = build_grid(window)
    cross_section = render_cross_section(
        grid,
        line_centers=line_centers,
        line_intensities=line_intensities,
        mass_da=metadata["mass_da"],
        temperature_k=case.temperature_k,
        pressure_torr=case.pressure_torr,
        gamma0=gamma0_value,
        n_exponent=n_exponent_value,
        line_cutoff=line_cutoff,
    )
    absorbance = cross_section_to_absorbance(cross_section, case)
    result = SpectrumResult(
        wavenumber=grid,
        values=absorbance,
        quantity="absorbance",
        metadata={
            "dataset_stem": DATASET_STEM,
            "transition_files": [path.name for path in transition_files],
            "line_cutoff": line_cutoff,
            "effective_wing": effective_wing,
            "intensity_threshold": intensity_threshold,
            "transition_stats": stats,
            "partition_function": partition_function,
            "gamma0": gamma0_value,
            "n_exponent": n_exponent_value,
        },
    )

    if output_dir is not None:
        ensure_directory(output_dir)
        stem = output_stem or (
            f"CH4_MM_abs_{int(window.wn_min)}_{int(window.wn_max)}"
            f"_T{case.temperature_k:g}K_P{case.pressure_torr:g}Torr"
        )
        result.csv_path = save_spectrum_csv(output_dir / f"{stem}.csv", result.wavenumber, result.values, "absorbance")
        result.html_path = save_spectrum_html(
            output_dir / f"{stem}.html",
            result.wavenumber,
            result.values,
            title=f"CH4 ExoMol MM Absorbance {window.wn_min:g}-{window.wn_max:g} cm^-1",
            trace_name="ExoMol MM",
            y_label="Absorbance",
        )
    return result


def symmetry_label(code: int) -> str:
    if code <= 0 or code > len(EXOMOL_SYMMETRIES):
        return "?"
    return EXOMOL_SYMMETRIES[code - 1]


def exomol_band_label(n1_value: int, n2_value: int, n3_value: int, n4_value: int, sym_code: int) -> str:
    return f"{n1_value} {n2_value} {n3_value} {n4_value} {symmetry_label(sym_code)}"


def hitran_style_band_label(n1_value: int, n2_value: int, n3_value: int, n4_value: int, sym_code: int) -> str:
    raw_symmetry = symmetry_label(sym_code)
    return f"{n1_value} {n2_value} {n3_value} {n4_value} {HITRAN_STYLE_MAP.get(raw_symmetry, raw_symmetry)}"


def scan_band_types(
    *,
    window: SpectralWindow,
    data_dir: Path | None = None,
    allow_other_modes: bool = False,
    output_dir: Path | None = None,
    output_stem: str | None = None,
    print_top_exact: int = 20,
) -> SummaryResult:
    resolved_data_dir = dataset_dir(data_dir)
    require_zero_other_modes = not allow_other_modes
    def_path = resolved_data_dir / f"{DATASET_STEM}.def"
    states_path = resolved_data_dir / f"{DATASET_STEM}.states.bz2"
    metadata = parse_def_file(def_path)
    energies, n1, n2, n3, n4, symmetry = load_scan_state_arrays(states_path, int(metadata["nstates"]))
    transition_files, _ = overlapping_transition_files(resolved_data_dir, wn_min=window.wn_min, wn_max=window.wn_max)

    category_counts: Counter[tuple[int, int]] = Counter()
    exact_counts: Counter[tuple[str, str]] = Counter()
    exact_to_hitran: dict[tuple[str, str], tuple[str, str]] = {}
    total_in_window = 0
    total_kept = 0

    for path in transition_files:
        for raw_line in iter_bz2_text_lines(path):
            parts = raw_line.split()
            if len(parts) != 3:
                continue
            upper_id = int(parts[0])
            lower_id = int(parts[1])
            wavenumber = energies[upper_id] - energies[lower_id]
            if wavenumber < window.wn_min or wavenumber > window.wn_max:
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
            total_kept += 1

    rows: list[dict[str, object]] = []
    for (lower_label, upper_label), count in exact_counts.most_common():
        lower_hitran, upper_hitran = exact_to_hitran[(lower_label, upper_label)]
        lower_q = int(lower_label.split()[2])
        upper_q = int(upper_label.split()[2])
        rows.append(
            {
                "category": f"nu3 {lower_q}->{upper_q}",
                "nu3_lower": lower_q,
                "nu3_upper": upper_q,
                "lower_label_exomol": lower_label,
                "upper_label_exomol": upper_label,
                "lower_label_hitran_style": lower_hitran,
                "upper_label_hitran_style": upper_hitran,
                "transition_count": count,
            }
        )

    result = SummaryResult(
        rows=rows,
        metadata={
            "window_transitions_found": total_in_window,
            "nu3_transitions_kept": total_kept,
            "distinct_exact_pairs": len(exact_counts),
            "category_counts": {f"{low}->{high}": count for (low, high), count in sorted(category_counts.items())},
            "top_exact_count": print_top_exact,
        },
    )

    if output_dir is not None:
        ensure_directory(output_dir)
        stem = output_stem or f"CH4_MM_nu3_band_types_{int(window.wn_min)}_{int(window.wn_max)}"
        result.csv_path = write_rows_csv(output_dir / f"{stem}_exact_pairs.csv", rows)

    return result


def sorted_nu3_band_dir(root_dir: Path | None = None) -> Path:
    paths = default_paths(root_dir=root_dir)
    return paths.root_dir / SORTED_NU3_FOLDER_NAME


def _parse_sorted_nu3_filename(path: Path) -> dict[str, object]:
    match = SORTED_NU3_FILENAME_PATTERN.match(path.name)
    if match is None:
        raise ValueError(f"Unrecognized sorted nu3 export filename: {path.name}")
    groups = match.groupdict()
    lower_mode = tuple(int(groups[f"lower_n{index}"]) for index in range(1, 5))
    upper_mode = tuple(int(groups[f"upper_n{index}"]) for index in range(1, 5))
    return {
        "path": path,
        "lower_mode": lower_mode,
        "upper_mode": upper_mode,
        "lower_sym": str(groups["lower_sym"]),
        "upper_sym": str(groups["upper_sym"]),
        "temperature_k": float(groups["temperature"]),
    }


def _progression_label(lower_mode: tuple[int, int, int, int], upper_mode: tuple[int, int, int, int]) -> str:
    return f"nu3 {lower_mode[2]}->{upper_mode[2]}"


def _progression_slug(lower_mode: tuple[int, int, int, int], upper_mode: tuple[int, int, int, int]) -> str:
    return f"nu3_{lower_mode[2]}_to_{upper_mode[2]}"


def _mode_label(mode: tuple[int, int, int, int]) -> str:
    return " ".join(str(value) for value in mode)


def _extract_j_value(local_label: str) -> int:
    match = J_VALUE_PATTERN.search(local_label.strip())
    if match is None:
        raise ValueError(f"Could not extract J from local label: {local_label!r}")
    return int(match.group(1))


def _iter_sorted_band_rows(path: Path):
    header: list[str] | None = None
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if header is None:
                header = line.split("\t")
                continue
            values = line.split("\t")
            if len(values) != len(header):
                continue
            yield dict(zip(header, values))


def _color_for_index(index: int, total: int) -> str:
    if total <= 1:
        return "#1f77b4"
    cmap = plt.get_cmap("turbo")
    return mcolors.to_hex(cmap(index / max(1, total - 1)), keep_alpha=False)


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
    y_key: str,
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
            y_values = np.asarray(trace[y_key], dtype=float)
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


def _label_candidates(
    traces: list[dict[str, object]],
    label_top_n: int,
    *,
    y_key: str = "intensity",
    peak_key: str = "peak_intensity",
    total_key: str = "total_intensity",
    required_labels: tuple[str, ...] = (),
) -> list[dict[str, object]]:
    if label_top_n <= 0 and not required_labels:
        return []
    ranked = sorted(
        traces,
        key=lambda trace: (
            float(trace[peak_key]),
            float(trace[total_key]),
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

    if not selected:
        return []

    candidates: list[dict[str, object]] = []
    for trace in selected:
        x_values = np.asarray(trace["wavenumber"], dtype=float)
        y_values = np.asarray(trace[y_key], dtype=float)
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


def _decimate_series(x_values: np.ndarray, y_values: np.ndarray, max_points: int) -> tuple[np.ndarray, np.ndarray]:
    if max_points <= 0 or len(x_values) <= max_points:
        return x_values, y_values
    indices = np.linspace(0, len(x_values) - 1, max_points, dtype=np.int64)
    return x_values[indices], y_values[indices]


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
        positioned.append(
            {
                **item,
                "label_x": anchor_x,
                "label_y": float(target_y),
            }
        )
    return positioned


def _save_progression_png(
    path: Path,
    *,
    progression_label: str,
    lower_mode: tuple[int, int, int, int],
    upper_mode: tuple[int, int, int, int],
    traces: list[dict[str, object]],
    labeled_traces: list[dict[str, object]],
    wn_min: float,
    wn_max: float,
) -> Path:
    ensure_directory(path.parent)
    figure, axis = plt.subplots(figsize=(16, 9), constrained_layout=True)

    for trace in traces:
        axis.plot(
            trace["wavenumber"],
            trace["intensity"],
            color=trace["color"],
            linewidth=1.0,
            alpha=0.95,
        )

    for item in labeled_traces:
        axis.plot(
            [item["peak_x"], item["label_x"]],
            [item["peak_y"], item["label_y"]],
            color=item["color"],
            linewidth=0.8,
            alpha=0.8,
        )
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
    axis.set_xlabel("Wavenumber (cm^-1)")
    axis.set_ylabel("Intensity (cm/molecule)")
    axis.set_title(
        f"{progression_label} intensity by J pair\n"
        f"{_mode_label(lower_mode)} -> {_mode_label(upper_mode)}"
    )
    axis.grid(True, alpha=0.25, linewidth=0.5)
    axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

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
    labeled_traces: list[dict[str, object]],
    wn_min: float,
    wn_max: float,
) -> Path:
    ensure_directory(path.parent)
    labeled_keys = {item["text"] for item in labeled_traces}
    figure = go.Figure()

    for trace in traces:
        figure.add_trace(
            go.Scattergl(
                x=trace["wavenumber"],
                y=trace["intensity"],
                mode="lines",
                name=str(trace["jpair_label"]),
                meta=str(trace["jpair_label"]),
                legendgroup=str(trace["jpair_label"]),
                showlegend=str(trace["jpair_label"]) in labeled_keys,
                line={"color": trace["color"], "width": 1.3},
                hovertemplate=(
                    "J pair: %{meta}<br>"
                    "Wavenumber: %{x:.6f} cm^-1<br>"
                    "Intensity: %{y:.3e} cm/molecule<br>"
                    "Points in trace: "
                    + str(trace["point_count"])
                    + "<extra></extra>"
                ),
            )
        )

    for item in labeled_traces:
        figure.add_annotation(
            x=item["label_x"],
            y=item["label_y"],
            text=item["text"],
            showarrow=True,
            arrowhead=0,
            arrowcolor=item["color"],
            arrowwidth=1.0,
            axref="x",
            ayref="y",
            ax=item["peak_x"],
            ay=item["peak_y"],
            font={"color": item["color"], "size": 11},
            bgcolor="rgba(255,255,255,0.7)",
        )

    figure.update_layout(
        title=(
            f"{progression_label} intensity by J pair"
            f"<br><sup>{_mode_label(lower_mode)} -> {_mode_label(upper_mode)}</sup>"
        ),
        template="plotly_white",
        xaxis_title="Wavenumber (cm^-1)",
        yaxis_title="Intensity (cm/molecule)",
        hovermode="x unified",
        width=1400,
        height=800,
        legend_title_text="Labeled J pairs",
    )
    figure.update_xaxes(range=[wn_min, wn_max])
    figure.write_html(path, include_plotlyjs="cdn")
    return path


def plot_sorted_nu3_intensity_progressions(
    *,
    input_dir: Path | None = None,
    output_dir: Path | None = None,
    wn_min: float = 2500.0,
    wn_max: float = 3500.0,
    label_top_n: int = 8,
) -> SummaryResult:
    if wn_max <= wn_min:
        raise ValueError("wn_max must be greater than wn_min")
    if label_top_n < 0:
        raise ValueError("label_top_n must be non-negative")

    resolved_input_dir = (input_dir or sorted_nu3_band_dir()).resolve()
    if not resolved_input_dir.exists():
        raise FileNotFoundError(f"Sorted ExoMol folder not found: {resolved_input_dir}")

    resolved_output_dir = ensure_directory(
        (output_dir or (default_paths().artifacts_dir / "exomol_sorted_nu3_intensity")).resolve()
    )

    grouped_files: dict[tuple[tuple[int, int, int, int], tuple[int, int, int, int]], list[dict[str, object]]] = defaultdict(list)
    for path in sorted(resolved_input_dir.glob("*.txt")):
        parsed = _parse_sorted_nu3_filename(path)
        key = (parsed["lower_mode"], parsed["upper_mode"])
        grouped_files[key].append(parsed)

    if not grouped_files:
        raise FileNotFoundError(f"No sorted ExoMol band text files found in {resolved_input_dir}")

    summary_rows: list[dict[str, object]] = []

    for (lower_mode, upper_mode), file_infos in sorted(grouped_files.items(), key=lambda item: (item[0][0][2], item[0][1][2])):
        progression_label = _progression_label(lower_mode, upper_mode)
        progression_slug = _progression_slug(lower_mode, upper_mode)
        grouped_rows: dict[tuple[int, int], dict[str, object]] = defaultdict(
            lambda: {"wavenumber": [], "intensity": [], "source_files": set()}
        )
        row_count = 0

        for file_info in file_infos:
            path = Path(file_info["path"])
            for row in _iter_sorted_band_rows(path):
                wavenumber = float(row["wavenumber_cm-1"])
                if wavenumber < wn_min or wavenumber > wn_max:
                    continue
                lower_j = _extract_j_value(row["local_lower_quanta"])
                upper_j = _extract_j_value(row["local_upper_quanta"])
                group = grouped_rows[(lower_j, upper_j)]
                group["wavenumber"].append(wavenumber)
                group["intensity"].append(float(row["line_intensity_cm_per_molecule"]))
                group["source_files"].add(path.name)
                row_count += 1

        traces: list[dict[str, object]] = []
        for index, ((lower_j, upper_j), payload) in enumerate(sorted(grouped_rows.items())):
            x_values = np.asarray(payload["wavenumber"], dtype=float)
            y_values = np.asarray(payload["intensity"], dtype=float)
            order = np.argsort(x_values)
            x_values = x_values[order]
            y_values = y_values[order]
            trace = {
                "trace_index": index + 1,
                "lower_j": lower_j,
                "upper_j": upper_j,
                "jpair_label": _format_jpair_label(lower_j, upper_j),
                "wavenumber": x_values,
                "intensity": y_values,
                "point_count": int(len(x_values)),
                "peak_intensity": float(np.max(y_values)) if len(y_values) else 0.0,
                "total_intensity": float(np.sum(y_values)) if len(y_values) else 0.0,
                "source_file_count": len(payload["source_files"]),
            }
            traces.append(trace)

        if not traces:
            continue

        for index, trace in enumerate(traces):
            trace["color"] = _color_for_index(index, len(traces))

        labeled_candidates = _label_candidates(traces, label_top_n)
        labeled_traces = _label_positions(labeled_candidates, wn_min, wn_max)
        labeled_names = {item["text"] for item in labeled_traces}

        mapping_rows: list[dict[str, object]] = []
        ranked_for_mapping = sorted(
            traces,
            key=lambda trace: (
                float(trace["peak_intensity"]),
                float(trace["total_intensity"]),
                -int(trace["lower_j"]),
                -int(trace["upper_j"]),
            ),
            reverse=True,
        )
        rank_lookup = {
            (int(trace["lower_j"]), int(trace["upper_j"])): rank
            for rank, trace in enumerate(ranked_for_mapping, start=1)
        }
        for trace in traces:
            mapping_rows.append(
                {
                    "trace_index": trace["trace_index"],
                    "strength_rank": rank_lookup[(int(trace["lower_j"]), int(trace["upper_j"]))],
                    "lower_j": trace["lower_j"],
                    "upper_j": trace["upper_j"],
                    "jpair_label": trace["jpair_label"],
                    "point_count": trace["point_count"],
                    "peak_intensity": f"{float(trace['peak_intensity']):.12e}",
                    "total_intensity": f"{float(trace['total_intensity']):.12e}",
                    "color_hex": trace["color"],
                    "source_file_count": trace["source_file_count"],
                    "labeled_on_figure": "yes" if trace["jpair_label"] in labeled_names else "no",
                }
            )

        png_path = _save_progression_png(
            resolved_output_dir / f"{progression_slug}_intensity.png",
            progression_label=progression_label,
            lower_mode=lower_mode,
            upper_mode=upper_mode,
            traces=traces,
            labeled_traces=labeled_traces,
            wn_min=wn_min,
            wn_max=wn_max,
        )
        html_path = _save_progression_html(
            resolved_output_dir / f"{progression_slug}_intensity.html",
            progression_label=progression_label,
            lower_mode=lower_mode,
            upper_mode=upper_mode,
            traces=traces,
            labeled_traces=labeled_traces,
            wn_min=wn_min,
            wn_max=wn_max,
        )
        mapping_path = write_rows_csv(
            resolved_output_dir / f"{progression_slug}_jpairs.csv",
            mapping_rows,
            fieldnames=[
                "trace_index",
                "strength_rank",
                "lower_j",
                "upper_j",
                "jpair_label",
                "point_count",
                "peak_intensity",
                "total_intensity",
                "color_hex",
                "source_file_count",
                "labeled_on_figure",
            ],
        )

        summary_rows.append(
            {
                "progression_label": progression_label,
                "lower_mode": _mode_label(lower_mode),
                "upper_mode": _mode_label(upper_mode),
                "source_file_count": len(file_infos),
                "row_count": row_count,
                "jpair_count": len(traces),
                "labeled_j_pairs": ", ".join(item["text"] for item in labeled_traces) if labeled_traces else "",
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
            "labeled_j_pairs",
            "png_file",
            "html_file",
            "mapping_csv",
        ],
    )

    report_lines = [
        "# ExoMol Sorted nu3 Intensity Progressions",
        "",
        f"- Input folder: `{resolved_input_dir}`",
        f"- Wavenumber window: `{wn_min:g}` to `{wn_max:g} cm^-1`",
        "- Y axis: `line_intensity_cm_per_molecule`",
        "- Curve grouping: full J pair `(lower J, upper J)`",
        f"- On-figure labels: strongest `{label_top_n}` J pairs per progression",
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
                f"- Rows plotted: `{row['row_count']}`",
                f"- J-pair curves: `{row['jpair_count']}`",
                f"- Labeled J pairs: `{row['labeled_j_pairs'] or 'none'}`",
                f"- Outputs: [PNG]({row['png_file']}), [HTML]({row['html_file']}), [J-pair CSV]({row['mapping_csv']})",
                "",
                f"![{row['progression_label']}]({row['png_file']})",
                "",
            ]
        )

    report_path = write_markdown(resolved_output_dir / "report.md", "\n".join(report_lines))
    result = SummaryResult(
        rows=summary_rows,
        csv_path=summary_csv_path,
        metadata={
            "input_dir": resolved_input_dir,
            "output_dir": resolved_output_dir,
            "report_path": report_path,
            "progression_count": len(summary_rows),
            "label_top_n": label_top_n,
            "wn_min": wn_min,
            "wn_max": wn_max,
        },
    )
    result.metadata["report_path"] = report_path
    return result


def _collect_sorted_progression_groups(
    input_dir: Path,
    *,
    wn_min: float,
    wn_max: float,
    min_line_intensity: float = 0.0,
) -> list[dict[str, object]]:
    grouped_files: dict[tuple[tuple[int, int, int, int], tuple[int, int, int, int]], list[dict[str, object]]] = defaultdict(list)
    for path in sorted(input_dir.glob("*.txt")):
        parsed = _parse_sorted_nu3_filename(path)
        key = (parsed["lower_mode"], parsed["upper_mode"])
        grouped_files[key].append(parsed)

    progression_groups: list[dict[str, object]] = []
    for (lower_mode, upper_mode), file_infos in sorted(grouped_files.items(), key=lambda item: (item[0][0][2], item[0][1][2])):
        grouped_rows: dict[tuple[int, int], dict[str, object]] = defaultdict(
            lambda: {"wavenumber": [], "intensity": [], "source_files": set()}
        )
        row_count = 0

        for file_info in file_infos:
            path = Path(file_info["path"])
            for row in _iter_sorted_band_rows(path):
                wavenumber = float(row["wavenumber_cm-1"])
                if wavenumber < wn_min or wavenumber > wn_max:
                    continue
                intensity = float(row["line_intensity_cm_per_molecule"])
                if intensity < min_line_intensity:
                    continue
                lower_j = _extract_j_value(row["local_lower_quanta"])
                upper_j = _extract_j_value(row["local_upper_quanta"])
                group = grouped_rows[(lower_j, upper_j)]
                group["wavenumber"].append(wavenumber)
                group["intensity"].append(intensity)
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


def _save_validation_png(
    path: Path,
    *,
    grid: np.ndarray,
    module_absorbance: np.ndarray,
    legacy_absorbance: np.ndarray,
    delta_absorbance: np.ndarray,
) -> Path:
    ensure_directory(path.parent)
    figure, axes = plt.subplots(2, 1, figsize=(15, 10), constrained_layout=True, sharex=True)

    axes[0].plot(grid, legacy_absorbance, color="#222222", linewidth=1.2, label="Legacy render")
    axes[0].plot(grid, module_absorbance, color="#1f77b4", linewidth=1.0, alpha=0.85, label="Module render")
    axes[0].set_ylabel("Absorbance")
    axes[0].set_title("Absorbance module validation: nu3 0->1")
    axes[0].grid(True, alpha=0.25, linewidth=0.5)
    axes[0].legend(loc="upper right")

    axes[1].plot(grid, delta_absorbance, color="#d62728", linewidth=1.0)
    axes[1].axhline(0.0, color="#666666", linewidth=0.8, linestyle="--")
    axes[1].set_xlabel("Wavenumber (cm^-1)")
    axes[1].set_ylabel("Module - Legacy")
    axes[1].grid(True, alpha=0.25, linewidth=0.5)

    figure.savefig(path, dpi=200)
    plt.close(figure)
    return path


def _save_validation_html(
    path: Path,
    *,
    grid: np.ndarray,
    module_absorbance: np.ndarray,
    legacy_absorbance: np.ndarray,
    delta_absorbance: np.ndarray,
    max_points: int,
) -> Path:
    ensure_directory(path.parent)
    grid_plot, module_plot = _decimate_series(grid, module_absorbance, max_points)
    _, legacy_plot = _decimate_series(grid, legacy_absorbance, max_points)
    _, delta_plot = _decimate_series(grid, delta_absorbance, max_points)

    figure = go.Figure()
    figure.add_trace(go.Scattergl(x=grid_plot, y=legacy_plot, mode="lines", name="Legacy render", line={"color": "#222222"}))
    figure.add_trace(go.Scattergl(x=grid_plot, y=module_plot, mode="lines", name="Module render", line={"color": "#1f77b4"}))
    figure.add_trace(go.Scattergl(x=grid_plot, y=delta_plot, mode="lines", name="Module - Legacy", line={"color": "#d62728"}, yaxis="y2"))
    figure.update_layout(
        title="Absorbance module validation: nu3 0->1",
        template="plotly_white",
        xaxis_title="Wavenumber (cm^-1)",
        yaxis_title="Absorbance",
        yaxis2={"title": "Module - Legacy", "overlaying": "y", "side": "right", "showgrid": False},
        width=1400,
        height=800,
    )
    figure.write_html(path, include_plotlyjs="cdn")
    return path


def validate_sorted_nu3_absorbance_module(
    *,
    input_dir: Path | None = None,
    output_dir: Path | None = None,
    wn_min: float = 2500.0,
    wn_max: float = 3500.0,
    wn_step: float = 0.25,
    pressure_torr: float = 3.0,
    mole_fraction: float = 0.008,
    path_length_cm: float = 100.0,
    line_cutoff: float = 0.5,
    min_line_intensity: float = 0.0,
    html_max_points: int = 5000,
) -> SummaryResult:
    resolved_input_dir = (input_dir or sorted_nu3_band_dir()).resolve()
    resolved_output_dir = ensure_directory(
        (output_dir or (default_paths().artifacts_dir / "exomol_sorted_nu3_absorbance" / "validation")).resolve()
    )
    case = GasCase(
        temperature_k=296.0,
        pressure_torr=pressure_torr,
        mole_fraction=mole_fraction,
        path_length_cm=path_length_cm,
    )
    window = SpectralWindow(wn_min=wn_min, wn_max=wn_max, wn_step=wn_step)

    progression_groups = _collect_sorted_progression_groups(
        resolved_input_dir,
        wn_min=wn_min,
        wn_max=wn_max,
        min_line_intensity=min_line_intensity,
    )
    target = next(
        (
            item
            for item in progression_groups
            if item["lower_mode"] == (0, 0, 0, 0) and item["upper_mode"] == (0, 0, 1, 0)
        ),
        None,
    )
    if target is None:
        raise RuntimeError("Could not find the sorted nu3 0->1 progression for validation")

    all_centers: list[float] = []
    all_intensities: list[float] = []
    for payload in target["grouped_rows"].values():
        all_centers.extend(payload["wavenumber"])
        all_intensities.extend(payload["intensity"])

    line_centers = np.asarray(all_centers, dtype=np.float64)
    line_intensities = np.asarray(all_intensities, dtype=np.float64)
    if line_centers.size == 0:
        raise RuntimeError("Validation progression 0->1 has no lines after filtering")

    metadata = parse_def_file(dataset_dir() / f"{DATASET_STEM}.def")
    grid = build_grid(window)
    module_absorbance = render_absorbance_on_grid(
        grid,
        case=case,
        line_centers=line_centers,
        line_intensities=line_intensities,
        mass_da=float(metadata["mass_da"]),
        gamma0=float(metadata["gamma0"]),
        n_exponent=float(metadata["n_exponent"]),
        line_cutoff=line_cutoff,
    )
    legacy_cross_section = render_cross_section(
        grid,
        line_centers=line_centers,
        line_intensities=line_intensities,
        mass_da=float(metadata["mass_da"]),
        temperature_k=case.temperature_k,
        pressure_torr=case.pressure_torr,
        gamma0=float(metadata["gamma0"]),
        n_exponent=float(metadata["n_exponent"]),
        line_cutoff=line_cutoff,
    )
    legacy_absorbance = cross_section_to_absorbance(legacy_cross_section, case)
    delta_absorbance = module_absorbance - legacy_absorbance

    legacy_peak = float(np.max(np.abs(legacy_absorbance)))
    signal_floor = max(1.0e-12, legacy_peak * 1.0e-8)
    signal_mask = np.abs(legacy_absorbance) > signal_floor
    relative_delta = np.zeros_like(delta_absorbance)
    if np.any(signal_mask):
        relative_delta[signal_mask] = np.abs(delta_absorbance[signal_mask]) / np.abs(legacy_absorbance[signal_mask])

    module_peak_index = int(np.argmax(module_absorbance))
    legacy_peak_index = int(np.argmax(legacy_absorbance))
    metrics_row = {
        "progression_label": "nu3 0->1",
        "line_count": int(line_centers.size),
        "grid_point_count": int(grid.size),
        "max_abs_difference": f"{float(np.max(np.abs(delta_absorbance))):.12e}",
        "mean_abs_difference": f"{float(np.mean(np.abs(delta_absorbance))):.12e}",
        "max_relative_difference": f"{float(np.max(relative_delta)):.12e}",
        "mean_relative_difference": f"{float(np.mean(relative_delta[signal_mask])):.12e}" if np.any(signal_mask) else "0.000000000000e+00",
        "legacy_peak_absorbance": f"{float(np.max(legacy_absorbance)):.12e}",
        "module_peak_absorbance": f"{float(np.max(module_absorbance)):.12e}",
        "legacy_peak_wavenumber_cm-1": f"{float(grid[legacy_peak_index]):.6f}",
        "module_peak_wavenumber_cm-1": f"{float(grid[module_peak_index]):.6f}",
        "peak_wavenumber_delta_cm-1": f"{abs(float(grid[module_peak_index]) - float(grid[legacy_peak_index])):.6e}",
    }

    png_path = _save_validation_png(
        resolved_output_dir / "validation_nu3_0_to_1.png",
        grid=grid,
        module_absorbance=module_absorbance,
        legacy_absorbance=legacy_absorbance,
        delta_absorbance=delta_absorbance,
    )
    html_path = _save_validation_html(
        resolved_output_dir / "validation_nu3_0_to_1.html",
        grid=grid,
        module_absorbance=module_absorbance,
        legacy_absorbance=legacy_absorbance,
        delta_absorbance=delta_absorbance,
        max_points=html_max_points,
    )
    csv_path = write_rows_csv(
        resolved_output_dir / "validation_metrics.csv",
        [metrics_row],
        fieldnames=list(metrics_row.keys()),
    )

    report_lines = [
        "# Absorbance Module Validation",
        "",
        "- Validation target: `nu3 0->1` from the sorted ExoMol folder",
        "- Comparison: new `research.absorbance` module vs existing `research.exomol.render_cross_section` path",
        "- Note: the sorted line intensities are exported at `296 K`, so the validation case keeps `T = 296 K`",
        f"- Window: `{wn_min:g}` to `{wn_max:g} cm^-1` with `step = {wn_step:g} cm^-1`",
        f"- Broadening cutoff: `{line_cutoff:g} cm^-1`",
        f"- Minimum line intensity kept: `{min_line_intensity:.3e} cm/molecule`",
        "",
        f"- Metrics CSV: [{csv_path.name}]({csv_path.name})",
        f"- Validation HTML: [{html_path.name}]({html_path.name})",
        "",
        "## Metrics",
        "",
    ]
    for key, value in metrics_row.items():
        report_lines.append(f"- {key}: `{value}`")
    report_lines.extend(["", f"![Validation overlay]({png_path.name})", ""])
    report_path = write_markdown(resolved_output_dir / "validation_nu3_0_to_1.md", "\n".join(report_lines))

    result = SummaryResult(
        rows=[metrics_row],
        csv_path=csv_path,
        metadata={
            "report_path": report_path,
            "png_path": png_path,
            "html_path": html_path,
            "window": window,
            "case": case,
            "line_cutoff": line_cutoff,
            "min_line_intensity": min_line_intensity,
        },
    )
    return result


def _save_absorbance_progression_png(
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
        f"{progression_label} absorbance by J pair, split by delta J\n"
        f"{_mode_label(lower_mode)} -> {_mode_label(upper_mode)}"
    )

    for axis, delta_j in zip(axes, ABSORBANCE_DELTA_J_VALUES):
        branch_traces = [trace for trace in traces if int(trace["delta_j"]) == delta_j and bool(trace["plotted_in_figure"])]
        for trace in branch_traces:
            axis.plot(
                trace["wavenumber"],
                trace["absorbance"],
                color=trace["color"],
                linewidth=1.0,
                alpha=0.95,
            )

        labeled_traces = labeled_traces_by_delta_j.get(delta_j, [])
        for item in labeled_traces:
            axis.plot(
                [item["peak_x"], item["label_x"]],
                [item["peak_y"], item["label_y"]],
                color=item["color"],
                linewidth=0.8,
                alpha=0.8,
            )
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
        y_min, y_max = compute_panel_y_limits(branch_traces, y_key="absorbance", label_items=labeled_traces)
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


def _save_absorbance_progression_html(
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
        title=(
            f"{progression_label} absorbance by J pair, split by delta J"
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


def plot_sorted_nu3_absorbance_progressions(
    *,
    input_dir: Path | None = None,
    output_dir: Path | None = None,
    wn_min: float = 2500.0,
    wn_max: float = 3500.0,
    wn_step: float = 0.25,
    pressure_torr: float = 3.0,
    mole_fraction: float = 0.008,
    path_length_cm: float = 100.0,
    line_cutoff: float = 0.5,
    min_line_intensity: float = 0.0,
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

    resolved_input_dir = (input_dir or sorted_nu3_band_dir()).resolve()
    if not resolved_input_dir.exists():
        raise FileNotFoundError(f"Sorted ExoMol folder not found: {resolved_input_dir}")

    resolved_output_dir = ensure_directory(
        (output_dir or (default_paths().artifacts_dir / "exomol_sorted_nu3_absorbance")).resolve()
    )
    validation_output_dir = ensure_directory(resolved_output_dir / "validation")
    effective_forced_j_pairs = _merge_forced_j_pairs(DEFAULT_FORCED_ABSORBANCE_J_PAIRS, forced_j_pairs)

    validation_result = validate_sorted_nu3_absorbance_module(
        input_dir=resolved_input_dir,
        output_dir=validation_output_dir,
        wn_min=wn_min,
        wn_max=wn_max,
        wn_step=wn_step,
        pressure_torr=pressure_torr,
        mole_fraction=mole_fraction,
        path_length_cm=path_length_cm,
        line_cutoff=line_cutoff,
        min_line_intensity=min_line_intensity,
        html_max_points=html_max_points,
    )

    progression_groups = _collect_sorted_progression_groups(
        resolved_input_dir,
        wn_min=wn_min,
        wn_max=wn_max,
        min_line_intensity=min_line_intensity,
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

    summary_rows: list[dict[str, object]] = []
    for progression in progression_groups:
        lower_mode = progression["lower_mode"]
        upper_mode = progression["upper_mode"]
        progression_label = str(progression["progression_label"])
        progression_slug = str(progression["progression_slug"])
        grouped_rows = progression["grouped_rows"]

        traces: list[dict[str, object]] = []
        for index, ((lower_j, upper_j), payload) in enumerate(sorted(grouped_rows.items())):
            line_centers = np.asarray(payload["wavenumber"], dtype=np.float64)
            line_intensities = np.asarray(payload["intensity"], dtype=np.float64)
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
            traces.append(
                {
                    "trace_index": index + 1,
                    "lower_j": lower_j,
                    "upper_j": upper_j,
                    "delta_j": _delta_j_value(lower_j, upper_j),
                    "branch_label": _delta_j_branch_label(_delta_j_value(lower_j, upper_j)),
                    "jpair_label": _format_jpair_label(lower_j, upper_j),
                    "wavenumber": grid,
                    "absorbance": absorbance,
                    "grid_point_count": int(grid.size),
                    "line_count": int(line_centers.size),
                    "peak_absorbance": float(np.max(absorbance)) if len(absorbance) else 0.0,
                    "integrated_absorbance": float(np.trapezoid(absorbance, grid)) if len(absorbance) else 0.0,
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
        png_path = _save_absorbance_progression_png(
            resolved_output_dir / f"{progression_slug}_absorbance.png",
            progression_label=progression_label,
            lower_mode=lower_mode,
            upper_mode=upper_mode,
            traces=traces,
            labeled_traces_by_delta_j=labeled_traces_by_delta_j,
            wn_min=wn_min,
            wn_max=wn_max,
        )
        html_path = _save_absorbance_progression_html(
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
                "source_file_count": len(progression["file_infos"]),
                "row_count": progression["row_count"],
                "jpair_count": len(traces),
                "grid_point_count": int(grid.size),
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

    validation_report_path = Path(validation_result.metadata["report_path"])
    validation_png_path = Path(validation_result.metadata["png_path"])
    validation_html_path = Path(validation_result.metadata["html_path"])
    validation_metrics = validation_result.rows[0]

    report_lines = [
        "# ExoMol Sorted nu3 Absorbance Progressions",
        "",
        f"- Input folder: `{resolved_input_dir}`",
        f"- Wavenumber window: `{wn_min:g}` to `{wn_max:g} cm^-1` with `step = {wn_step:g} cm^-1`",
        "- Y axis: `absorbance`",
        "- Curve grouping: full J pair `(lower J, upper J)`, split into `delta J = -1, 0, +1` panels",
        "- Line intensities come from the sorted `T296.0K` ExoMol exports, so this workflow keeps `T = 296 K`",
        f"- Pressure: `{pressure_torr:g} Torr`",
        f"- Mole fraction: `{mole_fraction:g}`",
        f"- Path length: `{path_length_cm:g} cm`",
        f"- Broadening cutoff: `{line_cutoff:g} cm^-1`",
        f"- Minimum line intensity kept: `{min_line_intensity:.3e} cm/molecule`",
        f"- On-figure labels: strongest `{label_top_n_per_delta_j}` J pairs per `delta J` panel, plus forced labels for `{', '.join(_format_jpair_label(lower_j, upper_j) for lower_j, upper_j in effective_forced_j_pairs)}` when present",
        f"- HTML traces are decimated to at most `{html_max_points}` points per J pair for responsiveness",
        f"- Summary CSV: [{summary_csv_path.name}]({summary_csv_path.name})",
        "",
        "## Validation",
        "",
        f"- Validation report: [validation/{validation_report_path.name}](validation/{validation_report_path.name})",
        f"- Validation HTML: [validation/{validation_html_path.name}](validation/{validation_html_path.name})",
        f"- Max absolute difference: `{validation_metrics['max_abs_difference']}`",
        f"- Mean absolute difference: `{validation_metrics['mean_abs_difference']}`",
        f"- Max relative difference: `{validation_metrics['max_relative_difference']}`",
        f"- Peak wavenumber delta: `{validation_metrics['peak_wavenumber_delta_cm-1']} cm^-1`",
        "",
        f"![Validation overlay](validation/{validation_png_path.name})",
        "",
    ]

    for row in summary_rows:
        report_lines.extend(
            [
                f"## {row['progression_label']}",
                "",
                f"- Modes: `{row['lower_mode']} -> {row['upper_mode']}`",
                f"- Files merged: `{row['source_file_count']}`",
                f"- Lines kept after filtering: `{row['row_count']}`",
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
    result = SummaryResult(
        rows=summary_rows,
        csv_path=summary_csv_path,
        metadata={
            "input_dir": resolved_input_dir,
            "output_dir": resolved_output_dir,
            "report_path": report_path,
            "validation_report_path": validation_report_path,
            "progression_count": len(summary_rows),
            "case": case,
            "window": window,
            "line_cutoff": line_cutoff,
            "min_line_intensity": min_line_intensity,
        },
    )
    return result


def build_hitran_table(*args, **kwargs):
    raise NotImplementedError(
        "build_hitran_table has not been ported into the first framework pass yet. "
        "Use scripts/build_exomol_ch4_mm_hitran_db.py as the archived reference."
    )


def extract_bands(*args, **kwargs):
    raise NotImplementedError(
        "ExoMol band-text extraction has not been ported into the first framework pass yet. "
        "Use scripts/extract_exomol_ch4_mm_pure_nu3_band_texts.py as the archived reference."
    )
