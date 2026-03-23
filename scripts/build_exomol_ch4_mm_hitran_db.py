"""
Build a HITRAN-style fixed-width database from the local ExoMol CH4 MM dataset.

This converter streams the ExoMol `.states/.trans/.pf/.def` files and writes a
HAPI-compatible pair:

- `hitran_db/{table_name}.data`
- `hitran_db/{table_name}.header`

What is mapped directly from ExoMol:
- `nu`: upper-state energy minus lower-state energy
- `sw`: LTE line intensity at `--reference-temperature-k`
- `a`: Einstein-A
- `elower`: lower-state energy
- `gp`, `gpp`: upper and lower state degeneracies
- `global_*_quanta`: compact CH4 vibrational labels built from `n1 n2 n3 n4`
  plus Td symmetry
- `local_*_quanta`: compact CH4 state labels built from `J`, `Grot`, and `irot`

What is not present in ExoMol and therefore filled with defaults/assumptions:
- `gamma_air`: ExoMol default Lorentz HWHM from the `.def` file unless
  overridden
- `gamma_self`: defaults to `gamma_air` unless overridden
- `n_air`: ExoMol default temperature exponent from the `.def` file
- `delta_air`: `0.0`
- `ierr`, `iref`, `line_mixing_flag`: placeholders

This script is specific to the CH4 `12C-1H4__MM` dataset currently stored in
this repo. It is designed to scale to all locally available transition chunks,
but the full build can still be very large and slow.
"""

from __future__ import annotations

import argparse
import bz2
import copy
import json
import math
import re
from pathlib import Path
from typing import Iterator

import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = ROOT_DIR / "exomol_db" / "CH4" / "12C-1H4" / "MM"
DEFAULT_OUTPUT_DIR = ROOT_DIR / "hitran_db"
DEFAULT_TABLE_NAME = "CH4_EXOMOL_MM_I1"
REFERENCE_TEMPERATURE_K = 296.0
SECOND_RADIATION_CONSTANT_CM_K = 1.438776877
LIGHT_SPEED_CM_S = 2.99792458e10

HITRAN_STYLE_MAP = {
    "A1": "1A1",
    "A2": "1A2",
    "E": "1E",
    "F1": "1F1",
    "F2": "1F2",
    "T1": "1F1",
    "T2": "1F2",
}

HITRAN_HEADER_TEMPLATE = {
    "table_name": "###",
    "number_of_rows": -1,
    "format": {
        "a": "%10.3E",
        "gamma_air": "%5.4f",
        "gp": "%7.1f",
        "local_iso_id": "%1d",
        "molec_id": "%2d",
        "sw": "%10.3E",
        "local_lower_quanta": "%15s",
        "local_upper_quanta": "%15s",
        "gpp": "%7.1f",
        "elower": "%10.4f",
        "n_air": "%4.2f",
        "delta_air": "%8.6f",
        "global_upper_quanta": "%15s",
        "iref": "%12s",
        "line_mixing_flag": "%1s",
        "ierr": "%6s",
        "nu": "%12.6f",
        "gamma_self": "%5.3f",
        "global_lower_quanta": "%15s",
    },
    "default": {
        "a": 0.0,
        "gamma_air": 0.0,
        "gp": "FFF",
        "local_iso_id": 0,
        "molec_id": 0,
        "sw": 0.0,
        "local_lower_quanta": "000",
        "local_upper_quanta": "000",
        "gpp": "FFF",
        "elower": 0.0,
        "n_air": 0.0,
        "delta_air": 0.0,
        "global_upper_quanta": "000",
        "iref": "EEE",
        "line_mixing_flag": "E",
        "ierr": "EEE",
        "nu": 0.0,
        "gamma_self": 0.0,
        "global_lower_quanta": "000",
    },
    "table_type": "column-fixed",
    "size_in_bytes": -1,
    "order": [
        "molec_id",
        "local_iso_id",
        "nu",
        "sw",
        "a",
        "gamma_air",
        "gamma_self",
        "elower",
        "n_air",
        "delta_air",
        "global_upper_quanta",
        "global_lower_quanta",
        "local_upper_quanta",
        "local_lower_quanta",
        "ierr",
        "iref",
        "line_mixing_flag",
        "gp",
        "gpp",
    ],
    "description": {
        "a": "Einstein A-coefficient in s-1",
        "gamma_air": "Air-broadened Lorentzian HWHM at p = 1 atm and T = 296 K; ExoMol default when unavailable",
        "gp": "Upper state degeneracy",
        "local_iso_id": "Local isotopologue ID",
        "molec_id": "HITRAN molecule ID",
        "sw": "Line intensity, multiplied by isotopologue abundance, at T = 296 K",
        "local_lower_quanta": "Compact ExoMol-derived lower-state label using J, Grot, and irot",
        "local_upper_quanta": "Compact ExoMol-derived upper-state label using J, Grot, and irot",
        "gpp": "Lower state degeneracy",
        "elower": "Lower-state energy",
        "n_air": "Temperature exponent for the air-broadened HWHM; ExoMol default when unavailable",
        "delta_air": "Pressure shift induced by air; set to 0.0 because ExoMol MM does not provide it here",
        "global_upper_quanta": "Compact ExoMol-derived upper-state vibrational label using n1 n2 n3 n4 and symmetry",
        "iref": "Placeholder reference identifier string",
        "line_mixing_flag": "Placeholder line-mixing flag",
        "ierr": "Placeholder uncertainty string",
        "nu": "Transition wavenumber",
        "gamma_self": "Self-broadened HWHM at 1 atm and 296 K; defaults to gamma_air here",
        "global_lower_quanta": "Compact ExoMol-derived lower-state vibrational label using n1 n2 n3 n4 and symmetry",
    },
    "position": {
        "molec_id": 0,
        "local_iso_id": 2,
        "nu": 3,
        "sw": 15,
        "a": 25,
        "gamma_air": 35,
        "gamma_self": 40,
        "elower": 45,
        "n_air": 55,
        "delta_air": 59,
        "global_upper_quanta": 67,
        "global_lower_quanta": 82,
        "local_upper_quanta": 97,
        "local_lower_quanta": 112,
        "ierr": 127,
        "iref": 133,
        "line_mixing_flag": 145,
        "gp": 146,
        "gpp": 153,
    },
    "extra": [],
    "extra_format": {},
    "extra_separator": ",",
}

FIELD_FORMAT_RE = re.compile(r"%(\d+)(?:\.\d+)?([A-Za-z])")
TRANSITION_FILE_RE_TEMPLATE = r"^{dataset}__(\d+)-(\d+)\.trans(?:\.bz2)?$"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a HITRAN-style HAPI table from the local ExoMol CH4 MM dataset.",
    )
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR, help="Directory containing the ExoMol files.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help="Directory for the HITRAN-style table.")
    parser.add_argument(
        "--dataset-stem",
        default=None,
        help="ExoMol dataset stem such as 12C-1H4__MM. Defaults to the only .def file found in --data-dir.",
    )
    parser.add_argument("--table-name", default=DEFAULT_TABLE_NAME, help="Output table stem inside --output-dir.")
    parser.add_argument("--molec-id", type=int, default=6, help="HITRAN molecule ID. Default 6 for CH4.")
    parser.add_argument("--local-iso-id", type=int, default=1, help="Local isotopologue ID.")
    parser.add_argument(
        "--reference-temperature-k",
        type=float,
        default=REFERENCE_TEMPERATURE_K,
        help="Reference temperature used for sw in Kelvin.",
    )
    parser.add_argument(
        "--intensity-threshold",
        type=float,
        default=0.0,
        help="Discard transitions weaker than this LTE intensity in cm/molecule.",
    )
    parser.add_argument("--abundance", type=float, default=1.0, help="Scale sw by this isotopologue abundance factor.")
    parser.add_argument("--wn-min", type=float, default=None, help="Optional minimum wavenumber in cm^-1.")
    parser.add_argument("--wn-max", type=float, default=None, help="Optional maximum wavenumber in cm^-1.")
    parser.add_argument(
        "--assume-ordered",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Assume transitions are approximately ordered by wavenumber within each file so scans can break early.",
    )
    parser.add_argument(
        "--enforce-ordering",
        action="store_true",
        default=False,
        help="Fail if transition wavenumbers decrease within a file by more than 1e-2 cm^-1.",
    )
    parser.add_argument(
        "--break-margin-cm",
        type=float,
        default=10.0,
        help="When --assume-ordered is on, stop a file after nu exceeds wn_max plus this margin.",
    )
    parser.add_argument(
        "--gamma-air",
        type=float,
        default=None,
        help="Override gamma_air. Defaults to the ExoMol .def value.",
    )
    parser.add_argument(
        "--gamma-self",
        type=float,
        default=None,
        help="Override gamma_self. Defaults to gamma_air because ExoMol MM does not provide a separate value here.",
    )
    parser.add_argument(
        "--delta-air",
        type=float,
        default=0.0,
        help="Air pressure shift to store in the HITRAN-style table. Default 0.0.",
    )
    parser.add_argument(
        "--max-transitions",
        type=int,
        default=None,
        help="Optional hard stop after parsing this many transition rows. Useful for smoke tests.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1_000_000,
        help="Print progress every N parsed transition rows.",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.reference_temperature_k <= 0.0:
        raise ValueError("--reference-temperature-k must be positive")
    if args.intensity_threshold < 0.0:
        raise ValueError("--intensity-threshold must be non-negative")
    if args.abundance < 0.0:
        raise ValueError("--abundance must be non-negative")
    if (args.wn_min is None) != (args.wn_max is None):
        raise ValueError("--wn-min and --wn-max must be provided together")
    if args.wn_min is not None and args.wn_max is not None and args.wn_max <= args.wn_min:
        raise ValueError("--wn-max must be greater than --wn-min")
    if args.break_margin_cm < 0.0:
        raise ValueError("--break-margin-cm must be non-negative")
    if args.max_transitions is not None and args.max_transitions <= 0:
        raise ValueError("--max-transitions must be positive")
    if args.progress_every <= 0:
        raise ValueError("--progress-every must be positive")


def autodetect_dataset_stem(data_dir: Path) -> str:
    def_paths = sorted(data_dir.glob("*.def"))
    if len(def_paths) != 1:
        raise RuntimeError(
            f"Expected exactly one .def file in {data_dir}, found {len(def_paths)}; pass --dataset-stem explicitly"
        )
    return def_paths[0].stem


def iter_text_lines(path: Path, chunk_size: int = 8 * 1024 * 1024) -> Iterator[str]:
    if path.suffix != ".bz2":
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                yield line.rstrip("\r\n")
        return

    decompressor = bz2.BZ2Decompressor()
    pending = b""
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            data = decompressor.decompress(chunk)
            if not data:
                continue
            pending += data
            lines = pending.split(b"\n")
            pending = lines.pop()
            for line in lines:
                yield line.decode("utf-8").rstrip("\r")
        if pending:
            yield pending.decode("utf-8").rstrip("\r")


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

            if comment.startswith("No. of states in .states file"):
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

    required_keys = {"nstates", "gamma0", "n_exponent", "has_uncertainty", "has_lifetime", "has_lande_g"}
    missing = required_keys - metadata.keys()
    if missing:
        raise RuntimeError(f"Missing required metadata in {def_path}: {sorted(missing)}")

    quantum_labels = metadata["quantum_labels"]
    if not quantum_labels:
        raise RuntimeError(f"No quantum labels found in {def_path}")

    metadata["ordered_quantum_labels"] = [quantum_labels[index] for index in sorted(quantum_labels)]
    metadata["ordered_auxiliary_titles"] = [
        metadata["auxiliary_titles"][index] for index in sorted(metadata["auxiliary_titles"])
    ]
    return metadata


def states_column_indices(def_metadata: dict[str, object]) -> dict[str, int]:
    base_columns = ["i", "E", "gtot", "J"]
    if bool(def_metadata["has_uncertainty"]):
        base_columns.append("unc")
    if bool(def_metadata["has_lifetime"]):
        base_columns.append("lifetime")
    if bool(def_metadata["has_lande_g"]):
        base_columns.append("lande_g")

    all_columns = base_columns + list(def_metadata["ordered_quantum_labels"]) + list(
        def_metadata["ordered_auxiliary_titles"]
    )
    return {name: index for index, name in enumerate(all_columns)}


def load_partition_function(pf_path: Path) -> tuple[np.ndarray, np.ndarray]:
    pf_table = np.loadtxt(pf_path)
    if pf_table.ndim != 2 or pf_table.shape[1] < 2:
        raise RuntimeError(f"Unexpected partition function format in {pf_path}")
    return pf_table[:, 0], pf_table[:, 1]


def interpolate_partition_function(
    temperatures: np.ndarray,
    partition_values: np.ndarray,
    temperature_k: float,
) -> float:
    if temperature_k < temperatures[0] or temperature_k > temperatures[-1]:
        raise ValueError(
            f"Temperature {temperature_k:g} K is outside the partition-function range "
            f"{temperatures[0]:g}-{temperatures[-1]:g} K"
        )
    return float(np.interp(temperature_k, temperatures, partition_values))


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


def to_hitran_symmetry(label: str) -> str:
    return HITRAN_STYLE_MAP.get(label.strip(), label.strip() or "000")


def format_j_value(value: str) -> str:
    j_value = float(value)
    if abs(j_value - round(j_value)) < 1.0e-9:
        return str(int(round(j_value)))
    return f"{j_value:g}"


def compact_global_quanta(parts: list[str], column_indices: dict[str, int]) -> str:
    return (
        f"{int(parts[column_indices['n1']])} "
        f"{int(parts[column_indices['n2']])} "
        f"{int(parts[column_indices['n3']])} "
        f"{int(parts[column_indices['n4']])} "
        f"{to_hitran_symmetry(parts[column_indices['Gtot']])}"
    )


def compact_local_quanta(parts: list[str], column_indices: dict[str, int]) -> str:
    if "Grot" not in column_indices or "irot" not in column_indices:
        return "000"
    j_text = format_j_value(parts[column_indices["J"]])
    grot = parts[column_indices["Grot"]].strip() or "000"
    irot = int(parts[column_indices["irot"]])
    return f"J{j_text} {grot} {irot}"


def load_state_arrays(
    states_path: Path,
    nstates: int,
    column_indices: dict[str, int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[str], list[str]]:
    required_labels = {"n1", "n2", "n3", "n4", "Gtot"}
    missing_labels = required_labels - column_indices.keys()
    if missing_labels:
        raise RuntimeError(f"Missing required state columns: {sorted(missing_labels)}")

    energies = np.empty(nstates + 1, dtype=np.float64)
    g_totals = np.empty(nstates + 1, dtype=np.float64)
    global_label_ids = np.zeros(nstates + 1, dtype=np.int32)
    local_label_ids = np.zeros(nstates + 1, dtype=np.int32)
    energies.fill(np.nan)
    g_totals.fill(np.nan)

    global_labels = ["000"]
    local_labels = ["000"]
    global_label_to_id = {"000": 0}
    local_label_to_id = {"000": 0}

    for line_number, raw_line in enumerate(iter_text_lines(states_path), start=1):
        parts = raw_line.split()
        if len(parts) <= max(column_indices.values()):
            continue

        state_id = int(parts[column_indices["i"]])
        energies[state_id] = float(parts[column_indices["E"]])
        g_totals[state_id] = float(parts[column_indices["gtot"]])

        global_label = compact_global_quanta(parts, column_indices)
        local_label = compact_local_quanta(parts, column_indices)

        global_label_id = global_label_to_id.get(global_label)
        if global_label_id is None:
            global_label_id = len(global_labels)
            global_labels.append(global_label)
            global_label_to_id[global_label] = global_label_id
        global_label_ids[state_id] = global_label_id

        local_label_id = local_label_to_id.get(local_label)
        if local_label_id is None:
            local_label_id = len(local_labels)
            local_labels.append(local_label)
            local_label_to_id[local_label] = local_label_id
        local_label_ids[state_id] = local_label_id

        if line_number % 1_000_000 == 0:
            print(f"loaded {line_number:,} states")

    if np.isnan(energies[1:]).any() or np.isnan(g_totals[1:]).any():
        raise RuntimeError(f"State table {states_path} did not fill all expected state IDs")

    return energies, g_totals, global_label_ids, local_label_ids, global_labels, local_labels


def parse_transition_range(path: Path, dataset_stem: str) -> tuple[float | None, float | None]:
    pattern = re.compile(TRANSITION_FILE_RE_TEMPLATE.format(dataset=re.escape(dataset_stem)))
    match = pattern.match(path.name)
    if match is None:
        return None, None
    return float(match.group(1)), float(match.group(2))


def discover_transition_files(data_dir: Path, dataset_stem: str) -> list[Path]:
    candidates = sorted(list(data_dir.glob(f"{dataset_stem}*.trans")) + list(data_dir.glob(f"{dataset_stem}*.trans.bz2")))
    unique_paths: list[Path] = []
    seen: set[Path] = set()
    for path in candidates:
        if path not in seen:
            unique_paths.append(path)
            seen.add(path)
    if not unique_paths:
        raise FileNotFoundError(f"No transition files found for {dataset_stem} in {data_dir}")
    return sorted(unique_paths, key=lambda path: (parse_transition_range(path, dataset_stem), path.name))


def select_transition_files(
    all_paths: list[Path],
    dataset_stem: str,
    wn_min: float | None,
    wn_max: float | None,
) -> list[Path]:
    if wn_min is None or wn_max is None:
        return all_paths

    selected: list[Path] = []
    for path in all_paths:
        start_cm, end_cm = parse_transition_range(path, dataset_stem)
        if start_cm is None or end_cm is None:
            selected.append(path)
            continue
        if end_cm <= wn_min:
            continue
        if start_cm >= wn_max:
            continue
        selected.append(path)
    return selected


def parse_fixed_width(format_string: str) -> tuple[int, str]:
    match = FIELD_FORMAT_RE.fullmatch(format_string)
    if match is None:
        raise ValueError(f"Unsupported fixed-width format: {format_string}")
    return int(match.group(1)), match.group(2)


def format_field(value: object, format_string: str) -> str:
    width, kind = parse_fixed_width(format_string)
    if kind.lower() == "s":
        text = str(value)
        if len(text) > width:
            text = text[:width]
        return text.rjust(width)

    text = format_string % value
    # Some local HITRAN tables use a five-character gamma field encoded as
    # `.0650` rather than `0.0650`, so mirror that storage convention.
    if len(text) == width + 1 and text.startswith("0"):
        text = text[1:]
    elif len(text) == width + 2 and text.startswith("-0"):
        text = "-" + text[2:]
    if len(text) > width:
        raise ValueError(f"Formatted value {text!r} exceeds width {width} for {format_string}")
    return text


def format_hitran_row(row: dict[str, object]) -> str:
    return "".join(
        format_field(row[name], HITRAN_HEADER_TEMPLATE["format"][name]) for name in HITRAN_HEADER_TEMPLATE["order"]
    )


def build_header(
    *,
    table_name: str,
    number_of_rows: int,
    size_in_bytes: int,
) -> dict[str, object]:
    header = copy.deepcopy(HITRAN_HEADER_TEMPLATE)
    header["table_name"] = table_name
    header["number_of_rows"] = number_of_rows
    header["size_in_bytes"] = size_in_bytes
    return header


def main() -> None:
    args = parse_args()
    validate_args(args)

    dataset_stem = args.dataset_stem or autodetect_dataset_stem(args.data_dir)
    def_path = args.data_dir / f"{dataset_stem}.def"
    pf_path = args.data_dir / f"{dataset_stem}.pf"
    states_path = args.data_dir / f"{dataset_stem}.states.bz2"
    if not states_path.exists():
        states_path = args.data_dir / f"{dataset_stem}.states"

    args.output_dir.mkdir(parents=True, exist_ok=True)
    data_path = args.output_dir / f"{args.table_name}.data"
    header_path = args.output_dir / f"{args.table_name}.header"

    def_metadata = parse_exomol_def(def_path)
    column_indices = states_column_indices(def_metadata)
    pf_temperatures, pf_values = load_partition_function(pf_path)
    partition_function = interpolate_partition_function(
        pf_temperatures,
        pf_values,
        args.reference_temperature_k,
    )

    gamma_air = float(def_metadata["gamma0"]) if args.gamma_air is None else args.gamma_air
    gamma_self = gamma_air if args.gamma_self is None else args.gamma_self
    n_air = float(def_metadata["n_exponent"])

    print(f"dataset stem: {dataset_stem}")
    print(f"loading states from {states_path.name}")
    energies, g_totals, global_label_ids, local_label_ids, global_labels, local_labels = load_state_arrays(
        states_path,
        int(def_metadata["nstates"]),
        column_indices,
    )

    all_transition_files = discover_transition_files(args.data_dir, dataset_stem)
    transition_files = select_transition_files(all_transition_files, dataset_stem, args.wn_min, args.wn_max)
    print(f"selected {len(transition_files)} transition file(s)")
    for path in transition_files:
        print(f"  {path.name}")

    stats = {
        "files": len(transition_files),
        "parsed": 0,
        "nonpositive": 0,
        "in_window": 0,
        "kept": 0,
    }
    with data_path.open("w", encoding="utf-8", newline="") as out_handle:
        for path in transition_files:
            print(f"scan {path.name}")
            last_nu = -1.0e300
            for raw_line in iter_text_lines(path):
                parts = raw_line.split()
                if len(parts) != 3:
                    continue

                upper_id = int(parts[0])
                lower_id = int(parts[1])
                a_value = float(parts[2])
                stats["parsed"] += 1

                wavenumber = float(energies[upper_id] - energies[lower_id])
                if wavenumber <= 0.0:
                    stats["nonpositive"] += 1
                    if args.max_transitions is not None and stats["parsed"] >= args.max_transitions:
                        break
                    continue

                if args.enforce_ordering and wavenumber < last_nu - 1.0e-2:
                    raise RuntimeError(
                        f"Transitions are not ordered by wavenumber in {path.name}: "
                        f"last={last_nu:.6f} current={wavenumber:.6f}"
                    )
                last_nu = wavenumber

                if args.wn_min is not None and wavenumber < args.wn_min:
                    if args.max_transitions is not None and stats["parsed"] >= args.max_transitions:
                        break
                    continue

                if args.wn_max is not None and wavenumber > args.wn_max:
                    if args.assume_ordered and wavenumber > args.wn_max + args.break_margin_cm:
                        break
                    if args.max_transitions is not None and stats["parsed"] >= args.max_transitions:
                        break
                    continue

                stats["in_window"] += 1

                sw = lte_line_intensity_cm_per_molecule(
                    wavenumber=wavenumber,
                    a_coefficient=a_value,
                    g_upper=float(g_totals[upper_id]),
                    lower_energy_cm=float(energies[lower_id]),
                    temperature_k=args.reference_temperature_k,
                    partition_function=partition_function,
                )
                sw *= args.abundance
                if sw < args.intensity_threshold:
                    if args.max_transitions is not None and stats["parsed"] >= args.max_transitions:
                        break
                    continue

                row = {
                    "molec_id": args.molec_id,
                    "local_iso_id": args.local_iso_id,
                    "nu": wavenumber,
                    "sw": sw,
                    "a": a_value,
                    "gamma_air": gamma_air,
                    "gamma_self": gamma_self,
                    "elower": float(energies[lower_id]),
                    "n_air": n_air,
                    "delta_air": args.delta_air,
                    "global_upper_quanta": global_labels[int(global_label_ids[upper_id])],
                    "global_lower_quanta": global_labels[int(global_label_ids[lower_id])],
                    "local_upper_quanta": local_labels[int(local_label_ids[upper_id])],
                    "local_lower_quanta": local_labels[int(local_label_ids[lower_id])],
                    "ierr": "EEE",
                    "iref": "EXOMOL-MM",
                    "line_mixing_flag": "N",
                    "gp": float(g_totals[upper_id]),
                    "gpp": float(g_totals[lower_id]),
                }
                out_handle.write(format_hitran_row(row))
                out_handle.write("\n")
                stats["kept"] += 1

                if stats["parsed"] % args.progress_every == 0:
                    print(
                        f"  parsed {stats['parsed']:,} transitions, "
                        f"window {stats['in_window']:,}, kept {stats['kept']:,}"
                    )

                if args.max_transitions is not None and stats["parsed"] >= args.max_transitions:
                    break

            if args.max_transitions is not None and stats["parsed"] >= args.max_transitions:
                print(f"stopped early after {stats['parsed']:,} parsed transitions because --max-transitions was reached")
                break

    if stats["kept"] == 0:
        raise RuntimeError("No transitions were written; adjust the current spectral or intensity filters.")

    header = build_header(
        table_name=args.table_name,
        number_of_rows=stats["kept"],
        size_in_bytes=data_path.stat().st_size,
    )
    header_path.write_text(json.dumps(header, indent=2), encoding="utf-8")

    print()
    print(f"Q({args.reference_temperature_k:g} K) = {partition_function:.8g}")
    print(f"parsed transitions: {stats['parsed']:,}")
    print(f"positive transitions in window: {stats['in_window']:,}")
    print(f"written rows: {stats['kept']:,}")
    print(f"skipped nonpositive transitions: {stats['nonpositive']:,}")
    print(f"saved {data_path}")
    print(f"saved {header_path}")


if __name__ == "__main__":
    main()
