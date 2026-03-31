"""
Extract pure-CH4-nu3 ExoMol MM transitions into per-band text files.

This script is intentionally modeled after `extract_ch4_nu3_band_texts.py`,
but it works directly from the ExoMol MM `.states/.trans/.pf/.def` files.

Inputs
------
- `--data-dir/{DATASET_STEM}.def`: ExoMol definition/metadata file used to
  discover the state-table column layout and required flags.
- `--data-dir/{DATASET_STEM}.pf`: partition-function table used to interpolate
  the partition function at `--reference-temperature-k`.
- `--data-dir/{DATASET_STEM}.states.bz2`: compressed state table containing
  energies, degeneracies, and vibrational/symmetry labels for each state ID.
- `--data-dir/{DATASET_STEM}__*.trans.bz2`: one or more compressed transition
  chunks containing `upper_id lower_id Einstein_A`; the script scans only the
  chunks overlapping `--wn-min/--wn-max` when a window is provided.
- Command-line filters: `--reference-temperature-k`,
  `--intensity-threshold`, `--wn-min`, `--wn-max`, and
  `--require-unit-step` control which transitions are exported.

Definition of "pure nu3" used here
----------------------------------
- lower and upper states must both satisfy `n1 = n2 = n4 = 0`
- the transition must be upward in `n3`
- optionally, `--require-unit-step` restricts the progression to `n3+1`

Important note about intensity
------------------------------
Unlike HITRAN, ExoMol does not provide a ready-made reference-temperature line
intensity in each transition row. This script computes one using the requested
`--reference-temperature-k` and writes that value to the exported text rows.

Outputs
-------
- `exomol_ch4_mm_pure_nu3_band_texts/*.txt`: one tab-delimited text file per
  exact ExoMol band pair
- `exomol_ch4_mm_pure_nu3_band_texts/exomol_pure_nu3_band_text_summary.csv`:
  manifest CSV with one row per exported band and its transition count

Band labeling used by this extractor
------------------------------------
The MM `.states` schema defines many ExoMol quantum labels plus auxiliary
titles. This extractor now preserves the full labeled MM state signature in the
exported band label:

- leading mode counters: `n1`, `n2`, `n3`, `n4`
- total symmetry: `Gtot`
- all remaining MM quantum labels from `.def`
- all auxiliary titles from `.def`

The only state fields not included in the grouping key are the standard base
columns such as `i`, `E`, `gtot`, `J`, `unc`, and `lifetime`, because those are
per-state bookkeeping/physical values rather than the named MM label set.
"""

from __future__ import annotations

import argparse
import bz2
import csv
import hashlib
import math
import re
from array import array
from bisect import bisect_left
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

ROOT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT_DIR / "exomol_db" / "CH4" / "12C-1H4" / "MM"
DATASET_STEM = "12C-1H4__MM"
SECOND_RADIATION_CONSTANT_CM_K = 1.438776877
LIGHT_SPEED_CM_S = 2.99792458e10


OUTPUT_DIR = ROOT_DIR / "exomol_ch4_mm_pure_nu3_band_texts"
SUMMARY_CSV_NAME = "exomol_pure_nu3_band_text_summary.csv"
MAX_OPEN_OUTPUTS = 128
DATA_COLUMNS = [
    "upper_id",
    "lower_id",
    "upper_energy_cm-1",
    "lower_energy_cm-1",
    "wavenumber_cm-1",
    "einstein_A_s-1",
    "line_intensity_cm_per_molecule",
    "local_upper_quanta",
    "local_lower_quanta",
]
DATA_HEADER = "\t".join(DATA_COLUMNS)
HITRAN_STYLE_MAP = {
    "A1": "1A1",
    "A2": "1A2",
    "E": "1E",
    "F1": "1F1",
    "F2": "1F2",
    "T1": "1F1",
    "T2": "1F2",
}


def iter_bz2_text_lines(path: Path, chunk_size: int = 8 * 1024 * 1024):
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


def transition_filename(start_cm: int) -> str:
    return f"{DATASET_STEM}__{start_cm:05d}-{start_cm + 100:05d}.trans.bz2"


def available_transition_starts(data_dir: Path) -> list[int]:
    starts: list[int] = []
    pattern = re.compile(rf"^{re.escape(DATASET_STEM)}__(\d{{5}})-(\d{{5}})\.trans\.bz2$")
    for path in data_dir.glob(f"{DATASET_STEM}__*.trans.bz2"):
        match = pattern.match(path.name)
        if match is None:
            continue
        start_cm = int(match.group(1))
        end_cm = int(match.group(2))
        if end_cm != start_cm + 100:
            continue
        starts.append(start_cm)
    if not starts:
        raise FileNotFoundError(f"No ExoMol transition files found in {data_dir}")
    return sorted(starts)


def overlapping_transition_files(data_dir: Path, wn_min: float, wn_max: float, wing: float) -> tuple[list[Path], float]:
    available_starts = available_transition_starts(data_dir)
    available_min = available_starts[0]
    available_max = available_starts[-1] + 100

    start_chunk = int(math.floor(max(0.0, wn_min - wing) / 100.0) * 100)
    end_chunk = int(math.ceil((wn_max + wing) / 100.0) * 100)

    clipped_start = max(start_chunk, available_min)
    clipped_end = min(end_chunk, available_max)
    if clipped_end <= clipped_start:
        raise FileNotFoundError(f"No available transition files overlap {wn_min:g}-{wn_max:g} cm^-1 in {data_dir}")

    effective_wing = wing
    if clipped_start > start_chunk:
        effective_wing = min(effective_wing, max(0.0, wn_min - clipped_start))
    if clipped_end < end_chunk:
        effective_wing = min(effective_wing, max(0.0, clipped_end - wn_max))

    files = []
    for start_cm in range(clipped_start, clipped_end, 100):
        path = data_dir / transition_filename(start_cm)
        if not path.exists():
            raise FileNotFoundError(f"Missing transition file: {path}")
        files.append(path)
    return files, effective_wing


def load_partition_function(pf_path: Path) -> tuple[list[float], list[float]]:
    temperatures: list[float] = []
    partition_values: list[float] = []
    with pf_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            parts = raw_line.split()
            if len(parts) < 2:
                continue
            temperatures.append(float(parts[0]))
            partition_values.append(float(parts[1]))

    if not temperatures:
        raise RuntimeError(f"Unexpected partition function format in {pf_path}")
    return temperatures, partition_values


def interpolate_partition_function(temperatures: list[float], partition_values: list[float], temperature_k: float) -> float:
    if temperature_k < temperatures[0] or temperature_k > temperatures[-1]:
        raise ValueError(
            f"Temperature {temperature_k:g} K is outside the partition-function range "
            f"{temperatures[0]:g}-{temperatures[-1]:g} K"
        )

    index = bisect_left(temperatures, temperature_k)
    if index < len(temperatures) and temperatures[index] == temperature_k:
        return partition_values[index]
    if index == 0 or index == len(temperatures):
        raise ValueError(f"Cannot interpolate partition function at {temperature_k:g} K")

    lower_temperature = temperatures[index - 1]
    upper_temperature = temperatures[index]
    lower_value = partition_values[index - 1]
    upper_value = partition_values[index]
    fraction = (temperature_k - lower_temperature) / (upper_temperature - lower_temperature)
    return lower_value + fraction * (upper_value - lower_value)


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


@dataclass(frozen=True)
class BandSignature:
    """Full labeled ExoMol MM state signature used to define an exported band.

    Fields kept explicitly for filtering/grouping:
    - `n1`: A1-symmetry normal-mode quantum number.
    - `n2`: E-symmetry normal-mode quantum number.
    - `n3`: F1-symmetry normal-mode quantum number; this is the `nu3` mode that
      the script tracks for pure-`nu3` progressions.
    - `n4`: F2-symmetry normal-mode quantum number.
    - `gtot_sym`: ExoMol quantum label `Gtot`, the total Td symmetry label of
      the state. This is not the numeric lower-case `gtot` degeneracy column.
    - `other_labels`: ordered `(name, value)` pairs covering every remaining MM
      quantum label and auxiliary title defined in `.def`, preserving the values
      from the `.states` row. This includes labels such as `P`, `Pnum`, `L2`,
      `L3`, `M3`, `L4`, `M4`, `Gvib`, `irot`, `Grot`, `ivib`, `Coef`, `v1`-`v9`,
      `SourceType`, and `Ecal`.

    Base state fields like `i`, `E`, `gtot`, `J`, `unc`, and `lifetime` are not
    stored here; those are standard state-table columns rather than MM labels,
    and including them would turn the export into a near one-state-per-band
    grouping instead of a label-defined band grouping.
    """
    n1: int
    n2: int
    n3: int
    n4: int
    gtot_sym: str
    other_labels: tuple[tuple[str, str], ...]


# ExoMol MM state labels defined in `12C-1H4__MM.def`:
# - `Gtot`: total Td symmetry of the full state.
# - `P`: polyad number.
# - `Pnum`: polyad counting number.
# - `n1`: A1-symmetry normal-mode quantum number.
# - `n2`: E-symmetry normal-mode quantum number.
# - `L2`: vibrational angular-momentum quantum number for mode 2.
# - `n3`: F1-symmetry normal-mode quantum number.
# - `L3`: vibrational angular-momentum quantum number for mode 3.
# - `M3`: multiplicity index for mode 3.
# - `n4`: F2-symmetry normal-mode quantum number.
# - `L4`: vibrational angular-momentum quantum number for mode 4.
# - `M4`: multiplicity index for mode 4.
# - `Gvib`: vibrational symmetry.
# - `irot`: index in the rotational symmetry block.
# - `Grot`: rotational symmetry.
# - `ivib`: index in the vibrational symmetry block.
# - `Coef`: coefficient of the largest basis-function contribution.
# - `v1`-`v9`: local-mode vibrational quantum numbers.
# - auxiliary `SourceType`: state provenance flag
#   (`Ma` MARVEL, `Ca` calculated, `EH` effective Hamiltonian,
#   `IE` isotopologue extrapolation).
# - auxiliary `Ecal`: calculated energy in cm^-1.
# This extractor keeps `n1`, `n2`, `n3`, `n4`, and `Gtot` explicitly, then
# stores every remaining MM quantum label and auxiliary title in
# `BandSignature.other_labels` in the order defined by `.def`.


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract pure-nu3 CH4 ExoMol MM transitions into per-band text files.",
    )
    # Runtime inputs: source dataset location plus export/intensity/window filters.
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR, help="Directory containing the ExoMol MM files.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Directory for exported text files.")
    parser.add_argument(
        "--reference-temperature-k",
        type=float,
        default=296.0,
        help="Reference temperature used to compute line intensity in cm/molecule.",
    )
    parser.add_argument(
        "--intensity-threshold",
        type=float,
        default=0.0,
        help="Discard transitions weaker than this reference-temperature line intensity.",
    )
    parser.add_argument("--wn-min", type=float, default=3000, help="Optional minimum wavenumber in cm^-1.")
    parser.add_argument("--wn-max", type=float, default=3000.001, help="Optional maximum wavenumber in cm^-1.")
    parser.add_argument(
        "--require-unit-step",
        action="store_true",
        default=False,
        help="Require a unit-step pure nu3 progression, e.g. 0->1 or 1->2.",
    )
    parser.add_argument(
        "--group-by",
        choices=("exact", "hitran-style"),
        default="exact",
        help="Group exported files by the exact MM state signature or directly by the coarse HITRAN-style band label.",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.reference_temperature_k <= 0.0:
        raise ValueError("--reference-temperature-k must be positive")
    if args.intensity_threshold < 0.0:
        raise ValueError("--intensity-threshold must be non-negative")
    if (args.wn_min is None) != (args.wn_max is None):
        raise ValueError("--wn-min and --wn-max must be provided together")
    if args.wn_min is not None and args.wn_max is not None and args.wn_max <= args.wn_min:
        raise ValueError("--wn-max must be greater than --wn-min")


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

    required_keys = {"nstates", "has_uncertainty", "has_lifetime", "has_lande_g"}
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


def signature_label_names(def_metadata: dict[str, object]) -> list[str]:
    excluded = {"Gtot", "n1", "n2", "n3", "n4"}
    ordered_names = list(def_metadata["ordered_quantum_labels"]) + list(def_metadata["ordered_auxiliary_titles"])
    return [name for name in ordered_names if name not in excluded]


def build_signature(
    parts: list[str],
    column_indices: dict[str, int],
    extra_label_names: list[str],
) -> BandSignature:
    return BandSignature(
        n1=int(parts[column_indices["n1"]]),
        n2=int(parts[column_indices["n2"]]),
        n3=int(parts[column_indices["n3"]]),
        n4=int(parts[column_indices["n4"]]),
        gtot_sym=parts[column_indices["Gtot"]],
        other_labels=tuple((name, parts[column_indices[name]]) for name in extra_label_names),
    )


def exomol_signature_label(signature: BandSignature) -> str:
    tail = " ".join(f"{name}={value}" for name, value in signature.other_labels)
    label = f"{signature.n1} {signature.n2} {signature.n3} {signature.n4} Gtot={signature.gtot_sym}"
    if tail:
        label += f" {tail}"
    return label


def exomol_coarse_signature_label(signature: BandSignature) -> str:
    return f"{signature.n1} {signature.n2} {signature.n3} {signature.n4} {signature.gtot_sym}"


def hitran_style_signature_label(signature: BandSignature) -> str:
    hitran_symmetry = HITRAN_STYLE_MAP.get(signature.gtot_sym, signature.gtot_sym)
    return f"{signature.n1} {signature.n2} {signature.n3} {signature.n4} {hitran_symmetry}"


def safe_label_fragment(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_") or "000"


def format_j_value(value: str) -> str:
    j_value = float(value)
    if abs(j_value - round(j_value)) < 1.0e-9:
        return str(int(round(j_value)))
    return f"{j_value:g}"


def compact_local_quanta(parts: list[str], column_indices: dict[str, int]) -> str:
    required_labels = {"Grot", "irot"}
    missing_labels = required_labels - column_indices.keys()
    if missing_labels:
        raise RuntimeError(f"Missing required state columns for local quanta: {sorted(missing_labels)}")

    j_text = format_j_value(parts[column_indices["J"]])
    grot = parts[column_indices["Grot"]].strip() or "000"
    irot = int(parts[column_indices["irot"]])
    return f"J{j_text} {grot} {irot}"


def float_array(length: int) -> array:
    return array("d", [math.nan]) * length


def int_array(length: int) -> array:
    return array("I", [0]) * length


def all_values_present(values: array) -> bool:
    for value in values[1:]:
        if math.isnan(value):
            return False
    return True


def signature_stem(signature: BandSignature) -> str:
    # Full MM labels are too verbose for a safe Windows path, so filenames keep
    # the core mode counters plus a stable digest of the full signature label.
    digest = hashlib.sha1(exomol_signature_label(signature).encode("utf-8")).hexdigest()[:12]
    return (
        f"{signature.n1}_{signature.n2}_{signature.n3}_{signature.n4}_"
        f"Gtot_{safe_label_fragment(signature.gtot_sym)}_{digest}"
    )


def output_path_for_band(
    output_dir: Path,
    lower_signature: BandSignature,
    upper_signature: BandSignature,
    reference_temperature_k: float,
) -> Path:
    filename = (
        "EXOMOL_CH4_MM_pure_nu3_"
        f"{signature_stem(lower_signature)}_to_{signature_stem(upper_signature)}"
        f"_T{reference_temperature_k:g}K.txt"
    )
    return output_dir / filename


def output_path_for_hitran_style_band(
    output_dir: Path,
    band_label_hitran_style: str,
    reference_temperature_k: float,
) -> Path:
    lower_label, upper_label = [part.strip() for part in band_label_hitran_style.split("->", maxsplit=1)]
    filename = (
        "EXOMOL_CH4_MM_pure_nu3_hitran_style_"
        f"{safe_label_fragment(lower_label)}_to_{safe_label_fragment(upper_label)}"
        f"_T{reference_temperature_k:g}K.txt"
    )
    return output_dir / filename


def band_label_pair(
    lower_signature: BandSignature,
    upper_signature: BandSignature,
    group_by: str,
) -> tuple[str, str]:
    if group_by == "hitran-style":
        band_label = f"{exomol_coarse_signature_label(lower_signature)} -> {exomol_coarse_signature_label(upper_signature)}"
    else:
        band_label = f"{exomol_signature_label(lower_signature)} -> {exomol_signature_label(upper_signature)}"
    band_label_hitran_style = (
        f"{hitran_style_signature_label(lower_signature)} -> {hitran_style_signature_label(upper_signature)}"
    )
    return band_label, band_label_hitran_style


def load_state_arrays(
    states_path: Path,
    nstates: int,
    column_indices: dict[str, int],
    extra_label_names: list[str],
) -> tuple[array, array, array, list[BandSignature | None], array, list[str]]:
    energies = float_array(nstates + 1)
    g_totals = float_array(nstates + 1)
    signature_ids = int_array(nstates + 1)
    local_label_ids = int_array(nstates + 1)

    signature_lookup: list[BandSignature | None] = [None]
    signature_to_id: dict[BandSignature, int] = {}
    local_labels = ["000"]
    local_label_to_id = {"000": 0}

    required_labels = {"n1", "n2", "n3", "n4", "Gtot"} | set(extra_label_names)
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

        local_label = compact_local_quanta(parts, column_indices)
        local_label_id = local_label_to_id.get(local_label)
        if local_label_id is None:
            local_label_id = len(local_labels)
            local_labels.append(local_label)
            local_label_to_id[local_label] = local_label_id
        local_label_ids[state_id] = local_label_id

        signature = build_signature(parts, column_indices, extra_label_names)
        if signature.n1 != 0 or signature.n2 != 0 or signature.n4 != 0:
            continue

        signature_id = signature_to_id.get(signature)
        if signature_id is None:
            signature_id = len(signature_lookup)
            signature_lookup.append(signature)
            signature_to_id[signature] = signature_id
        signature_ids[state_id] = signature_id

        if line_number % 1_000_000 == 0:
            print(f"loaded {line_number:,} states")

    if not all_values_present(energies) or not all_values_present(g_totals):
        raise RuntimeError(f"State table {states_path} did not fill all expected state IDs")

    return energies, g_totals, signature_ids, signature_lookup, local_label_ids, local_labels


def transition_files_for_args(args: argparse.Namespace) -> list[Path]:
    if args.wn_min is not None and args.wn_max is not None:
        files, _ = overlapping_transition_files(args.data_dir, args.wn_min, args.wn_max, 0.0)
        return files

    print("warning: no wavenumber window was provided; scanning all available transition chunks")
    return [args.data_dir / transition_filename(start_cm) for start_cm in available_transition_starts(args.data_dir)]


def in_selected_window(wavenumber: float, wn_min: float | None, wn_max: float | None) -> bool:
    if wn_min is None or wn_max is None:
        return True
    return wn_min <= wavenumber <= wn_max


def is_allowed_progression(
    lower_signature: BandSignature,
    upper_signature: BandSignature,
    *,
    require_unit_step: bool,
) -> bool:
    if upper_signature.n3 <= lower_signature.n3:
        return False
    if require_unit_step and upper_signature.n3 != lower_signature.n3 + 1:
        return False
    return True


def write_file_preamble(
    handle: TextIO,
    *,
    band_label: str,
    band_label_hitran_style: str,
    reference_temperature_k: float,
    wn_min: float | None,
    wn_max: float | None,
) -> None:
    handle.write("# ExoMol MM pure-nu3 band export\n")
    handle.write(f"# dataset_stem={DATASET_STEM}\n")
    handle.write(f"# reference_temperature_k={reference_temperature_k:g}\n")
    handle.write("# pure_nu3_filter=lower and upper states satisfy n1=n2=n4=0\n")
    handle.write(f"# band_label={band_label}\n")
    handle.write(f"# hitran_style_band_label={band_label_hitran_style}\n")
    if wn_min is not None and wn_max is not None:
        handle.write(f"# wavenumber_window_cm-1={wn_min:g} to {wn_max:g}\n")
    handle.write(f"{DATA_HEADER}\n")


def write_summary_csv(output_dir: Path, rows: list[dict[str, object]]) -> Path:
    summary_path = output_dir / SUMMARY_CSV_NAME
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "band_label",
                "band_label_hitran_style",
                "mode_pair",
                "line_count",
                "reference_temperature_k",
                "wn_min_cm-1",
                "wn_max_cm-1",
                "txt_path",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "band_label": row["band_label"],
                    "band_label_hitran_style": row["band_label_hitran_style"],
                    "mode_pair": row["mode_pair"],
                    "line_count": row["line_count"],
                    "reference_temperature_k": row["reference_temperature_k"],
                    "wn_min_cm-1": row["wn_min_cm-1"],
                    "wn_max_cm-1": row["wn_max_cm-1"],
                    "txt_path": row["txt_path"],
                }
            )
    return summary_path


def main() -> None:
    args = parse_args()
    validate_args(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    def_path = args.data_dir / f"{DATASET_STEM}.def"
    pf_path = args.data_dir / f"{DATASET_STEM}.pf"
    states_path = args.data_dir / f"{DATASET_STEM}.states.bz2"

    def_metadata = parse_exomol_def(def_path)
    column_indices = states_column_indices(def_metadata)
    extra_label_names = signature_label_names(def_metadata)
    pf_temperatures, pf_values = load_partition_function(pf_path)
    partition_function = interpolate_partition_function(
        pf_temperatures,
        pf_values,
        args.reference_temperature_k,
    )

    print("loading states")
    energies, g_totals, signature_ids, signature_lookup, local_label_ids, local_labels = load_state_arrays(
        states_path,
        int(def_metadata["nstates"]),
        column_indices,
        extra_label_names,
    )

    transition_files = transition_files_for_args(args)
    print(f"scanning {len(transition_files)} transition files")

    output_handles: OrderedDict[tuple[object, object], TextIO] = OrderedDict()
    band_rows: dict[tuple[object, object], dict[str, object]] = {}
    stats = defaultdict(int)
    stats["files"] = len(transition_files)

    try:
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
                    require_unit_step=args.require_unit_step,
                ):
                    continue
                stats["matching_states"] += 1

                wavenumber = float(energies[upper_id] - energies[lower_id])
                if not in_selected_window(wavenumber, args.wn_min, args.wn_max):
                    continue
                stats["in_window"] += 1

                intensity = lte_line_intensity_cm_per_molecule(
                    wavenumber=wavenumber,
                    a_coefficient=a_value,
                    g_upper=float(g_totals[upper_id]),
                    lower_energy_cm=float(energies[lower_id]),
                    temperature_k=args.reference_temperature_k,
                    partition_function=partition_function,
                )
                if intensity < args.intensity_threshold:
                    continue

                band_label, band_label_hitran_style = band_label_pair(lower_signature, upper_signature, args.group_by)
                if args.group_by == "hitran-style":
                    key = (band_label, band_label_hitran_style)
                else:
                    key = (lower_signature_id, upper_signature_id)
                handle = output_handles.get(key)
                if handle is None:
                    is_new_band = key not in band_rows
                    if is_new_band:
                        if args.group_by == "hitran-style":
                            out_path = output_path_for_hitran_style_band(
                                args.output_dir,
                                band_label_hitran_style,
                                args.reference_temperature_k,
                            )
                        else:
                            out_path = output_path_for_band(
                                args.output_dir,
                                lower_signature,
                                upper_signature,
                                args.reference_temperature_k,
                            )
                        band_rows[key] = {
                            "band_label": band_label,
                            "band_label_hitran_style": band_label_hitran_style,
                            "mode_pair": f"nu3 {lower_signature.n3}->{upper_signature.n3}",
                            "line_count": 0,
                            "reference_temperature_k": args.reference_temperature_k,
                            "wn_min_cm-1": "" if args.wn_min is None else args.wn_min,
                            "wn_max_cm-1": "" if args.wn_max is None else args.wn_max,
                            "txt_path": str(out_path),
                            "sort_lower_q": lower_signature.n3,
                            "sort_upper_q": upper_signature.n3,
                        }
                        open_mode = "w"
                    else:
                        out_path = Path(str(band_rows[key]["txt_path"]))
                        open_mode = "a"

                    if len(output_handles) >= MAX_OPEN_OUTPUTS:
                        _, old_handle = output_handles.popitem(last=False)
                        old_handle.close()

                    handle = out_path.open(open_mode, encoding="utf-8", newline="")
                    if is_new_band:
                        write_file_preamble(
                            handle,
                            band_label=str(band_rows[key]["band_label"]),
                            band_label_hitran_style=str(band_rows[key]["band_label_hitran_style"]),
                            reference_temperature_k=args.reference_temperature_k,
                            wn_min=args.wn_min,
                            wn_max=args.wn_max,
                        )
                    output_handles[key] = handle
                else:
                    output_handles.move_to_end(key)

                handle.write(
                    f"{upper_id}\t{lower_id}\t{float(energies[upper_id]):.10f}\t{float(energies[lower_id]):.10f}\t"
                    f"{wavenumber:.10f}\t{a_value:.10e}\t{intensity:.10e}\t"
                    f"{local_labels[int(local_label_ids[upper_id])]}\t{local_labels[int(local_label_ids[lower_id])]}\n"
                )
                band_rows[key]["line_count"] = int(band_rows[key]["line_count"]) + 1
                stats["kept"] += 1

                if stats["parsed"] % 1_000_000 == 0:
                    # Progress summary:
                    # - parsed: total transition rows read
                    # - matching_states: transitions whose states pass the pure-nu3 filter
                    # - in_window: matching transitions inside the selected wavenumber range
                    # - kept: transitions exported after all filters, including intensity
                    # - bands: distinct lower-signature -> upper-signature groups exported so far
                    print(
                        f"  parsed {stats['parsed']:,} transitions, "
                        f"matching {stats['matching_states']:,}, "
                        f"window {stats['in_window']:,}, kept {stats['kept']:,}, "
                        f"bands {len(band_rows):,}"
                    )
    finally:
        for handle in output_handles.values():
            handle.close()

    if not band_rows:
        raise RuntimeError("No pure-nu3 band files were produced for the current filters.")

    summary_rows = sorted(
        band_rows.values(),
        key=lambda row: (
            int(row["sort_upper_q"]) - int(row["sort_lower_q"]),
            int(row["sort_lower_q"]),
            str(row["band_label"]),
        ),
    )
    summary_path = write_summary_csv(args.output_dir, summary_rows)

    print()
    print(f"parsed transitions: {stats['parsed']:,}")
    print(f"matching pure-nu3 state pairs: {stats['matching_states']:,}")
    print(f"in selected window: {stats['in_window']:,}")
    print(f"exported transitions: {stats['kept']:,}")
    print(f"exported bands: {len(summary_rows):,}")
    print(f"saved {summary_path}")


if __name__ == "__main__":
    main()
