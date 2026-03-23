"""
Render a CH4 ExoMol MM absorbance spectrum directly from ExoMol .states/.trans data.

This is separate from the HITRAN/HAPI workflow used elsewhere in the repo.
It computes LTE line intensities from the ExoMol states/trans/partition-function
files and then applies a Voigt profile to produce an absorbance curve.

How the ExoMol inputs are handled
---------------------------------
- The ExoMol `.states.bz2` and `.trans.bz2` files are compressed text files.
  `iter_bz2_text_lines(...)` streams each file, decompresses it in memory with
  `bz2.BZ2Decompressor`, and yields plain-text lines to the caller; it does not
  convert the dataset into standalone `.txt` files on disk.
- The `.def` file is read for metadata such as isotopologue mass, default
  Lorentz half-width (`gamma0`), and the temperature exponent (`n_exponent`)
  by `parse_def_file(...)`.
- The `.pf` file is read to get the partition function at the requested gas
  temperature by `load_partition_function(...)` and
  `interpolate_partition_function(...)`.

How absorbance is computed here
-------------------------------
- This script does not call HAPI to compute coefficients or transmittance.
- Instead it reads ExoMol state energies, total degeneracies, and Einstein-A
  values directly from the ExoMol files via `load_state_arrays(...)` and
  `collect_relevant_transitions(...)`.
- For each kept transition it computes the LTE line intensity with
  `lte_line_intensity_cm_per_molecule(...)`.
- It then builds a Voigt-broadened cross section on the requested wavenumber
  grid with `render_cross_section(...)`.
- Finally it converts cross section to absorbance with Beer-Lambert scaling:
  `absorbance = cross_section * absorber_number_density * path_length`,
  where absorber number density is computed from `number_density_cm3(...)`
  times the CH4 mole fraction.

Where to set the case conditions
--------------------------------
- The main case settings are command-line arguments defined in `parse_args()`.
- Important ones are:
  `--temperature-k`, `--pressure-torr`, `--mole-fraction`, `--path-length-cm`,
  `--wn-min`, `--wn-max`, and `--wn-step`.
- The default values currently correspond to:
  `T=600 K`, `P=3 Torr`, `x=0.008`, `L=100 cm`, `3000-3010 cm^-1`.
- If you want ppm, convert to mole fraction before passing it in:
  `1000 ppm = 1000e-6 = 0.001`.

Notes
-----
- The pressure broadening uses the dataset defaults in the `.def` file unless
  overridden on the command line.
- The script streams transition files and only keeps transitions within the
  requested spectral window plus the line-wing cutoff.
- Large windows may still take a long time because the ExoMol MM files are huge.
"""

from __future__ import annotations

import argparse
import bz2
import csv
import math
import re
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from scipy.special import wofz


ROOT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT_DIR / "exomol_db" / "CH4" / "12C-1H4" / "MM"
OUTPUT_DIR = ROOT_DIR / "exomol_ch4_mm_plots"
DATASET_STEM = "12C-1H4__MM"

T_REF_K = 296.0
BAR_PER_ATM = 1.01325
AMU_TO_KG = 1.66053906660e-27
LIGHT_SPEED_CM_S = 2.99792458e10
LIGHT_SPEED_M_S = 2.99792458e8
BOLTZMANN_J_K = 1.380649e-23
SECOND_RADIATION_CONSTANT_CM_K = 1.438776877
SQRT_LN2 = math.sqrt(math.log(2.0))
SQRT_2LN2 = math.sqrt(2.0 * math.log(2.0))
SQRT_2PI = math.sqrt(2.0 * math.pi)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build CH4 MM absorbance CSV/HTML outputs from ExoMol states/trans files.",
    )
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR, help="Directory containing the ExoMol MM files.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Directory for CSV/HTML outputs.")
    parser.add_argument("--wn-min", type=float, default=3000.0, help="Minimum wavenumber in cm^-1.")
    parser.add_argument("--wn-max", type=float, default=3010.0, help="Maximum wavenumber in cm^-1.")
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
        "--show-top-lines",
        type=int,
        default=10,
        help="Print the strongest kept transitions for quick inspection.",
    )
    return parser.parse_args()


def iter_bz2_text_lines(path: Path, chunk_size: int = 8 * 1024 * 1024):
    decompressor = bz2.BZ2Decompressor()
    pending = b""

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break

            data = decompressor.decompress(chunk)
            if data:
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
        raise FileNotFoundError(
            f"No available transition files overlap {wn_min:g}-{wn_max:g} cm^-1 in {data_dir}"
        )

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

    required_keys = {"mass_da", "nstates", "gamma0", "n_exponent"}
    missing = required_keys - metadata.keys()
    if missing:
        raise RuntimeError(f"Missing required metadata in {def_path}: {sorted(missing)}")
    return metadata


def load_partition_function(pf_path: Path) -> tuple[np.ndarray, np.ndarray]:
    pf_table = np.loadtxt(pf_path)
    if pf_table.ndim != 2 or pf_table.shape[1] < 2:
        raise RuntimeError(f"Unexpected partition function format in {pf_path}")
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

    for line_number, raw_line in enumerate(iter_bz2_text_lines(states_path), start=1):
        parts = raw_line.split(None, 4)
        if len(parts) < 3:
            continue
        state_id = int(parts[0])
        energies[state_id] = float(parts[1])
        g_totals[state_id] = float(parts[2])

        if line_number % 1_000_000 == 0:
            print(f"loaded {line_number:,} states")

    if np.isnan(energies[1:]).any() or np.isnan(g_totals[1:]).any():
        raise RuntimeError(f"State table {states_path} did not fill all expected state IDs")

    return energies, g_totals


def lte_line_intensity_cm_per_molecule(
    wavenumber: float,
    a_coefficient: float,
    g_upper: float,
    lower_energy_cm: float,
    temperature_k: float,
    partition_function: float,
) -> float:
    """
    Convert one ExoMol transition into an LTE line intensity.

    Output meaning
    --------------
    Returns the integrated line intensity `S(T)` in units of `cm/molecule`
    at the requested temperature.

    This is not absorbance yet.
    It is the line strength for one transition before line broadening is applied.

    How it is used later
    --------------------
    - `collect_relevant_transitions(...)` computes one intensity value per kept
      transition with this function.
    - `render_cross_section(...)` multiplies that intensity by a normalized
      Voigt profile to spread the line over the spectral grid and form the
      cross section `sigma(nu)`.
    - `main()` then converts cross section to absorbance with:
      `absorbance = sigma(nu) * absorber_number_density * path_length`.

    Inputs
    ------
    - `wavenumber`: line center in `cm^-1`
    - `a_coefficient`: Einstein-A value from the ExoMol `.trans` row
    - `g_upper`: upper-state total degeneracy from the `.states` file
    - `lower_energy_cm`: lower-state energy in `cm^-1`
    - `temperature_k`: gas temperature in K
    - `partition_function`: `Q(T)` at the same temperature

    In short:
    this function maps ExoMol's `(A, state energies, degeneracy, Q(T))`
    into a temperature-dependent spectroscopic line strength.
    """
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


def doppler_hwhm_cm(wavenumbers: np.ndarray, temperature_k: float, mass_da: float) -> np.ndarray:
    mass_kg = mass_da * AMU_TO_KG
    factor = math.sqrt(2.0 * BOLTZMANN_J_K * temperature_k * math.log(2.0) / (mass_kg * LIGHT_SPEED_M_S**2))
    return wavenumbers * factor


def voigt_profile_cm(
    grid: np.ndarray,
    center_cm: float,
    doppler_hwhm_cm_value: float,
    lorentz_hwhm_cm_value: float,
) -> np.ndarray:
    sigma = max(doppler_hwhm_cm_value / SQRT_2LN2, 1.0e-12)
    z = ((grid - center_cm) + 1j * lorentz_hwhm_cm_value) / (sigma * math.sqrt(2.0))
    return np.real(wofz(z)) / (sigma * SQRT_2PI)


def number_density_cm3(pressure_torr: float, temperature_k: float) -> float:
    pressure_pa = pressure_torr * 133.32236842105263
    density_m3 = pressure_pa / (BOLTZMANN_J_K * temperature_k)
    return density_m3 / 1.0e6


def collect_relevant_transitions(
    transition_files: list[Path],
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
    stats = {
        "files": len(transition_files),
        "parsed": 0,
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

            if stats["parsed"] % 1_000_000 == 0:
                print(
                    f"  parsed {stats['parsed']:,} transitions, "
                    f"window {stats['in_window']:,}, kept {stats['kept']:,}"
                )

    return np.asarray(kept_wavenumbers), np.asarray(kept_intensities), stats


def render_cross_section(
    grid: np.ndarray,
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
    for index, (center, strength, alpha_d) in enumerate(zip(line_centers, line_intensities, doppler_hwhm), start=1):
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

        if index % 10_000 == 0:
            print(f"rendered {index:,} / {len(line_centers):,} kept transitions")

    return cross_section


def save_outputs(
    output_dir: Path,
    stem: str,
    grid: np.ndarray,
    absorbance: np.ndarray,
    metadata_lines: list[str],
) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"{stem}.csv"
    html_path = output_dir / f"{stem}.html"

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["wavenumber_cm-1", "absorbance"])
        for wavenumber, value in zip(grid, absorbance):
            writer.writerow([wavenumber, value])

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=grid, y=absorbance, mode="lines", name="ExoMol MM"))
    fig.update_layout(
        title="<br>".join(metadata_lines),
        xaxis_title="Wavenumber (cm^-1)",
        yaxis_title="Absorbance",
        template="plotly_white",
    )
    fig.update_xaxes(range=[float(grid[0]), float(grid[-1])])
    fig.write_html(html_path, include_plotlyjs="cdn")
    return csv_path, html_path


def print_strongest_lines(line_centers: np.ndarray, line_intensities: np.ndarray, count: int) -> None:
    if count <= 0 or len(line_centers) == 0:
        return
    order = np.argsort(line_intensities)[::-1][:count]
    print("strongest kept lines:")
    for rank, idx in enumerate(order, start=1):
        print(f"  {rank:2d}. nu={line_centers[idx]:10.6f} cm^-1  S={line_intensities[idx]:.3e} cm/molecule")


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


def main() -> None:
    args = parse_args()
    validate_args(args)

    def_path = args.data_dir / f"{DATASET_STEM}.def"
    pf_path = args.data_dir / f"{DATASET_STEM}.pf"
    states_path = args.data_dir / f"{DATASET_STEM}.states.bz2"

    metadata = parse_def_file(def_path)
    gamma0 = args.gamma0 if args.gamma0 is not None else metadata["gamma0"]
    n_exponent = args.n_exponent if args.n_exponent is not None else metadata["n_exponent"]

    pf_temperatures, pf_values = load_partition_function(pf_path)
    partition_function = interpolate_partition_function(pf_temperatures, pf_values, args.temperature_k)

    print("loading states")
    energies, g_totals = load_state_arrays(states_path, int(metadata["nstates"]))

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

    line_centers, line_intensities, stats = collect_relevant_transitions(
        transition_files=transition_files,
        energies=energies,
        g_totals=g_totals,
        wn_min=args.wn_min,
        wn_max=args.wn_max,
        wing=effective_wing,
        temperature_k=args.temperature_k,
        partition_function=partition_function,
        intensity_threshold=args.intensity_threshold,
    )
    if len(line_centers) == 0:
        raise RuntimeError("No transitions survived the current spectral/intensity filters")

    grid = np.arange(args.wn_min, args.wn_max + 0.5 * args.wn_step, args.wn_step, dtype=np.float64)
    cross_section = render_cross_section(
        grid=grid,
        line_centers=line_centers,
        line_intensities=line_intensities,
        mass_da=metadata["mass_da"],
        temperature_k=args.temperature_k,
        pressure_torr=args.pressure_torr,
        gamma0=gamma0,
        n_exponent=n_exponent,
        line_cutoff=args.line_cutoff,
    )

    absorber_density_cm3 = number_density_cm3(args.pressure_torr, args.temperature_k) * args.mole_fraction
    absorbance = cross_section * absorber_density_cm3 * args.path_length_cm

    stem = (
        f"CH4_MM_abs_{int(args.wn_min)}_{int(args.wn_max)}"
        f"_T{args.temperature_k:g}K_P{args.pressure_torr:g}Torr"
        f"_x{args.mole_fraction:g}_L{args.path_length_cm:g}cm"
    )
    metadata_lines = [
        f"CH4 ExoMol MM Absorbance {args.wn_min:g}-{args.wn_max:g} cm^-1",
        (
            f"T={args.temperature_k:g} K, P={args.pressure_torr:g} Torr, "
            f"x={args.mole_fraction:g}, L={args.path_length_cm:g} cm"
        ),
        (
            f"kept {stats['kept']:,} transitions above {args.intensity_threshold:.1e} cm/molecule, "
            f"gamma0={gamma0:g} cm^-1/bar, n={n_exponent:g}"
        ),
    ]
    csv_path, html_path = save_outputs(
        output_dir=args.output_dir,
        stem=stem,
        grid=grid,
        absorbance=absorbance,
        metadata_lines=metadata_lines,
    )

    print_strongest_lines(line_centers, line_intensities, args.show_top_lines)
    print(
        f"parsed {stats['parsed']:,} transitions across {stats['files']} file(s); "
        f"window {stats['in_window']:,}; kept {stats['kept']:,}"
    )
    print(f"partition function Q({args.temperature_k:g} K) = {partition_function:.6g}")
    print(f"saved {csv_path}")
    print(f"saved {html_path}")


if __name__ == "__main__":
    main()
