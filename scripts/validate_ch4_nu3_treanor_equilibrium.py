from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable, cast

import numpy as np

from _bootstrap import ensure_repo_root  # pyright: ignore[reportImplicitRelativeImport]

ROOT_DIR = ensure_repo_root()

from research.exomol import (
    collect_nu3_transitions_by_jpair as _collect_nu3_transitions_by_jpair,  # pyright: ignore[reportUnknownVariableType]
)
from research.nonlte import TransitionGroups, collect_nu3_transitions_nonlte


DEFAULT_DATA_DIR = Path("exomol_db/CH4/12C-1H4/MM")
T = 600.0
RELATIVE_ERROR_LIMIT = 1e-10
CollectByJPair = Callable[..., tuple[TransitionGroups, dict[str, object]]]
collect_nu3_transitions_by_jpair = cast(
    CollectByJPair, _collect_nu3_transitions_by_jpair
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate CH4 nu3 Treanor equilibrium identity end-to-end."
    )
    _ = parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help="Path to the ExoMol CH4 MM data directory.",
    )
    return parser.parse_args(argv)


def resolve_data_dir(data_dir: Path) -> Path:
    return data_dir if data_dir.is_absolute() else ROOT_DIR / data_dir


def compare_groups(
    grouped_lte: TransitionGroups,
    grouped_nonlte: TransitionGroups,
) -> tuple[int, float, list[str]]:
    compared_groups = 0
    max_relative_error = 0.0
    issues: list[str] = []

    for key, lte_group in grouped_lte.items():
        nonlte_group = grouped_nonlte.get(key)
        if nonlte_group is None:
            issues.append(f"Missing non-LTE group for key {key}")
            max_relative_error = float("inf")
            continue

        lte_values = np.asarray(lte_group["intensity"], dtype=float)
        nonlte_values = np.asarray(nonlte_group["intensity"], dtype=float)

        if lte_values.shape != nonlte_values.shape:
            issues.append(
                f"Intensity shape mismatch for key {key}: {lte_values.shape} vs {nonlte_values.shape}"
            )
            max_relative_error = float("inf")
            continue

        if lte_values.size == 0:
            compared_groups += 1
            continue

        relative_error = np.abs(nonlte_values - lte_values) / (
            np.abs(lte_values) + 1e-300
        )
        group_max = float(np.max(relative_error))
        max_relative_error = max(max_relative_error, group_max)
        compared_groups += 1

    return compared_groups, max_relative_error, issues


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    data_dir = resolve_data_dir(cast(Path, args.data_dir))

    if not data_dir.exists():
        print(f"Error: data directory does not exist: {data_dir}")
        return 1

    try:
        grouped_lte, _ = collect_nu3_transitions_by_jpair(
            data_dir=data_dir,
            temperature_k=T,
            wn_min=2500.0,
            wn_max=3500.0,
        )
        grouped_nonlte, _ = collect_nu3_transitions_nonlte(
            data_dir=data_dir,
            temperature_k=T,
            vibrational_temperature_k=T,
            wn_min=2500.0,
            wn_max=3500.0,
        )
    except FileNotFoundError as exc:
        print(f"Error reading ExoMol data: {exc}")
        return 1

    compared_groups, max_relative_error, issues = compare_groups(
        grouped_lte, grouped_nonlte
    )
    passed = (
        not issues and compared_groups > 0 and max_relative_error < RELATIVE_ERROR_LIMIT
    )
    status = "PASS" if passed else "FAIL"

    print(f"Groups compared: {compared_groups}")
    print(f"Max relative error: {max_relative_error:.3e}")
    if issues:
        for issue in issues:
            print(f"Issue: {issue}")
    print(status)

    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
