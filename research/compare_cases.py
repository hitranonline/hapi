"""Compare absorbance progressions between two temperature cases.

Reads raw data archives (.npz) produced by
``plot_combined_exomol_i1_absorbance_progressions`` with ``save_raw_data=True``
and generates 4-panel comparison figures with a 2-color scheme.

Usage pattern::

    from research.compare_cases import CaseData, load_raw_case, plot_dual_case_comparison

    case_a = load_raw_case(npz_path, meta_csv_path, label="T = 300 K", color="#1f77b4", temperature_k=300)
    case_b = load_raw_case(npz_path, meta_csv_path, label="T = 2500 K", color="#ff7f0e", temperature_k=2500)
    png = plot_dual_case_comparison(case_a, case_b, progression_label="nu3 1<-0", output_dir=Path("out"))
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .io import ensure_directory

_ABSORBANCE_DELTA_J_VALUES = (-1, 0, 1)


@dataclass(frozen=True)
class CaseData:
    """One temperature case worth of raw data for a single progression."""

    label: str
    color: str
    temperature_k: float
    grid: np.ndarray
    traces: list[dict]


def load_raw_case(
    npz_path: Path,
    meta_csv_path: Path,
    *,
    label: str,
    color: str,
    temperature_k: float,
) -> CaseData:
    """Load raw data from one .npz archive into a CaseData structure.

    Parameters
    ----------
    npz_path
        Path to the ``*_raw.npz`` file.
    meta_csv_path
        Path to the ``*_raw_meta.csv`` sidecar.
    label
        Human-readable label, e.g. ``"T = 300 K"``.
    color
        Hex color for all curves in this case.
    temperature_k
        Temperature in Kelvin.
    """
    archive = np.load(npz_path)
    grid = np.array(archive["wavenumber"])

    traces: list[dict] = []
    with meta_csv_path.open("r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            key = row["array_key"]
            traces.append(
                {
                    "lower_j": int(row["lower_j"]),
                    "upper_j": int(row["upper_j"]),
                    "delta_j": int(row["delta_j"]),
                    "jpair_label": row["jpair_label"],
                    "absorbance": np.array(archive[key]),
                }
            )

    archive.close()

    return CaseData(
        label=label,
        color=color,
        temperature_k=temperature_k,
        grid=grid,
        traces=traces,
    )


def _count_by_delta_j(traces: list[dict], delta_j: int) -> int:
    return sum(1 for t in traces if t["delta_j"] == delta_j)


def _panel_y_limits(curves: list[np.ndarray]) -> tuple[float, float]:
    """Compute y-axis limits with padding for a set of curves (streaming, no vstack)."""
    if not curves:
        return 0.0, 1.0
    lo = float(min(np.min(c) for c in curves))
    hi = float(max(np.max(c) for c in curves))
    if hi <= lo:
        return 0.0, max(1.0, hi)
    pad = (hi - lo) * 0.08
    return max(0.0, lo - pad * 0.5), hi + pad


def plot_dual_case_comparison(
    case_a: CaseData,
    case_b: CaseData,
    *,
    progression_label: str,
    output_dir: Path,
    wn_min: float = 2500.0,
    wn_max: float = 3500.0,
) -> Path:
    """Generate a 4-panel comparison PNG for one progression.

    Panels: P branch (dJ=-1), Q branch (dJ=0), R branch (dJ=+1), All overlaid.
    Each panel shows curves from both cases, with all curves from the same case
    in the same color (2 colors total).

    Returns the path to the generated PNG file.
    """
    if not np.array_equal(case_a.grid, case_b.grid):
        raise ValueError(
            "Grid mismatch between cases. "
            f"Case A has {case_a.grid.size} points, "
            f"Case B has {case_b.grid.size} points. "
            "Both cases must use the same wn_min, wn_max, and wn_step."
        )

    ensure_directory(output_dir)
    figure, axes = plt.subplots(
        4, 1, figsize=(16, 16), sharex=True, constrained_layout=True
    )

    cases = [case_a, case_b]
    figure.suptitle(
        f"{progression_label} absorbance comparison: {case_a.label} vs {case_b.label}"
    )

    wn_mask = (case_a.grid >= wn_min) & (case_a.grid <= wn_max)

    # Panels 1-3: individual branches
    for axis, delta_j in zip(axes[:3], _ABSORBANCE_DELTA_J_VALUES):
        branch_name = {-1: "P", 0: "Q", 1: "R"}[delta_j]
        counts = ", ".join(
            f"{c.label}: {_count_by_delta_j(c.traces, delta_j)} J pairs" for c in cases
        )
        axis.set_title(f"{branch_name} branch (dJ = {delta_j:+d}): {counts}")

        for case in cases:
            branch_traces = [t for t in case.traces if t["delta_j"] == delta_j]
            for trace in branch_traces:
                axis.plot(
                    case.grid,
                    trace["absorbance"],
                    color=case.color,
                    linewidth=0.8,
                    alpha=0.85,
                )

        axis.set_ylabel("Absorbance")
        axis.grid(alpha=0.25)
        all_curves = [
            t["absorbance"][wn_mask]
            for case in cases
            for t in case.traces
            if t["delta_j"] == delta_j
        ]
        if all_curves:
            axis.set_ylim(*_panel_y_limits(all_curves))

    # Panel 4: all overlaid
    overlay_axis = axes[3]
    overlay_axis.set_title("All branches overlaid")

    for case in cases:
        for trace in case.traces:
            overlay_axis.plot(
                case.grid,
                trace["absorbance"],
                color=case.color,
                linewidth=0.8,
                alpha=0.85,
            )
        # One legend entry per case (invisible line for legend only)
        overlay_axis.plot(
            [], [], color=case.color, linewidth=2.0,
            label=f"{case.label} ({len(case.traces)} J pairs)",
        )

    overlay_axis.set_xlabel("Wavenumber (cm^-1)")
    overlay_axis.set_ylabel("Absorbance")
    overlay_axis.legend(loc="upper right")
    overlay_axis.grid(alpha=0.25)

    all_curves = [
        t["absorbance"][wn_mask]
        for case in cases
        for t in case.traces
    ]
    if all_curves:
        overlay_axis.set_ylim(*_panel_y_limits(all_curves))

    slug = progression_label.replace(" ", "_").replace("<-", "_from_")
    png_path = output_dir / f"{slug}_comparison.png"
    figure.savefig(png_path, dpi=150)
    plt.close(figure)
    return png_path
