from __future__ import annotations

from pathlib import Path

from .io import write_markdown
from .models import ComparisonResult


def nearest_lines(
    exomol_rows: list[dict[str, float]],
    hitran_rows: list[dict[str, float]],
    *,
    tolerance: float,
) -> list[dict[str, float]]:
    if tolerance <= 0.0:
        raise ValueError("tolerance must be positive")

    matches: list[dict[str, float]] = []
    for exomol_row in exomol_rows:
        exomol_nu = float(exomol_row["wavenumber"])
        best_row: dict[str, float] | None = None
        best_delta: float | None = None

        for hitran_row in hitran_rows:
            hitran_nu = float(hitran_row["wavenumber"])
            delta = abs(exomol_nu - hitran_nu)
            if delta > tolerance:
                continue
            if best_delta is None or delta < best_delta:
                best_delta = delta
                best_row = hitran_row

        if best_row is not None and best_delta is not None:
            matches.append(
                {
                    "exomol_wavenumber": exomol_nu,
                    "hitran_wavenumber": float(best_row["wavenumber"]),
                    "exomol_intensity": float(exomol_row["intensity"]),
                    "hitran_intensity": float(best_row["intensity"]),
                    "delta_wavenumber": best_delta,
                    "delta_intensity": float(exomol_row["intensity"]) - float(best_row["intensity"]),
                }
            )
    return matches


def build_report(comparison: ComparisonResult, path: Path) -> Path:
    lines = ["# Comparison Report", ""]
    if comparison.band_label:
        lines.append(f"- ExoMol band: `{comparison.band_label}`")
    if comparison.matched_band_label:
        lines.append(f"- HITRAN band: `{comparison.matched_band_label}`")
    for key, value in comparison.metrics.items():
        lines.append(f"- {key}: {value}")
    lines.append("")
    return write_markdown(path, "\n".join(lines))


def compare_bands(*args, **kwargs):
    raise NotImplementedError(
        "compare_bands has not been ported into the first framework pass yet. "
        "Use scripts/report_and_plot_band_intensity_compare.py or "
        "scripts/compare_hitran_exomol_nearest_lines.py as archived references."
    )
