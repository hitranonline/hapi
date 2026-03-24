from __future__ import annotations

import re
from collections import defaultdict
from typing import Iterable

from .io import write_html_table, write_rows_csv
from .models import SummaryResult


def clean_quanta_label(value: str) -> str:
    text = value.strip()
    return text if text else "000"


def safe_label_fragment(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", clean_quanta_label(value)).strip("_") or "000"


def format_band_label(lower_label: str, upper_label: str) -> str:
    return f"{clean_quanta_label(lower_label)} -> {clean_quanta_label(upper_label)}"


def sort_bands(
    rows: Iterable[dict[str, object]],
    *,
    key: str = "band_label",
    reverse: bool = False,
) -> list[dict[str, object]]:
    return sorted(rows, key=lambda row: str(row.get(key, "")), reverse=reverse)


def merge_bands(
    rows: Iterable[dict[str, object]],
    *,
    key_fields: tuple[str, ...] = ("band_label",),
    sum_fields: tuple[str, ...] = ("line_count", "band_count"),
) -> list[dict[str, object]]:
    grouped: dict[tuple[object, ...], dict[str, object]] = {}

    for row in rows:
        key = tuple(row.get(field) for field in key_fields)
        if key not in grouped:
            grouped[key] = dict(row)
            for field in sum_fields:
                if field in grouped[key]:
                    grouped[key][field] = float(grouped[key][field])
            continue

        target = grouped[key]
        for field in sum_fields:
            if field in row and field in target:
                target[field] = float(target[field]) + float(row[field])

    merged_rows = list(grouped.values())
    return sort_bands(merged_rows)


def build_summary(
    rows: list[dict[str, object]],
    *,
    csv_path=None,
    html_path=None,
    title: str = "Band Summary",
) -> SummaryResult:
    ordered_rows = sort_bands(rows)
    result = SummaryResult(rows=ordered_rows)

    if csv_path is not None:
        result.csv_path = write_rows_csv(csv_path, ordered_rows)
    if html_path is not None:
        result.html_path = write_html_table(html_path, ordered_rows, title=title)

    return result
