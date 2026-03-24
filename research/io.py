from __future__ import annotations

import bz2
import csv
import html
from pathlib import Path
from typing import Iterable

from .models import DatasetPaths


ROOT_DIR = Path(__file__).resolve().parents[1]


def default_paths(root_dir: Path | None = None, artifacts_dir: Path | None = None) -> DatasetPaths:
    root = (root_dir or ROOT_DIR).resolve()
    return DatasetPaths(
        root_dir=root,
        hitran_db_dir=root / "hitran_db",
        exomol_db_dir=root / "exomol_db",
        artifacts_dir=(artifacts_dir or root / "artifacts").resolve(),
    )


def ensure_directory(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def iter_bz2_text_lines(path: Path, chunk_size: int = 8 * 1024 * 1024) -> Iterable[str]:
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


def write_rows_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> Path:
    ensure_directory(path.parent)
    if fieldnames is None:
        keys: list[str] = []
        seen: set[str] = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    seen.add(key)
                    keys.append(key)
        fieldnames = keys

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_html_table(path: Path, rows: list[dict[str, object]], title: str = "Table") -> Path:
    ensure_directory(path.parent)
    headers: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                headers.append(key)

    table_rows: list[str] = []
    for row in rows:
        cells = "".join(f"<td>{html.escape(str(row.get(header, '')))}</td>" for header in headers)
        table_rows.append(f"<tr>{cells}</tr>")

    header_html = "".join(f"<th>{html.escape(header)}</th>" for header in headers)
    body_html = "\n".join(table_rows)
    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
  <style>
    body {{ font-family: sans-serif; margin: 2rem; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #ccc; padding: 0.4rem 0.6rem; text-align: left; }}
    th {{ background: #f5f5f5; }}
  </style>
</head>
<body>
  <h1>{html.escape(title)}</h1>
  <table>
    <thead><tr>{header_html}</tr></thead>
    <tbody>
{body_html}
    </tbody>
  </table>
</body>
</html>
"""
    path.write_text(document, encoding="utf-8")
    return path


def write_markdown(path: Path, text: str) -> Path:
    ensure_directory(path.parent)
    path.write_text(text, encoding="utf-8")
    return path
