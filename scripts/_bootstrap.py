"""Shared bootstrap helpers for direct script execution from scripts/."""

from pathlib import Path
import sys


ROOT_DIR = Path(__file__).resolve().parents[1]


def ensure_repo_root() -> Path:
    root_text = str(ROOT_DIR)
    if root_text not in sys.path:
        sys.path.insert(0, root_text)
    return ROOT_DIR
