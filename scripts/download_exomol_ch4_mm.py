import argparse
import math
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


BASE_URL = "https://exomol.com/db/CH4/12C-1H4/MM"
MOLECULE_TAG = "12C-1H4__MM"
METADATA_FILES = (
    f"{MOLECULE_TAG}.def",
    f"{MOLECULE_TAG}.pf",
    f"{MOLECULE_TAG}.states.bz2",
)


def build_transition_name(start_cm: int) -> str:
    end_cm = start_cm + 100
    return f"{MOLECULE_TAG}__{start_cm:05d}-{end_cm:05d}.trans.bz2"


def iter_transition_files(wn_min: float, wn_max: float):
    if wn_max <= wn_min:
        raise ValueError("--wn-max must be greater than --wn-min")

    start_chunk = int(math.floor(wn_min / 100.0) * 100)
    end_chunk = int(math.ceil(wn_max / 100.0) * 100)
    for start_cm in range(start_chunk, end_chunk, 100):
        yield build_transition_name(start_cm)


def remote_size(url: str) -> int | None:
    request = Request(url, method="HEAD")
    with urlopen(request) as response:
        value = response.headers.get("Content-Length")
    return int(value) if value else None


def download_file(url: str, destination: Path, overwrite: bool, retries: int = 3) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)

    size_hint = None
    try:
        size_hint = remote_size(url)
    except Exception:
        size_hint = None

    if destination.exists() and not overwrite:
        if size_hint is not None and destination.stat().st_size == size_hint:
            print(f"skip  {destination.name} (already complete)")
            return
        if size_hint is None:
            print(f"skip  {destination.name} (already exists)")
            return

    for attempt in range(1, retries + 1):
        tmp_path = None
        try:
            with tempfile.NamedTemporaryFile(
                delete=False, dir=str(destination.parent), suffix=".part"
            ) as tmp_file:
                tmp_path = Path(tmp_file.name)
                request = Request(url, headers={"User-Agent": "hapi-exomol-downloader/1.0"})
                with urlopen(request) as response:
                    shutil.copyfileobj(response, tmp_file, length=8 * 1024 * 1024)
            os.replace(tmp_path, destination)
            print(f"saved {destination.name}")
            return
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            if tmp_path and tmp_path.exists():
                tmp_path.unlink()
            if attempt == retries:
                raise RuntimeError(f"failed to download {url}: {exc}") from exc
            wait_seconds = 2 ** (attempt - 1)
            print(f"retry {destination.name} after error: {exc} (attempt {attempt}/{retries})")
            time.sleep(wait_seconds)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download CH4 12C-1H4 MM ExoMol files into a local folder."
    )
    parser.add_argument(
        "--output-dir",
        default=str(Path("exomol_db") / "CH4" / "12C-1H4" / "MM"),
        help="Destination directory for downloaded files.",
    )
    parser.add_argument(
        "--wn-min",
        type=float,
        default=3200.0,
        help="Minimum wavenumber in cm^-1 for transition chunks.",
    )
    parser.add_argument(
        "--wn-max",
        type=float,
        default=3600.0,
        help="Maximum wavenumber in cm^-1 for transition chunks.",
    )
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Download only .def, .pf, and .states.bz2 files.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Redownload files even if matching files already exist locally.",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filenames = list(METADATA_FILES)
    if not args.metadata_only:
        filenames.extend(iter_transition_files(args.wn_min, args.wn_max))

    print(f"output: {output_dir}")
    for name in filenames:
        url = f"{BASE_URL}/{name}"
        download_file(url, output_dir / name, overwrite=args.overwrite)

    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
