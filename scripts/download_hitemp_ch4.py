"""
Download the official CH4 HITEMP 2020 file from HITRAN/HITEMP.

This script is intentionally similar in spirit to `download_exomol_ch4_mm.py`,
but HITEMP access differs in two important ways:

1. The methane file is published as a single compressed archive:
   `06_HITEMP2020.par.bz2`
2. HITRAN web login is typically required before the file can be downloaded.

Features
--------
- downloads the official CH4 HITEMP archive with retry logic and atomic writes
- optionally signs into the HITRAN website using your account credentials
- optionally extracts the `.bz2` archive into `.par`
- can extract only a wavenumber sub-range from the `.par.bz2` file

Credential handling
-------------------
You can provide HITRAN web-login credentials in either of these ways:

- command line:
  `python scripts/download_hitemp_ch4.py --email you@example.com`
  and you will be prompted for the password
- environment variables:
  `HITRAN_EMAIL` and `HITRAN_PASSWORD`

Optional API key handling
-------------------------
If you also want to send the API key in request headers, you can provide it in
either of these ways:

- edit `DEFAULT_HITRAN_API_KEY` below
- environment variable: `HITRAN_API_KEY`
- command line: `python scripts/download_hitemp_ch4.py --api-key YOUR_KEY`

Examples
--------
- Download the compressed archive only:
  `python scripts/download_hitemp_ch4.py --email you@example.com`

- Download and extract the full `.par` file:
  `python scripts/download_hitemp_ch4.py --email you@example.com --extract-par`

- Download and extract only 2500-3500 cm^-1:
  `python scripts/download_hitemp_ch4.py --email you@example.com --extract-par --numin 2500 --numax 3500`
"""

from __future__ import annotations

import argparse
import bz2
import getpass
import html
import json
import os
import shutil
import sys
import tempfile
import time
from dataclasses import dataclass
from html.parser import HTMLParser
from http.cookiejar import CookieJar
from pathlib import Path
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urljoin, urlparse
from urllib.request import HTTPCookieProcessor, Request, build_opener


ROOT_DIR = Path(__file__).resolve().parents[1]
LOGIN_URL = "https://hitran.org/login/"
DOWNLOAD_URL = "https://hitran.org/files/HITEMP/bzip2format/06_HITEMP2020.par.bz2"
DEFAULT_OUTPUT_DIR = ROOT_DIR / "hitemp_db" / "CH4"
HEADER_TEMPLATE_PATH = ROOT_DIR / "hitran_db" / "CH4_M6_I1.header"
DEFAULT_USER_AGENT = "hapi-hitemp-ch4-downloader/1.0"
DEFAULT_TIMEOUT_SECONDS = 60.0
NU_START = 3
NU_STOP = 15

# Optional local config block.
# Fill these in directly if you prefer editing the script over using CLI flags
# or environment variables. CLI flags and environment variables still override
# these values.
DEFAULT_HITRAN_EMAIL = ""
DEFAULT_HITRAN_PASSWORD = ""
DEFAULT_HITRAN_API_KEY = ""


@dataclass
class LoginForm:
    action: str
    method: str
    inputs: dict[str, str]
    username_field: str | None
    password_field: str | None


class LoginFormParser(HTMLParser):
    """Extract the first form that looks like a login form."""

    def __init__(self) -> None:
        super().__init__()
        self._in_form = False
        self._current_action = ""
        self._current_method = "post"
        self._current_inputs: dict[str, str] = {}
        self._current_username_field: str | None = None
        self._current_password_field: str | None = None
        self.forms: list[LoginForm] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attrs_dict = {key.lower(): value for key, value in attrs}
        if tag.lower() == "form":
            self._in_form = True
            self._current_action = attrs_dict.get("action") or ""
            self._current_method = (attrs_dict.get("method") or "post").lower()
            self._current_inputs = {}
            self._current_username_field = None
            self._current_password_field = None
            return

        if not self._in_form or tag.lower() != "input":
            return

        input_name = attrs_dict.get("name")
        if not input_name:
            return

        input_type = (attrs_dict.get("type") or "text").lower()
        input_value = attrs_dict.get("value") or ""
        self._current_inputs[input_name] = input_value

        if input_type in {"text", "email"} and self._current_username_field is None:
            self._current_username_field = input_name
        elif input_type == "password" and self._current_password_field is None:
            self._current_password_field = input_name

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() != "form" or not self._in_form:
            return

        self.forms.append(
            LoginForm(
                action=self._current_action,
                method=self._current_method,
                inputs=dict(self._current_inputs),
                username_field=self._current_username_field,
                password_field=self._current_password_field,
            )
        )
        self._in_form = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download official CH4 HITEMP 2020 data from HITRAN/HITEMP.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where the downloaded CH4 HITEMP file will be stored.",
    )
    parser.add_argument(
        "--filename",
        default=Path(urlparse(DOWNLOAD_URL).path).name,
        help="Filename to use for the downloaded `.par.bz2` archive.",
    )
    parser.add_argument(
        "--email",
        default=os.environ.get("HITRAN_EMAIL") or DEFAULT_HITRAN_EMAIL,
        help="HITRAN account email. Defaults to HITRAN_EMAIL if set.",
    )
    parser.add_argument(
        "--password",
        default=os.environ.get("HITRAN_PASSWORD") or DEFAULT_HITRAN_PASSWORD,
        help="HITRAN account password. Defaults to HITRAN_PASSWORD if set.",
    )
    parser.add_argument(
        "--api-key",
        default=os.environ.get("HITRAN_API_KEY") or DEFAULT_HITRAN_API_KEY,
        help="Optional HITRAN API key. Defaults to HITRAN_API_KEY if set.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Redownload the archive even if a matching local file already exists.",
    )
    parser.add_argument(
        "--extract-par",
        action="store_true",
        help="After download, extract the `.par.bz2` archive into `.par` text.",
    )
    parser.add_argument(
        "--write-header",
        action="store_true",
        help="Write a HAPI-style `.header` file for the resolved `.par` output path.",
    )
    parser.add_argument(
        "--par-output",
        type=Path,
        default=None,
        help="Explicit output path for the extracted `.par` file.",
    )
    parser.add_argument("--numin", type=float, default=2500, help="Optional minimum wavenumber in cm^-1.")
    parser.add_argument("--numax", type=float, default=3500, help="Optional maximum wavenumber in cm^-1.")
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=DEFAULT_TIMEOUT_SECONDS,
        help="Per-request timeout in seconds.",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.timeout_seconds <= 0:
        raise ValueError("--timeout-seconds must be positive")
    if args.email and not args.password:
        # Prompt later in main(), but allow the common "--email only" workflow.
        pass
    elif args.password and not args.email:
        raise ValueError("--password requires --email (or HITRAN_EMAIL)")
    if (args.numin is None) != (args.numax is None):
        raise ValueError("--numin and --numax must be provided together")
    if args.numin is not None and args.numax is not None and args.numax <= args.numin:
        raise ValueError("--numax must be greater than --numin")


def create_opener() -> tuple[object, CookieJar]:
    cookie_jar = CookieJar()
    opener = build_opener(HTTPCookieProcessor(cookie_jar))
    opener.addheaders = [("User-Agent", DEFAULT_USER_AGENT)]
    return opener, cookie_jar


def build_headers(api_key: str | None = None, extra_headers: dict[str, str] | None = None) -> dict[str, str]:
    headers = {"User-Agent": DEFAULT_USER_AGENT}
    if api_key:
        headers["X-API-Key"] = api_key
    if extra_headers:
        headers.update(extra_headers)
    return headers


def read_text_response(opener: object, url: str, timeout_seconds: float, api_key: str | None = None) -> str:
    request = Request(url, headers=build_headers(api_key=api_key))
    with opener.open(request, timeout=timeout_seconds) as response:
        charset = response.headers.get_content_charset() or "utf-8"
        return response.read().decode(charset, errors="replace")


def select_login_form(html_text: str) -> LoginForm:
    parser = LoginFormParser()
    parser.feed(html_text)
    for form in parser.forms:
        if form.username_field and form.password_field:
            return form
    raise RuntimeError("Could not locate a login form on the HITRAN login page.")


def form_action_url(base_url: str, action: str) -> str:
    if not action:
        return base_url
    return urljoin(base_url, html.unescape(action))


def login_to_hitran(
    opener: object,
    email: str,
    password: str,
    timeout_seconds: float,
    api_key: str | None = None,
) -> None:
    next_path = urlparse(DOWNLOAD_URL).path
    login_page_url = f"{LOGIN_URL}?next={next_path}"
    login_html = read_text_response(opener, login_page_url, timeout_seconds=timeout_seconds, api_key=api_key)
    login_form = select_login_form(login_html)

    payload = dict(login_form.inputs)
    payload[login_form.username_field] = email
    payload[login_form.password_field] = password

    encoded_payload = urlencode(payload).encode("utf-8")
    request = Request(
        form_action_url(login_page_url, login_form.action),
        data=encoded_payload,
        headers=build_headers(
            api_key=api_key,
            extra_headers={
                "Content-Type": "application/x-www-form-urlencoded",
                "Referer": login_page_url,
            },
        ),
    )
    with opener.open(request, timeout=timeout_seconds) as response:
        final_url = response.geturl()
        content_type = response.headers.get("Content-Type", "")
        body = response.read()

    if "/login" in final_url:
        response_text = body.decode("utf-8", errors="replace")
        if "User Login" in response_text or "Password:" in response_text:
            raise RuntimeError("HITRAN login failed. Check your email/password and try again.")

    if "text/html" in content_type and b"User Login" in body:
        raise RuntimeError("HITRAN login failed. The server returned the login page again.")


def remote_size(opener: object, url: str, timeout_seconds: float, api_key: str | None = None) -> int | None:
    request = Request(url, method="HEAD", headers=build_headers(api_key=api_key))
    with opener.open(request, timeout=timeout_seconds) as response:
        if "/login" in response.geturl():
            raise RuntimeError("Download requires login before file metadata can be checked.")
        value = response.headers.get("Content-Length")
    return int(value) if value else None


def ensure_parent_dir(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def looks_like_login_page(response_url: str, content_type: str, first_chunk: bytes) -> bool:
    if "/login" in response_url:
        return True
    if "text/html" not in content_type.lower():
        return False
    head = first_chunk.decode("utf-8", errors="ignore")
    return "User Login" in head or "Email Address:" in head


def download_file(
    opener: object,
    url: str,
    destination: Path,
    overwrite: bool,
    timeout_seconds: float,
    api_key: str | None = None,
    retries: int = 3,
) -> None:
    ensure_parent_dir(destination)

    size_hint = None
    try:
        size_hint = remote_size(opener, url, timeout_seconds=timeout_seconds, api_key=api_key)
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
        tmp_path: Path | None = None
        try:
            with tempfile.NamedTemporaryFile(delete=False, dir=str(destination.parent), suffix=".part") as tmp_file:
                tmp_path = Path(tmp_file.name)
                request = Request(url, headers=build_headers(api_key=api_key))
                with opener.open(request, timeout=timeout_seconds) as response:
                    first_chunk = response.read(8 * 1024)
                    if looks_like_login_page(
                        response_url=response.geturl(),
                        content_type=response.headers.get("Content-Type", ""),
                        first_chunk=first_chunk,
                    ):
                        raise RuntimeError(
                            "The download request was redirected to the HITRAN login page. "
                            "Provide valid credentials with --email/--password or the "
                            "HITRAN_EMAIL/HITRAN_PASSWORD environment variables."
                        )
                    tmp_file.write(first_chunk)
                    shutil.copyfileobj(response, tmp_file, length=8 * 1024 * 1024)
            os.replace(tmp_path, destination)
            print(f"saved {destination.name}")
            return
        except (HTTPError, URLError, TimeoutError, OSError, RuntimeError) as exc:
            if tmp_path and tmp_path.exists():
                tmp_path.unlink()
            if attempt == retries:
                raise RuntimeError(f"failed to download {url}: {exc}") from exc
            wait_seconds = 2 ** (attempt - 1)
            print(f"retry {destination.name} after error: {exc} (attempt {attempt}/{retries})")
            time.sleep(wait_seconds)


def default_par_output_path(archive_path: Path, numin: float | None, numax: float | None) -> Path:
    if archive_path.suffix != ".bz2":
        return archive_path.with_suffix(".par")

    base_name = archive_path.stem
    if numin is None or numax is None:
        return archive_path.with_suffix("")
    return archive_path.with_name(f"{base_name}_{numin:g}-{numax:g}.par")


def count_text_lines(path: Path) -> int:
    line_count = 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_count, _line in enumerate(handle, start=1):
            pass
    return line_count


def iter_filtered_par_lines(archive_path: Path, numin: float | None, numax: float | None) -> Iterable[str]:
    with bz2.open(archive_path, mode="rt", encoding="utf-8", errors="replace", newline="") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if len(line) < NU_STOP:
                raise RuntimeError(f"Unexpected short HITEMP line at {line_number}: {line!r}")
            try:
                nu = float(line[NU_START:NU_STOP])
            except ValueError as exc:
                raise RuntimeError(
                    f"Could not parse wavenumber from HITEMP line {line_number}: {line!r}"
                ) from exc

            if numin is not None and nu < numin:
                continue
            if numax is not None and nu > numax:
                break

            yield line


def extract_par_file(archive_path: Path, output_path: Path, numin: float | None, numax: float | None) -> int:
    ensure_parent_dir(output_path)
    kept_lines = 0
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        for line in iter_filtered_par_lines(archive_path, numin=numin, numax=numax):
            handle.write(line)
            handle.write("\n")
            kept_lines += 1
    print(f"extracted {kept_lines} lines to {output_path}")
    return kept_lines


def write_hapi_header(par_path: Path, line_count: int | None = None) -> Path:
    if not par_path.exists():
        raise FileNotFoundError(f"Cannot write header because the `.par` file does not exist: {par_path}")
    if not HEADER_TEMPLATE_PATH.exists():
        raise FileNotFoundError(f"Missing header template: {HEADER_TEMPLATE_PATH}")

    with HEADER_TEMPLATE_PATH.open("r", encoding="utf-8") as handle:
        header = json.load(handle)

    if line_count is None:
        line_count = count_text_lines(par_path)

    header["table_name"] = par_path.stem
    header["number_of_rows"] = line_count
    header["size_in_bytes"] = par_path.stat().st_size

    header_path = par_path.with_suffix(".header")
    with header_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(header, handle, indent=2)
        handle.write("\n")

    print(f"wrote header to {header_path}")
    return header_path


def main() -> int:
    args = parse_args()
    validate_args(args)

    archive_path = args.output_dir / args.filename
    par_output = args.par_output or default_par_output_path(archive_path, args.numin, args.numax)

    opener, _cookie_jar = create_opener()
    email = args.email
    password = args.password
    api_key = args.api_key

    if email and not password:
        password = getpass.getpass("HITRAN password: ")

    if email and password:
        print(f"login: {LOGIN_URL}")
        login_to_hitran(
            opener,
            email=email,
            password=password,
            timeout_seconds=args.timeout_seconds,
            api_key=api_key,
        )
        print("login successful")
    else:
        print("login: skipped (no credentials provided)")

    print(f"output: {archive_path}")
    download_file(
        opener,
        url=DOWNLOAD_URL,
        destination=archive_path,
        overwrite=args.overwrite,
        timeout_seconds=args.timeout_seconds,
        api_key=api_key,
    )

    if args.extract_par:
        line_count = extract_par_file(archive_path, output_path=par_output, numin=args.numin, numax=args.numax)
        if args.write_header:
            write_hapi_header(par_output, line_count=line_count)
    elif args.write_header:
        write_hapi_header(par_output)

    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
