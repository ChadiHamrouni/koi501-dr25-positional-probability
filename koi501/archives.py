"""Query layer for the services the analysis touches.

Responses are cached under ``results/.cache`` so a rerun is fast and gives the
same answer. Delete that directory to force a fresh fetch.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import http.client
import io
import json
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Iterable, Iterator

from .config import RESULTS

NASA_TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
VIZIER = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
XMATCH = "https://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync"

HEADERS = {"User-Agent": "koi501-reproduction (hamrounyshady@gmail.com)"}
CACHE = RESULTS / ".cache"


def _cached(key: str, fetch) -> str:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{hashlib.sha256(key.encode()).hexdigest()[:16]}.txt"
    if path.exists():
        return path.read_text(encoding="utf8")
    text = fetch()
    path.write_text(text, encoding="utf8")
    return text


def _open(req: urllib.request.Request, timeout: int, attempts: int = 4) -> bytes:
    """Fetch with retries; public archives drop long requests intermittently."""
    for attempt in range(attempts):
        try:
            return urllib.request.urlopen(req, timeout=timeout).read()
        except (OSError, http.client.HTTPException) as error:
            client_error = isinstance(error, urllib.error.HTTPError) and error.code < 500
            if client_error or attempt == attempts - 1:
                raise
            wait = 15 * 2 ** attempt
            print(f"  {req.full_url.split('?')[0]}: {error}; retrying in {wait} s")
            time.sleep(wait)
    raise AssertionError


def _download(url: str, attempts: int = 8, chunk: int = 1 << 20) -> bytes:
    """A large file read in chunks, resuming with a Range request when the
    connection drops part way (a single read() of a long response can stop short)."""
    data = bytearray()
    for attempt in range(attempts):
        headers = dict(HEADERS, **({"Range": f"bytes={len(data)}-"} if data else {}))
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=headers),
                                        timeout=900) as response:
                total = response.headers.get("Content-Length")
                expected = len(data) + int(total) if total else None
                while block := response.read(chunk):
                    data.extend(block)
            if expected is None or len(data) >= expected:
                return bytes(data)
        except (OSError, http.client.HTTPException) as error:
            print(f"  {url}: {error}; {len(data)} bytes so far, resuming")
        time.sleep(5 * (attempt + 1))
    raise OSError(f"{url}: incomplete after {attempts} attempts")


def _get(url: str, timeout: int = 600) -> str:
    req = urllib.request.Request(url, headers=HEADERS)
    return _open(req, timeout).decode("utf8", "replace")


def _post(url: str, data: bytes, headers: dict, timeout: int = 900) -> str:
    req = urllib.request.Request(url, data, headers={**HEADERS, **headers})
    return _open(req, timeout).decode("utf8", "replace")


def _rows(text: str) -> list[dict]:
    return list(csv.DictReader(io.StringIO(text)))


def nasa_tap(query: str) -> list[dict]:
    """NASA Exoplanet Archive TAP. Serves the DR25 KOI and TCE tables."""
    url = f"{NASA_TAP}?query={urllib.parse.quote(query)}&format=csv"
    return _rows(_cached(f"nasa:{query}", lambda: _get(url)))


def vizier_cone(source: str, columns: str, ra: float, dec: float,
                radius_as: float) -> list[dict]:
    """Cone search of a VizieR catalogue. Used for the KIC, 2MASS, Pan-STARRS."""
    params = [
        ("-source", source), ("-out", columns), ("-out.max", "500"),
        ("-mime", "csv"), ("-c", f"{ra} {dec:+f}"), ("-c.rs", str(radius_as)),
        ("-out.add", "_r"), ("-sort", "_r"),
    ]
    url = f"{VIZIER}?{urllib.parse.urlencode(params)}"
    key = f"vizier:{source}:{columns}:{ra}:{dec}:{radius_as}"
    text = _cached(key, lambda: _get(url))
    lines = [ln for ln in text.split("\n") if ln and not ln.startswith("#")]
    if len(lines) < 4:
        return []
    header = [c.strip() for c in lines[0].split(";")]
    out = []
    for line in lines[3:]:
        cells = line.split(";")
        if len(cells) == len(header):
            out.append({h: cells[i].strip() for i, h in enumerate(header)})
    return out


def vizier_all(source: str, columns: str, max_rows: int = 200_000) -> list[dict]:
    """Every row of a VizieR catalogue."""
    return vizier_by_id(source, columns, None, None, max_rows)


def vizier_by_id(source: str, columns: str, id_column: str | None,
                 value: int | str | None, max_rows: int = 50) -> list[dict]:
    """Rows of a VizieR catalogue selected by an identifier column."""
    params = [("-source", source), ("-out", columns),
              ("-mime", "csv"), ("-out.max", str(max_rows))]
    if id_column is not None:
        params.append((id_column, str(value)))
    url = f"{VIZIER}?{urllib.parse.urlencode(params)}"
    text = _cached(f"vizier-id:{source}:{columns}:{id_column}:{value}", lambda: _get(url))
    lines = [ln for ln in text.split("\n") if ln and not ln.startswith("#")]
    if len(lines) < 4:
        return []
    header = [c.strip() for c in lines[0].split(";")]
    return [{h: cells[i].strip() for i, h in enumerate(header)}
            for cells in (ln.split(";") for ln in lines[3:]) if len(cells) == len(header)]


def xmatch(rows: Iterable[dict], catalogue: str, radius_as: float,
           ra_col: str = "ra", dec_col: str = "dec",
           selection: str = "all") -> list[dict]:
    """Bulk positional cross-match through the CDS X-Match service.

    ``rows`` is uploaded as CSV; every column is echoed back alongside the
    matched catalogue columns and ``angDist``.
    """
    buf = io.StringIO()
    rows = list(rows)
    writer = csv.DictWriter(buf, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    payload = buf.getvalue()

    boundary = "----koi501xmatch"
    body = io.BytesIO()
    fields = {
        "request": "xmatch", "distMaxArcsec": str(radius_as),
        "RESPONSEFORMAT": "csv", "cat2": catalogue,
        "colRA1": ra_col, "colDec1": dec_col, "selection": selection,
    }
    for key, value in fields.items():
        body.write(f'--{boundary}\r\nContent-Disposition: form-data; '
                   f'name="{key}"\r\n\r\n{value}\r\n'.encode())
    body.write(f'--{boundary}\r\nContent-Disposition: form-data; name="cat1"; '
               f'filename="in.csv"\r\nContent-Type: text/csv\r\n\r\n'.encode())
    body.write(payload.encode())
    body.write(f"\r\n--{boundary}--\r\n".encode())

    key = f"xmatch:{catalogue}:{radius_as}:{selection}:{hashlib.sha256(payload.encode()).hexdigest()[:16]}"
    headers = {"Content-Type": f"multipart/form-data; boundary={boundary}"}
    return _rows(_cached(key, lambda: _post(XMATCH, body.getvalue(), headers)))


def mast_table(url: str) -> list[dict]:
    """A gzipped tab-separated bulk catalogue file from MAST."""
    fetch = lambda: gzip.decompress(_open(
        urllib.request.Request(url, headers=HEADERS), 600)).decode("utf8", "replace")
    return list(csv.DictReader(io.StringIO(_cached(f"mast:{url}", fetch)),
                               delimiter="\t"))


EXOFOP_FILE = "https://exofop.ipac.caltech.edu/tess/get_file.php?id={}"


def exofop_file(file_id: int, sha256: str) -> bytes:
    """A file from ExoFOP, cached, and checked against its expected checksum."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"exofop_{file_id}.bin"
    if not path.exists():
        req = urllib.request.Request(EXOFOP_FILE.format(file_id), headers=HEADERS)
        path.write_bytes(_open(req, 600))
    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if digest != sha256:
        path.unlink()
        raise ValueError(f"ExoFOP file {file_id}: checksum {digest} != expected {sha256}")
    return data


def arxiv_member(eprint: str, member: str, sha256: str) -> str:
    """One file from a pinned arXiv source tarball, cached, checked against its checksum."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"arxiv_{eprint}_{member.replace('/', '_')}"
    if not path.exists():
        tarball = _download(f"https://arxiv.org/e-print/{eprint}")
        with tarfile.open(fileobj=io.BytesIO(tarball), mode="r:*") as tar:
            path.write_bytes(tar.extractfile(member).read())
    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    if digest != sha256:
        path.unlink()
        raise ValueError(f"arXiv {eprint} {member}: checksum {digest} != expected {sha256}")
    return data.decode("utf8")


def kepler_fov(kic: int) -> list[dict]:
    """MAST Kepler field-of-view lookup for one KIC star."""
    url = ("https://archive.stsci.edu/kepler/kepler_fov/search.php?"
           f"kic_kepler_id={kic}&outputformat=JSON&action=Search")
    return json.loads(_cached(f"kepler_fov:{kic}", lambda: _get(url, 180)))


def validate(model, rows: Iterable[dict], rename: dict[str, str] | None = None
             ) -> Iterator:
    """Coerce raw catalogue rows through a pydantic model.

    ``rename`` maps catalogue column names onto model field names.
    """
    for row in rows:
        if rename:
            row = {rename.get(k, k): v for k, v in row.items()}
        yield model.model_validate(row)
