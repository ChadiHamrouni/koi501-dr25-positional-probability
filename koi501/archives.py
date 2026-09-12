"""Query layer for the four services the analysis touches.

Responses are cached under ``results/.cache`` so a rerun is fast and gives the
same answer. Delete that directory to force a fresh fetch.
"""

from __future__ import annotations

import csv
import hashlib
import io
import urllib.parse
import urllib.request
from typing import Iterable, Iterator

from .config import RESULTS

NASA_TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
VIZIER = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
XMATCH = "http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync"

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


def _get(url: str, timeout: int = 600) -> str:
    req = urllib.request.Request(url, headers=HEADERS)
    return urllib.request.urlopen(req, timeout=timeout).read().decode("utf8", "replace")


def _post(url: str, data: bytes, headers: dict, timeout: int = 900) -> str:
    req = urllib.request.Request(url, data, headers={**HEADERS, **headers})
    return urllib.request.urlopen(req, timeout=timeout).read().decode("utf8", "replace")


def _rows(text: str) -> list[dict]:
    return list(csv.DictReader(io.StringIO(text)))


def nasa_tap(query: str) -> list[dict]:
    """NASA Exoplanet Archive TAP. Serves the DR25 KOI and TCE tables."""
    url = f"{NASA_TAP}?query={urllib.parse.quote(query)}&format=csv"
    return _rows(_cached(f"nasa:{query}", lambda: _get(url)))


def gaia_tap(query: str) -> list[dict]:
    """Gaia archive TAP."""
    body = urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
    ).encode()
    return _rows(_cached(f"gaia:{query}", lambda: _post(GAIA_TAP, body, {})))


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

def validate(model, rows: Iterable[dict], rename: dict[str, str] | None = None
             ) -> Iterator:
    """Coerce raw catalogue rows through a pydantic model.

    ``rename`` maps catalogue column names onto model field names.
    """
    for row in rows:
        if rename:
            row = {rename.get(k, k): v for k, v in row.items()}
        yield model.model_validate(row)
