"""Shared helpers: small maths, writing results, printing them."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

from .config import RESULTS


def angsep(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    """Angular separation in arcsec."""
    r1, d1, r2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    cos = (math.sin(d1) * math.sin(d2)
           + math.cos(d1) * math.cos(d2) * math.cos(r1 - r2))
    return math.degrees(math.acos(min(1.0, max(-1.0, cos)))) * 3600.0


def flux(mag: float) -> float:
    """Relative flux from a magnitude."""
    return 10.0 ** (-0.4 * mag)


def write(name: str, payload: dict[str, Any]) -> Path:
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / f"{name}.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf8")
    return path


def read(name: str) -> dict[str, Any]:
    return json.loads((RESULTS / f"{name}.json").read_text(encoding="utf8"))


def write_table(name: str, title: str, columns: list[tuple[str, str, str]],
                rows: list[dict[str, Any]]) -> Path:
    """Write a CSV whose '#' header gives each column's unit and meaning.

    ``columns`` is a list of (key, unit, description); units use "-" for
    dimensionless values and text columns.
    """
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / f"{name}.csv"
    with path.open("w", newline="", encoding="utf8") as handle:
        handle.write(f"# {title}\n# {len(rows)} rows\n#\n")
        width = max(len(key) for key, _, _ in columns)
        for key, unit, description in columns:
            handle.write(f"# {key:<{width}}  [{unit}]  {description}\n")
        writer = csv.writer(handle)
        writer.writerow(key for key, _, _ in columns)
        for row in rows:
            writer.writerow("" if row.get(key) is None else row[key]
                            for key, _, _ in columns)
    return path


class Table:
    """Prints the key values a step produces, grouped under a heading."""

    def __init__(self, title: str, section: str) -> None:
        self.title = title
        self.section = section
        self.rows: list[tuple[str, Any]] = []

    def __call__(self, label: str, value: Any) -> None:
        self.rows.append((label, value))

    def show(self, result_file: str) -> None:
        width = max(len(label) for label, _ in self.rows)
        print()
        print(f"  {self.title}")
        print(f"  Paper: {self.section}")
        print()
        for label, value in self.rows:
            if isinstance(value, float):
                value = f"{value:.6g}"
            print(f"    {label:<{width}}   {value}")
        print()
        print(f"  Full output: results/{result_file}.json")
