"""The paper's own numbers, each tied to the output that produces it (claims.csv).

Each row names a results file, a dotted path into it, a format and the text the
paper prints. check() renders every value with its format and fails a row whose
text does not contain the rendered value, so a recomputation that moves a
printed number is caught here rather than by a reader.

Formats: a Python format spec (".3f", "d"); "pct0"/"pct1" for a fraction printed
as a LaTeX percentage; "tex1" for a LaTeX power of ten with one decimal ("tex1@-4" to fix the exponent);
"words" for a small integer written in English.
"""

from __future__ import annotations

import csv
import json
import math

from . import config

CLAIMS = config.ROOT / "claims.csv"
WORDS = ("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")


def value(file: str, path: str):
    node = json.loads((config.ROOT / file).read_text(encoding="utf8"))
    for key in path.split("."):
        node = node[int(key)] if isinstance(node, list) else node[key]
    return node


def render(x, fmt: str) -> str:
    if fmt in ("pct0", "pct1"):
        return f"{100 * x:.{fmt[-1]}f}\\%"
    if fmt.startswith("tex1@"):                      # fixed exponent, e.g. tex1@-4
        exponent = int(fmt[len("tex1@"):])
        return f"{x / 10 ** exponent:.1f}\\times10^{{{exponent}}}"
    if fmt == "tex1":
        exponent = math.floor(math.log10(abs(x)))
        mantissa = x / 10 ** exponent
        if round(mantissa, 1) >= 10:
            mantissa, exponent = mantissa / 10, exponent + 1
        return f"{mantissa:.1f}\\times10^{{{exponent}}}"
    if fmt == "words":
        return WORDS[int(x)]
    text = format(x, fmt)
    return text.replace(",", "{,}") if "," in fmt else text


def rows() -> list[dict]:
    with CLAIMS.open(encoding="utf8") as handle:
        return list(csv.DictReader(handle))


def check() -> list[str]:
    problems = []
    for row in rows():
        try:
            shown = render(value(row["file"], row["path"]), row["format"])
        except (OSError, KeyError, IndexError, ValueError) as error:
            problems.append(f"{row['section']} {row['claim']}: cannot read {row['file']} {row['path']} ({error})")
            continue
        text = row["text"].replace("{,}", "")
        if shown.replace("{,}", "") not in text and shown.lstrip("+") not in text:
            problems.append(f"{row['section']} {row['claim']}: output gives {shown}, paper prints {row['text']}")
    return problems
