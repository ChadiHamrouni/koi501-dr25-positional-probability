"""Step 99 - the consistency check: one source for every value.

Four checks, and the run stops if any fails.

1. Measured inputs. Every output that records the measured values it was made
   from (under "measured_inputs") is compared with the step output those values
   come from (koi501/inputs.py). An output made from a value that has since
   changed is stale.

2. Typed numbers. Every fixed value lives in koi501/config.py. This scans the
   programs in steps/, fit/, fpp/ and koi501/ for numeric literals and fails on
   any that is not a unit conversion or a mathematical constant (the ALLOWED set
   below) and is not on a line marked "# literal: <reason>" (plot layout,
   numerical safeguards, display rounding).

3. vespa's inputs. The committed input files of the vespa run are regenerated
   from the step outputs (fpp/vespa/make_inputs.py) and compared value by value.

4. The paper's numbers. Every row of claims.csv is rendered from its output and
   compared with the text the paper prints (koi501/claims.py).
"""

import ast
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import claims, config, inputs
from koi501.report import Table, write

#: Numbers that carry no assumption: identities, halves and doubles, the
#: magnitude definition (-2.5, 0.4, 10), percent, and unit conversions (time,
#: angle, hours to degrees, km to m).
ALLOWED = {0, 1, 2, -1, 0.5, -0.5, 10, 100, -2.5, 2.5, 0.4, -0.4,
           15, 24, 60, 1000, 1440, 3600, 86400}
#: Keyword arguments that only lay out a figure or a printout.
DISPLAY_KEYWORDS = {"lw", "linewidth", "fontsize", "labelsize", "s", "alpha", "zorder",
                    "dpi", "figsize", "xytext", "width", "height", "pad", "markersize",
                    "ms", "mew", "capsize", "elinewidth", "ncol", "handlelength", "bins"}
PROGRAMS = ("steps", "fit", "fpp", "koi501")
#: config.py is where the values belong; archives.py and fits_image.py hold
#: network timeouts and the FITS format, not analysis choices.
EXEMPT = {config.ROOT / "koi501" / name for name in ("config.py", "archives.py", "fits_image.py")}


def structural(tree: ast.AST) -> set[int]:
    """ids of literals that are positions, exponents, rounding digits or layout."""
    skip = set()

    def is_pi(node):
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow):
            node = node.left
        return isinstance(node, ast.Attribute) and node.attr == "pi"

    def mark(node):
        for sub in ast.walk(node):
            if isinstance(sub, ast.Constant):
                skip.add(id(sub))

    for node in ast.walk(tree):
        if isinstance(node, ast.Subscript):
            mark(node.slice)                                   # x[3], x[:, 1:4]
        elif isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow):
            mark(node.right)                                   # r ** 3, m ** (2 / 3)
        elif isinstance(node, ast.BinOp) and any(is_pi(s) for s in (node.left, node.right)):
            for side in (node.left, node.right):               # 4 / 3 * pi, 3 * pi, 4 * pi ** 2
                if not is_pi(side):
                    mark(side)
        elif isinstance(node, ast.Call):
            name = getattr(node.func, "id", getattr(node.func, "attr", ""))
            if name == "round" and len(node.args) > 1:
                mark(node.args[1])                             # round(x, 3)
            if name in ("range", "enumerate", "zip"):
                for arg in node.args:
                    if isinstance(arg, ast.Constant):
                        mark(arg)
            for kw in node.keywords:
                if kw.arg in DISPLAY_KEYWORDS:
                    mark(kw.value)
    return skip


def typed_numbers() -> list[str]:
    problems = []
    for folder in PROGRAMS:
        for path in sorted((config.ROOT / folder).rglob("*.py")):
            if path in EXEMPT:
                continue
            source = path.read_text(encoding="utf8")
            lines = source.splitlines()
            tree = ast.parse(source)
            skip = structural(tree)
            for node in ast.walk(tree):
                if not (isinstance(node, ast.Constant) and type(node.value) in (int, float)):
                    continue
                value = node.value
                if value in ALLOWED or id(node) in skip:
                    continue
                if "# literal:" in lines[node.lineno - 1]:
                    continue
                problems.append(f"{path.relative_to(config.ROOT).as_posix()}:{node.lineno}: "
                                f"{value!r}  | {lines[node.lineno - 1].strip()[:90]}")
    return problems


def vespa_inputs() -> list[str]:
    """vespa's committed input files against those regenerated from the step outputs."""
    spec = importlib.util.spec_from_file_location("make_inputs", config.ROOT / "fpp" / "vespa" / "make_inputs.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return [f"fpp/vespa/{name}: differs from the value regenerated from the step outputs"
            for name, text in module.files().items()
            if module.numbers((module.HERE / name).read_text(encoding="utf8")) != module.numbers(text)]


def main() -> bool:
    stale = inputs.check(stop=False)
    typed = typed_numbers()
    vespa = vespa_inputs()
    printed = claims.check()
    table = Table("One source for every value", "all")
    table("outputs recording their measured inputs", len(inputs.outputs_with_inputs()))
    table("  inputs that disagree with their source", len(stale))
    table("vespa input files that disagree with the step outputs", len(vespa))
    table("numbers typed into the programs outside config.py", len(typed))
    table(f"paper numbers (claims.csv, {len(claims.rows())}) the outputs no longer give", len(printed))
    table.show("99_check")
    for line in stale + vespa + typed + printed:
        print("  " + line)
    write("99_check", {"stale_inputs": stale, "vespa_inputs": vespa, "typed_numbers": typed, "claims": printed})
    return not stale and not typed and not vespa and not printed


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
