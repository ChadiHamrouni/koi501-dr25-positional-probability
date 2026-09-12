"""Run every step in order.

    python run_all.py            # all steps
    python run_all.py 02 04      # selected steps

Each step downloads what it needs, prints the values it computes with the
section of the paper they belong to, and writes its full output to
results/<step>.json. Step 8 draws the figures.
"""

from __future__ import annotations

import runpy
import sys
import time
from pathlib import Path

STEPS = sorted((Path(__file__).parent / "steps").glob("[0-9][0-9]_*.py"))


def main(argv: list[str]) -> int:
    wanted = [s for s in STEPS if not argv or s.name[:2] in argv]
    if not wanted:
        print(f"no step matches {argv}; available: "
              f"{', '.join(s.name[:2] for s in STEPS)}")
        return 2

    started = time.time()
    for step in wanted:
        rule = "=" * 72
        print(f"\n{rule}\nStep {step.name[:2]}  ({step.name})\n{rule}")
        try:
            runpy.run_path(str(step), run_name="__main__")
        except SystemExit as exit_:
            if exit_.code:
                return int(exit_.code)

    print("\n" + "=" * 72)
    print(f"Finished {len(wanted)} steps in {time.time() - started:.0f} s.")
    print("Values: results/*.json    Figures: figures/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
