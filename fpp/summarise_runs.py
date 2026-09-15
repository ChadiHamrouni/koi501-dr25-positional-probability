"""The numbers the letter quotes about the published false-positive runs.

vespa        FPP of fpp/vespa/results.txt; mean and standard deviation over the
             100 bootstrap resamplings (results_bootstrap.txt); the population
             size calcfpp generated (calcfpp.log); the radius of the isochrone
             fit vespa made to the adopted star (mist_starmodel_single.h5)
TRICERATOPS  the 20 published runs (results/paper_runs/triceratops_runs.csv):
             mean, median, largest, runs returning exactly zero, runs clearing
             both bars

    python fpp/summarise_runs.py     # needs pandas and PyTables for the .h5

Writes results/fpp_published_runs.json.
"""

from __future__ import annotations

import csv
import json
import re
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from koi501 import config  # noqa: E402

VESPA = HERE / "vespa"
OUT = ROOT / "results" / "fpp_published_runs.json"


def vespa_table(path: Path) -> list[dict]:
    lines = [ln.split() for ln in path.read_text(encoding="utf8").splitlines() if ln.strip()]
    return [dict(zip(lines[0], map(float, row))) for row in lines[1:]]


def main() -> None:
    import pandas as pd

    fpp = vespa_table(VESPA / "results.txt")[0]["FPP"]
    boot = [r["FPP"] for r in vespa_table(VESPA / "results_bootstrap.txt")]
    targets = re.findall(r"(\d+) Transiting planet systems generated \(target (\d+)\)",
                         (VESPA / "calcfpp.log").read_text(encoding="utf8"))
    radius = pd.read_hdf(VESPA / "mist_starmodel_single.h5", "samples")["radius_0_0"]

    with (ROOT / "results" / "paper_runs" / "triceratops_runs.csv").open(encoding="utf8") as handle:
        runs = list(csv.DictReader(handle))
    fpps = [float(r["FPP"]) for r in runs]
    nfpps = [float(r["NFPP"]) for r in runs]

    result = {
        "vespa": {"fpp": fpp, "bootstrap_n": len(boot), "bootstrap_mean": statistics.mean(boot),
                  "bootstrap_sd": statistics.stdev(boot),
                  "population_generated": int(targets[-1][0]), "population_target": int(targets[-1][1]),
                  "isochrone_radius_rsun": [float(radius.median()), float(radius.std())]},
        "triceratops_published": {
            "n_runs": len(runs), "fpp_mean": statistics.mean(fpps), "fpp_median": statistics.median(fpps),
            "fpp_max": max(fpps), "n_fpp_exactly_zero": sum(f == 0 for f in fpps),
            "nfpp_mean": statistics.mean(nfpps), "nfpp_max": max(nfpps),
            "n_pass_both_bars": sum(f < config.TRICERATOPS_FPP_BAR and n < config.TRICERATOPS_NFPP_BAR
                                    for f, n in zip(fpps, nfpps))},
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf8")
    v, t = result["vespa"], result["triceratops_published"]
    print(f"  vespa FPP {v['fpp']:.2e}; bootstrap sd {v['bootstrap_sd']:.2e} over {v['bootstrap_n']}; "
          f"population target {v['population_target']}; radius {v['isochrone_radius_rsun'][0]:.3f} "
          f"+/- {v['isochrone_radius_rsun'][1]:.3f}")
    print(f"  TRICERATOPS published: mean {t['fpp_mean']:.2e}, median {t['fpp_median']:.1e}, "
          f"max {t['fpp_max']:.4f}, {t['n_fpp_exactly_zero']} zeros, {t['n_pass_both_bars']} pass")


if __name__ == "__main__":
    main()
