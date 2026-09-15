"""Step 16 - how often something instrumental looks like this signal.

Thompson et al. (2018) ran the Kepler pipeline on light curves that contain no
transits by construction - inverted, and with the quarters scrambled - and
published the curated detections (VizieR J/ApJS/235/38, tables 1 and 2) with the
Robovetter's verdict on each. The instrumental share in a cell of period and
signal strength is the number of those detections the Robovetter passes, per
synthetic run, over the real detections in the same cell.

For the cell holding KOI-501.01 the count of synthetic detections passed is
zero, so the share is bounded, not measured. The upper limit on a Poisson count
of zero at confidence c is -ln(1 - c) events. It is given per synthetic run
(the two runs averaged, Thompson's convention) and without averaging (the
conservative reading), over two denominators: every DR25 detection (TCE) in the
cell, and the DR25 candidates in it.

Reproduces the instrumental term of Section 6 and Table 4.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config
from koi501.report import Table, write


def in_cell(period: float, mes: float) -> bool:
    lo, hi = config.INSTRUMENTAL_PERIOD_RANGE_D
    return lo < period <= hi and mes > config.INSTRUMENTAL_MES_MIN


def main() -> bool:
    synthetic = {}
    for run, table in (("inverted", "J/ApJS/235/38/table1"), ("scrambled", "J/ApJS/235/38/table2")):
        rows = archives.vizier_all(table, "Per,MES,Disp")
        cell = [r for r in rows if r["Per"] and r["MES"] and in_cell(float(r["Per"]), float(r["MES"]))]
        synthetic[run] = {"n_rows": len(rows), "n_in_cell": len(cell),
                          "n_passed_in_cell": sum(r["Disp"].strip() == "PC" for r in cell)}

    tces = archives.nasa_tap("select kepid,tce_plnt_num,tce_period,tce_max_mult_ev from q1_q17_dr25_tce")
    mes = {(r["kepid"], r["tce_plnt_num"]): float(r["tce_max_mult_ev"]) for r in tces}
    n_tce = sum(in_cell(float(r["tce_period"]), float(r["tce_max_mult_ev"])) for r in tces)
    kois = archives.nasa_tap("select kepid,kepoi_name,koi_tce_plnt_num,koi_period,koi_pdisposition "
                             "from q1_q17_dr25_koi where koi_pdisposition='CANDIDATE'")
    n_candidates = sum(in_cell(float(r["koi_period"]), mes[(r["kepid"], r["koi_tce_plnt_num"])])
                       for r in kois if (r["kepid"], r["koi_tce_plnt_num"]) in mes)

    passed = sum(s["n_passed_in_cell"] for s in synthetic.values())
    n_noise = sum(s["n_in_cell"] for s in synthetic.values())
    limit = -math.log(1 - config.INSTRUMENTAL_CONFIDENCE) if passed == 0 else None
    bounds = None
    if limit is not None:
        bounds = {f"per_run_over_{name}": limit / config.INSTRUMENTAL_RUNS / n
                  for name, n in (("tces", n_tce), ("candidates", n_candidates))}
        bounds.update({f"unaveraged_over_{name}": limit / n
                       for name, n in (("tces", n_tce), ("candidates", n_candidates))})

    payload = {
        "cell": {"period_range_d": config.INSTRUMENTAL_PERIOD_RANGE_D,
                 "mes_above": config.INSTRUMENTAL_MES_MIN},
        "synthetic": synthetic, "n_noise_in_cell": n_noise, "n_noise_passed": passed,
        "n_dr25_tces_in_cell": n_tce, "n_dr25_candidates_in_cell": n_candidates,
        "confidence": config.INSTRUMENTAL_CONFIDENCE, "poisson_upper_limit_events": limit,
        "share_upper_bounds": bounds,
        "point_estimate": 0.0 if passed == 0 else passed / config.INSTRUMENTAL_RUNS / n_candidates,
    }
    write("16_instrumental_share", payload)

    table = Table("How often the instrument makes a signal like this one", "Section 6; Table 4")
    table("synthetic detections in the cell (inverted + scrambled)",
          f"{synthetic['inverted']['n_in_cell']} + {synthetic['scrambled']['n_in_cell']}")
    table("  passed by the Robovetter", passed)
    table("real DR25 detections in the cell / candidates", f"{n_tce} / {n_candidates}")
    if bounds:
        table(f"{config.INSTRUMENTAL_CONFIDENCE:.0%} upper limit on the count", round(limit, 2))
        for key, value in bounds.items():
            table(f"  share bound, {key.replace('_', ' ')}", f"{value:.2%}")
    table.show("16_instrumental_share")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
