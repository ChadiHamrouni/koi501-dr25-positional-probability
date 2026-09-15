"""Step 0 - the DR25 measurements of KOI-501.01 that later steps and programs use.

The DR25 threshold-crossing event (TCE) row carries the transit fit the
pipeline's own vetting and the false-positive codes use: period, epoch,
duration, depth, radius ratio, detection statistic, weak-secondary search and
odd-even statistic. The DR25 object-of-interest (KOI) row carries the period to
more digits, the Kepler magnitude and the difference-image offsets on each axis.
Both rows are written whole, so every later use reads them from here
(koi501/inputs.py) rather than retyping them.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config
from koi501.report import Table, write

TCE_COLUMNS = ("kepid,tce_plnt_num,tce_period,tce_time0bk,tce_duration,tce_duration_err,"
               "tce_depth,tce_depth_err,tce_ror,tce_max_mult_ev,tce_num_transits,"
               "wst_depth,wst_depth_err,wst_robstat,tce_bin_oedp_stat,"
               "tce_dikco_msky,tce_dikco_msky_err,tce_dicco_msky,tce_dicco_msky_err,"
               "tce_fwm_srao,tce_fwm_srao_err,tce_fwm_sdeco,tce_fwm_sdeco_err,tce_fwm_stat")
KOI_COLUMNS = ("kepoi_name,kepid,koi_period,koi_time0bk,koi_kepmag,koi_depth,"
               "koi_dikco_mra,koi_dikco_mra_err,koi_dikco_mdec,koi_dikco_mdec_err,"
               "koi_dikco_msky,koi_dikco_msky_err,koi_dicco_msky,koi_dicco_msky_err,"
               "koi_fwm_stat_sig")


def numbers(row: dict) -> dict:
    out = {}
    for key, value in row.items():
        try:
            out[key] = int(value) if value.lstrip("-").isdigit() else float(value)
        except (ValueError, AttributeError):
            out[key] = value
    return out


def main() -> bool:
    tce = numbers(archives.nasa_tap(
        f"select {TCE_COLUMNS} from q1_q17_dr25_tce "
        f"where kepid={config.KEPID} and tce_plnt_num={config.TCE_PLANET_NUMBER}")[0])
    koi = numbers(archives.nasa_tap(
        f"select {KOI_COLUMNS} from q1_q17_dr25_koi where kepoi_name='{config.KOI}'")[0])
    write("00_dr25_target", {"tce": tce, "koi": koi})

    table = Table("The DR25 measurements of KOI-501.01", "Table 3; inputs to Sections 6 and 7")
    table("period (d), detection table / candidate table",
          f"{tce['tce_period']} / {koi['koi_period']}")
    table("epoch (BKJD), detection table", tce["tce_time0bk"])
    table("duration (h)", f"{tce['tce_duration']} +- {tce['tce_duration_err']}")
    table("depth (ppm)", f"{tce['tce_depth']} +- {tce['tce_depth_err']}")
    table("multiple event statistic / transits", f"{tce['tce_max_mult_ev']} / {tce['tce_num_transits']}")
    table("weak secondary (ppm) / robust statistic",
          f"{tce['wst_depth']} +- {tce['wst_depth_err']} / {tce['wst_robstat']}")
    table("odd-even statistic", tce["tce_bin_oedp_stat"])
    table.show("00_dr25_target")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
