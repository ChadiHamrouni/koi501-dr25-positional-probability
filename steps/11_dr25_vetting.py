"""Step 11 - the mission's own vetting statistics for KOI-501.01.

Every row of Table 4 marked DR25 comes from the DR25 threshold-crossing event
table or the DR25 KOI table, queried here rather than copied from the Data
Validation report.

Reproduces the DR25 rows of Table 4, Section 6.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config
from koi501.models import KOIRow, TCERow
from koi501.report import Table, write

TCE_COLUMNS = ("kepid,tce_plnt_num,tce_max_mult_ev,tce_num_transits,tce_depth,"
               "tce_depth_err,boot_fap,tce_bin_oedp_stat,wst_depth,wst_depth_err,"
               "wst_robstat,tce_maxmesd,tce_maxmes,tce_cap_stat,tce_hap_stat,"
               "tce_rb_tcount0,tce_rb_tcount1,tce_rb_tcount2,tce_rb_tcount3,"
               "tce_rb_tcount4,tce_fwm_stat")
KOI_COLUMNS = ("kepoi_name,kepid,koi_disposition,koi_score,koi_max_mult_ev,"
               "koi_fpflag_nt,koi_fpflag_ss,koi_fpflag_co,koi_fpflag_ec")


def main() -> bool:
    tce = next(archives.validate(TCERow, archives.nasa_tap(
        f"select {TCE_COLUMNS} from q1_q17_dr25_tce "
        f"where kepid={config.KEPID} and tce_plnt_num=1")))
    koi = next(archives.validate(KOIRow, archives.nasa_tap(
        f"select {KOI_COLUMNS} from q1_q17_dr25_koi "
        f"where kepoi_name='{config.KOI}'")))

    rolling = [tce.tce_rb_tcount0, tce.tce_rb_tcount1, tce.tce_rb_tcount2,
               tce.tce_rb_tcount3, tce.tce_rb_tcount4]
    payload = {
        "tce": tce.model_dump(), "koi": koi.model_dump(),
        "secondary_3sigma_limit_ppm": tce.wst_depth + 3 * tce.wst_depth_err,
        "rolling_band_severity_ge1": sum(rolling[1:]),
        "rolling_band_transits": sum(rolling),
        "core_over_halo": tce.tce_cap_stat / tce.tce_hap_stat,
        "robovetter_flags": koi.flags,
    }
    write("11_dr25_vetting", payload)

    table = Table("The mission's own vetting statistics", "Section 6, Table 4 (DR25 rows)")
    table("multiple event statistic", round(tce.tce_max_mult_ev, 2))
    table("transits", tce.tce_num_transits)
    table("bootstrap false-alarm probability", f"{tce.boot_fap:.2e}")
    table("rolling-band contamination, severity >= 1",
          f"{payload['rolling_band_severity_ge1']} of {payload['rolling_band_transits']}")
    table("weak secondary depth (ppm)", f"{tce.wst_depth} +- {tce.wst_depth_err}")
    table("weak secondary significance (sigma)", tce.wst_robstat)
    table("  at phase (days after transit)", tce.tce_maxmesd)
    table("secondary 3 sigma upper limit (ppm)", round(payload["secondary_3sigma_limit_ppm"], 1))
    table("odd-even depth statistic", tce.tce_bin_oedp_stat)
    table("optical ghost, core / halo statistic", f"{tce.tce_cap_stat} / {tce.tce_hap_stat}")
    table("Robovetter disposition score", koi.koi_score)
    table("Robovetter false-positive flags", koi.flags)
    table.show("11_dr25_vetting")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
