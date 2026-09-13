"""Step 10 - could another Kepler variable be putting its signal in these pixels?

A bright eclipsing binary elsewhere on the detector can contaminate a target
through the point-response wings, a reflection, or crosstalk, and the target then
inherits its period and epoch. This compares KOI-501.01 against every other DR25
object of interest and every Kirk et al. (2016) eclipsing binary under the two
criteria of Coughlin et al. (2014, Section III.2): a statistical period and epoch
match, and a physical path for the light.

Reproduces the ephemeris-matching paragraph of Section 6.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config
from koi501.models import EclipsingBinary, KeplerFieldPosition, KOIRow
from koi501.report import Table, angsep, write

SIGMA_P_MIN, SIGMA_T_MIN = 3.5, 2.0
DP_MAX = math.erfc(SIGMA_P_MIN / math.sqrt(2.0))
DT_MAX = math.erfc(SIGMA_T_MIN / math.sqrt(2.0))
FOV_CENTRE = (290.66667, 44.5)      # 19h22m40s, +44d30m00s
BKJD_FROM_KIRK = 54833.0            # Kirk et al. tabulate BJD - 2400000


def fractional_mismatch(a: float, b: float) -> float:
    x = (a - b) / a
    return abs(x - round(x))


def epoch_mismatch(period: float, epoch: float, other_epoch: float) -> float:
    x = (epoch - other_epoch) / period
    return abs(x - round(x))


def field_position(kic: int) -> KeplerFieldPosition:
    row = archives.kepler_fov(kic)[0]
    return KeplerFieldPosition.model_validate({
        "ra_sexagesimal": row["RA (J2000)"], "dec_sexagesimal": row["Dec (J2000)"],
        "kepmag": row.get("kepmag"),
        **{f"module_{i}": row.get(f"Module_{i}") for i in range(4)}})


def physical_path(a: KeplerFieldPosition, b: KeplerFieldPosition) -> dict:
    """Coughlin et al. (2014) criterion 2: at least one route for the light."""
    sep = angsep(a.ra, a.dec, b.ra, b.dec)
    brighter = min(m for m in (a.kepmag, b.kepmag, 20.0) if m is not None)
    d_max = 50.0 * math.sqrt(1e6 * 10 ** (-0.4 * brighter) + 1.0)
    anti = angsep(2 * FOV_CENTRE[0] - a.ra, 2 * FOV_CENTRE[1] - a.dec, b.ra, b.dec)
    shared = sorted(a.modules & b.modules)
    return {"separation_as": sep, "prf_reach_as": d_max, "prf_wings": sep < d_max,
            "antipodal_offset_as": anti, "antipodal_reflection": anti < 50.0,
            "shared_modules": shared, "same_module": bool(shared),
            "satisfied": sep < d_max or anti < 50.0 or bool(shared)}


def main() -> bool:
    kois = list(archives.validate(KOIRow, archives.nasa_tap(
        "select kepoi_name,kepid,koi_period,koi_time0bk from q1_q17_dr25_koi "
        "where koi_period is not null and koi_time0bk is not null")))
    target = next(k for k in kois if k.kepoi_name == config.KOI)
    p0, t0 = target.koi_period, target.koi_time0bk

    others = [(k.kepoi_name, k.kepid, k.koi_period, k.koi_time0bk)
              for k in kois if k.kepoi_name != config.KOI]
    ebs = list(archives.validate(
        EclipsingBinary,
        archives.vizier_all("J/AJ/151/68/catalog", "KIC,Per,BJD0"),
        {"KIC": "kic", "Per": "per", "BJD0": "bjd0"}))
    others += [(f"KIC {e.kic}", e.kic, e.per, e.bjd0 - BKJD_FROM_KIRK) for e in ebs]

    compared = []
    for name, kic, period, epoch in others:
        dp = min(fractional_mismatch(p0, period), fractional_mismatch(period, p0))
        dt = min(epoch_mismatch(p0, t0, epoch), epoch_mismatch(period, epoch, t0))
        compared.append({"name": name, "kic": kic, "period_d": period,
                         "dP_prime": dp, "dT_prime": dt})

    statistical = [c for c in compared if c["dP_prime"] < DP_MAX and c["dT_prime"] < DT_MAX]
    home = field_position(config.KEPID)
    for match in statistical:
        match["physical_path"] = physical_path(home, field_position(match["kic"]))
    full = [m for m in statistical if m["physical_path"]["satisfied"]]
    period_only = sorted((c for c in compared if c["dP_prime"] < DP_MAX),
                         key=lambda c: c["dP_prime"])

    payload = {
        "target": {"koi": config.KOI, "period_d": p0, "epoch_bkjd": t0,
                   "modules": sorted(home.modules)},
        "thresholds": {"sigma_P_min": SIGMA_P_MIN, "sigma_T_min": SIGMA_T_MIN,
                       "dP_prime_max": DP_MAX, "dT_prime_max": DT_MAX},
        "n_kois": len(kois) - 1, "n_eclipsing_binaries": len(ebs),
        "n_compared": len(compared),
        "n_period_only": len(period_only),
        "statistical_matches": statistical,
        "n_full_matches": len(full),
        "period_coincident": period_only[:12],
    }
    write("10_ephemeris_match", payload)

    table = Table("Ephemeris matches against every Kepler object and binary",
                  "Section 6, ephemeris matching")
    table("other DR25 objects of interest compared", payload["n_kois"])
    table("Kirk et al. 2016 eclipsing binaries compared", len(ebs))
    table("total compared", len(compared))
    table("period agrees (sigma_P > 3.5)", len(period_only))
    table("period and epoch agree (also sigma_T > 2.0)", len(statistical))
    for m in statistical:
        path = m["physical_path"]
        table(f"  {m['name']}: separation (deg)", round(path["separation_as"] / 3600, 1))
        table(f"  {m['name']}: shares a CCD module", path["same_module"])
        table(f"  {m['name']}: any physical path", path["satisfied"])
    table("matches under both Coughlin criteria", len(full))
    table.show("10_ephemeris_match")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
