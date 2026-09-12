"""Step 6 - where the light actually changed, from the sixteen difference images.

The per-quarter centroid offsets are transcribed from the mission's Data
Validation report (see data/README.md). This step averages them, places the two
favoured neighbours against that mean, and gives the result under five error
models so the headline figure is not resting on one choice.

Reproduces Section 3.
"""

import csv
import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config
from koi501.models import KICStar, QuarterCentroid
from koi501.report import Table, write

DV = config.DATA / "dv" / "dvr_quarterly_centroids.csv"
KIC_COLUMNS = "KIC,RAJ2000,DEJ2000,kepmag,rmag,Jmag"
KIC_RENAME = {"KIC": "kic", "RAJ2000": "ra", "DEJ2000": "dec", "Jmag": "jmag"}


def main() -> bool:
    quarters = [QuarterCentroid.model_validate(r)
                for r in csv.DictReader(DV.open())]
    n = len(quarters)
    ra = [q.d_ra for q in quarters]
    dec = [q.d_dec for q in quarters]

    mean = (statistics.mean(ra), statistics.mean(dec))
    scatter = (statistics.stdev(ra), statistics.stdev(dec))
    err_scatter = (scatter[0] / math.sqrt(n), scatter[1] / math.sqrt(n))

    formal = [(q.e_ra, q.e_dec) for q in quarters]
    err_formal = (1 / math.sqrt(sum(1 / e[0] ** 2 for e in formal)),
                  1 / math.sqrt(sum(1 / e[1] ** 2 for e in formal)))

    kic = {s.kic: s for s in archives.validate(
        KICStar, archives.vizier_cone("V/133/kic", KIC_COLUMNS,
                                      config.RA, config.DEC, 30.0), KIC_RENAME)}
    host = kic[config.KEPID]
    cos_dec = math.cos(math.radians(config.DEC))

    def offset_of(star: KICStar) -> tuple[float, float]:
        return ((star.ra - host.ra) * 3600 * cos_dec,
                (star.dec - host.dec) * 3600)

    def sigma(star_offset, err) -> float:
        return math.hypot((mean[0] - star_offset[0]) / err[0],
                          (mean[1] - star_offset[1]) / err[1])

    neighbours = {}
    for kid in config.NEIGHBOURS:
        pos = offset_of(kic[kid])
        neighbours[kid] = {
            "offset_ra_as": pos[0], "offset_dec_as": pos[1],
            "separation_as": math.hypot(*pos),
            "sigma_scatter_error": sigma(pos, err_scatter),
            "sigma_formal_error": sigma(pos, err_formal),
            "sigma_four_independent": sigma(pos, tuple(e * 2 for e in err_scatter)),
            "sigma_single_systematic": sigma(pos, scatter),
        }

    primary = offset_of(kic[config.NEIGHBOURS[0]])
    per_quarter = []
    for q in quarters:
        to_target = math.hypot(q.d_ra, q.d_dec)
        to_star = math.hypot(q.d_ra - primary[0], q.d_dec - primary[1])
        per_quarter.append({
            "quarter": q.quarter, "d_target_as": to_target,
            "d_neighbour_as": to_star,
            "sigma_from_neighbour": math.hypot((q.d_ra - primary[0]) / q.e_ra,
                                               (q.d_dec - primary[1]) / q.e_dec),
            "closer_to_target": to_target < to_star,
        })

    payload = {
        "n_quarters": n, "mean_offset_as": mean, "scatter_as": scatter,
        "error_on_mean_as": err_scatter, "error_formal_as": err_formal,
        "neighbours": neighbours, "per_quarter": per_quarter,
        "n_closer_to_target": sum(q["closer_to_target"] for q in per_quarter),
    }
    write("06_difference_images", payload)

    first = neighbours[config.NEIGHBOURS[0]]
    second = neighbours[config.NEIGHBOURS[1]]
    weakest = min(per_quarter, key=lambda q: q["sigma_from_neighbour"])
    leaning = min(per_quarter, key=lambda q: q["d_neighbour_as"])

    table = Table('Where the light changed, from the 16 difference images',
                  'Section 3')
    table('difference images used', n)
    table('mean offset, east-west (arcsec)', round(mean[0], 2))
    table('mean offset, north-south (arcsec)', round(mean[1], 2))
    table('scatter between images, east-west (arcsec)', round(scatter[0], 2))
    table('scatter between images, north-south (arcsec)', round(scatter[1], 2))
    table('uncertainty of the mean, east-west (arcsec)', round(err_scatter[0], 2))
    table('uncertainty of the mean, north-south (arcsec)', round(err_scatter[1], 2))
    table('KIC 4951867, distance from measured source (sigma)', round(first["sigma_scatter_error"], 1))
    table('KIC 4951861, distance from measured source (sigma)', round(second["sigma_scatter_error"], 1))
    table("  KIC 4951867, using the mission's own error bars", round(first["sigma_formal_error"], 1))
    table('  KIC 4951867, if only 4 images were independent', round(first["sigma_four_independent"], 1))
    table('  KIC 4951867, if all 16 shared one error', round(first["sigma_single_systematic"], 1))
    table('images placing the source nearer KOI-501.01', payload["n_closer_to_target"])
    table('image leaning most toward KIC 4951867 (sigma)', round(leaning["sigma_from_neighbour"], 1))
    table('weakest single image against KIC 4951867 (sigma)', round(weakest["sigma_from_neighbour"], 1))
    table.show("06_difference_images")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
