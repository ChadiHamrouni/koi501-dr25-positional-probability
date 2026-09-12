"""Step 4 - the archive's positional probability, its inputs, and its consequence.

The `koiapp` table has no API, so its rows come from CSVs exported by hand from
the archive's interactive table view (see data/README.md). This step validates
them, reproduces the correction of Section 2.4, and audits the twelve other
candidates that have an impossible-colour neighbour.
"""

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config
from koi501.models import KOIRow, PositionalProbability
from koi501.report import Table, read, write

KOIAPP = config.DATA / "koiapp"


def load_koiapp() -> list[PositionalProbability]:
    records = []
    for path in sorted(KOIAPP.glob("*.csv")):
        rows = list(csv.DictReader(
            ln for ln in path.read_text().splitlines() if not ln.startswith("#")))
        records.extend(PositionalProbability.model_validate(r) for r in rows)
    return records


def main() -> bool:
    records = load_koiapp()
    by_koi = {r.kepoi_name: r for r in records}
    target = by_koi[config.KOI]
    print(f"  {len(records)} positional-probability records from data/koiapp/")

    photometry = {s["kic"]: s for s in read("01_target_photometry")["stars"]}
    corrected = []
    for tag, kic in (("pp_1hi", 4951867), ("pp_2hi", 4951861)):
        kic_mag = getattr(target, f"{tag}_kepmag")
        depth = getattr(target, f"{tag}_mod_depth")
        gaia_g = photometry[kic]["gaia_G"]
        factor = 10.0 ** (0.4 * (gaia_g - kic_mag))
        corrected.append({
            "kic": kic, "catalogue_kepmag": kic_mag, "gaia_G": gaia_g,
            "magnitude_error": gaia_g - kic_mag, "flux_factor": factor,
            "tabulated_depth_ppm": depth,
            "corrected_depth_ppm": depth * factor,
            "corrected_depth_percent": depth * factor / 1e4,
            "exceeds_rejection_threshold": depth * factor > config.REJECTION_DEPTH_PPM,
        })

    kois = {k.kepoi_name: k for k in archives.validate(KOIRow, archives.nasa_tap(
        "select kepoi_name,kepid,koi_disposition from q1_q17_dr25_koi"))}
    offenders: dict[str, list[int]] = {}
    for row in read("02_colour_pathology")["neighbours"]:
        offenders.setdefault(row["koi"], []).append(row["kic"])

    audit = {"favours_the_bad_star": [], "host_favoured": [],
             "no_usable_value": [], "disfavoured_for_another_reason": []}
    for record in sorted(records, key=lambda r: r.kepoi_name):
        star = record.pp_1hi_starid or ""
        bad_stars = offenders.get(record.kepoi_name, [])
        if not record.usable:
            audit["no_usable_value"].append(record.kepoi_name)
        elif record.favours_host:
            audit["host_favoured"].append(record.kepoi_name)
        elif any(str(k) in star for k in bad_stars):
            audit["favours_the_bad_star"].append(record.kepoi_name)
        else:
            audit["disfavoured_for_another_reason"].append(record.kepoi_name)

    payload = {
        "target": target.model_dump(),
        "corrected_depths": corrected,
        "candidate_audit": audit,
        "records": [r.model_dump() for r in records],
    }
    write("04_positional_probability", payload)

    first, second = corrected
    table = Table("The archive's positional probability, and what it was built from",
                  'Sections 2.3, 2.4, 2.6')
    table("probability given to KOI-501.01's own star", target.pp_host_rel_prob)
    table("archive's quality score for that value", target.pp_host_prob_score)
    table('star the archive favours instead', target.pp_1hi_starid)
    table('  probability given to it', target.pp_1hi_rel_prob)
    table('  Kepler magnitude the archive used', target.pp_1hi_kepmag)
    table('second star, Kepler magnitude used', target.pp_2hi_kepmag)
    table('source of both magnitudes', f"{target.pp_1hi_prob_prov}/{target.pp_2hi_prob_prov}")
    table('eclipse needed on KIC 4951867, as tabulated (ppm)', first["tabulated_depth_ppm"])
    table('eclipse needed on KIC 4951861, as tabulated (ppm)', second["tabulated_depth_ppm"])
    table('eclipse needed on KIC 4951867, corrected (% of its light)', round(first["corrected_depth_percent"]))
    table('eclipse needed on KIC 4951861, corrected (% of its light)', round(second["corrected_depth_percent"]))
    table("KIC 4951861 above the archive's own rejection limit", second["exceeds_rejection_threshold"])
    table('candidates with a bad neighbour, records checked', len(records))
    table('  archive favours the bad neighbour', ", ".join(audit["favours_the_bad_star"]))
    table('  archive favours the host star', len(audit["host_favoured"]))
    table('  no usable value', len(audit["no_usable_value"]))
    table('  host disfavoured for another reason', len(audit["disfavoured_for_another_reason"]))
    table.show("04_positional_probability")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
