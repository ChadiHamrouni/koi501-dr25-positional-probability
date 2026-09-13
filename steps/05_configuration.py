"""Step 5 - when the defect actually costs an object its disposition.

A corrupted neighbour only matters if the scene model can plausibly hand it the
transit, which needs it close and catalogued bright. This step applies that
condition, compares the objects it selects using the difference-image centroid,
and shows the conclusion does not depend on where the two thresholds are put.

Reproduces Table 2 and the threshold grid in Section 2.6.
"""

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config
from koi501.models import KOIRow
from koi501.report import Table, flux, read, write, write_table

COLUMNS = ("kepoi_name,kepid,koi_disposition,koi_score,koi_kepmag,"
           "koi_dikco_msky,koi_dikco_msky_err,koi_max_mult_ev,"
           "koi_fpflag_nt,koi_fpflag_ss,koi_fpflag_co,koi_fpflag_ec")


def main() -> bool:
    kois = {k.kepoi_name: k for k in archives.validate(
        KOIRow, archives.nasa_tap(f"select {COLUMNS} from q1_q17_dr25_koi"))}

    pairs = []
    for row in read("02_colour_pathology")["neighbours"]:
        host = kois.get(row["koi"])
        if host is None or host.koi_kepmag is None or row["kepmag"] is None:
            continue
        host_flux, neighbour_flux = flux(host.koi_kepmag), flux(row["kepmag"])
        pairs.append({
            "koi": row["koi"], "kic": row["kic"], "sep_as": row["sep_as"],
            "share": neighbour_flux / (host_flux + neighbour_flux),
            "disposition": host.koi_disposition,
        })

    selected = [p for p in pairs
                if p["sep_as"] <= config.CONFIG_SEP_MAX
                and p["share"] >= config.CONFIG_SHARE_MIN]

    rows = []
    for p in sorted(selected, key=lambda x: x["sep_as"]):
        k = kois[p["koi"]]
        rows.append({**p, "score": k.koi_score, "flags": k.flags,
                      "offset_as": k.koi_dikco_msky,
                      "offset_err_as": k.koi_dikco_msky_err,
                      "offset_sigma": k.offset_sigma,
                      "mes": k.koi_max_mult_ev})

    grid = []
    for sep_max in config.CONFIG_GRID_SEP:
        for share_min in config.CONFIG_GRID_SHARE:
            hits = [p for p in pairs
                    if p["sep_as"] <= sep_max and p["share"] >= share_min]
            grid.append({
                "sep_max": sep_max, "share_min": share_min, "n": len(hits),
                "counts": dict(Counter(p["disposition"] for p in hits)),
                "candidates": sorted({p["koi"] for p in hits
                                      if p["disposition"] == "CANDIDATE"}),
            })

    only_target = [g for g in grid if g["candidates"] == [config.KOI]]
    empty = [g for g in grid if not g["candidates"]]

    payload = {"table": rows, "grid": grid, "all_pairs": pairs,
               "n_grid": len(grid), "n_only_target": len(only_target),
               "n_empty": len(empty),
               "n_other_candidate": len(grid) - len(only_target) - len(empty)}
    write("05_configuration", payload)
    write_table(
        "configuration",
        "Objects with an impossible-colour neighbour within 8 arcsec carrying "
        "at least 25% of the pair's KIC flux (Table 2)",
        [("koi", "-", "DR25 object of interest"),
         ("disposition", "-", "DR25 disposition"),
         ("kic", "-", "KIC identifier of the neighbour"),
         ("sep_as", "arcsec", "neighbour separation"),
         ("share", "-", "neighbour flux / (host + neighbour flux), KIC Kepler magnitudes"),
         ("offset_as", "arcsec", "DR25 koi_dikco_msky, difference-image offset from the KIC position"),
         ("offset_err_as", "arcsec", "DR25 koi_dikco_msky_err"),
         ("offset_sigma", "-", "offset_as / offset_err_as"),
         ("score", "-", "DR25 Robovetter disposition score"),
         ("flags", "-", "DR25 false-positive flags set: nt, ss, co, ec"),
         ("mes", "-", "DR25 maximum multiple event statistic")],
        rows)

    target = next(r for r in rows if r["koi"] == config.KOI)
    others = [r for r in rows if r["koi"] != config.KOI]
    displaced = [r for r in others if r["offset_sigma"] and r["offset_sigma"] > 2]

    table = Table('When the catalogue error can change a verdict',
                  'Section 2.6, Table 2')
    table('candidates with a bad neighbour within 8 arcsec and >= 25% of the light', len(rows))
    table('  of which FALSE POSITIVE', sum(1 for r in others if r["disposition"] == "FALSE POSITIVE"))
    table('  of which CANDIDATE', sum(1 for r in rows if r["disposition"] == "CANDIDATE"))
    table('KOI-501.01, neighbour distance (arcsec)', round(target["sep_as"], 2))
    table("KOI-501.01, neighbour's share of the light (%)", round(target["share"] * 100))
    table('KOI-501.01, difference-image offset (arcsec)', round(target["offset_as"], 2))
    table('KOI-501.01, difference-image offset (sigma)', round(target["offset_sigma"], 1))
    table('KOI-501.01, Robovetter score', target["score"])
    table('KOI-501.01, Robovetter flags', target["flags"])
    table('the other 8: displaced in the difference image', len(displaced))
    table('the other 8: Robovetter score above 0.001', sum(1 for r in others if r["score"] > 0.001))
    table('threshold grid (15 settings): KOI-501.01 only candidate', len(only_target))
    table('threshold grid: no candidate at all', len(empty))
    table('threshold grid: some other candidate', payload["n_other_candidate"])
    table.show("05_configuration")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
