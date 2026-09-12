"""Step 2 - how far the colour defect extends across the Kepler field.

The r - J test of Section 2.2 uses two KIC columns and nothing else, so it runs
on the whole catalogue. Every DR25 object of interest is cross-matched against
the KIC inside the blend radius and the same test applied.

Reproduces the incidence figures in Section 2.6.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config
from koi501.models import KICNeighbour, KOIRow
from koi501.report import Table, write

KIC_CAT = "vizier:V/133/kic"
XMATCH_COLUMNS = {"id": "koi", "KIC": "kic", "angDist": "sep_as",
                  "RA": "ra", "Dec": "dec", "Jmag": "jmag"}


def main() -> bool:
    kois = list(archives.validate(KOIRow, archives.nasa_tap(
        "select kepoi_name,kepid,ra,dec,koi_disposition,koi_kepmag "
        "from q1_q17_dr25_koi where ra is not null and dec is not null")))
    print(f"  {len(kois)} objects of interest, every one, no sampling")

    pairs = archives.xmatch(
        [{"id": k.kepoi_name, "ra": k.ra, "dec": k.dec} for k in kois],
        KIC_CAT, config.BLEND_RADIUS)
    print(f"  {len(pairs)} object-neighbour pairs inside "
          f"{config.BLEND_RADIUS:.0f} arcsec")

    host_of = {k.kepoi_name: k.kepid for k in kois}
    validated = archives.validate(KICNeighbour, pairs, XMATCH_COLUMNS)

    tested, impossible, suspicious = [], [], []
    for pair in validated:
        if not pair.testable:
            continue
        record = {
            "koi": pair.koi, "kic": pair.kic, "sep_as": pair.sep_as,
            "rmag": pair.rmag, "jmag": pair.jmag,
            "r_minus_J": pair.r_minus_j, "kepmag": pair.kepmag,
            "ra": pair.ra, "dec": pair.dec,
            "is_host": host_of.get(pair.koi) == pair.kic,
        }
        tested.append(record)
        if record["r_minus_J"] < config.COLOUR_IMPOSSIBLE:
            impossible.append(record)
        elif record["r_minus_J"] < config.COLOUR_SUSPICIOUS:
            suspicious.append(record)

    neighbours = [r for r in impossible if not r["is_host"]]
    hosts = [r for r in impossible if r["is_host"]]
    affected = sorted({r["koi"] for r in neighbours})

    disposition = {k.kepoi_name: k.koi_disposition for k in kois}
    counts: dict[str, int] = {}
    for koi in affected:
        counts[disposition.get(koi, "?")] = counts.get(disposition.get(koi, "?"), 0) + 1

    # Histogram of every tested colour, so the cut can be shown to sit in a
    # tail rather than through the bulk of the distribution.
    edges = [round(-6.0 + 0.25 * i, 2) for i in range(int(12 / 0.25) + 1)]
    histogram = [0] * (len(edges) - 1)
    for record in tested:
        colour = record["r_minus_J"]
        if edges[0] <= colour < edges[-1]:
            histogram[int((colour - edges[0]) / 0.25)] += 1

    payload = {
        "n_kois": len(kois),
        "n_pairs_tested": len(tested),
        "colour_histogram": {"edges": edges, "counts": histogram},
        "n_impossible_total": len(impossible),
        "n_impossible_neighbours": len(neighbours),
        "n_impossible_hosts": len(hosts),
        "n_suspicious": len(suspicious),
        "n_kois_affected": len(affected),
        "dispositions": counts,
        "koi501_recovered": any(r["koi"] == config.KOI for r in neighbours),
        "neighbours": sorted(neighbours, key=lambda r: r["r_minus_J"]),
    }
    write("02_colour_pathology", payload)

    table = Table('Impossible catalogue colours, across every Kepler candidate',
                  'Section 2.6')
    table('Kepler candidates tested', len(kois))
    table('neighbouring stars with both r and J', len(tested))
    table('neighbours with r - J < -0.5', len(neighbours))
    table('candidates with at least one such neighbour', len(affected))
    table('  of which FALSE POSITIVE', counts.get("FALSE POSITIVE"))
    table('  of which CONFIRMED', counts.get("CONFIRMED"))
    table('  of which CANDIDATE', counts.get("CANDIDATE"))
    table('  KOI-501.01 is one of them', payload["koi501_recovered"])
    table('host stars with the same defect', len(hosts))
    table.show("02_colour_pathology")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
