"""Step 3 - Gaia settles which of the two KIC columns is at fault.

The colour test says the catalogue contradicts itself. It does not say whether
r is too bright or J too faint. For a normal star Gaia's G sits within a few
tenths of r, so G - r near zero acquits the optical column and a large positive
G - r convicts it.

Reproduces the second half of Section 2.6.
"""

import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config
from koi501.models import GaiaMatch
from koi501.report import Table, read, write, write_table

GAIA_CAT = "vizier:I/355/gaiadr3"


def _convict(stars: dict[int, dict]) -> dict:
    """Match a set of KIC stars to Gaia and ask which column G agrees with."""
    matches = archives.xmatch(
        [{"kic": k, "ra": r["ra"], "dec": r["dec"]} for k, r in stars.items()],
        GAIA_CAT, config.GAIA_XMATCH_AS, selection="best")

    records = []
    for match in archives.validate(GaiaMatch, matches,
                                   {"Gmag": "gmag", "angDist": "sep_as"}):
        star = stars.get(match.kic)
        if star is None or match.gmag is None:
            continue
        records.append({
            "kic": match.kic, "kic_r": star["rmag"], "kic_J": star["jmag"],
            "r_minus_J": star["r_minus_J"], "gaia_G": match.gmag,
            "G_minus_r": match.gmag - star["rmag"],
            "G_minus_J": match.gmag - star["jmag"],
        })
    convicted = [r for r in records if r["G_minus_r"] > config.OPTICAL_COLUMN_WRONG_MAG]
    return {
        "n_stars": len(stars),
        "n_with_gaia": len(records),
        "n_optical_column_wrong": len(convicted),
        "median_G_minus_r": statistics.median(r["G_minus_r"] for r in records),
        "median_G_minus_J": statistics.median(r["G_minus_J"] for r in records),
        "records": sorted(records, key=lambda r: -r["G_minus_r"]),
    }


def main() -> bool:
    step2 = read("02_colour_pathology")
    neighbours = step2["neighbours"]
    print(f"  {len(neighbours)} offending neighbours")

    # The paper quotes the neighbour population. The host population is
    # reported alongside it because an earlier draft mixed the two.
    by_neighbour = _convict({r["kic"]: r for r in neighbours})
    print(f"  neighbours: {by_neighbour['n_with_gaia']} of "
          f"{by_neighbour['n_stars']} unique stars have a Gaia counterpart")

    payload = {"neighbours": by_neighbour}
    write("03_kic_gaia_comparison", payload)

    gaia = {r["kic"]: r for r in by_neighbour["records"]}
    rows = [{**r, "gaia_G": gaia.get(r["kic"], {}).get("gaia_G"),
             "G_minus_r": gaia.get(r["kic"], {}).get("G_minus_r")}
            for r in neighbours]
    write_table(
        "impossible_colour_neighbours",
        f"KIC neighbours within {config.BLEND_RADIUS:g} arcsec of a DR25 object of interest with "
        f"r - J < {config.COLOUR_IMPOSSIBLE} (Section 2.6)",
        [("koi", "-", "DR25 object of interest"),
         ("disposition", "-", "DR25 disposition of that object"),
         ("kic", "-", "KIC identifier of the neighbour"),
         ("ra", "deg", "KIC right ascension, J2000"),
         ("dec", "deg", "KIC declination, J2000"),
         ("sep_as", "arcsec", "separation from the object of interest"),
         ("kepmag", "mag", "KIC Kepler magnitude"),
         ("rmag", "mag", "KIC r magnitude"),
         ("jmag", "mag", "KIC J magnitude"),
         ("r_minus_J", "mag", "KIC r minus KIC J"),
         ("gaia_G", "mag", f"Gaia DR3 G of the best match within {config.GAIA_XMATCH_AS:g} arcsec; blank if none"),
         ("G_minus_r", "mag", "Gaia G minus KIC r; blank if no match")],
        rows)

    table = Table('KIC r against Gaia DR3 G for the impossible-colour neighbours',
                  'Section 2.6')
    table('distinct neighbouring stars', by_neighbour["n_stars"])
    table('  with a Gaia match', by_neighbour["n_with_gaia"])
    table(f'  Gaia more than {config.OPTICAL_COLUMN_WRONG_MAG:g} mag fainter than KIC r', by_neighbour["n_optical_column_wrong"])
    table('median Gaia G minus KIC r (mag)', round(by_neighbour["median_G_minus_r"], 2))
    table.show("03_kic_gaia_comparison")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
