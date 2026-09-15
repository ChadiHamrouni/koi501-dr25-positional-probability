"""Step 15 - earlier verdicts on KOI-501.01, and the population the letter compares it with.

Published values (read, not reproduced): the scores of Armstrong et al. (2021)
and ExoMiner (Valizadegan et al. 2022), from pinned arXiv source tables checked
against their SHA-256; ExoMiner V1.2 (Valizadegan et al. 2023) and the Berger et
al. (2018) evolutionary state, from VizieR.

Calculated here:
  implied prior   ExoMiner combines its score s with a prior pi by Bayes' rule,
                  s' = s pi / (s pi + (1 - s)(1 - pi)); inverting it for the two
                  published scores gives pi = s'(1 - s) / (s(1 - s') + s'(1 - s))
  population      DR25 objects of interest the NASA Exoplanet Archive lists as
                  CONFIRMED, matched to a named Kepler planet in pscomppars whose
                  mass provenance is neither Mass nor Msini, and how many of their
                  hosts Berger et al. (2018) classify as subgiants and giants.
                  The archive changes, so this count carries the date of the query.

Reproduces the verdicts paragraph of Section 1 and the host-star paragraph of Section 4.
"""

import csv
import io
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config
from koi501.report import Table, write


def arxiv_rows(name: str) -> list[dict]:
    eprint, member, sha256 = config.ARXIV_TABLES[name]
    return list(csv.DictReader(io.StringIO(archives.arxiv_member(eprint, member, sha256))))


def main() -> bool:
    armstrong = next(r for r in arxiv_rows("armstrong2021") if r["kepoi_name"] == config.KOI)
    exominer = next(r for r in arxiv_rows("valizadegan2022")
                    if int(float(r["target_id"])) == config.KEPID
                    and int(float(r["tce_plnt_num"])) == config.TCE_PLANET_NUMBER)
    s, s_prior = float(exominer["exominer"]), float(exominer["exominer_prior"])
    implied_prior = s_prior * (1 - s) / (s * (1 - s_prior) + s_prior * (1 - s))
    v12 = archives.vizier_by_id("J/AJ/166/28/table11", "KIC,Score1,Score2,Prob,RUWE",
                                "KIC", config.KEPID)[0]

    berger = {int(r["KIC"]): int(r["Evol"])
              for r in archives.vizier_all("J/ApJ/866/99/table1", "KIC,Evol")
              if r["KIC"] and r["Evol"] != ""}
    dr25 = {r["kepoi_name"]: int(r["kepid"])
            for r in archives.nasa_tap("select kepid,kepoi_name from q1_q17_dr25_koi")}
    cumulative = {r["kepoi_name"]: r for r in archives.nasa_tap(
        "select kepoi_name,kepler_name,koi_disposition from cumulative")}
    masses = {r["pl_name"]: r["pl_bmassprov"] for r in archives.nasa_tap(
        "select pl_name,pl_bmassprov from pscomppars where disc_facility like '%Kepler%'")}
    validated = [kepid for koi, kepid in dr25.items()
                 if cumulative.get(koi, {}).get("koi_disposition") == "CONFIRMED"
                 and cumulative[koi]["kepler_name"] in masses
                 and masses[cumulative[koi]["kepler_name"]] not in config.MEASURED_MASS_PROVENANCE]
    states = [berger.get(k) for k in validated]
    population = {
        "query_date": str(date.today()),
        "n_validated": len(validated),
        "n_with_berger_class": sum(x is not None for x in states),
        "n_subgiant_host": states.count(config.BERGER_SUBGIANT),
        "n_giant_host": states.count(config.BERGER_GIANT),
        "n_without_berger_class": states.count(None),
    }

    payload = {
        "armstrong2021": {"source": "arXiv:{} {}".format(*config.ARXIV_TABLES["armstrong2021"][:2]),
                          "gpc": float(armstrong["PP_GP_final"]), "rfc": float(armstrong["PP_RFC_final"]),
                          "vespa_fpp": float(armstrong["fpp_prob"]),
                          "evolutionary_state_flag": int(float(armstrong["Evolutionary_State_Flag"])),
                          "gaia_flag": int(float(armstrong["gaia_flag"])),
                          "pp_host_rel_prob": float(armstrong["pp_host_rel_prob"])},
        "valizadegan2022": {"source": "arXiv:{} {}".format(*config.ARXIV_TABLES["valizadegan2022"][:2]),
                            "exominer": s, "exominer_with_prior": s_prior,
                            "implied_prior": implied_prior},
        "valizadegan2023": {"vizier": "J/AJ/166/28/table11", "exominer_v12": float(v12["Score1"]),
                            "positional_probability_used": float(v12["Prob"])},
        "berger2018_host_evolutionary_state": berger.get(config.KEPID),
        "population": population,
    }
    write("15_literature", payload)

    table = Table("Earlier verdicts, and the population of validated Kepler planets",
                  "Sections 1 and 4")
    table("Armstrong et al. 2021, GP classifier score", round(payload["armstrong2021"]["gpc"], 4))
    table("ExoMiner score / with priors", f"{s:.3f} / {s_prior:.3f}")
    table("  prior those two scores imply", f"{implied_prior:.2e}")
    table("ExoMiner V1.2 score", round(payload["valizadegan2023"]["exominer_v12"], 3))
    table("Berger et al. 2018 state of the host (1 = subgiant)", payload["berger2018_host_evolutionary_state"])
    table(f"Kepler planets confirmed without a measured mass ({population['query_date']})",
          population["n_validated"])
    table("  hosts Berger et al. classify", population["n_with_berger_class"])
    table("  of which subgiants / giants", f"{population['n_subgiant_host']} / {population['n_giant_host']}")
    table.show("15_literature")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
