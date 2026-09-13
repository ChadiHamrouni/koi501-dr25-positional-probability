"""Step 9 - the host star's radius, and the planet radius it implies.

The radius comes from the 2MASS Ks magnitude and the Gaia distance, where the
extinction is a tenth of the optical value and the bolometric correction barely
depends on temperature, evaluated at the APOGEE DR17 temperature. The radius
ratio is the MCMC posterior of the transit fit in data/transit_fit/.

Reproduces Section 4 and Table 3.
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config
from koi501.models import (ApogeeStar, CatalogueRadius, GaiaDistance,
                           TwoMassPhotometry, WisePhotometry)
from koi501.report import Table, read, write

TRANSIT_FIT = config.DATA / "transit_fit" / "transit_fit.json"
SEED = 20260904
N = 400_000

MBOL_SUN = 4.74
TEFF_SUN = 5772.0
TEFF_WIDENED_K = 60.0       # APOGEE formal 33 K, widened for systematics
LOGG_WIDENED = 0.035        # APOGEE formal 0.026
G_SI, M_SUN_KG, R_SUN_M = 6.67430e-11, 1.98892e30, 6.957e8
R_EARTH_PER_R_SUN = 109.076


def bc_ks(teff):
    """Ks-band bolometric correction, a quadratic through the solar-metallicity
    table of Casagrande & VandenBerg (2018) over 5000-6500 K."""
    x = (teff - 5772.0) / 1000.0
    return 1.470 - 0.284 * x - 0.070 * x ** 2


def one(model, source, columns, rename, radius=3.0):
    rows = archives.vizier_cone(source, columns, config.RA, config.DEC, radius)
    return next(archives.validate(model, rows[:1], rename))


def by_kic(model, source, columns, rename):
    rows = archives.vizier_by_id(source, columns, "KIC", config.KEPID)
    return next(archives.validate(model, rows[:1], rename))


def spread(x):
    lo, med, hi = np.percentile(x, [15.865, 50.0, 84.135])
    return float(med), float(0.5 * (hi - lo))


def main() -> bool:
    tm = one(TwoMassPhotometry, "II/246/out", "Hmag,e_Hmag,Kmag,e_Kmag,Qflg",
             {"Hmag": "hmag", "e_Hmag": "e_hmag", "Kmag": "kmag",
              "e_Kmag": "e_kmag", "Qflg": "qflg"})
    wise = one(WisePhotometry, "II/328/allwise", "W2mag,e_W2mag",
               {"W2mag": "w2mag", "e_W2mag": "e_w2mag"})
    dist = one(GaiaDistance, "I/352/gedr3dis", "rgeo,b_rgeo,B_rgeo", {})
    apogee = one(ApogeeStar, "III/286/catalog", "Teff,e_Teff,logg,e_logg,SNR",
                 {"Teff": "teff", "e_Teff": "e_teff", "e_logg": "e_logg",
                  "SNR": "snr"})
    dr25 = by_kic(CatalogueRadius, "J/ApJS/229/30/catalog", "Rad,Teff,Mod",
                  {"Rad": "rad", "Teff": "teff"})
    berger = by_kic(CatalogueRadius, "J/AJ/159/280/table2", "Rad,Teff",
                    {"Rad": "rad", "Teff": "teff"})
    host_gaia = read("01_target_photometry")["host_gaia"]
    fit = json.loads(TRANSIT_FIT.read_text(encoding="utf8"))["runs"]["adopted"]

    rng = np.random.default_rng(SEED)
    d_pc = rng.normal(dist.rgeo, 0.5 * (dist.B_rgeo - dist.b_rgeo), N)
    ks = rng.normal(tm.kmag, tm.e_kmag, N)
    a_ks = np.clip(rng.normal(0.05, 0.03, N), 0.0, 0.25)
    teff = rng.normal(apogee.teff, TEFF_WIDENED_K, N)
    logg = rng.normal(apogee.logg, LOGG_WIDENED, N)

    m_ks = ks - 5 * np.log10(d_pc) + 5 - a_ks
    lum = 10 ** (-0.4 * (m_ks + bc_ks(teff) - MBOL_SUN))
    radius = np.sqrt(lum) * (TEFF_SUN / teff) ** 2
    mass = 10 ** logg / 100.0 * (radius * R_SUN_M) ** 2 / G_SI / M_SUN_KG
    density = mass / radius ** 3
    a_ks_rjce = 0.918 * ((tm.hmag - wise.w2mag) - 0.08)

    ror = fit["rp_over_rstar"]
    ror_draw = rng.normal(ror["median"], 0.5 * (ror["plus"] + ror["minus"]), N)
    rp_adopted = ror_draw * radius * R_EARTH_PER_R_SUN

    payload = {
        "inputs": {"twomass": tm.model_dump(), "allwise": wise.model_dump(),
                   "distance_pc": dist.model_dump(), "apogee": apogee.model_dump(),
                   "gaia_parallax_mas": [host_gaia["parallax"], host_gaia["parallax_error"]],
                   "dr25_radius_rsun": dr25.rad, "dr25_teff_k": dr25.teff,
                   "berger2020_radius_rsun": berger.rad,
                   "transit_fit": str(TRANSIT_FIT.relative_to(config.ROOT))},
        "a_ks_rjce": a_ks_rjce,
        "teff_k": spread(teff), "luminosity_lsun": spread(lum),
        "radius_rsun": spread(radius), "mass_msun": spread(mass),
        "density_solar": spread(density),
        "rp_over_rstar": ror,
        "planet_radius_earth": {
            "dr25_star": ror["median"] * dr25.rad * R_EARTH_PER_R_SUN,
            "berger2020_star": ror["median"] * berger.rad * R_EARTH_PER_R_SUN,
            "adopted_star": spread(rp_adopted),
        },
        "seed": SEED, "draws": N,
    }
    write("09_host_star", payload)

    rp = payload["planet_radius_earth"]
    table = Table("The host star's radius, and the planet radius it implies",
                  "Section 4, Table 3")
    table("APOGEE DR17 temperature (K)", f"{apogee.teff:.0f} +- {apogee.e_teff:.0f}")
    table("2MASS Ks (mag)", f"{tm.kmag} +- {tm.e_kmag}")
    table("Gaia geometric distance (pc)", round(dist.rgeo))
    table("DR25 catalogue radius (R_sun)", dr25.rad)
    table("Berger et al. 2020 radius (R_sun)", berger.rad)
    table("radius from Ks and distance (R_sun)", "%.3f +- %.3f" % payload["radius_rsun"])
    table("mass from APOGEE log g (M_sun)", "%.2f +- %.2f" % payload["mass_msun"])
    table("mean density (solar)", "%.3f +- %.3f" % payload["density_solar"])
    table("radius ratio from the transit fit",
          f"{ror['median']:.5f} +{ror['plus']:.5f} -{ror['minus']:.5f}")
    table("planet radius on the DR25 star (R_earth)", round(rp["dr25_star"], 2))
    table("planet radius on the Berger 2020 star (R_earth)", round(rp["berger2020_star"], 2))
    table("planet radius on the adopted star (R_earth)", "%.2f +- %.2f" % rp["adopted_star"])
    table.show("09_host_star")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
