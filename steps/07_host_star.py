"""Step 7 - the host star, the planet radius it implies, and the planet's orbit.

The radius comes from the 2MASS Ks magnitude and the Gaia distance, where the
extinction is a tenth of the optical value and the bolometric correction barely
depends on temperature, evaluated at the APOGEE DR17 temperature. The mass
follows from the APOGEE surface gravity and that radius. The radius ratio is the
MCMC posterior of fit/transit_fit.py in data/transit_fit/. The orbital distance
follows from Kepler's third law on the adopted mass alone, and the equilibrium
temperature and insolation from that distance.

Reproduces Section 4 and Table 3.
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config, inputs
from koi501.models import (ApogeeStar, CatalogueRadius, GaiaDistance,
                           TwoMassPhotometry, WisePhotometry)
from koi501.report import Table, read, write

TRANSIT_FIT = config.DATA / "transit_fit" / "transit_fit.json"


def bc_ks(teff):
    """Ks-band bolometric correction, quadratic in (Teff - Teff_sun) / 1000 K."""
    c0, c1, c2 = config.BC_KS_COEFFS
    x = (teff - config.TEFF_SUN_K) / 1000.0
    return c0 + c1 * x + c2 * x ** 2


def one(model, source, columns, rename):
    rows = archives.vizier_cone(source, columns, config.RA, config.DEC, config.POINT_MATCH_AS)
    return next(archives.validate(model, rows[:1], rename))


def by_kic(model, source, columns, rename):
    rows = archives.vizier_by_id(source, columns, "KIC", config.KEPID)
    return next(archives.validate(model, rows[:1], rename))


def spread(x):
    lo, med, hi = np.percentile(x, config.ONE_SIGMA_PERCENTILES)
    return float(med), float(0.5 * (hi - lo))


def main() -> bool:
    tm = one(TwoMassPhotometry, "II/246/out", "Jmag,e_Jmag,Hmag,e_Hmag,Kmag,e_Kmag,Qflg",
             {"Jmag": "jmag", "e_Jmag": "e_jmag", "Hmag": "hmag", "e_Hmag": "e_hmag",
              "Kmag": "kmag", "e_Kmag": "e_kmag", "Qflg": "qflg"})
    wise = one(WisePhotometry, "II/328/allwise", "W2mag,e_W2mag",
               {"W2mag": "w2mag", "e_W2mag": "e_w2mag"})
    dist = one(GaiaDistance, "I/352/gedr3dis", "rgeo,b_rgeo,B_rgeo", {})
    apogee = one(ApogeeStar, "III/286/catalog",
                 "Teff,e_Teff,logg,e_logg,[Fe/H],e_[Fe/H],SNR",
                 {"Teff": "teff", "e_Teff": "e_teff", "e_logg": "e_logg",
                  "[Fe/H]": "feh", "e_[Fe/H]": "e_feh", "SNR": "snr"})
    dr25 = by_kic(CatalogueRadius, "J/ApJS/229/30/catalog", "Rad,Teff,Mod",
                  {"Rad": "rad", "Teff": "teff"})
    berger = by_kic(CatalogueRadius, "J/AJ/159/280/table2", "Rad,Teff",
                    {"Rad": "rad", "Teff": "teff"})
    host_gaia = read("01_target_photometry")["host_gaia"]
    fit = json.loads(TRANSIT_FIT.read_text(encoding="utf8"))["runs"]["adopted"]
    period_d = inputs.get("transit.period_koi_d")

    n = config.HOST_DRAWS
    rng = np.random.default_rng(config.HOST_SEED)
    d_pc = rng.normal(dist.rgeo, 0.5 * (dist.B_rgeo - dist.b_rgeo), n)
    ks = rng.normal(tm.kmag, tm.e_kmag, n)
    ak_mean, ak_sd, ak_max = config.A_KS_PRIOR
    a_ks = np.clip(rng.normal(ak_mean, ak_sd, n), 0.0, ak_max)
    teff = rng.normal(apogee.teff, config.TEFF_WIDENED_K, n)
    logg = rng.normal(apogee.logg, config.LOGG_WIDENED, n)

    m_ks = ks - 5 * np.log10(d_pc) + 5 - a_ks  # literal: distance modulus
    lum = 10 ** (-0.4 * (m_ks + bc_ks(teff) - config.MBOL_SUN))
    radius = np.sqrt(lum) * (config.TEFF_SUN_K / teff) ** 2
    mass = 10 ** logg / 100.0 * (radius * config.R_SUN_M) ** 2 / config.G_SI / config.M_SUN_KG
    density = mass / radius ** 3
    rjce_scale, rjce_offset = config.RJCE
    a_ks_rjce = rjce_scale * ((tm.hmag - wise.w2mag) - rjce_offset)

    ror = fit["rp_over_rstar"]
    ror_draw = rng.normal(ror["median"], 0.5 * (ror["plus"] + ror["minus"]), n)
    rp_adopted = ror_draw * radius * config.R_EARTH_PER_R_SUN

    # the orbit on the adopted solution alone: Kepler's third law on its mass
    a_au = (mass * (period_d / config.YEAR_D) ** 2) ** (1.0 / 3.0)
    teq = teff * np.sqrt(radius * config.R_SUN_M / (2 * a_au * config.AU_M))
    insolation = lum / a_au ** 2

    def mean_std(x):
        return float(np.median(x)), float(np.std(x))

    payload = {
        "measured_inputs": inputs.record(["transit.period_koi_d"]),
        "inputs": {"twomass": tm.model_dump(), "allwise": wise.model_dump(),
                   "distance_pc": dist.model_dump(), "apogee": apogee.model_dump(),
                   "gaia_parallax_mas": [host_gaia["parallax"], host_gaia["parallax_error"]],
                   "dr25_radius_rsun": dr25.rad, "dr25_teff_k": dr25.teff,
                   "berger2020_radius_rsun": berger.rad,
                   "transit_fit": TRANSIT_FIT.relative_to(config.ROOT).as_posix()},
        "a_ks_rjce": a_ks_rjce,
        "teff_k": spread(teff), "luminosity_lsun": spread(lum),
        "radius_rsun": spread(radius), "mass_msun": spread(mass),
        "density_solar": spread(density),
        "rp_over_rstar": ror,
        "planet_radius_earth": {
            "dr25_star": ror["median"] * dr25.rad * config.R_EARTH_PER_R_SUN,
            "berger2020_star": ror["median"] * berger.rad * config.R_EARTH_PER_R_SUN,
            "adopted_star": spread(rp_adopted),
        },
        "planet_radius_rsun": spread(rp_adopted / config.R_EARTH_PER_R_SUN),
        # median and standard deviation, as the paper quotes the orbit
        "orbital_distance_au": mean_std(a_au),
        "equilibrium_temperature_k": mean_std(teq),
        "insolation_earth": mean_std(insolation),
        "distance_pc": dist.rgeo,
        "seed": config.HOST_SEED, "draws": n,
    }
    write("07_host_star", payload)

    rp = payload["planet_radius_earth"]
    table = Table("The host star, the planet radius it implies, and the orbit",
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
    table("orbital distance (AU)", "%.4f +- %.4f" % payload["orbital_distance_au"])
    table("equilibrium temperature (K)", "%.0f +- %.0f" % payload["equilibrium_temperature_k"])
    table("insolation (Earth)", "%.1f +- %.1f" % payload["insolation_earth"])
    table.show("07_host_star")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
