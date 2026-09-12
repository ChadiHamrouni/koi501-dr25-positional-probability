"""Step 1 - the photometry of the target and its two favoured neighbours.

Reproduces Table 1 and Section 2.2: the KIC magnitudes, the three later surveys
that disagree with them, and the r - J colour that no star has.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config
from koi501.models import (GaiaSource, KICStar, PanstarrsSource,
                           TwoMassSource)
from koi501.report import Table, angsep, flux, write

KIC_COLUMNS = "KIC,RAJ2000,DEJ2000,kepmag,rmag,Jmag"
KIC_RENAME = {"KIC": "kic", "RAJ2000": "ra", "DEJ2000": "dec", "Jmag": "jmag"}
SKY_RENAME = {"RAJ2000": "ra", "DEJ2000": "dec", "Jmag": "jmag"}


def main() -> bool:
    kic = {s.kic: s for s in archives.validate(
        KICStar,
        archives.vizier_cone("V/133/kic", KIC_COLUMNS,
                             config.RA, config.DEC, 30.0),
        KIC_RENAME)}

    gaia = list(archives.validate(GaiaSource, archives.gaia_tap(
        "SELECT source_id, ra, dec, phot_g_mean_mag, parallax, "
        "parallax_error, ruwe, astrometric_excess_noise, "
        "ipd_frac_multi_peak, non_single_star FROM gaiadr3.gaia_source "
        f"WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
        f"CIRCLE('ICRS',{config.RA},{config.DEC},30./3600.))")))

    twomass = list(archives.validate(
        TwoMassSource,
        archives.vizier_cone("II/246/out", "_r,RAJ2000,DEJ2000,Jmag",
                             config.RA, config.DEC, 30.0),
        SKY_RENAME))
    panstarrs = list(archives.validate(
        PanstarrsSource,
        archives.vizier_cone("II/349/ps1", "_r,RAJ2000,DEJ2000,rmag,e_rmag,Nr",
                             config.RA, config.DEC, 12.0),
        {**SKY_RENAME, "e_rmag": "rmag_error", "Nr": "n_detections"}))

    def nearest(sources, ra, dec):
        """The source closest to a position, and how far away it is."""
        best, best_sep = None, 1e9
        for source in sources:
            sep = angsep(source.ra, source.dec, ra, dec)
            if sep < best_sep:
                best, best_sep = source, sep
        return best, best_sep

    stars = []
    host = kic[config.KEPID]
    for kid in (config.KEPID, *config.NEIGHBOURS):
        star = kic[kid]
        g, g_sep = nearest(gaia, star.ra, star.dec)
        tm, tm_sep = nearest(twomass, star.ra, star.dec)
        ps, ps_sep = nearest(panstarrs, star.ra, star.dec)
        stars.append({
            "kic": kid,
            "sep_from_host_as": angsep(star.ra, star.dec, host.ra, host.dec),
            "kic_kepmag": star.kepmag,
            "kic_r": star.rmag,
            "kic_J": star.jmag,
            "r_minus_J": star.r_minus_j,
            "gaia_G": g.phot_g_mean_mag if g else None,
            "gaia_match_as": g_sep,
            "twomass_J": tm.jmag if tm else None,
            "twomass_match_as": tm_sep,
            "panstarrs_r": ps.rmag if ps else None,
            "panstarrs_match_as": ps_sep,
        })

    for s in stars:
        if s["gaia_G"] is not None and s["kic_kepmag"] is not None:
            s["kepmag_minus_G"] = s["kic_kepmag"] - s["gaia_G"]

    shares = {}
    for label, key in (("kic", "kic_kepmag"), ("gaia", "gaia_G")):
        f = [flux(s[key]) for s in stars]
        shares[label] = [x / sum(f) for x in f]

    brightest = min(
        (g for g in gaia
         if g.phot_g_mean_mag is not None
         and angsep(g.ra, g.dec, config.RA, config.DEC) > 1.0
         and angsep(g.ra, g.dec, config.RA, config.DEC) <= config.BLEND_RADIUS),
        key=lambda g: g.phot_g_mean_mag)

    host_gaia = min(gaia, key=lambda g: angsep(g.ra, g.dec, host.ra, host.dec))

    payload = {"stars": stars, "flux_shares": shares,
               "host_gaia": host_gaia.model_dump(),
               "brightest_other_gaia_within_25as": {
                   "G": brightest.phot_g_mean_mag,
                   "sep_as": angsep(brightest.ra, brightest.dec,
                                    config.RA, config.DEC)}}
    write("01_target_photometry", payload)

    host_row, n1, n2 = stars
    table = Table('Brightness of KOI-501.01 and its two neighbours, four catalogues',
                  'Section 2.2, Table 1')
    table('KOI-501.01, Kepler magnitude in the KIC', host_row["kic_kepmag"])
    table('KIC 4951867, Kepler magnitude in the KIC', n1["kic_kepmag"])
    table('KIC 4951861, Kepler magnitude in the KIC', n2["kic_kepmag"])
    table('KIC 4951867, distance from KOI-501.01 (arcsec)', round(n1["sep_from_host_as"], 2))
    table('KIC 4951861, distance from KOI-501.01 (arcsec)', round(n2["sep_from_host_as"], 2))
    table('KOI-501.01, r - J colour in the KIC', round(host_row["r_minus_J"], 2))
    table('KIC 4951867, r - J colour in the KIC', round(n1["r_minus_J"], 2))
    table('KIC 4951867, Gaia G magnitude', n1["gaia_G"])
    table('KIC 4951861, Gaia G magnitude', n2["gaia_G"])
    table('KIC 4951867, too bright in the KIC by (mag)', round(-n1["kepmag_minus_G"], 2))
    table('KIC 4951861, too bright in the KIC by (mag)', round(-n2["kepmag_minus_G"], 2))
    table("KOI-501.01, share of the three stars' light, KIC (%)", round(shares["kic"][0] * 100))
    table("KOI-501.01, share of the three stars' light, Gaia (%)", round(shares["gaia"][0] * 100))
    table('brightest other Gaia star within 25 arcsec (G)', round(payload["brightest_other_gaia_within_25as"]["G"], 2))
    table('KOI-501.01, Gaia RUWE', round(host_gaia.ruwe, 3))
    table('KOI-501.01, Gaia astrometric excess noise', host_gaia.astrometric_excess_noise)
    table('KOI-501.01, Gaia multi-peak fraction', host_gaia.ipd_frac_multi_peak)
    table('KOI-501.01, Gaia non-single-star flag', host_gaia.non_single_star)
    table.show("01_target_photometry")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
