"""Step 13 - the combined centroid offset and the positional probability from Gaia DR3.

Two measurements of where the light changed are combined here, one axis at a
time: DR25's difference-image offset from the catalogue position, and this
work's PRF fit to the sixteen difference images (data/centroid/prf_fit.json; the
fit itself is not reproduced, see README). They are combined two ways.

Section 3.1 adds the per-axis systematic to each input first, combines by
inverse variance, and quotes the length of the resulting offset.

Section 7, Equation 1, combines them by their formal errors and adds the
systematic afterwards. Which catalogued star best explains that position? Every
Gaia DR3 source within 20 arcsec is a candidate host. A source i carrying a
fraction f_i of the summed G-band flux is admitted only if it could produce the
observed depth when totally eclipsed (depth <= f_i). Each admitted source gets
the same prior weight:

    P_i = L_i / sum_j L_j
    L_i = exp[ -(x_c - x_i)^2 / (2 s_x^2) - (y_c - y_i)^2 / (2 s_y^2) ]
    s^2 = s_fit^2 + s_sys^2

The same calculation is repeated on the Section 3.1 centroid, whose errors
already hold the systematic, to show the choice does not matter.

Catalogued sources only: a star closer to the target than Gaia resolves is not
among them. That case belongs to the blend scenarios of the false-positive
calculations (fpp/).
"""

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config
from koi501.models import GaiaSource
from koi501.report import Table, flux, write

PRF_FIT = config.DATA / "centroid" / "prf_fit.json"
GAIA_RENAME = {"Source": "source_id", "RA_ICRS": "ra", "DE_ICRS": "dec",
               "Gmag": "phot_g_mean_mag"}


def combine(v1: float, e1: float, v2: float, e2: float) -> tuple[float, float]:
    """Inverse-variance mean of two measurements and its error."""
    w1, w2 = 1.0 / e1 ** 2, 1.0 / e2 ** 2
    return (v1 * w1 + v2 * w2) / (w1 + w2), 1.0 / math.sqrt(w1 + w2)


def score(sources: list[dict], depth: float, x: float, y: float,
          sx: float, sy: float) -> list[dict]:
    """Equation 1: each source's distance from the centroid and its probability."""
    rows = []
    for s in sources:
        sigma = math.hypot((x - s["dx_as"]) / sx, (y - s["dy_as"]) / sy)
        capable = depth <= s["flux_fraction"]
        rows.append({"source_id": s["source_id"], "capable": capable,
                     "sigma_from_centroid": sigma,
                     "likelihood": math.exp(-0.5 * sigma ** 2) if capable else 0.0})
    norm = sum(r["likelihood"] for r in rows)
    for r in rows:
        r["probability"] = r["likelihood"] / norm
    return rows


def summary(sources: list[dict], rows: list[dict]) -> dict:
    """Probability on the target (the nearest source) and the closest rival."""
    j = min(range(1, len(rows)), key=lambda i: rows[i]["sigma_from_centroid"])
    return {"p_on_target": rows[0]["probability"],
            "p_all_others": sum(r["probability"] for r in rows[1:]),
            "nearest_alternative": {
                "source_id": sources[j]["source_id"], "G": sources[j]["G"],
                "sep_as": sources[j]["sep_as"],
                "sigma_from_centroid": rows[j]["sigma_from_centroid"]}}


def main() -> bool:
    row = archives.nasa_tap(
        "select koi_dikco_mra,koi_dikco_mra_err,koi_dikco_mdec,koi_dikco_mdec_err,"
        f"koi_depth from q1_q17_dr25_koi where kepoi_name='{config.KOI}'")[0]
    dr25 = {k: float(v) for k, v in row.items()}
    depth = dr25["koi_depth"] * 1e-6
    sys_as = config.CENTROID_SYSTEMATIC_AS

    prf = json.loads(PRF_FIT.read_text(encoding="utf8"))
    # the fit's own scatter, split equally between the two axes
    prf_err = prf["sigma_as"] / math.sqrt(2.0)
    inputs = {"x": (dr25["koi_dikco_mra"], dr25["koi_dikco_mra_err"], prf["east_as"]),
              "y": (dr25["koi_dikco_mdec"], dr25["koi_dikco_mdec_err"], prf["north_as"])}

    # Section 3.1: systematic on each input, then combined as a vector
    vector = {axis: combine(d, math.hypot(e, sys_as), p, math.hypot(prf_err, sys_as))
              for axis, (d, e, p) in inputs.items()}
    (xv, exv), (yv, eyv) = vector["x"], vector["y"]
    offset = math.hypot(xv, yv)
    offset_err = math.hypot(xv * exv, yv * eyv) / offset

    # Equation 1: formal errors, systematic added afterwards
    formal = {axis: combine(d, e, p, prf_err) for axis, (d, e, p) in inputs.items()}
    (xc, exc), (yc, eyc) = formal["x"], formal["y"]
    sx, sy = math.hypot(exc, sys_as), math.hypot(eyc, sys_as)

    gaia = archives.validate(
        GaiaSource,
        archives.vizier_cone("I/355/gaiadr3", "Source,RA_ICRS,DE_ICRS,Gmag",
                             config.RA, config.DEC, 30.0),
        GAIA_RENAME)
    cos_dec = math.cos(math.radians(config.DEC))
    sources = []
    for s in gaia:
        if s.phot_g_mean_mag is None:
            continue
        dx = (s.ra - config.RA) * cos_dec * 3600.0
        dy = (s.dec - config.DEC) * 3600.0
        if math.hypot(dx, dy) <= config.POSPROB_RADIUS_AS:
            sources.append({"source_id": s.source_id, "G": s.phot_g_mean_mag,
                            "dx_as": dx, "dy_as": dy, "sep_as": math.hypot(dx, dy)})
    sources.sort(key=lambda s: s["sep_as"])
    total = sum(flux(s["G"]) for s in sources)
    for s in sources:
        s["flux_fraction"] = flux(s["G"]) / total

    eq1 = score(sources, depth, xc, yc, sx, sy)
    alt = score(sources, depth, xv, yv, exv, eyv)
    for s, a, b in zip(sources, eq1, alt):
        s["capable"] = a["capable"]
        s["equation_1"] = {k: a[k] for k in ("sigma_from_centroid", "probability")}
        s["total_error_weighting"] = {k: b[k] for k in ("sigma_from_centroid", "probability")}

    payload = {
        "dr25_offset_as": {"ra": [dr25["koi_dikco_mra"], dr25["koi_dikco_mra_err"]],
                           "dec": [dr25["koi_dikco_mdec"], dr25["koi_dikco_mdec_err"]]},
        "prf_fit": {k: prf[k] for k in ("east_as", "north_as", "sigma_as", "n_panels")},
        "systematic_per_axis_as": sys_as,
        "section_3_1": {"east_as": [xv, exv], "north_as": [yv, eyv],
                        "offset_as": offset, "offset_err_as": offset_err,
                        "offset_sigma": offset / offset_err},
        "equation_1": {"centroid_east_as": [xc, exc], "centroid_north_as": [yc, eyc],
                       "combined_error_as": [sx, sy], **summary(sources, eq1)},
        "total_error_weighting": {"centroid_east_as": [xv, exv],
                                  "centroid_north_as": [yv, eyv],
                                  **summary(sources, alt)},
        "depth_ppm": dr25["koi_depth"],
        "radius_as": config.POSPROB_RADIUS_AS,
        "n_sources": len(sources),
        "n_capable": sum(s["capable"] for s in sources),
        "includes_unresolved_sources": False,
        "sources": sources,
    }
    write("13_gaia_positional_probability", payload)

    main_result, other = payload["equation_1"], payload["total_error_weighting"]
    near, near_alt = main_result["nearest_alternative"], other["nearest_alternative"]
    table = Table("The combined centroid, and the positional probability from Gaia DR3",
                  "Section 3.1; Section 7, Equation 1")
    table("combined offset, systematic on each input (arcsec)",
          f"{offset:.3f} +/- {offset_err:.3f}")
    table("  in standard errors", round(offset / offset_err, 2))
    table("Gaia DR3 sources within 20 arcsec", len(sources))
    table("  bright enough for the depth if totally eclipsed", payload["n_capable"])
    table("Equation 1 centroid, east (arcsec)", f"{xc:+.3f} +/- {exc:.3f}")
    table("Equation 1 centroid, north (arcsec)", f"{yc:+.3f} +/- {eyc:.3f}")
    table("  systematic then added per axis (arcsec)", sys_as)
    table("probability on KOI-501.01",
          "> 0.999" if main_result["p_on_target"] > 0.999 else main_result["p_on_target"])
    table("nearest alternative, standard errors from centroid",
          f"{near['sigma_from_centroid']:.1f}")
    table("  its separation (arcsec) / G (mag)", f"{near['sep_as']:.2f} / {near['G']:.2f}")
    table("same with the Section 3.1 centroid: probability on target",
          "> 0.999" if other["p_on_target"] > 0.999 else other["p_on_target"])
    table("  nearest alternative, standard errors",
          f"{near_alt['sigma_from_centroid']:.1f}")
    table.show("13_gaia_positional_probability")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
