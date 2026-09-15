"""Step 9 - the dynamical mass limit from nineteen archival APOGEE spectra.

The transit fixes the orbital phase, so a circular Keplerian at the known
ephemeris has one free amplitude. Fitting it to the DR17 visit velocities bounds
the transiting body's mass without ever detecting it. A straight line through
the same velocities bounds the acceleration a wide bound companion would cause,
which sets how close a bound pair of stars could orbit.

Reproduces Section 5.1.
"""

import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config, inputs
from koi501.models import ApogeeVisit
from koi501.report import Table, write

MEASURED = ["transit.period_koi_d", "transit.epoch_bkjd", "star.mass_msun"]


def semi_amplitude(m2_msun: float, m1_msun: float, period_d: float) -> float:
    """Circular-orbit velocity semi-amplitude of the primary, m/s, sin i = 1."""
    p_s = period_d * config.DAY_S
    return ((2 * math.pi * config.G_SI / p_s) ** (1 / 3) * m2_msun * config.M_SUN_KG
            / ((m1_msun + m2_msun) * config.M_SUN_KG) ** (2 / 3))


def semi_amplitude_to_mass(k_ms: float, period_d: float, m_star_msun: float) -> float:
    """Companion mass in Jupiter masses from a circular-orbit semi-amplitude."""
    period_s = period_d * config.DAY_S
    m_star = m_star_msun * config.M_SUN_KG
    scale = (2 * math.pi * config.G_SI / period_s) ** (1 / 3)
    mass = k_ms * m_star ** (2 / 3) / scale
    for _ in range(config.MASS_FUNCTION_ITERATIONS):  # the mass function uses the total mass
        mass = k_ms * (m_star + mass) ** (2 / 3) / scale
    return mass / config.M_JUP_KG


def weighted_line(x, y, w):
    """Weighted least squares for y = a + b x: (a, b, variance of b, residuals)."""
    s0 = sum(w)
    s1 = sum(wi * xi for wi, xi in zip(w, x))
    s2 = sum(wi * xi * xi for wi, xi in zip(w, x))
    t0 = sum(wi * yi for wi, yi in zip(w, y))
    t1 = sum(wi * xi * yi for wi, xi, yi in zip(w, x, y))
    det = s0 * s2 - s1 * s1
    b = (s0 * t1 - s1 * t0) / det
    a = (s2 * t0 - s1 * t1) / det
    residual = [yi - (a + b * xi) for xi, yi in zip(x, y)]
    return a, b, s0 / det, residual


def main() -> bool:
    period = inputs.get("transit.period_koi_d")
    m_star = inputs.get("star.mass_msun")
    raw = archives.vizier_cone(
        "III/286/allvis", "JD,VHelio,e_RV,Ncomp,SNR,Date",
        config.RA, config.DEC, config.POINT_MATCH_AS)
    visits = list(archives.validate(
        ApogeeVisit, raw,
        {"JD": "jd", "VHelio": "vhelio_kms", "e_RV": "e_rv_kms",
         "Ncomp": "ncomp"}))

    jd = [v.jd for v in visits]
    rv = [v.vhelio_kms * 1000 for v in visits]
    err = [v.e_rv_kms * 1000 for v in visits]
    ncomp = [v.ncomp for v in visits if v.ncomp is not None]
    dates = sorted(r["Date"][:10] for r in raw if r.get("Date"))
    n = len(visits)

    weight = [1 / e ** 2 for e in err]
    gamma = sum(v * w for v, w in zip(rv, weight)) / sum(weight)
    gamma_err = 1 / math.sqrt(sum(weight))

    epoch = config.BKJD_OFFSET + inputs.get("transit.epoch_bkjd")
    phase = [((t - epoch) / period) % 1.0 for t in jd]

    # gamma + K on a fixed-phase circular orbit: velocity = gamma - K sin(2 pi phase)
    basis = [-math.sin(2 * math.pi * p) for p in phase]
    gamma_fit, k, k_var, residual = weighted_line(basis, rv, weight)
    chi2_red = sum((r / e) ** 2 for r, e in zip(residual, err)) / (n - 2)
    k_err = math.sqrt(k_var * chi2_red)
    limit = abs(k) + config.LIMIT_SIGMA * k_err

    ordered = sorted(phase)
    gaps = [b - a for a, b in zip(ordered, ordered[1:])] + [ordered[0] + 1 - ordered[-1]]

    # a straight line through the velocities: the acceleration of a wide companion
    t_mean = statistics.mean(jd)
    _, drift, drift_var, drift_res = weighted_line([t - t_mean for t in jd], rv, weight)
    drift_scale = max(1.0, math.sqrt(sum(w * r * r for w, r in zip(weight, drift_res)) / (n - 2)))
    drift_limit = (abs(drift) + config.LIMIT_SIGMA * math.sqrt(drift_var) * drift_scale) * config.YEAR_D
    acceleration = drift_limit / (config.YEAR_D * config.DAY_S)          # m s^-2
    pair_min_au = math.sqrt(config.G_SI * config.LIGHTEST_ECLIPSING_PAIR_MSUN * config.M_SUN_KG
                            / acceleration) / config.AU_M

    lo, hi = config.EB_COMPANION_RANGE_MSUN
    eb_ratio = [semi_amplitude(m, m_star, period) / limit for m in (lo, hi)]

    payload = {
        "measured_inputs": inputs.record(MEASURED),
        "n_visits": n, "first_visit": dates[0], "last_visit": dates[-1],
        "baseline_d": max(jd) - min(jd),
        "baseline_periods": (max(jd) - min(jd)) / period,
        "largest_phase_gap": max(gaps),
        "systemic_kms": gamma / 1000, "systemic_err_kms": gamma_err / 1000,
        "scatter_ms": statistics.stdev(rv),
        "error_range_ms": [min(err), max(err)],
        "semi_amplitude_ms": k, "semi_amplitude_err_ms": k_err,
        "chi2_reduced": chi2_red,
        "limit_3sigma_ms": limit,
        "companion_mass_limit_mjup": semi_amplitude_to_mass(limit, period, m_star),
        "single_component_every_visit": set(ncomp) == {1},
        "drift_ms_per_yr": drift * config.YEAR_D,
        "drift_limit_3sigma_ms_per_yr": drift_limit,
        "bound_pair_msun": config.LIGHTEST_ECLIPSING_PAIR_MSUN,
        "bound_pair_excluded_inside_au": pair_min_au,
        "eclipsing_binary_companion_msun": [lo, hi],
        "eclipsing_binary_k_over_limit": eb_ratio,
    }
    write("09_radial_velocities", payload)

    table = Table('Mass limit from the 19 APOGEE spectra', 'Section 5.1')
    table('spectra', n)
    table('first spectrum', payload["first_visit"])
    table('last spectrum', payload["last_visit"])
    table('time span, in orbits of the planet', round(payload["baseline_periods"]))
    table('largest gap in orbital coverage (fraction of an orbit)', round(payload["largest_phase_gap"], 2))
    table("star's mean velocity (km/s)", payload["systemic_kms"])
    table('scatter of the velocities (m/s)', round(payload["scatter_ms"]))
    table('smallest velocity error (m/s)', round(payload["error_range_ms"][0]))
    table('largest velocity error (m/s)', round(payload["error_range_ms"][1]))
    table('fitted wobble K (m/s)', round(k))
    table('  uncertainty on K (m/s)', round(k_err))
    table('largest wobble allowed, 3 sigma (m/s)', round(limit))
    table('largest mass allowed (Jupiter masses)', round(payload["companion_mass_limit_mjup"], 1))
    table('single star seen in every spectrum', payload["single_component_every_visit"])
    table('velocity drift allowed, 3 sigma (m/s per year)', round(drift_limit))
    table('lightest eclipsing pair excluded inside (AU)', round(pair_min_au, 1))
    table('eclipsing binary on the target, times the limit',
          f"{eb_ratio[0]:.0f} to {eb_ratio[1]:.0f}")
    table.show("09_radial_velocities")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
