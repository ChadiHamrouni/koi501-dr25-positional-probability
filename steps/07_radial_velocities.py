"""Step 7 - the dynamical mass limit from nineteen archival APOGEE spectra.

The transit fixes the orbital phase, so a circular Keplerian at the known
ephemeris has one free amplitude. Fitting it to the DR17 visit velocities bounds
the transiting body's mass without ever detecting it.

Reproduces Section 5.1.
"""

import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config
from koi501.models import ApogeeVisit
from koi501.report import Table, write

G = 6.674e-11
M_SUN = 1.989e30
M_JUP = 1.898e27
DAY = 86400.0
BKJD_OFFSET = 2454833.0


def semi_amplitude_to_mass(k_ms: float, period_d: float, m_star_msun: float) -> float:
    """Companion mass in Jupiter masses from a circular-orbit semi-amplitude."""
    period_s = period_d * DAY
    m_star = m_star_msun * M_SUN
    scale = (2 * math.pi * G / period_s) ** (1 / 3)
    mass = k_ms * m_star ** (2 / 3) / scale
    for _ in range(8):  # the mass function uses the total mass
        mass = k_ms * (m_star + mass) ** (2 / 3) / scale
    return mass / M_JUP


def main() -> bool:
    raw = archives.vizier_cone(
        "III/286/allvis", "JD,VHelio,e_RV,Ncomp,SNR,Date",
        config.RA, config.DEC, 3.0)
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

    epoch = BKJD_OFFSET + config.EPOCH_BKJD
    phase = [((t - epoch) / config.PERIOD_D) % 1.0 for t in jd]

    # least squares for gamma + K on a fixed-phase circular orbit
    basis = [(1.0, -math.sin(2 * math.pi * p)) for p in phase]
    a11 = sum(w * b[0] * b[0] for w, b in zip(weight, basis))
    a12 = sum(w * b[0] * b[1] for w, b in zip(weight, basis))
    a22 = sum(w * b[1] * b[1] for w, b in zip(weight, basis))
    b1 = sum(w * v * b[0] for w, v, b in zip(weight, rv, basis))
    b2 = sum(w * v * b[1] for w, v, b in zip(weight, rv, basis))
    det = a11 * a22 - a12 * a12
    k = (a11 * b2 - a12 * b1) / det
    k_var = a11 / det

    residual = [v - ((a22 * b1 - a12 * b2) / det + k * b[1])
                for v, b in zip(rv, basis)]
    chi2_red = sum((r / e) ** 2 for r, e in zip(residual, err)) / (n - 2)
    k_err = math.sqrt(k_var * chi2_red)
    limit = abs(k) + 3 * k_err

    ordered = sorted(phase)
    gaps = [b - a for a, b in zip(ordered, ordered[1:])] + [ordered[0] + 1 - ordered[-1]]

    payload = {
        "n_visits": n, "first_visit": dates[0], "last_visit": dates[-1],
        "baseline_d": max(jd) - min(jd),
        "baseline_periods": (max(jd) - min(jd)) / config.PERIOD_D,
        "largest_phase_gap": max(gaps),
        "systemic_kms": gamma / 1000, "systemic_err_kms": gamma_err / 1000,
        "scatter_ms": statistics.stdev(rv),
        "error_range_ms": [min(err), max(err)],
        "semi_amplitude_ms": k, "semi_amplitude_err_ms": k_err,
        "chi2_reduced": chi2_red,
        "limit_3sigma_ms": limit,
        "companion_mass_limit_mjup": semi_amplitude_to_mass(
            limit, config.PERIOD_D, config.STELLAR_MASS),
        "single_component_every_visit": set(ncomp) == {1},
    }
    write("07_radial_velocities", payload)

    table = Table('Mass limit from the 19 APOGEE spectra',
                  'Section 5.1')
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
    table.show("07_radial_velocities")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
