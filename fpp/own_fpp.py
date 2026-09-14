"""The third false-positive calculation of Section 7, described in Appendix A.

Seven scenarios are weighed by one Monte Carlo: a planet on the target (P); an
eclipsing binary on the target (EB); a hierarchical triple, a bound companion
eclipsed by its own companion (HEB); a background eclipsing binary within 15.9
arcsec (BEB); and each binary scenario detected at half its true period (EB2,
HEB2, BEB2). For every draw the scenario is built physically - which star is
eclipsed, by what, on what orbit - and the depth, duration, transit shape,
absent secondary eclipse and odd-even statistic it predicts are scored against
the measured values. The difference-image centroid weighs the on-target
scenarios against the background one. Priors are occurrence rates only: the
Raghavan et al. (2010) binary and triple fractions and period distribution, the
Kepler candidate rate per period window, and the TRILEGAL star density.

Run A uses the measured inputs and the constraint set of the calibrated
calculation. Run B multiplies the scenario weights by three factors from the
follow-up data: the APOGEE velocities against a binary on the target (a floor
of 1e-30 on a likelihood that underflows), the same velocities against a bound
pair inside 11.3 AU (which removes part of the triple prior), and the J-band
contrast curve against a background pair (which removes part of its sky area).

Each run adds a flat 0.06% for uncatalogued contaminating stars (Morton et al.
2016) and 0.14% for an instrumental origin. Twenty repeats of 2,000,000 draws.

Inputs
  fpp/data/trilegal_starfield.npz  Kepler_mag, Mact, logg of the TRILEGAL field
                                   simulated for vespa (one square degree)
  fpp/data/cumulative_koi.csv      the cumulative KOI table, 2026-08-22 snapshot
  results/12_contrast_curves.json  the UKIRT J contrast curve (step 12)

    pip install -r fpp/requirements.txt
    python fpp/own_fpp.py               # about N minutes
    python fpp/own_fpp.py some/dir      # write somewhere else

Writes results/fpp_own.json.
"""

from __future__ import annotations

import json
import sys
from math import erf, sqrt
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "results"

SEED = 20260830
RUN_SEED_OFFSET = 50100
REPEATS = 20
N = 2_000_000

G_CGS, R_SUN, M_SUN, R_EARTH, DAY = 6.674e-8, 6.957e10, 1.989e33, 6.371e8, 86400.0
RSUN_RE = R_SUN / R_EARTH
TEFF_SUN = 5772.0

# Raghavan et al. (2010): binary fraction, log-normal period distribution in
# log10(P/d), triple fraction
F_BINARY, LOGP_MU, LOGP_SIG, F_TRIPLE = 0.44, 5.03, 2.28, 0.11
N_TARGETS = 190_000                  # Kepler stars searched
BLEND_RADIUS_AS = 4.0 * 3.98         # four Kepler pixels
SEC_NSIGMA = 7.0                     # depth at which a secondary counts as detected
ECC_PLANET = (0.867, 3.03)           # Kipping (2013) Beta distribution
ECC_BINARY = (1.0, 2.0)
T_INT = 0.0204340                    # long-cadence integration, d
OE_SCALE = 2.43                      # spread of the odd-even statistic on confirmed planets
P_CONTAM_RESIDUAL = 0.0006           # Morton et al. (2016)
P_INSTRUMENT = 0.0014                # Thompson et al. (2018) curated inverted/scrambled sets
R_STELLAR_MIN = 0.086                # smallest hydrogen-burning star, R_sun
Q_MIN = 0.06                         # smallest companion mass ratio
F_RUWE = 0.33                        # credit for Gaia astrometry showing no companion
CEN_FLOOR = np.exp(-0.5 * 5.0 ** 2)  # never more than 5 sigma from the centroid

# KOI-501.01: DR25 transit and vetting measurements, and the host's temperature
TARGET = dict(P=24.79630, dur=8.4830, dur_err=0.4484, depth=575.3e-6,
              depth_err=15.47e-6, kepmag=14.612, wst=44.89, wst_rob=2.60,
              oe=0.1114, teff=5764.0)

# the spectroscopic stellar solution (Section 4)
R_MED, R_SD = 1.746, 0.060
LOGG_MED, LOGG_SD = 4.054, 0.035

# follow-up constraints of run B
EB_FLOOR = 1e-30                     # velocity likelihood floor for an EB on the target
A_MIN_BOUND_AU = 11.3                # bound stellar pairs excluded inside this
M_TRIPLE_MIN = 1.26 + 0.16           # host plus the lightest bound pair, M_sun


def cdfn(x, mu, s):
    return 0.5 * (1 + erf((x - mu) / (s * sqrt(2))))


def mass_radius(m):
    """Main-sequence radius in solar units, floored at the smallest real star."""
    return np.maximum(m ** 0.92, R_STELLAR_MIN)


def flux_ratio(m2, m1):
    """Kepler-band flux ratio of a main-sequence companion, L ~ M^4.5."""
    return (m2 / m1) ** 4.5


def ms_teff(m):
    """Main-sequence temperature from L = M^4.5 and R = M^0.92."""
    return TEFF_SUN * ((m ** 4.5) / mass_radius(m) ** 2) ** 0.25


def blocked(b, k):
    """Fraction of a uniform disc hidden by a disc of radius ratio k at
    separation b (in primary radii)."""
    b = np.asarray(b, float)
    k = np.clip(np.asarray(k, float), 1e-9, 1.0)
    shape = np.broadcast(b, k).shape
    bb, kk = np.broadcast_to(b, shape), np.broadcast_to(k, shape)
    full = bb <= 1.0 - kk
    out = np.where(full, kk ** 2, np.zeros(shape))
    part = (~full) & (bb < 1.0 + kk)
    if part.any():
        bp, kp = bb[part], kk[part]
        c1 = np.clip((bp ** 2 + 1.0 - kp ** 2) / (2 * bp), -1, 1)
        c2 = np.clip((bp ** 2 + kp ** 2 - 1.0) / (2 * bp * kp), -1, 1)
        tri = np.sqrt(np.clip((-bp + kp + 1) * (bp + kp - 1)
                              * (bp - kp + 1) * (bp + kp + 1), 0, None))
        out[part] = (np.arccos(c1) + kp ** 2 * np.arccos(c2) - 0.5 * tri) / np.pi
    return out


def centroid_terms(off, err, ap_area):
    """Centroid likelihood for a source on the target (offset zero) and for a
    blend anywhere inside the search area (uniform in area)."""
    on = max(np.exp(-0.5 * (off / err) ** 2) / (err * np.sqrt(2 * np.pi)), CEN_FLOOR)
    smax = np.sqrt(ap_area / np.pi)
    s = np.linspace(0, smax, 4000)
    kern = np.exp(-0.5 * ((off - s) / err) ** 2) / (err * np.sqrt(2 * np.pi))
    blend = max(float(np.trapezoid(kern * 2 * s / smax ** 2, s)), 1e-300)
    return on, blend


def score(k, aR, dil, P, t, rng, qsec=None, sec=None, ecc=ECC_PLANET,
          double=False, oe=None, teff1=None, teff2=None, shape=None):
    """Mean likelihood of the measured depth, duration, shape, odd-even
    statistic and secondary-eclipse search over isotropic, eccentric orbits."""
    n = len(k)
    if double:
        aR = aR * 2.0 ** (2.0 / 3.0)
        P = 2.0 * P
    cosi = rng.random(n)
    sini = np.sqrt(np.clip(1.0 - cosi ** 2, 1e-12, 1.0))
    e = np.clip(rng.beta(ecc[0], ecc[1], n), 0.0, 0.95)
    esw = e * np.sin(rng.uniform(0, 2 * np.pi, n))
    rfac = (1.0 - e ** 2) / (1.0 + esw)      # separation at transit / semi-major axis
    b = aR * cosi * rfac
    dep = blocked(b, k) * dil
    if double:
        q_ = np.clip(qsec, 1e-6, 1.0)
        dep_sec_ = dep * q_ ** 2.66
        oe_delta = np.abs(dep - dep_sec_)
        dep = 0.5 * (dep + dep_sec_)
    else:
        oe_delta = np.zeros_like(dep)
    # Winn (2010) eq. 16
    arg = np.sqrt(np.clip((1.0 + k) ** 2 - b ** 2, 0, None)) / (aR * sini)
    dur = ((P * 24.0 / np.pi) * np.arcsin(np.clip(arg, 0, 1))
           * np.sqrt(1.0 - e ** 2) / (1.0 + esw))
    ecl = b < (1.0 + k)

    if shape is not None:
        chord_t = np.sqrt(np.clip((1.0 + k) ** 2 - b ** 2, 0, None))
        chord_f = np.sqrt(np.clip((1.0 - k) ** 2 - b ** 2, 0, None))
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = np.where(chord_t > 0, (chord_t - chord_f) / (2 * chord_t), 0.5)
        tau_obs = np.clip(dur * frac, 0, None) + T_INT * 24.0
        ttau = np.where(tau_obs > 0, dur / tau_obs, 2.0)
        l_shape = np.exp(-0.5 * ((shape[0] - ttau) / shape[1]) ** 2)
    else:
        l_shape = 1.0

    ld = np.exp(-0.5 * ((dep - t["depth"]) / t["depth_err"]) ** 2) \
        / (t["depth_err"] * np.sqrt(2 * np.pi))
    lt = np.exp(-0.5 * ((dur - t["dur"]) / t["dur_err"]) ** 2) \
        / (t["dur_err"] * np.sqrt(2 * np.pi))
    v = np.where(ecl, ld * lt * l_shape, 0.0)
    if oe is not None and np.isfinite(oe):
        v = v * np.exp(-0.5 * ((oe - oe_delta / (2.0 * t["depth_err"])) / OE_SCALE) ** 2)
    if double:
        sec = None
    if sec is not None and np.isfinite(sec[1]) and sec[1] > 0:
        sec_obs, sec_sigma, detected = sec
        if qsec is None:
            dep_sec = 0.0
        elif teff1 is not None and teff2 is not None:
            dep_sec = dep * np.clip(teff2 / teff1, 0.0, None) ** 4
        else:
            dep_sec = dep * np.clip(qsec, 1e-6, 1.0) ** 2.66
        if detected:
            v = v * np.exp(-0.5 * ((sec_obs - dep_sec) / sec_sigma) ** 2)
        else:
            v = v * ndtr(SEC_NSIGMA - dep_sec / sec_sigma)
    return float(v.mean())


def stellar_samples(n, rng):
    """Radius from its measurement; mass from log g and that radius."""
    R = rng.normal(R_MED, R_SD, n)
    g = 10 ** rng.normal(LOGG_MED, LOGG_SD, n) / 100.0
    M = g * (R * 6.957e8) ** 2 / 6.67430e-11 / 1.98892e30
    return M, R


def evaluate(t, M1, R1, field, koi, rng, n=N):
    """Scenario weights (prior times likelihood) for one repeat."""
    P, depth = t["P"], t["depth"]
    rho1 = M1 * M_SUN / (4 / 3 * np.pi * (R1 * R_SUN) ** 3)
    aR1 = (G_CGS * rho1 * (P * DAY) ** 2 / (3 * np.pi)) ** (1 / 3)

    kmag, Mbg = field["Kepler_mag"], field["Mact"]
    Rbg = np.sqrt(G_CGS * Mbg * M_SUN / 10 ** field["logg"]) / R_SUN
    fbg = 10 ** (-0.4 * (kmag - t["kepmag"]))
    possible = (fbg / (1.0 + fbg)) > depth
    ap_area = np.pi * BLEND_RADIUS_AS ** 2
    n_bg = possible.sum() * ap_area / 3600.0 ** 2

    f_per = cdfn(np.log10(P * 1.5), LOGP_MU, LOGP_SIG) - cdfn(np.log10(P / 1.5), LOGP_MU, LOGP_SIG)
    f_per2 = (cdfn(np.log10(2 * P * 1.5), LOGP_MU, LOGP_SIG)
              - cdfn(np.log10(2 * P / 1.5), LOGP_MU, LOGP_SIG))
    a_cm = (G_CGS * M_SUN * (P * DAY) ** 2 / (4 * np.pi ** 2)) ** (1 / 3)
    candidates = koi.koi_disposition.isin(["CONFIRMED", "CANDIDATE"])
    n_win = int((candidates & koi.koi_period.between(P / 1.5, P * 1.5)).sum())
    prior = dict(P=n_win / (N_TARGETS * (R_SUN / a_cm)),
                 EB=F_BINARY * f_per, BEB=n_bg * F_BINARY * f_per, HEB=F_TRIPLE * f_per,
                 EB2=F_BINARY * f_per2, BEB2=n_bg * F_BINARY * f_per2, HEB2=F_TRIPLE * f_per2)

    sec = (t["wst"] * 1e-6, t["wst"] * 1e-6 / t["wst_rob"], t["wst_rob"] >= SEC_NSIGMA)
    oe, shape, teff = t["oe"], t["shape"], t["teff"]
    radii = koi.koi_prad[candidates & koi.koi_prad.between(0.3, 25.0)].to_numpy()

    # planet on the target
    Rp = rng.choice(radii, n) / RSUN_RE
    L_P = score(Rp / R1, aR1, 1.0, P, t, rng, sec=sec, oe=oe, shape=shape)

    # eclipsing binary on the target
    q2 = rng.uniform(Q_MIN, 1.0, n)
    M2 = q2 * M1
    R2, f2 = mass_radius(M2), flux_ratio(M2, M1)
    kw = dict(qsec=q2, ecc=ECC_BINARY, oe=oe, shape=shape, teff1=teff, teff2=ms_teff(M2))
    L_EB = score(R2 / R1, aR1, 1.0 / (1.0 + f2), P, t, rng, sec=sec, **kw)
    L_EB2 = score(R2 / R1, aR1, 1.0 / (1.0 + f2), P, t, rng, double=True, **kw)

    # hierarchical triple, diluted by the target's light
    qb = rng.uniform(0.10, 1.0, n)
    Mb = qb * M1
    Rb, fb = mass_radius(Mb), flux_ratio(Mb, M1)
    aRb = (G_CGS * (Mb * M_SUN / (4 / 3 * np.pi * (Rb * R_SUN) ** 3))
           * (P * DAY) ** 2 / (3 * np.pi)) ** (1 / 3)
    q3 = rng.uniform(Q_MIN, 1.0, n)
    R3 = mass_radius(q3 * Mb)
    kw = dict(qsec=q3, ecc=ECC_BINARY, oe=oe, shape=shape,
              teff1=ms_teff(Mb), teff2=ms_teff(q3 * Mb))
    L_HEB = score(R3 / Rb, aRb, fb / (1.0 + fb), P, t, rng, sec=sec, **kw)
    L_HEB2 = score(R3 / Rb, aRb, fb / (1.0 + fb), P, t, rng, double=True, **kw)

    # background eclipsing binary, drawn from the TRILEGAL stars bright enough
    j = rng.choice(np.flatnonzero(possible), n)
    Mg, Rg, fg = Mbg[j], Rbg[j], fbg[j]
    aRg = (G_CGS * (Mg * M_SUN / (4 / 3 * np.pi * (Rg * R_SUN) ** 3))
           * (P * DAY) ** 2 / (3 * np.pi)) ** (1 / 3)
    qg = rng.uniform(Q_MIN, 1.0, n)
    Rgc = mass_radius(qg * Mg)
    kw = dict(qsec=qg, ecc=ECC_BINARY, oe=oe, shape=shape,
              teff1=ms_teff(Mg), teff2=ms_teff(qg * Mg))
    L_BEB = score(Rgc / Rg, aRg, fg / (1.0 + fg), P, t, rng, sec=sec, **kw)
    L_BEB2 = score(Rgc / Rg, aRg, fg / (1.0 + fg), P, t, rng, double=True, **kw)

    on, blend = centroid_terms(t["cen"], t["cen_err"], ap_area)
    L = dict(P=L_P * on, EB=L_EB * on, BEB=L_BEB * blend, HEB=L_HEB * on * F_RUWE,
             EB2=L_EB2 * on, BEB2=L_BEB2 * blend, HEB2=L_HEB2 * on * F_RUWE)
    w = {h: prior[h] * L[h] for h in prior}
    tot = sum(w.values())
    return {"w": w, "prior": prior, "n_bg": float(n_bg), "ap_area": ap_area,
            "possible": possible, "p_on_target": 1 - (w["BEB"] + w["BEB2"]) / tot}


def with_allowances(w):
    """Astrophysical FPP from scenario weights, then the two flat allowances."""
    fpp = 1 - w["P"] / sum(w.values())
    fpp = 1 - (1 - fpp) * (1 - P_CONTAM_RESIDUAL)
    return fpp, 1 - (1 - fpp) * (1 - P_INSTRUMENT)


def triple_outer_fraction():
    """Share of the Raghavan period distribution beyond the velocity-drift limit."""
    p_yr = np.sqrt(A_MIN_BOUND_AU ** 3 / M_TRIPLE_MIN)
    return float(1.0 - cdfn(np.log10(p_yr * 365.25), LOGP_MU, LOGP_SIG))


def contrast_area_factor(field, possible, kepmag, ap_area):
    """Mean share of the search area where the J-band curve would not have
    detected each TRILEGAL star bright enough to matter. The curve is in J and
    the magnitudes are Kepler; a background star redder than the target has a
    smaller J contrast, so this overstates the surviving area."""
    cc = json.loads((RESULTS / "12_contrast_curves.json").read_text(encoding="utf8"))["J"]
    seps = np.array(cc["separation_as"])
    dmag = np.maximum.accumulate(np.array(cc["delta_mag_5sigma"]))
    iwa, smax = float(seps.min()), float(np.sqrt(ap_area / np.pi))
    dm = (field["Kepler_mag"] - kepmag)[possible]
    r_det = np.interp(dm, dmag, seps, left=iwa, right=smax)
    r_eff = np.clip(np.maximum(r_det, iwa), 0.0, smax)
    return float((np.pi * r_eff ** 2).sum() / (np.pi * smax ** 2) / max(possible.sum(), 1))


def measured_inputs() -> dict:
    """The centroid offset of step 13 (Section 3.1) and the transit shape of
    step 14, the latter as duration over integration-broadened ingress,
    T / (tau + T_int), at the median duration and the median, 16th and 84th
    percentile ingress, with half the 16-84 spread as its error."""
    offset = json.loads((RESULTS / "13_gaia_positional_probability.json").read_text(
        encoding="utf8"))["section_3_1"]
    boot = json.loads((RESULTS / "14_blend_bounds.json").read_text(encoding="utf8"))["bootstrap"]
    T = boot["duration_d_median"]
    tau16, tau50, tau84 = boot["ingress_d_p16_p50_p84"]
    shape = (T / (tau50 + T_INT), 0.5 * (T / (tau16 + T_INT) - T / (tau84 + T_INT)))
    return dict(cen=offset["offset_as"], cen_err=offset["offset_err_as"], shape=shape)


def main() -> int:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else RESULTS
    field = dict(np.load(HERE / "data" / "trilegal_starfield.npz"))
    koi = pd.read_csv(HERE / "data" / "cumulative_koi.csv", comment="#")
    t = dict(TARGET, **measured_inputs())

    run_a, run_b, first = [], [], None
    for k in range(REPEATS):
        rng = np.random.default_rng(SEED + RUN_SEED_OFFSET + k)
        M1, R1 = stellar_samples(N, rng)
        r = evaluate(t, M1, R1, field, koi, rng)
        if first is None:
            first = r
            factors = dict(P=1.0, EB=EB_FLOOR, EB2=EB_FLOOR,
                           HEB=triple_outer_fraction(), HEB2=triple_outer_fraction(),
                           BEB=contrast_area_factor(field, r["possible"], t["kepmag"],
                                                    r["ap_area"]))
            factors["BEB2"] = factors["BEB"]
        run_a.append(with_allowances(r["w"])[1])
        run_b.append(with_allowances({h: r["w"][h] * factors[h] for h in r["w"]})[1])
        print(f"  repeat {k + 1:2d}   run A {run_a[-1]:.4e}   run B {run_b[-1]:.4e}", flush=True)

    flat = 1 - (1 - P_CONTAM_RESIDUAL) * (1 - P_INSTRUMENT)

    def summary(values):
        v = np.array(values)
        mean = float(v.mean())
        return {"mean": mean, "sem": float(v.std(ddof=1) / np.sqrt(len(v))),
                "max": float(v.max()), "values": [float(x) for x in v],
                "astrophysical_part": 1 - (1 - mean) / (1 - flat)}

    payload = {
        "inputs": {"target": t, "stellar_radius_rsun": [R_MED, R_SD],
                   "stellar_logg": [LOGG_MED, LOGG_SD], "draws": N, "repeats": REPEATS,
                   "seeds": [SEED + RUN_SEED_OFFSET, SEED + RUN_SEED_OFFSET + REPEATS - 1],
                   "blend_radius_as": BLEND_RADIUS_AS, "raghavan": [F_BINARY, LOGP_MU,
                                                                    LOGP_SIG, F_TRIPLE],
                   "ruwe_credit": F_RUWE, "n_trilegal_stars": int(len(field["Kepler_mag"])),
                   "n_koi_rows": int(len(koi))},
        "flat_allowances": {"residual_contamination": P_CONTAM_RESIDUAL,
                            "instrumental": P_INSTRUMENT, "combined": flat},
        "priors": first["prior"], "expected_background_stars": first["n_bg"],
        "p_on_target": float(first["p_on_target"]),
        "run_b_factors": factors,
        "run_a": summary(run_a),
        "run_b": summary(run_b),
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "fpp_own.json").write_text(json.dumps(payload, indent=2), encoding="utf8")
    print(f"run A {payload['run_a']['mean']:.3e} +/- {payload['run_a']['sem']:.1e}   "
          f"run B {payload['run_b']['mean']:.3e} +/- {payload['run_b']['sem']:.1e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
