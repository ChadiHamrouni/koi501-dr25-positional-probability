"""The transit fit behind the radius ratio of Section 4.

Downloads the Kepler long-cadence PDCSAP light curve of KIC 4951877 from MAST,
fits a limb-darkened transit by MCMC, and writes the posterior that step 07 of
run_all.py reads.

Model: batman quadratic limb darkening, 29.4-minute exposures supersampled 15
times, circular orbit, period fixed at the DR25 value. Each transit window of
+-4 durations is divided by a straight line through its out-of-transit points.
Free parameters: Rp/R*, impact parameter b, stellar density rho (solar units),
epoch shift, Kipping q1 and q2, and a multiplicative flux-error scale.

Runs, on identical data:
  adopted   rho ~ N(step 07 density), the spectroscopic stellar solution (the one used)
  free_rho  rho ~ U(0.01, 5)
  fixed_ld  as adopted, with u1 and u2 held fixed (config.TRANSIT_FIT_FIXED_LD)

    pip install -r fit/requirements.txt
    python fit/transit_fit.py                 # all three runs, about an hour
    python fit/transit_fit.py adopted         # one run

Writes data/transit_fit/transit_fit.json and caches the light curve in
results/.cache/kic4951877_lc.npz.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import batman
import emcee
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from koi501 import config, inputs  # noqa: E402

LC_CACHE = ROOT / "results" / ".cache" / "kic4951877_lc.npz"
OUT = ROOT / "data" / "transit_fit" / "transit_fit.json"

KIC = config.KEPID
USED = ["transit.period_koi_d", "transit.epoch_bkjd", "transit.duration_h",
        "star.density_solar", "star.density_err_solar", "star.radius_rsun", "star.radius_err_rsun"]
P = inputs.get("transit.period_koi_d")
T0 = inputs.get("transit.epoch_bkjd")
DUR_D = inputs.get("transit.duration_h") / 24.0
CADENCE_D = config.LONG_CADENCE_D
RHO_ADOPTED = (inputs.get("star.density_solar"), inputs.get("star.density_err_solar"))
R_STAR = (inputs.get("star.radius_rsun"), inputs.get("star.radius_err_rsun"))
G_CGS = config.G_SI * 1000
LIMITS = config.TRANSIT_FIT_LIMITS
N_PARAMETERS = len(config.TRANSIT_FIT_START)

SEED = config.TRANSIT_FIT_SEED
NWALKERS, NBURN, NSTEPS = config.TRANSIT_FIT_SAMPLER
RUNS = {"adopted": ("adopted", False), "free_rho": ("free", False),
        "fixed_ld": ("adopted", True)}


def download_light_curve() -> dict[str, np.ndarray]:
    """Every long-cadence quarter, PDCSAP flux, each divided by its median."""
    if LC_CACHE.exists():
        return dict(np.load(LC_CACHE))
    import lightkurve as lk

    files = lk.search_lightcurve(f"KIC {KIC}", mission="Kepler", cadence="long").download_all(
        flux_column="pdcsap_flux", quality_bitmask="none")
    t, f, e, q = [], [], [], []
    for lc in files:
        flux = np.asarray(lc.flux.value, float)
        median = float(np.nanmedian(flux))
        if not np.isfinite(median) or median <= 0:
            continue
        t.append(np.asarray(lc.time.value, float))
        f.append(flux / median)
        e.append(np.asarray(lc.flux_err.value, float) / median)
        q.append(np.asarray(lc.quality.value, float))
    t, f, e, q = (np.concatenate(x) for x in (t, f, e, q))
    order = np.argsort(t)
    t, f, e, q = t[order], f[order], e[order], q[order]
    keep = np.isfinite(t) & np.isfinite(f) & np.isfinite(e)
    data = {"time": t[keep], "flux": f[keep].astype(np.float32),
            "flux_err": e[keep].astype(np.float32), "quality": q[keep].astype(np.int32),
            "n_quarters": np.int64(len(files))}
    LC_CACHE.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(LC_CACHE, **data)
    return data


def transit_windows(lc: dict[str, np.ndarray]):
    t, f = lc["time"], lc["flux"].astype(float)
    e, q = lc["flux_err"].astype(float), lc["quality"]
    ok = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (q == 0) & (f > 0)
    t, f, e = t[ok], f[ok], e[ok]
    xs, ys, es, n_windows, n_covered = [], [], [], 0, 0
    for n in range(int(np.floor((t.min() - T0) / P)) - 1, int(np.ceil((t.max() - T0) / P)) + 2):
        centre = T0 + n * P
        w = np.abs(t - centre) < config.TRANSIT_FIT_WINDOW_DURATIONS * DUR_D
        if w.sum() < config.TRANSIT_FIT_MIN_POINTS:
            continue
        dt, ff, ee = t[w] - centre, f[w], e[w]
        out = np.abs(dt) > config.OOT_DURATIONS * DUR_D
        if out.sum() < config.TRANSIT_FIT_MIN_OOT:
            continue
        base = np.polyval(np.polyfit(dt[out], ff[out], 1), dt)
        xs.append(dt)
        ys.append(ff / base)
        es.append(ee / base)
        n_windows += 1
        n_min, fraction = config.TRANSIT_FIT_COVER
        n_covered += int(np.sum(np.abs(dt) < fraction * DUR_D) >= n_min)
    return (np.concatenate(xs), np.concatenate(ys), np.concatenate(es),
            n_windows, n_covered)


def a_over_r(rho_solar):
    return (G_CGS * rho_solar * config.RHO_SUN_CGS * (P * 86400.0) ** 2 / (3 * np.pi)) ** (1.0 / 3.0)  # literal: Kepler's third law


class TransitLikelihood:
    def __init__(self, x, y, e):
        self.y, self.e = y, e
        self.params = batman.TransitParams()
        start = config.BATMAN_START
        self.params.t0, self.params.per, self.params.rp = 0.0, P, start["rp"]
        self.params.a, self.params.inc, self.params.ecc, self.params.w = start["a"], start["inc"], 0.0, 90.0  # literal: circular orbit
        self.params.u, self.params.limb_dark = list(start["u"]), "quadratic"
        self.model = batman.TransitModel(self.params, x, supersample_factor=config.TRANSIT_FIT_SUPERSAMPLE,
                                         exp_time=CADENCE_D)

    def light_curve(self, theta, fixed_ld):
        rp, b, rho, dt0, q1, q2, _ = theta
        aor = a_over_r(rho)
        self.params.t0, self.params.rp, self.params.a = dt0, rp, aor
        self.params.inc = np.degrees(np.arccos(b / aor))
        if fixed_ld:
            self.params.u = list(config.TRANSIT_FIT_FIXED_LD)
        else:
            sq = np.sqrt(q1)
            self.params.u = [2.0 * sq * q2, sq * (1.0 - 2.0 * q2)]
        return self.model.light_curve(self.params)

    def __call__(self, theta, prior, fixed_ld):
        rp, b, rho, dt0, q1, q2, s = theta
        if not (LIMITS["ror"][0] < rp < LIMITS["ror"][1] and 0.0 <= b < 1.0 + rp
                and LIMITS["rho"][0] < rho < LIMITS["rho"][1] and abs(dt0) < LIMITS["dt0"]
                and 0.0 < q1 < 1.0 and 0.0 < q2 < 1.0 and LIMITS["scale"][0] < s < LIMITS["scale"][1]):
            return -np.inf
        if b / a_over_r(rho) >= 1.0:
            return -np.inf
        log_prior = 0.0
        if prior == "adopted":
            log_prior = -0.5 * ((rho - RHO_ADOPTED[0]) / RHO_ADOPTED[1]) ** 2
        sigma = s * self.e
        resid = self.y - self.light_curve(theta, fixed_ld)
        return log_prior - 0.5 * np.sum((resid / sigma) ** 2) - np.sum(np.log(sigma))


def summarise(x):
    lo, med, hi = np.percentile(x, config.ONE_SIGMA_PERCENTILES)
    return {"median": float(med), "minus": float(med - lo), "plus": float(hi - med)}


def run(name, likelihood, rng):
    prior, fixed_ld = RUNS[name]
    start = np.array(config.TRANSIT_FIT_START)
    p0 = start + rng.normal(0, config.TRANSIT_FIT_START_SPREAD, (NWALKERS, N_PARAMETERS))
    p0[:, 1] = np.abs(p0[:, 1])
    sampler = emcee.EnsembleSampler(NWALKERS, N_PARAMETERS, likelihood, args=(prior, fixed_ld))
    sampler.random_state = np.random.RandomState(SEED).get_state()
    tic = time.time()
    state = sampler.run_mcmc(p0, NBURN, progress=False)
    sampler.reset()
    sampler.run_mcmc(state, NSTEPS, progress=False)
    chain = sampler.get_chain(flat=True)
    try:
        tau = sampler.get_autocorr_time(tol=0).tolist()
    except emcee.autocorr.AutocorrError:
        tau = None
    rp, b, rho = chain[:, 0], chain[:, 1], chain[:, 2]
    r_star = rng.normal(*R_STAR, len(rp))
    sq = np.sqrt(chain[:, 4])
    best = chain[np.argmax(sampler.get_log_prob(flat=True))]
    result = {
        "prior": prior, "fixed_ld": fixed_ld,
        "rp_over_rstar": summarise(rp),
        "impact_parameter": summarise(b),
        "stellar_density_solar": summarise(rho),
        "a_over_rstar": summarise(a_over_r(rho)),
        "epoch_shift_min": summarise(chain[:, 3] * 1440.0),
        "u1": None if fixed_ld else summarise(2 * sq * chain[:, 5]),
        "u2": None if fixed_ld else summarise(sq * (1 - 2 * chain[:, 5])),
        "error_scale": summarise(chain[:, 6]),
        "depth_rp2_ppm": summarise(rp ** 2 * config.PPM),
        "planet_radius_earth": summarise(rp * r_star * config.R_EARTH_PER_R_SUN),
        "acceptance_fraction": float(np.mean(sampler.acceptance_fraction)),
        "autocorr_steps": tau,
        "best_fit_rms_ppm": float(np.std(likelihood.y - likelihood.light_curve(best, fixed_ld)) * config.PPM),
        "seconds": time.time() - tic,
    }
    ror = result["rp_over_rstar"]
    print(f"{name:9s} Rp/R* = {ror['median']:.5f} +{ror['plus']:.5f} -{ror['minus']:.5f}   "
          f"planet {result['planet_radius_earth']['median']:.2f} R_earth   "
          f"({result['seconds']:.0f} s)", flush=True)
    return result


def main():
    names = sys.argv[1:] or list(RUNS)
    lc = download_light_curve()
    x, y, e, n_windows, n_covered = transit_windows(lc)
    likelihood = TransitLikelihood(x, y, e)
    results = json.loads(OUT.read_text()) if OUT.exists() else {"runs": {}}
    results.update({
        "data": {"source": f"MAST Kepler long-cadence PDCSAP, KIC {KIC}",
                 "n_quarters": int(lc["n_quarters"]), "n_points": int(len(x)),
                 "n_windows": n_windows, "n_windows_with_transit_points": n_covered,
                 "period_d": P, "t0_bkjd": T0, "window_durations": config.TRANSIT_FIT_WINDOW_DURATIONS,
                 "exposure_d": CADENCE_D, "supersample": config.TRANSIT_FIT_SUPERSAMPLE},
        "sampler": {"walkers": NWALKERS, "burn": NBURN, "steps": NSTEPS, "seed": SEED},
        "stellar_radius_rsun": list(R_STAR),
        "measured_inputs": inputs.record(USED),
    })
    for name in names:
        results["runs"][name] = run(name, likelihood, np.random.default_rng(SEED))
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(json.dumps(results, indent=2))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
