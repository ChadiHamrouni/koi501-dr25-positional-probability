"""The Bayes factor of Section 6: is a transit needed to explain the light curve?

Following Matesic et al. (2024), the photometry around each transit is modelled
twice, with the same noise model in both:

  TGP  a limb-darkened transit plus a Matern-3/2 Gaussian process and a jitter
  GP   the Gaussian process and jitter alone

The evidence of each is computed by nested sampling, and B = ln Z_TGP - ln Z_GP;
B >= 10 is decisive. Two controls follow: the transit priors narrowed to five
posterior standard deviations (a Bayes factor depends on prior volume), and the
same test at an epoch 0.371 periods away, with points near real transits
removed, where it should find nothing.

Data: the long-cadence PDCSAP light curve in data/lightcurve/ (quality flag 0).
Each window of +-2.5 durations is divided by a straight line through its
out-of-transit points and clipped at 5 robust standard deviations. The windows
are laid on one time axis 100 d apart; across that gap the kernel correlation is
exp(-100/rho), zero to machine precision, so one Gaussian-process call is exact.

Model: batman quadratic limb darkening held at u1 = 0.42, u2 = 0.24, 29.4-minute
exposures supersampled 7 times, circular orbit, period fixed at the DR25 value.
UltraNest, 400 live points, frac_remain 0.3. The sampler is not seeded, so a
rerun agrees with the published values within the quoted evidence errors.

    pip install -r fit/requirements.txt
    python fit/bayes_factor.py                      # about 40 minutes
    python fit/bayes_factor.py out.json             # write somewhere else

Writes results/bayes_factor.json by default. Sampler run directories go to
results/.cache/ultranest/.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import batman
import celerite2
import numpy as np
import ultranest
from celerite2 import terms

ROOT = Path(__file__).resolve().parent.parent
LIGHT_CURVE = ROOT / "data" / "lightcurve" / "kic4951877_pdcsap.npz"
RUNDIR = ROOT / "results" / ".cache" / "ultranest"
OUT = ROOT / "results" / "bayes_factor.json"

P = 24.7962801              # d, DR25
T0 = 145.529                # BKJD, DR25
DUR_D = 8.483 / 24.0        # DR25 duration
CADENCE_D = 29.4244 / 60 / 24
WINDOW = 2.5                # window half-width, in durations
LD_U1, LD_U2 = 0.42, 0.24   # quadratic, Kepler band, Teff 5764 K, log g 4.05
CONTROL_PHASE = 0.371       # control epoch offset, in periods

PRIORS_WIDE = {
    "t0": (-0.05, 0.05),          # d
    "log_rp": (-3.0, -1.0),       # log10 Rp/R*
    "aor": (5.0, 60.0),           # a/R*
    "b": (0.0, 0.99),
    "log_sigma": (-6.0, -2.0),    # log10 GP amplitude, relative flux
    "log_rho": (-1.5, 1.5),       # log10 GP timescale, d
    "log_jit": (-5.0, -3.0),      # log10 white-noise jitter, relative flux
}
TGP_NAMES = ["t0", "log_rp", "aor", "b", "log_sigma", "log_rho", "log_jit"]
GP_NAMES = ["log_sigma", "log_rho", "log_jit"]


def load_light_curve():
    d = np.load(LIGHT_CURVE)
    t, f, e, q = d["time"], d["flux"].astype(float), d["flux_err"].astype(float), d["quality"]
    m = (q == 0) & np.isfinite(t) & np.isfinite(f) & (f > 0) & np.isfinite(e)
    return t[m], f[m], e[m] / f[m]


def detrend_and_clip(x, y, ye):
    """Divide by a line through the out-of-transit points; clip 5 robust sigma."""
    oot = np.abs(x) > 0.75 * DUR_D
    if oot.sum() < 20:
        return None
    y = y / np.polyval(np.polyfit(x[oot], y[oot], 1), x)
    r = y - 1.0
    sd = 1.4826 * np.median(np.abs(r[oot] - np.median(r[oot])))
    keep = np.abs(r - np.median(r[oot])) < 5 * sd
    return x[keep], y[keep], ye[keep], keep


def transit_windows():
    t, f, e = load_light_curve()
    ep = np.round((t - T0) / P).astype(int)
    inwin = np.abs(((t - T0) / P - ep) * P) < WINDOW * DUR_D
    blocks = []
    for k in np.unique(ep[inwin]):
        s = inwin & (ep == k)
        if s.sum() < 40:
            continue
        done = detrend_and_clip(t[s] - (T0 + k * P), f[s], e[s])
        if done is None:
            continue
        x, y, ye, keep = done
        if keep.sum() < 40 or (np.abs(x) < 0.5 * DUR_D).sum() < 3:
            continue
        blocks.append((x, y, ye))
    return blocks


def control_windows(t0_alt):
    """The same extraction at an epoch with no transit, real transits removed."""
    t, f, e = load_light_curve()
    away = np.abs(((t - T0) / P + 0.5) % 1.0 - 0.5) * P > 1.5 * DUR_D
    t, f, e = t[away], f[away], e[away]
    ep = np.round((t - t0_alt) / P).astype(int)
    x_all = t - (t0_alt + ep * P)
    inwin = np.abs(x_all) < WINDOW * DUR_D
    blocks = []
    for k in np.unique(ep[inwin]):
        s = inwin & (ep == k)
        if s.sum() < 40:
            continue
        done = detrend_and_clip(x_all[s], f[s], e[s])
        if done is None:
            continue
        x, y, ye, keep = done
        if keep.sum() >= 40:
            blocks.append((x, y, ye))
    return blocks


def one_axis(blocks):
    """All windows on one time axis, 100 d apart, sorted."""
    x = np.concatenate([b[0] + i * 100.0 for i, b in enumerate(blocks)])
    y = np.concatenate([b[1] for b in blocks])
    e = np.concatenate([b[2] for b in blocks])
    order = np.argsort(x)
    return x[order], y[order], e[order], order


def gp_loglike(resid, log_sigma, log_rho, log_jit, x, e):
    try:
        gp = celerite2.GaussianProcess(
            terms.Matern32Term(sigma=10 ** log_sigma, rho=10 ** log_rho), mean=0.0)
        gp.compute(x, diag=e ** 2 + (10 ** log_jit) ** 2, check_sorted=False)
        v = gp.log_likelihood(resid)
    except Exception:
        return -1e30
    return v if np.isfinite(v) else -1e30


def tgp_likelihood(blocks):
    """Transit plus GP. The batman model is built once; only parameters change."""
    x, y, e, order = one_axis(blocks)
    pars = batman.TransitParams()
    pars.t0, pars.per, pars.rp, pars.a = 0.0, P, 0.02, 20.0
    pars.inc, pars.ecc, pars.w = 89.0, 0.0, 90.0
    pars.u, pars.limb_dark = [LD_U1, LD_U2], "quadratic"
    model = batman.TransitModel(pars, np.concatenate([b[0] for b in blocks]),
                                supersample_factor=7, exp_time=CADENCE_D)

    def loglike(p):
        t0, log_rp, aor, b, ls, lr, lj = p
        pars.t0, pars.rp, pars.a = t0, 10 ** log_rp, aor
        pars.inc = np.degrees(np.arccos(np.clip(b / aor, 0, 0.999)))
        try:
            flux = model.light_curve(pars)
        except Exception:
            return -1e30
        return gp_loglike(y - flux[order], ls, lr, lj, x, e)
    return loglike


def gp_likelihood(blocks):
    x, y, e, _ = one_axis(blocks)
    resid = y - 1.0
    return lambda p: gp_loglike(resid, p[0], p[1], p[2], x, e)


def sample(name, names, loglike, bounds):
    lo = np.array([bounds[n][0] for n in names])
    hi = np.array([bounds[n][1] for n in names])
    folder = RUNDIR / name
    if folder.exists():
        shutil.rmtree(folder)
    sampler = ultranest.ReactiveNestedSampler(
        names, loglike, lambda cube: lo + cube * (hi - lo), log_dir=str(folder),
        resume="overwrite", vectorized=False)
    tic = time.time()
    res = sampler.run(min_num_live_points=400, frac_remain=0.3,
                      show_status=False, viz_callback=None)
    print(f"  {name:12s} ln Z = {res['logz']:.2f} +/- {res['logzerr']:.2f}"
          f"   ({time.time() - tic:.0f} s)", flush=True)
    return res


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else OUT
    blocks = transit_windows()
    print(f"{len(blocks)} transit windows, {sum(len(b[0]) for b in blocks)} points")

    r_tgp = sample("tgp_wide", TGP_NAMES, tgp_likelihood(blocks), PRIORS_WIDE)
    r_gp = sample("gp_wide", GP_NAMES, gp_likelihood(blocks), PRIORS_WIDE)
    B = r_tgp["logz"] - r_gp["logz"]
    mean, sd = r_tgp["posterior"]["mean"], r_tgp["posterior"]["stdev"]

    narrow = dict(PRIORS_WIDE)
    for i, n in enumerate(TGP_NAMES[:4]):
        w = max(5 * sd[i], 1e-4)
        lo, hi = mean[i] - w, mean[i] + w
        if n == "b":
            lo, hi = max(0.0, lo), min(0.99, hi)
        if n == "aor":
            lo = max(3.0, lo)
        narrow[n] = (float(lo), float(hi))
    r_narrow = sample("tgp_narrow", TGP_NAMES, tgp_likelihood(blocks), narrow)

    control = control_windows(T0 + CONTROL_PHASE * P)
    print(f"control: {len(control)} windows, {sum(len(b[0]) for b in control)} points")
    r_tgp_c = sample("tgp_control", TGP_NAMES, tgp_likelihood(control), PRIORS_WIDE)
    r_gp_c = sample("gp_control", GP_NAMES, gp_likelihood(control), PRIORS_WIDE)

    results = {
        "data": {"n_windows": len(blocks), "n_points": int(sum(len(b[0]) for b in blocks)),
                 "window_durations": WINDOW, "period_d": P, "t0_bkjd": T0,
                 "limb_darkening": [LD_U1, LD_U2], "supersample": 7,
                 "live_points": 400, "frac_remain": 0.3},
        "wide": {"priors": PRIORS_WIDE,
                 "logz_tgp": float(r_tgp["logz"]), "logzerr_tgp": float(r_tgp["logzerr"]),
                 "logz_gp": float(r_gp["logz"]), "logzerr_gp": float(r_gp["logzerr"]),
                 "B": float(B),
                 "B_err": float(np.hypot(r_tgp["logzerr"], r_gp["logzerr"])),
                 "tgp_mean": {n: float(mean[i]) for i, n in enumerate(TGP_NAMES)},
                 "tgp_sd": {n: float(sd[i]) for i, n in enumerate(TGP_NAMES)},
                 "depth_ppm": float(10 ** (2 * mean[1]) * 1e6)},
        "narrow": {"priors": narrow, "logz_tgp": float(r_narrow["logz"]),
                   "B": float(r_narrow["logz"] - r_gp["logz"])},
        "control": {"epoch_offset_periods": CONTROL_PHASE,
                    "n_windows": len(control),
                    "logz_tgp": float(r_tgp_c["logz"]), "logz_gp": float(r_gp_c["logz"]),
                    "B": float(r_tgp_c["logz"] - r_gp_c["logz"])},
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2))
    print(f"B = {B:.1f} (narrow priors {results['narrow']['B']:.1f}, "
          f"control {results['control']['B']:.1f}); wrote {out}")


if __name__ == "__main__":
    main()
