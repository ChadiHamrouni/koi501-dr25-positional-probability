"""The light-curve rows of Table 4 marked "here", and the transit times.

On the long-cadence PDCSAP light curve (data/lightcurve/, quality flag 0), each
contiguous segment is divided by a polynomial fitted to its out-of-transit
points with iterative clipping. A box depth is the mean baseline minus the mean
in-transit flux. Then:

  wrong epochs     the box depth at 20,000 epochs drawn away from the transit;
                   the deepest is the "deepest false dip"
  single transits  each transit's own box depth; how many are deeper than zero
                   by more than their error
  secondary        the box depth at phase 0.5
  odd and even     the depth of odd- and even-numbered transits separately
  halves           the depth in the first and second half of the mission
  M0V host         the minimum eccentricity (Dawson & Johnson 2012) an M0V host
                   would need, from the stellar density a circular orbit implies
                   (the free-density transit fit, fit/transit_fit.py)
  transit times    each transit fitted for its centre alone, with the trapezoid
                   of step 14 held fixed; the true precision measured by
                   injecting the template into out-of-transit stretches; the
                   odd and even times compared
  crowding         CROWDSAP from the light-curve headers, per quarter (MAST)

    pip install -r fit/requirements.txt
    python fit/lightcurve_checks.py          # a few minutes

Writes results/lightcurve_checks.json.
"""

from __future__ import annotations

import csv
import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

from koi501 import config, inputs  # noqa: E402

LIGHT_CURVE = ROOT / "data" / "lightcurve" / "kic4951877_pdcsap.npz"
TRANSIT_FIT = ROOT / "data" / "transit_fit" / "transit_fit.json"
DV_QUARTERS = ROOT / "data" / "dv" / "dvr_quarterly_centroids.csv"
OUT = ROOT / "results" / "lightcurve_checks.json"

USED = ["transit.period_koi_d", "transit.epoch_bkjd", "transit.duration_h",
        "shape.ingress_d_p50", "shape.depth_ppm_p50"]
P = inputs.get("transit.period_koi_d")
T0 = inputs.get("transit.epoch_bkjd")
DUR_D = inputs.get("transit.duration_h") / 24.0
CADENCE_D = config.LONG_CADENCE_D
MIN_PER_D = 1440.0


def load():
    d = np.load(LIGHT_CURVE)
    t, f, e, q = d["time"], d["flux"].astype(float), d["flux_err"].astype(float), d["quality"]
    m = (q == 0) & np.isfinite(t) & np.isfinite(f) & (f > 0)
    return t[m], f[m], e[m]


def detrend(t, f, e):
    """Each segment divided by a clipped polynomial through its out-of-transit points."""
    seg = np.zeros(len(t), int)
    for i, g in enumerate(np.where(np.diff(t) > config.LC_SEGMENT_GAP_D)[0]):
        seg[g + 1:] = i + 1
    out = np.empty_like(f)
    ph = ((t - T0) / P + 0.5) % 1.0 - 0.5
    oot = np.abs(ph) > config.OOT_DURATIONS * DUR_D / P
    short, long_ = config.LC_DETREND_DEGREE
    for s in np.unique(seg):
        m = seg == s
        if m.sum() < config.LC_SEGMENT_MIN_POINTS:
            out[m] = np.nan
            continue
        x = t[m] - t[m].mean()
        x = x / max(x.std(), 1e-9)  # literal: technical, no division by zero
        deg = short if m.sum() < config.LC_DETREND_LONG_POINTS else long_
        X = np.vander(x, deg + 1)
        keep = oot[m].copy()
        for _ in range(config.LC_DETREND_ITERATIONS):
            if keep.sum() < deg + config.LC_DETREND_SPARE:
                break
            b, *_ = np.linalg.lstsq(X[keep], f[m][keep], rcond=None)
            r = f[m] - X @ b
            sd = config.MAD_TO_SIGMA * np.median(np.abs(r[keep] - np.median(r[keep])))
            keep = oot[m] & (np.abs(r - np.median(r[keep])) < config.LC_DETREND_CLIP_SIGMA * sd)
        if keep.sum() < deg + config.LC_DETREND_SPARE:
            out[m] = np.nan
            continue
        b, *_ = np.linalg.lstsq(X[keep], f[m][keep], rcond=None)
        out[m] = f[m] / (X @ b)
    ok = np.isfinite(out)
    return t[ok], out[ok], (e / f)[ok]


def box_depth(t, f, t0):
    """Baseline minus in-transit mean, ppm, and its error."""
    ph = ((t - t0) / P + 0.5) % 1.0 - 0.5
    hw = 0.5 * DUR_D / P
    i = np.abs(ph) < config.LC_BOX_IN * hw
    o = (np.abs(ph) > config.LC_BOX_OUT * hw) & (np.abs(ph) < config.LC_BOX_OUT_PHASE)
    if i.sum() < config.LC_BOX_MIN_IN or o.sum() < config.LC_BOX_MIN_OUT:
        return np.nan, np.nan
    d = (np.mean(f[o]) - np.mean(f[i])) * config.PPM
    s = np.hypot(np.std(f[i], ddof=1) / np.sqrt(i.sum()),
                 np.std(f[o], ddof=1) / np.sqrt(o.sum())) * config.PPM
    return float(d), float(s)


def compare(a, b):
    return {"first_ppm": a, "second_ppm": b,
            "sigma": float(abs(a[0] - b[0]) / np.hypot(a[1], b[1]))}


def depth_checks(rng):
    t, f, _ = detrend(*load())
    ph = ((t - T0) / P + 0.5) % 1.0 - 0.5
    hw = 0.5 * DUR_D / P
    ep = np.round((t - T0) / P).astype(int)
    depth = box_depth(t, f, T0)

    lo, hi = config.LC_NULL_PHASES
    nulls = np.array([box_depth(t, f, T0 + rng.uniform(lo, hi) * P)[0]
                      for _ in range(config.LC_NULL_EPOCHS)])
    nulls = nulls[np.isfinite(nulls)]

    single = []
    for k in np.unique(ep[np.abs(ph) < hw]):
        m = ep == k
        loc = m & (np.abs(ph) < config.LC_BOX_IN * hw)
        base = m & (np.abs(ph) > config.LC_BOX_OUT * hw) & (np.abs(ph) < config.LC_SINGLE_OUT_PHASE)
        if loc.sum() < config.LC_SINGLE_MIN_IN or base.sum() < config.LC_SINGLE_MIN_OUT:
            continue
        dd = (np.mean(f[base]) - np.mean(f[loc])) * config.PPM
        ss = np.std(f[base], ddof=1) * config.PPM * np.sqrt(1 / loc.sum() + 1 / base.sum())
        single.append((int(k), float(dd), float(ss)))
    dep = np.array([s[1] for s in single])
    err = np.array([s[2] for s in single])

    odd, even = ep % 2 == 1, ep % 2 == 0
    tmid = 0.5 * (t.min() + t.max())
    return {
        "n_points": int(len(t)),
        "depth_ppm": depth,
        "wrong_epochs": {"n": int(len(nulls)), "deepest_ppm": float(nulls.max()),
                         "mean_ppm": float(nulls.mean()), "sd_ppm": float(nulls.std()),
                         "n_as_deep_as_transit": int((nulls >= depth[0]).sum())},
        "single_transits": {"n": len(single), "n_deeper_than_zero_1sigma": int((dep > err).sum()),
                            "depths": [{"epoch": k, "depth_ppm": d, "err_ppm": s} for k, d, s in single]},
        "secondary_at_phase_half_ppm": box_depth(t, f, T0 + config.LC_SECONDARY_PHASE * P),
        "odd_even": compare(box_depth(t[odd], f[odd], T0), box_depth(t[even], f[even], T0)),
        "halves": {"split_bkjd": float(tmid),
                   **compare(box_depth(t[t < tmid], f[t < tmid], T0),
                             box_depth(t[t >= tmid], f[t >= tmid], T0))},
    }


def m0v_eccentricity() -> dict:
    """e_min = |1 - r^(2/3)| / (1 + r^(2/3)), r the circular-orbit density over the host's."""
    fit = json.loads(TRANSIT_FIT.read_text(encoding="utf8"))["runs"]["free_rho"]
    aor = fit["a_over_rstar"]["median"]
    rho_circ = (3 * np.pi * aor ** 3 / (config.G_SI * (P * config.DAY_S) ** 2)
                / (config.RHO_SUN_CGS * 1000))
    mass, radius = config.M0V_MASS_RADIUS
    r = rho_circ / (mass / radius ** 3)
    return {"a_over_rstar_free_density": aor, "rho_circ_solar": float(rho_circ),
            "m0v_density_solar": mass / radius ** 3,
            "e_min": float(abs(1 - r ** (2 / 3)) / (1 + r ** (2 / 3)))}


def template(x, tc, depth, tau):
    off = np.linspace(-0.5, 0.5, config.SHAPE_SUPERSAMPLE) * CADENCE_D
    acc = np.zeros((len(off), len(x)))
    flat = 0.5 * (DUR_D - 2 * tau)
    outer = flat + tau
    for i, o in enumerate(off):
        u = np.abs(x - tc + o)
        v = np.ones_like(u)
        v[u <= flat] = 1 - depth
        ramp = (u > flat) & (u < outer)
        v[ramp] = 1 - depth * (outer - u[ramp]) / tau
        acc[i] = v
    return acc.mean(axis=0)


def centre(x, y, sd, depth, tau):
    """Best centre from a chi-squared grid and a parabola through its minimum."""
    grid = np.linspace(-config.TIMING_GRID_D, config.TIMING_GRID_D, config.TIMING_GRID_SIZE)
    chi2 = np.array([np.sum((y - template(x, g, depth, tau)) ** 2) / sd ** 2 for g in grid])
    j, h = int(np.argmin(chi2)), config.TIMING_PARABOLA
    if j == 0 or j == len(grid) - 1:
        return None
    a, b, _ = np.polyfit(grid[j - h:j + h + 1], chi2[j - h:j + h + 1], 2)
    if a <= 0:
        return None
    return -b / (2 * a), 1.0 / np.sqrt(a), chi2[j]


def flatten(x, y, oot):
    y = y / np.polyval(np.polyfit(x[oot], y[oot], 1), x)
    return y, np.std(y[oot], ddof=1)


def transit_times(rng) -> dict:
    t, f, _ = load()
    depth = inputs.get("shape.depth_ppm_p50") / config.PPM
    tau = inputs.get("shape.ingress_d_p50")
    ep_all = np.round((t - T0) / P).astype(int)
    ph = (t - T0) / P - ep_all

    rows = []
    for k in np.unique(ep_all[np.abs(ph) < config.TIMING_SEARCH_PHASE]):
        s = (ep_all == k) & (np.abs(ph * P) < config.TIMING_WINDOW_D)
        if s.sum() < config.TIMING_MIN_POINTS:
            continue
        x, y = t[s] - (T0 + k * P), f[s]
        oot = np.abs(x) > config.OOT_DURATIONS * DUR_D
        if (oot.sum() < config.TIMING_MIN_OOT
                or (np.abs(x) < config.TIMING_IN_FRACTION * DUR_D).sum() < config.TIMING_MIN_IN):
            continue
        y, sd = flatten(x, y, oot)
        if not np.isfinite(sd) or sd <= 0:
            continue
        found = centre(x, y, sd, depth, tau)
        if found is None:
            continue
        tc, err, _ = found
        if abs(tc) > config.TIMING_MAX_OFFSET_D or not np.isfinite(err) or err > config.TIMING_MAX_ERROR_D:
            continue
        rows.append((int(k), float(tc), float(err)))

    ep = np.array([r[0] for r in rows], float)
    oc, er = np.array([r[1] for r in rows]), np.array([r[2] for r in rows])
    w = 1.0 / er ** 2
    A = np.vstack([np.ones_like(ep), ep]).T
    beta = np.linalg.inv(A.T @ (A * w[:, None])) @ (A.T @ (w * oc))
    resid = oc - A @ beta
    dof = len(rows) - 2
    scale = max(1.0, np.sqrt(float(np.sum(w * resid ** 2)) / dof))

    odd, even = ep.astype(int) % 2 == 1, ep.astype(int) % 2 == 0
    mo = np.sum(w[odd] * resid[odd]) / np.sum(w[odd])
    me = np.sum(w[even] * resid[even]) / np.sum(w[even])
    so, se = scale / np.sqrt(np.sum(w[odd])), scale / np.sqrt(np.sum(w[even]))

    recovered, tries = [], 0
    while len(recovered) < config.TIMING_INJECTIONS and tries < config.TIMING_TRIES:
        tries += 1
        c0 = rng.uniform(t.min() + 1.0, t.max() - 1.0)
        if np.abs(((c0 - T0) / P + 0.5) % 1.0 - 0.5) * P < config.TIMING_CLEARANCE * DUR_D:
            continue
        s2 = np.abs(t - c0) < config.TIMING_WINDOW_D
        if s2.sum() < config.TIMING_MIN_POINTS:
            continue
        x = t[s2] - c0
        true_off = rng.uniform(-config.TIMING_INJECTED_OFFSET_D, config.TIMING_INJECTED_OFFSET_D)
        y = f[s2] * template(x, true_off, depth, tau)
        oot = np.abs(x - true_off) > config.OOT_DURATIONS * DUR_D
        if oot.sum() < config.TIMING_MIN_OOT:
            continue
        y, sd = flatten(x, y, oot)
        if not np.isfinite(sd) or sd <= 0:
            continue
        found = centre(x, y, sd, depth, tau)
        if found is not None:
            recovered.append(found[0] - true_off)
    recovered = np.array(recovered)
    precision = float(np.std(recovered[np.abs(recovered) < config.TIMING_RECOVERY_CLIP_D], ddof=1))

    return {"n_transits": len(rows), "template": {"depth_ppm": depth * config.PPM, "ingress_min": tau * MIN_PER_D},
            "scatter_min": float(np.std(resid) * MIN_PER_D),
            "formal_precision_min": float(np.median(er) * MIN_PER_D),
            "n_injections": int(len(recovered)),
            "true_precision_min": precision * MIN_PER_D,
            "chi2_red_true_precision": float(np.sum((resid / precision) ** 2) / dof),
            "odd_even_sigma": float(abs(mo - me) / np.hypot(so, se)),
            "times": [{"epoch": k, "bkjd": T0 + k * P + tc, "err_d": e} for k, tc, e in rows]}


def crowding() -> dict:
    """CROWDSAP per quarter from the light-curve file headers (MAST)."""
    import lightkurve as lk
    files = lk.search_lightcurve(f"KIC {config.KEPID}", mission="Kepler",
                                 cadence="long").download_all()
    values = {int(lc.meta["QUARTER"]): float(lc.meta["CROWDSAP"]) for lc in files
              if lc.meta.get("CROWDSAP") is not None}
    # the quarters with a difference image, as the DV report lists them
    with DV_QUARTERS.open(encoding="utf8") as handle:
        imaged = {int(row["quarter"]) for row in csv.DictReader(handle)}
    in_images = [values[q] for q in sorted(imaged & set(values))]
    return {"per_quarter": {str(q): values[q] for q in sorted(values)},
            "range": [min(values.values()), max(values.values())],
            "difference_image_quarters": sorted(imaged),
            "range_difference_image_quarters": [min(in_images), max(in_images)]}


def main() -> None:
    result = {"measured_inputs": inputs.record(USED)}
    result.update(depth_checks(np.random.default_rng(config.LC_CHECK_SEED)))
    result["m0v_host"] = m0v_eccentricity()
    result["transit_times"] = transit_times(np.random.default_rng(config.LC_CHECK_SEED))
    result["crowdsap"] = crowding()
    OUT.write_text(json.dumps(result, indent=2), encoding="utf8")

    r, tt = result, result["transit_times"]
    print(f"  depth {r['depth_ppm'][0]:.0f} +/- {r['depth_ppm'][1]:.1f} ppm; deepest of "
          f"{r['wrong_epochs']['n']} wrong epochs {r['wrong_epochs']['deepest_ppm']:.0f} ppm")
    print(f"  transits deeper than zero: {r['single_transits']['n_deeper_than_zero_1sigma']} "
          f"of {r['single_transits']['n']}")
    print(f"  phase 0.5: {r['secondary_at_phase_half_ppm'][0]:+.1f} +/- {r['secondary_at_phase_half_ppm'][1]:.1f} ppm")
    for key in ("odd_even", "halves"):
        a, b = r[key]["first_ppm"], r[key]["second_ppm"]
        print(f"  {key}: {a[0]:.0f} +/- {a[1]:.0f} / {b[0]:.0f} +/- {b[1]:.0f} ppm ({r[key]['sigma']:.2f} sigma)")
    print(f"  M0V host e_min {r['m0v_host']['e_min']:.3f}")
    print(f"  times: {tt['n_transits']} transits, scatter {tt['scatter_min']:.1f} min vs precision "
          f"{tt['true_precision_min']:.1f} min, chi2_red {tt['chi2_red_true_precision']:.2f}, "
          f"odd-even {tt['odd_even_sigma']:.2f} sigma")
    print(f"  CROWDSAP {r['crowdsap']['range'][0]:.3f} to {r['crowdsap']['range'][1]:.3f}")
    print(f"  wrote {OUT.relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()
