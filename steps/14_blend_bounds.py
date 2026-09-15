"""Step 14 - how faint a blended star could be, and whether any catalogued star qualifies.

A star that dilutes a deeper eclipse into the observed 575 ppm dip is bounded
two ways.

Brightness: a totally eclipsed star contributes at most its own light, so it
can be at most -2.5 log10(depth) = 8.10 mag fainter than the target.

Shape: a diluted eclipse is intrinsically deeper, and a deeper eclipse has a
longer ingress relative to its duration. With total duration T and ingress tau,
the flat fraction x = (T - 2 tau) / T bounds the radius ratio at k = (1 - x) /
(1 + x), zero impact parameter being the most permissive case (Seager &
Mallen-Ornelas 2003), so the true depth is at most k^2. A trapezoid integrated
over the 29.4-minute long-cadence exposure is fitted to the flattened transits, then refitted
to 600 point-bootstrap resamples; the 95th percentile of the implied magnitude
bound is kept, because a larger permitted depth is the weaker bound.

Every Gaia DR3 source within 25 arcsec is then tested against both bounds and
the centroid: a source more than three standard errors of the combined offset
(step 13) from the target is excluded.

Reproduces "Blend scenarios" in Section 6.
"""

import sys
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives, config, inputs
from koi501.models import GaiaSource
from koi501.report import Table, angsep, read, write

LIGHT_CURVE = config.DATA / "lightcurve" / "kic4951877_pdcsap.npz"
GAIA_RENAME = {"Source": "source_id", "RA_ICRS": "ra", "DE_ICRS": "dec",
               "Gmag": "phot_g_mean_mag"}

USED = ["transit.period_koi_d", "transit.epoch_bkjd", "transit.duration_h", "transit.depth_ppm"]
PERIOD = inputs.get("transit.period_koi_d")
EPOCH = inputs.get("transit.epoch_bkjd")
DURATION = inputs.get("transit.duration_h") / 24.0
DEPTH = inputs.get("transit.depth_ppm") / config.PPM
CADENCE = config.LONG_CADENCE_D
BOUNDS = [list(b) for b in config.SHAPE_BOUNDS]


def flattened_transits() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Points within 0.06 in phase of each transit, divided by a straight line
    through that transit's out-of-transit points; one error per transit, the
    scatter of those points."""
    lc = np.load(LIGHT_CURVE)
    t, f, q = lc["time"], lc["flux"].astype(float), lc["quality"]
    good = (q == 0) & np.isfinite(t) & np.isfinite(f) & (f > 0)
    t, f = t[good], f[good]
    phase = ((t - EPOCH) / PERIOD + 0.5) % 1.0 - 0.5
    epoch = np.round((t - EPOCH) / PERIOD).astype(int)
    near = np.abs(phase) < config.SHAPE_PHASE_WINDOW
    xs, ys, es = [], [], []
    for k in np.unique(epoch[near]):
        window = near & (epoch == k)
        if window.sum() < config.SHAPE_MIN_POINTS:
            continue
        out = window & (np.abs(phase) > config.OOT_DURATIONS * DURATION / PERIOD)
        if out.sum() < config.SHAPE_MIN_OOT:
            continue
        centre = EPOCH + k * PERIOD
        x = t[window] - centre
        line = np.polyfit(t[out] - centre, f[out], 1)
        xs.append(x)
        ys.append(f[window] / np.polyval(line, x))
        es.append(np.full(window.sum(),
                          np.std(f[out] / np.polyval(line, t[out] - centre))))
    return np.concatenate(xs), np.concatenate(ys), np.concatenate(es)


def trapezoid(x, t0, depth, T, tau):
    """Trapezoid of total duration T and ingress tau, averaged over the exposure
    at nine points; tau is comparable to one exposure, so the averaging matters."""
    offsets = np.linspace(-0.5, 0.5, config.SHAPE_SUPERSAMPLE) * CADENCE
    y = np.zeros((len(offsets), len(x)))
    for i, o in enumerate(offsets):
        u = np.abs(x - t0 + o)
        flat = 0.5 * max(T - 2 * tau, 1e-6)  # literal: technical, no zero width
        outer = flat + tau
        v = np.ones_like(u)
        v[u <= flat] = 1 - depth
        ramp = (u > flat) & (u < outer)
        v[ramp] = 1 - depth * (outer - u[ramp]) / max(tau, 1e-6)  # literal: technical, no division by zero
        y[i] = v
    return y.mean(axis=0)


def max_true_depth(T, tau):
    """Largest intrinsic depth a transit of this shape allows (b = 0)."""
    x = np.clip((T - 2 * tau) / T, 0.0, 1.0)
    k = (1 - x) / (1 + x)
    return k ** 2


def percentiles(values, scale=1.0) -> list[float]:
    return [float(np.percentile(values, p) * scale) for p in config.SHAPE_SUMMARY_PERCENTILES]


def main() -> bool:
    # resamples where the covariance cannot be estimated still return a fit
    warnings.simplefilter("ignore", OptimizeWarning)
    x, y, e = flattened_transits()
    best, _ = curve_fit(trapezoid, x, y, p0=[0.0, DEPTH, DURATION, config.SHAPE_INGRESS_GUESS_D],
                        sigma=e, bounds=BOUNDS, maxfev=config.SHAPE_MAXFEV[0])

    rng = np.random.default_rng(config.SHAPE_SEED)
    index = np.arange(len(x))
    fits = []
    for _ in range(config.SHAPE_BOOTSTRAPS):
        s = rng.choice(index, size=len(index), replace=True)
        try:
            p, _ = curve_fit(trapezoid, x[s], y[s], p0=best, sigma=e[s],
                             bounds=BOUNDS, maxfev=config.SHAPE_MAXFEV[1])
        except (RuntimeError, ValueError):
            continue
        fits.append(p)
    fits = np.array(fits)
    depths, Ts, taus = fits[:, 1], fits[:, 2], fits[:, 3]

    # dilution only makes an eclipse shallower, so the observed depth is a floor
    dmax = np.maximum(max_true_depth(Ts, taus), DEPTH)
    dmag = -2.5 * np.log10(DEPTH / dmax)
    shape_bound = float(np.percentile(dmag, config.SHAPE_BOUND_PERCENTILE))
    brightness_bound = float(-2.5 * np.log10(DEPTH))

    offset = read("13_gaia_positional_probability")["section_3_1"]
    gaia = [g for g in archives.validate(
        GaiaSource,
        archives.vizier_cone("I/355/gaiadr3", "Source,RA_ICRS,DE_ICRS,Gmag",
                             config.RA, config.DEC, config.CONE_RADIUS_AS),
        GAIA_RENAME) if g.phot_g_mean_mag is not None]
    rows = sorted(({"source_id": g.source_id, "G": g.phot_g_mean_mag,
                    "sep_as": angsep(g.ra, g.dec, config.RA, config.DEC)} for g in gaia),
                  key=lambda r: r["sep_as"])
    rows = [r for r in rows if r["sep_as"] < config.BLEND_RADIUS]
    target, neighbours = rows[0], rows[1:]
    for r in neighbours:
        r["dmag"] = r["G"] - target["G"]
        r["bright_enough"] = r["dmag"] < brightness_bound
        r["shape_allows"] = r["dmag"] < shape_bound
        r["centroid_sigma"] = r["sep_as"] / offset["offset_err_as"]
        r["survives"] = (r["bright_enough"] and r["shape_allows"]
                         and r["centroid_sigma"] < config.CENTROID_EXCLUSION_SIGMA)

    payload = {
        "measured_inputs": inputs.record(USED),
        "light_curve": {"file": LIGHT_CURVE.relative_to(config.ROOT).as_posix(),
                        "n_points": int(len(x)), "period_d": PERIOD, "epoch_bkjd": EPOCH,
                        "exposure_d": CADENCE},
        "best_fit": {"epoch_offset_min": float(best[0] * 1440), "depth_ppm": float(best[1] * config.PPM),
                     "duration_h": float(best[2] * 24), "ingress_min": float(best[3] * 1440)},
        "bootstrap": {"requested": config.SHAPE_BOOTSTRAPS, "converged": int(len(fits)),
                      "seed": config.SHAPE_SEED,
                      "depth_ppm_p16_p50_p84": percentiles(depths, config.PPM),
                      "duration_h_p16_p50_p84": percentiles(Ts, 24),
                      "ingress_min_p16_p50_p84": percentiles(taus, 1440),
                      "T_over_tau_p16_p50_p84": percentiles(Ts / taus),
                      # unscaled, for fpp/own_fpp.py
                      "duration_d_median": float(np.median(Ts)),
                      "ingress_d_p16_p50_p84": percentiles(taus)},
        "max_true_depth_p95_ppm": float(np.percentile(dmax, config.SHAPE_BOUND_PERCENTILE) * config.PPM),
        "dmag_shape_bound": shape_bound,
        "dmag_shape_bound_median": float(np.median(dmag)),
        "dmag_brightness_bound": brightness_bound,
        "centroid_error_as": offset["offset_err_as"],
        "target": target,
        "blend_radius_as": config.BLEND_RADIUS,
        "n_neighbours_within_25as": len(neighbours),
        "n_survivors": sum(r["survives"] for r in neighbours),
        "neighbours": neighbours,
    }
    write("14_blend_bounds", payload)

    boot = payload["bootstrap"]
    table = Table("How faint a blended star could be", "Section 6, Blend scenarios")
    table("points near transit, after flattening", len(x))
    table("best-fit ingress (min) / duration (h)",
          f"{best[3] * 1440:.2f} / {best[2] * 24:.3f}")
    table("bootstrap fits converged", f"{len(fits)} of {config.SHAPE_BOOTSTRAPS}")
    table("  ingress, median [16th, 84th] (min)",
          "{1:.2f} [{0:.2f}, {2:.2f}]".format(*boot["ingress_min_p16_p50_p84"]))
    table("brightness bound (mag fainter than target)", f"{brightness_bound:.2f}")
    table("shape bound, 95th percentile (mag)", f"{shape_bound:.2f}")
    table(f"Gaia DR3 sources within {config.BLEND_RADIUS:g} arcsec besides the target", len(neighbours))
    table("  brighter than the shape bound", sum(r["shape_allows"] for r in neighbours))
    table(f"  within {config.CENTROID_EXCLUSION_SIGMA:g} standard errors of the centroid",
          sum(r["centroid_sigma"] < config.CENTROID_EXCLUSION_SIGMA for r in neighbours))
    table("  surviving all three", payload["n_survivors"])
    table.show("14_blend_bounds")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
