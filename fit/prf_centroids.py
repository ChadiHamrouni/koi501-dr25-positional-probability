"""Where the light changed: a PRF fit to the difference images, and its systematic.

For every quarterly difference image of KOI-501.01 (data/diffimages/kic4951877.npz)
the Kepler pixel response function (Bryson et al. 2010; the mission's 2011
calibration files) is fitted to the image: position, flux and a constant. The
fitted position is turned into an offset from the target's catalogue position
through the image's WCS, and the offsets of all good panels are combined by
their median, east and north separately, with a robust scatter over the root of
the panel count as the error on each axis.

The same fit is run on 243 control objects (data/diffimages/controls.csv), a
seeded random draw (seed 20260829) from the Kepler objects whose pixel data the
project had downloaded, restricted to DR25 difference-image offsets with errors
below 2 arcsec. Those that DR25 places on their target at better than one
standard error calibrate the systematic: their fitted offsets should follow a
Rayleigh distribution whose median is sqrt(2 ln 2) times the per-axis error. The
per-axis error that median implies, less the median formal per-axis error in
quadrature, is the per-axis systematic used throughout the paper.

    pip install -r fit/requirements.txt
    python fit/prf_centroids.py          # a few minutes; PRF files cached on first run

Writes results/prf_centroids.json and results/prf_centroids_controls.csv. The
difference images were built from the Kepler target pixel files by
fit/diff_images.py.
"""

from __future__ import annotations

import csv
import json
import math
import sys
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import scipy.interpolate
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

from koi501 import config  # noqa: E402

DIFFIMAGES = ROOT / "data" / "diffimages"
CONTROLS = DIFFIMAGES / "controls.csv"
PRF_DIR = ROOT / "results" / ".cache" / "kepler_prf"
PRF_URL = "https://archive.stsci.edu/missions/kepler/fpc/prf/"
OUT = ROOT / "results" / "prf_centroids.json"
OUT_CONTROLS = ROOT / "results" / "prf_centroids_controls.csv"
#: median of a Rayleigh distribution, in units of its per-axis sigma
RAYLEIGH_MEDIAN = math.sqrt(2.0 * math.log(2.0))

_calibration: dict = {}


def calibration(module: int, output: int):
    """The five reference PRF images of one detector output and their positions."""
    key = (module, output)
    if key not in _calibration:
        name = f"kplr{module:02d}.{output}_2011265_prf.fits"
        local = PRF_DIR / name
        if not local.exists():
            PRF_DIR.mkdir(parents=True, exist_ok=True)
            local.write_bytes(urllib.request.urlopen(PRF_URL + name, timeout=600).read())  # literal: technical
        with fits.open(local) as h:
            imgs = [np.array(x.data, float) for x in h[1:]]
            head = [x.header for x in h[1:]]
        _calibration[key] = (imgs, *(np.array([float(x[k]) for x in head])
                                     for k in ("CRVAL1P", "CRVAL2P", "CDELT1P", "CDELT2P")))
    return _calibration[key]


class KeplerPRF:
    """The PRF for a stamp of `shape` whose corner pixel is (column, row).

    The five reference images are blended by inverse distance to the stamp
    centre, normalised, splined, and evaluated on the stamp's pixel grid.
    """

    def __init__(self, shape, module, output, column, row):
        imgs, crval1p, crval2p, cdelt1p, cdelt2p = calibration(module, output)
        rowdim, coldim = shape
        ref_col, ref_row = column + 0.5 * coldim, row + 0.5 * rowdim
        prf = np.zeros_like(imgs[0], dtype=float)
        for i, img in enumerate(imgs):
            prf += img / max(np.hypot(ref_col - crval1p[i], ref_row - crval2p[i]),
                             1e-6)  # literal: technical, avoids division by zero
        prf /= np.nansum(prf) * cdelt1p[0] * cdelt2p[0]
        pc = (np.arange(0.5, prf.shape[1] + 0.5) - prf.shape[1] / 2) * cdelt1p[0]
        pr = (np.arange(0.5, prf.shape[0] + 0.5) - prf.shape[0] / 2) * cdelt2p[0]
        self.col_coord = np.arange(column + 1, column + coldim + 1)
        self.row_coord = np.arange(row + 1, row + rowdim + 1)
        self.interpolate = scipy.interpolate.RectBivariateSpline(pr, pc, prf)

    def evaluate(self, center_col, center_row, flux=1.0):
        return float(flux) * np.asarray(self.interpolate(self.row_coord - center_row,
                                                         self.col_coord - center_col))


def index_to_prf(col0, row0, ix, iy):
    return col0 + ix + config.PRF_INDEX_OFFSET, row0 + iy + config.PRF_INDEX_OFFSET


def fit_prf(image, prf, col0, row0):
    """Least-squares (column, row, flux, constant), seeded at the brightest pixel.

    Returns the position, the flux and the correlation of model with data.
    """
    ny, nx = image.shape
    data = np.nan_to_num(image, nan=0.0)
    iy, ix = np.unravel_index(int(np.argmax(data)), data.shape)
    seed_col, seed_row = index_to_prf(col0, row0, ix, iy)

    def model(p):
        return p[3] + p[2] * prf.evaluate(p[0], p[1], 1.0)

    def cost(p):
        return float(np.sum((model(p) - data) ** 2))

    amp = float(np.nanmax(data) - np.nanmedian(data))
    base = float(np.nanmedian(data))
    lo_c, lo_r = index_to_prf(col0, row0, 0, 0)
    hi_c, hi_r = index_to_prf(col0, row0, nx - 1, ny - 1)
    r = minimize(cost, np.array([seed_col, seed_row, max(amp, 1e-9), base]),  # literal: technical
                 method="L-BFGS-B",
                 bounds=[(lo_c, hi_c), (lo_r, hi_r), (0, None), (None, None)])
    with np.errstate(invalid="ignore"):
        q = float(np.corrcoef(model(r.x).ravel(), data.ravel())[0, 1])
    return float(r.x[0]), float(r.x[1]), float(r.x[2]), q if np.isfinite(q) else 0.0


def panels(path: Path):
    d = np.load(path, allow_pickle=True)
    for key in [k for k in d.files if k.startswith("meta_")]:
        suffix = key[len("meta_"):]
        err = d["err_" + suffix] if "err_" + suffix in d.files else None
        yield (np.asarray(d["d_plain_" + suffix], float),
               None if err is None else np.asarray(err, float),
               json.loads(str(d[key])))


def measure(path: Path, prf_cache: dict) -> dict | None:
    """Offset of the PRF-fitted source from the catalogue position, arcsec."""
    east, north, snrs, quality, skipped = [], [], [], [], 0
    for image, err, meta in panels(path):
        wcs = WCS(fits.Header.fromstring(meta["wcs"]), relax=True)
        a = np.nan_to_num(image)
        noise = np.nanmedian(err) if err is not None else None
        if noise is None or not np.isfinite(noise) or noise <= 0:
            noise = config.MAD_TO_SIGMA * np.median(np.abs(a - np.median(a))) or 1.0  # literal: technical
        if a.max() / noise < config.PRF_MIN_PANEL_SNR:
            skipped += 1
            continue
        key = (meta["module"], meta["output"], meta["column"], meta["row"], a.shape)
        if key not in prf_cache:
            prf_cache[key] = KeplerPRF(a.shape, *key[:4])
        col, row, _, q = fit_prf(a, prf_cache[key], meta["column"], meta["row"])
        if not np.isfinite(q) or q < config.PRF_MIN_QUALITY:
            skipped += 1
            continue
        ra, dec = wcs.all_pix2world([[col - meta["column"] - config.PRF_INDEX_OFFSET,
                                      row - meta["row"] - config.PRF_INDEX_OFFSET]], 0)[0]
        if not (np.isfinite(ra) and np.isfinite(dec)):
            skipped += 1
            continue
        east.append((float(ra) - meta["ra"]) * np.cos(np.radians(meta["dec"])) * 3600.0)
        north.append((float(dec) - meta["dec"]) * 3600.0)
        snrs.append(a.max() / noise)
        quality.append(q)

    if len(east) < config.PRF_MIN_PANELS:
        return None
    de, dn = np.array(east), np.array(north)
    e, n = float(np.median(de)), float(np.median(dn))
    se = float(config.MAD_TO_SIGMA * np.median(np.abs(de - e))) / np.sqrt(len(de))
    sn = float(config.MAD_TO_SIGMA * np.median(np.abs(dn - n))) / np.sqrt(len(dn))
    off = float(np.hypot(e, n))
    sig = float(np.hypot(e * se, n * sn) / off) if off > 0 else float(np.hypot(se, sn))
    return {"offset_as": off, "sigma_as": sig, "east_as": e, "north_as": n,
            "east_err_as": se, "north_err_as": sn, "n_panels": len(de),
            "n_skipped": skipped, "median_panel_snr": float(np.median(snrs)),
            "median_panel_quality": float(np.median(quality))}


def main() -> None:
    cache: dict = {}
    target = measure(DIFFIMAGES / f"kic{config.KEPID}.npz", cache)

    controls = []
    listed = list(csv.DictReader(CONTROLS.open(encoding="utf8")))
    for i, r in enumerate(listed, 1):
        m = measure(DIFFIMAGES / f"kic{r['kepid']}.npz", cache)
        if m is not None:
            controls.append({"kepid": int(r["kepid"]), **m,
                             "dr25_offset_as": float(r["dr25_offset_as"]),
                             "dr25_offset_err_as": float(r["dr25_offset_err_as"])})
        if i % 50 == 0:  # literal: display, progress
            print(f"  controls {i} / {len(listed)}", flush=True)

    on_target = [c for c in controls
                 if c["dr25_offset_as"] / c["dr25_offset_err_as"] < config.PRF_ON_TARGET_SIGMA]
    median_offset = float(np.median([c["offset_as"] for c in on_target]))
    total_per_axis = median_offset / RAYLEIGH_MEDIAN
    formal_per_axis = float(np.median([c["sigma_as"] for c in on_target])) / math.sqrt(2.0)
    systematic = math.sqrt(total_per_axis ** 2 - formal_per_axis ** 2)

    bright = [c for c in controls if c["median_panel_snr"] > config.PRF_COMPARISON_SNR]
    diff = np.array([c["offset_as"] - c["dr25_offset_as"] for c in bright])
    bias = float(np.median(diff))
    offsets = np.array([c["offset_as"] for c in controls])

    payload = {
        "target": {"kepid": config.KEPID, **target,
                   "controls_with_smaller_offset_percent":
                       float(100.0 * np.mean(offsets < target["offset_as"]))},
        "controls": {
            "n_listed": len(listed), "n_measured": len(controls),
            "n_on_target": len(on_target),
            "on_target_median_offset_as": median_offset,
            "total_per_axis_as": total_per_axis,
            "formal_per_axis_as": formal_per_axis,
            "n_bright": len(bright),
            "correlation_with_dr25": float(np.corrcoef(
                [c["offset_as"] for c in bright], [c["dr25_offset_as"] for c in bright])[0, 1]),
            "bias_minus_dr25_as": bias,
            "robust_scatter_minus_dr25_as": float(config.MAD_TO_SIGMA * np.median(np.abs(diff - bias))),
        },
        "systematic_per_axis_as": systematic,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf8")
    with OUT_CONTROLS.open("w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=list(controls[0]))
        w.writeheader()
        w.writerows(controls)

    t = payload["target"]
    print(f"\n  KOI-501.01: {t['n_panels']} panels, east {t['east_as']:+.3f}, north {t['north_as']:+.3f} arcsec")
    print(f"  offset {t['offset_as']:.3f} +/- {t['sigma_as']:.3f} arcsec")
    print(f"  controls measured {len(controls)} of {len(listed)}; on target {len(on_target)}")
    print(f"  on-target median offset {median_offset:.3f} -> per axis {total_per_axis:.3f}; "
          f"formal {formal_per_axis:.3f}")
    print(f"  systematic per axis {systematic:.4f} arcsec")
    print(f"  wrote {OUT.relative_to(ROOT).as_posix()} and {OUT_CONTROLS.relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()
