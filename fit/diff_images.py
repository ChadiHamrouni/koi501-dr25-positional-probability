"""The difference images: built from the Kepler target pixel files.

For every long-cadence quarter of a star and every DR25 detection (TCE) on it:
cadences within half a duration of mid-transit are "in"; cadences between 0.8
and 3 durations are "out". For each single transit with enough of both, the
mean out-of-transit frame of that transit minus its mean in-transit frame is
formed and weighted by the number of in-transit cadences; the weighted mean over
transits is the quarter's difference image. The out-of-transit mean over the
quarter is the direct image. Each pixel's error is a robust scatter of the
cadence-to-cadence differences, scaled by the numbers of in and out cadences.

    python fit/diff_images.py 4951877              # rebuild one star into results/.cache/diffimages/
    python fit/diff_images.py 4951877 --check      # and compare with data/diffimages/

The committed files (KOI-501.01 and the 243 controls of fit/prf_centroids.py)
were built with this procedure; the controls keep only the difference image,
its error and the metadata.
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

from koi501 import archives, config  # noqa: E402

COMMITTED = ROOT / "data" / "diffimages"
OUT = ROOT / "results" / ".cache" / "diffimages"


def mean_frame(cube, select):
    return np.nan_to_num(np.nanmean(cube[select], axis=0), nan=0.0)


def difference(cube, phase, epoch, inside, outside):
    """In-transit-count weighted mean over single transits of (out - in)."""
    total, weight, used = None, 0.0, 0
    for k in np.unique(epoch[inside]):
        in_k = inside & (epoch == k)
        if in_k.sum() < config.DIFF_MIN_IN_TRANSIT:
            continue
        out_k = np.where(outside & (epoch == k))[0]
        if out_k.size < config.DIFF_MIN_OUT_TRANSIT:
            continue
        w = float(in_k.sum())
        local = (mean_frame(cube, out_k) - mean_frame(cube, in_k)) * w
        total = local if total is None else total + local
        weight += w
        used += 1
    return (None, 0) if total is None else (total / weight, used)


def build(kepid: int) -> dict:
    import lightkurve as lk

    tces = archives.nasa_tap(
        "select tce_plnt_num,tce_period,tce_time0bk,tce_duration from q1_q17_dr25_tce "
        f"where kepid={kepid}")
    search = lk.search_targetpixelfile(f"KIC {kepid}", mission="Kepler", cadence="long")
    tpfs = [f for f in search[:config.DIFF_MAX_TPF].download_all()]
    payload = {}
    for tce in tces:
        period, t0 = float(tce["tce_period"]), float(tce["tce_time0bk"])
        duration = float(tce["tce_duration"]) / 24.0
        for q_index, tpf in enumerate(tpfs):
            t = np.asarray(tpf.time.value, float)
            cube = np.asarray(tpf.flux.value, float)
            usable = np.isfinite(cube).any(axis=0)
            if usable.sum() < config.DIFF_MIN_PIXELS:
                continue
            ok = np.isfinite(t) & (np.isfinite(cube[:, usable]).mean(axis=1) > config.DIFF_FINITE_FRACTION)
            t, cube = t[ok], cube[ok]
            if t.size < config.DIFF_MIN_CADENCES:
                continue
            phase = ((t - t0 + 0.5 * period) % period) - 0.5 * period
            a = np.abs(phase) / duration
            epoch = np.round((t - t0) / period).astype(int)
            lo, hi = config.DIFF_OUT_RANGE
            inside, outside = a <= config.DIFF_IN_FRACTION, (a > lo) & (a <= hi)
            if inside.sum() < config.DIFF_MIN_IN or outside.sum() < config.DIFF_MIN_OUT:
                continue
            diff, used = difference(cube, phase, epoch, inside, outside)
            if diff is None:
                continue
            steps = np.diff(np.nan_to_num(cube, nan=0.0), axis=0)
            scatter = config.MAD_TO_SIGMA * np.nanmedian(np.abs(steps), axis=0) / np.sqrt(2.0)
            err = scatter * np.sqrt(1.0 / inside.sum() + 1.0 / outside.sum())
            key = f"{tce['tce_plnt_num']}_{q_index}"
            header = tpf.hdu[0].header
            payload["direct_" + key] = mean_frame(cube, outside).astype(np.float32)
            payload["err_" + key] = np.maximum(err, config.DIFF_ERROR_FLOOR).astype(np.float32)
            payload["d_plain_" + key] = diff.astype(np.float32)
            payload["meta_" + key] = np.array(json.dumps(dict(
                tce=int(tce["tce_plnt_num"]), q=q_index, module=int(header["MODULE"]),
                output=int(header["OUTPUT"]), column=int(tpf.column), row=int(tpf.row),
                ra=float(tpf.ra), dec=float(tpf.dec), wcs=tpf.wcs.to_header().tostring(),
                n_transits=used)))
    return payload


def main() -> int:
    kepid = int(sys.argv[1])
    payload = build(kepid)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / f"kic{kepid}.npz", **payload)
    print(f"  kic{kepid}: {sum(k.startswith('meta_') for k in payload)} difference images")
    if "--check" not in sys.argv:
        return 0
    committed = np.load(COMMITTED / f"kic{kepid}.npz", allow_pickle=True)
    worst, missing = 0.0, []
    for key in committed.files:
        if key.startswith(("d_plain_", "err_", "direct_")):
            if key not in payload:
                missing.append(key)
                continue
            worst = max(worst, float(np.max(np.abs(committed[key] - payload[key]))))
    print(f"  committed arrays missing here: {len(missing)}; largest pixel difference {worst:.3g}")
    return 1 if missing or worst > 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())
