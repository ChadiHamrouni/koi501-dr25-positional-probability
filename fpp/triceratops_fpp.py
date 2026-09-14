"""TRICERATOPS, the second false-positive calculation of Section 7: twenty runs.

TRICERATOPS (Giacalone et al. 2021) weighs the transit scenarios - a planet on
the target, eclipsing binaries on the target, on a bound companion or on a
background star, and transits or eclipses on each resolved neighbour - using the
folded light curve, the stars around the target and a contrast curve. A single
run moves by factors of several on identical inputs, so twenty are made and the
mean is quoted.

Inputs
  light curve  Kepler long-cadence PDCSAP from MAST (lightkurve, default quality
               mask), stitched, flattened with a 401-point window with the
               transits masked, clipped at 4 sigma above, folded within 3
               durations of mid-transit and averaged into 100 bins; the flux
               error is the scatter of the out-of-transit bins
  star         R = 1.746 R_sun, Teff = 5764 K, log g = 4.054, hence
               M = 10^(log g - 4.438) R^2 = 1.259 M_sun; the Gaia DR3 parallax
               from results/01_target_photometry.json (step 01)
  transit      depth 575.3 ppm, period 24.7963 d, epoch 145.529 BKJD,
               exposure 29.4 min
  field        TESS Input Catalog stars within 10 pixels. Rows the catalogue
               marks as DUPLICATE or SPLIT copies of the target, or places at
               zero separation, are removed: each would otherwise be offered as
               a neighbour able to host the transit
  aperture     the mission's optimal aperture from the first quarter's target
               pixel file (left alone, the library assumes a 5x5 box)
  imaging      the UKIRT J contrast curve of step 12, as a running maximum out to
               10 arcsec (results/12_contrast_curves.json)
  population   fpp/data/4951877_TRILEGAL.csv, the TRILEGAL simulation of the
               field that the published runs used

Bars (Giacalone et al. 2021): FPP < 0.015 and NFPP < 0.001. They were set on
two-minute TESS light curves; Kepler's long cadence is 29.4 minutes.

The published runs each rebuilt the field from MAST and were not seeded. This
script builds the field once, gives each run its own copy, and seeds numpy's
generator before each run's draws, so its reruns repeat; they agree with the
published values statistically rather than digit for digit.

    pip install -r fpp/requirements.txt
    python fpp/triceratops_fpp.py                 # 20 runs, about 10 minutes
    python fpp/triceratops_fpp.py 3 some/dir      # 3 runs, written to some/dir

Writes results/fpp_triceratops.json and results/fpp_triceratops_runs.csv.
"""

from __future__ import annotations

import copy
import csv
import json
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
warnings.filterwarnings("ignore")

import _compat  # noqa: F401,E402  (must precede triceratops)
import lightkurve as lk  # noqa: E402
import numpy as np  # noqa: E402
import triceratops.triceratops as T  # noqa: E402
from astroquery.mast import Catalogs, Observations  # noqa: E402

ROOT = HERE.parent
RESULTS = ROOT / "results"
TRILEGAL = HERE / "data" / "4951877_TRILEGAL.csv"

KEPID = 4951877
PERIOD_D, EPOCH_BKJD = 24.7963, 145.529
DURATION_H, DEPTH_PPM = 8.483, 575.3
EXPTIME_D = 29.4 / 60.0 / 24.0
R_STAR, TEFF, LOGG = 1.746, 5764.0, 4.054
M_STAR = 10 ** (LOGG - 4.438) * R_STAR ** 2
SEARCH_RADIUS_PX = 10
N_RUNS, N_DRAWS, SEED = 20, 100_000, 20260904
FPP_BAR, NFPP_BAR = 0.015, 0.001


def contrast_curve(folder: Path) -> Path:
    """Step 12's J-band curve in the comma-separated form triceratops reads."""
    cc = json.loads((RESULTS / "12_contrast_curves.json").read_text(encoding="utf8"))["J"]
    sep = np.array(cc["separation_as"])
    dmag = np.maximum.accumulate(np.array(cc["delta_mag_5sigma"]))
    keep = sep <= 10.0
    path = folder / "ukirt_J_triceratops.txt"
    np.savetxt(path, np.column_stack([sep[keep], dmag[keep]]), fmt="%.4f,%.4f")
    return path


def retried(call, attempts: int = 5):
    """MAST's portal times out intermittently; wait and ask again."""
    for client in (Observations, Catalogs):
        client._portal_api_connection.TIMEOUT = 180    # astroquery's default is 600 s
    for attempt in range(attempts):
        try:
            return call()
        except OSError as error:            # includes timeouts and dropped connections
            if attempt == attempts - 1:
                raise
            print(f"  MAST: {str(error)[:80]}; retrying in {60 * 2 ** attempt} s", flush=True)
            time.sleep(60 * 2 ** attempt)


def folded_light_curve(window: float = 3.0, nbin: int = 100):
    lc = retried(lambda: lk.search_lightcurve(f"KIC {KEPID}", mission="Kepler", cadence="long")
                 .download_all(quality_bitmask="default")).stitch().remove_nans()
    dur = DURATION_H / 24.0
    phase = (lc.time.value - EPOCH_BKJD + 0.5 * PERIOD_D) % PERIOD_D - 0.5 * PERIOD_D
    lc = lc.flatten(window_length=401, mask=np.abs(phase) < dur).remove_outliers(
        sigma_upper=4, sigma_lower=20)
    phase = (lc.time.value - EPOCH_BKJD + 0.5 * PERIOD_D) % PERIOD_D - 0.5 * PERIOD_D
    near = np.abs(phase) < window * dur
    t, f = phase[near], np.asarray(lc.flux.value)[near]
    ok = np.isfinite(t) & np.isfinite(f)
    order = np.argsort(t[ok])
    t, f = t[ok][order], f[ok][order]
    edges = np.linspace(t.min(), t.max(), nbin + 1)
    idx = np.clip(np.digitize(t, edges) - 1, 0, nbin - 1)
    bt = np.array([t[idx == i].mean() if (idx == i).any() else np.nan for i in range(nbin)])
    bf = np.array([f[idx == i].mean() if (idx == i).any() else np.nan for i in range(nbin)])
    good = np.isfinite(bt) & np.isfinite(bf)
    bt, bf = bt[good], bf[good]
    oot = np.abs(bt) > 1.5 * dur
    return bt, bf / np.median(bf[oot]), float(np.std(bf[oot]))


def first_quarter() -> int:
    sr = retried(lambda: lk.search_targetpixelfile(f"KIC {KEPID}", mission="Kepler",
                                                   cadence="long"))
    quarters = {int(str(m).split()[-1]) for m in sr.table["mission"]
                if str(m).split() and str(m).split()[-1].isdigit()}
    if not quarters:
        raise RuntimeError(f"no target pixel files for KIC {KEPID}")
    return min(quarters)


def mission_aperture(quarter: int, target):
    """The pipeline's optimal aperture in triceratops' absolute pixel frame
    (tpf.column + column index); refused unless it contains the target."""
    tpf = retried(lambda: lk.search_targetpixelfile(
        f"KIC {KEPID}", mission="Kepler", cadence="long", quarter=quarter)[0].download())
    rows, cols = np.where(tpf.pipeline_mask)
    ap = np.array([cols + tpf.column, rows + tpf.row]).T
    centre = np.round(target.pix_coords[0][0])
    if not ((ap[:, 0] == centre[0]) & (ap[:, 1] == centre[1])).any():
        raise RuntimeError("the target is not inside the pipeline aperture")
    return [ap for _ in range(len(target.pix_coords))]


def field(quarter: int, parallax: float):
    """The target, its neighbours and the aperture, built once from MAST.

    Each run then works on its own copy, so the archive is queried once rather
    than twenty times; nothing random happens before calc_probs."""
    target = retried(lambda: T.target(ID=KEPID, sectors=np.array([quarter]), mission="Kepler",
                                      search_radius=SEARCH_RADIUS_PX, verify_ssl=False,
                                      trilegal_fname=str(TRILEGAL)))
    # the TESS Input Catalog carries no mass, radius or Teff for this star
    tid = target.stars["ID"].values[0]
    for param, value in (("mass", M_STAR), ("rad", R_STAR), ("Teff", TEFF)):
        try:
            target.update_star(tid, param, value)
        except Exception:
            target.stars.loc[0, param] = value
    target.stars.loc[0, "plx"] = parallax

    stars = target.stars
    duplicates = [str(stars["ID"].values[i]) for i in range(1, len(stars))
                  if str(stars["disposition"].values[i]) in ("DUPLICATE", "SPLIT")
                  or str(stars["duplicate_id"].values[i]) not in ("nan", "None", "")
                  or float(stars["sep (arcsec)"].values[i]) == 0.0]
    for star_id in duplicates:
        try:
            target.remove_star(int(float(star_id)))
        except Exception as error:
            print(f"  could not remove {star_id}: {error}")

    return target, duplicates, mission_aperture(quarter, target)


def one_run(k: int, lc, prepared, cc_file: Path) -> dict:
    t, f, err = lc
    base, duplicates, aperture = prepared
    target = copy.deepcopy(base)
    target.calc_depths(tdepth=DEPTH_PPM * 1e-6, all_ap_pixels=aperture)
    np.random.seed(SEED + k)
    target.calc_probs(time=t, flux_0=f, flux_err_0=err, P_orb=PERIOD_D,
                      contrast_curve_file=str(cc_file), filt="J",
                      exptime=EXPTIME_D, N=N_DRAWS, verbose=0)
    row = {"run": k, "seed": SEED + k, "FPP": float(target.FPP), "NFPP": float(target.NFPP),
           "n_stars": len(target.stars), "n_duplicates_removed": len(duplicates),
           "aperture_pixels": len(aperture[0]), "n_draws": N_DRAWS}
    print(f"  run {k:2d}   FPP {row['FPP']:.5f}   NFPP {row['NFPP']:.5f}   "
          f"{row['n_stars']} stars, {row['n_duplicates_removed']} duplicates removed, "
          f"{row['aperture_pixels']}-pixel aperture", flush=True)
    return row


def summarise(rows: list[dict]) -> dict:
    fpp = np.array([r["FPP"] for r in rows])
    nfpp = np.array([r["NFPP"] for r in rows])
    return {"n_runs": len(rows),
            "fpp_mean": float(fpp.mean()),
            "fpp_sem": float(fpp.std(ddof=1) / np.sqrt(len(fpp))) if len(fpp) > 1 else None,
            "fpp_median": float(np.median(fpp)), "fpp_max": float(fpp.max()),
            "n_fpp_exactly_zero": int((fpp == 0).sum()),
            "nfpp_mean": float(nfpp.mean()), "nfpp_max": float(nfpp.max()),
            "n_pass_both_bars": int(((fpp < FPP_BAR) & (nfpp < NFPP_BAR)).sum())}


def main() -> int:
    n_runs = int(sys.argv[1]) if len(sys.argv) > 1 else N_RUNS
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else RESULTS
    out.mkdir(parents=True, exist_ok=True)

    parallax = json.loads((RESULTS / "01_target_photometry.json").read_text(
        encoding="utf8"))["host_gaia"]["parallax"]
    cc_file = contrast_curve(out)
    lc = folded_light_curve()
    quarter = first_quarter()
    print(f"folded light curve: {len(lc[0])} bins, out-of-transit scatter "
          f"{lc[2] * 1e6:.0f} ppm; target pixel file from quarter {quarter}")

    prepared = field(quarter, parallax)
    rows = [one_run(k, lc, prepared, cc_file) for k in range(1, n_runs + 1)]
    with (out / "fpp_triceratops_runs.csv").open("w", newline="", encoding="utf8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = summarise(rows)
    payload = {
        "inputs": {"star": {"radius_rsun": R_STAR, "teff_k": TEFF, "logg": LOGG,
                            "mass_msun": M_STAR, "parallax_mas": parallax},
                   "transit": {"depth_ppm": DEPTH_PPM, "period_d": PERIOD_D,
                               "epoch_bkjd": EPOCH_BKJD, "duration_h": DURATION_H,
                               "exposure_min": EXPTIME_D * 1440},
                   "light_curve": {"bins": len(lc[0]), "flux_err": lc[2]},
                   "search_radius_px": SEARCH_RADIUS_PX, "draws_per_run": N_DRAWS,
                   "contrast_curve": {"band": "J", "file": cc_file.name},
                   "trilegal": TRILEGAL.relative_to(ROOT).as_posix(),
                   "quarter": quarter},
        "bars": {"fpp": FPP_BAR, "nfpp": NFPP_BAR},
        **summary,
    }
    (out / "fpp_triceratops.json").write_text(json.dumps(payload, indent=2), encoding="utf8")
    print(f"mean FPP {summary['fpp_mean']:.2e}, NFPP {summary['nfpp_mean']:.5f}; "
          f"{summary['n_pass_both_bars']} of {len(rows)} runs clear both bars")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
