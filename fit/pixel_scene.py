"""Kepler's own pixels: the fifth test of the catalogue magnitudes (Section 2.2).

For each of the sixteen quarters the mean out-of-transit image of the target's
pixels (data/figure_diffimage/diffimage_stamps.fits) is modelled with the Kepler
pixel response function (Bryson et al. 2010; lightkurve.prf.KeplerPRF, which
reads the mission's 2011 PRF calibration files from MAST). Pixels stored as
exactly zero lie outside the download mask and are not used.

1. Scene fit. Every Gaia DR3 source within 40 arcsec (data/figure_diffimage/
   gaia_neighbours.csv) is placed at its Gaia position through the stamp's WCS.
   The fluxes of the target, KIC 4951867 and KIC 4951861 are free, with one
   common positional shift and a constant background; the other sources are
   held at their Gaia G-band flux ratios to the target. The two neighbours'
   flux ratios are compared with those implied by the KIC Kepler magnitudes
   and by Gaia (step 01).
2. Light centre. The DV report gives, per quarter, the PRF-fit location of the
   target in the out-of-transit image and the KIC reference position
   (data/dv/dvr_quarterly_oot_centroids.csv). A single point source (position,
   flux, background) is fitted to two model scenes, all sources at their G
   fluxes or the target and the two neighbours at their KIC magnitudes, and to
   the stored image itself, and each fitted position is compared with the
   report's.

Every fit is made with uniform pixel weights and with photon-noise weights
(variance proportional to counts). The paper quotes the uniform-weight values.

    pip install -r fit/requirements.txt
    python fit/pixel_scene.py              # about two minutes
    python fit/pixel_scene.py out.json     # write somewhere else

Reads step 01's results/01_target_photometry.json. Writes results/pixel_scene.json;
the PRF files are cached in results/.cache/kepler_prf/.
"""

from __future__ import annotations

import csv
import json
import sys
import urllib.request
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from lightkurve.prf import KeplerPRF
from lightkurve.utils import module_output_to_channel
from scipy.optimize import least_squares

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent.parent
STAMPS = ROOT / "data" / "figure_diffimage" / "diffimage_stamps.fits"
GAIA = ROOT / "data" / "figure_diffimage" / "gaia_neighbours.csv"
DV = ROOT / "data" / "dv" / "dvr_quarterly_oot_centroids.csv"
PHOT = ROOT / "results" / "01_target_photometry.json"
PRF_DIR = ROOT / "results" / ".cache" / "kepler_prf"
PRF_URL = "https://archive.stsci.edu/missions/kepler/fpc/prf/"
OUT = ROOT / "results" / "pixel_scene.json"

KIC_RA_H, KIC_DEC = 19.89821400, 40.07587000      # the report's KIC reference
GAIA_RADIUS_AS = 40.0
NAMED = {2073562489151839872: 4951877,              # Gaia DR3 source -> KIC
         2073562489151838080: 4951867,
         2073562695310268032: 4951861}
WEIGHTINGS = ("uniform", "photon")


def local_prf_files():
    """lightkurve opens each PRF file from its URL once per extension and per
    stamp; read a local copy instead, downloading it the first time."""
    original = KeplerPRF._read_prf_calibration_file

    def read(self, path, ext):
        local = PRF_DIR / Path(path).name
        if not local.exists():
            PRF_DIR.mkdir(parents=True, exist_ok=True)
            local.write_bytes(urllib.request.urlopen(PRF_URL + local.name, timeout=600).read())
        return original(self, str(local), ext)
    KeplerPRF._read_prf_calibration_file = read


def dv_centroids() -> dict:
    out = {}
    for r in csv.DictReader(DV.open(encoding="utf8")):
        out[int(r["quarter"])] = {
            key: dict(row=float(r[f"{key}_row"]), col=float(r[f"{key}_col"]),
                      ra_h=float(r[f"{key}_ra_h"]), dec=float(r[f"{key}_dec"]))
            for key in ("oot", "kic")}
    return out


def sky_offset(ra_h, dec, ra0_h=KIC_RA_H, dec0=KIC_DEC):
    """Arcsec east and north of a reference, RA given in hours."""
    return ((ra_h - ra0_h) * 15.0 * np.cos(np.radians(dec0)) * 3600.0,
            (dec - dec0) * 3600.0)


def load_gaia() -> list[dict]:
    rows = []
    for line in GAIA.read_text(encoding="utf8").splitlines():
        if line.startswith("#") or line.startswith("source_id"):
            continue
        s = line.split(",")
        rows.append(dict(source_id=int(s[0]), ra=float(s[1]), dec=float(s[2]),
                         G=float(s[3]), sep_as=float(s[6])))
    return rows


class Stamp:
    def __init__(self, hdu):
        h = hdu.header
        self.image = np.asarray(hdu.data, float)
        self.shape = self.image.shape
        self.col0, self.row0 = int(h["COLUMN0"]), int(h["ROW0"])
        self.module, self.output = int(h["MODULE"]), int(h["OUTPUT"])
        self.wcs = WCS(h)
        self.prf = KeplerPRF(channel=module_output_to_channel(self.module, self.output),
                             shape=self.shape, column=self.col0, row=self.row0)
        # pixels outside the target's download mask are stored as exactly zero
        self.good = np.isfinite(self.image) & (self.image != 0)

    def pixel(self, ra, dec):
        """lightkurve PRF coordinates: pixel centres at column0 + i + 0.5."""
        x, y = self.wcs.world_to_pixel_values(ra, dec)
        return self.col0 + float(x) + 0.5, self.row0 + float(y) + 0.5

    def sky(self, col, row):
        c = self.wcs.pixel_to_world_values(col - self.col0 - 0.5, row - self.row0 - 0.5)
        ra, dec = float(c[0]), float(c[1])
        return ((ra - KIC_RA_H * 15.0) * np.cos(np.radians(KIC_DEC)) * 3600.0,
                (dec - KIC_DEC) * 3600.0)

    def model(self, stars):
        img = np.zeros(self.shape)
        for col, row, flux in stars:
            img += self.prf.evaluate(col, row, flux)
        return img

    def weights(self, image, weighting):
        if weighting == "uniform":
            return np.ones(int(self.good.sum()))
        return 1.0 / np.sqrt(np.clip(image[self.good], 1.0, None))

    def fit_single(self, image, weighting):
        """One point source plus a constant background on the measured pixels."""
        good, w = self.good, self.weights(image, weighting)
        r0, c0 = np.unravel_index(np.nanargmax(np.where(good, image, np.nan)), image.shape)
        p0 = [self.col0 + c0 + 0.5, self.row0 + r0 + 0.5, np.nansum(image[good]), 0.0]

        def resid(p):
            return (self.prf.evaluate(p[0], p[1], p[2]) + p[3] - image)[good] * w
        return least_squares(resid, p0, x_scale=[0.1, 0.1, p0[2] * 0.1, 10.0]).x


def main() -> int:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else OUT
    local_prf_files()
    dv = dv_centroids()
    hdus = [h for h in fits.open(STAMPS)[1:] if h.name.startswith("DIRECT")]
    gaia = load_gaia()
    phot = {s["kic"]: s for s in json.loads(PHOT.read_text(encoding="utf8"))["stars"]}
    g_target = next(g["G"] for g in gaia if NAMED.get(g["source_id"]) == 4951877)
    idx = {NAMED[g["source_id"]]: i for i, g in enumerate(gaia) if g["source_id"] in NAMED}
    if len(dv) != len(hdus):
        raise SystemExit("the report's quarters and the stored stamps do not pair up")

    per_quarter = []
    for hdu, q in zip(hdus, sorted(dv)):
        st = Stamp(hdu)
        d = dv[q]
        # stamps are stored in quarter order; check the pairing on the KIC reference pixel
        ref_row = st.row0 + hdu.header["CRPIX2"] - 1
        ref_col = st.col0 + hdu.header["CRPIX1"] - 1
        if abs(d["kic"]["row"] - ref_row) > 0.1 or abs(d["kic"]["col"] - ref_col) > 0.1:
            raise SystemExit(f"stamp {hdu.name} does not match quarter {q}")
        target_col, target_row = ref_col + 0.5, ref_row + 0.5

        def scene(kic_mags):
            stars = []
            for g in gaia:
                kic = NAMED.get(g["source_id"])
                mag = phot[kic]["kic_kepmag"] if (kic_mags and kic) else g["G"]
                ref = phot[4951877]["kic_kepmag"] if kic_mags else g_target
                col, row = st.pixel(g["ra"], g["dec"])
                stars.append((col, row, 10 ** (-0.4 * (mag - ref))))
            return stars

        tot = float(np.nansum(st.image[st.good]))
        scale = tot / float(np.sum(st.model(scene(False))[st.good]))
        synthetic = {"gaia_scene": st.model(scene(False)) * scale,
                     "kic_scene": st.model(scene(True)) * scale}
        base = scene(False)

        def scene_image(p):
            dx, dy, ft, f1, f2, bg = p
            img = np.full(st.shape, bg)
            for i, (col, row, ratio) in enumerate(base):
                flux = {idx[4951867]: f1, idx[4951861]: f2, idx[4951877]: ft}.get(i, ft * ratio)
                img += st.prf.evaluate(col + dx, row + dy, flux)
            return img

        rec = dict(quarter=q, module=st.module, output=st.output,
                   n_pixels=int(st.good.sum()),
                   observed_oot_minus_kic_as=list(sky_offset(d["oot"]["ra_h"], d["oot"]["dec"],
                                                             d["kic"]["ra_h"], d["kic"]["dec"])),
                   observed_oot_minus_kic_px=[d["oot"]["col"] - d["kic"]["col"],
                                              d["oot"]["row"] - d["kic"]["row"]])
        for weighting in WEIGHTINGS:
            fitted = {}
            for label, image in synthetic.items():
                p = st.fit_single(image, weighting)
                fitted[f"predicted_{label}_as"] = list(st.sky(p[0], p[1]))
            p = st.fit_single(st.image, weighting)
            fitted["real_stamp_single_prf_as"] = list(st.sky(p[0], p[1]))
            fitted["real_stamp_single_prf_px"] = [float(p[0] - target_col),
                                                  float(p[1] - target_row)]
            w = st.weights(st.image, weighting)
            sol = least_squares(lambda p: (scene_image(p) - st.image)[st.good] * w,
                                [0.0, 0.0, tot, 0.05 * tot, 0.05 * tot, 0.0],
                                x_scale=[0.05, 0.05, 0.1 * tot, 0.05 * tot, 0.05 * tot, 10.0])
            dx, dy, ft, f1, f2, bg = sol.x
            resid = (scene_image(sol.x) - st.image)[st.good]
            peak = float(np.max(st.image[st.good]))
            fitted["scene_fit"] = dict(shift_px=[float(dx), float(dy)],
                                       ratio_4951867=float(f1 / ft), ratio_4951861=float(f2 / ft),
                                       background=float(bg),
                                       max_abs_residual_over_peak=float(np.max(np.abs(resid)) / peak))
            rec[weighting] = fitted
        per_quarter.append(rec)
        print(f"  Q{q:02d} module {st.module:2d}.{st.output}  flux ratios "
              f"{rec['uniform']['scene_fit']['ratio_4951867']:.3f} "
              f"{rec['uniform']['scene_fit']['ratio_4951861']:.3f}", flush=True)

    def vectors(key, weighting=None):
        v = np.array([(r[weighting] if weighting else r)[key] for r in per_quarter])
        return v, np.hypot(v[:, 0], v[:, 1])

    obs, obs_d = vectors("observed_oot_minus_kic_as")
    obs_px, _ = vectors("observed_oot_minus_kic_px")
    catalogue = {
        "kic": {k: 10 ** (-0.4 * (phot[k]["kic_kepmag"] - phot[4951877]["kic_kepmag"]))
                for k in (4951867, 4951861)},
        "gaia": {k: 10 ** (-0.4 * (phot[k]["gaia_G"] - phot[4951877]["gaia_G"]))
                 for k in (4951867, 4951861)},
    }
    summary = dict(observed=dict(mean_as=obs.mean(axis=0).tolist(),
                                 median_distance_as=float(np.median(obs_d)),
                                 max_distance_as=float(obs_d.max())))
    for weighting in WEIGHTINGS:
        real, _ = vectors("real_stamp_single_prf_as", weighting)
        real_px, _ = vectors("real_stamp_single_prf_px", weighting)
        gaia_v, gaia_d = vectors("predicted_gaia_scene_as", weighting)
        kic_v, kic_d = vectors("predicted_kic_scene_as", weighting)
        r1 = np.array([r[weighting]["scene_fit"]["ratio_4951867"] for r in per_quarter])
        r2 = np.array([r[weighting]["scene_fit"]["ratio_4951861"] for r in per_quarter])
        worst = np.array([r[weighting]["scene_fit"]["max_abs_residual_over_peak"]
                          for r in per_quarter])
        summary[weighting] = dict(
            real_fit_minus_report=dict(
                rms_as=float(np.sqrt(np.mean(np.sum((real - obs) ** 2, axis=1)))),
                rms_px=float(np.sqrt(np.mean(np.sum((real_px - obs_px) ** 2, axis=1))))),
            predicted_gaia=dict(mean_as=gaia_v.mean(axis=0).tolist(),
                                median_distance_as=float(np.median(gaia_d))),
            predicted_kic=dict(mean_as=kic_v.mean(axis=0).tolist(),
                               median_distance_as=float(np.median(kic_d)),
                               min_distance_as=float(kic_d.min()),
                               quarters_above_largest_observed=int((kic_d > obs_d.max()).sum())),
            flux_ratio_4951867=dict(median=float(np.median(r1)), range=[float(r1.min()), float(r1.max())]),
            flux_ratio_4951861=dict(median=float(np.median(r2)), range=[float(r2.min()), float(r2.max())]),
            scene_worst_pixel_over_peak=dict(median=float(np.median(worst)), max=float(worst.max())))

    result = dict(n_quarters=len(per_quarter),
                  modules=sorted({r["module"] for r in per_quarter}),
                  n_gaia_sources=len(gaia), gaia_max_separation_as=max(g["sep_as"] for g in gaia),
                  gaia_radius_as=GAIA_RADIUS_AS,
                  stamp_pixels_measured=[min(r["n_pixels"] for r in per_quarter),
                                         max(r["n_pixels"] for r in per_quarter)],
                  catalogue_flux_ratio={s: {str(k): v for k, v in c.items()}
                                        for s, c in catalogue.items()},
                  summary=summary, per_quarter=per_quarter,
                  note=("pixel positions in the report are given to 0.01 pixel; "
                        "G stands in for the Kepler magnitude of Gaia sources"))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2), encoding="utf8")

    u = summary["uniform"]
    print(f"\n  KIC 4951867 carries {u['flux_ratio_4951867']['median']:.3f} of the target's flux "
          f"(KIC {catalogue['kic'][4951867]:.3f}, Gaia {catalogue['gaia'][4951867]:.3f})")
    print(f"  KIC 4951861 carries {u['flux_ratio_4951861']['median']:.3f} "
          f"(KIC {catalogue['kic'][4951861]:.3f}, Gaia {catalogue['gaia'][4951861]:.3f})")
    print(f"  light centre: report {summary['observed']['median_distance_as']:.3f} arcsec; "
          f"model with Gaia magnitudes {u['predicted_gaia']['median_distance_as']:.2f}, "
          f"with KIC magnitudes {u['predicted_kic']['median_distance_as']:.1f}")
    print(f"  wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
