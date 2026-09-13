"""Step 12 - how faint a companion each archival follow-up image could have shown.

For each frame: the target's centroid, the full width at half maximum of its
radial profile, the target's flux, and then, in annuli of increasing radius, the
scatter among photometric apertures placed around the target after subtracting
the target's own azimuthally averaged profile. Five times that scatter is the
faintest companion the frame could have revealed at that separation.

The frames are the Kepler Community Follow-up Observing Program images, fetched
from ExoFOP by file id and checked against their SHA-256 checksums.

Reproduces the contrast curves of Section 5.2.
"""

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from koi501 import archives
from koi501.fits_image import read_image
from koi501.models import FitsImageHeader
from koi501.report import Table, write, write_table

FRAMES = [
    # band, ExoFOP file id, file name, header seeing converted to arcsec, sha256
    ("U", 121811, "501Ic-SH20110625u.fits", 1.46,
     "803b84a02ced8a20c200f793b31467585d316e5ded775f5c36d3459cef2bb5aa"),
    ("B", 124500, "501Ic-SH20110625b.fits", 1.44,
     "f13a64db63a4b02433abda8fe176b47359ffa7f20a22190839925f3c7b9ae3f0"),
    ("V", 127189, "501Ic-SH20110625v.fits", 1.31,
     "e558ba2aab0c53ae863d342bc8da98bed760eb95108a28e7e4f21c48f8cdba0a"),
    ("J", 132444, "501Ic-dc20140206UKJ.fits", 0.66,
     "2ccdb1f7e393f3aa9755e37bb18cb060dcc4152b51395b3b52d13294e5857518"),
]
SIGMA = 5.0


def centroid(img, guess, box):
    y0, x0 = guess
    ys, xs = np.mgrid[0:img.shape[0], 0:img.shape[1]]
    inside = ((xs - x0) ** 2 + (ys - y0) ** 2) < box ** 2
    weight = np.clip(img - np.median(img), 0, None) * inside
    if weight.sum() <= 0:
        return guess
    return float((ys * weight).sum() / weight.sum()), float((xs * weight).sum() / weight.sum())


def radial_fwhm(img, yc, xc, scale):
    ys, xs = np.mgrid[0:img.shape[0], 0:img.shape[1]]
    r = np.hypot(ys - yc, xs - xc)
    background = np.median(img[r > 0.35 * min(img.shape)])
    profile = img - background
    peak = profile[int(round(yc)), int(round(xc))]
    radii = np.linspace(0, 12, 200)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # rings finer than a pixel can be empty
        values = np.array([np.mean(profile[(r >= a) & (r < a + 0.5)]) for a in radii])
    ok = np.isfinite(values)
    below = np.where(values[ok] < 0.5 * peak)[0]
    half_width = radii[ok][below[0]] if len(below) else np.nan
    return float(2 * half_width * scale), float(background)


def contrast_curve(img, header: FitsImageHeader, seeing_as: float) -> dict:
    scale = header.pixel_scale_as
    ny, nx = img.shape
    yc, xc = centroid(img, (ny / 2, nx / 2), box=max(4, int(6 / scale)))
    fwhm_as, background = radial_fwhm(img, yc, xc, scale)
    if not np.isfinite(fwhm_as) or fwhm_as <= 0:
        fwhm_as = seeing_as
    aperture_px = max(2.0, 0.5 * fwhm_as / scale)

    ys, xs = np.mgrid[0:ny, 0:nx]
    r = np.hypot(ys - yc, xs - xc)
    profile = img - background
    ring_edges = np.arange(0, r.max() + 1.0, 1.0)
    ring = np.digitize(r, ring_edges) - 1
    ring_median = np.array([np.median(profile[ring == i]) if (ring == i).sum() > 4 else 0.0
                            for i in range(len(ring_edges))])
    residual = profile - ring_median[np.clip(ring, 0, len(ring_median) - 1)]

    def aperture_flux(y0, x0, image):
        return float(image[np.hypot(ys - y0, xs - x0) <= aperture_px].sum())

    target_flux = aperture_flux(yc, xc, profile)
    separations, limits = [], []
    for sep in np.arange(max(fwhm_as, 2 * aperture_px * scale), 0.45 * nx * scale,
                         max(0.5 * fwhm_as, scale)):
        rr = sep / scale
        n_apertures = max(8, int(2 * np.pi * rr / (2 * aperture_px)))
        fluxes = []
        for theta in np.linspace(0, 2 * np.pi, n_apertures, endpoint=False):
            y1, x1 = yc + rr * np.sin(theta), xc + rr * np.cos(theta)
            if aperture_px < x1 < nx - aperture_px and aperture_px < y1 < ny - aperture_px:
                fluxes.append(aperture_flux(y1, x1, residual))
        if len(fluxes) < 6:
            continue
        fluxes = np.array(fluxes)
        scatter = 1.4826 * np.median(np.abs(fluxes - np.median(fluxes)))
        if scatter <= 0 or target_flux <= 0:
            continue
        separations.append(float(sep))
        limits.append(float(2.5 * np.log10(target_flux / (SIGMA * scatter))))

    return {"pixel_scale_as": scale, "fwhm_as": fwhm_as, "target_flux": target_flux,
            "inner_working_angle_as": max(fwhm_as, min(separations)),
            "separation_as": separations, "delta_mag_5sigma": limits}


def main() -> bool:
    curves, rows = {}, []
    for band, file_id, name, seeing, sha256 in FRAMES:
        img, cards = read_image(archives.exofop_file(file_id, sha256))
        header = FitsImageHeader.model_validate(cards)
        curve = contrast_curve(img, header, seeing)
        curves[band] = {"file": name, "exofop_file_id": file_id, **curve}
        rows += [{"band": band, "separation_as": round(s, 3), "delta_mag_5sigma": round(d, 3)}
                 for s, d in zip(curve["separation_as"], curve["delta_mag_5sigma"])]

    write("12_contrast_curves", curves)
    write_table("contrast_curves",
                "5-sigma contrast curves of the archival follow-up images (Section 5.2)",
                [("band", "-", "U, B, V (WIYN 0.9 m) or J (UKIRT WFCAM)"),
                 ("separation_as", "arcsec", "separation from the target"),
                 ("delta_mag_5sigma", "mag", "faintest companion detectable at 5 sigma")],
                rows)

    j = curves["J"]
    sep, dmag = np.array(j["separation_as"]), np.array(j["delta_mag_5sigma"])
    table = Table("How faint a companion each follow-up image could have shown",
                  "Section 5.2")
    table("UKIRT J, measured FWHM (arcsec)", round(j["fwhm_as"], 2))
    table("UKIRT J, inner working angle (arcsec)", round(j["inner_working_angle_as"], 2))
    near = int(np.argmin(np.abs(sep - 1.5)))
    table(f"UKIRT J, 5 sigma limit at {sep[near]:.2f} arcsec (mag)", round(float(dmag[near]), 2))
    table("UKIRT J, deepest limit (mag)", round(float(dmag.max()), 2))
    table("  reached at (arcsec)", round(float(sep[dmag.argmax()]), 2))
    for band in ("U", "B", "V"):
        table(f"{band}, deepest limit (mag)", round(max(curves[band]["delta_mag_5sigma"]), 2))
    table.show("12_contrast_curves")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
