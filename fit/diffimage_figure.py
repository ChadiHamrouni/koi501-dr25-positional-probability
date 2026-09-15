"""Figure 1: the difference images of KOI-501.01, and the data behind them.

Stage 1 (--data) writes data/figure_diffimage/ from the committed difference
images (data/diffimages/kic4951877.npz, built by fit/diff_images.py):

  diffimage_stamps.fits  per quarter the out-of-transit image (DIRECT_Qnn) and
                         the difference image (DIFF_Qnn), with module, output,
                         origin pixel, transit count and WCS
  gaia_neighbours.csv    Gaia DR3 sources within 40 arcsec, with offsets from
                         the target and pixel positions on the first stamp

Stage 2 (default) draws the figure from those files:

  row (a)  all quarters, each shifted so the target lands on the centre of a
           7 x 7 grid before averaging
  row (b)  one season, the quarters on module 6, which share a stamp shape and
           place the target on the same pixel, averaged without resampling

    python fit/diffimage_figure.py --data     # rebuild the figure data (queries Gaia)
    python fit/diffimage_figure.py            # draw figures/K00501.01_diffimage.png

pixel_scene.py reads the same stamps and Gaia list.
"""

from __future__ import annotations

import csv
import json
import math
import sys
import warnings
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.io import fits  # noqa: E402
from astropy.wcs import WCS  # noqa: E402
from scipy.ndimage import shift as nd_shift  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore", module="astropy")

from koi501 import archives, config  # noqa: E402

DIFFIMAGES = ROOT / "data" / "diffimages" / f"kic{config.KEPID}.npz"
DATA = ROOT / "data" / "figure_diffimage"
FIGURE = "K00501.01_diffimage.png"
GRID = config.FIGURE_GRID
CENTRE = (GRID - 1) / 2.0

plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.linewidth": 0.8,  # literal: display
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",  # literal: display
    "axes.labelsize": 9, "axes.titlesize": 9,  # literal: display
})


def build_data() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    d = np.load(DIFFIMAGES, allow_pickle=True)
    tags = sorted((k[len("meta_"):] for k in d.files if k.startswith("meta_")),
                  key=lambda t: int(t.split("_")[1]))
    hdus = [fits.PrimaryHDU()]
    hdus[0].header["OBJECT"] = f"{config.KOI} / KIC {config.KEPID}"
    hdus[0].header["COMMENT"] = "DIRECT = out-of-transit mean; DIFF = out minus in"
    shapes, stamps = Counter(), []
    for tag in tags:
        meta = json.loads(str(d["meta_" + tag]))
        direct = np.asarray(d["direct_" + tag], float)
        diff = np.asarray(d["d_plain_" + tag], float)
        header = fits.Header.fromstring(meta["wcs"])
        for key, value in (("QUARTER", meta["q"]), ("MODULE", meta["module"]),
                           ("OUTPUT", meta["output"]), ("COLUMN0", meta["column"]),
                           ("ROW0", meta["row"]), ("NTRANSIT", meta["n_transits"])):
            header[key] = value
        label = f"Q{meta['q']:02d}"
        hdus.append(fits.ImageHDU(direct, header, name=f"DIRECT_{label}"))
        hdus.append(fits.ImageHDU(diff, header, name=f"DIFF_{label}"))
        stamps.append((meta, direct, diff, header))
        shapes[diff.shape] += 1

    # the stack of stamps sharing the most common shape, kept for reference
    shape = shapes.most_common(1)[0][0]
    used = [s for s in stamps if s[2].shape == shape]
    stack_header = used[0][3].copy()
    stack_header["NSTACK"] = len(used)
    stack_header["QUARTERS"] = ",".join(str(s[0]["q"]) for s in used)
    stack_header["MODULES"] = ",".join(str(s[0]["module"]) for s in used)
    stack_header["COMMENT"] = "mean of stamps sharing the most common shape; WCS of the first"
    hdus.append(fits.ImageHDU(np.mean([s[1] for s in used], axis=0), stack_header, name="STACK_DIRECT"))
    hdus.append(fits.ImageHDU(np.mean([s[2] for s in used], axis=0), stack_header, name="STACK_DIFF"))
    fits.HDUList(hdus).writeto(DATA / "diffimage_stamps.fits", overwrite=True)

    wcs = WCS(stack_header)
    target_x, target_y = wcs.all_world2pix(config.RA, config.DEC, 0)
    gaia = archives.vizier_cone("I/355/gaiadr3", "Source,RA_ICRS,DE_ICRS,Gmag",
                                config.RA, config.DEC, config.PIXEL_GAIA_RADIUS_AS)
    rows = []
    for g in gaia:
        ra, dec = float(g["RA_ICRS"]), float(g["DE_ICRS"])
        dra = (ra - config.RA) * math.cos(math.radians(config.DEC)) * 3600.0
        ddec = (dec - config.DEC) * 3600.0
        sep = math.hypot(dra, ddec)
        if sep >= config.PIXEL_GAIA_RADIUS_AS or not g["Gmag"]:
            continue
        x, y = wcs.all_world2pix(ra, dec, 0)
        rows.append({"source_id": g["Source"], "ra": ra, "dec": dec, "phot_g_mean_mag": g["Gmag"],
                     "dra_as": round(dra, 3), "ddec_as": round(ddec, 3),  # literal: display
                     "sep_as": round(sep, 3), "x_pix": round(float(x), 3),  # literal: display
                     "y_pix": round(float(y), 3)})  # literal: display
    rows.sort(key=lambda r: r["sep_as"])
    with (DATA / "gaia_neighbours.csv").open("w", newline="", encoding="utf8") as handle:
        handle.write(f"# Gaia DR3 sources within {config.PIXEL_GAIA_RADIUS_AS:g} arcsec of "
                     f"KIC {config.KEPID} marked in the figure\n")
        handle.write("# ra, dec [deg]; phot_g_mean_mag [mag]; dra_as, ddec_as, sep_as [arcsec] "
                     "from the target; x_pix, y_pix [pixel, 0-based] on the STACK WCS\n")
        handle.write(f"# target at x_pix={float(target_x):.3f}, y_pix={float(target_y):.3f}\n")
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"  wrote {DATA.relative_to(ROOT).as_posix()}: {len(stamps)} quarters, {len(rows)} Gaia sources")


def stamps():
    with fits.open(DATA / "diffimage_stamps.fits") as hdul:
        for hdu in hdul:
            if hdu.name.startswith("DIFF_Q"):
                label = hdu.name[len("DIFF_"):]
                yield (hdu.header["QUARTER"], hdu.header["MODULE"],
                       hdul["DIRECT_" + label].data.astype(float),
                       hdu.data.astype(float), hdu.header)


def aligned_stack(selected):
    """Place each stamp on a GRID x GRID canvas with the target at CENTRE."""
    directs, diffs = [], []
    for _, _, direct, diff, header in selected:
        x, y = WCS(header).all_world2pix(config.RA, config.DEC, 0)
        out = []
        for image in (direct, diff):
            canvas = np.full((GRID, GRID), np.nan)
            ny, nx = image.shape
            oy, ox = int(round(CENTRE - y)), int(round(CENTRE - x))
            y0, x0 = max(oy, 0), max(ox, 0)
            y1, x1 = min(oy + ny, GRID), min(ox + nx, GRID)
            canvas[y0:y1, x0:x1] = image[y0 - oy:y1 - oy, x0 - ox:x1 - ox]
            residual = (CENTRE - y - oy, CENTRE - x - ox)
            valid = np.isfinite(canvas).astype(float)
            moved = nd_shift(np.nan_to_num(canvas), residual, order=1, mode="constant")
            weight = nd_shift(valid, residual, order=1, mode="constant")
            out.append(np.where(weight > 0.99, moved / np.maximum(weight, 1e-9), np.nan))  # literal: technical
        directs.append(out[0])
        diffs.append(out[1])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmean(directs, axis=0), np.nanmean(diffs, axis=0)


def neighbours():
    with (DATA / "gaia_neighbours.csv").open(encoding="utf8") as handle:
        rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    return [r for r in rows if float(r["sep_as"]) > config.FIGURE_TARGET_MATCH_AS]


def draw_row(axes, direct, diff, target_xy, gaia_xy, title):
    tx, ty = target_xy
    for ax, image, name in ((axes[0], direct, "direct image: all the light"),
                            (axes[1], diff, "difference image: the light that dims")):
        im = ax.imshow(np.ma.masked_invalid(image), origin="lower", cmap="magma",
                       interpolation="nearest")
        cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)  # literal: display
        cb.ax.tick_params(labelsize=6.5)
        h, w = image.shape
        for gx, gy, size in gaia_xy:
            if -0.5 < gx < w - 0.5 and -0.5 < gy < h - 0.5:
                ax.scatter(gx, gy, s=size, marker="o", facecolor="none",
                           edgecolor="#7ec8e3", lw=1.0, zorder=5)
        ax.plot(tx, ty, "+", color="#39ff14", ms=13, mew=2.0, zorder=6)
        peak = np.unravel_index(np.nanargmax(image), image.shape)
        ax.plot(peak[1], peak[0], "s", ms=15, mfc="none", mec="white", mew=0.8,
                ls="none", zorder=4)
        rows = np.where(np.isfinite(image).any(axis=1))[0]
        cols = np.where(np.isfinite(image).any(axis=0))[0]
        ax.set_xlim(cols[0] - 0.5, cols[-1] + 0.5)
        ax.set_ylim(rows[0] - 0.5, rows[-1] + 0.5)
        ax.set_title(name, fontsize=8.5)
        ax.set_xlabel("pixels")
    axes[0].set_ylabel(f"{title}\npixels")


def draw() -> None:
    all_stamps = list(stamps())
    stack_direct, stack_diff = aligned_stack(all_stamps)
    season = [s for s in all_stamps if s[1] == config.FIGURE_SEASON_MODULE]
    season_direct = np.mean([s[2] for s in season], axis=0)
    season_diff = np.mean([s[3] for s in season], axis=0)
    season_wcs = WCS(season[0][4])
    stx, sty = season_wcs.all_world2pix(config.RA, config.DEC, 0)

    gaia = neighbours()
    g_mag = np.array([float(r["phot_g_mean_mag"]) for r in gaia])
    sizes = 12 + 110 * 10 ** (-0.4 * (g_mag - g_mag.min()))  # literal: display, marker area by brightness
    # every stamp carries the same sky orientation, so a source's pixel offset
    # from the target is the same in all of them; the grid only moves the origin
    reference = WCS(all_stamps[0][4])
    rtx, rty = reference.all_world2pix(config.RA, config.DEC, 0)
    grid_xy, season_xy = [], []
    for r, s in zip(gaia, sizes):
        ra, dec = float(r["ra"]), float(r["dec"])
        gx, gy = reference.all_world2pix(ra, dec, 0)
        grid_xy.append((CENTRE + float(gx - rtx), CENTRE + float(gy - rty), s))
        sx, sy = season_wcs.all_world2pix(ra, dec, 0)
        season_xy.append((float(sx), float(sy), s))

    header = fits.Header()
    header["OBJECT"] = f"{config.KOI} / KIC {config.KEPID}"
    header["NSTACK"] = len(all_stamps)
    header["TARGETX"], header["TARGETY"] = CENTRE, CENTRE
    header["COMMENT"] = "all quarters shifted to put the target at (TARGETX, TARGETY), then averaged"
    module = config.FIGURE_SEASON_MODULE
    fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(stack_direct, header, name="ALIGNED_DIRECT"),
                  fits.ImageHDU(stack_diff, header, name="ALIGNED_DIFF"),
                  fits.ImageHDU(season_direct, season[0][4], name=f"SEASON_M{module}_DIRECT"),
                  fits.ImageHDU(season_diff, season[0][4], name=f"SEASON_M{module}_DIFF")]
                 ).writeto(DATA / "diffimage_stacks.fits", overwrite=True)

    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.2))
    draw_row(axes[0], stack_direct, stack_diff, (CENTRE, CENTRE), grid_xy,
             f"(a) all {len(all_stamps)} quarters, aligned")
    quarters = ", ".join(f"Q{s[0]}" for s in season)
    draw_row(axes[1], season_direct, season_diff, (float(stx), float(sty)), season_xy,
             f"(b) module {module} only ({quarters})")
    handles = [
        plt.Line2D([], [], marker="+", color="#2e8b1e", ms=9, mew=2, ls="none", label="the target star"),
        plt.Line2D([], [], marker="o", mfc="none", mec="#3a97bd", ms=7, ls="none",
                   label="Gaia DR3 sources, sized by brightness"),
        plt.Line2D([], [], marker="s", mfc="none", mec="#555555", ms=8, ls="none",
                   label="brightest pixel of each image"),
    ]
    fig.tight_layout(rect=(0, 0.04, 1, 1))  # literal: display
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=7, frameon=False)
    config.FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(config.FIGURES / FIGURE)
    plt.close(fig)

    for name, image, (tx, ty) in (("aligned direct", stack_direct, (CENTRE, CENTRE)),
                                  ("aligned diff", stack_diff, (CENTRE, CENTRE)),
                                  ("season direct", season_direct, (stx, sty)),
                                  ("season diff", season_diff, (stx, sty))):
        py, px = np.unravel_index(np.nanargmax(image), image.shape)
        print(f"  {name:15s} brightest pixel ({px}, {py}); target at ({float(tx):.2f}, {float(ty):.2f});"
              f" offset {np.hypot(px - tx, py - ty) * config.KEPLER_PIXEL_AS:.1f} arcsec")
    print(f"  wrote figures/{FIGURE} and {(DATA / 'diffimage_stacks.fits').relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    if "--data" in sys.argv:
        build_data()
    draw()
