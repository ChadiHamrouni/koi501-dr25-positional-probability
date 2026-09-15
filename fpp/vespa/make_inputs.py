"""Write vespa's input files from the step outputs, or check the committed ones.

    python fpp/vespa/make_inputs.py            # compare with the committed files
    python fpp/vespa/make_inputs.py --write    # overwrite them

star.ini   Teff, log g and [Fe/H] of step 07 with the widened errors of
           koi501/config.py; Gaia G and parallax (step 01); 2MASS J, H, Ks
fpp.ini    period and radius ratio from the DR25 detection table (step 00);
           maxrad = difference-image offset + 3 errors; secthresh = weak
           secondary depth + 3 x (depth / robust statistic); the velocity
           contrast constraint of config
lc.csv     the folded light curve fpp/triceratops_fpp.py writes
UKIRT_J.cc the J-band contrast curve of step 12 as a running maximum, out to
           config.CONTRAST_VESPA_MAX_AS

The check compares values, not layout: every number in the regenerated files
must equal the committed one at the precision written. A difference means the
committed vespa run (results.txt) was made from other inputs and must be redone
with run_vespa.sh.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT))

from koi501 import config, inputs  # noqa: E402

USED = ["transit.period_tce_d", "transit.ror", "transit.dikco_msky_as", "transit.dikco_msky_err_as",
        "transit.weak_secondary_ppm", "transit.weak_secondary_robstat", "transit.kepmag",
        "star.gaia_g", "star.parallax_mas", "star.parallax_err_mas", "star.j_mag", "star.j_err",
        "star.h_mag", "star.h_err", "star.ks_mag", "star.ks_err", "star.teff_k", "star.teff_err_k",
        "star.logg", "star.logg_err", "star.feh", "star.feh_err"]
FOLDED_LC = ROOT / "results" / "fpp_triceratops_lc.csv"


def sexagesimal(value: float, hours: bool) -> str:
    v = abs(value / 15.0 if hours else value)
    d = int(v)
    m = int((v - d) * 60)
    s = ((v - d) * 60 - m) * 60
    sign = "-" if value < 0 else ""
    return f"{sign}{d:02d}:{m:02d}:{s:05.2f}"


def files() -> dict[str, str]:
    g = inputs.get
    maxrad = g("transit.dikco_msky_as") + config.LIMIT_SIGMA * g("transit.dikco_msky_err_as")
    wst, rob = g("transit.weak_secondary_ppm"), g("transit.weak_secondary_robstat")
    secthresh = (wst + config.LIMIT_SIGMA * wst / rob) / config.PPM
    fpp = (f"name = {config.KOI.replace('.', '_')}\n"
           f"ra = {sexagesimal(config.RA, True)}\n"
           f"dec = {sexagesimal(config.DEC, False)}\n\n"
           f"period = {g('transit.period_tce_d'):.6f}   # days\n"
           f"rprs = {g('transit.ror'):.6f}\n"
           f"photfile = lc.csv\n"
           f"cadence = {config.LONG_CADENCE_D:.7f}   # Kepler long cadence, days\n"
           f"band = Kepler\n\n"
           f"[constraints]\n"
           f"maxrad = {maxrad:.3f}   # arcsec, DR25 tce_dikco_msky + 3 errors\n"
           f"secthresh = {secthresh:.8f}   # DR25 weak secondary + 3 x depth / robust statistic\n"
           f"ccfiles = UKIRT_J.cc   # step 12\n"
           f"vcc = {config.VESPA_VCC[0]}, {config.VESPA_VCC[1]}   # APOGEE DR17, single-lined\n")
    star = (f"# {config.KOI} / KIC {config.KEPID}, adopted spectroscopic solution\n"
            f"Kepler = {g('transit.kepmag'):.3f}\n"
            f"G = {g('star.gaia_g'):.3f}, {config.VESPA_G_ERROR}\n"
            f"parallax = {g('star.parallax_mas'):.4f}, {g('star.parallax_err_mas'):.4f}   # mas, Gaia DR3\n"
            f"J = {g('star.j_mag'):.3f}, {g('star.j_err'):.3f}\n"
            f"H = {g('star.h_mag'):.3f}, {g('star.h_err'):.3f}\n"
            f"K = {g('star.ks_mag'):.3f}, {g('star.ks_err'):.3f}\n"
            f"Teff = {g('star.teff_k'):.0f}, {g('star.teff_err_k'):.0f}   # APOGEE DR17 ASPCAP\n"
            f"logg = {g('star.logg'):.3f}, {g('star.logg_err'):.3f}\n"
            f"feh = {g('star.feh'):.3f}, {g('star.feh_err'):.2f}\n")
    cc = json.loads((config.RESULTS / "12_contrast_curves.json").read_text(encoding="utf8"))["J"]
    running, rows = float("-inf"), []
    for sep, dmag in zip(cc["separation_as"], cc["delta_mag_5sigma"]):
        running = max(running, dmag)
        if sep <= config.CONTRAST_VESPA_MAX_AS:
            rows.append(f"{sep:.4f}  {running:.4f}")
    out = {"fpp.ini": fpp, "star.ini": star, "UKIRT_J.cc": "\n".join(rows) + "\n"}
    if FOLDED_LC.exists():
        out["lc.csv"] = FOLDED_LC.read_text(encoding="utf8")
    return out


def numbers(text: str) -> list[float]:
    body = "\n".join(line.split("#")[0] for line in text.splitlines())
    return [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", body)]


def main() -> int:
    made = files()
    if "--write" in sys.argv:
        for name, text in made.items():
            (HERE / name).write_text(text, encoding="utf8", newline="\n")
        (HERE / "measured_inputs.json").write_text(
            json.dumps({"measured_inputs": inputs.record(USED)}, indent=2), encoding="utf8")
        print(f"wrote {', '.join(made)}")
        return 0
    bad = []
    for name, text in made.items():
        committed = numbers((HERE / name).read_text(encoding="utf8"))
        if committed != numbers(text):
            bad.append(name)
    missing = [] if "lc.csv" in made else ["lc.csv (run fpp/triceratops_fpp.py first)"]
    for name in made:
        print(f"  {name:11s} {'differs' if name in bad else 'matches'}")
    for name in missing:
        print(f"  not checked: {name}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
