#!/usr/bin/env bash
# vespa, the first false-positive calculation of Section 7, on the inputs in
# this folder. Run under Linux (the paper's run used WSL Ubuntu) in a Python 3.8
# environment with vespa 0.6, isochrones 1.2.2, numpy 1.19.5, scipy 1.5.4,
# pandas 0.25.3, tables 3.6.1, emcee 2.2.1, astropy 4.0.4, numba 0.53.1,
# configobj 5.0.9 and matplotlib 3.3.4.
#
# Inputs
#   star.ini    the adopted stellar solution: Teff 5764 +/- 60 K and log g
#               4.054 +/- 0.035 (APOGEE DR17), [Fe/H] 0.004, the Gaia DR3
#               parallax, and Gaia G and 2MASS J, H, Ks photometry
#   fpp.ini     period, radius ratio, 29.4-minute cadence; maxrad 1.440 arcsec
#               (the DR25 offset 0.675 plus three times its error 0.255);
#               secthresh 96.6 ppm (the DR25 weak secondary depth, 44.9 ppm,
#               plus three times that depth over its robust statistic, 2.60 -
#               looser than the 90.2 ppm limit of Table 4); the UKIRT J
#               contrast curve; the APOGEE velocity constraint (25 km/s, 3 mag)
#   lc.csv      the folded, binned light curve fpp/triceratops_fpp.py builds
#               (days from mid-transit, normalised flux, error), written to
#               eight decimals
#   UKIRT_J.cc  the J-band contrast curve of step 12 as a running maximum
#
# Outputs of the paper's run, committed here: starfit.log, calcfpp.log,
# results.txt (FPP 1.69e-4), results_bootstrap.txt (100 resamplings) and
# FPPsummary.png. The simulated populations (popset.h5, 139 MB) and the TRILEGAL
# field (starfield.h5, 142 MB) exceed GitHub's file limit and are not included.
# calcfpp regenerates both, as new random realisations, so a rerun agrees with
# the committed result within the bootstrap scatter rather than exactly.
set -euo pipefail
cd "$(dirname "$0")"

starfit --all .              # single, binary and triple stellar models (about 1 minute)
calcfpp .                    # populations, TRILEGAL query and FPP (about 10 minutes)
calcfpp --bootstrap 100 .    # 100 bootstrap resamplings of the cached populations
