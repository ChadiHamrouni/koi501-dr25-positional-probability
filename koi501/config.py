"""Constants and paths. Every threshold used in the paper is defined here."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

KOI = "K00501.01"
KEPID = 4951877
RA = 298.473210
DEC = 40.075871

#: The two neighbours the DR25 positional probability favours over the host.
NEIGHBOURS = (4951867, 4951861)

#: Blend search radius, arcsec. Armstrong et al. (2021) adopt the same.
BLEND_RADIUS = 25.0

#: r - J below this is not a stellar colour; even a white dwarf stays above it.
#: Normal stars run from about +0.3 (hot) to +4 and redder.
COLOUR_IMPOSSIBLE = -0.5
COLOUR_SUSPICIOUS = 0.0

#: A corrupted neighbour can only take the transit if the scene model can
#: plausibly place it there: close, and carrying real light.
CONFIG_SEP_MAX = 8.0
CONFIG_SHARE_MIN = 0.25

#: Grid used to show the two thresholds above are not what drives the result.
CONFIG_GRID_SEP = (6.0, 8.0, 10.0, 12.0, 15.0)
CONFIG_GRID_SHARE = (0.15, 0.25, 0.35)

#: Depth above which Bryson & Morton (2017) reject a star outright, ppm.
REJECTION_DEPTH_PPM = 3.0e6

#: Adopted stellar mass, solar units.
STELLAR_MASS = 1.26

#: Orbital period, days, and transit epoch, BKJD.
PERIOD_D = 24.7963
EPOCH_BKJD = 145.529
