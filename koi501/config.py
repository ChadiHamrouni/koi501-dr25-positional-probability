"""Every fixed value the analysis uses, in one place.

What belongs here: thresholds and choices, assumptions taken from the
literature, algorithm settings that change a result (windows, grids, clipping,
priors, bounds, sample sizes, seeds), physical constants, and identifiers.

What does not: measured values. The DR25 transit parameters, the stellar
solution, the centroid, the transit shape and every other measurement are
written by a step into results/ and read from there through koi501/inputs.py, so
a measurement exists in exactly one place. steps/99_check.py fails the run if a
program carries a numeric literal that is not either here, mathematical (0, 1,
2, 3, 4, 0.5, 10 and 0.4 in magnitude arithmetic, unit conversions) or marked
on its line as display or technical.
"""

import math
from pathlib import Path

# ---------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
CACHE = RESULTS / ".cache"

# ---------------------------------------------------------------- the target
KOI = "K00501.01"
KEPID = 4951877
TCE_PLANET_NUMBER = 1
#: KIC catalogue position, deg: the centre of every cone search.
RA = 298.473210
DEC = 40.075871
#: The two neighbours the DR25 positional probability favours over the host.
NEIGHBOURS = (4951867, 4951861)

# ---------------------------------------------------------------- constants
G_SI = 6.67430e-11             # m^3 kg^-1 s^-2, CODATA 2018
M_SUN_KG = 1.98892e30
R_SUN_M = 6.957e8              # IAU 2015 nominal
R_EARTH_M = 6.3781e6           # IAU 2015 nominal, equatorial
M_JUP_KG = 1.898e27
AU_M = 1.495979e11
DAY_S = 86400.0
YEAR_D = 365.25
MINUTES_PER_DAY = 1440.0
#: Parts per million in one.
PPM = 1.0e6
TEFF_SUN_K = 5772.0            # IAU 2015 nominal
MBOL_SUN = 4.74                # IAU 2015
R_EARTH_PER_R_SUN = R_SUN_M / R_EARTH_M
#: Solar mean density, g cm^-3, and surface gravity, log10 cgs.
RHO_SUN_CGS = M_SUN_KG * 1e3 / (4.0 / 3.0 * math.pi * (R_SUN_M * 100.0) ** 3)
LOGG_SUN_CGS = math.log10(G_SI * M_SUN_KG / R_SUN_M ** 2 * 100.0)
#: Normal-distribution percentiles for a one-sigma interval.
ONE_SIGMA_PERCENTILES = (15.865, 50.0, 84.135)
#: Median absolute deviation to standard deviation, for a Gaussian.
MAD_TO_SIGMA = 1.4826

# ---------------------------------------------------------------- Kepler
KEPLER_PIXEL_AS = 3.98
#: Long-cadence integration, d (1765.5 s).
LONG_CADENCE_D = 0.0204340
#: Kepler barycentric Julian date offset: BKJD = BJD - 2454833.
BKJD_OFFSET = 2454833.0

# ---------------------------------------------------------------- queries
#: Radius of the catalogue cones around the target, arcsec.
CONE_RADIUS_AS = 30.0
PANSTARRS_CONE_AS = 12.0
#: Radius for picking a single star's own catalogue row, arcsec.
POINT_MATCH_AS = 3.0
#: Best-match radius for the KIC-to-Gaia cross-match of step 03, arcsec.
GAIA_XMATCH_AS = 2.0
#: Gaia G fainter than KIC r by more than this convicts the optical column, mag.
OPTICAL_COLUMN_WRONG_MAG = 1.0

# ---------------------------------------------------------------- steps 01-06
#: Blend search radius, arcsec. Armstrong et al. (2021) adopt the same.
BLEND_RADIUS = 25.0
#: r - J below this is not a stellar colour; even a white dwarf stays above it.
#: Normal stars run from about +0.3 (hot) to +4 and redder.
COLOUR_IMPOSSIBLE = -0.5
COLOUR_SUSPICIOUS = 0.0
#: Histogram of r - J drawn by step 08: range and bin width, mag.
COLOUR_HIST_RANGE = (-6.0, 6.0)
COLOUR_HIST_BIN = 0.25
#: A corrupted neighbour can only take the transit if the scene model can
#: plausibly place it there: close, and carrying real light.
CONFIG_SEP_MAX = 8.0
CONFIG_SHARE_MIN = 0.25
#: Grid used to show the two thresholds above are not what drives the result.
CONFIG_GRID_SEP = (6.0, 8.0, 10.0, 12.0, 15.0)
CONFIG_GRID_SHARE = (0.15, 0.25, 0.35)
#: Robovetter score at or below which step 05 counts an object as failed.
ROBOVETTER_FAIL_SCORE = 0.001
#: A difference-image offset above this many sigma counts as displaced (step 05).
DISPLACED_SIGMA = 2.0
#: Depth above which Bryson & Morton (2017) reject a star outright, ppm.
REJECTION_DEPTH_PPM = 3.0e6
#: The significance multiple used for upper limits and exclusions.
LIMIT_SIGMA = 3.0
#: Step 06 worst case: only this many of the sixteen quarters independent.
WORST_CASE_INDEPENDENT_QUARTERS = 4

# ---------------------------------------------------------------- step 07, host star
HOST_SEED = 20260904
HOST_DRAWS = 400_000
#: APOGEE formal errors widened for systematics: 33 K to 60 K, 0.026 to 0.035.
TEFF_WIDENED_K = 60.0
LOGG_WIDENED = 0.035
#: [Fe/H] uncertainty adopted for vespa's isochrone fit (APOGEE formal 0.007).
FEH_WIDENED = 0.10
#: Ks-band extinction prior, mag: normal(mean, sd) clipped to [0, max].
A_KS_PRIOR = (0.05, 0.03, 0.25)
#: RJCE extinction (Majewski et al. 2011): A_Ks = 0.918 (H - W2 - 0.08).
RJCE = (0.918, 0.08)
#: Ks bolometric correction, quadratic in (Teff - 5772)/1000 through the
#: solar-metallicity table of Casagrande & VandenBerg (2018).
BC_KS_COEFFS = (1.470, -0.284, -0.070)

# ---------------------------------------------------------------- step 09, velocities
#: Iterations solving the mass function for the total mass.
MASS_FUNCTION_ITERATIONS = 8
#: The lightest pair of stars that could eclipse (two 0.08 M_sun stars), M_sun.
LIGHTEST_ECLIPSING_PAIR_MSUN = 0.16
#: Companions an eclipsing binary on the target could have, M_sun: from the
#: smallest hydrogen-burning star to a solar twin.
EB_COMPANION_RANGE_MSUN = (0.08, 1.0)

# ---------------------------------------------------------------- step 10, ephemeris
#: Coughlin et al. (2014) statistical thresholds.
EPHEMERIS_SIGMA_P, EPHEMERIS_SIGMA_T = 3.5, 2.0
#: Kepler field-of-view centre (19h22m40s, +44d30m), deg.
KEPLER_FOV_CENTRE = (290.66667, 44.5)
#: Kirk et al. (2016) tabulate BJD - 2400000; subtracting this gives BKJD.
KIRK_BJD_OFFSET = 54833.0
#: Coughlin et al. (2014) physical path: PRF reach and antipodal tolerance, arcsec.
PRF_REACH_AS, ANTIPODAL_AS = 50.0, 50.0
#: Coughlin et al. (2014) reach formula d = 50 sqrt(scale 10^(-0.4 Kp) + 1).
PRF_REACH_FLUX_SCALE = 1.0e6
#: Magnitude assumed when neither object has one.
DEFAULT_KEPMAG = 20.0

# ---------------------------------------------------------------- step 12, contrast curves
#: ExoFOP follow-up frames: band, file id, name, header seeing in arcsec, sha256.
CONTRAST_FRAMES = (
    ("U", 121811, "501Ic-SH20110625u.fits", 1.46,
     "803b84a02ced8a20c200f793b31467585d316e5ded775f5c36d3459cef2bb5aa"),
    ("B", 124500, "501Ic-SH20110625b.fits", 1.44,
     "f13a64db63a4b02433abda8fe176b47359ffa7f20a22190839925f3c7b9ae3f0"),
    ("V", 127189, "501Ic-SH20110625v.fits", 1.31,
     "e558ba2aab0c53ae863d342bc8da98bed760eb95108a28e7e4f21c48f8cdba0a"),
    ("J", 132444, "501Ic-dc20140206UKJ.fits", 0.66,
     "2ccdb1f7e393f3aa9755e37bb18cb060dcc4152b51395b3b52d13294e5857518"),
)
CONTRAST_SIGMA = 5.0
#: Pixels beyond this fraction of the frame's short side estimate the background.
CONTRAST_BACKGROUND_FRACTION = 0.35
#: Radial profile for the FWHM: maximum radius, samples, ring width (pixels).
CONTRAST_PROFILE = (12.0, 200, 0.5)
#: Centroid box: at least this many pixels, or this many arcsec.
CONTRAST_CENTROID_BOX = (4, 6.0)
#: Rings with at most this many pixels give no median.
CONTRAST_MIN_RING_PIXELS = 4
#: Separations reach this fraction of the frame width.
CONTRAST_OUTER_FRACTION = 0.45
#: Apertures per ring, at least; fewer measured apertures than this give no limit.
CONTRAST_MIN_APERTURES = 8
CONTRAST_MIN_FLUXES = 6
#: Separation at which the text quotes the J limit, arcsec.
CONTRAST_QUOTE_SEP_AS = 1.5
#: TRICERATOPS reads the curve out to this separation, arcsec.
CONTRAST_TRICERATOPS_MAX_AS = 10.0

# ---------------------------------------------------------------- step 13, positional probability
POSPROB_RADIUS_AS = 20.0
#: Probability above which the target is reported as "> 0.999".
POSPROB_REPORT_LEVEL = 0.999

# ---------------------------------------------------------------- step 14, blend bounds
SHAPE_BOOTSTRAPS = 600
SHAPE_SEED = 20260904
#: Points within this phase of each transit are fitted.
SHAPE_PHASE_WINDOW = 0.06
#: A transit window needs this many points, and this many out of transit.
SHAPE_MIN_POINTS, SHAPE_MIN_OOT = 20, 10
#: Out of transit means beyond this many durations from mid-transit.
OOT_DURATIONS = 0.75
#: Exposure supersampling of the trapezoid.
SHAPE_SUPERSAMPLE = 9
#: Trapezoid (epoch offset d, depth, duration d, ingress d): lower and upper bounds,
#: the ingress starting guess, and evaluation limits for the fit and the bootstrap.
SHAPE_BOUNDS = ((-0.05, 1e-5, 0.1, 1e-4), (0.05, 0.05, 0.6, 0.25))
SHAPE_INGRESS_GUESS_D = 0.05
SHAPE_MAXFEV = (40000, 20000)
#: The shape bound is this percentile of the bootstrap, the permissive side.
SHAPE_BOUND_PERCENTILE = 95.0
#: Percentiles summarising the bootstrap.
SHAPE_SUMMARY_PERCENTILES = (16, 50, 84)
#: A catalogued star further than this many standard errors of the combined
#: difference-image offset from the target cannot host the signal.
CENTROID_EXCLUSION_SIGMA = 3.0

# ---------------------------------------------------------------- step 15, literature
#: Tables published only in arXiv sources: pinned version, file, SHA-256 when first read.
ARXIV_TABLES = {
    "armstrong2021": ("2008.10516v1", "TableA1.csv",
                      "cf5a342f55a06a7333f9fb01a35234696ece02a8fdce0f36748f8c2ce0b5aa03"),
    "valizadegan2022": ("2111.10009v2", "Tables/table_13.csv",
                        "4f25d4f656736539cad04fc943bdaa072b57f025fc98c72e252ab7a36e2c9e84"),
}
#: Berger et al. (2018) evolutionary-state flag: subgiant, red giant.
BERGER_SUBGIANT, BERGER_GIANT = 1, 2
#: Mass provenances that count as a measured mass in pscomppars.
MEASURED_MASS_PROVENANCE = ("Mass", "Msini")

# ---------------------------------------------------------------- step 16, instrumental share
#: The cell of Thompson et al. (2018) holding KOI-501.01: period range (d) and
#: the multiple event statistic above which it counts.
INSTRUMENTAL_PERIOD_RANGE_D, INSTRUMENTAL_MES_MIN = (10.0, 40.0), 30.0
#: One-sided confidence of the upper limit on a count of zero.
INSTRUMENTAL_CONFIDENCE = 0.95
#: Synthetic runs averaged into one realisation of the noise (inverted, scrambled).
INSTRUMENTAL_RUNS = 2
#: The bound the paper quotes: the widest, each synthetic run counted as its own
#: realisation of the noise and the DR25 candidates as the denominator.
INSTRUMENTAL_QUOTED_BOUND = "unaveraged_over_candidates"

# ---------------------------------------------------------------- fit/transit_fit.py
TRANSIT_FIT_SEED = 20260913
TRANSIT_FIT_SAMPLER = (40, 3000, 6000)          # walkers, burn-in, steps
#: Windows of this many durations either side of each transit, at least this
#: many points and this many out of transit; a window "covers" the transit with
#: this many points inside this fraction of a duration.
TRANSIT_FIT_WINDOW_DURATIONS = 4.0
TRANSIT_FIT_MIN_POINTS, TRANSIT_FIT_MIN_OOT = 20, 12
TRANSIT_FIT_COVER = (5, 0.4)
TRANSIT_FIT_SUPERSAMPLE = 15
#: Sensitivity run with the limb darkening held fixed (u1, u2).
TRANSIT_FIT_FIXED_LD = (0.36, 0.29)
#: Parameter limits: Rp/R*, density (solar), |epoch shift| d, error scale.
TRANSIT_FIT_LIMITS = dict(ror=(0.005, 0.1), rho=(0.01, 5.0), dt0=0.05, scale=(0.3, 3.0))
#: Walker start (Rp/R*, b, rho, dt0, q1, q2, scale) and its spread.
TRANSIT_FIT_START = (0.0216, 0.3, 0.236, 0.0, 0.35, 0.3, 1.0)
TRANSIT_FIT_START_SPREAD = (5e-4, 0.05, 0.01, 1e-3, 0.03, 0.03, 0.02)
#: batman starting geometry before the first likelihood call.
BATMAN_START = dict(rp=0.02, a=22.0, inc=89.9, u=(0.4, 0.2))

# ---------------------------------------------------------------- fit/bayes_factor.py
BAYES_SEED = 20260915
BAYES_WINDOW_DURATIONS = 2.5
BAYES_MIN_POINTS, BAYES_MIN_OOT, BAYES_MIN_IN = 40, 20, 3
BAYES_CLIP_SIGMA = 5.0
#: Quadratic limb darkening for Teff 5764 K, log g 4.05, [Fe/H] 0, Kepler band.
BAYES_LD = (0.42, 0.24)
BAYES_SUPERSAMPLE = 7
BAYES_LIVE_POINTS, BAYES_FRAC_REMAIN = 400, 0.3
#: Control epoch offset, periods; points within this many durations of a real
#: transit are removed from the control.
BAYES_CONTROL_PHASE = 0.371
BAYES_CONTROL_EXCLUSION_DURATIONS = 1.5
#: Windows are laid on one time axis this far apart, d.
BAYES_AXIS_GAP_D = 100.0
BAYES_PRIORS = {
    "t0": (-0.05, 0.05),          # d
    "log_rp": (-3.0, -1.0),       # log10 Rp/R*
    "aor": (5.0, 60.0),           # a/R*
    "b": (0.0, 0.99),
    "log_sigma": (-6.0, -2.0),    # log10 GP amplitude, relative flux
    "log_rho": (-1.5, 1.5),       # log10 GP timescale, d
    "log_jit": (-5.0, -3.0),      # log10 white-noise jitter, relative flux
}
#: Narrow priors: this many posterior standard deviations, at least this wide,
#: a/R* at least this.
BAYES_NARROW_SIGMA, BAYES_NARROW_MIN_WIDTH, BAYES_NARROW_MIN_AOR = 5.0, 1e-4, 3.0
#: Log-likelihood returned for a failed evaluation.
BAYES_BAD_LOGL = -1e30
BATMAN_BAYES_START = dict(rp=0.02, a=20.0, inc=89.0)

# ---------------------------------------------------------------- fit/lightcurve_checks.py
LC_CHECK_SEED = 20260904
#: Detrending: a gap longer than this (d) starts a segment; segments need this
#: many points; polynomial degree (short, long) switching at this many points;
#: clipping iterations and threshold (robust sigma); points needed beyond the degree.
LC_SEGMENT_GAP_D, LC_SEGMENT_MIN_POINTS = 0.5, 40
LC_DETREND_DEGREE, LC_DETREND_LONG_POINTS = (3, 5), 400
LC_DETREND_ITERATIONS, LC_DETREND_CLIP_SIGMA, LC_DETREND_SPARE = 6, 4.0, 6
#: Box depth: in-transit within this fraction of the half-duration; baseline
#: from this many half-durations out to this phase; minimum points in and out.
LC_BOX_IN, LC_BOX_OUT, LC_BOX_OUT_PHASE = 0.7, 1.6, 0.3
LC_BOX_MIN_IN, LC_BOX_MIN_OUT = 5, 50
#: The null distribution: this many wrong epochs, uniform in this phase range.
LC_NULL_EPOCHS, LC_NULL_PHASES = 20_000, (0.1, 0.9)
#: Single transits: baseline out to this phase; minimum points in and out.
LC_SINGLE_OUT_PHASE, LC_SINGLE_MIN_IN, LC_SINGLE_MIN_OUT = 0.20, 4, 30
#: Circular-orbit secondary eclipse phase.
LC_SECONDARY_PHASE = 0.5
#: An M0V dwarf, mass and radius in solar units (Pecaut & Mamajek 2013).
M0V_MASS_RADIUS = (0.5, 0.47)
#: Transit times: window half-width (d), minimum points, transits searched
#: within this phase; timing grid half-width (d) and size; parabola half-width
#: (grid steps); largest accepted offset and error (d).
TIMING_WINDOW_D, TIMING_MIN_POINTS, TIMING_SEARCH_PHASE = 0.9, 25, 0.06
TIMING_GRID_D, TIMING_GRID_SIZE, TIMING_PARABOLA = 0.25, 501, 3
TIMING_MAX_OFFSET_D, TIMING_MAX_ERROR_D = 0.2, 0.2
TIMING_MIN_OOT, TIMING_MIN_IN, TIMING_IN_FRACTION = 12, 3, 0.4
#: Injection and recovery: injections, tries, clearance from real transits
#: (durations), injected offset range (d), recovered outliers beyond (d).
TIMING_INJECTIONS, TIMING_TRIES, TIMING_CLEARANCE = 300, 4000, 1.2
TIMING_INJECTED_OFFSET_D, TIMING_RECOVERY_CLIP_D = 0.05, 0.15

# ---------------------------------------------------------------- fit/diff_images.py
#: Target pixel files read per star (long cadence, in archive order).
DIFF_MAX_TPF = 17
#: A stamp needs this many usable pixels; a cadence is kept when this fraction
#: of them is finite; a quarter needs this many cadences.
DIFF_MIN_PIXELS, DIFF_FINITE_FRACTION, DIFF_MIN_CADENCES = 9, 0.9, 50
#: In transit: within this fraction of a duration of mid-transit; comparison
#: frames between these multiples of the duration.
DIFF_IN_FRACTION, DIFF_OUT_RANGE = 0.5, (0.8, 3.0)
#: Minimum in- and out-of-transit cadences per quarter, and per single transit.
DIFF_MIN_IN, DIFF_MIN_OUT = 3, 12
DIFF_MIN_IN_TRANSIT, DIFF_MIN_OUT_TRANSIT = 2, 4
#: Error floor of a difference-image pixel.
DIFF_ERROR_FLOOR = 1e-9

# ---------------------------------------------------------------- fit/prf_centroids.py
#: A difference-image panel is fitted only above this peak signal-to-noise, and
#: kept only if the fitted PRF correlates with the data at least this well.
PRF_MIN_PANEL_SNR = 4.0
PRF_MIN_QUALITY = 0.5
#: An object needs this many good panels.
PRF_MIN_PANELS = 3
#: A PRF centred at column0 + c peaks at array index c - 1 (verified by round trip).
PRF_INDEX_OFFSET = 1.0
#: Controls DR25 places on their target at better than this many standard errors
#: calibrate the systematic.
PRF_ON_TARGET_SIGMA = 1.0
#: Controls above this median panel signal-to-noise give the comparison with DR25.
PRF_COMPARISON_SNR = 6.0

# ---------------------------------------------------------------- fit/diffimage_figure.py
#: Aligned stack: a square grid this many pixels across, target at its centre.
FIGURE_GRID = 7
#: The one-season panel: the CCD module whose quarters share a stamp shape.
FIGURE_SEASON_MODULE = 6
#: A Gaia source this close to the target (arcsec) is the target itself.
FIGURE_TARGET_MATCH_AS = 0.5

# ---------------------------------------------------------------- fit/pixel_scene.py
PIXEL_GAIA_RADIUS_AS = 40.0
#: A stamp pairs with a report quarter when their KIC reference pixels agree to this.
PIXEL_PAIRING_TOLERANCE_PX = 0.1
#: least_squares x_scale: position (px), flux fraction of the total, background.
PIXEL_XSCALE = (0.1, 0.1, 10.0)
#: Scene fit: starting neighbour flux as a fraction of the total, position scale.
PIXEL_SCENE_START_FRACTION = 0.05
PIXEL_SCENE_XSCALE_PX = 0.05
PIXEL_WEIGHTINGS = ("uniform", "photon")

# ---------------------------------------------------------------- fpp/triceratops_fpp.py
TRICERATOPS_RUNS, TRICERATOPS_DRAWS, TRICERATOPS_SEED = 20, 100_000, 20260904
TRICERATOPS_SEARCH_RADIUS_PX = 10
#: Giacalone et al. (2021) validation bars.
TRICERATOPS_FPP_BAR, TRICERATOPS_NFPP_BAR = 0.015, 0.001
#: Light-curve fold: window in durations, bins, flatten window (cadences),
#: outlier clip (upper, lower sigma), out-of-transit bins beyond this many durations.
TRICERATOPS_FOLD = dict(window=3.0, bins=100, flatten=401, clip=(4, 20), oot=1.5)

# ---------------------------------------------------------------- fpp/own_fpp.py
OWN_SEED, OWN_RUN_SEED_OFFSET, OWN_REPEATS, OWN_DRAWS = 20260830, 50100, 20, 2_000_000
#: Raghavan et al. (2010): binary fraction, log-normal period distribution
#: in log10(P/d), triple fraction.
RAGHAVAN = dict(binary=0.44, logp_mu=5.03, logp_sigma=2.28, triple=0.11)
#: Kepler stars searched.
KEPLER_TARGETS = 190_000
#: Background blends within four Kepler pixels.
OWN_BLEND_RADIUS_PIXELS = 4.0
#: Depth, in standard errors, at which a secondary counts as detected.
SECONDARY_DETECTION_SIGMA = 7.0
#: Eccentricity Beta distributions: Kipping (2013) for planets; broader for binaries.
ECC_PLANET, ECC_BINARY, ECC_MAX = (0.867, 3.03), (1.0, 2.0), 0.95
#: Spread of the DR25 odd-even statistic among confirmed planets of similar
#: strength and period (3% exceed 5.27).
ODD_EVEN_SCALE = 2.43
#: Flat allowances: uncatalogued contaminating stars (Morton et al. 2016) and an
#: instrumental origin (curated inverted and scrambled runs, Thompson et al. 2018).
P_CONTAM_RESIDUAL, P_INSTRUMENT = 0.0006, 0.0014
#: Smallest hydrogen-burning star, R_sun (Dieterich et al. 2014).
R_STELLAR_MIN = 0.086
#: Companion mass ratios: smallest considered, and smallest bound triple member.
Q_MIN, Q_MIN_TRIPLE = 0.06, 0.10
#: Credit on the triple scenarios for Gaia astrometry showing no companion.
RUWE_CREDIT = 0.33
#: The centroid likelihood is never taken below this many sigma.
CENTROID_FLOOR_SIGMA = 5.0
#: Main-sequence relations: R = M^a, L = M^b; secondary depth ratio q^c.
MS_RADIUS_EXPONENT, MS_LUMINOSITY_EXPONENT, MS_SECONDARY_EXPONENT = 0.92, 4.5, 2.66
#: Period window for the priors: P/f to P*f.
PERIOD_WINDOW_FACTOR = 1.5
#: Candidate radii sampled for the planet scenario, R_earth.
PLANET_RADIUS_RANGE = (0.3, 25.0)
#: Numerical settings: centroid integration grid, floors.
OWN_CENTROID_GRID = 4000
OWN_FLOORS = dict(k=1e-9, sini=1e-12, q=1e-6, blend=1e-300)
#: Run B: likelihood floor for an eclipsing binary on the target.
EB_FLOOR = 1e-30

# ---------------------------------------------------------------- fpp/vespa
#: Velocity contrast constraint from the single-lined APOGEE spectra: a second
#: star within this velocity separation (km/s) and magnitude difference.
VESPA_VCC = (25.0, 3.0)
#: Gaia G uncertainty adopted for the isochrone fit, mag.
VESPA_G_ERROR = 0.02
#: vespa reads the contrast curve out to this separation, arcsec.
CONTRAST_VESPA_MAX_AS = 12.0
