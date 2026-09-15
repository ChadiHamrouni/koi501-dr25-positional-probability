# KOI-501.01: DR25 positional probability and validation

Code and data for *A Catalogue Error in the DR25 Positional Probability for
KOI-501.01, and a Validated Sub-Saturn*.

Morton et al. (2016) computed a false-positive probability of
(2.0 ± 0.7) × 10⁻⁴ for KOI-501.01 but did not validate it, because the DR25
astrophysical positional probability of Bryson & Morton (2017) gives the host
a relative probability of 3.0 × 10⁻⁴ and gives 0.950 to KIC 4951867, 6.7 arcsec
away. The paper shows that this value was computed from Kepler Input Catalog
(KIC) magnitudes that are 2.96 and 2.03 mag too bright for the two neighbours,
that the defect recurs across the DR25 objects of interest but changes a verdict
only for KOI-501.01, that the Data Validation difference images place the
transit source on the target, and that APOGEE DR17 velocities bound the
transiting body below 3.3 M_Jup.

This repository recomputes those results from the public archives, together
with the fit of Kepler's own pixels in Section 2.2, the positional probability
of Section 7, the blend bounds and Bayes factor of Section 6, and the three
false-positive probabilities of Section 7.

## Running

```bash
python -m venv .venv && . .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python run_all.py            # all steps
python run_all.py 04 06      # selected steps
```

Queries are cached in `results/.cache/`; with the cache present the full run
takes under a minute, most of it the 600 transit-shape fits of step 14, and a
first run spends a few minutes downloading. Requests to the archives
are retried automatically if a connection drops.

The calculations that need heavier packages sit outside `run_all.py`, each
listed with its requirements:

| Program | Paper | Run time | Requirements | Output |
|---|---|---|---|---|
| `fit/diff_images.py` | §3, the difference images | minutes per star (downloads) | `fit/requirements.txt` | `results/.cache/diffimages/` (compare with `data/diffimages/`) |
| `fit/prf_centroids.py` | §3.1, PRF fit and centroid systematic | a few minutes | `fit/requirements.txt` | `results/prf_centroids.json`, `results/prf_centroids_controls.csv` |
| `fit/diffimage_figure.py` | Figure 1 | under a minute | `fit/requirements.txt` | `figures/K00501.01_diffimage.png`, `data/figure_diffimage/` |
| `fit/transit_fit.py` | §4, radius ratio | about 1 hour | `fit/requirements.txt` | `data/transit_fit/transit_fit.json` |
| `fit/pixel_scene.py` | §2.2, Kepler's own pixels | under a minute | `fit/requirements.txt` | `results/pixel_scene.json` |
| `fit/lightcurve_checks.py` | §6, Table 4 rows marked "here" | a few minutes | `fit/requirements.txt` | `results/lightcurve_checks.json` |
| `fit/bayes_factor.py` | §6, Table 4, Bayes factor | 40 to 70 minutes | `fit/requirements.txt` | `results/bayes_factor.json` |
| `fpp/vespa/run_vespa.sh` | §7, vespa | about 15 minutes | Python 3.8 and vespa 0.6, listed in the script | `fpp/vespa/results*.txt` |
| `fpp/vespa/make_inputs.py` | §7, vespa inputs | seconds | none | compares or rewrites `fpp/vespa/*.ini`, `lc.csv`, `UKIRT_J.cc` |
| `fpp/triceratops_fpp.py` | §7, TRICERATOPS | about 15 minutes | `fpp/requirements.txt` | `results/fpp_triceratops*` |
| `fpp/own_fpp.py` | §7 and Appendix A | about 10 minutes | `fpp/requirements.txt` | `results/fpp_own.json` |
| `fpp/summarise_runs.py` | §7, the published runs | seconds | pandas, PyTables | `results/fpp_published_runs.json` |

Order: steps 00-07 and 09-12; `fit/prf_centroids.py`; steps 13-16; then
`fit/transit_fit.py` (reads step 07), `fit/lightcurve_checks.py` (reads the
transit fit and step 14) and the rest; step 99 last.

## One source for every value

Every fixed value (a threshold, a seed, a physical constant, a prior) is
defined once, in `koi501/config.py`. Every measured value (a period, a stellar
radius, a centroid, a velocity limit) is read from the output of the step that
fetched or computed it, through `koi501/inputs.py`, and each output records the
values it was made from under `measured_inputs`.

Step 99 (`steps/99_check.py`) stops the run if

1. a recorded input no longer equals the value its source step now writes;
2. vespa's committed input files differ from those regenerated from the step
   outputs (`fpp/vespa/make_inputs.py`);
3. any program outside `config.py` contains a typed number that is not a unit
   conversion, a mathematical constant, an array position or a line marked
   `# literal: <reason>` (plot layout, numerical safeguards);
4. a number the paper prints, listed with its source in `claims.csv`, is not what
   the output now gives (`koi501/claims.py`).

Numbers taken from published tables or archives are quoted, not reproduced;
every number the paper derives is produced by a program here.

## Requirements

| | |
|---|---|
| Language | Python 3.12.7; exact package versions in `requirements.txt` (steps), `fit/requirements.txt` and `fpp/requirements.txt` |
| Tested on | Windows 11 Pro (10.0.26200); vespa under Ubuntu on WSL |
| Network | HTTPS access to exoplanetarchive.ipac.caltech.edu, vizier.cds.unistra.fr, cdsxmatch.u-strasbg.fr, archive.stsci.edu, exofop.ipac.caltech.edu; `fpp/triceratops_fpp.py` also uses mast.stsci.edu |
| Disk | about 11 MB of cached queries, 1 MB of outputs |
| Inputs | the archives above and the committed files in `data/` and `fpp/` |
| Outputs | `results/*.json`, `results/*.csv`, `figures/colour_pathology.{pdf,png}`, `figures/K00501.01_diffimage.png` |

`run_all.py` executes the scripts in `steps/` in order. Each one queries its
sources, validates every returned record, performs one part of the analysis,
prints the quantities the paper quotes under the relevant section heading, and
writes its full output to `results/<step>.json`. The script reports values
only; it does not compare them with the manuscript.

## Claims and where they are reproduced

| Paper | Claim | Step or program | Data |
|---|---|---|---|
| §2.2, Table 1 | KIC gives KIC 4951867 Kp = 14.830 and KIC 4951861 Kp = 15.485; Gaia DR3 gives G = 17.787 and 17.513 | 01 | KIC, Gaia DR3 |
| §2.2 | KIC 4951867 has r − J = −1.49 in the KIC's own columns; the host has +1.17 | 01 | KIC |
| §2.2 | 2MASS J and Pan-STARRS DR1 r agree with Gaia, not with the KIC; the host's KIC magnitude is correct to 0.02 mag | 01 | 2MASS, Pan-STARRS DR1 |
| §2.2 | No Gaia source brighter than G = 16.23 other than the host within 25 arcsec | 01 | Gaia DR3 |
| §2.2 | Kepler's own pixels: fitted on the sixteen out-of-transit images with the 60 Gaia DR3 sources within 40 arcsec, KIC 4951867 carries 0.052 and KIC 4951861 0.088 of the target's flux (0.050 and 0.091 with photon-noise weights), against 0.818 and 0.448 from the KIC and 0.053 and 0.068 from Gaia; four CCD modules; no pixel off by more than 2.7% of the peak. The report's out-of-transit light centre lies 0.187 arcsec from the catalogue position; a single source fitted to the model scene lands 0.15 arcsec away with Gaia magnitudes and 2.0 arcsec with KIC magnitudes, beyond the largest observed offset, 0.260 arcsec, in all 16 quarters; fitted to the stored images it reproduces the report to 0.19 arcsec rms | `fit/pixel_scene.py` | `data/figure_diffimage/`, `data/dv/`, step 01 |
| §4 | Host RUWE 0.970, astrometric excess noise 0, multi-peak fraction 0, no non-single-star entry | 01 | Gaia DR3 |
| §2.3 | The archive's positional-probability record lists Kp = 14.830 and 15.480 with provenance `KIC`, and required depths of 69,203 and 580,986 ppm | 04 | `data/koiapp/` |
| §2.4 | With Gaia magnitudes the required eclipses become 105% and 378%; the second exceeds the 3 × 10⁶ ppm rejection limit | 04 | `data/koiapp/`, Gaia DR3 |
| §2.3, §2.4 | An earlier release of the table gives the host 2.9 × 10⁻⁴ and KIC 4951867 0.96 from the same KIC magnitudes; corrected eclipses 89% and 226% | 04 | MAST bulk file, Gaia DR3 |
| §2.6 | Of 29,673 KIC neighbours within 25 arcsec of the 8,054 DR25 objects of interest, 135 beside 124 objects have r − J < −0.5; 35 hosts have the same defect | 02 | DR25 KOI table, KIC |
| §2.6 | Dispositions of the 124: 80 FALSE POSITIVE, 31 CONFIRMED, 13 CANDIDATE | 02 | DR25 KOI table |
| §2.6 | The 135 are 106 stars; 92 have a Gaia counterpart and all 92 have G more than 1 mag fainter than KIC r (median G − r = +3.15) | 03 | KIC, Gaia DR3 |
| §2.6 | Of the twelve other CANDIDATE objects, none has its positional probability taken by the impossible-colour neighbour | 04 | `data/koiapp/` |
| §2.6, Table 2 | Nine neighbours lie within 8 arcsec and carry ≥ 25% of the pair's light; eight hosts are FALSE POSITIVE, KOI-501.01 is the ninth | 05 | DR25 KOI table |
| §2.6 | Over 15 threshold combinations, 12 return KOI-501.01 as the only candidate, 3 return none, none returns another | 05 | DR25 KOI table |
| §3 | Mean difference-image offset (+0.28, −0.47) arcsec, error on the mean 0.22 arcsec per axis | 06 | `data/dv/` |
| §3 | KIC 4951867 lies 30.7σ and KIC 4951861 38.6σ from the source; 99.8σ, 30.5σ (the report's robust-mean errors), 15.4σ and 7.7σ under the alternative error models; all 16 quarters closer to the target, the weakest still 12.5σ from KIC 4951867 | 06 | `data/dv/`, KIC |
| §3.1 | This work's PRF fit to the sixteen difference images; the per-axis systematic, 0.321 arcsec, from 73 on-target controls among 243 | `fit/prf_centroids.py` | `data/diffimages/` |
| §3.1 | This work's fit combined with DR25's offset as vectors: 0.299 ± 0.270 arcsec, 1.11σ | 13 | DR25 KOI table, `results/prf_centroids.json` |
| §3.1 | Flux-weighted centroid: 1.037 ± 0.340 arcsec, 3.05σ; 2.22σ with the per-axis systematic | 11 | DR25 TCE table, `results/prf_centroids.json` |
| Figure 1 | Difference images, all quarters aligned and the module-6 season | `fit/diffimage_figure.py` | `data/diffimages/`, Gaia DR3 |
| §5.1 | 19 APOGEE visits; K = 12 ± 63 m s⁻¹; 3σ limit 199 m s⁻¹, i.e. M < 3.3 M_Jup; drift limit and the bound pairs and eclipsing binaries it excludes | 09 | APOGEE DR17, step 07 |
| §4, Table 3 | Host radius 1.746 ± 0.061 R☉ from 2MASS Ks and the Gaia distance at the APOGEE temperature; mass 1.26 ± 0.13 M☉; density 0.236 ± 0.021 ρ☉; orbital distance, equilibrium temperature and insolation | 07 | 2MASS, AllWISE, Gaia DR3 distances, APOGEE DR17 |
| §4, Table 3 | Rp/R★ = 0.0218 ± 0.0005; planet 2.72, 4.06 and 4.16 ± 0.17 R⊕ on the DR25, Berger et al. 2020 and adopted stellar radii | 07 | `data/transit_fit/`, Mathur et al. 2017, Berger et al. 2020 |
| §5.2 | UKIRT J contrast curve: FWHM 0.70 arcsec, inner working angle 0.80 arcsec, 5σ ΔJ = 5.1 mag at 1.5 arcsec and 7.21 mag at 2.9 arcsec; U, B, V reach 5.68, 5.88 and 6.70 mag | 12 | ExoFOP follow-up images |
| §6 | 8053 objects of interest and 2876 eclipsing binaries compared; one statistical match (K01278.01), no physical path, no match under Coughlin et al. 2014 | 10 | DR25 KOI table, Kirk et al. 2016, MAST field of view |
| §6, Table 4 | DR25 rows: MES 33.06, bootstrap false-alarm probability 1.4 × 10⁻²³⁸, rolling bands 0 of 52, weak secondary 44.9 ± 15.1 ppm, odd–even 0.111, ghost core/halo 17.96/6.23, score 0.998, no flags | 11 | DR25 TCE and KOI tables |
| §6, Table 4 | A blend can be at most 4.41 mag fainter than the target by the transit shape, 8.10 by brightness; none of the 26 other Gaia DR3 sources within 25 arcsec survives both bounds and the centroid | 14 | `data/lightcurve/`, Gaia DR3, step 13 |
| §6, Table 4 | Transit plus Gaussian process against the process alone: ΔlnZ = 392.0 ± 0.5; 397.3 with priors narrowed to five posterior standard deviations (397.1 in the unseeded run of `results/paper_runs/`); −0.9 at an epoch with no transit | `fit/bayes_factor.py` | `data/lightcurve/` |
| §7, Eq. 1 | 22 Gaia DR3 sources within 20 arcsec, all bright enough for the depth; centroid (+0.225 ± 0.143, −0.108 ± 0.145) arcsec plus 0.321 arcsec per axis; P > 0.999 on the target; nearest alternative 19.0 standard errors away, G = 20.83 at 6.91 arcsec; 24.5 on the §3.1 centroid | 13 | Gaia DR3, DR25 KOI table, `results/prf_centroids.json` |
| §1, §4 | ExoMiner's two scores imply a prior of 2.6 × 10⁻⁴; of the Kepler planets confirmed without a measured mass whose hosts Berger et al. 2018 classify, 2080, 252 orbit subgiants | 15 | arXiv tables, VizieR, NASA Exoplanet Archive |
| §6, Table 4 | Instrumental share: no synthetic detection of Thompson et al. 2018 passes in the cell (0 of 22; 1229 real detections), so the share is given as a 95% upper bound | 16 | VizieR, DR25 TCE and KOI tables |
| §6, Table 4 | Deepest false dip in 20,000 wrong epochs 73 ppm against 516; 50 of 50 transits deeper than zero; −0.1 ± 15.0 ppm at phase 0.5; odd/even 520 ± 22 / 512 ± 21 ppm and first/second half 512 ± 22 / 520 ± 20 ppm (0.26σ); transit times scatter 31 min against 40 min precision (χ²ν 0.62), odd and even times 1.07σ apart; an M0V host would need e > 0.76; CROWDSAP 0.922–0.993 | `fit/lightcurve_checks.py` | `data/lightcurve/`, `data/transit_fit/`, MAST |
| §7 | vespa: FPP 1.7 × 10⁻⁴, standard deviation 0.9 × 10⁻⁴ over 100 bootstrap resamplings; populations of 20,000; isochrone radius 1.76 R☉ | `fpp/vespa/`, `fpp/summarise_runs.py` | committed inputs; TRILEGAL |
| §7 | TRICERATOPS, 20 runs of 10⁵ draws: mean FPP 7.2 × 10⁻⁴, NFPP 0.00000; five runs exactly zero, median 2 × 10⁻¹², largest 0.0143 (`results/paper_runs/`). Seeded rerun: mean 9.6 × 10⁻⁴, NFPP 0.00000, largest 0.0127, all 20 runs below both bars | `fpp/triceratops_fpp.py` | MAST, `fpp/data/`, steps 01 and 12 |
| §7, App. A | This work's calculation: run A 2.2 × 10⁻³, run B 2.1 × 10⁻³, of which 2.00 × 10⁻³ are the flat allowances | `fpp/own_fpp.py` | `fpp/data/`, steps 12–14 |

Step 08 draws `figures/colour_pathology.pdf`, a supporting figure of the step 02
and 05 results that does not appear in the paper. The left panel is the r − J
distribution of all tested neighbours with the −0.5 cut; the right panel places
each impossible-colour neighbour by separation and light share, with the
configuration of Table 2 shaded.

Not reproduced here:

- the choice of the 243 control objects of `fit/prf_centroids.py`, a seeded
  random draw (seed 20260829) from the stars whose pixel data the analysis
  pipeline had downloaded; the list and their difference images are committed;
- the TRILEGAL simulations of the field, which are committed as outputs
  (`fpp/data/`) rather than regenerated;
- vespa's stellar models and simulated populations, which `fpp/vespa/run_vespa.sh`
  regenerates as new random realisations (the paper's populations are 139 MB and
  are not included).

Every sampler is seeded (`koi501/config.py`), so a rerun in the same
environment repeats the committed outputs. Step 14's bootstrap fits differ in
the fifth significant figure between numpy 2.1 and numpy 2.5, which reaches
`fpp/own_fpp.py` at one part in 10,000. The paper's TRICERATOPS runs and its
first Bayes-factor run were not seeded; they are kept in `results/paper_runs/`
and agree with the seeded reruns within their scatter.

## Sources

| Source | Table or catalogue | Access | Steps |
|---|---|---|---|
| NASA Exoplanet Archive | `q1_q17_dr25_koi` (DR25 objects of interest) | TAP | 02, 04, 05, 13 |
| NASA Exoplanet Archive | `koiapp` (DR25 positional probabilities) | manual export, committed | 04 |
| MAST | `kepler_koiapp.txt.gz` (earlier positional-probability release) | bulk file | 04 |
| Kepler DV report | per-quarter difference-image centroids, KIC 4951877 | manual extraction, committed | 06 |
| Kepler DV report | per-quarter out-of-transit and KIC reference centroids, KIC 4951877 | manual extraction, committed | `fit/pixel_scene.py` |
| MAST | Kepler PRF calibration files (2011) for the four CCD modules | download on first use, cached in `results/.cache/` | `fit/pixel_scene.py` |
| CDS VizieR | Gaia DR3 `I/355/gaiadr3` (the CDS copy of `gaia_source`) | cone search | 01, 13, 14 |
| CDS VizieR | KIC `V/133/kic`, 2MASS `II/246`, Pan-STARRS DR1 `II/349`, APOGEE DR17 `III/286/allvis` | cone search | 01, 06, 07 |
| CDS X-Match | DR25 objects of interest × KIC; flagged neighbours × Gaia DR3 `I/355/gaiadr3` | bulk cross-match | 02, 03 |
| CDS VizieR | 2MASS `II/246`, AllWISE `II/328`, Bailer-Jones distances `I/352`, APOGEE DR17 `III/286/catalog`, Mathur et al. 2017 `J/ApJS/229/30`, Berger et al. 2020 `J/AJ/159/280`, Kirk et al. 2016 `J/AJ/151/68` | cone or identifier query | 09, 10 |
| NASA Exoplanet Archive | `q1_q17_dr25_tce` (DR25 threshold-crossing events) | TAP | 11 |
| MAST | Kepler field-of-view service (CCD modules per season) | JSON query | 10 |
| ExoFOP | Kepler follow-up images of KIC 4951877 (TIC 170740547): WIYN 0.9 m U, B, V (file ids 121811, 124500, 127189) and UKIRT WFCAM J (132444) | file download, SHA-256 checked | 12 |
| MAST | Kepler long-cadence PDCSAP light curve of KIC 4951877, 17 quarters | lightkurve download by `fit/transit_fit.py`; the download is committed as `data/lightcurve/` | 14, `fit/` |
| Transit fit | MCMC posterior of the radius ratio, written by `fit/transit_fit.py` | committed, `data/transit_fit/` | 07 |
| MAST | Kepler target pixel files of KIC 4951877 | lightkurve | `fit/diff_images.py` |
| arXiv | source tables of Armstrong et al. 2021 (2008.10516v1) and Valizadegan et al. 2022 (2111.10009v2) | pinned version, SHA-256 checked | 15 |
| CDS VizieR | Thompson et al. 2018 `J/ApJS/235/38` tables 1-2; Berger et al. 2018 `J/ApJ/866/99`; Valizadegan et al. 2023 `J/AJ/166/28` | full table or identifier query | 15, 16 |
| NASA Exoplanet Archive | `cumulative`, `pscomppars` | TAP | 15 |
| MAST | TESS Input Catalog around the target; Kepler light curves and quarter-1 target pixel file | astroquery and lightkurve | `fpp/triceratops_fpp.py` |
| TRILEGAL | simulated star populations along this line of sight, one for TRICERATOPS and one for vespa (three columns kept) | committed, `fpp/data/` | `fpp/` |
| NASA Exoplanet Archive | `cumulative` KOI table, downloaded 2026-08-22 (four columns kept) | committed, `fpp/data/` | `fpp/own_fpp.py` |

### Committed data

Inputs not served by any query interface are committed in `data/`.

**`data/koiapp/` — DR25 positional probabilities (13 CSV files).** The DR25
`koiapp` table is not exposed by the archive's TAP service or its
`nph-nstedAPI` endpoint (both return `"koiapp" is not a valid table`); it is
available through the interactive viewer. MAST's bulk file
`https://archive.stsci.edu/pub/kepler/catalogs/kepler_koiapp.txt.gz` is an
earlier release with different values (for KOI-501.01 a 655 ppm depth rather
than 575 ppm); step 04 downloads it and reports it separately, but it is not
the DR25 table. To reproduce an export:

1. Open <https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=koiapp>.
2. Choose **Select Columns → All Columns**. The default column set omits
   `pp_1hi_kepmag`, `pp_2hi_kepmag`, `pp_1hi_mod_depth`, `pp_2hi_mod_depth` and
   the provenance flags, which are the inputs examined in §2.3.
3. Filter `kepoi_name_display` on the object, e.g. `K00501.01`.
4. Choose **Download Table → CSV**. The file header records the constraint
   and download time and is left unedited.

`koiapp_2026.09.07_07.39.04.csv` is KOI-501.01. The other twelve are K00696.01,
K01552.01, K03156.01, K03469.01, K05071.01, K05129.01, K05157.01, K05179.01,
K05440.01, K06557.01, K06969.01 and K08180.01: every object that step 02 finds
with disposition CANDIDATE and an impossible-colour neighbour. They support
the §2.6 statement that no other candidate is affected.

**`data/dv/dvr_quarterly_centroids.csv` — difference-image centroids (16 rows).**
The archive tabulates only the combined offset (`koi_dikco_msky` =
0.675 ± 0.255 arcsec). The per-quarter values appear only in the DV report, a
copy of which is in `docs/`:

1. Download
   `https://exoplanetarchive.ipac.caltech.edu/data/KeplerData/004/004951/004951877/dv/kplr004951877-20160209194854_dvr.pdf`
   (the path given by `koi_datalink_dvr`).
2. Extract the text, e.g. `pdftotext -layout <file>.pdf out.txt`.
3. For each quarter, under *PRF Fit of the Difference Image*, read the `Offset`
   row of *Mean offset from the KIC RA and Dec*: RA and Dec offsets in arcsec
   with their errors.

The report's multi-quarter summary on the same pages (robust weighted mean
+0.3309 ± 0.225 and −0.5879 ± 0.235 arcsec) is the cross-check quoted in §3.

**`data/dv/dvr_quarterly_oot_centroids.csv` — out-of-transit centroids (16 rows).**
From the same pages of the same report: for each quarter, under *PRF Fit of the
Difference Image*, the `Out of Transit Image Centroid` row of *Offset from the
PRF fit to the out of transit image* and the `KIC Reference Centroid` row of
*Offset from the KIC RA and Dec converted to pixels via motion polynomials*,
each as row and column [pixel, 0-based, to 0.01] and RA [hours] and Dec [deg].

**`data/lightcurve/kic4951877_pdcsap.npz` — the Kepler light curve.** The file
`fit/transit_fit.py` writes to `results/.cache/` when it downloads all 17
long-cadence quarters from MAST (lightkurve, `quality_bitmask="none"`), each
quarter's PDCSAP flux and error divided by its median. Committed so that step 14
and `fit/bayes_factor.py` run without lightkurve. Arrays `time` [BKJD], `flux`
and `flux_err` [-], `quality` [Kepler bit flags], `n_quarters`; 64,797 cadences.

**`data/diffimages/` — the difference images (244 files).** `kic4951877.npz`
holds, per quarter, the out-of-transit image (`direct_*`), the difference image
(`d_plain_*`) and its error (`err_*`) [e⁻ s⁻¹ per pixel], and metadata (`meta_*`:
quarter, module, output, stamp origin, WCS, transit count), as built by
`fit/diff_images.py`, which rebuilds it identically from MAST. The other 243
files are the controls of `fit/prf_centroids.py`, with the difference image,
error and metadata only; `controls.csv` lists them with their DR25
difference-image offsets [arcsec].

**`data/dv/dvr_multiquarter_mean.csv` — the report's robust weighted mean.** The
multi-quarter offsets from the KIC position and from the out-of-transit
centroid, with errors [arcsec], from the same DV report pages; step 06 uses the
errors of the first.

**`fpp/data/` — inputs of the false-positive programs.** `4951877_TRILEGAL.csv`
is the TRILEGAL simulation the paper's TRICERATOPS runs used, as TRICERATOPS
saved it. `trilegal_starfield.npz` keeps `Kepler_mag` [mag], `Mact` [M☉] and
`logg` [cgs] of the 772,038 stars of the TRILEGAL simulation vespa made for this
field (`starfield.h5`, 142 MB), in their original order. `cumulative_koi.csv`
keeps `kepoi_name`, `koi_disposition`, `koi_period` [d] and `koi_prad` [R⊕] of the
archive's cumulative table as downloaded on 2026-08-22, in its original row
order; `fpp/own_fpp.py` samples planet radii from it, so the order matters.

**`fpp/vespa/` — the vespa run.** Inputs (`fpp.ini`, `star.ini`, `lc.csv`,
`UKIRT_J.cc`) and the outputs of the paper's run (`starfit.log`, `calcfpp.log`,
`results.txt`, `results_bootstrap.txt`, `FPPsummary.png`); `run_vespa.sh`
documents each input and the environment.

**`results/paper_runs/` — the unseeded runs quoted in the paper.**
`bayes_factor.json` is the nested-sampling run behind §6 and Table 4, and
`triceratops_runs.csv` the twenty TRICERATOPS runs behind §7 (their local
contrast-curve path column removed). Both were made by the analysis pipeline's
scripts, of which `fit/bayes_factor.py` and `fpp/triceratops_fpp.py` are ports;
the ports' own reruns are written beside the step outputs in `results/`.

SHA-256 checksums:

```
62d2408d0a5d370f37b13ffb9ec11ee34e9a0de7f6a1f36a6064ebac9a27d864  data/lightcurve/kic4951877_pdcsap.npz
270969f9a72f67d0b7bb3a8738434dafe2768666d8e2c6ad19d69acbdf88ce9e  fpp/data/4951877_TRILEGAL.csv
4e4c7f6c8ef4fd56aca5a9eee23966cfe252c6d91f7accefd75ef0d24b03e732  fpp/data/cumulative_koi.csv
70bb8da9c67c65fc3e413521a919098b1ea802d9dc3391229cf840bdb062c396  fpp/data/trilegal_starfield.npz
```

## Files

### Inputs (`data/`)

| File | Format | Columns used [unit] | Paper |
|---|---|---|---|
| `koiapp/koiapp_*.csv` (13) | CSV, archive export; `#` lines are the archive's own header | `kepid`, `kepoi_name` [-]; `pp_koi_depth` [ppm]; `pp_host_rel_prob`, `pp_host_prob_score`, `pp_1hi_rel_prob`, `pp_2hi_rel_prob` [-]; `pp_1hi_kepmag`, `pp_2hi_kepmag` [mag]; `pp_1hi_mod_depth`, `pp_2hi_mod_depth` [ppm]; `pp_1hi_ra`, `pp_1hi_dec`, `pp_2hi_ra`, `pp_2hi_dec` [deg]; `*_prob_prov`, `*_starid` [-] | §2.3, §2.4, §2.6 |
| `dv/dvr_quarterly_centroids.csv` | CSV | `quarter` [-]; `d_ra`, `d_dec` [arcsec], offset of the difference-image source from the KIC position as tabulated in the DV report; `e_ra`, `e_dec` [arcsec], the report's 1σ errors | §3 |
| `dv/dvr_quarterly_oot_centroids.csv` | CSV | `quarter` [-]; `oot_row`, `oot_col`, `kic_row`, `kic_col` [pixel]; `oot_ra_h`, `kic_ra_h` [hours]; `oot_dec`, `kic_dec` [deg]: the report's out-of-transit PRF centroid and KIC reference position | §2.2 |

### Transit fit (`data/transit_fit/transit_fit.json`)

The MCMC posterior step 07 takes the radius ratio from, written by
`fit/transit_fit.py`. That script is not part of `run_all.py` because it needs
heavier packages and runs for about an hour:

```bash
pip install -r fit/requirements.txt
python fit/transit_fit.py            # all three runs
```

It downloads the light curve from MAST itself. Model: batman 2.5.3
quadratic limb darkening, 29.4-minute exposures supersampled 15 times, circular
orbit, period fixed at the DR25 value; each transit window of ±4 durations
divided by a straight line through its out-of-transit points. Free: Rp/R★,
impact parameter, stellar density, epoch shift, both limb-darkening
coefficients (Kipping q1, q2), and a flux-error scale. Sampled with emcee 3.1.6,
40 walkers, 3000 burn-in and 6000 steps, seed 20260913, on the DR25 PDCSAP light
curve (quality flag 0; 50 windows containing transit points). Three runs:
`adopted` (density prior from step 07, the one used), `free_rho` (density
free) and `fixed_ld` (limb darkening held at `config.TRANSIT_FIT_FIXED_LD`). The
ephemeris, duration, density and stellar radius are read from steps 00 and 07
and recorded under `measured_inputs`.

### Data behind Figure 1 (`data/figure_diffimage/`)

| File | Format | Contents [unit] |
|---|---|---|
| `diffimage_stamps.fits` | FITS, one image pair per quarter | `DIRECT_Qnn`: out-of-transit mean; `DIFF_Qnn`: out-of-transit minus in-transit mean [e⁻ s⁻¹ per pixel]. Each header carries the stamp's WCS, `QUARTER`, `MODULE`, `OUTPUT`, `COLUMN0`, `ROW0` [pixel] and `NTRANSIT`. Built from the archived Kepler target pixel files |
| `diffimage_stacks.fits` | FITS | `ALIGNED_DIRECT`, `ALIGNED_DIFF`: all 16 quarters shifted to put the target at pixel (`TARGETX`, `TARGETY`) and averaged (Figure 1a); `SEASON_M6_DIRECT`, `SEASON_M6_DIFF`: the four module-6 quarters averaged without resampling (Figure 1b) |
| `gaia_neighbours.csv` | CSV, `#` header gives units | Gaia DR3 sources within 40 arcsec: `ra`, `dec` [deg]; `phot_g_mean_mag` [mag]; offsets from the target [arcsec]; pixel position on the Q0 stamp's WCS |

The stamps share one sky orientation, but the target falls on a different pixel
on each of the four CCD modules, up to one pixel apart, which is why Figure 1a
aligns the quarters before averaging.

Column definitions for `koiapp` are in the archive's
[positional probabilities documentation](https://exoplanetarchive.ipac.caltech.edu/docs/API_koiapp_columns.html).

### Outputs (`results/`)

| File | Format | Contents [unit] | Paper |
|---|---|---|---|
| `01_target_photometry.json` | JSON | per star: KIC `kic_kepmag`, `kic_r`, `kic_J`, Gaia `gaia_G`, 2MASS `twomass_J`, Pan-STARRS `panstarrs_r` [mag]; `sep_from_host_as`, `*_match_as` [arcsec]; `flux_shares` [fraction]; host Gaia astrometry (`parallax` [mas], `ruwe` [-]) | §2.2, Table 1, §4 |
| `02_colour_pathology.json` | JSON | counts [-]; `colour_histogram` edges [mag] and counts [-]; one record per impossible-colour neighbour (see the CSV below) | §2.6 |
| `03_kic_gaia_comparison.json` | JSON | per neighbour star: `kic_r`, `kic_J`, `gaia_G`, `G_minus_r`, `G_minus_J` [mag]; counts [-] | §2.6 |
| `04_positional_probability.json` | JSON | validated `koiapp` records (units as in the input table); `corrected_depths`: `magnitude_error` [mag], `flux_factor` [-], depths [ppm] and [%]; `earlier_release` (MAST) in the same form; `candidate_audit` [KOI names] | §2.3, §2.4, §2.6 |
| `05_configuration.json` | JSON | Table 2 rows (see the CSV below); `grid`: `sep_max` [arcsec], `share_min` [fraction], counts [-]; `all_pairs` | §2.6, Table 2 |
| `06_difference_images.json` | JSON | `mean_offset_as`, `scatter_as`, `error_on_mean_as`, `error_formal_as` [arcsec, RA and Dec]; per neighbour `sigma_*` [σ] under each error model; `per_quarter` distances [arcsec] and significance [σ] | §3 |
| `09_radial_velocities.json` | JSON | `systemic_kms` [km s⁻¹]; `scatter_ms`, `error_range_ms`, `semi_amplitude_ms`, `semi_amplitude_err_ms`, `limit_3sigma_ms` [m s⁻¹]; `baseline_d` [d]; `largest_phase_gap` [orbital phase]; `companion_mass_limit_mjup` [M_Jup]; `drift_limit_3sigma_ms_per_yr` [m s⁻¹ yr⁻¹]; `bound_pair_excluded_inside_au` [AU]; `eclipsing_binary_k_over_limit` [-] | §5.1 |
| `00_dr25_target.json` | JSON | the DR25 TCE and KOI rows of the target: period [d], epoch [BKJD], duration [h], depth [ppm], radius ratio, MES, weak secondary [ppm], centroid offsets [arcsec], Kepler magnitude [mag] | read by every program |
| `07_host_star.json` | JSON | catalogue inputs; `radius_rsun`, `mass_msun`, `density_solar`, `luminosity_lsun`, `teff_k` as [median, 1σ]; `planet_radius_earth` on three stellar radii [R⊕]; `orbital_distance_au` [AU], `equilibrium_temperature_k` [K], `insolation_earth` [S⊕] | §4, Table 3 |
| `10_ephemeris_match.json` | JSON | thresholds; counts [-]; statistical matches with `dP_prime`, `dT_prime` [-] and their physical-path test (`separation_as` [arcsec], shared CCD modules) | §6 |
| `12_contrast_curves.json` | JSON | per band: `pixel_scale_as`, `fwhm_as`, `inner_working_angle_as`, `separation_as` [arcsec]; `delta_mag_5sigma` [mag]; ExoFOP file id | §5.2 |
| `pixel_scene.json` | JSON | written by `fit/pixel_scene.py`: per quarter and for uniform and photon-noise weights, the scene-fit flux ratios of KIC 4951867 and KIC 4951861 to the target [-], common shift [pixel], background [e⁻ s⁻¹], worst residual over peak [-]; the report's out-of-transit centroid minus the KIC position and the single-source positions fitted to the model scenes and to the image [arcsec; pixel]; summaries, and the flux ratios the KIC and Gaia magnitudes imply | §2.2 |
| `contrast_curves.csv` | CSV, `#` header gives units | the four curves: band, separation [arcsec], 5σ limit [mag] | §5.2 |
| `11_dr25_vetting.json` | JSON | DR25 TCE statistics: depths [ppm], `boot_fap` [-], `tce_maxmesd` [d], rolling-band counts [-], ghost core/halo statistics [-]; KOI score and flags; `flux_weighted_centroid` offsets [arcsec] and significance [σ] with and without the systematic | §3.1, §6, Table 4 |
| `15_literature.json` | JSON | published scores of Armstrong et al. 2021, ExoMiner and ExoMiner V1.2 [-]; the implied prior [-]; the population counts with the query date | §1, §4 |
| `16_instrumental_share.json` | JSON | synthetic detections and those passed in the cell [-]; real TCEs and candidates in the cell [-]; the 95% upper limit [events] and the share bounds [fraction] | §6, Table 4 |
| `99_check.json` | JSON | stale inputs, vespa input mismatches, typed numbers and paper numbers the outputs no longer give (all empty on a clean run) | all |
| `prf_centroids.json`, `prf_centroids_controls.csv` | JSON, CSV | written by `fit/prf_centroids.py`: the target's east and north offsets, errors and total [arcsec], panels [-]; the controls' statistics and `systematic_per_axis_as` [arcsec]; per control the same fit and DR25's offset [arcsec] | §3.1 |
| `lightcurve_checks.json` | JSON | written by `fit/lightcurve_checks.py`: box depths [ppm]; wrong-epoch statistics; single-transit depths [ppm]; odd/even and halves [ppm, σ]; M0V minimum eccentricity [-]; transit times [BKJD, d] with scatter and precision [min]; CROWDSAP per quarter [-] | §2.6, §6, Table 4 |
| `fpp_published_runs.json` | JSON | written by `fpp/summarise_runs.py`: vespa FPP, bootstrap mean and standard deviation [-], population [-], isochrone radius [R☉]; the published TRICERATOPS runs' mean, median, maximum [-], zeros and passes [-] | §7 |
| `fpp_triceratops_lc.csv` | CSV | written by `fpp/triceratops_fpp.py`: the folded light curve [d from mid-transit, relative flux, error], vespa's `lc.csv` | §7 |
| `impossible_colour_neighbours.csv` | CSV, `#` header gives unit and meaning of every column | the 135 neighbours: KOI, disposition, KIC ID, position [deg], separation [arcsec], KIC and Gaia magnitudes [mag] | §2.6 |
| `configuration.csv` | CSV, `#` header as above | the nine objects of Table 2: separation [arcsec], light share [fraction], DR25 centroid offset [arcsec] and significance [σ], score, flags, MES | Table 2 |
| `13_gaia_positional_probability.json` | JSON | inputs: `dr25_offset_as`, `prf_fit` [arcsec] and `systematic_per_axis_as` from `prf_centroids.json`; `section_3_1`: combined `east_as`, `north_as` [value, error], `offset_as`, `offset_err_as` [arcsec], `offset_sigma` [σ]; `equation_1` and `total_error_weighting`: centroid [arcsec], `p_on_target`, `p_all_others` [-], `nearest_alternative` (`sep_as` [arcsec], `G` [mag], `sigma_from_centroid` [σ]); per source `G` [mag], offsets `dx_as`, `dy_as`, `sep_as` [arcsec], `flux_fraction` [-], `capable`, distance [σ] and probability under both weightings | §3.1, §7 |
| `14_blend_bounds.json` | JSON | `best_fit`: epoch offset, ingress [min], depth [ppm], duration [h]; `bootstrap`: 16th, 50th and 84th percentiles of depth [ppm], duration [h], ingress [min], T/τ [-], and the unscaled duration and ingress [d]; `max_true_depth_p95_ppm` [ppm]; `dmag_shape_bound`, `dmag_brightness_bound` [mag]; per Gaia neighbour `G`, `dmag` [mag], `sep_as` [arcsec], `centroid_sigma` [σ] and the three tests | §6, Table 4 |
| `bayes_factor.json` | JSON | written by `fit/bayes_factor.py`: `wide`, `narrow`, `control` runs with their priors, `logz_*` and `B` [natural log], transit posterior means and standard deviations (`t0` [d], `log_rp` [log10], `aor` [-], `b` [-], GP `log_sigma` [log10 flux], `log_rho` [log10 d], `log_jit` [log10 flux]) | §6, Table 4 |
| `fpp_triceratops.json`, `fpp_triceratops_runs.csv`, `ukirt_J_triceratops.txt` | JSON, CSV, text | written by `fpp/triceratops_fpp.py`: inputs; per run `FPP`, `NFPP` [-], seed, stars, duplicates removed, aperture pixels; mean, standard error, median, maximum, exact zeros, runs clearing both bars; the contrast curve as TRICERATOPS reads it, separation [arcsec] and contrast [mag] | §7 |
| `fpp_own.json` | JSON | written by `fpp/own_fpp.py`: inputs; scenario priors [-]; expected background stars [-]; run B factors [-]; `run_a` and `run_b`: mean, standard error, maximum and the twenty values [-], and the astrophysical part after removing the flat allowances [-] | §7, Appendix A |
| `paper_runs/` | JSON, CSV | the unseeded runs quoted in the paper (see Committed data) | §6, §7 |

`figures/K00501.01_diffimage.png` (Figure 1) is drawn by `fit/diffimage_figure.py`.
`figures/colour_pathology.{pdf,png}` is drawn by step 08 from
`02_colour_pathology.json` and `05_configuration.json`; it is a supporting
figure and does not appear in the paper.

## Validation

Every record is parsed into a [pydantic](https://docs.pydantic.dev) model
(`koi501/models.py`) on arrival: KIC stars, Gaia sources, 2MASS and Pan-STARRS
matches, cross-match pairs, DR25 objects of interest, positional-probability
records, quarterly centroids and APOGEE visits. Magnitudes and coordinates must
be numeric, probabilities must lie in [0, 1], uncertainties must be positive,
and blank cells become `None`. A change in an archive's schema raises an error
rather than altering a result.

All thresholds and constants are defined in `koi501/config.py`. The separation (8 arcsec) and
light-share (25%) limits of Table 2 were set after KOI-501.01's position was
known; step 05 therefore also evaluates the 15-point grid reported in §2.6. No
subsampling is used anywhere: step 02 covers all 8,054 objects of interest.

## Layout

```
run_all.py    runs the steps in order
koi501/       config.py    thresholds and constants
              models.py    pydantic record models
              archives.py  TAP, VizieR, X-Match, MAST and ExoFOP clients with an on-disk cache
              report.py    shared functions and printed tables
              fits_image.py  reader for the single-image FITS follow-up frames
              inputs.py    measured values read from the step outputs, and the check
              claims.py    the paper's numbers against the outputs (claims.csv)
steps/        00_dr25_target.py
              01_target_photometry.py
              02_colour_pathology.py
              03_kic_gaia_comparison.py
              04_positional_probability.py
              05_configuration.py
              06_difference_images.py
              07_host_star.py
              08_figures.py
              09_radial_velocities.py
              10_ephemeris_match.py
              11_dr25_vetting.py
              12_contrast_curves.py
              13_gaia_positional_probability.py
              14_blend_bounds.py
              15_literature.py
              16_instrumental_share.py
              99_check.py
data/         committed datasets (see above)
docs/         the NASA reports values are transcribed or quoted from, with checksums
fit/          diff_images.py     the difference images, from the target pixel files
              prf_centroids.py   the PRF fit to them and the centroid systematic
              diffimage_figure.py  Figure 1
              transit_fit.py     the transit fit of Section 4
              pixel_scene.py     Kepler's own pixels, Section 2.2
              lightcurve_checks.py  the Table 4 light-curve rows and transit times
              bayes_factor.py    the Bayes factor of Section 6
              requirements.txt
fpp/          triceratops_fpp.py the TRICERATOPS runs of Section 7
              own_fpp.py         this work's calculation, Section 7 and Appendix A
              summarise_runs.py  the published vespa and TRICERATOPS runs
              _compat.py         numpy and scipy renames triceratops still imports
              vespa/             the vespa run: inputs, outputs, run_vespa.sh, make_inputs.py
              data/              TRILEGAL fields and the cumulative KOI snapshot
              requirements.txt
claims.csv    the paper's own numbers and the outputs that produce them
results/      output of each step and of the programs above; paper_runs/
figures/      output of step 08 and fit/diffimage_figure.py
```

## Citation and licence

Please cite the paper if you use this code. Data are from the NASA Exoplanet
Archive (operated by Caltech under contract with NASA), the ESA *Gaia* mission,
2MASS, Pan-STARRS1, SDSS-IV APOGEE, and the CDS VizieR and X-Match services.
Code is MIT-licensed; data retain the terms of their original providers.
