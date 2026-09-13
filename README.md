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

This repository recomputes those results from the public archives.

## Running

```bash
python -m venv .venv && . .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python run_all.py            # all steps
python run_all.py 04 06      # selected steps
```

Queries are cached in `results/.cache/`; with the cache present the full run
takes under 30 seconds, and a first run spends a few minutes downloading.
Requests to the archives are retried automatically if a connection drops.

## Requirements

| | |
|---|---|
| Language | Python 3.12.7; exact package versions in `requirements.txt` |
| Tested on | Windows 11 Pro (10.0.26200) |
| Network | HTTPS access to exoplanetarchive.ipac.caltech.edu, gea.esac.esa.int, vizier.cds.unistra.fr, cdsxmatch.u-strasbg.fr, archive.stsci.edu |
| Disk | about 11 MB of cached queries, 0.5 MB of outputs |
| Inputs | the archives above and the committed files in `data/` |
| Outputs | `results/*.json`, `results/*.csv`, `figures/colour_pathology.{pdf,png}` |

`run_all.py` executes the scripts in `steps/` in order. Each one queries its
sources, validates every returned record, performs one part of the analysis,
prints the quantities the paper quotes under the relevant section heading, and
writes its full output to `results/<step>.json`. The script reports values
only; it does not compare them with the manuscript.

## Claims and where they are reproduced

| Paper | Claim | Step | Data |
|---|---|---|---|
| §2.2, Table 1 | KIC gives KIC 4951867 Kp = 14.830 and KIC 4951861 Kp = 15.485; Gaia DR3 gives G = 17.787 and 17.513 | 01 | KIC, Gaia DR3 |
| §2.2 | KIC 4951867 has r − J = −1.49 in the KIC's own columns; the host has +1.17 | 01 | KIC |
| §2.2 | 2MASS J and Pan-STARRS DR1 r agree with Gaia, not with the KIC; the host's KIC magnitude is correct to 0.02 mag | 01 | 2MASS, Pan-STARRS DR1 |
| §2.2 | No Gaia source brighter than G = 16.23 other than the host within 25 arcsec | 01 | Gaia DR3 |
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
| §3 | KIC 4951867 lies 30.7σ and KIC 4951861 38.6σ from the source; 99.8σ, 15.4σ and 7.7σ under the alternative error models; all 16 quarters closer to the target | 06 | `data/dv/`, KIC |
| §5.1 | 19 APOGEE visits; K = 12 ± 63 m s⁻¹; 3σ limit 199 m s⁻¹, i.e. M < 3.3 M_Jup | 07 | APOGEE DR17 |

Step 08 draws `figures/colour_pathology.pdf`, a supporting figure of the step 02
and 05 results that does not appear in the paper. The left panel is the r − J
distribution of all tested neighbours with the −0.5 cut; the right panel places
each impossible-colour neighbour by separation and light share, with the
configuration of Table 2 shaded.

Not reproduced here: the stellar radius and planet radius (§4), the contrast
curves (§5.2), the vetting diagnostics (§6) and the false-positive probability
(§7).

## Sources

| Source | Table or catalogue | Access | Steps |
|---|---|---|---|
| NASA Exoplanet Archive | `q1_q17_dr25_koi` (DR25 objects of interest) | TAP | 02, 04, 05 |
| NASA Exoplanet Archive | `koiapp` (DR25 positional probabilities) | manual export, committed | 04 |
| MAST | `kepler_koiapp.txt.gz` (earlier positional-probability release) | bulk file | 04 |
| Kepler DV report | per-quarter difference-image centroids, KIC 4951877 | manual extraction, committed | 06 |
| ESA Gaia archive | Gaia DR3 `gaia_source` | TAP | 01 |
| CDS VizieR | KIC `V/133/kic`, 2MASS `II/246`, Pan-STARRS DR1 `II/349`, APOGEE DR17 `III/286/allvis` | cone search | 01, 06, 07 |
| CDS X-Match | DR25 objects of interest × KIC; flagged neighbours × Gaia DR3 `I/355/gaiadr3` | bulk cross-match | 02, 03 |

### Committed data

Two datasets are not served by any query interface and are committed in `data/`.

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

## Files

### Inputs (`data/`)

| File | Format | Columns used [unit] | Paper |
|---|---|---|---|
| `koiapp/koiapp_*.csv` (13) | CSV, archive export; `#` lines are the archive's own header | `kepid`, `kepoi_name` [-]; `pp_koi_depth` [ppm]; `pp_host_rel_prob`, `pp_host_prob_score`, `pp_1hi_rel_prob`, `pp_2hi_rel_prob` [-]; `pp_1hi_kepmag`, `pp_2hi_kepmag` [mag]; `pp_1hi_mod_depth`, `pp_2hi_mod_depth` [ppm]; `pp_1hi_ra`, `pp_1hi_dec`, `pp_2hi_ra`, `pp_2hi_dec` [deg]; `*_prob_prov`, `*_starid` [-] | §2.3, §2.4, §2.6 |
| `dv/dvr_quarterly_centroids.csv` | CSV | `quarter` [-]; `d_ra`, `d_dec` [arcsec], offset of the difference-image source from the KIC position as tabulated in the DV report; `e_ra`, `e_dec` [arcsec], the report's 1σ errors | §3 |

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
| `07_radial_velocities.json` | JSON | `systemic_kms` [km s⁻¹]; `scatter_ms`, `error_range_ms`, `semi_amplitude_ms`, `semi_amplitude_err_ms`, `limit_3sigma_ms` [m s⁻¹]; `baseline_d` [d]; `largest_phase_gap` [orbital phase]; `companion_mass_limit_mjup` [M_Jup] | §5.1 |
| `impossible_colour_neighbours.csv` | CSV, `#` header gives unit and meaning of every column | the 135 neighbours: KOI, disposition, KIC ID, position [deg], separation [arcsec], KIC and Gaia magnitudes [mag] | §2.6 |
| `configuration.csv` | CSV, `#` header as above | the nine objects of Table 2: separation [arcsec], light share [fraction], DR25 centroid offset [arcsec] and significance [σ], score, flags, MES | Table 2 |

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

All thresholds are defined in `koi501/config.py`. The separation (8 arcsec) and
light-share (25%) limits of Table 2 were set after KOI-501.01's position was
known; step 05 therefore also evaluates the 15-point grid reported in §2.6. No
subsampling is used anywhere: step 02 covers all 8,054 objects of interest.

## Layout

```
run_all.py    runs the steps in order
koi501/       config.py    thresholds and constants
              models.py    pydantic record models
              archives.py  TAP, VizieR, X-Match and MAST clients with an on-disk cache
              report.py    shared functions and printed tables
steps/        01_target_photometry.py
              02_colour_pathology.py
              03_kic_gaia_comparison.py
              04_positional_probability.py
              05_configuration.py
              06_difference_images.py
              07_radial_velocities.py
              08_figures.py
data/         committed datasets (see above)
docs/         the NASA reports values are transcribed or quoted from, with checksums
results/      output of each step
figures/      output of step 08
```

## Citation and licence

Please cite the paper if you use this code. Data are from the NASA Exoplanet
Archive (operated by Caltech under contract with NASA), the ESA *Gaia* mission,
2MASS, Pan-STARRS1, SDSS-IV APOGEE, and the CDS VizieR and X-Match services.
Code is MIT-licensed; data retain the terms of their original providers.
