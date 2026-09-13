# Reports

NASA Kepler mission documents from which the paper takes values or quotations.
They are PDFs with no query interface, so copies are kept here exactly as
downloaded. Each can be checked against its source with the SHA-256 checksum.

| File | Document | Used for | Source |
|---|---|---|---|
| `kplr004951877-20160209194854_dvr.pdf` | Kepler Data Validation report, KIC 4951877, Q1–Q17 DR25 (SOC 9.3, 2016) | the 16 per-quarter difference-image offsets in `data/dv/dvr_quarterly_centroids.csv` (§3); the robust weighted mean and confusion radius quoted in §3 | [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/data/KeplerData/004/004951/004951877/dv/kplr004951877-20160209194854_dvr.pdf), path given by `koi_datalink_dvr` in the DR25 KOI table |
| `kplr004951877_q1_q16_tcert.pdf` | Q1–Q16 vetting (TCERT) report, KIC 4951877 (SOC 9.1, 2013) | the flux-weighted centroid offset of 0.14 ± 0.33 arcsec (§3) and the model-shift secondary test in Table 4 (pages 5–6) | [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/data/KeplerData/004/004951/004951877/tcert/kplr004951877_q1_q16_tcert.pdf) |
| `KSCI-19108-001.pdf` | Bryson & Morton (2017), *Planet Reliability Metrics: Astrophysical Positional Probabilities for Data Release 25* | how the DR25 positional probability is computed, its input catalogue, the quality score, and the 3 × 10⁶ ppm rejection limit (§2.1, §2.4, §2.5) | [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/docs/KSCI-19108-001.pdf) |
| `KSCI-19092-001.pdf` | Bryson & Morton (2015), *Planet Reliability Metrics: Astrophysical Positional Probabilities* | documentation of the DR24 release of the table distributed by MAST (§2.3) | [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/docs/target_probability_KSCI.pdf) |

## Checksums (SHA-256)

```
2d6c0b9f93727c9c7b809466f8a683dead2f04b7935f94a8f861d51b01e992a8  kplr004951877-20160209194854_dvr.pdf
26761f0b8820c723d15a2fd7b7bd355282da2b4761c20b4461c318f9e38f2451  kplr004951877_q1_q16_tcert.pdf
fdaa0732a33ac62d72bdd3eabf0946b076176d3d1e93104f55d58266517dfc0a  KSCI-19108-001.pdf
3e86c97b122f7793b3f1aeef1e2f79f088837bfd29b46df82e7fc8083c6286c1  KSCI-19092-001.pdf
```

Verify with `sha256sum -c` on Linux or macOS, or `Get-FileHash -Algorithm SHA256 <file>` in
PowerShell.
