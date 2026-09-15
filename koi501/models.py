"""Pydantic models for every external record the analysis consumes.

Archive and catalogue queries return untyped CSV with blank cells for missing
values. Validating on ingest means a schema change upstream fails loudly here
rather than silently producing a wrong number downstream.
"""

from __future__ import annotations

from typing import Annotated, Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator

Mag = Annotated[float, Field(gt=-5, lt=30)]
Arcsec = Annotated[float, Field(ge=0)]
Prob = Annotated[float, Field(ge=0, le=1)]


class Strict(BaseModel):
    model_config = ConfigDict(extra="ignore", frozen=True)

    @field_validator("*", mode="before")
    @classmethod
    def _blank_to_none(cls, v):
        if isinstance(v, str) and v.strip() == "":
            return None
        return v


class KOIRow(Strict):
    """A row of the DR25 Kepler Objects of Interest table."""

    kepoi_name: str
    kepid: int
    ra: Optional[float] = None
    dec: Optional[float] = None
    koi_disposition: Optional[str] = None
    koi_score: Optional[Prob] = None
    koi_kepmag: Optional[Mag] = None
    koi_dikco_msky: Optional[Arcsec] = None
    koi_dikco_msky_err: Optional[float] = None
    koi_max_mult_ev: Optional[float] = None
    koi_period: Optional[float] = Field(default=None, gt=0)
    koi_time0bk: Optional[float] = None
    koi_fpflag_nt: Optional[int] = None
    koi_fpflag_ss: Optional[int] = None
    koi_fpflag_co: Optional[int] = None
    koi_fpflag_ec: Optional[int] = None

    @property
    def flags(self) -> str:
        named = (("nt", self.koi_fpflag_nt), ("ss", self.koi_fpflag_ss),
                 ("co", self.koi_fpflag_co), ("ec", self.koi_fpflag_ec))
        return " ".join(n for n, v in named if v == 1) or "none"

    @property
    def offset_sigma(self) -> Optional[float]:
        if self.koi_dikco_msky and self.koi_dikco_msky_err:
            return self.koi_dikco_msky / self.koi_dikco_msky_err
        return None


class TCERow(Strict):
    """The DR25 threshold-crossing event table's vetting statistics."""

    kepid: int
    tce_plnt_num: int
    tce_max_mult_ev: float
    tce_num_transits: int
    tce_depth: float
    tce_depth_err: float = Field(gt=0)
    boot_fap: float = Field(ge=0, le=1)
    tce_bin_oedp_stat: float
    wst_depth: float
    wst_depth_err: float = Field(gt=0)
    wst_robstat: float
    tce_maxmesd: float
    tce_maxmes: float
    tce_cap_stat: float
    tce_hap_stat: float
    tce_rb_tcount0: int
    tce_rb_tcount1: int
    tce_rb_tcount2: int
    tce_rb_tcount3: int
    tce_rb_tcount4: int
    tce_fwm_stat: float


class KICStar(Strict):
    """A Kepler Input Catalog entry, as served by VizieR V/133/kic."""

    kic: int
    ra: Optional[float] = None
    dec: Optional[float] = None
    kepmag: Optional[Mag] = None
    rmag: Optional[Mag] = None
    jmag: Optional[Mag] = None

    @property
    def r_minus_j(self) -> Optional[float]:
        if self.rmag is not None and self.jmag is not None:
            return self.rmag - self.jmag
        return None


class TwoMassSource(Strict):
    """A 2MASS point source (VizieR II/246)."""

    ra: float
    dec: float
    jmag: Optional[Mag] = None


class TwoMassPhotometry(Strict):
    """2MASS J, H and Ks photometry of one source (VizieR II/246)."""

    jmag: Mag
    e_jmag: float = Field(gt=0)
    hmag: Mag
    e_hmag: float = Field(gt=0)
    kmag: Mag
    e_kmag: float = Field(gt=0)
    qflg: str


class WisePhotometry(Strict):
    """AllWISE W2 photometry of one source (VizieR II/328)."""

    w2mag: Mag
    e_w2mag: float = Field(gt=0)


class GaiaDistance(Strict):
    """Bailer-Jones et al. (2021) geometric distance (VizieR I/352)."""

    rgeo: float = Field(gt=0)
    b_rgeo: float = Field(gt=0)
    B_rgeo: float = Field(gt=0)


class ApogeeStar(Strict):
    """APOGEE DR17 ASPCAP parameters of one star (VizieR III/286/catalog)."""

    teff: float = Field(gt=2000, lt=10000)  # literal: schema sanity range
    e_teff: float = Field(gt=0)
    logg: float
    e_logg: float = Field(gt=0)
    feh: float
    e_feh: float = Field(gt=0)
    snr: float = Field(gt=0)


class CatalogueRadius(Strict):
    """A stellar radius from a published catalogue, with its source."""

    rad: float = Field(gt=0)
    teff: Optional[float] = None


class PanstarrsSource(Strict):
    """A Pan-STARRS DR1 source (VizieR II/349)."""

    ra: float
    dec: float
    rmag: Optional[Mag] = None
    rmag_error: Optional[float] = None
    n_detections: Optional[int] = None


class GaiaSource(Strict):
    """A Gaia DR3 source."""

    source_id: int
    ra: float
    dec: float
    phot_g_mean_mag: Optional[Mag] = None
    parallax: Optional[float] = None
    parallax_error: Optional[float] = None
    ruwe: Optional[float] = None
    astrometric_excess_noise: Optional[float] = None
    ipd_frac_multi_peak: Optional[int] = None
    non_single_star: Optional[int] = None


class KICNeighbour(Strict):
    """One object-to-KIC pair returned by the CDS X-Match service.

    The uploaded columns come back alongside the catalogue's, so this carries
    both the object it belongs to and the neighbour's photometry.
    """

    koi: str
    kic: Optional[int] = None
    sep_as: Optional[Arcsec] = None
    ra: Optional[float] = None
    dec: Optional[float] = None
    kepmag: Optional[Mag] = None
    rmag: Optional[Mag] = None
    jmag: Optional[Mag] = None

    @property
    def testable(self) -> bool:
        """True when the colour test can be applied to this pair."""
        return None not in (self.kic, self.rmag, self.jmag, self.sep_as)

    @property
    def r_minus_j(self) -> Optional[float]:
        if self.rmag is not None and self.jmag is not None:
            return self.rmag - self.jmag
        return None


class GaiaMatch(Strict):
    """One KIC-star-to-Gaia pair returned by the CDS X-Match service."""

    kic: int
    gmag: Optional[Mag] = None
    sep_as: Optional[Arcsec] = None


class FitsImageHeader(Strict):
    """The header keywords a 2-D follow-up image must carry."""

    model_config = ConfigDict(extra="ignore", frozen=True, populate_by_name=True)

    simple: bool = Field(alias="SIMPLE")
    bitpix: int = Field(alias="BITPIX")
    naxis: int = Field(alias="NAXIS", ge=2, le=2)
    naxis1: int = Field(alias="NAXIS1", gt=0)
    naxis2: int = Field(alias="NAXIS2", gt=0)
    bzero: float = Field(default=0.0, alias="BZERO")
    bscale: float = Field(default=1.0, alias="BSCALE")
    cd1_1: float = Field(alias="CD1_1")
    cd1_2: float = Field(alias="CD1_2")
    cd2_1: float = Field(alias="CD2_1")
    cd2_2: float = Field(alias="CD2_2")

    @property
    def pixel_scale_as(self) -> float:
        """Arcsec per pixel from the linear part of the CD matrix."""
        det = self.cd1_1 * self.cd2_2 - self.cd1_2 * self.cd2_1
        return abs(det) ** 0.5 * 3600.0


class EclipsingBinary(Strict):
    """A Kepler eclipsing binary (Kirk et al. 2016, VizieR J/AJ/151/68)."""

    kic: int
    per: float = Field(gt=0)
    bjd0: float


class KeplerFieldPosition(Strict):
    """MAST Kepler field-of-view record: position, magnitude, CCD modules."""

    ra_sexagesimal: str
    dec_sexagesimal: str
    kepmag: Optional[Mag] = None
    module_0: Optional[int] = None
    module_1: Optional[int] = None
    module_2: Optional[int] = None
    module_3: Optional[int] = None

    @property
    def ra(self) -> float:
        h, m, s = (float(v) for v in self.ra_sexagesimal.split())
        return 15.0 * (h + m / 60.0 + s / 3600.0)

    @property
    def dec(self) -> float:
        text = self.dec_sexagesimal.strip()
        sign = -1.0 if text.startswith("-") else 1.0
        d, m, s = (float(v) for v in text.lstrip("+-").split())
        return sign * (d + m / 60.0 + s / 3600.0)

    @property
    def modules(self) -> set[int]:
        return {m for m in (self.module_0, self.module_1, self.module_2,
                            self.module_3) if m}


class PositionalProbability(Strict):
    """A row of a positional-probability table (`koiapp`).

    The DR25 rows are CSVs exported from the NASA Exoplanet Archive's table
    viewer. MAST's bulk file is an earlier release and has no KOI name column.
    """

    kepid: int
    kepoi_name: Optional[str] = None
    pp_koi_depth: Optional[float] = None
    pp_host_rel_prob: Optional[Prob] = None
    pp_host_prob_score: Optional[Prob] = None
    pp_host_prob_prov: Optional[str] = None
    pp_1hi_starid: Optional[str] = None
    pp_1hi_kepmag: Optional[Mag] = None
    pp_1hi_rel_prob: Optional[Prob] = None
    pp_1hi_mod_depth: Optional[float] = None
    pp_1hi_prob_prov: Optional[str] = None
    pp_1hi_ra: Optional[float] = None
    pp_1hi_dec: Optional[float] = None
    pp_2hi_starid: Optional[str] = None
    pp_2hi_kepmag: Optional[Mag] = None
    pp_2hi_rel_prob: Optional[Prob] = None
    pp_2hi_mod_depth: Optional[float] = None
    pp_2hi_prob_prov: Optional[str] = None
    pp_2hi_ra: Optional[float] = None
    pp_2hi_dec: Optional[float] = None

    @property
    def favours_host(self) -> bool:
        return self.pp_host_rel_prob is not None and self.pp_host_rel_prob >= 0.5

    @property
    def usable(self) -> bool:
        return (self.pp_host_prob_score or 0) > 0 and self.pp_host_prob_prov != "FAILED"


class QuarterCentroid(Strict):
    """One quarter's difference-image centroid offset from the DV report."""

    quarter: int
    d_ra: float
    e_ra: float = Field(gt=0)
    d_dec: float
    e_dec: float = Field(gt=0)


class ApogeeVisit(Strict):
    """One APOGEE DR17 visit."""

    jd: float
    vhelio_kms: float
    e_rv_kms: float = Field(gt=0)
    ncomp: Optional[int] = None
