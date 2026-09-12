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


class PositionalProbability(Strict):
    """A row of the archive's DR25 positional-probability table (`koiapp`).

    This table has no API. Rows come from CSVs exported by hand from the
    archive's interactive table view; see data/README.md.
    """

    kepid: int
    kepoi_name: str
    pp_host_rel_prob: Optional[Prob] = None
    pp_host_prob_score: Optional[Prob] = None
    pp_host_prob_prov: Optional[str] = None
    pp_1hi_starid: Optional[str] = None
    pp_1hi_kepmag: Optional[Mag] = None
    pp_1hi_rel_prob: Optional[Prob] = None
    pp_1hi_mod_depth: Optional[float] = None
    pp_1hi_prob_prov: Optional[str] = None
    pp_2hi_starid: Optional[str] = None
    pp_2hi_kepmag: Optional[Mag] = None
    pp_2hi_rel_prob: Optional[Prob] = None
    pp_2hi_mod_depth: Optional[float] = None
    pp_2hi_prob_prov: Optional[str] = None

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
