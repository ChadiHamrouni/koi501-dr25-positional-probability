"""Measured inputs, read from the step outputs in results/.

A measurement lives in one place: the results file of the step that made or
fetched it. Programs ask for it here by a dotted name, record what they used in
their own output under "inputs", and check() compares every recorded input with
its source, so an output made from values that have since changed stops the run.

Standard library only, so every program can import it whatever its environment.
"""

from __future__ import annotations

import json
from typing import Any, Callable

from . import config


def _read(name: str) -> dict:
    path = config.RESULTS / f"{name}.json"
    if not path.exists():
        raise SystemExit(f"inputs: results/{name}.json is missing; run the step that writes it")
    return json.loads(path.read_text(encoding="utf8"))


def _transit() -> dict:
    d = _read("00_dr25_target")
    koi, tce = d["koi"], d["tce"]
    return {
        # the detection table (TCE): the fit the false-positive codes use
        "period_tce_d": tce["tce_period"], "epoch_bkjd": tce["tce_time0bk"],
        "duration_h": tce["tce_duration"], "duration_err_h": tce["tce_duration_err"],
        "depth_ppm": tce["tce_depth"], "depth_err_ppm": tce["tce_depth_err"],
        "ror": tce["tce_ror"], "mes": tce["tce_max_mult_ev"],
        "n_transits": tce["tce_num_transits"],
        "weak_secondary_ppm": tce["wst_depth"],
        "weak_secondary_err_ppm": tce["wst_depth_err"],
        "weak_secondary_robstat": tce["wst_robstat"],
        "odd_even_stat": tce["tce_bin_oedp_stat"],
        "dikco_msky_as": tce["tce_dikco_msky"], "dikco_msky_err_as": tce["tce_dikco_msky_err"],
        # the candidate table (KOI): the period to more digits, used for folding
        "period_koi_d": koi["koi_period"], "kepmag": koi["koi_kepmag"],
        "dikco_mra_as": koi["koi_dikco_mra"], "dikco_mra_err_as": koi["koi_dikco_mra_err"],
        "dikco_mdec_as": koi["koi_dikco_mdec"], "dikco_mdec_err_as": koi["koi_dikco_mdec_err"],
    }


def _star() -> dict:
    host = _read("07_host_star")
    phot = _read("01_target_photometry")
    apogee = host["inputs"]["apogee"]
    twomass = host["inputs"]["twomass"]
    return {
        "teff_k": apogee["teff"], "teff_err_k": config.TEFF_WIDENED_K,
        "logg": apogee["logg"], "logg_err": config.LOGG_WIDENED,
        "feh": apogee["feh"], "feh_err": config.FEH_WIDENED,
        "radius_rsun": host["radius_rsun"][0], "radius_err_rsun": host["radius_rsun"][1],
        "mass_msun": host["mass_msun"][0], "mass_err_msun": host["mass_msun"][1],
        "density_solar": host["density_solar"][0], "density_err_solar": host["density_solar"][1],
        "parallax_mas": phot["host_gaia"]["parallax"],
        "parallax_err_mas": phot["host_gaia"]["parallax_error"],
        "gaia_g": phot["host_gaia"]["phot_g_mean_mag"],
        "j_mag": twomass["jmag"], "j_err": twomass["e_jmag"],
        "h_mag": twomass["hmag"], "h_err": twomass["e_hmag"],
        "ks_mag": twomass["kmag"], "ks_err": twomass["e_kmag"],
    }


def _centroid() -> dict:
    s = _read("13_gaia_positional_probability")["section_3_1"]
    return {"offset_as": s["offset_as"], "offset_err_as": s["offset_err_as"],
            "systematic_per_axis_as": _read("prf_centroids")["systematic_per_axis_as"]}


def _shape() -> dict:
    b = _read("14_blend_bounds")["bootstrap"]
    return {"duration_d_median": b["duration_d_median"],
            "ingress_d_p16": b["ingress_d_p16_p50_p84"][0],
            "ingress_d_p50": b["ingress_d_p16_p50_p84"][1],
            "ingress_d_p84": b["ingress_d_p16_p50_p84"][2],
            "depth_ppm_p50": b["depth_ppm_p16_p50_p84"][1]}


def _velocities() -> dict:
    v = _read("09_radial_velocities")
    return {"k_limit_ms": v["limit_3sigma_ms"],
            "bound_pair_min_au": v["bound_pair_excluded_inside_au"]}


def _contrast() -> dict:
    j = _read("12_contrast_curves")["J"]
    return {"j_separation_as": j["separation_as"], "j_delta_mag": j["delta_mag_5sigma"]}


def _instrumental() -> dict:
    s = _read("16_instrumental_share")
    return {"n_passed": s["n_noise_passed"], "n_synthetic": s["n_noise_in_cell"],
            "upper_bound": s["share_upper_bounds"][config.INSTRUMENTAL_QUOTED_BOUND]}


NAMESPACES: dict[str, Callable[[], dict]] = {
    "transit": _transit, "star": _star, "centroid": _centroid, "shape": _shape,
    "velocities": _velocities, "contrast": _contrast, "instrumental": _instrumental,
}
_cache: dict[str, dict] = {}


def get(name: str) -> Any:
    """One measured value by its dotted name, e.g. "star.radius_rsun"."""
    space, key = name.split(".", 1)
    if space not in _cache:
        _cache[space] = NAMESPACES[space]()
    return _cache[space][key]


def record(names: list[str]) -> dict:
    """The values of several inputs, to store under "inputs" in an output."""
    return {n: get(n) for n in names}


def outputs_with_inputs() -> list:
    """Every JSON output that records the measured inputs it was made from."""
    candidates = sorted(config.RESULTS.glob("*.json")) + [
        config.DATA / "transit_fit" / "transit_fit.json"]
    out = []
    for path in candidates:
        try:
            if path.exists() and "measured_inputs" in json.loads(path.read_text(encoding="utf8")):
                out.append(path)
        except ValueError:
            continue
    return out


def check(stop: bool = True) -> list[str]:
    """Compare every recorded input with its current source value.

    An output stores the inputs it was made from under "measured_inputs". If any
    of them no longer equals the value its source step now writes, the output is
    stale: it has to be remade before it can be trusted, and the run stops.
    """
    problems = []
    _cache.clear()
    for path in outputs_with_inputs():
        stored = json.loads(path.read_text(encoding="utf8"))["measured_inputs"]
        for name, value in stored.items():
            try:
                current = get(name)
            except SystemExit as missing:
                problems.append(f"{path.name}: {name}: {missing}")
                continue
            if current != value:
                problems.append(f"{path.name}: {name} was {value}, source now gives {current}")
    if problems and stop:
        raise SystemExit("inputs disagree with their sources:\n  " + "\n  ".join(problems))
    return problems
