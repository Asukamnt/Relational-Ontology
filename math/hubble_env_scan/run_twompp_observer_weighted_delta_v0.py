"""
2M++ observer-weighted delta_g* (v0)
===================================

Goal
----
Define an *observer-centered* effective overdensity proxy

  delta_eff = (∫ delta_g*(x) w(x) d^3x) / (∫ w(x) d^3x),

where the weight function w(r, n) approximates the distance-ladder / Hubble-flow
estimator geometry (radial window + sky coverage).

This is meant to operationalize the "observer environment" quantity in a way that
is closer to how H0 is actually inferred (weighted by a tracer / survey geometry),
and to compare it to the delta required to explain the full Planck→SH0ES shift.

Inputs
------
- Pantheon+SH0ES public distance table (for an empirical radial window and a
  crude latitude-only sky selection proxy).
- 2M++ reconstructed density field (Carrick+2015): delta_g* on a 257^3 cube.
- Optional: read beta / H0 values from results_v1.json (repo output).

Outputs
-------
- twompp_observer_weighted_delta_results_v0.json
- twompp_observer_weighted_delta_plot_v0.png

Run (repo root)
---------------
python math/hubble_env_scan/run_twompp_observer_weighted_delta_v0.py

Notes
-----
- delta_g* is a luminosity-weighted galaxy overdensity proxy (not delta_m).
- The latitude model here is axisymmetric (depends only on |b|), and is a
  conservative proxy for a complex survey footprint.
- Keep console output ASCII-friendly (Windows terminals).
"""

from __future__ import annotations

import argparse
import json
import math
import urllib.request
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d

import twompp_density_lookup_v0 as twompp


PANTHEONPLUS_URL = (
    "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
    "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat"
)

C_KMS = 299792.458

DEFAULT_RESULTS_V1 = Path(__file__).parent / "results_v1.json"

OUTPUT_JSON = Path(__file__).parent / "twompp_observer_weighted_delta_results_v0.json"
OUTPUT_PLOT = Path(__file__).parent / "twompp_observer_weighted_delta_plot_v0.png"


def _cache_dir() -> Path:
    d = Path(__file__).parent / "_cache"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _download_if_missing(url: str, path: Path) -> None:
    if path.exists():
        return
    print(f"Downloading Pantheon+ distance file -> {path}")
    print(f"Source: {url}")
    urllib.request.urlretrieve(url, path)  # noqa: S310 (intentional external download)


def _json_sanitize(obj: Any):
    if isinstance(obj, (float, np.floating)):
        v = float(obj)
        return v if math.isfinite(v) else None
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, dict):
        return {k: _json_sanitize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_sanitize(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _json_sanitize(obj.tolist())
    return obj


@dataclass
class PpRow:
    cid: str
    idsurvey: int
    zhd: float
    zcmb: float
    ra_deg: float
    dec_deg: float
    is_calibrator: int
    used_in_sh0es_hf: int


def read_pantheonplus_dat(path: Path) -> List[PpRow]:
    with path.open("r", encoding="utf-8") as f:
        header = f.readline().split()
        idx = {name: i for i, name in enumerate(header)}

        required = ["CID", "IDSURVEY", "zHD", "zCMB", "RA", "DEC", "IS_CALIBRATOR", "USED_IN_SH0ES_HF"]
        missing = [c for c in required if c not in idx]
        if missing:
            raise ValueError(f"Missing columns in Pantheon+ file: {missing}")

        rows: List[PpRow] = []
        for ln in f:
            if not ln.strip():
                continue
            parts = ln.split()
            try:
                cid = parts[idx["CID"]]
                idsurvey = int(float(parts[idx["IDSURVEY"]]))
                zhd = float(parts[idx["zHD"]])
                zcmb = float(parts[idx["zCMB"]])
                ra = float(parts[idx["RA"]])
                dec = float(parts[idx["DEC"]])
                is_cal = int(float(parts[idx["IS_CALIBRATOR"]]))
                used_hf = int(float(parts[idx["USED_IN_SH0ES_HF"]]))
            except Exception:
                continue
            if not (math.isfinite(zhd) and math.isfinite(zcmb) and zhd > 0 and zcmb > 0):
                continue
            if not (math.isfinite(ra) and math.isfinite(dec)):
                continue
            rows.append(
                PpRow(
                    cid=str(cid),
                    idsurvey=int(idsurvey),
                    zhd=float(zhd),
                    zcmb=float(zcmb),
                    ra_deg=float(ra),
                    dec_deg=float(dec),
                    is_calibrator=int(is_cal),
                    used_in_sh0es_hf=int(used_hf),
                )
            )
    return rows


def dedup_by_cid(rows: List[PpRow]) -> List[PpRow]:
    # Keep first occurrence (Pantheon+ has per-CID duplicates rarely; for selection geometry it doesn't matter much).
    seen = set()
    out: List[PpRow] = []
    for r in rows:
        if r.cid in seen:
            continue
        seen.add(r.cid)
        out.append(r)
    return out


def z_to_r_mpc_h(z: np.ndarray) -> np.ndarray:
    # Low-z approximation consistent with other scripts: r [Mpc/h] = (c z) / 100.
    return (C_KMS * np.asarray(z, dtype=float)) / 100.0


def build_radial_weight_from_sample(
    r_mpc_h: np.ndarray,
    *,
    rmin: float,
    rmax: float,
    bin_width: float,
    smooth_sigma_bins: float,
) -> Dict[str, Any]:
    r = np.asarray(r_mpc_h, dtype=float)
    r = r[np.isfinite(r)]
    r = r[(r >= float(rmin)) & (r <= float(rmax))]
    if r.size < 50:
        raise ValueError("Not enough Pantheon+ rows after r cuts to build radial weights.")

    # Histogram counts N_i over r; then use w(r) ~ N_i / r_mid^2 so that
    # volume integral weights shells approximately by object counts N_i.
    edges = np.arange(float(rmin), float(rmax) + float(bin_width), float(bin_width), dtype=float)
    if edges.size < 6:
        raise ValueError("Too few radial bins; increase r-range or reduce bin_width.")

    counts, _ = np.histogram(r, bins=edges)
    mids = 0.5 * (edges[:-1] + edges[1:])

    # Avoid divide-by-zero; rmin should keep us safely away from 0.
    w = counts.astype(float) / np.maximum(mids * mids, 1.0)

    if float(smooth_sigma_bins) > 0:
        w = gaussian_filter1d(w, sigma=float(smooth_sigma_bins), mode="nearest")

    w = np.clip(w, 0.0, None)
    if not np.any(w > 0):
        raise ValueError("Radial weight is all zeros after smoothing/clipping.")

    # Normalize to a convenient scale (mean=1 over the supported range, weighted by bin widths).
    mean_w = float(np.average(w, weights=np.ones_like(w)))
    if mean_w > 0:
        w = w / mean_w

    return {
        "r_edges_mpc_h": edges,
        "r_mids_mpc_h": mids,
        "counts": counts,
        "w_r_mids": w,
        "rmin_mpc_h": float(rmin),
        "rmax_mpc_h": float(rmax),
        "bin_width_mpc_h": float(bin_width),
        "smooth_sigma_bins": float(smooth_sigma_bins),
        "definition": "w_r(mid) ~ N(mid)/mid^2, smoothed, then normalized to mean=1; used as w(r) in delta_eff integral",
    }


def build_latitude_weight_from_sample(
    b_deg: np.ndarray,
    *,
    n_bins: int,
    smooth_sigma_bins: float,
    clip_min: float,
    clip_max: float,
) -> Dict[str, Any]:
    b = np.asarray(b_deg, dtype=float)
    b = b[np.isfinite(b)]
    absb = np.abs(b)
    absb = absb[(absb >= 0.0) & (absb <= 90.0)]
    if absb.size < 50:
        raise ValueError("Not enough rows to build latitude weights.")

    edges = np.linspace(0.0, 90.0, int(n_bins) + 1, dtype=float)
    mids = 0.5 * (edges[:-1] + edges[1:])
    counts, _ = np.histogram(absb, bins=edges)

    # Isotropic expectation for |b| in [0,90]: pdf(b) = cos(b) with b in radians (normalized).
    # Bin expected fraction = sin(b_hi)-sin(b_lo).
    lo = np.deg2rad(edges[:-1])
    hi = np.deg2rad(edges[1:])
    expected_frac = np.sin(hi) - np.sin(lo)  # sums to 1
    expected = expected_frac * float(np.sum(counts))

    eps = 1e-12
    ratio = (counts.astype(float) + eps) / (expected + eps)

    if float(smooth_sigma_bins) > 0:
        ratio = gaussian_filter1d(ratio, sigma=float(smooth_sigma_bins), mode="nearest")

    ratio = np.clip(ratio, float(clip_min), float(clip_max))

    # Normalize so that the isotropic-weighted mean is 1.
    norm = float(np.sum(ratio * expected_frac))
    if norm > 0:
        ratio = ratio / norm

    return {
        "absb_edges_deg": edges,
        "absb_mids_deg": mids,
        "counts": counts,
        "expected_counts_isotropic": expected,
        "w_lat_mids": ratio,
        "n_bins": int(n_bins),
        "smooth_sigma_bins": float(smooth_sigma_bins),
        "clip_min": float(clip_min),
        "clip_max": float(clip_max),
        "definition": "axisymmetric latitude weight w_lat(|b|) from sample / isotropic, smoothed, clipped, isotropic-mean normalized",
    }


def _interp_on_mids(x: np.ndarray, mids: np.ndarray, values: np.ndarray) -> np.ndarray:
    # Piecewise-linear interpolation with zero outside range.
    x = np.asarray(x, dtype=float)
    mids = np.asarray(mids, dtype=float)
    values = np.asarray(values, dtype=float)
    return np.interp(x, mids, values, left=0.0, right=0.0)


def _axis_coords() -> np.ndarray:
    idx = np.arange(twompp.GRID_N, dtype=np.float32)
    return (idx - twompp.GRID_CENTER_INDEX) * twompp.GRID_SPACING_MPC_H


def _ensure_twompp_file_for_sigma(total_sigma_mpc_h: float) -> Path:
    total_sigma_mpc_h = float(total_sigma_mpc_h)
    if abs(total_sigma_mpc_h - twompp.BASE_SMOOTHING_GAUSSIAN_MPC_H) < 1e-9:
        base = twompp._base_density_path()
        twompp._download_if_missing(twompp.TWOMPP_DENSITY_URL, base)
        return base

    p = twompp._smoothed_density_path(total_sigma_mpc_h)
    if p.exists():
        return p
    # Compute and save.
    _ = twompp._load_density_field(total_sigma_mpc_h=total_sigma_mpc_h)
    return p


def read_beta_and_h0_from_results_v1(path: Path, *, env_id_local: str) -> Dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    beta = float(data["metrics"]["mc_beta_direct"]["p50"]) if data.get("metrics", {}).get("mc_beta_direct") else float(data["fit_direct"]["beta"])
    beta_p16 = float(data["metrics"]["mc_beta_direct"]["p16"]) if data.get("metrics", {}).get("mc_beta_direct") else None
    beta_p84 = float(data["metrics"]["mc_beta_direct"]["p84"]) if data.get("metrics", {}).get("mc_beta_direct") else None

    H0_global = float(data["metadata"]["H0_anchor"]) if data.get("metadata", {}).get("H0_anchor") is not None else float(data["fit_direct"]["H0_global"])

    # Find the local point by env_id in fit points.
    H0_local = None
    for p in data.get("fit_direct", {}).get("fit_points", []):
        if str(p.get("env_id")) == str(env_id_local):
            H0_local = float(p.get("H0"))
            break
    if H0_local is None:
        # Fallback: search predictions_mc list.
        for p in data.get("predictions_mc", []):
            if str(p.get("env_id")) == str(env_id_local):
                H0_local = float(p.get("H0"))
                break
    if H0_local is None:
        raise ValueError(f"Could not find env_id={env_id_local} in results_v1.json.")

    dH_over_H = (H0_local / H0_global) - 1.0
    delta_required = -dH_over_H / beta if beta != 0 else None
    delta_required_p16 = -dH_over_H / beta_p84 if (beta_p84 is not None and beta_p84 != 0) else None
    delta_required_p84 = -dH_over_H / beta_p16 if (beta_p16 is not None and beta_p16 != 0) else None

    return {
        "results_v1_path": str(path),
        "env_id_local": str(env_id_local),
        "H0_global": float(H0_global),
        "H0_local": float(H0_local),
        "dH_over_H": float(dH_over_H),
        "beta_p50": float(beta),
        "beta_p16": float(beta_p16) if beta_p16 is not None else None,
        "beta_p84": float(beta_p84) if beta_p84 is not None else None,
        "delta_required_p50": float(delta_required) if delta_required is not None else None,
        "delta_required_p16": float(delta_required_p16) if delta_required_p16 is not None else None,
        "delta_required_p84": float(delta_required_p84) if delta_required_p84 is not None else None,
        "note": "delta_required is the delta (in the same model equation) needed to explain H0_local/H0_global given beta. This is not delta_g*.",
    }


def compute_delta_eff_for_models(
    *,
    radial_weight: Dict[str, Any],
    lat_weight: Optional[Dict[str, Any]],
    b_cut_deg: float,
    r_shell_min: float,
    r_shell_max: float,
    smoothing_sigmas_mpc_h: List[float],
) -> Dict[str, Any]:
    xs = _axis_coords().astype(np.float32)
    ys = xs
    N = int(twompp.GRID_N)

    r_mids = np.asarray(radial_weight["r_mids_mpc_h"], dtype=float)
    w_r_mids = np.asarray(radial_weight["w_r_mids"], dtype=float)

    if lat_weight is not None:
        absb_mids = np.asarray(lat_weight["absb_mids_deg"], dtype=float)
        w_lat_mids = np.asarray(lat_weight["w_lat_mids"], dtype=float)
    else:
        absb_mids = None
        w_lat_mids = None

    # Accumulators: denom depends only on weights; numerator depends on density field (sigma).
    model_keys = ["uniform_shell", "sh0es_radial_fullsky", "sh0es_radial_galcut", "sh0es_radial_latweight"]

    denoms: Dict[str, float] = {k: 0.0 for k in model_keys}
    numerators: Dict[str, Dict[str, float]] = {k: {} for k in model_keys}
    numer2: Dict[str, Dict[str, float]] = {k: {} for k in model_keys}  # for weighted variance

    density_maps: Dict[str, np.ndarray] = {}
    sigma_paths: Dict[str, str] = {}
    for s in smoothing_sigmas_mpc_h:
        p = _ensure_twompp_file_for_sigma(float(s))
        density_maps[str(float(s))] = np.load(p, mmap_mode="r")  # type: ignore
        sigma_paths[str(float(s))] = str(p)
        for k in model_keys:
            numerators[k][str(float(s))] = 0.0
            numer2[k][str(float(s))] = 0.0

    # Loop over z-slices
    for kk in range(N):
        z = float(xs[kk])
        z2 = np.float32(z * z)

        # r(x,y; z)
        r2 = (xs[:, None] ** 2) + (ys[None, :] ** 2) + z2
        r = np.sqrt(r2, dtype=np.float32)

        # Galactic latitude b = asin(z/r)
        with np.errstate(divide="ignore", invalid="ignore"):
            brad = np.arcsin(np.where(r > 0, z / r, 0.0), dtype=np.float32)
        b_deg = (brad * (180.0 / np.pi)).astype(np.float32, copy=False)
        absb = np.abs(b_deg).astype(np.float32, copy=False)

        # Common masks / weights
        shell_mask = (r >= float(r_shell_min)) & (r < float(r_shell_max))
        w_shell = shell_mask.astype(np.float32, copy=False)  # 0/1

        wr = _interp_on_mids(r.astype(float, copy=False), r_mids, w_r_mids).astype(np.float32, copy=False)
        # Impose the same r-range support as the histogram
        wr = np.where(shell_mask, wr, 0.0).astype(np.float32, copy=False)

        bmask = (absb >= float(b_cut_deg))
        w_galcut = (wr * bmask.astype(np.float32, copy=False)).astype(np.float32, copy=False)

        if absb_mids is not None and w_lat_mids is not None:
            wlat = _interp_on_mids(absb.astype(float, copy=False), absb_mids, w_lat_mids).astype(np.float32, copy=False)
            w_lat = (wr * wlat).astype(np.float32, copy=False)
        else:
            w_lat = np.zeros_like(wr)

        weight_slices = {
            "uniform_shell": w_shell,
            "sh0es_radial_fullsky": wr,
            "sh0es_radial_galcut": w_galcut,
            "sh0es_radial_latweight": w_lat,
        }

        # Update denominators (once; independent of sigma).
        for mk in model_keys:
            denoms[mk] += float(np.sum(weight_slices[mk]))

        # Update numerators for each sigma.
        for s in smoothing_sigmas_mpc_h:
            skey = str(float(s))
            dens = density_maps[skey][:, :, kk].astype(np.float32, copy=False)
            for mk in model_keys:
                w = weight_slices[mk]
                numerators[mk][skey] += float(np.sum(dens * w))
                numer2[mk][skey] += float(np.sum((dens * dens) * w))

    # Finalize means/stds
    out: Dict[str, Any] = {
        "sigma_paths": sigma_paths,
        "models": {},
    }
    for mk in model_keys:
        model_out: Dict[str, Any] = {
            "denom_sumw": float(denoms[mk]),
            "by_sigma": {},
        }
        for s in smoothing_sigmas_mpc_h:
            skey = str(float(s))
            denom = float(denoms[mk])
            if denom <= 0:
                mean = None
                std = None
            else:
                mean = float(numerators[mk][skey] / denom)
                # Var_w = E[w x^2]/E[w] - mean^2
                ex2 = float(numer2[mk][skey] / denom)
                var = max(0.0, ex2 - mean * mean)
                std = float(math.sqrt(var))
            model_out["by_sigma"][skey] = {"delta_g_star_eff_mean": mean, "delta_g_star_eff_std": std}
        out["models"][mk] = model_out
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="2M++ observer-weighted delta_g* (v0)")
    ap.add_argument("--results-v1", type=str, default=str(DEFAULT_RESULTS_V1), help="Path to results_v1.json (for beta/H0).")
    ap.add_argument("--env-id-local", type=str, default="sh0es_local", help="Local env_id to compare against (in results_v1.json).")

    ap.add_argument("--data-file", type=str, default="", help="Optional local Pantheon+SH0ES.dat path; if empty, download.")
    ap.add_argument("--z-min", type=float, default=0.015)
    ap.add_argument("--z-max", type=float, default=0.060)
    ap.add_argument("--z-field", type=str, default="zCMB", choices=["zCMB", "zHD"])
    ap.add_argument(
        "--only-sh0es-hf",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use only USED_IN_SH0ES_HF==1 (recommended).",
    )
    ap.add_argument(
        "--exclude-calibrators",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Exclude IS_CALIBRATOR==1 (recommended).",
    )
    ap.add_argument(
        "--dedup-cid",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Deduplicate by CID (recommended).",
    )

    ap.add_argument("--r-shell-min", type=float, default=40.0, help="Observer-centered rmin for integration support [Mpc/h].")
    ap.add_argument("--r-shell-max", type=float, default=200.0, help="Observer-centered rmax for integration support [Mpc/h].")
    ap.add_argument("--r-bin-width", type=float, default=5.0, help="Radial bin width for building w_r(r) [Mpc/h].")
    ap.add_argument("--r-smooth-sigma-bins", type=float, default=1.0, help="Gaussian smoothing sigma (in bins) for w_r.")

    ap.add_argument("--b-cut-deg", type=float, default=15.0, help="Galactic latitude cut |b|>=b_cut for one model [deg].")
    ap.add_argument("--lat-bins", type=int, default=18, help="Bins for latitude weight on |b| in [0,90].")
    ap.add_argument("--lat-smooth-sigma-bins", type=float, default=1.0, help="Gaussian smoothing sigma (in bins) for w_lat.")
    ap.add_argument("--lat-clip-min", type=float, default=0.0, help="Min clip for w_lat.")
    ap.add_argument("--lat-clip-max", type=float, default=5.0, help="Max clip for w_lat.")

    ap.add_argument("--smooth-list", type=str, default="4,20,50,100", help="Comma-separated 2M++ total smoothing sigmas [Mpc/h].")
    ap.add_argument("--output-json", type=str, default=str(OUTPUT_JSON))
    ap.add_argument("--output-plot", type=str, default=str(OUTPUT_PLOT))
    args = ap.parse_args()

    results_v1_path = Path(str(args.results_v1))
    ref = read_beta_and_h0_from_results_v1(results_v1_path, env_id_local=str(args.env_id_local))

    # Pantheon+ input
    if str(args.data_file).strip():
        data_path = Path(str(args.data_file))
    else:
        data_path = _cache_dir() / "Pantheon+SH0ES.dat"
        _download_if_missing(PANTHEONPLUS_URL, data_path)

    rows = read_pantheonplus_dat(data_path)
    n_raw = len(rows)

    if bool(args.exclude_calibrators):
        rows = [r for r in rows if int(r.is_calibrator) == 0]
    if bool(args.only_sh0es_hf):
        rows = [r for r in rows if int(r.used_in_sh0es_hf) == 1]
    if bool(args.dedup_cid):
        rows = dedup_by_cid(rows)

    # z filter for geometry window, and enforce 2M++ extent via z_max choice.
    z_min = float(args.z_min)
    z_max = float(args.z_max)
    if str(args.z_field) == "zHD":
        rows = [r for r in rows if z_min <= float(r.zhd) <= z_max]
        z_arr = np.array([r.zhd for r in rows], dtype=float)
    else:
        rows = [r for r in rows if z_min <= float(r.zcmb) <= z_max]
        z_arr = np.array([r.zcmb for r in rows], dtype=float)

    if rows:
        # Convert RA/Dec -> Galactic latitude b (deg) for selection footprint proxy.
        b_list: List[float] = []
        for r in rows:
            _l, b = twompp._radec_to_galactic_deg(float(r.ra_deg), float(r.dec_deg))
            b_list.append(float(b))
        b_arr = np.array(b_list, dtype=float)
    else:
        b_arr = np.array([], dtype=float)

    if rows is None or len(rows) < 50:
        raise SystemExit("Not enough Pantheon+ rows after filtering; relax z-range or flags.")

    r_arr = z_to_r_mpc_h(z_arr)

    radial_weight = build_radial_weight_from_sample(
        r_arr,
        rmin=float(args.r_shell_min),
        rmax=float(args.r_shell_max),
        bin_width=float(args.r_bin_width),
        smooth_sigma_bins=float(args.r_smooth_sigma_bins),
    )

    lat_weight = build_latitude_weight_from_sample(
        b_arr,
        n_bins=int(args.lat_bins),
        smooth_sigma_bins=float(args.lat_smooth_sigma_bins),
        clip_min=float(args.lat_clip_min),
        clip_max=float(args.lat_clip_max),
    )

    # Smoothing list
    smoothing_sigmas: List[float] = []
    for t in str(args.smooth_list).split(","):
        tt = t.strip()
        if not tt:
            continue
        try:
            smoothing_sigmas.append(float(tt))
        except Exception:
            continue
    if not smoothing_sigmas:
        smoothing_sigmas = [4.0, 20.0, 50.0, 100.0]

    # Compute delta_eff
    model_results = compute_delta_eff_for_models(
        radial_weight=radial_weight,
        lat_weight=lat_weight,
        b_cut_deg=float(args.b_cut_deg),
        r_shell_min=float(args.r_shell_min),
        r_shell_max=float(args.r_shell_max),
        smoothing_sigmas_mpc_h=smoothing_sigmas,
    )

    # Attach prediction mapping using beta (note: delta_g* is not delta_m; this is a diagnostic).
    beta = float(ref["beta_p50"])
    H0_global = float(ref["H0_global"])
    for mk, md in model_results["models"].items():
        for skey, rec in md["by_sigma"].items():
            d = rec.get("delta_g_star_eff_mean")
            if d is None or beta == 0:
                rec["predicted_dH_over_H_using_delta_g_star"] = None
                rec["predicted_H0_using_delta_g_star"] = None
            else:
                pred = -beta * float(d)
                rec["predicted_dH_over_H_using_delta_g_star"] = float(pred)
                rec["predicted_H0_using_delta_g_star"] = float(H0_global * (1.0 + pred))

    out: Dict[str, Any] = {
        "metadata": {"timestamp": datetime.now().isoformat(), "script": Path(__file__).name},
        "inputs": {
            "pantheonplus_url": PANTHEONPLUS_URL,
            "pantheonplus_file": str(data_path),
            "results_v1": ref,
            "filters": {
                "z_field": str(args.z_field),
                "z_min": float(args.z_min),
                "z_max": float(args.z_max),
                "only_sh0es_hf": bool(args.only_sh0es_hf),
                "exclude_calibrators": bool(args.exclude_calibrators),
                "dedup_cid": bool(args.dedup_cid),
            },
            "counts": {"pantheonplus_rows_raw": int(n_raw), "pantheonplus_rows_used": int(len(rows))},
            "weight_config": {
                "r_shell_min": float(args.r_shell_min),
                "r_shell_max": float(args.r_shell_max),
                "b_cut_deg": float(args.b_cut_deg),
            },
            "twompp": {
                "grid_N": int(twompp.GRID_N),
                "grid_spacing_mpc_h": float(twompp.GRID_SPACING_MPC_H),
                "grid_extent_mpc_h": float(twompp.GRID_HALF_EXTENT_MPC_H),
                "base_smoothing_mpc_h": float(twompp.BASE_SMOOTHING_GAUSSIAN_MPC_H),
            },
        },
        "weights": {"radial": radial_weight, "latitude_axisymmetric": lat_weight},
        "results": model_results,
        "outputs": {"json": str(Path(args.output_json)), "plot": str(Path(args.output_plot))},
        "notes": [
            "delta_eff is computed as a 3D weighted volume average of 2M++ delta_g* around the origin.",
            "w_r(r) is derived from Pantheon+ SH0ES-HF sample radial counts (as a window proxy).",
            "w_lat(|b|) is an axisymmetric latitude-only proxy for sky coverage; real footprint is more complex.",
            "Predicted dH/H uses the model beta but treats delta_g* as delta (diagnostic only).",
        ],
    }

    # Write JSON
    Path(args.output_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_json).write_text(json.dumps(_json_sanitize(out), indent=2, ensure_ascii=False), encoding="utf-8")

    # Plot
    sigmas = [float(s) for s in smoothing_sigmas]
    sigmas_sorted = sorted(sigmas)
    x = np.array(sigmas_sorted, dtype=float)

    req = ref.get("delta_required_p50")
    req16 = ref.get("delta_required_p16")
    req84 = ref.get("delta_required_p84")
    obs_dH = float(ref["dH_over_H"])

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), sharex=False)
    ax0, ax1 = axes

    styles = {
        "uniform_shell": {"label": f"uniform shell [{args.r_shell_min:.0f},{args.r_shell_max:.0f}] full-sky", "color": "#444444", "marker": "o"},
        "sh0es_radial_fullsky": {"label": "w_r(r) from SH0ES-HF (full-sky)", "color": "#1f77b4", "marker": "o"},
        "sh0es_radial_galcut": {"label": f"w_r(r) + |b|>={args.b_cut_deg:.0f} deg", "color": "#ff7f0e", "marker": "o"},
        "sh0es_radial_latweight": {"label": "w_r(r) * w_lat(|b|) (axisym)", "color": "#2ca02c", "marker": "o"},
    }

    for mk, st in styles.items():
        y = []
        ypred = []
        for s in sigmas_sorted:
            skey = str(float(s))
            rec = out["results"]["models"][mk]["by_sigma"][skey]
            y.append(rec["delta_g_star_eff_mean"])
            ypred.append(rec["predicted_dH_over_H_using_delta_g_star"])
        ax0.plot(x, y, marker=st["marker"], color=st["color"], label=st["label"], linewidth=1.6)
        ax1.plot(x, ypred, marker=st["marker"], color=st["color"], label=st["label"], linewidth=1.6)

    if req is not None:
        ax0.axhline(float(req), color="#111111", linestyle="--", linewidth=1.2, label="delta_required (p50)")
        if req16 is not None and req84 is not None:
            ax0.fill_between(x, float(req16), float(req84), color="#111111", alpha=0.10, label="delta_required (p16-p84)")

    ax1.axhline(obs_dH, color="#111111", linestyle="--", linewidth=1.2, label="observed dH/H (SH0ES vs Planck)")

    ax0.set_title("2M++ delta_g* effective observer overdensity")
    ax0.set_xlabel("2M++ total smoothing sigma (Mpc/h)")
    ax0.set_ylabel("delta_g*_eff (weighted volume mean)")
    ax0.grid(alpha=0.25)
    ax0.legend(fontsize=8)

    ax1.set_title("Predicted dH/H using beta and delta_g* (diagnostic)")
    ax1.set_xlabel("2M++ total smoothing sigma (Mpc/h)")
    ax1.set_ylabel("predicted dH/H = -beta * delta_g*_eff")
    ax1.grid(alpha=0.25)
    ax1.legend(fontsize=8)

    fig.tight_layout()
    Path(args.output_plot).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(args.output_plot), dpi=160)
    plt.close(fig)

    print("Observer-weighted delta_g* completed.")
    print(f"Pantheon+ rows used: {len(rows)} (raw={n_raw})")
    print(f"Output JSON: {args.output_json}")
    print(f"Output plot: {args.output_plot}")


if __name__ == "__main__":
    main()

