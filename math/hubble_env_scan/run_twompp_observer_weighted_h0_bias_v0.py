"""
2M++ observer-weighted H0 bias from peculiar velocity field (v0)
==============================================================

Goal
----
Use a *gravitational-response proxy* (peculiar velocity field) to estimate the
observer-window-induced bias on an H0 fit, under a given window w(r, n).

This is motivated by the "hidden relations" route:
- If luminous tracers (delta_g*) are not the right proxy for the relational/mass
  environment, then a velocity-field proxy (which responds to total gravitating
  structure) is a better observable handle.

We compute a simple weighted least-squares slope estimator:

  v_obs = 100 * r + v_r,   with r in (Mpc/h), v in (km/s),

  H_fit_tilde = sum(w * r * v_obs) / sum(w * r^2)
              = 100 + sum(w * r * v_r) / sum(w * r^2),

and report the fractional bias:

  deltaH_over_H = (H_fit_tilde - 100) / 100.

Inputs
------
- Pantheon+SH0ES distance table (to build an empirical w(r) and a crude sky proxy)
- 2M++ velocity field (Carrick+2015) in numpy format: twompp_velocity.npy
- Optional: results_v1.json (for the observed Planck→SH0ES dH/H reference)

Outputs
-------
- twompp_observer_weighted_h0_bias_results_v0.json
- twompp_observer_weighted_h0_bias_plot_v0.png

Run (repo root)
---------------
python math/hubble_env_scan/run_twompp_observer_weighted_h0_bias_v0.py

Notes
-----
- The velocity field is expected in Galactic Cartesian coordinates on the same
  257^3 cube as twompp_density.npy, spanning [-200,200] Mpc/h.
- This is still a proxy (model-based reconstruction). It is nonetheless closer
  to "gravitational mass response" than raw luminous density.
- Keep console output ASCII-friendly (Windows terminals).
"""

from __future__ import annotations

import argparse
import json
import math
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

import run_twompp_observer_weighted_delta_v0 as wdelta
import twompp_density_lookup_v0 as twompp


TWOMPP_VELOCITY_URL = "https://cosmicflows.iap.fr/assets/data/twompp_velocity.npy"

DEFAULT_RESULTS_V1 = Path(__file__).parent / "results_v1.json"

OUTPUT_JSON = Path(__file__).parent / "twompp_observer_weighted_h0_bias_results_v0.json"
OUTPUT_PLOT = Path(__file__).parent / "twompp_observer_weighted_h0_bias_plot_v0.png"


def _cache_dir() -> Path:
    d = Path(__file__).parent / "_cache"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _base_velocity_path() -> Path:
    return _cache_dir() / "twompp_velocity.npy"


def _download_if_missing(url: str, path: Path) -> None:
    if path.exists():
        return
    print(f"Downloading 2M++ velocity field -> {path}")
    print(f"Source: {url}")
    urllib.request.urlretrieve(url, path)  # noqa: S310 (intentional external download)


def _load_velocity_field_mmap() -> np.ndarray:
    path = _base_velocity_path()
    _download_if_missing(TWOMPP_VELOCITY_URL, path)
    arr = np.load(path, mmap_mode="r")
    return arr


def _axis_coords() -> np.ndarray:
    idx = np.arange(twompp.GRID_N, dtype=np.float32)
    return (idx - twompp.GRID_CENTER_INDEX) * twompp.GRID_SPACING_MPC_H


def _vel_components_slice(vfield: np.ndarray, k: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    # Support both common layouts:
    # - (N,N,N,3) with last axis [vx,vy,vz]
    # - (3,N,N,N) with first axis [vx,vy,vz]
    if vfield.ndim != 4:
        raise ValueError(f"Unexpected velocity field ndim={vfield.ndim}; expected 4.")

    N = int(twompp.GRID_N)
    if vfield.shape[:3] == (N, N, N) and vfield.shape[3] == 3:
        vx = np.asarray(vfield[:, :, k, 0])
        vy = np.asarray(vfield[:, :, k, 1])
        vz = np.asarray(vfield[:, :, k, 2])
        return vx, vy, vz

    if vfield.shape[0] == 3 and vfield.shape[1:] == (N, N, N):
        vx = np.asarray(vfield[0, :, :, k])
        vy = np.asarray(vfield[1, :, :, k])
        vz = np.asarray(vfield[2, :, :, k])
        return vx, vy, vz

    raise ValueError(f"Unexpected velocity field shape={vfield.shape}; cannot infer component layout.")


def compute_h0_bias(
    *,
    vfield: np.ndarray,
    radial_weight: Dict[str, Any],
    lat_weight: Dict[str, Any],
    b_cut_deg: float,
    r_shell_min: float,
    r_shell_max: float,
) -> Dict[str, Any]:
    xs = _axis_coords().astype(np.float32)
    ys = xs
    N = int(twompp.GRID_N)

    r_mids = np.asarray(radial_weight["r_mids_mpc_h"], dtype=float)
    w_r_mids = np.asarray(radial_weight["w_r_mids"], dtype=float)

    absb_mids = np.asarray(lat_weight["absb_mids_deg"], dtype=float)
    w_lat_mids = np.asarray(lat_weight["w_lat_mids"], dtype=float)

    model_keys = ["uniform_shell", "sh0es_radial_fullsky", "sh0es_radial_galcut", "sh0es_radial_latweight"]

    # Regression accumulators:
    # S_rr = sum(w_obj * r^2), S_rvr = sum(w_obj * r * v_r)
    S_rr: Dict[str, float] = {k: 0.0 for k in model_keys}
    S_rvr: Dict[str, float] = {k: 0.0 for k in model_keys}

    # Optional diagnostic: mean(v_r/r) under object measure
    S_w: Dict[str, float] = {k: 0.0 for k in model_keys}
    S_vr_over_r: Dict[str, float] = {k: 0.0 for k in model_keys}

    for kk in range(N):
        z = float(xs[kk])
        z2 = np.float32(z * z)

        r2 = (xs[:, None] ** 2) + (ys[None, :] ** 2) + z2
        r = np.sqrt(r2, dtype=np.float32)

        with np.errstate(divide="ignore", invalid="ignore"):
            brad = np.arcsin(np.where(r > 0, z / r, 0.0), dtype=np.float32)
        b_deg = (brad * (180.0 / np.pi)).astype(np.float32, copy=False)
        absb = np.abs(b_deg).astype(np.float32, copy=False)

        shell_mask = (r >= float(r_shell_min)) & (r < float(r_shell_max))
        w_shell = shell_mask.astype(np.float32, copy=False)

        wr = wdelta._interp_on_mids(r.astype(float, copy=False), r_mids, w_r_mids).astype(np.float32, copy=False)
        wr = np.where(shell_mask, wr, 0.0).astype(np.float32, copy=False)

        bmask = (absb >= float(b_cut_deg))
        w_galcut = (wr * bmask.astype(np.float32, copy=False)).astype(np.float32, copy=False)

        wlat = wdelta._interp_on_mids(absb.astype(float, copy=False), absb_mids, w_lat_mids).astype(np.float32, copy=False)
        w_lat = (wr * wlat).astype(np.float32, copy=False)

        weight_slices = {
            "uniform_shell": w_shell,
            "sh0es_radial_fullsky": wr,
            "sh0es_radial_galcut": w_galcut,
            "sh0es_radial_latweight": w_lat,
        }

        # Velocity radial component at slice
        vx, vy, vz = _vel_components_slice(vfield, kk)
        vx = vx.astype(np.float32, copy=False)
        vy = vy.astype(np.float32, copy=False)
        vz = vz.astype(np.float32, copy=False)

        # v_r = (v·xhat) with xhat = x/r
        # Broadcast x,y coords
        x = xs[:, None].astype(np.float32, copy=False)
        y = ys[None, :].astype(np.float32, copy=False)
        zz = np.float32(z)
        with np.errstate(divide="ignore", invalid="ignore"):
            inv_r = np.where(r > 0, 1.0 / r, 0.0).astype(np.float32, copy=False)
        vr = (vx * x + vy * y + vz * zz) * inv_r

        r_f = r.astype(np.float32, copy=False)
        r2_f = (r_f * r_f).astype(np.float32, copy=False)

        for mk in model_keys:
            w = weight_slices[mk]
            # Object measure sums
            S_w[mk] += float(np.sum(w))
            S_vr_over_r[mk] += float(np.sum(w * (vr * inv_r)))  # vr/r

            # Regression sums
            S_rr[mk] += float(np.sum(w * r2_f))
            S_rvr[mk] += float(np.sum(w * (r_f * vr)))

    results: Dict[str, Any] = {}
    for mk in model_keys:
        denom = float(S_rr[mk])
        if denom <= 0:
            results[mk] = {"error": "zero_denom"}
            continue

        deltaH_tilde = float(S_rvr[mk] / denom)  # km/s/(Mpc/h)
        H_fit_tilde = 100.0 + deltaH_tilde
        deltaH_over_H = float(deltaH_tilde / 100.0)

        mean_vr_over_r = float(S_vr_over_r[mk] / S_w[mk]) if float(S_w[mk]) > 0 else None

        results[mk] = {
            "S_rr": float(S_rr[mk]),
            "S_rvr": float(S_rvr[mk]),
            "deltaH_tilde_kms_per_Mpch": float(deltaH_tilde),
            "H_fit_tilde_kms_per_Mpch": float(H_fit_tilde),
            "deltaH_over_H": float(deltaH_over_H),
            "mean_vr_over_r_kms_per_Mpch2": float(mean_vr_over_r) if mean_vr_over_r is not None else None,
            "note_units": "r in Mpc/h, v in km/s. H_tilde is in km/s per (Mpc/h). Fractional bias is deltaH_over_H=(H_fit_tilde-100)/100.",
        }
    return results


def main() -> None:
    ap = argparse.ArgumentParser(description="2M++ observer-weighted H0 bias from velocity field (v0)")
    ap.add_argument("--results-v1", type=str, default=str(DEFAULT_RESULTS_V1), help="Path to results_v1.json (for observed dH/H reference).")
    ap.add_argument("--env-id-local", type=str, default="sh0es_local", help="Local env_id in results_v1.json for reference dH/H.")

    ap.add_argument("--data-file", type=str, default="", help="Optional local Pantheon+SH0ES.dat path; if empty, download (via the other script).")
    ap.add_argument("--z-min", type=float, default=0.015)
    ap.add_argument("--z-max", type=float, default=0.060)
    ap.add_argument("--z-field", type=str, default="zCMB", choices=["zCMB", "zHD"])
    ap.add_argument("--only-sh0es-hf", action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument("--exclude-calibrators", action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument("--dedup-cid", action=argparse.BooleanOptionalAction, default=True)

    ap.add_argument("--r-shell-min", type=float, default=40.0)
    ap.add_argument("--r-shell-max", type=float, default=200.0)
    ap.add_argument("--r-bin-width", type=float, default=5.0)
    ap.add_argument("--r-smooth-sigma-bins", type=float, default=1.0)

    ap.add_argument("--b-cut-deg", type=float, default=15.0)
    ap.add_argument("--lat-bins", type=int, default=18)
    ap.add_argument("--lat-smooth-sigma-bins", type=float, default=1.0)
    ap.add_argument("--lat-clip-min", type=float, default=0.0)
    ap.add_argument("--lat-clip-max", type=float, default=5.0)

    ap.add_argument("--output-json", type=str, default=str(OUTPUT_JSON))
    ap.add_argument("--output-plot", type=str, default=str(OUTPUT_PLOT))
    args = ap.parse_args()

    ref = wdelta.read_beta_and_h0_from_results_v1(Path(str(args.results_v1)), env_id_local=str(args.env_id_local))

    # Build weights from Pantheon+ subset (reuse helper script functions)
    if str(args.data_file).strip():
        data_path = Path(str(args.data_file))
    else:
        data_path = wdelta._cache_dir() / "Pantheon+SH0ES.dat"
        wdelta._download_if_missing(wdelta.PANTHEONPLUS_URL, data_path)

    rows = wdelta.read_pantheonplus_dat(data_path)
    n_raw = len(rows)
    if bool(args.exclude_calibrators):
        rows = [r for r in rows if int(r.is_calibrator) == 0]
    if bool(args.only_sh0es_hf):
        rows = [r for r in rows if int(r.used_in_sh0es_hf) == 1]
    if bool(args.dedup_cid):
        rows = wdelta.dedup_by_cid(rows)

    z_min = float(args.z_min)
    z_max = float(args.z_max)
    if str(args.z_field) == "zHD":
        rows = [r for r in rows if z_min <= float(r.zhd) <= z_max]
        z_arr = np.array([r.zhd for r in rows], dtype=float)
    else:
        rows = [r for r in rows if z_min <= float(r.zcmb) <= z_max]
        z_arr = np.array([r.zcmb for r in rows], dtype=float)

    if len(rows) < 50:
        raise SystemExit("Not enough Pantheon+ rows after filtering; relax z-range or flags.")

    b_list: List[float] = []
    for r in rows:
        _l, b = twompp._radec_to_galactic_deg(float(r.ra_deg), float(r.dec_deg))
        b_list.append(float(b))
    b_arr = np.array(b_list, dtype=float)

    r_arr = wdelta.z_to_r_mpc_h(z_arr)

    radial_weight = wdelta.build_radial_weight_from_sample(
        r_arr,
        rmin=float(args.r_shell_min),
        rmax=float(args.r_shell_max),
        bin_width=float(args.r_bin_width),
        smooth_sigma_bins=float(args.r_smooth_sigma_bins),
    )
    lat_weight = wdelta.build_latitude_weight_from_sample(
        b_arr,
        n_bins=int(args.lat_bins),
        smooth_sigma_bins=float(args.lat_smooth_sigma_bins),
        clip_min=float(args.lat_clip_min),
        clip_max=float(args.lat_clip_max),
    )

    vfield = _load_velocity_field_mmap()

    bias = compute_h0_bias(
        vfield=vfield,
        radial_weight=radial_weight,
        lat_weight=lat_weight,
        b_cut_deg=float(args.b_cut_deg),
        r_shell_min=float(args.r_shell_min),
        r_shell_max=float(args.r_shell_max),
    )

    # Attach comparisons to observed dH/H, and map to an implied delta (diagnostic) using beta.
    obs_dh = float(ref["dH_over_H"])
    beta = float(ref["beta_p50"])
    H0_global = float(ref["H0_global"])
    for mk, rec in bias.items():
        if rec.get("deltaH_over_H") is None:
            continue
        dh = float(rec["deltaH_over_H"])
        rec["predicted_H0_kms_per_Mpc"] = float(H0_global * (1.0 + dh))
        rec["residual_dH_over_H_needed"] = float(obs_dh - dh)
        rec["implied_delta_eff_if_using_beta"] = float(-dh / beta) if beta != 0 else None

    out: Dict[str, Any] = {
        "metadata": {"timestamp": datetime.now().isoformat(), "script": Path(__file__).name},
        "inputs": {
            "results_v1_ref": ref,
            "pantheonplus_file": str(data_path),
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
            "twompp_velocity_url": TWOMPP_VELOCITY_URL,
            "twompp_velocity_file": str(_base_velocity_path()),
            "twompp_grid": {
                "N": int(twompp.GRID_N),
                "spacing_mpc_h": float(twompp.GRID_SPACING_MPC_H),
                "extent_mpc_h": float(twompp.GRID_HALF_EXTENT_MPC_H),
            },
            "assumption": "Baseline Hubble flow in Mpc/h units is v=100*r (so H_tilde_true=100).",
        },
        "weights": {"radial": radial_weight, "latitude_axisymmetric": lat_weight},
        "results": {"h0_bias": bias},
        "outputs": {"json": str(Path(args.output_json)), "plot": str(Path(args.output_plot))},
        "notes": [
            "This uses the 2M++ reconstructed peculiar velocity field as a gravitational-response proxy.",
            "deltaH/H is computed from a weighted least-squares slope estimator for v_obs vs r.",
            "This is a diagnostic; it does not by itself validate any specific relational mechanism.",
        ],
    }

    Path(args.output_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_json).write_text(json.dumps(wdelta._json_sanitize(out), indent=2, ensure_ascii=False), encoding="utf-8")

    # Plot
    labels = {
        "uniform_shell": "uniform shell [40,200] full-sky",
        "sh0es_radial_fullsky": "w_r(r) from SH0ES-HF (full-sky)",
        "sh0es_radial_galcut": "w_r(r) + |b|>=15 deg",
        "sh0es_radial_latweight": "w_r(r) * w_lat(|b|) (axisym)",
    }
    colors = {
        "uniform_shell": "#444444",
        "sh0es_radial_fullsky": "#1f77b4",
        "sh0es_radial_galcut": "#ff7f0e",
        "sh0es_radial_latweight": "#2ca02c",
    }

    obs = float(ref["dH_over_H"])

    fig, ax = plt.subplots(1, 1, figsize=(8.4, 4.6))
    for mk, lab in labels.items():
        rec = bias.get(mk, {})
        if "deltaH_over_H" not in rec:
            continue
        ax.bar(lab, float(rec["deltaH_over_H"]), color=colors.get(mk, "#999999"), alpha=0.85)

    ax.axhline(0.0, color="#888888", linestyle="--", linewidth=1.0)
    ax.axhline(obs, color="#111111", linestyle="--", linewidth=1.2, label="observed dH/H (SH0ES vs Planck)")
    ax.set_ylabel("deltaH/H from 2M++ velocity field (window-weighted)")
    ax.set_title("Observer-window H0 bias proxy from 2M++ peculiar velocities")
    ax.grid(axis="y", alpha=0.25)
    ax.tick_params(axis="x", labelrotation=20)
    ax.legend(fontsize=9)
    fig.tight_layout()
    Path(args.output_plot).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(args.output_plot), dpi=160)
    plt.close(fig)

    print("Observer-weighted H0 bias from 2M++ velocity field completed.")
    print(f"Pantheon+ rows used: {len(rows)} (raw={n_raw})")
    print(f"Output JSON: {args.output_json}")
    print(f"Output plot: {args.output_plot}")


if __name__ == "__main__":
    main()

