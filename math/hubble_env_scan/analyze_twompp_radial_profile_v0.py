"""
Radial profile of the 2M++ density field around the observer (origin).

This script answers a key v2 sanity question:
  Does the public 2M++ reconstruction show an underdensity (or overdensity)
  around the origin on the relevant distance-ladder scales?

We compute:
  - shell mean δ_g* in radial bins (volume-weighted, grid-cell average)
  - cumulative mean δ_g* within a sphere of radius r

Input is the 2M++ delta_g* field (Carrick+2015), optionally with extra Gaussian
smoothing to reach a target total sigma (Mpc/h).
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

import matplotlib.pyplot as plt
import numpy as np

from twompp_density_lookup_v0 import (
    GRID_CENTER_INDEX,
    GRID_HALF_EXTENT_MPC_H,
    GRID_N,
    GRID_SPACING_MPC_H,
    _load_density_field,
)


HERE = Path(__file__).parent
DEFAULT_OUT_JSON = HERE / "twompp_radial_profile_results_v0.json"
DEFAULT_OUT_PLOT = HERE / "twompp_radial_profile_plot_v0.png"


def _axis_coords() -> np.ndarray:
    idx = np.arange(GRID_N, dtype=np.float32)
    return (idx - GRID_CENTER_INDEX) * float(GRID_SPACING_MPC_H)


def _radius_grid() -> np.ndarray:
    xs = _axis_coords()
    r2 = (xs[:, None, None] ** 2) + (xs[None, :, None] ** 2) + (xs[None, None, :] ** 2)
    return np.sqrt(r2, dtype=np.float32)


def compute_radial_profile(
    *,
    smoothing_sigma_mpc_h: float,
    bin_width_mpc_h: float,
    max_r_mpc_h: float,
) -> Dict[str, Any]:
    sig = float(smoothing_sigma_mpc_h)
    bw = float(bin_width_mpc_h)
    rmax = float(max_r_mpc_h)

    if bw <= 0:
        raise ValueError("bin_width_mpc_h must be > 0")
    if rmax <= 0 or rmax > GRID_HALF_EXTENT_MPC_H:
        raise ValueError(f"max_r_mpc_h must be in (0,{GRID_HALF_EXTENT_MPC_H}]")

    density = _load_density_field(total_sigma_mpc_h=sig).astype(np.float32, copy=False)
    r = _radius_grid()

    # Flatten and filter to r<rmax
    rf = r.reshape(-1)
    df = density.reshape(-1)
    m = rf < rmax
    rf = rf[m]
    df = df[m]

    n_bins = int(np.ceil(rmax / bw))
    # Uniform bins: bin index by floor(r/bw)
    bin_idx = np.floor(rf / bw).astype(np.int32)
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)

    counts = np.bincount(bin_idx, minlength=n_bins).astype(np.int64)
    sums = np.bincount(bin_idx, weights=df, minlength=n_bins).astype(np.float64)
    sums2 = np.bincount(bin_idx, weights=(df.astype(np.float64) ** 2), minlength=n_bins).astype(np.float64)

    with np.errstate(invalid="ignore", divide="ignore"):
        means = sums / counts
        vars_ = (sums2 / counts) - (means**2)
        stds = np.sqrt(np.maximum(vars_, 0.0))

    # Bin centers
    r_edges = np.arange(0.0, n_bins * bw + 1e-9, bw, dtype=np.float64)
    r_mids = (r_edges[:-1] + r_edges[1:]) / 2.0

    # Cumulative within sphere of radius r_edge[k]
    cum_counts = np.cumsum(counts).astype(np.int64)
    cum_sums = np.cumsum(sums).astype(np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        cum_means = cum_sums / cum_counts

    shell_rows: List[Dict[str, Any]] = []
    for i in range(n_bins):
        shell_rows.append(
            {
                "bin": int(i),
                "rmin_mpc_h": float(r_edges[i]),
                "rmax_mpc_h": float(r_edges[i + 1]),
                "rmid_mpc_h": float(r_mids[i]),
                "n_cells": int(counts[i]),
                "delta_g_star_mean": float(means[i]) if np.isfinite(means[i]) else None,
                "delta_g_star_std": float(stds[i]) if np.isfinite(stds[i]) else None,
            }
        )

    cum_rows: List[Dict[str, Any]] = []
    for i in range(n_bins):
        cum_rows.append(
            {
                "rmax_mpc_h": float(r_edges[i + 1]),
                "n_cells": int(cum_counts[i]),
                "delta_g_star_mean": float(cum_means[i]) if np.isfinite(cum_means[i]) else None,
            }
        )

    return {
        "config": {"smoothing_sigma_mpc_h": sig, "bin_width_mpc_h": bw, "max_r_mpc_h": rmax},
        "shell_profile": shell_rows,
        "cumulative_profile": cum_rows,
        "global_field_stats": {
            "mean": float(np.mean(density)),
            "std": float(np.std(density)),
            "min": float(np.min(density)),
            "max": float(np.max(density)),
        },
    }


def plot_profile(results: Dict[str, Any], *, out_plot: Path) -> None:
    shell = results["shell_profile"]
    cum = results["cumulative_profile"]

    rmid = np.array([r["rmid_mpc_h"] for r in shell], dtype=float)
    smean = np.array([r["delta_g_star_mean"] if r["delta_g_star_mean"] is not None else np.nan for r in shell], dtype=float)
    sstd = np.array([r["delta_g_star_std"] if r["delta_g_star_std"] is not None else np.nan for r in shell], dtype=float)

    rcum = np.array([r["rmax_mpc_h"] for r in cum], dtype=float)
    cmean = np.array([r["delta_g_star_mean"] if r["delta_g_star_mean"] is not None else np.nan for r in cum], dtype=float)

    sig = results.get("config", {}).get("smoothing_sigma_mpc_h")
    bw = results.get("config", {}).get("bin_width_mpc_h")

    fig, axes = plt.subplots(2, 1, figsize=(9.6, 6.8), sharex=True, gridspec_kw={"height_ratios": [1, 1]})

    axes[0].plot(rmid, smean, color="#1f77b4", linewidth=1.5)
    axes[0].fill_between(rmid, smean - sstd, smean + sstd, color="#1f77b4", alpha=0.18, label="shell mean ± 1σ (field scatter)")
    axes[0].axhline(0.0, color="#999999", linestyle="--", linewidth=0.9)
    axes[0].set_ylabel("Shell mean δ_g*")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.2)

    axes[1].plot(rcum, cmean, color="#d62728", linewidth=1.6, label="cumulative mean (<r)")
    axes[1].axhline(0.0, color="#999999", linestyle="--", linewidth=0.9)
    axes[1].set_xlabel("Radius r (Mpc/h)")
    axes[1].set_ylabel("Cumulative mean δ_g*")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.2)

    fig.suptitle(f"2M++ δ_g* radial profile (σ={sig} Mpc/h, bin={bw} Mpc/h)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_plot, dpi=150)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="2M++ radial profile around origin (v0).")
    ap.add_argument("--smooth", type=float, default=4.0, help="Total Gaussian smoothing σ [Mpc/h] (>=4)")
    ap.add_argument("--bin-width", type=float, default=10.0, help="Radial bin width [Mpc/h]")
    ap.add_argument("--max-r", type=float, default=200.0, help="Max radius [Mpc/h] (<=200)")
    ap.add_argument("--output-json", type=str, default=str(DEFAULT_OUT_JSON))
    ap.add_argument("--output-plot", type=str, default=str(DEFAULT_OUT_PLOT))
    args = ap.parse_args()

    res = compute_radial_profile(smoothing_sigma_mpc_h=args.smooth, bin_width_mpc_h=args.bin_width, max_r_mpc_h=args.max_r)
    res["metadata"] = {
        "timestamp": datetime.now().isoformat(),
        "script": Path(__file__).name,
        "grid": {"N": GRID_N, "spacing_mpc_h": GRID_SPACING_MPC_H, "extent_mpc_h": GRID_HALF_EXTENT_MPC_H},
    }

    out_json = Path(args.output_json)
    out_plot = Path(args.output_plot)
    out_json.parent.mkdir(parents=True, exist_ok=True)

    with out_json.open("w", encoding="utf-8") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)

    plot_profile(res, out_plot=out_plot)

    print("2M++ radial profile completed.")
    print(f"Output plot: {out_plot}")
    print(f"Output JSON: {out_json}")


if __name__ == "__main__":
    main()

