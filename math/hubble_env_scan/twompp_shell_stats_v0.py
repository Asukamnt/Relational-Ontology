"""
2M++ density shell statistics (v0)
=================================

Goal
----
Compute an "observer-based environment" scalar from a public density field by
averaging the reconstructed density contrast within a spherical shell around
the Local Group (origin).

This matches the kind of question we care about for H0 measurements:
"What is the mean overdensity around the observer on the distance scales that
enter the distance ladder / Hubble-flow fit?"

Notes
-----
- Uses the 2M++ luminosity-weighted density contrast delta_g* on a 257^3 cube.
- Optionally applies extra Gaussian smoothing to reach a target sigma (Mpc/h).
- The 2M++ published product is already smoothed at sigma=4 Mpc/h.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from twompp_density_lookup_v0 import (
    GRID_CENTER_INDEX,
    GRID_HALF_EXTENT_MPC_H,
    GRID_N,
    GRID_SPACING_MPC_H,
    lookup_delta_g_star,  # for consistency reference
    _load_density_field,
)


OUTPUT_JSON = Path(__file__).parent / "twompp_shell_stats_v0.json"


def _axis_coords() -> np.ndarray:
    # X = (i-128)*400/256; grid spacing is 400/256
    idx = np.arange(GRID_N, dtype=np.float32)
    return (idx - GRID_CENTER_INDEX) * GRID_SPACING_MPC_H


def shell_mean_delta_g_star(rmin: float, rmax: float, smoothing_gaussian_mpc_h: float) -> dict:
    rmin = float(rmin)
    rmax = float(rmax)
    if rmin < 0 or rmax <= 0 or rmax <= rmin:
        raise ValueError("Require 0 <= rmin < rmax.")
    if rmax > GRID_HALF_EXTENT_MPC_H:
        raise ValueError(f"Require rmax <= {GRID_HALF_EXTENT_MPC_H} Mpc/h so the full sphere fits in the cube.")

    density = _load_density_field(total_sigma_mpc_h=float(smoothing_gaussian_mpc_h))

    xs = _axis_coords()
    ys = xs
    zs = xs

    # Radius^2 via broadcasting (float32 to reduce memory)
    r2 = (xs[:, None, None] ** 2) + (ys[None, :, None] ** 2) + (zs[None, None, :] ** 2)
    r = np.sqrt(r2, dtype=np.float32)
    mask = (r >= rmin) & (r < rmax)

    vals = density[mask]
    mean = float(np.mean(vals)) if vals.size else float("nan")
    std = float(np.std(vals)) if vals.size else float("nan")

    # Reference: point value at the origin (arbitrary l,b when r=0)
    origin_point = lookup_delta_g_star(l_deg=0.0, b_deg=0.0, r_mpc_h=0.0, smoothing_gaussian_mpc_h=float(smoothing_gaussian_mpc_h))

    return {
        "rmin_mpc_h": rmin,
        "rmax_mpc_h": rmax,
        "smoothing_gaussian_mpc_h": float(smoothing_gaussian_mpc_h),
        "delta_g_star_mean": mean,
        "delta_g_star_std": std,
        "n_cells": int(vals.size),
        "origin_delta_g_star": origin_point["delta_g_star"],
        "grid": {
            "N": GRID_N,
            "spacing_mpc_h": GRID_SPACING_MPC_H,
            "extent_mpc_h": GRID_HALF_EXTENT_MPC_H,
        },
        "notes": [
            "delta_g* is luminosity-weighted galaxy density contrast from 2M++ (Carrick et al. 2015).",
            "Shell mean is a simple volume average over grid cells whose centers fall in [rmin,rmax).",
            "This is not directly delta_m unless a bias model is assumed.",
        ],
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rmin", type=float, required=True, help="Inner radius of shell [Mpc/h]")
    ap.add_argument("--rmax", type=float, required=True, help="Outer radius of shell [Mpc/h] (<=200)")
    ap.add_argument("--smooth", type=float, default=4.0, help="Total Gaussian smoothing sigma [Mpc/h] (>=4)")
    args = ap.parse_args()

    out = shell_mean_delta_g_star(args.rmin, args.rmax, args.smooth)

    with OUTPUT_JSON.open("w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    print(f"Wrote: {OUTPUT_JSON}")
    print(json.dumps(out, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

