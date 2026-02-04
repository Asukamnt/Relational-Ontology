"""
2M++ overdense observer candidates scan (v0)
============================================

Motivation
----------
Our Hubble-tension model is primarily about **observer environment**. However, we currently lack
clean "direct" overdense observer points (the desired “4th bullet”).

This script produces a *proxy* shortlist of **strongly overdense candidate observer locations**
within the 2M++ cube, to guide follow-up searches for independent H0 measurements in overdense regions.

Method
------
Using the public 2M++ reconstructed density field (Carrick et al. 2015), we:

1) Load (optionally smoothed) delta_g* field.
2) Find the top-K grid cells by delta_g* (after enforcing a minimum separation).
3) For each candidate center, compute "observer-shell" averaged delta_g* in selected shells,
   approximating the environment an observer placed there would have.

Outputs
-------
- twompp_overdense_observer_candidates_results_v0.json
- twompp_overdense_observer_candidates_plot_v0.png

Run (repo root)
---------------
python math/hubble_env_scan/run_twompp_overdense_observer_candidates_v0.py
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

import twompp_density_lookup_v0 as twompp


OUTPUT_JSON = Path(__file__).parent / "twompp_overdense_observer_candidates_results_v0.json"
OUTPUT_PLOT = Path(__file__).parent / "twompp_overdense_observer_candidates_plot_v0.png"


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


def _grid_to_xyz(i: int, j: int, k: int) -> Tuple[float, float, float]:
    # twompp._xyz_to_grid_index uses: i_float = x/spacing + center
    # so inverse is: x = (i - center) * spacing
    x = (float(i) - float(twompp.GRID_CENTER_INDEX)) * float(twompp.GRID_SPACING_MPC_H)
    y = (float(j) - float(twompp.GRID_CENTER_INDEX)) * float(twompp.GRID_SPACING_MPC_H)
    z = (float(k) - float(twompp.GRID_CENTER_INDEX)) * float(twompp.GRID_SPACING_MPC_H)
    return x, y, z


def _xyz_to_grid(x: float, y: float, z: float) -> Tuple[int, int, int]:
    i = twompp._xyz_to_grid_index(float(x))
    j = twompp._xyz_to_grid_index(float(y))
    k = twompp._xyz_to_grid_index(float(z))
    return int(i), int(j), int(k)


def _sphere_mask_offsets(rmin: float, rmax: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Offsets (di,dj,dk) for grid points within shell [rmin,rmax] (in Mpc/h)."""
    rmin = float(rmin)
    rmax = float(rmax)
    if rmax <= 0 or rmax < rmin:
        raise ValueError("Invalid shell radii")
    step = float(twompp.GRID_SPACING_MPC_H)
    n = int(math.ceil(rmax / step))
    grid = np.arange(-n, n + 1, dtype=np.int16)
    di, dj, dk = np.meshgrid(grid, grid, grid, indexing="ij")
    rr = np.sqrt(((di.astype(np.float32) * step) ** 2) + ((dj.astype(np.float32) * step) ** 2) + ((dk.astype(np.float32) * step) ** 2))
    m = (rr >= rmin) & (rr <= rmax)
    return di[m].astype(np.int16), dj[m].astype(np.int16), dk[m].astype(np.int16)


@dataclass
class Candidate:
    rank: int
    i: int
    j: int
    k: int
    x_mpc_h: float
    y_mpc_h: float
    z_mpc_h: float
    delta_g_star_center: float
    shells: Dict[str, Any]


def _within_cube(i: int, j: int, k: int) -> bool:
    n = int(twompp.GRID_N)
    return 0 <= i < n and 0 <= j < n and 0 <= k < n


def _shell_mean(
    density: np.ndarray,
    *,
    center_ijk: Tuple[int, int, int],
    di: np.ndarray,
    dj: np.ndarray,
    dk: np.ndarray,
) -> Tuple[Optional[float], Optional[int]]:
    ci, cj, ck = center_ijk
    ii = ci + di
    jj = cj + dj
    kk = ck + dk
    m = (ii >= 0) & (ii < density.shape[0]) & (jj >= 0) & (jj < density.shape[1]) & (kk >= 0) & (kk < density.shape[2])
    if not np.any(m):
        return None, None
    vals = density[ii[m], jj[m], kk[m]]
    return float(np.mean(vals)), int(vals.size)


def pick_topk_with_separation(
    density: np.ndarray,
    *,
    k: int,
    min_sep_mpc_h: float,
) -> List[Tuple[int, int, int, float]]:
    """Pick top-K grid cells by density value, enforcing minimum spatial separation."""
    k = int(k)
    if k <= 0:
        return []

    flat = density.reshape(-1)
    # Get many top candidates first, then filter by separation.
    # NOTE: For smoothed fields (e.g. σ=50–100), the very top cells can be concentrated
    # in a small region. We therefore inspect a much larger pool before applying the
    # separation filter.
    n_try = min(flat.size, max(k * 20000, 200000))
    top_idx = np.argpartition(flat, -n_try)[-n_try:]
    top_idx = top_idx[np.argsort(flat[top_idx])[::-1]]

    chosen: List[Tuple[int, int, int, float]] = []
    chosen_xyz: List[Tuple[float, float, float]] = []
    min_sep2 = float(min_sep_mpc_h) ** 2

    for idx in top_idx:
        i = int(idx // (twompp.GRID_N * twompp.GRID_N))
        rem = int(idx % (twompp.GRID_N * twompp.GRID_N))
        j = int(rem // twompp.GRID_N)
        k2 = int(rem % twompp.GRID_N)
        v = float(density[i, j, k2])
        x, y, z = _grid_to_xyz(i, j, k2)

        ok = True
        for (x2, y2, z2) in chosen_xyz:
            if (x - x2) ** 2 + (y - y2) ** 2 + (z - z2) ** 2 < min_sep2:
                ok = False
                break
        if not ok:
            continue
        chosen.append((i, j, k2, v))
        chosen_xyz.append((x, y, z))
        if len(chosen) >= k:
            break

    return chosen


def _xyz_to_lb_r(x: float, y: float, z: float) -> Tuple[float, float, float]:
    r = float(math.sqrt(x * x + y * y + z * z))
    if r <= 0:
        return 0.0, 0.0, 0.0
    l = math.degrees(math.atan2(y, x)) % 360.0
    b = math.degrees(math.asin(max(-1.0, min(1.0, z / r))))
    return float(l), float(b), float(r)


def _galactic_to_radec_deg(l_deg: float, b_deg: float) -> Tuple[float, float]:
    """
    Galactic (l,b) -> Equatorial (RA,Dec) (J2000) using inverse of the rotation in twompp.
    """
    l = math.radians(float(l_deg))
    b = math.radians(float(b_deg))
    x_g = math.cos(b) * math.cos(l)
    y_g = math.cos(b) * math.sin(l)
    z_g = math.sin(b)

    # Equatorial -> Galactic matrix (R); inverse is transpose.
    r11, r12, r13 = -0.0548755604, -0.8734370902, -0.4838350155
    r21, r22, r23 = 0.4941094279, -0.4448296300, 0.7469822445
    r31, r32, r33 = -0.8676661490, -0.1980763734, 0.4559837762

    x_eq = r11 * x_g + r21 * y_g + r31 * z_g
    y_eq = r12 * x_g + r22 * y_g + r32 * z_g
    z_eq = r13 * x_g + r23 * y_g + r33 * z_g

    ra = math.degrees(math.atan2(y_eq, x_eq)) % 360.0
    dec = math.degrees(math.asin(max(-1.0, min(1.0, z_eq))))
    return float(ra), float(dec)


def main() -> None:
    ap = argparse.ArgumentParser(description="2M++ overdense observer candidates scan (v0)")
    ap.add_argument("--smooth", type=float, default=20.0, help="Total Gaussian smoothing sigma for delta_g* [Mpc/h] (>=4)")
    ap.add_argument("--topk", type=int, default=50, help="Number of candidate centers to report")
    ap.add_argument("--min-sep", type=float, default=20.0, help="Minimum separation between candidates [Mpc/h]")
    ap.add_argument(
        "--shells",
        type=str,
        default="0-75,30-75,0-150,40-150",
        help="Comma-separated shell ranges rmin-rmax [Mpc/h]. Note: large rmax can be truncated near cube edges.",
    )
    ap.add_argument("--output-json", type=str, default=str(OUTPUT_JSON))
    ap.add_argument("--output-plot", type=str, default=str(OUTPUT_PLOT))
    args = ap.parse_args()

    smooth = float(args.smooth)
    if smooth < 4.0:
        raise SystemExit("smooth must be >= 4.0 Mpc/h")

    density = twompp._load_density_field(total_sigma_mpc_h=smooth)

    # Parse shell list
    shell_specs: List[Tuple[float, float]] = []
    for part in str(args.shells).split(","):
        t = part.strip()
        if not t:
            continue
        if "-" not in t:
            raise SystemExit(f"Bad shell spec '{t}', expected 'rmin-rmax'")
        a, b = t.split("-", 1)
        rmin = float(a.strip())
        rmax = float(b.strip())
        shell_specs.append((rmin, rmax))
    if not shell_specs:
        shell_specs = [(0.0, 67.0), (0.0, 75.0), (0.0, 150.0), (0.0, 200.0), (40.0, 200.0)]

    # Precompute offset masks per shell for speed
    shell_offsets: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for (rmin, rmax) in shell_specs:
        key = f"{rmin:g}-{rmax:g}"
        shell_offsets[key] = _sphere_mask_offsets(rmin, rmax)

    # Pick candidates
    picked = pick_topk_with_separation(density, k=int(args.topk), min_sep_mpc_h=float(args.min_sep))

    candidates: List[Candidate] = []
    for idx, (i, j, k, v) in enumerate(picked, start=1):
        x, y, z = _grid_to_xyz(i, j, k)
        l_deg, b_deg, r_mpc_h = _xyz_to_lb_r(x, y, z)
        ra_deg, dec_deg = _galactic_to_radec_deg(l_deg, b_deg) if r_mpc_h > 0 else (0.0, 0.0)
        shells_out: Dict[str, Any] = {}
        for key, (di, dj, dk) in shell_offsets.items():
            mean, n = _shell_mean(density, center_ijk=(i, j, k), di=di, dj=dj, dk=dk)
            n_full = int(di.size)
            cov = (float(n) / float(n_full)) if (n is not None and n_full > 0) else None
            shells_out[key] = {"mean_delta_g_star": mean, "n_cells": n, "n_cells_full": n_full, "coverage": cov}
        candidates.append(
            Candidate(
                rank=int(idx),
                i=int(i),
                j=int(j),
                k=int(k),
                x_mpc_h=float(x),
                y_mpc_h=float(y),
                z_mpc_h=float(z),
                delta_g_star_center=float(v),
                shells=shells_out,
            )
        )

    # Summaries for plot
    center_vals = np.array([c.delta_g_star_center for c in candidates], dtype=float) if candidates else np.array([], dtype=float)
    shell_keys = list(shell_offsets.keys())
    shell_means: Dict[str, List[float]] = {k: [] for k in shell_keys}
    for c in candidates:
        for sk in shell_keys:
            shell_means[sk].append(float(c.shells[sk]["mean_delta_g_star"]) if c.shells[sk]["mean_delta_g_star"] is not None else float("nan"))

    # Plot: center delta and a few shell means vs rank
    fig, ax = plt.subplots(figsize=(10.0, 6.0))
    ranks = np.arange(1, len(candidates) + 1, dtype=int)
    if candidates:
        ax.plot(ranks, center_vals, "o-", color="#111111", linewidth=1.2, markersize=4, label="center δ_g*")
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
        for ci, sk in enumerate(shell_keys[:5]):
            y = np.array(shell_means[sk], dtype=float)
            ax.plot(ranks, y, "o--", linewidth=1.0, markersize=3, color=colors[ci % len(colors)], label=f"shell mean δ_g* ({sk})")
    ax.set_xlabel("candidate rank (by center δ_g*)")
    ax.set_ylabel("δ_g*")
    ax.set_title(f"2M++ overdense observer candidates (σ={smooth:g} Mpc/h; min_sep={float(args.min_sep):g} Mpc/h)")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(Path(args.output_plot), dpi=150)
    plt.close(fig)

    out: Dict[str, Any] = {
        "metadata": {"timestamp": datetime.now().isoformat(), "script": Path(__file__).name},
        "inputs": {
            "smoothing_gaussian_mpc_h": float(smooth),
            "topk": int(args.topk),
            "min_sep_mpc_h": float(args.min_sep),
            "shells": [f"{a:g}-{b:g}" for (a, b) in shell_specs],
        },
        "notes": [
            "This is a proxy shortlist using 2M++ delta_g* (galaxy overdensity).",
            "Shell means approximate 'observer environment' around candidate centers.",
            "Periodic-BC FFT smoothing is a limitation at large sigma.",
        ],
        "candidates": [
            {
                "rank": c.rank,
                "grid_index": {"i": c.i, "j": c.j, "k": c.k},
                "xyz_mpc_h": {"x": c.x_mpc_h, "y": c.y_mpc_h, "z": c.z_mpc_h},
                "lbr_deg_mpc_h": (
                    lambda lb: {"l_deg": float(lb[0]), "b_deg": float(lb[1]), "r_mpc_h": float(lb[2])}
                )(lb_r := _xyz_to_lb_r(c.x_mpc_h, c.y_mpc_h, c.z_mpc_h)),
                "radec_deg": (
                    lambda rd: {"ra_deg": float(rd[0]), "dec_deg": float(rd[1])}
                )(_galactic_to_radec_deg(lb_r[0], lb_r[1]) if lb_r[2] > 0 else (0.0, 0.0)),
                "delta_g_star_center": c.delta_g_star_center,
                "shells": c.shells,
            }
            for c in candidates
        ],
        "outputs": {"json": str(Path(args.output_json)), "plot": str(Path(args.output_plot))},
        "reference": "Carrick et al. 2015 (2M++ density field; https://cosmicflows.iap.fr/download/)",
    }

    Path(args.output_json).write_text(json.dumps(_json_sanitize(out), indent=2, ensure_ascii=False), encoding="utf-8")
    print("2M++ overdense observer candidates scan completed.")
    print(f"Candidates: {len(candidates)}")
    print(f"Output JSON: {args.output_json}")
    print(f"Output plot: {args.output_plot}")


if __name__ == "__main__":
    main()

