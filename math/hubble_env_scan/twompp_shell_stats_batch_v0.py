"""
Batch shell statistics for the 2M++ density field (v0)
=====================================================

Reads `twompp_shell_definitions_v0.csv` and produces:
- `twompp_shell_stats_grid_v0.csv`
- `twompp_shell_stats_grid_v0.json`

Each row corresponds to (shell_id, smoothing_sigma) and reports:
- mean/std of delta_g* inside the shell around the Local Group (origin)

Important: delta_g* is a luminosity-weighted galaxy density contrast (Carrick+2015).
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from twompp_density_lookup_v0 import (
    BASE_SMOOTHING_GAUSSIAN_MPC_H,
    GRID_CENTER_INDEX,
    GRID_HALF_EXTENT_MPC_H,
    GRID_N,
    GRID_SPACING_MPC_H,
    _load_density_field,
    lookup_delta_g_star,
)


SHELLS_CSV = Path(__file__).parent / "twompp_shell_definitions_v0.csv"
OUTPUT_CSV = Path(__file__).parent / "twompp_shell_stats_grid_v0.csv"
OUTPUT_JSON = Path(__file__).parent / "twompp_shell_stats_grid_v0.json"


def _parse_float(value: str) -> Optional[float]:
    if value is None:
        return None
    t = value.strip()
    if not t:
        return None
    try:
        return float(t)
    except ValueError:
        return None


def read_shells(path: Path) -> List[Dict]:
    rows: List[Dict] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            shell_id = (r.get("shell_id") or "").strip()
            rmin = _parse_float(r.get("rmin_mpc_h", ""))
            rmax = _parse_float(r.get("rmax_mpc_h", ""))
            if not shell_id or rmin is None or rmax is None:
                continue
            rows.append(
                {
                    "shell_id": shell_id,
                    "rmin_mpc_h": float(rmin),
                    "rmax_mpc_h": float(rmax),
                    "notes": (r.get("notes") or "").strip(),
                }
            )
    return rows


def _axis_coords() -> np.ndarray:
    idx = np.arange(GRID_N, dtype=np.float32)
    return (idx - GRID_CENTER_INDEX) * GRID_SPACING_MPC_H


def _radius_sq_grid() -> np.ndarray:
    xs = _axis_coords()
    # r^2 via broadcasting (float32)
    return (xs[:, None, None] ** 2) + (xs[None, :, None] ** 2) + (xs[None, None, :] ** 2)


def compute_shell_stats(
    density: np.ndarray, r2: np.ndarray, rmin: float, rmax: float
) -> Tuple[float, float, int]:
    rmin2 = float(rmin) * float(rmin)
    rmax2 = float(rmax) * float(rmax)
    mask = (r2 >= rmin2) & (r2 < rmax2)
    vals = density[mask]
    if vals.size == 0:
        return float("nan"), float("nan"), 0
    return float(np.mean(vals)), float(np.std(vals)), int(vals.size)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--smooth",
        type=str,
        default="4,20,50,100",
        help="Comma-separated total Gaussian smoothing sigmas [Mpc/h] (>=4)",
    )
    args = ap.parse_args()

    smooth_list = []
    for part in (args.smooth or "").split(","):
        part = part.strip()
        if not part:
            continue
        smooth_list.append(float(part))

    shells = read_shells(SHELLS_CSV)
    if not shells:
        raise SystemExit(f"No shells found in {SHELLS_CSV}")

    # Quick sanity: shells must fit in cube
    for s in shells:
        if s["rmax_mpc_h"] > GRID_HALF_EXTENT_MPC_H:
            raise SystemExit(f"Shell {s['shell_id']} has rmax>{GRID_HALF_EXTENT_MPC_H} Mpc/h (not supported).")

    r2 = _radius_sq_grid()

    out_rows: List[Dict] = []

    for sigma in smooth_list:
        if sigma < BASE_SMOOTHING_GAUSSIAN_MPC_H:
            raise SystemExit(
                f"Requested smooth={sigma} Mpc/h < base smoothing {BASE_SMOOTHING_GAUSSIAN_MPC_H} Mpc/h"
            )

        density = _load_density_field(total_sigma_mpc_h=float(sigma))
        origin = lookup_delta_g_star(
            l_deg=0.0, b_deg=0.0, r_mpc_h=0.0, smoothing_gaussian_mpc_h=float(sigma)
        )

        for s in shells:
            mean, std, n = compute_shell_stats(
                density=density,
                r2=r2,
                rmin=float(s["rmin_mpc_h"]),
                rmax=float(s["rmax_mpc_h"]),
            )
            out_rows.append(
                {
                    "shell_id": s["shell_id"],
                    "rmin_mpc_h": float(s["rmin_mpc_h"]),
                    "rmax_mpc_h": float(s["rmax_mpc_h"]),
                    "smoothing_gaussian_mpc_h": float(sigma),
                    "delta_g_star_mean": mean,
                    "delta_g_star_std": std,
                    "n_cells": n,
                    "origin_delta_g_star": float(origin["delta_g_star"]),
                    "notes": s["notes"],
                    "reference": "Carrick et al. 2015 (2M++ density field; https://cosmicflows.iap.fr/download/)",
                }
            )

    with OUTPUT_JSON.open("w", encoding="utf-8") as f:
        json.dump(
            {
                "shells_file": str(SHELLS_CSV),
                "grid": {
                    "N": GRID_N,
                    "spacing_mpc_h": GRID_SPACING_MPC_H,
                    "extent_mpc_h": GRID_HALF_EXTENT_MPC_H,
                },
                "rows": out_rows,
            },
            f,
            indent=2,
            ensure_ascii=False,
        )

    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "shell_id",
            "rmin_mpc_h",
            "rmax_mpc_h",
            "smoothing_gaussian_mpc_h",
            "delta_g_star_mean",
            "delta_g_star_std",
            "n_cells",
            "origin_delta_g_star",
            "reference",
            "notes",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in out_rows:
            w.writerow({k: r.get(k) for k in fieldnames})

    print(f"Wrote: {OUTPUT_CSV}")
    print(f"Wrote: {OUTPUT_JSON}")


if __name__ == "__main__":
    main()

