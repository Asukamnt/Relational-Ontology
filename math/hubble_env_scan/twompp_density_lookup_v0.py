"""
2M++ density-field lookup helper (v0)
====================================

Purpose
-------
Provide a reproducible way to assign an *observer-based* local overdensity proxy
from a public reconstructed density field (2M++).

This is an optional helper for the "observer-based delta" track:
- The 2M++ product provides a luminosity-weighted galaxy density contrast delta_g*
  on a 257^3 Cartesian grid in Galactic coordinates, smoothed with a Gaussian of
  scale 4 Mpc/h (Carrick et al. 2015).
- Under a bias assumption, one may relate delta_m ~ delta_g / b, but b is model-dependent.

Data source
-----------
Download page: https://cosmicflows.iap.fr/download/
Direct numpy file: https://cosmicflows.iap.fr/assets/data/twompp_density.npy

Usage
-----
1) Lookup by Galactic coordinates (l,b,r):
   python twompp_density_lookup_v0.py --l 58.1 --b 88.0 --r 70

2) Lookup by Equatorial coordinates (ra,dec,r) (J2000):
   python twompp_density_lookup_v0.py --ra 194.95 --dec 27.98 --r 70

Notes
-----
- Distances are in comoving Mpc/h (h=H0/100). The 2M++ grid spans [-200,200] Mpc/h.
- This script downloads a relatively large (~O(100MB)) numpy file into a local cache folder.
  The cache is under: math/hubble_env_scan/_cache/
"""

from __future__ import annotations

import argparse
import math
import urllib.request
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np


TWOMPP_DENSITY_URL = "https://cosmicflows.iap.fr/assets/data/twompp_density.npy"

GRID_N = 257
GRID_HALF_EXTENT_MPC_H = 200.0
GRID_SPACING_MPC_H = 400.0 / 256.0  # 1.5625 Mpc/h
GRID_CENTER_INDEX = 128
BASE_SMOOTHING_GAUSSIAN_MPC_H = 4.0

# In-memory cache: key is total Gaussian smoothing (Mpc/h)
_DENSITY_CACHE: Dict[float, np.ndarray] = {}


def _cache_dir() -> Path:
    cache_dir = Path(__file__).parent / "_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def _base_density_path() -> Path:
    return _cache_dir() / "twompp_density.npy"


def _smoothed_density_path(total_sigma_mpc_h: float) -> Path:
    # Keep filename stable and filesystem-friendly.
    tag = f"{total_sigma_mpc_h:.3f}".replace(".", "p")
    return _cache_dir() / f"twompp_density_smoothed_sigma_{tag}_mpch.npy"


def _download_if_missing(url: str, path: Path) -> None:
    if path.exists():
        return
    print(f"Downloading 2M++ density field to: {path}")
    print(f"Source: {url}")
    urllib.request.urlretrieve(url, path)  # noqa: S310 (intentional external download)


def _load_density_field(total_sigma_mpc_h: float = BASE_SMOOTHING_GAUSSIAN_MPC_H) -> np.ndarray:
    """Load (and optionally further smooth) the 2M++ density field.

    The published field is already smoothed with a Gaussian of sigma=4 Mpc/h.
    If total_sigma_mpc_h > 4, we apply an additional Gaussian smoothing so that the
    resulting field has (approximately) total sigma=total_sigma_mpc_h.

    Smoothing is done via FFT with periodic boundary conditions (a limitation).
    """
    total_sigma_mpc_h = float(total_sigma_mpc_h)
    if total_sigma_mpc_h <= 0:
        raise ValueError("total_sigma_mpc_h must be positive")

    # Quantize cache key to avoid float-key surprises.
    cache_key = round(total_sigma_mpc_h, 6)
    if cache_key in _DENSITY_CACHE:
        return _DENSITY_CACHE[cache_key]

    base_path = _base_density_path()
    _download_if_missing(TWOMPP_DENSITY_URL, base_path)

    density = np.load(base_path)
    if density.shape[:3] != (GRID_N, GRID_N, GRID_N):
        raise ValueError(f"Unexpected density shape {density.shape}; expected {(GRID_N, GRID_N, GRID_N)}")

    if abs(total_sigma_mpc_h - BASE_SMOOTHING_GAUSSIAN_MPC_H) < 1e-9:
        _DENSITY_CACHE[cache_key] = density
        return density

    if total_sigma_mpc_h < BASE_SMOOTHING_GAUSSIAN_MPC_H:
        raise ValueError(
            f"Requested sigma={total_sigma_mpc_h} Mpc/h < base smoothing {BASE_SMOOTHING_GAUSSIAN_MPC_H} Mpc/h"
        )

    add_sigma = math.sqrt(max(total_sigma_mpc_h * total_sigma_mpc_h - BASE_SMOOTHING_GAUSSIAN_MPC_H**2, 0.0))

    smoothed_path = _smoothed_density_path(total_sigma_mpc_h)
    if smoothed_path.exists():
        smoothed = np.load(smoothed_path)
        if smoothed.shape[:3] != (GRID_N, GRID_N, GRID_N):
            raise ValueError(f"Unexpected smoothed shape {smoothed.shape}; expected {(GRID_N, GRID_N, GRID_N)}")
        _DENSITY_CACHE[cache_key] = smoothed
        return smoothed

    print(
        f"Computing additional Gaussian smoothing: base={BASE_SMOOTHING_GAUSSIAN_MPC_H} -> total={total_sigma_mpc_h} Mpc/h "
        f"(additional sigma≈{add_sigma:.3f} Mpc/h)"
    )

    # Build k-grid (rad / (Mpc/h))
    freqs = np.fft.fftfreq(GRID_N, d=GRID_SPACING_MPC_H).astype(np.float32)
    ks = (2.0 * np.pi * freqs).astype(np.float32)
    kx, ky, kz = np.meshgrid(ks, ks, ks, indexing="ij")
    k2 = kx * kx + ky * ky + kz * kz
    kernel = np.exp(-0.5 * k2 * (add_sigma * add_sigma)).astype(np.float32)

    # FFT smooth
    ft = np.fft.fftn(density)
    smoothed = np.fft.ifftn(ft * kernel).real.astype(density.dtype, copy=False)

    np.save(smoothed_path, smoothed)
    _DENSITY_CACHE[cache_key] = smoothed
    return smoothed


def _radec_to_galactic_deg(ra_deg: float, dec_deg: float) -> Tuple[float, float]:
    """Convert Equatorial (J2000) RA/Dec in degrees to Galactic l/b in degrees.

    Uses the standard J2000 rotation matrix (same as astropy).
    """
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)

    # Unit vector in equatorial Cartesian coordinates
    x_eq = math.cos(dec) * math.cos(ra)
    y_eq = math.cos(dec) * math.sin(ra)
    z_eq = math.sin(dec)

    # Equatorial -> Galactic rotation matrix (J2000)
    r11, r12, r13 = -0.0548755604, -0.8734370902, -0.4838350155
    r21, r22, r23 = 0.4941094279, -0.4448296300, 0.7469822445
    r31, r32, r33 = -0.8676661490, -0.1980763734, 0.4559837762

    x_g = r11 * x_eq + r12 * y_eq + r13 * z_eq
    y_g = r21 * x_eq + r22 * y_eq + r23 * z_eq
    z_g = r31 * x_eq + r32 * y_eq + r33 * z_eq

    l = math.degrees(math.atan2(y_g, x_g)) % 360.0
    b = math.degrees(math.asin(max(-1.0, min(1.0, z_g))))
    return l, b


def _galactic_spherical_to_cartesian_mpc_h(l_deg: float, b_deg: float, r_mpc_h: float) -> Tuple[float, float, float]:
    l = math.radians(l_deg)
    b = math.radians(b_deg)
    x = r_mpc_h * math.cos(b) * math.cos(l)
    y = r_mpc_h * math.cos(b) * math.sin(l)
    z = r_mpc_h * math.sin(b)
    return x, y, z


def _xyz_to_grid_index(x_mpc_h: float) -> int:
    # Inverse of: X = (i-128)*400/256
    i_float = (x_mpc_h / GRID_SPACING_MPC_H) + GRID_CENTER_INDEX
    i = int(round(i_float))
    return max(0, min(GRID_N - 1, i))


def lookup_delta_g_star(
    l_deg: float,
    b_deg: float,
    r_mpc_h: float,
    *,
    smoothing_gaussian_mpc_h: float = BASE_SMOOTHING_GAUSSIAN_MPC_H,
) -> dict:
    """Nearest-cell lookup of delta_g* in the 2M++ density field.

    Parameters
    ----------
    smoothing_gaussian_mpc_h:
        Total Gaussian smoothing sigma (Mpc/h). The published product is sigma=4 Mpc/h.
        If larger values are requested, we apply an additional Gaussian smoothing (FFT).
    """
    if abs(r_mpc_h) > GRID_HALF_EXTENT_MPC_H:
        raise ValueError(f"r={r_mpc_h} Mpc/h outside 2M++ grid extent ±{GRID_HALF_EXTENT_MPC_H} Mpc/h")

    x, y, z = _galactic_spherical_to_cartesian_mpc_h(l_deg, b_deg, r_mpc_h)
    if any(abs(v) > GRID_HALF_EXTENT_MPC_H for v in (x, y, z)):
        raise ValueError("Position outside 2M++ cube (±200 Mpc/h).")

    i = _xyz_to_grid_index(x)
    j = _xyz_to_grid_index(y)
    k = _xyz_to_grid_index(z)

    density = _load_density_field(total_sigma_mpc_h=smoothing_gaussian_mpc_h)

    delta_g = float(density[i, j, k])

    return {
        "l_deg": float(l_deg),
        "b_deg": float(b_deg),
        "r_mpc_h": float(r_mpc_h),
        "x_mpc_h": float(x),
        "y_mpc_h": float(y),
        "z_mpc_h": float(z),
        "grid_index": {"i": i, "j": j, "k": k},
        "delta_g_star": delta_g,
        "smoothing_gaussian_mpc_h": float(smoothing_gaussian_mpc_h),
        "reference": "Carrick et al. 2015 (2M++ density field; https://cosmicflows.iap.fr/download/)",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--l", type=float, default=None, help="Galactic longitude l [deg]")
    ap.add_argument("--b", type=float, default=None, help="Galactic latitude b [deg]")
    ap.add_argument("--ra", type=float, default=None, help="Equatorial RA (J2000) [deg]")
    ap.add_argument("--dec", type=float, default=None, help="Equatorial Dec (J2000) [deg]")
    ap.add_argument("--r", type=float, required=True, help="Comoving distance r [Mpc/h]")
    ap.add_argument("--smooth", type=float, default=BASE_SMOOTHING_GAUSSIAN_MPC_H, help="Total Gaussian smoothing sigma [Mpc/h] (>=4)")
    ap.add_argument("--bias", type=float, default=None, help="Optional galaxy bias b to compute delta_m~delta_g/b")

    args = ap.parse_args()

    if args.l is None or args.b is None:
        if args.ra is None or args.dec is None:
            raise SystemExit("Provide either (--l,--b) or (--ra,--dec).")
        l_deg, b_deg = _radec_to_galactic_deg(args.ra, args.dec)
    else:
        l_deg, b_deg = args.l, args.b

    out = lookup_delta_g_star(l_deg=l_deg, b_deg=b_deg, r_mpc_h=args.r, smoothing_gaussian_mpc_h=args.smooth)

    if args.bias is not None and args.bias != 0:
        out["bias_b"] = float(args.bias)
        out["delta_m_assuming_bias"] = float(out["delta_g_star"] / args.bias)

    import json

    print(json.dumps(out, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

