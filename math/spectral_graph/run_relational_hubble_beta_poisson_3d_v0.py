"""
Relational Hubble beta from Poisson + entropic flow (toy) v0
===========================================================

Goal
----
Provide a *first-principles within-repo* derivation chain for the sign and
order-of-magnitude of the environmental Hubble bias:

  H_obs(delta) = H_global * (1 - beta * delta)

This script stitches together three primitives that already exist in this repo:

1) Geometry operator: graph Laplacian on a 3D torus grid
   (same construction style as run_alpha_from_coulomb_3d_v0.py).

2) "Potential" from relations: Poisson / Green operator on the Laplacian
     L phi = - delta_rel
   In the repo language, you can view phi as the *linear-response* proxy of the
   diffusion-entropy potential S = -log K_ii(t) (v0 coarse-grained limit).

3) Entropic flow: v ~ -grad(phi), and H0 regression bias from a shell:
     deltaH/H = sum w(r) * r * v_r / sum w(r) * r^2
   This matches the estimator used in hubble_env_scan velocity-bias scripts.

Notes / Scope
-------------
This is a toy. It demonstrates:
- monotonic sign: overdense -> deltaH/H < 0, void -> deltaH/H > 0
- near-linear scaling: deltaH/H ~ -beta * delta_eff

It does NOT claim to fix cosmological normalization (growth factor f, units).
Those can be attached later by relating the entropic-flow time scale to the
emergent causal time (dark_energy_time_arrow_toy) or other dynamics.

Run (repo root)
--------------
python math/spectral_graph/run_relational_hubble_beta_poisson_3d_v0.py

Outputs
-------
- relational_hubble_beta_poisson_3d_results_v0.json
- relational_hubble_beta_poisson_3d_plot_v0.png
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu


OUT_JSON = Path(__file__).parent / "relational_hubble_beta_poisson_3d_results_v0.json"
OUT_PLOT = Path(__file__).parent / "relational_hubble_beta_poisson_3d_plot_v0.png"


def _json_sanitize(obj: Any):
    if isinstance(obj, (float, np.floating)):
        v = float(obj)
        return v if math.isfinite(v) else None
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, dict):
        return {str(k): _json_sanitize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_sanitize(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _json_sanitize(obj.tolist())
    return obj


def _idx(x: int, y: int, z: int, n: int) -> int:
    return int(x) * int(n) * int(n) + int(y) * int(n) + int(z)


def _build_laplacian_torus_3d(n: int, w: float) -> csr_matrix:
    n = int(n)
    if n < 6:
        raise ValueError("n must be >= 6.")

    N = n * n * n
    rows: List[int] = []
    cols: List[int] = []
    data: List[float] = []

    w = float(w)
    for x in range(n):
        for y in range(n):
            for z in range(n):
                i = _idx(x, y, z, n)
                for dx, dy, dz in (
                    (1, 0, 0),
                    (-1, 0, 0),
                    (0, 1, 0),
                    (0, -1, 0),
                    (0, 0, 1),
                    (0, 0, -1),
                ):
                    nx = (x + dx) % n
                    ny = (y + dy) % n
                    nz = (z + dz) % n
                    j = _idx(nx, ny, nz, n)
                    rows.append(i)
                    cols.append(j)
                    data.append(w)

    A = sparse.coo_matrix((data, (rows, cols)), shape=(N, N), dtype=float).tocsr()
    deg = np.asarray(A.sum(axis=1)).ravel()
    L = sparse.diags(deg, offsets=0, format="csr") - A
    return L


def _factor_poisson(L: csr_matrix, ref_index: int) -> Tuple[splu, int]:
    N = int(L.shape[0])
    ref = int(ref_index)
    if ref < 0 or ref >= N:
        raise ValueError("ref_index out of range.")
    A = L.tolil()
    # Gauge fix: enforce phi(ref)=0 by replacing the row with identity.
    A.rows[ref] = [ref]
    A.data[ref] = [1.0]
    lu = splu(A.tocsc())
    return lu, ref


def _solve_phi(lu: splu, ref: int, rhs: np.ndarray) -> np.ndarray:
    b = np.asarray(rhs, dtype=float).copy()
    b[int(ref)] = 0.0
    phi = lu.solve(b)
    return np.asarray(phi, dtype=float)


def _signed_min_image(d: np.ndarray, n: int) -> np.ndarray:
    """
    Signed minimal image displacement for periodic box:
      d in Z (can be negative), maps to [-n/2, n/2] with sign.
    """
    n = int(n)
    dd = np.asarray(d, dtype=int)
    # bring to [-n, n)
    dd = ((dd + n) % (2 * n)) - n
    # now choose minimal magnitude representative
    # If n is even, prefer the negative direction for the exactly-half case (rare; doesn't matter much).
    half = n // 2
    dd = np.where(dd > half, dd - 2 * n, dd)
    dd = np.where(dd < -half, dd + 2 * n, dd)
    return dd.astype(float)


def make_compensated_gaussian_delta(
    n: int,
    *,
    delta0: float,
    sigma: float,
    center: Tuple[int, int, int],
) -> np.ndarray:
    """
    delta(x) = -delta0 * exp(-r^2/(2*sigma^2))  minus its global mean (compensation),
    so sum delta ≈ 0 (needed for physical Poisson on a torus).
    """
    n = int(n)
    cx, cy, cz = (int(center[0]), int(center[1]), int(center[2]))
    xs = np.arange(n, dtype=int)
    ys = np.arange(n, dtype=int)
    zs = np.arange(n, dtype=int)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    dx = _signed_min_image(X - cx, n)
    dy = _signed_min_image(Y - cy, n)
    dz = _signed_min_image(Z - cz, n)
    r2 = dx * dx + dy * dy + dz * dz
    sig2 = float(sigma) ** 2
    if sig2 <= 0:
        raise ValueError("sigma must be > 0")
    core = -float(delta0) * np.exp(-0.5 * r2 / sig2)
    core = core - float(np.mean(core))
    return core.reshape(-1)


def make_compensated_tophat_delta(
    n: int,
    *,
    delta_in: float,
    R: float,
    center: Tuple[int, int, int],
) -> np.ndarray:
    """
    Spherical top-hat contrast with exact zero mean on the discrete torus:
      delta = delta_in   for r <= R
      delta = delta_out  for r >  R
    where delta_out is chosen so that mean(delta)=0.
    """
    n = int(n)
    if float(R) <= 0:
        raise ValueError("R must be > 0")

    cx, cy, cz = (int(center[0]), int(center[1]), int(center[2]))
    xs = np.arange(n, dtype=int)
    ys = np.arange(n, dtype=int)
    zs = np.arange(n, dtype=int)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    dx = _signed_min_image(X - cx, n)
    dy = _signed_min_image(Y - cy, n)
    dz = _signed_min_image(Z - cz, n)
    r = np.sqrt(dx * dx + dy * dy + dz * dz)

    m_in = r <= float(R)
    n_tot = int(n * n * n)
    n_in = int(np.sum(m_in))
    n_out = int(n_tot - n_in)
    if n_in <= 0 or n_out <= 0:
        raise ValueError("Invalid R: sphere is empty or fills the whole box.")

    delta_in = float(delta_in)
    delta_out = -delta_in * float(n_in) / float(n_out)
    delta = np.where(m_in, delta_in, delta_out).astype(float)
    # Numerical check (debug-level): mean should be ~0.
    return delta.reshape(-1)

def grad_periodic_3d(phi_3d: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Central-difference gradient on periodic cubic grid, lattice spacing=1."""
    phi = np.asarray(phi_3d, dtype=float)
    dpx = (np.roll(phi, -1, axis=0) - np.roll(phi, 1, axis=0)) * 0.5
    dpy = (np.roll(phi, -1, axis=1) - np.roll(phi, 1, axis=1)) * 0.5
    dpz = (np.roll(phi, -1, axis=2) - np.roll(phi, 1, axis=2)) * 0.5
    return dpx, dpy, dpz


def shell_weight_uniform(r: np.ndarray, rmin: float, rmax: float) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    w = np.zeros_like(r, dtype=float)
    m = (r >= float(rmin)) & (r <= float(rmax)) & np.isfinite(r)
    w[m] = 1.0
    return w


def compute_observer_shell_stats(
    *,
    n: int,
    delta_flat: np.ndarray,
    v_flat_xyz: Tuple[np.ndarray, np.ndarray, np.ndarray],
    center: Tuple[int, int, int],
    rmin: float,
    rmax: float,
) -> Dict[str, float]:
    """
    Compute:
      - delta_eff: uniform-shell average of delta
      - deltaH_over_H: weighted regression bias estimator using v_r
          deltaH/H = sum w r v_r / sum w r^2   (baseline H_true=1 in lattice units)
    """
    n = int(n)
    cx, cy, cz = (int(center[0]), int(center[1]), int(center[2]))

    xs = np.arange(n, dtype=int)
    ys = np.arange(n, dtype=int)
    zs = np.arange(n, dtype=int)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    dx = _signed_min_image(X - cx, n)
    dy = _signed_min_image(Y - cy, n)
    dz = _signed_min_image(Z - cz, n)
    r = np.sqrt(dx * dx + dy * dy + dz * dz)
    r_flat = r.reshape(-1)

    # Unit vectors (avoid r=0)
    eps = 1e-12
    invr = 1.0 / np.maximum(r_flat, eps)
    rx = (dx.reshape(-1) * invr).astype(float)
    ry = (dy.reshape(-1) * invr).astype(float)
    rz = (dz.reshape(-1) * invr).astype(float)

    vx, vy, vz = v_flat_xyz
    vx = np.asarray(vx, dtype=float).reshape(-1)
    vy = np.asarray(vy, dtype=float).reshape(-1)
    vz = np.asarray(vz, dtype=float).reshape(-1)
    vr = vx * rx + vy * ry + vz * rz

    w = shell_weight_uniform(r_flat, rmin=float(rmin), rmax=float(rmax))
    # Exclude r=0 from all sums
    w = np.where(r_flat < eps, 0.0, w)

    sumw = float(np.sum(w))
    if sumw <= 0:
        raise ValueError("empty shell: choose a wider rmin/rmax.")

    delta = np.asarray(delta_flat, dtype=float).reshape(-1)
    delta_eff = float(np.sum(w * delta) / sumw)

    S_rr = float(np.sum(w * (r_flat**2)))
    S_rvr = float(np.sum(w * r_flat * vr))
    if S_rr <= 0:
        raise ValueError("S_rr <= 0 (unexpected).")

    deltaH_over_H = float(S_rvr / S_rr)  # baseline H_true=1

    beta_eff = None
    if abs(delta_eff) > 1e-12:
        beta_eff = float(-deltaH_over_H / delta_eff)

    return {
        "delta_eff": delta_eff,
        "deltaH_over_H": deltaH_over_H,
        "beta_eff": beta_eff if beta_eff is not None else float("nan"),
        "sumw": sumw,
        "S_rr": S_rr,
        "S_rvr": S_rvr,
    }


@dataclass(frozen=True)
class Config:
    n: int = 18
    w_edge: float = 1.0
    ref_index: int = 0

    # Density profile on the 3D torus (mean-subtracted / compensated)
    profile: str = "tophat"  # "tophat" or "gaussian"
    sigma: float = 3.0  # gaussian width (if profile=gaussian)
    R_void: float = 8.0  # top-hat radius (if profile=tophat)
    delta0_list: Tuple[float, ...] = (0.10, 0.20, 0.30, 0.40, 0.50)

    # Shell for the H0 regression proxy (lattice units)
    rmin: float = 3.0
    rmax: float = 7.0

    # Optional velocity scale factor (dimensionless). Keep default 1 for "pure" toy.
    v_scale: float = 1.0

    seed: int = 0  # reserved (kept for schema stability)


def fit_beta(delta_eff: np.ndarray, deltaH_over_H: np.ndarray) -> Dict[str, float]:
    x = np.asarray(delta_eff, dtype=float).reshape(-1)
    y = np.asarray(deltaH_over_H, dtype=float).reshape(-1)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if x.size < 2:
        return {"beta": float("nan"), "intercept": float("nan"), "r2": float("nan"), "pearson_r": float("nan")}

    # Fit y = a*x + b ; expect a ~ -beta
    A = np.stack([x, np.ones_like(x)], axis=1)
    coef, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    a = float(coef[0])
    b = float(coef[1])
    y_hat = a * x + b
    ss_res = float(np.sum((y - y_hat) ** 2))
    ss_tot = float(np.sum((y - float(np.mean(y))) ** 2)) + 1e-18
    r2 = float(1.0 - ss_res / ss_tot)

    # Pearson r
    xr = x - float(np.mean(x))
    yr = y - float(np.mean(y))
    denom = float(np.sqrt(np.sum(xr * xr) * np.sum(yr * yr))) + 1e-18
    pr = float(np.sum(xr * yr) / denom)

    return {"beta": float(-a), "intercept": b, "r2": r2, "pearson_r": pr}


def main() -> None:
    p = argparse.ArgumentParser(description="Toy derivation: beta from Poisson + entropic flow on 3D torus (v0).")
    p.add_argument("--n", type=int, default=Config.n)
    p.add_argument("--w-edge", type=float, default=Config.w_edge)
    p.add_argument("--profile", type=str, default=Config.profile, choices=["tophat", "gaussian"])
    p.add_argument("--sigma", type=float, default=Config.sigma)
    p.add_argument("--R-void", type=float, default=Config.R_void, help="top-hat void radius (lattice units)")
    p.add_argument("--delta0-list", type=str, default="0.10,0.20,0.30,0.40,0.50")
    p.add_argument("--rmin", type=float, default=Config.rmin)
    p.add_argument("--rmax", type=float, default=Config.rmax)
    p.add_argument("--v-scale", type=float, default=Config.v_scale)
    args = p.parse_args()

    n = int(args.n)
    if n < 6:
        raise ValueError("--n must be >= 6")
    if float(args.w_edge) <= 0:
        raise ValueError("--w-edge must be > 0")
    if str(args.profile).strip().lower() == "gaussian" and float(args.sigma) <= 0:
        raise ValueError("--sigma must be > 0 for profile=gaussian")
    if str(args.profile).strip().lower() == "tophat" and float(args.R_void) <= 0:
        raise ValueError("--R-void must be > 0 for profile=tophat")
    if float(args.rmax) <= float(args.rmin):
        raise ValueError("--rmax must be > --rmin")

    delta0_list = [float(x) for x in str(args.delta0_list).split(",") if str(x).strip()]
    if not delta0_list:
        raise ValueError("--delta0-list must be non-empty")

    cfg = Config(
        n=int(n),
        w_edge=float(args.w_edge),
        profile=str(args.profile).strip().lower(),
        sigma=float(args.sigma),
        R_void=float(args.R_void),
        delta0_list=tuple(delta0_list),
        rmin=float(args.rmin),
        rmax=float(args.rmax),
        v_scale=float(args.v_scale),
    )

    # Observer at box center
    center = (cfg.n // 2, cfg.n // 2, cfg.n // 2)

    print("=== relational hubble beta poisson 3d toy (v0) ===")
    print(
        f"n={cfg.n}, profile={cfg.profile}, sigma={cfg.sigma}, R_void={cfg.R_void}, "
        f"shell=[{cfg.rmin},{cfg.rmax}], v_scale={cfg.v_scale}"
    )

    # Build Laplacian and factor once.
    L = _build_laplacian_torus_3d(cfg.n, float(cfg.w_edge))
    lu, ref = _factor_poisson(L, int(cfg.ref_index))

    rows: List[Dict[str, float]] = []

    for delta0 in cfg.delta0_list:
        # Use a compensated "void" with amplitude delta0 (delta is negative at center).
        if str(cfg.profile) == "gaussian":
            delta = make_compensated_gaussian_delta(cfg.n, delta0=float(delta0), sigma=float(cfg.sigma), center=center)
        else:
            delta = make_compensated_tophat_delta(cfg.n, delta_in=-float(delta0), R=float(cfg.R_void), center=center)
        # Graph Laplacian L = D - A corresponds to -∇^2 in the continuum limit.
        # For a gravity-like sign convention (overdensity attracts / void repels),
        # we solve:  L phi = - delta   (so that ∇^2 phi = delta).
        phi = _solve_phi(lu, ref, -delta)
        phi_3d = phi.reshape((cfg.n, cfg.n, cfg.n))

        dpx, dpy, dpz = grad_periodic_3d(phi_3d)
        # Entropic flow (toy): v = -grad(phi)
        vx = (-float(cfg.v_scale) * dpx).reshape(-1)
        vy = (-float(cfg.v_scale) * dpy).reshape(-1)
        vz = (-float(cfg.v_scale) * dpz).reshape(-1)

        stats = compute_observer_shell_stats(
            n=cfg.n,
            delta_flat=delta,
            v_flat_xyz=(vx, vy, vz),
            center=center,
            rmin=float(cfg.rmin),
            rmax=float(cfg.rmax),
        )

        rows.append(
            {
                "delta0": float(delta0),
                "delta_eff": float(stats["delta_eff"]),
                "deltaH_over_H": float(stats["deltaH_over_H"]),
                "beta_eff": float(stats["beta_eff"]),
            }
        )

        print(
            f"delta0={float(delta0):.3f}  delta_eff={float(stats['delta_eff']): .4f}  "
            f"deltaH/H={float(stats['deltaH_over_H']): .4f}  beta_eff={float(stats['beta_eff']): .4f}"
        )

    delta_eff = np.array([r["delta_eff"] for r in rows], dtype=float)
    dH = np.array([r["deltaH_over_H"] for r in rows], dtype=float)
    fit = fit_beta(delta_eff, dH)

    out: Dict[str, Any] = {
        "metadata": {"timestamp": datetime.now().isoformat(), "script": str(Path(__file__).name)},
        "config": _json_sanitize(asdict(cfg)),
        "model": {
            "definitions": {
                    "delta_rel": "node scalar field on 3D torus (profile=tophat or gaussian), constructed to have mean~0",
                    "potential_phi": "Poisson solve on graph Laplacian: L phi = -delta_rel (gauge fixed by phi(ref)=0); graph L≈-∇^2 so ∇^2 phi = delta_rel",
                "velocity": "entropic flow proxy v = -v_scale * grad(phi) (central difference, periodic)",
                "hubble_bias_estimator": "deltaH/H = sum w r v_r / sum w r^2 with uniform shell weight",
            },
            "notes": [
                "This toy focuses on sign + near-linearity, not cosmological normalization.",
                "To map to beta=f/3 in standard linear theory, attach a dynamics-derived scaling between v and delta.",
            ],
        },
        "series": rows,
        "fit_deltaH_over_H_vs_delta_eff": fit,
        "outputs": {"json": str(OUT_JSON), "plot": str(OUT_PLOT)},
    }

    with OUT_JSON.open("w", encoding="utf-8") as f:
        json.dump(_json_sanitize(out), f, indent=2, ensure_ascii=False)

    # Plot: deltaH/H vs delta_eff + fit line.
    fig, ax = plt.subplots(1, 1, figsize=(7.6, 5.2))
    ax.axhline(0.0, color="#888888", linewidth=1.0, alpha=0.6)
    ax.axvline(0.0, color="#888888", linewidth=1.0, alpha=0.6)
    ax.scatter(delta_eff, dH, s=50, color="#1f77b4", label="toy points")

    if np.isfinite(fit.get("beta", float("nan"))):
        xs = np.linspace(float(np.min(delta_eff)), float(np.max(delta_eff)), 200)
        ys = -(float(fit["beta"])) * xs + float(fit["intercept"])
        ax.plot(xs, ys, color="#d62728", linewidth=2.0, label=f"fit: beta={float(fit['beta']):.3f}, r={float(fit['pearson_r']):.3f}")

    ax.set_xlabel("delta_eff (shell mean of delta_rel)")
    ax.set_ylabel("deltaH/H (proxy from v_r regression)")
    ax.set_title("Relational Hubble bias toy: Poisson + entropic flow (3D torus)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT_PLOT, dpi=150)
    plt.close(fig)

    print(f"Output JSON: {OUT_JSON}")
    print(f"Output plot: {OUT_PLOT}")
    print("DONE.")


if __name__ == "__main__":
    main()

