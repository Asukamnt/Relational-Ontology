"""
CF2 cone proxy bins (v0)
========================

Build a small *proxy-track* dataset (delta_g* vs H0) from Cosmicflows-2 (CF2) within
a sky cone, using the public 2M++ reconstructed density field for delta_g*.

Why this exists
---------------
We want to "try Shapley" (or any overdense direction) without hand-curating many
papers. This script:

- selects CF2 objects (galaxy groups by default) in a cone around (l,b)
- computes delta_g* from 2M++ at the object position
- computes H0 ~ V/Dist (optionally with a CF2-provided flow correction PVa)
- bins objects by delta_g* quantiles and outputs (delta_mean, H0_mean) points

Important caveat
----------------
This is *not* the same quantity as "observer environment bias" in the main
direct-delta model. Use these points as exploratory / proxy evidence only.

Outputs
-------
1) bins CSV in hubble_env_data_v1 schema (appendable):
   hubble_env_cf2_cone_{cone_name}_proxy_bins_sigma{smooth}_v0.csv
2) optional merged dataset (base v1 + bins):
   hubble_env_data_v1_cf2_cone_{cone_name}_proxy_sigma{smooth}_v0.csv
3) summary JSON + plot:
   cf2_cone_{cone_name}_proxy_bins_results_sigma{smooth}_v0.json
   cf2_cone_{cone_name}_proxy_bins_plot_sigma{smooth}_v0.png

Run (repo root)
---------------
python math/hubble_env_scan/run_cf2_cone_proxy_bins_v0.py --cone-name shapley --smooth 20
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

import build_cf2_proxy_bins_v0 as cf2
import run_hubble_env_scan_v1 as scan
from twompp_density_lookup_v0 import lookup_delta_g_star


@dataclass(frozen=True)
class ConeConfig:
    name: str
    l_deg: float
    b_deg: float
    radius_deg: float


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


def _angular_sep_deg(l1_deg: float, b1_deg: float, l2_deg: float, b2_deg: float) -> float:
    """Great-circle distance on the sphere (deg) for Galactic coordinates."""
    l1 = math.radians(float(l1_deg))
    b1 = math.radians(float(b1_deg))
    l2 = math.radians(float(l2_deg))
    b2 = math.radians(float(b2_deg))
    cosang = math.sin(b1) * math.sin(b2) + math.cos(b1) * math.cos(b2) * math.cos(l1 - l2)
    cosang = max(-1.0, min(1.0, cosang))
    return math.degrees(math.acos(cosang))


def _fit_line(d: np.ndarray, y: np.ndarray, yerr: np.ndarray) -> Dict[str, Any]:
    """Fit y = a - b*d (weighted), return a,b and derived beta=b/a."""
    d = np.asarray(d, dtype=float)
    y = np.asarray(y, dtype=float)
    yerr = np.asarray(yerr, dtype=float)
    m = np.isfinite(d) & np.isfinite(y) & np.isfinite(yerr) & (yerr > 0)
    d = d[m]
    y = y[m]
    yerr = yerr[m]
    if d.size < 2 or np.unique(d).size < 2:
        return {"error": "insufficient_points", "n_points": int(d.size)}
    a, b, cov, chi2, chi2_dof, corr = scan.fit_ab_weighted(d, y, yerr)
    beta = (b / a) if np.isfinite(a) and a != 0 else float("nan")
    return {
        "H0_global_like_a": float(a),
        "b": float(b),
        "beta_like": float(beta),
        "chi2": chi2,
        "chi2_dof": chi2_dof,
        "corr_delta_H0": corr,
        "cov_ab": cov.tolist() if cov is not None else None,
        "n_points": int(d.size),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Build CF2 cone proxy bins (v0)")
    ap.add_argument("--cone-name", type=str, default="shapley")
    # Default to the Shapley "concentration" center used in twompp_targets_delta_*.csv
    ap.add_argument("--center-l", type=float, default=311.531118363705, help="Cone center l (deg)")
    ap.add_argument("--center-b", type=float, default=32.31073097097159, help="Cone center b (deg)")
    ap.add_argument("--radius-deg", type=float, default=30.0, help="Cone radius (deg)")
    ap.add_argument("--smooth", type=float, default=20.0, help="Total Gaussian smoothing for 2M++ [Mpc/h] (>=4)")
    ap.add_argument("--table", type=str, default="table2", choices=["table1", "table2"])
    ap.add_argument("--out-max", type=int, default=12000, help="VizieR row cap (CF2 has 8315 rows)")
    ap.add_argument("--force-download", action="store_true", help="Redownload CF2 even if cached file exists")

    ap.add_argument("--dist-min", type=float, default=60.0, help="Min CF2 distance [Mpc]")
    ap.add_argument("--dist-max", type=float, default=250.0, help="Max CF2 distance [Mpc]")
    ap.add_argument("--err-max", type=float, default=0.25, help="Max CF2 fractional distance error")
    ap.add_argument("--vmin", type=float, default=6000.0, help="Min velocity [km/s] (after optional flow correction)")
    ap.add_argument("--vmax", type=float, default=20000.0, help="Max velocity [km/s]")
    ap.add_argument("--peculiar-vel-err", type=float, default=250.0, help="Extra velocity error [km/s] for H0_err")
    ap.add_argument(
        "--velocity-mode",
        type=str,
        default="hubble_flow",
        choices=["raw", "hubble_flow"],
        help="raw: use Vcmba. hubble_flow: use (Vcmba - PVa) when available (table2 only).",
    )
    ap.add_argument(
        "--h-ref",
        type=float,
        default=0.6736,
        help="Reference h=H0/100 used to convert CF2 Dist[Mpc] -> r[Mpc/h] for 2M++ lookup.",
    )
    ap.add_argument(
        "--r-from",
        type=str,
        default="dist",
        choices=["dist", "vcmba"],
        help="How to set r (Mpc/h) for 2M++ lookup: dist (Dist*h_ref) or vcmba (V/100).",
    )
    ap.add_argument("--group-base-dist-frac-err", type=float, default=0.20)
    ap.add_argument("--group-min-n1", type=int, default=2)

    ap.add_argument("--n-bins", type=int, default=5, help="Quantile bins in delta_g* (within cone)")
    ap.add_argument("--min-per-bin", type=int, default=15, help="Minimum objects per bin")

    ap.add_argument(
        "--base-data-file",
        type=str,
        default=str(Path(__file__).parent / "hubble_env_data_v1.csv"),
        help="Base v1 dataset to merge with (optional).",
    )
    ap.add_argument("--out-bins-file", type=str, default="")
    ap.add_argument("--out-merged-file", type=str, default="")
    ap.add_argument("--out-results-json", type=str, default="")
    ap.add_argument("--out-plot", type=str, default="")
    args = ap.parse_args()

    cone = ConeConfig(
        name=str(args.cone_name).strip(),
        l_deg=float(args.center_l),
        b_deg=float(args.center_b),
        radius_deg=float(args.radius_deg),
    )
    smooth = float(args.smooth)
    if smooth < 4.0:
        raise ValueError("smooth must be >= 4.0 (2M++ base smoothing)")

    # Resolve outputs
    out_tag = f"{cone.name}_sigma{smooth:.0f}".replace(".", "p")
    out_bins = Path(args.out_bins_file) if str(args.out_bins_file).strip() else Path(__file__).parent / f"hubble_env_cf2_cone_{out_tag}_proxy_bins_v0.csv"
    out_merged = Path(args.out_merged_file) if str(args.out_merged_file).strip() else Path(__file__).parent / f"hubble_env_data_v1_cf2_cone_{out_tag}_proxy_v0.csv"
    out_json = Path(args.out_results_json) if str(args.out_results_json).strip() else Path(__file__).parent / f"cf2_cone_{out_tag}_proxy_bins_results_v0.json"
    out_plot = Path(args.out_plot) if str(args.out_plot).strip() else Path(__file__).parent / f"cf2_cone_{out_tag}_proxy_bins_plot_v0.png"

    # Download CF2
    table = cf2.CF2_TABLE2 if str(args.table).strip() == "table2" else cf2.CF2_TABLE1
    cache_path = cf2._cache_dir() / f"cf2_{Path(table).name}_asu_tsv_outmax_{int(args.out_max)}.tsv"
    url = cf2._build_cf2_request_url(table=table, out_max=int(args.out_max))
    if bool(args.force_download) and cache_path.exists():
        cache_path.unlink()
    cf2._download_if_missing(url, cache_path)

    raw = cf2._read_cf2_asu_tsv(cache_path, table=table, group_base_dist_frac_err=float(args.group_base_dist_frac_err))

    # Select and compute per-object proxy values
    n_raw = len(raw)
    n_after_filters = 0
    deltas: List[float] = []
    h0s: List[float] = []
    h0_errs: List[float] = []
    angs: List[float] = []

    for r in raw:
        if r.dist_mpc < float(args.dist_min) or r.dist_mpc > float(args.dist_max):
            continue

        # Flow-corrected velocity selection (table2 only, if available)
        v_sel = float(r.vcmba_kms)
        if str(args.velocity_mode) == "hubble_flow" and table.endswith("/table2") and r.pva_kms is not None:
            v_sel = float(r.vcmba_kms) - float(r.pva_kms)

        if v_sel < float(args.vmin) or v_sel > float(args.vmax):
            continue
        if table.endswith("/table2") and (r.n1 is None or int(r.n1) < int(args.group_min_n1)):
            continue
        if r.dist_frac_err is not None and (r.dist_frac_err <= 0 or r.dist_frac_err > float(args.err_max)):
            continue

        ang = _angular_sep_deg(float(r.glon_deg), float(r.glat_deg), cone.l_deg, cone.b_deg)
        if ang > float(cone.radius_deg):
            continue

        # 2M++ lookup radius (comoving Mpc/h)
        if str(args.r_from) == "vcmba":
            r_mpc_h = float(r.vcmba_kms) / 100.0
        else:
            r_mpc_h = float(r.dist_mpc) * float(args.h_ref)
        if r_mpc_h <= 0 or r_mpc_h > 200.0:
            continue

        try:
            lookup = lookup_delta_g_star(
                l_deg=float(r.glon_deg),
                b_deg=float(r.glat_deg),
                r_mpc_h=float(r_mpc_h),
                smoothing_gaussian_mpc_h=float(smooth),
            )
            delta_g = float(lookup["delta_g_star"])
        except Exception:
            continue
        if not math.isfinite(delta_g):
            continue

        # H0 computation
        v_h0 = float(r.vcmba_kms)
        if str(args.velocity_mode) == "hubble_flow" and table.endswith("/table2") and r.pva_kms is not None:
            v_h0 = float(r.vcmba_kms) - float(r.pva_kms)
        if not math.isfinite(v_h0) or v_h0 <= 0:
            continue

        H0 = float(v_h0) / float(r.dist_mpc)
        sigma_v = float(args.peculiar_vel_err)
        if r.vel_disp_sigma_kms is not None and math.isfinite(float(r.vel_disp_sigma_kms)):
            sigma_v = math.sqrt(sigma_v * sigma_v + float(r.vel_disp_sigma_kms) * float(r.vel_disp_sigma_kms))
        dist_frac_err = float(r.dist_frac_err) if r.dist_frac_err is not None else float(args.err_max)
        H0_err = math.sqrt((H0 * dist_frac_err) ** 2 + (sigma_v / float(r.dist_mpc)) ** 2)
        if not (math.isfinite(H0) and math.isfinite(H0_err) and H0_err > 0):
            continue

        n_after_filters += 1
        deltas.append(delta_g)
        h0s.append(H0)
        h0_errs.append(H0_err)
        angs.append(float(ang))

    deltas_np = np.asarray(deltas, dtype=float)
    h0_np = np.asarray(h0s, dtype=float)
    h0e_np = np.asarray(h0_errs, dtype=float)

    # Build quantile bins (fallback to a single aggregate point if too few objects)
    if deltas_np.size == 0:
        raise RuntimeError("No CF2 objects survived filters in the cone. Try increasing cone radius or relaxing cuts.")
    if deltas_np.size < int(args.n_bins) * max(1, int(args.min_per_bin)):
        print(
            f"Warning: filtered rows={int(deltas_np.size)} may be too small for stable binning "
            f"(n_bins={int(args.n_bins)}, min_per_bin={int(args.min_per_bin)})."
        )

    qs = np.linspace(0.0, 1.0, int(args.n_bins) + 1)
    edges = np.quantile(deltas_np, qs) if deltas_np.size else np.array([])
    edges = np.unique(edges)

    bin_rows: List[Dict[str, object]] = []
    bin_debug: List[Dict[str, Any]] = []
    total_used = 0

    if edges.size >= 3:
        for bi in range(edges.size - 1):
            lo = float(edges[bi])
            hi = float(edges[bi + 1])
            if bi == edges.size - 2:
                m = (deltas_np >= lo) & (deltas_np <= hi)
            else:
                m = (deltas_np >= lo) & (deltas_np < hi)
            idx = np.where(m)[0]
            n = int(idx.size)
            if n < int(args.min_per_bin):
                continue
            total_used += n

            d = deltas_np[idx]
            y = h0_np[idx]
            ye = h0e_np[idx]
            w = 1.0 / (ye * ye)

            d_mean = cf2._weighted_mean(d, w)
            y_mean = cf2._weighted_mean(y, w)
            y_se = cf2._weighted_se_mean(w)
            if not (math.isfinite(d_mean) and math.isfinite(y_mean) and math.isfinite(y_se)):
                continue

            env_id = f"cf2_cone_{cone.name}_bin_{bi+1:02d}_sigma{smooth:.0f}"
            env_name = f"CF2 cone[{cone.name}] bin {bi+1} (δ_g*∈[{lo:.3f},{hi:.3f}], n={n})"
            bin_rows.append(
                {
                    "env_id": env_id,
                    "env_name": env_name,
                    "category": cf2._assign_category(float(d_mean)),
                    "data_kind": "observed",
                    "use_for_fit": True,
                    "is_global_anchor": False,
                    "delta": float(d_mean),
                    "delta_err": None,
                    "H0": float(y_mean),
                    "H0_err": float(y_se),
                    "redshift_range": "local",
                    "scale_mpc": f"Cone {cone.radius_deg:.0f}deg @ (l,b)=({cone.l_deg:.1f},{cone.b_deg:.1f}); CF2 {Path(table).name}; smooth={smooth:.0f} Mpc/h",
                    "source": "Cosmicflows-2 (Tully et al. 2013; VizieR J/AJ/146/86) + 2M++ (Carrick et al. 2015)",
                    "notes": (
                        f"Within-cone quantile bin. H0=V/Dist with velocity_mode={str(args.velocity_mode)} "
                        f"and peculiar_vel_err={float(args.peculiar_vel_err):.0f} km/s."
                    ),
                    "delta_kind": "delta_g_star_2mpp",
                    "delta_scale_mpc_h": float(smooth),
                    "delta_source": str(cache_path),
                }
            )
            bin_debug.append(
                {"bin_index": bi + 1, "lo": lo, "hi": hi, "n": n, "delta_mean": d_mean, "H0_mean": y_mean, "H0_se": y_se}
            )

    if not bin_rows:
        # Fallback: single weighted-mean point for overlay (not a fit point).
        w_all = 1.0 / np.square(h0e_np)
        d_mean = cf2._weighted_mean(deltas_np, w_all)
        y_mean = cf2._weighted_mean(h0_np, w_all)
        y_se = cf2._weighted_se_mean(w_all)
        total_used = int(deltas_np.size)
        bin_debug.append(
            {
                "bin_index": 1,
                "lo": float(np.min(deltas_np)),
                "hi": float(np.max(deltas_np)),
                "n": int(deltas_np.size),
                "delta_mean": d_mean,
                "H0_mean": y_mean,
                "H0_se": y_se,
                "note": "fallback_single_aggregate (not enough samples for stable quantile bins)",
            }
        )
        bin_rows = [
            {
                "env_id": f"cf2_cone_{cone.name}_aggregate_sigma{smooth:.0f}",
                "env_name": f"CF2 cone[{cone.name}] aggregate (n={int(deltas_np.size)})",
                "category": cf2._assign_category(float(d_mean)),
                "data_kind": "observed",
                "use_for_fit": False,
                "is_global_anchor": False,
                "delta": float(d_mean),
                "delta_err": None,
                "H0": float(y_mean),
                "H0_err": float(y_se),
                "redshift_range": "local",
                "scale_mpc": f"Cone {cone.radius_deg:.0f}deg @ (l,b)=({cone.l_deg:.1f},{cone.b_deg:.1f}); CF2 {Path(table).name}; smooth={smooth:.0f} Mpc/h",
                "source": "Cosmicflows-2 (Tully et al. 2013; VizieR J/AJ/146/86) + 2M++ (Carrick et al. 2015)",
                "notes": (
                    f"Fallback aggregate point (not for fit). H0=V/Dist with velocity_mode={str(args.velocity_mode)} "
                    f"and peculiar_vel_err={float(args.peculiar_vel_err):.0f} km/s."
                ),
                "delta_kind": "delta_g_star_2mpp",
                "delta_scale_mpc_h": float(smooth),
                "delta_source": str(cache_path),
            }
        ]

    # Write bins-only CSV
    cf2._write_v1_rows(out_bins, bin_rows)

    # Write merged dataset
    base_path = Path(args.base_data_file)
    cf2._merge_v1_dataset(base_path, bin_rows, out_merged)

    # Fit bins trend (proxy-only, for direction check; only meaningful if bins are fit points)
    d_fit = np.array([float(r["delta"]) for r in bin_rows], dtype=float)
    y_fit = np.array([float(r["H0"]) for r in bin_rows], dtype=float)
    ye_fit = np.array([float(r["H0_err"]) for r in bin_rows], dtype=float)
    fit = _fit_line(d_fit, y_fit, ye_fit) if sum(bool(r.get("use_for_fit")) for r in bin_rows) >= 2 else {"error": "not_applicable", "n_points": int(len(bin_rows))}

    # Plot
    fig, ax = plt.subplots(figsize=(9.2, 6.1))
    ax.scatter(deltas_np, h0_np, s=10, alpha=0.18, color="#888888", label=f"CF2 objects in cone (n={int(deltas_np.size)})")
    ax.errorbar(d_fit, y_fit, yerr=ye_fit, fmt="o", color="#1f77b4", capsize=3, label=f"binned means (n={len(bin_rows)})")
    if not fit.get("error"):
        a = float(fit["H0_global_like_a"])
        beta = float(fit["beta_like"])
        xs = np.linspace(float(np.min(d_fit)) - 0.05, float(np.max(d_fit)) + 0.05, 200)
        ys = a * (1.0 - beta * xs)
        ax.plot(xs, ys, color="#111111", linewidth=1.6, label=f"fit: a={a:.2f}, beta_like={beta:.3f}")
    ax.axvline(0.0, color="#999999", linestyle="--", linewidth=0.9)
    ax.set_xlabel("δ_g* (2M++ proxy)")
    ax.set_ylabel("H0 (km/s/Mpc)  [proxy from CF2]")
    ax.set_title(f"CF2 cone proxy bins ({cone.name}; σ={smooth:.0f} Mpc/h)")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=9)
    fig.tight_layout()
    out_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_plot, dpi=150)
    plt.close(fig)

    results: Dict[str, Any] = {
        "metadata": {
            "timestamp": datetime.now().isoformat(),
            "script": Path(__file__).name,
        },
        "config": {
            "cone": {"name": cone.name, "l_deg": cone.l_deg, "b_deg": cone.b_deg, "radius_deg": cone.radius_deg},
            "smooth_mpc_h": smooth,
            "table": str(Path(table).name),
            "velocity_mode": str(args.velocity_mode),
            "r_from": str(args.r_from),
            "h_ref": float(args.h_ref),
            "dist_min_mpc": float(args.dist_min),
            "dist_max_mpc": float(args.dist_max),
            "vmin_kms": float(args.vmin),
            "vmax_kms": float(args.vmax),
            "err_max": float(args.err_max),
            "peculiar_vel_err_kms": float(args.peculiar_vel_err),
            "n_bins": int(args.n_bins),
            "min_per_bin": int(args.min_per_bin),
        },
        "counts": {
            "n_cf2_rows": int(n_raw),
            "n_selected": int(n_after_filters),
            "n_used_in_bins": int(total_used),
            "n_bins_output": int(len(bin_rows)),
        },
        "bin_edges": edges.tolist(),
        "bins": bin_debug,
        "fit_bins": fit,
        "outputs": {"bins_csv": str(out_bins), "merged_csv": str(out_merged), "json": str(out_json), "plot": str(out_plot)},
        "notes": [
            "This proxy track is exploratory and not equivalent to observer-environment bias.",
            "If beta_like is negative here, it usually indicates residual flow/systematic effects dominating H0=V/Dist.",
        ],
    }
    out_json.write_text(json.dumps(_json_sanitize(results), indent=2, ensure_ascii=False), encoding="utf-8")

    print("CF2 cone proxy bins completed.")
    print(f"Selected objects: {int(deltas_np.size)} (from raw {n_raw})")
    print(f"Bins written: {out_bins}")
    print(f"Merged dataset: {out_merged}")
    print(f"Results JSON: {out_json}")
    print(f"Plot: {out_plot}")
    if fit.get("error"):
        print(f"Fit: error={fit.get('error')}")
    else:
        print(f"Fit beta_like={fit.get('beta_like')}")


if __name__ == "__main__":
    main()

