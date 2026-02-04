"""
Pantheon+ low-z environment correlation (v0)
===========================================

Goal
----
Use real low-redshift SN Ia distances (Pantheon+ / SH0ES distance moduli) and the
public 2M++ reconstructed density field (delta_g*) to test whether there is an
observable trend between local density proxy and inferred H0.

This is an *observational proxy test*; it is not identical to the "observer
environment bias" mechanism. However, it provides a data-driven check:

  Do low-z inferred H0 values show a monotonic negative correlation with delta_g* ?

Data source
-----------
Pantheon+ DataRelease (public):
  https://github.com/PantheonPlusSH0ES/DataRelease
File:
  Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES.dat

We download the .dat file into math/hubble_env_scan/_cache/ if missing.

Outputs
-------
- pantheonplus_env_correlation_results_v0.json
- pantheonplus_env_correlation_plot_v0.png

Run (repo root)
---------------
python math/hubble_env_scan/run_pantheonplus_env_correlation_v0.py
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

import run_hubble_env_scan_v1 as scan
import twompp_density_lookup_v0 as twompp


PANTHEONPLUS_URL = (
    "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
    "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat"
)

C_KMS = 299792.458

OUTPUT_JSON = Path(__file__).parent / "pantheonplus_env_correlation_results_v0.json"
OUTPUT_PLOT = Path(__file__).parent / "pantheonplus_env_correlation_plot_v0.png"


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


def _parse_float(x: str) -> Optional[float]:
    t = str(x).strip()
    if not t or t.lower() == "nan":
        return None
    try:
        return float(t)
    except Exception:
        return None


@dataclass
class PpRow:
    cid: str
    idsurvey: int
    zhd: float
    zcmb: float
    mu_sh0es: float
    mu_sh0es_err_diag: float
    ra_deg: float
    dec_deg: float
    is_calibrator: int
    used_in_sh0es_hf: int


def read_pantheonplus_dat(path: Path) -> List[PpRow]:
    with path.open("r", encoding="utf-8") as f:
        header = f.readline().split()
        idx = {name: i for i, name in enumerate(header)}

        required = ["CID", "IDSURVEY", "zHD", "zCMB", "MU_SH0ES", "MU_SH0ES_ERR_DIAG", "RA", "DEC", "IS_CALIBRATOR", "USED_IN_SH0ES_HF"]
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
                mu = float(parts[idx["MU_SH0ES"]])
                mu_err = float(parts[idx["MU_SH0ES_ERR_DIAG"]])
                ra = float(parts[idx["RA"]])
                dec = float(parts[idx["DEC"]])
                is_cal = int(float(parts[idx["IS_CALIBRATOR"]]))
                used_hf = int(float(parts[idx["USED_IN_SH0ES_HF"]]))
            except Exception:
                continue
            # MU_SH0ES is -9 for some entries; treat as missing.
            if not math.isfinite(mu) or mu <= 0:
                continue
            if not math.isfinite(mu_err) or mu_err <= 0:
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
                    mu_sh0es=float(mu),
                    mu_sh0es_err_diag=float(mu_err),
                    ra_deg=float(ra),
                    dec_deg=float(dec),
                    is_calibrator=int(is_cal),
                    used_in_sh0es_hf=int(used_hf),
                )
            )
    return rows


def dedup_by_cid(rows: List[PpRow]) -> List[PpRow]:
    best: Dict[str, PpRow] = {}
    for r in rows:
        prev = best.get(r.cid)
        if prev is None or float(r.mu_sh0es_err_diag) < float(prev.mu_sh0es_err_diag):
            best[r.cid] = r
    return list(best.values())


def mu_to_distance_mpc(mu: np.ndarray) -> np.ndarray:
    mu = np.asarray(mu, dtype=float)
    return np.power(10.0, (mu - 25.0) / 5.0)


def mu_err_to_frac_dist_err(mu_err: np.ndarray) -> np.ndarray:
    # D ~ 10^(mu/5) => sigma_D/D = ln(10)/5 * sigma_mu
    mu_err = np.asarray(mu_err, dtype=float)
    return (math.log(10.0) / 5.0) * mu_err


def compute_delta_g_star(
    rows: List[PpRow],
    *,
    smoothing_mpc_h: float,
    r_from: str,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Compute delta_g* for each SN row."""
    smoothing_mpc_h = float(smoothing_mpc_h)
    if smoothing_mpc_h < 4.0:
        raise ValueError("smoothing_mpc_h must be >= 4.0 (2M++ base smoothing)")

    deltas = np.full(len(rows), np.nan, dtype=float)
    n_ok = 0
    n_oob = 0

    for i, r in enumerate(rows):
        # Convert RA/Dec -> Galactic l/b
        l_deg, b_deg = twompp._radec_to_galactic_deg(float(r.ra_deg), float(r.dec_deg))
        # Set comoving r (Mpc/h) from redshift velocity
        if r_from == "zHD":
            r_mpc_h = (C_KMS * float(r.zhd)) / 100.0
        else:
            r_mpc_h = (C_KMS * float(r.zcmb)) / 100.0
        try:
            out = twompp.lookup_delta_g_star(
                l_deg=float(l_deg),
                b_deg=float(b_deg),
                r_mpc_h=float(r_mpc_h),
                smoothing_gaussian_mpc_h=float(smoothing_mpc_h),
            )
            deltas[i] = float(out["delta_g_star"])
            n_ok += 1
        except Exception:
            n_oob += 1
            deltas[i] = float("nan")

    meta = {"n_ok": int(n_ok), "n_failed": int(n_oob), "smoothing_mpc_h": float(smoothing_mpc_h), "r_from": str(r_from)}
    return deltas, meta


def make_bins(x: np.ndarray, y: np.ndarray, yerr: np.ndarray, *, n_bins: int) -> List[Dict[str, Any]]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    yerr = np.asarray(yerr, dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(yerr) & (yerr > 0)
    x = x[m]
    y = y[m]
    yerr = yerr[m]
    if x.size < max(20, n_bins * 10):
        n_bins = max(3, min(int(n_bins), max(3, int(x.size // 10))))
    qs = np.linspace(0.0, 1.0, int(n_bins) + 1)
    edges = np.unique(np.quantile(x, qs))
    bins: List[Dict[str, Any]] = []
    for bi in range(edges.size - 1):
        lo = float(edges[bi])
        hi = float(edges[bi + 1])
        sel = (x >= lo) & (x <= hi) if bi == edges.size - 2 else (x >= lo) & (x < hi)
        if int(np.sum(sel)) < 10:
            continue
        xx = x[sel]
        yy = y[sel]
        ye = yerr[sel]
        w = 1.0 / (ye * ye)
        x_mean = float(np.sum(w * xx) / np.sum(w))
        y_mean = float(np.sum(w * yy) / np.sum(w))
        y_se = float(math.sqrt(1.0 / float(np.sum(w))))
        bins.append({"lo": lo, "hi": hi, "n": int(xx.size), "x_mean": x_mean, "y_mean": y_mean, "y_se": y_se})
    return bins


def main() -> None:
    ap = argparse.ArgumentParser(description="Pantheon+ low-z environment correlation (v0)")
    ap.add_argument("--data-file", type=str, default="", help="Optional local Pantheon+SH0ES.dat path; if empty, download")
    ap.add_argument("--z-min", type=float, default=0.015)
    ap.add_argument("--z-max", type=float, default=0.060)
    ap.add_argument("--z-field", type=str, default="zCMB", choices=["zCMB", "zHD"])
    ap.add_argument("--exclude-calibrators", action="store_true", help="Exclude IS_CALIBRATOR==1 (recommended)")
    ap.add_argument("--only-sh0es-hf", action="store_true", help="Keep only USED_IN_SH0ES_HF==1")
    ap.add_argument("--dedup-cid", action="store_true", help="Deduplicate by CID (keep smallest MU_SH0ES_ERR_DIAG)")
    ap.add_argument("--smooth-list", type=str, default="20,50,100", help="Comma-separated 2M++ smoothing sigmas (Mpc/h)")
    ap.add_argument("--r-from", type=str, default="zCMB", choices=["zCMB", "zHD"], help="Use zCMB or zHD to set r for 2M++ lookup")
    ap.add_argument("--n-bins", type=int, default=8, help="Quantile bins for plotting binned means")
    ap.add_argument("--output-json", type=str, default=str(OUTPUT_JSON))
    ap.add_argument("--output-plot", type=str, default=str(OUTPUT_PLOT))
    args = ap.parse_args()

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

    z_field = str(args.z_field)
    if z_field == "zHD":
        rows = [r for r in rows if float(args.z_min) <= float(r.zhd) <= float(args.z_max)]
    else:
        rows = [r for r in rows if float(args.z_min) <= float(r.zcmb) <= float(args.z_max)]

    if not rows:
        raise SystemExit("No rows after filtering. Try adjusting z-range or flags.")

    # Distances and per-SN inferred H0
    mu = np.array([r.mu_sh0es for r in rows], dtype=float)
    mu_err = np.array([r.mu_sh0es_err_diag for r in rows], dtype=float)
    D_mpc = mu_to_distance_mpc(mu)
    fracD = mu_err_to_frac_dist_err(mu_err)

    v_zhd = C_KMS * np.array([r.zhd for r in rows], dtype=float)
    v_zcmb = C_KMS * np.array([r.zcmb for r in rows], dtype=float)

    H0_zhd = v_zhd / D_mpc
    H0_zcmb = v_zcmb / D_mpc

    # For plotting: uncertainty on H0 from distance-modulus diag only
    H0_zhd_err = np.abs(H0_zhd) * fracD
    H0_zcmb_err = np.abs(H0_zcmb) * fracD

    # Smoothing list
    smooth_list: List[float] = []
    for t in str(args.smooth_list).split(","):
        v = _parse_float(t)
        if v is not None:
            smooth_list.append(float(v))
    smooth_list = sorted(set(smooth_list)) if smooth_list else [20.0, 50.0, 100.0]

    # Compute delta_g* for each smoothing scale
    per_smooth: Dict[str, Any] = {}
    for s in smooth_list:
        d, meta = compute_delta_g_star(rows, smoothing_mpc_h=float(s), r_from=str(args.r_from))

        # Correlations and linear fits (proxy only)
        m = np.isfinite(d) & np.isfinite(H0_zhd) & np.isfinite(H0_zhd_err) & (H0_zhd_err > 0)
        n_use = int(np.sum(m))
        d_use = d[m]

        def _corr(a: np.ndarray, b: np.ndarray) -> Optional[float]:
            if a.size < 3 or np.std(a) <= 0 or np.std(b) <= 0:
                return None
            return float(np.corrcoef(a, b)[0, 1])

        fit_zhd = scan.fit_ab_weighted(d_use, H0_zhd[m], H0_zhd_err[m])
        fit_zcmb = scan.fit_ab_weighted(d_use, H0_zcmb[m], H0_zcmb_err[m])

        # Derive beta-like from the (a,b) fit: H0 = a - b*delta => beta=b/a
        def _beta_like(a: float, b: float) -> Optional[float]:
            if not (math.isfinite(a) and math.isfinite(b)) or a == 0:
                return None
            return float(b / a)

        a_zhd, b_zhd, cov_zhd, chi2_zhd, chi2dof_zhd, corr_zhd = fit_zhd
        a_zcmb, b_zcmb, cov_zcmb, chi2_zcmb, chi2dof_zcmb, corr_zcmb = fit_zcmb

        per_smooth[str(int(round(s)))] = {
            "delta_meta": meta,
            "n_used": n_use,
            "corr": {
                "delta_vs_H0_zhd": _corr(d_use, H0_zhd[m]),
                "delta_vs_H0_zcmb": _corr(d_use, H0_zcmb[m]),
            },
            "fit": {
                "zhd": {
                    "a": float(a_zhd),
                    "b": float(b_zhd),
                    "beta_like": _beta_like(float(a_zhd), float(b_zhd)),
                    "chi2": _json_sanitize(chi2_zhd),
                    "chi2_dof": _json_sanitize(chi2dof_zhd),
                    "corr_delta_H0": _json_sanitize(corr_zhd),
                    "cov_ab": cov_zhd.tolist() if cov_zhd is not None else None,
                },
                "zcmb": {
                    "a": float(a_zcmb),
                    "b": float(b_zcmb),
                    "beta_like": _beta_like(float(a_zcmb), float(b_zcmb)),
                    "chi2": _json_sanitize(chi2_zcmb),
                    "chi2_dof": _json_sanitize(chi2dof_zcmb),
                    "corr_delta_H0": _json_sanitize(corr_zcmb),
                    "cov_ab": cov_zcmb.tolist() if cov_zcmb is not None else None,
                },
            },
            "bins": {
                "zhd": make_bins(d_use, H0_zhd[m], H0_zhd_err[m], n_bins=int(args.n_bins)),
                "zcmb": make_bins(d_use, H0_zcmb[m], H0_zcmb_err[m], n_bins=int(args.n_bins)),
            },
        }

    # Plot: 2 columns (zHD, zCMB) × up to 3 smoothing rows (20/50/100)
    n_rows = min(3, len(smooth_list))
    n_cols = 2
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(11.5, 3.6 * n_rows), sharex=False, sharey=False)
    if n_rows == 1:
        axes = np.array([axes])  # type: ignore

    for ri, s in enumerate(smooth_list[:n_rows]):
        key = str(int(round(s)))
        d, _meta = compute_delta_g_star(rows, smoothing_mpc_h=float(s), r_from=str(args.r_from))
        m = np.isfinite(d) & np.isfinite(H0_zhd) & np.isfinite(H0_zhd_err) & (H0_zhd_err > 0)
        d_use = d[m]

        for ci, (label, y, ye) in enumerate(
            [
                ("zHD", H0_zhd[m], H0_zhd_err[m]),
                ("zCMB", H0_zcmb[m], H0_zcmb_err[m]),
            ]
        ):
            ax = axes[ri, ci]
            ax.scatter(d_use, y, s=10, alpha=0.20, color="#777777", label=f"SN (n={int(d_use.size)})")

            # binned means
            bins = per_smooth[key]["bins"][label.lower()]
            if bins:
                bx = np.array([b["x_mean"] for b in bins], dtype=float)
                by = np.array([b["y_mean"] for b in bins], dtype=float)
                be = np.array([b["y_se"] for b in bins], dtype=float)
                ax.errorbar(bx, by, yerr=be, fmt="o", capsize=3, color="#1f77b4", label=f"binned (n={len(bins)})")

            fit = per_smooth[key]["fit"][label.lower()]
            if fit.get("beta_like") is not None and fit.get("a") is not None:
                a = float(fit["a"])
                beta_like = float(fit["beta_like"])
                xs = np.linspace(float(np.nanmin(d_use)) - 0.05, float(np.nanmax(d_use)) + 0.05, 250)
                ys = a * (1.0 - beta_like * xs)
                ax.plot(xs, ys, color="#111111", linewidth=1.4, label=f"fit: beta_like={beta_like:.3f}")

            ax.axvline(0.0, color="#999999", linestyle="--", linewidth=0.9)
            ax.set_xlabel(f"δ_g* (2M++; σ={float(s):.0f} Mpc/h)")
            ax.set_ylabel(f"H0 from {label} / MU_SH0ES")
            ax.set_title(f"Pantheon+ low-z ({label}); σ={float(s):.0f} Mpc/h")
            ax.grid(alpha=0.2)
            ax.legend(fontsize=9)

    fig.tight_layout()
    Path(args.output_plot).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(args.output_plot), dpi=150)
    plt.close(fig)

    out: Dict[str, Any] = {
        "metadata": {"timestamp": datetime.now().isoformat(), "script": Path(__file__).name},
        "inputs": {"pantheonplus_url": PANTHEONPLUS_URL, "data_file": str(data_path)},
        "filters": {
            "z_field": z_field,
            "z_min": float(args.z_min),
            "z_max": float(args.z_max),
            "exclude_calibrators": bool(args.exclude_calibrators),
            "only_sh0es_hf": bool(args.only_sh0es_hf),
            "dedup_cid": bool(args.dedup_cid),
            "r_from_for_2mpp": str(args.r_from),
        },
        "counts": {"n_raw": int(n_raw), "n_after_filters": int(len(rows))},
        "summary": {
            "H0_stats": {
                "zhd": {
                    "mean": float(np.mean(H0_zhd)),
                    "std": float(np.std(H0_zhd, ddof=1)) if H0_zhd.size >= 2 else 0.0,
                    "median": float(np.median(H0_zhd)),
                },
                "zcmb": {
                    "mean": float(np.mean(H0_zcmb)),
                    "std": float(np.std(H0_zcmb, ddof=1)) if H0_zcmb.size >= 2 else 0.0,
                    "median": float(np.median(H0_zcmb)),
                },
            }
        },
        "per_smoothing": per_smooth,
        "outputs": {"json": str(Path(args.output_json)), "plot": str(Path(args.output_plot))},
        "notes": [
            "H0 is inferred as (c*z)/D where D is from MU_SH0ES (SH0ES-calibrated distance modulus).",
            "This is a proxy test; it does not isolate observer-environment bias uniquely.",
            "zHD includes VPEC corrections; zCMB does not. Comparing them helps diagnose flow-correction impact.",
        ],
    }

    Path(args.output_json).write_text(json.dumps(_json_sanitize(out), indent=2, ensure_ascii=False), encoding="utf-8")
    print("Pantheon+ environment correlation completed.")
    print(f"Rows used: {len(rows)} (raw={n_raw})")
    print(f"Output: {args.output_json}")
    print(f"Plot: {args.output_plot}")


if __name__ == "__main__":
    main()

