"""
Hubble tension environment scan (v1)
====================================

Goal
----
Extend the v0 environment scan by adding:

- Monte-Carlo uncertainty propagation (delta_err, H0_err) -> uncertainty bands for (H0_global, beta)
- Predictive intervals for forecast / inferred points
- Optional overlay of overdense-environment candidate points with delta proxies (2M++)

Model
-----
    H0(delta) = H0_global * (1 - beta * delta)

Equivalent linear form:
    H0 = a - b * delta, where a=H0_global and b=a*beta.

Run (repo root)
--------------
python math/hubble_env_scan/run_hubble_env_scan_v1.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np


DATA_FILE = Path(__file__).parent / "hubble_env_data_v0.csv"
OVERDENSE_CANDIDATES_FILE = Path(__file__).parent / "overdense_env_candidates_v0.csv"

OUTPUT_JSON = Path(__file__).parent / "results_v1.json"
OUTPUT_PLOT = Path(__file__).parent / "plot_v1.png"


@dataclass
class EnvPoint:
    env_id: str
    env_name: str
    category: str
    data_kind: str
    use_for_fit: bool
    is_global_anchor: bool
    delta: Optional[float]
    delta_err: Optional[float]
    H0: Optional[float]
    H0_err: Optional[float]
    redshift_range: str
    scale_mpc: str
    source: str
    notes: str
    delta_kind: Optional[str] = None
    delta_scale_mpc_h: Optional[float] = None
    delta_source: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "env_id": self.env_id,
            "env_name": self.env_name,
            "category": self.category,
            "data_kind": self.data_kind,
            "use_for_fit": self.use_for_fit,
            "is_global_anchor": self.is_global_anchor,
            "delta": self.delta,
            "delta_err": self.delta_err,
            "H0": self.H0,
            "H0_err": self.H0_err,
            "redshift_range": self.redshift_range,
            "scale_mpc": self.scale_mpc,
            "source": self.source,
            "notes": self.notes,
            "delta_kind": self.delta_kind,
            "delta_scale_mpc_h": self.delta_scale_mpc_h,
            "delta_source": self.delta_source,
        }


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


def _parse_float(value: str) -> Optional[float]:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _parse_bool(value: str) -> bool:
    text = (value or "").strip().lower()
    return text in {"true", "1", "yes", "y"}


def read_dataset(path: Path) -> List[EnvPoint]:
    points: List[EnvPoint] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            points.append(
                EnvPoint(
                    env_id=(row.get("env_id", "") or "").strip(),
                    env_name=(row.get("env_name", "") or "").strip(),
                    category=(row.get("category", "") or "").strip(),
                    data_kind=(row.get("data_kind", "") or "").strip(),
                    use_for_fit=_parse_bool(row.get("use_for_fit", "")),
                    is_global_anchor=_parse_bool(row.get("is_global_anchor", "")),
                    delta=_parse_float(row.get("delta", "")),
                    delta_err=_parse_float(row.get("delta_err", "")),
                    H0=_parse_float(row.get("H0", "")),
                    H0_err=_parse_float(row.get("H0_err", "")),
                    redshift_range=(row.get("redshift_range", "") or "").strip(),
                    scale_mpc=(row.get("scale_mpc", "") or "").strip(),
                    source=(row.get("source", "") or "").strip(),
                    notes=(row.get("notes", "") or "").strip(),
                    delta_kind=(row.get("delta_kind", "") or "").strip() or None,
                    delta_scale_mpc_h=_parse_float(row.get("delta_scale_mpc_h", "") or ""),
                    delta_source=(row.get("delta_source", "") or "").strip() or None,
                )
            )
    return points


def read_overdense_candidates(path: Path) -> List[Dict[str, Any]]:
    """
    Load optional "overdense environment candidates" table.

    The table has a different schema. We use it only for plotting overlays
    (not for fitting) because delta_proxy is not necessarily on the same scale
    as literature delta used in hubble_env_data_v0.csv.
    """
    if not path.exists():
        return []
    rows: List[Dict[str, Any]] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            delta = _parse_float(r.get("delta_proxy", "") or "")
            H0 = _parse_float(r.get("H0", "") or "")
            if delta is None or H0 is None:
                continue
            rows.append(
                {
                    "env_id": (r.get("env_id") or "").strip(),
                    "env_name": (r.get("env_name") or "").strip(),
                    "category": (r.get("category") or "").strip(),
                    "measurement_type": (r.get("measurement_type") or "").strip(),
                    "delta_proxy": delta,
                    "delta_proxy_kind": (r.get("delta_proxy_kind") or "").strip(),
                    "delta_proxy_smoothing_mpc_h": _parse_float(r.get("delta_proxy_smoothing_mpc_h") or ""),
                    "H0": H0,
                    "H0_err": _parse_float(r.get("H0_err") or ""),
                    "reference": (r.get("reference") or "").strip(),
                    "notes": (r.get("notes") or "").strip(),
                    "delta_proxy_source": (r.get("delta_proxy_source") or "").strip(),
                }
            )
    return rows


def get_global_anchor(points: List[EnvPoint]) -> Tuple[float, str]:
    anchors = [p for p in points if p.is_global_anchor and p.H0 is not None]
    if anchors:
        return float(anchors[0].H0), anchors[0].env_id
    raise ValueError("No global anchor with H0 found in dataset.")


def _get_fit_points(points: List[EnvPoint]) -> List[EnvPoint]:
    return [
        p
        for p in points
        if p.data_kind == "observed"
        and p.use_for_fit
        and p.delta is not None
        and p.H0 is not None
    ]


def _get_fit_points_direct(points: List[EnvPoint]) -> List[EnvPoint]:
    return [
        p
        for p in points
        if p.data_kind == "observed"
        and p.use_for_fit
        and p.delta is not None
        and p.H0 is not None
        and (p.delta_kind is None or p.delta_kind == "direct")
    ]


def _get_fit_points_direct_including_inferred(points: List[EnvPoint]) -> List[EnvPoint]:
    return [
        p
        for p in points
        if p.data_kind in {"observed", "inferred"}
        and p.use_for_fit
        and p.delta is not None
        and p.H0 is not None
        and (p.delta_kind is None or p.delta_kind == "direct")
    ]


def _get_fit_points_gavas(points: List[EnvPoint]) -> List[EnvPoint]:
    """Select Gavas-2024 Fig.3 inferred points (delta_den within 10 Mpc/h).

    These points are kept separate from the "direct" fit because their delta
    definition/scale differs from literature void deltas in the main dataset.
    """
    return [
        p
        for p in points
        if p.data_kind == "inferred"
        and p.use_for_fit
        and p.delta is not None
        and p.H0 is not None
        and p.delta_kind == "gavas_delta_den_10mpch"
    ]


def _get_fit_points_proxy(
    points: List[EnvPoint],
    *,
    scale_min: Optional[float],
    scale_max: Optional[float],
) -> List[EnvPoint]:
    out: List[EnvPoint] = []
    for p in points:
        if p.data_kind != "observed" or p.delta is None or p.H0 is None:
            continue
        # Respect the dataset "use_for_fit" flag for proxy-mode fits too.
        # Otherwise, "overlay-only" proxy points can silently leak into the fit.
        if not p.use_for_fit:
            continue
        if p.delta_kind != "delta_g_star_2mpp":
            continue
        if p.delta_scale_mpc_h is None:
            continue
        if scale_min is not None and float(p.delta_scale_mpc_h) < float(scale_min):
            continue
        if scale_max is not None and float(p.delta_scale_mpc_h) > float(scale_max):
            continue
        out.append(p)
    return out


def _compute_plot_xs(points: List[EnvPoint], *, n: int = 220) -> np.ndarray:
    deltas_all = [p.delta for p in points if p.delta is not None]
    if deltas_all:
        delta_min = float(min(deltas_all)) - 0.15
        delta_max = float(max(deltas_all)) + 0.15
    else:
        delta_min, delta_max = -0.8, 0.8
    return np.linspace(delta_min, delta_max, int(n))


def _fit_is_valid_for_plot(fit: Optional[Dict[str, Any]], *, min_points: int) -> bool:
    if not isinstance(fit, dict):
        return False
    if fit.get("error"):
        return False
    n_points = fit.get("n_points", 0)
    try:
        if int(n_points) < int(min_points):
            return False
    except Exception:
        return False
    try:
        H0_global = float(fit.get("H0_global"))
        beta = float(fit.get("beta"))
    except Exception:
        return False
    return bool(np.isfinite(H0_global) and np.isfinite(beta))


def fit_ab_weighted(
    deltas: np.ndarray,
    H0s: np.ndarray,
    sigmas: np.ndarray,
) -> Tuple[float, float, Optional[np.ndarray], Optional[float], Optional[float], Optional[float]]:
    """
    Weighted least squares fit: H0 = a - b*delta.
    Returns (a, b, cov_ab, chi2, chi2_dof, corr_delta_H0)
    """
    deltas = np.asarray(deltas, dtype=float)
    H0s = np.asarray(H0s, dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)

    weights = 1.0 / np.square(sigmas)
    X = np.column_stack([np.ones_like(deltas), -deltas])
    XtW = X.T * weights
    XtWX = XtW @ X
    XtWy = XtW @ H0s

    try:
        theta = np.linalg.solve(XtWX, XtWy)
    except np.linalg.LinAlgError:
        return float("nan"), float("nan"), None, None, None, None

    a = float(theta[0])
    b = float(theta[1])

    model = X @ theta
    residuals = H0s - model
    chi2 = float(np.sum(weights * np.square(residuals)))
    dof = max(len(H0s) - 2, 0)
    chi2_dof = (chi2 / dof) if dof > 0 else None

    cov = None
    try:
        cov = np.linalg.inv(XtWX)
        if chi2_dof is not None:
            cov = cov * chi2_dof
    except np.linalg.LinAlgError:
        cov = None

    corr = None
    if len(H0s) >= 2 and np.std(deltas) > 0 and np.std(H0s) > 0:
        corr = float(np.corrcoef(deltas, H0s)[0, 1])

    return a, b, cov, chi2, chi2_dof, corr


def fit_H0_global_and_beta(
    points: List[EnvPoint],
    *,
    missing_H0_err: Optional[float] = None,
    weight_mode: str = "h0_err",
    fit_points: Optional[List[EnvPoint]] = None,
) -> Dict[str, Any]:
    """
    Fit (H0_global, beta) using weighted least squares.

    weight_mode:
      - "h0_err": weights = 1/H0_err^2 (missing -> missing_H0_err or 1)
      - "equal":  equal weights
    """
    fit_points = fit_points if fit_points is not None else _get_fit_points(points)
    if len(fit_points) < 2:
        return {"n_points": len(fit_points), "error": "insufficient_points"}

    deltas = np.array([float(p.delta) for p in fit_points], dtype=float)
    H0s = np.array([float(p.H0) for p in fit_points], dtype=float)

    if np.unique(deltas).size < 2:
        return {"n_points": len(fit_points), "error": "insufficient_delta_diversity"}

    # default missing sigma: median of present errors
    present = [float(p.H0_err) for p in fit_points if p.H0_err is not None and np.isfinite(p.H0_err)]
    default_missing = float(np.median(np.asarray(present, dtype=float))) if present else 1.0
    missing_val = float(missing_H0_err) if missing_H0_err is not None else default_missing

    if weight_mode == "equal":
        sigmas = np.ones_like(H0s, dtype=float)
    else:
        sigmas = np.array(
            [float(p.H0_err) if p.H0_err is not None else float(missing_val) for p in fit_points],
            dtype=float,
        )
        sigmas = np.maximum(sigmas, 1e-9)

    a, b, cov, chi2, chi2_dof, corr = fit_ab_weighted(deltas, H0s, sigmas)

    beta = (b / a) if np.isfinite(a) and a != 0 else float("nan")
    a_sigma = None
    b_sigma = None
    beta_sigma = None

    if cov is not None:
        a_sigma = math.sqrt(float(cov[0, 0])) if float(cov[0, 0]) >= 0 else None
        b_sigma = math.sqrt(float(cov[1, 1])) if float(cov[1, 1]) >= 0 else None
        if np.isfinite(a) and a != 0 and np.isfinite(b):
            d_beta_da = -b / (a * a)
            d_beta_db = 1.0 / a
            var_beta = (
                d_beta_da * d_beta_da * float(cov[0, 0])
                + 2.0 * d_beta_da * d_beta_db * float(cov[0, 1])
                + d_beta_db * d_beta_db * float(cov[1, 1])
            )
            beta_sigma = math.sqrt(var_beta) if var_beta >= 0 else None

    return {
        "H0_global": a,
        "H0_global_sigma": a_sigma,
        "b": b,
        "b_sigma": b_sigma,
        "beta": beta,
        "beta_sigma": beta_sigma,
        "chi2": chi2,
        "dof": max(len(fit_points) - 2, 0),
        "chi2_dof": chi2_dof,
        "corr_delta_H0": corr,
        "n_points": len(fit_points),
        "weight_mode": weight_mode,
        "missing_H0_err_used": missing_val,
        "fit_points": [p.to_dict() for p in fit_points],
        "sigmas_used": sigmas.tolist(),
        "cov_ab": cov.tolist() if cov is not None else None,
    }


def fit_beta_fixed_anchor(
    points: List[EnvPoint],
    *,
    H0_anchor: float,
    missing_H0_err: Optional[float] = None,
    weight_mode: str = "h0_err",
    fit_points: Optional[List[EnvPoint]] = None,
) -> Dict[str, Any]:
    """
    Fit only beta, keeping H0_global fixed to the provided anchor.
    """
    fit_points = fit_points if fit_points is not None else _get_fit_points(points)
    if len(fit_points) < 2:
        return {"n_points": len(fit_points), "error": "insufficient_points"}

    deltas = np.array([float(p.delta) for p in fit_points], dtype=float)
    H0s = np.array([float(p.H0) for p in fit_points], dtype=float)

    if np.unique(deltas).size < 2:
        return {"n_points": len(fit_points), "error": "insufficient_delta_diversity"}

    present = [float(p.H0_err) for p in fit_points if p.H0_err is not None and np.isfinite(p.H0_err)]
    default_missing = float(np.median(np.asarray(present, dtype=float))) if present else 1.0
    missing_val = float(missing_H0_err) if missing_H0_err is not None else default_missing

    if weight_mode == "equal":
        sigmas = np.ones_like(H0s, dtype=float)
    else:
        sigmas = np.array(
            [float(p.H0_err) if p.H0_err is not None else float(missing_val) for p in fit_points],
            dtype=float,
        )
        sigmas = np.maximum(sigmas, 1e-9)

    w = 1.0 / np.square(sigmas)
    a = float(H0_anchor)
    # minimize sum w*(H0 - (a - b*delta))^2 -> b = sum w*delta*(a-H0) / sum w*delta^2
    denom = float(np.sum(w * deltas * deltas))
    if denom <= 0:
        return {"n_points": len(fit_points), "error": "degenerate_weights"}
    b = float(np.sum(w * deltas * (a - H0s)) / denom)
    beta = (b / a) if a != 0 else float("nan")

    # crude uncertainty for b (assuming a fixed): Var(b) ≈ 1 / sum(w*delta^2) * chi2_dof
    model = a - b * deltas
    resid = H0s - model
    chi2 = float(np.sum(w * resid * resid))
    dof = max(len(fit_points) - 1, 0)
    chi2_dof = (chi2 / dof) if dof > 0 else None
    b_sigma = None
    beta_sigma = None
    if chi2_dof is not None:
        var_b = (1.0 / denom) * chi2_dof
        b_sigma = math.sqrt(var_b) if var_b >= 0 else None
        beta_sigma = (b_sigma / abs(a)) if b_sigma is not None and a != 0 else None

    return {
        "H0_global_fixed": a,
        "b": b,
        "b_sigma": b_sigma,
        "beta": beta,
        "beta_sigma": beta_sigma,
        "chi2": chi2,
        "dof": dof,
        "chi2_dof": chi2_dof,
        "n_points": len(fit_points),
        "weight_mode": weight_mode,
        "missing_H0_err_used": missing_val,
        "fit_points": [p.to_dict() for p in fit_points],
        "sigmas_used": sigmas.tolist(),
    }


def _quantiles(x: np.ndarray, qs: Tuple[float, ...] = (0.16, 0.5, 0.84)) -> Dict[str, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {"p16": float("nan"), "p50": float("nan"), "p84": float("nan"), "mean": float("nan"), "std": float("nan")}
    return {
        "p16": float(np.quantile(x, qs[0])),
        "p50": float(np.quantile(x, qs[1])),
        "p84": float(np.quantile(x, qs[2])),
        "mean": float(np.mean(x)),
        "std": float(np.std(x, ddof=1)) if x.size >= 2 else 0.0,
    }


def mc_uncertainty_propagation(
    points: List[EnvPoint],
    *,
    n_samples: int,
    seed: int,
    missing_H0_err: Optional[float] = None,
    missing_delta_err: float = 0.0,
    include_delta_err: bool = True,
    weight_mode: str = "h0_err",
    fix_anchor: bool = False,
    H0_anchor: Optional[float] = None,
    fit_points: Optional[List[EnvPoint]] = None,
) -> Dict[str, Any]:
    """
    Monte-Carlo propagate measurement errors -> distribution of (H0_global, beta).
    """
    rng = np.random.default_rng(int(seed))
    fit_points = fit_points if fit_points is not None else _get_fit_points(points)
    if len(fit_points) < 2:
        return {"n_points": len(fit_points), "error": "insufficient_points"}
    deltas0 = np.array([float(p.delta) for p in fit_points], dtype=float)
    H0s0 = np.array([float(p.H0) for p in fit_points], dtype=float)

    if np.unique(deltas0).size < 2:
        return {"n_points": len(fit_points), "error": "insufficient_delta_diversity"}

    present_H0_err = [float(p.H0_err) for p in fit_points if p.H0_err is not None and np.isfinite(p.H0_err)]
    default_missing_H0 = float(np.median(np.asarray(present_H0_err, dtype=float))) if present_H0_err else 1.0
    missing_H0 = float(missing_H0_err) if missing_H0_err is not None else default_missing_H0

    H0_err_use = np.array(
        [float(p.H0_err) if p.H0_err is not None else float(missing_H0) for p in fit_points],
        dtype=float,
    )
    H0_err_use = np.maximum(H0_err_use, 1e-9)

    delta_err_use = np.array(
        [
            float(p.delta_err)
            if (p.delta_err is not None and np.isfinite(p.delta_err))
            else float(missing_delta_err)
            for p in fit_points
        ],
        dtype=float,
    )
    delta_err_use = np.maximum(delta_err_use, 0.0)

    # Pre-draw noise
    eps_H0 = rng.normal(loc=0.0, scale=1.0, size=(int(n_samples), len(fit_points)))
    eps_delta = rng.normal(loc=0.0, scale=1.0, size=(int(n_samples), len(fit_points))) if include_delta_err else None

    a_s = np.full(int(n_samples), np.nan, dtype=float)
    b_s = np.full(int(n_samples), np.nan, dtype=float)
    beta_s = np.full(int(n_samples), np.nan, dtype=float)

    for i in range(int(n_samples)):
        H0s = H0s0 + eps_H0[i] * H0_err_use
        if include_delta_err and eps_delta is not None:
            deltas = deltas0 + eps_delta[i] * delta_err_use
        else:
            deltas = deltas0

        if weight_mode == "equal":
            sigmas = np.ones_like(H0s, dtype=float)
        else:
            sigmas = H0_err_use

        if fix_anchor:
            a = float(H0_anchor) if H0_anchor is not None else float("nan")
            if not np.isfinite(a) or a == 0:
                continue
            w = 1.0 / np.square(sigmas)
            denom = float(np.sum(w * deltas * deltas))
            if denom <= 0:
                continue
            b = float(np.sum(w * deltas * (a - H0s)) / denom)
        else:
            a, b, _cov, _chi2, _chi2_dof, _corr = fit_ab_weighted(deltas, H0s, sigmas)

        if not np.isfinite(a) or not np.isfinite(b) or a == 0:
            continue

        a_s[i] = float(a)
        b_s[i] = float(b)
        beta_s[i] = float(b / a)

    return {
        "n_samples": int(n_samples),
        "seed": int(seed),
        "missing_H0_err_used": float(missing_H0),
        "missing_delta_err_used": float(missing_delta_err),
        "include_delta_err": bool(include_delta_err),
        "weight_mode": weight_mode,
        "fix_anchor": bool(fix_anchor),
        "H0_anchor_used": float(H0_anchor) if H0_anchor is not None else None,
        "summary": {
            "H0_global": _quantiles(a_s),
            "b": _quantiles(b_s),
            "beta": _quantiles(beta_s),
        },
        # Keep samples for plotting bands (but avoid dumping huge arrays to JSON)
        "samples_for_plot": {
            "H0_global": a_s.tolist(),
            "beta": beta_s.tolist(),
        },
    }


def predict_band(xs: np.ndarray, a_samples: np.ndarray, beta_samples: np.ndarray) -> Dict[str, Any]:
    xs = np.asarray(xs, dtype=float)
    a = np.asarray(a_samples, dtype=float)
    b = np.asarray(beta_samples, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    a = a[m]
    b = b[m]
    if a.size == 0:
        return {"p16": None, "p50": None, "p84": None}
    ys = a[:, None] * (1.0 - b[:, None] * xs[None, :])
    return {
        "p16": np.quantile(ys, 0.16, axis=0),
        "p50": np.quantile(ys, 0.50, axis=0),
        "p84": np.quantile(ys, 0.84, axis=0),
    }


def plot_results(
    points: List[EnvPoint],
    fit: Dict[str, Any],
    *,
    band: Optional[Dict[str, Any]] = None,
    overdense_candidates: Optional[List[Dict[str, Any]]] = None,
    fit_proxy: Optional[Dict[str, Any]] = None,
    band_proxy: Optional[Dict[str, Any]] = None,
    fit_direct_inferred: Optional[Dict[str, Any]] = None,
    fit_gavas: Optional[Dict[str, Any]] = None,
    plot_forecasts: bool = False,
):
    colors = {"void": "#1f77b4", "cluster": "#d62728", "global": "#2ca02c", "supercluster": "#9467bd"}

    fig, ax = plt.subplots(figsize=(9.2, 6.1))

    # Observed points
    for p in points:
        if p.data_kind != "observed" or p.delta is None or p.H0 is None:
            continue
        color = colors.get(p.category, "#333333")
        marker = "D" if p.delta_kind == "delta_g_star_2mpp" else "o"
        alpha = 0.65 if p.delta_kind == "delta_g_star_2mpp" else 0.9
        ax.errorbar(
            p.delta,
            p.H0,
            xerr=p.delta_err if p.delta_err is not None else None,
            yerr=p.H0_err if p.H0_err is not None else None,
            fmt=marker,
            color=color,
            ecolor=color,
            capsize=3,
            alpha=alpha,
        )

    # Inferred points (optional; never part of the default fit unless explicitly enabled)
    for p in points:
        if p.data_kind != "inferred" or p.delta is None or p.H0 is None:
            continue
        color = colors.get(p.category, "#333333")
        ax.errorbar(
            p.delta,
            p.H0,
            xerr=p.delta_err if p.delta_err is not None else None,
            yerr=p.H0_err if p.H0_err is not None else None,
            fmt="^",
            mfc="none",
            mec=color,
            ecolor=color,
            capsize=3,
            alpha=0.55,
            label="inferred (not in default fit)",
        )

    # Optional overlay: overdense candidates with delta proxies (2M++)
    if overdense_candidates:
        xs = np.array([float(r["delta_proxy"]) for r in overdense_candidates], dtype=float)
        ys = np.array([float(r["H0"]) for r in overdense_candidates], dtype=float)
        ax.scatter(
            xs,
            ys,
            marker="D",
            s=38,
            facecolors="none",
            edgecolors="#444444",
            alpha=0.45,
            label="overdense candidates (delta_proxy, not fit)",
        )

    xs = _compute_plot_xs(points, n=220)

    H0_global = float(fit.get("H0_global", float("nan")))
    beta = float(fit.get("beta", float("nan")))
    ys = H0_global * (1.0 - beta * xs)

    # Uncertainty band (from MC propagation)
    if band is not None and band.get("p16") is not None:
        ax.fill_between(xs, band["p16"], band["p84"], color="#999999", alpha=0.25, label="MC 68% band")
        ax.plot(xs, band["p50"], color="#111111", linewidth=1.2, alpha=0.7, label="MC median")

    ax.plot(xs, ys, color="#111111", linewidth=1.8, label=f"fit (primary): H0_global={H0_global:.2f}, beta={beta:.3f}")

    # Forecast points: show where the *model* predicts at specified δ
    #
    # Note: forecast points are *not* fit points. They are displayed only when the
    # primary fit corresponds to the same δ-definition (typically direct δ).
    if plot_forecasts:
        forecast_label_added = False
        for p in points:
            if p.data_kind != "forecast" or p.delta is None:
                continue
            x0 = float(p.delta)
            label = "forecast (model @ δ; MC 68%)" if not forecast_label_added else "_nolegend_"
            if band is not None and band.get("p16") is not None and band.get("p50") is not None and band.get("p84") is not None:
                y16 = float(np.interp(x0, xs, np.asarray(band["p16"], dtype=float)))
                y50 = float(np.interp(x0, xs, np.asarray(band["p50"], dtype=float)))
                y84 = float(np.interp(x0, xs, np.asarray(band["p84"], dtype=float)))
                ax.errorbar(
                    x0,
                    y50,
                    yerr=[[y50 - y16], [y84 - y50]],
                    fmt="*",
                    color="#ff7f0e",
                    ecolor="#ff7f0e",
                    capsize=3,
                    alpha=0.70,
                    markersize=9,
                    label=label,
                )
            else:
                y0 = H0_global * (1.0 - beta * x0)
                ax.scatter(
                    [x0],
                    [y0],
                    marker="*",
                    s=90,
                    color="#ff7f0e",
                    alpha=0.70,
                    label=label,
                )
            forecast_label_added = True

    if _fit_is_valid_for_plot(fit_direct_inferred, min_points=4):
        H0_global_i = float(fit_direct_inferred.get("H0_global", float("nan")))
        beta_i = float(fit_direct_inferred.get("beta", float("nan")))
        ys_i = H0_global_i * (1.0 - beta_i * xs)
        ax.plot(
            xs,
            ys_i,
            color="#111111",
            linewidth=1.3,
            linestyle="--",
            alpha=0.8,
            label=f"fit (direct+inferred): H0_global={H0_global_i:.2f}, beta={beta_i:.3f}",
        )

    # Optional: fit only the Gavas-2024 inferred points (kept separate from direct fit)
    if _fit_is_valid_for_plot(fit_gavas, min_points=3):
        H0_global_g = float(fit_gavas.get("H0_global", float("nan")))
        beta_g = float(fit_gavas.get("beta", float("nan")))
        ys_g = H0_global_g * (1.0 - beta_g * xs)
        ax.plot(
            xs,
            ys_g,
            color="#8b0000",
            linewidth=1.5,
            linestyle=":",
            alpha=0.85,
            label=f"fit (Gavas sim): H0_global={H0_global_g:.2f}, beta={beta_g:.3f}",
        )

    if _fit_is_valid_for_plot(fit_proxy, min_points=3):
        H0_global_p = float(fit_proxy.get("H0_global", float("nan")))
        beta_p = float(fit_proxy.get("beta", float("nan")))
        ys_p = H0_global_p * (1.0 - beta_p * xs)
        ax.plot(xs, ys_p, color="#444444", linewidth=1.4, linestyle="--", label=f"fit (proxy): H0_global={H0_global_p:.2f}, beta={beta_p:.3f}")
    if band_proxy is not None and band_proxy.get("p16") is not None and _fit_is_valid_for_plot(fit_proxy, min_points=3):
        ax.fill_between(xs, band_proxy["p16"], band_proxy["p84"], color="#777777", alpha=0.18, label="proxy 68% band")
    ax.axvline(0.0, color="#999999", linestyle="--", linewidth=0.9)

    ax.set_xlabel("Density contrast δ")
    ax.set_ylabel("H0 (km/s/Mpc)")
    ax.set_title("Hubble tension: environment scan (v1)")
    ax.grid(alpha=0.2)

    # Legend cleanup
    handles, labels = ax.get_legend_handles_labels()
    seen = set()
    hh, ll = [], []
    for h, l in zip(handles, labels):
        if l not in seen:
            seen.add(l)
            hh.append(h)
            ll.append(l)
    if hh:
        ax.legend(hh, ll, fontsize=9)

    fig.tight_layout()
    fig.savefig(OUTPUT_PLOT, dpi=150)
    plt.close(fig)


def main():
    global OUTPUT_JSON, OUTPUT_PLOT
    ap = argparse.ArgumentParser(description="Hubble environment scan (v1): MC uncertainty bands.")
    ap.add_argument("--data-file", type=str, default=str(DATA_FILE))
    ap.add_argument("--output-json", type=str, default=str(OUTPUT_JSON))
    ap.add_argument("--output-plot", type=str, default=str(OUTPUT_PLOT))
    ap.add_argument("--n-mc", type=int, default=5000, help="Number of MC samples for uncertainty propagation")
    ap.add_argument("--seed", type=int, default=0, help="Random seed for MC")
    ap.add_argument(
        "--missing-h0-err",
        type=float,
        default=None,
        help="Impute missing H0_err in fit points (default: median of present H0_err)",
    )
    ap.add_argument("--missing-delta-err", type=float, default=0.0, help="Impute missing delta_err (default: 0)")
    ap.add_argument("--no-delta-err", action="store_true", help="Do not sample delta_err in MC (only sample H0_err)")
    ap.add_argument("--weight-mode", type=str, default="h0_err", choices=["h0_err", "equal"])
    ap.add_argument("--no-overdense-overlay", action="store_true", help="Do not overlay overdense candidate points")
    ap.add_argument("--fit-mode", type=str, default="dual", choices=["direct", "proxy", "dual"])
    ap.add_argument("--proxy-scale-min", type=float, default=50.0)
    ap.add_argument("--proxy-scale-max", type=float, default=100.0)
    ap.add_argument(
        "--include-inferred-fit",
        action="store_true",
        help="Also fit direct points with data_kind=inferred (must still have use_for_fit=true).",
    )
    ap.add_argument(
        "--include-gavas-fit",
        action="store_true",
        help="Also fit inferred points with delta_kind=gavas_delta_den_10mpch (kept separate from direct/proxy fits).",
    )
    args = ap.parse_args()

    # Allow custom outputs without changing code.
    OUTPUT_JSON = Path(str(args.output_json))
    OUTPUT_PLOT = Path(str(args.output_plot))

    data_file = Path(args.data_file)
    points = read_dataset(data_file)
    H0_anchor, anchor_id = get_global_anchor(points)

    fit_points_direct = _get_fit_points_direct(points)
    fit_points_direct_inferred = _get_fit_points_direct_including_inferred(points) if bool(args.include_inferred_fit) else []
    fit_points_gavas = _get_fit_points_gavas(points) if bool(args.include_gavas_fit) else []
    fit_points_proxy = _get_fit_points_proxy(
        points,
        scale_min=float(args.proxy_scale_min) if args.proxy_scale_min is not None else None,
        scale_max=float(args.proxy_scale_max) if args.proxy_scale_max is not None else None,
    )

    fit_direct: Dict[str, Any] = {}
    fit_fixed_direct: Dict[str, Any] = {}
    mc_direct: Dict[str, Any] = {}
    band_direct = None
    xs_plot = _compute_plot_xs(points, n=220)
    fit_direct_inferred: Dict[str, Any] = {}
    fit_gavas: Dict[str, Any] = {}

    if args.fit_mode in {"direct", "dual"}:
        fit_direct = fit_H0_global_and_beta(
            points,
            missing_H0_err=args.missing_h0_err,
            weight_mode=str(args.weight_mode),
            fit_points=fit_points_direct,
        )
        fit_fixed_direct = fit_beta_fixed_anchor(
            points,
            H0_anchor=float(H0_anchor),
            missing_H0_err=args.missing_h0_err,
            weight_mode=str(args.weight_mode),
            fit_points=fit_points_direct,
        )
        mc_direct = mc_uncertainty_propagation(
            points,
            n_samples=int(args.n_mc),
            seed=int(args.seed),
            missing_H0_err=args.missing_h0_err,
            missing_delta_err=float(args.missing_delta_err),
            include_delta_err=(not bool(args.no_delta_err)),
            weight_mode=str(args.weight_mode),
            fix_anchor=False,
            fit_points=fit_points_direct,
        )
        if fit_points_direct:
            a_samples = np.asarray(mc_direct.get("samples_for_plot", {}).get("H0_global", []), dtype=float)
            beta_samples = np.asarray(mc_direct.get("samples_for_plot", {}).get("beta", []), dtype=float)
            band_direct = predict_band(xs_plot, a_samples, beta_samples)

        if fit_points_direct_inferred:
            fit_direct_inferred = fit_H0_global_and_beta(
                points,
                missing_H0_err=args.missing_h0_err,
                weight_mode=str(args.weight_mode),
                fit_points=fit_points_direct_inferred,
            )

    if fit_points_gavas:
        fit_gavas = fit_H0_global_and_beta(
            points,
            missing_H0_err=args.missing_h0_err,
            weight_mode=str(args.weight_mode),
            fit_points=fit_points_gavas,
        )

    fit_proxy: Dict[str, Any] = {}
    fit_fixed_proxy: Dict[str, Any] = {}
    mc_proxy: Dict[str, Any] = {}
    band_proxy = None

    if args.fit_mode in {"proxy", "dual"}:
        fit_proxy = fit_H0_global_and_beta(
            points,
            missing_H0_err=args.missing_h0_err,
            weight_mode=str(args.weight_mode),
            fit_points=fit_points_proxy,
        )
        fit_fixed_proxy = fit_beta_fixed_anchor(
            points,
            H0_anchor=float(H0_anchor),
            missing_H0_err=args.missing_h0_err,
            weight_mode=str(args.weight_mode),
            fit_points=fit_points_proxy,
        )
        # MC bands for the proxy fit are only meaningful when the proxy dataset
        # has enough degrees of freedom (>=3 points).
        if len(fit_points_proxy) >= 3:
            mc_proxy = mc_uncertainty_propagation(
                points,
                n_samples=int(args.n_mc),
                seed=int(args.seed),
                missing_H0_err=args.missing_h0_err,
                missing_delta_err=float(args.missing_delta_err),
                include_delta_err=(not bool(args.no_delta_err)),
                weight_mode=str(args.weight_mode),
                fix_anchor=False,
                fit_points=fit_points_proxy,
            )
            a_samples = np.asarray(mc_proxy.get("samples_for_plot", {}).get("H0_global", []), dtype=float)
            beta_samples = np.asarray(mc_proxy.get("samples_for_plot", {}).get("beta", []), dtype=float)
            band_proxy = predict_band(xs_plot, a_samples, beta_samples)

    overdense = [] if bool(args.no_overdense_overlay) else read_overdense_candidates(OVERDENSE_CANDIDATES_FILE)

    # Plot
    plot_results(
        points,
        fit_direct if args.fit_mode in {"direct", "dual"} else fit_proxy,
        band=band_direct,
        overdense_candidates=overdense,
        fit_proxy=fit_proxy if args.fit_mode in {"proxy", "dual"} else None,
        band_proxy=band_proxy,
        fit_direct_inferred=fit_direct_inferred if fit_direct_inferred else None,
        fit_gavas=fit_gavas if fit_gavas else None,
        plot_forecasts=bool(args.fit_mode in {"direct", "dual"}),
    )

    # Predictive intervals for all points with delta
    def _predict_from_samples(a_samples: np.ndarray, beta_samples: np.ndarray) -> List[Dict[str, Any]]:
        preds: List[Dict[str, Any]] = []
        m = np.isfinite(a_samples) & np.isfinite(beta_samples)
        a_ok = a_samples[m]
        b_ok = beta_samples[m]
        for p in points:
            if p.delta is None:
                continue
            delta = float(p.delta)
            if a_ok.size:
                H0_draws = a_ok * (1.0 - b_ok * delta)
                q = _quantiles(H0_draws)
            else:
                q = {"p16": float("nan"), "p50": float("nan"), "p84": float("nan"), "mean": float("nan"), "std": float("nan")}
            preds.append({**p.to_dict(), "model_H0_mc": q})
        return preds

    preds_direct: List[Dict[str, Any]] = []
    preds_proxy: List[Dict[str, Any]] = []
    if args.fit_mode in {"direct", "dual"}:
        a_samples = np.asarray(mc_direct.get("samples_for_plot", {}).get("H0_global", []), dtype=float)
        beta_samples = np.asarray(mc_direct.get("samples_for_plot", {}).get("beta", []), dtype=float)
        preds_direct = _predict_from_samples(a_samples, beta_samples)
    if args.fit_mode in {"proxy", "dual"}:
        a_samples = np.asarray(mc_proxy.get("samples_for_plot", {}).get("H0_global", []), dtype=float)
        beta_samples = np.asarray(mc_proxy.get("samples_for_plot", {}).get("beta", []), dtype=float)
        preds_proxy = _predict_from_samples(a_samples, beta_samples)

    results: Dict[str, Any] = {
        "metadata": {
            "timestamp": datetime.now().isoformat(),
            "script": Path(__file__).name,
            "data_file": str(data_file),
            "anchor_env_id": anchor_id,
            "H0_anchor": float(H0_anchor),
            "overdense_candidates_file": str(OVERDENSE_CANDIDATES_FILE),
        },
        "config": {
            "data_file": str(data_file),
            "overdense_candidates_file": str(OVERDENSE_CANDIDATES_FILE),
            "n_mc": int(args.n_mc),
            "seed": int(args.seed),
            "missing_h0_err": float(args.missing_h0_err) if args.missing_h0_err is not None else None,
            "missing_delta_err": float(args.missing_delta_err),
            "include_delta_err": (not bool(args.no_delta_err)),
            "weight_mode": str(args.weight_mode),
            "no_overdense_overlay": bool(args.no_overdense_overlay),
            "fit_mode": str(args.fit_mode),
            "proxy_scale_min": float(args.proxy_scale_min),
            "proxy_scale_max": float(args.proxy_scale_max),
            "include_inferred_fit": bool(args.include_inferred_fit),
        },
        "metrics": {
            "fit_direct": {
                "H0_global": fit_direct.get("H0_global"),
                "beta": fit_direct.get("beta"),
                "beta_sigma": fit_direct.get("beta_sigma"),
                "chi2_dof": fit_direct.get("chi2_dof"),
                "n_points": fit_direct.get("n_points"),
            },
            "fit_direct_including_inferred": {
                "H0_global": fit_direct_inferred.get("H0_global") if fit_direct_inferred else None,
                "beta": fit_direct_inferred.get("beta") if fit_direct_inferred else None,
                "beta_sigma": fit_direct_inferred.get("beta_sigma") if fit_direct_inferred else None,
                "chi2_dof": fit_direct_inferred.get("chi2_dof") if fit_direct_inferred else None,
                "n_points": fit_direct_inferred.get("n_points") if fit_direct_inferred else 0,
            },
            "fit_gavas_inferred": {
                "H0_global": fit_gavas.get("H0_global") if fit_gavas else None,
                "beta": fit_gavas.get("beta") if fit_gavas else None,
                "beta_sigma": fit_gavas.get("beta_sigma") if fit_gavas else None,
                "chi2_dof": fit_gavas.get("chi2_dof") if fit_gavas else None,
                "n_points": fit_gavas.get("n_points") if fit_gavas else 0,
            },
            "fit_proxy": {
                "H0_global": fit_proxy.get("H0_global"),
                "beta": fit_proxy.get("beta"),
                "beta_sigma": fit_proxy.get("beta_sigma"),
                "chi2_dof": fit_proxy.get("chi2_dof"),
                "n_points": fit_proxy.get("n_points"),
            },
            "mc_beta_direct": mc_direct.get("summary", {}).get("beta"),
            "mc_beta_proxy": mc_proxy.get("summary", {}).get("beta"),
        },
        "fit": fit_direct if args.fit_mode in {"direct", "dual"} else fit_proxy,
        "fit_direct": fit_direct,
        "fit_with_proxy": fit_proxy,
        "fit_fixed_anchor": fit_fixed_direct if args.fit_mode in {"direct", "dual"} else fit_fixed_proxy,
        "fit_fixed_anchor_proxy": fit_fixed_proxy,
        "mc_uncertainty": {
            **{k: v for k, v in mc_direct.items() if k != "samples_for_plot"},
            "note": "MC samples are not stored in JSON (except summarized quantiles).",
        },
        "mc_uncertainty_proxy": {
            **{k: v for k, v in mc_proxy.items() if k != "samples_for_plot"},
            "note": "MC samples are not stored in JSON (except summarized quantiles).",
        },
        "predictions_mc_direct": preds_direct,
        "predictions_mc_proxy": preds_proxy,
        "predictions_mc": preds_direct if preds_direct else preds_proxy,
        "overdense_candidates_overlay_n": int(len(overdense)),
        "outputs": {"plot": str(OUTPUT_PLOT), "json": str(OUTPUT_JSON)},
    }

    with OUTPUT_JSON.open("w", encoding="utf-8") as f:
        json.dump(_json_sanitize(results), f, indent=2, ensure_ascii=False)

    print("Environment scan v1 completed.")
    print(f"Anchor: {anchor_id}, H0_anchor={H0_anchor}")
    if "H0_global" in fit_direct and "beta" in fit_direct:
        print(f"Fit (direct): H0_global={fit_direct['H0_global']:.4f}, beta={fit_direct['beta']:.4f}, n={fit_direct.get('n_points')}")
    if "H0_global" in fit_proxy and "beta" in fit_proxy:
        print(f"Fit (proxy): H0_global={fit_proxy['H0_global']:.4f}, beta={fit_proxy['beta']:.4f}, n={fit_proxy.get('n_points')}")
    if "beta" in fit_fixed_direct:
        print(f"Fit (fixed anchor, direct): beta={fit_fixed_direct['beta']:.4f}, n={fit_fixed_direct.get('n_points')}")
    if "beta" in fit_fixed_proxy:
        print(f"Fit (fixed anchor, proxy): beta={fit_fixed_proxy['beta']:.4f}, n={fit_fixed_proxy.get('n_points')}")
    if mc_direct.get("summary"):
        print(f"MC beta (direct): {mc_direct['summary'].get('beta')}")
    if mc_proxy.get("summary"):
        print(f"MC beta (proxy): {mc_proxy['summary'].get('beta')}")
    print(f"Output: {OUTPUT_JSON}")
    print(f"Plot: {OUTPUT_PLOT}")


if __name__ == "__main__":
    main()

