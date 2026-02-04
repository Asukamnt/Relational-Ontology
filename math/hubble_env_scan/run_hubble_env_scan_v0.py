"""
Hubble tension environment scan (v0)
====================================

Goal:
- Fit (H0_global, beta) in H0(delta) = H0_global * (1 - beta * delta)
- Equivalent linear form: H0 = a - b * delta, where a=H0_global and b=H0_global*beta
- Produce JSON summary + plot for reproducibility

Run:
python run_hubble_env_scan_v0.py
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt


DATA_FILE = Path(__file__).parent / "hubble_env_data_v0.csv"
OUTPUT_JSON = Path(__file__).parent / "results_v0.json"
OUTPUT_PLOT = Path(__file__).parent / "plot_v0.png"
WOJTAK2014_TABLE1_FILE = Path(__file__).parent / "observer_env_bias_wojtak2014_table1.csv"
ODDERSKOV2015_TABLE2_FILE = Path(__file__).parent / "observer_env_bias_odderskov2015_table2.csv"

# Planck-2018-like baseline for simple theory comparisons (can be overridden later in v1).
OMEGA_M0_DEFAULT = 0.315


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

    def to_dict(self) -> Dict:
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
        }


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


def _parse_int(value: str) -> Optional[int]:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def read_dataset(path: Path) -> List[EnvPoint]:
    points: List[EnvPoint] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            points.append(
                EnvPoint(
                    env_id=row.get("env_id", "").strip(),
                    env_name=row.get("env_name", "").strip(),
                    category=row.get("category", "").strip(),
                    data_kind=row.get("data_kind", "").strip(),
                    use_for_fit=_parse_bool(row.get("use_for_fit", "")),
                    is_global_anchor=_parse_bool(row.get("is_global_anchor", "")),
                    delta=_parse_float(row.get("delta", "")),
                    delta_err=_parse_float(row.get("delta_err", "")),
                    H0=_parse_float(row.get("H0", "")),
                    H0_err=_parse_float(row.get("H0_err", "")),
                    redshift_range=row.get("redshift_range", "").strip(),
                    scale_mpc=row.get("scale_mpc", "").strip(),
                    source=row.get("source", "").strip(),
                    notes=row.get("notes", "").strip(),
                )
            )
    return points


def read_wojtak2014_table1(path: Path) -> List[Dict]:
    """Read a small, human-entered summary table (Wojtak+2014 Table 1).

    This file is used as an external, simulation-based consistency check for the
    *direction* and *scale dependence* of environment-driven Hloc biases.
    """
    if not path.exists():
        return []

    rows: List[Dict] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(
                {
                    "observer_scheme": (r.get("observer_scheme") or "").strip(),
                    "tracer_halo_mass_range": (r.get("tracer_halo_mass_range") or "").strip(),
                    "reference_frame": (r.get("reference_frame") or "").strip(),
                    "rmax_mpc_h": _parse_int(r.get("rmax_mpc_h", "")),
                    "mu_percent": _parse_float(r.get("mu_percent", "")),
                    "sigma_percent": _parse_float(r.get("sigma_percent", "")),
                    "reference": (r.get("reference") or "").strip(),
                    "notes": (r.get("notes") or "").strip(),
                }
            )
    return rows


def beta_linear_theory(Omega_m0: float = OMEGA_M0_DEFAULT, gamma: float = 0.55) -> Dict:
    """Very simple linear-theory baseline.

    In linear theory, velocity divergence relates to density contrast as
    ∇·v ~ -a H f δ with f ≈ Ω_m^gamma. For an isotropic local expansion rate,
    δH/H ~ -(f/3) δ (top-hat / large-scale limit), hence beta ≈ f/3.

    This is not a replacement for data fitting; it is a sanity-check baseline.
    """
    Omega_m0 = float(Omega_m0)
    gamma = float(gamma)
    f = Omega_m0 ** gamma if Omega_m0 > 0 else float("nan")
    beta = f / 3.0 if np.isfinite(f) else float("nan")
    return {
        "Omega_m0": Omega_m0,
        "gamma": gamma,
        "f_growth": f,
        "beta": beta,
        "derivation": "linear theory: deltaH/H ~ -(f/3) delta  => beta ~ f/3",
        "notes": "Order-of-magnitude baseline; depends on delta smoothing / estimator details.",
    }


def beta_turner1992_effective(Omega_m0: float = OMEGA_M0_DEFAULT) -> Dict:
    """Turner, Cen & Ostriker (1992) empirical 'local→global' correction baseline.

    They reported an approximate relation:
      (deltaH0/H0) ≈ -0.6 * delta_n_gal * Omega_m^0.4
    which is not directly identical to mass-density delta and is bias/definition dependent.

    If one *formally* compares to deltaH/H = -beta*delta, then an effective coefficient is
      beta_eff ≈ 0.6 * Omega_m^0.4    (for delta interpreted as galaxy-number overdensity).
    """
    Omega_m0 = float(Omega_m0)
    beta_eff = 0.6 * (Omega_m0 ** 0.4) if Omega_m0 > 0 else float("nan")
    return {
        "Omega_m0": Omega_m0,
        "beta_eff": beta_eff,
        "mapping": "deltaH/H = -beta*delta  ~  -0.6*delta_n_gal*Omega_m^0.4",
        "notes": "Applies to galaxy-number overdensity with bias/systematics; use as rough external reference only.",
        "reference": "Turner, Cen & Ostriker 1992 (AJ 103, 1427; DOI:10.1086/116156)",
    }


def read_odderskov2015_table2(path: Path) -> List[Dict]:
    """Read Odderskov et al. (2015; arXiv:1407.7364) Table 2 (long-form CSV)."""
    if not path.exists():
        return []

    rows: List[Dict] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(
                {
                    "analysis_id": (r.get("analysis_id") or "").strip(),
                    "analysis_desc": (r.get("analysis_desc") or "").strip(),
                    "rmax_mpc_h": _parse_int(r.get("rmax_mpc_h", "")),
                    "mu_percent": _parse_float(r.get("mu_percent", "")),
                    "sigma_percent": _parse_float(r.get("sigma_percent", "")),
                    "reference": (r.get("reference") or "").strip(),
                    "notes": (r.get("notes") or "").strip(),
                }
            )
    return rows


def _enrich_bias_rows_with_implied_delta(
    rows: List[Dict],
    beta_fit: Optional[float],
    beta_linear: Optional[float],
) -> List[Dict]:
    """Add derived fields to external bias tables.

    We map reported mean bias mu_percent (=100*(Hloc/H0-1)) to a fractional deviation:
      deltaH_over_H = mu_percent/100

    Under the model deltaH/H ≈ -beta*delta, we can define an *effective* delta:
      delta_eff = -(deltaH_over_H)/beta

    This is not a measurement of delta for those simulations; it is a convenient way to
    compare scales and signs across literature under a shared parametric form.
    """
    out: List[Dict] = []
    beta_fit_f = float(beta_fit) if beta_fit is not None and np.isfinite(beta_fit) else None
    beta_lin_f = (
        float(beta_linear) if beta_linear is not None and np.isfinite(beta_linear) else None
    )

    for r in rows:
        mu = r.get("mu_percent", None)
        mu_f = float(mu) if mu is not None and np.isfinite(mu) else None
        deltaH_over_H = (mu_f / 100.0) if mu_f is not None else None

        implied_fit = None
        if deltaH_over_H is not None and beta_fit_f not in (None, 0.0):
            implied_fit = -deltaH_over_H / beta_fit_f

        implied_lin = None
        if deltaH_over_H is not None and beta_lin_f not in (None, 0.0):
            implied_lin = -deltaH_over_H / beta_lin_f

        out.append(
            {
                **r,
                "deltaH_over_H": deltaH_over_H,
                "implied_delta_beta_fit": implied_fit,
                "implied_delta_beta_linear": implied_lin,
            }
        )

    return out


def get_global_anchor(points: List[EnvPoint]) -> Tuple[float, str]:
    anchors = [p for p in points if p.is_global_anchor and p.H0 is not None]
    if anchors:
        return anchors[0].H0, anchors[0].env_id
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


def fit_H0_global_and_beta(points: List[EnvPoint]) -> Dict:
    fit_points = _get_fit_points(points)

    if len(fit_points) < 2:
        return {
            "H0_global": float("nan"),
            "H0_global_sigma": None,
            "b": float("nan"),
            "b_sigma": None,
            "beta": float("nan"),
            "beta_sigma": None,
            "chi2": None,
            "dof": 0,
            "chi2_dof": None,
            "corr_delta_H0": None,
            "n_points": len(fit_points),
            "fit_points": fit_points,
            "weights": [],
            "residuals": [],
            "cov_ab": None,
        }

    deltas = np.array([p.delta for p in fit_points], dtype=float)
    H0s = np.array([p.H0 for p in fit_points], dtype=float)

    # Need at least two distinct delta values to identify both a and b in H0 = a - b*delta.
    if np.unique(deltas).size < 2:
        return {
            "H0_global": float("nan"),
            "H0_global_sigma": None,
            "b": float("nan"),
            "b_sigma": None,
            "beta": float("nan"),
            "beta_sigma": None,
            "chi2": None,
            "dof": 0,
            "chi2_dof": None,
            "corr_delta_H0": None,
            "n_points": len(fit_points),
            "fit_points": fit_points,
            "weights": [],
            "residuals": [],
            "cov_ab": None,
            "error": "insufficient_delta_diversity",
        }
    sigmas = np.array(
        [p.H0_err if p.H0_err is not None else 1.0 for p in fit_points],
        dtype=float,
    )
    weights = 1.0 / np.square(sigmas)

    # Weighted least squares for: H0 = a - b * delta
    # Design matrix: [1, -delta]
    X = np.column_stack([np.ones_like(deltas), -deltas])
    XtW = X.T * weights  # shape (2, n)
    XtWX = XtW @ X       # shape (2, 2)
    XtWy = XtW @ H0s     # shape (2,)

    try:
        theta = np.linalg.solve(XtWX, XtWy)
    except np.linalg.LinAlgError:
        theta = np.array([float("nan"), float("nan")], dtype=float)

    a = float(theta[0])
    b = float(theta[1])

    model = X @ theta
    residuals = H0s - model
    chi2 = float(np.sum(weights * np.square(residuals)))
    dof = max(len(fit_points) - 2, 0)
    chi2_dof = chi2 / dof if dof > 0 else None

    cov = None
    a_sigma = None
    b_sigma = None
    beta_sigma = None

    try:
        cov = np.linalg.inv(XtWX)
        if chi2_dof is not None:
            cov = cov * chi2_dof
        a_sigma = math.sqrt(float(cov[0, 0])) if float(cov[0, 0]) >= 0 else None
        b_sigma = math.sqrt(float(cov[1, 1])) if float(cov[1, 1]) >= 0 else None
    except np.linalg.LinAlgError:
        cov = None

    beta = (b / a) if np.isfinite(a) and a != 0 else float("nan")

    # Delta-method uncertainty for beta=b/a
    if cov is not None and np.isfinite(a) and a != 0 and np.isfinite(b):
        d_beta_da = -b / (a * a)
        d_beta_db = 1.0 / a
        var_beta = (
            d_beta_da * d_beta_da * float(cov[0, 0])
            + 2.0 * d_beta_da * d_beta_db * float(cov[0, 1])
            + d_beta_db * d_beta_db * float(cov[1, 1])
        )
        beta_sigma = math.sqrt(var_beta) if var_beta >= 0 else None

    corr = None
    if len(fit_points) >= 2:
        if np.std(deltas) > 0 and np.std(H0s) > 0:
            corr = float(np.corrcoef(deltas, H0s)[0, 1])

    return {
        "H0_global": a,
        "H0_global_sigma": a_sigma,
        "b": b,
        "b_sigma": b_sigma,
        "beta": beta,
        "beta_sigma": beta_sigma,
        "chi2": chi2,
        "dof": dof,
        "chi2_dof": chi2_dof,
        "corr_delta_H0": corr,
        "n_points": len(fit_points),
        "fit_points": fit_points,
        "weights": weights.tolist(),
        "residuals": residuals.tolist(),
        "cov_ab": cov.tolist() if cov is not None else None,
    }


def leave_one_out(points: List[EnvPoint]) -> List[Dict]:
    fit_points = _get_fit_points(points)
    results: List[Dict] = []
    for i, removed in enumerate(fit_points):
        subset = [p for j, p in enumerate(fit_points) if j != i]
        if len(subset) < 2:
            results.append(
                {
                    "removed_env_id": removed.env_id,
                    "H0_global": None,
                    "beta": None,
                    "chi2_dof": None,
                    "n_points": len(subset),
                }
            )
            continue
        sub_fit = fit_H0_global_and_beta(subset)
        results.append(
            {
                "removed_env_id": removed.env_id,
                "H0_global": sub_fit["H0_global"],
                "beta": sub_fit["beta"],
                "chi2_dof": sub_fit["chi2_dof"],
                "n_points": sub_fit["n_points"],
            }
        )
    return results


def build_results(points: List[EnvPoint], H0_global_anchor: float, anchor_id: str) -> Dict:
    fit = fit_H0_global_and_beta(points)
    H0_global = fit["H0_global"]
    beta_fit = fit["beta"]

    observed = []
    for p in points:
        if p.data_kind != "observed":
            continue
        if p.delta is None or p.H0 is None:
            continue
        model_H0 = H0_global * (1.0 - fit["beta"] * p.delta)
        observed.append(
            {
                **p.to_dict(),
                "model_H0": model_H0,
                "residual": p.H0 - model_H0,
            }
        )

    unplotted = []
    for p in points:
        if p.data_kind != "observed":
            continue
        if p.H0 is None:
            continue
        if p.delta is None:
            unplotted.append({**p.to_dict(), "reason": "missing delta"})

    forecasts = []
    for p in points:
        if p.data_kind != "forecast":
            continue
        if p.delta is None:
            continue
        model_H0 = H0_global * (1.0 - fit["beta"] * p.delta)
        forecasts.append({**p.to_dict(), "model_H0": model_H0})

    inferred = []
    for p in points:
        if p.data_kind != "inferred":
            continue
        if p.delta is None or p.H0 is None:
            continue
        model_H0 = H0_global * (1.0 - fit["beta"] * p.delta)
        inferred.append({**p.to_dict(), "model_H0": model_H0})

    # Baselines / external references (added to JSON for reproducibility)
    baseline_linear = beta_linear_theory()
    wojtak_rows = read_wojtak2014_table1(WOJTAK2014_TABLE1_FILE)
    odderskov_rows = read_odderskov2015_table2(ODDERSKOV2015_TABLE2_FILE)

    beta_lin = baseline_linear.get("beta")

    wojtak_rows_enriched = _enrich_bias_rows_with_implied_delta(
        wojtak_rows, beta_fit=beta_fit, beta_linear=beta_lin
    )
    odderskov_rows_enriched = _enrich_bias_rows_with_implied_delta(
        odderskov_rows, beta_fit=beta_fit, beta_linear=beta_lin
    )

    return {
        "metadata": {
            "data_file": str(DATA_FILE),
            "anchor_env_id": anchor_id,
            "external_files": {
                "wojtak2014_table1": str(WOJTAK2014_TABLE1_FILE),
                "odderskov2015_table2": str(ODDERSKOV2015_TABLE2_FILE),
            },
        },
        "fit": {
            "H0_global_anchor": H0_global_anchor,
            "H0_global": H0_global,
            "H0_global_sigma": fit["H0_global_sigma"],
            "b": fit["b"],
            "b_sigma": fit["b_sigma"],
            "beta": fit["beta"],
            "beta_sigma": fit["beta_sigma"],
            "chi2": fit["chi2"],
            "dof": fit["dof"],
            "chi2_dof": fit["chi2_dof"],
            "corr_delta_H0": fit["corr_delta_H0"],
            "n_points": fit["n_points"],
            "weights": "1/sigma^2 (missing sigma -> 1)",
            "model": "H0 = H0_global * (1 - beta * delta)  (fit via H0 = a - b*delta, beta=b/a)",
        },
        "baselines": {
            "linear_theory_beta": baseline_linear,
            "turner1992_beta_eff": beta_turner1992_effective(),
            "wojtak2014_table1": wojtak_rows_enriched,
            "odderskov2015_table2": odderskov_rows_enriched,
            "notes": {
                "implied_delta_fields": (
                    "External bias tables include deltaH_over_H (=mu_percent/100) and implied_delta "
                    "computed via deltaH/H ≈ -beta*delta, for beta=beta_fit and beta=beta_linear."
                )
            },
        },
        "leave_one_out": leave_one_out(points),
        "observed": observed,
    "unplotted": unplotted,
        "forecasts": forecasts,
        "inferred": inferred,
    }


def plot_results(points: List[EnvPoint], H0_global: float, beta: float):
    colors = {
        "void": "#1f77b4",
        "cluster": "#d62728",
        "global": "#2ca02c",
    }

    fig, ax = plt.subplots(figsize=(8, 5.5))

    # Observed points
    for p in points:
        if p.data_kind != "observed" or p.delta is None or p.H0 is None:
            continue
        color = colors.get(p.category, "#333333")
        ax.errorbar(
            p.delta,
            p.H0,
            xerr=p.delta_err if p.delta_err is not None else None,
            yerr=p.H0_err if p.H0_err is not None else None,
            fmt="o",
            color=color,
            ecolor=color,
            capsize=3,
            label=f"{p.category} (obs)" if p.category else "observed",
        )

    # Forecast points
    for p in points:
        if p.data_kind != "forecast" or p.delta is None or p.H0 is None:
            continue
        color = colors.get(p.category, "#333333")
        ax.scatter(
            p.delta,
            p.H0,
            marker="s",
            facecolors="none",
            edgecolors=color,
            label=f"{p.category} (forecast)" if p.category else "forecast",
        )

    # Inferred points (model-based or literature-derived)
    for p in points:
        if p.data_kind != "inferred" or p.delta is None or p.H0 is None:
            continue
        color = colors.get(p.category, "#333333")
        ax.scatter(
            p.delta,
            p.H0,
            marker="^",
            facecolors="none",
            edgecolors=color,
            label=f"{p.category} (inferred)" if p.category else "inferred",
        )

    deltas_all = [p.delta for p in points if p.delta is not None]
    if deltas_all:
        delta_min = min(deltas_all) - 0.1
        delta_max = max(deltas_all) + 0.1
    else:
        delta_min, delta_max = -0.8, 0.8

    xs = np.linspace(delta_min, delta_max, 200)
    ys = H0_global * (1.0 - beta * xs)
    ax.plot(xs, ys, color="#111111", linewidth=1.5, label=f"fit: H0_global={H0_global:.2f}, beta={beta:.3f}")

    ax.axvline(0.0, color="#999999", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Density contrast delta")
    ax.set_ylabel("H0 (km/s/Mpc)")
    ax.set_title("Hubble tension: environment scan (v0)")
    ax.grid(alpha=0.2)

    # Deduplicate legend entries
    handles, labels = ax.get_legend_handles_labels()
    seen = set()
    new_handles = []
    new_labels = []
    for h, l in zip(handles, labels):
        if l not in seen:
            seen.add(l)
            new_handles.append(h)
            new_labels.append(l)
    ax.legend(new_handles, new_labels, fontsize=8)

    fig.tight_layout()
    fig.savefig(OUTPUT_PLOT, dpi=150)
    plt.close(fig)


def main():
    points = read_dataset(DATA_FILE)
    H0_global_anchor, anchor_id = get_global_anchor(points)
    results = build_results(points, H0_global_anchor, anchor_id)

    def _json_sanitize(obj):
        # Convert NaN/Inf (not valid JSON) to None recursively.
        if isinstance(obj, (float, np.floating)):
            v = float(obj)
            return v if math.isfinite(v) else None
        if isinstance(obj, dict):
            return {k: _json_sanitize(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_json_sanitize(v) for v in obj]
        return obj

    with OUTPUT_JSON.open("w", encoding="utf-8") as f:
        json.dump(_json_sanitize(results), f, indent=2, ensure_ascii=False)

    plot_results(points, results["fit"]["H0_global"], results["fit"]["beta"])

    print("Environment scan completed.")
    print(f"Anchor (measured): {anchor_id}, H0_global={H0_global_anchor}")
    print(f"Fit: H0_global={results['fit']['H0_global']:.4f}, beta={results['fit']['beta']:.4f} (n={results['fit']['n_points']})")
    print(f"Output: {OUTPUT_JSON}")
    print(f"Plot: {OUTPUT_PLOT}")


if __name__ == "__main__":
    main()
