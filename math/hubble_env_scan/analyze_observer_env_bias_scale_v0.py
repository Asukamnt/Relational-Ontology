"""
Analyze observer-environment bias vs scale from literature tables.

This script turns the curated CSV summaries into a quick "scale dependence" plot:
- Wojtak et al. 2014 Table 1 (H_loc bias for different observer selections)
- Odderskov et al. 2015 Table 2 (H_loc bias under different observer/sky setups)

All y-values are reported as mu_percent = 100 * (H_loc/H_bg - 1).
The sigma_percent columns are the reported scatter (not the error on the mean).
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).parent

DEFAULT_WOJTAK = HERE / "observer_env_bias_wojtak2014_table1.csv"
DEFAULT_ODDERSKOV = HERE / "observer_env_bias_odderskov2015_table2.csv"

DEFAULT_OUT_JSON = HERE / "observer_env_bias_scale_results_v0.json"
DEFAULT_OUT_PLOT = HERE / "observer_env_bias_scale_plot_v0.png"


def _parse_float(v: str) -> Optional[float]:
    t = (v or "").strip()
    if not t:
        return None
    try:
        return float(t)
    except Exception:
        return None


@dataclass
class BiasRow:
    label: str
    rmax_mpc_h: float
    mu_percent: float
    sigma_percent: float
    reference: str
    notes: str


def _read_csv_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _build_wojtak_curves(rows: List[Dict[str, str]]) -> Dict[str, List[BiasRow]]:
    """
    Focus on Wojtak+2014 Table 1 rows that cleanly express the 'observer environment' effect.

    We plot the CMB-frame measurements, with tracer haloes log10M>13, for:
      - random_in_space (baseline)
      - random_in_haloes (overdense observers)
      - void_centres (underdense observers)
    """
    wanted_schemes = {"random_in_space", "random_in_haloes", "void_centres"}
    out: Dict[str, List[BiasRow]] = {}
    for r in rows:
        scheme = (r.get("observer_scheme") or "").strip()
        frame = (r.get("reference_frame") or "").strip()
        tracer = (r.get("tracer_halo_mass_range") or "").strip()
        if scheme not in wanted_schemes:
            continue
        if frame != "CMB":
            continue
        if tracer != "log10Mhalo>13":
            continue
        rmax = _parse_float(r.get("rmax_mpc_h") or "")
        mu = _parse_float(r.get("mu_percent") or "")
        sig = _parse_float(r.get("sigma_percent") or "")
        if rmax is None or mu is None or sig is None:
            continue
        out.setdefault(scheme, []).append(
            BiasRow(
                label=scheme,
                rmax_mpc_h=float(rmax),
                mu_percent=float(mu),
                sigma_percent=float(sig),
                reference=(r.get("reference") or "").strip(),
                notes=(r.get("notes") or "").strip(),
            )
        )
    for k in list(out.keys()):
        out[k] = sorted(out[k], key=lambda rr: rr.rmax_mpc_h)
    return out


def _build_odderskov_curves(rows: List[Dict[str, str]]) -> Dict[str, List[BiasRow]]:
    """
    Select a small, high-signal subset of Odderskov+2015 Table 2:
      - A.0 Random positions in space
      - A.1 Random positions in halos
      - A.2 Local Group-like halos
    """
    wanted = {
        "A.0": "Random in space",
        "A.1": "Random in halos",
        "A.2": "Local-Group-like halos",
    }
    out: Dict[str, List[BiasRow]] = {}
    for r in rows:
        aid = (r.get("analysis_id") or "").strip()
        if aid not in wanted:
            continue
        rmax = _parse_float(r.get("rmax_mpc_h") or "")
        mu = _parse_float(r.get("mu_percent") or "")
        sig = _parse_float(r.get("sigma_percent") or "")
        if rmax is None or mu is None or sig is None:
            continue
        out.setdefault(aid, []).append(
            BiasRow(
                label=wanted[aid],
                rmax_mpc_h=float(rmax),
                mu_percent=float(mu),
                sigma_percent=float(sig),
                reference=(r.get("reference") or "").strip(),
                notes=(r.get("notes") or "").strip(),
            )
        )
    for k in list(out.keys()):
        out[k] = sorted(out[k], key=lambda rr: rr.rmax_mpc_h)
    return out


def _rows_to_jsonable(rows: List[BiasRow]) -> List[Dict[str, Any]]:
    return [
        {
            "label": r.label,
            "rmax_mpc_h": r.rmax_mpc_h,
            "mu_percent": r.mu_percent,
            "sigma_percent": r.sigma_percent,
            "reference": r.reference,
            "notes": r.notes,
        }
        for r in rows
    ]


def _plot_curves(ax: Any, curves: Dict[str, List[BiasRow]], *, title: str, labels_map: Optional[Dict[str, str]] = None) -> None:
    palette = {
        "random_in_space": "#444444",
        "random_in_haloes": "#d62728",
        "void_centres": "#1f77b4",
        "A.0": "#444444",
        "A.1": "#d62728",
        "A.2": "#ff7f0e",
    }
    markers = {
        "random_in_space": "o",
        "random_in_haloes": "s",
        "void_centres": "^",
        "A.0": "o",
        "A.1": "s",
        "A.2": "D",
    }
    for key, rows in curves.items():
        x = np.array([r.rmax_mpc_h for r in rows], dtype=float)
        y = np.array([r.mu_percent for r in rows], dtype=float)
        yerr = np.array([r.sigma_percent for r in rows], dtype=float)
        name = (labels_map or {}).get(key, rows[0].label if rows else key)
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            marker=markers.get(key, "o"),
            linewidth=1.4,
            markersize=5.5,
            capsize=3,
            color=palette.get(key, "#333333"),
            alpha=0.9,
            label=name,
        )
    ax.axhline(0.0, color="#999999", linestyle="--", linewidth=0.9)
    ax.set_title(title)
    ax.set_xlabel(r"$r_{\max}$ (Mpc/h)")
    ax.set_ylabel(r"$\mu$ [%]  (100×(H_loc/H_bg−1))")
    ax.grid(alpha=0.2)


def main() -> None:
    ap = argparse.ArgumentParser(description="Analyze observer-environment H_loc bias vs scale (v0).")
    ap.add_argument("--wojtak-csv", type=str, default=str(DEFAULT_WOJTAK))
    ap.add_argument("--odderskov-csv", type=str, default=str(DEFAULT_ODDERSKOV))
    ap.add_argument("--output-json", type=str, default=str(DEFAULT_OUT_JSON))
    ap.add_argument("--output-plot", type=str, default=str(DEFAULT_OUT_PLOT))
    args = ap.parse_args()

    wojtak_path = Path(args.wojtak_csv)
    odd_path = Path(args.odderskov_csv)
    out_json = Path(args.output_json)
    out_plot = Path(args.output_plot)

    wojtak_rows = _read_csv_rows(wojtak_path) if wojtak_path.exists() else []
    odd_rows = _read_csv_rows(odd_path) if odd_path.exists() else []

    wojtak_curves = _build_wojtak_curves(wojtak_rows)
    odd_curves = _build_odderskov_curves(odd_rows)

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6), sharey=True)
    _plot_curves(
        axes[0],
        wojtak_curves,
        title="Wojtak+2014 Table 1 (CMB frame; log10M>13)",
        labels_map={
            "random_in_space": "Random in space",
            "random_in_haloes": "Observers in haloes (overdense)",
            "void_centres": "Observers in void centres (underdense)",
        },
    )
    _plot_curves(
        axes[1],
        odd_curves,
        title="Odderskov+2015 Table 2 (selected analyses)",
        labels_map={
            "A.0": "Random in space",
            "A.1": "Random in haloes",
            "A.2": "Local-Group-like haloes",
        },
    )
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=8)
    fig.suptitle("Observer-environment bias vs scale (reported scatter as error bars)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_plot, dpi=150)
    plt.close(fig)

    results: Dict[str, Any] = {
        "metadata": {
            "timestamp": datetime.now().isoformat(),
            "script": Path(__file__).name,
            "wojtak_csv": str(wojtak_path),
            "odderskov_csv": str(odd_path),
        },
        "wojtak_curves": {k: _rows_to_jsonable(v) for k, v in wojtak_curves.items()},
        "odderskov_curves": {k: _rows_to_jsonable(v) for k, v in odd_curves.items()},
        "outputs": {"json": str(out_json), "plot": str(out_plot)},
        "note": "sigma_percent is the reported scatter of H_loc across observers, not the error on the mean.",
    }

    out_json.parent.mkdir(parents=True, exist_ok=True)
    with out_json.open("w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print("Observer-bias scale analysis completed.")
    print(f"Output JSON: {out_json}")
    print(f"Output plot: {out_plot}")


if __name__ == "__main__":
    main()

