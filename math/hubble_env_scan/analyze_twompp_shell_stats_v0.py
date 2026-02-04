"""
Analyze 2M++ shell statistics grid (v0).

Input:
  - twompp_shell_stats_grid_v0.csv

Outputs:
  - twompp_shell_stats_plot_v0.png
  - twompp_shell_stats_summary_v0.json

The goal is to make the "observer environment" proxy δ_g* visible as a function
of smoothing scale and shell definition, so the model v2 can explicitly carry
its delta_kind/scale dependence (rather than mixing δ definitions).
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
DEFAULT_INPUT = HERE / "twompp_shell_stats_grid_v0.csv"
DEFAULT_OUT_PLOT = HERE / "twompp_shell_stats_plot_v0.png"
DEFAULT_OUT_JSON = HERE / "twompp_shell_stats_summary_v0.json"


def _parse_float(v: str) -> Optional[float]:
    t = (v or "").strip()
    if not t:
        return None
    try:
        return float(t)
    except Exception:
        return None


@dataclass
class Row:
    shell_id: str
    rmin_mpc_h: float
    rmax_mpc_h: float
    smoothing_gaussian_mpc_h: float
    delta_g_star_mean: float
    delta_g_star_std: float
    n_cells: int
    origin_delta_g_star: float
    notes: str


def read_rows(path: Path) -> List[Row]:
    rows: List[Row] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            shell_id = (r.get("shell_id") or "").strip()
            rmin = _parse_float(r.get("rmin_mpc_h") or "")
            rmax = _parse_float(r.get("rmax_mpc_h") or "")
            sigma = _parse_float(r.get("smoothing_gaussian_mpc_h") or "")
            mean = _parse_float(r.get("delta_g_star_mean") or "")
            std = _parse_float(r.get("delta_g_star_std") or "")
            n_cells = _parse_float(r.get("n_cells") or "")
            origin = _parse_float(r.get("origin_delta_g_star") or "")
            if not shell_id or rmin is None or rmax is None or sigma is None or mean is None or std is None or n_cells is None or origin is None:
                continue
            rows.append(
                Row(
                    shell_id=shell_id,
                    rmin_mpc_h=float(rmin),
                    rmax_mpc_h=float(rmax),
                    smoothing_gaussian_mpc_h=float(sigma),
                    delta_g_star_mean=float(mean),
                    delta_g_star_std=float(std),
                    n_cells=int(n_cells),
                    origin_delta_g_star=float(origin),
                    notes=(r.get("notes") or "").strip(),
                )
            )
    return rows


def _group_by_shell(rows: List[Row]) -> Dict[str, List[Row]]:
    out: Dict[str, List[Row]] = {}
    for r in rows:
        out.setdefault(r.shell_id, []).append(r)
    for k in list(out.keys()):
        out[k] = sorted(out[k], key=lambda rr: rr.smoothing_gaussian_mpc_h)
    return out


def _select_shells(shell_ids: List[str], rows_by_shell: Dict[str, List[Row]]) -> List[str]:
    # Keep only those present; preserve user order.
    out = []
    seen = set()
    for sid in shell_ids:
        if sid in rows_by_shell and sid not in seen:
            out.append(sid)
            seen.add(sid)
    # Fallback: if user passed nothing valid, pick a small informative default.
    if not out:
        for sid in [
            "shell_0_67",
            "shell_30_67",
            "shell_0_75",
            "shell_30_75",
            "shell_0_150",
            "shell_40_150",
            "shell_40_200",
        ]:
            if sid in rows_by_shell and sid not in seen:
                out.append(sid)
                seen.add(sid)
    return out


def _nice_shell_label(rows: List[Row]) -> str:
    if not rows:
        return ""
    r0 = rows[0]
    return f"{r0.shell_id}  ({int(r0.rmin_mpc_h)}–{int(r0.rmax_mpc_h)} Mpc/h)"


def main() -> None:
    ap = argparse.ArgumentParser(description="Analyze 2M++ shell stats grid (v0).")
    ap.add_argument("--input", type=str, default=str(DEFAULT_INPUT))
    ap.add_argument("--output-plot", type=str, default=str(DEFAULT_OUT_PLOT))
    ap.add_argument("--output-json", type=str, default=str(DEFAULT_OUT_JSON))
    ap.add_argument(
        "--shells",
        type=str,
        default="shell_0_67,shell_30_67,shell_0_75,shell_30_75,shell_0_150,shell_40_150,shell_40_200",
        help="Comma-separated shell_ids to plot (default: common rmax=67/75/150 plus 40-200).",
    )
    args = ap.parse_args()

    in_path = Path(args.input)
    out_plot = Path(args.output_plot)
    out_json = Path(args.output_json)

    rows = read_rows(in_path)
    if not rows:
        raise SystemExit(f"No rows read from {in_path}")

    rows_by_shell = _group_by_shell(rows)
    shell_ids = [s.strip() for s in (args.shells or "").split(",") if s.strip()]
    sel = _select_shells(shell_ids, rows_by_shell)

    # Plot: δ_g* mean vs smoothing sigma for selected shells.
    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    for sid in sel:
        rr = rows_by_shell[sid]
        x = np.array([r.smoothing_gaussian_mpc_h for r in rr], dtype=float)
        y = np.array([r.delta_g_star_mean for r in rr], dtype=float)
        ax.plot(x, y, marker="o", linewidth=1.6, markersize=4.8, alpha=0.9, label=_nice_shell_label(rr))

    ax.axhline(0.0, color="#999999", linestyle="--", linewidth=0.9)
    ax.set_xlabel("Total Gaussian smoothing σ (Mpc/h)")
    ax.set_ylabel("Shell mean δ_g* (2M++)")
    ax.set_title("2M++ observer shell mean δ_g* vs smoothing (selected shells)")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8, ncol=1)
    fig.tight_layout()
    fig.savefig(out_plot, dpi=150)
    plt.close(fig)

    # JSON summary for docs (selected shells only)
    summary: Dict[str, Any] = {
        "metadata": {
            "timestamp": datetime.now().isoformat(),
            "script": Path(__file__).name,
            "input": str(in_path),
        },
        "selected_shells": sel,
        "shells": {},
        "outputs": {"plot": str(out_plot), "json": str(out_json)},
        "note": "delta_g* is luminosity-weighted galaxy density contrast from 2M++ (Carrick+2015).",
    }
    for sid in sel:
        rr = rows_by_shell[sid]
        summary["shells"][sid] = {
            "rmin_mpc_h": rr[0].rmin_mpc_h,
            "rmax_mpc_h": rr[0].rmax_mpc_h,
            "notes": rr[0].notes,
            "rows": [
                {
                    "smoothing_gaussian_mpc_h": r.smoothing_gaussian_mpc_h,
                    "delta_g_star_mean": r.delta_g_star_mean,
                    "delta_g_star_std": r.delta_g_star_std,
                    "n_cells": r.n_cells,
                    "origin_delta_g_star": r.origin_delta_g_star,
                }
                for r in rr
            ],
        }

    out_json.parent.mkdir(parents=True, exist_ok=True)
    with out_json.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    print("2M++ shell stats analysis completed.")
    print(f"Output plot: {out_plot}")
    print(f"Output JSON: {out_json}")


if __name__ == "__main__":
    main()

