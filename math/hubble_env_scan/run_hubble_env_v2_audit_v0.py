"""
Hubble tension (v2) dataset audit + scale mapping
=================================================

This script is a "v2 hardening" utility:

- Audits the v1 dataset (`hubble_env_data_v1.csv`) by track:
  - direct   : literature δ (primary fit)
  - sim      : simulation δ_den (kept separate; reference only)
  - proxy    : 2M++ δ_g* (kept separate; overlay / exploratory)

- Produces a *reproducible* report that records:
  - which points participate in which fit track
  - what δ definition/scale each point uses
  - suggested 2M++ "observer-shell" proxy mapping for points whose `scale_mpc`
    is expressed as an explicit distance range (e.g. "40-300").

Outputs (in this folder):
- hubble_env_v2_audit_results_v0.json
- hubble_env_v2_audit_plot_v0.png
- hubble_env_v2_audit_report_v0.md

Run (repo root):
    python math/hubble_env_scan/run_hubble_env_v2_audit_v0.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

import run_hubble_env_scan_v1 as scan


DATA_FILE = Path(__file__).parent / "hubble_env_data_v1.csv"
SHELL_DEF_FILE = Path(__file__).parent / "twompp_shell_definitions_v0.csv"
SHELL_STATS_FILE = Path(__file__).parent / "twompp_shell_stats_grid_v0.csv"

OUTPUT_JSON = Path(__file__).parent / "hubble_env_v2_audit_results_v0.json"
OUTPUT_PLOT = Path(__file__).parent / "hubble_env_v2_audit_plot_v0.png"
OUTPUT_MD = Path(__file__).parent / "hubble_env_v2_audit_report_v0.md"


@dataclass(frozen=True)
class ShellDef:
    shell_id: str
    rmin_mpc_h: float
    rmax_mpc_h: float
    notes: str


@dataclass(frozen=True)
class ShellStats:
    shell_id: str
    rmin_mpc_h: float
    rmax_mpc_h: float
    smoothing_gaussian_mpc_h: float
    delta_g_star_mean: float
    delta_g_star_std: float
    n_cells: int
    origin_delta_g_star: float
    reference: str
    notes: str


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


def _parse_float(s: str) -> Optional[float]:
    if s is None:
        return None
    t = str(s).strip()
    if not t:
        return None
    try:
        return float(t)
    except ValueError:
        return None


def read_shell_defs(path: Path) -> List[ShellDef]:
    shells: List[ShellDef] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            shell_id = (row.get("shell_id") or "").strip()
            rmin = _parse_float(row.get("rmin_mpc_h") or "")
            rmax = _parse_float(row.get("rmax_mpc_h") or "")
            if not shell_id or rmin is None or rmax is None:
                continue
            shells.append(
                ShellDef(
                    shell_id=shell_id,
                    rmin_mpc_h=float(rmin),
                    rmax_mpc_h=float(rmax),
                    notes=(row.get("notes") or "").strip(),
                )
            )
    return shells


def read_shell_stats(path: Path) -> List[ShellStats]:
    rows: List[ShellStats] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            shell_id = (r.get("shell_id") or "").strip()
            if not shell_id:
                continue
            rmin = _parse_float(r.get("rmin_mpc_h") or "")
            rmax = _parse_float(r.get("rmax_mpc_h") or "")
            smooth = _parse_float(r.get("smoothing_gaussian_mpc_h") or "")
            mean = _parse_float(r.get("delta_g_star_mean") or "")
            std = _parse_float(r.get("delta_g_star_std") or "")
            n_cells = _parse_float(r.get("n_cells") or "")
            origin = _parse_float(r.get("origin_delta_g_star") or "")
            if None in (rmin, rmax, smooth, mean, std, n_cells, origin):
                continue
            rows.append(
                ShellStats(
                    shell_id=shell_id,
                    rmin_mpc_h=float(rmin),
                    rmax_mpc_h=float(rmax),
                    smoothing_gaussian_mpc_h=float(smooth),
                    delta_g_star_mean=float(mean),
                    delta_g_star_std=float(std),
                    n_cells=int(float(n_cells)),
                    origin_delta_g_star=float(origin),
                    reference=(r.get("reference") or "").strip(),
                    notes=(r.get("notes") or "").strip(),
                )
            )
    return rows


def _track_from_point(p: scan.EnvPoint) -> str:
    if p.delta is None:
        return "none"
    if p.delta_kind is None or p.delta_kind == "direct":
        return "direct"
    if p.delta_kind == "delta_g_star_2mpp":
        return "proxy"
    if p.delta_kind == "gavas_delta_den_10mpch":
        return "sim"
    return "other"


def _role_from_point(p: scan.EnvPoint) -> str:
    track = _track_from_point(p)
    if p.data_kind == "forecast":
        return "forecast_marker"
    if p.data_kind == "observed" and p.use_for_fit and track == "direct":
        return "primary_fit"
    if p.data_kind == "inferred" and p.use_for_fit and track == "direct":
        return "optional_fit_direct"
    if p.data_kind == "inferred" and p.use_for_fit and track == "sim":
        return "reference_fit_sim"
    if p.data_kind == "observed" and track == "proxy":
        return "proxy_observed_overlay"
    return "overlay_only"


_RANGE_RE = re.compile(r"^\s*(?P<a>\d+(?:\.\d+)?)\s*-\s*(?P<b>\d+(?:\.\d+)?)\s*$")


def _suggest_shell_id_from_scale(scale_mpc: str, shells: List[ShellDef], *, rmax_cap: float = 200.0) -> Tuple[Optional[str], str]:
    """
    Suggest a 2M++ shell_id from a human scale string (e.g. "40-300").

    Returns (shell_id, note). If no suggestion is possible, returns (None, reason).

    Notes:
    - We interpret the input range as Mpc or Mpc/h only at the level of *coarse matching*.
    - We cap rmax to the maximal sphere fully inside the 2M++ cube (200 Mpc/h).
    """
    text = (scale_mpc or "").strip()
    if not text:
        return None, "scale_mpc_empty"
    if text.lower() == "global":
        return None, "scale_mpc_global"

    m = _RANGE_RE.match(text)
    if not m:
        return None, "scale_mpc_not_range"

    rmin = float(m.group("a"))
    rmax = float(m.group("b"))
    if rmax < rmin:
        rmin, rmax = rmax, rmin

    note_parts: List[str] = []
    if rmax > float(rmax_cap):
        note_parts.append(f"rmax_capped_to_{rmax_cap:g}")
        rmax = float(rmax_cap)

    # Prefer exact match first.
    for s in shells:
        if abs(s.rmin_mpc_h - rmin) < 1e-9 and abs(s.rmax_mpc_h - rmax) < 1e-9:
            return s.shell_id, (";".join(note_parts) if note_parts else "exact_match")

    # Otherwise, pick the closest by L1 distance in (rmin,rmax).
    best = None
    best_dist = float("inf")
    for s in shells:
        dist = abs(s.rmin_mpc_h - rmin) + abs(s.rmax_mpc_h - rmax)
        if dist < best_dist:
            best_dist = dist
            best = s
    if best is None:
        return None, "no_shell_defs"
    note_parts.append(f"approx_match_dist={best_dist:.3g}")
    return best.shell_id, ";".join(note_parts)


def _lookup_shell_stats(
    stats: Dict[Tuple[str, float], ShellStats],
    *,
    shell_id: str,
    smoothing: float,
) -> Optional[ShellStats]:
    # Use a small tolerance in case of float formatting differences.
    for (sid, sm), v in stats.items():
        if sid != shell_id:
            continue
        if abs(float(sm) - float(smoothing)) < 1e-9:
            return v
    return None


def build_markdown_report(
    *,
    points_audit: List[Dict[str, Any]],
    counts: Dict[str, Any],
    shell_notes: List[str],
) -> str:
    lines: List[str] = []
    lines.append("# 哈勃张力（模型 v2）点集口径审计报告（v0）")
    lines.append("")
    lines.append(f"_生成时间：{datetime.now().isoformat(timespec='seconds')}_")
    lines.append("")
    lines.append("## 摘要")
    lines.append("")
    lines.append(f"- 总行数：**{counts['n_total']}**")
    lines.append(f"- 观测：{counts['n_observed']}；预测：{counts['n_forecast']}；推断：{counts['n_inferred']}")
    lines.append(f"- 默认主线拟合点（direct+observed+use_for_fit）：**{counts['n_primary_fit']}**")
    lines.append("")
    lines.append("按 δ 口径轨道（track）计数：")
    for k, v in counts["by_track"].items():
        lines.append(f"- **{k}**：{v}")
    lines.append("")
    lines.append("按用途角色（role）计数：")
    for k, v in counts["by_role"].items():
        lines.append(f"- **{k}**：{v}")
    lines.append("")
    lines.append("## 关键口径提醒（v2）")
    lines.append("")
    lines.append("- **direct**：文献/模型给出的 δ（主线拟合用这个口径）。")
    lines.append("- **sim**：模拟定义的 \u03b4_den（例如 10 Mpc/h 内），只做外部参照；不要与 direct 混合拟合。")
    lines.append("- **proxy**：2M++ 的 \u03b4_g*（带平滑尺度），只做公共参照/overlay；不要与 direct 混合拟合。")
    lines.append("")
    lines.append("## 点集清单（按 role）")
    lines.append("")

    def _emit(role: str):
        rows = [p for p in points_audit if p.get("role") == role]
        if not rows:
            return
        lines.append(f"### {role}")
        lines.append("")
        for p in rows:
            env_id = p.get("env_id")
            delta = p.get("delta")
            H0 = p.get("H0")
            dk = p.get("delta_kind") or "direct"
            extra = []
            if p.get("delta_scale_mpc_h") is not None:
                extra.append(f"scale={p.get('delta_scale_mpc_h')} Mpc/h")
            if p.get("suggested_shell_id"):
                extra.append(f"shell={p.get('suggested_shell_id')}")
            lines.append(f"- `{env_id}`：δ={delta}, H0={H0}（delta_kind={dk}；" + ("；".join(extra) if extra else "无额外标注") + "）")
        lines.append("")

    for role in [
        "primary_fit",
        "optional_fit_direct",
        "reference_fit_sim",
        "proxy_observed_overlay",
        "forecast_marker",
        "overlay_only",
    ]:
        _emit(role)

    if shell_notes:
        lines.append("## 壳层映射说明")
        lines.append("")
        for s in shell_notes:
            lines.append(f"- {s}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main():
    ap = argparse.ArgumentParser(description="Hubble v2 audit: dataset tracks + shell mapping")
    ap.add_argument("--data-file", type=str, default=str(DATA_FILE))
    ap.add_argument("--shell-def-file", type=str, default=str(SHELL_DEF_FILE))
    ap.add_argument("--shell-stats-file", type=str, default=str(SHELL_STATS_FILE))
    ap.add_argument("--output-json", type=str, default=str(OUTPUT_JSON))
    ap.add_argument("--output-plot", type=str, default=str(OUTPUT_PLOT))
    ap.add_argument("--output-md", type=str, default=str(OUTPUT_MD))
    ap.add_argument("--rmax-cap", type=float, default=200.0, help="Cap rmax when mapping scale ranges to 2M++ shells")
    ap.add_argument(
        "--plot-smoothing",
        type=str,
        default="50,100",
        help="Comma-separated smoothing scales (Mpc/h) to plot for delta_observer shell means",
    )
    args = ap.parse_args()

    data_file = Path(args.data_file)
    shells = read_shell_defs(Path(args.shell_def_file))
    shell_stats_rows = read_shell_stats(Path(args.shell_stats_file))
    shell_stats: Dict[Tuple[str, float], ShellStats] = {(r.shell_id, float(r.smoothing_gaussian_mpc_h)): r for r in shell_stats_rows}

    points = scan.read_dataset(data_file)

    # Audit per point.
    points_audit: List[Dict[str, Any]] = []
    shell_notes: List[str] = []

    for p in points:
        track = _track_from_point(p)
        role = _role_from_point(p)
        suggested_shell_id, shell_note = _suggest_shell_id_from_scale(p.scale_mpc, shells, rmax_cap=float(args.rmax_cap))
        suggested_shell_stats: Dict[str, Any] = {}
        if suggested_shell_id:
            # Provide a small, fixed set of proxy lookups for readability.
            for smooth in (50.0, 100.0):
                ss = _lookup_shell_stats(shell_stats, shell_id=suggested_shell_id, smoothing=smooth)
                if ss is not None:
                    suggested_shell_stats[str(int(smooth))] = {
                        "smoothing_gaussian_mpc_h": ss.smoothing_gaussian_mpc_h,
                        "delta_g_star_mean": ss.delta_g_star_mean,
                        "delta_g_star_std": ss.delta_g_star_std,
                        "n_cells": ss.n_cells,
                    }
        if suggested_shell_id and shell_note:
            shell_notes.append(f"`{p.env_id}`：scale_mpc='{p.scale_mpc}' → `{suggested_shell_id}`（{shell_note}）")

        points_audit.append(
            {
                **p.to_dict(),
                "track": track,
                "role": role,
                "suggested_shell_id": suggested_shell_id,
                "suggested_shell_note": shell_note,
                "suggested_shell_stats": suggested_shell_stats,
            }
        )

    # Counts
    by_track: Dict[str, int] = {}
    by_role: Dict[str, int] = {}
    for p in points_audit:
        by_track[p["track"]] = by_track.get(p["track"], 0) + 1
        by_role[p["role"]] = by_role.get(p["role"], 0) + 1

    counts = {
        "n_total": len(points_audit),
        "n_observed": sum(1 for p in points_audit if p["data_kind"] == "observed"),
        "n_forecast": sum(1 for p in points_audit if p["data_kind"] == "forecast"),
        "n_inferred": sum(1 for p in points_audit if p["data_kind"] == "inferred"),
        "n_primary_fit": sum(1 for p in points_audit if p["role"] == "primary_fit"),
        "by_track": dict(sorted(by_track.items(), key=lambda kv: kv[0])),
        "by_role": dict(sorted(by_role.items(), key=lambda kv: kv[0])),
    }

    # Plot: delta_g* mean vs rmax for selected shells and smoothing scales.
    smoothing_list: List[float] = []
    for t in str(args.plot_smoothing).split(","):
        v = _parse_float(t)
        if v is not None:
            smoothing_list.append(float(v))
    smoothing_list = sorted(set(smoothing_list)) if smoothing_list else [50.0, 100.0]

    # Candidate shell families: rmin in {0, 30, 40} with various rmax
    families: Dict[str, List[ShellDef]] = {
        "rmin=0": [s for s in shells if abs(s.rmin_mpc_h - 0.0) < 1e-9],
        "rmin=30": [s for s in shells if abs(s.rmin_mpc_h - 30.0) < 1e-9],
        "rmin=40": [s for s in shells if abs(s.rmin_mpc_h - 40.0) < 1e-9],
    }

    fig, ax = plt.subplots(figsize=(9.4, 6.1))
    color_cycle = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    style_cycle = ["-", "--", ":", "-."]

    curve_data: Dict[str, Any] = {}
    for fi, (fam_label, fam_shells) in enumerate(families.items()):
        fam_shells = sorted(fam_shells, key=lambda s: s.rmax_mpc_h)
        for si, smooth in enumerate(smoothing_list):
            xs: List[float] = []
            ys: List[float] = []
            for sh in fam_shells:
                ss = _lookup_shell_stats(shell_stats, shell_id=sh.shell_id, smoothing=float(smooth))
                if ss is None:
                    continue
                xs.append(float(ss.rmax_mpc_h))
                ys.append(float(ss.delta_g_star_mean))
            if len(xs) < 2:
                continue
            label = f"{fam_label}, σ={smooth:g}"
            c = color_cycle[fi % len(color_cycle)]
            ls = style_cycle[si % len(style_cycle)]
            ax.plot(xs, ys, linestyle=ls, linewidth=1.8, color=c, alpha=0.9, label=label)
            curve_data[label] = {"rmax_mpc_h": xs, "delta_g_star_mean": ys}

    # Mark the commonly referenced shell used as a KBC-like proxy in docs.
    for mark_shell in ["shell_40_200", "shell_0_75", "shell_0_67", "shell_0_150"]:
        sd = next((s for s in shells if s.shell_id == mark_shell), None)
        if sd is None:
            continue
        ax.axvline(float(sd.rmax_mpc_h), color="#999999", linewidth=0.9, alpha=0.35)
        ax.text(
            float(sd.rmax_mpc_h) + 1.5,
            ax.get_ylim()[0] + 0.01,
            mark_shell.replace("shell_", ""),
            rotation=90,
            fontsize=8,
            color="#666666",
            va="bottom",
        )

    ax.set_xlabel("r_max (Mpc/h)  [2M++ shell outer radius]")
    ax.set_ylabel("mean δ_g* within shell")
    ax.set_title("2M++ observer-shell proxy: δ_g* mean vs r_max (v0)")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(Path(args.output_plot), dpi=150)
    plt.close(fig)

    # Save JSON + Markdown report.
    report_md = build_markdown_report(points_audit=points_audit, counts=counts, shell_notes=shell_notes)
    Path(args.output_md).write_text(report_md, encoding="utf-8")

    out: Dict[str, Any] = {
        "metadata": {
            "timestamp": datetime.now().isoformat(),
            "script": Path(__file__).name,
        },
        "inputs": {
            "data_file": str(data_file),
            "shell_def_file": str(Path(args.shell_def_file)),
            "shell_stats_file": str(Path(args.shell_stats_file)),
            "rmax_cap": float(args.rmax_cap),
            "plot_smoothing": smoothing_list,
        },
        "counts": counts,
        "points_audit": points_audit,
        "plot_data": curve_data,
        "outputs": {
            "json": str(Path(args.output_json)),
            "plot": str(Path(args.output_plot)),
            "report_md": str(Path(args.output_md)),
        },
        "notes": [
            "This audit does NOT refit beta; it only classifies points and records δ definitions/scale metadata.",
            "Suggested shell mapping is only a coarse proxy (2M++ limited to 200 Mpc/h). Do not treat it as a replacement for literature δ.",
        ],
    }

    Path(args.output_json).write_text(json.dumps(_json_sanitize(out), indent=2, ensure_ascii=False), encoding="utf-8")

    print("Hubble v2 audit completed.")
    print(f"Output JSON : {args.output_json}")
    print(f"Output plot : {args.output_plot}")
    print(f"Output MD   : {args.output_md}")


if __name__ == "__main__":
    main()

