"""
Batch 2M++ density lookup for named targets (v0)
===============================================

Reads `twompp_targets_v0.csv` and produces:
- `twompp_targets_delta_v0.json`
- `twompp_targets_delta_v0.csv`

This is a helper to start "delta-ifying" concrete sky locations (observer positions)
under a common public density-field reference (2M++).
"""

from __future__ import annotations

import csv
import json
import argparse
from pathlib import Path
from typing import Dict, List, Optional

from twompp_density_lookup_v0 import _radec_to_galactic_deg, lookup_delta_g_star


TARGETS_CSV = Path(__file__).parent / "twompp_targets_v0.csv"


def _output_paths(smooth_sigma_mpc_h: float) -> tuple[Path, Path]:
    tag = f"{smooth_sigma_mpc_h:.3f}".replace(".", "p")
    out_json = Path(__file__).parent / f"twompp_targets_delta_sigma_{tag}_mpch.json"
    out_csv = Path(__file__).parent / f"twompp_targets_delta_sigma_{tag}_mpch.csv"
    return out_json, out_csv


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


def read_targets(path: Path) -> List[Dict]:
    rows: List[Dict] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(
                {
                    "target_id": (r.get("target_id") or "").strip(),
                    "target_name": (r.get("target_name") or "").strip(),
                    "coord_system": (r.get("coord_system") or "").strip().lower(),
                    "ra_deg": _parse_float(r.get("ra_deg", "")),
                    "dec_deg": _parse_float(r.get("dec_deg", "")),
                    "l_deg": _parse_float(r.get("l_deg", "")),
                    "b_deg": _parse_float(r.get("b_deg", "")),
                    "r_mpc_h": _parse_float(r.get("r_mpc_h", "")),
                    "reference": (r.get("reference") or "").strip(),
                    "notes": (r.get("notes") or "").strip(),
                }
            )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--smooth", type=float, default=4.0, help="Total Gaussian smoothing sigma [Mpc/h] (>=4)")
    args = ap.parse_args()

    targets = read_targets(TARGETS_CSV)
    out_rows: List[Dict] = []

    for t in targets:
        if t["r_mpc_h"] is None:
            raise ValueError(f"Missing r_mpc_h for target_id={t['target_id']}")

        coord = t["coord_system"]
        if coord == "galactic":
            if t["l_deg"] is None or t["b_deg"] is None:
                raise ValueError(f"Missing l/b for target_id={t['target_id']}")
            l_deg, b_deg = float(t["l_deg"]), float(t["b_deg"])
        elif coord == "equatorial":
            if t["ra_deg"] is None or t["dec_deg"] is None:
                raise ValueError(f"Missing ra/dec for target_id={t['target_id']}")
            l_deg, b_deg = _radec_to_galactic_deg(float(t["ra_deg"]), float(t["dec_deg"]))
        else:
            raise ValueError(f"Unknown coord_system='{coord}' for target_id={t['target_id']}")

        lookup = lookup_delta_g_star(
            l_deg=l_deg,
            b_deg=b_deg,
            r_mpc_h=float(t["r_mpc_h"]),
            smoothing_gaussian_mpc_h=float(args.smooth),
        )
        out_rows.append({**t, "galactic_l_deg": l_deg, "galactic_b_deg": b_deg, **lookup})

    out_json, out_csv = _output_paths(float(args.smooth))

    with out_json.open("w", encoding="utf-8") as f:
        json.dump(out_rows, f, indent=2, ensure_ascii=False)

    # Flat CSV export (subset)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "target_id",
            "target_name",
            "coord_system",
            "galactic_l_deg",
            "galactic_b_deg",
            "r_mpc_h",
            "x_mpc_h",
            "y_mpc_h",
            "z_mpc_h",
            "delta_g_star",
            "smoothing_gaussian_mpc_h",
            "reference",
            "notes",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in out_rows:
            w.writerow({k: r.get(k) for k in fieldnames})

    print(f"Wrote: {out_json}")
    print(f"Wrote: {out_csv}")


if __name__ == "__main__":
    main()

