"""
Hubble environment dataset validator (v0)
=========================================

Goal
----
Before expanding the dataset, we want a repeatable "sanity report" that answers:

- Are required columns present?
- Are env_id values unique?
- How many observed points exist? How many are used for fit?
- Do fit points have sufficient delta diversity to identify (H0_global, beta)?
- What errors are missing (delta_err / H0_err)?

Run (repo root)
--------------
python math/hubble_env_scan/validate_hubble_env_data_v0.py

Outputs
-------
- dataset_summary_v0.json
"""

from __future__ import annotations

import json
import math
import csv
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

import run_hubble_env_scan_v1 as scan


DATA_FILE = Path(__file__).parent / "hubble_env_data_v0.csv"
OUTPUT_JSON = Path(__file__).parent / "dataset_summary_v0.json"


REQUIRED_COLUMNS = [
    "env_id",
    "env_name",
    "category",
    "data_kind",
    "use_for_fit",
    "is_global_anchor",
    "delta",
    "delta_err",
    "H0",
    "H0_err",
    "redshift_range",
    "scale_mpc",
    "source",
    "notes",
]

OPTIONAL_COLUMNS = [
    "delta_kind",
    "delta_scale_mpc_h",
    "delta_source",
]


def _json_sanitize(obj: Any):
    if isinstance(obj, (float, np.floating)):
        v = float(obj)
        return v if math.isfinite(v) else None
    if isinstance(obj, dict):
        return {k: _json_sanitize(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_json_sanitize(v) for v in obj]
    return obj


def main():
    import argparse

    ap = argparse.ArgumentParser(description="Validate hubble environment dataset")
    ap.add_argument("--data-file", type=str, default=str(DATA_FILE))
    ap.add_argument("--output-json", type=str, default=None)
    args = ap.parse_args()

    data_file = Path(args.data_file)
    if args.output_json:
        output_json = Path(args.output_json)
    else:
        suffix = "v1" if "v1" in data_file.stem else "v0"
        output_json = data_file.parent / f"dataset_summary_{suffix}.json"

    # Read raw CSV header to validate schema.
    with data_file.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []

    missing_cols = [c for c in REQUIRED_COLUMNS if c not in fieldnames]
    missing_optional = [c for c in OPTIONAL_COLUMNS if c not in fieldnames]
    extra_cols = [c for c in fieldnames if c not in REQUIRED_COLUMNS and c not in OPTIONAL_COLUMNS]

    # Parse points.
    points = scan.read_dataset(data_file)

    env_ids = [p.env_id for p in points if p.env_id]
    duplicates = [env_id for env_id, c in Counter(env_ids).items() if c > 1]

    anchors = [p for p in points if p.is_global_anchor and p.H0 is not None]
    anchor_ids = [p.env_id for p in anchors]

    observed = [p for p in points if p.data_kind == "observed"]
    fit_points = [p for p in observed if p.use_for_fit and p.delta is not None and p.H0 is not None]

    deltas = [p.delta for p in fit_points if p.delta is not None]
    unique_deltas = sorted(set(float(d) for d in deltas))

    missing = {
        "delta_err_missing": [p.env_id for p in fit_points if p.delta_err is None],
        "H0_err_missing": [p.env_id for p in fit_points if p.H0_err is None],
    }

    ok_for_two_param_fit = len(unique_deltas) >= 2 and len(fit_points) >= 2

    summary: Dict[str, Any] = {
        "data_file": str(data_file),
        "schema": {
            "required_columns": REQUIRED_COLUMNS,
            "optional_columns": OPTIONAL_COLUMNS,
            "present_columns": fieldnames,
            "missing_columns": missing_cols,
            "missing_optional_columns": missing_optional,
            "extra_columns": extra_cols,
        },
        "counts": {
            "n_rows_total": len(points),
            "n_observed": len(observed),
            "n_forecast": sum(1 for p in points if p.data_kind == "forecast"),
            "n_inferred": sum(1 for p in points if p.data_kind == "inferred"),
            "n_fit_points": len(fit_points),
        },
        "env_id": {
            "n_unique": len(set(env_ids)),
            "duplicates": duplicates,
        },
        "anchors": {
            "n_anchors": len(anchor_ids),
            "anchor_ids": anchor_ids,
        },
        "fit_identifiability": {
            "unique_deltas": unique_deltas,
            "n_unique_deltas": len(unique_deltas),
            "ok_for_two_param_fit": ok_for_two_param_fit,
            "note": "Need >=2 distinct delta values to identify (H0_global, beta) in H0=a-b*delta.",
        },
        "categories_observed": dict(Counter(p.category for p in observed)),
        "missing_errors_in_fit_points": missing,
        "delta_kind_counts": dict(Counter(p.delta_kind or "direct" for p in points if p.delta is not None)),
        "proxy_scales_mpc_h": sorted(
            set(float(p.delta_scale_mpc_h) for p in points if p.delta_scale_mpc_h is not None)
        ),
    }

    with output_json.open("w", encoding="utf-8") as f:
        json.dump(_json_sanitize(summary), f, indent=2, ensure_ascii=False)

    print("Dataset validation completed.")
    print(f"Output: {output_json}")
    print(f"Fit points: {len(fit_points)}, unique deltas: {len(unique_deltas)} -> ok_for_two_param_fit={ok_for_two_param_fit}")
    if duplicates:
        print("WARNING: duplicate env_id:", duplicates)
    if len(anchor_ids) != 1:
        print("WARNING: expected exactly 1 anchor; found:", anchor_ids)
    if missing_cols:
        print("WARNING: missing required columns:", missing_cols)
    if missing_optional:
        print("NOTE: missing optional columns:", missing_optional)


if __name__ == "__main__":
    main()

