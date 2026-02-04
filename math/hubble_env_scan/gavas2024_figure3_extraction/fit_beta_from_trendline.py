#!/usr/bin/env python3
"""
fit_beta_from_trendline.py

Fit beta from manually extracted trendline points from Gavas et al. 2024 Figure 3.

Model: delta_H = -beta * delta_den + intercept
"""

import csv
import json
from pathlib import Path
from scipy import stats
import numpy as np

SCRIPT_DIR = Path(__file__).parent
OUR_BETA_V0 = 0.183
OUR_BETA_LINEAR = 0.177

def load_csv(filepath):
    with open(filepath, "r") as f:
        reader = csv.DictReader(f)
        return list(reader)

def fit_beta(delta_H, delta_den):
    """Fit delta_H = slope * delta_den + intercept, return beta = -slope"""
    slope, intercept, r_value, p_value, std_err = stats.linregress(delta_den, delta_H)
    return {
        "slope": slope,
        "intercept": intercept,
        "beta": -slope,
        "beta_std_err": std_err,
        "r_squared": r_value ** 2,
        "pearson_r": r_value,
    }

def main():
    files = {
        "150mpc": SCRIPT_DIR / "gavas2024_fig3_trendline_150mpc.csv",
        "500mpc": SCRIPT_DIR / "gavas2024_fig3_trendline_500mpc.csv",
        "1000mpc": SCRIPT_DIR / "gavas2024_fig3_trendline_1000mpc.csv",
    }
    
    results = {"boxes": {}, "combined": None}
    all_dH, all_dd = [], []
    
    print("=" * 60)
    print("Gavas et al. 2024 Trendline Beta Fit")
    print("=" * 60)
    
    for name, fpath in files.items():
        data = load_csv(fpath)
        dH = np.array([float(r["delta_H"]) for r in data])
        dd = np.array([float(r["delta_den"]) for r in data])
        
        fit = fit_beta(dH, dd)
        results["boxes"][name] = fit
        
        print(f"\n{name}:")
        print(f"  n_points: {len(dH)}")
        print(f"  beta = {fit['beta']:.4f} +/- {fit['beta_std_err']:.4f}")
        print(f"  R^2 = {fit['r_squared']:.3f}")
        
        all_dH.extend(dH)
        all_dd.extend(dd)
    
    # Combined fit
    all_dH = np.array(all_dH)
    all_dd = np.array(all_dd)
    combined = fit_beta(all_dH, all_dd)
    results["combined"] = combined
    
    print(f"\nCombined (all boxes):")
    print(f"  n_points: {len(all_dH)}")
    print(f"  beta = {combined['beta']:.4f} +/- {combined['beta_std_err']:.4f}")
    print(f"  R^2 = {combined['r_squared']:.3f}")
    
    # Comparison
    print("\n" + "=" * 60)
    print("Comparison with our model")
    print("=" * 60)
    print(f"  Gavas trendline beta:  {combined['beta']:.4f}")
    print(f"  Our v0 fitted beta:    {OUR_BETA_V0:.4f}")
    print(f"  Our linear theory:     {OUR_BETA_LINEAR:.4f}")
    
    diff_v0 = (combined["beta"] - OUR_BETA_V0) / OUR_BETA_V0 * 100
    diff_lin = (combined["beta"] - OUR_BETA_LINEAR) / OUR_BETA_LINEAR * 100
    print(f"\n  Difference from v0:    {diff_v0:+.1f}%")
    print(f"  Difference from linear:{diff_lin:+.1f}%")
    
    if abs(diff_v0) < 50:
        print("\n  [OK] Gavas simulation beta is within 50% of our v0 fit!")
    
    # Save
    output_path = SCRIPT_DIR / "gavas2024_trendline_fit_results.json"
    results["comparison"] = {
        "our_beta_v0": OUR_BETA_V0,
        "our_beta_linear": OUR_BETA_LINEAR,
        "diff_from_v0_pct": diff_v0,
        "diff_from_linear_pct": diff_lin,
    }
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {output_path}")

if __name__ == "__main__":
    main()
