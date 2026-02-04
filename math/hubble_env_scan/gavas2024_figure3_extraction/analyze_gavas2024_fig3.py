#!/usr/bin/env python3
"""
analyze_gavas2024_fig3.py

分析从 Gavas et al. 2024 Figure 3 提取的 δ_den vs δ_H 数据，
拟合 β 并与我们的 v0 结果对比。

用法：
    python analyze_gavas2024_fig3.py

依赖：
    numpy, scipy, matplotlib
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json
import csv
import sys

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


# ============================================================
# 常量
# ============================================================

SCRIPT_DIR = Path(__file__).parent
DATA_FILES = {
    "150": SCRIPT_DIR / "gavas2024_fig3_150mpc.csv",
    "500": SCRIPT_DIR / "gavas2024_fig3_500mpc.csv",
    "1000": SCRIPT_DIR / "gavas2024_fig3_1000mpc.csv",
}
OUTPUT_JSON = SCRIPT_DIR / "gavas2024_fig3_analysis.json"
OUTPUT_PLOT = SCRIPT_DIR / "gavas2024_fig3_fit.png"

# 我们 v0 拟合的 β 参考值
BETA_V0_FIT = 0.183  # 从 results_v0.json


# ============================================================
# 数据读取
# ============================================================

def read_extracted_data(filepath: Path) -> Optional[List[Dict]]:
    """
    读取从 WebPlotDigitizer 导出的 CSV 文件。
    
    期望列：delta_H, delta_den（或 X, Y）
    """
    if not filepath.exists():
        return None
    
    rows = []
    with open(filepath, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                # 支持多种列名格式
                delta_H = float(row.get("delta_H") or row.get("X") or row.get("x", 0))
                delta_den = float(row.get("delta_den") or row.get("Y") or row.get("y", 0))
                rows.append({"delta_H": delta_H, "delta_den": delta_den})
            except (ValueError, KeyError):
                continue
    return rows if rows else None


# ============================================================
# 拟合
# ============================================================

def fit_beta_from_gavas(
    delta_H_arr: np.ndarray,
    delta_den_arr: np.ndarray,
) -> Dict:
    """
    从 Gavas 的 δ_den vs δ_H 散点拟合 β。
    
    模型：δ_H = -β * δ_den + c
    （允许截距 c 来吸收系统偏移）
    
    返回：拟合结果字典
    """
    # 线性回归：delta_H = slope * delta_den + intercept
    slope, intercept, r_value, p_value, std_err = stats.linregress(delta_den_arr, delta_H_arr)
    
    # β = -slope（因为 δ_H = -β δ）
    beta = -slope
    beta_err = std_err  # slope 的标准误差
    
    # Pearson 和 Spearman 相关系数
    pearson_r, pearson_p = stats.pearsonr(delta_den_arr, delta_H_arr)
    spearman_r, spearman_p = stats.spearmanr(delta_den_arr, delta_H_arr)
    
    return {
        "n_points": len(delta_H_arr),
        "slope": float(slope),
        "intercept": float(intercept),
        "beta": float(beta),
        "beta_std_err": float(beta_err),
        "r_squared": float(r_value ** 2),
        "pearson_r": float(pearson_r),
        "pearson_p": float(pearson_p),
        "spearman_r": float(spearman_r),
        "spearman_p": float(spearman_p),
    }


# ============================================================
# 主流程
# ============================================================

def main():
    results = {
        "source": "Gavas et al. 2024 (arXiv:2407.10139) Figure 3",
        "model": "delta_H = -beta * delta_den + intercept",
        "reference_beta_v0": BETA_V0_FIT,
        "box_results": {},
        "combined": None,
    }
    
    all_delta_H = []
    all_delta_den = []
    
    print("=" * 60)
    print("Gavas et al. 2024 Figure 3 分析")
    print("=" * 60)
    
    # 逐盒子分析
    for box_size, filepath in DATA_FILES.items():
        print(f"\n--- {box_size} Mpc/h 盒子 ---")
        data = read_extracted_data(filepath)
        
        if data is None:
            print(f"  [跳过] 文件不存在或为空: {filepath}")
            results["box_results"][box_size] = {"status": "no_data", "file": str(filepath)}
            continue
        
        delta_H = np.array([d["delta_H"] for d in data])
        delta_den = np.array([d["delta_den"] for d in data])
        
        print(f"  数据点数: {len(data)}")
        print(f"  δ_H 范围: [{delta_H.min():.3f}, {delta_H.max():.3f}]")
        print(f"  δ_den 范围: [{delta_den.min():.3f}, {delta_den.max():.3f}]")
        
        fit_result = fit_beta_from_gavas(delta_H, delta_den)
        results["box_results"][box_size] = fit_result
        
        print(f"  拟合 β: {fit_result['beta']:.4f} ± {fit_result['beta_std_err']:.4f}")
        print(f"  Pearson r: {fit_result['pearson_r']:.3f} (p={fit_result['pearson_p']:.2e})")
        print(f"  Spearman r: {fit_result['spearman_r']:.3f}")
        
        # 与 v0 对比
        diff_pct = (fit_result["beta"] - BETA_V0_FIT) / BETA_V0_FIT * 100
        print(f"  与 v0 β 差异: {diff_pct:+.1f}%")
        
        all_delta_H.extend(delta_H)
        all_delta_den.extend(delta_den)
    
    # 合并分析
    if all_delta_H:
        print(f"\n--- 合并所有盒子 ---")
        all_delta_H = np.array(all_delta_H)
        all_delta_den = np.array(all_delta_den)
        
        combined_fit = fit_beta_from_gavas(all_delta_H, all_delta_den)
        results["combined"] = combined_fit
        
        print(f"  总数据点数: {combined_fit['n_points']}")
        print(f"  合并 β: {combined_fit['beta']:.4f} ± {combined_fit['beta_std_err']:.4f}")
        print(f"  Pearson r: {combined_fit['pearson_r']:.3f}")
        
        # 绘图
        plot_fit(all_delta_den, all_delta_H, combined_fit)
    
    # 保存结果
    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\n结果已保存到: {OUTPUT_JSON}")
    
    # 总结
    print("\n" + "=" * 60)
    print("总结")
    print("=" * 60)
    
    if results["combined"]:
        beta_gavas = results["combined"]["beta"]
        print(f"Gavas 模拟拟合 β: {beta_gavas:.4f}")
        print(f"我们 v0 拟合 β:   {BETA_V0_FIT:.4f}")
        print(f"差异: {(beta_gavas - BETA_V0_FIT) / BETA_V0_FIT * 100:+.1f}%")
        
        if abs(beta_gavas - BETA_V0_FIT) / BETA_V0_FIT < 0.5:
            print("\n✓ 模拟数据与我们的观测拟合在量级上一致！")
        else:
            print("\n⚠ 模拟数据与我们的观测拟合存在显著差异，需进一步调查。")
    else:
        print("⚠ 没有找到提取的数据文件。请先完成步骤 1-3。")


def plot_fit(delta_den: np.ndarray, delta_H: np.ndarray, fit_result: Dict):
    """绘制散点图和拟合线"""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # 散点
    ax.scatter(delta_den, delta_H, alpha=0.5, s=10, label="Gavas et al. data")
    
    # 拟合线
    x_fit = np.linspace(delta_den.min(), delta_den.max(), 100)
    y_fit = fit_result["slope"] * x_fit + fit_result["intercept"]
    ax.plot(x_fit, y_fit, "r-", lw=2, label=f"Fit: β = {fit_result['beta']:.3f}")
    
    # v0 参考线
    y_v0 = -BETA_V0_FIT * x_fit
    ax.plot(x_fit, y_v0, "g--", lw=2, alpha=0.7, label=f"v0 fit: β = {BETA_V0_FIT:.3f}")
    
    ax.axhline(0, color="gray", lw=0.5, ls=":")
    ax.axvline(0, color="gray", lw=0.5, ls=":")
    
    ax.set_xlabel(r"$\delta_{den}$ (local overdensity, 10 Mpc/h)", fontsize=12)
    ax.set_ylabel(r"$\delta_H = (H_L - H_G) / H_G$", fontsize=12)
    ax.set_title("Gavas et al. 2024 Figure 3: Simulation vs Our Fit", fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    fig.tight_layout()
    fig.savefig(OUTPUT_PLOT, dpi=150)
    print(f"图像已保存到: {OUTPUT_PLOT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
