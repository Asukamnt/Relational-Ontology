# 复现：哈勃张力（环境修正）v0/v1

本页给出“最短复现路径”，用于对外分享（例如给 arXiv endorsers / 合作者）。

---

## 环境

- Python 3.10+（推荐 3.11）

安装依赖（在仓库根目录）：

```bash
pip install -r requirements.txt
```

---

## 复现 1：环境扫描（v1，含 MC 误差传播）

运行：

```bash
python math/hubble_env_scan/run_hubble_env_scan_v1.py
```

期望输出：

- `math/hubble_env_scan/results_v1.json`
- `math/hubble_env_scan/plot_v1.png`

关键数值（写在输出 JSON 里）：

- \(\beta\approx 0.171\)
- MC 16/50/84：0.139 / 0.169 / 0.201

---

## 复现 2：关系论一阶推导 toy（Poisson + entropic flow）

运行：

```bash
python math/spectral_graph/run_relational_hubble_beta_poisson_3d_v0.py
```

期望输出：

- `math/spectral_graph/relational_hubble_beta_poisson_3d_results_v0.json`
- `math/spectral_graph/relational_hubble_beta_poisson_3d_plot_v0.png`

在默认参数 `v_scale=1.0` 下，toy 会给出 \(\beta_{\mathrm{eff}}\approx 1/3\)，用于验证符号与线性关系；将其映射到 \(\beta\approx0.17\) 需引入动力学尺度（见论文/文档解释）。

---

## 文档入口

- 预印本草稿：`papers/hubble_tension_environment_v0.md`
- arXiv LaTeX 稿（英文，推荐用于投稿）：`papers/hubble_tension_environment_arxiv_v0.tex`
- arXiv 上传包生成器：`python papers/build_arxiv_hubble_tension_environment.py`
- 机制与证据链：`docs/physics/hubble_tension.md`
- 可证伪预测清单：`docs/verification/predictions.md`
- 管线说明：`math/hubble_env_scan/README.md`

