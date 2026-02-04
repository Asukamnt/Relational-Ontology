# 哈勃张力：环境扫描管线（v0 / v1）

目标：把“密度环境影响局部 \(H_0\)”这件事，从单点吻合升级为**可复现的环境曲线拟合**。

核心模型（v0/v1）：

```
H0(delta) = H0_ref * (1 - beta * delta)
```

- 等价线性形式：`H0 = a - b * delta`，其中 `a=H0_ref`、`b=H0_ref*beta`
- v0 拟合采用**加权最小二乘**先拟合 `(a,b)`，再由 `beta=b/a` 得到 \(\beta\) 及其误差（delta method）

- `delta`：密度对比度 \(\delta=(\rho_{\text{local}}-\rho_{\text{mean}})/\rho_{\text{mean}}\)
- `H0_ref`：参考锚点（通常用 CMB / Planck 的推断值近似；在关系论口径下对应 \(\delta\approx 0\) 的“平均观察者参考值”，不是上帝视角）
- `beta`：密度-膨胀耦合系数（需要由多环境点集拟合）

---

## 文件

- `hubble_env_data_v0.csv`：v0 点集（观测点 + 预测点 + 文献推断点）
- `hubble_env_data_v1.csv`：v1 点集（新增代理 δ 字段与双轨拟合）
- `run_hubble_env_scan_v0.py`：读表 → 拟合 → 出图 → 导出 JSON
- `results_v0.json`：输出结果（运行脚本生成）
- `plot_v0.png`：输出图（运行脚本生成）
- `run_hubble_env_scan_v1.py`：v1：在 v0 基础上增加 **MC 误差传播**（δ_err/H0_err）与 **预测区间**，并可选 overlay `overdense_env_candidates_v0.csv`
- `results_v1.json`：v1 输出（包含拟合结果 + MC 分位数 + 每个点的预测区间）
- `plot_v1.png`：v1 输出图（含 68% MC band）
- `overdense_env_candidates_v0.csv`：P0 过密结构候选数据（目前多为“校准链/环境效应”层面的数；对部分条目补充了 `delta_proxy`，例如 Coma 的 2M++ \(\delta_g^\*\)@20 Mpc/h）
- `observer_env_bias_wojtak2014_table1.csv`：模拟中“过密观测者偏低 Hloc、空洞观测者偏高 Hloc”的定量摘要（Table 1）
- `gavas2024_figure3_extraction/`：从 Gavas et al. 2024 Figure 3 提取的趋势线数据及 beta 拟合（beta ~ 0.10，与 v0 同一量级）
- `observer_env_bias_odderskov2015_table2.csv`：模拟中“观测者选择/天空覆盖/光锥效应”对 \(H_{\mathrm{loc}}\) 偏差与方差的定量摘要（Table 2）
- `analyze_observer_env_bias_scale_v0.py`：把 Wojtak+2014 / Odderskov+2015 的摘要表整理成“偏差随尺度变化”的对照图（输出：`observer_env_bias_scale_results_v0.json`、`observer_env_bias_scale_plot_v0.png`）
- `analyze_twompp_shell_stats_v0.py`：把 2M++ 壳层均值 \(\delta_g^\*\)（观测者环境 proxy）随平滑尺度的变化画出来（输出：`twompp_shell_stats_summary_v0.json`、`twompp_shell_stats_plot_v0.png`）
- `analyze_twompp_radial_profile_v0.py`：计算 2M++ 在原点周围的径向壳层均值/球内累积均值 \(\delta_g^\*(r)\)，用于检查是否存在“类 KBC”尺度的系统性欠密（输出：`twompp_radial_profile_results_v0.json`、`twompp_radial_profile_plot_v0.png`）
- `run_twompp_observer_weighted_delta_v0.py`：把“距离梯子/样本几何”写成可计算的 \(w(r,\hat n)\)（径向窗口 + 天区权重 proxy），并用 2M++ 计算观测者有效 \(\delta_g^\*\)（输出：`twompp_observer_weighted_delta_results_v0.json`、`twompp_observer_weighted_delta_plot_v0.png`）
- `run_twompp_observer_weighted_h0_bias_v0.py`：把“隐藏关系/质量场响应”换成速度场 proxy：用 2M++ 速度场在同一 \(w(r,\hat n)\) 下估计窗口诱导的 H0 回归偏差（输出：`twompp_observer_weighted_h0_bias_results_v0.json`、`twompp_observer_weighted_h0_bias_plot_v0.png`）
- `run_hubble_env_v2_audit_v0.py`：模型 v2 的“口径审计 + scale→shell 映射”（输出：`hubble_env_v2_audit_results_v0.json`、`hubble_env_v2_audit_plot_v0.png`、`hubble_env_v2_audit_report_v0.md`）
- `run_cf2_cone_proxy_bins_v0.py`：在指定天空锥形区域（如 Shapley 方向）从 CF2 构造 \(\delta_g^\*\)–\(H_0\) 的 proxy 分箱点（输出：`cf2_cone_*_proxy_bins_results_v0.json`、`cf2_cone_*_proxy_bins_plot_v0.png`、`hubble_env_cf2_cone_*_proxy_bins_v0.csv`）
- `run_pantheonplus_env_correlation_v0.py`：Pantheon+ 低红移（距离模量 MU_SH0ES）× 2M++ \(\delta_g^\*\) 的环境相关检验（输出：`pantheonplus_env_correlation_results_v0.json`、`pantheonplus_env_correlation_plot_v0.png`）
- `run_twompp_overdense_observer_candidates_v0.py`：2M++：扫描“强过密的候选观测者位置”（proxy），并给出壳层均值 \(\delta_g^\*\)（输出：`twompp_overdense_observer_candidates_results_v0.json`、`twompp_overdense_observer_candidates_plot_v0.png`）

---

## 数据口径（CSV 字段）

| 字段 | 含义 |
|---|---|
| `env_id` | 行唯一标识 |
| `env_name` | 环境名称 |
| `category` | `void` / `cluster` / `global` 等 |
| `data_kind` | `observed` / `forecast` / `inferred` |
| `use_for_fit` | 是否参与拟合（预测点必须为 `false`） |
| `is_global_anchor` | 是否为参考锚点（历史字段名；语义为 `H0_ref`） |
| `delta`, `delta_err` | \(\delta\) 与误差（可缺省） |
| `H0`, `H0_err` | \(H_0\) 与误差（可缺省） |
| `redshift_range` | 红移范围或 `local` |
| `scale_mpc` | 平滑尺度（Mpc） |
| `source` | 文献/数据来源 |
| `notes` | 假设与口径说明 |
| `delta_kind` | `direct` / `delta_g_star_2mpp` |
| `delta_scale_mpc_h` | 代理 δ 的平滑尺度（Mpc/h） |
| `delta_source` | 代理 δ 来源文件 |

**v0 处理规则**：
- 缺失 `H0_err` 时，默认等权（weight=1）。
- `forecast` / `inferred` 行**不参与拟合**（只用于对照与展示）。
- 若拟合点集中 \(\delta\) 取值不足两类（例如全都相同），则无法同时识别 \(H0_\mathrm{global}\) 与 \(\beta\)，脚本会返回 `null` 并提示需要扩点。

**v1 新增规则（用于“更精确定量”）**：
- v1 会对拟合点做 **MC 误差传播**：默认同时采样 `H0_err` 与 `delta_err`，输出 \(\beta\) 的 16/50/84 分位数（近似 68% 区间）。
- 拟合点中缺失 `H0_err` 时，v1 默认用“已有误差的中位数”做插补（可用 `--missing-h0-err` 覆盖）。
- v1 另外提供“固定锚点”的单参数拟合：把 \(H0_\mathrm{global}\) 固定为 CMB 锚点，只拟合 \(\beta\)（见 `results_v1.json` 的 `fit_fixed_anchor`）。
- v1 默认会把 `overdense_env_candidates_v0.csv`（带 `delta_proxy` 的点）以**灰色菱形**叠加在图上，但不参与拟合（因为 `delta_proxy` 口径/尺度不一定与文献 δ 一致）。
- v1 支持**双轨拟合**：`direct`（文献 δ）与 `proxy`（2M++ δ_g*），用 `--fit-mode direct|proxy|dual` 选择，并可用 `--proxy-scale-min/--proxy-scale-max` 做尺度筛选。
  - `proxy` 轨只使用 `delta_kind=delta_g_star_2mpp` 且 `use_for_fit=true` 的行（避免“仅用于展示的 proxy 点”误入拟合）。
  - `proxy` 轨的 MC band/拟合线默认要求拟合点 **≥3**（否则自由度不足，容易产生误导性的带宽/外推）。
  - `--include-inferred-fit` 可选：把 `data_kind=inferred` 且 `use_for_fit=true` 的 direct 点纳入一个“对照拟合”（图中会额外画一条虚线）。默认主拟合仍只用 `observed` 点。
  - `--include-gavas-fit` 可选：对 `data_kind=inferred` 且 `delta_kind=gavas_delta_den_10mpch` 的点做**单独拟合**（图中会额外画一条红色点划线）。这组点来自 Gavas et al. 2024 Figure 3（10 Mpc/h 内过密度定义），**不建议与主线 direct δ 混合**，主要用于“方向/量级”外部参照。

---

## 运行方式

依赖（推荐在仓库根目录安装）：

```
pip install -r requirements.txt
```

运行（在仓库根目录执行）：

```
python math/hubble_env_scan/run_hubble_env_scan_v0.py
```

输出：
- `results_v0.json`
- `plot_v0.png`

v1（推荐：带误差条/预测区间）：

```
python math/hubble_env_scan/run_hubble_env_scan_v1.py
python math/hubble_env_scan/run_hubble_env_scan_v1.py --data-file math/hubble_env_scan/hubble_env_data_v1.csv --fit-mode dual --proxy-scale-min 50 --proxy-scale-max 100
python math/hubble_env_scan/run_hubble_env_scan_v1.py --data-file math/hubble_env_scan/hubble_env_data_v1.csv --fit-mode direct --include-inferred-fit
```

输出：
- `results_v1.json`
- `plot_v1.png`

模型 v2 口径审计（推荐在扩点/写文档前跑一次）：

```
python math/hubble_env_scan/run_hubble_env_v2_audit_v0.py
```

输出：
- `hubble_env_v2_audit_results_v0.json`
- `hubble_env_v2_audit_plot_v0.png`
- `hubble_env_v2_audit_report_v0.md`

其中 `results_v0.json` 里除了 `fit`（拟合结果）外，还包含 `baselines`：
- `baselines.linear_theory_beta`：线性理论的 \(\beta\approx f/3\) 基线（默认 \(\Omega_m=0.315\)）
- `baselines.turner1992_beta_eff`：Turner+1992 的经验修正系数（注意口径为星系数密度过/欠密，带偏置依赖）
- `baselines.wojtak2014_table1`：Wojtak+2014 Table 1 的“观测者环境→\(H_{\mathrm{loc}}\) 偏差”摘要表
- `baselines.odderskov2015_table2`：Odderskov+2015 Table 2 的“观测者选择/天空覆盖/光锥效应→\(H_{\mathrm{loc}}\) 偏差”摘要表

外部表行里还会附带三个派生字段（便于“观测者口径 δ 化”的快速对照）：
- `deltaH_over_H`：\(\delta H/H = \mu/100\)
- `implied_delta_beta_fit`：用本次拟合得到的 \(\beta\) 把 \(\delta H/H\) 映射成“有效 \(\delta\)”：\(\delta_{\mathrm{eff}}=-(\delta H/H)/\beta_{\mathrm{fit}}\)
- `implied_delta_beta_linear`：用线性理论 \(\beta_{\mathrm{lin}}\approx f/3\) 做同样映射

数据自检（推荐在扩点前先跑一遍）：

```
python math/hubble_env_scan/validate_hubble_env_data_v0.py
python math/hubble_env_scan/validate_hubble_env_data_v0.py --data-file math/hubble_env_scan/hubble_env_data_v1.csv
```

输出：`dataset_summary_v0.json`

---

## 观测者口径 δ（可选工具）

如果你想把“过密/欠密”从分类推进到**可重复计算的 \(\delta_{\mathrm{observer}}\)**，可以用 2M++ 重建密度场做一个公共参照（注意它输出的是**光度加权星系密度对比度** \(\delta_g^\*\)，不是严格的质量密度 \(\delta_m\)）：

- 脚本：`twompp_density_lookup_v0.py`
- 批量脚本：`twompp_density_batch_lookup_v0.py`（读取 `twompp_targets_v0.csv`）
- 壳层统计：`twompp_shell_stats_v0.py`（给定 `rmin/rmax` 直接算“观测者周围”的壳层平均 \(\delta_g^\*\)）
- 壳层批量：`twompp_shell_stats_batch_v0.py` + `twompp_shell_definitions_v0.csv`（输出 `twompp_shell_stats_grid_v0.csv/json`）
- 数据来源：`cosmicflows.iap.fr` 的 2M++ density field（会自动下载 `.npy` 到 `math/hubble_env_scan/_cache/`）
- 平滑尺度：高斯 \(4\,\mathrm{Mpc}/h\)

示例：

```
python math/hubble_env_scan/twompp_density_lookup_v0.py --ra 194.95 --dec 27.98 --r 70
```

批量运行（会生成 `twompp_targets_delta_sigma_*.json/csv`）：

```
python math/hubble_env_scan/twompp_density_batch_lookup_v0.py
```

如果要看更大尺度（更接近距离梯子/环境扫描口径）的 \(\delta\)，可增加总平滑尺度参数（单位 Mpc/h）：

```
python math/hubble_env_scan/twompp_density_batch_lookup_v0.py --smooth 50
python math/hubble_env_scan/twompp_density_batch_lookup_v0.py --smooth 100
```

如果想把“观测者环境”按距离范围压成一个标量（更贴近“distance ladder 使用的距离壳层”），可用壳层统计脚本：

```
python math/hubble_env_scan/twompp_shell_stats_v0.py --rmin 27 --rmax 200 --smooth 4
```

批量生成多个壳层与多个平滑尺度的表：

```
python math/hubble_env_scan/twompp_shell_stats_batch_v0.py --smooth 4,20,50,100
```

常用对照壳层（已写入 `twompp_shell_definitions_v0.csv`）包括：
- `0–67`, `30–67`：对齐 Odderskov+2015 常用的 \(r_{\max}=67\) 与 \(r_{\min}=30\)
- `0–75`, `30–75`：对齐 Wojtak+2014 的 \(r_{\max}=75\)
- `0–150`, `40–150`：对齐常用的 \(r_{\max}=150\)（并给出一个近似的 KBC-like 下限）

> 说明：`--r` 需要是 comoving Mpc/h；若想把 \(\delta_g^\*\) 近似映射到 \(\delta_m\)，可加 `--bias b`（\(\delta_m\approx\delta_g^\*/b\)），但偏置 \(b\) 本身依赖样本与尺度。

---

## 证伪条件（v0 口径）

- 若多环境拟合得到 \(\beta \le 0\) 且显著性 \(\ge 2\sigma\)，则密度-膨胀耦合方向需要被修正。
- 若在 \(\delta > 0\) 的过密环境中测得 \(H_0\) 系统性高于 \(H0_{\text{global}}\)（\(\ge 2\sigma\)），则该模型失效。
- 若残差呈显著非线性且一阶修正无法吸收，则需要升级模型（非线性或分尺度）。

---

## v0 局限

- 观测点集很小（基本只有“全局锚点 + 本地空洞”），当前只用于**流程验证**。
- 部分误差项缺失，只能等权处理。
- v0 的 `forecast` 行仍沿用旧文档中的 \(\beta=0.15\) 直觉值，仅作展示（历史遗留）。
- v1 的 `forecast` 行已对齐当前 direct 拟合（当前 \(\beta\approx 0.171\)），并在图中用星号标出“模型在指定 \(\delta\) 处的预测”（带 MC 68% 区间）。

---

## v1 下一步

- 扩展更多“明确给出 \(\delta\) 与 \(H_0\)”的环境点（并同步误差项）。
- 把 \(\delta\) 的定义按平滑尺度分层，做敏感性分析（\(\beta\) 随尺度变化）。
- 加入红移分箱：拟合 \(H_0(z,\delta)\) 而不是单一 \(H_0(\delta)\)。
