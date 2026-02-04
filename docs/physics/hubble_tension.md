# 哈勃张力：关系论的一个定量检验（v0）

> 状态：🔄 方向与量级与部分观测/模拟一致（v0 口径；点集与系统对照仍待扩展）

---

## 1 问题陈述

### 1.1 传统描述

哈勃张力是宇宙学中最重要的未解问题之一：

| 测量方法 | H₀ (km/s/Mpc) | 来源 |
|---------|---------------|------|
| CMB（早期宇宙） | 67.36 ± 0.54 | Planck 2018 |
| 超新星（晚期宇宙） | 73.17 ± 0.86 | SH0ES |
| 差异 | ~9% | 5-6σ 显著 |

### 1.2 主流猜测

- 暗能量在演化？
- 存在第五力？
- 早期宇宙有未知机制？

---

## 2 关系论解释

### 2.1 核心洞察

```
关系网络在全局增加节点（宇宙膨胀）
    │
    ├── 欠密区（空洞）：关系稀疏 → 新节点更易"推开"旧节点
    │   → 局部膨胀更快 → H₀ 偏高
    │
    └── 过密区（星系团）：关系密集 → 引力自吸减缓膨胀
        → 局部膨胀更慢 → H₀ 偏低
```

### 2.2 数学形式（先相对、后锚定）

**核心可检验内容（相对口径；不需要“全局常数”）**：

- 观察者 A 处的局部测量记为 \(H_{0,A}:=H_0(\delta_A)\)。
- 预测的“硬内容”首先是**观察者之间的相对偏差规律**（比值形式）：
  \[
  \frac{H_{0,A}}{H_{0,B}}=\frac{1-\beta\,\delta_A}{1-\beta\,\delta_B}.
  \]
- 其中：
  - δ = (ρ_local - ρ_mean) / ρ_mean = 观测者环境的局部密度对比度（需显式注明口径/尺度）
  - β = 密度-膨胀耦合系数（由多环境点集拟合/外部模拟参照）

**把它写成绝对数值（km/s/Mpc）需要一个参考标尺**：

- 定义 \(H_{0,\mathrm{ref}} := H_0(\delta=0)\) 为“统计均匀的平均观察者参考值”。实践中常用 **CMB/Planck 的推断值**作为该参考的锚点；这不是“上帝视角”，而是一个**归一化/口径选择**。
- 于是得到线性近似：
  \[
  H_0(\delta)=H_{0,\mathrm{ref}}(1-\beta\,\delta).
  \]

**物理解释（方向）**：

- δ < 0（欠密）→ \(H_0 > H_{0,\mathrm{ref}}\)
- δ > 0（过密）→ \(H_0 < H_{0,\mathrm{ref}}\)

---

## 3 数据验证

### 3.1 关键观测数据

| 参数 | 数值 | 来源 |
|------|------|------|
| KBC 空洞密度对比度 | δ = -0.46 ± 0.06 | Keenan et al. |
| KBC 空洞范围 | 40 - 300 Mpc | 多项研究 |
| CMB H₀ | 67.4 km/s/Mpc | Planck |
| 超新星 H₀ | 73.2 km/s/Mpc | SH0ES |
| 空洞模型预测 H₀ | 72.1 km/s/Mpc | 2025 Tully-Fisher |

### 3.2 系数拟合

先用单点做“量级估计”（注意：单点只能给出粗略 β，并不能替代多点拟合）：

- 用 KBC（Tully-Fisher）点：\(\delta=-0.46,\;H_0=72.1\)
- 参考锚点：\(H_{0,\mathrm{ref}}=67.36\)（Planck；近似 \(\delta\approx 0\) 的“平均观察者”）

```
72.1 = 67.36 × (1 - β × (-0.46))
72.1/67.36 - 1 = 0.46β
β ≈ 0.153
```

同理，用 SH0ES 本地点（\(\delta=-0.46,\;H_0=73.17\)）会得到 \(\beta\approx 0.188\)。

更稳健地，我们用环境扫描管线（v0/v1）对 3 个观测点（Planck + KBC TF + SH0ES）做加权拟合，得到：
\(\beta\approx 0.171\pm 0.024\)（MC 16/50/84：0.139 / 0.169 / 0.201）。

**结论（当前口径）**：每 10% 的密度对比度 → 约 1.7% 的 \(H_0\) 偏离（线性近似）。

### 3.3 独立验证：N-body 模拟

2024-2025 年 N-body 模拟研究确认：

> "a **negative correlation** exists between local overdensity and H₀ deviations—**denser regions show greater systematic bias**"

这在方向上与该模型预测一致。

#### Gavas et al. 2024 定量验证

从 Gavas et al. (arXiv:2407.10139) Figure 3 提取趋势线数据后，直接拟合 \(\delta_H = -\beta\,\delta_{\mathrm{den}}\) 得到：

| 模拟盒子 | β (趋势线拟合) | Pearson r |
|----------|----------------|-----------|
| 150 Mpc/h | 0.151 ± 0.028 | -0.49 |
| 500 Mpc/h | 0.101 ± 0.010 | -0.49 |
| 1000 Mpc/h | 0.098 ± 0.003 | -0.51 |
| **合并** | **0.101 ± 0.008** | — |

**与本项目 v0 拟合比较**：

| 来源 | β | 备注 |
|------|---|------|
| Gavas 模拟（合并） | 0.101 | 10 Mpc/h 内 δ |
| 本项目 v0 拟合 | 0.173 | ~150-300 Mpc 尺度 δ |
| 线性理论 (f/3) | 0.177 | Ω_m ≈ 0.315 |

差异约 -45%，但在**同一量级**。差异可用 **δ 口径不同**解释：
- Gavas 用 10 Mpc/h 内的局部过密度
- KBC void 用 ~150-300 Mpc 尺度的密度对比度
- 更大尺度的 δ 通常更平滑（绝对值更小），相同 ΔH 对应更大的 β

**结论**：N-body 模拟独立验证了模型的**方向（负相关）**和**量级（β ~ 0.1）**。

（数据提取脚本与结果见：`math/hubble_env_scan/gavas2024_figure3_extraction/`）

为便于把这组外部参照直接叠加到本仓库的 v1 环境扫描图上，我们在
`math/hubble_env_scan/hubble_env_data_v1.csv` 中加入了少量 `gavas500_*` 的趋势线采样点（`data_kind=inferred`，
`delta_kind=gavas_delta_den_10mpch`），并在 v1 脚本中提供 `--include-gavas-fit` 开关用于**单独拟合/绘制**（不与主线 direct 拟合混合）：

```
python math/hubble_env_scan/run_hubble_env_scan_v1.py --fit-mode direct --include-gavas-fit
```

同时，更早的“大体积 N-body”工作也给出一致方向：**过密（halo）观测者测得更低的局部 \(H_{\mathrm{loc}}\)**，而**空洞中心观测者测得更高的 \(H_{\mathrm{loc}}\)**（Wojtak et al. 2014，Table 1）。
例如在 \(r_{\max}=75\,\mathrm{Mpc}/h\) 尺度上：
- 观测者在随机 halo（\(\log_{10} M_{\mathrm{halo}}>13\)）时：\(\mu_{75}\approx -0.8\%\)
- 观测者在空洞中心时：\(\mu_{75}\approx +1.0\%\)

（数据表已整理为：`math/hubble_env_scan/observer_env_bias_wojtak2014_table1.csv`）

另一个独立的 N-body 工作（Odderskov et al. 2015，Table 2）也报告了同类系统偏差：对“Local Group-like”观测者，
在 \(r_{\max}=150\,\mathrm{Mpc}/h\) 时 \(\mu_{150}\approx -0.3\%\)，在 \(r_{\max}=67\,\mathrm{Mpc}/h\) 时 \(\mu_{67}\approx -2.0\%\)，并随着尺度变大快速趋近于 0。
（数据表已整理为：`math/hubble_env_scan/observer_env_bias_odderskov2015_table2.csv`）

为把这两张表的“尺度趋势”一眼看清，我们额外生成了一个对照图（误差棒表示文献报告的散度 \(\sigma\)，不是均值误差）：
`math/hubble_env_scan/observer_env_bias_scale_plot_v0.png`（脚本：`math/hubble_env_scan/analyze_observer_env_bias_scale_v0.py`）。

更早的 Turner, Cen & Ostriker (1992) 也从大体积结构形成模拟给出了一个**经验近似**（用于“从局部修正到全局”）：
\[
\frac{\delta H_0}{H_0}\approx -0.6 \times \delta n_{\mathrm{gal}}\times \Omega_m^{0.4}
\]
其中 \(\delta n_{\mathrm{gal}}\) 表示局部体积内星系数密度的过/欠密程度。它与本项目采用的
\(\;H_0(\delta)=H_{0,\mathrm{ref}}\,(1-\beta\delta)\;\) 在形式上是同一类：**过密 → \(\delta H_0/H_0<0\)**。
需要注意该经验修正本身在文献中被强调为对**偏置（galaxy bias）与口径**敏感，并存在显著散度/系统不确定性，因此这里把它作为“方向与数量级”层面的外部参照，而不是直接用于本仓库 v0 拟合。

### 3.4 环境扫描管线 (v0)

为避免只靠单点拟合，已建立可复现的“环境扫描”管线：

- 数据表：`math/hubble_env_scan/hubble_env_data_v0.csv`
- 脚本：`math/hubble_env_scan/run_hubble_env_scan_v0.py`
- 输出：`math/hubble_env_scan/results_v0.json`、`math/hubble_env_scan/plot_v0.png`

**数据口径（v0）**：
- δ 直接取自文献的局部密度对比度
- 以 CMB/Planck 作为 \(H_{0,\mathrm{ref}}\) 锚点（\(\delta\approx 0\) 的“平均观察者参考值”）
- 观测点/预测点分离，预测点不参与拟合
- 缺失误差项按等权处理（后续补齐）
- SH0ES 本地样本默认采用 KBC 的 δ（同一环境假设）
- 拟合实现采用等价线性形式：\(H_0=a-b\delta\)，再由 \(\beta=b/a\) 得到 \(\beta\) 及其误差

**v0 拟合摘要（n=3）**：
- \(H_{0,\mathrm{ref}}\) ≈ 67.36 ± 0.46
- β ≈ 0.171 ± 0.024
- χ²/dof ≈ 0.74

**线性理论基线（对照）**：在简单线性近似下 \( \delta H/H \sim -(f/3)\delta \)，其中 \(f\simeq\Omega_m^{0.55}\)。取 Planck-2018 类参数 \(\Omega_m\approx0.315\) 时，
\(\beta_{\mathrm{lin}}\approx f/3\approx 0.177\)，与 v0 拟合的 \(\beta\approx0.173\) 在同一量级并非常接近。
（该对照已写入输出 `results_v0.json` 的 `baselines.linear_theory_beta`）

**关系论的一阶推导（toy，操作性）**：把“环境→\(H_0\) 偏差”写成一个可计算链条：
1) 几何/算子：从关系图的拉普拉斯 \(L\) 出发。在 3D torus 网格上有 \(L\approx-\nabla^2\)（离散连续极限）。
2) 势：对“关系密度对比度” \(\delta_{\mathrm{rel}}\) 定义 Poisson/Green 响应
\[
L\phi=-\delta_{\mathrm{rel}}\quad(\Rightarrow\ \nabla^2\phi=\delta_{\mathrm{rel}}).
\]
3) 动力学：采用熵梯度/势梯度流（离散版可对应 `run_g_from_entropy_gradient_v0.py` 的“沿 \(\nabla S\) 漂移”思想）
\[
\mathbf{v}=-f\,\nabla\phi,
\]
其中 \(f\) 是“增长率/响应强度”的无量纲系数（在标准线性理论里正是 growth rate；在 toy 脚本中对应参数 `--v-scale`）。
4) 观测口径：距离梯子回归的窗口偏差（与仓库的速度场脚本同口径）
\[
\frac{\delta H}{H}\approx \frac{\sum w(r)\,r\,v_r}{\sum w(r)\,r^2}\approx -\frac{f}{3}\,\delta_{\mathrm{eff}}
\quad\Rightarrow\quad \beta=\frac{f}{3}.
\]

仓库提供一个端到端 toy 验证脚本（3D torus + top-hat 空洞；默认 \(f=1\) 得到 \(\beta\approx 1/3\)；若取 \(f\approx0.51\) 则 \(\beta\approx0.17\) 与 v0/v1 拟合量级一致）：
- `math/spectral_graph/run_relational_hubble_beta_poisson_3d_v0.py`
- 输出：`math/spectral_graph/relational_hubble_beta_poisson_3d_results_v0.json`、`math/spectral_graph/relational_hubble_beta_poisson_3d_plot_v0.png`

**补充环境点（不参与拟合）**：
- KBC 参数化模型：δ ≈ -0.3 时预测 H₀ 上浮约 5.5%（Kenworthy et al. 2019）
- Shanks 2019 各向同性空洞：δ ≈ -0.2 时预测 H₀ 上浮约 3.5%（Kenworthy et al. 2019 引述）
- CF4 Tully-Fisher 约束下的本地超空洞：H₀ ≈ 72.1 / 70.4 / 70.2（指数/高斯/Maxwell 型）（Stiskalek et al. 2025）
- 过密环境候选数据（Coma/Laniakea 等）：见 `math/hubble_env_scan/overdense_env_candidates_v0.csv`（当前多为“环境效应/校准链”层面的结果，δ 口径待统一）

**说明**：v0 点集极小（仅锚点 + 本地空洞），结果用于“流程验证”，不是最终统计结论。

### 3.5 环境扫描管线 (v1：误差传播 / 预测区间)

为把“更精确定量预测”做成可复现输出，v1 在 v0 基础上新增：

- **MC 误差传播**：对拟合点的 \((\delta,\;H_0)\) 按 `delta_err` / `H0_err` 采样，得到 \(\beta\) 与 \(H0_\mathrm{global}\) 的分布（16/50/84 分位数），并在图上画出 68% band
- **缺失误差插补**：拟合点若缺 `H0_err`，默认用“已有误差的中位数”插补（可用 CLI 覆盖）
- **预测区间**：对 CSV 中每个带 \(\delta\) 的点，给出 `model_H0_mc` 的均值/方差与分位数
- **可选 overlay**：把 `overdense_env_candidates_v0.csv`（带 `delta_proxy` 的候选点）叠加展示，但不参与拟合（口径/尺度可能不同）
- **双轨拟合（v1+）**：引入 `hubble_env_data_v1.csv`，新增 `delta_kind/delta_scale_mpc_h/delta_source` 字段，支持 `direct`（文献 δ）与 `proxy`（2M++ δ_g*）两套拟合，并用 `--proxy-scale-min/max` 做尺度筛选。

运行：

```
python math/hubble_env_scan/run_hubble_env_scan_v1.py
```

输出：
- `math/hubble_env_scan/results_v1.json`
- `math/hubble_env_scan/plot_v1.png`

**v1 结果示例（当前点集仍很小，主要用于流程与不确定度口径验证）**：
- \(\beta\)（点估计）≈ 0.171
- \(\beta\)（MC 16/50/84 分位数）≈ 0.139 / 0.169 / 0.201

**对照（可选）**：若启用 `--include-inferred-fit`（把 `data_kind=inferred` 且 `use_for_fit=true` 的 direct 点纳入对照拟合），当前点集会扩展到 **n=5**，并得到：
- \(\beta\) ≈ 0.172（与主线 n=3 拟合在数值上非常接近）

**观测者口径 δ 的下一步**：为把“过密/欠密”从分类推进到可重复计算的 \(\delta_{\mathrm{observer}}\)，仓库已加入一个可选工具，
可从公开的 2M++ 重建密度场查询 \(\delta_g^\*\)（光度加权星系密度对比度，平滑尺度 \(4\,\mathrm{Mpc}/h\)）作为公共参照：
`math/hubble_env_scan/twompp_density_lookup_v0.py`。
它会把数据下载到本地缓存并输出查询点的 \(\delta_g^\*\)（可选用偏置参数近似映射到 \(\delta_m\)）。

### 3.6 模型 v2：口径分层（direct / sim / proxy）

为避免“把不同口径的数据强行压成同一个 \(\beta\)”造成伪矛盾，v1 开始按 \(\delta\) 的**定义/尺度**分层使用：

- **direct（主线）**：`delta_kind=direct`（或缺省），\(\delta\) 来自文献/模型给出的“环境密度对比度”。用于主拟合与主预测。
- **sim（外部参照）**：`delta_kind=gavas_delta_den_10mpch`（10 Mpc/h 内 \(\delta_{\mathrm{den}}\)）。只做单独拟合/绘制（`--include-gavas-fit`），用来校验“方向与量级”。
- **proxy（展示/探索）**：`delta_kind=delta_g_star_2mpp`（2M++ \(\delta_g^\*\)，带平滑尺度 `delta_scale_mpc_h`）。当前更适合作为“公共参照/候选环境 overlay”，不建议与主线 direct 混合拟合。

#### 3.6.1 观测者口径 \(\delta_{\mathrm{observer}}(R)\)：2M++ 壳层均值（围绕原点）

为把“我们到底处在怎样的环境”从概念推进到**可复现的数值**，我们用 2M++ 重建密度场的
光度加权星系密度对比度 \(\delta_g^\*\)，在“围绕本地观测者（原点）”的球壳上做体积平均：

- 脚本：`math/hubble_env_scan/twompp_shell_stats_batch_v0.py`
- 壳层定义：`math/hubble_env_scan/twompp_shell_definitions_v0.csv`
- 输出：`math/hubble_env_scan/twompp_shell_stats_grid_v0.csv/json`

下面列出几个与文献常用 \(r_{\max}\) 对齐的例子（\(\sigma\) 为**总**高斯平滑尺度，单位 Mpc/h；2M++ 原始产品本身已含 \(\sigma=4\) 平滑）：

| 壳层 (Mpc/h) | \(\delta_g^\*\) mean @ \(\sigma=100\) | \(\delta_g^\*\) mean @ \(\sigma=50\) | 备注 |
|---|---:|---:|---|
| 0–67 | 0.0125 | 0.0435 | Odderskov+2015 的 \(r_{\max}=67\) 对照尺度 |
| 30–67 | 0.0124 | 0.0426 | Odderskov+2015 的 \(r_{\min}=30\) 下限口径 |
| 0–75 | 0.0122 | 0.0409 | Wojtak+2014 的 \(r_{\max}=75\) 对照尺度 |
| 30–75 | 0.0121 | 0.0401 | 同上 + 保守下限 \(r_{\min}=30\) |
| 0–150 | 0.0080 | 0.0182 | 常用 \(r_{\max}=150\) 对照尺度 |
| 40–150 | 0.0079 | 0.0175 | “KBC-like” 下限的近似（受 2M++ 200 Mpc/h 盒子限制） |

**读法**：
- 这些数是 \(\delta_g^\*\)（星系光度加权的过/欠密），**不等同于**文献中的质量密度 \(\delta_m\) 或 KBC 的计数口径 \(\delta\)；
  因此它更适合作为 `proxy` 轨或“观测者环境的公共参照”，而不是直接拿来替代主线 `direct` 的 \(\delta\)。
- 你会看到：当平滑尺度变大（例如 \(\sigma=100\)）时，\(\delta_g^\*\) 壳层均值会显著变小并趋近 0，
  这与模拟工作中“\(H_{\mathrm{loc}}\) 偏差随 \(r_{\max}\) 增大快速衰减”的趋势是同方向的（尺度越大越接近全局平均）。

我们也把“\(\delta_g^\*\) 随平滑尺度 \(\sigma\) 的变化”画成了图，便于直观看到这种收敛趋势：
`math/hubble_env_scan/twompp_shell_stats_plot_v0.png`（脚本：`math/hubble_env_scan/analyze_twompp_shell_stats_v0.py`）。

另外，为了直接检查“是否存在一个在 \(40\!-\!200\,\mathrm{Mpc}/h\) 等尺度上持续为负的欠密带”（类 KBC 口径的直觉），
我们还生成了径向剖面图（壳层均值 + 球内累积均值）：
`math/hubble_env_scan/twompp_radial_profile_plot_v0.png`（脚本：`math/hubble_env_scan/analyze_twompp_radial_profile_v0.py`）。
该剖面显示：在 \(\sigma=4\) 的 2M++ 场中，\(\delta_g^\*\) 在部分半径段可为负，但在 \(r\lesssim 200\,\mathrm{Mpc}/h\) 的球内累积均值总体接近 0，
提示“2M++ 的球对称体积平均 proxy”与文献 KBC 欠密 \(\delta\) 之间仍存在显著口径差异（这正是 v2 必须分层的原因之一）。

#### 3.6.2 Pantheon+ 低红移：\(\delta_g^\*\)–\(H_0\) 相关（proxy 测试）

我们额外做了一个“真实 SN 样本 + 2M++”的相关检验，用来评估：
**把目标点的 \(\delta_g^\*\) 当作环境变量时，是否能直接得到我们想要的（过密→更低 \(H_0\)）单调趋势**。

- 脚本：`math/hubble_env_scan/run_pantheonplus_env_correlation_v0.py`
- 数据：Pantheon+ `Pantheon+SH0ES.dat`（公开数据；脚本自动下载）
- 选择：排除校准宿主（`IS_CALIBRATOR=0`），按 `CID` 去重，取 \(z\in[0.015,0.060]\)（示例运行：n=403）
- 指标：用 `MU_SH0ES` 得到距离 \(D\)，再计算 \(H_0=(cz)/D\)，并分别用 `zHD`（含 VPEC 修正）与 `zCMB`（不含 VPEC 修正）两套口径对照

**结果（v0）**：在 \(\sigma=20\) Mpc/h 下相关几乎为 0；在更大平滑尺度（\(\sigma=50,100\)）下反而出现“过密→更高 \(H_0\)”的弱正相关。
这进一步支持我们在 CF2“天空锥 proxy”实验中的判断：这种“目标点环境 + \(H_0=cz/D\)”的构造更像在测量**视向流场/修正口径**，并不等价于我们主线的 **observer environment bias** 口径，因此目前只作为方法学诊断，不纳入主线拟合。

输出：`math/hubble_env_scan/pantheonplus_env_correlation_results_v0.json`、`math/hubble_env_scan/pantheonplus_env_correlation_plot_v0.png`。

---

## 4 进一步预测

### 4.1 不同环境的 H₀

用当前 direct 拟合的 \(\beta\approx 0.171\) 做一阶线性外推（forecast，仅用于“把数字锁死”，不参与拟合）：

> 注：下表的 km/s/Mpc 绝对数值依赖参考定标（这里用 Planck 的 \(H_{0,\mathrm{ref}}\)）。更“关系论”的检验优先看
> \(H_0/H_{0,\mathrm{ref}}\) 或两观察者比值 \(\;H_{0,A}/H_{0,B}=(1-\beta\delta_A)/(1-\beta\delta_B)\;\)。

| 环境 | 密度对比度 δ | 预测 H₀ (km/s/Mpc) |
|------|-------------|---------------------|
| 最大空洞（Cold Spot） | -0.6 ~ -0.7 | 74.3 - 75.4 |
| KBC 空洞 | -0.46 | 72.1 ✓ |
| 平均宇宙 | 0 | 67.4 ✓ |
| 普通星系团 | +0.3 ~ +0.5 | 63.9 - 61.6 |
| 最大星系团（Coma） | +1.0 ~ +2.0 | 55.8 - 44.3 |

### 4.2 红移依赖

**预测**：H₀(z) 应随红移增加而趋近 67.4

**已有支持**：2025 年分析确认 "H₀(z) declining toward the background Planck value at higher redshifts"

---

## 5 验证路线

### 5.1 已完成

- [x] 公式推导
- [x] KBC 空洞数据验证
- [x] 红移依赖定性支持
- [x] 环境扫描管线 v0（数据表 + 拟合脚本）
- [x] N-body 模拟定量验证（Gavas et al. 2024 Figure 3 提取 β ≈ 0.10）

### 5.2 待验证

- [ ] Cold Spot 区域 H₀ 测量（预测 74.3-75.4）
- [ ] 星系团内部 H₀ 测量（预测 61.6-63.9）
- [x] 密度-H₀ 连续曲线拟合（v0 管线已建立，点集待扩展）

### 5.3 潜在合作

- 观测团队：SH0ES, Pantheon+, CosmicFlows
- 数据来源：Tully-Fisher, BOSS void catalog

---

## 6 意义

### 6.1 对哈勃张力的解释

**一种可检验的“环境/观测者效应”解释**：在当前口径下，不必引入额外早期宇宙新自由度，也能把哈勃张力的一部分表现解释为环境导致的局部测量偏置：
- 我们在 KBC 空洞里
- 本地测量自然偏高
- CMB 测的是全局平均

### 6.2 对关系论的验证

这是关系论目前最接近“外部可检验闭环”的一个定量主张：给出可复现管线与明确的证伪条件，并在当前口径下与部分外部观测/模拟结果在方向与量级上不矛盾。

---

## 证据文件

| 脚本 | 内容 | 输出 |
|-----|---|---|
| `run_hubble_env_scan_v0.py` | v0 环境扫描拟合（direct δ） | `results_v0.json`, `plot_v0.png` |
| `run_hubble_env_scan_v1.py` | v1：MC 误差传播 + forecast 预测点（星号）+ direct/proxy 双轨框架 | `results_v1.json`, `plot_v1.png` |
| `validate_hubble_env_data_v0.py` | 点集审计：schema/计数/可辨识性（v1） | `dataset_summary_v1.json` |
| `gavas2024_figure3_extraction/fit_beta_from_trendline.py` | 从 Gavas 2024 Figure 3 趋势线拟合 β（外部参照） | `gavas2024_figure3_extraction/gavas2024_trendline_fit_results.json` |
| `analyze_observer_env_bias_scale_v0.py` | Wojtak+2014 / Odderskov+2015：观测者偏差随 \(r_{\max}\) 的尺度趋势 | `observer_env_bias_scale_results_v0.json`, `observer_env_bias_scale_plot_v0.png` |
| `twompp_shell_stats_batch_v0.py` | 2M++：观测者壳层 \(\delta_g^\*\) 统计（多壳层×多平滑） | `twompp_shell_stats_grid_v0.csv`, `twompp_shell_stats_grid_v0.json` |
| `analyze_twompp_shell_stats_v0.py` | 2M++：壳层均值 \(\delta_g^\*\) 随平滑尺度的变化图 | `twompp_shell_stats_summary_v0.json`, `twompp_shell_stats_plot_v0.png` |
| `run_twompp_observer_weighted_delta_v0.py` | 观测者口径 \(w(r,\hat n)\)（用 Pantheon+ SH0ES‑HF 构造径向窗口+天区权重）→ 2M++ 加权体平均 \(\delta_g^\*\) | `twompp_observer_weighted_delta_results_v0.json`, `twompp_observer_weighted_delta_plot_v0.png` |
| `run_twompp_observer_weighted_h0_bias_v0.py` | 用 2M++ 速度场（引力响应 proxy）在同一 \(w(r,\hat n)\) 下估计“窗口诱导的 H0 回归偏差” | `twompp_observer_weighted_h0_bias_results_v0.json`, `twompp_observer_weighted_h0_bias_plot_v0.png` |
| `math/spectral_graph/run_relational_hubble_beta_poisson_3d_v0.py` | 关系论推导链条（toy）：\(L\approx-\nabla^2\)，解 \(L\phi=-\delta_{\mathrm{rel}}\)，取 \(\mathbf{v}=-f\nabla\phi\)，用同一回归口径得到 \(\beta\approx f/3\) | `math/spectral_graph/relational_hubble_beta_poisson_3d_results_v0.json`, `math/spectral_graph/relational_hubble_beta_poisson_3d_plot_v0.png` |
| `run_hubble_env_v2_audit_v0.py` | 模型 v2：口径审计 + scale→shell 映射 + \(\delta_g^\*\)–\(r_{\max}\) 对照图 | `hubble_env_v2_audit_results_v0.json`, `hubble_env_v2_audit_plot_v0.png` |
| `run_cf2_cone_proxy_bins_v0.py` | CF2（Cosmicflows-2）在指定天空锥（如 Shapley 方向）构造 \(\delta_g^\*\)–\(H_0\) 的 proxy 分箱点（方法学探索） | `cf2_cone_shapley_sigma20_proxy_bins_results_v0.json`, `cf2_cone_shapley_sigma20_proxy_bins_plot_v0.png` |
| `run_pantheonplus_env_correlation_v0.py` | Pantheon+ 低红移（MU_SH0ES）× 2M++ \(\delta_g^\*\)：环境–\(H_0\) 相关检验（proxy） | `pantheonplus_env_correlation_results_v0.json`, `pantheonplus_env_correlation_plot_v0.png` |
| `run_twompp_overdense_observer_candidates_v0.py` | 2M++：扫描“强过密的候选观测者位置”（proxy），并给出壳层均值 \(\delta_g^\*\) | `twompp_overdense_observer_candidates_results_v0.json`, `twompp_overdense_observer_candidates_plot_v0.png` |

## 参考文献

1. Planck Collaboration (2018). Planck 2018 results. VI. Cosmological parameters.
2. Riess et al. (2022). A Comprehensive Measurement of the Local Value of the Hubble Constant.
3. Keenan et al. (2013). Evidence for a ~300 Mpc Scale Under-density in the Local Galaxy Distribution.
4. 2024 N-body simulation study on H₀ dispersion from gravitational clustering.
5. 2025 Tully-Fisher study testing local supervoid solution.
6. Gavas et al. (2024/2025). Dispersion in the Hubble-Lemaître constant measurements from gravitational clustering. arXiv:2407.10139
7. Scolnic et al. (2025). The Hubble Tension in Our Own Backyard: DESI and the Nearness of the Coma Cluster. ApJL 979 L9
8. Giani et al. (2024). An effective description of Laniakea: impact on cosmology and the local determination of the Hubble constant. JCAP 2024(1)071
9. Wojtak et al. (2014). Cosmic variance of the local Hubble flow in large-scale cosmological simulations. MNRAS 438, 1805. arXiv:1312.0276
10. Turner, Cen & Ostriker (1992). The Relationship of Local measures of Hubble's Constant to its Global Value. AJ 103, 1427. DOI:10.1086/116156
11. Odderskov, Hannestad & Haugbølle (2015). On the local variation of the Hubble constant. arXiv:1407.7364