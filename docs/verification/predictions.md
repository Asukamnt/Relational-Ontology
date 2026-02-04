# 可验证预测清单

---

## 状态说明

| 符号 | 含义 |
|------|------|
| ✅ | 有独立外部证据支持（方向/量级一致；仍可能受口径与系统误差影响） |
| 🔄 | 初步支持，待扩点与系统对照（整体结论尚不稳固） |
| ⏳ | 待验证 |
| ❓ | 需要新实验/观测 |

注：同一主题可能包含多条“子预测”（下表逐条标状态）；即使某条子预测已有 ✅，主题整体仍可能因点集/口径不足而保持 🔄。

**统一验证格式**：见 [`docs/verification/pipeline.md`](pipeline.md)。

---

## 宇宙学预测

### 1. 哈勃张力（环境修正） 🔄

**公式（direct / 文献 δ 口径；当前点集较小）**：
\[
H_0(\delta)=H_{0,\mathrm{ref}}(1-\beta\,\delta)
\]

- \(H_{0,\mathrm{ref}}\approx 67.36\)（Planck；作为 \(\delta\approx 0\) 的参考锚点/口径选择）
- \(\beta\approx 0.171\)（direct 拟合；MC 16/50/84：0.139 / 0.169 / 0.201）

**外部参照（不混合口径）**：Gavas 2024（\(\delta_{\mathrm{den}}\) within 10 Mpc/h）给出 \(\beta\approx 0.085\)（方向一致、量级更小；\(\delta\) 定义/尺度不同，单独对照即可）。

**模型 v2（口径分层）的一个“硬检验点”**：观测者环境的 proxy \(\delta_{\mathrm{observer}}(R)\) 可以从 2M++ 的 \(\delta_g^\*\) 壳层均值可重复计算（见 `math/hubble_env_scan/twompp_shell_stats_grid_v0.csv`，以及图 `math/hubble_env_scan/twompp_shell_stats_plot_v0.png`）。
它显示：当尺度/平滑增大时，\(\delta_g^\*\) 的壳层均值快速趋近于 0——因此若某个 \(H_0\) 测量口径主要依赖“大尺度壳层”，环境效应应快速衰减、\(H_0\) 也应更接近 CMB 值；反之，若测量更依赖近邻小尺度，则环境项更可能显著。

| 预测 | 状态 |
|------|------|
| KBC 空洞 \(H_0\) ≈ 72–73（对照：72.1±0.9 / 73.17±0.86） | 🔄 |
| 密度-H₀ 负相关 | ✅ |
| H₀(z) 趋近 CMB 值 | 🔄 |
| 密度-H₀ 连续曲线拟合（v0 管线） | 🔄 |
| Cold Spot \(H_0\) ≈ 74.3–75.4（取 \(\delta\approx-0.6\sim-0.7\) 外推） | ⏳ |
| 典型星系团 \(H_0\) ≈ 61.6–63.9（取 \(\delta\approx0.3\sim0.5\) 外推） | ⏳ |

#### 封存预测（2026-02-02，Einstein-style：先写死再等验证）

为避免“事后调参”，在 **direct / 文献 \(\delta\)** 口径下，我们把当下采用的参数与由此推出的数值**封存**如下（未来出现新点集时，可另起一条新版本预测，但不得回溯改写本条）：

- 固定：\(H_{0,\mathrm{ref}}=67.36\)（Planck 2018；作为 \(\delta\approx 0\) 的“平均观察者参考值”，不是上帝视角）
- 固定：\(\beta=0.171\)（direct 拟合点估计；见 `math/hubble_env_scan/results_v1.json`）
- 计算：\(H_0(\delta)=H_{0,\mathrm{ref}}(1-\beta\delta)\)

> 备注（关系论口径）：上表给出的 km/s/Mpc 绝对数值依赖一个参考定标（这里选 Planck 作为 \(H_{0,\mathrm{ref}}\)）。
> 模型的“硬内容”首先是无量纲的相对偏差律 \(H_0/H_{0,\mathrm{ref}}=1-\beta\delta\)，或更一般的两观察者比值
> \(\;H_{0,A}/H_{0,B}=(1-\beta\delta_A)/(1-\beta\delta_B)\;\)。

| 环境（示例） | \(\delta\)（direct） | 封存预测 \(H_0\) (km/s/Mpc) |
|---|---:|---:|
| Cold Spot 级空洞 | -0.70 | 75.42 |
| Cold Spot 级空洞 | -0.60 | 74.27 |
| 典型星系团 | +0.30 | 63.90 |
| 典型星系团 | +0.50 | 61.60 |
| Coma 级强过密（下限） | +1.00 | 55.84 |
| Coma 级强过密（上限） | +2.00 | 44.32 |

> 备注：这条“封存预测”要检验的是 **observer-environment** 口径（观测者所处的大尺度环境 \(\delta_{\mathrm{obs}}\)），不等价于“我们站在原点去测某个目标所在环境”的 target-proxy 相关。

**说明**：已建立 v0 “环境扫描”管线（`math/hubble_env_scan/`），并引入 v1 数据表与 `direct/proxy` 双轨拟合框架（2M++ δ_g* 代理口径）。目前外部证据对“方向（过密→更低 H₀，欠密→更高 H₀）”有支持，但点集与尺度口径仍偏少/不统一，且 proxy 拟合样本量有限，因此整体状态标记为 🔄。

**证伪条件（v0 口径）**：
- 若多环境拟合得到 β ≤ 0 且显著性 ≥ 2σ，则密度-膨胀耦合方向需被修正。
- 若在同一 \(\delta\) 定义/尺度下，出现系统性违反单调性：\(\delta_A<\delta_B\) 但测得 \(H_{0,A}\le H_{0,B}\)（≥2σ），则该模型失效。
- 若在 δ > 0 的过密环境中测得 H₀ 系统性高于 \(H_{0,\mathrm{ref}}\)（≥2σ），则该模型失效。
- 若残差呈显著非线性且一阶修正无法吸收，则需要升级模型（非线性或分尺度）。

**证伪条件（v2：尺度/口径依赖）**：
- 若在明确的“同一 \(\delta\) 定义 + 同一尺度”下，独立点集给出 \(\beta\le 0\)（≥2σ），则 v2 失效。
- 若 \(\delta_{\mathrm{observer}}(R)\) 已趋近 0（大尺度壳层/大平滑），但对应口径的 \(H_0\) 仍稳定显著偏离 \(H_{0,\mathrm{ref}}\)（排除系统误差后），则“环境项主导”的解释需要被放弃或大改。

### 2. 暗物质 🔄

**核心主张**：暗物质不是"物质"，是"隐藏的关系/边"。

引力 = 向高 degree（连接度）方向的漂移。如果星系外围的旋转速度比预期快，说明那里的有效 degree 比可见物质显示的高——存在"暗边"。

| 预测 | 状态 |
|------|------|
| 引力 = 向高 degree 漂移（toy 验证） | ✅ |
| 星系旋转曲线可从关系密度推导 | 🔄 |
| 暗物质晕 = 隐藏的长程弱连接 | ⏳ |
| 引力透镜效应与关系拓扑相关 | 🔄（toy） |

**数值验证（v0 toy）**：
- `run_dark_matter_entropy_vs_ds_v0.py`：比较 S, d_s, degree 三个场，发现 **degree 是唯一能把 walkers 完全吸向核心的场**
- `run_dark_matter_rotation_curve_v0.py`：core-halo 图的"旋转曲线"
- `run_dark_matter_escape_velocity_v0.py`：逃逸概率测试
- `run_dark_matter_lensing_v0.py`：透镜/时间延迟（CTQW + CAP）：core-only vs core+halo 的传播偏折 \(\Delta y\) 与峰到达延迟 \(\Delta t\)

**证伪条件**：
- 若高 degree 区域不产生向心漂移（在足够大的图上），则"引力=连接度梯度"失效
- 若旋转曲线无法用 effective degree 分布拟合，则需要真正的新粒子

#### 封存预测（2026-02-02，Einstein-style：暗物质涌现的相变条件）

从 v3 扫描的 215 个参数组合中（见 `math/spectral_graph/run_dark_matter_phase_boundary_v0.py`），我们提取出**平坦旋转曲线涌现的必要条件**，用三个无量纲参数描述：

| 无量纲参数 | 定义 | 临界区间 | 物理含义 |
|-----------|------|----------|----------|
| 密度对比 | \(\displaystyle\frac{k_{\mathrm{core}}}{k_{\mathrm{halo}}}\) | **2 – 4** | 核心要比光环致密 2–4 倍 |
| 耦合强度 | \(p_{\mathrm{core\text{-}halo}} \times n_{\mathrm{core}}\) | **0.07 – 0.49** | 每个 halo 节点从 core 获得的期望连接数 |
| 光环填充 | \(\displaystyle\frac{n_{\mathrm{halo}}}{n_{\mathrm{core}}}\) | **1 – 7.5** | 光环规模要比核心大 1–7.5 倍 |

**关键约束（最强判别力）**：
\[
\boxed{0.07 \lesssim p \cdot n_{\mathrm{core}} \lesssim 0.5}
\]

- **太弱**（< 0.07）→ 核心-光环"断联"，无引力传递，曲线不平坦
- **太强**（> 0.5）→ 核心-光环"融合"，丧失"暗晕"势阱效应

**预测（可检验）**：
1. 满足上述条件的星系，应具有平坦旋转曲线（"有暗物质"）
2. 不满足上述条件的星系，应**不**具有平坦旋转曲线
3. 若星系演化跨越临界面（如并合/瓦解），其旋转曲线形态应随之改变

**与 Tully-Fisher 关系的联动**：在满足相变条件的结构族内部，存在 \(M \propto v_{\mathrm{flat}}^\alpha\) 的幂律关系，toy 模型给出 \(\alpha \approx 3.9\)（观测值 \(\alpha \approx 3.94\)）。

**证伪条件（相变条件版）**：
- 若存在大量违反上述边界、但仍有平坦旋转曲线的星系，则该模型失效
- 若在满足条件的区间内，TF 指数 \(\alpha\) 系统性偏离 4（如 \(\alpha < 2\) 或 \(\alpha > 6\)），则拓扑解释需要大改

### 3. 暗能量 🔄

**预测**：暗能量是关系网络节点增加（熵增）的表现

| 预测 | 状态 |
|------|------|
| 膨胀加速与节点增加速率相关 | 🔄 |
| 不需要宇宙常数作为独立实体 | ⏳ |

---

## 量子物理预测

### 4. 量子纠缠 🔄

**预测**：纠缠粒子 d_R ≈ 0，d_M 很大（投影假象）

| 预测 | 状态 |
|------|------|
| 与 ER=EPR 猜想一致 | 🔄 |
| 纠缠强度与关系强度正相关 | ⏳ |
| Bell/CHSH 违反（S>2）且满足无信号（no-signalling） | 🔄（toy） |

**当前证据（toy，兼容性检查）**：
- `math/tensor_network/mera_demo_v0.py`：树距离 vs 纠缠代理 **r=-0.837**（“纠缠强→投影距离近”）
- `math/tensor_network/run_quantum_chsh_decoherence_v0.py`：**CHSH S=2.828≈2√2**，且无信号误差 max Δ≈0
- `math/tensor_network/run_entanglement_hidden_relations_distance_test_v0.py`：在固定 \(d_{\mathrm{geo}}\) 下扫描隐藏边强度，观察相关性主要随 \(d_{\mathrm{rel}}\) 变化；并在扫描 \(d_{\mathrm{geo}}\) 时对比“无隐藏关系”与“有隐藏关系”的衰减差异（相关性代理，toy）
- `math/tensor_network/run_entanglement_hidden_relations_chsh_bridge_v0.py`：桥接 demo：用 \(C(d_{\mathrm{rel}})=e^{-\alpha d_{\mathrm{rel}}}\) 把关系距离映射为 concurrence，再计算最大 CHSH \(S_{\max}\) 并显式 no-signalling 自检（toy）
- `math/tensor_network/run_entanglement_hidden_relations_chsh_bridge_v1.py`：桥接 v1：用传播子相关读出 \(p\) 并用最大熵 Werner 映射得到 \(S_{\max}\)（toy）
- `math/tensor_network/run_entanglement_hidden_relations_chsh_emergence_dynamics_v0.py`：涌现动力学：隐藏边按 \(|\rho|\) 自洽生长/淘汰，最强长程对出现 \(S_{\max}>2\)（toy）
- `math/tensor_network/run_entanglement_hidden_relations_chsh_emergence_dynamics_v1.py`：涌现动力学 v1：局部提议-强化-衰减（无全局 top‑K），best 长程对出现 \(S_{\max}>2\)（toy）
- `math/tensor_network/run_entanglement_hidden_relations_chsh_emergence_dynamics_v2.py`：涌现动力学 v2：候选集合由随机游走发现（无固定远距候选集）；默认点未越过 \(S=2\)（负结果记录，toy）
- `math/tensor_network/run_entanglement_hidden_relations_chsh_emergence_dynamics_v3.py`：涌现动力学 v3：随机游走候选 + 局部预算（mutual top‑K）+ 远近权衡，best \(S_{\max}>2\)（toy）
- `math/tensor_network/run_entanglement_hidden_relations_operational_chsh_v0.py`：Operational CHSH：不借用量子态映射，直接抽样/计数得到 \(S>2\)（含 base-only 与 \(\lambda=0\) 对照，且 no‑signalling 自检，toy）
- `math/tensor_network/run_entanglement_hidden_relations_operational_chsh_v1.py`：Operational CHSH v1：测量=松弛过程（Metropolis），纯计数得到 \(S>2\)（含对照与 no‑signalling 自检，toy）
- `math/tensor_network/run_entanglement_hidden_relations_operational_chsh_v2.py`：Operational CHSH v2：小子网 + “相位/复幅度中介”松弛；交叉项由中介闭合涌现（减少显式 \(J_{xy}\)），纯计数得到 \(S>2\)（含对照与 no‑signalling 自检，toy）
- `math/tensor_network/run_entanglement_hidden_relations_operational_chsh_v3.py`：Operational CHSH v3：中介通道升级为链/环小子网；交叉项由中介 Green 函数涌现，并比较 chain/ring 的稳健性（toy）
- `math/tensor_network/run_entanglement_hidden_relations_operational_chsh_v4.py`：Operational CHSH v4：中介通道加入 \(\phi^4\) 饱和非线性，并用对称翻转混合扇区后仍可 \(S>2\) 且 no‑signalling 自检通过（toy）
- `math/tensor_network/run_entanglement_hidden_relations_operational_chsh_v5.py`：Operational CHSH v5：中介通道改为角变量（XY/紧致约束），MCMC 松弛 + 纯计数得到 \(S>2\)（含对照与 no‑signalling 自检，toy）

> 说明：以上是“量子信息论口径的兼容性 suite”，不是“从关系动力学推导量子规则”的最终答案。

### 5. 波函数坍缩 🔄（toy）

**预测**：观测 = 建立关系 = 锁定可能性分布

| 预测 | 状态 |
|------|------|
| 测量扰动与关系建立强度相关 | ⏳ |
| 退相干：相干项 \(|\rho_{01}|\) 随环境规模 \(N\) 指数衰减（观测=建立关系的量化代理） | 🔄（toy） |
| 单次结果 + Born 频率：CTQW + 吸收屏幕（CAP）的 click 统计由 norm-loss 决定，且频率匹配 \(|\psi|^2\) 的吸收积分 | 🔄（toy） |
| which-way 退相干：若环境能局域区分两条路径（相位标记/关系建立），则双缝干涉项 \(I(y)\) 被压平，屏幕/点击分布趋近 \(P_1+P_2\) | 🔄（toy） |

**说明**：目前我们已在图上的量子行走（CTQW）里做出一个最小闭合：单位演化 + 非幺正“关系建立”（CAP）→ 单次 click + 频率律（Born-like）。  
但这仍是 toy：它尚未回答“为何在更一般情形下必然选择某个指针基/测量基”的最终版本，也未把量子规则（希尔伯特空间公设）完全从关系公理推导出来。

**证据（toy）**：
- `math/spectral_graph/run_graph_quantum_measurement_born_v0.py`：CAP（不可逆关系建立）→ 单次 click + Born‑like 频率
- `math/spectral_graph/run_graph_quantum_which_way_decoherence_v0.py`：局域 which‑way 标记（相位噪声）→ 干涉项被压平
- `math/spectral_graph/run_graph_quantum_measurement_which_way_clicks_v0.py`：把两者接起来做 click 分布对照：相干 `rel_L1≈0.423` → which‑way `rel_L1≈0.194`（趋近 \(P_1+P_2\)）

### 6. 不确定性原理 ⏳

**预测**：位置和动量是同一关系的两个投影面

| 预测 | 状态 |
|------|------|
| 不确定性来自投影，不是本体模糊 | 框架性解释 |

---

## 粒子物理预测

### 7. 三代费米子 🔄（toy）

**预测（toy 口径）**：三代 = 同一“电荷/拓扑缺陷”附近可稳定支持的 **3 个局域低模**（谱激发）；本征值 \(\lambda\) 作为“质量”代理。

| 子预测 | 状态 |
|------|------|
| 单缺陷出现三重局域低模平台（\(N_{\mathrm{localized}}=3\)） | 🔄（toy） |
| 若存在第四阶稳定模 → 第四代 | ❓ |
| 第四代夸克质量 8-20 TeV | 需下一代对撞机 |

**证据（toy）**：
- `math/spectral_graph/run_relational_link_variable_three_generations_v0.py`：扫描缺陷强度 \(\phi_0\)，统计缺陷附近的局域低模数量；在默认扫描中多处出现 \(N_{\mathrm{localized}}=3\) 平台。代表点：\(\phi_0\approx 1.689\) 给出三模 \(\lambda\approx 0.839,\ 0.969,\ 0.973\)（见 JSON/PNG）。
- `math/spectral_graph/run_relational_link_variable_su2u1_doublet_generations_v0.py`：把态提升为 SU(2) doublet（\(\mathbb{C}^2\)）后，在同一缺陷参数点得到 **6=3×2** 个局域候选，且 base 下谱两重简并（pair gap ~1e-15）；加入局域 Higgs/VEV 分裂项（\(\sigma_3\)）后简并被打破（见 JSON/PNG）。
- `math/spectral_graph/run_relational_link_variable_su2u1_doublet_generations_v1.py`：对 \((m_{\mathrm{split}}, r_H)\) 做小扫描，找到 split 后仍保持 **6=3×2** 的稳定区域，并自动选择代表点（例如 `m_split=0.30, r_higgs=3.0`，`ipr_over_uniform_min=4.0`）用于展示“分裂但对象仍在”（见 JSON/PNG）。
- `math/spectral_graph/run_relational_link_variable_su2u1_doublet_generations_v2.py`：把 \(U_{ij}\) 升级为非平庸 SU(2) 连接（局域 holonomy lump），对照显示 doublet 的精确两重简并被抬起（pair gap 从 \(\sim10^{-15}\) 到 \(\sim10^{-2}\)），从而不再是“两份拷贝”（见 JSON/PNG）。
- `math/spectral_graph/run_relational_link_variable_su2u1_doublet_generations_v3.py`：在固定非平庸 SU(2) 连接下扫描 \((m_{\mathrm{split}}, r_H)\)，找到 **SU2+Higgs 仍保持 6=3×2** 的区域与代表点（例如 `m_split=0.20, r_higgs=5.0`，`ipr_over_uniform_min=4.0`；见 JSON/PNG）。
- `math/spectral_graph/run_relational_link_variable_su2u1_doublet_generations_v4.py`：固定 Higgs 点（`m_split=0.20, r_higgs=5.0`）扫描 `su2_angle0` 并加密 1.1 附近；在离散网格上 ok6 仅在少数点出现（例如 `su2_angle0=1.10, 1.16`），提示相区边界可能很窄或对阈值敏感（见 JSON/PNG）。
- `math/spectral_graph/run_relational_link_variable_su2u1_doublet_generations_v5.py`：固定 `m_split=0.20` 扫 `(su2_angle0 × r_higgs)` 的 2D 相图切片；观测到 **4/55** 个 ok6 cell（例如 `(1.10,5.0),(1.00,2.0),(0.80,6.0),(1.20,6.0)`），说明 ok6 区域存在但在粗网格上呈碎片状，需要局部加密提取连续边界（见 JSON/PNG）。
- `math/spectral_graph/run_relational_link_variable_su2u1_doublet_generations_v6.py`：对 v5 的关键区域做 **zoom 加密**（`su2_angle0∈[1.08,1.22]`、`r_higgs∈[4.6,6.2]`），在 135 个 cell 中出现 **35** 个 ok6，并形成最大连通块 **21** 个 cell（相区不再是孤立点，可提取边界/面积占比；见 JSON/PNG）。
- `math/spectral_graph/run_relational_link_variable_su2u1_doublet_dirac_kahler_v0.py`：把同一 SU(2)×U(1) 背景下的节点协变拉普拉斯 \(L\) 以 Dirac–Kähler 方式“开方”成 \(D\)：构造 covariant incidence \(B\) 并验证 \(L=B^\*B\)（`rel≈1.7e-17`），且 \(\{\gamma_5,D\}=0\)（数值上为 0）；局域候选仍为 `n_candidates=6`，并呈现 \(\pm\sqrt{\lambda}\) 的粒子/反粒子成对谱（见 JSON/PNG）。
- `math/spectral_graph/run_relational_link_variable_su2u1_chiral_doublet_singlet_v0.py`：**弱相互作用手征性（toy）**：在同一 U(1) 磁通团背景上构造 LEFT(SU2×U1 doublet) 与 RIGHT(U1 singlet) 的协变拉普拉斯；在同一 IPR+峰距判据下得到 `LEFT=6`、`RIGHT=3`，且 SU(2) holonomy lump 只会在 LEFT 上把 pair-gap 中位数从 \(\sim10^{-15}\) 抬升到 \(\sim10^{-2}\)（RIGHT 不受影响；见 JSON/PNG）。
- `math/spectral_graph/run_relational_link_variable_su2u1_chiral_yukawa_v0.py`：**Yukawa/VEV（toy）**：在 LEFT(SU2×U1 doublet) 与 RIGHT(U1 singlet) 上加入局域 Yukawa 混合 \(H=\begin{pmatrix}L_L&Y\\Y^\*&L_R\end{pmatrix}\)，且 \(Y\) 只耦合 LEFT 的 down 分量（Higgs VEV 方向，支持于半径 `r_higgs`）；默认点（取 `su2_angle0=0` 隔离机制）下得到 `n_candidates=9`，其中 3 个几乎纯 `LEFT_up`（`w_right≈0`），6 个为 `LEFT_down` 与 `RIGHT` 近似 50/50 混合（`w_right≈0.5`；见 JSON/PNG）。
- `math/spectral_graph/run_relational_link_variable_su2u1_chiral_yukawa_v1.py`：**协变 Higgs 方向（toy）**：用 SU(2) 并行运输构造节点依赖 Higgs 双态方向 \(h(x)\)，并将 split 与 Yukawa 写成 \(V_H=m(2|h\rangle\langle h|-I)\)、\(Y=yh(x)\)，避免“固定 down 分量”在 SU(2) 非平庸连接下的基依赖（见 JSON/PNG）。
- `math/spectral_graph/run_relational_link_variable_su2u1_chiral_yukawa_covhiggs_scan_v4.py`：**SU(2) 开启下的 Yukawa 相区（toy）**：在 `su2_angle0=1.1` 固定下对 \((m_{\mathrm{split}}\times yukawa_m)\) 做 zoom 扫描，出现连通 ok 区（`ok cells=2/25`，例如 `m_split≈0.64–0.65`、`yukawa_m≈0.50`），在 `w_right` 阈值口径下得到 `3(neutral-like)+6(charged-like)`（见 JSON/PNG）。

### 8. 缪子 g-2 异常 ⏳

**预测**：高阶关系振荡产生附加相位

| 预测 | 状态 |
|------|------|
| 异常与关系密度环境相关 | ⏳ |

---

## 宏观物理预测

### 9. 黑洞奇点 ⏳

**预测**：奇点 = 连接密度超阈值，投影失效

| 预测 | 状态 |
|------|------|
| 信息不丢失，退出时空表示 | 与全息原理一致 |

### 9.1 黑洞振铃（ringdown）尺度律 🔄

**预测（可证伪口径）**：若“有限时间逃逸概率”定义出的视界代理与真实黑洞半径在结构上对应，则 ringdown 的主频应满足近似标度
\[
\omega \propto \frac{1}{R_h}\quad(\text{等价于 } \omega R_h \approx \text{const})
\]
这对应 GR 中“ringdown 频率 \(\propto 1/M\)”的最基本维度分析结论（忽略自旋等额外无量纲参数时）。

**当前进展（toy）**：
- **谱法（非厄米 QNM proxy）**：`math/spectral_graph/run_self_gravitating_ringdown_qnm_scan_v0.py`
  - 生成自引力几何 \(L(w)\)
  - 用 \(P_{escape}(T)\) 得到视界 mask 与 \(R_h^{area}=\sqrt{A_h/\pi}\)
  - 在吸收边界下构造非厄米生成元 \(A=(-i\gamma)L-\kappa_{out}\Pi_{out}-\kappa_{hor}\Pi_{hor}\)
  - 从特征值 \(\lambda=\sigma+i\omega\) 抽取 “QNM proxy” 的 \(\omega,\tau,Q\)，并拟合 \(\omega\)–\(R_h\) 的 log-log 斜率
- **波方程口径（更接近 ringdown 的“波动”定义）**：`math/spectral_graph/run_self_gravitating_ringdown_wavefit_scan_v0.py`
  - 在自引力几何上演化阻尼波方程 \(u_{tt}=-c^2Lu-\Gamma u_t\)
  - 在探针点拟合 \(\omega,\tau\)
  - 同时输出两种半径代理：
    - \(R_h^{area}=\sqrt{A_h/\pi}\)（节点面积口径）
    - \(R_h^{proper}\)：以 \(\ell_{ij}=1/\sqrt{w_{ij}}\) 的最短路距离定义（取中心到视界边界的中位数距离）
- **控制验证（解析可解，校验提取管线）**：`math/spectral_graph/run_wavefit_control_ring_lattice_v0.py`
  - 在 1D 环格点上解析地知道 \(\omega_k=c\sqrt{\lambda_k}\)，并且对 \(k=1\) 有 \(\omega\sim 1/N\)
  - 数值结果（见 `wavefit_control_ring_lattice_results_v0.json`）给出拟合斜率约 \(p\approx -1.00\)，说明“波方程+拟合”的 ω 提取本身是可靠的

**状态说明**：目前该标度在 toy 中对“视界口径/吸收层/模式选择/拟合窗口”仍敏感，尚不足以宣称已经稳定复现 \(-1\) 斜率，因此标记为 ⏳（待稳健化与对照校准）。

**v0 现象提示（用于自检，而非最终结论）**：
- 以当前脚本默认口径：
  - 谱法（`self_gravitating_ringdown_qnm_scan_results_v0.json`）的 \(\omega\)–\(R_h\) 斜率仍可能显著偏离 \(-1\)。
  - 波方程拟合（`self_gravitating_ringdown_wavefit_scan_results_v0.json`）在一组默认 toy 参数下得到：
    - \(\omega\)–\(R_h^{area}\) 的斜率约 \(p\approx 0.04\)
    - \(\omega\)–\(R_h^{proper}\) 的斜率约 \(p\approx -0.11\)
  这说明：当前模型下 ringdown 主频对“视界半径代理”的依赖仍偏弱/未呈现 GR 期望的 \(-1\) 标度，需要继续升级口径或动力学。
- **有限盒效应诊断**：当采用“散射式激发”（`--excite-mode scatter`）时，若外吸收层不够“无反射”，拟合到的 ω 可能更像由“有效腔长”控制（例如与 \(L_{\text{cavity}}=R_{\text{outer}}-R_h^{proper}\) 呈近似幂律），提示观测到的“振铃”可能来自有限域/边界反射而非黑洞势垒的 QNM。
- **边界敏感性扫描（固定几何）**：`math/spectral_graph/run_self_gravitating_ringdown_boundary_sensitivity_scan_v0.py`（结果：`self_gravitating_ringdown_boundary_sensitivity_results_v0.json` / `..._plot_v0.png`）
  - 固定同一套自引力几何（同一 mass_total/seed/alpha_w），只扫描外吸收层参数（outer_layer/profile/gamma_outer）。
  - 由于 \(\omega\) 提取对“时间窗口/拟合起点”敏感，扫描中我们改用**统一口径**（源/探针相对中心固定；固定早期时间窗），以避免“窗长变化导致 FFT 分辨率变化”的伪漂移。
  - 在该口径下：`shell` 的 \(\omega\) 对 outer_layer 的漂移较小（典型 ω_rel_span≈0.11–0.13），`scatter` 的漂移更明显（典型 ω_rel_span≈0.24–0.27），并出现 \(\omega\) 随 \(L_{\text{cavity}}\) 的幂律趋势（盒子效应更强）。
  - 这提示：若要把“ringdown=QNM”作为可证伪对应，必须先把“外边界实现”做成对 \(\omega\) **不主导**的设置（更强的 PML、或更大的域、或更严格的窗口定义）。
  
- **推远外边界（但保持核心几何不变）**：`math/spectral_graph/run_self_gravitating_ringdown_padding_domain_scan_v0.py`
  - 为避免“视界 proxy 随域大小漂移”的混淆，我们先在 \(N_0\times N_0\) 上生成一次自引力几何与 horizon mask，然后把它嵌入更大的平坦网格 \(N\times N\)，仅仅把边界推远。
  - 在固定早期时间窗的口径下，`shell` 的 \(\omega\) 在 \(N=30\to 48\) 上几乎不变（示例：约 0.805），说明该 ω 更像“局域响应/几何固有尺度”；而 `scatter` 更接近“入射脉冲频谱”，需要另行定义“散射后的 ringdown 窗口”才能讨论腔体/共振。
  - **散射后的 ringdown 窗口（可操作口径）**：对 `scatter`，把 ringdown 定义为“探针信号的**第一个显著峰值**之后的一段固定长度窗口”（例如 10s），即 `--fit-start-mode first_peak --fit-max-duration 10`。
    - 在该口径下（示例：baseN=30，padding 到 \(N=48,60,72\)，且 \(L_{\text{cavity}}\) 明显变大），拟合到的 \(\omega\) 仍基本保持不变（约 0.688），这表明**早期 ringdown** 可以做到对外边界不敏感；而“更长窗口/更晚时刻”才更容易被腔体低频成分主导。

- **padded 质量扫描（用稳定窗口测试尺度律）**：`math/spectral_graph/run_self_gravitating_ringdown_mass_scan_padded_v0.py`
  - 以 `scatter` 的“first_peak + 10s”窗口为 ringdown 口径，在固定 baseN=30、padN=72 下扫描 mass_total（多 seed）。
  - 当前 toy 参数下得到的 \(\omega\)–\(R_h^{proper}\) 拟合斜率仍为**正**（示例：\(p\approx +1.76\)），与 GR 期望 \(-1\) 不一致，提示“该 ω proxy”更像由权重/耦合强度（局域刚度）控制，而非单纯由视界半径控制。

- **势垒项升级（把“曲率散射势”显式写入波动方程）**：`math/spectral_graph/run_self_gravitating_ringdown_potential_beta_scan_v0.py`
  - 动机：GR 的 ringdown/QNM 在有效 1D 径向方程里由 Regge–Wheeler / Zerilli 势垒控制；仅用 \(L(w)\) + 阻尼时，toy 中的 \(\omega\) 往往更像“局域刚度”而不是“视界尺度”。
  - 做法（toy）：在波方程里加入对角势
    \[
    M(w)u_{tt}+\Gamma u_t + c^2Lu + \beta V(w)u=0,\quad V(w):=|\nabla\log w|^2
    \]
    其中 \(V\) 在外部区域取值并做 percentile 归一化；\(\beta\) 为强度系数。
  - 现象（padded 域、seeds=0..2、扫描 \(\beta=0..6\)）：\(\beta\) 增大时，\(\omega\)–\(R_h^{proper}\) 的 slope 会从接近 0 转为显著负值，出现“更像 QNM 势垒控制”的标度趋势。一次扫描（`self_gravitating_ringdown_potential_beta_scan_results_v0.json`）的 slope_ok(proper) 示例如下：
    - β=0：p≈+0.12
    - β=1：p≈−0.50
    - **β=2：p≈−1.17（在该扫描点集中最接近 −1）**
    - β=3：p≈−1.49
    - β=4：p≈−1.63
    - β=5：p≈−1.46
    - β=6：p≈−1.78
  - 说明：默认 `--fit-method lsq` 在部分点会出现 `lsq_failed`（此时仍保留 FFT 的 \(\omega\)），脚本用 **ω>0 的全部点** 计算 `slope_ok`，并同时报告 `fit_ok` 比例以便诊断。
  - 进一步细扫（padded 域、seeds=0..2、β=1.4..2.2 step=0.1，见 `self_gravitating_ringdown_potential_beta_scan_fine_1p4_2p2_step0p1_results_v0.json`）：
    - **β=1.7：p≈−1.021（abs error≈0.021）**（当前这组 toy 设置下“最接近 −1”的点）
  - **重要备注（已推进）**：我们新增了一个把 \(\beta\) 从“扫描调参”推进到“结构可计算”的口径（`auto_dim`），把 **\(\beta\)** 绑定到“几何的有效维数”上，从而在保持 \(V\) 定义不变的前提下，给出一个独立约束。

- **β 的结构约束（`auto_dim`：由“体积增长维数”给出无量纲系数）**：
  - 代码升级：
    - `math/spectral_graph/run_self_gravitating_ringdown_wavefit_scan_v0.py`
      - 新增 `--potential-beta-mode {manual,auto,auto_dim}`
      - 在输出 JSON 中记录 `potential.beta_eff`（实际进入方程的系数）以及 `potential.dim_growth`（维数估计详情）
    - `math/spectral_graph/run_self_gravitating_ringdown_mass_scan_padded_v0.py` 同步暴露该开关
  - 定义（toy，但**可计算**）：
    - 先在 proper 距离上做“体积增长”拟合：
      \[
      N(r):=\#\{i\in\text{mask}\mid d(i)\le r\}\ \propto\ r^{d_{\mathrm{eff}}}
      \]
      这里的 mask 取“视界外且不在吸收层内”的区域（并剔除最外一圈以避开边界伪影）。
    - 取径向化（radial reduction）对应的系数：
      \[
      f_{\mathrm{dim}}:=\frac{d_{\mathrm{eff}}-1}{2}
      \]
    - 当 \(V\) 采用 percentile 归一化（例如 p95）时，令脚本内部的有效系数为：
      \[
      \beta_{\mathrm{eff}}:=c_R^2\cdot \beta_0\cdot \text{scale}(V)\cdot f_{\mathrm{dim}}
      \]
      其中 `scale(V)` 是该几何上 \(V\) 在 mask 内的 p95（已记录在 `potential.scale`），\(\beta_0\) 为一个无量纲“基准系数”（`--potential-beta 0` 时取默认 1）。
      这样，方程中对 **未归一化的 \(V_{\mathrm{raw}}\)** 的实际系数变成 \(c_R^2\beta_0 f_{\mathrm{dim}}\)，从而把“归一化带来的尺度任意性”吸收掉，并引入一个由几何自洽给出的维数因子。
  - 结果（padded 质量扫描，scatter 的稳定窗口，seeds=0..2）：
    - 运行命令（示例）：
      - `python math/spectral_graph/run_self_gravitating_ringdown_mass_scan_padded_v0.py --seed-list 0,1,2 --mass-total-list 0.6,0.8,1.0,1.2,1.4 --pad-N 64 --dt 0.04 --steps 6000 --excite-mode scatter --scatter-offset 12 --scatter-probe-offset 12 --first-peak-min-frac 0.005 --fit-max-duration 10.0 --fit-end-mode full --fit-method lsq --lsq-min-cycles 1.0 --potential-mode gradlogw_sq --potential-beta-mode auto_dim --potential-beta 0 --potential-norm p95`
    - 得到 \(\omega\)–\(R_h^{proper}\) 的 log-log 斜率：**p≈−0.988（目标 −1）**
    - 同一批记录中典型统计量（来自输出 JSON 的汇总）：
      - \(d_{\mathrm{eff}}\) 均值约 **2.33**
      - \(f_{\mathrm{dim}}\) 均值约 **0.666（接近 2/3）**
      - \(\beta_{\mathrm{eff}}\) 均值约 **1.74**（与此前“细扫得到的 β≈1.7”一致）
  - 解释：
    - 这一步不再把 \(\beta\) 当作纯调参，而是把它视为“从几何有效维数（体积增长）到 1D 径向有效方程”的结构系数（toy 中可计算），因此比单纯扫描更接近“预测口径”。

**证伪条件（内部一致性）**：
- 若在固定口径（固定 \(T\)、阈值、吸收层参数）下，多 seed 扫描仍无法得到稳定的 \(\omega R_h\) 平台（或斜率显著偏离 \(-1\) 且不可通过合理的模式选择自洽修复），则“该视界定义 + 该 QNM 提取”的黑洞对应关系需要被修正。

### 10. 真空灾变 ⏳

**预测**：零点能取决于可激活边比例，非绝对总数

| 预测 | 状态 |
|------|------|
| 极端曲率区零点能跃升 | ⏳ |

---

## 实验室可测预测

### 11. 声致发光 ⏳

**预测**：气泡塌缩 = 关系网拓扑相变

| 预测 | 状态 |
|------|------|
| 电解质浓度改变光强 | 实验室可测 |

---

## 认知验证

### 12. 智能 = 同构 ⏳

**预测**：以关系网络为核心的认知系统，应能展现“智能 = 同构”所描述的跨域迁移与泛化行为。

| 预测 | 状态 |
|------|------|
| 概念图 = 内部关系网络 | ⏳ |
| 学习 = 同构映射 | ⏳ |
| ΔU/cost 驱动有效 | ⏳ |

---

## 关系网络动力学 / 新物理预测

### 13. 维度涨落 (Dimension Fluctuation) ⏳
**预测**：微观尺度下，$d_s$ 并非固定为 3，而是随采样尺度剧烈波动。
*   **现象**：高能散射截面出现“维度泄露”特征。
*   **状态**：⏳ 待计算具体的涨落谱。

### 14. 暗物质晕的粒度 (Halo Granularity) ⏳
**预测**：由于网络离散性，暗物质晕并非平滑流体，而是具有微观泡沫结构。
*   **现象**：引力透镜观测应存在微小焦散线（Caustics）。
*   **状态**：⏳ 待模拟精细结构。

### 15. 光速各向异性 (Anisotropy of c) ⏳
**预测**：光速是网络连接的统计平均，有限网络中必然存在微小涨落与方向依赖。
*   **现象**：极高精度干涉仪可测得方向性噪音。
*   **状态**：⏳ 待计算涨落量级。

### 16. 纠缠引力 (Entanglement Gravity) ⏳
**预测**：纠缠态改变了网络拓扑距离，从而改变局部熵力（有效引力）。
*   **现象**：纠缠态物质与非纠缠态物质的有效引力质量存在极微小差异。
*   **状态**：⏳ 待设计精密称重方案。

---

## 优先级排序

| 优先级 | 预测 | 理由 |
|--------|------|------|
| **P0** | 哈勃张力更多环境 | 已有支持，可扩展 |
| **P0** | 声致发光实验 | 实验室可做 |
| **P1** | 暗物质旋转曲线 | 数据丰富 |
| **P2** | 认知闭环验证 | 需独立实验与可复现管线 |
| **P2** | 粒子物理预测 | 需大型设备 |

---

## 框架自洽 / 数值 toy 验证（非观测预测）

这些不直接对应天文/实验室观测点，但用于检验“关系网络动力学”在最小模型下能否退化出熟悉结构：

- **图波动方程**：`math/spectral_graph/run_graph_wave_equation_v0.py`
- **图波动方程收敛性扫描（dt↓/n↑）**：`math/spectral_graph/run_graph_wave_equation_convergence_scan_v1.py`
- **图量子双缝（干涉）**：`math/spectral_graph/run_graph_quantum_double_slit_v0.py`
- **图薛定谔收敛性扫描（Crank-Nicolson，dt↓/n↑）**：`math/spectral_graph/run_graph_schrodinger_convergence_scan_v1.py`
- **图上的不确定性风格权衡**：`math/spectral_graph/run_graph_uncertainty_v0.py`
- **非均匀权重导致偏折（透镜类比）**：`math/spectral_graph/run_graph_geometric_lensing_v0.py`
- **自引力闭环（\(\rho\rightarrow w\rightarrow L\rightarrow S\rightarrow \rho\)）**：`math/spectral_graph/run_self_gravitating_network_v0.py`
- **自引力闭环相图扫描（alpha_w × gamma_S）**：`math/spectral_graph/run_self_gravitating_network_phase_scan_v1.py`
- **关系孤子 / 稳定子结构（物体=稳定模式）**：`math/spectral_graph/run_relational_soliton_v0.py`
- **关系孤子稳健性扫描（toy）**：`math/spectral_graph/run_relational_soliton_scan_v1.py`
- **关系孤子消融实验（toy）**：`math/spectral_graph/run_relational_soliton_ablation_v1.py`
- **关系孤子稳定性四件套（toy：持存/恢复/谱稳健/收敛）**：`math/spectral_graph/run_relational_soliton_stability_suite_v1.py`
- **关系孤子相互作用（并合/保持/解体，toy）**：`math/spectral_graph/run_relational_soliton_interaction_v0.py`
- **关系孤子相互作用多 seed 统计（toy）**：`math/spectral_graph/run_relational_soliton_interaction_scan_v1.py`
- **关系孤子相互作用阈值（cross_gamma 扫描，toy）**：`math/spectral_graph/run_relational_soliton_interaction_cross_gamma_scan_v1.py`
- **物体属性候选（质量/电荷 proxy，toy）**：`math/spectral_graph/run_relational_object_properties_v0.py`

- **量子硬指标（toy：CHSH + 无信号 + 退相干）**：`math/tensor_network/run_quantum_chsh_decoherence_v0.py`

### 基本常数内禀性验证

- **光速 c 由因果时间定义**：`math/spectral_graph/run_c_from_causal_time_v0.py`
- **光锥 c 缩放验证**：`math/spectral_graph/run_c_lightcone_scaling_v0.py`
- **c 从熵度量涌现（固定点迭代）**：`math/spectral_graph/run_c_entropy_metric_fixed_point_v0.py`
- **G 从熵梯度涌现**：`math/spectral_graph/run_g_from_entropy_gradient_v0.py`
- **G = 向高 degree 漂移**：`math/spectral_graph/run_dark_matter_entropy_vs_ds_v0.py`
- **暗物质旋转曲线 toy**：`math/spectral_graph/run_dark_matter_rotation_curve_v0.py`
- **暗物质逃逸速度 toy**：`math/spectral_graph/run_dark_matter_escape_velocity_v0.py`
- **精细结构常数 α 的操作性测量（toy）**：`math/spectral_graph/run_alpha_from_coulomb_3d_v0.py`（用 3D torus Poisson 的 \(\Delta U(r)\sim -q^2/(4\pi r)\) 系数定义 α_eff；当前仅完成“定义与测量口径”，未解释为何是 1/137）
- **α 的尺度流/屏蔽（toy）**：`math/spectral_graph/run_alpha_screening_hidden_edges_scan_v1.py`（在 3D torus 上加入“隐藏长程弱边”作为额外关系，α_eff 在中尺度出现近似平台，并随隐藏边密度被屏蔽；在某些 toy 参数下平台可落在 \(\sim 1/137\) 量级——这提供了“数值可检验的路线”，但仍不是严格推导）
- **α 的“几何约束屏蔽扫描”（toy）**：`math/spectral_graph/run_alpha_screening_geo_constraint_scan_v0.py`（把几何维度用 base torus 的 \(d_s^{geo}\approx3\) 约束住；把隐藏长程弱边只作为 Poisson 操作子的额外关系用于屏蔽，从而在 \(d_s^{geo}\approx3\) 的前提下扫描 \(w_{\mathrm{extra}}\times k_{\mathrm{extra}}\) 寻找 \(\alpha_{\mathrm{eff}}\) 的 plateau \(\approx 1/137\)）
- **α 的“隐藏关系自洽涌现→屏蔽”（toy）**：`math/spectral_graph/run_alpha_emergent_hidden_edges_v0.py`（隐藏长程弱边不作为外参扫描，而是按传播子相关 \(|\rho_{ij}|\) 的“提议-强化-衰减-预算（mutual top‑K）”规则自洽生长；在固定可见几何 \(d_s^{geo}\approx3\) 的前提下，用最终 \(A_{geo}+A_{hidden}\) 的 Poisson/Coulomb 测量 \(\alpha_{\mathrm{eff}}\)；在部分参数点可得到 \(\alpha_{\mathrm{eff}}\approx 1/137\)）
- **α 的“自洽涌现 + 几何约束屏蔽扫描”（toy）**：`math/spectral_graph/run_alpha_emergent_geo_constraint_scan_v0.py`（把“隐藏边自洽生长”接到“几何维度约束屏蔽”路线：可见几何仍用 base torus 的 \(d_s^{geo}\approx3\) 约束；隐藏边按传播子相关 \(|\rho_{ij}|\) 的局部生态规则涌现，同时在更新里加入各向同性/小世界化惩罚项；最终用 \(A_{geo}+A_{hidden}\) 的 Poisson/Coulomb 测量 \(\alpha_{\mathrm{eff}}\)，并用中尺度 plateau/拟合共同打分）

#### 证据文件

| 脚本 | 内容 | 输出 |
|-----|---|---|
| `run_alpha_from_coulomb_3d_v0.py` | 3D torus 上用 Poisson/Coulomb 系数定义并测量 \(\alpha_{\mathrm{eff}}\)（口径管线） | `alpha_from_coulomb_3d_results_v0.json`, `alpha_from_coulomb_3d_plot_v0.png` |
| `run_alpha_screening_hidden_edges_scan_v1.py` | 加入隐藏长程弱边后 \(\alpha_{\mathrm{eff}}\) 的尺度流/屏蔽扫描（含 \(d_s\) 诊断） | `alpha_screening_hidden_edges_scan_results_v1.json`, `alpha_screening_hidden_edges_scan_plot_v1.png` |
| `run_alpha_screening_geo_constraint_scan_v0.py` | **几何维度约束**：用可见几何 \(A_{geo}\) 计算 \(d_s^{geo}\approx3\)；用 \(A_{geo}+A_{hidden}\) 的 Poisson 测量 \(\alpha_{\mathrm{eff}}\) 并扫描屏蔽强度，寻找 plateau \(\approx 1/137\) | `alpha_screening_geo_constraint_scan_results_v0.json`, `alpha_screening_geo_constraint_scan_plot_v0.png` |
| `run_alpha_emergent_hidden_edges_v0.py` | **隐藏关系自洽涌现**：隐藏边按 \(|\rho_{ij}|\) 的局部规则自洽生长（含预算/竞争），再对最终全关系图做 Poisson/Coulomb 测量 \(\alpha_{\mathrm{eff}}\) | `alpha_emergent_hidden_edges_results_v0.json`, `alpha_emergent_hidden_edges_plot_v0.png` |
| `run_alpha_emergent_geo_constraint_scan_v0.py` | **自洽涌现 + 几何约束屏蔽扫描**：隐藏边按 \(|\rho_{ij}|\) 涌现，并加入各向同性/小世界化惩罚；在 \(d_s^{geo}\approx3\) 下扫描 \((\mu, K, w_{\max}, w_{\mathrm{prune}})\) 寻找 \(\alpha_{\mathrm{eff}}\approx 1/137\) 的区域 | `alpha_emergent_geo_constraint_scan_results_v0.json`, `alpha_emergent_geo_constraint_scan_plot_v0.png` |
| `run_alpha_emergent_geo_constraint_scan_v0.py`（tag: `n10_refine_K18`） | **细扫**：固定 \(K=18\)，在 \(\mu\) 与 \(w_{\max}\) 小网格中细扫 | `alpha_emergent_geo_constraint_scan_n10_refine_K18_results_v0.json`, `alpha_emergent_geo_constraint_scan_n10_refine_K18_plot_v0.png` |
| `run_alpha_emergent_geo_constraint_scan_v0.py`（tag: `n10_fine_v1`） | **细扫**：在 \(\mu\in[0.11,0.13]\)、\(w_{\max}\in[0.28,0.32]\)、\(w_{\mathrm{prune}}\in[0.002,0.003]\) 网格内细扫 | `alpha_emergent_geo_constraint_scan_n10_fine_v1_results_v0.json`, `alpha_emergent_geo_constraint_scan_n10_fine_v1_plot_v0.png` |
| `run_alpha_emergent_geo_constraint_scan_v0.py`（tag: `n10_locate_muK_wmax_v1`） | **定位扫**：\(\mu\in\{0.06,0.08,0.10\}\)、\(K\in\{30,40,50\}\)、\(w_{\max}\in\{0.20,0.25,0.30,0.35\}\)（2 seeds） | `alpha_emergent_geo_constraint_scan_n10_locate_muK_wmax_v1_results_v0.json`, `alpha_emergent_geo_constraint_scan_n10_locate_muK_wmax_v1_plot_v0.png` |
| `run_alpha_emergent_geo_constraint_scan_v0.py`（tag: `n10_mu0p08_wmax_verif_s8`） | **复核**：固定 \(\mu=0.08,K=30,w_{\mathrm{prune}}=0.003\)，扫 \(w_{\max}\in\{0.28,0.30,0.32\}\)（8 seeds） | `alpha_emergent_geo_constraint_scan_n10_mu0p08_wmax_verif_s8_results_v0.json`, `alpha_emergent_geo_constraint_scan_n10_mu0p08_wmax_verif_s8_plot_v0.png` |
| `run_alpha_emergent_geo_constraint_scan_v0.py`（tag: `n10_extend_wmax_v1`） | **扩 \(w_{\max}\) 试探**：\(\mu\in\{0.06,0.08,0.10\}\)、\(w_{\max}\in\{0.30,0.50,0.70\}\)（2 seeds） | `alpha_emergent_geo_constraint_scan_n10_extend_wmax_v1_results_v0.json`, `alpha_emergent_geo_constraint_scan_n10_extend_wmax_v1_plot_v0.png` |
| `run_alpha_emergent_geo_constraint_scan_v0.py`（tag: `n10_globalprop_smoke`） | **全局候选提议（长程候选混入）烟雾测试**：启用 `--proposal-global-frac`（2 seeds） | `alpha_emergent_geo_constraint_scan_n10_globalprop_smoke_results_v0.json`, `alpha_emergent_geo_constraint_scan_n10_globalprop_smoke_plot_v0.png` |
| `run_alpha_emergent_geo_constraint_scan_v0.py`（tag: `n10_prop0p2_wmax_zoom_v1`） | **新旋钮：提议覆盖率**：固定 \(\mu=0.08,K=30,w_{\mathrm{prune}}=0.003\)，设置 `propose_rate=0.2`，细扫 \(w_{\max}\in\{0.31,0.32,0.33\}\)（4 seeds） | `alpha_emergent_geo_constraint_scan_n10_prop0p2_wmax_zoom_v1_results_v0.json`, `alpha_emergent_geo_constraint_scan_n10_prop0p2_wmax_zoom_v1_plot_v0.png` |
| `run_alpha_emergent_geo_constraint_scan_v0.py`（tag: `n10_prop0p2_wmax0p32_verif_s10`） | **复核（快速测量口径）**：`propose_rate=0.2`，固定 \(w_{\max}=0.32\)（10 seeds，`n_pairs_per_r=16`） | `alpha_emergent_geo_constraint_scan_n10_prop0p2_wmax0p32_verif_s10_results_v0.json`, `alpha_emergent_geo_constraint_scan_n10_prop0p2_wmax0p32_verif_s10_plot_v0.png` |
| `run_alpha_emergent_geo_constraint_scan_v0.py`（tag: `n10_prop0p2_wmax0p32_verif_s10_npairs64`） | **复核（更低噪声测量）**：同上，但 `n_pairs_per_r=64`（显示测量噪声/口径敏感性） | `alpha_emergent_geo_constraint_scan_n10_prop0p2_wmax0p32_verif_s10_npairs64_results_v0.json`, `alpha_emergent_geo_constraint_scan_n10_prop0p2_wmax0p32_verif_s10_npairs64_plot_v0.png` |
| `run_alpha_emergent_geo_constraint_scan_v0.py`（tag: `n10_prop0p2_wmax0p26_verif_s10_npairs64`） | **当前更稳的证据点候选**：`propose_rate=0.2`，固定 \(w_{\max}=0.26\)（10 seeds，`n_pairs_per_r=64`） | `alpha_emergent_geo_constraint_scan_n10_prop0p2_wmax0p26_verif_s10_npairs64_results_v0.json`, `alpha_emergent_geo_constraint_scan_n10_prop0p2_wmax0p26_verif_s10_npairs64_plot_v0.png` |
| `run_alpha_emergent_geo_constraint_scan_v0.py`（tag: `n10_prop0p2_wmax_zoom_v2_npairs64`） | **定位 \(w_{\max}\) 的甜点区**：`propose_rate=0.2`，扫 \(w_{\max}\in\{0.25,0.26,0.27\}\)（6 seeds，`n_pairs_per_r=64`） | `alpha_emergent_geo_constraint_scan_n10_prop0p2_wmax_zoom_v2_npairs64_results_v0.json`, `alpha_emergent_geo_constraint_scan_n10_prop0p2_wmax_zoom_v2_npairs64_plot_v0.png` |
| `run_alpha_emergent_geo_constraint_scan_v0.py`（tag: `n10_prop0p2_wmax0p255_verif_s10_npairs64`） | **反例/波动记录**：同口径下 \(w_{\max}=0.255\)（10 seeds，`n_pairs_per_r=64`）波动较大 | `alpha_emergent_geo_constraint_scan_n10_prop0p2_wmax0p255_verif_s10_npairs64_results_v0.json`, `alpha_emergent_geo_constraint_scan_n10_prop0p2_wmax0p255_verif_s10_npairs64_plot_v0.png` |
| `run_alpha_emergent_geo_constraint_scan_v0.py`（tag: `n10_wmax_refine_mu0p08_prop0p2_npairs64_v2`） | **加密 \(w_{\max}\)（固定 \(\mu=0.08\)）**：`propose_rate=0.2`，扫 \(w_{\max}\in\{0.259,0.260,0.261\}\)（4 seeds，`n_pairs_per_r=64`） | `alpha_emergent_geo_constraint_scan_n10_wmax_refine_mu0p08_prop0p2_npairs64_v2_results_v0.json`, `alpha_emergent_geo_constraint_scan_n10_wmax_refine_mu0p08_prop0p2_npairs64_v2_plot_v0.png` |
| `run_alpha_emergent_geo_constraint_scan_v0.py`（tag: `n10_prop0p2_wmax0p259_verif_s20_npairs64_v2`） | **更高 seed 复核**：固定 `propose_rate=0.2`、\(w_{\max}=0.259\)（20 seeds，`n_pairs_per_r=64`） | `alpha_emergent_geo_constraint_scan_n10_prop0p2_wmax0p259_verif_s20_npairs64_v2_results_v0.json`, `alpha_emergent_geo_constraint_scan_n10_prop0p2_wmax0p259_verif_s20_npairs64_v2_plot_v0.png` |