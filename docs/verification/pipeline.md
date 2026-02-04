# 验证管线模板（统一格式）

目标：把“框架性解释”变成**可复现、可扩点、可证伪**的研究管线。

本仓库优先采用同一种套路来写每个可验证主张：

```
主张（公式/方向/假设）
  ↓
明确数据口径（变量定义、尺度、误差项）
  ↓
可复现数据表（CSV/JSON）
  ↓
可复现脚本（fit + sanity-check + plot + export）
  ↓
输出（results.json + plot.png）
  ↓
证伪条件（失败即改模型）
```

---

## 1. 最小文件结构（推荐）

为每个“可验证命题”新建一个目录（放在 `math/` 或 `docs/physics/` 对应子目录下均可，但建议 `math/` 放计算管线）：

```
<topic>/
├── README.md                   # 目的、口径、运行方式、局限、下一步
├── data_v0.csv                 # 明确点集（观测/预测/推断分开标注）
├── run_<topic>_v0.py            # 读取数据 → 拟合/统计 → 出图 → 导出 JSON
├── results_v0.json              # 机器可读结果（由脚本生成）
└── plot_v0.png                  # 图（由脚本生成）
```

---

## 2. README 必须包含的 6 件事

- **问题与主张**：一句话 + 公式（变量含义）
- **数据口径**：δ/尺度/红移段/选择规则/误差缺失如何处理
- **拟合/统计方法**：最小二乘、贝叶斯、LOO 等（至少给出 chi²/dof 或等价指标）
- **依赖与运行方式**：推荐统一用仓库根目录的 `requirements.txt`（`pip install -r requirements.txt`），并给出从根目录可直接运行的命令
- **输出说明**：`results_v0.json` 的字段含义
- **证伪条件**：出现什么观测模式就必须推翻/修正模型
- **v0 局限 & v1 路线**：缺点与下一步扩点策略

---

## 3. 输出与版本管理（强制）

为保证“可复现 + 可审核”，输出与版本管理采用统一规则：

- **输出命名**：`results_v*.json` 与 `plot_v*.png`
- **应提交到 git**：
  - `results_v*.json`（机器可读指标与配置）
  - `plot_v*.png`（可视化证据）
  - 数据表（CSV/JSON）与脚本
- **应忽略（不进 git）**：
  - 下载缓存/临时文件：`_cache/`、`*.npy`、`*.npz`
  - 临时演示截图（如终端渲染图）
  - 其他一次性/可重建中间产物

这些规则在 `/.gitignore` 中统一维护。

---

## 4. results JSON 最小结构（统一字段）

每个管线的 `results_v*.json` 至少包含以下四块：

- `metadata`：时间戳、脚本名、数据源等
- `config`：关键参数（运行配置）
- `metrics` 或 `final`：核心指标摘要
- `outputs`：输出文件路径（json/plot）

**最小示例：**

```json
{
  "metadata": {
    "timestamp": "2026-02-01T12:00:00",
    "script": "run_example_v1.py",
    "data_file": "path/to/data_v1.csv"
  },
  "config": {
    "seed": 42,
    "n_steps": 200
  },
  "metrics": {
    "beta": 0.16,
    "chi2_dof": 0.9
  },
  "outputs": {
    "json": "path/to/results_v1.json",
    "plot": "path/to/plot_v1.png"
  },
  "results": [
    { "item": "detail-1" }
  ]
}
```

> 说明：`metrics`/`final` 用于“摘要指标”；`results` 用于详细记录（可选但推荐保留）。

---

## 5. 数据表（CSV）推荐字段

不同主题字段会不同，但建议至少有：

| 字段 | 含义 |
|---|---|
| `id` | 行唯一标识 |
| `category` | 环境/类型（void/cluster/global…） |
| `data_kind` | `observed` / `forecast` / `inferred` |
| `use_for_fit` | 是否参与拟合（预测点必须为 false） |
| `x`, `x_err` | 自变量与误差（如 δ） |
| `y`, `y_err` | 因变量与误差（如 H0） |
| `source` | 文献/数据来源 |
| `notes` | 假设与口径说明 |

---

## 6. 两个现成范例

- **哈勃张力环境扫描**：`math/hubble_env_scan/`
  - 数据表 + 拟合 + LOOCV + JSON/图输出 + 证伪条件
- **维度涌现（谱维度）**：`math/spectral_graph/`
  - 图生成 + 谱估计 + 对照组校准 + JSON/图输出
- **黑洞观测桥接（方法论）**：`docs/verification/black_hole_observational_bridge.md`
  - 把“同构直觉”落到“前向模型 + 反演 + 多观测约束 + 证伪条件”的可复现套路
- **已对齐输出结构**：
  - `math/hubble_env_scan/run_hubble_env_scan_v1.py`
  - `math/spectral_graph/run_dimension_emergence_v0.py`

---

## 7. 最重要的原则：证伪优先

“能解释很多现象”不是验证；**能明确失败的条件**才是理论。

写任何预测时，请先写：

- 零假设是什么？
- 什么数据模式会让模型必然失败？
- 失败后该升级哪一层（口径、模型形式、或者框架本身）？

