# 哈勃张力（模型 v2）点集口径审计报告（v0）

_生成时间：2026-02-02T20:30:54_

## 摘要

- 总行数：**24**
- 观测：9；预测：6；推断：9
- 默认主线拟合点（direct+observed+use_for_fit）：**3**

按 δ 口径轨道（track）计数：
- **direct**：14
- **none**：2
- **proxy**：4
- **sim**：4

按用途角色（role）计数：
- **forecast_marker**：6
- **optional_fit_direct**：2
- **overlay_only**：5
- **primary_fit**：3
- **proxy_observed_overlay**：4
- **reference_fit_sim**：4

## 关键口径提醒（v2）

- **direct**：文献/模型给出的 δ（主线拟合用这个口径）。
- **sim**：模拟定义的 δ_den（例如 10 Mpc/h 内），只做外部参照；不要与 direct 混合拟合。
- **proxy**：2M++ 的 δ_g*（带平滑尺度），只做公共参照/overlay；不要与 direct 混合拟合。

## 点集清单（按 role）

### primary_fit

- `cmb_global`：δ=0.0, H0=67.36（delta_kind=direct；无额外标注）
- `kbc_void_tf`：δ=-0.46, H0=72.1（delta_kind=direct；shell=shell_40_200）
- `sh0es_local`：δ=-0.46, H0=73.17（delta_kind=direct；shell=shell_40_200）

### optional_fit_direct

- `kbc_void_model`：δ=-0.3, H0=71.06（delta_kind=direct；无额外标注）
- `shanks_void_model`：δ=-0.2, H0=69.72（delta_kind=direct；无额外标注）

### reference_fit_sim

- `gavas500_denden_0p7`：δ=0.7, H0=66.01（delta_kind=gavas_delta_den_10mpch；scale=10.0 Mpc/h）
- `gavas500_denden_1p4`：δ=1.4, H0=59.28（delta_kind=gavas_delta_den_10mpch；scale=10.0 Mpc/h）
- `gavas500_denden_2p2`：δ=2.2, H0=55.24（delta_kind=gavas_delta_den_10mpch；scale=10.0 Mpc/h）
- `gavas500_denden_3p6`：δ=3.6, H0=48.5（delta_kind=gavas_delta_den_10mpch；scale=10.0 Mpc/h）

### proxy_observed_overlay

- `coma_cz_over_D_2mpp20`：δ=0.4322764151332613, H0=74.42（delta_kind=delta_g_star_2mpp；scale=20.0 Mpc/h）
- `coma_cz_over_D_2mpp50`：δ=0.07677854526028692, H0=74.42（delta_kind=delta_g_star_2mpp；scale=50.0 Mpc/h）
- `coma_cz_over_D_2mpp100`：δ=0.012647045225691635, H0=74.42（delta_kind=delta_g_star_2mpp；scale=100.0 Mpc/h）
- `coma_benisty2025_hubbleflow_2mpp20`：δ=0.4322764151332613, H0=73.1（delta_kind=delta_g_star_2mpp；scale=20.0 Mpc/h）

### forecast_marker

- `cold_spot_min`：δ=-0.7, H0=75.4（delta_kind=direct；无额外标注）
- `cold_spot_max`：δ=-0.6, H0=74.3（delta_kind=direct；无额外标注）
- `cluster_min`：δ=0.3, H0=63.9（delta_kind=direct；无额外标注）
- `cluster_max`：δ=0.5, H0=61.6（delta_kind=direct；无额外标注）
- `coma_min`：δ=1.0, H0=55.8（delta_kind=direct；无额外标注）
- `coma_max`：δ=2.0, H0=44.3（delta_kind=direct；无额外标注）

### overlay_only

- `cf4_void_exp`：δ=-0.46, H0=72.08（delta_kind=direct；无额外标注）
- `cf4_void_gauss`：δ=-0.46, H0=70.4（delta_kind=direct；无额外标注）
- `cf4_void_mb`：δ=-0.46, H0=70.18（delta_kind=direct；无额外标注）
- `coma_cz_over_D`：δ=None, H0=74.42（delta_kind=direct；无额外标注）
- `coma_desi_fp`：δ=None, H0=76.5（delta_kind=direct；无额外标注）

## 壳层映射说明

- `kbc_void_tf`：scale_mpc='40-300' → `shell_40_200`（rmax_capped_to_200）
- `sh0es_local`：scale_mpc='40-300' → `shell_40_200`（rmax_capped_to_200）
