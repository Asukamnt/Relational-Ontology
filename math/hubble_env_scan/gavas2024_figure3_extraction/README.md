# Gavas et al. 2024 Figure 3 数据提取

## 目的

从 [Gavas et al. 2024 (arXiv:2407.10139)](https://arxiv.org/abs/2407.10139) 的 Figure 3 提取 $\delta_{den}$ vs $\delta_H$ 散点数据，用于独立拟合 $\beta$。

## Figure 3 说明

Figure 3 包含 6 个子图，分为 3 行（对应 150/500/1000 Mpc/h 模拟盒子）：
- **左列**：10 Mpc/h 厚切片中观测者的空间分布（红=H₀偏低，绿=H₀偏高）
- **右列**：$\delta_{den}$（10 Mpc/h 内局部过密度）vs $\delta_H$（H₀ 偏差）的散点图 ← **这是我们需要的**

## 步骤 1：下载图像

从 arXiv HTML 版本直接下载右列散点图（Figure 3 右栏）：

**直接下载链接**：
- [150 Mpc/h 盒子](https://arxiv.org/html/2407.10139v4/extracted/6143171/images/rdmk150_1024halo_r_5p.png)
- [500 Mpc/h 盒子](https://arxiv.org/html/2407.10139v4/extracted/6143171/images/rdmk500_1024halo_r_5p.png)
- [1000 Mpc/h 盒子](https://arxiv.org/html/2407.10139v4/extracted/6143171/images/rdmk1000_1024halo_r_5p.png)

将它们保存到本目录下的 `images/` 文件夹：
```
images/
  rdmk150_1024halo_r_5p.png
  rdmk500_1024halo_r_5p.png
  rdmk1000_1024halo_r_5p.png
```

## 步骤 2：使用 WebPlotDigitizer 提取数据

推荐使用免费在线工具 [WebPlotDigitizer](https://apps.automeris.io/wpd/)：

1. 上传图像
2. 选择 "2D (X-Y) Plot"
3. 校准坐标轴（根据图中的刻度）
4. 手动或自动提取散点坐标
5. 导出为 CSV

### 坐标轴参考（从论文读取）

根据 Figure 3 右列：
- **X 轴**：$\delta_H = (H_L - H_G)/H_G$，范围约 [-0.5, 0.7]
- **Y 轴**：$\delta_{den}$（局部过密度），范围约 [-1, 30]（对数刻度可能）

## 步骤 3：保存提取的数据

将提取的数据保存为 CSV，格式如下：

```csv
delta_H,delta_den,source_box_mpc_h,halo_type
-0.05,2.3,150,all
0.02,-0.5,150,all
...
```

保存到本目录下：
- `gavas2024_fig3_150mpc.csv`
- `gavas2024_fig3_500mpc.csv`
- `gavas2024_fig3_1000mpc.csv`

## 步骤 4：运行分析脚本

提取完成后，运行：

```bash
python analyze_gavas2024_fig3.py
```

## 论文报告的统计量（备用验证）

如果无法精确提取所有散点，可以用论文报告的相关系数来验证方向：

| Box Size | Pearson (all) | Spearman (all) | Pearson (MWh) | Spearman (MWh) |
|----------|---------------|----------------|---------------|----------------|
| 150 Mpc/h | ~-0.27 to -0.35 | ~-0.28 to -0.34 | similar | similar |
| 500 Mpc/h | similar | similar | similar | similar |
| 1000 Mpc/h | similar | similar | similar | similar |

（MWh = Milky Way-sized halos）

## 与我们模型的对应关系

我们的模型：$H_0(\delta) = H_{0,global} \times (1 - \beta \delta)$

等价于：$\delta_H = (H_L - H_G)/H_G = -\beta \delta_{den}$

因此：
- 从 Gavas 散点图的 **斜率** 可以直接读出 $-\beta$
- 负相关 → $\beta > 0$（符合预期）

## 参考文献

- Gavas, S., Bagla, J. S., & Khandai, N. (2024). *Dispersion in the Hubble-Lemaître constant measurements from gravitational clustering*. arXiv:2407.10139
