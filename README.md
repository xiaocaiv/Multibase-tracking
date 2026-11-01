# Multibase-tracking

毕业课题 ——《基于目标噪声特征关联的多基地定位跟踪技术研究》

当前仓库实现的是课题的基础场景：**双基地机动平台协同跟踪双目标**，纯方位角量测，UKF + JPDA + 三角测量初始化，并以 LOFAR/DEMON 舰船噪声频谱特征辅助关联。

## 目录

```
双站机动跟踪多目标/
├── main.m                      主脚本: 仿真循环 + 评价 + 出图
├── init_parameters.m           全部参数 (仿真/UKF/JPDA/特征/GOSPA)
├── generate_truth.m            真值轨迹 + 频谱模板
├── generate_measurements.m     方位角量测 + 杂波
├── generate_sigma_points.m     UKF: sigma 点
├── ukf_predict.m               UKF: 预测
├── ukf_measurement_predict.m   UKF: 量测预测
├── measurement_bearing.m       方位角量测方程
├── feasible_events.m           JPDA: 可行联合事件
├── jpda_associate.m            JPDA: 边缘关联概率
├── jpda_process_platform.m     单平台 JPDA + UKF 更新流程
├── triangulate_init.m          双平台三角测量起批
├── manage_tracks.m             航迹起始 / 删除
├── platform_ecm_optimize.m     平台机动优化
├── gospa_metric.m              GOSPA 指标 (单时刻, 含分解)
├── evaluate_gospa.m            GOSPA 逐帧评估与汇总
├── test_gospa.m                GOSPA 正确性自检
├── plot_results.m              可视化 + 指标打印
├── sweep_clutter.m             杂波密度蒙特卡洛扫描
├── compare_ship_spectra.m      舰船频谱对比工具
└── *.wav                       舰船噪声样本 (仿真输入数据)
```

## 运行

需要 MATLAB R2024a。所有文件在同一个扁平目录下，需加入 MATLAB 路径。

```bash
# 完整仿真 (含 RMSE 与 GOSPA 汇总、7 张图)
matlab -batch "main"

# GOSPA 实现自检 (819 项断言)
matlab -batch "test_gospa"

# 杂波密度蒙特卡洛扫描 (50 次 × 5 个杂波密度, 约 8 分钟)
matlab -batch "sweep_clutter"
```

依赖：Statistics and Machine Learning Toolbox（`chi2inv`、`poissrnd`）；Signal Processing Toolbox（`spectrogram`、`hilbert`、`butter`、`filtfilt`，仅当 `use_features=true` 且导入 wav 时需要）。

## 评价指标

仓库同时提供两类指标，二者互补，建议在论文中同时报告。

| 指标 | 含义 | 局限 / 特点 |
|---|---|---|
| **RMSE** | 目标位置、速度的均方根误差 | 只统计"已正确关联"的航迹，**忽略漏检目标与虚假航迹**，数值偏乐观 |
| **GOSPA** | 广义最优子模式指派距离 | 同时惩罚定位误差、漏检、虚警，更严格；可分解为定位分量与基数分量 |

### GOSPA 定义

```
GOSPA = [ min_γ ( Σ_(i,j)∈γ d(x_i, y_j)^p + (c^p/α)(n + m − 2|γ|) ) ]^(1/p)
```

- `X`：真值点集（n 个目标）；`Y`：估计点集（m 条确认航迹）
- `γ`：X 到 Y 的**部分单射**（每个点最多匹配一次），在所有可能指派中取最优
- 未被匹配的点（漏检或虚警）每个贡献 `c^p/α`
- 可加分解：`total^p = loc^p + card^p`，用于区分"跟踪不准"与"跟踪漏/多"

### 参数取值

| 参数 | 默认 | 说明 |
|---|---|---|
| `gospa_c` | 100 m | 截断距离。建议取典型位置误差的 3~5 倍；本场景位置误差约 10~20 m、目标间距 ≥ 569 m，故取 100 m |
| `gospa_p` | 2 | 距离阶数，p=2 时 GOSPA 满足度量公理 |
| `gospa_alpha` | 2 | 基数分配系数，漏检与虚警各承担 `c^p` 的一半 |
| `gospa_use_confirmed` | true | 是否只用已确认航迹参与评估 |
| `gospa_conv_step` | 20 | 计算"收敛段均值"的起始步（首帧必然全漏检，不宜计入） |

`evaluate_gospa` 返回逐帧序列与汇总统计，论文中通常引用 `mean_total_conv`（收敛段均值）。

## 默认场景结果

`rng(42)`、`use_features = false`：

| 指标 | 数值 |
|---|---|
| 位置 RMSE（目标1 / 目标2） | 9.38 m / 17.26 m |
| 速度 RMSE — Vx（目标1 / 目标2） | 1.10 m/s / 1.57 m/s |
| 速度 RMSE — Vy（目标1 / 目标2） | 0.71 m/s / 1.23 m/s |
| GOSPA（收敛段均值 / 末帧） | 119.00 m / 173.87 m |
| 平均漏检目标数 | 0.07 个/帧 |
| 平均虚假航迹数 | 2.74 条/帧 |

两者差距主要来自**虚假航迹**：RMSE 看不到它，GOSPA 会直接计入成本。这也是引入 GOSPA 的主要理由。

## 参考文献

- A. S. Rahmani, M. Vihola, "A Two-Level Implementation of GOSPA Metric for Multi-Target Tracking Performance Evaluation," IEEE FUSION, 2019.
