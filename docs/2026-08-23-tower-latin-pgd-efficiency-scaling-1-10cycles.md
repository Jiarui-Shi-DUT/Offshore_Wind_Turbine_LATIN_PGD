# 海上风机塔筒 LATIN-PGD 计算效率与规模扩展阶段总结：1–10 周期 wall-time 验证

**日期：** 2026-08-23  
**项目：** Offshore_Wind_Turbine_LATIN_PGD  
**分支：** `feature/offshore-wind-turbine-tower-fatigue`  
**阶段：** FOM-3D  
**当前状态：** 1、2、5、10 周期基础 wall-time scaling 已完成；成本来源诊断尚未开展

---

## 1. 本文档的起点与范围

本文档承接上一份正式阶段总结：

`docs/2026-08-23-tower-matched-asymmetric-latin-multicycle-validation-1-10cycles.md`

上一份总结已经正式记录：

> **FOM-3C：海上风机塔筒 LATIN-PGD 非对称循环 1–10 周期匹配精度验证 PASS。**

其对应 Git checkpoint 为：

`9822737  docs: summarize matched asymmetric LATIN multicycle validation`

在该 checkpoint 中已经确认：

- 1、2、5、10 周期 matched asymmetric LATIN-PGD 均稳定收敛；
- 10 周期应力和弹性应变全历史误差约为 `0.14%`；
- 10 周期主要内部变量全历史误差约为 `1%`；
- cycle 10 最大塑性应变误差约为 `0.734%`；
- cycle 10 最大损伤误差约为 `0.688%`；
- 时间点数由 41 增长到 401 时，PGD rank 仅由 11 增长到 21；
- 1–10 周期范围内没有观察到明显的误差累积失控。

但是上一阶段尚未回答一个关键问题：

> **低 rank 是否真正转化为了 wall-time 层面的计算效率优势。**

因此，本阶段从 FOM-3C 的精度验证转入 FOM-3D 的计算效率和规模扩展验证。

本文档只总结上一份正式 Markdown 之后完成的工作，重点包括：

1. 计时口径和公平比较原则；
2. 单周期计时脚本建立与 API 校核；
3. 1、2、5、10 周期的三次独立 wall-time 重复试验；
4. FOM 与 LATIN-PGD 的规模扩展趋势；
5. PGD rank 与实际 wall time 之间出现的关键差异；
6. 当前实现中效率 crossover 的初步定位；
7. 为什么此时不应直接继续到 100 周期；
8. 下一阶段 wall-time 成本来源诊断的研究入口。

本阶段仍然没有修改任何 `latin/` 核心代码。

---

## 2. 为什么在 10 周期精度验证之后没有直接计算 100 周期

上一阶段已经得到一个表面上十分有利的低秩现象：

| 周期数 | Nt | PGD rank |
|---:|---:|---:|
| 1 | 41 | 11 |
| 2 | 81 | 13 |
| 5 | 201 | 17 |
| 10 | 401 | 21 |

从 1 周期增加到 10 周期：

- `Nt` 从 41 增加到 401，约增加 `9.78 倍`；
- PGD rank 从 11 增加到 21，仅增加 `1.91 倍`。

如果只观察 rank，很容易形成一个过早判断：

> 时间域增加近 10 倍，而 PGD rank 只增加不到 2 倍，因此 LATIN-PGD 的计算成本也应明显优于 FOM。

但这种判断并不充分。

LATIN-PGD 的实际 wall time 不只由最终 rank 决定，还包括：

- Local stage；
- search-direction construction；
- fixed-basis temporal update；
- Trial A evaluation；
- PGD enrichment fixed-point；
- Trial B evaluation；
- Global finishing / reconstruction；
- outer LATIN iterations；
- 不同 trial 中重复进行的全时间域运算。

因此，本阶段决定先回答：

> **当前代码实现的实际 wall time 如何随周期数和 Nt 增长。**

在这个问题没有回答之前，直接运行 100 周期只会得到一个更大的总时间，却无法解释时间花在了哪里。

所以本阶段采用的原则是：

> **先做 1、2、5、10 周期受控 scaling，再决定是否进入 100 周期。**

---

## 3. 计时 benchmark 与上一阶段保持完全一致

本阶段没有重新定义新的物理 benchmark，而是直接沿用已经完成精度验证的 matched asymmetric benchmark。

### 3.1 荷载合同

- `Fmax = +1.0 MN`
- `Fmin = -0.5 MN`
- `R = -0.5`
- `Fmean = +0.25 MN`
- 周期 `T = 10 s`
- 每周期 40 个增量
- 每周期完整闭区间 41 个时间点

周期路径为：

`+0.25 -> +1.00 -> +0.25 -> -0.50 -> +0.25 MN`

FOM 继续使用独立的：

`0 -> +0.25 MN`

弹性预加载。

LATIN 继续从已验证等价的 `+0.25 MN` 弹性初始化状态直接开始周期历史。

### 3.2 空间离散合同

所有 timing benchmark 均保持：

- 梁柱单元数：10；
- 每单元 Gauss 点数：2；
- 每 Gauss 点环向纤维数：16；
- 径向纤维层数：1；
- 总材料点数：`Nq = 320`；
- 相同 `MaterialParameters`；
- 相同塔筒几何；
- 相同 `ViscoplasticDamageTowerSystem2D`；
- 相同 equilibrium operator；
- 相同 `stress_to_force_factor`。

### 3.3 LATIN 求解合同

继续采用上一阶段已经验证的参数：

- `spatial_strategy="residual_ls"`
- `tolerance=1.0e-5`
- `mode_significance_tolerance=0.0`
- `acceptance_tolerance=0.0`
- `max_fixed_point_iterations=200`

因此，本阶段 wall-time 变化不能归因于 benchmark 或求解参数改变。

---

## 4. 计时口径的确定

### 4.1 为什么不能只比较 LATIN solver 时间与 FOM 总时间

如果 LATIN 只计：

`solve_tower_latin_pgd(...)`

而 FOM 计：

`run_nonlinear_asymmetric_analysis(...)`

那么二者的计时边界并不完全一致。

这会人为遗漏 LATIN 的：

- geometry 和 mesh 构建；
- tower material-point system 构建；
- equilibrium operator 构建；
- load vectors 构建；
- elastic initialization。

因此，本阶段将 LATIN 拆为：

- `LATIN setup wall time`
- `LATIN solver wall time`
- `LATIN total wall time = setup + solver`

正式用于当前阶段对比的是：

`FOM analysis wall time / LATIN total wall time`

定义：

$$ S_{\mathrm{pilot}}=\frac{t_{\mathrm{FOM}}}{t_{\mathrm{LATIN,total}}} $$

其中：

- `S_pilot` 大于 1：当前测量中 LATIN 更快；
- `S_pilot` 等于 1：两者 wall time 接近；
- `S_pilot` 小于 1：当前测量中 FOM 更快。

### 4.2 FOM 的计时边界

FOM 计时包围：

`run_nonlinear_asymmetric_analysis(...)`

该函数当前包括：

- 塔筒 mesh 与 system 构建；
- 非对称工况的 10 步弹性预加载；
- 后续全部周期非线性 load-step 求解；
- 位移、应变、应力、状态历史快照；
- response 组装。

因此该时间代表当前正式 FOM driver 的实际 wall time。

### 4.3 LATIN 的计时边界

LATIN setup 包括：

- `create_tower_geometry(...)`
- `create_uniform_vertical_tower_mesh(...)`
- `ViscoplasticDamageTowerSystem2D(...)`
- `build_tower_equilibrium_operator(...)`
- 周期 `load_vectors`
- `compute_tower_elastic_initialization(...)`

LATIN solver 计时包围：

`solve_tower_latin_pgd(...)`

因此当前阶段使用：

`LATIN total = LATIN setup + LATIN solver`

### 4.4 当前计时口径的性质

这套口径适合回答：

> **当前两个正式代码路径在同一台机器上的实际 wall-time 相对量级。**

但它仍然不是最终论文级性能 benchmark，因为尚未系统记录：

- CPU 型号；
- 核数和线程数；
- BLAS 实现；
- Python / NumPy 版本；
- 内存峰值；
- CPU frequency / thermal state；
- FOM 与 LATIN 执行顺序是否产生系统偏差。

因此当前所有 ratio 均应理解为：

> **当前实现、当前机器、当前 benchmark 下的受控数值诊断结果。**

---

## 5. 单周期计时脚本的建立与 API 校核

本阶段首先建立：

`tower_asymmetric_1cycle_efficiency_pilot.py`

随后经过 API 校核形成：

`tower_asymmetric_1cycle_efficiency_pilot_v2.py`

最终可完整运行的版本为：

`tower_asymmetric_1cycle_efficiency_pilot_v3.py`

这些均为本地 diagnostic，不修改 `latin/` 核心代码。

### 5.1 第一个 API 问题：初始化平衡残差属性

最初脚本错误使用：

`initialization.max_free_equilibrium_residual`

当前 `TowerElasticInitialization` 实际公开属性为：

`initialization.maximum_free_equilibrium_residual`

因此最初两次运行都在 FOM 已完成、LATIN setup 已完成之后，于打印初始化残差时退出。

这两次 FOM 时间分别约为：

- `37.069276 s`
- `36.471076 s`

由于同一进程内 LATIN solver 没有完成，因此这两次运行没有形成成对的 FOM/LATIN timing sample。

所以：

> **这两次失败运行没有纳入后续统计。**

### 5.2 第二个 API 问题：LATIN state 时间点属性

脚本最初计划使用：

`initialization.state.n_time_steps`

当前 `LatinStateTower` 正确属性为：

`initialization.state.n_time`

该问题在重新运行前同步修正。

### 5.3 第三个 API 问题：PGD basis rank 属性

v2 已经成功完成：

- FOM；
- LATIN setup；
- LATIN solver。

而且已经打印：

- `converged = True`
- `iterations = 18`
- `trial evaluations = 29`

但在打印 PGD rank 时错误使用：

`latin.basis.rank`

当前 `PGDBasisTower` 的正确 reduced-basis dimension 属性为：

`latin.basis.n_modes`

因此该次运行仍然没有打印最终 wall-time summary。

### 5.4 同步校核的结果属性

为了避免再次在求解完成后因为 diagnostic 输出失败，本阶段同时校核并统一采用：

- `latin.termination_reason.value`
- `latin.basis.n_modes`
- `latin.total_modes_added`
- `latin.final_indicator`
- `latin.attempted_iterations`

最终形成的 v3 脚本可以完整结束并输出全部 timing 数据。

### 5.5 这一轮脚本调试的性质

需要明确：

> 以上问题都只是外部 diagnostic 脚本调用了错误的公开属性名。

它们不是：

- FOM 求解错误；
- LATIN 求解错误；
- material model 错误；
- convergence 错误；
- PGD basis 错误。

因此没有修改任何核心算法。

---

## 6. FOM-3D-1：1-cycle wall-time baseline

### 6.1 数值一致性

三个有效的 1-cycle timing run 中，LATIN 数值行为完全一致：

- `Nt = 41`
- `Nq = 320`
- `termination_reason = converged`
- `converged = True`
- `failure_reason = None`
- `iterations = 18`
- `attempted = 18`
- `trial evaluations = 29`
- `PGD rank = 11`
- `modes added = 11`
- `final xi = 7.918424536257e-06`

这与上一阶段的单周期 matched asymmetric 精度验证完全一致。

因此：

> **加入 wall-time 计时没有改变求解路径和收敛结果。**

### 6.2 三次独立运行

| Run | FOM wall time / s | LATIN setup / s | LATIN solver / s | LATIN total / s | FOM/LATIN |
|---:|---:|---:|---:|---:|---:|
| 1 | 36.713106 | 0.134500 | 23.415871 | 23.550371 | 1.558918 |
| 2 | 36.697423 | 0.141240 | 23.300499 | 23.441739 | 1.565474 |
| 3 | 38.015724 | 0.150231 | 24.036705 | 24.186936 | 1.571746 |

平均值：

- FOM：`37.142084 s`
- LATIN setup：`0.141990 s`
- LATIN solver：`23.584358 s`
- LATIN total：`23.726349 s`
- 样本 ratio 平均：`1.565379`

### 6.3 重复性

| 指标 | 标准差 | CV |
|---|---:|---:|
| FOM wall time | 0.756635 s | 2.04% |
| LATIN total | 0.402561 s | 1.70% |
| FOM/LATIN ratio | 0.006415 | 0.41% |

ratio 的离散远小于单独 wall time 的离散。

### 6.4 单周期阶段判断

按平均时间：

- LATIN total 约为 FOM 的 `63.9%`；
- 当前测量下 LATIN wall time 约减少 `36.1%`。

因此可以记录：

> **FOM-3D-1：1-cycle timing baseline repeatability PASS。**

但该结论不应写成普适的“LATIN 加速 1.56 倍”。

正确表述应是：

> 在当前 benchmark、当前实现和当前运行环境下，单周期 LATIN total wall time 稳定低于 FOM wall time，当前平均 FOM/LATIN ratio 约为 1.565。

---

## 7. 统一 scaling 脚本

在完成单周期 pilot 后，为避免分别维护 2、5、10 周期脚本，本阶段建立统一 diagnostic：

`tower_asymmetric_efficiency_scaling_pilot.py`

通过参数：

`--cycles`

选择：

- 1 cycle；
- 2 cycles；
- 5 cycles；
- 10 cycles。

例如：

`python tower_asymmetric_efficiency_scaling_pilot.py --cycles 5`

统一脚本保持：

- 相同 benchmark；
- 相同 timing boundary；
- 相同 LATIN 参数；
- 相同输出字段；
- 相同 solver diagnostics。

---

## 8. FOM-3D-2：2-cycle wall-time baseline

### 8.1 数值一致性

三次 2-cycle timing run 均得到：

- `Nt = 81`
- `Nq = 320`
- `converged = True`
- `iterations = 23`
- `attempted = 23`
- `trial evaluations = 36`
- `PGD rank = 13`
- `modes added = 13`
- `final xi = 8.893567635885e-06`

这与上一阶段 2-cycle 精度验证完全一致。

### 8.2 三次独立运行

| Run | FOM wall time / s | LATIN setup / s | LATIN solver / s | LATIN total / s | FOM/LATIN |
|---:|---:|---:|---:|---:|---:|
| 1 | 69.958600 | 0.256368 | 59.110212 | 59.366580 | 1.178417 |
| 2 | 70.336209 | 0.258324 | 59.898090 | 60.156414 | 1.169222 |
| 3 | 67.621059 | 0.237522 | 57.306900 | 57.544422 | 1.175111 |

平均值：

- FOM：`69.305289 s`
- LATIN setup：`0.250738 s`
- LATIN solver：`58.771734 s`
- LATIN total：`59.022472 s`
- 样本 ratio 平均：`1.174250`

### 8.3 重复性

| 指标 | 标准差 | CV |
|---|---:|---:|
| FOM wall time | 1.470755 s | 2.12% |
| LATIN total | 1.339565 s | 2.27% |
| FOM/LATIN ratio | 0.004658 | 0.40% |

### 8.4 两周期阶段判断

按平均时间：

- LATIN total 约为 FOM 的 `85.2%`；
- 当前测量下 LATIN wall time 约减少 `14.8%`。

因此：

> **FOM-3D-2：2-cycle timing baseline PASS。**

同时，ratio 已从 1-cycle 的约 `1.565` 降低到约 `1.174`。

---

## 9. FOM-3D-3：5-cycle wall-time baseline

### 9.1 数值一致性

三次 5-cycle timing run 均得到：

- `Nt = 201`
- `Nq = 320`
- `converged = True`
- `iterations = 33`
- `attempted = 33`
- `trial evaluations = 50`
- `PGD rank = 17`
- `modes added = 17`
- `final xi = 9.560582155456e-06`

### 9.2 三次独立运行

| Run | FOM wall time / s | LATIN setup / s | LATIN solver / s | LATIN total / s | FOM/LATIN |
|---:|---:|---:|---:|---:|---:|
| 1 | 157.182482 | 0.534686 | 200.070046 | 200.604732 | 0.783543 |
| 2 | 162.886189 | 0.567093 | 208.220647 | 208.787740 | 0.780152 |
| 3 | 163.749973 | 0.539853 | 206.159463 | 206.699316 | 0.792213 |

平均值：

- FOM：`161.272881 s`
- LATIN setup：`0.547211 s`
- LATIN solver：`204.816719 s`
- LATIN total：`205.363929 s`
- 样本 ratio 平均：`0.785303`

### 9.3 重复性

| 指标 | 标准差 | CV |
|---|---:|---:|
| FOM wall time | 3.568621 s | 2.21% |
| LATIN total | 4.251805 s | 2.07% |
| FOM/LATIN ratio | 0.006220 | 0.79% |

### 9.4 五周期阶段判断

按平均时间：

- LATIN total 约为 FOM 的 `127.3%`；
- LATIN 比 FOM 慢约 `27.3%`。

因此：

> **FOM-3D-3：5-cycle timing baseline PASS。**

并且：

> **5-cycle 已经稳定出现 LATIN wall time 高于 FOM 的结果。**

---

## 10. FOM-3D-4：10-cycle wall-time baseline

### 10.1 数值一致性

三次 10-cycle timing run 均得到：

- `Nt = 401`
- `Nq = 320`
- `converged = True`
- `iterations = 39`
- `attempted = 39`
- `trial evaluations = 60`
- `PGD rank = 21`
- `modes added = 21`
- `final xi = 8.941607234831e-06`

### 10.2 三次独立运行

| Run | FOM wall time / s | LATIN setup / s | LATIN solver / s | LATIN total / s | FOM/LATIN |
|---:|---:|---:|---:|---:|---:|
| 1 | 307.974583 | 0.990150 | 462.872319 | 463.862469 | 0.663935 |
| 2 | 317.447711 | 1.031023 | 476.065041 | 477.096064 | 0.665375 |
| 3 | 298.250990 | 1.029030 | 460.097022 | 461.126052 | 0.646788 |

平均值：

- FOM：`307.891095 s`
- LATIN setup：`1.016734 s`
- LATIN solver：`466.344794 s`
- LATIN total：`467.361528 s`
- 样本 ratio 平均：`0.658699`

### 10.3 重复性

| 指标 | 标准差 | CV |
|---|---:|---:|
| FOM wall time | 9.598633 s | 3.12% |
| LATIN total | 8.540661 s | 1.83% |
| FOM/LATIN ratio | 0.010341 | 1.57% |

### 10.4 十周期阶段判断

按平均时间：

- LATIN total 约为 FOM 的 `151.8%`；
- LATIN 比 FOM 慢约 `51.8%`。

因此：

> **FOM-3D-4：10-cycle timing baseline PASS。**

---

## 11. 1、2、5、10 周期完整 scaling 汇总

### 11.1 求解规模与算法状态

| 周期数 | Nt | PGD rank | LATIN iterations | trial evaluations | final xi |
|---:|---:|---:|---:|---:|---:|
| 1 | 41 | 11 | 18 | 29 | `7.918424536257e-06` |
| 2 | 81 | 13 | 23 | 36 | `8.893567635885e-06` |
| 5 | 201 | 17 | 33 | 50 | `9.560582155456e-06` |
| 10 | 401 | 21 | 39 | 60 | `8.941607234831e-06` |

全部规模点均稳定收敛。

### 11.2 平均 wall-time

| 周期数 | Nt | FOM 平均 / s | LATIN setup 平均 / s | LATIN solver 平均 / s | LATIN total 平均 / s |
|---:|---:|---:|---:|---:|---:|
| 1 | 41 | 37.142084 | 0.141990 | 23.584358 | 23.726349 |
| 2 | 81 | 69.305289 | 0.250738 | 58.771734 | 59.022472 |
| 5 | 201 | 161.272881 | 0.547211 | 204.816719 | 205.363929 |
| 10 | 401 | 307.891095 | 1.016734 | 466.344794 | 467.361528 |

### 11.3 平均 FOM/LATIN ratio

| 周期数 | 平均 FOM/LATIN ratio | 当前 wall-time 关系 |
|---:|---:|---|
| 1 | 1.565379 | LATIN 更快 |
| 2 | 1.174250 | LATIN 更快，但优势明显缩小 |
| 5 | 0.785303 | FOM 更快 |
| 10 | 0.658699 | FOM 明显更快 |

---

## 12. 关键发现一：当前测试点将效率 crossover 夹在 2–5 周期之间

在当前离散测试点上：

- 2-cycle：ratio 约 `1.174`
- 5-cycle：ratio 约 `0.785`

因此可以确认：

> **当前 benchmark 的实际 wall-time 优势在 2-cycle 与 5-cycle 两个已测试规模点之间发生反转。**

需要特别注意：

> 当前数据只能说明 crossover 被 2 和 5 两个测试点夹住，不能据此声称精确 crossover cycle 已经确定。

如果未来需要精确定位，可增加 3-cycle 与 4-cycle。

但当前更重要的问题不是精确 crossover cycle，而是：

> **为什么 LATIN 的成本增长速度会超过 FOM。**

---

## 13. 关键发现二：低 PGD rank 不等于低 wall time

从 1 cycle 到 10 cycles：

- `Nt`：41 -> 401，增长约 `9.78 倍`
- PGD rank：11 -> 21，增长约 `1.91 倍`
- LATIN iterations：18 -> 39，增长约 `2.17 倍`
- trial evaluations：29 -> 60，增长约 `2.07 倍`

但是平均 wall time：

- FOM：`37.142084 -> 307.891095 s`
- LATIN：`23.726349 -> 467.361528 s`

对应：

- FOM wall time 增长约 `8.29 倍`
- LATIN wall time 增长约 `19.70 倍`

因此：

> **低秩表示能力已经得到数值支持，但低秩表示尚未转化为当前实现的 wall-time 优势。**

这是本阶段最重要的研究认识之一。

---

## 14. 当前区间的经验 scaling exponent

为了对 1-cycle 到 10-cycle 的总体增长速度做描述性比较，可用：

`time ~ Nt^p`

对首尾两点做经验指数估计：

$$ p=\frac{\ln(t_{10}/t_1)}{\ln(N_{t,10}/N_{t,1})} $$

得到约：

- FOM：`p ≈ 0.93`
- LATIN：`p ≈ 1.31`

因此在当前 41–401 个时间点的区间内：

> FOM wall time 接近线性增长，而当前 LATIN implementation 的 wall time 增长快于线性趋势。

但这个 `p` 只能作为描述性指标，不能视为严格算法复杂度。

不能写成：

> LATIN-PGD 理论复杂度就是 `O(Nt^1.31)`。

---

## 15. 关键发现三：LATIN setup 不是主要瓶颈

| 周期数 | setup / s | LATIN total / s | setup 占 total |
|---:|---:|---:|---:|
| 1 | 0.141990 | 23.726349 | 约 0.60% |
| 2 | 0.250738 | 59.022472 | 约 0.42% |
| 5 | 0.547211 | 205.363929 | 约 0.27% |
| 10 | 1.016734 | 467.361528 | 约 0.22% |

因此：

> **当前效率问题几乎全部位于 `solve_tower_latin_pgd(...)` 内部，而不是 problem setup。**

没有必要优先优化：

- geometry 构建；
- mesh 构建；
- equilibrium operator 构建；
- load-vector 构建；
- elastic initialization。

真正需要拆解的是 LATIN solver 内部成本。

---

## 16. wall-time 结果必须与精度结果联合解释

当前 wall time 显示：

- 1 cycle：LATIN 更快；
- 2 cycles：LATIN 略快；
- 5 cycles：LATIN 更慢；
- 10 cycles：LATIN 明显更慢。

但 5-cycle 和 10-cycle 的 LATIN 结果仍具有良好精度。

### 16.1 5-cycle 已验证精度

Global/FOM 全历史主要误差约为：

- total strain：`0.211%`
- elastic strain：`0.146%`
- stress：`0.146%`
- eps_p：`1.080%`
- alpha：`1.059%`
- r_bar：`1.080%`
- damage：`1.273%`

### 16.2 10-cycle 已验证精度

Global/FOM 全历史主要误差约为：

- total strain：`0.297%`
- elastic strain：`0.138%`
- stress：`0.138%`
- eps_p：`0.999%`
- alpha：`0.927%`
- r_bar：`1.006%`
- damage：`0.945%`

因此当前最准确的综合结论是：

> **LATIN-PGD 在 1–10 周期范围内保持良好精度、稳定收敛和较低 PGD rank，但当前代码实现的 wall-time scaling 没有保持相对于 FOM 的效率优势。**

---

## 17. 为什么当前结果并不否定 LATIN-PGD 方法本身

当前测量的是：

`current implementation wall time`

而不是：

`LATIN-PGD theoretical minimum cost`

当前 tower solver 是 transactional whole-history nonlinear implementation。

每个 outer LATIN iteration 可能包含：

- Local constitutive update；
- search-direction calculation；
- fixed-basis time-function update；
- Trial A；
- saturation decision；
- enrichment；
- enrichment fixed-point iterations；
- Trial B；
- convergence indicator；
- accepted-state reconstruction。

这些步骤中有很多操作会在整个 `Nt × Nq` 历史上重复进行。

因此即使最终 PGD rank 较低，也可能因为：

- 全历史数组被多次扫描；
- trial evaluation 重复；
- fixed-point 内部重复组装；
- reduced temporal system 随 Nt 变大；
- Python 层循环；
- 数组复制和 immutable value reconstruction；

使 wall time 快速增加。

所以目前应提出：

> **当前 tower LATIN-PGD implementation 的实现复杂度需要进一步诊断和优化。**

而不是：

> **LATIN-PGD 方法不具有降阶价值。**

---

## 18. 当前不能推出的结论

### 18.1 不能宣称普适 speedup

不能写：

> LATIN-PGD 相对于 FOM 加速 1.56 倍。

因为 1.56 只出现在当前单周期 benchmark。

### 18.2 不能宣称 rank 决定计算复杂度

`rank << Nt` 不是 speedup 的充分条件。

### 18.3 不能宣称 100 周期一定更慢

当前 1–10 周期 trend 表明继续扩展可能十分昂贵，但没有实际运行 100-cycle LATIN。

### 18.4 不能宣称 crossover 精确位于某一周期

当前只能说已测试点将 crossover 夹在 2 与 5 周期之间。

### 18.5 不能把经验指数写成理论复杂度

`p ≈ 0.93` 与 `p ≈ 1.31` 只是当前区间的经验描述。

---

## 19. 当前可以正式确认的结论

本阶段已经有充分数值依据确认：

1. FOM-3D 的计时口径已经建立；
2. LATIN setup 与 solver 时间已经分开计量；
3. 1、2、5、10 周期均完成三次独立 timing run；
4. 四个规模点的 convergence 与上一阶段完全一致；
5. 1-cycle 与 2-cycle 下，当前 LATIN total wall time 低于 FOM；
6. 5-cycle 与 10-cycle 下，当前 LATIN total wall time 高于 FOM；
7. 当前测试点将 wall-time crossover 夹在 2 与 5 周期之间；
8. FOM wall time 从 1 到 10 周期增长约 8.29 倍；
9. LATIN wall time 从 1 到 10 周期增长约 19.70 倍；
10. PGD rank 同期只由 11 增加到 21；
11. 低 rank 与实际 wall-time speedup 在当前实现中已经明显脱钩；
12. LATIN setup 占 total 的比例始终低于约 0.6%，不是主要瓶颈；
13. 当前主要效率问题位于 `solve_tower_latin_pgd(...)` 内部；
14. 此时直接跳到 100-cycle 的研究价值低于先做 solver 内部成本分解。

---

## 20. FOM-3D 当前阶段判定

可以正式记录：

> **FOM-3D 基础效率 scaling：1、2、5、10-cycle wall-time benchmark PASS。**

这里的 PASS 指：

- timing protocol 成功执行；
- 三次重复性满足诊断需要；
- 所有 LATIN solve 正常收敛；
- scaling trend 已经清晰建立。

但：

> **FOM-3D 整体仍然 OPEN。**

原因是尚未完成：

- LATIN solver 内部成本来源分解；
- memory scaling；
- 100-cycle LATIN 实际 wall time；
- 时间多尺度或高周加速；
- 优化后是否恢复 scaling 优势。

---

## 21. 下一阶段：FOM-3D-5 wall-time 成本来源诊断

下一阶段不应立即运行 100 周期。

优先任务是：

> **对 `solve_tower_latin_pgd(...)` 的内部调用链做 timing-point 审计。**

建议优先拆分：

1. Local stage；
2. search-direction construction；
3. frozen global data preparation；
4. fixed-basis temporal update；
5. Trial A evaluation；
6. saturation decision；
7. PGD enrichment；
8. enrichment fixed-point；
9. Trial B evaluation；
10. Global finishing / candidate construction；
11. convergence indicator；
12. state / basis copy 与 reconstruction。

### 21.1 第一阶段只做源码审计

在没有明确 timing point 前：

> 不应直接在 `latin/` 核心文件中到处加入 `perf_counter()`。

先通过当前源码调用链确认：

- 每个 outer iteration 调用了哪些函数；
- Trial A 和 Trial B 分别调用了什么；
- enrichment fixed-point 内哪些函数被重复执行；
- 哪些函数包含全 `Nt × Nq` 扫描；
- 哪些 immutable copy 可能造成大量数组复制。

### 21.2 第二阶段建立非侵入式 profiling diagnostic

优先考虑：

- 外部 wrapper；
- monkey patch；
- `cProfile`；
- 函数级累计时间；
- 调用次数；
- 单次平均耗时。

原则仍然是：

> **先诊断，再决定是否修改 core implementation。**

### 21.3 后续需要回答的问题

- LATIN 从 1 到 10 周期增加的 solver 时间主要增长在哪里；
- 是 Local stage 主导，还是 enrichment 主导；
- fixed-basis temporal update 是否随 `Nt × rank` 快速增长；
- Trial A / Trial B 是否造成重复成本；
- fixed-point iteration 数是否随时间域增长；
- Python loops、数组复制或线性代数求解是否为主要热点；
- 哪些成本属于 LATIN 理论必需；
- 哪些成本只是当前实现造成的额外开销。

---

## 22. 与原论文理论主线的关系

本阶段没有修改 Bhattacharyya 等人 LATIN-PGD 的理论基准。

仍然保持此前的来源层级。

### 22.1 原论文明确给出的内容

包括：

- LATIN；
- PGD separated representation；
- fixed-basis temporal update；
- enrichment；
- search direction；
- convergence indicator；
- 时间方向 DG0 说明。

### 22.2 当前项目推导得到的内容

包括：

- fiber beam-column 到 canonical material-point coordinate；
- tower equilibrium operator；
- `q ↔ (element, Gauss point, fiber)`；
- 塔筒相容和平衡表达。

### 22.3 当前数值验证得到的内容

包括：

- 1–10 周期 matched asymmetric 精度；
- 1–10 周期 PGD rank growth；
- 1–10 周期 wall-time scaling；
- 2–5 周期之间出现的 wall-time crossover 区间；
- 当前实现中低 rank 与 wall-time speedup 脱钩。

### 22.4 当前工程实现选择

包括：

- `residual_ls`；
- transactional Trial A / Trial B；
- 当前 convergence tolerance；
- 当前 wall-time diagnostic；
- 当前 timing boundary。

这些工程实现选择不能直接等同于原论文理论复杂度。

---

## 23. 当前仍然保留的理论与实现开放问题

### 23.1 原论文 DG0 与当前 Global BE-like 的关系

原论文明确指出时间多变量微分问题使用 discontinuous Galerkin method of order zero。

当前 Global 内部变量更新具有 BE-like 特征。

两者严格等价性仍未完成数学证明。

### 23.2 Local RK4 与 Global BE-like 的差异

当前：

- Local 为逐材料点 RK4 constitutive integration；
- Global 对部分内部变量进行 BE-like reconstruction；
- damage 主要继承 Local。

该差异已经在精度阶段被诊断，但当前没有依据为了效率测试而修改。

### 23.3 residual-LS 的理论来源

`residual_ls` 是当前 tower implementation 的工程扩展。

它不能描述为原论文 Eq. (65)–(71) 的直接代数恒等改写。

### 23.4 高周疲劳仍未真正解决

当前 10 周期只是 multicycle validation。

工程高周疲劳所需的：

- 很大 cycle count；
- time multiscale；
- cycle jump；
- 周期压缩；

仍未建立。

---

## 24. 本阶段使用的 diagnostic 文件

本阶段新增或使用的本地 timing diagnostic 包括：

- `tower_asymmetric_1cycle_efficiency_pilot.py`
- `tower_asymmetric_1cycle_efficiency_pilot_v2.py`
- `tower_asymmetric_1cycle_efficiency_pilot_v3.py`
- `tower_asymmetric_efficiency_scaling_pilot.py`

其中：

- 前两个版本主要用于 API 调试；
- `v3` 是完整单周期 timing pilot；
- `tower_asymmetric_efficiency_scaling_pilot.py` 是统一 1/2/5/10 周期 scaling diagnostic。

这些文件的最终 Git 状态在提交本总结前应通过：

`git status --short`

重新确认。

除非后续决定将统一 scaling diagnostic 正式纳入仓库，否则不应自动把所有临时 pilot 版本一起提交。

---

## 25. 本阶段最重要的研究认识

本阶段真正得到的不是一个简单的“LATIN 快或者慢”的结论，而是把此前两个容易混淆的问题分开。

### 25.1 表示层面

当前 tower response 在 1–10 周期 whole-history 下具有明显低秩特征：

- `Nt` 增长约 9.78 倍；
- rank 只增长约 1.91 倍。

因此：

> **PGD representation 在当前问题中具有较好的低秩表达能力。**

### 25.2 实现层面

当前 transactional LATIN-PGD solver 的 wall time：

- 增长速度快于 FOM；
- 5 周期以后已经失去当前 wall-time 优势；
- 10 周期明显慢于 FOM。

因此：

> **当前实现尚未把低秩表达能力转化为实际计算效率优势。**

### 25.3 后续研究的真正方向

后续工作的重点不应该是盲目增加 cycle count，而应该是：

> **寻找低秩表示已经存在、但 wall-time 优势没有实现的原因。**

---

## 26. 下一次工作开始时应恢复的 checkpoint

下一次继续工作时，应从以下结论开始：

> **FOM-3C 精度验证已经关闭。**

> **FOM-3D 的 1、2、5、10 周期基础 wall-time scaling 已完成。**

当前最关键的综合数据为：

| cycles | Nt | rank | iterations | trials | FOM / s | LATIN / s | FOM/LATIN |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 41 | 11 | 18 | 29 | 37.142 | 23.726 | 1.565 |
| 2 | 81 | 13 | 23 | 36 | 69.305 | 59.022 | 1.174 |
| 5 | 201 | 17 | 33 | 50 | 161.273 | 205.364 | 0.785 |
| 10 | 401 | 21 | 39 | 60 | 307.891 | 467.362 | 0.659 |

下一项正式任务应为：

> **FOM-3D-5：LATIN solver wall-time 成本来源诊断。**

第一步应当是：

> **源码级 timing-point / call-chain 审计，不修改核心算法。**

在完成该诊断以前：

> **暂不直接运行 100-cycle LATIN whole-history。**

---

## 27. 当前阶段最终判定

本阶段可以正式记录：

> **FOM-3D 基础效率与规模扩展验证：1、2、5、10 周期 wall-time scaling PASS。**

但必须同时记录：

> **当前 LATIN-PGD implementation 的 wall-time 优势没有随时间域扩大而保持。**

具体表现为：

- 1-cycle：LATIN 更快；
- 2-cycle：LATIN 仍更快，但优势显著缩小；
- 5-cycle：FOM 已经更快；
- 10-cycle：LATIN 平均耗时约为 FOM 的 1.52 倍。

因此整个 FOM-3D 仍然没有关闭。

下一阶段的核心问题已经从：

> LATIN-PGD 能不能准确计算多循环塔筒响应？

转变为：

> **为什么已经表现出低 PGD rank 的 whole-history LATIN-PGD，在当前实现中 wall time 反而比 FOM 增长得更快？**

这将作为下一阶段代码审计、profiling 和算法效率改进的直接入口。
