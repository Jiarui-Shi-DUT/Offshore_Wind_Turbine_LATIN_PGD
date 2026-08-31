# OPT-3 阶段总结：四阶龙格-库塔局部阶段调用扁平化与外力插值复用优化

**日期：2026-08-31**

**仓库：** `Jiarui-Shi-DUT/Offshore_Wind_Turbine_LATIN_PGD`  
**分支：** `perf/tower-local-stage-optimization`  
**OPT-3 代码提交：** `2a88b68` — `perf: flatten RK4 local-stage hot path`

---

## 1. 阶段定位

本阶段是在 OPT-1 和 OPT-2 完成之后，对塔筒 LATIN-PGD 10 周期算例重新进行函数级耗时分析，并据此确定的第三轮实现级性能优化。

本阶段不改变 LATIN-PGD 的数学算法，不改变 PGD 富集策略，不改变局部本构演化方程，也不改变四阶龙格-库塔时间积分公式。其目标仅为削减局部本构阶段中高频 Python 函数调用和重复外力插值产生的实现开销。

本阶段正式定义为：

> **OPT-3：四阶龙格-库塔阶段调用扁平化与外力插值复用优化及多周期等价性验证**

最终判定：

> **PASS**

---

## 2. 为什么在 OPT-3 前重新做性能定位

OPT-1 和 OPT-2 已经显著改变了原始代码的耗时结构，因此不能继续直接使用优化前的热点比例判断第三轮优化对象。

优化前 10 周期 LATIN-PGD 正常运行总时间约为：

$$
T_{\mathrm{LATIN,original}} \approx 467.362\ \mathrm{s}
$$

OPT-1 后约为：

$$
T_{\mathrm{LATIN,OPT1}} = 223.160690\ \mathrm{s}
$$

OPT-2 后约为：

$$
T_{\mathrm{LATIN,OPT2}} = 174.717540\ \mathrm{s}
$$

因此在进入 OPT-3 前，重新使用 `cProfile` 对当前 OPT-2 版本进行 10 周期函数级耗时分析。

运行命令：

```powershell
python -m cProfile -o tower_10cycle_profile_opt2.pstats .\tower_asymmetric_efficiency_scaling_pilot.py --cycles 10
```

性能分析运行本身引入了显著额外开销，因此其绝对时间不作为正式效率结果，只用于比较函数之间的相对耗时层级和调用次数。

---

## 3. OPT-2 后重新性能分析结果

性能分析运行中：

- FOM：`452.165537 s`
- LATIN setup：`1.483000 s`
- LATIN solver：`224.160156 s`
- LATIN total：`225.643156 s`

数值路径保持：

- `iterations = 39`
- `attempted = 39`
- `trial evaluations = 60`
- `PGD rank = 21`
- `modes added = 21`
- `final xi = 8.941607234831e-06`

这说明性能分析并未改变求解路径。

### 3.1 当前 LATIN 顶层热点

`solve_tower_local_stage`：

- 调用次数：39
- 累计时间：`199.196 s`

LATIN 求解器总性能分析时间：

$$
224.155\ \mathrm{s}
$$

因此当前局部本构阶段约占：

$$
\frac{199.196}{224.155}\times100\%
\approx 88.9\%
$$

虽然相比优化前约 95.6% 已明显下降，但局部本构阶段仍是第一瓶颈。

### 3.2 当前局部阶段内部热点

核心函数统计：

| 函数 | 调用次数 | 自身时间 / s | 累计时间 / s |
|---|---:|---:|---:|
| `solve_tower_local_stage` | 39 | 22.546 | 199.196 |
| `_integrate_one_local_step` | 4,992,000 | 39.261 | 158.149 |
| `rate` | 19,968,000 | 6.518 | 100.664 |
| `_local_state_rate` | 19,968,000 | 18.712 | 94.145 |
| `local_rates_from_forces` | 24,972,480 | 50.852 | 75.900 |
| `_positive_part` | 49,944,960 | 7.012 | 12.874 |
| `_interpolate` | 79,872,000 | 8.571 | 8.571 |
| `isotropic_force_from_transformed_force` | 29,976,960 | 7.965 | 7.965 |
| `_clip_damage_scalar` | 34,968,960 | 4.216 | 4.216 |
| `unilateral_elastic_strain` | 5,004,480 | 1.804 | 2.410 |

其中 `_integrate_one_local_step` 的累计时间约占整个 LATIN solver 性能分析时间：

$$
\frac{158.149}{224.155}\times100\%
\approx 70.6\%
$$

说明在 OPT-1、OPT-2 完成之后，**四阶龙格-库塔材料状态时间推进仍然是最主要的剩余热点。**

---

## 4. OPT-3 的优化对象

原有热路径为：

```text
_integrate_one_local_step
    ↓
rate
    ↓
_local_state_rate
    ↓
_interpolate × 4
    ↓
local_rates_from_forces
```

10 周期下：

- `_integrate_one_local_step`：4,992,000 次
- `rate`：19,968,000 次
- `_local_state_rate`：19,968,000 次
- `_interpolate`：79,872,000 次

对于每一个材料点时间步，四阶龙格-库塔具有四个阶段：

$$
k_1=f(z_n,t_n)
$$

$$
k_2=f\left(z_n+\frac{\Delta t}{2}k_1,t_n+\frac{\Delta t}{2}\right)
$$

$$
k_3=f\left(z_n+\frac{\Delta t}{2}k_2,t_n+\frac{\Delta t}{2}\right)
$$

$$
k_4=f(z_n+\Delta t k_3,t_{n+1})
$$

最终：

$$
z_{n+1}
=
z_n+
\frac{\Delta t}{6}
\left(k_1+2k_2+2k_3+k_4\right)
$$

在给定外力历史线性插值的条件下，四个阶段实际上只对应三个不同的外力时刻：

$$
0,\qquad \frac12,\qquad 1
$$

其中 $k_2$ 和 $k_3$ 共用同一个中点外力状态。

原实现中，每个 RK4 阶段都通过 `_local_state_rate()` 重新调用 4 次 `_interpolate()`。因此每个材料时间步共有：

$$
4\times4=16
$$

次插值调用。

10 周期总插值调用数：

$$
4,992,000\times16
=
79,872,000
$$

与性能分析计数完全一致。

---

## 5. OPT-3 实际代码修改

本阶段只修改：

```text
latin/local_stage.py
```

核心修改位于：

```python
_integrate_one_local_step()
```

### 5.1 调用路径扁平化

原实现通过内部函数 `rate()` 再调用 `_local_state_rate()`。

OPT-3 后，RK4 四个阶段直接调用：

```python
local_rates_from_forces(...)
```

从热路径中绕过：

- `rate()`
- `_local_state_rate()`

但这两个函数暂时保留在文件中，没有在本轮删除，以避免将性能优化与代码清理混为一次修改。

### 5.2 外力插值复用

每一个材料时间步预先计算：

- 起点外力状态；
- 中点外力状态；
- 终点外力状态。

其中中点由 $k_2$ 和 $k_3$ 共用。

为了尽可能保持原有浮点计算顺序，起点和终点仍保持原插值表达式形式，而没有直接简单替换为 `old` 和 `new`。

### 5.3 明确没有修改的内容

本阶段没有修改：

- LATIN 外层迭代；
- PGD 富集算法；
- 空间最小二乘求解；
- 时间函数更新；
- 局部本构方程；
- 材料参数；
- RK4 四阶段定义；
- RK4 权重；
- 损伤演化方程；
- 塑性、运动硬化和各向同性硬化演化方程。

因此本阶段属于严格的实现级优化。

---

## 6. 回归测试

### 6.1 局部阶段单元测试

运行：

```powershell
python -m pytest .\tests\test_latin_local_stage.py -q
```

结果：

```text
2 passed, 10 warnings in 5.39s
```

### 6.2 全仓库回归测试

运行：

```powershell
python -m pytest -q
```

结果：

```text
313 passed, 10 warnings in 102.29s
```

10 条 warning 均为 Matplotlib / setuptools 的 `distutils` 弃用提示，与 OPT-3 无关。

---

## 7. 1 周期数值与性能验证

运行：

```powershell
python .\tower_asymmetric_efficiency_scaling_pilot.py --cycles 1
```

结果：

- `Nt = 41`
- FOM wall time：`35.086349 s`
- FOM max Newton iterations：4
- FOM max top displacement：`1.173015956970e+00 m`
- LATIN setup：`0.129788 s`
- initial equilibrium residual：`1.715324872498e-06 N`
- iterations：18
- attempted：18
- trial evaluations：29
- PGD rank：11
- modes added：11
- final xi：`7.918424536257e-06`
- LATIN solver：`8.661120 s`
- LATIN total：`8.790908 s`
- FOM/LATIN：`3.991209`

与 OPT-2 相比：

$$
9.433297\ \mathrm{s}
\rightarrow
8.790908\ \mathrm{s}
$$

节省：

$$
0.642389\ \mathrm{s}
$$

相对降低：

$$
\frac{9.433297-8.790908}{9.433297}\times100\%
\approx 6.81\%
$$

增量加速约：

$$
1.073\times
$$

数值路径与 OPT-2 完全一致。

---

## 8. 10 周期数值与性能验证

运行：

```powershell
python .\tower_asymmetric_efficiency_scaling_pilot.py --cycles 10
```

结果：

- `Nt = 401`
- FOM wall time：`295.177001 s`
- FOM max Newton iterations：4
- FOM max top displacement：`1.426595421862e+00 m`
- LATIN setup：`1.012338 s`
- initial equilibrium residual：`2.073353735670e-06 N`
- iterations：39
- attempted：39
- trial evaluations：60
- PGD rank：21
- modes added：21
- final xi：`8.941607234831e-06`
- LATIN solver：`161.851950 s`
- LATIN total：`162.864288 s`
- FOM/LATIN：`1.812411`

数值路径与 OPT-2 完全一致。

### 8.1 相对 OPT-2 的增量收益

OPT-2：

$$
174.717540\ \mathrm{s}
$$

OPT-3：

$$
162.864288\ \mathrm{s}
$$

节省：

$$
174.717540-162.864288
=
11.853252\ \mathrm{s}
$$

相对降低：

$$
\frac{11.853252}{174.717540}\times100\%
\approx 6.78\%
$$

增量加速约：

$$
1.073\times
$$

### 8.2 机器状态参考

OPT-2 的 10 周期 FOM 时间：

$$
294.307999\ \mathrm{s}
$$

OPT-3 的 10 周期 FOM 时间：

$$
295.177001\ \mathrm{s}
$$

两者仅相差约 0.3%，说明两次基准运行的机器状态非常接近。OPT-3 中 LATIN 约 6.8% 的下降明显超过该波动量级，因此可以认为本轮优化具有真实且可重复验证的性能收益。

---

## 9. 三轮实现级优化的累计结果

| 阶段 | 10 周期 LATIN 总时间 / s | 相对上一阶段 |
|---|---:|---:|
| 原始实现 | ≈467.362 | — |
| OPT-1 | 223.160690 | ↓约52.3% |
| OPT-2 | 174.717540 | ↓约21.7% |
| OPT-3 | **162.864288** | **↓约6.8%** |

从原始实现到 OPT-3：

$$
467.362
\rightarrow
162.864288\ \mathrm{s}
$$

累计时间降低约：

$$
\frac{467.362-162.864288}{467.362}\times100\%
\approx 65.2\%
$$

累计加速约：

$$
\frac{467.362}{162.864288}
\approx 2.87\times
$$

---

## 10. 当前 10 周期 FOM / LATIN 对比

当前 OPT-3 基准：

$$
T_{\mathrm{FOM}}=295.177001\ \mathrm{s}
$$

$$
T_{\mathrm{LATIN}}=162.864288\ \mathrm{s}
$$

因此：

$$
\frac{T_{\mathrm{FOM}}}{T_{\mathrm{LATIN}}}
=
1.812411
$$

当前 LATIN-PGD 实际计算时间约为 FOM 的：

$$
\frac{162.864288}{295.177001}\times100\%
\approx55.2\%
$$

需要强调：这一结论只针对当前塔筒模型、当前离散规模、当前 Python 实现和当前收敛参数。

当前 FOM 自身在函数级性能分析中仍存在大量标量 NumPy `clip` 等实现开销，因此该结果不能直接解释为 LATIN-PGD 与 FOM 两种数学方法在最优实现条件下的理论效率比。

---

## 11. 科研认识

### 11.1 OPT-3 再次证明原始效率交叉不能被解释为算法固有边界

原始 Python 实现曾在 2 到 5 周期之间出现 FOM/LATIN 计算效率交叉。

经过 OPT-1、OPT-2 和 OPT-3 三轮均不改变数学算法的实现级优化后，10 周期 LATIN 总时间已经从约 467 s 降至约 163 s，并显著快于当前 FOM。

因此此前观察到的效率交叉主要受原始实现开销影响，不能直接视为 LATIN-PGD 方法本身的固有效率边界。

### 11.2 局部材料历史积分仍然是当前主要成本来源

OPT-2 后重新性能分析显示，局部本构阶段仍占 LATIN solver 约 88.9%，其中 RK4 材料状态时间推进约占整个 LATIN solver 性能分析时间 70.6%。

因此 Local 仍是当前主要优化对象。

### 11.3 实现级微优化的边际收益正在递减

三轮优化的相对收益呈现：

$$
\mathrm{OPT1}\; \gt \;\mathrm{OPT2}\; \gt \;\mathrm{OPT3}
$$

OPT-3 相对 OPT-2 的收益约为 6.8%，明显低于前两轮。

这说明当前代码已经逐步离开最明显的标量 NumPy 调度和固定量重复计算瓶颈。若希望获得下一阶段更显著的性能提升，仅继续消除个别 Python 小函数调用的潜在收益可能逐渐有限。

---

## 12. 当前已经确认的事实

本阶段可以确认：

1. OPT-2 后局部本构阶段仍是当前 LATIN-PGD 第一性能瓶颈。
2. 10 周期下局部本构阶段在带性能分析器运行中约占 LATIN solver 的 88.9%。
3. `_integrate_one_local_step` 仍是局部阶段最大热点。
4. 原实现存在约 1997 万次 `rate()` 调用、1997 万次 `_local_state_rate()` 调用和约 7987 万次 `_interpolate()` 调用。
5. OPT-3 在 RK4 热路径中绕过了 `rate()` 和 `_local_state_rate()`，并复用了 RK4 中点外力状态。
6. RK4 数学公式和本构演化方程未改变。
7. 局部阶段测试为 `2 passed`。
8. 全仓库测试为 `313 passed`。
9. 1 周期数值路径完全保持，LATIN 总时间由 OPT-2 的约 9.433 s 降至约 8.791 s。
10. 10 周期数值路径完全保持，LATIN 总时间由 OPT-2 的约 174.718 s 降至约 162.864 s。
11. OPT-3 相对 OPT-2 的 10 周期时间降低约 6.8%。
12. 从原始实现到 OPT-3，10 周期 LATIN 总时间累计降低约 65.2%，约获得 2.87 倍加速。

---

## 13. 当前仍不能声称的结论

本阶段仍不能证明：

1. 当前局部阶段已经达到最优实现；
2. 当前 LATIN-PGD 实现已经达到理论最优效率；
3. LATIN-PGD 对任意长周期问题都必然快于 FOM；
4. 当前 FOM/LATIN 比值能够代表两种数学方法本身的理论效率；
5. OPT-4 的材料点向量化一定会获得显著加速；
6. 进一步性能提升一定需要修改算法；
7. 当前 Python 解释器层面的所有局部本构开销都已被消除。

---

## 14. 下一阶段建议

原计划中的 OPT-4 是：

> **材料点方向的批量化 / 向量化。**

理论基础是：对于同一时间位置，不同材料点的局部本构更新在空间上相互独立，而同一个材料点沿时间方向仍必须保持历史顺序。

因此可能的组织形式是：

$$
\boxed{\text{时间方向顺序推进}\; + \;\text{材料点方向批量计算}}
$$

当前塔筒共有：

$$
N_q=320
$$

个材料点。

如果进入 OPT-4，应保持：

- 单材料点时间历史顺序不变；
- RK4 数学公式不变；
- 本构演化关系不变；

仅改变多个材料点在同一时间步上的计算组织方式。

不过在正式进入 OPT-4 前，应先完成 OPT-3 的代码提交、阶段文档提交和远端同步，形成稳定 checkpoint。

---

## 15. Git 状态

OPT-3 代码已经独立提交：

```text
2a88b68  perf: flatten RK4 local-stage hot path
```

当前分支：

```text
perf/tower-local-stage-optimization
```

代码提交后本地分支相对远端领先 1 个 commit。

当前仍存在未跟踪的性能分析和诊断文件，包括：

```text
tower_10cycle_profile.pstats
tower_10cycle_profile_opt2.pstats
tower_asymmetric_10cycle_latin_diagnostic.py
tower_asymmetric_1cycle_efficiency_pilot.py
tower_asymmetric_1cycle_efficiency_pilot_v2.py
tower_asymmetric_1cycle_efficiency_pilot_v3.py
tower_asymmetric_1cycle_latin_diagnostic.py
tower_asymmetric_1cycle_latin_tight_diagnostic.py
tower_asymmetric_2cycle_latin_diagnostic.py
tower_asymmetric_2cycle_latin_diagnostic_fixed.py
tower_asymmetric_5cycle_latin_diagnostic.py
tower_asymmetric_efficiency_scaling_pilot.py
tower_fresh_local_fom_global_diagnostic.py
tower_local_fom_global_diagnostic.py
```

这些文件不属于 OPT-3 正式提交内容，后续 Git 操作仍应避免使用：

```powershell
git add .
```

---

## 16. 阶段结论

本阶段重新测量了 OPT-2 后的性能热点，确认局部本构阶段仍占 LATIN-PGD 主要计算成本，并据此实施了严格受控的四阶龙格-库塔热路径优化。

OPT-3 未改变 LATIN、PGD、本构或 RK4 数学结构，仅通过调用路径扁平化和中点外力插值复用削减高频 Python 调度开销。

1 周期和 10 周期均保持完全相同的迭代次数、试算次数、PGD 秩、模式数和最终收敛指标；全仓库 313 项测试全部通过。

10 周期 LATIN 总时间由：

$$
174.717540\ \mathrm{s}
$$

进一步下降至：

$$
162.864288\ \mathrm{s}
$$

相对 OPT-2 降低约 6.8%。

从原始实现约 467.362 s 到当前约 162.864 s，累计时间降低约 65.2%，累计加速约 2.87 倍。

因此：

> **OPT-3：四阶龙格-库塔阶段调用扁平化与外力插值复用优化及多周期等价性验证 — PASS**
