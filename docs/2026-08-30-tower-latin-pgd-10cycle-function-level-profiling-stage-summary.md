# 海上风机塔筒 LATIN-PGD 10 周期函数级耗时剖析阶段总结

**日期：** 2026-08-30  
**项目：** Offshore_Wind_Turbine_LATIN_PGD  
**分支：** `feature/offshore-wind-turbine-tower-fatigue`  
**建议阶段编号：** FOM-3D-6  
**当前状态：** 10 周期 `cProfile` 函数级耗时剖析已完成；当前实现的主要耗时链已定位；低侵入式真实 wall-time 分项确认尚未开始  
**核心原则：** 本阶段只做诊断，不修改任何 `latin/` 核心算法

---

## 1. 本阶段的起点与研究问题

上一阶段已经完成 1、2、5、10 周期 matched asymmetric benchmark 的精度与总 wall-time scaling，并确认：

- 10 周期时间点数为 `Nt = 401`；
- 总材料点数为 `Nq = 320`；
- LATIN-PGD 稳定收敛；
- 外层 LATIN iteration 为 `39`；
- trial evaluation 为 `60`；
- 最终 PGD rank 为 `21`；
- 10 周期未加 profiler 的正式 LATIN total wall time 均值约为 `467.36 s`；
- 10 周期 FOM wall time 均值约为 `307.89 s`；
- 因而当前实现下 10 周期 LATIN-PGD 已经慢于 FOM；
- 源码审计已经确认：39 次 outer LATIN iteration 中，每次必有一次 Trial A fixed-basis temporal update；21 次 successful enrichment 又各带来一次 enlarged-basis full temporal re-optimization，因此共有 `39 + 21 = 60` 次 full-basis temporal update；
- 但在上一阶段结束时，仍不知道总 wall time 主要花在 Local stage、temporal update、PGD enrichment、residual-LS / LSMR，还是其他模块。

因此，本阶段的核心问题被严格定义为：

> **在完全保持 10 周期 benchmark、算法参数和求解路径不变的条件下，对当前 LATIN-PGD Python 实现进行函数级 profiling，定量回答各模块到底消耗多少时间，并识别真正的第一耗时瓶颈。**

本阶段尤其希望验证前期的几个主要怀疑：

1. 60 次 full-basis temporal update 是否是 10 周期效率恶化的主要来源；
2. PGD enrichment 中的 fixed-point / residual-LS / LSMR 是否是主要来源；
3. Local nonlinear stage 是否存在显著的 whole-history 重复积分成本；
4. 当前耗时问题究竟更偏向算法层重复工作，还是实现层 Python / NumPy 标量调用开销；
5. 哪个方向才值得进入下一阶段的低侵入式 wall-time 验证和后续优化。

---

## 2. Profiling benchmark 与原 10 周期正式 benchmark 保持一致

本阶段没有重新定义新的物理问题，而是直接复用已经完成精度和效率验证的：

`tower_asymmetric_efficiency_scaling_pilot.py --cycles 10`

10 周期 benchmark 保持：

- `cycles = 10`
- `Nt = 401`
- `Nq = 320`
- `Fmax = +1.0 MN`
- `Fmin = -0.5 MN`
- `R = -0.5`
- 每周期 40 个增量
- `spatial_strategy = "residual_ls"`
- `tolerance = 1.0e-5`
- 其余 LATIN、PGD、enrichment 参数保持原正式 benchmark 不变

Profiling 命令为：

```powershell
python -m cProfile -o tower_10cycle_profile.pstats .\tower_asymmetric_efficiency_scaling_pilot.py --cycles 10
```

生成：

`tower_10cycle_profile.pstats`

随后使用 `pstats` 逐级查询：

- 顶层累计耗时；
- `solve_tower_latin_pgd()` 的直接 callees；
- `solve_tower_local_stage()` 的直接 callees；
- `_integrate_one_local_step()` 的直接 callees；
- `_local_state_rate()` 的直接 callees；
- `local_rates_from_forces()` 的直接 callees；
- `enrich_tower_pgd_basis_once()` 的直接 callees；
- `_raw_fixed_point()` 的直接 callees。

因此，本阶段不是根据代码结构“估计”时间，而是基于实际 10 周期运行的函数调用计数和累计时间进行诊断。

---

## 3. Profiling 没有改变 LATIN 的数值求解路径

`cProfile` 运行得到：

- `termination_reason = converged`
- `converged = True`
- `failure_reason = None`
- `iterations = 39`
- `attempted = 39`
- `trial evaluations = 60`
- `PGD rank = 21`
- `modes added = 21`
- `final xi = 8.941607234831e-06`

这与未加 profiler 的正式 10 周期 benchmark 的求解路径一致。

因此可以确认：

> **profiling 没有改变 LATIN-PGD 的迭代次数、trial 数、最终 rank 或最终收敛指标。**

这使得本次 `.pstats` 可以用于回答“当前这条求解路径的计算时间花在哪里”。

---

## 4. 必须严格区分：profiling 时间不能直接替代正式 wall time

`cProfile` 运行得到：

- FOM wall time：`473.474468 s`
- LATIN setup：`1.471594 s`
- LATIN solver：`699.096638 s`
- LATIN total：`700.568232 s`

而此前未加 profiler 的正式 10 周期均值约为：

- FOM：`307.891095 s`
- LATIN total：`467.361528 s`

因此 profiler 明显引入了额外开销。

粗略比较：

$$ \frac{700.568232}{467.361528}\approx1.50 $$

$$ \frac{473.474468}{307.891095}\approx1.54 $$

也就是说，本次函数级 profiling 使 FOM 与 LATIN 的 wall time 都大幅增加。

同时，`pstats` 记录：

- 总函数调用数约 `2.615 × 10^9`
- primitive calls 约 `2.477 × 10^9`

这说明对于当前这种包含数千万级 Python / NumPy 小函数调用的程序，`cProfile` 的跟踪开销不可忽略。

因此，本阶段所有百分比都必须按照下面的证据边界解释：

> **函数级 profiling 可用于定位当前程序的主要耗时链和相对量级，但不能把 `699 s` 或其中的分项秒数直接当作未经 profiler 的最终生产 wall time。**

尤其对于大量标量函数调用的 Local stage，`cProfile` 可能放大其绝对成本。因此，本阶段之后仍需要一次低侵入式 `perf_counter` 分项计时来确认真实 wall-time 占比。

---

## 5. 第一层结果：10 周期 LATIN-PGD 的绝对主耗时是 Local stage

`solve_tower_latin_pgd()` 的直接 callees 为：

| 顶层模块 | 调用次数 | `cumtime / s` | 占 profiler 下 LATIN solver |
|---|---:|---:|---:|
| `solve_tower_local_stage` | 39 | `668.605` | `95.64%` |
| `enrich_tower_pgd_basis_once` | 21 | `19.037` | `2.72%` |
| `compute_tower_descent_search_directions` | 39 | `6.122` | `0.88%` |
| Trial-A `update_tower_pgd_time_functions` | 39 | `2.423` | `0.35%` |
| `evaluate_tower_trial` | 60 | `1.091` | `0.16%` |
| `build_unrelaxed_candidate` | 60 | `1.071` | `0.15%` |
| `prepare_frozen_global_data` | 39 | `0.194` | `0.03%` |
| solver 自身控制、复制、封装等剩余开销 | — | 约 `0.35` | 约 `0.05%` |

其中：

$$ \frac{668.605}{699.091}\approx95.64\% $$

也就是说，在本次 profiled run 中：

> **Local nonlinear stage 一项就占了 LATIN solver 约 95.6% 的累计时间。**

除 Local stage 外，其他全部 LATIN-PGD 顶层模块合计只有：

$$ 699.091-668.605\approx30.486\ \mathrm{s} $$

约占：

$$ 4.36\% $$

这使本阶段的主要认识发生了明显变化：

> **当前 10 周期实现的第一耗时矛盾，不是 PGD rank、LSMR 或 full-basis temporal update，而是每个 outer LATIN iteration 都重新执行的完整 Local constitutive history integration。**

---

## 6. Local stage 的调用次数与当前算法结构严格闭合

### 6.1 每次 Local stage 都完整遍历时间历史和全部材料点

当前塔筒 Local stage 的实现逻辑为：

- 不同材料点 `q` 在空间上相互独立；
- 每一个材料点内部的历史变量需要沿时间顺序推进；
- 内部变量 `[eps_p, alpha, r_bar, D]` 使用与已验证一维算例相同的经典 RK4 material-point integrator；
- 对于每次 outer LATIN iteration，Local stage 都从该 iteration 的 global baseline 出发，对完整时间历史重新做一次 local projection。

10 周期中：

- `Nt = 401`
- 真正的时间增量数为 `400`
- `Nq = 320`
- outer LATIN iteration = `39`

因此理论上应该执行：

$$ 39\times400\times320=4,992,000 $$

次 material-point time-step integration。

Profiler 实测：

`_integrate_one_local_step` 调用次数正好为：

`4,992,000`

这与算法结构严格一致。

因此可以把这一点提升为确定事实：

> **当前 10 周期 LATIN 实际执行了 4,992,000 次材料点 RK4 时间步积分。**

---

## 7. RK4 四阶段进一步将本构状态速率评价放大到近两千万次

每一个 `_integrate_one_local_step()` 使用经典 RK4：

- `k1`
- `k2`
- `k3`
- `k4`

因此每个材料时间步需要 4 次 `_local_state_rate()`。

理论调用数：

$$ 4\times4,992,000=19,968,000 $$

Profiler 实测：

- `rate`：`19,968,000`
- `_local_state_rate`：`19,968,000`

完全一致。

因此，10 周期 LATIN 中仅 RK4 中间阶段就执行了：

> **19,968,000 次本构状态速率评价。**

---

## 8. 为什么 `local_rates_from_forces` 达到 24,972,480 次

除了 RK4 四阶段内部的：

`19,968,000`

次 `local_rates_from_forces()`，当前 `tower_local_stage.py` 还会在每个正式时间点重新计算并存储：

- `plastic_strain_rate`
- `alpha_rate`
- `r_bar_rate`
- `damage_rate`

每次 Local stage 的正式时间点数为：

$$ 401\times320=128,320 $$

39 次 Local stage：

$$ 39\times401\times320=5,004,480 $$

因此总调用数：

$$ 19,968,000+5,004,480=24,972,480 $$

Profiler 实测：

`local_rates_from_forces = 24,972,480`

完全闭合。

因此：

> **当前 10 周期 LATIN 共执行了 24,972,480 次材料演化速率公式评价。**

---

## 9. 其他 Local 调用计数同样与公式完全一致

### 9.1 单边弹性应变

`unilateral_elastic_strain` 每个正式 `(t,q)` 点计算一次，因此：

$$ 39\times401\times320=5,004,480 $$

Profiler：

`5,004,480`

### 9.2 RK4 内部固定力历史插值

每一次 `_local_state_rate` 需要插值：

- stress
- beta
- R_bar
- energy release rate

共 4 次 `_interpolate()`。

因此：

$$ 19,968,000\times4=79,872,000 $$

Profiler：

`79,872,000`

再次完全一致。

这一系列严格对应说明：

> 本次 profiler 中最主要的巨量函数调用并不是异常重复、死循环或错误 retry，而是当前 Local stage 的设计与实现自然产生的计算量。

---

## 10. Local stage 第一层耗时分解

`solve_tower_local_stage()` 的 profiled cumulative time：

`668.605 s`

直接 callees 为：

| Local 子模块 | 调用次数 | `cumtime / s` | 占 Local stage |
|---|---:|---:|---:|
| `_integrate_one_local_step` | 4,992,000 | `499.921` | `74.77%` |
| 时间点最终 `local_rates_from_forces` | 5,004,480 | `86.640` | `12.96%` |
| `unilateral_elastic_strain` | 5,004,480 | `55.085` | `8.24%` |
| Local 主循环、索引、赋值及其他自身开销 | — | 约 `26.7` | 约 `4.0%` |
| 材料序列、field copy 等 | — | 很小 | 很小 |

前三项合计：

$$ 499.921+86.640+55.085=641.646\ \mathrm{s} $$

约占 Local stage：

$$ \frac{641.646}{668.605}\approx96\% $$

因此可以确认：

> **Local stage 的绝大多数时间确实花在材料本构历史更新，而不是 state object 构造、数组复制或外围 bookkeeping。**

---

## 11. RK4 单步积分内部的耗时结构

`_integrate_one_local_step()` 总 cumulative time：

`499.921 s`

其直接 callees：

| RK4 单步内部 | 调用次数 | `cumtime / s` | 占该模块 |
|---|---:|---:|---:|
| 4 个 RK4 阶段的 `rate()` | 19,968,000 | `383.946` | `76.8%` |
| 每步结束后的 `np.clip` 损伤限幅 | 4,992,000 | `49.318` | `9.9%` |
| `np.all(...)` 有限性检查 | 4,992,000 | `18.305` | `3.7%` |
| `_integrate_one_local_step` 自身数组运算等 | — | `45.402` | `9.1%` |
| `np.asarray` | 9,984,000 | `1.226` | `0.25%` |
| ndarray `copy` | 4,992,000 | `1.724` | `0.34%` |

因此：

> **RK4 本身的第一大成本是四个阶段反复调用本构状态速率，而不是 RK4 加权组合公式本身。**

`rate()` 的 `383.946 s` 单独已经约占 profiler 下整个 LATIN solver：

$$ \frac{383.946}{699.091}\approx54.9\% $$

---

## 12. `_local_state_rate`：真正的大头不是历史插值，而是材料演化公式

`_local_state_rate()`：

- 调用次数：`19,968,000`
- cumulative time：`376.287 s`

其内部：

| 子项 | 调用次数 | `cumtime / s` |
|---|---:|---:|
| `local_rates_from_forces` | 19,968,000 | `337.178` |
| `_interpolate` | 79,872,000 | `9.366` |
| `np.asarray` | 19,968,000 | `7.539` |

因此：

$$ \frac{337.178}{376.287}\approx89.6\% $$

这说明：

> **四个固定力历史的线性插值并不是 Local stage 的主要瓶颈；真正昂贵的是每个 RK4 stage 中对材料演化方程本身的重复评价。**

---

## 13. `local_rates_from_forces` 是目前最重要的实现级热点

所有调用路径合计：

- `local_rates_from_forces` 调用次数：`24,972,480`
- cumulative time：`423.818 s`

其中来自：

- RK4 stage：`19,968,000` 次，`337.178 s`
- 正式时间点最终速率重算：`5,004,480` 次，`86.640 s`

严格满足：

$$ 337.178+86.640=423.818\ \mathrm{s} $$

其内部主要耗时：

| `local_rates_from_forces` 内部 | 调用次数 | `cumtime / s` |
|---|---:|---:|
| 标量 `np.clip(damage, ...)` | 24,972,480 | `258.060` |
| 函数自身公式与 Python 运算 | — | `96.776` |
| `isotropic_force_from_transformed_force` | 24,972,480 | `25.351` |
| `_positive_part` | 49,944,960 | `17.294` |
| `np.finfo(float).eps` / `getlimits` | 24,972,480 | `10.941` |
| `k_viscoplastic` property | 24,972,480 | `5.789` |
| `Y0` property | 24,972,480 | `6.210` |
| built-in `abs` | 49,944,960 | `3.398` |

最显著的是：

> **24,972,480 次对单个 damage 标量调用 `np.clip`，在 profiled run 中累计达到 258.060 s。**

当前源码等价于反复执行：

```python
damage_safe = float(
    np.clip(
        damage,
        0.0,
        material.damage_upper_bound,
    )
)
```

这里每次 `np.clip` 处理的不是大数组，而是单个 scalar。

因此本阶段识别出一个明确的实现级现象：

> **当前 Local constitutive kernel 中存在大量标量级 NumPy 调用。NumPy 的通用 dispatch / wrapping 成本在单次调用中很小，但在 10^7 量级重复后成为显著热点。**

类似地：

```python
np.finfo(float).eps
```

在每次本构速率评价中都重新访问一次，累计达到近 2500 万次。

必须强调：

> 这并不说明 `np.clip` 或 `np.finfo` 在一般数组计算中低效，而是说明“在 Python 热循环内对 scalar 重复调用通用 NumPy API”是当前实现不利的使用模式。

---

## 14. Local stage 的问题具有“算法层 + 实现层”两层叠加性质

本阶段不能把 Local stage 的高成本简单归结为“Python 慢”，也不能只归结为“RK4 贵”。

更准确的结构是：

### 14.1 算法/实现路径层面的重复历史成本

每一个 outer LATIN iteration 都重新执行完整 local projection：

$$ 39\times320\times400=4,992,000 $$

次材料时间步积分。

也就是说，PGD 的低秩表示主要作用在 Global correction / separated representation 一侧，但当前 Local nonlinear constitutive stage 仍然在完整 `(t,q)` 网格上反复计算。

### 14.2 单次材料积分的 Python 实现成本

每一个材料时间步内部又采用：

- RK4 四阶段；
- 多次 Python 函数调用；
- scalar `np.clip`；
- scalar property access；
- scalar `np.finfo`；
- 小数组 `np.asarray` / `np.all`。

因此当前高成本是：

$$ \text{完整历史反复积分}\times\text{高频 scalar Python/NumPy constitutive kernel} $$

共同造成的。

这一区分非常重要，因为后续优化可以有不同层级：

- 只优化 scalar kernel：不改变算法；
- 向材料点维度向量化/编译：改变实现方式，不改变本构方程；
- 减少 Local stage 的完整历史重复次数：属于更高层算法结构优化，需要更严格的理论与数值论证。

本阶段只完成诊断，不进行上述任何修改。

---

## 15. PGD enrichment 并不是当前第一耗时瓶颈

`enrich_tower_pgd_basis_once()`：

- 调用次数：`21`
- cumulative time：`19.037 s`
- 约占 profiler 下 LATIN solver：`2.72%`

其内部：

| enrichment 子项 | 调用次数 | `cumtime / s` |
|---|---:|---:|
| `_raw_fixed_point` | 21 | `18.009` |
| enlarged-basis full `update_tower_pgd_time_functions` | 21 | `0.925` |
| `_post_fixed_point_transform` | 21 | `0.066` |
| input validation | 21 | `0.028` |
| residual norm 等 | — | 很小 |

因此约：

$$ \frac{18.009}{19.037}\approx94.6\% $$

的 enrichment 时间都在新 mode 的 raw fixed-point。

---

## 16. 21 个 PGD 模态实际上经历了 523 次空间—时间 fixed-point 更新

`_raw_fixed_point()` 内部 profiler 记录：

- `solve_tower_residual_ls_spatial`：`523` 次
- `_temporal_solve`：`544` 次
- `_pair_change`：`523` 次
- `_inloop_orthogonalize_spatial`：`544` 次
- `apply_spatial`：`544` 次
- seed：`21` 次

由于每个 fixed-point iteration 都包含一次 residual-LS spatial solve，因此：

> **当前 21 次 enrichment 合计执行了 523 次空间—时间 fixed-point iteration。**

平均每个新 PGD mode：

$$ \frac{523}{21}\approx24.9 $$

次 fixed-point。

而 temporal solve 为：

$$ 544=523+21 $$

其含义与代码结构一致：

- 每个 enrichment 根据初始 spatial seed 先求一次 temporal function：21 次；
- 每个 fixed-point spatial update 后再求一次 temporal function：523 次。

因此一共：

`544`

次新模态 temporal solve。

---

## 17. Enrichment 内部真正最贵的是 residual-LS / LSMR，但其全局占比仍然很小

`_raw_fixed_point()` 的主要时间：

| 子项 | 调用次数 | `cumtime / s` |
|---|---:|---:|
| `solve_tower_residual_ls_spatial` | 523 | `12.236` |
| `_temporal_solve` | 544 | `2.521` |
| first import / `_find_and_load` | 1 | `2.194` |
| `_pair_change` | 523 | `0.874` |
| in-loop spatial orthogonalization | 544 | `0.095` |
| equilibrium `apply_spatial` | 544 | `0.041` |

前期 profile 还确认：

- `solve_tower_residual_ls_spatial`：`523` 次，`12.236 s`
- 其中 `scipy.sparse.linalg.lsmr`：`523` 次，`11.857 s`

因此 residual-LS spatial solve 内：

$$ \frac{11.857}{12.236}\approx96.9\% $$

的累计时间位于 LSMR。

但相对于整个 profiled LATIN solver：

$$ \frac{11.857}{699.091}\approx1.70\% $$

所以可以确认：

> **LSMR 是 enrichment 内部的主要空间求解成本，但不是整个 10 周期 LATIN-PGD 的主要 wall-time 瓶颈。**

另外，`_raw_fixed_point` 中 `2.194 s` 的 `_find_and_load` 是一次性 Python / SciPy import loading，应与重复算法成本区分。它不能解释为每次 fixed-point 都存在的 2.2 s 开销。

---

## 18. 60 次 full-basis temporal update 的实际 profiled 成本很小

此前源码审计确认：

- 每个 outer LATIN iteration 有一次 Trial-A fixed-basis full temporal update：39 次；
- 每次 successful enrichment 后，对 enlarged basis 再进行一次 full temporal re-optimization：21 次；

因此：

$$ 39+21=60 $$

次 full-basis temporal update。

本次 profiling 中：

- 39 次 Trial-A `update_tower_pgd_time_functions`：`2.423 s`
- 21 次 enrichment 后 enlarged-basis `update_tower_pgd_time_functions`：`0.925 s`

合计：

$$ 2.423+0.925=3.348\ \mathrm{s} $$

平均每次：

$$ \frac{3.348}{60}\approx0.0558\ \mathrm{s} $$

相对于 profiled LATIN solver：

$$ \frac{3.348}{699.091}\approx0.48\% $$

因此，本阶段能够明确修正前期的主要怀疑：

> **对于当前 10 周期、401 时间点、最终 rank 21 的实现，60 次 full-basis temporal update 并不是总 wall time 的主要来源。**

这也意味着此前讨论的：

- `320 × m` least-squares 是否直接压缩为 `m × m`；
- temporal update 是否采用另一种小型直接求解；
- full-history temporal update 是否成为首要优化目标；

虽然仍具有数学和实现研究价值，但从当前 10 周期总耗时诊断看，都不应列为第一优先级。

---

## 19. 本阶段对前期成本假设的重新排序

### 19.1 已被 profiling 明确弱化的假设

#### 假设 A：反复 full-basis temporal update 是第一大耗时源

当前数据：

- 60 次合计约 `3.348 s`
- profiled LATIN solver 占比约 `0.48%`

因此目前不支持其作为第一瓶颈。

#### 假设 B：PGD enrichment / LSMR 是第一大耗时源

当前数据：

- 21 次 enrichment：`19.037 s`
- 523 次 LSMR：`11.857 s`
- 占整个 profiled LATIN solver 分别约 `2.72%` 和 `1.70%`

因此目前也不支持其作为第一瓶颈。

### 19.2 当前得到最强支持的假设

#### 假设 C：Local whole-history constitutive projection 是第一瓶颈

当前数据：

- 39 次 Local stage：`668.605 s`
- profiled LATIN solver 占比：`95.64%`
- 4,992,000 次 RK4 material time-step integration
- 24,972,480 次 `local_rates_from_forces`
- 大量 scalar NumPy / Python hot-loop 调用

因此，在当前 profiled implementation 中，该假设得到非常强的支持。

---

## 20. 对“低秩为什么没有转化为低 wall time”的科研认识更新

此前 1–10 周期结果已经确认：

- 时间点数从 41 增长到 401；
- PGD rank 仅从 11 增长到 21；
- 因而结构响应在空间—时间表示意义上保持较好的低秩可压缩性。

但是本阶段进一步揭示：

> **低秩可表示性并不自动意味着整个 LATIN-PGD 求解过程低成本。**

当前实现中：

1. Global correction / PGD basis 可以低秩表达；
2. 但 Local nonlinear constitutive update 仍在完整 `(time × material point)` 网格上计算；
3. 每个 outer LATIN iteration 又重新执行一次完整 Local history；
4. 因而低 rank 并没有压缩当前 Local stage 的主要工作量。

可以将当前核心矛盾概括为：

$$ \boxed{\text{Global response is low-rank, but Local constitutive history is still full-order and repeatedly recomputed.}} $$

用中文更准确地说：

> **当前实现已经证明塔筒响应具有低秩表达能力，但 Local 本构历史积分仍保持全时间、全材料点、逐 LATIN 迭代重复计算，因此“低秩”尚未转化为“低总计算成本”。**

这一认识比简单归因于“PGD rank 仍然太高”或“LSMR 太慢”更符合当前实际证据。

---

## 21. 必须保持的证据边界

### 21.1 已经确认的事实

- profiled 10-cycle LATIN 数值路径与正式 benchmark 一致；
- 39 次 Local stage；
- 4,992,000 次 RK4 material time-step integration；
- 19,968,000 次 RK4 stage state-rate evaluation；
- 24,972,480 次 `local_rates_from_forces`；
- 21 次 enrichment；
- 523 次 residual-LS spatial solve / fixed-point spatial update；
- 544 次 new-mode temporal solve；
- 523 次 LSMR；
- 60 次 full-basis temporal update；
- `cProfile` 下 Local stage 占 LATIN solver 约 95.6%；
- `cProfile` 下 enrichment 约 2.7%；
- `cProfile` 下 60 次 full-basis temporal update 合计约 0.48%；
- `cProfile` 下 LSMR 约 1.70%。

### 21.2 目前只能作为“profiled run”结论，不能直接外推的内容

以下绝对百分比仍需低侵入式计时确认：

- 未加 profiler 时 Local 是否仍然精确占 95.6%；
- 未加 profiler 时 scalar `np.clip` 是否仍然占约 36.9% 的 LATIN solver；
- 未加 profiler 时各子函数的绝对秒数是否保持同样比例。

原因是 `cProfile` 对高频 Python 函数调用存在明显扰动。

### 21.3 尚不能做出的结论

当前不能声称：

- Local stage 的高成本是 LATIN-PGD 理论本身不可避免的；
- 原论文作者的实现也存在同样的 scalar NumPy 开销；
- 原论文采用的 DG0 temporal update 可以解决当前 Local bottleneck；
- 改用 whole-time temporal solve 会显著改善当前总效率；
- 直接优化 LSMR 会显著缩短总 wall time；
- 只修改 `np.clip` 就能恢复 LATIN 对 FOM 的速度优势；
- 当前 profiling 已经等价于最终论文级 performance benchmark。

这些都需要额外数值证据。

---

## 22. 与原论文和当前工程实现的区分

本阶段的耗时结论属于：

> **当前 `Offshore_Wind_Turbine_LATIN_PGD` 塔筒 Python implementation 的源码级和实测 profiling 结论。**

不能把它直接表述为原论文 LATIN-PGD 方法的普遍复杂度结论。

尤其需要保持以下区分：

- 原论文明确包含 Local / Global LATIN 结构和 PGD reduction；
- 当前塔筒实现使用经过一维验证迁移而来的 RK4 local material-point integrator；
- 当前 temporal update 为项目中的 backward-Euler-like weighted least-squares 实现，而非已经证明与原论文 DG0 严格等价的实现；
- 当前 spatial enrichment 使用 `residual_ls + matrix-free LSMR`，属于当前项目实现选择；
- 当前 profiling 反映的是这些工程选择共同作用下的实际 Python 运行成本。

因此，后续汇报中最稳妥的表述是：

> **“在当前塔筒 LATIN-PGD Python 实现中，10 周期 profiling 显示 Local 本构历史积分是绝对主耗时模块；PGD enrichment、LSMR 和 full-basis temporal update 的占比明显较小。该结论反映当前实现，而不是对 LATIN-PGD 理论方法本身的普遍效率判定。”**

---

## 23. 当前最值得关注的优化层级

基于本阶段 profiling，后续优化问题应该按层级区分，而不是混在一起。

### 23.1 实现级、低风险优化候选

不改变理论公式和算法路径，例如：

- 避免 hot loop 内对 scalar 重复调用通用 NumPy API；
- 将常数如 machine epsilon 预先缓存；
- 减少不必要的 scalar property access；
- 对 q-point 维度做向量化；
- 考察 Numba / compiled constitutive kernel；
- 将有限性检查、clip 等安全逻辑以等价但更低开销的方式实现。

这些属于性能工程问题。

### 23.2 算法级、高价值但高风险问题

例如：

- Local stage 是否必须在每个 outer LATIN iteration 从头完整积分 whole history；
- 是否能够利用 previous Local solution / warm-start；
- 是否能在 LATIN framework 下减少重复 constitutive history integration；
- 是否能进行时间多尺度、周期跳跃或其他高周疲劳专用处理。

这些问题可能带来更大加速，但必须先证明不改变理论正确性和数值精度。

### 23.3 当前不应列为第一优先级的优化

根据本次 10 周期 profile：

- 仅优化 full-basis temporal LS；
- 仅将 `320 × m` temporal least-squares 改写为直接 `m × m` normal system；
- 仅优化 LSMR；
- 仅优化 PGD orthogonalization。

这些可能有局部收益，但目前不是总 wall time 的主要控制项。

---

## 24. 下一阶段必须先做低侵入式 wall-time 确认

由于 `cProfile` 明显增加函数调用开销，本阶段不能直接以 `95.6%` 作为最终 production wall-time 占比。

因此下一阶段最合理的任务是：

> **FOM-3D-7：10 周期 LATIN-PGD 低侵入式顶层 wall-time 分项确认。**

原则：

- 不修改核心算法；
- 不改变数值路径；
- 不运行新的 100-cycle LATIN；
- 只在少数顶层模块外部使用 `time.perf_counter()`；
- 避免对数百万级材料函数逐调用打点；
- 重点记录：
  - Local stage 总时间；
  - search direction；
  - frozen global data；
  - Trial-A full temporal update；
  - enrichment 总时间；
  - Trial construction / evaluation；
  - solver 总时间；
- 最终比较低侵入式分项与本次 `cProfile` 排名是否一致。

下一阶段的 PASS 条件建议为：

1. 数值求解路径仍保持 `39 iterations / 60 trials / rank 21 / final xi` 一致；
2. 各顶层分项 wall time 可近似加回 LATIN solver 总 wall time；
3. Local stage 是否仍为第一耗时模块得到独立确认；
4. 若 Local 仍绝对占优，再进入 Local kernel 的低风险优化研究；
5. 若低侵入式结果与 `cProfile` 比例显著不同，则先解释 profiler bias，不直接做算法优化。

---

## 25. 本阶段的正式结论

本阶段 10 周期函数级 profiling 得到以下核心结论。

1. 当前 LATIN-PGD 10 周期数值路径保持稳定：39 次 outer LATIN iteration、60 次 trial、最终 rank 21、最终 `xi = 8.941607234831e-06`。

2. `cProfile` 显著增加绝对 wall time，因此本阶段主要用于定位相对热点，而不是替代正式 timing benchmark。

3. 在 profiled run 中，39 次 `solve_tower_local_stage()` 累计 `668.605 s`，约占 LATIN solver `95.64%`，远高于所有其他模块。

4. Local stage 的高成本与算法结构严格一致：39 次 outer LATIN iteration × 320 个材料点 × 400 个时间增量 = `4,992,000` 次 RK4 material time-step integration。

5. RK4 四阶段进一步产生 `19,968,000` 次 `_local_state_rate()`，再加正式时间点速率重算后，`local_rates_from_forces()` 总调用达到 `24,972,480` 次。

6. `local_rates_from_forces()` 是当前最显著的 constitutive kernel hotspot；其内部 profiled 成本中，scalar `np.clip`、`np.finfo` 等高频 NumPy 调用占据明显比例，表明当前实现存在典型的“scalar NumPy in Python hot loop”性能问题。

7. 21 次 PGD enrichment 总计仅约 `19.037 s`，其中 523 次 residual-LS spatial solve / LSMR 是 enrichment 的主要内部成本，但 LSMR 相对于整个 LATIN solver 仅约 `1.70%`。

8. 60 次 full-basis temporal update 在 profiled run 中合计仅约 `3.348 s`，约占 LATIN solver `0.48%`。因此，反复 temporal full-basis re-optimization 不是当前 10 周期实现的第一耗时来源。

9. 本阶段最重要的科研认识是：当前塔筒响应确实具有较好的低秩可表示性，但 Local nonlinear constitutive history 仍然在完整 `(time × material point)` 网格上、随每个 LATIN outer iteration 重复计算，因此低 rank 尚未转化为低总 wall time。

10. 下一阶段不应立即修改核心算法，也不应继续 100-cycle LATIN，而应先用低侵入式 `perf_counter` 顶层计时独立确认 Local stage 在未加 profiler 条件下是否仍然占绝对主导。

---

## 26. 当前阶段状态

可以将本阶段定义为：

> **FOM-3D-6：10 周期 LATIN-PGD 函数级 profiling 与主要成本链定位 —— PASS。**

当前最重要的已确认成本链为：

$$ \boxed{\text{outer LATIN} \rightarrow \text{whole-history Local stage} \rightarrow \text{RK4 material integration} \rightarrow \text{constitutive state-rate evaluation}} $$

当前最主要的待验证结论为：

> **Local stage 是否在低侵入式真实 wall-time 测量中仍然保持第一主导项。**

在完成该验证以前，不对核心 Local algorithm、temporal solver、PGD enrichment 或 LSMR 做生产代码优化。

---

## 27. 建议 Git 记录

本阶段建议新增正式阶段总结：

`docs/2026-08-30-tower-latin-pgd-10cycle-function-level-profiling-stage-summary.md`

建议在确认文档内容无误后单独提交，例如：

```text
docs: summarize 10-cycle LATIN function-level profiling
```

本阶段 profiling 数据文件：

`tower_10cycle_profile.pstats`

属于本地诊断产物。是否纳入 Git 应根据仓库对二进制/大体积 diagnostic 文件的管理原则另行决定；阶段结论不依赖将 `.pstats` 提交到仓库。

