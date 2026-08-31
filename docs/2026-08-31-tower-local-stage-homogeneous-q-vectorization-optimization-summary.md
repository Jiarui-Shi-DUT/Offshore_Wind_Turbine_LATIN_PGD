# OPT-4 阶段总结：同质材料点方向批量化局部本构积分与逐点参考路径等价性验证

**日期：2026-08-31**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**仓库：`Jiarui-Shi-DUT/Offshore_Wind_Turbine_LATIN_PGD`**  
**当前分支：`perf/tower-local-stage-optimization`**  
**OPT-4 代码提交：`423a9c4`，提交信息：`perf: vectorize homogeneous tower local stage`**

---

# 1. 阶段定位

OPT-4 是在 OPT-1、OPT-2 和 OPT-3 完成之后，对塔筒 LATIN-PGD 局部本构阶段开展的第四轮实现级性能优化。

本阶段不修改 LATIN-PGD 外层算法，不修改 PGD 富集策略，不修改搜索方向，不修改材料本构方程，也不修改四阶龙格-库塔时间积分公式。

本阶段只改变局部本构阶段的计算组织方式：

- 时间方向继续严格顺序推进；
- 同质材料条件下，不同材料点沿材料点方向进行 NumPy 批量计算；
- 异质材料条件下，继续使用原逐材料点路径；
- 原逐点实现被保留为参考与后备路径。

本阶段正式定义为：

> **OPT-4：同质材料点方向批量化局部本构积分与逐点参考路径等价性验证**

最终判定：

> **PASS**

---

# 2. 为什么 OPT-4 具有明确的数值基础

当前塔筒中的标准材料点索引为：

$$ q \leftrightarrow (\mathrm{element},\mathrm{Gauss\ point},\mathrm{fiber}) $$

塔筒所有局部材料历史采用统一二维存储：

$$ (N_t,N_q) $$

其中 $N_t$ 为时间点数，$N_q$ 为材料点数。

当前正式塔筒 benchmark 使用：

- 10 个梁柱单元；
- 每单元 2 个高斯点；
- 每高斯点 16 个周向纤维；
- 1 个径向层。

因此材料点总数为：

$$ N_q = 10 \times 2 \times 16 \times 1 = 320 $$

在 LATIN 局部阶段中，当前采用的上升搜索方向使热力学力历史在局部投影过程中保持给定。

因此不同材料点之间没有局部本构耦合。

对固定时间位置 $t_n$，不同材料点状态：

$$ z_{1,n},z_{2,n},\ldots,z_{N_q,n} $$

可以彼此独立地计算。

但是对任意单个材料点 $q$，内部变量具有历史依赖：

$$ z_{q,n} \longrightarrow z_{q,n+1} $$

所以时间方向不能被打乱，也不能把同一材料点的全部时间步直接并行化。

OPT-4 的核心计算结构因此为：

$$ \boxed{\mathrm{time\ sequential} + \mathrm{q\ direction\ batched}} $$

对应中文含义是：

> **时间方向保持顺序，材料点方向进行批量计算。**

---

# 3. OPT-4 前的局部阶段执行结构

OPT-3 完成后，`solve_tower_local_stage()` 的主体仍然是外层遍历材料点、内层遍历时间步。

概念结构为：

```text
for q in material_points:
    initialize internal state of q
    compute initial rates of q
    compute initial elastic strain of q

    for step in time_steps:
        advance q by RK4
        store internal state
        recompute local rates
        recompute elastic strain
```

也就是先完整计算第 1 个材料点的全部时间历程，再计算第 2 个材料点，以此类推。

对于 10 周期正式算例：

- `Nt = 401`；
- 时间区间数为 400；
- `Nq = 320`；
- LATIN 外层迭代次数为 39。

因此单材料点材料状态时间推进总次数为：

$$ 39 \times 400 \times 320 = 4{,}992{,}000 $$

即约 499 万次材料状态时间推进。

这一计数与此前函数级耗时分析中 `_integrate_one_local_step()` 的实际调用次数严格一致。

---

# 4. OPT-4 的核心重构

同质材料情况下，OPT-4 将执行顺序重构为：

```text
for step in time_steps:
    advance all q-points together
```

即每一个时间步只进入一次 Python 层循环，然后把全部 320 个材料点的状态作为数组统一送入 NumPy 运算。

时间顺序仍然保持：

$$ \mathbf Z_n \longrightarrow \mathbf Z_{n+1} $$

其中批量内部状态写为：

$$ \mathbf Z_n \in \mathbb R^{4\times N_q} $$

四行分别对应：

$$ [\varepsilon^p,\alpha,\bar r,D]^T $$

因此 OPT-4 不是改变材料本构，也不是跳过材料点，而是把原来的大量 Python 标量调度转换为数组批量计算。

在 10 周期、39 次 LATIN 迭代下，Python 层材料状态推进循环的概念数量由：

$$ 39 \times 400 \times 320 = 4{,}992{,}000 $$

变为：

$$ 39 \times 400 = 15{,}600 $$

对应循环数量比例为：

$$ \frac{4{,}992{,}000}{15{,}600} = 320 $$

需要特别强调：

> **这里的 320 倍只是 Python 层循环粒度的变化，不是物理计算量减少 320 倍，也不是实际运行时间必然加速 320 倍。**

全部 320 个材料点的本构运算仍然存在，只是被 NumPy 以更合适的数组粒度执行。

---

# 5. 为什么当前正式 benchmark 可以进入同质材料快速路径

正式效率 benchmark 在公共输入中只创建一次：

```python
material = MaterialParameters()
```

随后同一个材料对象被传给：

- FOM；
- 塔筒系统；
- LATIN 初始化；
- `solve_tower_latin_pgd()`。

因此当前正式 benchmark 的 320 个材料点采用完全相同的一组材料参数。

这使得当前 benchmark 可以直接进入同质材料快速路径。

但是原有 API 允许每个材料点使用不同的 `MaterialParameters`。

因此 OPT-4 不能简单地删除逐点算法。

---

# 6. 双路径设计

OPT-4 最终采用同质材料快速路径和异质材料后备路径并存的结构。

概念结构为：

```text
solve_tower_local_stage
        |
        |-- material parameters are homogeneous
        |       |
        |       -> batched homogeneous path
        |
        |-- material parameters are heterogeneous
                |
                -> validated pointwise path
```

新增材料一致性判断：

```python
_materials_are_homogeneous(...)
```

同质材料快速路径：

```python
_solve_tower_local_stage_homogeneous(...)
```

原逐材料点路径被完整保留为：

```python
_solve_tower_local_stage_pointwise(...)
```

这一设计具有两个目的。

第一，当前正式同质钢材塔筒可以获得批量化收益。

第二，原来允许不同材料点采用不同材料参数的通用能力没有被破坏。

因此原逐点路径同时承担：

- 异质材料情况下的正式后备路径；
- 批量路径的参考实现。

---

# 7. 向量化本构速率计算

OPT-4 新增：

```python
_vectorized_local_rates_from_forces(...)
```

该函数是原标量函数：

```python
local_rates_from_forces(...)
```

在材料点方向上的数组版本。

数学公式保持不变。

## 7.1 损伤限幅

原标量逻辑与批量逻辑均满足：

$$ D_{\mathrm{safe}} = \min(\max(D,0),D_{\max}) $$

OPT-4 批量路径使用 `np.clip()` 一次处理全部材料点。

这一点与 OPT-1 并不矛盾。

OPT-1 中 `np.clip()` 每次只处理一个 Python 标量，因此通用数组调度成本过高。

OPT-4 中 `np.clip()` 一次处理长度为 $N_q$ 的数组，属于 NumPy 适合的使用粒度。

## 7.2 相对有效应力

保持：

$$ \sigma_{\mathrm{rel}} = \frac{\sigma}{1-D_{\mathrm{safe}}} - \beta $$

## 7.3 各向同性硬化力

保持中间变量：

$$ q_R = \frac{\bar R\sqrt{\gamma}}{2R_{\infty}} $$

保持物理硬化力：

$$ R = R_{\infty}q_R(2-q_R) $$

## 7.4 屈服函数

保持：

$$ f = |\sigma_{\mathrm{rel}}| + \frac{a\beta^2}{2C} - R - \sigma_y $$

## 7.5 Norton 黏塑性乘子

保持：

$$ \dot p = K^{-n}\langle f\rangle_+^n $$

代码中继续使用 OPT-2 已缓存的：

```python
material.k_viscoplastic
```

因此数学参数没有变化。

## 7.6 流动方向

原标量实现中，当相对有效应力的绝对值处于机器精度范围内时，流动方向取零，否则取符号函数。

批量形式保持同一规则：

$$ \mathrm{dir} = 0 \quad \mathrm{if}\quad |\sigma_{\mathrm{rel}}|\le \varepsilon_{\mathrm{machine}} $$

其余情况：

$$ \mathrm{dir} = \operatorname{sign}(\sigma_{\mathrm{rel}}) $$

## 7.7 塑性、硬化与损伤速率

塑性应变率、运动硬化速率、各向同性硬化速率和损伤速率的数学形式均保持不变。

OPT-4 只将逐点标量计算改为长度为 $N_q$ 的数组计算。

---

# 8. 向量化单边弹性关系

OPT-4 新增：

```python
_vectorized_unilateral_elastic_strain(...)
```

拉应力条件下继续采用：

$$ \varepsilon^e = \frac{\sigma}{E(1-D)} $$

压应力条件下继续采用：

$$ \varepsilon^e = \frac{\sigma}{E(1-hD)} $$

因此材料模型中拉压不同的损伤弹性恢复关系没有改变。

---

# 9. 批量四阶龙格-库塔推进

OPT-4 新增：

```python
_integrate_homogeneous_material_points_one_step(...)
```

其输入状态由单材料点向量：

$$ z_n = [\varepsilon^p,\alpha,\bar r,D]^T $$

改为所有同质材料点共同组成的矩阵：

$$ \mathbf Z_n \in \mathbb R^{4\times N_q} $$

但四阶龙格-库塔公式完全保持。

第一阶段：

$$ K_1 = F(\mathbf Z_n,t_n) $$

第二阶段：

$$ K_2 = F(\mathbf Z_n+\frac{\Delta t}{2}K_1,t_n+\frac{\Delta t}{2}) $$

第三阶段：

$$ K_3 = F(\mathbf Z_n+\frac{\Delta t}{2}K_2,t_n+\frac{\Delta t}{2}) $$

第四阶段：

$$ K_4 = F(\mathbf Z_n+\Delta tK_3,t_{n+1}) $$

最终更新：

$$ \mathbf Z_{n+1} = \mathbf Z_n+\frac{\Delta t}{6}(K_1+2K_2+2K_3+K_4) $$

OPT-3 已经建立的三个不同外力采样位置继续保持：

$$ 0,\frac{1}{2},1 $$

其中 $K_2$ 和 $K_3$ 继续复用同一个中点外力状态。

因此：

> **OPT-4 不改变 RK4 数学格式，只改变 RK4 在材料点方向上的执行粒度。**

---

# 10. 本阶段正式代码修改

OPT-4 正式代码提交为：

```text
423a9c4
perf: vectorize homogeneous tower local stage
```

提交涉及：

```text
latin/tower_local_stage.py
tests/test_tower_local_stage_opt4_equivalence.py
```

Git 提交统计为：

```text
2 files changed, 596 insertions(+), 30 deletions(-)
create mode 100644 tests/test_tower_local_stage_opt4_equivalence.py
```

其中：

- `latin/tower_local_stage.py`：增加同质材料批量路径和异质材料逐点后备路径；
- `tests/test_tower_local_stage_opt4_equivalence.py`：增加 OPT-4 专门等价性测试。

---

# 11. OPT-4 专门等价性测试的设计

在正式修改批量路径之前，先建立独立测试：

```text
tests/test_tower_local_stage_opt4_equivalence.py
```

测试不是只使用零初始状态，而是专门包含：

- 非均匀时间步；
- 拉压应力反转；
- 不同材料点的非零初始塑性应变；
- 不同材料点的非零初始运动硬化变量；
- 不同材料点的非零初始各向同性硬化变量；
- 不同材料点的非零初始损伤；
- 全部材料点使用同一个 `MaterialParameters()`。

参考侧采用已经验证的一维逐点局部阶段：

```python
solve_local_stage(...)
```

塔筒侧采用：

```python
solve_tower_local_stage(...)
```

对 `LatinStateTower.MATERIAL_FIELD_NAMES` 中的全部材料场逐项比较。

比较容差为：

```python
rtol = 0.0
atol = 1.0e-14
```

---

# 12. OPT-4 修改前的测试基线

在尚未修改 `tower_local_stage.py` 时，专门等价性测试结果为：

```text
1 passed in 0.12s
```

这一步证明新增测试本身在 OPT-3 稳定实现上有效。

原有塔筒局部阶段相关测试基线为：

```text
5 passed in 0.15s
```

---

# 13. OPT-4 修改后的专门等价性测试

完成材料点批量化后，再次运行：

```powershell
python -m pytest .\tests\test_tower_local_stage_opt4_equivalence.py -q
```

结果：

```text
1 passed in 0.12s
```

因此在当前专门测试覆盖范围内，OPT-4 批量路径与已验证逐点参考路径在绝对容差 $10^{-14}$ 下保持一致。

可正式表述为：

> **在当前同质材料、非均匀时间步、拉压反转和非零初始内部变量测试中，OPT-4 批量局部阶段与逐点参考局部阶段保持机器精度量级的一致性。**

---

# 14. 原有塔筒局部阶段测试

运行：

```powershell
python -m pytest .\tests\test_tower_local_search_bridge.py .\tests\test_tower_local_search_bridge_integration.py -q
```

OPT-4 修改前：

```text
5 passed in 0.15s
```

OPT-4 修改后：

```text
5 passed in 0.14s
```

原有测试覆盖：

- 塔筒局部阶段与一维逐点实现的一致性；
- 不同材料点使用不同材料参数；
- 单材料对象广播；
- 输入状态不可变；
- 热力学力历史保持固定；
- 返回状态只读；
- 实际梁单元、高斯点和纤维材料点网格。

因此异质材料逐点后备路径没有发生退化。

---

# 15. 全仓库回归测试

运行：

```powershell
python -m pytest -q
```

结果：

```text
314 passed, 10 warnings in 102.74s
```

OPT-3 稳定版本原有测试数为 313。

增加的第 314 个测试正是：

```text
tests/test_tower_local_stage_opt4_equivalence.py
```

10 条 warning 仍然来自 Matplotlib 和 setuptools 的 `distutils` 弃用提示，与 OPT-4 无关。

因此：

> **OPT-4 没有破坏当前仓库已有功能。**

---

# 16. 1 周期正式数值验证

运行：

```powershell
python .\tower_asymmetric_efficiency_scaling_pilot.py --cycles 1
```

OPT-4 实测结果：

```text
cycles                    = 1
Nt                        = 41
Nq                        = 320
FOM analysis wall time    = 34.956793 s
LATIN setup wall time     = 0.127857 s
LATIN solver wall time    = 2.119564 s
LATIN total wall time     = 2.247421 s
sample FOM/LATIN ratio    = 15.554181
```

完整数值路径：

```text
termination_reason = converged
converged          = True
iterations         = 18
attempted          = 18
trial evaluations  = 29
PGD rank           = 11
modes added        = 11
final xi           = 7.918424536257e-06
```

其他关键量：

```text
FOM max Newton iterations = 4
FOM max |top displacement|= 1.173015956970e+00 m
LATIN initial eq residual = 1.715324872498e-06 N
```

这些收敛路径指标与 OPT-3 完全一致。

---

# 17. 1 周期相对 OPT-3 的性能提升

OPT-3 的 1 周期 LATIN solver 为：

$$ T_{\mathrm{solver,OPT3}} = 8.661120\ \mathrm{s} $$

OPT-4 为：

$$ T_{\mathrm{solver,OPT4}} = 2.119564\ \mathrm{s} $$

solver 耗时降低：

$$ \frac{8.661120-2.119564}{8.661120}\times100\% \approx 75.53\% $$

solver 增量加速约：

$$ \frac{8.661120}{2.119564} \approx 4.09 $$

OPT-3 的 1 周期 LATIN total 为：

$$ T_{\mathrm{total,OPT3}} = 8.790908\ \mathrm{s} $$

OPT-4 为：

$$ T_{\mathrm{total,OPT4}} = 2.247421\ \mathrm{s} $$

总时间降低约：

$$ 74.43\% $$

总时间增量加速约：

$$ 3.91 $$

---

# 18. 10 周期正式数值验证

运行：

```powershell
python .\tower_asymmetric_efficiency_scaling_pilot.py --cycles 10
```

OPT-4 实测结果：

```text
cycles                    = 10
Nt                        = 401
Nq                        = 320
FOM analysis wall time    = 293.417296 s
LATIN setup wall time     = 0.998693 s
LATIN solver wall time    = 31.701352 s
LATIN total wall time     = 32.700045 s
sample FOM/LATIN ratio    = 8.972994
```

完整数值路径：

```text
termination_reason = converged
converged          = True
iterations         = 39
attempted          = 39
trial evaluations  = 60
PGD rank           = 21
modes added        = 21
final xi           = 8.941607234831e-06
```

其他关键量：

```text
FOM max Newton iterations = 4
FOM max |top displacement|= 1.426595421862e+00 m
LATIN initial eq residual = 2.073353735670e-06 N
```

这些数值路径指标与 OPT-3 完全一致。

---

# 19. 10 周期相对 OPT-3 的性能提升

OPT-3 的 10 周期 LATIN solver 为：

$$ T_{\mathrm{solver,OPT3}} = 161.851950\ \mathrm{s} $$

OPT-4 为：

$$ T_{\mathrm{solver,OPT4}} = 31.701352\ \mathrm{s} $$

solver 耗时降低：

$$ \frac{161.851950-31.701352}{161.851950}\times100\% \approx 80.41\% $$

solver 增量加速约：

$$ \frac{161.851950}{31.701352} \approx 5.11 $$

OPT-3 的 10 周期 LATIN total 为：

$$ T_{\mathrm{total,OPT3}} = 162.864288\ \mathrm{s} $$

OPT-4 为：

$$ T_{\mathrm{total,OPT4}} = 32.700045\ \mathrm{s} $$

总时间降低约：

$$ 79.92\% $$

总时间增量加速约：

$$ 4.98 $$

---

# 20. OPT-1 至 OPT-4 的累计优化轨迹

10 周期 LATIN 总时间依次为：

```text
Original : 467.362 s
OPT-1    : 223.160690 s
OPT-2    : 174.717540 s
OPT-3    : 162.864288 s
OPT-4    : 32.700045 s
```

从原始实现到 OPT-4：

$$ 467.362\ \mathrm{s} \longrightarrow 32.700045\ \mathrm{s} $$

累计耗时降低约：

$$ 93.00\% $$

累计加速约：

$$ \frac{467.362}{32.700045} \approx 14.29 $$

因此当前可以确认：

> **在不改变 LATIN-PGD 数学算法的前提下，四轮实现级优化已经使当前 10 周期 LATIN 总时间降低约 93%，累计加速约 14.29 倍。**

---

# 21. 如何重新理解原先的效率交叉

原始实现下，10 周期结果约为：

$$ T_{\mathrm{FOM,original}} \approx 307.891\ \mathrm{s} $$

$$ T_{\mathrm{LATIN,original}} \approx 467.362\ \mathrm{s} $$

当时 LATIN-PGD 比 FOM 更慢，并且 1 至 10 周期效率缩放结果曾显示效率交叉位于 2 至 5 周期之间。

OPT-1 至 OPT-4 均没有修改 LATIN-PGD 数学算法。

但是 OPT-4 后，当前 10 周期结果变为：

$$ T_{\mathrm{FOM,current}} = 293.417296\ \mathrm{s} $$

$$ T_{\mathrm{LATIN,current}} = 32.700045\ \mathrm{s} $$

当前两套具体 Python 实现的实测时间比为：

$$ \frac{T_{\mathrm{FOM,current}}}{T_{\mathrm{LATIN,current}}} = 8.972994 $$

因此此前的效率交叉不能再被解释为 LATIN-PGD 方法本身的固有效率边界。

更严谨的科研表述为：

> **原始 Python 实现中曾在 2 至 5 周期附近出现 FOM 与 LATIN-PGD 的效率交叉，但经过四轮严格保持数学算法不变的实现级优化后，10 周期 LATIN-PGD 已显著快于当前 FOM。因此原先的效率交叉主要受到实现开销污染，不能直接视为 LATIN-PGD 方法的理论效率边界。**

---

# 22. 为什么当前 FOM/LATIN 比值不能直接解释为方法理论加速比

当前 10 周期实测：

$$ \frac{T_{\mathrm{FOM}}}{T_{\mathrm{LATIN}}} \approx 8.97 $$

这一结果是真实的当前实现实测值。

但是此前函数级耗时分析已经表明，当前 FOM 自身也存在明显的 Python 标量实现开销，例如高频标量损伤处理等。

因此当前 8.97 倍只能解释为：

> **当前两套具体 Python 实现之间的实测效率比。**

不能直接宣称：

> **LATIN-PGD 这一数值方法理论上必然比 FOM 快约 9 倍。**

如果后续需要开展严格的方法级效率比较，需要对 FOM 建立更公平的实现基线。

---

# 23. OPT-4 带来的第一项科研认识：局部阶段具有明确的空间独立性

OPT-4 的高收益不是偶然的代码技巧，而是利用了局部阶段本身的结构。

同一个材料点沿时间方向具有历史依赖，因此时间必须顺序推进。

但是在给定 LATIN 全局状态后，不同材料点之间的局部本构更新彼此独立。

因此当前局部阶段具有明确结构：

$$ \boxed{\mathrm{time\ sequential} + \mathrm{space\ independent}} $$

这一结构意味着材料点方向天然适合：

- NumPy 向量化；
- 编译优化；
- 多核并行；
- GPU 或其他批量执行后端。

OPT-4 当前只验证了 NumPy 批量化这一层。

---

# 24. OPT-4 带来的第二项科研认识：算法复杂度与 Python 实现耗时必须严格区分

OPT-4 相对 OPT-3 在 10 周期下使 LATIN total 从约 163 s 降到约 33 s。

期间数学算法没有变化。

因此可以确认当前性能损失中曾存在大量实现层开销，包括：

- 逐材料点 Python 循环；
- 高频标量函数调用；
- Python 与 NumPy 之间的小粒度切换；
- 标量数组调度；
- 重复的小对象构造。

所以后续分析必须坚持：

$$ \boxed{\mathrm{algorithmic\ complexity} \neq \mathrm{current\ Python\ wall\ time}} $$

这也是为什么任何“方法快慢”的科研结论都不能只由未优化代码的一次计时直接给出。

---

# 25. OPT-4 带来的第三项科研认识：NumPy 本身不是快慢问题，关键是使用粒度

OPT-1 曾证明：

> 对一个标量频繁调用通用 `np.clip()` 很低效。

OPT-4 则证明：

> 对 320 个材料点组成的数组统一调用 NumPy 运算可以获得显著收益。

两者并不矛盾。

真正的区别是数据粒度：

```text
one scalar per NumPy call
    -> high dispatch overhead

320 q-points per NumPy call
    -> suitable vectorized granularity
```

因此更准确的认识是：

> **NumPy 是否高效取决于一次调用处理的数据规模，而不是简单判断某个 NumPy 函数本身快或慢。**

---

# 26. 当前已经证明的内容

截至 OPT-4，可以正式确认以下事实。

1. 当前正式 1 周期和 10 周期 benchmark 使用同一组材料参数作用于全部 320 个材料点。
2. 同质材料情况下，材料点方向可以批量计算，同时保持时间方向严格顺序推进。
3. 专门等价性测试在 `rtol=0.0`、`atol=1.0e-14` 下通过。
4. 原有异质材料逐点路径继续通过一维参考对照。
5. 原有塔筒局部阶段相关测试 5 项全部通过。
6. 全仓库回归测试达到 `314 passed, 10 warnings`。
7. 1 周期 LATIN-PGD 收敛路径与 OPT-3 完全一致。
8. 10 周期 LATIN-PGD 收敛路径与 OPT-3 完全一致。
9. 10 周期 LATIN total 从 `162.864288 s` 降到 `32.700045 s`。
10. OPT-4 相对 OPT-3 的 10 周期总时间增量加速约为 4.98 倍。
11. 从最初实现到 OPT-4 的 10 周期累计加速约为 14.29 倍。

---

# 27. 当前仍然没有证明的内容

以下结论目前仍不能宣称。

1. 不能宣称 OPT-4 对所有材料模型都具有约 5 倍加速。
2. 不能宣称异质材料塔筒也具有与同质材料相同的批量化收益。
3. 不能把当前约 8.97 倍 FOM/LATIN 比值解释为两种数值方法的理论效率比。
4. 不能宣称当前局部阶段已经达到性能最优。
5. 不能宣称 LATIN-PGD 在任意更长周期问题中都会保持相同加速比例。
6. 不能宣称未来材料点数增加时收益一定严格按比例扩大。
7. 不能宣称下一轮优化一定还能获得与 OPT-4 同数量级的性能提升。
8. 不能因为当前 benchmark 为同质材料，就删除原有异质材料接口和后备路径。

---

# 28. OPT-4 的正式验证证据链

当前 OPT-4 的验证证据按层级可以整理为：

```text
Level 1
专门同质材料等价性测试
1 passed in 0.12 s
absolute tolerance = 1e-14

Level 2
原有塔筒局部阶段与真实材料点网格测试
5 passed in 0.14 s

Level 3
全仓库回归
314 passed, 10 warnings in 102.74 s

Level 4
1-cycle full LATIN-PGD benchmark
18 iterations
29 trial evaluations
rank = 11
xi = 7.918424536257e-06
LATIN total = 2.247421 s

Level 5
10-cycle full LATIN-PGD benchmark
39 iterations
60 trial evaluations
rank = 21
xi = 8.941607234831e-06
LATIN total = 32.700045 s
```

因此 OPT-4 的 PASS 不是只依靠性能计时，而是同时建立在局部等价性、模块回归、全仓库回归和完整多周期求解路径一致性之上。

---

# 29. 当前优化链

当前性能优化分支上的主要稳定节点为：

```text
900ddf2  docs: summarize 10-cycle LATIN function-level profiling

c27af15  perf: optimize scalar damage clipping in local stage
55d578e  docs: summarize scalar clipping optimization

e2d15b9  perf: cache material derived constants
2f0c15a  docs: summarize derived constant caching optimization

2a88b68  perf: flatten RK4 local-stage hot path
f9507b5  docs: summarize RK4 hot-path optimization

423a9c4  perf: vectorize homogeneous tower local stage
```

OPT-4 文档提交将在本总结人工检查通过后单独完成。

---

# 30. 下一阶段为什么不应直接照旧计划进入 OPT-5

OPT-4 已经把 10 周期 LATIN total 从约 163 s 大幅压缩到约 33 s。

因此 OPT-4 之前的性能热点比例已经失效。

原来占据绝大多数耗时的局部本构阶段在 OPT-4 后必然显著下降。

与此同时，原本次要的模块可能成为新的主要开销，例如：

- PGD 富集；
- 空间最小二乘求解；
- 搜索方向计算；
- 时间函数更新；
- 状态对象构造；
- 大数组复制；
- 其他 Python 层事务逻辑。

因此下一阶段不应直接假定“继续优化 Local 就一定最有效”。

---

# 31. 下一阶段推荐动作

OPT-4 文档归档完成后，下一步推荐重新进行 10 周期函数级耗时分析。

目标是建立新的性能结构：

```text
OPT-4 stable checkpoint
        |
        -> fresh 10-cycle function-level profiling
        |
        -> identify new first bottleneck
        |
        -> decide whether OPT-5 is justified
```

新的函数级耗时分析重点应回答：

- `solve_tower_local_stage()` 当前占 LATIN solver 的比例是多少；
- PGD enrichment 当前占比是多少；
- LSMR 当前占比是多少；
- 搜索方向计算当前占比是多少；
- 时间函数更新当前占比是多少；
- 状态构造与数组复制是否已经成为显著开销。

只有完成这一轮性能结构刷新后，才应决定下一优化对象。

---

# 32. 正式阶段结论

OPT-4 完成了同质材料塔筒局部本构阶段的材料点方向批量化。

批量化过程中：

- 时间历史依赖保持；
- 四阶龙格-库塔数学格式保持；
- 材料本构方程保持；
- LATIN 外层算法保持；
- PGD 富集策略保持；
- 异质材料逐点后备路径保持。

专门等价性测试、原有塔筒局部阶段测试、全仓库回归测试以及 1 周期和 10 周期正式 benchmark 均通过。

10 周期 LATIN 总时间由 OPT-3 的：

$$ 162.864288\ \mathrm{s} $$

降低至 OPT-4 的：

$$ 32.700045\ \mathrm{s} $$

相对降低约：

$$ 79.92\% $$

增量加速约：

$$ 4.98 $$

从最初实现到 OPT-4：

$$ 467.362\ \mathrm{s} \longrightarrow 32.700045\ \mathrm{s} $$

累计降低约：

$$ 93.00\% $$

累计加速约：

$$ 14.29 $$

10 周期求解路径继续保持：

```text
iterations         = 39
attempted          = 39
trial evaluations  = 60
PGD rank           = 21
modes added        = 21
final xi           = 8.941607234831e-06
```

因此正式判定：

> **OPT-4：同质材料点方向批量化局部本构积分与逐点参考路径等价性验证 — PASS**

当前稳定代码 checkpoint：

```text
423a9c4
perf: vectorize homogeneous tower local stage
```

下一步应先重新开展 OPT-4 后的 10 周期函数级耗时分析，再决定是否进入下一轮优化。

