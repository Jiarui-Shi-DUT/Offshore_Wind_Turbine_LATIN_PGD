# FOM-OPT3：塔筒全阶模型 RK4 状态有限性检查标量化优化

## 1. 阶段定位

本阶段继续服务于塔筒全阶模型公平效率基线重建，不改变全阶模型的控制方程、材料本构、时间离散、Newton 求解框架或截面数值切线算法。

FOM-OPT3 的目标非常明确：

> 将 RK4 单步积分末端针对 4 分量内部变量状态向量的 NumPy 有限性检查，改写为标量 `math.isfinite()` 检查，减少高频小数组上的 NumPy 调度开销。

因此，本优化属于纯实现层优化，而不是算法层或模型层修改。

---

## 2. FOM-OPT2 后重新 profiling 的发现

FOM-OPT2 生产版本 10 周期 profiling 得到：

| 指标 | FOM-OPT2 profiling 结果 |
| --- | ---: |
| profiling 总耗时 | 41.339223 s |
| 函数调用总数 | 75,034,584 |
| `_numerical_tangent` 累计时间 | 30.028 s |
| `rk4_step` 累计时间 | 26.004 s |
| `_scalar_state_rate_components` 累计时间 | 13.327 s |
| `evaluate_state` 累计时间 | 5.814 s |
| `np.all(np.isfinite(...))` 相关累计时间 | 约 4.769 s |

数值签名保持为：

- 存储时间点：401；
- 节点位移张量：`(401, 33)`；
- 纤维应变张量：`(401, 10, 2, 16)`；
- 纤维应力张量：`(401, 10, 2, 16)`；
- 纤维状态张量：`(401, 10, 2, 16, 4)`；
- 最大 Newton 迭代数：4；
- 最大残差范数：`8.761245319024e-03`；
- 最终最大塑性应变绝对值：`4.954213909559e-04`；
- 最终最大损伤：`8.634555499717e-03`；
- 临界位置：`(0, 0, 4)`。

这次 profiling 表明，FOM-OPT2 已显著压低 RK4 内部的小数组构造和函数调度成本，但 RK4 尾部对长度为 4 的状态向量执行 `np.all(np.isfinite(new_state))` 仍然形成了可观的高频 NumPy 调度开销。

---

## 3. 为什么先不修改中心差分数值切线

FOM-OPT2 后，`_numerical_tangent` 仍然是累计时间最大的上层函数。

当前每次截面试算保持：

$$ N_{\mathrm{fiber\ integration}} = 5 N_{\mathrm{section\ trial}} $$

即：

- 1 次真实截面响应积分；
- 2 次轴向应变中心差分扰动；
- 2 次曲率中心差分扰动。

因此，约 80% 的材料积分由数值切线扰动产生。

但是，中心差分切线属于全阶求解算法结构的一部分。如果直接减少扰动次数、改用解析切线或改变切线构造方式，就不再是与 LATIN-PGD 中标量化、常数缓存、小数组消除同一类别的纯实现优化。

所以本阶段继续遵循以下原则：

> 在修改中心差分切线之前，优先消除其内部仍然存在的、数学上严格等价的高频实现开销。

---

## 4. FOM-OPT3 的具体修改

FOM-OPT2 中，RK4 完成四个阶段后先构造：

```python
new_state = np.array([...], dtype=np.float64)
```

随后执行：

```python
if not np.all(np.isfinite(new_state)):
    raise FloatingPointError(...)
```

FOM-OPT3 改为先得到四个标量：

```python
new_eps_p
new_alpha
new_r_bar
new_damage
```

并执行：

```python
if not (
    math.isfinite(new_eps_p)
    and math.isfinite(new_alpha)
    and math.isfinite(new_r_bar)
    and math.isfinite(new_damage)
):
    raise FloatingPointError(...)
```

检查通过后再构造最终长度为 4 的 `np.array`。

经典 RK4 更新关系保持不变：

$$ \vec{z}_{n+1} = \vec{z}_n + \frac{\Delta t}{6}\left(\vec{k}_1 + 2\vec{k}_2 + 2\vec{k}_3 + \vec{k}_4\right) $$

其中各阶段状态速率、本构方程、损伤截断、时间步长和应变插值均未改变。

---

## 5. 单周期内存 A/B 验证

单周期基准中：

| 指标 | FOM-OPT2 基线 A | FOM-OPT3 候选 B |
| --- | ---: | ---: |
| wall-time | 3.568918 s | 3.114895 s |
| A/B | 1.145759 | |
| 时间降低 | | 12.722% |

完整场比较：

| 场量 | 相对二范数误差 | 最大绝对误差 |
| --- | ---: | ---: |
| 节点位移 | 0 | 0 |
| 纤维应变 | 0 | 0 |
| 纤维应力 | 0 | 0 |
| 纤维内部变量 | 0 | 0 |

同时：

- Newton 迭代序列完全一致；
- 残差范数逐项完全一致；
- 临界位置均为 `(0, 0, 4)`；
- 最大 Newton 迭代数均为 4；
- 最终最大损伤均为 `1.964038879640e-03`；
- 最终最大塑性应变绝对值均为 `1.387953504904e-04`。

结论：单周期 A/B 验证 PASS。

---

## 6. 10 周期内存 A/B 验证

10 周期匹配塔筒算例中：

| 指标 | FOM-OPT2 基线 A | FOM-OPT3 候选 B |
| --- | ---: | ---: |
| wall-time | 29.795887 s | 26.154806 s |
| A/B | 1.139213 | |
| 时间降低 | | 12.220% |

完整场仍保持严格一致：

| 场量 | 相对二范数误差 | 最大绝对误差 |
| --- | ---: | ---: |
| 节点位移 | 0 | 0 |
| 纤维应变 | 0 | 0 |
| 纤维应力 | 0 | 0 |
| 纤维内部变量 | 0 | 0 |

其他关键量：

- Newton 迭代序列完全一致；
- 残差范数逐项完全一致；
- 临界位置均为 `(0, 0, 4)`；
- 最大 Newton 迭代数均为 4；
- 最终最大损伤均为 `8.634555499717e-03`；
- 最终最大塑性应变绝对值均为 `4.954213909559e-04`。

结论：10 周期 A/B 验证 PASS。

---

## 7. 生产版本落地与直接验证

生产代码仅修改：

```text
material/viscoplastic_damage_1d.py
```

修改内容仅包括：

1. 新增 `import math`；
2. 将 RK4 尾部 4 分量状态的 NumPy 有限性检查改写为标量 `math.isfinite()`；
3. 在有限性检查通过后再构造最终 `np.array`。

生产版 10 周期直接运行得到：

```text
stored time points          = 401
nodal displacement shape   = (401, 33)
fiber strain shape         = (401, 10, 2, 16)
fiber stress shape         = (401, 10, 2, 16)
fiber state shape          = (401, 10, 2, 16, 4)
max Newton iterations      = 4
max residual norm          = 8.761245319024e-03
final max |plastic strain| = 4.954213909559e-04
final max damage           = 8.634555499717e-03
critical location          = (0, 0, 4)
```

单次非 profiling 运行时间为：

```text
26.020437 s
```

该时间与 10 周期候选 B 的 `26.154806 s` 量级一致，说明优化收益已经稳定落入生产代码。

但是该值仍然只是一次生产一致性验证，不作为最终 FOM 与 LATIN-PGD 公平效率比较中的正式 wall-time。

---

## 8. 回归测试

生产代码落地后依次完成：

| 测试范围 | 结果 |
| --- | --- |
| 一维材料模型 | 3 passed |
| 黏塑性纤维截面 | 7 passed |
| 黏塑性塔筒系统 | 7 passed |
| 非对称循环塔筒专项 | 17 passed, 10 warnings |
| 完整测试集 | 318 passed, 10 warnings |

完整测试集耗时：

```text
67.48 s
```

另外，在 Windows + conda 环境下直接执行 `pytest` 时曾出现仓库根目录没有进入模块搜索路径的问题，表现为：

```text
ModuleNotFoundError: No module named 'material'
```

改用：

```powershell
python -m pytest
```

后测试正常执行。因此后续本项目回归测试统一优先使用 `python -m pytest`。

---

## 9. FOM-OPT3 的科研含义

FOM-OPT3 再次说明，在 Python/NumPy 实现的高频材料积分内核中：

> 对只有几个标量分量的小状态向量反复调用 NumPy 通用接口，可能产生与真实浮点计算同量级甚至更高的调度成本。

因此，FOM-OPT1、FOM-OPT2 和 FOM-OPT3 共同完成了一条清晰的实现优化链：

- FOM-OPT1：标量损伤截断；
- FOM-OPT2：RK4 与材料速率标量化；
- FOM-OPT3：RK4 状态有限性检查标量化。

这些修改均未减少自由度、单元、积分点、纤维数量、时间增量或材料内部变量，也未引入降阶近似。

因此，它们应被归类为：

> 全阶模型实现成熟度提升，而不是模型降阶收益。

---

## 10. 对公平效率比较的影响

当前仍不能使用旧 FOM 与优化后 LATIN-PGD 的时间比作为 LATIN-PGD 的算法加速比。

最终应区分：

$$ S_{\mathrm{FOM,impl}} = \frac{t_{\mathrm{FOM,old}}}{t_{\mathrm{FOM,opt}}} $$

$$ S_{\mathrm{LATIN,impl}} = \frac{t_{\mathrm{LATIN,old}}}{t_{\mathrm{LATIN,opt}}} $$

真正用于方法比较的公平加速比仍定义为：

$$ S_{\mathrm{fair}} = S_{\mathrm{method}} = \frac{t_{\mathrm{FOM,opt}}}{t_{\mathrm{LATIN,opt}}} $$

其中 `FOM,opt` 必须是在与 LATIN-PGD 相近实现成熟度下得到的优化全阶模型，而不是早期包含大量 Python/NumPy 额外开销的旧实现。

---

## 11. 当前 checkpoint

FOM-OPT3 生产代码提交：

```text
c2a1406  perf: optimize tower FOM RK4 finite-state validation
```

当前阶段结论：

> FOM-OPT3 已通过单周期 A/B、10 周期 A/B、生产 10 周期数值签名验证以及完整 318 项回归测试，可以作为新的稳定生产 checkpoint。

---

## 12. 下一步

下一步不应立即修改中心差分数值切线。

首先应针对 FOM-OPT3 生产版本重新执行 10 周期 profiling，使用新的独立诊断文件，例如：

```text
tower_fom_10cycle_profile_opt3.py
tower_fom_10cycle_profile_opt3.pstats
```

然后重新判断：

- RK4 内核还剩多少纯实现层热点；
- `evaluate_state`、`safe_damage`、应变插值闭包或响应对象构造是否仍值得优化；
- 中心差分数值切线在消除这些实现开销之后是否真正成为不可回避的主导成本。

只有完成新的 profiling 后，才能决定是否存在 FOM-OPT4，以及 FOM 优化何时应当停止并冻结为最终 `FOM,opt`。
