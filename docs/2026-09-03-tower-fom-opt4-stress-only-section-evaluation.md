# FOM-OPT4：塔筒全阶模型纤维截面应力专用求值优化

## 1. 阶段定位

本阶段继续服务于塔筒全阶模型公平效率基线重建。

FOM-OPT4 不改变：

- 全阶有限元离散；
- 纤维数量、Gauss 积分点数量和时间增量；
- 一维黏塑性损伤本构；
- RK4 时间积分；
- 截面中心差分数值切线；
- Newton 迭代；
- 收敛容限；
- 非对称循环加载历史。

本阶段只针对纤维截面中的一个高频实现冗余：

> 当调用点最终只需要应力时，不再调用返回完整材料诊断量的 `evaluate_state()`，而直接调用已有的 `stress_from_state()`。

因此，FOM-OPT4 属于纯实现层优化，不属于模型降阶或求解算法改变。

---

## 2. FOM-OPT3 后重新 profiling

FOM-OPT3 生产版本完成后，对 10 周期匹配非对称塔筒全阶模型重新进行 profiling。

关键结果为：

| 指标 | FOM-OPT3 profiling |
| --- | ---: |
| profiling 总耗时 | 35.908981 s |
| 函数调用总数 | 68,365,784 |
| `_numerical_tangent` 累计时间 | 25.691 s |
| `rk4_step` 累计时间 | 20.694 s |
| `_scalar_state_rate_components` 累计时间 | 13.271 s |
| `evaluate_state` 累计时间 | 5.635 s |
| `safe_damage` 累计时间 | 2.037 s |

与 FOM-OPT2 profiling 相比：

- profiling 总耗时由 `41.339223 s` 降至 `35.908981 s`；
- 降低约 `13.14%`；
- 函数调用数由约 `75.03 million` 降至约 `68.37 million`；
- FOM-OPT3 针对 `np.all(np.isfinite(...))` 的优化热点已经基本消失。

这说明 FOM-OPT3 的优化方向成立，同时新的实现热点已经暴露出来。

---

## 3. 新热点：只为应力调用完整 `evaluate_state()`

profiling 中：

```text
evaluate_state
1,798,720 calls
cumtime = 5.635 s
```

进一步检查纤维截面生产源码发现，在两个高频位置：

1. `_integrate_fibers()`；
2. `_response_at_known_state()`；

调用形式均为：

```python
evaluate_state(...)[0]
```

也就是说，调用点只使用返回值中的第一个分量，即应力。

但是 `evaluate_state()` 除应力外，还会继续计算：

- 随动硬化相关量；
- 各向同性硬化相关量；
- 有效相对应力；
- 屈服函数；
- 损伤能量释放率。

这些结果在上述两个纤维截面调用点中均未被使用。

因此形成了一个明确的纯实现层冗余：

> 为获得一个应力标量，重复计算整套材料诊断量。

---

## 4. FOM-OPT4 的优化思路

材料模块中已经存在独立的：

```python
stress_from_state(...)
```

该函数根据：

- 总应变；
- 塑性应变；
- 损伤变量；
- 材料参数；

按照与原 `evaluate_state()` 相同的单轴损伤弹性关系计算应力。

因此，FOM-OPT4 不引入新的应力公式，只改变调用路径：

原路径：

```text
fiber section
    → evaluate_state()
        → stress
        → hardening
        → yield function
        → damage energy release rate
    → only use [0]
```

优化后：

```text
fiber section
    → stress_from_state()
    → directly obtain stress
```

这意味着删除的是未被消费的诊断量计算，而不是删除物理量或改变材料状态演化。

---

## 5. 单周期 A/B pilot

单周期匹配塔筒算例中：

| 指标 | FOM-OPT3 基线 A | FOM-OPT4 候选 B |
| --- | ---: | ---: |
| wall-time | 3.102259 s | 2.759861 s |
| A/B | 1.124063 | |
| 时间降低 | | 11.037% |

完整场比较：

| 场量 | 相对二范数误差 | 最大绝对误差 |
| --- | ---: | ---: |
| 节点位移 | 0 | 0 |
| 纤维应变 | 0 | 0 |
| 纤维应力 | 0 | 0 |
| 纤维内部变量 | 0 | 0 |

求解器和状态签名：

- Newton 迭代序列完全一致；
- 残差范数逐项完全一致；
- 临界位置均为 `(0, 0, 4)`；
- 最大 Newton 迭代数均为 4；
- 最终最大损伤均为 `1.964038879640e-03`；
- 最终最大塑性应变绝对值均为 `1.387953504904e-04`。

结论：

> FOM-OPT4 单周期候选 PASS。

---

## 6. 10 周期 A/B 验证

10 周期匹配非对称塔筒算例中：

| 指标 | FOM-OPT3 基线 A | FOM-OPT4 候选 B |
| --- | ---: | ---: |
| wall-time | 26.025971 s | 23.401038 s |
| A/B | 1.112172 | |
| 时间降低 | | 10.086% |

完整场仍保持严格一致：

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
- 最终最大损伤均为 `8.634555499717e-03`；
- 最终最大塑性应变绝对值均为 `4.954213909559e-04`。

结论：

> FOM-OPT4 10 周期候选 PASS。

---

## 7. 生产代码修改

生产代码仅修改：

```text
fem/viscoplastic_fiber_section.py
```

实际修改包括三处：

1. 将材料模块导入由 `evaluate_state` 改为 `stress_from_state`；
2. `_integrate_fibers()` 中，以 `stress_from_state()` 直接计算已完成 RK4 更新后的纤维应力；
3. `_response_at_known_state()` 中，以相同方式直接计算已知状态下的纤维应力。

未修改：

- RK4；
- 材料状态更新；
- 损伤截断；
- 截面合力积分；
- 中心差分切线；
- 单元刚度；
- Newton 迭代；
- 状态 commit / rollback。

---

## 8. 回归测试

生产代码落地后依次完成：

| 测试范围 | 结果 |
| --- | --- |
| 黏塑性纤维截面 | 7 passed |
| 黏塑性塔筒系统 | 7 passed |
| 非对称循环塔筒专项 | 17 passed, 10 warnings |
| 完整测试集 | 318 passed, 10 warnings |

完整测试集耗时：

```text
64.98 s
```

所有测试均使用：

```powershell
python -m pytest
```

---

## 9. 生产版本 10 周期直接验证

当前生产源码直接运行，不使用 monkey patch。

得到：

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

单次非 profiling wall-time：

```text
22.721388 s
```

该时间与 10 周期候选 B 的 `23.401038 s` 处于一致量级，说明 FOM-OPT4 的性能收益已经稳定进入生产代码。

但是：

> `22.721388 s` 仍只是单次生产一致性验证，不作为最终 FOM 与 LATIN-PGD 公平效率比较的正式时间。

---

## 10. 科研含义

FOM-OPT4 进一步说明，当前塔筒全阶模型早期较高的运行时间中，存在相当比例的 Python 层重复计算成本。

FOM-OPT4 没有减少任何物理自由度，也没有改变数值求解方法，却仍然获得约 10% 的额外 10 周期运行时间下降。

因此，在最终进行 LATIN-PGD 与 FOM 的效率比较时，必须区分：

- 全阶模型实现成熟度提升；
- LATIN-PGD 实现成熟度提升；
- 真正由 LATIN-PGD 模型降阶带来的方法收益。

不能把旧 FOM 中的实现冗余计入 LATIN-PGD 的算法加速比。

---

## 11. 当前 FOM 优化链

截至 FOM-OPT4：

- FOM-OPT1：标量损伤截断优化；
- FOM-OPT2：RK4 / 材料速率标量化；
- FOM-OPT3：RK4 状态有限性检查标量化；
- FOM-OPT4：纤维截面应力专用求值优化。

这四项均属于纯实现层优化。

从历史 10 周期旧 FOM 的约 `300 s` 量级，到当前单次生产验证约 `22.7 s`，已经可以确认：

> 早期 FOM 的大量耗时来自 Python / NumPy 实现开销，旧 FOM 不能再作为 LATIN-PGD 最终算法效率对照基线。

但正式实现加速比仍应通过统一机器、统一边界、预热和多次重复运行后再报告。

---

## 12. 公平效率比较仍采用三层指标

全阶模型实现优化收益：

$$ S_{\mathrm{FOM,impl}} = \frac{t_{\mathrm{FOM,old}}}{t_{\mathrm{FOM,opt}}} $$

LATIN-PGD 实现优化收益：

$$ S_{\mathrm{LATIN,impl}} = \frac{t_{\mathrm{LATIN,old}}}{t_{\mathrm{LATIN,opt}}} $$

最终真正用于方法比较的公平加速比：

$$ S_{\mathrm{fair}} = S_{\mathrm{method}} = \frac{t_{\mathrm{FOM,opt}}}{t_{\mathrm{LATIN,opt}}} $$

当前仍不应根据单次 `22.721388 s` 与 LATIN-PGD OPT6 的单次时间直接形成最终加速结论。

---

## 13. 当前 checkpoint

FOM-OPT4 生产代码提交：

```text
9357e32  perf: optimize tower FOM stress-only section evaluation
```

当前阶段可以判定：

> FOM-OPT4 已通过单周期 A/B、10 周期 A/B、生产 10 周期数值签名验证以及完整 318 项仓库回归，可以作为新的稳定生产 checkpoint。

---

## 14. 下一步

FOM-OPT4 完成后，下一步仍然应先 profiling，而不是直接修改中心差分数值切线。

建议建立：

```text
tower_fom_10cycle_profile_opt4.py
tower_fom_10cycle_profile_opt4.pstats
```

重新回答以下问题：

- `evaluate_state` 热点是否按预期明显下降；
- RK4 / `_scalar_state_rate_components` 是否仍存在可消除的纯实现层开销；
- 纤维应变插值闭包是否值得优化；
- 截面响应对象构造和状态复制是否开始成为主要成本；
- 中心差分数值切线在去除这些实现开销后是否真正成为不可回避的主导成本。

只有新的 profiling 表明纯实现层热点已经基本耗尽，才应考虑停止 FOM 优化，冻结最终 `FOM,opt`，进入 1 / 2 / 5 / 10 周期正式公平基准。
