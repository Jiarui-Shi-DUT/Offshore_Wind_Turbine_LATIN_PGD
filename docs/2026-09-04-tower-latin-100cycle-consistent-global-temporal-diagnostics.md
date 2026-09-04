# 海上风机塔筒 LATIN-PGD 100 周期一致全时间域时间更新诊断阶段总结

**日期：** 2026-09-04
**项目：** Offshore_Wind_Turbine_LATIN_PGD
**分支：** `perf/tower-fom-optimization`
**上一正式阶段总结：** `docs/2026-09-03-tower-latin-100cycle-long-horizon-temporal-diagnostics.md`
**上一正式 checkpoint：** `5b45f9f3db18e53e23a1a39565b93acbf94427ef`
**当前阶段性质：** 100 周期长时间域时间离散一致性诊断，不涉及生产代码修改
**当前状态：** 已进一步确认 100 周期下原有富集失败与“新模态时间求解”和“全部模态时间重优化”采用不一致时间处理有关；将二者统一为 whole-time global-BE 后，原有 `full_residual_benefit_insufficient` 富集失败消失，但完整求解在 50 次外层迭代后仍未达到 LATIN 收敛阈值，并表现出持续富集、PGD 阶数快速增长和计算成本显著上升。下一阶段应停止盲目增加迭代次数，转而严格回到原论文 Eq. (58)–(72) 与 DG0 时间离散重新推导。

---

## 1. 本阶段从哪里开始

上一阶段已经完成了 100 周期长时间域的第一轮系统诊断，并将结果正式总结在：

`docs/2026-09-03-tower-latin-100cycle-long-horizon-temporal-diagnostics.md`

对应稳定 checkpoint：

`5b45f9f3db18e53e23a1a39565b93acbf94427ef`

上一阶段已经得到以下重要认识：

1. 当前生产版 LATIN-PGD 在 100 周期下不能正常完成收敛。
2. 失败发生在第 5 个 PGD 模态富集阶段。
3. 失败原因是 `full_residual_benefit_insufficient`。
4. 失败模态固定点已经收敛。
5. 空间新颖度、时间显著度、正交性和场不变性均正常。
6. 原生产时间函数更新得到的第 5 个模态会使完整残差明显恶化。
7. 对同一个第 5 个空间模态，如果改用 whole-time global-BE 时间函数，则完整残差可以下降。
8. 但是，如果从求解开始就只把“新模态时间半步”改成 whole-time global-BE，而后续“全部模态时间重优化”仍保留 production sequential 更新，完整求解仍然会失败，而且会更早在 rank 1 附近失败。

因此，上一阶段结束时留下了一个非常明确的问题：

> 如果 whole-time global-BE 能让原生产轨迹中的第 5 个空间模态成为下降方向，那么为什么从求解开始整体替换新模态时间半步以后，完整求解反而在第 2 个模态附近更早失败？

本阶段的全部工作就是围绕这个问题继续展开。

---

## 2. 本阶段研究边界

### 2.1 没有修改生产代码

本阶段所有试验均通过本地诊断脚本进行，没有修改正式生产源码。因此，本阶段得到的是对当前实现机制的诊断证据，而不是已经完成并验证的新生产算法。

### 2.2 没有继续跑 1000 周期

1000 周期仍然暂停。100 周期已经暴露出新的算法问题，在这些问题没有得到理论澄清前，直接推进到 1000 周期既缺乏科学意义，也会产生巨大的时间和内存成本。

### 2.3 没有声称 whole-time global-BE 等价于原论文 DG0

本阶段使用的 whole-time global-BE 只是为了诊断当前时间更新机制而构造的全时间域后向欧拉最小二乘形式。目前只能称为诊断性 formulation，不能称为原论文 Eq. (72) 的严格离散，也不能称为原论文 DG0 的精确实现。原论文一致性仍然需要下一阶段重新推导。

---

## 3. 本阶段统一数值问题

本阶段继续采用与上一阶段完全相同的海上风机塔筒问题。

### 3.1 空间离散

- 梁柱单元：10；
- 每单元 Gauss 点：2；
- 每 Gauss 点环向纤维：16；
- 径向层数：1；
- 总材料点数：320；
- 总自由度：33；
- 自由自由度：30；
- 相同塔筒几何；
- 相同材料参数；
- 相同纤维截面离散。

### 3.2 循环荷载

- 最大塔顶水平力：`+1.0 MN`；
- 最小塔顶水平力：`-0.5 MN`；
- 荷载比：`R = -0.5`；
- 平均荷载：`+0.25 MN`；
- 荷载幅值：`0.75 MN`；
- 周期：`10 s`；
- 每周期：40 个时间增量。

100 周期时：

$$N_t=40\times100+1=4001.$$

### 3.3 当前正式 LATIN-PGD 主要参数

```text
spatial_strategy = "residual_ls"
tolerance = 1.0e-5
fixed_point_tolerance = 1.0e-5
max_fixed_point_iterations = 200
mode_significance_tolerance = 0.0
acceptance_tolerance = 0.0
```

本阶段没有通过放宽容限绕过失败。

---

## 4. 上一阶段 global-BE 完整求解为什么需要继续诊断

上一阶段已经做过一个重要试验：从求解开始，将新模态固定点中的 temporal solve 改成 whole-time global-BE。

该试验结果为：

```text
setup_time_s             ≈ 0.971
solver_time_s            ≈ 23.514
total_time_s             ≈ 24.484

converged                = False
termination_reason       = enrichment_failed
failure_reason           = full_residual_benefit_insufficient

attempted_iterations     = 4
committed_iterations     = 3
trial_evaluations        = 5

accepted_pgd_rank        = 1
accepted_indicator       ≈ 1.7042e-02
total_modes_added        = 1
```

也就是说，把 whole-time global-BE 从一开始用于“新模态时间半步”，并没有直接让完整 100 周期求解成功，而且失败位置从原生产算法中的第 5 个模态提前到了第 2 个模态附近。

这说明，仅仅知道“原生产第 5 空间模态在 global-BE 时间坐标下有下降潜力”还不够，因为一旦改变时间半步，整个 PGD 交替搜索轨迹都会改变。

下一步必须具体检查：新 global-BE 轨迹中的第 2 个模态到底为什么失败。

---

## 5. 第一项新诊断：global-BE 完整轨迹中的第 2 模态为什么失败

本阶段首先构造：

`tower_latin_100cycle_global_be_new_failure_diagnostic.py`

该脚本保持 residual-LS 空间半步、外层 LATIN 逻辑、饱和判据、接受判据、材料模型、网格与荷载不变，只额外捕获新的失败富集，并分析固定点收敛、空间新颖度、时间显著度、正交误差、场不变性以及富集前后残差变化。

---

## 6. global-BE 新轨迹总体结果复现

新诊断成功复现了上一阶段 global-BE full-solver 的早期失败：

```text
setup_time_s                 = 0.993230
solver_time_s                = 23.871006
total_time_s                 = 24.864236
external_elapsed_s           = 24.874924

converged                    = False
termination_reason           = enrichment_failed
failure_reason               = full_residual_benefit_insufficient

attempted_iterations         = 4
committed_iterations         = 3
trial_evaluations            = 5

accepted_pgd_rank            = 1
accepted_indicator           = 1.704194979217e-02
total_modes_added            = 1
```

说明诊断没有改变该轨迹的主要失败行为。

---

## 7. 新失败对应的 Trial-A

最后一个 Trial-A：

```text
indicator   = 0.014213313725739357
saturation  = 0.09050111079085522
converged   = False
```

当前富集阈值：

```text
zeta_enrich = 0.1
```

由于 `0.09050 < 0.1`，算法判断现有 rank 1 基底还不够，需要增加第 2 个 PGD 模态。这里的富集触发逻辑是正常的。

---

## 8. 新第 2 模态固定点已经正常收敛

失败富集：

```text
accepted                               = False
failure_reason                         = full_residual_benefit_insufficient
fixed_point_iterations                 = 6
fixed_point_converged                  = True
```

固定点历史：

```text
[1.28964378e-01
 7.90698734e-03
 1.25922450e-03
 2.01936744e-04
 3.23069613e-05
 6.15191639e-06]
```

最终变化量 `6.1519e-06` 小于 `fixed_point_tolerance = 1.0e-5`。

因此，新 global-BE 轨迹中的第 2 模态失败，也不是因为固定点没有收敛。

---

## 9. 新第 2 模态空间与时间质量都正常

诊断得到：

```text
orthogonal_scale                       = 1.0
spatial_novelty                        = 1.0
temporal_significance                  = 0.6582366910753831
orthogonality_error                    = 2.7304230315292536e-16
plastic_field_invariance_error         = 1.3314004711364718e-18
plastic_rate_field_invariance_error    = 2.036728564984244e-19
stress_field_invariance_error          = 8.301024423194015e-17
```

这说明：

- 空间新颖度为 1；
- 时间显著度约为 0.658；
- 正交误差接近机器精度；
- 场不变性误差也极小。

因此，新第 2 模态不是一个接近已有基底的退化模态，也不是一个几乎没有时间贡献的无意义模态。

---

## 10. 本阶段最关键的新发现之一：第 2 原始模态其实非常有效

本阶段将失败过程拆成三个阶段。

定义：

- A：富集前，只有已经接受的 rank 1 基底；
- B：加入固定点得到的原始第 2 模态以后；
- C：随后执行 production 的全部模态时间函数重优化以后。

得到：

```text
A residual_before                   = 5.103561333469e-01
B residual_raw_appended             = 2.883168625467e-01
C residual_after_temporal_reopt     = 5.642874050011e-01
```

对应：

```text
B vs A raw_mode_benefit             = +4.350673114166e-01
C vs A final_benefit                = -1.056737993929e-01
C vs B reopt_improvement            = -9.571779465716e-01
```

换成百分比：

- A → B：残差下降约 43.51%；
- B → C：残差恶化约 95.72%；
- A → C：最终比富集前还恶化约 10.57%。

这与上一阶段原生产第 5 模态失败形成了完全不同的结构。

---

## 11. 原生产第 5 模态与 global-BE 新第 2 模态的失败机制不同

### 11.1 原生产第 5 模态

上一阶段结果：

```text
A ≈ 0.04777
B ≈ 0.08431
C ≈ 0.05018
```

原始第 5 模态一加入就使残差恶化约 76.5%；后续全部模态时间重优化把大部分恶化修回来，但最终仍比富集前差约 5.05%。因此，原生产第 5 模态的主要问题是：新模态自己的 sequential 时间函数不是完整残差下降方向。

### 11.2 global-BE 新第 2 模态

本阶段结果：

```text
A ≈ 0.51036
B ≈ 0.28832
C ≈ 0.56429
```

新第 2 模态本身非常有效，原始加入后完整残差下降约 43.5%；真正破坏结果的是后续 production 的 all-mode sequential 时间重优化。

因此，新第 2 模态的主要问题是：global-BE 新模态时间半步和 production sequential all-mode 时间重优化之间存在明显不一致。

这是本阶段第一个重要新结论。

---

## 12. 一个非常重要的实现事实：之前所谓“global-BE full solver”其实仍然是混合时间处理

上一阶段那个“global-BE full solver”诊断实际上只替换了新模态固定点中的 `_temporal_solve`。

但是在新模态固定点结束以后，代码仍然调用 `update_tower_pgd_time_functions(...)` 对全部已有模态和新模态的时间函数重新求解，而这一部分仍然是当前 production 的 sequential backward-Euler fixed-basis update。

因此上一阶段的 global-BE full solver 实际上是：

```text
新模态时间半步：
whole-time global-BE

全部模态时间重优化：
production sequential BE
```

也就是说，它并不是一个时间离散完全一致的 whole-time global-BE 求解器。

这一点在本阶段被明确识别出来。

---

## 13. 第二项新诊断：两个模态一起做 whole-time global-BE 时间重优化

为了验证“混合时间处理是否是第 2 模态失败的直接原因”，本阶段进一步构造：

`tower_latin_100cycle_two_mode_joint_global_temporal_diagnostic.py`

目标是比较四种状态：

- A：富集前 rank 1；
- B：加入第 2 模态，保持该模态的 global-BE 时间函数；
- C：production sequential all-mode 时间重优化；
- D：对相同两个空间模态执行 whole-time global-BE 联合时间重优化。

如果 D 真的是在同一个 BE 残差定义下对两个时间坐标做联合最小化，则：

$$D\le B.$$

并且由于 C 也是相同两个空间模态上的另一组可行时间坐标，理论上还应有：

$$D\le C.$$

因此，`D <= B` 和 `D <= C` 可以作为诊断脚本内部非常强的数学自检。

---

## 14. 联合 global-BE 时间最小二乘的离散结构

对固定的多个空间模态，记空间塑性应变矩阵为 $P$，空间应力矩阵为 $S$，时间坐标向量为 $\lambda_n$，时间步长为 $\Delta t_n$，完整 forcing 为 $f_n$。

本阶段诊断采用后向欧拉残差：

$$r_n=P\frac{\lambda_n-\lambda_{n-1}}{\Delta t_n}-H_{\sigma,n}S\lambda_n-f_n.$$

令：

$$A_n=\frac{P}{\Delta t_n}-H_{\sigma,n}S.$$

以及：

$$B_n=-\frac{P}{\Delta t_n}.$$

则：

$$r_n=A_n\lambda_n+B_n\lambda_{n-1}-f_n.$$

对整个时间域求加权最小二乘：

$$\min_{\lambda_1,\ldots,\lambda_{N_t-1}}\sum_{n=1}^{N_t-1}\|r_n\|_{W_n}^2.$$

由于每个时间步只耦合 $\lambda_n$ 与 $\lambda_{n-1}$，所以正规方程形成对称块三对角系统。

对于本阶段 rank 2、100 周期：

```text
unknown temporal coordinates = 4000 × 2 = 8000
normal matrix shape           = (8000, 8000)
```

这在诊断层面仍然可直接求解。

---

## 15. 第一版联合 global-BE 脚本出现了自检失败

第一版运行结果：

```text
A = 5.103561333469e-01
B = 2.883168625467e-01
C = 5.642874050011e-01
D = 7.793562197338e-01
```

并出现：

```text
SELF-CHECK: FAIL  D > B
```

甚至 `D > C`。

这在数学上是不合理的。因为如果 D 真的是在相同空间基底上、对相同完整 BE 残差进行联合全时间域最小化，那么 B 和 C 都是 D 所在可行空间中的候选时间坐标，因此 D 不应该比 B 和 C 更差。

这说明第一版联合 global-BE 诊断脚本自身存在代数或实现错误。这一轮结果不能用于解释 LATIN-PGD 算法。

---

## 16. 第一版联合 global-BE 诊断脚本的具体错误

重新检查正规方程后发现：

残差定义为：

$$r=A\lambda-f.$$

因此最小化：

$$\|A\lambda-f\|_W^2.$$

正规方程应为：

$$A^TWA\lambda=A^TWf.$$

所以 forcing 投影在右端项应为正号。

第一版脚本误写成了负号。对于块时间系统，也就是把原本应为：

```text
rhs = + A^T W f + B^T W f
```

错误写成：

```text
rhs = - A^T W f - B^T W f
```

因此第一版求解的是错误的最小二乘问题。

这是一个诊断脚本错误，不是生产算法错误。

---

## 17. 为什么必须正式记录这个脚本错误

这一点值得在阶段总结中保留，而不是删除。

第一，它说明本阶段坚持了先建立数学自检、再解释数值结果。如果没有 `D <= B` 这类自检，很容易错误地把脚本错误解释成“global-BE 联合时间更新更差”。

第二，它说明当前所有新结论都必须经过相同残差定义下的可行域检查。

第三，这为后续 Eq. (58)–(72) 的严格推导提供了一个重要经验：时间离散中的正负号、forcing 定义、defect 定义和残差定义必须统一，不能仅凭形式相似直接拼接。

---

## 18. 修正后的 V2 联合 global-BE 诊断

修正脚本：

`tower_latin_100cycle_two_mode_joint_global_temporal_diagnostic_v2.py`

主要修正：

```text
rhs = + forcing projection
```

同时增强自检为：

```text
D <= B
and
D <= C
```

修正后运行得到：

```text
A  before enrichment                     = 5.103561333469e-01
B  raw second mode, global-BE temporal    = 2.883168625467e-01
C  production sequential all-mode update  = 5.642874050011e-01
D  JOINT global-BE all-mode update        = 2.789879315282e-01
```

---

## 19. V2 自检正式通过

输出：

```text
SELF-CHECK: PASS  D <= B and D <= C
```

对应：

```text
D = 0.27899
B = 0.28832
C = 0.56429
```

因此：

$$D<B<C.$$

更重要的是：

$$D<A.$$

说明联合 whole-time global-BE 时间重优化不仅没有破坏第 2 模态的收益，而且比原始 B 进一步降低了完整残差。

---

## 20. 四阶段残差的定量比较

修正后：

```text
A = 0.5103561333469
B = 0.2883168625467
C = 0.5642874050011
D = 0.2789879315282
```

对应：

```text
B vs A benefit     = +0.4350673114166
C vs A benefit     = -0.1056737993929
D vs A benefit     = +0.4533465686038
D vs B improvement = +0.0323565223904
D vs C improvement = +0.5055924887644
```

即：

- B 相对 A 降低约 43.51%；
- C 相对 A 恶化约 10.57%；
- D 相对 A 降低约 45.33%；
- D 相对 B 进一步改善约 3.24%；
- D 相对 C 改善约 50.56%。

这是本阶段第二个核心新发现。

---

## 21. 由 V2 可以支持的直接结论

对于这一条具体 global-BE 轨迹中的第 2 模态，可以直接说：

> 第 2 个空间模态本身是有效的。

可以直接说：

> whole-time global-BE 得到的新模态时间函数是完整 BE 残差的下降方向。

可以直接说：

> 对相同两个空间模态进行一致的 whole-time global-BE 联合时间重优化，可以保持并进一步改善该下降效果。

可以直接说：

> 之前第 2 模态的拒绝不是因为第 2 空间模态无效，而是因为后续切回 production sequential all-mode 时间重优化后破坏了已有下降。

---

## 22. 目前不能从 V2 直接推出什么

即使 V2 通过，也仍不能直接说：

- global-BE 就是原论文正确的 DG0；
- 原论文 Eq. (72) 必须整体时间耦合求解；
- production sequential 时间更新在所有情况下都是错误的；
- LATIN-PGD 在 100 周期下已经解决。

V2 只证明：在当前构造的 whole-time BE 残差定义下，新模态与 all-mode 时间坐标如果采用一致的全时间域最小二乘处理，可以避免前面观察到的第 2 模态人为破坏。

---

## 23. 第三项新诊断：构造时间处理完全一致的 100 周期 full solver

在 V2 之后，问题自然变为：如果我们不再混合两套时间更新，而是让 100 周期完整求解中所有关键 PGD 时间更新都采用同一个 whole-time global-BE 定义，原来的富集失败是否会消失？

因此构造：

`tower_latin_100cycle_consistent_global_be_full_solver_diagnostic.py`

该脚本在当前 Python 进程中临时替换三个位置：

1. 新模态固定点中的 temporal solve；
2. Trial-A 固定基时间函数更新；
3. 富集完成后的 all-mode 时间函数重优化。

三者全部统一使用 whole-time global-BE temporal treatment。

与此同时保持 residual-LS 空间半步、外层 LATIN 事务式 Trial-A / Trial-B、饱和逻辑、接受逻辑、材料模型、搜索方向、网格、荷载以及所有正式容限。

这仍然是诊断，不修改生产代码。

---

## 24. 一致 global-BE full solver 的最重要结果

100 周期运行：

```text
setup_time_s                 = 1.017535
solver_time_s                = 841.349486
total_time_s                 = 842.367020
external_elapsed_s           = 842.376921

converged                    = False
termination_reason           = max_iterations
failure_reason               = None

attempted_iterations         = 50
committed_iterations         = 50
trial_evaluations            = 96

accepted_pgd_rank            = 46
accepted_indicator           = 1.082333446290e-03
total_modes_added            = 46
```

这是本阶段最重要的 full-solver 结果。

---

## 25. 原来的富集失败机制已经消失

原生产算法 100 周期：

```text
termination_reason = enrichment_failed
failure_reason     = full_residual_benefit_insufficient
```

只成功接受 4 个模态。

而一致 global-BE：

```text
termination_reason = max_iterations
failure_reason     = None
```

并且成功接受 46 个模态。

在 50 次外层迭代内，没有再次出现 `full_residual_benefit_insufficient`。

因此可以明确说：

> 将新模态时间求解、Trial-A 固定基时间更新和 all-mode 时间重优化统一为相同 whole-time global-BE 定义后，原先 100 周期下的受控富集拒绝机制消失了。

这是本阶段第三个核心发现。

---

## 26. 为什么这个结果比“单个模态改善”更强

上一阶段只证明：对原生产轨迹中的第 5 个固定空间模态，global-BE 时间函数可以让完整残差下降。这仍然是一个局部、冻结空间方向的结果。

本阶段现在进一步证明：在一个完整 100 周期 LATIN-PGD 求解过程中，将关键时间更新保持一致以后，算法能够连续通过 46 次 PGD 富集，而不再触发原来的 full residual acceptance failure。

这说明：时间更新不一致不是一个偶然发生在单个模态上的数值细节，而是会真正影响完整长时间域 PGD 富集轨迹的鲁棒性。

---

## 27. 但一致 global-BE 仍然没有让 100 周期正式收敛

必须同样强调：

```text
converged = False
```

原因不是富集失败，而是：

```text
termination_reason = max_iterations
```

当前最多允许 50 个外层 LATIN 迭代。

在第 50 次正式提交以后：

```text
accepted_indicator = 1.082333446290e-03
```

而正式目标：

```text
tolerance = 1.0e-5
```

因此仍相差约两个数量级。

所以：一致 global-BE 解决的是当前观察到的长时间域富集鲁棒性问题，但没有自动解决完整 LATIN 收敛效率问题。

---

## 28. LATIN 指标在后期呈稳定下降，而不是振荡或爆炸

最后 15 个 accepted indicator：

```text
[0.00153539
 0.00149449
 0.00145523
 0.00141754
 0.00138133
 0.00134651
 0.00131298
 0.00128065
 0.00124946
 0.00121933
 0.00119015
 0.00116193
 0.00113459
 0.00110807
 0.00108233]
```

这一段具有明显特点：持续下降，没有明显大幅回升，没有爆炸，也没有停滞到机器精度，但下降速度较慢。

因此目前更准确的描述是：求解轨迹是稳定的，但收敛效率偏低，而不是算法数值不稳定。

---

## 29. 饱和指标揭示了为什么 PGD 阶数一直增长

最后 15 个饱和指标：

```text
[0.01246359
 0.01232978
 0.01232846
 0.01220712
 0.01220378
 0.01210020
 0.01210872
 0.01200082
 0.01200094
 0.01190436
 0.01190434
 0.01182237
 0.01182252
 0.01175095
 0.01175118]
```

当前富集阈值：

```text
zeta_enrich = 0.1
```

而后期：

```text
zeta ≈ 0.012
```

明显低于 `0.1`。

因此当前饱和控制持续判断：更新已有时间函数带来的改善不足，需要继续增加新的 PGD 模态。

这与最后的 commit history 完全一致。

---

## 30. 后期几乎每一步都在提交 Trial-B

最后 15 个 commit：

```text
('B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B')
```

这意味着后期几乎每一个 LATIN 外层迭代都依赖新增一个 PGD 模态继续前进。

因此 100 周期一致 global-BE 轨迹表现出持续富集型收敛，而不是固定低阶基底上的快速时间函数修正收敛。

---

## 31. PGD 阶数从 1–10 周期到 100 周期的变化

正式 1、2、5、10 周期公平 benchmark 中：

```text
1 cycle   rank = 11
2 cycles  rank = 13
5 cycles  rank = 17
10 cycles rank = 21
```

100 周期一致 global-BE 在未收敛时已经达到：

```text
rank = 46
```

但必须严格注意：不能说“100 周期最终需要 rank 46”。因为 `xi ≈ 1.08e-3` 仍远高于 `1e-5`。

正确表述是：

> 在 100 周期一致 global-BE 诊断中，当算法运行到第 50 次外层迭代、接受 46 个 PGD 模态时，LATIN 指标仍为约 `1.08e-3`，因此最终收敛所需阶数尚未知。

---

## 32. 当前对低秩性的认识必须更加谨慎

早期 1–10 周期结果证明塔筒问题存在明显的低秩表示能力。

但是 100 周期结果提醒我们：长时间域下“存在低秩结构”和“当前算法能够高效找到低秩结构”是两个不同问题。

100 周期 rank 持续增长可能来自多种原因：

1. 真实空间-时间响应需要更多模态；
2. 当前时间离散不够贴合原理论；
3. residual-LS 空间半步和 global temporal half-step 组合不是理想的交替极小；
4. 饱和指标在长时间域下触发过于积极；
5. whole-history 残差范数随时间域增长改变了各模态的相对贡献；
6. 当前时间更新虽然一致，但不一定是原论文最优的 PGD 时间更新结构。

目前不能把这些原因区分开。

---

## 33. 一致 global-BE 最后一个富集模态仍然是正常下降方向

最后一次富集：

```text
accepted                     = True
failure_reason               = None

fixed_point_iterations       = 23
fixed_point_converged        = True

spatial_novelty              = 1.000000000000e+00
temporal_significance        = 2.711231285761e-04

residual_norm_before         = 2.882045646621e-04
residual_norm_after          = 2.757497574278e-04
residual_benefit             = 4.321516298274e-02
```

即最后一个新模态仍使完整残差下降约 4.32%。

因此，第 50 次停止不是因为新模态无效、固定点失效、残差下降失败、空间线性相关或时间求解崩溃，而只是因为达到 `max_iterations = 50`。

---

## 34. 最后一个固定点仍然正常收敛

最后一个新模态固定点：

```text
[3.76152107e-01
 1.36986965e-01
 7.33748821e-02
 5.10242728e-02
 3.97789820e-02
 3.24861357e-02
 2.69281369e-02
 2.23654144e-02
 1.85246621e-02
 1.52789383e-02
 1.25469681e-02
 1.02621593e-02
 8.36405839e-03
 2.71362573e-02
 6.35943918e-03
 1.72975292e-03
 8.45229866e-04
 1.50776208e-04
 4.45043159e-05
 8.70261702e-05
 6.41609126e-05
 1.46942428e-05
 5.83811860e-06]
```

最终 `5.838e-06 < 1e-5`，所以固定点仍然通过。

虽然中间有局部回升，但最终仍然收敛。因此，当前主要矛盾已经从“单个新模态固定点是否能收敛”转移到“需要多少次持续富集才能让外层 LATIN 指标达到正式容限”。

---

## 35. 计算时间暴露出新的工程问题

一致 global-BE 100 周期：

```text
solver_time_s = 841.349486
```

约为 14.0 min，而且 `converged = False`。

这意味着目前没有理由直接把最大迭代数从 50 粗暴提高到 100、200 或更高。

即使继续运行，也可能只是消耗几十分钟到数小时，让 PGD 阶数继续增长，却仍然无法解释为什么长时间域下富集如此频繁，也无法确认这种 global-BE 时间更新是否忠实于原论文。

所以继续长跑的科研信息增益已经明显低于理论回溯。

---

## 36. 目前不能把 841 s 当作正式效率 benchmark

必须强调，`841.349486 s` 不是一个可与 FOM 正式比较的 LATIN 时间。

原因包括：

1. 当前求解没有收敛；
2. 当前使用的是诊断性 whole-time global-BE；
3. 当前并非生产实现；
4. 当前没有公平重复次数；
5. 当前算法每次 all-mode 时间更新需要求解更大的全时间域系统；
6. 当前运行只用于机制诊断。

因此，不能用 841 s 计算 FOM/LATIN 公平效率比，也不能据此说 LATIN 比 FOM 慢多少倍。

---

## 37. 当前阶段已经得到的三层证据

### 37.1 第一层：单模态固定空间方向证据

上一阶段已经证明：原生产第 5 空间模态在 sequential temporal solve 下不是下降方向，但在 whole-time global-BE 时间坐标下可以成为下降方向。

### 37.2 第二层：相同两个空间模态上的时间重优化证据

本阶段 V2 证明：global-BE 第 2 模态加入后残差下降，且对两个模态进行一致的 joint global-BE all-mode 时间重优化可以进一步下降。

### 37.3 第三层：完整求解器轨迹证据

本阶段 consistent full solver 证明：将所有关键 PGD 时间更新统一为 whole-time global-BE 后，原生产版在 100 周期下出现的 `full_residual_benefit_insufficient` 富集失败在前 50 次外层迭代中消失。

三层证据相互一致。

---

## 38. 当前可以认为已经较强支持的机制判断

当前可以较强支持：

> 100 周期下，当前生产实现中的 sequential temporal treatment 与完整长时间域残差下降之间存在明显局限。

当前可以较强支持：

> 只把新模态时间半步改成 whole-time global-BE，但 all-mode 时间重优化仍使用 sequential 更新，会产生内部不一致。

当前可以较强支持：

> 这种时间处理不一致可以直接导致一个原本有效的新空间模态在后续 all-mode 时间重优化以后被破坏，并最终触发 `full_residual_benefit_insufficient`。

当前可以较强支持：

> 当关键时间更新统一为 whole-time global-BE 后，原观察到的富集失败机制会消失。

---

## 39. 当前仍然没有证明的核心问题

### 39.1 whole-time global-BE 是否等价于原论文 Eq. (72)

没有证明。

### 39.2 原论文 DG0 是否应形成 whole-time 耦合系统

没有证明。

### 39.3 production sequential 更新是否理论上错误

没有证明。它可能只是当前工程离散选择，在短时间域足够有效，但在长时间域逐渐偏离完整残差最小化。现在不能直接称其“理论错误”。

### 39.4 100 周期最终需要多少 PGD 模态

没有证明。

### 39.5 100 周期 consistent global-BE 最终是否一定能收敛到 `1e-5`

没有证明。

### 39.6 PGD 阶数持续增长的根本原因

没有证明。

### 39.7 当前饱和阈值 `0.1` 是否适用于 100 周期

没有证明。

---

## 40. 为什么下一步不应该直接调饱和阈值

看到 `zeta ≈ 0.012` 而 `zeta_enrich = 0.1`，一个直接想法是把富集阈值调小，是不是就能少加模态。

当前不应该马上这样做。原因是饱和判据的数学意义必须和原论文算法结构一起理解。

如果不了解 Eq. (58)–(72) 和原论文饱和策略之间的关系，仅凭 100 周期结果调阈值，很可能只是用参数补偿 formulation 问题。这种做法不适合作为下一阶段主线。

---

## 41. 为什么下一步也不应该直接把 max_iterations 调大

把 50 改成 100 甚至 200，确实可能让 indicator 进一步下降。

但这种运行只能回答“继续跑会不会继续下降”，却不能回答“为什么要这么多模态”，也不能回答“当前 whole-time global-BE 是否忠实于原论文”，更不能回答“当前时间离散到底应该是 sequential 还是 global coupled”。

所以从科研效率看，理论推导的优先级已经高于继续长算例。

---

## 42. 下一阶段必须重新回到原论文的核心公式

下一阶段建议围绕：

- Eq. (58)
- Eq. (59)
- Eq. (60)
- Eq. (70)
- Eq. (71)
- Eq. (72)

以及 DG0 时间离散重新完整推导。

核心目标不是再提出一个新算法，而是回答：原论文到底如何从连续时间变分形式得到时间函数更新方程？

特别需要确认：

1. 时间试函数是什么；
2. 时间离散单元是什么；
3. DG0 中时间函数在哪个时间区间内为常数；
4. 时间导数如何处理；
5. 跳跃项如何出现；
6. 相邻时间 slab 如何耦合；
7. Eq. (72) 是局部顺序递推，还是全时间域代数系统的一部分；
8. 原文在 implementation 层面是否允许逐步求解；
9. 当前 production sequential BE 与论文形式哪里一致；
10. 当前 whole-time global-BE 与论文形式哪里一致、哪里不一致。

---

## 43. 下一阶段理论推导的核心判定问题

下一阶段最核心的问题可以浓缩成一句话：

> 原论文 Eq. (72) 在 DG0 时间离散以后，时间函数更新究竟应当表现为逐时间步顺序求解，还是表现为全时间域耦合的离散系统？

这是现在最值得优先回答的问题。

因为它直接关系到：

- 100 周期富集鲁棒性；
- 时间方向残差下降；
- PGD 阶数增长；
- 长时间域效率；
- 原论文一致性；
- 是否应修改 production 时间更新。

---

## 44. 需要特别比较的三套时间处理

下一阶段理论推导时，应把三套东西并排比较。

### 44.1 原论文理论时间更新

来源：Bhattacharyya 等人的 LATIN-PGD 循环损伤原论文。必须严格依据原文。

### 44.2 当前 production sequential 时间更新

当前工程实现位于 `latin/tower_pgd_time_update.py`，其特点是固定空间基底以后，在每个时间步顺序求当前时间坐标，并使用上一步时间坐标作为已知量。

### 44.3 本阶段 diagnostic whole-time global-BE

特点是把所有时间坐标作为联合未知量，对同一个后向欧拉完整时间残差做全时间域加权最小二乘。

三者不能混为一谈。

---

## 45. 当前不同残差范数定义也值得下一阶段重新核对

当前实现中存在一个容易被忽略的细节。

新模态接受中的 `_be_residual_norm` 使用 slab 右端点的 BE 残差和 `dt` 加权。

而 `update_tower_pgd_time_functions` 中的报告范数使用 trapezoidal time integration。

本阶段 whole-time global-BE 诊断为了与 enrichment acceptance 保持一致，主要针对 `_be_residual_norm` 构造全时间域最小二乘。

这意味着下一阶段需要检查：原论文中 PGD 时间更新、饱和指标和富集接受是否应该基于完全相同的时间积分定义。

目前不能断言这里就是问题，但它已经成为值得正式核对的 formulation consistency 项。

---

## 46. 本阶段关于“时间一致性”的最通俗理解

可以把 PGD 的空间模态理解为“往哪个方向修正”，把时间函数理解为“每个时刻沿这个方向修正多少”。

原生产第 5 模态失败时，空间方向本身不一定坏，但 sequential 时间分配使其在完整 100 周期历史上成为非下降方向。

后续 global-BE 第 2 模态诊断又说明，即使给新模态安排了一套好的全时间域时间分配，如果随后再用另一套 sequential 规则把所有模态时间函数重新安排一次，也可能把刚刚得到的好结果破坏掉。

因此，本阶段最直观的认识是：长时间域下，PGD 时间函数求解不仅要“能算出来”，还要在新模态求解、已有基底更新和富集后重优化之间保持同一个残差定义和同一种离散逻辑。

---

## 47. 本阶段与原生产第 5 模态失败的完整因果链

```text
原生产 100 周期
        ↓
第 5 模态 fixed point 收敛
        ↓
原始第 5 模态使 full residual 恶化
        ↓
富集被拒绝
        ↓
同一空间模态改 global-BE 时间函数
        ↓
full residual 可以下降
        ↓
说明该空间方向存在下降潜力
        ↓
从求解开始替换新模态 temporal solve
        ↓
第 2 模态附近再次失败
        ↓
进一步拆分 A/B/C
        ↓
发现第 2 模态本身使 residual 降低 43.5%
        ↓
production all-mode sequential update 把结果破坏
        ↓
joint global-BE all-mode V2
        ↓
residual 进一步降低
        ↓
构造 consistent global-BE full solver
        ↓
原 full_residual_benefit_insufficient 消失
        ↓
算法连续富集到 rank 46
        ↓
但 50 次外层迭代仍未达到 xi = 1e-5
        ↓
新的主问题变成：
长时间域低秩效率与原论文时间离散一致性
```

---

## 48. 本阶段哪些属于实测事实

以下均为本阶段实际运行得到的数值事实。

### 48.1 global-BE 新第 2 模态

```text
fixed_point_converged = True
fixed_point_iterations = 6
spatial_novelty = 1.0
temporal_significance ≈ 0.658
```

### 48.2 A/B/C

```text
A = 0.5103561333469
B = 0.2883168625467
C = 0.5642874050011
```

### 48.3 V2 joint global-BE

```text
D = 0.2789879315282
```

且 `D <= B`、`D <= C`。

### 48.4 consistent global-BE full solver

```text
50 attempted iterations
50 committed iterations
46 accepted PGD modes
xi ≈ 1.0823e-3
termination = max_iterations
failure_reason = None
```

### 48.5 最后一次富集

```text
accepted = True
fixed_point_converged = True
residual_benefit ≈ +4.32%
```

这些属于已测事实。

---

## 49. 本阶段哪些属于基于事实的机制解释

以下是当前证据强烈支持的解释，但仍属于机制判断。

### 49.1 原 global-BE 第 2 模态失败来自时间处理混用

依据：B 明显优于 A，C 明显劣于 B，D 又明显优于 C，而且 consistent full solver 不再出现相同富集失败。

因此可以较强判断：混合 global-BE new-mode update 与 sequential all-mode update 是该新失败的直接原因。

### 49.2 production sequential 时间更新在长时间域上存在局限

依据：原生产第 5 模态 sequential 时间函数非下降，同空间方向 global-BE 可下降，all-mode sequential update 又会破坏 global-BE 第 2 模态。

因此可以较强判断：当前 sequential 时间处理至少在 100 周期长时间域中不能稳定保证完整残差下降。

但这仍然不等价于“它违反原论文”。

---

## 50. 本阶段哪些结论仍然不能成立

以下说法目前不能写进论文结论或对导师作为确定事实表述：

- “原论文 Eq. (72) 应使用 whole-time global-BE”；
- “当前 production sequential BE 写错了”；
- “100 周期需要 46 个 PGD 模态”；
- “consistent global-BE 一定会最终收敛”；
- “100 周期低秩性消失了”；
- “LATIN-PGD 不适合长时间域疲劳”。

当前只能说：当前塔筒 implementation 在 100 周期下暴露出时间离散一致性和持续富集效率问题。

---

## 51. 与 1–10 周期结果的关系

前面 1、2、5、10 周期正式公平 benchmark 仍然有效。它们属于当前正式 production implementation 在短到中等时间域下的实测效率结果。

本阶段没有推翻这些结果。

新的认识是：当时间域扩展到 100 周期后，之前在 1–10 周期中不明显的时间方向 formulation 问题开始暴露。

这意味着短时间域验证通过不能自动推出长时间域算法结构一定合理。

---

## 52. 与一维三材料杆验证的关系

此前一维三材料杆上的 LATIN-PGD 已经成功复现并验证。本阶段不应把塔筒 100 周期问题解释为“一维验证无效”。

更合理的理解是：

- 一维杆验证了基本 LATIN-PGD 实现链条；
- 塔筒 1–10 周期验证了空间扩展和短时间域低秩表示；
- 塔筒 100 周期进一步暴露长时间域时间函数更新问题。

这体现的是数值问题复杂度逐层提高以后，新的 formulation 限制逐步显现。

---

## 53. 对当前研究路线的影响

本阶段结果帮助我们避免了两条风险较大的路线。

第一条错误路线是看到 100 周期失败就直接放宽 acceptance tolerance，这会掩盖真正的时间方向问题。

第二条错误路线是看到 global-BE 单模态有效就立即把整个生产代码替换成 global-BE。本阶段 consistent full solver 说明，即使鲁棒性问题消失，仍然存在持续富集和效率问题。

所以现在最合理的路线是：先回到原论文严格推导，再决定 production 时间更新应该如何修改。

---

## 54. 当前推荐的下一阶段工作顺序

下一阶段建议严格按以下顺序推进。

### 第一优先级

重新阅读原论文 Eq. (58)–(72) 周围完整上下文。

### 第二优先级

从连续 LATIN-PGD 弱式或变分形式重新推导时间函数方程。

### 第三优先级

严格推导 DG0 时间离散。

### 第四优先级

把原论文离散形式与 production sequential BE、diagnostic whole-time global-BE 逐项对照。

### 第五优先级

只有在理论关系澄清以后，才决定是否修改生产时间更新、是否需要新的一致 all-mode update、是否重新评估 saturation criterion、是否继续 100 周期正式收敛测试。

### 第六优先级

100 周期稳定后，再设计 1000 周期内存安全 benchmark。

---

## 55. 下一阶段暂时不做的事情

暂时不：

- 跑 1000 周期；
- 把 `max_iterations` 直接调到 100 或 200；
- 调低 `zeta_enrich`；
- 放宽 `acceptance_tolerance`；
- 把 diagnostic global-BE 合入 production；
- 做新的 FOM/LATIN 公平效率 benchmark；
- 对 841 s 做正式效率解读；
- 根据 rank 46 外推最终阶数。

---

## 56. 本阶段新增本地诊断脚本

本阶段新增的重要本地诊断脚本包括：

```text
tower_latin_100cycle_global_be_new_failure_diagnostic.py
tower_latin_100cycle_two_mode_joint_global_temporal_diagnostic.py
tower_latin_100cycle_two_mode_joint_global_temporal_diagnostic_v2.py
tower_latin_100cycle_consistent_global_be_full_solver_diagnostic.py
```

其中：

### `tower_latin_100cycle_global_be_new_failure_diagnostic.py`

作用：诊断 global-BE full-solver 新轨迹中第 2 模态失败。

### `tower_latin_100cycle_two_mode_joint_global_temporal_diagnostic.py`

作用：第一版两个模态 joint global-BE 时间重优化。

注意：该版本正规方程 forcing RHS 符号有误，结果不可用于算法解释。

### `tower_latin_100cycle_two_mode_joint_global_temporal_diagnostic_v2.py`

作用：修正 forcing RHS 符号并通过 `D <= B`、`D <= C` 自检。这是当前有效的 joint global-BE 诊断版本。

### `tower_latin_100cycle_consistent_global_be_full_solver_diagnostic.py`

作用：在完整 100 周期求解轨迹中，使新模态时间半步、Trial-A 固定基时间更新和 all-mode 时间重优化统一使用 whole-time global-BE。

所有这些脚本当前仍应保持本地 untracked diagnostic，不要直接提交为生产代码。

---

## 57. 本阶段推荐形成的正式 checkpoint

建议本阶段只提交本阶段总结 Markdown。

不提交：

- 临时诊断脚本；
- profiling 文件；
- CSV/JSON 临时输出；
- patch；
- build 脚本。

建议文档路径：

`docs/2026-09-04-tower-latin-100cycle-consistent-global-temporal-diagnostics.md`

建议提交信息：

```text
docs: summarize consistent global temporal diagnostics
```

---

## 58. 本阶段核心结论摘要

如果把本阶段压缩成最重要的六条，可以概括为：

1. global-BE 新轨迹中的第 2 个空间模态本身是有效的，原始加入后完整残差下降约 43.5%。
2. production sequential all-mode 时间重优化会把这一有效结果破坏，使残差最终比富集前恶化约 10.6%。
3. 对相同两个空间模态采用 joint whole-time global-BE 时间重优化后，残差进一步下降到约 0.279，自检 `D <= B` 和 `D <= C` 通过。
4. 因此，之前 global-BE full-solver 的第 2 模态失败直接暴露了“global new-mode temporal solve + sequential all-mode reoptimization”的时间处理不一致。
5. 将关键 PGD 时间更新统一为 whole-time global-BE 后，原 100 周期 `full_residual_benefit_insufficient` 富集失败在前 50 次外层迭代中消失，算法连续接受 46 个模态。
6. 但是求解仍未达到 `1e-5` LATIN 容限，且阶数持续增长、计算成本显著提高，因此下一阶段应停止盲目长跑，严格回到原论文 Eq. (58)–(72) 与 DG0 时间离散重新推导。

---

## 59. 当前最准确的科研判断

截至本阶段，最准确、最克制的表述是：

> 当前塔筒 LATIN-PGD implementation 在 100 周期长时间域下暴露出明显的时间更新一致性问题。原生产实现中的 sequential 时间函数处理不能稳定保证新增模态在完整时间历史上的残差下降；只将新模态时间半步替换为 whole-time global-BE 又会与后续 sequential all-mode 时间重优化形成内部不一致。通过将关键时间更新统一为同一个 whole-time global-BE 残差定义，可以消除此前观察到的 `full_residual_benefit_insufficient` 富集失败，并使完整求解持续稳定富集。然而，该诊断版本在 50 次外层迭代后仍未达到正式 LATIN 收敛阈值，且 PGD 阶数增长到 46，说明长时间域下仍存在低秩效率、饱和控制或时间离散 formulation 方面尚未澄清的问题。由于 whole-time global-BE 目前只是诊断性工程离散，不能视为原论文 DG0 的等价实现，下一步必须回到 Eq. (58)–(72) 与 DG0 原始推导，重新建立 paper-fidelity 基准以后再决定生产算法修改方向。

---

## 60. 下一阶段唯一优先问题

下一阶段首先回答：

> 原论文 Eq. (72) 在 DG0 时间离散下到底应该如何离散、组装和求解？

并由此继续回答：

> 当前 production sequential 时间更新、当前 diagnostic whole-time global-BE 与原论文之间分别是什么关系？

在这个问题解决前：

> 不继续进行 100 周期长跑，不推进 1000 周期，不修改生产核心代码。
