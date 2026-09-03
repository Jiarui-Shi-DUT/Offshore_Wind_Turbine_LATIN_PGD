# 海上风机塔筒 LATIN-PGD 100 周期长时间域失败与时间函数诊断阶段总结

**日期：** 2026-09-03
**项目：** Offshore_Wind_Turbine_LATIN_PGD
**分支：** `perf/tower-fom-optimization`
**上一正式阶段总结：** `docs/2026-09-03-tower-fom-latin-optimized-fair-efficiency-benchmark.md`
**上一正式 checkpoint：** `0b061cb6bc301eaa0e510120215e52ab018e6f9f`
**当前阶段性质：** 长时间域鲁棒性诊断，不涉及生产代码修改
**当前状态：** 100 周期 LATIN-PGD 首次出现受控富集失败；已完成多轮针对性诊断并基本定位到“长时间域下新模态时间函数求解与完整残差下降之间的不一致”，但正式算法修正尚未完成

---

## 1. 本阶段为什么发生

上一阶段已经完成了 FOM 与 LATIN-PGD 的实现级优化，并在相同海上风机塔筒问题上完成了 1、2、5、10 周期的正式公平效率比较。

上一阶段得到的核心效率结果为：

| 周期数 | FOM,opt 中位时间 / s | LATIN-PGD,opt 总时间中位数 / s | 公平效率比 $S_{\mathrm{fair}}$ |
|---:|---:|---:|---:|
| 1 | 2.739699 | 1.454202 | 1.883987 |
| 2 | 5.071664 | 2.876301 | 1.763259 |
| 5 | 11.738477 | 9.852742 | 1.191392 |
| 10 | 22.851421 | 24.985172 | 0.914599 |

其中：

$$S_{\mathrm{fair}}=\frac{\operatorname{median}(t_{\mathrm{FOM,opt}})}{\operatorname{median}(t_{\mathrm{LATIN,opt,total}})}.$$

这说明：

- 1、2、5 周期下 LATIN-PGD 更快；
- 10 周期下 FOM 略快；
- 当前效率交叉区间位于 5 到 10 周期之间；
- PGD 阶数仍然很低，但这种低秩性并没有在 10 周期时完全转化为 wall-time 优势。

上一阶段因此留下了两个重要问题：

1. 当时间域进一步增长到 100 周期、1000 周期时，这种趋势是否会进一步放大；
2. 当前 LATIN-PGD 是否仍然能够保持稳定收敛，还是会出现新的长时间域算法问题。

基于这个问题，本阶段没有立即继续做 1 到 10 周期 profiling，而是先增加一个更长的时间尺度锚点：

> 100 周期。

---

## 2. 本阶段统一数值问题保持不变

本阶段所有诊断均基于当前稳定生产实现及相同塔筒问题。

### 2.1 空间离散

- 梁柱单元：10；
- 每单元 Gauss 点：2；
- 每 Gauss 点环向纤维：16；
- 径向层数：1；
- 总材料点数：320；
- 自由度：33，总自由自由度为 30；
- 相同塔筒几何；
- 相同材料参数；
- 相同纤维截面离散。

### 2.2 循环荷载

- 最大塔顶水平力：`+1.0 MN`；
- 最小塔顶水平力：`-0.5 MN`；
- 荷载比：`R = -0.5`；
- 平均荷载：`+0.25 MN`；
- 荷载幅值：`0.75 MN`；
- 周期：`10 s`；
- 每周期：40 个时间增量。

对于 $n$ 个周期：

$$N_t=40n+1.$$

因此 100 周期时：

$$N_t=4001.$$

### 2.3 LATIN-PGD 参数

本阶段继续采用正式 benchmark 中的稳定参数：

```text
spatial_strategy = "residual_ls"
tolerance = 1.0e-5
fixed_point_tolerance = 1.0e-5
max_fixed_point_iterations = 200
mode_significance_tolerance = 0.0
acceptance_tolerance = 0.0
```

特别强调：

- 本阶段没有通过修改容限强行使 100 周期通过；
- 没有把 `acceptance_tolerance` 改成负数；
- 没有降低 LATIN 收敛要求；
- 没有修改生产代码；
- 所有额外试验均为本地诊断脚本。

---

## 3. 第一项长时间域试验：100 周期公平 benchmark 试跑

执行：

```powershell
python tower_fom_latin_fair_benchmark.py --cycles 100 --repeats 1 --output-prefix tower_fair_100cycle_pilot
```

环境：

```text
timestamp_local          = 2026-09-03 19:20:45
git_branch               = perf/tower-fom-optimization
git_commit               = 0b061cb6bc301eaa0e510120215e52ab018e6f9f
python                   = 3.8.20
numpy                    = 1.23.5
platform                 = Windows-10-10.0.26200-SP0
logical_cpu_count        = 20
```

正式 100 周期试跑前完成：

```text
1-cycle FOM warm-up
1-cycle LATIN-PGD warm-up
```

随后：

```text
100-cycle FOM
100-cycle LATIN-PGD
```

### 3.1 关键结果

100 周期 LATIN-PGD 没有完成正常收敛，而是在 PGD 富集阶段终止：

```text
RuntimeError:
LATIN did not converge:
enrichment_failed / full_residual_benefit_insufficient
```

这意味着：

> 100 周期问题并不是简单地“计算时间变长”，而是首次暴露了一个长时间域 PGD 富集鲁棒性问题。

### 3.2 这一失败不是崩溃，而是受控失败

当前外层 LATIN-PGD 求解器采用事务式 Trial-A / Trial-B 结构。

当某次新模态富集失败时：

- 已接受的持久状态不被破坏；
- 失败模态不会被强行加入 PGD 基；
- 求解器返回最后一个有效 persistent snapshot；
- `termination_reason = enrichment_failed`；
- 同时记录具体 `failure_reason`。

因此：

> 本次 100 周期失败是算法保护机制主动拒绝一个无效候选模态，而不是程序异常崩溃。

这一点非常重要。

---

## 4. 第一轮失败诊断：100 周期到底在哪一步失败

使用专门诊断脚本重新运行 100 周期 LATIN-PGD，并打印失败瞬间的完整求解状态。

### 4.1 总体结果

```text
setup_time_s             = 1.000707
solver_time_s            = 102.211783
total_time_s             = 103.212489

converged                = False
termination_reason       = enrichment_failed
failure_reason           = full_residual_benefit_insufficient

attempted_iterations     = 10
committed_iterations     = 9
trial_evaluations        = 14

accepted_pgd_rank        = 4
accepted_indicator       = 6.554750722475e-03
total_modes_added        = 4
```

这说明：

- 前 9 次持久状态更新均已成功；
- 第 10 次 LATIN 外层迭代打开；
- 当前已经成功接受 4 个 PGD 模态；
- 当前接受状态的 LATIN 指标约为 `6.55e-3`；
- 第 10 次迭代要求继续富集，但第 5 个模态被拒绝。

### 4.2 已接受 LATIN 指标历史

```text
[0.06150361
 0.01765225
 0.01397141
 0.01981406
 0.01424584
 0.01098389
 0.00951909
 0.00769659
 0.00655475]
```

注意：

> LATIN 指标并不是严格单调下降。

例如：

```text
0.01397141 → 0.01981406
```

曾经出现一次回升。

这与当前事务式 LATIN 外层控制一致：

- 并不是每次 Trial 都要求比上一次更小；
- 有效 Trial-B 可以原子提交；
- 真正的最终收敛由 LATIN 指标与饱和控制共同决定。

### 4.3 第 10 次 Trial-A

第 10 次外层迭代中：

```text
baseline xi = 0.00655475
Trial-A xi  = 0.00579630
```

对应饱和指标：

$$\zeta=\frac{\xi_{\mathrm{previous}}-\xi_{\mathrm{current}}}{\xi_{\mathrm{previous}}+\xi_{\mathrm{current}}}.$$

得到：

```text
zeta = 0.06140817
```

当前富集阈值为：

```text
zeta_enrich = 0.1
```

因此：

```text
0.0614 < 0.1
```

算法判断：

> 仅更新已有 4 个 PGD 模态还不够，需要增加第 5 个空间-时间模态。

---

## 5. 第 5 个 PGD 模态的固定点其实已经正常收敛

失败模态的固定点历史为：

```text
[3.92483574e-01
 4.04432570e-01
 6.37269744e-01
 2.52311597e-01
 8.76645859e-02
 1.37175919e-02
 2.88354559e-03
 5.73956952e-04
 1.15995490e-04
 2.33724940e-05
 4.71235668e-06]
```

最终：

```text
fixed_point_iterations = 11
fixed_point_converged  = True
```

当前固定点收敛阈值为：

```text
1.0e-5
```

最后固定点变化量：

```text
4.71235668e-06
```

因此：

> 第 5 个新模态并不是因为固定点没有收敛而失败。

这是本阶段第一个重要排除项。

---

## 6. 第 5 个模态也没有表现出明显的线性相关或正交化退化

失败富集诊断：

```text
spatial_novelty                      = 1.0
temporal_significance                = 0.17373662099752832
orthogonality_error                  = 8.901105358732973e-16
plastic_field_invariance_error       = 8.515541232132033e-18
plastic_rate_field_invariance_error  = 3.520974100650146e-20
stress_field_invariance_error        = 1.560347329636251e-14
```

这些结果说明：

### 6.1 空间新颖度

```text
spatial_novelty = 1.0
```

说明该候选空间模态并没有表现为与已有 4 个空间基底接近线性相关。

### 6.2 时间显著度

```text
temporal_significance ≈ 0.174
```

不是一个趋近于零的无意义时间函数。

### 6.3 正交性

```text
orthogonality_error ≈ 8.9e-16
```

接近机器精度。

### 6.4 场不变性

正交化与坐标变换前后：

- 塑性应变修正场；
- 塑性应变率修正场；
- 应力修正场；

均保持到了约 `1e-14` 或更小的相对误差。

因此：

> 本次失败不能简单归因于 Gram-Schmidt 正交化错误，也不能归因于新模态与旧模态严重重复。

---

## 7. 真正触发失败的是完整机械残差变坏

第 5 个候选模态在最终接受判据下得到：

```text
residual_norm_before = 0.04776694562884108
residual_norm_after  = 0.05017696487787036
```

定义完整残差收益：

$$\eta_{\mathrm{res}}=1-\frac{\lVert R_{\mathrm{after}}\rVert}{\lVert R_{\mathrm{before}}\rVert}.$$

得到：

```text
residual_benefit = -0.050453702184678484
```

即：

> 加入第 5 个候选模态并完成当前全部后处理后，完整机械残差反而增加约 5.05%。

当前：

```text
acceptance_tolerance = 0
```

源码判据为：

```text
benefit <= acceptance_tolerance
```

则拒绝富集。

因此本次拒绝是合理的：

```text
-0.05045 <= 0
```

而且这不是 `1e-12` 量级的数值噪声，而是约 5% 的明确恶化。

---

## 8. 第一个重要科研认识：固定点收敛不等于完整残差下降

本阶段通过 100 周期问题首次清楚观察到：

> 一个 PGD 新模态可以满足固定点收敛，但仍然不能保证它降低最终完整机械残差。

当前 residual-LS 富集固定点的收敛条件主要检查：

- 相邻空间-时间模态对的变化；
- 原始固定点映射变化；
- 加速后变化。

而最终新模态是否真正“有用”，是在固定点结束以后才通过完整机械残差收益检查。

因此：

```text
fixed-point converged
```

与：

```text
full residual decreased
```

是两个不同层次的判据。

在 1 到 10 周期范围内，这种区别没有导致整体失败；

在 100 周期下，这种区别第一次成为决定性因素。

---

## 9. Aitken 是否是 100 周期失败的原因

由于当前稳定生产版本采用 OPT6 受保护 Aitken 加速，因此首先需要排除：

> 100 周期失败是否是 OPT6 长时间域副作用。

诊断方式：

- 不修改生产代码；
- 当前进程临时禁止 `_aitken_tail_ready(...)` 返回 True；
- 即完全关闭 residual-LS 固定点中的 Aitken；
- 其他参数全部保持不变；
- 重新跑 100 周期。

### 9.1 无 Aitken 结果

```text
setup_time_s             = 0.992768
solver_time_s            = 102.433402
total_time_s             = 103.426170

converged                = False
termination_reason       = enrichment_failed
failure_reason           = full_residual_benefit_insufficient

attempted_iterations     = 10
committed_iterations     = 9
trial_evaluations        = 14
accepted_pgd_rank        = 4
accepted_indicator       = 6.554750261969e-03
```

最后失败模态：

```text
fixed_point_iterations   = 11
fixed_point_converged    = True
spatial_novelty          = 1.0
temporal_significance    = 0.17373681354221882
residual_norm_before     = 0.04776692517485954
residual_norm_after      = 0.050176925416246715
residual_benefit         = -0.05045332586439111
```

### 9.2 与生产 OPT6 比较

| 指标 | 开启 Aitken | 关闭 Aitken |
|---|---:|---:|
| attempted iterations | 10 | 10 |
| committed iterations | 9 | 9 |
| PGD rank | 4 | 4 |
| fixed-point iterations | 11 | 11 |
| accepted xi | 0.0065547507 | 0.0065547503 |
| residual benefit | -5.045370% | -5.045333% |

两者几乎完全重合。

因此可以明确：

> 当前 100 周期失败不是 OPT6 Aitken 加速造成的。

这使得 OPT6 在当前问题上的嫌疑基本排除。

---

## 10. 第二轮定位：到底是“第 5 模态本身”还是“全部时间函数重优化”导致残差恶化

当前新模态富集流程大致为：

```text
已有 m 个 PGD 模态
→ residual-LS 固定点得到第 m+1 个原始模态
→ 正交化与坐标变换
→ 将全部 m+1 个模态的时间函数重新优化
→ 计算完整机械残差收益
→ 接受或拒绝
```

因此需要区分：

1. 第 5 个原始模态一加入就已经有害；
2. 第 5 个原始模态其实有益，是后续全部时间函数重优化把它变坏。

专门诊断比较三个阶段。

### 10.1 阶段 A：加入第 5 模态之前

```text
A residual_before = 4.776694562884e-02
```

### 10.2 阶段 B：加入已固定点收敛的第 5 原始模态，但不做全部时间函数重优化

```text
B residual_raw_appended = 8.431447845015e-02
```

相对于 A：

```text
raw_mode_benefit = -0.7651218293354
```

即：

> 原始第 5 模态一加入，完整机械残差就恶化约 76.51%。

### 10.3 阶段 C：对 5 个模态全部时间函数重新优化以后

```text
C residual_after_temporal_reopt = 5.017696487787e-02
```

相对于 A：

```text
final_benefit = -0.05045370218468
```

即最终仍比 A 差约 5.05%。

但是 C 相对于 B：

```text
reopt_change = 0.4048831730894
```

即：

> 全部时间函数重优化实际上把 B 的大幅恶化修复了约 40.49%。

因此可以排除：

> 不是全部时间函数重优化把一个原本有用的新模态变坏。

相反：

> 原始第 5 模态在当前时间函数下本身已经是完整残差非下降方向，而全部时间函数重优化只是在尽力补救。

---

## 11. 全部时间函数重优化本身也没有明显病态

第 5 模态加入后，5 模态时间函数重新优化诊断：

```text
candidate_relative_residual       = 0.4180425371327
condition_max                     = 2.360384983392
condition_median                  = 1.037478502045
```

最大条件数约 2.36，中位条件数约 1.04。

因此：

> 当前 5 模态 reduced temporal system 并不存在明显的矩阵病态问题。

这一结果进一步排除了：

- 数值条件数爆炸；
- 5 模态时间基严重相关；
- `np.linalg.lstsq` 在病态系统中失真；

作为主要原因的可能性。

---

## 12. 核心问题进一步缩小：空间模态有问题，还是该空间模态对应的时间函数有问题

前述诊断说明：

- 第 5 模态固定点收敛；
- 空间新颖度良好；
- 正交性良好；
- 全部时间函数重优化没有制造问题；
- 但是使用生产 `_temporal_solve(...)` 得到的第 5 原始空间-时间模态会使完整机械残差恶化。

于是问题集中到：

> 对“同一个第 5 空间模态”，当前顺序时间函数是否真的适合整个 100 周期时间域。

---

## 13. 当前生产时间半步的结构

当前 `latin/tower_pgd_enrichment.py` 中 `_temporal_solve(...)` 采用顺序后向欧拉形式。

其核心结构是：

```text
lambda_0
→ lambda_1
→ lambda_2
→ ...
→ lambda_Nt-1
```

对于每个时间步，只根据：

- 上一步时间系数；
- 当前空间模态；
- 当前搜索方向；
- 当前缺陷；

求解当前时间系数。

因此它具有明显的因果顺序结构：

> 每个时间步一旦确定，不会因为后面更远时间的信息再回头整体调整。

这是目前从一维杆复现中继承并经过短周期验证的工程实现。

---

## 14. whole-time global BE 诊断的目的

为了回答：

> 第 5 个空间模态本身到底有没有下降潜力？

本阶段构造了一个**仅用于诊断**的 whole-time global BE temporal least-squares 解。

它不改变第 5 个空间模态，只改变时间函数求解方式。

设第 5 个固定空间模态为：

```text
p(q)
s(q)
```

当前生产顺序时间函数为：

```text
lambda_seq(t)
```

诊断用全时间域时间函数为：

```text
lambda_global(t)
```

诊断目标是：

> 在同一个后向欧拉完整机械残差定义下，把整个 100 周期时间域中的所有时间系数同时优化。

需要强调：

> 这个 whole-time global BE 只是诊断工具，当前不能等同于原论文 DG0 的精确时间离散。

---

## 15. 同一个第 5 空间模态：顺序时间函数与全时间域时间函数的决定性对比

100 周期诊断结果：

```text
A residual_before             = 4.776694562884e-02
B production_sequential_after = 8.431447845015e-02
C diagnostic_global_BE_after  = 4.374795543223e-02
```

### 15.1 当前生产顺序时间函数

```text
B vs A benefit = -0.7651218293354
```

即：

> 使用当前顺序时间函数，第 5 空间模态使残差恶化约 76.51%。

### 15.2 同一个空间模态 + 全时间域 BE 时间函数

```text
C vs A benefit = 0.08413747506145
```

即：

> 完全相同的第 5 空间模态，如果换成 whole-time global BE 时间函数，反而可以把完整机械残差降低约 8.41%。

### 15.3 global BE 相对于顺序时间函数的改善

```text
C vs B improvement = 0.4811335344012
```

即：

> 全时间域时间函数相对于生产顺序时间函数，将该模态对应的残差降低约 48.11%。

这是本阶段目前最重要的数值证据。

---

## 16. 由此可以明确什么

这项对照试验证明：

> 第 5 个空间模态本身不是一个“天然无用”的空间方向。

因为对于完全相同的空间模态，存在一个时间函数，使得它可以使完整机械残差下降约 8.41%。

因此：

> 100 周期原始失败不能简单归咎为 residual-LS 找到了错误空间方向。

更准确的说法是：

> 当前 residual-LS 找到的第 5 空间模态具有下降潜力，但生产顺序时间更新没有为它找到能够在整个 100 周期时间域中发挥这种下降潜力的时间函数。

---

## 17. 一个非常重要的通俗理解

可以把一个 PGD 模态理解为：

```text
空间模态 = 往哪个方向修正
时间函数 = 每个时刻修正多少
```

100 周期当前失败意味着：

> “往哪个方向修正”其实可以是对的，但“100 个周期里每个时刻修正多少”安排得不够好。

所以：

```text
空间方向可以
+
时间分配不合适
=
完整残差仍然变坏
```

而 whole-time global BE 诊断说明：

```text
完全相同空间方向
+
重新从整个100周期一起安排时间分配
=
完整残差下降
```

---

## 18. 这一结论与此前 fully-reversed 1 周期问题不能混为一谈

项目此前在 fully-reversed 单周期问题中也曾研究过：

- PGD 固定点振荡；
- spatial half-step formulation；
- whole-time BE 时间最小化；
- residual-LS 空间求解。

当时的重要发现是：

> 单纯改成 whole-time BE 时间求解，并不能消除当时的 period-3 / fixed-point pathology。

因此此前 fully-reversed 单周期问题的主要矛盾后来指向：

> 空间半步的离散一致性与 residual-LS formulation。

而本次 100 周期 asymmetric 问题的诊断却显示：

> 同一个空间模态在顺序时间函数下是非下降方向，在 whole-time 时间函数下却可以成为下降方向。

所以目前应该严格区分：

### 18.1 fully-reversed 单周期问题

主要暴露：

> spatial half-step formulation / fixed-point map 的问题。

### 18.2 asymmetric 100 周期问题

目前主要暴露：

> 长时间域下 sequential temporal update 无法充分发挥一个有效空间方向的问题。

这两者是不同的数值机制。

---

## 19. 对 temporal amplitude “相对误差”的一个纠正

global temporal diagnostic 输出过：

```text
temporal_amplitude_relL2_diff = 4.128825249494e-03
```

这个量不能直接解释成：

> 两个时间函数只差 0.41%。

原因是诊断脚本为了避免分母过小，使用了：

```text
max(1.0, ||lambda_global||)
```

而时间函数本身约为 `1e-4` 量级。

因此这个量实际上更接近：

> 绝对 L2 差异被 1 归一化。

真正更直观的幅值对比为：

```text
max_abs_sequential_lambda = 1.086196404677e-04
max_abs_global_lambda     = 3.791513360069e-05
```

两者最大幅值差异明显。

这一点不改变残差结论，但后续正式文档与论文描述中不能把 `4.13e-3` 当成真正的相对误差百分比。

---

## 20. 为什么还不能直接把 whole-time global BE 作为正式生产算法

虽然 whole-time global BE 对失败第 5 模态非常有效，但它目前仍然只能作为：

> 因果诊断工具。

原因包括：

### 20.1 它不是当前原论文 fidelity 路线已经严格证明的时间离散

当前项目一直坚持：

- 原论文是理论基准；
- 工程实现与原论文内容必须严格区分；
- 当前后向欧拉 tower 实现是从一维杆复现中继承并验证的离散；
- 不能因为某个替代算法数值效果好，就直接声称它等价于原论文 Eq. (72) / DG0。

### 20.2 whole-time global BE 改变整个 PGD 固定点轨迹

PGD 新模态是空间函数与时间函数交替构造的：

```text
time function
→ spatial function
→ time function
→ spatial function
→ ...
```

因此如果从第 1 个模态开始就更换时间半步：

> 后续空间模态也会跟着改变。

所以：

> “对原生产轨迹中的第 5 空间模态，global BE 更好”

并不自动等价于：

> “把整个算法从第一步开始改成 global BE 后，整个 100 周期求解一定会更好”。

这一点必须单独验证。

---

## 21. 完整求解器级对照：从第一个模态开始使用 global BE temporal half-step

为了进行因果验证，本阶段继续做了一个完整求解器级诊断：

- 生产源码不修改；
- 当前 Python 进程中临时替换新模态 `_temporal_solve(...)`；
- 使用 whole-time global BE temporal half-step；
- residual-LS 空间半步保持不变；
- 外层 Trial-A / Trial-B 保持不变；
- 饱和规则保持不变；
- 完整残差收益保护保持不变；
- 100 周期重新从头求解。

### 21.1 结果

```text
setup_time_s             = 0.970597
solver_time_s            = 23.513718
total_time_s             = 24.484315

converged                = False
termination_reason       = enrichment_failed
failure_reason           = full_residual_benefit_insufficient

attempted_iterations     = 4
committed_iterations     = 3
trial_evaluations        = 5

accepted_pgd_rank        = 1
accepted_indicator       = 1.704194979217e-02
total_modes_added        = 1
```

### 21.2 与原生产算法的 100 周期失败位置不同

原生产算法：

```text
attempted iterations = 10
committed iterations = 9
accepted PGD rank     = 4
accepted xi           = 0.00655475
第5模态失败
```

global BE temporal half-step 从头替换后：

```text
attempted iterations = 4
committed iterations = 3
accepted PGD rank     = 1
accepted xi           = 0.01704195
第2模态附近失败
```

也就是说：

> global BE 从算法第一步开始介入以后，整个 PGD 搜索轨迹被明显改变。

---

## 22. global BE 完整求解器试验并没有推翻前面的单模态诊断

这两个试验回答的是不同问题。

### 22.1 单模态 fixed-spatial diagnostic 回答

问题：

> 原生产算法在第 10 次迭代找到的那个第 5 空间模态，本身有没有下降潜力？

答案：

> 有。

证据：

```text
same spatial fifth mode
sequential temporal solve: residual worsens 76.51%
global temporal solve: residual improves 8.41%
```

### 22.2 global BE full-solver diagnostic 回答

问题：

> 如果从第一个 PGD 模态开始就把 temporal half-step 全部换成 global BE，100 周期整体算法能否直接恢复正常？

答案：

> 不能直接恢复。

并且：

> 算法搜索轨迹改变，在只有 1 个已接受模态时就进入新的富集失败。

所以正确结论是：

> global BE 证明了原第 5 空间模态有潜在价值，但 global BE 本身尚未成为一个可直接替代当前生产 temporal half-step 的稳定完整算法。

---

## 23. 当前已经可以排除的原因

经过本阶段连续诊断，目前对于原生产 100 周期失败，可以较有把握排除以下原因。

### 23.1 不是简单运行时间过长

失败是：

```text
enrichment_failed
```

而不是：

- timeout；
- MemoryError；
- Python 异常崩溃；
- Newton 不收敛；
- FOM 失败。

### 23.2 不是 Aitken 加速导致

开启和关闭 Aitken 得到几乎完全一致的：

- 失败位置；
- xi；
- PGD rank；
- fixed-point history；
- residual benefit。

### 23.3 不是 fixed-point 本身未收敛

第 5 模态：

```text
11 次 fixed-point
最后 chi ≈ 4.7e-6
```

满足 `1e-5`。

### 23.4 不是新空间模态明显线性相关

```text
spatial_novelty = 1.0
```

### 23.5 不是 Gram-Schmidt 变换破坏物理场

场不变性误差约：

```text
1e-14 或更小
```

### 23.6 不是 all-mode temporal reoptimization 把好模态变坏

实际上：

```text
raw appended residual = 0.0843145
reoptimized residual  = 0.0501770
```

all-mode temporal reoptimization 是在修复，而不是制造失败。

### 23.7 不是 5 模态 reduced temporal system 病态

条件数：

```text
max    ≈ 2.36
median ≈ 1.04
```

十分健康。

---

## 24. 当前最重要的“已证实事实”

以下结论属于本阶段**实测事实**。

### 24.1 事实 A

当前稳定 LATIN-PGD 在：

```text
1、2、5、10 周期
```

均能正常收敛。

### 24.2 事实 B

当前稳定 LATIN-PGD 在 100 周期首次出现：

```text
enrichment_failed / full_residual_benefit_insufficient
```

### 24.3 事实 C

失败时：

```text
accepted rank = 4
accepted xi   = 6.55475e-3
```

并不是已经形成高阶 PGD 基后才发生。

### 24.4 事实 D

失败的第 5 个模态固定点已经正常收敛。

### 24.5 事实 E

失败的第 5 空间模态：

- 空间新颖度良好；
- 正交性良好；
- 时间显著度非零；
- 场不变性正常。

### 24.6 事实 F

在生产顺序时间函数下，该第 5 模态使完整机械残差恶化约 76.51%。

### 24.7 事实 G

all-mode temporal reoptimization 将该恶化大幅修复，但最终仍比富集前差约 5.05%，因此新模态被正确拒绝。

### 24.8 事实 H

关闭 Aitken 后，失败几乎完全复现，因此 Aitken 不是根因。

### 24.9 事实 I

保持完全相同的第 5 空间模态，若改用诊断用 whole-time global BE 时间函数：

```text
完整机械残差可降低约 8.41%
```

### 24.10 事实 J

如果从整个算法第一步开始就把新模态 temporal half-step 全部替换为 global BE：

- 原失败轨迹改变；
- 只接受 1 个模态；
- 第 2 模态附近再次出现 `full_residual_benefit_insufficient`；
- 整体 100 周期仍未收敛。

---

## 25. 当前“强证据支持，但仍属于机制解释”的认识

以下内容目前不是原论文直接结论，也不是形式化数学证明，但已经得到较强数值证据支持。

### 25.1 当前生产 sequential temporal update 存在长时间域局限

对原生产第 5 空间模态：

```text
sequential temporal function → residual worsens
global whole-time temporal function → residual decreases
```

因此强烈说明：

> 当前顺序 temporal half-step 在 100 周期长时间域下没有找到该空间模态对应的 whole-time residual descent temporal coordinate。

### 25.2 长时间域下“固定点自洽”与“全时间域残差下降”发生分离

第 5 模态：

```text
fixed-point converged
```

但：

```text
full residual worsened
```

因此：

> 当前 fixed-point map 的自洽性不再自动保证 whole-history mechanical residual descent。

### 25.3 当前 100 周期问题更接近时间方向的 formulation 问题，而不是单纯空间方向失效

原因是：

> 相同空间模态在更合适的 whole-time temporal coordinate 下可以成为下降方向。

但是这一判断仍需要后续与原论文 Eq. (72) / DG0 formulation 对照。

---

## 26. 当前仍然没有被证明的事情

必须严格避免把以下内容提前写成结论。

### 26.1 还不能说原论文 Eq. (72) 本身在 100 周期失效

当前失败的是：

> 我们当前 tower-v1 的 sequential BE implementation。

它不等于：

> 原论文 DG0 时间离散的完整数学实现。

因此不能说：

> “原论文 LATIN-PGD 不适用于 100 周期”。

### 26.2 还不能说 whole-time global BE 就是正确正式修正

它目前只证明：

> 对某个固定空间模态，它能够找到更好的时间函数。

但 full-solver 从头替换后仍然失败。

### 26.3 还不能说 residual-LS 空间方向完全没有问题

对于原第 5 模态，它具有下降潜力；

但 global BE full-solver 改变轨迹后，第 2 模态附近又失败。

因此：

> 空间半步与时间半步之间的耦合仍然需要进一步诊断。

### 26.4 还不能给出 1000 周期结论

100 周期尚未稳定求解，因此：

> 现在没有科学依据直接进入 1000 周期 benchmark。

### 26.5 还不能给出 100 周期 FOM/LATIN 公平速度比

原 100 周期 fair benchmark 中 LATIN 后续失败，因此整套 pair validation 没有完成。

即使 FOM 在前面已经执行，也没有形成当前认可的完整正式 paired timing record。

因此：

> 100 周期效率比较尚未成立。

---

## 27. 为什么当前不应该直接跑 1000 周期

本阶段开始时原计划：

```text
100 周期
→ 1000 周期
```

用于观察长期 scaling。

但是 100 周期已经发现：

```text
算法鲁棒性问题
```

而不是单纯：

```text
性能问题
```

因此如果现在直接跑 1000 周期：

- 很可能在更早或类似位置再次失败；
- 无法获得有意义的长期效率数据；
- 反而会混淆算法失败与复杂度增长；
- 还会带来明显内存压力。

所以研究顺序已经合理调整为：

```text
先解决或解释100周期富集问题
→ 再恢复100周期完整精度/效率 benchmark
→ 再决定是否进入1000周期
```

---

## 28. 对 1000 周期内存问题的阶段认识

当前 whole-history LATIN-PGD 状态场尺寸约为：

```text
Nt × Nq
```

1000 周期：

```text
Nt = 40001
Nq = 320
```

单个 `float64` 场：

```text
40001 × 320 × 8 bytes
```

约为 100 MB 量级。

而当前 LATIN 状态包含多个：

- 应力；
- 总应变；
- 塑性应变；
- 累积塑性变量；
- 损伤；
- 搜索方向；
- 局部状态；
- baseline / trial；
- forcing / residual；

等 whole-history 数组。

因此：

> 1000 周期不仅是 CPU 时间问题，还很可能成为 whole-history 内存复杂度问题。

但这一点目前尚未正式 profiling。

---

## 29. 当前研究问题已经从“效率”转向“长时间域算法一致性”

上一阶段的核心问题是：

> 优化后的 FOM 与优化后的 LATIN-PGD 到底谁更快。

本阶段之后，优先级发生变化。

当前更重要的问题变成：

> 为什么当前 tower LATIN-PGD 在 100 周期时，新模态的固定点已经收敛，但对应的 sequential temporal coordinate 却不能保证完整机械残差下降。

进一步来说：

> 当前 sequential BE temporal half-step 与 residual-LS spatial half-step 在长时间域下是否仍然构成一致的交替最小化/固定点结构。

这个问题目前比继续做 wall-time profiling 更重要。

---

## 30. 与原论文 fidelity 的关系

本项目始终要求严格区分：

1. 原论文明确内容；
2. 我们严格推导的内容；
3. 一维杆中已经验证的内容；
4. tower implementation 的工程选择；
5. 本阶段诊断用替代算法。

### 30.1 当前 production sequential BE

属于：

> 从一维杆复现中继承并经过验证的工程离散实现。

### 30.2 residual-LS spatial half-step

属于：

> 为解决 tower paper-Galerkin 空间半步实际数值问题而形成的工程增强 formulation。

### 30.3 whole-time global BE temporal solve

属于：

> 本阶段诊断工具。

不能声称为：

> 原论文 Eq. (72) DG0 的精确复现。

### 30.4 下一步 paper-fidelity 目标

应该重新回到原论文时间离散，特别是：

> Eq. (72) 与 DG0 时间离散到底对应怎样的离散代数系统。

只有在这一点严格厘清后，才能判断：

- 当前 sequential BE 是否只是一个短周期近似；
- whole-time coupling 是否更接近论文原始离散；
- 长时间域失败是否来自时间离散 fidelity 偏差。

---

## 31. 本阶段最重要的学术认识

如果要向导师用一句话概括本阶段的新发现，可以表述为：

> 当海上风机塔筒 LATIN-PGD 从 10 周期扩展到 100 周期时，首次出现了 PGD 新模态固定点已收敛但完整机械残差反而增加的富集失效；进一步诊断表明，该失败空间模态本身具有下降潜力，而当前逐时间步顺序后向欧拉时间函数未能在整个长时间域中充分发挥该空间模态的残差降低能力，说明当前 tower LATIN-PGD 的长时间域瓶颈已经从单纯计算效率问题转向空间-时间 PGD 交替更新的一致性问题。

这个表述目前有实测证据支持，但最后一句“空间-时间一致性问题”仍应作为：

> 强证据支持的机制判断，

而不是已经完成形式数学证明的结论。

---

## 32. 本阶段对“低秩”的认识也需要更新

此前 10 周期时：

```text
Nt = 401
PGD rank = 21
```

低秩性非常明显。

100 周期失败时：

```text
Nt = 4001
accepted PGD rank = 4
```

不能简单把：

```text
rank = 4
```

解释成：

> 100 周期问题只需要 4 阶。

因为算法是在第 5 模态富集时失败，中途终止。

因此：

> 4 阶只是失败前最后一个可接受 PGD 基规模，不是 100 周期最终收敛阶数。

这一点后续汇报时必须明确。

---

## 33. 为什么 100 周期 solver time 约 102 s 也不能直接用于正式效率比较

生产算法 100 周期失败诊断中：

```text
solver_time ≈ 102 s
```

这个时间只代表：

> 算法运行到第 10 次外层迭代、第 5 模态富集失败位置所消耗的时间。

它不是：

> 100 周期收敛 LATIN-PGD 的完整 wall-time。

因此不能拿：

```text
102 s
```

直接和 FOM 100 周期时间做方法效率比较。

同样：

global BE full-solver diagnostic 中：

```text
solver_time ≈ 23.5 s
```

也不是更快的 100 周期收敛算法，因为它在：

```text
rank = 1
```

时就提前失败。

所以：

> 失败运行时间不能作为公平性能指标。

---

## 34. 当前 global BE full-solver 失败的意义

global BE full-solver 试验说明：

> “修正某一个失败模态的时间函数”与“替换整个 PGD 固定点中的 temporal half-step”不是同一件事。

因为 PGD 是交替构造：

```text
时间函数改变
→ 下一空间模态改变
→ 再下一时间函数改变
→ 整条富集轨迹改变
```

因此从第一个模态开始改 temporal half-step：

> 会改变整个 PGD 搜索路径。

这也是为什么：

```text
原生产算法：第5模态失败
global BE full solver：第2模态附近失败
```

并不构成逻辑矛盾。

---

## 35. 下一步应优先做什么

在正式修改算法之前，下一步应该先完成：

> global BE full-solver 新失败点的详细诊断。

需要回答：

1. 第 2 个候选模态是否 fixed-point 收敛；
2. 空间新颖度是多少；
3. 时间显著度是多少；
4. residual before / after 分别是多少；
5. 原始模态加入前后残差如何；
6. all-mode temporal reoptimization 是修复还是恶化；
7. 新失败是否与原生产第 5 模态属于同一种机制；
8. global BE 时间半步是否改变了空间 fixed-point map 的收敛支路。

只有得到这些信息以后，才决定：

- 是否继续沿 global BE 诊断；
- 是否回到原论文 DG0 时间离散推导；
- 是否需要构造新的 temporal half-step；
- 是否需要重新设计固定点收敛/富集收益的一致性判据。

---

## 36. 推荐下一阶段研究顺序

### 阶段 A：诊断 global BE full-solver 的新失败

只做诊断，不修改生产源码。

目标：

> 比较原生产“第 5 模态失败”与 global BE full-solver“第 2 模态失败”的机制是否相同。

### 阶段 B：回到原论文 Eq. (72) / DG0

重新严格推导：

- 原论文时间函数试探空间；
- DG0 离散；
- 时间 slab 上的未知量定义；
- Eq. (72) 的离散形式；
- 是否形成全时间域耦合；
- 是否可以顺序求解；
- 当前 BE 版本与论文 DG0 的差别。

### 阶段 C：建立 paper-fidelity temporal diagnostic

如果论文 DG0 明确导出不同于当前 sequential BE 的系统：

> 先做本地诊断实现，而不是直接改生产代码。

### 阶段 D：重新验证 1 / 2 / 5 / 10 / 100 周期

必须重新检查：

- convergence；
- PGD rank；
- xi；
- FOM matched accuracy；
- fixed-point iterations；
- enrichment benefit；
- wall-time。

### 阶段 E：只有 100 周期稳定后，再考虑 1000 周期

1000 周期应同时记录：

- wall-time；
- peak memory；
- PGD rank；
- outer iterations；
- fixed-point iterations；
- accuracy；
- internal state stability。

---

## 37. 当前不建议做的事情

### 37.1 不建议直接放宽 acceptance_tolerance

例如：

```text
acceptance_tolerance < 0
```

强行接受残差恶化模态。

这会破坏当前完整残差保护机制。

### 37.2 不建议直接增加 max_fixed_point_iterations

第 5 模态已经 fixed-point converged。

增加次数不能解决：

```text
residual benefit < 0
```

的根本问题。

### 37.3 不建议直接修改 Aitken

已经实测排除。

### 37.4 不建议直接把 global BE 写入 production

full-solver 诊断已经证明：

> 它从头介入后并不能直接恢复 100 周期收敛。

### 37.5 不建议现在跑 1000 周期

100 周期鲁棒性尚未解决。

### 37.6 不建议现在继续做正式公平 speedup 图

因为 100 周期 LATIN 尚未形成完整收敛运行。

---

## 38. 本阶段使用的本地诊断脚本

以下脚本均属于本地诊断，目前不应使用 `git add .` 将其全部加入版本库。

### 38.1 100 周期公平试跑

```text
tower_fom_latin_fair_benchmark.py
```

### 38.2 100 周期富集失败详细诊断

```text
tower_latin_100cycle_enrichment_failure_diagnostic.py
```

### 38.3 100 周期关闭 Aitken 对照

```text
tower_latin_100cycle_no_aitken_diagnostic.py
```

### 38.4 第 5 模态三阶段残差诊断

```text
tower_latin_100cycle_enrichment_stage_diagnostic.py
```

### 38.5 同空间模态 sequential 与 global temporal solve 对照

```text
tower_latin_100cycle_global_temporal_ls_diagnostic.py
```

### 38.6 global BE temporal half-step 完整求解器对照

```text
tower_latin_100cycle_global_temporal_full_solver_diagnostic.py
```

这些脚本目前主要用于：

> 保留科学诊断证据，

而不是生产 API。

---

## 39. 当前 Git checkpoint 的意义

当前稳定远端 checkpoint：

```text
0b061cb6bc301eaa0e510120215e52ab018e6f9f
```

包含：

> 1 到 10 周期 FOM,opt 与 LATIN-PGD,opt 正式公平效率总结。

本阶段 100 周期工作发生在该 checkpoint 之后。

重要的是：

> 本阶段尚未修改生产源码。

因此：

```text
0b061cb
```

仍然是当前生产代码与正式 1 到 10 周期 benchmark 的稳定基准。

本阶段的新内容主要是：

- 长时间域试跑；
- failure diagnostics；
- alternative temporal diagnostics；
- 新的科研认识。

---

## 40. 建议本阶段新的正式 checkpoint

建议将本 Markdown 阶段总结作为新的独立 docs commit。

建议文件名：

```text
docs/2026-09-03-tower-latin-100cycle-long-horizon-temporal-diagnostics.md
```

建议 commit message：

```text
docs: summarize tower LATIN 100-cycle temporal diagnostics
```

该 commit 应：

- 只包含本总结 Markdown；
- 不包含本地诊断脚本；
- 不包含 CSV / JSON；
- 不包含 `.pstats`；
- 不包含临时 patch；
- 不包含其他 untracked 文件。

---

## 41. 最终阶段结论

截至本阶段，可以形成以下结论。

### 41.1 关于 100 周期可计算性

当前稳定 tower LATIN-PGD：

> 尚不能完成 100 周期收敛求解。

失败形式为：

```text
enrichment_failed / full_residual_benefit_insufficient
```

### 41.2 关于失败位置

原生产算法在：

```text
第10次外层尝试
已接受4个PGD模态
accepted xi ≈ 6.55e-3
```

时尝试第 5 模态失败。

### 41.3 关于固定点

第 5 模态：

> 固定点正常收敛。

因此不是固定点迭代次数不足。

### 41.4 关于空间模态质量

第 5 空间模态：

> 具有良好的空间新颖度、正交性和时间显著度。

并且：

> 存在一个 whole-time temporal coordinate 能使它降低完整残差约 8.41%。

因此不能简单判定为空间方向本身无效。

### 41.5 关于当前 sequential temporal update

当前生产顺序时间函数下，同一空间模态使残差恶化约 76.51%。

这为：

> sequential temporal update 的长时间域局限

提供了强数值证据。

### 41.6 关于 Aitken

关闭 Aitken 后结果几乎完全一致。

因此：

> Aitken 不是 100 周期失败根因。

### 41.7 关于 all-mode temporal reoptimization

它不是问题制造者。

相反：

> 它将原始第 5 模态造成的大幅残差恶化显著修复，只是仍不足以达到净下降。

### 41.8 关于 whole-time global BE

对固定失败空间模态：

> 有效。

作为完整 temporal half-step 从第一模态开始替换：

> 仍未使 100 周期整体收敛，并改变了 PGD 搜索轨迹。

因此：

> 目前只能作为诊断工具，不能直接升级为正式生产算法。

### 41.9 关于下一阶段

当前最合理的下一步不是：

```text
1000-cycle benchmark
```

也不是：

```text
直接修改生产代码
```

而是：

> 先诊断 global BE full-solver 新出现的第 2 模态富集失败，再回到原论文 Eq. (72) / DG0 时间离散进行 paper-fidelity 推导与对照。

---

## 42. 一句话总结

> 1 到 10 周期阶段主要回答了“优化后的 LATIN-PGD 到底快不快”；100 周期阶段首次暴露出一个更本质的问题：当前 PGD 新模态可以在固定点意义下收敛，但其顺序时间函数在长时间域中不再保证完整机械残差下降，而同一个空间模态在 whole-time 时间优化下又能够恢复下降潜力，这表明下一阶段的核心研究对象已经从性能优化转向长时间域下 PGD 空间-时间交替更新与原论文时间离散的一致性问题。
