# 海上风机塔筒 LATIN-PGD 内部验证闭环与 FOM 前置判定

日期：2026-08-20  
当前分支：`feature/offshore-wind-turbine-tower-fatigue`  
进入 I-6 阶段时的基线提交：`369d424ab4a16d3d1e2f0b60973d73a10215169d`

> 本文档用于系统总结海上风机塔筒 LATIN-PGD 方法从塔筒平衡算子、PGD 富集、residual-LS 空间更新，到完整 outer LATIN 收敛验证的全过程，并给出是否具备进入全阶模型 FOM 对比阶段的明确判定。
>
> 为避免 PyCharm Markdown 预览对行内 LaTeX 支持不稳定，本文档不再使用 `\(...\)`、`$...$` 或 `\[...\]` 等行内数学语法。所有公式均采用普通 Unicode 数学符号或代码块形式书写，以保证在 PyCharm 中稳定显示。

---

# 1. 当前阶段的核心问题

本阶段不是直接回答“LATIN-PGD 是否已经正确复现 FOM”，而是先回答更基础的问题：

> 当前面向海上风机塔筒建立的 LATIN-PGD 计算链条，是否已经在内部数值层面保持平衡性、残差下降、PGD 增秩、事务一致性以及 outer LATIN 收敛，从而具备进入 FOM 对比验证的资格？

为回答这一问题，当前内部验证被拆分为 I-0 至 I-6 七个阶段：

1. I-0：冻结可追溯的代码基线；
2. I-1：验证塔筒平衡算子；
3. I-2：定位第四模态固定点失稳，并验证 residual-LS 空间更新；
4. I-3：验证 matrix-free residual-LS 与稠密诊断形式的数值等价性；
5. I-4：将 residual-LS 正式接入 PGD enrichment；
6. I-5：验证完整 outer LATIN-PGD 是否能够继续迭代并最终收敛；
7. I-6：进行进入 FOM 前的总回归与方法边界闭环。

当前 I-6 已完成 31 个联合回归测试，全部通过。

---

# 2. 必须长期保留的理论边界

这一部分非常重要。后续无论是论文、阶段总结、代码说明，还是与 FOM 的对比，都必须明确区分下面三类内容。

## 2.1 原论文明确给出的 LATIN-PGD 内容

原论文明确给出了以下理论框架：

- LATIN 的局部阶段和全局阶段交替；
- 空间函数与时间函数分离的 PGD 表示；
- 固定空间基时的时间函数更新；
- 新空间模态和新时间函数的富集；
- 空间固定点与时间固定点交替更新；
- 原论文 Eq. (65)–(71) 对应的空间 Galerkin 更新；
- 原论文 Eq. (72) 对应的时间残差最小化思想；
- 富集后的 Gram-Schmidt 正交化；
- 富集后重新修正时间坐标。

因此：

> 后续涉及 Eq. (58)–(72) 的理论解释时，必须以论文原式为理论基准。

但同时也必须避免把当前工程实现中新增的 residual-LS 路径误写成原论文 Eq. (65)–(71) 的直接等价形式。

---

## 2.2 塔筒有限元离散是对原理论的结构化翻译

当前塔筒采用纤维梁柱离散。为了把原论文中的全局可容许条件映射到塔筒有限元体系，我们构造了塔筒平衡算子。

对于一个空间塑性应变型修正模态 p，通过平衡投影得到应变修正 ε：

```text
ε = E p
```

对应的应力修正模态 s 为：

```text
s = C0 (ε - p)
```

其中：

- p：塑性应变型空间函数；
- ε：通过塔筒平衡算子得到的应变修正；
- s：对应的应力修正模态；
- C0：弹性刚度；
- E：塔筒平衡投影算子。

对于这一类“齐次修正应力模态”，自由自由度平衡条件应满足：

```text
Hᵀ M s ≈ 0
```

其中：

- H：结构位移到材料点应变的映射；
- M：材料点积分权重对应的度量；
- Hᵀ M s：该应力修正模态在自由自由度上的平衡残差。

需要特别强调：

> `Hᵀ M s ≈ 0` 针对的是由平衡算子生成的齐次修正应力模态，而不是包含外荷载作用的总物理应力。

对于受真实外荷载作用的塔筒，总体平衡应理解为：

```text
内部力 - 外部荷载 ≈ 0
```

而不是简单要求总应力场对应的 `Hᵀ M s` 为零。

---

## 2.3 residual-LS 不是原论文空间 Galerkin 方程的代数改写

这是当前阶段最需要澄清的方法学边界。

对于固定的时间函数 λₙ 及其时间导数 λ̇ₙ，当前 residual-LS 空间更新使用如下机械残差：

```text
rₙ(p) = λ̇ₙ p - λₙ Hσ,n Aσ p + Δₙ
```

其中：

- p：待求空间塑性应变模态；
- Aσ p：由 p 通过塔筒平衡算子生成的应力模态；
- Hσ,n：第 n 个时间节点对应的搜索方向参数；
- Δₙ：当前 frozen global stage 中的已知残差项。

定义：

```text
Bₙ = λ̇ₙ I - λₙ Hσ,n Aσ
```

则有：

```text
rₙ(p) = Bₙ p + Δₙ
```

residual-LS 的目标是直接最小化加权机械残差：

```text
Jₕ(p)
= 1/2 · Σₙ wₙ · rₙ(p)ᵀ · M · Hσ,n⁻¹ · rₙ(p)
```

这一最小二乘形式会自然包含三类时间权重结构：

```text
λ̇² Hσ⁻¹
λ λ̇
λ² Hσ
```

而原论文空间 Galerkin 更新中自然出现的主要结构是：

```text
〈Hσ λ²〉
〈λ λ̇〉
〈Δ λ〉
```

两者虽然都服务于“求一个新的空间模态”，但离散变分结构并不相同。

因此后续应使用下面这些表述：

- residual-LS 空间富集；
- 基于机械残差最小二乘的空间更新；
- 塔筒 LATIN-PGD 的 residual-LS enrichment strategy；
- 经过内部验证的工程化扩展。

不应写成：

> residual-LS 是原论文 Eq. (65)–(71) 的另一种代数写法。

---

# 3. 当前时间离散与 residual-LS 离散方式

当前成功的 residual-LS 路径沿用了诊断阶段已经验证有效的一组离散选择。

## 3.1 时间更新

当前时间函数更新采用 backward Euler。

这意味着时间导数在离散层面采用前后时间节点的后向差分关系。

## 3.2 residual-LS 空间积分

当前 residual-LS 空间最小二乘采用节点型梯形求积，并且包含初始时间节点 t₀。

这一点与当前 paper-Galerkin 空间更新使用的时间离散习惯并不完全相同。

因此现阶段应当把它理解为：

> 一种已经通过当前内部 benchmark 验证的有效离散组合。

而不能提前声称：

> 它已经被严格证明为原论文时空弱式的唯一一致离散。

这一问题应在后续 FOM-4 离散收敛阶段继续检查，包括：

- 时间步长细化；
- 梯形求积与其他求积方式对比；
- 是否包含 t₀ 的敏感性；
- backward Euler 与空间 residual-LS 离散之间的一致性。

---

# 4. I-0：代码基线冻结

## 4.1 目的

I-0 的作用不是做数值验证，而是先建立可靠的版本边界。

在开始第四模态问题诊断之前，需要确保：

- 当前代码可追溯；
- 后续所有诊断都可以与一个固定版本进行对比；
- 一旦出现错误，可以明确判断问题来自哪一阶段修改；
- 不把多个大阶段混在一个 commit 中。

## 4.2 阶段结论

I-0 完成后，后续的 residual-LS 诊断、matrix-free 重构、enrichment 接入和 outer solver 接入都分别保留了独立 checkpoint。

结论：

```text
I-0：PASS
```

---

# 5. I-1：塔筒平衡条件验证

## 5.1 核心问题

如果 PGD 新增空间模态本身都不能满足结构平衡，那么后续固定点收敛和残差下降都没有物理意义。

因此 I-1 首先检查塔筒平衡算子。

## 5.2 验证内容

验证包含两类：

### A. 算子级验证

对任意空间源项 p：

1. 通过平衡算子生成 ε；
2. 构造应力模态 s；
3. 计算自由自由度平衡残差 `Hᵀ M s`；
4. 检查该残差是否达到数值精度。

### B. 真实塔筒离散验证

在真实塔筒的：

- 梁单元网格；
- Gauss 积分点；
- 环向纤维；
- 材料点度量；

条件下重复平衡检查。

## 5.3 代表性结果

第四模态诊断中，不同候选空间更新得到的应力模态平衡残差均达到接近机器精度的水平。

例如 residual-LS 路径中代表性的相对平衡残差约为：

```text
6.4 × 10⁻¹⁶ ～ 8.0 × 10⁻¹⁶
```

这说明第四模态 period-3 失稳并不是由“空间应力模态不平衡”造成的。

## 5.4 阶段结论

```text
I-1：PASS
```

---

# 6. I-2：第四模态固定点失稳诊断

这是整个内部验证过程中最关键的一步之一。

## 6.1 原始问题

在完整 outer LATIN-PGD 过程中，前三个 PGD 模态可以正常建立。

当算法进入第四模态 enrichment 时，原有 paper-Galerkin 空间更新与当前时间更新组合无法达到固定点收敛。

最初可能的解释包括：

- 固定点最大迭代次数不足；
- 模态正交性不足；
- 新模态与旧 basis 太接近；
- 固定点初始 seed 不合理；
- 简单 damping 不足；
- 时间离散导致数值振荡；
- 空间 Galerkin 离散与当前时间更新存在不兼容。

经过逐项诊断后，前面几类简单原因均没有成为主要解释。

## 6.2 period-3 轨道

将固定点迭代次数提高后发现，误差不是缓慢下降，而是稳定地在三个状态之间循环。

即固定点指标表现出类似：

```text
χₖ
→ 0.577...
→ 0.658...
→ 0.603...
→ 0.577...
→ ...
```

并且 lag-3 差异很小。

这说明问题本质不是“还没迭代够”，而是：

> 固定点映射进入了稳定的 period-3 轨道。

---

## 6.3 A/B/C 对比诊断

为了定位问题来源，构造三组对比。

### Case A

```text
paper-Galerkin 空间更新
+ 当前 backward-Euler 时间更新
+ 不增加新的模态内正交化
```

结果：

```text
未收敛
200 次固定点迭代后仍处于 period-3
最终 χ ≈ 0.6035
```

### Case B

```text
paper-Galerkin 空间更新
+ 当前 backward-Euler 时间更新
+ 增加与 1D 验证一致的模态内 M-正交化
```

结果：

```text
仍未收敛
200 次固定点迭代后仍处于 period-3
最终 χ ≈ 0.6063
```

说明：

> 单纯增加正交化并不能解决第四模态失稳。

### Case C

```text
residual-LS 空间更新
+ 相同 backward-Euler 时间更新
+ 相同模态内 M-正交化
```

结果：

```text
fixed-point converged = True
fixed-point iterations = 27
final χ ≈ 7.55 × 10⁻⁷
```

与此同时：

```text
raw-mode residual benefit ≈ 26.44%
maximum basis overlap ≈ 5.4 × 10⁻¹⁶
equilibrium relative residual ≈ 6.4 × 10⁻¹⁶
```

## 6.4 对这一结果的解释

B 和 C 之间唯一的核心变化是空间 half-step。

因此当前最有力的证据是：

> 第四模态 period-3 失稳与空间更新形式密切相关，而不是由 basis orthogonality 单独引起。

同时需要保持谨慎：

> 由于 residual-LS 还沿用了与 1D Case C 相同的梯形节点求积，因此不能把成功完全归因于“变分形式改变”这一单一因素。

这一点必须在后续离散化敏感性研究中继续验证。

## 6.5 阶段结论

```text
I-2：PASS
```

---

# 7. I-3：matrix-free residual-LS 等价性验证

## 7.1 为什么必须 matrix-free

诊断阶段可以显式构造应力映射：

```text
s = Aσ p
```

但真实塔筒中材料点数量会随着：

- 单元数量；
- Gauss 点数量；
- 环向纤维数量；
- 径向纤维数量；

快速增长。

显式构造 Aσ 会带来严重的内存和计算开销。

因此必须把 residual-LS 改写成 matrix-free 形式。

---

## 7.2 matrix-free 最小二乘形式

定义：

```text
Bₙ p = λ̇ₙ p - λₙ Hσ,n Aσ p
```

进一步定义加权算子：

```text
Dₙ = diag( sqrt(wₙ M / Hσ,n) )
```

则每个时间节点对应：

```text
Cₙ p = Dₙ Bₙ p
```

右端项：

```text
bₙ = -Dₙ Δₙ
```

最终求解：

```text
minimize || C p - b ||₂
```

使用 SciPy LSMR 求解，不再构造稠密矩阵 C。

---

## 7.3 等价性检查

I-3 分别比较：

1. matrix-free forward action；
2. 显式 dense forward action；
3. matrix-free adjoint action；
4. 显式 dense adjoint action；
5. adjoint identity；
6. LSMR solution；
7. dense stacked least-squares solution；
8. 最终 fixed-point pair。

代表性结果：

```text
RHS relative error                  = 0
forward relative error              ≈ 2.10 × 10⁻¹⁶
adjoint relative error              ≈ 6.46 × 10⁻¹⁶
adjoint identity error              ≈ 1.74 × 10⁻¹⁷
LSMR vs dense spatial raw error     ≈ 5.36 × 10⁻¹²
full fixed-point spatial error      ≈ 1.07 × 10⁻¹¹
full fixed-point stress error       ≈ 9.39 × 10⁻¹²
```

LSMR 本身：

```text
converged = True
iterations = 26
estimated condition ≈ 1.13 × 10³
```

full fixed-point：

```text
converged = True
iterations = 27
final χ ≈ 7.55 × 10⁻⁷
```

## 7.4 阶段结论

matrix-free residual-LS 与 dense 诊断形式在当前 benchmark 下保持了非常高的数值一致性。

结论：

```text
I-3：STRONG PASS
```

---

# 8. I-4：residual-LS 正式接入 PGD enrichment

I-3 通过以后，residual-LS 才从“诊断工具”进入正式生产路径。

---

## 8.1 I-4A：独立 residual-LS spatial solver

新增核心模块：

```text
latin/tower_residual_ls_spatial.py
```

其主要职责包括：

- 构造 matrix-free weighted least-squares operator；
- 构造右端项；
- 调用 LSMR；
- 返回 spatial plastic strain mode；
- 通过塔筒平衡算子生成 spatial stress mode；
- 输出收敛状态和 LSMR 诊断信息。

同时保持：

```text
不显式构造 dense Aσ
不显式构造 normal equations
```

这一步为后续真实塔筒规模扩展提供了必要基础。

---

## 8.2 I-4B：接入 enrichment strategy

在：

```text
latin/tower_pgd_enrichment.py
```

中保留两条空间更新路径：

```text
paper_galerkin
residual_ls
```

并增加：

```python
spatial_strategy="paper_galerkin"
```

默认仍然是 `paper_galerkin`。

这意味着：

> residual-LS 在 I-4 阶段被接入，但没有偷偷替换原有默认算法。

---

## 8.3 第四模态生产路径 A/B 对比

在完全相同的已知第四模态 benchmark 上：

### A. paper_galerkin

```text
accepted = False
failure reason = fixed_point_not_converged
fixed-point iterations = 200
final χ ≈ 0.603487
```

即成功复现此前失稳。

### B. residual_ls

```text
accepted = True
failure reason = None
fixed-point converged = True
fixed-point iterations = 27
final χ ≈ 7.55 × 10⁻⁷
rank before = 3
rank after = 4
```

完整机械残差：

```text
before ≈ 4.25218 × 10⁻³
after  ≈ 3.04559 × 10⁻³
benefit ≈ 28.38%
```

同时：

```text
orthogonality error ≈ 5.83 × 10⁻¹⁶
equilibrium relative residual ≈ 7.96 × 10⁻¹⁶
```

并且事务审计确认：

```text
Trial-A basis unchanged = True
Trial-A residual unchanged = True
```

说明 enrichment 的接受或失败不会污染已经存在的 Trial-A。

## 8.4 阶段结论

```text
I-4：STRONG PASS
```

对应 checkpoint：

```text
ac44065
feat: integrate residual-LS into tower PGD enrichment
```

---

# 9. I-5：完整 outer LATIN-PGD 收敛验证

I-4 只能证明：

> 第四模态本身可以通过 residual-LS 富集。

但还不能证明：

> 完整 outer LATIN 能继续向后运行并最终收敛。

因此 I-5 将 `spatial_strategy` 正式暴露到：

```text
solve_tower_latin_pgd(...)
```

---

## 9.1 outer solver strategy plumbing

新增参数：

```python
spatial_strategy="paper_galerkin"
```

并从：

```text
solve_tower_latin_pgd
```

传递到：

```text
enrich_tower_pgd_basis_once
```

仍然保持：

```text
默认 = paper_galerkin
```

只有显式指定：

```python
spatial_strategy="residual_ls"
```

才使用 residual-LS。

---

## 9.2 solver-level regression

为避免只依靠真实 benchmark 判断接口正确性，还增加了单元测试验证：

- residual_ls 是否真正被传递到 enrichment；
- 非法 strategy 是否在 outer solver 层被拒绝；
- 原 Trial-A / Trial-B transaction 行为是否保持。

最终相关 regression suite：

```text
Ran 23 tests

OK
```

---

## 9.3 完整 nonlinear reversed-loading benchmark

采用此前已知 benchmark：

```text
水平荷载幅值：1.0 MN
周期：10
循环数：1
每循环增量：40
时间节点 Nt：41
梁单元：10
每单元 Gauss 点：2
每 Gauss 点环向纤维：16
总材料点 Nq：320
```

材料进入明显非线性区：

```text
elastic-init max stress ≈ 116.53 MPa
yield stress ≈ 80 MPa
```

因此这不是一个仅处于弹性阶段的 trivial test。

---

## 9.4 outer LATIN 最终结果

显式运行：

```python
spatial_strategy="residual_ls"
```

得到：

```text
termination reason          = converged
solver converged            = True
failure reason              = None
committed iterations        = 12
attempted iterations        = 12
trial evaluations           = 19
final PGD rank              = 7
total modes added           = 7
```

最终 LATIN indicator：

```text
ξfinal = 9.424149148753 × 10⁻⁵
```

第一次 committed indicator：

```text
ξ1 = 1.596571271 × 10⁻²
```

两者比值：

```text
ξfinal / ξ1 = 5.902742534452 × 10⁻³
```

即最终 committed indicator 约为第一次 committed indicator 的：

```text
0.5903%
```

---

## 9.5 committed history 的意义

12 次 committed iterations 中：

```text
committed decreases = 11
committed increases = 0
committed unchanged = 0
```

也就是说：

> 每一次真正被接受并写入 persistent baseline 的 LATIN state，其 indicator 都比前一次更小。

这比“最终勉强低于阈值”更有意义，因为它说明：

- outer LATIN 没有出现新的周期振荡；
- residual-LS 没有破坏 transaction；
- rank growth 没有造成 indicator 反复回升；
- accepted state 的全局趋势稳定下降。

---

## 9.6 A/B transaction 历史

commit kinds：

```text
B A A B A B A B B A B B
```

说明完整求解过程中既有：

```text
Trial-A commit
```

也有：

```text
Trial-B enrichment commit
```

并不是通过“关闭 enrichment”才收敛。

modes added per commit：

```text
1 0 0 1 0 1 0 1 1 0 1 1
```

最终：

```text
rank = 7
```

说明 residual-LS enrichment 在 outer LATIN 中被多次实际调用，而不是只在第四模态调用一次。

---

## 9.7 最后一次 enrichment 质量

最后一次 enrichment：

```text
accepted = True
failure reason = None
fixed-point converged = True
fixed-point iterations = 16
final χ ≈ 8.05 × 10⁻⁷
residual benefit ≈ 28.68%
orthogonality error ≈ 7.83 × 10⁻¹⁶
```

这说明即使已经进入较高 rank，residual-LS 仍然没有明显退化。

---

## 9.8 I-5 结论

```text
I-5：STRONG PASS
```

对应 checkpoint：

```text
369d424
feat: integrate residual-LS into tower LATIN-PGD outer solver
```

远端完整 SHA：

```text
369d424ab4a16d3d1e2f0b60973d73a10215169d
```

---

# 10. I-6：进入 FOM 前的联合闭环

I-6 不再继续修改算法，而是回答：

> 当前内部验证链条是否已经完整到足以进入 FOM？

---

## 10.1 I-6A 联合 regression

联合运行：

```text
tests.test_tower_equilibrium_operator
tests.test_tower_equilibrium_operator_integration
tests.test_tower_residual_ls_spatial
tests.test_tower_pgd_enrichment
tests.test_tower_pgd_enrichment_residual_ls
tests.test_tower_latin_pgd_solver
tests.test_tower_latin_pgd_solver_integration
```

结果：

```text
Ran 31 tests in 0.371s

OK
```

覆盖：

- equilibrium；
- real-tower equilibrium integration；
- residual-LS matrix-free solver；
- PGD enrichment；
- residual-LS enrichment；
- outer solver transaction；
- outer solver integration。

---

# 11. FOM entry gate

进入 FOM 之前，定义七项内部门槛。

| 编号 | 门槛 | 具体要求 | 当前状态 |
|---|---|---|---|
| G1 | 平衡性 | 齐次修正应力模态满足自由自由度平衡 | PASS |
| G2 | inner fixed-point | residual-LS 富集满足 complete-pair fixed-point 收敛 | PASS |
| G3 | 模态残差收益 | accepted rank growth 带来正的完整机械残差下降 | PASS |
| G4 | basis 质量 | 新模态具有 novelty，且正交化后场保持一致 | PASS |
| G5 | transaction | enrichment 失败或拒绝不会污染 Trial-A 或 persistent baseline | PASS |
| G6 | outer LATIN | nonlinear benchmark 达到 LATIN 收敛阈值 | PASS |
| G7 | 总回归 | equilibrium、spatial、enrichment、outer solver 联合测试通过 | PASS |

因此当前判定：

```text
内部 FOM 前置门槛：PASS
```

更准确地说：

> 当前塔筒 LATIN-PGD 已经通过“内部算法一致性资格检查”，可以进入 LATIN-PGD 与 FOM 的直接精度验证。

但这绝不等于：

```text
LATIN-PGD 已被 FOM 验证
```

---

# 12. 当前能够下结论的内容

根据 I-0 至 I-6，目前可以较有把握地说：

1. 塔筒平衡算子没有发现明显代数错误；
2. 第四模态 failure 不是由简单的自由自由度失衡造成；
3. 第四模态 failure 不是单纯由 basis orthogonality 不足造成；
4. paper-Galerkin 空间更新与当前离散时间更新组合在第四模态产生稳定 period-3；
5. residual-LS 空间更新能够解除该 period-3；
6. matrix-free residual-LS 与 dense diagnostic residual-LS 在当前 benchmark 下数值等价；
7. residual-LS 可被安全接入 PGD enrichment；
8. enrichment 的接受和拒绝不会破坏 Trial-A transaction；
9. residual-LS 可以驱动 PGD rank 从 3 增长到 4，并继续增长到 7；
10. 完整 outer LATIN-PGD 在 nonlinear reversed-loading benchmark 上达到绝对收敛阈值；
11. committed LATIN indicator 全程单调下降；
12. 当前 31 个内部 regression tests 全部通过。

---

# 13. 当前还不能下结论的内容

以下问题仍然没有被证明：

1. LATIN-PGD 与 FOM 的位移误差是否足够小；
2. LATIN-PGD 与 FOM 的应力误差是否足够小；
3. 塑性应变和损伤变量是否能够准确重现；
4. hysteresis loop 面积是否一致；
5. PGD rank = 7 是否已经足够；
6. 当前 Nt = 41 是否已经达到时间离散收敛；
7. 10 个梁单元是否已经达到网格收敛；
8. 16 个环向纤维是否已经达到截面离散收敛；
9. residual-LS 的梯形节点求积是否对结果存在敏感性；
10. 一循环 benchmark 的成功能否推广到多循环；
11. 多循环能否进一步推广到高周疲劳；
12. 当前方法在真实海上风机疲劳荷载谱下的计算效率优势究竟有多大。

因此下一阶段必须转向 FOM，而不是继续仅靠 LATIN-PGD 自身的内部指标证明正确性。

---

# 14. FOM-1：第一阶段直接精度验证

第一阶段 FOM 对比必须严格控制变量。

LATIN-PGD 与 FOM 必须使用完全相同的：

- 塔筒几何；
- 梁单元数量；
- 单元节点位置；
- Gauss 积分点；
- 环向纤维数量；
- 径向纤维数量；
- 材料参数；
- 初始状态；
- 荷载历史；
- 时间节点；
- 边界条件。

只有这样，二者差异才能主要解释为：

```text
PGD / LATIN 降阶与迭代误差
```

而不是离散模型本身不同。

---

# 15. FOM-1 建议优先比较的响应量

## 15.1 全局响应

优先：

- 塔顶水平位移时程；
- 塔底反力时程；
- 塔底弯矩时程；
- 全局 load-displacement hysteresis。

## 15.2 截面响应

建议选择：

- 塔底截面；
- 一个中部截面；
- 一个较高位置截面。

比较：

- 轴力；
- 弯矩；
- 截面曲率；
- 关键 Gauss 点响应。

## 15.3 材料点响应

至少选择：

- 最大拉应力纤维；
- 最大压应力纤维；
- 代表性中性轴附近纤维；
- 损伤最明显的材料点。

比较：

- stress history；
- plastic strain history；
- hardening / internal variable；
- damage history。

---

# 16. FOM-2：定量误差指标

不能只依靠曲线“看起来接近”。

## 16.1 位移相对误差

建议定义：

```text
eu = ||uPGD - uFOM|| / ||uFOM||
```

## 16.2 应力相对误差

```text
eσ = ||σPGD - σFOM|| / ||σFOM||
```

## 16.3 塑性应变误差

```text
ep = ||pPGD - pFOM|| / ||pFOM||
```

## 16.4 损伤误差

```text
ed = ||dPGD - dFOM|| / ||dFOM||
```

若 FOM 对应量接近零，则需要避免直接用相对误差，应增加绝对误差或加小正则项。

## 16.5 峰值误差

例如：

```text
epeak = |ypeak,PGD - ypeak,FOM| / |ypeak,FOM|
```

## 16.6 滞回能误差

对一个完整循环：

```text
Ehys = ∮ F du
```

再比较：

```text
eEhys = |Ehys,PGD - Ehys,FOM| / |Ehys,FOM|
```

这一指标对于循环塑性和疲劳问题尤其重要，因为即使峰值接近，滞回耗能也可能不同。

---

# 17. FOM-3：PGD rank 收敛

FOM 对比完成后，需要专门检查：

```text
rank = 1
rank = 2
rank = 3
...
```

对应误差如何变化。

核心问题是：

> 当 rank 增加时，PGD 误差是否系统下降？

如果 rank 增加后误差不再下降，则可能说明：

- PGD truncation 已经不是主要误差源；
- 时间离散误差占主导；
- 空间有限元误差占主导；
- residual-LS 离散误差占主导；
- outer LATIN tolerance 已成为主要限制。

---

# 18. FOM-4：离散化收敛

在确认 LATIN-PGD 与 FOM 基本一致后，再分别进行离散化研究。

## 18.1 时间步长

例如逐步比较：

```text
40 increments / cycle
80 increments / cycle
160 increments / cycle
```

观察：

- FOM 本身是否收敛；
- LATIN-PGD 是否同步收敛；
- residual-LS 结果对时间网格是否敏感。

## 18.2 梁单元网格

例如：

```text
10 elements
20 elements
40 elements
```

## 18.3 Gauss 点

例如：

```text
2 Gauss points
3 Gauss points
4 Gauss points
```

## 18.4 环向纤维

例如：

```text
16 fibers
32 fibers
64 fibers
```

## 18.5 residual-LS 求积

需要专门检查：

- 当前梯形节点求积；
- 是否包含 t₀；
- 其他合理时间求积方式。

这一步将直接回答：

> 当前 residual-LS 成功是否对某一特定离散组合高度敏感？

---

# 19. FOM-5：向多循环和高周疲劳推进

只有在单循环 nonlinear benchmark 通过 FOM 以后，才建议继续：

1. 5 cycles；
2. 10 cycles；
3. 20 cycles；
4. 100 cycles；
5. 更长疲劳时间窗；
6. 最终进入高周疲劳。

原因是：

> 如果一循环状态变量都不能与 FOM 保持一致，那么增加循环数只会累积误差，而不能证明方法正确。

---

# 20. 对当前第四模态问题的最终认识

经过 I-1 至 I-5 的逐层排查，现在对第四模态问题可以给出比最初更明确的判断。

## 20.1 已排除的主要解释

当前证据不支持把主要原因归结为：

- 最大迭代次数不足；
- 单纯 basis orthogonality 不足；
- 简单 constant damping 不足；
- 空间应力模态不满足自由自由度平衡；
- 新模态完全没有残差收益；
- rank 已经物理饱和。

## 20.2 当前最合理解释

目前最合理的解释是：

> 原论文导出的空间 Galerkin half-step 与当前项目采用的离散 temporal half-step，在第四模态问题上没有表现出足够稳定的离散固定点兼容性，因此形成了 period-3 轨道。

但这一点仍应称为：

```text
当前数值证据支持的工作假设
```

而不是严格数学定理。

## 20.3 residual-LS 的作用

residual-LS 不再通过原空间 Galerkin 投影关系间接更新空间模态，而是直接最小化当前离散机械残差。

因此它对当前 discrete BE residual 具有更直接的下降目标。

当前 benchmark 中的事实是：

```text
paper_galerkin → period-3
residual_ls    → fixed-point convergence
```

并且 residual-LS 进一步成功驱动完整 outer LATIN 收敛。

---

# 21. 当前阶段最终判定

综合 I-0 至 I-6：

```text
平衡算子             PASS
residual-LS 有效性   PASS
matrix-free 等价性   STRONG PASS
PGD enrichment       STRONG PASS
outer LATIN          STRONG PASS
31-test closure       PASS
```

因此：

```text
========================================
内部 FOM 前置资格：PASS
========================================
```

这意味着当前研究重点应正式从：

```text
“如何让塔筒 LATIN-PGD 跑起来”
```

转向：

```text
“塔筒 LATIN-PGD 相对于 FOM 到底有多准确，以及这种精度和效率如何随 rank 与离散参数变化”
```

这也是下一阶段最有科研价值的问题。

---

# 22. 下一阶段建议的严格顺序

建议继续保持目前“一步一验证”的原则：

```text
FOM-0：建立与 LATIN-PGD 完全一致的 FOM 基线
↓
FOM-1：单循环同网格直接响应对比
↓
FOM-2：建立位移、应力、塑性、损伤和滞回能误差指标
↓
FOM-3：研究 PGD rank 收敛
↓
FOM-4：研究时间 / 网格 / 纤维离散收敛
↓
FOM-5：推进多循环
↓
FOM-6：推进高周疲劳与真实海上风机载荷
```

不建议跳过单循环 FOM 验证，直接进入高周疲劳。

---

# 23. 当前项目状态一句话总结

> 当前海上风机塔筒 LATIN-PGD 已经完成从平衡算子、residual-LS 空间富集、matrix-free 实现、PGD 增秩、Trial-A / Trial-B transaction 到完整 outer LATIN 收敛的内部闭环；下一阶段应正式进入与对应全阶纤维梁柱模型 FOM 的定量精度和收敛性验证。
