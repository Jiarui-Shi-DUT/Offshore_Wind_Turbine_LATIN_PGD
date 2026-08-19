# Tower LATIN-PGD 迁移核心困难阶段总结：空间—时间子问题的离散一致性

**日期：2026-08-19**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**仓库：`Jiarui-Shi-DUT/Offshore_Wind_Turbine_LATIN_PGD`**  
**分支：`feature/offshore-wind-turbine-tower-fatigue`**  
**当前理论主线：Bhattacharyya et al. 原论文 $x-t$ LATIN-PGD → fiber beam-column offshore wind turbine tower**  
**本阶段性质：概念澄清与问题重新聚焦，不修改 `latin/` 核心代码。**  
**前序阶段 checkpoint：`f3a113f` — `test: isolate tower fourth-mode spatial half-step mismatch`**

---

# 1. 本阶段为什么需要重新总结

在将 LATIN-PGD 方法由一维三材料杆案例迁移至海上风机塔筒的过程中，已经完成了较多理论推导、塔筒结构离散、PGD enrichment、fixed-point diagnostics 和数值排查。

随着诊断逐步深入，问题已经从最初较宽泛的：

```text
“为什么一维杆可以收敛，而塔筒第 4 个 PGD mode 不收敛？”
```

逐渐收缩为一个更明确的核心问题：

```text
“在引入塔筒结构平衡约束后，
PGD 的 spatial half-step 与 temporal half-step
是否仍然处在同一套彼此一致的离散求解逻辑中？”
```

因此，本阶段的目的不是继续增加新的公式或新的 numerical probe，而是把目前已经得到的结论重新整理成一个直观、稳定、可作为后续研究基线的认识框架。

---

# 2. 最核心结论

目前整个 tower LATIN-PGD 迁移过程中，最核心的困难可以概括为：

> **真正困难的不是把一维杆替换成风机塔筒，而是在引入塔筒结构平衡约束后，如何仍然保持 PGD spatial subproblem 与 temporal subproblem 之间的数值一致性，使其 alternating iteration 能够稳定收敛。**

更具体地说：

```text
validated 1D implementation
│
├─ spatial half-step:
│  direct weighted residual minimisation
│
└─ temporal half-step:
   backward-Euler weighted residual minimisation

→ 两个 half-step 数值上彼此协调
→ enrichment fixed point 可以稳定收敛
```

而当前 tower v1 为了更忠实原论文，引入了：

```text
tower implementation
│
├─ spatial half-step:
│  paper-derived Eq. (70)-(71) Galerkin formulation
│
└─ temporal half-step:
   validated 1D backward-Euler residual minimisation

→ 两个 half-step 各自都有理论/数值依据
→ 但它们组合后的离散一致性尚未被证明
→ 第 4 PGD mode 出现稳定 period-3 orbit
```

因此目前问题的本质不是：

```text
“tower finite element 不工作”
```

也不是：

```text
“LATIN-PGD 本身无法用于塔筒”
```

而是：

```text
“当前 tower spatial update 与 temporal update
在离散层面可能不是一对彼此匹配的 alternating operators。”
```

---

# 3. 用最直观的方式理解 LATIN-PGD enrichment

在 PGD enrichment 中，每增加一个新的 separated mode，本质上都需要反复确定两个东西：

```text
空间模态 p(x)
+
时间函数 lambda(t)
```

可以把它理解为两个相互依赖的子问题：

```text
给定空间模态 p
        ↓
求时间函数 lambda
        ↓
给定新的 lambda
        ↓
重新求空间模态 p
        ↓
再重新求 lambda
        ↓
……
直到二者互相一致
```

理想情况下：

```text
p^(0)
  ↓
lambda^(0)
  ↓
p^(1)
  ↓
lambda^(1)
  ↓
p^(2)
  ↓
lambda^(2)
  ↓
……
  ↓
(p*, lambda*)
```

最终得到一个 ordinary fixed point，即 complete-pair change：

$$ \chi_{\rm fp}^{(k)}\rightarrow0 $$

---

# 4. 为什么 validated 1D 三材料杆能够工作

当前已经成功复现的一维三材料杆 LATIN-PGD implementation 中，spatial half-step 和 temporal half-step 在数值上采用了高度一致的思想。

## 4.1 1D spatial half-step

固定时间函数：

$$ \lambda_n,\qquad\dot\lambda_n $$

构造 mechanical residual：

$$ r_n(p)=\dot\lambda_n p-D_{H,n}\lambda_nA_\sigma p+\Delta_n $$

然后直接求：

$$ p=\operatorname*{arg\,min}_pJ_h(p) $$

其中：

$$ J_h(p)=\frac12\sum_nw_n r_n(p)^TMD_{H,n}^{-1}r_n(p) $$

也就是说：

> **空间模态的目标是让当前机械残差尽可能小。**

## 4.2 1D temporal half-step

固定 spatial pair：

$$ p,\qquad s=A_\sigma p $$

current validated 1D implementation 使用 backward-Euler temporal update。

对于 $n\ge1$：

$$ \dot\lambda_n=\frac{\lambda_n-\lambda_{n-1}}{\Delta t_n} $$

定义：

$$ g_n=\frac{p}{\Delta t_n}-D_{H,n}s $$

以及：

$$ b_n=\frac{p}{\Delta t_n}\lambda_{n-1}-\Delta_n $$

则当前时间步 residual 为：

$$ r_n=g_n\lambda_n-b_n $$

并选择：

$$ \lambda_n=\operatorname*{arg\,min}_{\lambda_n}\frac12r_n^TMD_{H,n}^{-1}r_n $$

也就是说：

> **时间函数的目标同样是让当前机械残差尽可能小。**

---

# 5. 为什么 1D 的 spatial 与 temporal half-step 容易彼此协调

虽然 current temporal solve 是 causal backward-Euler one-step minimisation，而不是严格的 whole-time global minimisation，但 spatial 和 temporal 两个 half-step 的数值目标具有高度一致性：

```text
spatial:
reduce weighted mechanical residual

temporal:
reduce weighted mechanical residual
```

因此可以直观理解为：

```text
两个人虽然分别负责空间和时间，
但都在用相同的“误差标准”修正自己的答案。
```

所以：

```text
空间修正
↓
时间修正
↓
空间修正
↓
时间修正
```

能够较自然地朝同一个方向靠近。

这也是 current validated 1D enrichment 能稳定产生新 PGD modes 的重要数值背景之一。

---

# 6. 从 1D 迁移到 tower 后发生了什么变化

风机塔筒与一维杆相比，需要显式处理结构平衡和兼容条件。

对于 tower fiber beam-column discretisation：

```text
element
↓
Gauss point
↓
section fiber
↓
material point q=(e,g,f)
```

需要建立：

- material-point metric $M$；
- compatibility operator $H$；
- reference modulus $C_0$；
- structural equilibrium projection；
- free structural DOFs；
- fiber section integration；
- tower-level admissibility。

因此 spatial problem 已经不再只是简单的一维材料点向量问题。

---

# 7. Tower structural equilibrium operator

当前 tower reference equilibrium mapping 为：

$$ \mathcal E_{\rm tower}=H(H^TMC_0H)^{-1}H^TMC_0 $$

给定 source strain：

$$ p $$

compatible strain 为：

$$ \varepsilon=\mathcal E_{\rm tower}p $$

stress correction 为：

$$ s=C_0(\varepsilon-p) $$

并满足：

$$ H^TMs=0 $$

因此 tower spatial mode 必须同时满足：

```text
material-point correction
+
fiber compatibility
+
global structural equilibrium
```

这一扩展是从 1D 到 tower 的必要结构变化，但目前没有证据说明 reference equilibrium operator 本身就是 period-3 的错误来源。

---

# 8. 为了更忠实原论文，我们改变了 spatial half-step

Bhattacharyya et al. 原论文在 enrichment 中采用 hybrid strategy：

```text
spatial functions
→ Galerkin formulation

time functions
→ residual minimisation
```

因此在 tower v1 推导中，我们没有继续直接照搬 1D spatial residual least-squares，而是尽量回到原论文 Eq. (65)-(71)。

这导致 tower spatial half-step 变成：

```text
fixed lambda
↓
construct temporal contractions
↓
construct W
↓
solve structural Galerkin equilibrium
↓
recover spatial plastic mode p
```

---

# 9. Tower Eq. (70)-(71) spatial half-step

固定：

$$ \lambda(t),\qquad\dot\lambda(t) $$

current tower implementation 定义：

$$ a_h=\sum_{n=1}^{N}\Delta t_n\lambda_n\dot\lambda_n $$

$$ A_q=\sum_{n=1}^{N}\Delta t_nH_{\sigma,nq}\lambda_n^2 $$

$$ \bar\delta_q=\sum_{n=1}^{N}\Delta t_n\Delta_{nq}\lambda_n $$

并构造：

$$ W_q^{-1}=A_q+\frac{a_h}{C_{0,q}} $$

然后求解：

$$ H^TMD_WH\bar{\tilde U}=-H^TMD_W\bar\delta $$

得到：

$$ \bar{\tilde\varepsilon}=H\bar{\tilde U} $$

以及：

$$ \bar\sigma=D_W(\bar{\tilde\varepsilon}+\bar\delta) $$

最终恢复：

$$ p=\frac{\bar{\tilde\varepsilon}}{a_h}-C_0^{-1}\bar\sigma $$

这就是 current tower Eq. (70)-(71) spatial half-step。

---

# 10. Eq. (70)-(71) 目前并没有被证明是错的

经过重新审计，目前已经确认：

- temporal Galerkin projection 的结构与原论文一致；
- $\lambda^2H_\sigma$ contraction 一致；
- $\lambda\dot\lambda$ contraction 一致；
- $\lambda\Delta$ contraction 一致；
- $W^{-1}$ 结构一致；
- Eq. (70) FE assembly 一致；
- RHS sign 一致；
- Eq. (71) plastic spatial mode recovery 一致；
- normalization 后重新通过 equilibrium operator 计算 stress 在数学上等价。

因此当前应该明确：

$$ \boxed{\text{没有证据表明 Eq. (70)-(71) 的理论推导或 tower FE implementation 本身存在明显代数错误}} $$

这点非常重要。

---

# 11. 真正发生变化的是“空间—时间求解组合”

current tower v1 采用：

```text
spatial:
paper-derived Galerkin Eq. (70)-(71)

temporal:
validated 1D backward-Euler residual minimisation
```

也就是：

$$ \boxed{\text{paper-style spatial Galerkin}+\text{project-BE temporal minimisation}} $$

而 validated 1D 更接近：

$$ \boxed{\text{residual-minimising spatial solve}+\text{project-BE temporal minimisation}} $$

因此，从一维杆迁移到塔筒过程中，真正改变的并不只是结构维度，而是：

> **spatial half-step 的 variational formulation 发生了改变。**

---

# 12. 为什么这种变化可能导致不收敛

可以把两个 half-step 想象成两个合作的人。

在 validated 1D 中：

```text
A：我根据 residual minimisation 求最好的空间模态。
B：我也根据 residual minimisation 求最好的时间函数。
```

两人的判定标准相近，因此容易逐步趋向共同答案。

而 current tower 中：

```text
A：我根据 Eq. (70)-(71) Galerkin 条件求空间模态。
B：我根据 BE residual minimisation 求时间函数。
```

这两个规则各自都有依据，但是：

> **它们是否在 current discrete formulation 下具有同一个稳定 fixed point，目前尚未被证明。**

因此可能出现：

```text
A认为当前最合适的是状态 1
↓
B修正后得到状态 2
↓
A再修正得到状态 3
↓
B再修正又把系统推回状态 1
```

于是形成：

```text
state A
↓
state B
↓
state C
↓
state A
↓
state B
↓
state C
```

这正是当前第 4 mode 中观察到的 period-3 orbit 的直观含义。

---

# 13. 第 4 mode 实际发生了什么

当前 reversed tower benchmark 中：

```text
10 elements
× 2 Gauss points
× 16 fibers
= 320 material points
```

在第 4 个 enrichment mode 上，current tower raw fixed-point map 没有趋向普通 fixed point，而是：

$$ \chi^{(k)}\approx0.577,\ 0.658,\ 0.603,\ 0.577,\ 0.658,\ 0.603,\ldots $$

同时：

$$ d(z^k,z^{k-3})\rightarrow O(10^{-6}) $$

经过 in-loop orthogonalisation 后甚至：

$$ d(z^k,z^{k-3})\rightarrow O(10^{-8}) $$

所以它不是随机振荡，而是一个非常稳定的 period-3 orbit。

---

# 14. 为什么 period-3 比“收敛慢”更严重

如果 fixed-point iteration 只是收敛很慢，那么增加迭代次数理论上仍可能解决问题。

但 period-3 不一样。

period-3 意味着：

```text
iteration 1
→ iteration 2
→ iteration 3
→ iteration 1
```

也就是说：

> **current alternating map 本身正在把解送进一个循环，而不是送向一个 fixed point。**

因此简单增加：

```text
120 iterations
→ 200
→ 400
```

不会从根本上解决问题。

---

# 15. 到目前为止已经排除了哪些可能原因

围绕第 4 mode，我们已经逐项检查了多个可能原因。

## 15.1 iteration cap

增加 fixed-point iterations 后仍然进入稳定 period-3。

因此：

```text
iteration cap
×
```

## 15.2 constant under-relaxation

测试多组 spatial under-relaxation 后，虽然 relaxed step 可以被人为压小，但 raw fixed-point defect 没有消失。

因此：

```text
simple damping
×
```

## 15.3 sequential BE temporal solve

whole-time coupled BE diagnostic 仍然保持明显周期行为。

因此：

```text
sequential BE scope
×
```

至少不是单独主因。

## 15.4 mode insignificance / basis saturation

对 period-3 三个 phase 进行 diagnostic post-processing 后发现：

- spatial novelty 很高；
- temporal significance 明显；
- residual reduction 仍约为 $10\%$–$19\%$。

因此：

```text
basis saturation
×

fourth mode insignificant
×
```

## 15.5 seed sensitivity

对 top-10 residual-energy rows 作为 deterministic initial seed 进行了 sweep。

结果：

```text
all 10 seeds
→ nonconvergent
→ collapse to same three phase families
```

因此，在当前 residual-driven seed family 中：

```text
seed / basin sensitivity
×
```

不是主要解释。

## 15.6 in-loop orthogonalisation

加入 validated 1D-style in-loop orthogonalisation 后：

$$ \max_j|\langle p,p_j\rangle_M|\approx5.55\times10^{-16} $$

说明 candidate 已经在 machine precision 下与 existing basis 正交。

但仍然：

```text
converged = False
```

并继续保持 period-3。

因此：

```text
lack of in-loop orthogonalisation
×
```

也不是主因。

---

# 16. 第一次真正解除 period-3 的变化

随后进行了关键 A/B/C diagnostic。

## A：current tower map

```text
tower Eq.(70)-(71)
+
no in-loop orthogonalisation
+
current BE temporal solve
```

结果：

```text
period-3
nonconvergent
```

## B：加入 1D-style in-loop orthogonalisation

```text
tower Eq.(70)-(71)
+
1D-style in-loop orthogonalisation
+
current BE temporal solve
```

结果：

```text
period-3
nonconvergent
```

## C：将 spatial half-step 替换为 literal 1D-style direct weighted LS

```text
literal 1D-style direct weighted-LS spatial solve
+
same in-loop orthogonalisation
+
same BE temporal solve
```

结果：

```text
converged = True
iterations = 27
```

最终：

$$ \chi_{\rm fp}=7.5524\times10^{-7}\lt10^{-6} $$

最后若干步为：

```text
8.5866e-05
5.5837e-05
3.6310e-05
2.3613e-05
1.5355e-05
9.9860e-06
6.4940e-06
4.2230e-06
2.7460e-06
1.7860e-06
1.1610e-06
7.5524e-07
```

表现出清楚的 contraction。

---

# 17. 为什么这个结果如此关键

这是 current diagnostics 中：

> **第一个被隔离出来、能够把 period-3 改成 ordinary fixed-point convergence 的算法差异。**

因此当前诊断链可以写成：

```text
iteration cap
    ×

simple damping
    ×

sequential-BE scope
    ×

mode insignificance
    ×

basis saturation
    ×

tested seed family
    ×

lack of in-loop orthogonalisation
    ×

spatial half-step formulation
    ✓
```

因此 spatial half-step formulation 成为当前最需要审计的核心。

---

# 18. 但是 Case C 不能简单理解为“修好了 Eq. (70)-(71)”

Case C 解的是一个不同的离散 variational problem。

它定义：

$$ B_n=\dot\lambda_nI-\lambda_nD_{H,n}A_\sigma $$

并最小化：

$$ J_h(p)=\frac12\sum_nw_n(B_np+\Delta_n)^TMD_{H,n}^{-1}(B_np+\Delta_n) $$

stationarity condition 为：

$$ \left(\sum_nw_nB_n^TMD_{H,n}^{-1}B_n\right)p=-\sum_nw_nB_n^TMD_{H,n}^{-1}\Delta_n $$

这就是 literal 1D-style direct weighted residual minimisation。

---

# 19. Case C 与 paper Eq. (70)-(71) 从哪里开始不同

展开 Case C 的 normal equation 会出现：

$$ \dot\lambda_n^2MD_{H,n}^{-1} $$

$$ -\lambda_n\dot\lambda_nMA_\sigma $$

$$ -\lambda_n\dot\lambda_nA_\sigma^TM $$

以及：

$$ \lambda_n^2A_\sigma^TMD_{H,n}A_\sigma $$

因此 Case C spatial equation 中包含：

```text
dot(lambda)^2 H_sigma^-1
lambda dot(lambda)
lambda^2 H_sigma
```

而原论文 spatial Galerkin equation 是通过 kinematic admissibility 与 temporal Galerkin projection 构造，只自然产生：

$$ \langle H_\sigma\lambda^2\rangle $$

$$ \langle\lambda\dot\lambda\rangle $$

$$ \langle\Delta\lambda\rangle $$

不会产生：

$$ \langle\dot\lambda^2H_\sigma^{-1}\rangle $$

所以：

$$ \boxed{\text{B 与 C 是不同 variational formulation，不是同一个方程的两种写法}} $$

---

# 20. 当前真正的核心矛盾

原论文 continuous formulation 是：

```text
spatial Galerkin
+
temporal minimisation
```

两者属于同一套 LATIN-PGD 构造。

但是 current tower v1 实际采用：

```text
paper-derived spatial Galerkin
+
project-specific causal BE temporal update
```

而 current temporal side 并不是 paper-exact DG0 implementation。

因此目前尚未证明：

> **paper spatial Galerkin half-step 与 project BE temporal half-step 在离散层面仍然保持原有的一致性。**

这就是当前真正需要解决的问题。

---

# 21. 为什么不能简单把问题归结为“塔筒太复杂”

塔筒确实比一维杆复杂得多。

从：

```text
1D rod
```

变成：

```text
2D beam-column tower
```

还增加了：

```text
fiber section
Gauss points
many material points
global equilibrium
compatibility matrix H
material metric M
reference projection
trial / commit / revert
```

但是目前这些基础部分已经能够工作。

例如：

- tower FOM 已建立；
- fiber material-point layout 已建立；
- equilibrium operator 已建立；
- local stage 已建立；
- search direction 已建立；
- Trial-A / Trial-B 已建立；
- PGD basis transaction 已建立；
- 前 3 个 new modes 能够加入；
- current LATIN-PGD 可以推进到第 4 enrichment failure。

因此现在的瓶颈已经不是：

```text
结构复杂度本身
```

而是：

```text
结构复杂度引入以后，
spatial PGD subproblem 被重新定义，
从而打破了原 1D numerical pair 的协调关系。
```

---

# 22. 最通俗的一句话

如果把整个问题讲给没有看过所有公式的人，可以说：

> 一维杆里，空间和时间两个求解器像两个使用同一把尺子的人，所以可以越改越接近；到了塔筒里，为了忠实原论文，我们给空间求解器换了一把 Galerkin 的尺子，但时间求解器还在使用原来的一维 BE 残差尺子。现在第 4 个模态说明，这两把尺子组合后可能无法保证双方最终对同一个答案达成一致。

---

# 23. 当前 period-3 的直观含义

current tower 第 4 mode 不是：

```text
“算不出来”
```

而更接近：

```text
“空间子问题和时间子问题各自都能给出答案，
但双方反复交换答案以后无法停在一个共同 fixed point。”
```

即：

```text
A
↓
B
↓
C
↓
A
↓
B
↓
C
```

因此：

$$ \mathcal F^3(z)\approx z $$

但：

$$ \mathcal F(z)\ne z $$

这就是 period-3 fixed-point pathology 的最直观表达。

---

# 24. 当前最可能的工作假设

目前最有解释力的工作假设是：

> **current tower mixed discrete map 的 spatial Galerkin half-step 与 causal BE temporal half-step 并非同一个 discrete objective 或同一个 fully consistent discrete weak formulation 的两个 coordinate updates，因此不必具有 contraction 性质。**

必须强调：

```text
这是工作假设
不是最终数学证明
```

它得到以下事实支持：

```text
Eq.(70)-(71) + BE
→ period-3

direct weighted-LS spatial + same BE
→ ordinary convergence
```

但还需要进一步 discrete consistency audit。

---

# 25. 目前仍未完全隔离的另一个差异

Case C 与 tower Eq. (70)-(71) 不只是 variational formulation 不同。

它们的 temporal quadrature 也不同。

current tower Eq. (70)-(71) contractions 使用：

$$ \sum_{n=1}^{N}\Delta t_n(\cdot)_n $$

属于 current right-endpoint BE-consistent treatment。

而 diagnostic C 为了忠实 1D spatial implementation，使用：

```text
trapezoidal nodal quadrature
+
includes n=0
```

因此 B → C 同时改变了：

```text
A. spatial variational formulation
B. temporal quadrature / t0 treatment
```

所以：

> **目前不能把 C 的成功 100% 归因于 Galerkin → residual-LS 这一项。**

这也是下一阶段需要进一步收缩的问题。

---

# 26. 当前不应该做什么

在 current understanding 下，不应该立即：

- 删除 Eq. (70)-(71)；
- 直接宣称原论文 Galerkin formulation 不适用于 tower；
- 直接把 diagnostic C 写进 production core；
- 直接接受 unconverged period-3 phase；
- 放宽 fixed-point convergence gate；
- 使用 cycle-aware acceptance 绕过普通 fixed point；
- 引入 Anderson acceleration；
- 引入 Aitken acceleration；
- 引入大量 damping；
- 直接跳到 high-cycle tower loading；
- 引入 cycle-phase PGD；
- 引入 multi-time-scale PGD。

因为 current core theoretical issue 尚未闭合。

---

# 27. 下一阶段真正应该做什么

后续工作应集中于：

> **spatial-temporal discrete consistency audit**

第一步：

```text
把 current causal BE temporal Eq. (72)
写成明确的 discrete optimality / weak equation
```

第二步：

```text
把 current tower Eq. (65)-(71)
写成完全对应的 discrete Galerkin equation
```

第三步：

```text
检查二者是否来自同一套 discrete formulation
```

第四步：

```text
区分：
variational-form mismatch
vs
quadrature mismatch
vs
t0 treatment mismatch
```

---

# 28. 后续可能出现的两条路线

## 路线 A：paper fidelity first

如果审计证明：

```text
paper Eq.(70)-(71)
本身应保留
```

那么需要恢复与其更一致的 temporal discretisation。

这可能最终要求重新审查：

```text
paper DG0 temporal treatment
```

此时 DG0 将重新成为必要问题。

## 路线 B：tower v1 migration stability first

如果目标优先是：

```text
先把 validated 1D numerical strategy
稳定推广到 tower
```

那么可以考虑：

```text
generalise 1D-style direct weighted residual spatial solve
to tower material points
```

并明确标注：

```text
tower-v1 engineering discretisation
```

而不是声称：

```text
paper-exact Eq.(70)-(71)
```

随后再用 paper-Galerkin formulation 作为 fidelity comparison。

---

# 29. 当前最重要的 source-layer 边界

后续必须继续区分四个 source layers。

## Layer 1：原论文明确

```text
LATIN local/global structure
x-t PGD
spatial Galerkin formulation
temporal residual minimisation
Eq.(61)-(72)
post-basis orthogonalisation
```

## Layer 2：严格推导

```text
tower H / M formulation
Eq.(70)-(71) FE translation
tower structural equilibrium mapping
```

## Layer 3：validated 1D implementation

```text
backward-Euler temporal solve
direct weighted-LS spatial solve
in-loop orthogonalisation
tested convergence behaviour
```

## Layer 4：tower v1 engineering choices

```text
current BE temporal discretisation
complete-pair graph-norm criterion
transaction semantics
diagnostic probes
temporary algorithm comparisons
```

后续所有修改都必须明确说明属于哪一层。

---

# 30. 目前可以正式冻结的核心结论

**结论 1**  
从一维三材料杆迁移到 offshore wind turbine tower 的最核心困难，不是结构维度增加本身，而是引入 structural equilibrium 后 spatial PGD subproblem 的形式发生变化。

**结论 2**  
validated 1D implementation 中 spatial 与 temporal half-step 都以 weighted mechanical residual reduction 为核心，因此具有较好的 numerical compatibility。

**结论 3**  
current tower v1 使用 paper-derived Eq. (70)-(71) spatial Galerkin half-step，同时保留 validated 1D backward-Euler temporal residual minimisation。

**结论 4**  
这两个 half-step 各自有依据，但它们组合后的 discrete consistency 尚未被证明。

**结论 5**  
current tower 第 4 enrichment mode 出现稳定 period-3 orbit，说明 alternating map 并没有趋向 ordinary fixed point。

**结论 6**  
iteration cap、simple damping、whole-time BE scope、mode insignificance、basis saturation、tested seed family 和 lack of in-loop orthogonalisation 都已经不能解释这一 failure。

**结论 7**  
将 spatial half-step 暂时替换为 literal 1D-style direct weighted residual minimisation 后，在保持同一 BE temporal solve 的情况下，第 4 mode 在 27 iterations 内恢复 ordinary convergence。

**结论 8**  
因此 spatial-half-step formulation 是 current diagnostics 中第一个被隔离出来、能够解除 period-3 的关键算法差异。

**结论 9**  
当前没有证据证明 tower Eq. (70)-(71) 本身存在 algebraic implementation error。

**结论 10**  
当前最需要研究的是 paper-style spatial Galerkin 与 project-BE temporal update 的离散一致性，而不是继续扩大模型复杂度。

---

# 31. 博士论文式概括

可以将当前最核心困难正式表述为：

> **The central difficulty in extending the validated 1D LATIN-PGD implementation to the offshore wind-turbine tower is not the increase in structural dimensionality itself, but the preservation of numerical consistency between the spatial and temporal PGD subproblems after structural equilibrium is introduced.**

中文：

> **将已验证的一维 LATIN-PGD 方法推广至海上风机塔筒时，核心困难并非结构维度增加本身，而是在引入塔筒结构平衡约束后，如何保持 PGD 空间子问题与时间子问题之间的数值一致性，使其交替迭代能够稳定收敛。**

---

# 32. 当前研究路线的最简逻辑图

```text
validated 1D bar
│
├─ spatial:
│  weighted residual minimisation
│
└─ temporal:
   BE weighted residual minimisation
        │
        ↓
   stable alternating iteration
        │
        ↓
      convergence


tower migration
│
├─ spatial:
│  paper Eq.(70)-(71) Galerkin
│
└─ temporal:
   validated 1D BE update
        │
        ↓
   discrete compatibility unknown
        │
        ↓
   fourth-mode period-3


diagnostic C
│
├─ spatial:
│  1D-style direct weighted LS
│
└─ temporal:
   same BE update
        │
        ↓
   ordinary contraction
        │
        ↓
   convergence in 27 iterations
```

---

# 33. 下一阶段任务

下一阶段建议标题：

```text
Tower LATIN-PGD spatial-temporal discrete consistency audit:
paper Galerkin spatial update
vs
project backward-Euler temporal update
```

下一阶段只回答：

```text
current Eq.(70)-(71) discrete spatial Galerkin
与
current causal BE Eq.(72)
到底是否属于同一套一致的 discrete LATIN-PGD formulation
```

在此之前：

```text
do not modify core
do not bypass fixed-point convergence
do not accept period-3 phase
do not introduce higher-level ROM extensions
```

---

# 34. 本阶段定位

本阶段不是新的数值算法修改，而是对前序所有 diagnostics 的概念收敛。

它的主要价值是把此前多个局部问题重新统一为一个更清晰的研究主线：

$$ \boxed{\text{从“为什么第 4 mode 不收敛？”转变为“如何保持 spatial-temporal discrete consistency？”}} $$

这应作为后续 tower LATIN-PGD 研究的 current core problem。
