# Tower LATIN-PGD 迁移核心困难：通俗版阶段总结

**日期：2026-08-19**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**用途：用于快速理解当前研究到底卡在哪里，以及为什么下一步要做 spatial-temporal discrete consistency audit。**

**本版补充：新增“为什么 1D 当初采用 residual-LS 而不是 paper spatial Galerkin”的历史背景，并明确区分现有代码可确认事实与对当时工程选择原因的合理重建。**

---

# 1. 一句话先说清楚

目前从一维三材料杆迁移到海上风机塔筒时，最核心的困难不是：

```text
“塔筒太复杂，所以 LATIN-PGD 算不动。”
```

而是：

> **一维程序里，空间模态和时间函数的求解方式彼此很匹配；迁移到塔筒以后，为了更忠实原论文，我们改变了空间模态的求法，但时间函数仍然暂时沿用一维程序的 BE 求法，于是这两个 half-step 在离散层面可能不再完全匹配。**

这就是当前最核心的问题。

---

# 2. LATIN-PGD enrichment 到底在做什么

每增加一个新的 PGD mode，都要反复求两个东西：

```text
空间模态 p(x)
+
时间函数 lambda(t)
```

最直观地理解就是：

```text
先给一个空间形状
      ↓
根据这个空间形状求时间函数
      ↓
再根据新的时间函数修正空间形状
      ↓
再修正时间函数
      ↓
……
直到空间和时间彼此一致
```

如果一切正常，应该出现：

```text
空间 1
↓
时间 1
↓
空间 2
↓
时间 2
↓
空间 3
↓
时间 3
↓
……
↓
最终稳定答案
```

也就是 fixed-point convergence。

---

# 3. 为什么一维三材料杆案例能成功

在已经成功复现的一维程序里，空间和时间两个子问题虽然算法形式不完全一样，但目标很接近：

```text
空间求解：
找一个空间模态，让 mechanical residual 尽可能小

时间求解：
找一个时间函数，让 mechanical residual 尽可能小
```

可以把它理解成：

> **空间求解器和时间求解器都在使用同一种“误差尺子”。**

所以它们反复交换答案以后，通常能够逐渐靠近同一个结果。

因此一维程序表现为：

```text
空间
↓
时间
↓
空间
↓
时间
↓
逐渐稳定
```

---

# 4. 到了风机塔筒以后，什么东西变了

塔筒当然比一维杆复杂得多。

从一维杆变成塔筒以后，需要考虑：

```text
beam-column elements
Gauss points
section fibers
many material points
compatibility
global equilibrium
structural DOFs
fiber integration
```

也就是说，空间模态不再只是一个简单的材料点向量。

它必须满足：

```text
材料点响应
+
截面纤维兼容
+
整个塔筒的结构平衡
```

所以我们建立了：

```text
compatibility matrix H
material-point metric M
reference equilibrium operator
fiber beam-column tower
```

这些都是必须做的，而且目前基本已经建立成功。

---

# 5. 真正关键的改变发生在 spatial half-step

为了尽量忠实 Bhattacharyya 原论文，我们在塔筒里没有继续完全照搬一维程序的 spatial residual minimisation。

而是采用了原论文 Eq. (70)-(71) 的 Galerkin spatial formulation。

所以当前 tower v1 实际上是：

```text
空间：
paper Eq.(70)-(71) Galerkin

时间：
validated 1D backward-Euler residual minimisation
```

也就是：

```text
paper-style spatial solver
+
1D-style temporal solver
```

这两个方法各自都有依据。

真正的问题是：

> **它们组合在一起以后，是否还是一对彼此协调的 alternating solvers？**

目前第 4 mode 告诉我们：至少在当前 benchmark 下，不一定。

---

# 6. 为什么一维三材料杆当初没有直接使用 paper spatial Galerkin

这里有一个现在必须补充清楚的历史背景。

当初做一维三材料杆复现时，并不是因为我们认为原论文不应该使用 spatial Galerkin，也不是因为一维问题不能做 Galerkin。

原论文 enrichment 的基本路线本来就是：

```text
spatial functions
→ Galerkin formulation

temporal functions
→ minimisation
```

但真正落实到一维 `pgd_enrichment.py` 的 numerical implementation 时，我们采用的是：

```text
residual-driven weighted least-squares spatial solve
```

而不是 paper Eq. (65)-(71) 的严格 Galerkin spatial solve。

所以需要明确：

> **一维三材料杆的成功复现，在 enrichment spatial half-step 这一小块，并不是 paper-exact reproduction，而是一套 residual-driven engineering implementation。**

---

# 7. 为什么当时 residual-LS 很自然

一维杆的 equilibrium structure 非常简单。

给定 spatial plastic mode：

$$ p $$

其 equilibrated stress correction 可以直接写成：

$$ s=A_\sigma p $$

因此给定 temporal function 后，mechanical residual 很容易写成：

$$ r=Bp+\Delta $$

于是 spatial problem 可以非常直接地理解成：

```text
找一个 p
↓
让 mechanical residual 尽可能小
```

也就是：

$$ \min_p \|Bp+\Delta\|^2 $$

这样做的现实优点很明显：

```text
公式直接
代码容易实现
容易调试
容易检查 residual 是否下降
容易判断新 mode 是否真正有效
和 temporal residual minimisation 的思路很接近
```

因此，在当时“先把一维 LATIN-PGD 算例稳定跑通”的目标下，residual-LS 是一条非常务实的 numerical route。

---

# 8. 为什么这个差异在一维阶段没有暴露成大问题

一维三材料杆本身非常简单。

它只有轴向变形，global equilibrium 很强，spatial correction 的自由度也远低于 tower。

所以即使：

```text
paper:
spatial Galerkin

1D implementation:
residual-LS spatial solve
```

两种方法在这个简单 benchmark 上都可能给出非常有效的 spatial directions。

最终我们看到的是：

```text
位移响应合理
应力历史合理
损伤演化合理
PGD enrichment 能工作
LATIN 能收敛
```

因此，从 numerical reproduction 的角度，一维案例是成功的。

但如果现在提高标准，问：

```text
“一维 enrichment spatial half-step
是不是严格按 paper Eq.(65)-(71) 实现？”
```

答案应该是：

```text
不是。
```

更准确的说法是：

> **一维三材料杆成功验证的是 LATIN-PGD 的主要算法结构，以及一套可工作的 residual-LS spatial + BE temporal numerical implementation，而不是 enrichment spatial half-step 的逐式 paper-exact reproduction。**

---

# 9. 为什么这个历史差异到了 tower 才突然重要

迁移到 tower 时，我们提高了对 paper fidelity 的要求。

我们的基本想法是：

```text
既然要真正把 LATIN-PGD 推广到海上风机塔筒
↓
spatial part 应尽可能回到原论文 Eq.(65)-(71)
```

所以 tower spatial half-step 改成了：

```text
paper Eq.(70)-(71) Galerkin
```

但是 temporal half-step 因为 paper-exact DG0 的全部离散细节尚未完全恢复，我们暂时继续沿用：

```text
validated 1D backward-Euler temporal solve
```

于是实际发生的并不是简单的：

```text
1D bar
↓
换成 tower geometry
```

而是：

```text
1D implementation
│
├─ spatial:
│  residual-LS
│
└─ temporal:
   BE
        ↓

tower implementation
│
├─ spatial:
│  改回 paper Galerkin
│
└─ temporal:
   仍然沿用 BE
```

也就是说：

> **我们在迁移结构的同时，也改变了 PGD spatial algorithm。**

这一点在早期没有充分暴露出来。

---

# 10. 这也解释了为什么 Case C 特别有意义

Case C 做的事情，本质上是：

```text
tower geometry / tower equilibrium
保持不变

但是 spatial half-step
重新换回 1D-style residual-LS
```

同时 temporal solve 仍然保持原来的 BE。

于是：

```text
1D-style spatial
+
same BE temporal
→ 重新收敛
```

这说明：

> **一维程序真正被数值验证过的，是“residual-LS spatial + BE temporal”这一对 numerical half-steps。**

Case C 进一步证明：这对 numerical strategy 推广到 tower material-point system 后，至少对于当前第 4 mode，也可以恢复 ordinary convergence。

因此，一维工程实现现在不是一个需要否定的“旧方案”，反而是当前诊断 tower 问题时非常重要的 numerical reference。

---

# 11. 这是不是说明一维复现当初做错了

不能这样理解。

更准确的评价是：

> **一维复现是一套成功的 engineering reproduction，但 enrichment spatial half-step 不是 paper-exact Galerkin reproduction。**

它完成了当时最重要的任务：

```text
LATIN local/global structure 跑通
PGD enrichment 跑通
三材料杆主要响应复现
稳定 spatial/temporal modes 获得
一套可工作的 numerical strategy 被验证
```

所以它并不是“错误实现”。

真正需要修正的是以后对它的表述。

以后不宜笼统说：

```text
“一维三材料杆完全严格复现了原论文 LATIN-PGD”
```

而应更准确地说：

> **一维三材料杆成功复现了 LATIN-PGD 的主要算法结构与 benchmark numerical response，但 enrichment spatial half-step 采用 residual-driven weighted least-squares engineering implementation，而不是 paper Eq. (65)-(71) 的严格 Galerkin implementation。**

---

# 12. 这个历史背景如何重新解释当前核心困难

现在可以把整个迁移过程看得更清楚：

```text
原论文
│
├─ spatial:
│  Galerkin
│
└─ temporal:
   minimisation
```

一维 working implementation：

```text
1D
│
├─ spatial:
│  residual-LS
│
└─ temporal:
   BE residual minimisation
        │
        ↓
      收敛
```

当前 tower implementation：

```text
tower
│
├─ spatial:
│  paper Galerkin
│
└─ temporal:
   1D BE
        │
        ↓
   第 4 mode period-3
```

Case C：

```text
tower
│
├─ spatial:
│  1D-style residual-LS
│
└─ temporal:
   same BE
        │
        ↓
   27 iterations convergence
```

因此现在最值得研究的问题更加明确：

> **period-3 不仅是在暴露 tower 的结构复杂性，也是在暴露“1D engineering implementation”和“paper-exact spatial formulation”之间此前一直存在、但在简单一维 benchmark 上没有显现出来的差异。**

这里还需要保持一个来源边界：

```text
可以确认的事实：
1D 实际代码采用 residual-driven weighted LS；
原论文 spatial enrichment 采用 Galerkin。

合理重建的历史原因：
residual-LS 更容易实现和调试；
与当时“先稳定复现 1D benchmark”的目标一致；
论文没有给出全部离散实现细节。
```

后面这一组原因属于基于项目历史和代码结构的合理解释，不应写成原论文明确规定，或当时已有逐条书面决策记录。

---

# 13. 第 4 mode 到底发生了什么

正常情况下应该是：

```text
A
↓
B
↓
C
↓
D
↓
……
↓
stable fixed point
```

但当前塔筒第 4 mode 实际表现为：

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
↓
……
```

也就是说，它没有收敛到一个答案，而是在三个答案之间循环。

这就是我们反复看到的：

```text
period-3
```

对应的 convergence indicator 也一直近似循环：

```text
0.577
→ 0.658
→ 0.603
→ 0.577
→ 0.658
→ 0.603
→ ...
```

所以这个问题不是：

```text
“再多算几十次就好了”
```

而是：

> **当前 alternating map 本身把解送进了一个三周期，而不是 fixed point。**

---

# 14. 为什么我们现在认为问题集中在 spatial half-step

因为前面已经排除了很多可能原因。

我们怀疑过：

```text
是不是迭代次数不够？
```

不是。

我们怀疑过：

```text
是不是需要 damping？
```

没有解决。

我们怀疑过：

```text
是不是 seed 选错了？
```

多个 seed 都进入相同的三周期家族。

我们怀疑过：

```text
是不是第 4 mode 本身没有意义？
```

不是，第 4 mode 仍然能明显降低 residual。

我们怀疑过：

```text
是不是已有 3 个 mode 已经把空间占满了？
```

不是，新 mode 仍然有很高 novelty。

我们怀疑过：

```text
是不是缺少 1D-style in-loop orthogonalisation？
```

也不是。

正交化以后 basis overlap 已经降到机器精度，但 period-3 仍然存在。

---

# 15. 最关键的 A/B/C 诊断

真正重要的是下面这组比较。

## A：当前 tower 方法

```text
Eq.(70)-(71) spatial Galerkin
+
current BE temporal solve
```

结果：

```text
不收敛
→ period-3
```

---

## B：加上 1D-style in-loop orthogonalisation

```text
Eq.(70)-(71) spatial Galerkin
+
1D-style orthogonalisation
+
current BE temporal solve
```

结果：

```text
还是不收敛
→ period-3
```

所以问题不是 orthogonalisation。

---

## C：把 spatial half-step 换回 1D-style residual minimisation

```text
1D-style direct weighted-LS spatial solve
+
same orthogonalisation
+
same BE temporal solve
```

结果：

```text
27 iterations
→ converged
```

最终：

$$ \chi_{\rm fp}=7.55\times10^{-7}\lt10^{-6} $$

而且后期是很清楚地逐步下降。

所以这是目前最关键的证据：

> **只改变 spatial half-step 的求法，就能把 period-3 变回普通收敛。**

---

# 16. 这是不是说明 Eq. (70)-(71) 错了

不是。

这是目前最需要避免的误解。

我们已经重新检查过：

```text
paper equations
↓
tower derivation
↓
H / M formulation
↓
Eq.(70) FE system
↓
Eq.(71) recovery
↓
current code
```

目前没有发现明显的代数错误。

所以当前更合理的理解不是：

```text
“Eq.(70)-(71) 写错了。”
```

而是：

> **Eq.(70)-(71) 本来属于原论文整套 Galerkin + temporal minimisation framework；我们现在只保留了空间 Galerkin 这一半，而时间部分暂时用了自己的一维 BE 离散，因此这两个离散 half-step 组合后可能不再完全协调。**

---

# 17. 最通俗的比喻

可以把 spatial solver 和 temporal solver 看成两个合作的人。

在一维杆中：

```text
空间求解器：
我按 residual minimisation 判断什么答案最好。

时间求解器：
我也按 residual minimisation 判断什么答案最好。
```

两个人基本使用同一把尺子。

所以：

```text
A给B一个答案
↓
B修正
↓
再还给A
↓
A再修正
↓
……
↓
最终达成一致
```

---

# 18. 到了塔筒以后

现在变成：

```text
空间求解器：
我按 paper Galerkin 规则判断。

时间求解器：
我按 1D backward-Euler residual 规则判断。
```

两个人的规则都可能是对的。

但如果两把尺子不完全一致，就可能出现：

```text
空间说：最合适的是 A
↓
时间说：那我改成 B
↓
空间说：按我的规则应该是 C
↓
时间说：那又应该回到 A
```

于是：

```text
A → B → C → A
```

这就是 period-3 最直观的解释。

---

# 19. 所以真正的难点不是“塔筒太复杂”

塔筒复杂当然是真的。

但是目前 tower 的基础部分已经做出来了：

```text
fiber beam-column tower
material-point layout
H matrix
M metric
reference equilibrium
local stage
search directions
Trial A / Trial B
PGD basis
transaction / rollback
```

而且前 3 个 enrichment modes 也能成功生成。

所以现在不能再简单说：

```text
“因为塔筒复杂，所以不收敛。”
```

真正准确的说法应该是：

> **塔筒结构平衡的引入改变了 spatial PGD subproblem 的形式，而新的 spatial solver 与当前 temporal solver 是否仍然离散一致，是目前真正没有闭合的问题。**

---

# 20. 当前核心困难的正式表述

如果用一句学术语言概括：

> **The central difficulty in extending the validated 1D LATIN-PGD implementation to the offshore wind-turbine tower is not the increase in structural dimensionality itself, but the preservation of numerical consistency between the spatial and temporal PGD subproblems after structural equilibrium is introduced.**

中文可以写成：

> **将已验证的一维 LATIN-PGD 方法推广至海上风机塔筒时，核心困难并非结构维度增加本身，而是在引入塔筒结构平衡约束后，如何保持 PGD 空间子问题与时间子问题之间的数值一致性，使其交替迭代能够稳定收敛。**

---

# 21. 现在最应该记住的逻辑图

```text
validated 1D bar
│
├─ spatial:
│  residual minimisation
│
└─ temporal:
   BE residual minimisation
        │
        ↓
   两者数值上比较协调
        │
        ↓
      convergence
```

而 tower 当前是：

```text
tower
│
├─ spatial:
│  paper Eq.(70)-(71) Galerkin
│
└─ temporal:
   validated 1D BE update
        │
        ↓
   discrete compatibility 尚未证明
        │
        ↓
   fourth-mode period-3
```

而 diagnostic C 是：

```text
tower diagnostic C
│
├─ spatial:
│  1D-style residual minimisation
│
└─ temporal:
   same BE update
        │
        ↓
   ordinary contraction
        │
        ↓
   converged in 27 iterations
```

---

# 22. 所以下一步真正要研究什么

下一步不是继续做更大的 tower model，也不是马上加更多循环。

真正应该先回答：

> **怎样让 tower spatial solver 和 temporal solver 重新处在同一套一致的离散 LATIN-PGD 逻辑中？**

具体就是检查：

```text
paper Eq.(70)-(71) spatial Galerkin
```

与：

```text
current causal BE Eq.(72) temporal solve
```

到底是不是离散一致的。

如果不是，就要决定：

```text
路线 A：
保留 paper Galerkin
→ 恢复与之更一致的 temporal discretisation

或

路线 B：
第一阶段优先稳定迁移
→ 采用 validated 1D-style spatial residual minimisation
→ 明确标注为 tower-v1 engineering discretisation
```

---

# 23. 当前暂时不要做的事情

在这个问题闭合之前，不应该急着：

```text
直接删除 Eq.(70)-(71)
直接接受 period-3 phase
放宽 fixed-point convergence gate
大量增加 damping
引入 Anderson / Aitken
进入 high-cycle simulation
引入 cycle-phase PGD
引入 multi-time-scale PGD
```

因为这些都会绕开当前真正的核心问题。

---

# 24. 当前阶段最终结论

最简单地说：

> **一维杆成功，是因为空间和时间两个求解器比较“合拍”；塔筒第 4 mode 出问题，是因为我们为了忠实原论文改变了空间求解方式，但时间求解仍沿用一维 BE 方法，这两个 half-step 现在可能“不合拍”。**

目前最强的数值证据是：

```text
paper-style spatial + BE temporal
→ period-3

1D-style spatial + same BE temporal
→ convergence
```

因此，当前整个 tower LATIN-PGD 研究的核心已经从：

```text
“为什么第 4 mode 不收敛？”
```

收缩为：

```text
“如何恢复 spatial 与 temporal PGD subproblems
在离散层面的数值一致性？”
```

这就是下一阶段最重要的研究问题。
