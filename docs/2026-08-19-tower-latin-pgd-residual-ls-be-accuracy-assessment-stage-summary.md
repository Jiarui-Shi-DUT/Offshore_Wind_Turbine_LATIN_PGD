# Tower LATIN-PGD：residual-LS spatial + BE temporal 的准确性判断阶段总结

**日期：2026-08-19**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**研究对象：将 LATIN-PGD 由一维三材料杆推广至 fiber beam-column offshore wind turbine tower**  
**本阶段主题：如果 tower 继续采用 residual-LS spatial + BE temporal，这套方法是否仍然可以得到准确结果？**  
**本阶段性质：方法路线判断与验证标准冻结，不修改 `latin/` 核心代码。**

---

# 1. 本阶段最核心的问题

当前已经发现：

```text
paper-style Eq.(70)-(71) spatial Galerkin
+
current BE temporal solve
→ fourth-mode period-3
```

而 diagnostic Case C 表明：

```text
1D-style residual-LS spatial solve
+
same BE temporal solve
→ ordinary convergence
```

因此自然产生一个重要问题：

> **如果在风机塔筒算例中不再坚持 paper Eq. (70)-(71) spatial Galerkin，而是继续采用已经在 1D 中验证过的 residual-LS spatial + BE temporal，那么最终得到的 tower response 还会不会准确？**

本阶段的核心结论是：

> **完全有可能准确。是否准确，不能由“是不是完全照搬原论文公式”来判断，而应该由它是否收敛到正确的 tower full-order physical solution 来判断。**

---

# 2. 必须先区分两种“准确”

这是整个问题里最重要的概念区分。

第一种是：

```text
algorithmic fidelity
```

即：

```text
是不是严格照原论文 Eq.(65)-(72) 的算法实现？
```

第二种是：

```text
physical / numerical accuracy
```

即：

```text
最终塔筒位移、应力、塑性应变、损伤等
是不是逼近正确的 full-order solution？
```

这两种准确不是一回事。

因此：

```text
不是 paper-exact
≠
物理结果不准确
```

更直接地说：

```text
是不是和原论文算法一模一样？
→ 不一定

是不是仍然可能把 tower response 算得很准？
→ 完全可能
```

---

# 3. 为什么“和论文不完全一样”不等于“结果不准确”

对于同一个非线性有限元问题，可以有不同的 numerical solvers，例如：

```text
Newton-Raphson
Modified Newton
arc-length
```

它们并不完全一样，但只要最终都收敛到同一个平衡解，结果都可以准确。

LATIN-PGD 也是同样的道理。

原论文选择的是：

```text
spatial:
Galerkin formulation

temporal:
minimisation
```

这是一套有效的 PGD enrichment strategy。

但它并不意味着：

```text
只有 Galerkin spatial
才能得到正确 tower response
```

如果 residual-LS spatial 也能够：

```text
满足结构平衡
降低机械残差
稳定生成 PGD modes
最终逼近 full-order solution
```

那么它同样可以是一套准确的 reduced-order formulation。

---

# 4. residual-LS spatial 并没有丢掉 tower structural equilibrium

这是一个非常容易产生的误解。

如果 spatial half-step 不再使用 Eq. (70)-(71)，并不意味着：

```text
不管结构平衡了
```

Case C 中仍然保留 tower equilibrium mapping：

$$ s=A_\sigma p $$

这里的 $A_\sigma$ 不是任意矩阵。

它来自已经建立好的 tower reference equilibrium operator。

给定 spatial plastic mode：

$$ p $$

程序仍然通过 tower equilibrium problem 得到：

$$ s=A_\sigma p $$

并满足：

$$ H^TMs=0 $$

因此：

$$ \boxed{\text{tower structural equilibrium 仍然被保留}} $$

改变的只是：

> **在所有满足 tower equilibrium mapping 的 candidate spatial modes 中，用什么准则选出最合适的那个 $p$。**

---

# 5. Galerkin 和 residual-LS 的区别到底是什么

可以把它们理解成两种不同的“选最佳空间模态”的规则。

## 5.1 paper Galerkin

给定 temporal function 后，通过 Galerkin weak condition 选取 spatial mode。

直观上可以理解为：

```text
让误差对规定的 test space 满足正交条件
```

---

# 6. residual-LS

residual-LS 的思路更直接：

```text
直接找一个 spatial mode
让当前 mechanical residual 的加权平方和最小
```

即：

$$ p=\operatorname*{arg\,min}_p J_h(p) $$

例如当前 diagnostic 中：

$$ J_h(p)=\frac12\sum_nw_n r_n(p)^TMD_{H,n}^{-1}r_n(p) $$

所以 residual-LS 的目标非常清楚：

> **每增加一个 spatial mode，都直接寻找一个最能降低当前 mechanical residual 的方向。**

这本身是一个有明确数学意义的 numerical criterion，并不是随意拟合。

---

# 7. 为什么 Case C 是一个很积极的信号

Case C 中保持不变的是：

```text
tower geometry
tower material-point layout
tower equilibrium operator
current Trial-A state
rank-3 basis
residual seed
H_sigma
M
BE temporal solve
fixed-point convergence criterion
```

主要改变的是：

```text
spatial half-step
```

从：

```text
paper Eq.(70)-(71) Galerkin
```

改为：

```text
1D-style direct weighted residual-LS
```

结果：

```text
period-3
↓
ordinary contraction
↓
converged in 27 iterations
```

最终：

$$ \chi_{\rm fp}=7.5524\times10^{-7}\lt10^{-6} $$

因此 Case C 至少证明：

> **residual-LS spatial 与 current BE temporal solver 在当前 tower fourth-mode benchmark 中具有更好的 numerical compatibility。**

---

# 8. 但“第 4 mode 收敛”还不能直接等于“整个 tower solution 准确”

这是必须严格保持的边界。

Case C 当前证明的是：

```text
inner enrichment fixed point
可以收敛
```

它还没有自动证明：

```text
整个 LATIN-PGD tower solution
已经逼近 full-order physical solution
```

这两者不是同一个层次。

所以现在不能直接写：

```text
residual-LS + BE
已经被证明是准确 tower method
```

目前只能写：

> **它已经表现出良好的内部 fixed-point convergence，但最终 solution accuracy 还必须通过 full-order benchmark 验证。**

---

# 9. 真正判断“准确”的裁判是谁

真正的 reference 不应该是 Eq. (70)-(71) 本身。

真正的裁判应该是：

$$ \boxed{\text{tower Full-Order Model (FOM)}} $$

也就是不用 PGD reduction，而是正常逐时间步执行 nonlinear tower solve：

```text
time step 1
→ nonlinear equilibrium solve

time step 2
→ nonlinear equilibrium solve

...

time step Nt
→ nonlinear equilibrium solve
```

这套 full-order tower result 应作为 reference solution。

---

# 10. LATIN-PGD 应该和 FOM 比较什么

如果最终采用：

```text
residual-LS spatial
+
BE temporal
```

那么至少要和 FOM 比较：

```text
tower-top displacement history
element / section force history
fiber stress history
fiber plastic-strain history
damage history
global equilibrium residual
LATIN residual
```

理想状态应满足：

$$ u_{\rm PGD}\approx u_{\rm FOM} $$

$$ \sigma_{\rm PGD}\approx\sigma_{\rm FOM} $$

$$ \varepsilon^p_{\rm PGD}\approx\varepsilon^p_{\rm FOM} $$

$$ D_{\rm PGD}\approx D_{\rm FOM} $$

如果这些都随着 LATIN iteration 和 PGD rank 增加而逼近 FOM，那么这套方法就是数值上准确的。

---

# 11. 为什么不能用“是不是 Galerkin”直接判断准确性

假设 FOM 得到 tower-top maximum displacement：

```text
0.1250 m
```

如果 paper-Galerkin LATIN-PGD 得到：

```text
0.1248 m
```

而 residual-LS LATIN-PGD 得到：

```text
0.1251 m
```

那么不能因为第二种方法没有使用 Eq. (70)-(71)，就说：

```text
0.1251 m 不准确
```

真正需要看的应该是：

```text
哪一种方法更稳定地逼近 FOM
```

而不是：

```text
哪一种方法中间的 basis generation
更像原论文
```

---

# 12. residual-LS 可能具有的现实优势

对于 tower fiber system：

```text
many material points
strong structural coupling
distributed nonlinearities
complex stress redistribution
```

直接围绕 mechanical residual 构造 spatial mode，有一个很自然的 numerical advantage：

> **每个新 mode 都直接针对当前尚未消除的 global mechanical error。**

也就是说：

```text
当前 residual 在哪里大
↓
新 spatial mode 就优先解决哪里
```

从 reduced-order modelling 的角度，这是一种合理而直接的 enrichment philosophy。

目前 Case C 已经表明：

```text
这种 residual-driven spatial correction
在当前 fourth-mode benchmark 中
比 mixed Galerkin + BE map 更容易形成 contraction
```

但目前仍然只能把这看成一个积极数值信号，而不是最终普适结论。

---

# 13. 如果正式采用 residual-LS + BE，最大的变化是什么

最大的变化不是：

```text
方法一定不准确
```

而是：

```text
它不能继续被称为 paper-exact Eq.(65)-(72) implementation
```

如果最终 tower v1 正式采用：

```text
residual-LS spatial
+
BE temporal
```

那么在论文和文档里应该明确把它描述成：

```text
tower-specific LATIN-PGD discretisation
```

或者：

```text
residual-based LATIN-PGD enrichment strategy
```

而不是：

```text
direct implementation of paper Eq.(65)-(72)
```

---

# 14. 这反而可能形成 tower-specific 方法贡献

当前博士研究目标不是简单地逐字复制原论文。

真正目标是：

> **将 LATIN-PGD 推广到由 fiber beam-column elements 组成的 offshore wind turbine tower fatigue problem。**

如果后续能够证明：

```text
paper-style spatial Galerkin + current BE
→ fixed-point pathology

residual-LS spatial + BE
→ stable convergence

同时：
residual-LS + BE
→ accurately reproduces FOM
```

那么这就不再只是一个“临时修补”。

它可以被发展成：

> **针对 fiber beam-column tower 的 residual-based LATIN-PGD enrichment strategy。**

这种工作具有明确的方法扩展意义。

---

# 15. 但目前不能因为 Case C 成功就立即修改 production core

当前应保持谨慎。

Case C 只证明：

```text
fourth-mode fixed point
恢复收敛
```

下一步还必须确认：

```text
rank 4
rank 5
rank 6
...
```

是否仍然保持稳定。

同时还要确认：

```text
outer LATIN iteration
是否稳定下降
```

以及：

```text
最终 physical response
是否逼近 FOM
```

所以当前不应该立即：

```text
删除 Eq.(70)-(71)
把 Case C diagnostic 直接写进 production core
声明 residual-LS 已经最终胜出
```

---

# 16. 如果走 residual-LS + BE 路线，需要通过哪四个核心验证

如果这四个验证全部通过，那么 residual-LS + BE 就具备成为正式 tower LATIN-PGD formulation 的基础。

---

# 17. 验证 1：结构平衡

需要确认 PGD correction 和最终 response 都满足：

$$ H^TM\sigma\approx0 $$

也就是：

```text
global tower equilibrium
不能因为 reduced representation 而被破坏
```

当前 tower equilibrium operator 已经为这一点提供了基础。

---

# 18. 验证 2：LATIN / mechanical residual 必须真正下降

随着 enrichment rank 增加，应表现为：

$$ \|R\|_{r+1}\lt\|R\|_r $$

并最终趋于所设 tolerance。

不仅要看：

```text
inner fixed-point converged
```

还要看：

```text
new mode 是否真正降低 Trial-A mechanical residual
```

以及：

```text
outer LATIN convergence indicator 是否下降
```

---

# 19. 验证 3：PGD solution 必须逼近 FOM

这是最重要的 accuracy criterion。

例如 displacement error：

$$ e_u=\frac{\|u_{\rm PGD}-u_{\rm FOM}\|}{\|u_{\rm FOM}\|} $$

stress error：

$$ e_\sigma=\frac{\|\sigma_{\rm PGD}-\sigma_{\rm FOM}\|}{\|\sigma_{\rm FOM}\|} $$

plastic strain error：

$$ e_{\varepsilon^p}=\frac{\|\varepsilon^p_{\rm PGD}-\varepsilon^p_{\rm FOM}\|}{\|\varepsilon^p_{\rm FOM}\|} $$

damage error：

$$ e_D=\frac{\|D_{\rm PGD}-D_{\rm FOM}\|}{\|D_{\rm FOM}\|} $$

随着 PGD rank 和 LATIN iteration 增加，这些误差应进入稳定且足够小的范围。

---

# 20. 验证 4：离散收敛性

还需要检查结果是否对 numerical discretisation 具有合理稳定性。

至少应考虑：

```text
time-step refinement
mesh refinement
fiber-number refinement
PGD-rank refinement
```

即：

```text
dt ↓
mesh finer
fiber count ↑
rank ↑
```

以后结果应趋于稳定。

否则即使某一个 coarse benchmark 看起来与 FOM 接近，也不能说明方法已经可靠。

---

# 21. 当前可以形成两条明确研究路线

## 路线 A：paper fidelity first

坚持：

```text
paper Galerkin spatial
```

那么下一步必须研究：

```text
它真正应该匹配什么 temporal discretisation
```

这很可能最终需要重新恢复：

```text
paper DG0 temporal treatment
```

这条路线的优点是：

```text
理论忠实度最高
```

缺点是：

```text
当前 Galerkin + BE 已出现 period-3
DG0 exact discrete details 尚未完全恢复
```

---

# 22. 路线 B：tower migration stability first

采用：

```text
residual-LS spatial
+
validated BE temporal
```

然后通过严格的 FOM benchmark 和 convergence study 验证。

如果验证通过，那么可以把它定义为：

```text
tower-specific residual-based LATIN-PGD
```

或者：

```text
engineering LATIN-PGD discretisation for fiber beam-column towers
```

这条路线的优点是：

```text
继承已验证 1D numerical pair
当前 Case C 已显示稳定 contraction
便于建立统一 residual-based enrichment logic
```

---

# 23. 当前最合理的研究态度

目前不能写：

```text
residual-LS 一定比 Galerkin 好
```

也不能写：

```text
Galekin 一定不适用于 tower
```

当前最准确的结论是：

> **residual-LS spatial + BE temporal 已经成为一个非常有竞争力的 tower LATIN-PGD candidate formulation，但其最终准确性仍需通过 full-order comparison、rank convergence、LATIN convergence 和 discretisation convergence 系统验证。**

---

# 24. 当前应如何理解 Case C 的地位

Case C 目前的角色应该定义为：

```text
不是最终 production algorithm

而是：

第一个解除 fourth-mode period-3 的
causally isolated numerical candidate
```

它证明了：

```text
tower geometry 本身
不是 fourth-mode failure 的充分原因
```

也证明了：

```text
current BE temporal solve
并非必然导致 nonconvergence
```

因为同一个 BE temporal solve 在 residual-LS spatial 下可以收敛。

所以问题被进一步集中到了：

```text
spatial formulation
+
spatial-temporal discrete compatibility
```

---

# 25. 最通俗的结论

可以用一句大白话理解：

> **换成 residual-LS spatial 并不是把塔筒物理问题“简化错了”，而只是换了一种挑选 PGD 空间模态的规则。塔筒的结构平衡、本构关系和 material-point system 仍然保留。只要最终能稳定收敛，并且和 full-order tower solution 对得上，这套方法就可以是准确的。**

---

# 26. 最终阶段结论

当前可以冻结以下认识。

**结论 1**

$$ \boxed{\text{paper fidelity 与 physical accuracy 不是同一概念}} $$

---

**结论 2**

residual-LS spatial + BE temporal 不是 paper-exact Eq. (65)-(72) implementation。

---

**结论 3**

但 residual-LS spatial 仍然通过 tower equilibrium operator 保留结构平衡，因此不能简单理解为“失去物理约束”。

---

**结论 4**

Case C 已经证明这套 spatial-temporal pair 在当前 fourth-mode benchmark 中可以恢复 ordinary fixed-point convergence。

---

**结论 5**

Case C 的成功不能单独证明 whole tower solution accuracy。

---

**结论 6**

真正的 accuracy criterion 应是：

$$ \boxed{\text{LATIN-PGD solution}\rightarrow\text{tower FOM solution}} $$

---

**结论 7**

如果 residual-LS + BE 同时通过：

```text
equilibrium verification
residual convergence
FOM accuracy comparison
discretisation convergence
```

那么它完全可以成为一套准确的 tower LATIN-PGD formulation。

---

**结论 8**

如果最终正式采用，应明确把它定义为：

```text
tower-specific residual-based LATIN-PGD discretisation
```

而不是声称为 paper-exact Galerkin implementation。

---

# 27. 下一阶段建议

在继续修改 production core 之前，建议先完成：

```text
residual-LS + BE accuracy validation plan
```

重点设计一个最小但完整的验证矩阵：

```text
FOM
vs
LATIN-PGD residual-LS + BE
```

至少比较：

```text
tower-top displacement
representative fiber stress
plastic strain
damage
mechanical residual
LATIN indicator
PGD rank convergence
```

只有当这一层闭合以后，才值得决定是否将 residual-LS spatial 从 diagnostic candidate 升级为 tower v1 正式方法。
