# Tower LATIN-PGD Eq. (65)–(69) 空间许可条件与空间子问题阶段总结

**日期：2026-08-17**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**当前研究路线：原论文 $x-t$ LATIN-PGD → fiber beam-column offshore wind turbine tower**  
**阶段范围：原论文 Eq. (65)–(69)**  
**上一阶段衔接：`2026-08-16-tower-latin-pgd-eq61-64-enrichment-stage-summary.md`**  
**下一阶段：Eq. (70) displacement-like spatial FE problem**

---

# 1. 本阶段定位

上一阶段已经完成 Eq. (61)–(64) 的逐式推导，并明确了 new-mode enrichment 的基本逻辑：

$$ \boxed{ \text{fixed-basis temporal update} \rightarrow \text{saturation check} \rightarrow \text{new-mode enrichment} } $$

在 enrichment 阶段，新的 PGD pair 写为：

$$ \boxed{ \Delta\dot{\varepsilon}^{p} = \dot{\lambda}(t)\bar{\varepsilon}^{p}(x), \qquad \Delta\sigma' = \lambda(t)\bar{\sigma}(x) } $$

其中：

- $\lambda(t)$：new temporal function；
- $\bar{\varepsilon}^{p}(x)$：new plastic-strain spatial mode；
- $\bar{\sigma}(x)$：由 spatial mode 诱导的 equilibrated stress spatial mode。

上一阶段 Eq. (64) 给出：

$$ \boxed{ \Delta\dot{\varepsilon}' = \Delta\dot{\varepsilon}^{p} + \mathbb C^{-1}\Delta\dot{\sigma}' } \tag{64} $$

再结合 Eq. (61)：

$$ \Delta\dot{\varepsilon}^{p} - H_\sigma\Delta\sigma' + \bar{\Delta} = 0 $$

以及 Eq. (63)：

$$ \Delta\sigma' = \lambda\bar{\sigma}, \qquad \Delta\dot{\sigma}' = \dot{\lambda}\bar{\sigma} $$

可得到：

$$ \boxed{ \Delta\dot{\varepsilon}' = H_\sigma\lambda\bar{\sigma} + \dot{\lambda}\mathbb C^{-1}\bar{\sigma} - \bar{\Delta} } $$

本阶段 Eq. (65)–(69) 的任务，就是把这一 material-point / local-search-direction relation 重新放回 global admissibility manifold，并逐步转化成一个可以由标准有限元空间离散求解的 spatial problem。

本阶段最核心的逻辑链为：

$$ \boxed{ \text{Eq. (65) space-time kinematic admissibility} \rightarrow \text{Eq. (66) temporal Galerkin projection} \rightarrow \text{Eq. (67) effective spatial operator} \rightarrow \text{Eq. (68) compatible auxiliary strain field} \rightarrow \text{Eq. (69) equilibrated stress field} } $$

---

# 2. 与原论文整体 LATIN 架构的关系

LATIN 将问题划分为：

$$ \boxed{ A = \text{global admissibility manifold} } $$

与：

$$ \boxed{ \Gamma = \text{local constitutive/evolution manifold} } $$

其中 global admissibility 主要包括：

1. kinematic admissibility；
2. static admissibility。

因此 new-mode enrichment 不能只满足 Eq. (61) 的 descent search direction，也不能只在 material points 上拟合 remaining defect。

新的 spatial mode 最终必须同时满足：

$$ \boxed{ \text{kinematic compatibility} + \text{static equilibrium} } $$

Eq. (65)–(69) 正是在完成这一过程。

这也是 LATIN-PGD enrichment 与单纯 snapshot SVD/POD enrichment 的根本区别之一。

---

# 3. Eq. (65)：kinematic admissibility in rate form

原论文 Eq. (65) 可写为：

$$ \boxed{ \int_{[0,T]\times\Omega} \Delta\dot{\varepsilon}' : \sigma^* \, d\Omega\,dt = 0, \qquad \forall \sigma^* \in S_0 } \tag{65} $$

这里：

- $\Delta\dot{\varepsilon}'$：当前 enrichment correction 对应的 total strain-rate correction；
- $\sigma^*$：test stress field；
- $S_0$：zero-load statically admissible stress space。

Eq. (65) 的本质不是新的 constitutive equation，也不是新的 equilibrium equation，而是：

$$ \boxed{ \Delta\dot{\varepsilon}' \text{ 必须属于 homogeneous compatible strain-rate correction space} } $$

---

# 4. 为什么 kinematic admissibility 要用 stress field 测试

论文不是直接通过强式：

$$ \Delta\dot{\varepsilon}' = \nabla^{sym}\Delta\dot{u} $$

实施 compatibility，而是使用 admissibility spaces 的对偶 weak characterization。

对 homogeneous correction，有：

$$ \Delta u = 0 \qquad \text{on prescribed-displacement boundary} $$

因此：

$$ \Delta u \in U_0 $$

对应的 compatible strain correction 属于：

$$ E_0 = \{ \varepsilon(u) \, | \, u \in U_0 \} $$

而 kinematically admissible strain space 与 zero-load statically admissible stress space 之间满足 work-conjugate orthogonality：

$$ \boxed{ E_0 = S_0^\perp } $$

这里的 pairing 为：

$$ \langle \varepsilon,\sigma \rangle = \int_\Omega \varepsilon:\sigma\,d\Omega $$

因此 Eq. (65) 用 $\sigma^*\in S_0$ 来测试 $\Delta\dot{\varepsilon}'$，并不是“用应力检查平衡”，而是利用 strain-stress duality 检查 compatibility。

---

# 5. Eq. (65) 右端为什么为零

这是本阶段首先澄清的重要问题。

真实结构的总场满足非齐次边界条件与外荷载，因此原始 admissibility weak form 通常有非零右端。

但是当前求解的是 enrichment correction。

若 baseline state 与 new state 都满足同一 prescribed displacement：

$$ u^{up} = u_d, \qquad u^{new} = u_d \qquad \text{on } \partial_1\Omega $$

则：

$$ \Delta u = u^{new} - u^{up} = 0 \qquad \text{on } \partial_1\Omega $$

因此 correction 是 homogeneous correction：

$$ \boxed{ \Delta u \in U_0 } $$

这使对应的 weak boundary term 消失，从而得到 Eq. (65) 的零右端。

因此必须记住：

$$ \boxed{ \text{Eq. (65) RHS}=0 \text{ 不是因为真实结构没有外荷载，而是因为当前求解的是 homogeneous correction} } $$

---

# 6. Eq. (65) 的 strong-form 直观解释

若：

$$ \Delta\dot{\varepsilon}' = \nabla^{sym}\Delta\dot{u} $$

则：

$$ \int_\Omega \sigma^*:\Delta\dot{\varepsilon}'\,d\Omega = \int_\Omega \sigma^*:\nabla^{sym}\Delta\dot{u}\,d\Omega $$

经分部积分：

$$ \int_\Omega \sigma^*:\nabla^{sym}\Delta\dot{u}\,d\Omega = -\int_\Omega (\nabla\cdot\sigma^*)\cdot\Delta\dot{u}\,d\Omega + \int_{\partial\Omega}(\sigma^*n)\cdot\Delta\dot{u}\,dS $$

因为 $\sigma^*\in S_0$：

$$ \nabla\cdot\sigma^* = 0 $$

在 Neumann 边界：

$$ \sigma^*n = 0 $$

在 Dirichlet 边界：

$$ \Delta\dot{u} = 0 $$

因此：

$$ \boxed{ \int_\Omega \sigma^*:\Delta\dot{\varepsilon}'\,d\Omega = 0 } $$

这与 Eq. (65) 完全一致。

---

# 7. 为什么 Eq. (65) 对整个时间区间积分

LATIN global iteration 面向整个：

$$ [0,T]\times\Omega $$

space-time domain。

因此 Eq. (65) 不是逐时间步 Newton equilibrium，而是：

$$ \boxed{ \text{entire-history space-time admissibility condition} } $$

这与原论文 $x-t$ PGD 架构完全一致。

当前阶段仍不引入：

- cycle-phase separation；
- $n-\tau-x$ PGD；
- cycle jump；
- temporal homogenization。

---

# 8. Eq. (65) 中的 $S_0$ 到底是什么

$S_0$ 是 zero-load statically admissible stress vector space。

其 weak definition 可以理解为：

$$ \boxed{ \int_\Omega \sigma^*:\varepsilon(u^*)\,d\Omega = 0, \qquad \forall u^*\in U_0 } $$

从 strong form 直观理解：

$$ \nabla\cdot\sigma^* = 0 $$

并且在 Neumann boundary：

$$ \sigma^*n = 0 $$

但在 prescribed-displacement boundary 上允许 non-zero reaction traction。

因此：

$$ \boxed{ \sigma^* \in S_0 \text{ 并不意味着所有边界 traction 都为零} } $$

---

# 9. Eq. (65) 与 Eq. (69) 的职责不能混淆

本阶段明确区分：

Eq. (65)：

$$ \boxed{ \text{kinematic admissibility of } \Delta\dot{\varepsilon}' } $$

它使用：

$$ \sigma^*\in S_0 $$

作为 test field。

Eq. (69)：

$$ \boxed{ \text{static admissibility of } \bar{\sigma} } $$

它使用：

$$ u^*\in U_0 $$

作为 test field。

二者构成经典虚功对偶关系。

---

# 10. Eq. (66) 的起点：固定 temporal function 求 spatial mode

new PGD pair 为：

$$ \lambda(t)\bar{\varepsilon}^{p}(x) $$

因为 $\lambda(t)$ 与 $\bar{\varepsilon}^{p}(x)$ 双线性耦合，原论文采用 fixed-point alternating enrichment：

$$ \boxed{ \lambda^{(k)} \rightarrow \bar{\varepsilon}^{p,(k+1)} \rightarrow \lambda^{(k+1)} \rightarrow \cdots } $$

Eq. (66) 属于 spatial half-step。

因此此时：

$$ \boxed{ \lambda(t) \text{ 暂时固定并视为已知} } $$

需要求的是新的 spatial mode。

---

# 11. Eq. (65) 中选择 separable test field

在 spatial subproblem 中选择：

$$ \boxed{ \sigma^*(x,t) = \lambda(t)\bar{\sigma}^*(x) } $$

其中：

$$ \bar{\sigma}^*(x)\in S_0 $$

这意味着 trial temporal function 与 test temporal function 采用同一个当前固定的 $\lambda(t)$。

因此这一操作是 temporal Galerkin projection：

$$ \boxed{ \text{trial temporal function} = \text{test temporal function} = \lambda(t) } $$

其目的是把 $(x,t)$ space-time weak problem 压缩为纯 $x$ spatial weak problem。

---

# 12. Eq. (66) 的逐项推导

由 Eq. (61)–(64)：

$$ \Delta\dot{\varepsilon}' = H_\sigma\lambda\bar{\sigma} + \dot{\lambda}\mathbb C^{-1}\bar{\sigma} - \bar{\Delta} $$

代入 Eq. (65)，并取：

$$ \sigma^* = \lambda\bar{\sigma}^* $$

得到：

$$ \int_0^T\int_\Omega \left( H_\sigma\lambda\bar{\sigma} + \dot{\lambda}\mathbb C^{-1}\bar{\sigma} - \bar{\Delta} \right) : \left( \lambda\bar{\sigma}^* \right)\,d\Omega\,dt = 0 $$

逐项展开后得到三个 temporal products：

$$ H_\sigma\lambda \times \lambda \rightarrow H_\sigma\lambda^2 $$

$$ \dot{\lambda} \times \lambda \rightarrow \lambda\dot{\lambda} $$

$$ \bar{\Delta} \times \lambda \rightarrow \bar{\Delta}\lambda $$

因此：

$$ \boxed{ \int_\Omega \left[ \left\langle H_\sigma\lambda^2 \right\rangle\bar{\sigma} + \left\langle \lambda\dot{\lambda} \right\rangle\mathbb C^{-1}\bar{\sigma} - \left\langle \bar{\Delta}\lambda \right\rangle \right] : \bar{\sigma}^*\,d\Omega = 0 } \tag{66} $$

其中：

$$ \langle a\rangle = \int_0^T a(t)\,dt $$

---

# 13. Eq. (66) 三项的来源

第一项：

$$ \boxed{ \left\langle H_\sigma\lambda^2 \right\rangle\bar{\sigma} } $$

来源于：

$$ H_\sigma\lambda\bar{\sigma} : \lambda\bar{\sigma}^* $$

其中 $\lambda^2$ 是 trial temporal mode 与 test temporal mode 相乘得到的 temporal Galerkin weight。

第二项：

$$ \boxed{ \left\langle \lambda\dot{\lambda} \right\rangle\mathbb C^{-1}\bar{\sigma} } $$

来源于：

$$ \Delta\dot{\sigma}' = \dot{\lambda}\bar{\sigma} $$

与 test temporal function $\lambda$ 的乘积。

第三项：

$$ \boxed{ -\left\langle \bar{\Delta}\lambda \right\rangle } $$

是 remaining descent-search-direction defect 沿当前 temporal direction $\lambda(t)$ 的 projection。

---

# 14. $\langle H_\sigma\lambda^2\rangle$ 一般不是全结构常数

因为：

$$ H_\sigma = H_\sigma(x,t) $$

所以：

$$ \left\langle H_\sigma\lambda^2 \right\rangle = \int_0^T H_\sigma(x,t)\lambda^2(t)\,dt $$

一般仍然依赖空间位置：

$$ \boxed{ \left\langle H_\sigma\lambda^2 \right\rangle = A(x) } $$

对 tower：

$$ x \rightarrow q=(e,g,f) $$

因此：

$$ \boxed{ A_q = \int_0^T H_{\sigma,q}(t)\lambda^2(t)\,dt } $$

不能把它错误地实现成一个全塔统一 scalar。

---

# 15. $\langle\lambda\dot{\lambda}\rangle$ 的边界性质

有：

$$ \lambda\dot{\lambda} = \frac{1}{2}\frac{d}{dt}\lambda^2 $$

因此：

$$ \boxed{ \left\langle\lambda\dot{\lambda}\right\rangle = \frac{1}{2}\left[\lambda^2(T)-\lambda^2(0)\right] } $$

但本阶段明确：

$$ \boxed{ \text{不能因为研究循环荷载就预先删除这一项} } $$

只有在明确满足：

$$ \lambda(T)=\lambda(0) $$

或其他使两端平方相同的特殊条件时，该积分才为零。

PGD temporal function $\lambda(t)$ 并不等同于外荷载周期函数，因此原论文保留此项，我们的 tower migration 也应保留。

---

# 16. $\langle\bar{\Delta}\lambda\rangle$ 的算法意义

定义：

$$ \boxed{ d(x) = \left\langle\bar{\Delta}\lambda\right\rangle = \int_0^T \bar{\Delta}(x,t)\lambda(t)\,dt } $$

这是将原：

$$ \bar{\Delta}(x,t) $$

沿 temporal direction $\lambda(t)$ contraction 后得到的 pure spatial forcing。

因此：

$$ \boxed{ \bar{\Delta}(x,t) \rightarrow d(x) } $$

从离散矩阵角度，如果：

$$ \bar{\Delta}\in\mathbb R^{N_q\times N_t} $$

则 temporal contraction 后：

$$ d\in\mathbb R^{N_q} $$

这是 Eq. (66) 将整个 space-time enrichment problem 压缩为 spatial problem 的关键一步。

---

# 17. Eq. (66) 的真正未知量

在 spatial fixed-point half-step 中：

- $\lambda(t)$ 已知；
- $H_\sigma(x,t)$ 已知；
- $\bar{\Delta}(x,t)$ 已知；
- $\mathbb C$ 已知。

因此 Eq. (66) 中真正需要求的是：

$$ \boxed{ \bar{\sigma}(x) } $$

但它不能是任意 material-point stress field，而必须属于：

$$ \boxed{ S_0 } $$

所以 Eq. (66) 是一个受 static admissibility constraint 限制的 spatial weak problem。

---

# 18. Eq. (66) 不能逐 material point 强制括号为零

不能从：

$$ \int_\Omega [\cdots]:\bar{\sigma}^*\,d\Omega = 0 \qquad \forall\bar{\sigma}^*\in S_0 $$

直接推出：

$$ [\cdots](x)=0 \qquad \forall x $$

原因是：

$$ \bar{\sigma}^* $$

不是任意 test field，而受：

$$ \bar{\sigma}^*\in S_0 $$

限制。

因此 Eq. (66) 是：

$$ \boxed{ \text{weak orthogonality condition} } $$

而不是：

$$ \boxed{ \text{pointwise constitutive equation} } $$

这解释了为什么不能逐 fiber 直接计算 $\bar{\sigma}_q$ 后就认为 spatial mode 已合法。

---

# 19. Eq. (66) 与 Eq. (59) 的区别

Eq. (59)：

$$ \boxed{ \text{known spatial basis} \rightarrow \text{solve temporal coefficients} } $$

Eq. (66)：

$$ \boxed{ \text{fixed temporal function} \rightarrow \text{solve new spatial mode} } $$

二者对应 fixed-point alternating enrichment 的两个不同 half-step。

因此不能把 Eq. (59) 的 temporal least-squares solve 与 Eq. (66) 的 spatial Galerkin solve 混为同一个 reduced system。

---

# 20. Eq. (67)：定义 effective spatial operator

Eq. (66) 中前两项都乘同一个 $\bar{\sigma}$，因此定义：

$$ \boxed{ W^{-1} = \left\langle H_\sigma\lambda^2 \right\rangle + \left\langle \dot{\lambda}\lambda \right\rangle\mathbb C^{-1} } $$

以及：

$$ \boxed{ \bar{\delta} = \left\langle \bar{\Delta}\lambda \right\rangle } $$

再定义：

$$ \boxed{ \bar{\tilde{\varepsilon}} = W^{-1}\bar{\sigma} - \bar{\delta} } \tag{67} $$

Eq. (67) 本身不引入新的 physics，而是对 Eq. (66) 的 pure spatial quantities 做重新组织。

---

# 21. $W^{-1}$ 的准确含义

$W^{-1}$ 由两部分组成：

$$ \left\langle H_\sigma\lambda^2 \right\rangle $$

对应 temporally projected LATIN descent search-direction contribution；

$$ \left\langle \dot{\lambda}\lambda \right\rangle\mathbb C^{-1} $$

对应 temporally projected elastic compliance contribution。

因此：

$$ \boxed{ W^{-1} = \text{effective compliance-like spatial operator} } $$

必须明确：

$$ \boxed{ W^{-1}\neq\mathbb C^{-1} } $$

也不能把：

$$ W $$

理解成：

- damaged tangent stiffness；
- material constitutive tangent；
- tower global stiffness matrix。

$W$ 仍首先是 continuum/material-point level 的 spatial operator。

---

# 22. 为什么记为 $W^{-1}$

Eq. (66) 中该 operator 作用在 stress spatial mode $\bar{\sigma}$ 上，并产生 strain-like quantity，因此从映射关系看具有 compliance-like 角色：

$$ \bar{\sigma} \xrightarrow{W^{-1}} \text{strain-like field} $$

因此论文自然记为：

$$ W^{-1} $$

相应的：

$$ W $$

在后续 Eq. (70) 中将表现为 stiffness-like operator。

---

# 23. $\bar{\delta}$ 的准确含义

定义：

$$ \boxed{ \bar{\delta}(x) = \int_0^T \bar{\Delta}(x,t)\lambda(t)\,dt } $$

它是：

$$ \boxed{ \text{temporal projection of the remaining LATIN search-direction defect} } $$

不是：

- Newton residual；
- equilibrium residual；
- plastic strain；
- independent damage mode；
- new constitutive variable。

它只是 Eq. (66) spatial problem 中的 known forcing/source field。

---

# 24. $\bar{\tilde{\varepsilon}}$ 在 Eq. (67) 时的身份

Eq. (67)：

$$ \bar{\tilde{\varepsilon}} = W^{-1}\bar{\sigma} - \bar{\delta} $$

此时首先只能称为：

$$ \boxed{ \text{auxiliary spatial strain-like field} } $$

它不是：

$$ \Delta\varepsilon' $$

不是：

$$ \Delta\dot{\varepsilon}' $$

不是：

$$ \bar{\varepsilon}^{e} $$

也不是最终 new PGD plastic spatial mode：

$$ \bar{\varepsilon}^{p} $$

最终 plastic spatial mode 要到后面的 Eq. (71) 才从 spatial solution 中恢复。

---

# 25. Eq. (67) 的代数反解形式

由：

$$ \bar{\tilde{\varepsilon}} = W^{-1}\bar{\sigma} - \bar{\delta} $$

得到：

$$ W^{-1}\bar{\sigma} = \bar{\tilde{\varepsilon}} + \bar{\delta} $$

因此：

$$ \boxed{ \bar{\sigma} = W\left( \bar{\tilde{\varepsilon}} + \bar{\delta} \right) } $$

这个形式非常重要，因为它已经呈现出：

$$ \boxed{ \text{stress} = \text{effective stiffness} \times \left( \text{compatible strain candidate} + \text{known source field} \right) } $$

的结构。

这正是 Eq. (70) 能转换成 standard FE weak form 的直接准备。

---

# 26. Eq. (68)：Eq. (66) 在 Eq. (67) 变量下的重写

将 Eq. (67) 代回 Eq. (66)，直接得到：

$$ \boxed{ \int_\Omega \bar{\tilde{\varepsilon}}:\bar{\sigma}^*\,d\Omega = 0, \qquad \forall\bar{\sigma}^*\in S_0 } \tag{68} $$

因此：

$$ \boxed{ \text{Eq. (68) 本身不是新 constitutive equation，而是 Eq. (66)+(67) 的直接重写} } $$

但 Eq. (68) 的意义非常关键，因为它赋予 $\bar{\tilde{\varepsilon}}$ kinematic admissibility。

---

# 27. 为什么 Eq. (68) 足以说明 $\bar{\tilde{\varepsilon}}$ compatible

Eq. (68) 表示：

$$ \boxed{ \bar{\tilde{\varepsilon}} \perp S_0 } $$

根据论文前面采用的 homogeneous kinematic-admissibility characterization：

$$ \boxed{ E_0 = S_0^\perp } $$

因此：

$$ \boxed{ \bar{\tilde{\varepsilon}}\in E_0 } $$

也就是说存在：

$$ \boxed{ \bar{\tilde{u}}\in U_0 } $$

使：

$$ \boxed{ \bar{\tilde{\varepsilon}} = \varepsilon(\bar{\tilde{u}}) } $$

所以 $\bar{\tilde{\varepsilon}}$ 的身份在 Eq. (68) 前后有一个重要变化：

Eq. (67) 时：

$$ \boxed{ \bar{\tilde{\varepsilon}} = \text{auxiliary strain-like field} } $$

Eq. (68) 后：

$$ \boxed{ \bar{\tilde{\varepsilon}} = \text{kinematically admissible auxiliary spatial strain field} } $$

---

# 28. Eq. (68) 只证明 displacement representation 的存在性

Eq. (68) 说明：

$$ \exists\bar{\tilde{u}}\in U_0 \quad \text{s.t.} \quad \bar{\tilde{\varepsilon}}=\varepsilon(\bar{\tilde{u}}) $$

但此时还没有实际求出 $\bar{\tilde{u}}$。

因此：

$$ \boxed{ \text{Eq. (68) establishes admissibility/existence, not yet the displacement solve} } $$

真正的 displacement-like FE equation 需要将 Eq. (67)、Eq. (68) 与 Eq. (69) 结合，得到下一阶段 Eq. (70)。

---

# 29. Eq. (69)：static admissibility to zero

原论文 Eq. (69)：

$$ \boxed{ \int_\Omega \bar{\sigma}:\varepsilon(\bar{u}^*)\,d\Omega = 0, \qquad \forall\bar{u}^*\in U_0 } \tag{69} $$

这一式负责：

$$ \boxed{ \bar{\sigma}\in S_0 } $$

即 new stress spatial mode 必须是 zero-load statically admissible stress field。

Eq. (68) 与 Eq. (69) 因此构成一对：

$$ \boxed{ \begin{aligned} \text{Eq. (68)} &: \bar{\tilde{\varepsilon}}\in E_0, \\ \text{Eq. (69)} &: \bar{\sigma}\in S_0. \end{aligned} } $$

这两式共同关闭 new spatial mode 的 global admissibility。

---

# 30. Eq. (69) 右端为什么同样为零

Eq. (69) 处理的也是 enrichment correction mode，而不是 total structural field。

真实结构外荷载已经由当前 global state 与 elastic initialization 等承担。

new enrichment stress spatial mode 只能作为 homogeneous correction：

$$ \boxed{ \text{no additional external load is attached to the new stress mode} } $$

因此：

$$ \boxed{ \text{Eq. (69) RHS}=0 } $$

不表示：

$$ \boxed{ \text{tower has no wind / wave / gravity loading} } $$

只表示：

$$ \boxed{ \text{new spatial correction mode must be self-equilibrated under zero incremental external loading} } $$

---

# 31. Eq. (69) 的 strong-form 解释

对：

$$ \int_\Omega \bar{\sigma}:\nabla^{sym}\bar{u}^*\,d\Omega = 0 $$

分部积分：

$$ -\int_\Omega(\nabla\cdot\bar{\sigma})\cdot\bar{u}^*\,d\Omega + \int_{\partial\Omega}(\bar{\sigma}n)\cdot\bar{u}^*\,dS = 0 $$

因此对应：

$$ \boxed{ \nabla\cdot\bar{\sigma} = 0 \quad \text{in }\Omega } $$

以及 correction traction：

$$ \boxed{ \bar{\sigma}n = 0 \quad \text{on Neumann boundary} } $$

在 prescribed-displacement boundary 上：

$$ \bar{u}^*=0 $$

因此允许 non-zero reaction traction。

---

# 32. 对固定基础 tower 的物理意义

对于 cantilever offshore wind turbine tower：

- 基底为 constrained DOFs；
- 塔身和塔顶为 free structural DOFs。

new PGD stress mode 可以引起：

- base shear correction；
- base bending-moment correction；
- base axial reaction correction。

这些都是允许的，因为 constrained DOFs 可以产生 reaction。

但是在 free DOFs 上必须满足：

$$ \boxed{ \text{zero structural residual force} } $$

这就是 Eq. (69) 对 fiber beam-column tower 最直观的物理解释。

---

# 33. Tower kinematic mapping

此前 tower spatial derivation 已建立：

$$ \vec{\varepsilon} = H\vec{U} $$

其中：

- $\vec{U}$：free structural DOF vector；
- $H$：nodal DOF → fiber material-point strain mapping；
- $\vec{\varepsilon}$：flatten 后的 fiber material-point strain vector。

对 coarse tower：

$$ N_e=10,\qquad N_g=2,\qquad N_f=16 $$

因此：

$$ \boxed{ N_q=N_eN_gN_f=320 } $$

所以：

$$ \vec{\varepsilon}\in\mathbb R^{320} $$

---

# 34. Tower quadrature metric

定义：

$$ M = \operatorname{diag}(v_q) $$

其中：

$$ v_q = A_{egf}w_gJ_e $$

连续内积：

$$ \int_\Omega a(x)b(x)\,d\Omega $$

离散为：

$$ \vec{a}^{\,T}M\vec{b} $$

因此 Eq. (68) 和 Eq. (69) 都可直接迁移到 fiber material-point quadrature space。

---

# 35. Tower Eq. (68) 的离散形式

Eq. (68)：

$$ \int_\Omega \bar{\tilde{\varepsilon}}:\bar{\sigma}^*\,d\Omega = 0 $$

离散后：

$$ \boxed{ (\vec{\bar{\sigma}}^{\,*})^T M \vec{\bar{\tilde{\varepsilon}}} = 0, \qquad \forall\vec{\bar{\sigma}}^{\,*}\in S_{0,\mathrm{tower}} } $$

而 zero-load statically admissible stress space：

$$ \boxed{ S_{0,\mathrm{tower}} = \ker(H^TM) } $$

因此：

$$ \vec{\bar{\tilde{\varepsilon}}} \perp_M \ker(H^TM) $$

---

# 36. Tower Eq. (68) 的线性代数证明

compatible fiber-strain space 为：

$$ \boxed{ E_{0,\mathrm{tower}} = \operatorname{Range}(H) } $$

因为：

$$ \vec{\varepsilon}=H\vec{U} $$

另一方面：

$$ H^TM\vec{\sigma}^*=0 $$

说明：

$$ \vec{\sigma}^* \perp_M \operatorname{Range}(H) $$

因此：

$$ \boxed{ S_{0,\mathrm{tower}} = \operatorname{Range}(H)^{\perp_M} } $$

Eq. (68) 又给出：

$$ \vec{\bar{\tilde{\varepsilon}}} \perp_M S_{0,\mathrm{tower}} $$

所以：

$$ \vec{\bar{\tilde{\varepsilon}}} \in \left[\operatorname{Range}(H)^{\perp_M}\right]^{\perp_M} $$

有限维空间中：

$$ \left(V^\perp\right)^\perp = V $$

于是：

$$ \boxed{ \vec{\bar{\tilde{\varepsilon}}}\in\operatorname{Range}(H) } $$

因此存在：

$$ \boxed{ \vec{\bar{\tilde{U}}} } $$

使：

$$ \boxed{ \vec{\bar{\tilde{\varepsilon}}} = H\vec{\bar{\tilde{U}}} } $$

这是原论文 Eq. (68) 在 tower 上的精确离散对应。

---

# 37. Tower Eq. (69) 的离散形式

任意 virtual displacement：

$$ \vec{U}^* $$

产生 virtual fiber strain：

$$ \vec{\varepsilon}^* = H\vec{U}^* $$

因此 Eq. (69)：

$$ \int_\Omega \bar{\sigma}:\varepsilon(\bar{u}^*)\,d\Omega = 0 $$

离散为：

$$ (\vec{\bar{\sigma}})^T M H\vec{U}^* = 0 $$

转置：

$$ (\vec{U}^*)^T H^TM\vec{\bar{\sigma}} = 0 $$

对所有 free virtual DOFs 任意成立，因此：

$$ \boxed{ H^TM\vec{\bar{\sigma}} = 0 } $$

这就是 tower zero-load static admissibility。

---

# 38. Eq. (69) 中可以真正定义 structural equilibrium residual

这里需要与前面的 $\bar{\Delta}$ 明确区分。

remaining LATIN defect：

$$ \bar{\Delta}_q(t) $$

是：

$$ \boxed{ \text{material-point descent-search-direction defect} } $$

而：

$$ H^TM\vec{\bar{\sigma}} $$

确实是 structural equilibrium residual。

因此可以定义：

$$ \boxed{ \vec{r}_{eq} = H^TM\vec{\bar{\sigma}} } $$

并要求：

$$ \boxed{ \vec{r}_{eq,F}=0 } $$

这两个“residual”不是同一个量，也不在同一个物理层级。

---

# 39. Tower Eq. (67) 的 material-point 形式

对 scalar fiber material：

$$ C_0=E_0 $$

定义：

$$ A_q = \left\langle H_{\sigma,q}\lambda^2 \right\rangle $$

$$ b = \left\langle \dot{\lambda}\lambda \right\rangle $$

$$ \bar{\delta}_q = \left\langle \bar{\Delta}_q\lambda \right\rangle $$

则：

$$ \boxed{ W_q^{-1} = A_q + \frac{b}{E_0} } $$

并：

$$ \boxed{ \bar{\tilde{\varepsilon}}_q = W_q^{-1}\bar{\sigma}_q - \bar{\delta}_q } $$

若 $W_q^{-1}$ 可逆，则：

$$ \boxed{ W_q = \frac{1}{A_q+b/E_0} } $$

因此每个 fiber material point 一般拥有自己的：

$$ W_q $$

不能把它实现成一个全塔统一 stiffness scalar。

---

# 40. Tower Eq. (67) 的向量形式

定义：

$$ D_{W^{-1}} = \operatorname{diag}(W_1^{-1},\ldots,W_{N_q}^{-1}) $$

以及：

$$ \vec{\bar{\delta}} = [\bar{\delta}_1,\ldots,\bar{\delta}_{N_q}]^T $$

则：

$$ \boxed{ \vec{\bar{\tilde{\varepsilon}}} = D_{W^{-1}}\vec{\bar{\sigma}} - \vec{\bar{\delta}} } $$

反解：

$$ \boxed{ \vec{\bar{\sigma}} = D_W\left(\vec{\bar{\tilde{\varepsilon}}}+\vec{\bar{\delta}}\right) } $$

其中：

$$ D_W = \operatorname{diag}(W_1,\ldots,W_{N_q}) $$

再利用 Eq. (68)：

$$ \vec{\bar{\tilde{\varepsilon}}}=H\vec{\bar{\tilde{U}}} $$

可预先得到：

$$ \boxed{ \vec{\bar{\sigma}} = D_W\left(H\vec{\bar{\tilde{U}}}+\vec{\bar{\delta}}\right) } $$

这一式正是下一阶段 Eq. (70) 的直接起点。

---

# 41. Eq. (68) 与 Eq. (69) 在 tower 上形成完整对偶结构

Eq. (68)：

$$ \boxed{ \vec{\bar{\tilde{\varepsilon}}} = H\vec{\bar{\tilde{U}}} } $$

负责：

$$ \boxed{ \text{kinematic compatibility} } $$

Eq. (69)：

$$ \boxed{ H^TM\vec{\bar{\sigma}} = 0 } $$

负责：

$$ \boxed{ \text{static equilibrium} } $$

所以 new enrichment spatial solution 必须同时满足：

$$ \boxed{ \begin{cases} \vec{\bar{\tilde{\varepsilon}}}=H\vec{\bar{\tilde{U}}}, \\ H^TM\vec{\bar{\sigma}}=0. \end{cases} } $$

这就是 original continuum admissibility 在 fiber beam-column tower 上的离散实现。

---

# 42. 为什么不能逐 fiber 独立求解

如果逐 fiber 直接写：

$$ W_q^{-1}\bar{\sigma}_q - \bar{\delta}_q = 0 $$

则：

$$ \bar{\sigma}_q = W_q\bar{\delta}_q $$

虽然每个 material point 都有一个 local solution，但一般：

$$ H^TM\vec{\bar{\sigma}}\neq0 $$

说明该 stress field 在结构层面不平衡。

同样，如果任意构造：

$$ \vec{\bar{\tilde{\varepsilon}}}\in\mathbb R^{N_q} $$

也不保证：

$$ \vec{\bar{\tilde{\varepsilon}}}\in\operatorname{Range}(H) $$

可能违反：

- beam kinematics；
- section plane-section assumption；
- element interpolation；
- nodal continuity；
- global displacement compatibility。

所以必须同时保留 Eq. (68) 与 Eq. (69)。

---

# 43. $H$ 与 $H_\sigma$ 必须严格区分

当前推导中同时出现两个完全不同的符号：

LATIN descent search direction：

$$ \boxed{ H_\sigma } $$

Tower kinematic mapping：

$$ \boxed{ H } $$

前者作用于 constitutive/search-direction relation；

后者完成：

$$ \vec{U}\rightarrow\vec{\varepsilon} $$

的结构运动学映射。

后续代码实现中建议采用清晰的变量名区分，例如：

```text
H_sigma   -> LATIN descent search-direction operator
H_kin     -> tower DOF-to-fiber-strain kinematic mapping
```

数学文档可以沿用当前符号，但代码层必须避免混淆。

---

# 44. Eq. (47)–(53) 中的 $C_0$-based equilibrium operator 与当前 $W$-based problem 的关系

此前已建立 reference tower equilibrium projection：

$$ \mathcal E_{\mathrm{tower}} = H(H^TMC_0H)^{-1}H^TMC_0 $$

以及：

$$ \sigma^{eq} = C_0(\mathcal E_{\mathrm{tower}}-I)q $$

其 reference global matrix 为：

$$ \boxed{ K^0 = H^TMC_0H } $$

但 Eq. (67)–(69) 中 spatial enrichment problem 出现的是：

$$ \boxed{ W_q } $$

而不是固定：

$$ C_0=E_0 $$

因此本阶段形成的重要判断是：

> 原来的 tower equilibrium operator 的**运动学映射、quadrature metric、global assembly 和 free-DOF equilibrium architecture**很可能可以复用，但其 material-point operator 在 Eq. (70) spatial enrichment solve 中不能未经推导就继续固定使用 $C_0$。

下一阶段需要严格检查是否应形成：

$$ \boxed{ K_W = H^TMD_WH } $$

这样的 effective spatial stiffness matrix。

这是 Eq. (70) 必须解决的核心问题之一。

---

# 45. 当前不能直接把 $\mathcal E_{\mathrm{tower}}$ 当成 Eq. (68) 的 $W$-projection

以前：

$$ \mathcal E_{\mathrm{tower}} $$

是基于：

$$ C_0 $$

metric 构造的 reference projection。

而当前 Eq. (67) 引入：

$$ W(x) $$

作为 temporally projected effective spatial operator。

虽然 compatible subspace 都是：

$$ \operatorname{Range}(H) $$

但不同 operator 对应的 projection metric 未必相同。

因此：

$$ \boxed{ \text{same compatible subspace} \neq \text{same projection operator} } $$

在 Eq. (70) 推导完成之前，不应直接把旧的 $C_0$-based projection 当成当前 enrichment spatial solve 的最终 operator。

---

# 46. 本阶段形成的完整 Eq. (65)–(69) 逻辑链

首先：

$$ \Delta\dot{\varepsilon}' = H_\sigma\lambda\bar{\sigma} + \dot{\lambda}\mathbb C^{-1}\bar{\sigma} - \bar{\Delta} $$

Eq. (65)：

$$ \boxed{ \int_{T\times\Omega}\Delta\dot{\varepsilon}':\sigma^*=0 } $$

要求 space-time kinematic admissibility。

取：

$$ \sigma^*=\lambda\bar{\sigma}^* $$

完成 temporal Galerkin projection。

Eq. (66)：

$$ \boxed{ \int_\Omega \left[ \langle H_\sigma\lambda^2\rangle\bar{\sigma} + \langle\lambda\dot{\lambda}\rangle C^{-1}\bar{\sigma} - \langle\bar{\Delta}\lambda\rangle \right]:\bar{\sigma}^*\,d\Omega=0 } $$

定义：

$$ W^{-1}=\langle H_\sigma\lambda^2\rangle+\langle\dot{\lambda}\lambda\rangle C^{-1} $$

$$ \bar{\delta}=\langle\bar{\Delta}\lambda\rangle $$

Eq. (67)：

$$ \boxed{ \bar{\tilde{\varepsilon}}=W^{-1}\bar{\sigma}-\bar{\delta} } $$

Eq. (68)：

$$ \boxed{ \bar{\tilde{\varepsilon}}\perp S_0 \Rightarrow \bar{\tilde{\varepsilon}}\in E_0 } $$

因此：

$$ \bar{\tilde{\varepsilon}}=\varepsilon(\bar{\tilde{u}}) $$

Eq. (69)：

$$ \boxed{ \bar{\sigma}\in S_0 } $$

因此 spatial enrichment mode 同时满足：

$$ \boxed{ \text{kinematic compatibility} + \text{static equilibrium} } $$

---

# 47. Tower migration 的完整对应链

连续 spatial coordinate：

$$ x $$

离散为：

$$ q=(e,g,f) $$

time-space defect：

$$ \bar{\Delta}(x,t) $$

变成：

$$ \bar{\Delta}_q(t) $$

temporal projection：

$$ \boxed{ \bar{\delta}_q=\int_0^T\bar{\Delta}_q(t)\lambda(t)\,dt } $$

effective compliance：

$$ \boxed{ W_q^{-1}=\int_0^T H_{\sigma,q}(t)\lambda^2(t)\,dt+\frac{1}{E_0}\int_0^T\dot{\lambda}(t)\lambda(t)\,dt } $$

auxiliary strain：

$$ \boxed{ \vec{\bar{\tilde{\varepsilon}}}=D_{W^{-1}}\vec{\bar{\sigma}}-\vec{\bar{\delta}} } $$

kinematic admissibility：

$$ \boxed{ \vec{\bar{\tilde{\varepsilon}}}=H\vec{\bar{\tilde{U}}} } $$

static admissibility：

$$ \boxed{ H^TM\vec{\bar{\sigma}}=0 } $$

stress relation：

$$ \boxed{ \vec{\bar{\sigma}}=D_W\left(H\vec{\bar{\tilde{U}}}+\vec{\bar{\delta}}\right) } $$

这些关系已经为 Eq. (70) 的 finite-element assembly 完整铺平道路。

---

# 48. 与原论文 $x-t$ 路线的一致性

到 Eq. (69) 为止，tower migration 仍然没有要求改变原论文的核心 PGD decomposition：

$$ \boxed{ \Delta\varepsilon^p(x,t)=\lambda(t)\bar{\varepsilon}^p(x) } $$

只是连续空间 $x$ 被 fiber beam-column discretization 表示为：

$$ \boxed{ x \rightarrow (e,g,f) } $$

因此目前仍然支持原定路线：

$$ \boxed{ \text{original }x-t\text{ LATIN-PGD} \rightarrow \text{fiber beam-column tower }x-t\text{ LATIN-PGD} } $$

不需要在当前阶段引入额外 PGD dimensions。

---

# 49. 本阶段明确不能做的事情

1. 不把 Eq. (65) 当成 Newton equilibrium equation；
2. 不把 $\bar{\Delta}$ 当成 structural force residual；
3. 不把 Eq. (66) 按 fiber 逐点强制为零；
4. 不把 $\bar{\delta}$ 当成 plastic strain；
5. 不把 $\bar{\tilde{\varepsilon}}$ 当成 final PGD plastic spatial mode；
6. 不把 $W^{-1}$ 当成原始 material compliance $C^{-1}$；
7. 不把 $W$ 当成 tower global stiffness matrix；
8. 不因为是 cyclic loading 就删除 $\langle\lambda\dot{\lambda}\rangle$；
9. 不把 $S_0$ 理解成所有 boundary traction 均为零；
10. 不把 constrained-DOF reactions 错误强制为零；
11. 不把旧 $C_0$-based equilibrium projection 未经推导直接用于新的 $W$-based enrichment problem；
12. 不把 Eq. (68) 简化成“任意 material-point strain vector”；
13. 不把 Eq. (69) 忽略，否则 new stress spatial mode 一般不满足 structural equilibrium；
14. 不在 Eq. (69) 前提前跳到 nodal solve；
15. 不在本阶段引入 $n-\tau-x$ 或 cycle-jump formulation。

---

# 50. 本阶段已经解决的问题

通过 Eq. (65)–(69) 的逐式推导，目前已经明确：

1. Eq. (65) 为什么是 kinematic admissibility；
2. 为什么 kinematic admissibility 用 $\sigma^*\in S_0$ 测试；
3. 为什么 Eq. (65) 右端为零；
4. 为什么 correction 的 homogeneous boundary condition 是关键；
5. 为什么 LATIN 在整个 $[0,T]\times\Omega$ 上施加 admissibility；
6. $S_0$ 的准确含义；
7. 为什么 constrained boundary 可以产生 reaction；
8. 为什么 spatial subproblem 中取 $\sigma^*=\lambda\bar{\sigma}^*$；
9. Eq. (66) 中 $\lambda^2$ 的来源；
10. Eq. (66) 中 $\lambda\dot{\lambda}$ 的来源；
11. Eq. (66) 中 $\bar{\Delta}\lambda$ 的来源；
12. 为什么 Eq. (66) 是 temporal Galerkin projection；
13. 为什么 Eq. (66) 不能逐 material point 强制为零；
14. Eq. (66) 与 Eq. (59) 的区别；
15. $W^{-1}$ 为什么是 effective compliance-like operator；
16. $W^{-1}$ 为什么不是 $C^{-1}$；
17. $W$ 为什么不是 global tower stiffness；
18. $\bar{\delta}$ 的准确含义；
19. $\bar{\tilde{\varepsilon}}$ 在 Eq. (67) 时为什么只是 strain-like field；
20. Eq. (68) 为什么能够证明 $\bar{\tilde{\varepsilon}}$ compatible；
21. 为什么 Eq. (68) 只建立 displacement representation 的存在性；
22. Eq. (69) 为什么是 static admissibility；
23. Eq. (69) 右端为什么为零；
24. 为什么 free DOFs 必须平衡而 constrained DOFs 可以有 reaction；
25. tower Eq. (68) 为什么对应 $\operatorname{Range}(H)$；
26. tower Eq. (69) 为什么对应 $H^TM\bar{\sigma}=0$；
27. $\bar{\Delta}$ 与 structural equilibrium residual 的区别；
28. 为什么之前建立的 $H$、$M$ 和 equilibrium architecture 可以继续复用；
29. 为什么 $C_0$-based projection 与当前 $W$-based spatial problem 仍需进一步区分；
30. 为什么 Eq. (70) 将是当前阶段真正进入 nodal spatial solve 的关键。

---

# 51. 当前尚未解决的问题

下一阶段 Eq. (70) 必须重点解决：

1. 如何从 Eq. (67)–(69) 严格得到 displacement-like weak form；
2. 为什么最终 unknown 可以选为 $\bar{\tilde{u}}$；
3. Eq. (70) 左端 stiffness-like bilinear form 的准确结构；
4. Eq. (70) 右端 projected-defect/source term 的准确符号；
5. 对 tower 是否得到：

$$ K_W = H^TMD_WH $$

6. $D_W$ 是否需要每次 enrichment fixed-point spatial iteration 重新构造；
7. 如何处理 $W_q^{-1}$ 接近零或非正定的数值问题；
8. 原 1D reproduction 中对 $H_\sigma$、$H_\beta$、$H_{\bar R}$ 的正定性经验如何迁移；
9. 旧 `equilibrium_operator.py` 是直接复用，还是只复用 assembly skeleton；
10. Eq. (70) 解出的 $\bar{\tilde U}$ 如何进一步恢复 $\bar{\sigma}$；
11. 再如何从 Eq. (71) 恢复真正的 new plastic-strain spatial mode $\bar{\varepsilon}^p$；
12. enrichment fixed-point spatial solve 与后续 temporal residual minimization 如何衔接。

---

# 52. 下一阶段工作入口

下一阶段严格从 Eq. (70) 开始，不直接跳到 Eq. (71)–(72)。

当前已经具备：

$$ \vec{\bar{\tilde{\varepsilon}}}=H\vec{\bar{\tilde U}} $$

$$ \vec{\bar{\sigma}}=D_W\left(H\vec{\bar{\tilde U}}+\vec{\bar{\delta}}\right) $$

以及：

$$ H^TM\vec{\bar{\sigma}}=0 $$

因此下一阶段的直接代入起点为：

$$ H^TM D_W\left(H\vec{\bar{\tilde U}}+\vec{\bar{\delta}}\right)=0 $$

形式上可整理为：

$$ \boxed{ H^TM D_W H\vec{\bar{\tilde U}} = -H^TM D_W\vec{\bar{\delta}} } $$

但这一式在本阶段只作为 Eq. (70) 推导的**入口提示**，尚未作为最终 tower Eq. (70) 正式确认。

下一阶段需要回到原论文 Eq. (70) 的连续 weak form，从原论文符号和符号方向重新逐项核对后，再确认 tower matrix form。

---

# 53. 阶段结论

本阶段的核心成果是完成了从 space-time enrichment correction 到 global-admissible spatial problem 的完整逻辑闭环：

$$ \boxed{ \bar{\Delta}(x,t) \xrightarrow{\lambda(t)\text{ temporal projection}} \bar{\delta}(x),\,W(x) \xrightarrow{\text{Eq. (67)}} \bar{\tilde{\varepsilon}}(x),\,\bar{\sigma}(x) \xrightarrow{\text{Eq. (68)}} \text{kinematic admissibility} \xrightarrow{\text{Eq. (69)}} \text{static admissibility} } $$

对 fiber beam-column tower，该结构直接对应：

$$ \boxed{ \bar{\Delta}_q(t) \rightarrow \bar{\delta}_q,\,W_q \rightarrow \vec{\bar{\tilde{\varepsilon}}}=H\vec{\bar{\tilde U}} \rightarrow H^TM\vec{\bar{\sigma}}=0 } $$

因此截至 Eq. (69)，原论文的 $x-t$ LATIN-PGD enrichment architecture 仍可以在不改变理论主线的条件下迁移到 tower fiber material-point space。

当前真正需要解决的下一关键问题已经非常明确：

$$ \boxed{ \text{如何把 }W_q\text{、}H\text{、}M\text{ 组装成 Eq. (70) 的 tower displacement-like spatial FE solve} } $$

这将决定我们之前的 fixed-reference $C_0$ equilibrium operator 在 new-mode enrichment 中究竟应该如何复用和扩展。
