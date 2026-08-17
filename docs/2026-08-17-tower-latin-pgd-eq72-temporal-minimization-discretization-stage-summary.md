# Tower LATIN-PGD Eq. (72) temporal minimisation、DG0 讨论与 1D 离散继承阶段总结

**日期：2026-08-17**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**当前研究路线：Bhattacharyya et al. 原论文 $x-t$ LATIN-PGD → 2D fiber beam-column offshore wind turbine tower**  
**阶段范围：原论文 Eq. (72) 的连续最小化结构、temporal weak form、连续强形式、DG0 信息边界、与现有 1D backward-Euler temporal update 的对应关系及后续研究决策**  
**上一阶段衔接：`2026-08-17-tower-latin-pgd-eq70-71-spatial-solve-plastic-mode-stage-summary.md`**  
**下一阶段：基于已验证 1D temporal discretization，推导 single-new-mode tower Eq. (72) 的离散标量更新方程**

---

# 1. 本阶段定位

上一阶段已经完成原论文 Eq. (70)–(71) 的 tower spatial half-step。

固定当前 temporal guess $\lambda^{(k)}(t)$ 后，已经得到：

$$ \boxed{ \lambda^{(k)}(t) \rightarrow W^{(k)},\bar\delta^{(k)} \rightarrow \bar{\tilde U}^{(k+1)} \rightarrow \bar{\tilde\varepsilon}^{(k+1)},\bar\sigma^{(k+1)} \rightarrow \bar\varepsilon^{p,(k+1)} } $$

其中 tower Eq. (70) 已离散为：

$$ \boxed{ H^TMD_WH\vec{\bar{\tilde U}}=-H^TMD_W\vec{\bar\delta} } $$

Eq. (71) 已恢复真正的 plastic-strain spatial mode：

$$ \boxed{ \vec{\bar\varepsilon}^{\,p}=\frac{1}{a}H\vec{\bar{\tilde U}}-C_0^{-1}D_W\left(H\vec{\bar{\tilde U}}+\vec{\bar\delta}\right) } $$

其中：

$$ \boxed{ a=\left\langle\dot\lambda\lambda\right\rangle } $$

因此 Eq. (71) 结束后，new-mode enrichment 仍未闭合。

当前缺失的是：

$$ \boxed{ \bar\varepsilon^{p,(k+1)} \rightarrow \lambda^{(k+1)}(t) } $$

这正是原论文 Eq. (72) 的任务。

本阶段没有修改代码，而是完成 Eq. (72) 的数学解释、tower 映射、时间离散信息边界和下一步离散策略决策。

---

# 2. 原论文 Eq. (72) 的理论地位

原论文在 Eq. (71) 得到新的 plastic-strain spatial function 后，明确说明：

> Then the time function $\lambda_{m+1}$ is calculated using a minimisation technique similar to the update stage.

因此 Eq. (72) 与 Eq. (58)–(59) 的 old-basis temporal update 属于同一类 mechanical-residual minimisation，只是作用对象不同：

- Eq. (58)–(59)：已有 $m$ 个 spatial modes 固定，联合更新其 temporal coefficients；
- Eq. (72)：new enrichment spatial mode 固定，只更新新 pair 的 temporal function。

所以 Eq. (72) 不是新的材料方程、平衡方程或动力学方程，而是 new PGD pair 的 temporal half-step。

---

# 3. Eq. (72) 必须回到 Eq. (61) 理解

原论文 enrichment search-direction relation 为：

$$ \boxed{ \Delta\dot\varepsilon^p_{i+1}-H_\sigma\Delta\sigma'_{i+1}+\bar\Delta_{i+1}=0 } \tag{61} $$

new rank-one pair 写为：

$$ \boxed{ \Delta\dot\varepsilon^p_{i+1}=\dot\lambda_{m+1}(t)\bar\varepsilon^p_{m+1}(x) } $$

以及：

$$ \boxed{ \Delta\sigma'_{i+1}=\lambda_{m+1}(t)\bar\sigma_{m+1}(x) } $$

其中：

$$ \boxed{ \bar\sigma_{m+1}=\mathcal C\bar\varepsilon^p_{m+1} } $$

这里 $\mathcal C$ 表示前面 Eq. (52) 所建立的 plastic-strain spatial mode 到 equilibrated stress spatial mode 的线性 operator，本阶段使用 $\mathcal C$ 只是为了避免与 reference elastic tensor $C_0$ 的记号混淆。

代入 Eq. (61)：

$$ \boxed{ \dot\lambda_{m+1}\bar\varepsilon^p_{m+1}-H_\sigma\lambda_{m+1}\bar\sigma_{m+1}+\bar\Delta_{i+1}=0 } $$

因此 Eq. (72) 的 minimisation 对象直接就是 Eq. (61) 的剩余 mechanical search-direction defect。

---

# 4. Eq. (72) 的 residual 定义

固定 Eq. (71) 得到的 spatial pair 后，定义：

$$ \boxed{ r(x,t;\lambda)=\dot\lambda(t)\bar\varepsilon^p(x)-H_\sigma(x,t)\lambda(t)\bar\sigma(x)+\bar\Delta(x,t) } $$

于是原论文 Eq. (72) 可写为：

$$ \boxed{ \lambda=\arg\min_\lambda\left\|r(\lambda)\right\|_{H_\sigma^{-1}} } \tag{72} $$

其核心含义是：

$$ \boxed{ \text{在固定 new spatial mode 后，寻找最合适的 temporal amplitude，使 Eq. (61) 的 residual 尽可能小} } $$

---

# 5. Eq. (72) 最小化的不是哪些量

Eq. (72) 不是在最小化：

- ordinary Newton residual；
- tower nodal equilibrium residual；
- displacement error；
- stress error 本身；
- damage error；
- external wind-load mismatch；
- POD/SVD reconstruction error；
- total structural potential energy。

它最小化的是：

$$ \boxed{ \text{LATIN descent search-direction mechanical residual} } $$

这一区分对未来 tower implementation 很重要，避免将 Eq. (72) 误写成普通 structural load-step equation。

---

# 6. 为什么使用 $H_\sigma^{-1}$ metric

对于 residual field $r(x,t)$，原论文采用与 Eq. (59) 相同的 search-direction metric：

$$ \boxed{ \left\|r\right\|_{H_\sigma^{-1}}^2=\int_0^T\int_\Omega r:H_\sigma^{-1}r\,d\Omega\,dt } $$

因此 $H_\sigma^{-1}$ 不是随意选择的数值权重，而是与 LATIN descent search direction 相匹配的 mechanical metric。

这与当前 1D reproduction 中采用的 $H_\sigma^{-1}$ weighted residual norm 在结构上保持一致。

---

# 7. 为什么 Eq. (71) 后必须重新更新 $\lambda$

Eq. (70)–(71) 的 spatial half-step 使用一个暂时固定的 temporal guess：

$$ \boxed{ \lambda^{(k)}(t) } $$

在此条件下求得：

$$ \boxed{ \bar\varepsilon^{p,(k+1)}(x) } $$

这个 spatial function 是在当前 $\lambda^{(k)}$ 条件下得到的最佳空间方向。

但是 spatial function 更新后，旧 temporal function 一般不再是新 spatial mode 对应的最佳 temporal amplitude，因此必须再求：

$$ \boxed{ \bar\varepsilon^{p,(k+1)} \rightarrow \lambda^{(k+1)} } $$

于是 new-mode fixed point 为：

$$ \boxed{ \lambda^{(k)} \rightarrow \bar\varepsilon^{p,(k+1)} \rightarrow \lambda^{(k+1)} \rightarrow \bar\varepsilon^{p,(k+2)} \rightarrow \cdots } $$

因此 Eq. (70)–(71) 与 Eq. (72) 必须作为一对 alternating half-steps 理解。

---

# 8. 原论文 enrichment 的 hybrid strategy

原论文明确将 new-mode enrichment 描述为 hybrid strategy：

- spatial function：Galerkin formulation；
- temporal function：mechanical residual minimisation。

所以：

$$ \boxed{ \text{Eq. (65)–(71)}=\text{spatial Galerkin half-step} } $$

$$ \boxed{ \text{Eq. (72)}=\text{temporal minimisation half-step} } $$

这也再次确认：未来 tower implementation 不应机械复制当前 1D `_solve_spatial_function()` 的 fixed-temporal least-squares spatial solve，而应继续优先采用已经推导完成的 explicit $W$-based Eq. (70)–(71) route。

---

# 9. Eq. (72) 三个 residual contributions 的物理/算法含义

将 residual 写成：

$$ \boxed{ r=\underbrace{\dot\lambda\bar\varepsilon^p}_{\text{new plastic-rate correction}}-\underbrace{H_\sigma\lambda\bar\sigma}_{\text{new stress contribution in search direction}}+\underbrace{\bar\Delta}_{\text{remaining shifted LATIN defect}} } $$

第一项：

$$ \boxed{ \dot\lambda\bar\varepsilon^p=\Delta\dot\varepsilon^p } $$

表示 new rank-one pair 对 plastic-strain-rate correction 的贡献。

第二项：

$$ \boxed{ H_\sigma\lambda\bar\sigma } $$

表示 corresponding equilibrated stress correction 在 LATIN descent relation 中的贡献。

第三项：

$$ \boxed{ \bar\Delta } $$

来自 update stage 后已有 PGD basis 尚未消除的 shifted defect。

所以 Eq. (72) 的目标可以概括为：

$$ \boxed{ \text{选择 }\lambda(t)\text{，使 new mode 最大程度抵消 }\bar\Delta(x,t) } $$

---

# 10. Eq. (72) 中不直接出现真实风载并不意味着载荷被忽略

Eq. (72) 中没有直接出现：

$$ F_{\mathrm{wind}}(t),\qquad F_{\mathrm{wave}}(t) $$

原因是物理荷载已经通过 elastic initialisation、current LATIN state、local stage state 和 updated global state 等信息进入 shifted defect：

$$ \boxed{ \bar\Delta(x,t) } $$

因此 Eq. (72) 解决的是：

$$ \boxed{ \text{当前 LATIN iteration 中的 correction problem} } $$

而不是重新求解一遍物理外载问题。

未来 tower 代码中不应把 Eq. (72) 的 forcing 写成 tower top load history，而应从 current shifted defect 构造。

---

# 11. 将 Eq. (72) 写成目标泛函

最小化范数与最小化范数平方等价，因此定义：

$$ \boxed{ J[\lambda]=\frac12\int_0^T\int_\Omega r:H_\sigma^{-1}r\,d\Omega\,dt } $$

在当前 temporal half-step 中：

- $\bar\varepsilon^p(x)$ 固定；
- $\bar\sigma(x)$ 固定；
- $H_\sigma(x,t)$ 已知；
- $\bar\Delta(x,t)$ 已知；
- 唯一 unknown 为 $\lambda(t)$。

因此 Eq. (72) 是 pure temporal optimisation problem。

---

# 12. 对 temporal function 做任意 variation

令：

$$ \boxed{ \lambda(t)\rightarrow\lambda(t)+\epsilon\eta(t) } $$

则：

$$ \boxed{ \dot\lambda(t)\rightarrow\dot\lambda(t)+\epsilon\dot\eta(t) } $$

因此 residual variation 为：

$$ \boxed{ \delta r=\dot\eta\bar\varepsilon^p-H_\sigma\eta\bar\sigma } $$

这里 $H_\sigma$ 在当前 temporal minimisation 中视为已知 search-direction field，不对 $\lambda$ 做变分。

---

# 13. Eq. (72) 的第一变分

由：

$$ \boxed{ J=\frac12\int r:H_\sigma^{-1}r } $$

可得：

$$ \boxed{ \delta J=\int_0^T\int_\Omega\left(\dot\eta\bar\varepsilon^p-H_\sigma\eta\bar\sigma\right):H_\sigma^{-1}r\,d\Omega\,dt } $$

最优 temporal function 满足：

$$ \boxed{ \delta J=0,\qquad\forall\eta } $$

因此得到 temporal weak optimality condition：

$$ \boxed{ \int_0^T\int_\Omega\left(\dot\eta\bar\varepsilon^p-H_\sigma\eta\bar\sigma\right):H_\sigma^{-1}\left(\dot\lambda\bar\varepsilon^p-H_\sigma\lambda\bar\sigma+\bar\Delta\right)d\Omega\,dt=0 } $$

这一步是从 Eq. (72) 目标泛函直接得到的数学推导，不是原论文另行给出的新公式。

---

# 14. 空间 contraction 的第一组量

由于 $H_\sigma$ 为对称 search-direction operator，有：

$$ \boxed{ (H_\sigma\bar\sigma):H_\sigma^{-1}r=\bar\sigma:r } $$

因此 temporal weak form 可写成：

$$ \boxed{ \int_0^T\left[\dot\eta\int_\Omega\bar\varepsilon^p:H_\sigma^{-1}r\,d\Omega-\eta\int_\Omega\bar\sigma:r\,d\Omega\right]dt=0 } $$

这说明 fixed spatial mode 后，所有空间信息可以先被积分 contraction，再形成纯时间问题。

---

# 15. 定义 temporal coefficient $A(t)$

定义：

$$ \boxed{ A(t)=\int_\Omega\bar\varepsilon^p:H_\sigma^{-1}(x,t)\bar\varepsilon^p\,d\Omega } $$

$A(t)$ 衡量 plastic spatial mode 在 $H_\sigma^{-1}$ metric 下的权重。

它控制 temporal rate $\dot\lambda$ 在 Eq. (72) objective 中的贡献。

---

# 16. 定义 coupling coefficient $B$

定义：

$$ \boxed{ B=\int_\Omega\bar\varepsilon^p:\bar\sigma\,d\Omega } $$

在当前 single-mode temporal half-step 中，$\bar\varepsilon^p$ 和 $\bar\sigma$ 都是固定 spatial functions，所以 $B$ 不随时间变化。

$B$ 表示 plastic-strain spatial mode 与 corresponding equilibrated stress spatial mode 之间的 spatial coupling。

---

# 17. 定义 temporal coefficient $C(t)$

定义：

$$ \boxed{ C(t)=\int_\Omega\bar\sigma:H_\sigma(x,t)\bar\sigma\,d\Omega } $$

它衡量 stress spatial mode 在 $H_\sigma$ metric 下的大小，控制 temporal amplitude $\lambda(t)$ 的 stress-side contribution。

---

# 18. 定义 defect projections $d(t)$ 与 $e(t)$

定义：

$$ \boxed{ d(t)=\int_\Omega\bar\varepsilon^p:H_\sigma^{-1}\bar\Delta\,d\Omega } $$

以及：

$$ \boxed{ e(t)=\int_\Omega\bar\sigma:\bar\Delta\,d\Omega } $$

这两个量是 remaining shifted LATIN defect 在当前 new spatial pair 上的 temporal projections。

因此真正驱动 Eq. (72) temporal update 的不是直接外荷载，而是：

$$ \boxed{ d(t),\qquad e(t) } $$

---

# 19. Eq. (72) 的纯时间 weak form

将 residual 展开并使用上述定义，可得到：

$$ \boxed{ \int_0^T\left[\dot\eta\left(A\dot\lambda-B\lambda+d\right)-\eta\left(B\dot\lambda-C\lambda+e\right)\right]dt=0,\qquad\forall\eta } $$

这一式是本阶段最重要的 continuous temporal weak form。

它表明原本的 $(x,t)$ space-time minimisation 在 fixed spatial mode 条件下已经 contraction 为一个仅依赖：

$$ \boxed{ A(t),\quad B,\quad C(t),\quad d(t),\quad e(t) } $$

的 temporal reduced problem。

---

# 20. 为什么原论文说 temporal minimisation 会产生 differential equation

定义：

$$ \boxed{ Q(t)=A(t)\dot\lambda(t)-B\lambda(t)+d(t) } $$

以及：

$$ \boxed{ R(t)=B\dot\lambda(t)-C(t)\lambda(t)+e(t) } $$

则 weak form 为：

$$ \boxed{ \int_0^T\left[\dot\eta Q-\eta R\right]dt=0 } $$

对第一项做时间分部积分：

$$ \boxed{ \int_0^T\dot\eta Q\,dt=\left[\eta Q\right]_0^T-\int_0^T\eta\dot Q\,dt } $$

因此：

$$ \boxed{ \left[\eta Q\right]_0^T-\int_0^T\eta\left(\dot Q+R\right)dt=0 } $$

如果只看时间区间内部，由 $\eta$ 的任意性得到：

$$ \boxed{ \dot Q+R=0 } $$

这解释了为什么 Eq. (72) continuous minimisation 会产生 temporal differential equation。

---

# 21. Continuous interior strong form

代入 $Q$ 和 $R$：

$$ \boxed{ \frac{d}{dt}\left(A\dot\lambda-B\lambda+d\right)+B\dot\lambda-C\lambda+e=0 } $$

由于当前 $B$ 不随时间变化：

$$ \boxed{ \frac{d}{dt}(-B\lambda)=-B\dot\lambda } $$

与外部的 $+B\dot\lambda$ 正好抵消，因此 interior equation 可写为：

$$ \boxed{ A(t)\ddot\lambda(t)+\dot A(t)\dot\lambda(t)-C(t)\lambda(t)+\dot d(t)+e(t)=0 } $$

或：

$$ \boxed{ \frac{d}{dt}\left[A(t)\dot\lambda(t)+d(t)\right]-C(t)\lambda(t)+e(t)=0 } $$

这说明 Eq. (72) 的 continuous Euler–Lagrange equation 一般包含 $\ddot\lambda$。

---

# 22. 为什么出现 $\ddot\lambda$ 并不奇怪

Eq. (72) 的 objective 显式依赖：

$$ \boxed{ \lambda,\qquad\dot\lambda } $$

因此连续 calculus of variations 的 Euler–Lagrange equation 会对 $\partial J/\partial\dot\lambda$ 再做一次时间导数，产生 $\ddot\lambda$。

所以不能把 Eq. (72) 的 continuous problem 直接理解为普通的一阶 initial-value ODE。

---

# 23. Temporal endpoint condition 目前不能自行补全

分部积分后仍有：

$$ \boxed{ \left[\eta Q\right]_0^T } $$

若在普通连续变分问题中假定两个端点 variation 都自由，则会自然得到相应 natural endpoint conditions。

但是原论文 Eq. (72) 附近没有给出：

- $\lambda(0)$ 是否 prescribed；
- $\lambda(T)$ 是否 free；
- temporal traces 如何定义；
- DG0 jump terms 如何写；
- first time slab 如何处理；
- endpoint conditions 如何与 DG formulation 配合。

因此当前不能未经作者离散细节证明就直接把：

$$ \lambda(0)=0 $$

或：

$$ Q(0)=Q(T)=0 $$

写入 tower algorithm。

---

# 24. 为什么也不能简单由 correction initial condition 推出单个 mode 的 $\lambda(0)=0$

虽然最终 LATIN correction 必须与整体 initial admissibility 相容，但 Eq. (72) 中的 $\lambda_{m+1}$ 是一个单独的 enrichment mode temporal coefficient。

new spatial mode 后续还要：

- 与已有 spatial basis 正交化；
- 进行 scaling/normalisation；
- 重新更新全部 temporal functions；
- 可能因 modified temporal norm 太小而拒绝。

因此：

$$ \boxed{ \text{最终 correction 的初始可容许性} \not\Rightarrow \text{无需推导即可单独指定 new mode 的 }\lambda_{m+1}(0) } $$

在 DG temporal representation 下还涉及 left/right trace convention，所以这一问题必须留给明确的时间离散 formulation。

---

# 25. 原论文关于 DG0 能确认到什么程度

原论文在 Eq. (59) 后明确说明：

- temporal minimisation gives a multi-variable differential equation；
- 该 differential equation 使用 discontinuous Galerkin method of order zero 求解；
- 更详细的 DG minimisation treatment 指向早期 LATIN 文献。

Eq. (72) 又明确说 new-mode temporal function 使用与 update stage similar minimisation technique。

因此可以确认：

$$ \boxed{ \text{Eq. (59)/(72) temporal minimisation} \rightarrow \text{DG0 treatment in the original formulation} } $$

但当前原论文正文没有给出 exact DG0 algebra。

---

# 26. 当前阶段不能声称已经恢复 original DG0 的 exact formulation

截至本阶段，我们仍未从当前可直接核实的原始论文正文中得到：

- exact temporal trial space；
- exact test space；
- left/right trace definitions；
- jump sign convention；
- jump bilinear form；
- endpoint treatment；
- first-slab algebra；
- final discrete matrix equation。

因此未来阶段总结和论文写作中不能写：

$$ \boxed{ \text{current backward Euler implementation is exactly the original DG0} } $$

除非后续找到了作者更完整的离散公式并完成逐项证明。

---

# 27. Generic DG0 understanding 与 original-paper evidence 必须区分

从一般 discontinuous Galerkin 时间离散的数学概念可理解：

- zero-order temporal approximation 通常意味着 time slab 内采用 piecewise constant representation；
- 相邻 time slabs 之间允许 discontinuity；
- temporal evolution 可通过 weak derivative 和 interface/jump contributions 传递。

但这些属于 generic DG interpretation。

当前原论文正文只明确说使用 DG0，没有把上述 exact algebra 写出来。

因此项目内必须区分：

$$ \boxed{ \text{paper-confirmed information} } $$

和：

$$ \boxed{ \text{generic DG interpretation / our mathematical inference} } $$

不能把后者写成作者原式。

---

# 28. Tower Eq. (72) 的 material-point residual

当前 tower spatial coordinate 已经离散为：

$$ \boxed{ q=(e,g,f) } $$

定义：

$$ \boxed{ \vec p=\vec{\bar\varepsilon}^{\,p}_{m+1} } $$

$$ \boxed{ \vec s=\vec{\bar\sigma}_{m+1} } $$

以及 shifted defect：

$$ \boxed{ \vec\Delta(t)=\vec{\bar\Delta}_{i+1}(t) } $$

search-direction diagonal operator：

$$ \boxed{ D_H(t)=\operatorname{diag}(H_{\sigma,q}(t)) } $$

则 tower Eq. (72) residual 为：

$$ \boxed{ \vec r(t)=\dot\lambda(t)\vec p-D_H(t)\lambda(t)\vec s+\vec\Delta(t) } $$

这保持了原论文 Eq. (72) 的完整 $x-t$ structure，只把 continuum spatial coordinate 替换为 fiber material-point discretisation。

---

# 29. Tower quadrature metric

此前已经定义：

$$ \boxed{ M=\operatorname{diag}(v_q) } $$

其中：

$$ \boxed{ v_q=A_{egf}w_gJ_e } $$

因此 continuum spatial integral：

$$ \int_\Omega a:b\,d\Omega $$

离散后由 fiber quadrature weighted inner product 表示。

Eq. (72) objective 因而变成：

$$ \boxed{ J[\lambda]=\frac12\int_0^T\vec r(t)^TMD_H^{-1}(t)\vec r(t)\,dt } $$

这一步不需要任何新的 cycle-phase-space separation。

---

# 30. Tower 中的 temporal coefficients

continuum coefficient $A(t)$ 对应：

$$ \boxed{ A(t)=\vec p^{\,T}MD_H^{-1}(t)\vec p } $$

coupling coefficient：

$$ \boxed{ B=\vec p^{\,T}M\vec s } $$

stress coefficient：

$$ \boxed{ C(t)=\vec s^{\,T}MD_H(t)\vec s } $$

defect projection：

$$ \boxed{ d(t)=\vec p^{\,T}MD_H^{-1}(t)\vec\Delta(t) } $$

以及：

$$ \boxed{ e(t)=\vec s^{\,T}M\vec\Delta(t) } $$

因此 Eq. (72) 固定 spatial mode 后仍然只是 pure temporal reduced problem。

---

# 31. Eq. (72) 不需要每个时间点重新求 tower FE system

Eq. (70) 是 spatial FE solve：

$$ \boxed{ K_W\vec{\bar{\tilde U}}=\vec f_\delta } $$

但 Eq. (72) 固定 spatial mode 后，所有 tower spatial information 已经 contraction 为：

$$ \boxed{ A(t),\ B,\ C(t),\ d(t),\ e(t) } $$

所以 Eq. (72) 不应在每个时间点重新调用一个完整 tower structural FE solve。

这是 original $x-t$ PGD 降维逻辑在 tower 中得以保持的重要体现。

---

# 32. 当前 1D `pgd_time_update.py` 的实际离散结构

现有 1D implementation 固定 spatial basis 后，使用 residual：

$$ \boxed{ r=P\dot{\lambda}-H_\sigma S\lambda-f } $$

其中：

- $P$：spatial plastic-strain basis matrix；
- $S$：spatial equilibrated-stress basis matrix；
- $f$：known plastic forcing；
- temporal coefficients 为 unknown。

当前代码对 temporal derivative 采用 backward difference：

$$ \boxed{ \dot\lambda_n\approx\frac{\lambda_n-\lambda_{n-1}}{\Delta t_n} } $$

随后在每个 discrete time level 上求 $H_\sigma^{-1}$ weighted least-squares problem。

---

# 33. 当前 1D temporal update 与 Eq. (72) 一致的部分

现有 1D `pgd_time_update.py` 与原论文 Eq. (72) 在以下核心结构上保持一致：

1. fixed spatial mode / basis；
2. unknown 为 temporal coefficient；
3. mechanical residual 中包含 plastic-rate term；
4. mechanical residual 中包含 $-H_\sigma$ weighted stress correction；
5. residual 使用 $H_\sigma^{-1}$ metric；
6. shifted defect 可通过 forcing sign convention 对应；
7. temporal update 完成后可重新构造 complete residual field。

所以它不是与 Eq. (72) 无关的经验算法，而是对同一 Eq. (61)/(72) residual structure 的一种离散处理。

---

# 34. Eq. (72) 与当前 1D forcing 的 sign convention

原论文 Eq. (72) residual：

$$ \boxed{ r=\dot\lambda\bar\varepsilon^p-H_\sigma\lambda\bar\sigma+\bar\Delta } $$

当前 1D code residual：

$$ \boxed{ r=P\dot\lambda-H_\sigma S\lambda-f } $$

两者对应条件为：

$$ \boxed{ f=-\bar\Delta } $$

当前 `pgd_enrichment.py` 在 new-mode temporal update 中确实将当前 residual 以负号作为 forcing 传入 `update_pgd_time_functions()`。

因此 new-mode temporal update 的 sign structure 与 Eq. (72) 是一致的。

---

# 35. 当前 1D backward-Euler route 与 original DG0 尚未证明严格等价

现有 1D implementation 的实际顺序是：

$$ \boxed{ \text{先离散 }\dot\lambda\text{ 为 backward difference} \rightarrow \text{再做 discrete weighted least squares} } $$

而从 Eq. (72) continuous functional 出发的理论路线是：

$$ \boxed{ \text{continuous minimisation} \rightarrow \text{temporal weak/differential problem} \rightarrow \text{original paper says DG0} } $$

由于 original paper 正文没有给出 exact DG0 jump algebra，目前还不能证明这两个离散路径严格产生完全相同的 algebraic system。

因此项目中必须避免写：

$$ \boxed{ \text{backward Euler}=\text{original DG0 exactly} } $$

---

# 36. 对当前 `pgd_time_update.py` docstring 的理解需要更谨慎

当前代码说明中将 backward-Euler treatment 描述为 original zero-order discontinuous Galerkin treatment 的 discrete counterpart。

现阶段更严谨的研究表述应是：

$$ \boxed{ \text{a backward-Euler-discretised }H_\sigma^{-1}\text{-weighted temporal residual minimisation consistent with Eq. (72)} } $$

而不是未经进一步证明就声称已经严格恢复作者 DG0 algebra。

当前阶段不立即修改 docstring，先将理论边界记录清楚，未来若正式重构 tower temporal solver 时再统一处理代码说明。

---

# 37. 为什么不能因为 DG0 细节未恢复就阻塞 tower 研究

当前已经确定的原论文核心理论结构包括：

$$ \boxed{ \bar\Delta } $$

$$ \boxed{ H_\sigma } $$

$$ \boxed{ \bar\varepsilon^p } $$

$$ \boxed{ \bar\sigma } $$

$$ \boxed{ H_\sigma^{-1}\text{ residual minimisation} } $$

这些决定了 Eq. (72) 的物理和数学核心。

DG0 只是 temporal discretization 层的进一步具体选择。

与此同时，1D three-material bar reproduction 已经用当前 backward-Euler temporal update 得到稳定、准确的复现结果。

因此若作者 DG0 exact algebra 最终仍无法恢复，不应因此停止 tower migration。

---

# 38. 当前形成的正式研究决策

如果后续能够找到 original DG0 的完整 discrete formulation，则：

1. 逐项比较 original DG0 与 current backward-Euler implementation；
2. 检查两者是否 algebraically equivalent 或仅 structurally similar；
3. 再决定 tower temporal solver 是否需要调整。

如果 original DG0 exact algebra 仍无法可靠恢复，则第一版 tower implementation 正式采用：

$$ \boxed{ \text{original-paper Eq. (72) continuous structure}+\text{validated 1D backward-Euler temporal discretization} } $$

即：

- 理论 residual 不改；
- $H_\sigma^{-1}$ metric 不改；
- fixed-spatial / updated-temporal architecture 不改；
- 只在时间离散层沿用已经在 1D benchmark 验证过的 backward-Euler weighted least-squares implementation。

---

# 39. 这一研究决策不等于“放弃原论文”

必须区分：

$$ \boxed{ \text{continuous algorithmic structure} } $$

和：

$$ \boxed{ \text{specific temporal discretization} } $$

当前决定保持 Bhattacharyya Eq. (61)–(72) 的 continuous LATIN-PGD architecture，只在作者正文没有给足离散细节的地方，采用已验证的项目内 numerical discretization。

因此该路线仍然属于：

$$ \boxed{ x-t\text{ LATIN-PGD tower migration} } $$

而不是改成新的 PGD theory。

---

# 40. Tower 第一版 Eq. (72) 的预定离散起点

固定 Eq. (71) 得到的：

$$ \boxed{ \vec p=\vec{\bar\varepsilon}^{\,p}_{m+1},\qquad\vec s=\vec{\bar\sigma}_{m+1} } $$

在 discrete time level $t_n$，tower Eq. (72) residual 为：

$$ \boxed{ \vec r_n=\dot\lambda_n\vec p-D_{H,n}\lambda_n\vec s+\vec\Delta_n } $$

若采用 1D 已验证的 backward-Euler treatment：

$$ \boxed{ \dot\lambda_n\approx\frac{\lambda_n-\lambda_{n-1}}{\Delta t_n} } $$

则：

$$ \boxed{ \vec r_n=\left(\frac{\vec p}{\Delta t_n}-D_{H,n}\vec s\right)\lambda_n-\frac{\vec p}{\Delta t_n}\lambda_{n-1}+\vec\Delta_n } $$

下一阶段将从这一式严格推导 single-new-mode tower scalar weighted least-squares update，不在本阶段继续展开。

---

# 41. 下一阶段需要重点检查的 sign convention

当前应始终保持：

$$ \boxed{ r_n=\text{new-mode correction}+\bar\Delta_n } $$

如果沿用 `pgd_time_update.py` 的 forcing-style interface，则需要：

$$ \boxed{ f_n=-\bar\Delta_n } $$

未来 tower code 不能因为变量命名变化而把 Eq. (72) 的 $+\bar\Delta$ 错写成相反符号。

---

# 42. 下一阶段需要重点检查的 spatial weighting

1D bar 中 spatial integration weight 为 element volume：

$$ \boxed{ v_e=A L_e } $$

Tower 中必须替换为 fiber quadrature metric：

$$ \boxed{ v_q=A_{egf}w_gJ_e } $$

因此 discrete Eq. (72) 的 weighted least-squares metric 应基于：

$$ \boxed{ MD_{H,n}^{-1} } $$

而不是 ordinary Euclidean norm。

这也是 1D temporal algorithm 迁移至 tower 时最主要的 spatial generalisation。

---

# 43. 下一阶段需要重点检查的 temporal initial treatment

由于 original DG0 exact endpoint treatment 尚未恢复，第一版 backward-Euler tower discretization 需要明确：

1. $\lambda_0$ 如何定义；
2. initial temporal rate 如何处理；
3. 是否完全继承现有 `pgd_time_update.py` 对 first time point 的处理；
4. new-mode enrichment 与 old-basis update 在 $t=0$ 是否应采用同一 convention；
5. 该 convention 是否保持 correction initial admissibility。

这部分必须通过 1D 代码行为和 tower consistency tests 验证，不能只靠形式推理。

---

# 44. 下一阶段需要继承的 1D numerical safeguards

即使 Eq. (72) continuous theory 不变，tower implementation 仍应继承已经在 1D reproduction 中证明有价值的 numerical controls，包括：

- strictly positive $H_\sigma$；
- weighted least-squares conditioning monitor；
- finite-value checks；
- minimum spatial norm；
- fixed-point maximum iteration limit；
- pair-change convergence criterion；
- candidate mode residual-reduction check；
- spatial normalisation 后 temporal rescaling；
- accepted new mode 后重新更新全部 temporal coefficients；
- overall LATIN convergence 与 PGD saturation 分离。

这些属于 validated implementation experience，不应误写成 Eq. (72) 原论文公式本身。

---

# 45. Eq. (72) 与 Eq. (70)–(71) 的完整 fixed-point 关系

截至本阶段，new-mode alternating enrichment 已经在理论上形成清楚结构：

$$ \boxed{ \lambda^{(k)} \overset{\text{Eq. (70)–(71)}}{\longrightarrow} \bar\varepsilon^{p,(k+1)} \overset{\text{Eq. (72)}}{\longrightarrow} \lambda^{(k+1)} } $$

重复：

$$ \boxed{ \lambda^{(k)}\rightarrow\bar\varepsilon^{p,(k+1)}\rightarrow\lambda^{(k+1)}\rightarrow\bar\varepsilon^{p,(k+2)}\rightarrow\cdots } $$

其共同目标始终是降低：

$$ \boxed{ \left\|\dot\lambda\bar\varepsilon^p-H_\sigma\lambda\bar\sigma+\bar\Delta\right\|_{H_\sigma^{-1}} } $$

---

# 46. Eq. (72) 后仍然还需要的 algorithmic finishing operations

即使 new pair fixed point 收敛，原论文 enrichment 仍未完全结束。

后续还需要：

1. new spatial function 与已有 spatial basis 做 Gram–Schmidt orthonormalisation；
2. spatial normalisation 后相应修正 temporal amplitude，保持 rank-one product 不变；
3. 更新 previously existing temporal functions；
4. 检查 modified new temporal function norm 是否 insignificant；
5. 必要时 rejection 当前 new pair；
6. new pair 接受后才正式进入 global-stage hardening/damage completion。

这些步骤后续应在 Eq. (72) 离散闭合后继续处理。

---

# 47. 当前阶段没有引入的内容

本阶段没有引入：

- cycle–phase–space / $n-\tau-x$ PGD；
- multi-time-scale PGD；
- temporal homogenisation；
- cycle jump；
- stochastic wind-wave loading；
- simultaneous PGD separation of all internal variables；
- damage PGD basis；
- new material model；
- new tower element formulation。

继续保持：

$$ \boxed{ x\rightarrow(e,g,f),\qquad\Delta\varepsilon^p(x,t)=\lambda(t)\bar\varepsilon^p(x) } $$

这一 original $x-t$ architecture。

---

# 48. 本阶段明确不能做的事情

1. 不把 Eq. (72) 当成 structural Newton equation；
2. 不把 $\bar\Delta$ 当成直接 wind load；
3. 不删除 $H_\sigma^{-1}$ weighting；
4. 不把 Eq. (72) 改成普通 Euclidean least squares；
5. 不假设 $\lambda(0)=0$ 而不检查 temporal discretization convention；
6. 不把 continuous strong form 的 natural endpoint condition 未经证明写成 original paper condition；
7. 不声称已经恢复 exact original DG0 algebra；
8. 不声称 backward Euler 与 original DG0 严格等价；
9. 不因为 DG0 细节缺失而重新设计整个 temporal algorithm；
10. 不机械复制当前 1D spatial least-squares half-step 到 tower；
11. 不在 Eq. (72) 尚未离散闭合前开始大范围代码重构；
12. 不提前引入多时间尺度扩展。

---

# 49. 本阶段已经解决的问题

本阶段已经明确：

1. Eq. (72) 的直接来源是 Eq. (61)；
2. Eq. (72) 最小化的是 mechanical LATIN search-direction residual；
3. $H_\sigma^{-1}$ 是与 descent search direction 相匹配的 metric；
4. 为什么 Eq. (71) 后必须重新更新 temporal function；
5. 为什么 Eq. (70)–(71) 与 Eq. (72) 构成 alternating fixed point；
6. 为什么 enrichment 是 spatial Galerkin + temporal minimisation 的 hybrid strategy；
7. Eq. (72) 三个 residual contributions 的准确含义；
8. 为什么真实 tower load 不直接出现在 Eq. (72)；
9. Eq. (72) objective functional 的准确形式；
10. Eq. (72) 的 first variation；
11. continuous temporal weak optimality condition；
12. spatial contraction 后 $A,B,C,d,e$ 五类 temporal coefficients；
13. continuous interior Euler–Lagrange differential equation；
14. 为什么该 continuous equation 一般含 $\ddot\lambda$；
15. 为什么 endpoint conditions 当前不能自行补全；
16. original paper 能确认到 DG0，但正文未给出 exact algebra；
17. generic DG interpretation 与 paper-confirmed information 必须区分；
18. tower material-point Eq. (72) residual；
19. tower quadrature metric 对 Eq. (72) 的作用；
20. tower Eq. (72) 为什么不需要 repeated structural FE solves；
21. 当前 1D `pgd_time_update.py` 的 backward-Euler residual-minimisation structure；
22. 1D forcing 与 original $+\bar\Delta$ 的 sign correspondence；
23. 当前 1D implementation 与 Eq. (72) 在 residual 和 metric 上高度一致；
24. backward Euler 与 original DG0 尚未证明严格等价；
25. 若 original DG0 细节无法恢复，可以采用 validated 1D temporal discretization；
26. 这一选择不改变 Eq. (72) continuous theory；
27. 下一阶段 tower discrete Eq. (72) 应从 backward-Euler single-mode equation 开始；
28. 后续仍需检查 first-time treatment、conditioning、fixed-point 和 basis rescaling。

---

# 50. 当前尚未解决的问题

下一阶段仍需解决：

1. fixed $\vec p$、$\vec s$ 后，single-new-mode backward-Euler tower Eq. (72) 的 exact scalar weighted least-squares equation；
2. 该 scalar equation 的 numerator / denominator 如何由 $(e,g,f)$ quadrature contraction 得到；
3. $\lambda_n$ 与 $\lambda_{n-1}$ 的递推关系；
4. first time point 的具体处理；
5. 当前 1D `update_pgd_time_functions()` 的 first-step treatment 是否直接适用于 tower；
6. denominator / conditioning safeguard；
7. temporal residual reduction check；
8. Eq. (72) 更新后是否立即返回 Eq. (70)–(71) 形成下一 fixed-point iteration；
9. pair convergence 应比较 $\lambda$、$\bar\varepsilon^p$ 还是完整 separated correction；
10. spatial orthonormalisation 与 temporal rescaling 应在 fixed-point convergence 前还是后进行；
11. final new pair 接受后如何重新更新 all temporal functions；
12. original DG0 exact algebra 如果后续找到，如何与当前 backward-Euler route 做严格对比。

---

# 51. 下一阶段的明确工作入口

下一阶段不重新讨论 Eq. (72) 的 continuous origin，而直接继承本阶段结论。

从：

$$ \boxed{ \vec r_n=\left(\frac{\vec p}{\Delta t_n}-D_{H,n}\vec s\right)\lambda_n-\frac{\vec p}{\Delta t_n}\lambda_{n-1}+\vec\Delta_n } $$

开始。

定义：

$$ \boxed{ \vec g_n=\frac{\vec p}{\Delta t_n}-D_{H,n}\vec s } $$

然后在 tower metric：

$$ \boxed{ MD_{H,n}^{-1} } $$

下，对 scalar unknown $\lambda_n$ 做 weighted least-squares minimisation。

下一阶段只需要把该 minimisation 的 first-order condition、scalar denominator、scalar RHS 逐项推出来，并说明每一项的 $(e,g,f)$ 来源。

在完成这一离散推导之前，仍不开始代码修改。

---

# 52. 阶段结论

截至本阶段，原论文 Eq. (72) 的 continuous mathematical structure 已经基本澄清。

核心 residual：

$$ \boxed{ r=\dot\lambda\bar\varepsilon^p-H_\sigma\lambda\bar\sigma+\bar\Delta } $$

核心 objective：

$$ \boxed{ \lambda=\arg\min_\lambda\|r\|_{H_\sigma^{-1}} } $$

其 continuous weak form 已得到：

$$ \boxed{ \int_0^T\left[\dot\eta(A\dot\lambda-B\lambda+d)-\eta(B\dot\lambda-C\lambda+e)\right]dt=0 } $$

其 interior strong form 为：

$$ \boxed{ A\ddot\lambda+\dot A\dot\lambda-C\lambda+\dot d+e=0 } $$

Tower 中保持原 $x-t$ architecture：

$$ \boxed{ x\rightarrow(e,g,f) } $$

并有：

$$ \boxed{ \vec r(t)=\dot\lambda(t)\vec p-D_H(t)\lambda(t)\vec s+\vec\Delta(t) } $$

同时已经明确：原论文确认 temporal problem 使用 DG0，但正文没有给出足以恢复 exact jump/endpoint algebra 的完整离散公式。

因此当前正式研究决策是：

$$ \boxed{ \text{若 exact original DG0 无法恢复，则采用 original Eq. (72) continuous structure + validated 1D backward-Euler temporal discretization} } $$

该决策的目标是：在不改变 Bhattacharyya $x-t$ LATIN-PGD 核心理论结构的前提下，优先继承已经通过 1D three-material bar reproduction 验证的数值经验，避免为了形式上追求 DG0 而重新引入未经验证的 temporal algorithm。

因此，当前真正的下一步已经非常明确：

$$ \boxed{ \text{derive the single-new-mode backward-Euler tower Eq. (72) scalar update} } $$

完成该步后，Eq. (70)–(72) 的 tower new-mode alternating enrichment 才能在离散层面真正闭合。
