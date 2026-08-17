# Tower LATIN-PGD Eq. (70)–(71) 空间有限元求解与塑性空间模态恢复阶段总结

**日期：2026-08-17**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**当前研究路线：原论文 $x-t$ LATIN-PGD → fiber beam-column offshore wind turbine tower**  
**阶段范围：原论文 Eq. (70)–(71)**  
**上一阶段衔接：`2026-08-17-tower-latin-pgd-eq65-69-spatial-admissibility-stage-summary.md`**  
**下一阶段：Eq. (72) fixed-spatial-mode temporal update**

---

# 1. 本阶段定位

上一阶段 Eq. (65)–(69) 已经完成了 new-mode enrichment spatial subproblem 的 global admissibility 构造。

已得到：

$$ \boxed{ \bar{\tilde{\varepsilon}} = W^{-1}\bar{\sigma} - \bar{\delta} } \tag{67} $$

以及：

$$ \boxed{ \bar{\tilde{\varepsilon}} \in E_0 } $$

因此存在：

$$ \boxed{ \bar{\tilde{\varepsilon}} = \varepsilon(\bar{\tilde{u}}), \qquad \bar{\tilde{u}}\in U_0 } $$

同时：

$$ \boxed{ \bar{\sigma}\in S_0 } $$

即：

$$ \boxed{ \int_\Omega \bar{\sigma}:\varepsilon(\bar{u}^*)\,d\Omega = 0, \qquad \forall\bar{u}^*\in U_0 } \tag{69} $$

本阶段 Eq. (70)–(71) 的任务是：

1. 将 Eq. (67)–(69) 转化为一个 standard displacement-like spatial FE problem；
2. 求得 auxiliary displacement-like spatial field $\bar{\tilde u}$；
3. 恢复 compatible auxiliary strain field $\bar{\tilde\varepsilon}$；
4. 恢复 equilibrated stress spatial mode $\bar\sigma$；
5. 最终恢复真正要加入 PGD basis 的 plastic-strain spatial mode $\bar\varepsilon^p_{m+1}$。

因此，本阶段完成的是：

$$ \boxed{ \text{global-admissible spatial problem} \rightarrow \text{FE solve} \rightarrow \text{new plastic PGD spatial mode} } $$

---

# 2. 上一阶段已经确定的 $W^{-1}$ 与 $\bar\delta$

Eq. (67) 中：

$$ \boxed{ W^{-1} = \left\langle H_\sigma\lambda^2 \right\rangle + \left\langle \dot{\lambda}\lambda \right\rangle \mathbb C^{-1} } $$

以及：

$$ \boxed{ \bar{\delta} = \left\langle \bar{\Delta}\lambda \right\rangle } $$

其中：

- $H_\sigma$：LATIN descent search-direction operator；
- $\lambda(t)$：当前 fixed-point spatial half-step 中暂时固定的 temporal function；
- $\mathbb C$：reference elastic operator；
- $\bar\Delta(x,t)$：remaining descent-search-direction defect；
- $\bar\delta(x)$：$\bar\Delta$ 沿当前 temporal function $\lambda(t)$ 的 temporal projection。

必须继续保持：

$$ \boxed{ W^{-1}\neq\mathbb C^{-1} } $$

同时：

$$ \boxed{ W \neq \text{tower global stiffness matrix} } $$

$W$ 首先是当前 enrichment spatial solve 中的 effective stiffness-like material-point/continuum spatial operator。

---

# 3. Eq. (70) 的直接代数起点

Eq. (67)：

$$ \bar{\tilde{\varepsilon}} = W^{-1}\bar{\sigma} - \bar{\delta} $$

移项：

$$ W^{-1}\bar{\sigma} = \bar{\tilde{\varepsilon}} + \bar{\delta} $$

因此：

$$ \boxed{ \bar{\sigma} = W\left(\bar{\tilde{\varepsilon}}+\bar{\delta}\right) } $$

而 Eq. (68) 已经证明：

$$ \boxed{ \bar{\tilde{\varepsilon}} = \varepsilon(\bar{\tilde u}) } $$

所以：

$$ \boxed{ \bar{\sigma} = W\left[\varepsilon(\bar{\tilde u})+\bar{\delta}\right] } $$

---

# 4. 将 Eq. (67)–(68) 代入 Eq. (69)

Eq. (69)：

$$ \int_\Omega \bar{\sigma}:\varepsilon(\bar u^*)\,d\Omega = 0 $$

代入：

$$ \bar{\sigma} = W\left[\varepsilon(\bar{\tilde u})+\bar{\delta}\right] $$

得到：

$$ \int_\Omega W\left[\varepsilon(\bar{\tilde u})+\bar{\delta}\right]:\varepsilon(\bar u^*)\,d\Omega = 0 $$

展开：

$$ \int_\Omega W\varepsilon(\bar{\tilde u}):\varepsilon(\bar u^*)\,d\Omega + \int_\Omega W\bar{\delta}:\varepsilon(\bar u^*)\,d\Omega = 0 $$

将 known source term 移至右端：

$$ \boxed{ \int_\Omega W\varepsilon(\bar{\tilde u}):\varepsilon(\bar u^*)\,d\Omega = -\int_\Omega W\bar{\delta}:\varepsilon(\bar u^*)\,d\Omega, \qquad \forall\bar u^*\in U_0 } \tag{70} $$

这就是 Eq. (70) 的核心结构。

---

# 5. Eq. (70) 没有引入新的 constitutive assumption

Eq. (70) 完全来自：

$$ \boxed{ \text{Eq. (67)} + \text{Eq. (68)} + \text{Eq. (69)} } $$

因此它不是新的材料本构方程，也不是新的物理平衡定律。

其作用是：

$$ \boxed{ \text{stress-space constrained problem} \rightarrow \text{displacement-space FE problem} } $$

---

# 6. 为什么 Eq. (70) 是“经典有限元形式”

标准线弹性 weak form 具有：

$$ \int_\Omega \mathbb C\varepsilon(u):\varepsilon(u^*)\,d\Omega = \text{external virtual work} $$

Eq. (70) 左端：

$$ \int_\Omega W\varepsilon(\bar{\tilde u}):\varepsilon(\bar u^*)\,d\Omega $$

与标准有限元形式完全同构，只是：

$$ \boxed{ \mathbb C \rightarrow W } $$

因此 $W$ 在 Eq. (70) 中承担 effective stiffness-like operator 的角色。

---

# 7. $W$ 的准确身份

$W$ 不是：

- reference elasticity $\mathbb C$；
- damaged tangent stiffness；
- algorithmic material tangent；
- tower global stiffness；
- accepted PGD basis property。

$W$ 是：

$$ \boxed{ \text{当前 fixed temporal function 下的 enrichment spatial effective stiffness-like operator} } $$

它由：

$$ W^{-1} = \left\langle H_\sigma\lambda^2 \right\rangle + \left\langle \dot{\lambda}\lambda \right\rangle\mathbb C^{-1} $$

定义。

因此 $W$ 随：

- 当前 temporal function $\lambda(t)$；
- 当前 search-direction field $H_\sigma(x,t)$；
- 当前 LATIN state；

变化。

---

# 8. Eq. (70) 右端的准确含义

右端为：

$$ \boxed{ -\int_\Omega W\bar{\delta}:\varepsilon(\bar u^*)\,d\Omega } $$

其中：

$$ \bar{\delta} = \left\langle\bar{\Delta}\lambda\right\rangle $$

所以右端来源于：

$$ \boxed{ \text{remaining LATIN search-direction defect 的 temporal projection} } $$

经 $W$ 转换后形成 equivalent spatial forcing。

因此它不是：

- wind load；
- wave load；
- gravity load；
- tower top load；
- Newton force residual。

更准确地称为：

$$ \boxed{ \text{projected-defect-induced equivalent spatial forcing} } $$

---

# 9. Eq. (70) 右端负号的来源

负号完全来自 Eq. (67) 的代数移项。

由：

$$ \bar{\sigma} = W(\bar{\tilde{\varepsilon}}+\bar{\delta}) $$

代入 equilibrium：

$$ \int_\Omega W\bar{\tilde{\varepsilon}}:\varepsilon(\bar u^*)\,d\Omega + \int_\Omega W\bar{\delta}:\varepsilon(\bar u^*)\,d\Omega = 0 $$

所以：

$$ \boxed{ \int_\Omega W\bar{\tilde{\varepsilon}}:\varepsilon(\bar u^*)\,d\Omega = -\int_\Omega W\bar{\delta}:\varepsilon(\bar u^*)\,d\Omega } $$

未来代码实现时不能凭经验把右端改成正号。

---

# 10. $\bar{\tilde u}$ 的准确含义

Eq. (70) 中求解：

$$ \boxed{ \bar{\tilde u} } $$

但它不是：

- 实际 tower displacement history；
- LATIN total displacement；
- physical vibration mode；
- final PGD plastic spatial mode；
- direct deformation snapshot。

它是：

$$ \boxed{ \text{为构造 new enrichment spatial mode 而引入的 auxiliary displacement-like spatial field} } $$

其作用是自动保证：

$$ \bar{\tilde{\varepsilon}} = \varepsilon(\bar{\tilde u}) $$

满足 kinematic admissibility。

---

# 11. 为什么不直接用 $\bar\sigma$ 做 unknown

若直接求 stress spatial mode $\bar\sigma$，必须显式满足：

$$ \bar\sigma\in S_0 $$

通常意味着需要：

- stress-space basis；
- equilibrium constraints；
- Lagrange multiplier；
- mixed formulation。

通过 Eq. (68) 引入：

$$ \bar{\tilde{\varepsilon}} = \varepsilon(\bar{\tilde u}) $$

则可将 constrained stress problem 转换成 conventional displacement-based FE problem。

因此 Eq. (70) 的主要数值价值是：

$$ \boxed{ \text{constrained stress solve} \rightarrow \text{standard displacement solve} } $$

---

# 12. Tower 的运动学映射

前面已经建立：

$$ \boxed{ \vec{\varepsilon} = H\vec U } $$

其中：

- $\vec U$：tower free structural DOF vector；
- $H$：nodal DOF → fiber material-point strain mapping；
- $\vec\varepsilon$：flatten 后 fiber strain vector。

对于当前 coarse tower：

$$ N_e=10,\qquad N_g=2,\qquad N_f=16 $$

因此：

$$ \boxed{ N_q=N_eN_gN_f=320 } $$

所以：

$$ \vec\varepsilon\in\mathbb R^{320} $$

---

# 13. Tower quadrature metric

定义：

$$ \boxed{ M=\operatorname{diag}(v_q) } $$

其中：

$$ \boxed{ v_q=A_{egf}w_gJ_e } $$

连续空间积分：

$$ \int_\Omega a(x)b(x)\,d\Omega $$

离散为：

$$ \vec a^{\,T}M\vec b $$

因此 Eq. (70) 可以直接迁移到 fiber quadrature space。

---

# 14. Tower 中 $W$ 的离散形式

当前每个 fiber material point 为 scalar axial material。

因此：

$$ C_0=E_0 $$

定义：

$$ A_q=\left\langle H_{\sigma,q}\lambda^2\right\rangle $$

以及：

$$ a=\left\langle\dot{\lambda}\lambda\right\rangle $$

则：

$$ \boxed{ W_q^{-1}=A_q+\frac{a}{E_0} } $$

若可逆：

$$ \boxed{ W_q=\frac{1}{A_q+a/E_0} } $$

定义：

$$ \boxed{ D_W=\operatorname{diag}(W_1,\ldots,W_{N_q}) } $$

---

# 15. Tower Eq. (70) 左端离散

有：

$$ \vec{\bar{\tilde{\varepsilon}}}=H\vec{\bar{\tilde U}} $$

任意 virtual field：

$$ \vec{\varepsilon}^*=H\vec U^* $$

Eq. (70) 左端离散为：

$$ (H\vec U^*)^TMD_W(H\vec{\bar{\tilde U}}) $$

整理：

$$ (\vec U^*)^TH^TMD_WH\vec{\bar{\tilde U}} $$

因此定义：

$$ \boxed{ K_W=H^TMD_WH } $$

这是当前 tower enrichment spatial half-step 的 effective structural stiffness matrix。

---

# 16. Tower Eq. (70) 右端离散

连续右端：

$$ -\int_\Omega W\bar\delta:\varepsilon(\bar u^*)\,d\Omega $$

离散后：

$$ -(H\vec U^*)^TMD_W\vec{\bar\delta} $$

整理：

$$ -(\vec U^*)^TH^TMD_W\vec{\bar\delta} $$

因此定义：

$$ \boxed{ \vec f_\delta=-H^TMD_W\vec{\bar\delta} } $$

所以 tower Eq. (70) matrix form 为：

$$ \boxed{ H^TMD_WH\vec{\bar{\tilde U}}=-H^TMD_W\vec{\bar\delta} } $$

即：

$$ \boxed{ K_W\vec{\bar{\tilde U}}=\vec f_\delta } $$

---

# 17. Eq. (70) 现在正式确认了此前的 $K_W$ 猜测

上一阶段只根据 Eq. (67)–(69) 预判可能出现：

$$ K_W=H^TMD_WH $$

本阶段通过 Eq. (70) 的严格 weak-form 离散正式确认：

$$ \boxed{ K_W=H^TMD_WH } $$

因此，这已经不是 heuristic guess，而是 tower Eq. (70) 的直接离散结果。

---

# 18. 与旧 $C_0$-based reference stiffness 的关系

此前 tower reference operator：

$$ \boxed{ K^0=H^TMC_0H } $$

当前 enrichment Eq. (70)：

$$ \boxed{ K_W=H^TMD_WH } $$

两者具有完全相同的 assembly skeleton：

$$ \boxed{ H^TM(\cdot)H } $$

区别只在 material-point operator。

reference problem：

$$ \boxed{ C_0 } $$

enrichment spatial problem：

$$ \boxed{ D_W } $$

因此形成明确代码设计判断：

> 旧 tower equilibrium/global assembly 框架可以高度复用，但不能在 Eq. (70) 中继续把 material-point operator 固定为 $C_0$。

---

# 19. 代码层未来应抽象的 operator 结构

概念上可将现有 fixed-reference assembly：

```text
K0 = H.T @ M @ C0 @ H
```

推广为：

```text
K(A) = H.T @ M @ A @ H
```

其中：

```text
A = C0
```

用于 reference equilibrium operator；

而：

```text
A = D_W
```

用于 Eq. (70) new-mode enrichment spatial solve。

当前阶段只形成设计结论，不立即修改代码。

---

# 20. 实际实现不应显式建立 dense $D_W$

数学上：

$$ D_W=\operatorname{diag}(W_q) $$

但实际代码不应建立 dense $N_q\times N_q$ diagonal matrix。

只需存：

$$ \boxed{ \vec W=[W_1,\ldots,W_{N_q}]^T } $$

然后对 $H$ 每一行乘：

$$ v_qW_q $$

从而：

$$ K_W=H^T[(v_qW_q)H_q] $$

同理：

$$ \vec f_\delta=-H^T[v_qW_q\bar\delta_q] $$

这样既符合公式，也更高效。

---

# 21. $W_q$ 在 fixed-point 中是否固定

在一次 spatial half-step 中：

$$ \lambda^{(k)}(t) $$

暂时固定，因此当前：

$$ W_q^{(k)} $$

也是固定的。

但是更新 temporal function 后：

$$ \lambda^{(k)}\rightarrow\lambda^{(k+1)} $$

一般：

$$ W_q^{(k+1)}\neq W_q^{(k)} $$

因此：

$$ \boxed{ K_W^{(k+1)} \text{ 一般需要重新组装} } $$

这与可以长期缓存的 reference stiffness $K^0$ 不同。

---

# 22. Eq. (70) 的计算量意义

虽然每次 enrichment spatial half-step 需要重新构造 $W_q$ 和 $K_W$，但它仍然只是：

$$ \boxed{ \text{one pure spatial solve} } $$

而不是：

$$ \boxed{ N_t \text{ 个 full FE solves} } $$

整个时间历史先被 contraction 为：

$$ W_q,\qquad\bar\delta_q $$

然后仅求：

$$ K_W\bar{\tilde U}=f_\delta $$

这保持了 original $x-t$ PGD 的核心降维逻辑。

---

# 23. $W_q^{-1}$ 的可逆性要求

因为 Eq. (67) 定义 $W^{-1}$，Eq. (70) 又需要 $W$，因此至少要求：

$$ \boxed{ W_q^{-1}\neq0 } $$

对 scalar fiber case。

更强的数值要求是：

$$ \boxed{ W_q>0 } $$

从而保证 $K_W$ 具有稳定的 stiffness-like 性质。

---

# 24. $H_\sigma$ 正定性与 Eq. (70) 可解性的关系

有：

$$ W_q^{-1}=\left\langle H_{\sigma,q}\lambda^2\right\rangle+\frac{a}{E_0} $$

其中：

$$ a=\left\langle\dot\lambda\lambda\right\rangle $$

如果：

$$ H_{\sigma,q}>0 $$

则：

$$ \left\langle H_{\sigma,q}\lambda^2\right\rangle>0 $$

通常为 $W_q^{-1}$ 提供主要正贡献。

因此形成数值链：

$$ \boxed{ H_\sigma\text{ positive} \rightarrow W\text{ stable} \rightarrow K_W\text{ solvable} } $$

这进一步解释了之前 1D reproduction 中为什么必须对 search-direction operator 的正值性进行保护。

---

# 25. $\langle\lambda\dot\lambda\rangle$ 在 $W^{-1}$ 中的作用

定义：

$$ a=\left\langle\lambda\dot\lambda\right\rangle $$

且：

$$ \boxed{ a=\frac12[\lambda^2(T)-\lambda^2(0)] } $$

因此 $a$ 可以为：

- 正；
- 零；
- 负。

所以 $W^{-1}$ 的正值性不能只依靠 $a/E_0$。

主要仍依赖：

$$ \left\langle H_\sigma\lambda^2\right\rangle $$

的正贡献。

---

# 26. 若 $D_W>0$，tower $K_W$ 的期望性质

若：

$$ M>0 $$

且：

$$ D_W>0 $$

并且 fixed-base boundary conditions 已移除 rigid-body modes，则：

$$ K_W=H^TMD_WH $$

在 free-DOF space 上应具有 standard symmetric positive-definite-like structure。

因此可以使用 conventional linear solver。

这一点属于 tower discrete-form 的数学推论，后续应通过 unit tests 验证，而不是只凭理论假定。

---

# 27. Eq. (70) solve 后如何恢复 spatial fields

求得：

$$ \boxed{ \vec{\bar{\tilde U}} } $$

后首先恢复：

$$ \boxed{ \vec{\bar{\tilde{\varepsilon}}}=H\vec{\bar{\tilde U}} } $$

然后：

$$ \boxed{ \vec{\bar{\sigma}}=D_W\left(\vec{\bar{\tilde{\varepsilon}}}+\vec{\bar\delta}\right) } $$

即：

$$ \boxed{ \vec{\bar{\sigma}}=D_W\left(H\vec{\bar{\tilde U}}+\vec{\bar\delta}\right) } $$

---

# 28. Eq. (70) 为什么自动保证 static equilibrium

由：

$$ K_W\vec{\bar{\tilde U}}=-H^TMD_W\vec{\bar\delta} $$

即：

$$ H^TMD_WH\vec{\bar{\tilde U}}+H^TMD_W\vec{\bar\delta}=0 $$

提取：

$$ H^TM D_W\left(H\vec{\bar{\tilde U}}+\vec{\bar\delta}\right)=0 $$

而：

$$ D_W\left(H\vec{\bar{\tilde U}}+\vec{\bar\delta}\right)=\vec{\bar\sigma} $$

因此：

$$ \boxed{ H^TM\vec{\bar\sigma}=0 } $$

所以 Eq. (70) solve 后恢复出的 stress spatial mode 自动满足 Eq. (69)。

---

# 29. Eq. (70) 为什么自动保证 compatibility

因为：

$$ \vec{\bar{\tilde{\varepsilon}}}=H\vec{\bar{\tilde U}} $$

所以：

$$ \boxed{ \vec{\bar{\tilde{\varepsilon}}}\in\operatorname{Range}(H) } $$

即自动满足 Eq. (68)。

因此 Eq. (70) 一次 displacement-like solve 同时关闭：

$$ \boxed{ \text{kinematic compatibility} + \text{static equilibrium} } $$

---

# 30. Eq. (70) 与 Newton solve 的区别

普通 Newton：

$$ K_{\mathrm{tan}}\Delta U=F_{\mathrm{ext}}-F_{\mathrm{int}} $$

Eq. (70)：

$$ K_W\bar{\tilde U}=-H^TMD_W\bar\delta $$

二者形式相似，但 RHS 完全不同。

Eq. (70) RHS 来自：

$$ \bar\delta=\left\langle\bar\Delta\lambda\right\rangle $$

因此：

$$ \boxed{ f_\delta\neq F_{\mathrm{ext}}-F_{\mathrm{int}} } $$

代码命名时不能使用：

```text
newton_residual
```

更合理的概念名为：

```text
projected_defect_load
```

或：

```text
enrichment_spatial_rhs
```

---

# 31. Eq. (70) 的“source-strain-like”解释

由：

$$ \bar{\sigma}=W(\bar{\tilde\varepsilon}+\bar\delta) $$

形式上类似：

$$ \sigma=C(\varepsilon-\varepsilon^*) $$

如果定义：

$$ \varepsilon^*=-\bar\delta $$

则：

$$ \bar{\sigma}=W(\bar{\tilde\varepsilon}-\varepsilon^*) $$

因此 Eq. (70) 可类比为：

$$ \boxed{ \text{effective stiffness} + \text{known source/eigenstrain-like field} } $$

但为了避免物理误解，当前仍建议使用：

$$ \boxed{ \text{projected-defect source field} } $$

而不是直接称 $\bar\delta$ 为真实 eigenstrain。

---

# 32. Eq. (70) 尚未得到最终 PGD plastic spatial mode

Eq. (70) 得到的是：

$$ \bar{\tilde u} $$

继而：

$$ \bar{\tilde\varepsilon},\qquad\bar\sigma $$

但最终需要加入 PGD basis 的量仍是：

$$ \boxed{ \bar\varepsilon^p_{m+1} } $$

因此 Eq. (70) 只是：

$$ \boxed{ \text{spatial FE solve} } $$

而不是：

$$ \boxed{ \text{final enrichment mode completion} } $$

最终 plastic-strain spatial mode 由 Eq. (71) 恢复。

---

# 33. Eq. (71) 的原始目标

Eq. (71) 的作用是：

$$ \boxed{ \bar{\tilde\varepsilon},\bar\sigma,\lambda \rightarrow \bar\varepsilon^p_{m+1} } $$

为简化记号，定义：

$$ \boxed{ a=\left\langle\dot\lambda\lambda\right\rangle } $$

Eq. (71) 可写成：

$$ \boxed{ \bar{\varepsilon}^{p}=\frac{1}{a}\left[\bar{\tilde{\varepsilon}}-a\mathbb C^{-1}W(\bar{\tilde{\varepsilon}}+\bar\delta)\right] } \tag{71} $$

利用：

$$ \bar\sigma=W(\bar{\tilde\varepsilon}+\bar\delta) $$

可改写为：

$$ \boxed{ \bar{\varepsilon}^{p}=\frac{1}{a}\bar{\tilde{\varepsilon}}-\mathbb C^{-1}\bar{\sigma} } \tag{71'} $$

Eq. (71') 是理解 Eq. (71) 最清楚的形式。

---

# 34. Eq. (71) 的推导必须回到 Eq. (61)

Eq. (61)：

$$ \boxed{ \Delta\dot\varepsilon^p-H_\sigma\Delta\sigma'+\bar\Delta=0 } $$

new pair：

$$ \Delta\dot\varepsilon^p=\dot\lambda\bar\varepsilon^p $$

$$ \Delta\sigma'=\lambda\bar\sigma $$

代入：

$$ \dot\lambda\bar\varepsilon^p-H_\sigma\lambda\bar\sigma+\bar\Delta=0 $$

两边乘 temporal test function：

$$ \lambda(t) $$

得到：

$$ \lambda\dot\lambda\bar\varepsilon^p-H_\sigma\lambda^2\bar\sigma+\bar\Delta\lambda=0 $$

---

# 35. 对时间积分得到 plastic spatial relation

因为 $\bar\varepsilon^p$ 与 $\bar\sigma$ 不依赖 $t$，所以：

$$ \left\langle\lambda\dot\lambda\right\rangle\bar\varepsilon^p-\left\langle H_\sigma\lambda^2\right\rangle\bar\sigma+\left\langle\bar\Delta\lambda\right\rangle=0 $$

即：

$$ \boxed{ a\bar\varepsilon^p=\left\langle H_\sigma\lambda^2\right\rangle\bar\sigma-\bar\delta } \tag{A} $$

因此 Eq. (71) 的 plastic spatial mode 本质上来自：

$$ \boxed{ \text{Eq. (61) 沿 }\lambda(t)\text{ 的 temporal Galerkin projection} } $$

---

# 36. Eq. (67) 与 Eq. (A) 的结合

Eq. (67) 展开为：

$$ \bar{\tilde\varepsilon}=\left[\left\langle H_\sigma\lambda^2\right\rangle+a\mathbb C^{-1}\right]\bar\sigma-\bar\delta $$

重新分组：

$$ \bar{\tilde\varepsilon}=\left[\left\langle H_\sigma\lambda^2\right\rangle\bar\sigma-\bar\delta\right]+a\mathbb C^{-1}\bar\sigma $$

根据 Eq. (A)：

$$ \left\langle H_\sigma\lambda^2\right\rangle\bar\sigma-\bar\delta=a\bar\varepsilon^p $$

因此：

$$ \boxed{ \bar{\tilde\varepsilon}=a\bar\varepsilon^p+a\mathbb C^{-1}\bar\sigma } $$

即：

$$ \boxed{ \bar{\tilde\varepsilon}=a\left(\bar\varepsilon^p+\mathbb C^{-1}\bar\sigma\right) } \tag{B} $$

---

# 37. Eq. (71) 中 $1/a$ 的来源

由 Eq. (B)：

$$ \frac{1}{a}\bar{\tilde\varepsilon}=\bar\varepsilon^p+\mathbb C^{-1}\bar\sigma $$

所以：

$$ \boxed{ \bar\varepsilon^p=\frac{1}{a}\bar{\tilde\varepsilon}-\mathbb C^{-1}\bar\sigma } $$

因此：

$$ \boxed{ \frac{1}{\left\langle\dot\lambda\lambda\right\rangle} } $$

不是人为 normalization。

它来自 plastic-rate term：

$$ \dot\lambda\bar\varepsilon^p $$

与 test temporal function：

$$ \lambda $$

相乘并做时间积分后得到的：

$$ \left\langle\lambda\dot\lambda\right\rangle\bar\varepsilon^p $$

因此分母是 temporal Galerkin projection 的自然结果。

---

# 38. Eq. (71) 的物理内核其实是 strain partition

定义：

$$ \boxed{ \bar\varepsilon=\frac{1}{a}\bar{\tilde\varepsilon} } $$

则 Eq. (71')：

$$ \bar\varepsilon^p=\bar\varepsilon-\mathbb C^{-1}\bar\sigma $$

而：

$$ \mathbb C^{-1}\bar\sigma=\bar\varepsilon^e $$

所以：

$$ \boxed{ \bar\varepsilon=\bar\varepsilon^e+\bar\varepsilon^p } $$

因此 Eq. (71) 看似复杂，其核心仍然是经典关系：

$$ \boxed{ \text{plastic strain}=\text{total strain}-\text{elastic strain} } $$

---

# 39. $\bar{\tilde\varepsilon}$ 在 Eq. (71) 后的更精确解释

此前 Eq. (67) 时只称：

$$ \bar{\tilde\varepsilon} $$

为 auxiliary strain-like field。

Eq. (68) 后可确认其 compatible。

Eq. (71) 后进一步得到：

$$ \boxed{ \bar{\tilde\varepsilon}=a\bar\varepsilon } $$

其中 $\bar\varepsilon$ 是与 final new PGD pair 对应的 compatible total-strain spatial mode。

因此更精确地说：

$$ \boxed{ \bar{\tilde\varepsilon}=\text{temporally scaled compatible total-strain spatial field} } $$

---

# 40. Eq. (71) 后重新回到 canonical $C$-based relation

由：

$$ \bar\varepsilon=\frac{1}{a}\bar{\tilde\varepsilon} $$

以及：

$$ \bar\varepsilon^p=\bar\varepsilon-\mathbb C^{-1}\bar\sigma $$

可得：

$$ \boxed{ \bar\sigma=\mathbb C(\bar\varepsilon-\bar\varepsilon^p) } $$

这说明：

> $W$ 只是 spatial enrichment solve 中的临时 effective operator；new mode 生成完成后，其 canonical stress/strain/plastic-strain relation 重新回到 reference elasticity $\mathbb C$。

---

# 41. 这与早先 Eq. (47)–(53) 的 tower spatial relation 一致

此前 tower mode 已定义：

$$ \bar\varepsilon^p \rightarrow \bar\varepsilon=\mathcal E_{\mathrm{tower}}\bar\varepsilon^p $$

以及：

$$ \bar\sigma=C_0(\bar\varepsilon-\bar\varepsilon^p) $$

Eq. (71) 说明 new-mode enrichment 最终恢复出的：

$$ \bar\varepsilon^p $$

与：

$$ \bar\sigma $$

仍然满足同一 canonical $C_0$-based relation。

因此 Eq. (70) 中 $W$-based solve 并没有改变最终 PGD basis 的物理定义。

---

# 42. $W$ 不需要成为 accepted PGD basis 的永久数据

accepted PGD mode 真正需要保存的是：

$$ \boxed{ \bar\varepsilon_j^p,\quad\bar\sigma_j,\quad\lambda_j(t),\quad\dot\lambda_j(t) } $$

而：

$$ W_q $$

只属于当前 enrichment fixed-point spatial half-step。

当：

$$ \lambda^{(k)}\rightarrow\lambda^{(k+1)} $$

时：

$$ W_q^{(k)}\rightarrow W_q^{(k+1)} $$

因此：

$$ \boxed{ W \text{ 是临时算法量，不是 permanent basis property} } $$

---

# 43. Tower Eq. (71) 的离散形式

定义：

$$ a=\left\langle\dot\lambda\lambda\right\rangle $$

Eq. (70) 已得到：

$$ \vec{\bar{\tilde\varepsilon}}=H\vec{\bar{\tilde U}} $$

以及：

$$ \vec{\bar\sigma}=D_W\left(H\vec{\bar{\tilde U}}+\vec{\bar\delta}\right) $$

所以：

$$ \boxed{ \vec{\bar\varepsilon}^{\,p}=\frac{1}{a}\vec{\bar{\tilde\varepsilon}}-C_0^{-1}\vec{\bar\sigma} } $$

即：

$$ \boxed{ \vec{\bar\varepsilon}^{\,p}=\frac{1}{a}H\vec{\bar{\tilde U}}-C_0^{-1}D_W\left(H\vec{\bar{\tilde U}}+\vec{\bar\delta}\right) } $$

对统一 scalar steel modulus：

$$ C_0=E_0I $$

因此：

$$ \boxed{ \vec{\bar\varepsilon}^{\,p}=\frac{1}{a}H\vec{\bar{\tilde U}}-\frac{1}{E_0}D_W\left(H\vec{\bar{\tilde U}}+\vec{\bar\delta}\right) } $$

---

# 44. Eq. (71) 最终返回 fiber material-point space

当前 coarse tower：

$$ N_q=320 $$

所以：

$$ \boxed{ \vec{\bar\varepsilon}^{\,p}_{m+1}\in\mathbb R^{320} } $$

Eq. (70) 在 nodal structural DOF space 求：

$$ \vec{\bar{\tilde U}} $$

Eq. (71) 再返回 fiber material-point plastic-strain space。

因此完整映射为：

$$ \boxed{ \text{fiber space-time defect} \rightarrow \text{nodal auxiliary solve} \rightarrow \text{fiber plastic spatial mode} } $$

这与项目当前坚持的 PGD 主 spatial variable 位于 fiber plastic-strain space 完全一致。

---

# 45. Tower Eq. (70)–(71) 的完整 spatial half-step

固定当前：

$$ \lambda(t) $$

第一步，计算：

$$ \boxed{ a=\left\langle\dot\lambda\lambda\right\rangle } $$

第二步，每个 fiber 计算：

$$ \boxed{ W_q^{-1}=\left\langle H_{\sigma,q}\lambda^2\right\rangle+\frac{a}{E_0} } $$

第三步：

$$ \boxed{ \bar\delta_q=\left\langle\bar\Delta_q\lambda\right\rangle } $$

第四步，组装：

$$ \boxed{ K_W=H^TMD_WH } $$

第五步，求：

$$ \boxed{ K_W\vec{\bar{\tilde U}}=-H^TMD_W\vec{\bar\delta} } $$

第六步：

$$ \boxed{ \vec{\bar{\tilde\varepsilon}}=H\vec{\bar{\tilde U}} } $$

第七步：

$$ \boxed{ \vec{\bar\sigma}=D_W\left(\vec{\bar{\tilde\varepsilon}}+\vec{\bar\delta}\right) } $$

第八步，Eq. (71)：

$$ \boxed{ \vec{\bar\varepsilon}^{\,p}=\frac{1}{a}\vec{\bar{\tilde\varepsilon}}-C_0^{-1}\vec{\bar\sigma} } $$

至此得到真正的 new plastic-strain spatial mode。

---

# 46. Eq. (70)–(71) 的三类 consistency checks

未来实现后，至少应验证以下三个关系。

第一，reference stress relation：

$$ \boxed{ \vec{\bar\sigma}\approx C_0\left(\vec{\bar\varepsilon}-\vec{\bar\varepsilon}^{\,p}\right) } $$

其中：

$$ \vec{\bar\varepsilon}=\frac{1}{a}\vec{\bar{\tilde\varepsilon}} $$

第二，static equilibrium：

$$ \boxed{ H^TM\vec{\bar\sigma}\approx0 } $$

第三，kinematic compatibility：

$$ \boxed{ \vec{\bar{\tilde\varepsilon}}=H\vec{\bar{\tilde U}} } $$

或者等价检查：

$$ \vec{\bar\varepsilon}\in\operatorname{Range}(H) $$

这三项同时通过，才能较强地证明 Eq. (70)–(71) tower implementation 正确。

---

# 47. Eq. (71) 的关键数值风险：$a\approx0$

因为：

$$ \bar\varepsilon^p\propto\frac{1}{a} $$

所以需要：

$$ \boxed{ a=\left\langle\dot\lambda\lambda\right\rangle\neq0 } $$

又因为：

$$ a=\frac12[\lambda^2(T)-\lambda^2(0)] $$

如果：

$$ |\lambda(T)|\approx|\lambda(0)| $$

则可能：

$$ |a|\ll1 $$

导致 Eq. (71) 数值放大甚至病态。

---

# 48. 循环荷载并不自动意味着 $a=0$

必须区分：

$$ \text{external load history} $$

与：

$$ \text{PGD temporal amplitude }\lambda(t) $$

即使：

$$ F(T)=F(0) $$

也不意味着：

$$ \lambda(T)=\lambda(0) $$

尤其 plastic correction 具有历史性和不可逆性。

因此不能因为“完成完整循环”就预先删除：

$$ a=\left\langle\dot\lambda\lambda\right\rangle $$

或认为它必然为零。

---

# 49. Eq. (71) 需要 denominator safeguard

未来代码实现至少需要检查：

$$ \boxed{ |a|>a_{\min} } $$

若：

$$ |a|\ll1 $$

不应简单使用：

```text
a = max(a, a_min)
```

因为这会改变原论文 Eq. (71) 的数学关系。

更合理的候选策略是：

- reject current temporal candidate；
- reinitialise temporal mode；
- rescale candidate pair；
- restart fixed-point enrichment；
- 采用与原论文/1D solver consistent 的 failure handling。

具体 safeguard 需要在 Eq. (72) 和 fixed-point stopping strategy 推导后再确定。

---

# 50. PGD pair 的 scaling indeterminacy

rank-one pair：

$$ \lambda(t)\bar\varepsilon^p(x) $$

具有 scaling freedom：

$$ \lambda^*=c\lambda $$

$$ \bar\varepsilon^{p*}=\frac{1}{c}\bar\varepsilon^p $$

因此：

$$ \boxed{ \lambda^*\bar\varepsilon^{p*}=\lambda\bar\varepsilon^p } $$

单独的 temporal/spatial amplitudes 并不唯一。

真正有物理/算法意义的是其乘积。

---

# 51. Eq. (71) 与 scaling indeterminacy 相容

若：

$$ \lambda^*=c\lambda $$

则：

$$ a^*=c^2a $$

同时：

$$ W^{*-1}=c^2W^{-1} $$

所以：

$$ W^*=\frac{1}{c^2}W $$

且：

$$ \bar\delta^*=c\bar\delta $$

相应可得到：

$$ \bar{\tilde\varepsilon}^*=c\bar{\tilde\varepsilon} $$

$$ \bar\sigma^*=\frac{1}{c}\bar\sigma $$

最终：

$$ \bar\varepsilon^{p*}=\frac{1}{c}\bar\varepsilon^p $$

所以：

$$ \boxed{ \lambda^*\bar\varepsilon^{p*}=\lambda\bar\varepsilon^p } $$

这说明 Eq. (70)–(71) 与 PGD scaling freedom 自洽。

---

# 52. 对 tower spatial orthonormalisation 的意义

Eq. (71) 得到 new raw spatial mode：

$$ \vec{\bar\varepsilon}_{m+1}^p $$

之后，仍需与已有：

$$ \vec{\bar\varepsilon}_1^p,\ldots,\vec{\bar\varepsilon}_m^p $$

做 weighted orthogonalisation。

tower spatial inner product 应基于：

$$ M=\operatorname{diag}(v_q) $$

若 normalization：

$$ \bar\varepsilon^p\rightarrow\frac{1}{s}\bar\varepsilon^p $$

则 temporal function 必须：

$$ \lambda\rightarrow s\lambda $$

以保持：

$$ \lambda\bar\varepsilon^p $$

不变。

---

# 53. 当前一维代码与原论文 Eq. (66)–(71) 路线存在实现差异

当前一维 `pgd_enrichment.py` 已形成成熟 fixed-point enrichment framework，但其 spatial half-step 更接近：

$$ \boxed{ \text{fixed-temporal weighted residual minimisation} } $$

而原论文逐式路线为：

$$ \boxed{ \text{Eq. (66)} \rightarrow W \rightarrow \text{Eq. (70) FE solve} \rightarrow \text{Eq. (71) recovery} } $$

因此，未来 tower implementation 不应简单机械复制现有 1D spatial-enrichment function。

需要专门比较：

1. 两种 formulation 是否在特定 assumptions 下等价；
2. 当前 1D solver 的 numerical safeguards 哪些可以复用；
3. 哪些部分必须改为原论文 explicit $W$-based route；
4. 如何尽可能保持与原论文 Eq. (66)–(72) 一致。

---

# 54. 本阶段明确不能做的事情

1. 不把 $W$ 当成 physical tangent stiffness；
2. 不把 $K_W$ 当成 ordinary Newton tangent matrix；
3. 不把 $f_\delta$ 当成真实外载或 Newton residual；
4. 不把 $\bar{\tilde U}$ 当成真实 tower displacement mode；
5. 不把 Eq. (70) 解出的 $\bar{\tilde\varepsilon}$ 直接存入 PGD basis；
6. 不把 $\bar{\tilde\varepsilon}$ 与 $\bar\varepsilon^p$ 混淆；
7. 不忽略 Eq. (71) 的 $1/a$；
8. 不因为 cyclic loading 就设 $a=0$；
9. 不在 $a\approx0$ 时简单硬截断 denominator；
10. 不把 $W$ 永久存为 accepted basis property；
11. 不显式建立 dense $D_W$；
12. 不机械复制当前 1D spatial solve 而跳过 Eq. (70)–(71) 的 original-paper structure；
13. 不在 Eq. (71) 后立即进入代码，尚需完成 Eq. (72) temporal half-step；
14. 不提前引入 $n-\tau-x$ 或 cycle-jump formulation。

---

# 55. 本阶段已经解决的问题

通过 Eq. (70)–(71) 的逐式推导，目前已经明确：

1. Eq. (70) 如何由 Eq. (67)–(69) 得到；
2. Eq. (70) 为什么是 displacement-like spatial FE problem；
3. 为什么 $W$ 在 Eq. (70) 中扮演 effective stiffness-like operator；
4. 为什么 $W$ 仍不是 material tangent stiffness；
5. Eq. (70) 右端为什么是负号；
6. Eq. (70) RHS 的 projected-defect forcing 含义；
7. $\bar{\tilde U}$ 的准确身份；
8. 为什么 displacement formulation 比 direct stress constrained solve 更适合 FE implementation；
9. tower Eq. (70) 为什么得到 $K_W=H^TMD_WH$；
10. tower Eq. (70) RHS 为什么是 $-H^TMD_W\bar\delta$；
11. 为什么旧 $C_0$-based assembly skeleton 可以复用；
12. 为什么 material-point operator 必须从 $C_0$ 推广到 $D_W$；
13. 为什么 $D_W$ 不应以 dense matrix 存储；
14. 为什么 $W_q$ 会随 fixed-point temporal update 变化；
15. 为什么 Eq. (70) 仍然只是一次 pure spatial solve；
16. 为什么 $H_\sigma$ 正定性影响 $W$ 与 $K_W$ 的稳定性；
17. Eq. (70) solve 后如何恢复 $\bar{\tilde\varepsilon}$；
18. Eq. (70) solve 后如何恢复 $\bar\sigma$；
19. 为什么恢复的 $\bar\sigma$ 自动满足 static equilibrium；
20. 为什么 $\bar{\tilde\varepsilon}$ 自动 compatible；
21. Eq. (71) 为什么必须回到 Eq. (61) 的 temporal projection；
22. Eq. (71) 中 $1/\langle\dot\lambda\lambda\rangle$ 的来源；
23. Eq. (71) 的物理内核为什么仍是 total-elastic-plastic strain partition；
24. $\bar{\tilde\varepsilon}$ 在 Eq. (71) 后更准确的 scaled-total-strain interpretation；
25. 为什么 final mode 再次满足 canonical $C_0$-based relation；
26. 为什么 $W$ 不需要成为 permanent PGD basis data；
27. tower Eq. (71) 的 fiber-level vector form；
28. Eq. (71) 最终为什么回到 $\mathbb R^{320}$ fiber material-point space；
29. Eq. (70)–(71) 的完整 tower spatial half-step；
30. 应建立哪些 consistency checks；
31. 为什么 $a\approx0$ 是关键 numerical risk；
32. 为什么 cyclic external loading 不意味着 $a=0$；
33. denominator safeguard 的原则；
34. Eq. (71) 与 PGD scaling indeterminacy 的相容性；
35. spatial orthonormalisation 后 temporal rescaling 的必要性；
36. 当前 1D code 与 original-paper explicit $W$-based route 的实现差异。

---

# 56. 当前尚未解决的问题

下一阶段 Eq. (72) 必须重点解决：

1. 固定 $\bar\varepsilon^p_{m+1}$ 后如何重新求 $\lambda_{m+1}(t)$；
2. Eq. (72) 的 objective/residual 到底从哪一个 search-direction relation 推导；
3. Eq. (72) 是否属于 weighted least-squares / residual minimisation；
4. weighting 是否为 $H_\sigma^{-1}$；
5. Eq. (72) 与现有 `pgd_time_update.py` 的数学结构是否一致；
6. temporal discretization 如何处理 $\dot\lambda$；
7. temporal boundary/initial conditions 如何施加；
8. $a=\langle\dot\lambda\lambda\rangle$ 如何在 fixed-point iterations 中保持数值可用；
9. fixed-point pair convergence 应采用什么 criterion；
10. spatial mode normalization 与 temporal function update 应在什么时机进行；
11. new pair fixed-point convergence 后如何进入 basis orthonormalisation；
12. Eq. (72) 后是否需要重新计算 Eq. (70)–(71) 直到 pair 收敛；
13. 原论文 Eq. (72) 与 current 1D enrichment implementation 的差异；
14. 哪些 1D convergence safeguards 可以直接迁移到 tower。

---

# 57. 下一阶段工作入口

下一阶段严格从 Eq. (72) 开始。

当前 spatial half-step 已闭合为：

$$ \boxed{ \lambda^{(k)}(t) \rightarrow W^{(k)},\bar\delta^{(k)} \rightarrow \bar{\tilde U}^{(k+1)} \rightarrow \bar{\tilde\varepsilon}^{(k+1)},\bar\sigma^{(k+1)} \rightarrow \bar\varepsilon^{p,(k+1)} } $$

下一步需要解决：

$$ \boxed{ \bar\varepsilon^{p,(k+1)} \rightarrow \lambda^{(k+1)}(t) } $$

从而构成完整 fixed-point enrichment：

$$ \boxed{ \lambda^{(k)} \rightarrow \bar\varepsilon^{p,(k+1)} \rightarrow \lambda^{(k+1)} \rightarrow \cdots } $$

在 Eq. (72) 完成前，暂不进入代码修改。

---

# 58. 阶段结论

本阶段完成了原论文 new-mode enrichment spatial half-step 中最关键的两步。

Eq. (70)：

$$ \boxed{ \text{temporally projected LATIN defect} \rightarrow \text{effective spatial FE solve} } $$

在 tower 中对应：

$$ \boxed{ H^TMD_WH\vec{\bar{\tilde U}}=-H^TMD_W\vec{\bar\delta} } $$

Eq. (71)：

$$ \boxed{ \text{auxiliary compatible strain + equilibrated stress} \rightarrow \text{plastic-strain spatial mode} } $$

在 tower 中对应：

$$ \boxed{ \vec{\bar\varepsilon}^{\,p}=\frac{1}{a}H\vec{\bar{\tilde U}}-C_0^{-1}D_W\left(H\vec{\bar{\tilde U}}+\vec{\bar\delta}\right) } $$

因此，截至 Eq. (71)，原论文 new-mode enrichment 的 fixed-temporal spatial half-step 已经完整闭合。

对于当前 tower migration，仍然无需改变 original $x-t$ PGD architecture：

$$ \boxed{ x \rightarrow (e,g,f), \qquad \Delta\varepsilon^p(x,t)=\lambda(t)\bar\varepsilon^p(x) } $$

当前真正缺失的最后一部分，是 Eq. (72)：

$$ \boxed{ \text{fixed spatial mode} \rightarrow \text{updated temporal function} } $$

只有完成 Eq. (72)，new rank-one PGD pair 的 alternating fixed-point enrichment 才能在理论上完全闭合。
