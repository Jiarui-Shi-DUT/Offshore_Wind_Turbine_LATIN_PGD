# 海上风机塔筒 LATIN-PGD：Eq. (60) 饱和判据、Eq. (76)–(77) LATIN 范数及一维成熟控制逻辑迁移阶段总结

**日期：2026-08-16**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**当前分支：`feature/offshore-wind-turbine-tower-fatigue`**  
**前置总结：`docs/2026-08-15-tower-latin-pgd-global-stage-temporal-damage-hardening-stage-summary.md`**  
**前一稳定提交：`dcf091a` — `docs: summarize tower LATIN-PGD global stage derivation`**  
**本阶段范围：原论文 Eq. (60) saturation criterion、Eq. (76)–(77) LATIN relative indicator 与 mechanical norm 的物理解释和 tower fiber-beam 离散；同时回查一维三材料杆中已经形成的 adaptive PGD / LATIN convergence 控制逻辑，明确哪些问题已经解决并应直接迁移。**  
**下一阶段：原论文 Eq. (61) 开始的 new PGD mode enrichment。**

---

# 1. 本文档的阶段定位

上一份阶段总结已经完成 tower LATIN global stage 的理论闭合，包括：

- Eq. (58)–(59) fixed spatial basis temporal update；
- fiber-level $H_\sigma$；
- paper-separated forcing；
- damage residual branch；
- hardening variables；
- $D,\dot D,Y$；
- plastic / damage / hardening 三条 global-stage branch；
- one fixed reference tower equilibrium operator；
- current one-dimensional code 与 tower v1 formulation 的继承关系。

在此基础上，本阶段继续处理：

$$
\boxed{\text{Eq. (60) saturation criterion}}
$$

以及它依赖的：

$$
\boxed{\text{Eq. (76)–(77) LATIN indicator and mechanical norm}}
$$

。

本阶段最初计划进一步重新讨论：

$$
\zeta_{\mathrm{enrich}}\quad\text{与}\quad\zeta_{\mathrm{end}}
$$

的嵌套控制逻辑。

但通过回查当前一维三材料杆的 `pgd_solver.py`、`pgd_saturation.py` 和 `iteration_control.py`，确认该控制问题此前已经形成了经过调试的成熟方案。

因此本阶段最终结论是：

> **Eq. (60) 及其 tower norm migration 已经解决；adaptive control 的逻辑主体不需要在 tower 阶段重新发明，而应继承一维杆中已经验证的两层控制架构。**

---

# 2. 本阶段核心问题

本阶段实际回答了以下问题：

1. Eq. (60) 中 $\zeta$ 到底衡量什么；
2. $\xi_i$ 与 $\xi_{i+1}$ 到底是什么；
3. $\zeta$ 是不是 PGD truncation error；
4. $\zeta$ 与 Eq. (59) reduced residual 的关系；
5. $\zeta=0.1$ 的定量含义；
6. $\zeta<0$ 应如何理解；
7. Eq. (76)–(77) 中的 LATIN mechanical state 包含哪些变量；
8. 为什么 $D,\dot D,Y$ 不显式进入 Eq. (77)；
9. Eq. (77) 如何从 continuum/bar 空间积分迁移为 fiber-beam material-point quadrature；
10. tower convergence norm 为什么不能改成节点位移、截面力或弯矩范数；
11. time integration 是否需要在 tower 阶段更改；
12. material-point quadrature volume 如何统一服务 Eq. (59) 和 Eq. (77)；
13. $\zeta_{\mathrm{enrich}}$、$\zeta_{\mathrm{end}}$ 的冲突是否已经在一维杆中解决；
14. tower v1 应继承什么控制逻辑；
15. Eq. (60) 阶段何时可以结束并进入 Eq. (61)。

---

# 3. Eq. (60) 在 hybrid PGD 中的位置

原论文的 hybrid PGD global-stage strategy 可以概括为：

```text
existing m spatial modes
        ↓
keep spatial basis fixed
        ↓
Eq. (58)–(59)
update temporal coefficients
        ↓
construct current LATIN global candidate
        ↓
evaluate LATIN indicator xi
        ↓
Eq. (60) saturation parameter zeta
        ↓
judge basis adequacy
        ↓
sufficient improvement → advance LATIN
insufficient improvement → enrich PGD basis
```

因此 Eq. (60) 并不是独立存在的 stopping formula。

它处于：

$$
\boxed{\text{fixed-basis temporal reuse}}
$$

和：

$$
\boxed{\text{new-mode enrichment}}
$$

之间。

---

# 4. 原论文 Eq. (60)

定义：

$$
\boxed{\zeta_i=\frac{\xi_i-\xi_{i+1}}{\xi_i+\xi_{i+1}}}
$$

。

原论文基本意图：

- 若当前 iteration 的 LATIN indicator 得到足够改善，则已有 PGD basis 仍然有效；
- 若改善不足，则认为当前 basis 已经不足，需要新增一个 space-time pair。

因此 $\zeta$ 是一个 adaptive PGD control parameter。

---

# 5. $\xi$ 的真正含义

原论文 Eq. (76)：

$$
\boxed{\xi_i=\frac{\|\hat s^p_{i+1/2}-s^p_{i+1}\|}{\|\hat s^p_{i+1/2}\|+\|s^p_{i+1}\|}}
$$

其中：

- $\hat s^p_{i+1/2}$：local-stage mechanical state；
- $s^p_{i+1}$：global-stage mechanical state。

所以：

$$
\boxed{\xi=\text{relative distance between local and global LATIN states}}
$$

。

当：

$$
\xi\to0
$$

意味着：

$$
\hat s^p\approx s^p
$$

也就是两个 LATIN admissibility manifolds 逐渐相交。

因此：

$$
\boxed{\xi=\text{LATIN nonlinear convergence indicator}}
$$

。

---

# 6. Eq. (60) 中 $\xi_i$ 和 $\xi_{i+1}$ 的精确理解

一个关键修正是：

> 不应把 Eq. (60) 简化表述为“Eq. (59) reduced residual 更新前后的两个 residual”。

Eq. (60) 比较的是 successive LATIN indicators：

$$
\xi_i\longrightarrow\xi_{i+1}
$$

。

因此：

$$
\boxed{\zeta_i=\text{successive LATIN-indicator relative improvement}}
$$

而不是：

$$
\boxed{\zeta_i=\text{PGD projection error}}
$$

。

---

# 7. $\zeta$ 为什么又可以用来判断 basis adequacy

在当前 hybrid strategy 中，global candidate 的主要 reduced approximation 来自：

$$
\sum_{j=1}^{m}\lambda_j(t)\bar\varepsilon_j^p(x)
$$

。

若固定已有 spatial basis，只更新时间函数以后：

$$
\xi_{i+1}\ll\xi_i
$$

则说明当前 reduced subspace 仍能明显改善 LATIN solution。

因此无需新增 spatial mode。

若：

$$
\xi_{i+1}\approx\xi_i
$$

则说明仅靠已有 spatial basis 已无法继续有效改善当前 global stage。

于是作者以：

$$
\zeta\text{ small}
$$

作为：

$$
\boxed{\text{current PGD basis insufficient / saturated}}
$$

的代理判据。

所以更准确的表述是：

> **$\zeta$ 不是直接测量 PGD truncation error，而是通过 LATIN convergence improvement 间接判断 reduced basis adequacy。**

---

# 8. $\xi$ 与 $\zeta$ 分别回答什么问题

$\xi$ 回答：

$$
\boxed{\text{当前 local/global states 相距多远？}}
$$

$\zeta$ 回答：

$$
\boxed{\text{这个距离相对于上一 accepted LATIN state 改善得有多快？}}
$$

因此逻辑链：

$$
\boxed{\xi\longrightarrow\text{successive improvement}\longrightarrow\zeta\longrightarrow\text{basis decision}}
$$

。

---

# 9. $\zeta$ 的数学变换

定义：

$$
r_i=\frac{\xi_{i+1}}{\xi_i}
$$

则：

$$
\boxed{\zeta_i=\frac{1-r_i}{1+r_i}}
$$

反解：

$$
\boxed{r_i=\frac{1-\zeta_i}{1+\zeta_i}}
$$

。

---

# 10. $\zeta>0$

若：

$$
\xi_{i+1}<\xi_i
$$

则：

$$
0<r_i<1
$$

所以：

$$
\boxed{0<\zeta_i<1}
$$

。

$\zeta$ 越大，说明 LATIN indicator 改善越明显。

---

# 11. $\zeta=0$

若：

$$
\xi_{i+1}=\xi_i
$$

则：

$$
\boxed{\zeta_i=0}
$$

。

表示：

$$
\boxed{\text{no LATIN progress}}
$$

。

---

# 12. $\zeta<0$

若：

$$
\xi_{i+1}>\xi_i
$$

则：

$$
\boxed{\zeta_i<0}
$$

。

这表示：

$$
\boxed{\text{LATIN indicator worsened}}
$$

。

因此 negative $\zeta$ 不能解释为 convergence saturation。

在当前一维实现中，negative $\zeta$ 被视为 basis insufficient / enrichment-needed，而不是 stopping success。

---

# 13. $\zeta$ 的取值范围

若：

$$
\xi_i\ge0,\qquad\xi_{i+1}\ge0
$$

且分母非零，则：

$$
\boxed{-1\le\zeta_i\le1}
$$

。

极限：

$$
\xi_{i+1}=0\Rightarrow\zeta_i=1
$$

表示一步把 LATIN indicator 降到零。

若：

$$
\xi_{i+1}\gg\xi_i
$$

则：

$$
\zeta_i\to-1
$$

表示严重恶化。

---

# 14. 为什么使用 symmetric relative change

论文不是使用：

$$
\frac{\xi_i-\xi_{i+1}}{\xi_i}
$$

而是：

$$
\frac{\xi_i-\xi_{i+1}}{\xi_i+\xi_{i+1}}
$$

。

其数学优点包括：

- symmetric normalisation；
- 不只把旧值作为 reference；
- 自然有界于 $[-1,1]$；
- 对不同绝对尺度的 $\xi$ 更容易比较。

因此可以把 $\zeta$ 理解为：

$$
\boxed{\text{symmetric relative LATIN-indicator change}}
$$

。

---

# 15. $\zeta_{\mathrm{enrich}}=0.1$ 的定量含义

由：

$$
r=\frac{1-\zeta}{1+\zeta}
$$

代入：

$$
\zeta=0.1
$$

得到：

$$
r=\frac{0.9}{1.1}\approx0.81818
$$

。

因此：

$$
\boxed{\zeta>0.1\Longleftrightarrow\xi_{i+1}\lesssim0.818\xi_i}
$$

即大致要求当前 iteration 带来超过：

$$
\boxed{18.2\%}
$$

的 relative LATIN-indicator reduction。

---

# 16. 数值例：显著改善

设：

$$
\xi_i=0.10,\qquad\xi_{i+1}=0.05
$$

则：

$$
\zeta=\frac{0.10-0.05}{0.10+0.05}=0.333
$$

。

因为：

$$
0.333>0.1
$$

所以已有 basis 的改善足够明显。

---

# 17. 数值例：改善不足

设：

$$
\xi_i=0.10,\qquad\xi_{i+1}=0.09
$$

则：

$$
\zeta=\frac{0.01}{0.19}\approx0.0526
$$

。

虽然 $\xi$ 下降，但：

$$
0.0526<0.1
$$

说明 improvement insufficient。

---

# 18. 数值例：indicator 恶化

设：

$$
\xi_i=0.10,\qquad\xi_{i+1}=0.101
$$

则：

$$
\zeta\approx-0.00498
$$

。

这意味着：

$$
\boxed{\text{current reduced approximation worsened the LATIN indicator}}
$$

。

---

# 19. Eq. (59) residual 与 Eq. (60) $\zeta$ 的层级区别

Eq. (59) 的 reduced residual：

$$
r_{\mathrm{PGD}}=P\Delta\dot\lambda-D_HS\Delta\lambda-f
$$

回答：

> 在 fixed spatial basis 中，当前 temporal coefficients 对 plastic search-direction equation 拟合得怎么样？

所以它属于：

$$
\boxed{\text{inner reduced global-stage problem}}
$$

。

Eq. (60)：

$$
\zeta=\frac{\xi_i-\xi_{i+1}}{\xi_i+\xi_{i+1}}
$$

回答：

> 当前整个 LATIN global candidate 相对于上一 accepted state 是否取得足够进展？

所以它属于：

$$
\boxed{\text{outer adaptive PGD control}}
$$

。

二者不能混为一个指标。

---

# 20. Tower Eq. (60) 公式本身是否需要修改

结论：

$$
\boxed{\text{No}}
$$

。

Tower 直接使用：

$$
\boxed{\zeta_i^{\mathrm{tower}}=\frac{\xi_i^{\mathrm{tower}}-\xi_{i+1}^{\mathrm{tower}}}{\xi_i^{\mathrm{tower}}+\xi_{i+1}^{\mathrm{tower}}}}
$$

。

真正需要迁移的是：

$$
\boxed{\xi^{\mathrm{tower}}}
$$

中的 mechanical norm。

---

# 21. Eq. (76) 的 tower form

定义：

$$
N_{\mathrm{diff}}=\|\hat s^p_{i+1/2}-s^p_{i+1}\|_{\mathrm{tower}}
$$

$$
N_{\mathrm{local}}=\|\hat s^p_{i+1/2}\|_{\mathrm{tower}}
$$

$$
N_{\mathrm{global}}=\|s^p_{i+1}\|_{\mathrm{tower}}
$$

则：

$$
\boxed{\xi_{i+1}^{\mathrm{tower}}=\frac{N_{\mathrm{diff}}}{N_{\mathrm{local}}+N_{\mathrm{global}}}}
$$

。

所以关键变成：

$$
\boxed{\|\cdot\|_{\mathrm{tower}}}
$$

如何定义。

---

# 22. Eq. (77) 的 mechanical state

进入原论文 mechanical norm 的 state fields 可以展开为：

$$
\boxed{\sigma,\beta,\bar R,\dot\varepsilon^p,\varepsilon^e,\dot\alpha,\dot{\bar r}}
$$

。

其中：

$$
X=(\alpha,\bar r)
$$

$$
Z=(\beta,\bar R)
$$

。

当前 `latin/iteration_control.py` 也正是通过 `_plastic_fields()` 返回这七个 fields。

---

# 23. Eq. (77) 的 continuum structure

原论文 mechanical norm：

$$
\boxed{\|s^p\|^2=\int_0^T\int_\Omega\left[\sigma:H_\sigma\sigma+Z:H_ZZ+\dot\varepsilon^p:H_\sigma^{-1}\dot\varepsilon^p+\varepsilon^e:C_0\varepsilon^e+\dot X:H_Z^{-1}\dot X\right]d\Omega\,dt}
$$

。

Diagonal search direction：

$$
H_Z=\operatorname{diag}(H_\beta,H_{\bar R})
$$

。

---

# 24. 为什么这些项成对出现

Plastic conjugate pair：

$$
\boxed{\sigma\leftrightarrow\dot\varepsilon^p}
$$

对应：

$$
\sigma:H_\sigma\sigma
$$

与：

$$
\dot\varepsilon^p:H_\sigma^{-1}\dot\varepsilon^p
$$

。

Hardening conjugate pair：

$$
\boxed{Z\leftrightarrow\dot X}
$$

对应：

$$
Z:H_ZZ
$$

与：

$$
\dot X:H_Z^{-1}\dot X
$$

。

Elastic metric：

$$
\boxed{\varepsilon^e:C_0\varepsilon^e}
$$

。

---

# 25. 为什么 $D,\dot D,Y$ 不显式进入 Eq. (77)

完整 LATIN state 中存在：

$$
D,\dot D,Y
$$

。

但当前 original formulation：

$$
\boxed{b^-=0}
$$

。

因此 damage conjugate pair 没有形成类似：

$$
Yb^-Y
$$

以及：

$$
\dot D(b^-)^{-1}\dot D
$$

的正定 norm terms。

所以：

$$
\boxed{D,\dot D,Y\text{ do not appear as independent Eq. (77) components}}
$$

。

但 damage 并未被忽略，因为它会影响：

- $\sigma$；
- $\varepsilon^e$；
- $H_\sigma$；
- $H_\beta$；
- $H_{\bar R}$。

---

# 26. Fiber material 中 tensor norm 的 scalar reduction

当前 tower fiber material law 是 scalar axial model。

因此：

$$
\sigma:H_\sigma\sigma\longrightarrow H_\sigma\sigma^2
$$

。

Plastic-rate term：

$$
\dot\varepsilon^p:H_\sigma^{-1}\dot\varepsilon^p\longrightarrow\frac{(\dot\varepsilon^p)^2}{H_\sigma}
$$

。

Elastic term：

$$
\varepsilon^e:C_0\varepsilon^e\longrightarrow E_0(\varepsilon^e)^2
$$

。

---

# 27. Hardening scalar reduction

由于：

$$
Z=\begin{bmatrix}\beta\\\bar R\end{bmatrix}
$$

所以：

$$
\boxed{Z^TH_ZZ=H_\beta\beta^2+H_{\bar R}\bar R^2}
$$

。

又因为：

$$
\dot X=\begin{bmatrix}\dot\alpha\\\dot{\bar r}\end{bmatrix}
$$

所以：

$$
\boxed{\dot X^TH_Z^{-1}\dot X=\frac{\dot\alpha^2}{H_\beta}+\frac{\dot{\bar r}^2}{H_{\bar R}}}
$$

。

---

# 28. 单个 tower fiber material point 的 norm density

定义：

$$
\mathcal N_q(t)
$$

则：

$$
\boxed{\mathcal N_q=H_{\sigma,q}\sigma_q^2+H_{\beta,q}\beta_q^2+H_{\bar R,q}\bar R_q^2+\frac{(\dot\varepsilon_q^p)^2}{H_{\sigma,q}}+E_{0,q}(\varepsilon_q^e)^2+\frac{\dot\alpha_q^2}{H_{\beta,q}}+\frac{\dot{\bar r}_q^2}{H_{\bar R,q}}}
$$

。

这就是 Eq. (77) 在 scalar fiber material-point space 中的直接 integrand。

---

# 29. 为什么 elastic metric 使用 $E_0$

Tower v1 保持：

$$
\boxed{C_0=E_0}
$$

作为 reference elastic metric。

因此 Eq. (77) 使用：

$$
\boxed{E_0(\varepsilon^e)^2}
$$

而不是 damaged modulus。

原因与 $\mathcal E_{\mathrm{tower}}$ 相同：global admissibility / LATIN metric 使用 fixed positive reference elasticity，而不是随 damage 退化的 current tangent。

---

# 30. Damage 仍然通过 elastic strain 进入 norm

例如 tension：

$$
\varepsilon^e=\frac{\sigma}{E_0(1-D)}
$$

则：

$$
E_0(\varepsilon^e)^2=\frac{\sigma^2}{E_0(1-D)^2}
$$

。

所以：

$$
\boxed{\text{fixed metric does not mean damage has no effect}}
$$

。

---

# 31. 从 continuum 空间积分到 fiber-beam quadrature

Continuum：

$$
\int_\Omega(\cdot)d\Omega
$$

对于 beam：

$$
d\Omega=dA\,ds
$$

因此：

$$
\int_\Omega(\cdot)d\Omega=\int_0^L\int_A(\cdot)dA\,ds
$$

。

Fiber beam-column 已分别离散 section integration 与 beam-axis integration。

---

# 32. Section fiber quadrature

在 beam integration point $(e,g)$：

$$
\boxed{\int_{A_{eg}}\mathcal N\,dA\approx\sum_fA_{egf}\mathcal N_{egf}}
$$

。

这里权重是：

$$
A_{egf}
$$

。

Fiber coordinate $y_{egf}$ 已经通过：

$$
\varepsilon_{egf}=\varepsilon_{0,eg}-y_{egf}\kappa_{eg}
$$

进入材料状态，不再作为 norm integration weight。

---

# 33. Beam-axis quadrature

Element $e$：

$$
\boxed{\int_{L_e}(\cdot)ds\approx\sum_gw_gJ_e(\cdot)_{eg}}
$$

。

若标准 parent coordinate：

$$
\xi\in[-1,1]
$$

则通常：

$$
J_e=\frac{L_e}{2}
$$

。

---

# 34. Tower material-point volume weight

结合两层 quadrature：

$$
\boxed{v_{egf}=A_{egf}w_gJ_e}
$$

。

量纲：

$$
[v_{egf}]=L^3
$$

。

所以它是每个 discrete fiber material point 所代表的 physical volume weight。

---

# 35. Tower Eq. (77)

最终：

$$
\boxed{\|s^p\|_{\mathrm{tower}}^2=\int_0^T\sum_{e,g,f}v_{egf}\mathcal N_{egf}(t)\,dt}
$$

其中：

$$
\boxed{v_{egf}=A_{egf}w_gJ_e}
$$

。

Flatten：

$$
q=q(e,g,f)
$$

后：

$$
\boxed{\|s^p\|_{\mathrm{tower}}^2=\int_0^T\sum_{q=1}^{N_q}v_q\mathcal N_q(t)\,dt}
$$

。

---

# 36. 当前 coarse tower 尺寸

当前：

$$
N_e=10
$$

$$
N_g=2
$$

$$
N_f=16
$$

所以：

$$
\boxed{N_q=320}
$$

。

因此：

```text
state field            : (N_t, 320)
search-direction field : (N_t, 320)
volume weights         : (320,)
```

。

---

# 37. Local/global difference norm

定义：

$$
\delta\sigma=\hat\sigma-\sigma
$$

$$
\delta\beta=\hat\beta-\beta
$$

$$
\delta\bar R=\hat{\bar R}-\bar R
$$

$$
\delta\dot\varepsilon^p=\hat{\dot\varepsilon}^p-\dot\varepsilon^p
$$

$$
\delta\varepsilon^e=\hat\varepsilon^e-\varepsilon^e
$$

$$
\delta\dot\alpha=\hat{\dot\alpha}-\dot\alpha
$$

$$
\delta\dot{\bar r}=\hat{\dot{\bar r}}-\dot{\bar r}
$$

。

Then：

$$
\boxed{\|\hat s^p-s^p\|_{\mathrm{tower}}^2=\int_0^T\sum_qv_q\left[H_{\sigma,q}(\delta\sigma_q)^2+H_{\beta,q}(\delta\beta_q)^2+H_{\bar R,q}(\delta\bar R_q)^2+\frac{(\delta\dot\varepsilon_q^p)^2}{H_{\sigma,q}}+E_0(\delta\varepsilon_q^e)^2+\frac{(\delta\dot\alpha_q)^2}{H_{\beta,q}}+\frac{(\delta\dot{\bar r}_q)^2}{H_{\bar R,q}}\right]dt}
$$

。

---

# 38. 同一 LATIN iteration 必须使用同一套 search directions

计算：

$$
\|\hat s^p\|
$$

$$
\|s^p\|
$$

以及：

$$
\|\hat s^p-s^p\|
$$

时，必须使用相同的：

$$
H_{\sigma,i+1/2},H_{\beta,i+1/2},H_{\bar R,i+1/2}
$$

。

当前一维 `relative_latin_indicator()` 已经严格遵守这一点。

---

# 39. Tower Eq. (76) 完整闭合

$$
\boxed{\xi_{i+1}^{\mathrm{tower}}=\frac{\|\hat s^p_{i+1/2}-s^p_{i+1}\|_{\mathrm{tower}}}{\|\hat s^p_{i+1/2}\|_{\mathrm{tower}}+\|s^p_{i+1}\|_{\mathrm{tower}}}}
$$

。

至此：

$$
\boxed{\xi^{\mathrm{tower}}}
$$

已经可以直接计算。

---

# 40. 为什么不能使用 section force norm

不能把 Eq. (77) 改成：

$$
N^2+M^2
$$

。

因为 $N,M$ 是 fiber stresses 的 section resultants，会发生正负相消。

特别是在 pure bending 中：

$$
N\approx0
$$

但 tensile / compressive fiber stresses 很大。

因此 Eq. (77) 必须保留 fiber-level material-state integration。

---

# 41. 为什么不能使用 nodal displacement norm

Tower-top displacement 可以作为工程响应 validation quantity。

但 Eq. (76) 的 $\xi$ 是：

$$
\boxed{\text{LATIN constitutive state-space convergence indicator}}
$$

。

所以 $u_{\mathrm{top}}$ 不能代替 $\xi$。

---

# 42. Time integration

当前一维 `iteration_control.py` 使用：

```text
space_integral = density @ element_volumes
value = np.trapz(space_integral, x=time)
```

。

因此 tower v1 最直接迁移是：

$$
I_n=\sum_qv_q\mathcal N_{nq}
$$

然后：

$$
\boxed{\|s^p\|^2\approx\operatorname{trapz}(I_n,t_n)}
$$

。

---

# 43. 是否现在改成严格 DG0 rectangular quadrature

本阶段决定：

$$
\boxed{\text{No}}
$$

。

原因：

- 当前目标是 spatial migration；
- 一维复现已经使用 trapezoidal diagnostic integration；
- 不应在 tower migration 同时改变 time norm quadrature；
- DG0 与 backward-Euler / trapz 的严格一致性可以作为以后单独理论问题。

因此 tower v1 首先继承 current one-dimensional time integration。

---

# 44. Eq. (59) 与 Eq. (77) 的 quadrature consistency

Eq. (59) temporal reduced problem 使用：

$$
\boxed{v_q=A_{egf}w_gJ_e}
$$

Eq. (77) LATIN norm 同样使用：

$$
\boxed{v_q=A_{egf}w_gJ_e}
$$

。

因此未来应建立唯一：

```text
material_point_volumes[q]
```

同时供：

```text
pgd_time_update
latin_indicator
reduced_residual_norm
```

使用。

---

# 45. Geometry sanity check

建议未来检查：

$$
\boxed{\sum_qv_q\approx V_{\mathrm{tower}}}
$$

。

Thin-walled tower 连续体积：

$$
V_{\mathrm{tower}}=\int_0^LA(s)\,ds
$$

近似：

$$
A(s)\approx\pi D(s)t(s)
$$

因此：

$$
\boxed{\sum_{e,g,f}A_{egf}w_gJ_e\approx\int_0^L\pi D(s)t(s)\,ds}
$$

。

---

# 46. Taper 对 fiber area 的影响

Tower section 随高度变化：

$$
D=D(s),\qquad t=t(s)
$$

。

因此严格应该使用：

$$
\boxed{A_{egf}}
$$

而不是假定所有 integration points 共用 $A_f$。

---

# 47. Tower norm 的矩阵实现

定义：

$$
\mathcal N\in\mathbb R^{N_t\times N_q}
$$

和：

$$
v=[v_1,\ldots,v_{N_q}]^T
$$

。

每个时间点：

$$
\boxed{I_n=\mathcal N_{n,:}v}
$$

然后：

$$
\boxed{\|s^p\|^2=\operatorname{trapz}(I,t)}
$$

。

这与一维当前实现 `density @ element_volumes` 完全同构。

---

# 48. 一维到 tower 的本质代码变化

一维：

$$
\boxed{v_e=AL_e}
$$

Tower：

$$
\boxed{v_q=A_{egf}w_gJ_e}
$$

。

因此 `iteration_control.py` 的理论核心无需重写。

本质变化只有：

```text
second dimension:
    n_elements
        ↓
    n_material_points
```

和：

```text
element_volumes
        ↓
material_point_volumes
```

。

---

# 49. Damage diagnostic 的定位

虽然原论文 Eq. (77) 不显式加入 damage field，future tower robustness 可额外记录：

$$
\eta_D=\frac{\|\hat D-D\|_V}{\|\hat D\|_V+\|D\|_V+\epsilon}
$$

作为：

$$
\boxed{\text{diagnostic only}}
$$

。

但不能用它替代 $\xi$。

---

# 50. Eq. (60) 的 tower complete chain

```text
local state hat_s
        ↓
global state s
        ↓
fiber mechanical norm
        ↓
xi_tower
        ↓
compare with previous accepted xi
        ↓
zeta_tower
        ↓
adaptive PGD basis decision
```

数学上：

$$
\boxed{\hat s^p,s^p\longrightarrow\|\cdot\|_{\mathrm{tower}}\longrightarrow\xi^{\mathrm{tower}}\longrightarrow\zeta^{\mathrm{tower}}}
$$

。

---

# 51. 回查一维三材料杆：adaptive control 是否已解决

答案：

$$
\boxed{\text{Yes}}
$$

。

当前一维 `pgd_solver.py` 已经明确：saturation parameter 只用于判断 current PGD basis 是否需要 enrichment，而不是 nonlinear LATIN solve 的唯一 stopping criterion。

因此一维复现已经把：

$$
\boxed{\text{basis saturation}}
$$

与：

$$
\boxed{\text{nonlinear convergence}}
$$

分开。

---

# 52. 一维成熟方案中的 $\xi_i$ reference

当前 solver 明确规定：

$$
\boxed{\xi_i=\text{previous accepted LATIN indicator}}
$$

。

在当前 LATIN global stage 内，如果连续 enrichment：

$$
m\to m+1\to m+2\to\cdots
$$

baseline：

$$
\xi_i
$$

保持不变。

每个 trial candidate 产生：

$$
\xi_{i+1}^{\mathrm{trial}}
$$

再计算：

$$
\boxed{\zeta^{\mathrm{trial}}=\frac{\xi_i-\xi_{i+1}^{\mathrm{trial}}}{\xi_i+\xi_{i+1}^{\mathrm{trial}}}}
$$

。

---

# 53. First LATIN iteration 的 baseline

当前一维 solver：

$$
\boxed{\xi_0=1.0}
$$

作为 natural normalized initial reference。

并且空 spatial basis 无法做 time-only update，因此：

> **第一次 global stage 先强制生成第一个 PGD pair，再评估 $\zeta$。**

这一逻辑 tower v1 也应继承。

---

# 54. 一维 solver 的第一层：absolute LATIN convergence

当前优先判断：

$$
\boxed{\xi_{i+1}\le\xi_{\mathrm{tol}}}
$$

默认：

$$
\boxed{\xi_{\mathrm{tol}}=10^{-4}}
$$

。

若满足，则：

$$
\boxed{\text{complete nonlinear LATIN solve converged}}
$$

直接结束。

这比 $\zeta$ enrichment decision 优先。

---

# 55. 一维 solver 的第二层：basis improvement

若：

$$
\xi_{i+1}>\xi_{\mathrm{tol}}
$$

再计算 $\zeta$。

当：

$$
\boxed{\zeta>\zeta_{\mathrm{enrich}}}
$$

默认：

$$
\zeta_{\mathrm{enrich}}=0.1
$$

则：

$$
\boxed{\text{accept current global candidate and advance to next LATIN iteration}}
$$

。

---

# 56. $\zeta\le0.1$ 时的处理

当前 mature one-dimensional solver 已经不再使用：

> small $\zeta$ → immediately stop nonlinear solve

这种机械逻辑。

相反：

$$
\boxed{\zeta\le0.1}
$$

原则上表示：

$$
\boxed{\text{current basis improvement insufficient}}
$$

所以优先考虑 enrichment。

Negative $\zeta$ 同样属于 inadequate basis / worsened candidate，而不是 successful saturation。

---

# 57. Reduced residual protection

一维调试中发现一个重要 corner case：

$$
\zeta\text{ requests enrichment}
$$

但：

$$
\boxed{r_{\mathrm{red}}\le r_{\mathrm{tol}}}
$$

。

这意味着当前 reduced global problem 已经求到 prescribed tolerance。

此时已经没有足够 residual 可以生成有意义的新 PGD mode。

若仍强行 enrichment，容易产生：

- nearly zero spatial mode；
- tiny temporal mode；
- noise mode；
- false enrichment failure。

因此当前 solver 采用：

$$
\boxed{r_{\mathrm{red}}\le r_{\mathrm{tol}}\Rightarrow\text{accept candidate and advance LATIN}}
$$

。

---

# 58. 为什么 reduced residual 与 $\zeta$ 必须联合判断

$\zeta$ 衡量：

$$
\boxed{\text{outer LATIN improvement}}
$$

而 reduced residual 衡量：

$$
\boxed{\text{inner reduced global problem solvability with current basis}}
$$

。

因此可能出现：

```text
zeta small
but
reduced residual already tiny
```

此时 outer nonlinear progress 小，不一定意味着还能通过增加 PGD mode 继续改善当前 linearised global stage。

---

# 59. Persistent stagnation criterion

当前一维 solver 还独立定义 nonlinear stagnation control。

默认：

$$
\boxed{\xi_{\mathrm{stag}}=10^{-3}}
$$

$$
\boxed{\Delta\xi_{\mathrm{abs,tol}}=10^{-6}}
$$

连续要求：

$$
\boxed{N_{\mathrm{stag}}=3}
$$

轮 accepted LATIN iterations。

也就是如果：

$$
\xi_{i+1}\le10^{-3}
$$

且：

$$
|\xi_{i+1}-\xi_i|\le10^{-6}
$$

连续三轮，则可以：

$$
\boxed{\text{terminate as persistent nonlinear stagnation}}
$$

。

---

# 60. 为什么 persistent stagnation 比直接用 $\zeta_{\mathrm{end}}$ 更合理

论文给出：

$$
\zeta_{\mathrm{enrich}}=0.1
$$

以及：

$$
\zeta_{\mathrm{end}}\sim10^{-4}
$$

容易形成逻辑嵌套。

因为：

$$
\zeta<10^{-4}
$$

当然也满足：

$$
\zeta<0.1
$$

。

一维成熟实现把整体 nonlinear stopping 改为：

$$
\boxed{\text{absolute xi convergence + persistent accepted-state stagnation}}
$$

因此避免了仅凭一个非常小 $\zeta$ 直接结束 nonlinear solver 的歧义。

---

# 61. `pgd_saturation.py` 与 outer solver 的角色区别

当前 `pgd_saturation.py` 仍定义：

```text
ADVANCE_LATIN
ENRICH_BASIS
STOP_SATURATED
```

并保留：

```text
enrichment_tolerance = 0.1
stopping_tolerance = 1e-4
```

。

但 `pgd_solver.py` 的成熟 outer logic 已明确：

> **small saturation value no longer terminates the nonlinear solve by itself.**

也就是说：

$$
\boxed{\text{final solver semantics are determined by pgd_solver.py, not by the enum label alone}}
$$

。

这一点在未来 tower implementation 中必须保持清楚。

---

# 62. 当前一维 adaptive architecture 的两层结构

第一层：

$$
\boxed{\text{PGD adaptive basis control}}
$$

主要看：

$$
\zeta
$$

加：

$$
r_{\mathrm{red}}
$$

。

它回答：

> 当前 global stage 是否需要继续增加 spatial mode？

第二层：

$$
\boxed{\text{LATIN nonlinear convergence control}}
$$

主要看：

$$
\xi
$$

以及 persistent stagnation。

它回答：

> 整个 nonlinear LATIN solve 是否结束？

---

# 63. 一维成熟流程的简化伪代码

```text
previous accepted global state
previous accepted indicator xi_i
        ↓
local stage
        ↓
search directions
        ↓
global stage with current PGD basis
        ↓
relax global candidate
        ↓
compute xi_trial
        ↓
if xi_trial <= xi_tol:
    accept
    terminate CONVERGED

else:
    compute zeta(xi_i, xi_trial)

    if zeta sufficiently large:
        accept
        advance LATIN

    else:
        if reduced residual <= reduced_tol:
            accept
            advance LATIN
        else:
            enrich basis by one mode
            recompute trial candidate
            keep xi_i fixed
            repeat

after accepted iteration:
    update persistent stagnation counter
    if persistent stagnation reached:
        terminate STAGNATED
```

。

---

# 64. 其他一维 adaptive safeguards

当前 outer solver 还具有：

- `max_iterations`；
- `max_enrichments_per_iteration`；
- `ENRICHMENT_FAILED`；
- `MAX_ITERATIONS`；
- `minimum_spatial_norm`；
- `acceptance_tolerance`；
- `rcond`；
- relaxed candidate evaluation；
- fixed-point tolerance；
- non-finite checks。

这些属于成熟的一维 numerical-control layer。

Tower v1 不应随意删掉。

---

# 65. Relaxation 的位置

当前一维流程不是直接对 raw global candidate 计算 $\xi$。

而是先：

$$
\boxed{s_{i+1}=\mu\breve s_{i+1}+(1-\mu)s_i}
$$

默认：

$$
\boxed{\mu=0.8}
$$

然后对 relaxed candidate 计算 $\xi_{i+1}$。

因此 tower saturation chain 也应保持：

```text
raw global-stage candidate
        ↓
LATIN relaxation
        ↓
xi
        ↓
zeta
```

。

---

# 66. 为什么 tower 不需要重新研究 $\zeta_{\mathrm{end}}$ 逻辑

因为一维杆已经把论文较模糊的 stopping / enrichment 关系转化为一个更清晰的 adaptive architecture：

$$
\boxed{\text{zeta controls basis adequacy}}
$$

而：

$$
\boxed{\text{xi + persistent stagnation control nonlinear termination}}
$$

。

Tower 阶段重新设计这一层会导致：

- 同时改变 spatial discretisation 和 outer algorithm；
- 失去与一维 validated reproduction 的可比性；
- 增加不必要的 debug dimensions。

因此：

$$
\boxed{\text{reuse the mature 1D control logic}}
$$

是当前最稳妥路线。

---

# 67. Tower v1 对一维 adaptive control 的继承内容

直接继承：

1. $\xi_0=1$；
2. first pair forced when basis empty；
3. relaxed candidate；
4. Eq. (76) relative indicator；
5. Eq. (60) $\zeta$；
6. $\zeta_{\mathrm{enrich}}=0.1$ 作为 initial reference threshold；
7. absolute $\xi$ tolerance；
8. reduced residual protection；
9. enrichment loop with fixed baseline $\xi_i$；
10. max enrichments per LATIN iteration；
11. persistent stagnation；
12. max LATIN iterations；
13. enrichment failure safeguards。

---

# 68. Tower 阶段真正需要修改的部分

主要不是 control logic，而是：

$$
\boxed{\text{state discretisation and spatial integration}}
$$

。

一维：

```text
field[t, element]
element_volumes = A * L_e
```

Tower：

```text
field[t, q]
material_point_volumes[q] = A_egf * w_g * J_e
```

。

因此：

$$
\boxed{\text{control flow inherited, norm discretisation migrated}}
$$

。

---

# 69. Eq. (60) 阶段的最终理论 closure

现在完整链条：

$$
\boxed{\hat s^p,s^p}
$$

通过：

$$
\boxed{\|\cdot\|_{\mathrm{tower}}}
$$

得到：

$$
\boxed{\xi^{\mathrm{tower}}}
$$

通过 successive accepted indicator comparison：

$$
\boxed{\zeta^{\mathrm{tower}}}
$$

再结合：

$$
\boxed{r_{\mathrm{red}}}
$$

完成 basis adequacy control。

同时 nonlinear termination 独立由：

$$
\boxed{\xi\text{ absolute tolerance + persistent stagnation}}
$$

管理。

所以：

$$
\boxed{\text{Eq. (60) stage is closed}}
$$

。

---

# 70. 本阶段最重要的方法学结论

不能把以下三件事混为一谈：

### A. Reduced residual

$$
r_{\mathrm{red}}
$$

表示：

$$
\boxed{\text{current reduced global problem residual}}
$$

。

### B. LATIN indicator

$$
\xi
$$

表示：

$$
\boxed{\text{local/global state distance}}
$$

。

### C. Saturation parameter

$$
\zeta
$$

表示：

$$
\boxed{\text{successive LATIN-indicator relative improvement}}
$$

。

对应关系：

```text
r_red
    → inner reduced solve quality

xi
    → nonlinear LATIN state convergence

zeta
    → basis adequacy / adaptive enrichment signal
```

。

---

# 71. 本阶段与上一阶段的关系

上一阶段完成：

$$
\boxed{\text{tower global stage state construction}}
$$

。

本阶段完成：

$$
\boxed{\text{tower global stage quality measurement and adaptive basis decision}}
$$

。

因此算法链现在是：

```text
local stage
    ↓
global stage
    ↓
relaxation
    ↓
Eq. (76) xi
    ↓
Eq. (60) zeta
    ↓
basis adequate?
    ↓
yes → next LATIN iteration
no  → enrichment
```

。

---

# 72. 现在为什么可以进入 Eq. (61)

进入 Eq. (61) 的 prerequisite 已经齐全：

1. 已知 current global stage 如何计算；
2. 已知 existing basis 如何做 temporal update；
3. 已知 reduced residual；
4. 已知 LATIN indicator；
5. 已知 saturation parameter；
6. 已知何时判断 basis insufficient；
7. 已知一维 mature control loop 如何触发 enrichment。

因此现在才真正具备进入：

$$
\boxed{\text{Eq. (61) new separated pair}}
$$

的完整上下文。

---

# 73. 下一阶段不应同时做什么

进入 Eq. (61) 后，仍应坚持：

- 不立即编码；
- 不同时推完 Eq. (61)–(72) 所有细节；
- 不同时改变 outer solver；
- 不重新设计 saturation control；
- 不引入 $n-\tau-x$；
- 不把 hardening / damage 加入 PGD unknown set。

下一步只解决：

> **原论文 Eq. (61) 中“新增一个 PGD separated pair”具体增加的是什么，以及这一 pair 在 tower fiber material-point space 中的变量对应和数据结构是什么。**

---

# 74. 对未来代码接口的建议

未来可将一维：

```text
relative_latin_indicator(
    local_state,
    global_state,
    directions,
    mesh,
    area,
    materials
)
```

逐步泛化为：

```text
relative_latin_indicator(
    local_state,
    global_state,
    directions,
    material_point_volumes,
    elastic_moduli
)
```

。

这样 convergence module 不需要了解 bar/tower mesh 拓扑，只需要：

$$
\boxed{\text{material-point states + positive integration weights}}
$$

。

---

# 75. Material-point-volume abstraction

建议未来统一：

```text
material_point_volumes : (N_q,)
```

满足：

$$
v_q>0
$$

。

同时用于：

```text
pgd_time_update
latin mechanical norm
reduced residual norm
future enrichment inner products
```

。

这可能成为 tower LATIN-PGD 中一个重要的通用 discretisation abstraction。

---

# 76. Tower norm regression tests

未来至少建立：

### Test 1

若所有 fields 为零：

$$
\|s^p\|=0
$$

。

### Test 2

若所有 fields 乘以常数 $a$，且 directions 固定：

$$
\|as^p\|=|a|\|s^p\|
$$

。

### Test 3

Difference state 为零：

$$
\xi=0
$$

。

### Test 4

所有 material-point volume weights 必须 positive。

### Test 5

总 material-point volume 与 tower geometric volume一致。

### Test 6

Bar special case 下 generalized norm 应退化到：

$$
AL_e
$$

权重形式。

---

# 77. Saturation regression tests

### Case A

$$
\xi_i=1,\qquad\xi_{i+1}=0.5
$$

应有：

$$
\zeta=\frac13
$$

。

### Case B

$$
\xi_i=\xi_{i+1}
$$

应有：

$$
\zeta=0
$$

。

### Case C

$$
\xi_{i+1}>\xi_i
$$

应有：

$$
\zeta<0
$$

并进入 enrichment-needed semantics。

### Case D

Both indicators numerical zero 时避免除零。

### Case E

NaN / inf indicators 明确拒绝。

---

# 78. Adaptive-control regression tests

未来 tower solver 应复用或移植以下行为测试：

1. empty basis forces first pair；
2. fixed baseline $\xi_i$ during same-stage enrichment；
3. $\xi\le\xi_{\mathrm{tol}}$ terminates before enrichment；
4. large $\zeta$ advances LATIN；
5. small $\zeta$ + large reduced residual enriches；
6. small $\zeta$ + already-small reduced residual accepts candidate；
7. negative $\zeta$ does not count as nonlinear convergence；
8. max enrichments only limits current global stage；
9. persistent stagnation requires consecutive accepted iterations；
10. max LATIN iterations terminates cleanly。

---

# 79. 本阶段不改代码的原因

虽然 Eq. (60) 与 tower norm 已经闭合，但本阶段仍不建议立刻改代码。

原因：

- Eq. (61)–(72) 的 enrichment data interface 尚未冻结；
- basis semantics 尚需与 new-mode addition 一起确认；
- temporal amplitude 是 total 还是 correction 的问题仍需在 enrichment 中统一；
- 如果现在先 generalized `iteration_control.py`，后面可能再次改 state/basis interface。

所以仍保持：

$$
\boxed{\text{theory first, code after enrichment interfaces freeze}}
$$

。

---

# 80. 本阶段最终结论汇总

1. Eq. (60) 公式本身无需为 tower 修改。
2. $\zeta$ 是 successive LATIN-indicator symmetric relative change。
3. $\zeta$ 不是 PGD truncation error。
4. $\zeta$ 不是 Eq. (59) reduced residual。
5. $\xi$ 是 local/global LATIN state relative distance。
6. $\zeta$ 通过 $\xi$ 的改善速度间接判断 basis adequacy。
7. $\zeta>0$ 表示 LATIN indicator 改善。
8. $\zeta=0$ 表示无进展。
9. $\zeta<0$ 表示 indicator 恶化，不能理解为成功收敛。
10. $\zeta=0.1$ 约对应 $\xi$ 至少下降 18.2%。
11. Eq. (76) 的 tower form保持不变。
12. 真正需要迁移的是 Eq. (77) mechanical norm。
13. Tower norm仍然是 material-point thermodynamic state norm。
14. 不应改成 nodal displacement norm。
15. 不应改成 section $N$–$M$ norm。
16. Eq. (77) 进入七个 mechanical fields。
17. $D,\dot D,Y$ 不显式进入 Eq. (77)。
18. Damage 通过 stress、elastic strain、search directions 间接进入 norm。
19. Fiber scalar norm density 已完整给出。
20. Tower spatial quadrature weight 为 $A_{egf}w_gJ_e$。
21. Eq. (59) 和 Eq. (77) 使用同一套 material-point volume weights。
22. 当前 320-point tower 可直接采用 $(N_t,320)$ norm arrays。
23. Time norm integration tower v1 先继承当前一维 `trapz`。
24. 同一 LATIN iteration 的 local/global/difference norm必须使用同一 search directions。
25. 当前一维三材料杆已经解决 adaptive control 的主要歧义。
26. $\zeta$ 主要负责 basis adequacy。
27. Nonlinear termination 不应只靠 small $\zeta$。
28. Absolute $\xi$ tolerance 是优先 stopping condition。
29. Persistent stagnation 是另一 nonlinear stopping condition。
30. Reduced residual 用于避免无意义 enrichment。
31. Same global stage 的 enrichment trials 始终使用固定 previous accepted $\xi_i$。
32. Empty basis 时先强制生成第一 mode。
33. Relaxed candidate 后再计算 $\xi$。
34. Tower v1 应继承 mature one-dimensional control architecture。
35. Tower 阶段只修改 norm discretisation，不重新发明 outer adaptive flow。
36. Eq. (60) 阶段已经闭合。
37. 下一步可以正式进入 Eq. (61) new PGD mode enrichment。

---

# 81. 一句话阶段定位

> **本阶段完成了原论文 Eq. (60) saturation criterion 与 Eq. (76)–(77) LATIN convergence norm 向 fiber beam-column offshore wind turbine tower 的完整迁移，并通过回查一维三材料杆 solver 确认 basis saturation 与 nonlinear convergence 的控制冲突此前已经被成熟地拆分为“$\zeta$ + reduced residual 的 PGD adaptive control”和“absolute $\xi$ + persistent stagnation 的 LATIN nonlinear control”；因此 Eq. (60) 阶段无需再重新设计控制逻辑，下一阶段可以直接进入 Eq. (61) new separated PGD pair。**
