# 海上风机塔筒 LATIN-PGD：固定空间基时间更新、损伤残差、硬化变量与完整 Global Stage 阶段总结

**日期：2026-08-15**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**当前分支：`feature/offshore-wind-turbine-tower-fatigue`**  
**前置总结：`docs/2026-08-13-tower-latin-pgd-eq47-53-spatial-derivation-summary.md`**  
**本阶段范围：自 Eq. (47)–(53) 空间问题总结之后，继续完成 Eq. (58)–(59) fixed spatial basis temporal update、search-direction fields、plastic forcing、plastic/damage 分裂、damage residual branch、hardening global update、damage-rate 与 energy-release-rate update，并形成完整 tower LATIN global-stage skeleton。**  
**下一阶段：原论文 Eq. (60) saturation criterion，然后进入 Eq. (61)–(72) new-mode enrichment。**

---

# 1. 本文档定位

上一份阶段总结已经完成原论文 Eq. (47)–(53) 从 bar/continuum 到 fiber beam-column offshore wind turbine tower 的空间问题迁移。

上一阶段建立了：

$$
\bar\varepsilon^p\longrightarrow\mathcal E_{\mathrm{tower}}\bar\varepsilon^p\longrightarrow C_0(\mathcal E_{\mathrm{tower}}-I)\bar\varepsilon^p
$$

并明确：

- tower PGD 的主空间未知量是 fiber-level plastic-strain spatial mode；
- compatible strain spatial mode 由 equilibrium projection 生成；
- stress spatial mode 不是独立未知量；
- reference elastic equilibrium operator 使用固定 $E_0$；
- 每个 stress spatial mode 在自由结构 DOFs 上满足离散弱平衡；
- 当前路线继续保持原论文的 $x-t$ 两变量 PGD，而不立即扩展为 $n-\tau-x$。

本阶段继续回答：

> **当 tower spatial modes 已经存在后，如何真正完成一个 LATIN global stage？**

---

# 2. 本阶段的核心问题

本阶段依次处理：

1. fixed spatial basis 下 temporal functions 如何更新；
2. tower fiber-level $H_\sigma$ 如何构造；
3. temporal reduced equation 的 forcing 如何定义；
4. 原论文 paper-separated forcing 与当前一维 repo coupled forcing 的差别；
5. tower v1 应采用哪一种 forcing；
6. damage residual branch 如何迁移；
7. hardening variables 如何更新；
8. $D,\dot D,Y$ 如何更新；
9. plastic、damage、hardening 三条 branch 如何重组成完整 global state；
10. 哪些现有一维代码可以继承，哪些暂不直接继承。

---

# 3. 当前总体策略再次确认

当前仍严格采用：

$$
\boxed{\text{original }x-t\text{ LATIN-PGD}\longrightarrow\text{fiber beam-column tower }x-t\text{ LATIN-PGD}}
$$

暂不做：

$$
x-t\longrightarrow n-\tau-x
$$

。

当前第一目标是：

> **确认原论文 LATIN-PGD 的核心 mathematical architecture 在 offshore wind turbine tower 的 fiber beam-column discretisation 上能够完整闭合。**

因此继续坚持：

- 不同时修改 PGD 理论与结构离散；
- 不把所有 constitutive variables 都 PGD 化；
- 不因为 offline SVD 显示低秩就立刻增加 online reduced variables；
- 先复现原论文 plastic/damage split；
- 先保持 diagonal descent search-direction approximation；
- 将 numerical safeguards 与 theory modifications 分开。

---

# 4. Tower material-point notation

连续空间意义下：

$$
x_{\mathrm{tower}}\equiv(s,y)
$$

有限元离散后：

$$
(s,y)\longrightarrow(e,g,f)
$$

其中：

- $e$：beam element；
- $g$：beam Gauss point；
- $f$：section fiber。

为了 reduced linear algebra，继续 flatten：

$$
q=q(e,g,f)
$$

并定义：

$$
q=1,\ldots,N_q
$$

其中：

$$
N_q=N_eN_gN_f
$$

当前 coarse tower model：

$$
N_e=10,\qquad N_g=2,\qquad N_f=16
$$

所以：

$$
\boxed{N_q=320}
$$

。

---

# 5. 前置 equilibrium projection

定义 fiber kinematic matrix：

$$
\varepsilon=HU_F
$$

定义 quadrature volume matrix：

$$
M=\operatorname{diag}(A_{egf}w_gJ_e)
$$

定义 reference elastic operator：

$$
C_0=E_0I
$$

则：

$$
K^0_{FF}=H^TMC_0H
$$

对于任意 source strain field $q$：

$$
\boxed{\mathcal E_{\mathrm{tower}}=H(H^TMC_0H)^{-1}H^TMC_0}
$$

compatible strain：

$$
\boxed{\varepsilon^{\mathrm{comp}}=\mathcal E_{\mathrm{tower}}q}
$$

equilibrated stress：

$$
\boxed{\sigma^{\mathrm{eq}}=C_0(\mathcal E_{\mathrm{tower}}-I)q}
$$

。

本阶段进一步确认：这个 operator 不仅用于 plastic source strain，也用于 damage residual source strain。

---

# 6. Eq. (58)–(59) 的任务

原论文 hybrid PGD strategy 在每一个 LATIN global stage 中首先：

> **保持已有 spatial PGD basis 不变，只重新更新 temporal functions。**

只有当更新 temporal functions 后 reduced mechanical residual 仍然过大，才进入 new-mode enrichment。

因此 Eq. (58)–(59) 的任务是：

$$
\boxed{\text{fixed spatial basis}\longrightarrow\text{update temporal coefficients}}
$$

。

---

# 7. Fixed spatial basis 的 tower matrix form

已有 $m$ 个 plastic-strain spatial modes：

$$
\bar\varepsilon^p_1,\ldots,\bar\varepsilon^p_m
$$

定义：

$$
\boxed{P=\begin{bmatrix}\bar\varepsilon^p_1&\cdots&\bar\varepsilon^p_m\end{bmatrix}}
$$

其中：

$$
P\in\mathbb R^{N_q\times m}
$$

。

每个 plastic mode 对应：

$$
\bar\sigma_j=C_0(\mathcal E_{\mathrm{tower}}-I)\bar\varepsilon^p_j
$$

定义：

$$
\boxed{S=\begin{bmatrix}\bar\sigma_1&\cdots&\bar\sigma_m\end{bmatrix}}
$$

其中：

$$
S\in\mathbb R^{N_q\times m}
$$

。

---

# 8. Temporal coefficient vector

定义：

$$
\Delta\lambda(t)=\begin{bmatrix}\Delta\lambda_1(t)&\cdots&\Delta\lambda_m(t)\end{bmatrix}^T
$$

以及：

$$
\Delta\dot\lambda(t)=\begin{bmatrix}\Delta\dot\lambda_1(t)&\cdots&\Delta\dot\lambda_m(t)\end{bmatrix}^T
$$

。

Tower Eq. (58)：

$$
\boxed{\Delta\dot\varepsilon^p(t)=P\Delta\dot\lambda(t)}
$$

$$
\boxed{\Delta\sigma'(t)=S\Delta\lambda(t)}
$$

。

此时：

- $P$ 固定；
- $S$ 固定；
- reduced unknown 只有 $\Delta\lambda(t)$。

---

# 9. 为什么 temporal update 不需要重新求 tower equilibrium

每个 stress mode 已满足：

$$
H^TM\bar\sigma_j=0
$$

所以：

$$
\Delta\sigma'(t)=\sum_{j=1}^{m}\Delta\lambda_j(t)\bar\sigma_j
$$

仍满足：

$$
\boxed{H^TM\Delta\sigma'(t)=0}
$$

。

因此 temporal coefficients 改变时，不需要：

- 重新构造 compatible spatial mode；
- 重新求 reference tower equilibrium；
- 重新做 tower nonlinear Newton。

---

# 10. Tower temporal residual

定义：

$$
D_H(t)=\operatorname{diag}\left(H_{\sigma,1}(t),\ldots,H_{\sigma,N_q}(t)\right)
$$

则：

$$
\boxed{r(t)=P\Delta\dot\lambda(t)-D_H(t)S\Delta\lambda(t)-f(t)}
$$

。

这里：

- $H_\sigma(x,t)$：known search-direction field；
- $f(x,t)$：known forcing field；
- $\Delta\lambda(t)$：reduced unknown。

---

# 11. Eq. (59) 的 tower weighted norm

定义每个 fiber material point 的 quadrature volume：

$$
\boxed{v_{egf}=A_{egf}w_gJ_e}
$$

flatten 后：

$$
V=\operatorname{diag}(v_1,\ldots,v_{N_q})
$$

定义：

$$
\boxed{W_H(t)=VD_H(t)^{-1}}
$$

则 temporal problem：

$$
\boxed{\Delta\lambda=\arg\min\int_0^T r(t)^TW_H(t)r(t)\,dt}
$$

。

这不是普通 Euclidean fitting，而是 LATIN mechanical-search-direction weighted projection。

---

# 12. 为什么 tower quadrature weights 必须进入 temporal norm

Tower 是变截面结构，不同 fiber material points 对总体 mechanical residual 的权重不同。

影响包括：

- fiber area；
- beam Gauss weight；
- element Jacobian；
- tower taper；
- element length；
- section geometry。

因此不能只使用：

$$
\sum_qr_q^2
$$

而应使用：

$$
\sum_qv_q\frac{r_q^2}{H_{\sigma,q}}
$$

。

---

# 13. 当前数据尺寸

当前：

$$
N_q=320
$$

若已有 $m$ 个 modes：

| Quantity | Shape |
|---|---:|
| $P$ | $320\times m$ |
| $S$ | $320\times m$ |
| $H_\sigma$ | $N_t\times320$ |
| $f$ | $N_t\times320$ |
| $\Delta\lambda$ | $N_t\times m$ |
| reduced matrix | $m\times m$ |

例如：

$$
m=5
$$

时每个时间点只需求解：

$$
5\times5
$$

规模的问题。

---

# 14. 时间离散

当前一维 `pgd_time_update.py` 使用：

$$
\boxed{\Delta\dot\lambda_n\approx\frac{\Delta\lambda_n-\Delta\lambda_{n-1}}{\Delta t_n}}
$$

。

于是：

$$
r_n=P\frac{\Delta\lambda_n-\Delta\lambda_{n-1}}{\Delta t_n}-D_{H,n}S\Delta\lambda_n-f_n
$$

整理为：

$$
\boxed{r_n=A_n\Delta\lambda_n-b_n}
$$

其中：

$$
\boxed{A_n=\frac{P}{\Delta t_n}-D_{H,n}S}
$$

$$
\boxed{b_n=f_n+\frac{P}{\Delta t_n}\Delta\lambda_{n-1}}
$$

。

---

# 15. 每一时间点的 weighted least squares

$$
\boxed{\Delta\lambda_n=\arg\min\left\|A_n\Delta\lambda_n-b_n\right\|_{W_{H,n}}^2}
$$

normal equation：

$$
\boxed{A_n^TW_{H,n}A_n\Delta\lambda_n=A_n^TW_{H,n}b_n}
$$

。

实际实现仍建议使用 weighted least-squares / QR / SVD，而不是显式 normal equation。

当前一维代码使用 `numpy.linalg.lstsq` 并保存 condition-history，这一策略可继承。

---

# 16. Reduced matrix 展开

$$
\boxed{K_{\mathrm{red},n}=\frac{1}{\Delta t_n^2}P^TVD_{H,n}^{-1}P-\frac{1}{\Delta t_n}(P^TVS+S^TVP)+S^TVD_{H,n}S}
$$

。

三部分分别来自：

1. plastic-strain-rate contribution；
2. rate-stress coupling；
3. search-direction stress contribution。

---

# 17. $\mathcal E_{\mathrm{tower}}$ 与 $H_\sigma$ 的角色区别

必须保持：

$$
\boxed{\mathcal E_{\mathrm{tower}}\text{ fixed}}
$$

因为：

$$
\bar\sigma_j=C_0(\mathcal E_{\mathrm{tower}}-I)\bar\varepsilon^p_j
$$

必须是纯 spatial mode。

但是：

$$
\boxed{H_\sigma(x,t)\text{ may vary}}
$$

因为它只是当前 LATIN iteration 中的 known weighting / descent coefficient field。

两者不矛盾。

---

# 18. Fiber-level $H_\sigma$

定义 effective relative stress：

$$
\boxed{\tau=\frac{\sigma}{1-D}-\beta}
$$

yield function：

$$
\boxed{f^p=|\tau|+\frac{a}{2C}\beta^2-R-\sigma_y}
$$

Norton multiplier：

$$
\boxed{\dot p=k\langle f^p\rangle_+^n}
$$

其中：

$$
k=K^{-n}
$$

。

Plastic strain rate：

$$
\boxed{\dot\varepsilon^p=k\langle f^p\rangle_+^n\frac{\operatorname{sign}(\tau)}{1-D}}
$$

。

---

# 19. $H_\sigma$ 的 raw tangent

在 smooth plastic region：

$$
f^p>0
$$

且：

$$
\tau\neq0
$$

有：

$$
\frac{\partial\tau}{\partial\sigma}=\frac{1}{1-D}
$$

所以：

$$
\boxed{H_\sigma^{\mathrm{raw}}=kn\frac{\langle f^p\rangle_+^{n-1}}{(1-D)^2}}
$$

。

---

# 20. $H_\sigma$ regularization

Elastic regime 下：

$$
f^p\le0
$$

会导致：

$$
H_\sigma^{\mathrm{raw}}=0
$$

但 Eq. (59) 需要 $H_\sigma^{-1}$。

因此：

$$
\boxed{H_\sigma=H_\sigma^{\mathrm{raw}}+\frac{\zeta_{\mathrm{reg}}}{E_0}}
$$

当前：

$$
\zeta_{\mathrm{reg}}=0.15
$$

。

所以未来 tower implementation 必须保证：

$$
\boxed{H_{\sigma,q,n}>0}
$$

。

---

# 21. $H_\sigma$ 的 tower data flow

对每一个 local-stage fiber state：

```text
hat_sigma
hat_beta
hat_R_bar
hat_D
```

先恢复 physical isotropic hardening force：

$$
R=R(\bar R)
$$

然后：

$$
\tau=\frac{\hat\sigma}{1-\hat D}-\hat\beta
$$

再计算：

$$
f^p
$$

进而得到：

$$
H_\sigma^{\mathrm{raw}}
$$

最后：

$$
H_\sigma=H_\sigma^{\mathrm{raw}}+\frac{\zeta_{\mathrm{reg}}}{E_0}
$$

。

因此：

$$
\boxed{H_\sigma\in\mathbb R^{N_t\times320}}
$$

。

量纲：

$$
\boxed{[H_\sigma]=(\text{stress}\cdot T)^{-1}}
$$

所以：

$$
H_\sigma\Delta\sigma
$$

与：

$$
\Delta\dot\varepsilon^p
$$

量纲一致。

---

# 22. Plastic forcing：从 total-stress search relation 出发

原论文 total-stress form：

$$
\boxed{\Delta\dot\varepsilon^p-H_\sigma\Delta\sigma+\Delta\bar\varepsilon=0}
$$

其中：

$$
\boxed{\Delta\bar\varepsilon=H_\sigma(\hat\sigma-\sigma_i)-(\hat{\dot\varepsilon}^p-\dot\varepsilon_i^p)}
$$

。

由于：

$$
\Delta\sigma=\Delta\sigma'+\Delta\tilde\sigma
$$

必须明确 plastic/damage split 之后 temporal reduced equation 中保留哪些项。

---

# 23. Paper-separated forcing

原论文在 plastic reduced equation 中只保留：

$$
\Delta\sigma'
$$

于是：

$$
\boxed{\Delta\dot\varepsilon^p-H_\sigma\Delta\sigma'+\Delta\bar\varepsilon\approx0}
$$

因此：

$$
\boxed{f^{\mathrm{sep}}=-\Delta\bar\varepsilon}
$$

也就是：

$$
\boxed{f^{\mathrm{sep}}=\hat{\dot\varepsilon}^p-\dot\varepsilon_i^p-H_\sigma(\hat\sigma-\sigma_i)}
$$

。

这里以后统一称为：

$$
\boxed{\text{paper-separated formulation}}
$$

。

---

# 24. 当前 ascent choice 下 forcing 的简化

当前 local stage 采用：

$$
(B^+)^{-1}=0
$$

因此：

$$
\boxed{\hat\sigma=\sigma_i}
$$

于是：

$$
\boxed{f^{\mathrm{sep}}=\hat{\dot\varepsilon}^p-\dot\varepsilon_i^p}
$$

。

这是 tower v1 当前正式采用的 plastic forcing。

---

# 25. Current repo coupled forcing

当前一维 `latin/global_stage.py` 与 `latin/pgd_global_stage.py` 实际使用：

$$
\boxed{f^{\mathrm{coupled}}=f^{\mathrm{sep}}+H_\sigma\Delta\tilde\sigma}
$$

。

这一形式可以从：

$$
\Delta\dot\varepsilon^p-H_\sigma(\Delta\sigma'+\Delta\tilde\sigma)+\Delta\bar\varepsilon=0
$$

直接移项得到。

因此它并不是“错误公式”，而是保留了 total-stress search relation 中 damage stress coupling 的 formulation。

---

# 26. 两种 forcing 的正式关系

$$
\boxed{f^{\mathrm{coupled}}-f^{\mathrm{sep}}=H_\sigma\Delta\tilde\sigma}
$$

。

因此二者差别完全由 damage-dependent stress correction 决定。

对于 low-damage metallic fatigue：

$$
\Delta\tilde\sigma
$$

通常应相对较小，因此原论文采用 plastic/damage separation 具有物理基础。

---

# 27. Tower v1 的 formulation 决策

本阶段正式决定：

$$
\boxed{\text{Tower v1 uses paper-separated forcing}}
$$

也就是：

$$
\boxed{f^{\mathrm{tower,v1}}=f^{\mathrm{sep}}}
$$

当前 ascent choice：

$$
\boxed{f^{\mathrm{tower,v1}}=\hat{\dot\varepsilon}^p-\dot\varepsilon_i^p}
$$

。

主要理由：

1. 当前目标是先迁移 original method；
2. 不希望同时改变 plastic/damage coupling；
3. 若 tower 结果出现差异，需要能够区分 structural migration 与 formulation modification；
4. paper-separated formulation 与 Eq. (41) → Eq. (53) → Eq. (58)–(59) 的原论文理论链更一致。

---

# 28. 对当前一维 repo 的处理

当前一维 reproduction 已经验证，并采用 coupled forcing。

因此本阶段决定：

> **不回改现有一维代码。**

而是：

- 明确记录 formulation difference；
- tower v1 使用 paper-separated forcing；
- 后续做 A/B sensitivity comparison。

未来比较：

```text
Variant A:
    paper-separated forcing

Variant B:
    coupled-total-stress forcing
```

观察：

- tower-top displacement；
- critical fiber stress；
- plastic strain；
- damage；
- LATIN iterations；
- PGD mode count；
- reduced residual；
- runtime。

---

# 29. Damage residual branch：符号约定

原论文 residual stress 使用 $\Delta R$，但当前 material model 已有：

$$
R
$$

作为 isotropic hardening force。

为避免混淆，本阶段统一把 residual stress 写成：

$$
\boxed{\Delta R^{\mathrm{res}}}
$$

程序命名建议：

```text
residual_stress
```

而 hardening 继续保留：

```text
R
R_bar
r_bar
```

。

---

# 30. Residual stress

对 tower fiber material point $q=(e,g,f)$：

$$
\boxed{\Delta R_q^{\mathrm{res}}=(\sigma_{i,q}-\hat\sigma_q)-E_0(\varepsilon_{i,q}^e-\hat\varepsilon_q^e)}
$$

。

其物理意义：

> **当前 global state 与 local state 之间，reference Hooke law 无法解释的 damaged-elastic mismatch。**

它不是新的 constitutive internal variable。

---

# 31. 无 damage 时 residual branch 自动消失

若：

$$
\sigma=E_0\varepsilon^e
$$

则：

$$
\sigma_i-\hat\sigma=E_0(\varepsilon_i^e-\hat\varepsilon^e)
$$

所以：

$$
\boxed{\Delta R^{\mathrm{res}}=0}
$$

进而：

$$
\boxed{\Delta\varepsilon^R=0}
$$

以及：

$$
\boxed{\Delta\tilde\sigma=0}
$$

。

因此该 branch 确实专门表示 nonlinear damaged elasticity。

---

# 32. 当前 ascent choice 下 residual stress 的简化

因为：

$$
\hat\sigma=\sigma_i
$$

所以：

$$
\boxed{\Delta R^{\mathrm{res}}=E_0(\hat\varepsilon^e-\varepsilon_i^e)}
$$

。

这说明即使 stress 不变，只要：

$$
\hat D\neq D_i
$$

damaged compliance 改变，就会有：

$$
\hat\varepsilon^e\neq\varepsilon_i^e
$$

因此 damage residual correction 仍然存在。

---

# 33. Residual strain

定义：

$$
\boxed{\Delta\varepsilon^R=C_0^{-1}\Delta R^{\mathrm{res}}}
$$

在 scalar fiber material 中：

$$
\boxed{\Delta\varepsilon_q^R=\frac{\Delta R_q^{\mathrm{res}}}{E_0}}
$$

。

当前 ascent choice：

$$
\boxed{\Delta\varepsilon^R=\hat\varepsilon^e-\varepsilon_i^e}
$$

。

需要强调：

> $\Delta\varepsilon^R$ 不是 damage increment，而是 known equivalent residual source strain field。

---

# 34. Damage branch 复用同一个 $\mathcal E_{\mathrm{tower}}$

Damage-dependent correction：

$$
\boxed{\Delta\tilde\sigma=C_0(\Delta\tilde\varepsilon-\Delta\varepsilon^R)}
$$

满足 static admissibility。

因此：

$$
\boxed{\Delta\tilde\varepsilon=\mathcal E_{\mathrm{tower}}\Delta\varepsilon^R}
$$

$$
\boxed{\Delta\tilde\sigma=C_0(\mathcal E_{\mathrm{tower}}-I)\Delta\varepsilon^R}
$$

。

所以不需要两套 spatial operators。

建议只保留：

```text
tower_equilibrium_operator(source_strain)
```

统一处理 plastic source strain 与 damage residual source strain。

---

# 35. Damage branch 的 fiber → section → element → tower 离散

定义：

$$
a_{egf}=\begin{bmatrix}1\\-y_{egf}\end{bmatrix}
$$

等效 section resultant：

$$
\boxed{\Delta r^R_{eg}=E_0\sum_fA_{egf}a_{egf}\Delta\varepsilon^R_{egf}}
$$

其中：

$$
\Delta r^R_{eg}=\begin{bmatrix}\Delta N^R_{eg}\\\Delta M^R_{eg}\end{bmatrix}
$$

。

Element equivalent load：

$$
\boxed{\Delta f_e^{R,\ell}=\sum_gw_gJ_eB_{eg}^T\Delta r^R_{eg}}
$$

坐标转换：

$$
\boxed{\Delta f_e^R=T_e^T\Delta f_e^{R,\ell}}
$$

装配：

$$
\boxed{\Delta F^R=\operatorname{Assembly}_e(\Delta f_e^R)}
$$

。

---

# 36. Damage branch 的 reference elastic solve

自由 DOFs：

$$
\boxed{K^0_{FF}\Delta\tilde U_F=\Delta F_F^R}
$$

然后：

$$
\boxed{\Delta\tilde\varepsilon=H\Delta\tilde U_F}
$$

。

由于：

$$
K^0_{FF}
$$

与 plastic spatial projection 完全相同，可以：

```text
assemble K0 once
factorise K0_FF once
reuse for many RHS
```

。

---

# 37. Damage correction 的 equilibrium property

由：

$$
\Delta\tilde\varepsilon=\mathcal E_{\mathrm{tower}}\Delta\varepsilon^R
$$

得到：

$$
H^TMC_0(\Delta\tilde\varepsilon-\Delta\varepsilon^R)=0
$$

因此：

$$
\boxed{H^TM\Delta\tilde\sigma=0}
$$

。

Plastic branch 同样有：

$$
H^TM\Delta\sigma'=0
$$

所以：

$$
\boxed{H^TM(\Delta\sigma'+\Delta\tilde\sigma)=0}
$$

。

---

# 38. Mechanical field recombination

Plastic strain：

$$
\boxed{\varepsilon_{i+1}^p=\varepsilon_i^p+\Delta\varepsilon^p}
$$

Stress：

$$
\boxed{\sigma_{i+1}=\sigma_i+\Delta\sigma'+\Delta\tilde\sigma}
$$

Total strain correction：

$$
\boxed{\Delta\varepsilon=\Delta\varepsilon'+\Delta\tilde\varepsilon}
$$

Elastic strain：

$$
\boxed{\varepsilon_{i+1}^e=\varepsilon_i^e+\Delta\varepsilon'-\Delta\varepsilon^p+\Delta\tilde\varepsilon}
$$

其中：

$$
\boxed{\Delta\varepsilon'=\mathcal E_{\mathrm{tower}}\Delta\varepsilon^p}
$$

。

---

# 39. Damage branch 不做 PGD

Tower v1 明确：

$$
\boxed{\text{damage residual branch is not PGD-reduced}}
$$

。

理由：

- $\Delta\varepsilon^R(x,t)$ 已知；
- 它不是待求的新 separated field；
- 只需 fixed reference operator 的 linear projection；
- 保持原论文 formulation；
- computationally 只是 fixed matrix multiple-RHS solve。

所以：

$$
\boxed{\text{plastic branch reduced, damage branch full projection}}
$$

。

---

# 40. Hardening variables

定义：

$$
X=(\alpha,\bar r)
$$

与：

$$
Z=(\beta,\bar R)
$$

。

Global state laws：

$$
\boxed{\beta=C\alpha}
$$

以及：

$$
\boxed{\bar R=R_\infty\bar r}
$$

。

Partial-normal transformation 的重要作用之一，就是将 isotropic hardening relation 转换为这一 linear state law。

---

# 41. Hardening descent equations

Diagonal approximation 下：

$$
H_Z=\operatorname{diag}(H_\beta,H_{\bar R})
$$

Kinematic hardening：

$$
\boxed{\dot\alpha-\hat{\dot\alpha}+H_\beta(\beta-\hat\beta)=0}
$$

所以：

$$
\boxed{\dot\alpha+H_\beta C\alpha=\hat{\dot\alpha}+H_\beta\hat\beta}
$$

。

Isotropic hardening：

$$
\boxed{\dot{\bar r}-\hat{\dot{\bar r}}+H_{\bar R}(\bar R-\hat{\bar R})=0}
$$

所以：

$$
\boxed{\dot{\bar r}+H_{\bar R}R_\infty\bar r=\hat{\dot{\bar r}}+H_{\bar R}\hat{\bar R}}
$$

。

---

# 42. Hardening 的 backward-Euler update

Kinematic hardening：

$$
\boxed{\alpha_n=\frac{\alpha_{n-1}+\Delta t_n(\hat{\dot\alpha}_n+H_{\beta,n}\hat\beta_n)}{1+\Delta t_nH_{\beta,n}C}}
$$

$$
\boxed{\dot\alpha_n=\frac{\alpha_n-\alpha_{n-1}}{\Delta t_n}}
$$

$$
\boxed{\beta_n=C\alpha_n}
$$

。

Isotropic hardening：

$$
\boxed{\bar r_n=\frac{\bar r_{n-1}+\Delta t_n(\hat{\dot{\bar r}}_n+H_{\bar R,n}\hat{\bar R}_n)}{1+\Delta t_nH_{\bar R,n}R_\infty}}
$$

$$
\boxed{\dot{\bar r}_n=\frac{\bar r_n-\bar r_{n-1}}{\Delta t_n}}
$$

$$
\boxed{\bar R_n=R_\infty\bar r_n}
$$

。

---

# 43. 为什么 hardening 不需要 tower equilibrium

这些变量：

$$
\alpha,\beta,\bar r,\bar R
$$

是 fiber-local internal variables。

它们不需要：

$$
B_{eg},\quad T_e,\quad K^0,\quad\mathcal E_{\mathrm{tower}}
$$

。

Tower geometry 只通过 local stress history 间接影响 hardening evolution。

---

# 44. $H_\beta$

定义：

$$
s=\operatorname{sign}(\tau)
$$

以及：

$$
c_t=kn\langle f^p\rangle_+^{n-1}
$$

。

Yield function 对 $\beta$ 的导数：

$$
\boxed{g_\beta=-s+\frac{a}{C}\beta}
$$

。

Raw diagonal Hessian：

$$
\boxed{H_\beta^{\mathrm{raw}}=c_tg_\beta^2+\dot p\frac{a}{C}}
$$

加入 regularization：

$$
\boxed{H_\beta=c_t\left(-s+\frac{a\beta}{C}\right)^2+\dot p\frac{a}{C}+\frac{\zeta_{\mathrm{reg}}}{C}}
$$

。

---

# 45. $H_{\bar R}$

定义：

$$
q_R=\frac{\sqrt{\gamma}\bar R}{2R_\infty}
$$

以及：

$$
\chi=1-\frac{\sqrt{\gamma}\bar R}{2R_\infty}
$$

。

Raw diagonal Hessian：

$$
\boxed{H_{\bar R}^{\mathrm{raw}}=c_t\gamma\chi^2+\dot p\frac{\gamma}{2R_\infty}}
$$

加入 regularization：

$$
\boxed{H_{\bar R}=c_t\gamma\left(1-\frac{\sqrt{\gamma}\bar R}{2R_\infty}\right)^2+\dot p\frac{\gamma}{2R_\infty}+\frac{\zeta_{\mathrm{reg}}}{R_\infty}}
$$

。

---

# 46. Hardening 不做 PGD 的理由

即使 offline SVD 显示：

$$
\alpha(x,t)
$$

或：

$$
\bar r(x,t)
$$

具有低秩结构，也不意味着 online solver 必须 reduced。

第一版判断标准是：

$$
\boxed{\text{computational bottleneck}}
$$

。

Hardening update 只是：

$$
O(N_tN_q)
$$

个 scalar local operations，不涉及 global structural solve。

因此 tower v1 不对 hardening fields 单独做 PGD。

---

# 47. Same-stage 解耦与 LATIN-iteration-level 耦合

由于采用 diagonal search direction：

$$
H_{\sigma\beta}=H_{\sigma\bar R}=H_{\beta\bar R}=0
$$

所以同一 global stage 内：

```text
plastic PGD temporal solve
```

与：

```text
hardening local linear ODEs
```

解耦。

但当前 iteration 得到的：

$$
\beta_{i+1},\bar R_{i+1}
$$

会进入下一轮 local stage，重新影响：

$$
f^p,\dot\varepsilon^p,H_\sigma
$$

。

因此它们是：

$$
\boxed{\text{decoupled within a global stage, coupled across LATIN iterations}}
$$

。

---

# 48. Damage evolution law

对每一个 fiber：

$$
\boxed{\dot D=k_d\langle Y-Y_0\rangle_+^{n_d}}
$$

其中：

$$
\boxed{Y_0=\frac{\sigma_y^2}{2E_0}}
$$

。

Damage evolution 仍然是 fiber-local material history。

---

# 49. Global-stage damage search direction

原论文采用：

$$
\boxed{b^-=0}
$$

所以：

$$
\boxed{\dot D_{i+1}=\hat{\dot D}_{i+1/2}}
$$

。

准确理解：

> global stage 不再沿 damage evolution direction 对 local-stage damage-rate history 做额外 correction。

不能误解为：

$$
D_{i+1}=D_i
$$

。

---

# 50. Damage history 的恢复

连续形式：

$$
\boxed{D_{i+1}(t)=D_0+\int_0^t\dot D_{i+1}(\tau)\,d\tau}
$$

当前一维 global-stage implementation 使用：

$$
\boxed{D_n=D_{n-1}+\Delta t_n\dot D_n}
$$

并施加：

$$
0\le D\le D_{\max}<1
$$

。

---

# 51. Damage 时间离散一致性问题

当前 local stage 使用 RK4 积分 internal variables。

Global stage 当前一维代码则：

1. 复制 $\hat{\dot D}$；
2. 再以 backward-Euler style 重建 $D$。

因此可能出现：

$$
D_{i+1}(t_n)\neq\hat D_{i+1/2}(t_n)
$$

的离散差异。

Tower v1 编码前需专门决定：

```text
Option A:
copy local damage_rate and re-integrate D

Option B:
directly inherit local integrated damage history
```

理论核心仍是：

$$
b^-=0
$$

即 global stage 不修改 local damage evolution。

---

# 52. Energy-release rate $Y$

Global state 必须满足：

$$
Y=Y(\sigma,D)
$$

。

当前 1D fiber unilateral damage specialization：

当：

$$
\sigma\ge0
$$

有：

$$
\boxed{Y=\frac{\sigma^2}{2E_0(1-D)^2}}
$$

当：

$$
\sigma<0
$$

有：

$$
\boxed{Y=\frac{h\sigma^2}{2E_0(1-hD)^2}}
$$

。

其中：

$$
0<h<1
$$

描述 compression closure / unilateral effect。

---

# 53. 为什么 $Y$ 必须使用最终 global stress 重算

Local stage 的：

$$
\hat Y
$$

对应：

$$
\hat\sigma,\hat D
$$

。

而 global mechanical correction 完成后：

$$
\sigma_{i+1}=\sigma_i+\Delta\sigma'+\Delta\tilde\sigma
$$

一般：

$$
\sigma_{i+1}\neq\hat\sigma
$$

。

因此不能简单设置：

$$
Y_{i+1}=\hat Y
$$

而必须：

$$
\boxed{Y_{i+1}=Y(\sigma_{i+1},D_{i+1})}
$$

。

---

# 54. 完整 tower LATIN global state

原论文 LATIN state 可写：

$$
s=\{\dot\varepsilon^p,\varepsilon^e,\dot X,\dot D,\sigma,Z,Y\}
$$

其中：

$$
X=(\alpha,\bar r)
$$

$$
Z=(\beta,\bar R)
$$

。

Tower fiber-level：

$$
\boxed{s_q(t)=\{\dot\varepsilon_q^p,\varepsilon_q^e,\dot\alpha_q,\dot{\bar r}_q,\dot D_q,\sigma_q,\beta_q,\bar R_q,Y_q\}}
$$

并为数值积分附加保存：

$$
\varepsilon_q^p,\alpha_q,\bar r_q,D_q
$$

。

---

# 55. 建议的 tower state arrays

当前：

$$
N_q=320
$$

建议：

```text
plastic_strain         : (N_t, 320)
plastic_strain_rate    : (N_t, 320)

elastic_strain         : (N_t, 320)
stress                 : (N_t, 320)

alpha                  : (N_t, 320)
alpha_rate             : (N_t, 320)
beta                   : (N_t, 320)

r_bar                  : (N_t, 320)
r_bar_rate             : (N_t, 320)
R_bar                  : (N_t, 320)

damage                 : (N_t, 320)
damage_rate            : (N_t, 320)
energy_release_rate    : (N_t, 320)
```

。

---

# 56. 完整 global-stage 三条 branch

| Branch | 核心变量 | 求解方式 |
|---|---|---|
| Structural plastic branch | $\Delta\varepsilon^p,\Delta\sigma'$ | LATIN-PGD reduced temporal solve |
| Damage residual mechanical branch | $\Delta\varepsilon^R,\Delta\tilde\sigma$ | full reference-equilibrium projection |
| Hardening/damage local branch | $\alpha,\beta,\bar r,\bar R,D,\dot D,Y$ | fiber-local ODE / algebraic state update |

这是当前 tower v1 最重要的结构结论之一。

---

# 57. 完整 global-stage algorithm skeleton

```text
INPUT:
    previous global state s_i
    current local state hat_s_(i+1/2)
    current PGD spatial basis
    reference tower operator

STEP A:
    compute search directions
        H_sigma
        H_beta
        H_R_bar
        b_damage = 0

STEP B:
    plastic PGD branch
        build paper-separated forcing
        update existing temporal functions
        obtain Delta_eps_p
        obtain Delta_sigma'
        obtain Delta_eps'

STEP C:
    damage residual branch
        compute residual stress
        compute Delta_eps_R
        apply E_tower
        obtain Delta_eps_tilde
        obtain Delta_sigma_tilde

STEP D:
    recombine mechanical fields
        sigma_(i+1)
        eps_e_(i+1)
        eps_p_(i+1)

STEP E:
    hardening global update
        alpha
        beta
        r_bar
        R_bar

STEP F:
    damage update
        dot_D_(i+1) = hat_dot_D
        recover D_(i+1)

STEP G:
    energy-release-rate state law
        Y_(i+1) = Y(sigma_(i+1),D_(i+1))

OUTPUT:
    new global state s_(i+1) in A
```

---

# 58. 与当前一维 repo 的继承关系

可直接继承的核心思想来自：

```text
latin/pgd_time_update.py
latin/search_directions.py
latin/global_stage.py
latin/pgd_global_stage.py
latin/local_stage.py
latin/equilibrium_operator.py
latin/pgd_basis.py
material/viscoplastic_damage_1d.py
```

可以继承：

- fixed-basis temporal update；
- weighted least squares；
- condition-number diagnostics；
- diagonal $H_\sigma,H_\beta,H_{\bar R}$；
- regularization；
- hardening backward-Euler update；
- damage residual projection；
- stress/elastic-strain recombination；
- source-strain equilibrium-operator abstraction。

---

# 59. 当前一维 repo 中不直接继承到 tower v1 的部分

最重要的是：

$$
\boxed{+H_\sigma\Delta\tilde\sigma}
$$

进入 plastic forcing 的 coupled formulation。

Current 1D：

$$
f^{\mathrm{coupled}}=f^{\mathrm{sep}}+H_\sigma\Delta\tilde\sigma
$$

Tower v1：

$$
\boxed{f^{\mathrm{tower,v1}}=f^{\mathrm{sep}}}
$$

。

这一差异未来代码必须显式标记。

---

# 60. Current repo 中应继续继承的 damage 处理

虽然 forcing coupling 暂不继承，但仍应继承：

1. residual strain calculation；
2. common equilibrium projection；
3. damage stress correction；
4. final stress recombination；
5. final elastic-strain recombination；
6. final $Y(\sigma,D)$ post-processing。

所以：

> **不把 damage stress 回灌到 plastic forcing，不等于忽略 damage correction。**

---

# 61. Future `tower_equilibrium_operator.py`

建议未来建立独立模块：

```text
latin/tower_equilibrium_operator.py
```

输入：

```text
source_strain[q]
or
source_strain[t,q]
```

内部：

```text
fiber source strain
    ↓
section equivalent resultants
    ↓
element equivalent loads
    ↓
global RHS
    ↓
fixed K0_FF solve
    ↓
compatible fiber strain
    ↓
equilibrated fiber stress
```

输出建议：

```text
source_strain
compatible_strain
stress
displacement
base_reaction
```

。

---

# 62. Reference stiffness 实现原则

Tower v1 应：

1. 用 $E_0$ 构造 reference section stiffness；
2. assembly $K^0$；
3. apply tower-base constraints；
4. 得到 $K^0_{FF}$；
5. factorize once；
6. reuse for all source-strain projections。

必须避免：

```text
current damaged algorithmic tangent
```

直接替代：

```text
reference elastic K0
```

。

---

# 63. Future search-direction arrays

建议：

```text
H_sigma   : (N_t, N_q)
H_beta    : (N_t, N_q)
H_R_bar   : (N_t, N_q)
b_damage  : (N_t, N_q)
```

其中：

$$
\boxed{b_{\mathrm{damage}}=0}
$$

。

所有 search directions 都由 current local-stage state 计算。

---

# 64. Future PGD basis arrays

单个 mode：

```text
spatial_plastic_strain : (N_q,)
spatial_stress         : (N_q,)
temporal_amplitude     : (N_t,)
temporal_rate          : (N_t,)
```

多个 modes：

```text
P          : (N_q, m)
S          : (N_q, m)
Lambda     : (N_t, m)
Lambda_dot : (N_t, m)
```

。

---

# 65. Physics index 与 reduced index

理论中继续写：

$$
(e,g,f)
$$

以保留物理意义。

Reduced solver 中使用：

$$
q
$$

更方便。

因此：

$$
\boxed{\text{physics uses }(e,g,f),\qquad\text{reduced algebra uses }q}
$$

。

---

# 66. Unit convention

Material stress 使用 MPa。

因此：

$$
E_0(\varepsilon-\varepsilon^p)
$$

得到 MPa stress。

只有在 stress integration 成：

$$
N,M
$$

以及 nodal forces 时才进行 MPa → Pa / N unit conversion。

所以：

> 不应将 $10^6$ conversion factor 写进 stress mode 本身。

应在 force/stiffness integration layer 统一处理。

---

# 67. 当前应继续保留的 numerical safeguards

Tower v1 应继承：

- $H_\sigma>0$；
- $H_\beta>0$；
- $H_{\bar R}>0$；
- search-direction regularization；
- reduced least-squares condition monitoring；
- `rcond`；
- non-finite checks；
- damage upper bound；
- minimum spatial norm；
- small-mode rejection；
- residual stagnation detection；
- LATIN convergence 与 reduced residual 不混淆。

这些属于 numerical robustness，不是 theory changes。

---

# 68. Equilibrium-operator validation checklist

未来应检查：

1. zero source strain gives zero correction；
2. linear scaling；
3. superposition；
4. free-DOF residual：
$$
H^TM\sigma^{\mathrm{eq}}\approx0
$$
5. fixed-base reactions允许非零；
6. repeated RHS 复用同一 factorisation；
7. plastic source 与 damage residual source 均能通过同一 operator。

---

# 69. Temporal-update validation checklist

应测试：

1. one-mode synthetic recovery；
2. two-mode synthetic recovery；
3. constant $H_\sigma$；
4. variable $H_\sigma(t,q)$；
5. weighted residual decreases；
6. coefficient histories finite；
7. condition-history stable；
8. regularization prevents singular weights；
9. quadrature-volume weights correct；
10. reconstructed stress remains equilibrated。

---

# 70. Forcing validation checklist

分别验证：

$$
f^{\mathrm{sep}}
$$

与：

$$
f^{\mathrm{coupled}}
$$

。

在 small-damage regime，应检查：

$$
H_\sigma\Delta\tilde\sigma
$$

相对主 forcing 的量级。

这将直接检验 paper-separated approximation 在 tower steel fatigue 中是否成立。

---

# 71. Damage residual validation checklist

测试：

1. no-damage limit：
$$
\Delta\varepsilon^R\approx0
$$
2. fixed stress + changed damage：
$$
\Delta\varepsilon^R\neq0
$$
3. damage correction free-DOF equilibrium；
4. plastic + damage recombination；
5. local/global elastic-state consistency；
6. same reference operator reusable。

---

# 72. Hardening validation checklist

对所有 fibers / times 检查：

$$
\boxed{\beta=C\alpha}
$$

以及：

$$
\boxed{\bar R=R_\infty\bar r}
$$

。

同时：

$$
H_\beta>0
$$

$$
H_{\bar R}>0
$$

。

---

# 73. Damage-state validation checklist

检查：

$$
\dot D_{i+1}=\hat{\dot D}
$$

以及：

$$
0\le D\le D_{\max}
$$

。

最终：

$$
Y=Y(\sigma,D)
$$

必须严格成立。

尤其测试 tension/compression switching。

---

# 74. 尚未最终决定的实现问题 1：damage history inheritance

仍需在编码前判断：

```text
Option A:
copy local damage_rate and re-integrate D

Option B:
directly inherit local integrated D history
```

。

理论要求的核心是：

$$
b^-=0
$$

即 global stage 不修改 local damage evolution。

---

# 75. 尚未最终决定的实现问题 2：temporal coefficient semantics

原论文 Eq. (58) 使用：

$$
\Delta\lambda_j
$$

当前一维 code 使用：

```text
temporal_amplitude
```

。

未来必须明确：

- basis 保存 total amplitude 还是 correction amplitude；
- global stage 是 overwrite 还是 accumulate；
- enrichment 后 re-optimisation 的 coefficient semantic。

该问题建议等 Eq. (60)–(72) 完整推导后统一冻结。

---

# 76. 尚未最终决定的实现问题 3：DG0 与 backward Euler

当前一维代码把 backward-Euler type derivative 作为原论文 DG0 的工程实现。

Tower v1 可先保持这一离散。

但未来论文表述应谨慎：

> 使用与一维 reproduction 一致的 first-order piecewise-constant / backward-Euler temporal treatment，作为原论文 DG0 strategy 的实现方式。

不应在没有进一步离散推导时声称两者所有细节完全等价。

---

# 77. 已关闭的问题

本阶段已经明确：

### Eq. (58)–(59) 是否需要新 PGD 理论？

$$
\boxed{\text{No}}
$$

### $H_\sigma$ 是否为 structural matrix？

$$
\boxed{\text{No}}
$$

它是 fiber-local scalar field。

### Damage correction 是否进入 tower v1 plastic forcing？

$$
\boxed{\text{No}}
$$

### Damage correction 是否忽略？

$$
\boxed{\text{No}}
$$

### Damage branch 是否需要新 spatial operator？

$$
\boxed{\text{No}}
$$

### Hardening 是否 PGD reduced？

$$
\boxed{\text{No for tower v1}}
$$

### $D,\dot D,Y$ 是否 PGD reduced？

$$
\boxed{\text{No for tower v1}}
$$

。

---

# 78. 本阶段最重要的 formulation decision

$$
\boxed{\text{Tower v1 keeps the original plastic/damage separation as clean as possible}}
$$

具体：

```text
plastic correction
    → PGD reduced

damage elastic residual
    → separate reference projection

hardening variables
    → fiber-local linear global-stage ODEs

damage evolution
    → local history retained with b_minus = 0
```

。

---

# 79. Offline low-rank analysis 的角色

此前 tower FOM 的 SVD/HOSVD 仍然有效，但角色限定为：

$$
\boxed{\text{offline evidence of compressibility}}
$$

而不是：

$$
\boxed{\text{proof that all constitutive fields should be PGD unknowns}}
$$

。

因此当前没有因为 damage / hardening 低秩就扩大 PGD variable set。

---

# 80. 当前完整 theoretical chain

上一阶段：

$$
\bar\varepsilon_j^p
$$

通过：

$$
\mathcal E_{\mathrm{tower}}
$$

得到：

$$
\bar\varepsilon_j
$$

和：

$$
\bar\sigma_j
$$

。

本阶段：

$$
P,S
$$

固定后求：

$$
\Delta\lambda_j(t)
$$

得到：

$$
\Delta\varepsilon^p,\Delta\sigma'
$$

。

同时：

$$
\Delta\varepsilon^R
$$

通过同一个：

$$
\mathcal E_{\mathrm{tower}}
$$

得到：

$$
\Delta\tilde\sigma
$$

。

Hardening / damage local updates 再补齐：

$$
\alpha,\beta,\bar r,\bar R,D,\dot D,Y
$$

最终得到：

$$
\boxed{s_{i+1}\in A}
$$

。

因此：

$$
\boxed{\text{tower LATIN global stage is theoretically closed}}
$$

。

---

# 81. 为什么现在还不能称为完整 tower LATIN-PGD solver

仍缺：

1. Eq. (60) saturation criterion；
2. 判断是否需要 enrichment；
3. Eq. (61) new separated pair；
4. alternating update；
5. new spatial-mode solve；
6. new temporal-mode solve；
7. normalization；
8. orthogonalisation；
9. mode acceptance；
10. enrichment stopping；
11. outer LATIN convergence；
12. tower-level numerical validation。

因此：

$$
\boxed{\text{global-stage foundation complete}\neq\text{complete tower LATIN-PGD solver}}
$$

。

---

# 82. 下一阶段：Eq. (60) saturation criterion

下一步应该首先解决：

> **已有 spatial basis 在重新优化 temporal functions 后，何时认为它已经无法继续充分表示新的 LATIN correction？**

也就是：

$$
\boxed{\text{Eq. (60) saturation criterion}}
$$

。

它决定：

```text
basis sufficient
    → continue LATIN with reused basis

basis saturated
    → add new PGD mode
```

。

---

# 83. 为什么不能跳过 Eq. (60) 直接 enrichment

如果每个 LATIN iteration 都无条件增加 mode：

- basis 快速膨胀；
- reduced-order benefit 下降；
- modes 可能重复；
- conditioning 变差；
- mode count 难解释；
- 不符合原论文 hybrid PGD strategy。

所以 Eq. (60) 是：

$$
\boxed{\text{temporal reuse}\longrightarrow\text{new-mode enrichment}}
$$

之间的 gate。

---

# 84. 下一阶段建议路线

严格按：

```text
Step 1:
    Eq. (60) saturation criterion

Step 2:
    Eq. (61) new separated pair

Step 3:
    Eq. (62)–(66) temporal part

Step 4:
    Eq. (67)–(72) spatial part

Step 5:
    tower fiber weighting and normalisation

Step 6:
    tower enrichment acceptance and safeguards

Step 7:
    freeze code interfaces

Step 8:
    begin implementation
```

。

---

# 85. 当前 Git 状态说明

上一份正式总结：

```text
docs/2026-08-13-tower-latin-pgd-eq47-53-spatial-derivation-summary.md
```

已提交并推送到：

```text
feature/offshore-wind-turbine-tower-fatigue
```

对应 commit：

```text
988122a
docs: derive tower LATIN-PGD spatial operator for Eqs 47-53
```

。

自该总结之后，本阶段主要完成理论推导与 formulation decisions，尚未因此修改 tower solver 代码。

---

# 86. 本阶段相关 repo 文件

```text
latin/pgd_time_update.py
latin/search_directions.py
latin/global_stage.py
latin/pgd_global_stage.py
latin/local_stage.py
latin/equilibrium_operator.py
latin/pgd_basis.py

material/viscoplastic_damage_1d.py
```

未来 tower implementation 对接：

```text
fem/fiber_section.py
fem/viscoplastic_fiber_section.py
fem/beam_column_2d.py
fem/viscoplastic_beam_column_2d.py
fem/viscoplastic_tower_system_2d.py
```

。

---

# 87. 本阶段最终结论汇总

1. Eq. (58)–(59) 可直接迁移到 tower fiber material-point space。
2. Existing spatial basis 用 $P$ 表示。
3. Equilibrated stress basis 用 $S$ 表示。
4. Fixed-basis update 中 $P,S$ 保持不变。
5. Reduced unknown 是 $\Delta\lambda(t)$。
6. Tower spatial integration weight 为 $A_fw_gJ_e$。
7. Eq. (59) 是 $H_\sigma^{-1}$ weighted least-squares projection。
8. 当前 backward-Euler temporal treatment 可作为 tower v1 起点。
9. $H_\sigma$ 是 fiber-local scalar tangent field。
10. $H_\sigma$ 可以 space-time varying。
11. $\mathcal E_{\mathrm{tower}}$ 必须 fixed reference elastic。
12. $H_\sigma$ 必须 regularized。
13. Paper-separated forcing 与 current repo coupled forcing 存在明确差异。
14. Tower v1 采用 paper-separated forcing。
15. 当前一维 coupled forcing 暂不回改。
16. Coupled forcing 保留为 future sensitivity variant。
17. Damage correction 不被忽略。
18. Damage residual branch 独立计算。
19. $\Delta\varepsilon^R$ 是 known residual source strain，不是 damage increment。
20. Damage branch 复用 $\mathcal E_{\mathrm{tower}}$。
21. Damage branch 不 PGD reduced。
22. Plastic 与 damage stress corrections 各自 weakly equilibrated。
23. 二者叠加仍满足 equilibrium。
24. Hardening variables 不需要 tower equilibrium。
25. Hardening variables 不 PGD reduced。
26. $\alpha,\beta$ 由 local linear ODE update。
27. $\bar r,\bar R$ 由 local linear ODE update。
28. $H_\beta$ 与 $H_{\bar R}$ 使用 diagonal Hessian approximation。
29. 两者均 regularized。
30. Hardening 与 plastic PGD 在 same global stage 内解耦。
31. Hardening 与 plastic evolution 在 LATIN iterations 间重新耦合。
32. Damage search direction 使用 $b^-=0$。
33. Global stage 接受 local damage-rate history。
34. $D$ 由 damage-rate history恢复或直接继承 local integrated history，具体离散方式待定。
35. $Y$ 必须使用 final global stress 与 damage 重算。
36. $D,\dot D,Y$ 不 PGD reduced。
37. Tower global state 可统一存为 $(N_t,N_q)$ fiber material-point arrays。
38. 当前 coarse tower 中 $N_q=320$。
39. 一个完整 tower LATIN global stage 的理论框架已经闭合。
40. 下一步应先处理 Eq. (60) saturation criterion，而不是直接进入 enrichment。

---

# 88. 一句话阶段定位

> **自 Eq. (47)–(53) 空间算子总结之后，本阶段已经把 fixed-basis temporal update、fiber-level search directions、plastic forcing、plastic/damage formulation choice、damage residual projection、hardening global update 与 damage thermodynamic update 全部迁移到 fiber beam-column tower，并由此形成了一个理论上闭合的 tower LATIN global stage；下一步可以正式进入 Eq. (60) spatial-basis saturation criterion。**
