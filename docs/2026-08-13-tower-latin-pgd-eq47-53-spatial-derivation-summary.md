# 海上风机塔筒 LATIN-PGD：原论文 Eq. (47)–(53) 空间问题迁移与推导总结

**日期：2026-08-13**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**阶段范围：原论文 Eq. (47)–(53)**  
**阶段目标：在尽量保持原论文 x-t LATIN-PGD 结构不变的前提下，将其空间问题从 bar/continuum 迁移到 fiber beam-column offshore wind turbine tower**  
**下一阶段：Eq. (58)–(59) 固定 spatial basis 下的 temporal update**

---

# 1. 本阶段的核心任务

当前路线不立即引入 n-τ-x 三变量 PGD，而是首先确认原论文的 x-t LATIN-PGD 是否可以直接应用到纤维梁柱海上风机塔筒。

本阶段集中解决四个问题：

1. 原论文 Eq. (47) 中的 plastic-strain spatial mode 在塔筒中是什么；
2. 原论文 Eq. (50)–(51) 的 equilibrium operator 如何在 fiber beam-column tower 中离散；
3. 原论文 Eq. (52)–(53) 的 stress correction operator 如何迁移；
4. 为什么由该 operator 生成的 stress spatial mode 自动满足塔筒自由 DOF 上的离散弱平衡。

本阶段暂不进入：

- Eq. (58)–(59) temporal update；
- new mode alternating solver；
- enrichment Eq. (61)–(72)；
- damage correction 的完整 tower-PGD 化；
- n-τ-x 三变量扩展。

---

# 2. 原论文 Eq. (47)–(53) 的基本逻辑

原论文首先对 plastic correction 引入 separated mode：

$$
\Delta\varepsilon^p_{i+1}(x,t)=\lambda^p(t)\bar\varepsilon^p(x)
$$

对应：

$$
\Delta\dot\varepsilon^p_{i+1}(x,t)=\dot\lambda^p(t)\bar\varepsilon^p(x)
$$

随后通过空间平衡求取 compatible strain mode：

$$
\bar\varepsilon=\mathcal E\bar\varepsilon^p
$$

再构造 stress correction：

$$
\Delta\sigma'=\lambda(t)\mathbb C(\mathcal E-\mathbb I)\bar\varepsilon^p
$$

定义：

$$
\bar{\mathbb C}=\mathbb C(\mathcal E-\mathbb I)
$$

于是原论文 Eq. (53) 可以写成：

$$
\Delta\dot\varepsilon^p_{i+1}(x,t)=\dot\lambda(t)\bar\varepsilon^p(x)
$$

$$
\Delta\sigma'_{i+1}(x,t)=\lambda(t)\bar{\mathbb C}\bar\varepsilon^p(x)
$$

核心链条是：

$$
\boxed{\bar\varepsilon^p\longrightarrow\mathcal E\bar\varepsilon^p\longrightarrow\mathbb C(\mathcal E-\mathbb I)\bar\varepsilon^p}
$$

因此主空间未知量是 $\bar\varepsilon^p$；compatible strain 和 stress mode 都是由空间算子派生出来的。

---

# 3. 塔筒中的空间坐标

从数学上，原论文中的 $x$ 是物理空间坐标，而不是有限元索引。

对于当前二维 Euler-Bernoulli 纤维梁柱塔筒，连续意义下可写为：

$$
x_{\mathrm{tower}}\equiv(s,y)
$$

其中：

- $s$：沿塔筒轴线的空间坐标；
- $y$：截面内 fiber 相对于中性轴的位置。

有限元离散后：

$$
(s,y)\longrightarrow(e,g,f)
$$

其中：

- $e$：beam element index；
- $g$：beam Gauss integration point；
- $f$：section fiber index。

因此更严谨的关系是：

$$
x\xrightarrow{\mathrm{FE\ discretisation}}(e,g,f)
$$

而不是简单把 $x$ 与 $(e,g,f)$ 当成完全相同的数学对象。

---

# 4. 当前 fiber beam-column 塔筒的运动学

第 $e$ 个二维 Euler-Bernoulli 梁单元局部自由度为：

$$
\mathbf d_e^\ell=[u_1,v_1,\theta_1,u_2,v_2,\theta_2]^T
$$

第 $g$ 个 beam Gauss point 的 generalized strain 为：

$$
\eta_{eg}=\begin{bmatrix}\varepsilon_{0,eg}\\\kappa_{eg}\end{bmatrix}
$$

并满足：

$$
\eta_{eg}=\mathbf B_{eg}\mathbf d_e^\ell
$$

fiber strain 采用：

$$
\varepsilon_{egf}=\varepsilon_{0,eg}-y_{egf}\kappa_{eg}
$$

定义：

$$
\mathbf a_{egf}=\begin{bmatrix}1\\-y_{egf}\end{bmatrix}
$$

则：

$$
\varepsilon_{egf}=\mathbf a_{egf}^T\eta_{eg}
$$

若 $\mathbf P_e$ 为单元 DOF 提取矩阵，$\mathbf T_e$ 为 global-to-local transformation，则：

$$
\mathbf d_e^\ell=\mathbf T_e\mathbf P_e\mathbf U
$$

所以：

$$
\varepsilon_{egf}=\mathbf a_{egf}^T\mathbf B_{eg}\mathbf T_e\mathbf P_e\mathbf U
$$

塔筒运动学链条为：

$$
\boxed{\mathbf U\longrightarrow\mathbf d_e^\ell\longrightarrow\eta_{eg}\longrightarrow\varepsilon_{egf}}
$$

---

# 5. Eq. (47) 在塔筒中的定义

连续形式：

$$
\Delta\varepsilon^p(s,y,t)=\lambda(t)\bar\varepsilon^p(s,y)
$$

离散形式：

$$
\boxed{\Delta\varepsilon^p_{egf}(t)=\lambda(t)\bar\varepsilon^p_{egf}}
$$

对应 plastic-strain rate：

$$
\boxed{\Delta\dot\varepsilon^p_{egf}(t)=\dot\lambda(t)\bar\varepsilon^p_{egf}}
$$

因此：

> **$\bar\varepsilon^p_{egf}$ 是定义在 tower fiber material points 上的一维轴向塑性应变修正空间模态。**

它不是 nodal displacement mode、tower vibration mode、curvature mode、moment mode 或 damage mode。

---

# 6. 一个 tower PGD spatial mode 的离散维数

若：

- beam elements 数为 $N_e$；
- 每个 element 的 beam Gauss points 数为 $N_g$；
- 每个 section 的 fibers 数为 $N_f$；

则：

$$
N_{\mathrm{fiberGP}}=N_eN_gN_f
$$

一个 spatial mode：

$$
\bar{\varepsilon}^p_j\in\mathbb R^{N_{\mathrm{fiberGP}}}
$$

当前 100-cycle 粗模型采用：

$$
N_e=10,\qquad N_g=2,\qquad N_f=16
$$

因此：

$$
N_{\mathrm{fiberGP}}=320
$$

即：

$$
\boxed{\bar{\varepsilon}^p_j\in\mathbb R^{320}}
$$

这 320 个量是 material-point spatial field values，不是 320 个新的结构自由度。

---

# 7. 为什么 plastic-strain spatial mode 本身不需要兼容

一般情况下：

$$
\bar{\varepsilon}^p\notin\operatorname{Range}(\mathbf H)
$$

也就是说，它通常不能直接由某一个 beam nodal displacement vector 产生。

这是合理的，因为 plastic strain 是材料点内部变量，而不是结构总应变。

原论文 Eq. (50) 的作用就是：

> 给定一个 material-point plastic-strain field，求一个满足结构运动学和整体平衡的 compatible total-strain correction。

所以必须区分：

$$
\bar{\varepsilon}^p\neq\bar{\varepsilon}
$$

其中：

- $\bar{\varepsilon}^p$：plastic material-point mode；
- $\bar{\varepsilon}$：beam-compatible total-strain mode。

---

# 8. 原论文 Eq. (50) 的 tower weak form

原论文 Eq. (50) 可理解为：

$$
\int_\Omega\mathbb C(\bar\varepsilon-\bar\varepsilon^p):\varepsilon(\bar u^*)\,d\Omega=0
$$

对 fiber beam-column：

$$
d\Omega=dA\,ds
$$

采用 reference elastic modulus $E_0$：

$$
\sum_e\int_{L_e}\int_{A(s)}E_0(\bar\varepsilon-\bar\varepsilon^p)\delta\varepsilon\,dA\,ds=0
$$

截面采用 fiber quadrature，梁长采用 Gauss quadrature 后：

$$
\boxed{\sum_{e,g,f}w_gJ_eA_{egf}E_0(\bar\varepsilon_{egf}-\bar\varepsilon^p_{egf})\delta\varepsilon_{egf}=0}
$$

这就是原论文 Eq. (50) 在当前 tower fiber discretisation 上的直接对应。

---

# 9. fiber → section：等效 plastic section resultants

定义：

$$
\mathbf a_{egf}=\begin{bmatrix}1\\-y_{egf}\end{bmatrix}
$$

参考弹性截面矩阵：

$$
\boxed{\mathbf D^0_{eg}=E_0\sum_fA_{egf}\mathbf a_{egf}\mathbf a_{egf}^T}
$$

展开：

$$
\mathbf D^0_{eg}=E_0\begin{bmatrix}\sum_fA_f&-\sum_fA_fy_f\\-\sum_fA_fy_f&\sum_fA_fy_f^2\end{bmatrix}
$$

对于关于中性轴对称的圆环截面：

$$
\sum_fA_fy_f\approx0
$$

所以：

$$
\mathbf D^0_{eg}\approx\begin{bmatrix}E_0A_{eg}&0\\0&E_0I_{eg}\end{bmatrix}
$$

plastic spatial mode 对应的等效截面广义力定义为：

$$
\boxed{\bar{\mathbf r}^p_{eg}=E_0\sum_fA_{egf}\mathbf a_{egf}\bar\varepsilon^p_{egf}}
$$

其中：

$$
\bar{\mathbf r}^p_{eg}=\begin{bmatrix}\bar N^p_{eg}\\\bar M^p_{eg}\end{bmatrix}
$$

具体：

$$
\bar N^p_{eg}=E_0\sum_fA_f\bar\varepsilon^p_{egf}
$$

$$
\bar M^p_{eg}=-E_0\sum_fA_fy_f\bar\varepsilon^p_{egf}
$$

因此：

$$
\boxed{\bar\varepsilon^p_{egf}\longrightarrow(\bar N^p_{eg},\bar M^p_{eg})}
$$

---

# 10. section → element → tower

等效单元局部 load-like vector：

$$
\boxed{\bar{\mathbf f}^{p,\ell}_e=\sum_gw_gJ_e\mathbf B_{eg}^T\bar{\mathbf r}^p_{eg}}
$$

转换至全局：

$$
\boxed{\bar{\mathbf f}^p_e=\mathbf T_e^T\bar{\mathbf f}^{p,\ell}_e}
$$

全局装配：

$$
\boxed{\bar{\mathbf F}^p=\operatorname{Assembly}_e(\bar{\mathbf f}^p_e)}
$$

于是得到：

$$
\boxed{\bar\varepsilon^p_{egf}\longrightarrow\bar{\mathbf r}^p_{eg}\longrightarrow\bar{\mathbf f}^p_e\longrightarrow\bar{\mathbf F}^p}
$$

也就是：

$$
\boxed{\text{fiber}\longrightarrow\text{section}\longrightarrow\text{element}\longrightarrow\text{tower}}
$$

---

# 11. tower reference elastic stiffness

第 $e$ 个 element 的 reference local stiffness：

$$
\mathbf K^{0,\ell}_e=\sum_gw_gJ_e\mathbf B_{eg}^T\mathbf D^0_{eg}\mathbf B_{eg}
$$

global form：

$$
\mathbf K^0_e=\mathbf T_e^T\mathbf K^{0,\ell}_e\mathbf T_e
$$

装配后：

$$
\boxed{\mathbf K^0=\operatorname{Assembly}_e(\mathbf K^0_e)}
$$

因此 tower Eq. (50) 可写成：

$$
\boxed{\mathbf K^0\bar{\mathbf U}=\bar{\mathbf F}^p}
$$

施加固定塔底边界条件后：

$$
\boxed{\mathbf K^0_{FF}\bar{\mathbf U}_F=\bar{\mathbf F}^p_F}
$$

其中 $F$ 表示 free DOFs。

---

# 12. Eq. (50) 的物理解释

$\bar{\mathbf F}^p$ 不是风荷载，也不是外部 correction load。

它是由：

$$
\bar{\varepsilon}^p
$$

诱导出来的等效 plastic eigenstrain load。

因此 Eq. (50) 的物理含义是：

$$
\boxed{\text{plastic eigenstrain}\longrightarrow\text{self-equilibrated structural correction}}
$$

完整链条：

$$
\boxed{\bar{\varepsilon}^p\longrightarrow\bar{\mathbf F}^p\longrightarrow\bar{\mathbf U}\longrightarrow\bar{\varepsilon}}
$$

---

# 13. Eq. (51)：tower equilibrium projection operator

定义全局 fiber kinematic matrix：

$$
\boxed{\varepsilon_f=\mathbf H\mathbf U_F}
$$

第 $(e,g,f)$ 行：

$$
\boxed{\mathbf H_{egf}=\mathbf a_{egf}^T\mathbf B_{eg}\mathbf T_e\mathbf P_{e,F}}
$$

定义几何积分权重矩阵：

$$
\boxed{\mathbf M=\operatorname{diag}(A_{egf}w_gJ_e)}
$$

定义 reference material operator：

$$
\boxed{\mathbf C_0=E_0\mathbf I}
$$

则：

$$
\boxed{\mathbf K^0_{FF}=\mathbf H^T\mathbf M\mathbf C_0\mathbf H}
$$

以及：

$$
\boxed{\bar{\mathbf F}^p_F=\mathbf H^T\mathbf M\mathbf C_0\bar{\varepsilon}^p}
$$

因此：

$$
\bar{\mathbf U}_F=(\mathbf H^T\mathbf M\mathbf C_0\mathbf H)^{-1}\mathbf H^T\mathbf M\mathbf C_0\bar{\varepsilon}^p
$$

compatible fiber strain：

$$
\bar{\varepsilon}=\mathbf H\bar{\mathbf U}_F
$$

所以：

$$
\boxed{\bar{\varepsilon}=\mathcal E_{\mathrm{tower}}\bar{\varepsilon}^p}
$$

其中：

$$
\boxed{\mathcal E_{\mathrm{tower}}=\mathbf H(\mathbf H^T\mathbf M\mathbf C_0\mathbf H)^{-1}\mathbf H^T\mathbf M\mathbf C_0}
$$

这就是原论文 Eq. (51) 的 tower counterpart。

---

# 14. equilibrium operator 的投影意义

$\mathcal E_{\mathrm{tower}}$ 把任意 fiber material-point field：

$$
\bar{\varepsilon}^p
$$

映射到：

$$
\operatorname{Range}(\mathbf H)
$$

中的 compatible strain field。

因此可理解为：

$$
\boxed{\text{fiber material-point strain space}\longrightarrow\text{beam-compatible strain space}}
$$

的 weighted equilibrium/compatibility projection。

权重由：

$$
\mathbf M\mathbf C_0
$$

决定。

---

# 15. 为什么必须首先使用 reference elastic modulus

第一版应基于：

$$
\boxed{E_0}
$$

构造 $\mathbf C_0$ 与 $\mathbf K^0$，而不是使用 FOM 当前 nonlinear numerical tangent。

若采用 current tangent：

$$
\mathbf K=\mathbf K(t)
$$

则：

$$
\mathcal E=\mathcal E(t)
$$

也会变成 time-dependent operator。

此时原论文中：

$$
\mathcal E[\lambda(t)\bar\varepsilon^p]
$$

不能再简单提取 $\lambda(t)$。

因此会破坏 Eq. (47)–(53) 的 x-t separation。

所以第一版必须保持：

$$
\boxed{\mathbf K^0\text{ 为固定的 reference elastic spatial operator}}
$$

damage、hardening 和 nonlinear elastic degradation 仍由 LATIN 的其他 correction/local-update 环节处理。

---

# 16. Eq. (52)：tower stress spatial mode

每个 fiber 上的 elastic strain spatial correction：

$$
\bar\varepsilon^e_{egf}=\bar\varepsilon_{egf}-\bar\varepsilon^p_{egf}
$$

因此：

$$
\boxed{\bar\sigma_{egf}=E_0(\bar\varepsilon_{egf}-\bar\varepsilon^p_{egf})}
$$

向量形式：

$$
\bar{\sigma}=\mathbf C_0(\bar{\varepsilon}-\bar{\varepsilon}^p)
$$

代入 Eq. (51)：

$$
\boxed{\bar{\sigma}=\mathbf C_0(\mathcal E_{\mathrm{tower}}-\mathbf I)\bar{\varepsilon}^p}
$$

定义 tower stress correction operator：

$$
\boxed{\bar{\mathbf C}_{\mathrm{tower}}=\mathbf C_0(\mathcal E_{\mathrm{tower}}-\mathbf I)}
$$

于是：

$$
\boxed{\bar{\sigma}=\bar{\mathbf C}_{\mathrm{tower}}\bar{\varepsilon}^p}
$$

---

# 17. stress spatial mode 不是独立 PGD unknown

一个 PGD spatial mode 中，真正的主空间未知量仍然是：

$$
\boxed{\bar{\varepsilon}^p_j}
$$

由它生成：

$$
\bar{\varepsilon}_j=\mathcal E_{\mathrm{tower}}\bar{\varepsilon}^p_j
$$

以及：

$$
\bar{\sigma}_j=\mathbf C_0(\mathcal E_{\mathrm{tower}}-\mathbf I)\bar{\varepsilon}^p_j
$$

因此：

$$
\boxed{\bar{\varepsilon}^p_j\longrightarrow\bar{\varepsilon}_j\longrightarrow\bar{\sigma}_j}
$$

如果把 $\bar{\sigma}_j$ 作为独立 spatial unknown，就不能自动保证结构平衡。

---

# 18. 为什么同一个 temporal amplitude 可以继续使用

若：

$$
\Delta\varepsilon^p(t)=\lambda(t)\bar{\varepsilon}^p
$$

由于 $\mathcal E_{\mathrm{tower}}$ 为固定纯空间线性算子：

$$
\Delta\varepsilon(t)=\mathcal E_{\mathrm{tower}}\Delta\varepsilon^p(t)
$$

所以：

$$
\Delta\varepsilon(t)=\lambda(t)\mathcal E_{\mathrm{tower}}\bar{\varepsilon}^p
$$

即：

$$
\boxed{\Delta\varepsilon(t)=\lambda(t)\bar{\varepsilon}}
$$

进一步：

$$
\Delta\sigma'(t)=\mathbf C_0[\Delta\varepsilon(t)-\Delta\varepsilon^p(t)]
$$

所以：

$$
\boxed{\Delta\sigma'(t)=\lambda(t)\bar{\sigma}}
$$

这就是原论文 Eq. (52)–(53) 中 displacement/strain/stress correction 共用同一 temporal amplitude 的原因。

---

# 19. Eq. (53) 的 tower counterpart

第 $j$ 个 PGD mode：

$$
\boxed{\Delta\dot{\varepsilon}^p_j(t)=\dot\lambda_j(t)\bar{\varepsilon}^p_j}
$$

$$
\boxed{\Delta\sigma'_j(t)=\lambda_j(t)\bar{\sigma}_j}
$$

其中：

$$
\boxed{\bar{\sigma}_j=\mathbf C_0(\mathcal E_{\mathrm{tower}}-\mathbf I)\bar{\varepsilon}^p_j}
$$

已有 $m$ 个 modes：

$$
\boxed{\Delta\dot{\varepsilon}^p(t)=\sum_{j=1}^{m}\dot\lambda_j(t)\bar{\varepsilon}^p_j}
$$

$$
\boxed{\Delta\sigma'(t)=\sum_{j=1}^{m}\lambda_j(t)\bar{\sigma}_j}
$$

因此当前仍是原论文的 x-t PGD，而不是三变量 PGD。

---

# 20. 为什么 stress spatial mode 自动满足 tower equilibrium

tower Eq. (50)：

$$
\mathbf H^T\mathbf M\mathbf C_0(\bar{\varepsilon}-\bar{\varepsilon}^p)=\mathbf 0
$$

定义：

$$
\bar{\sigma}=\mathbf C_0(\bar{\varepsilon}-\bar{\varepsilon}^p)
$$

立即得到：

$$
\boxed{\mathbf H^T\mathbf M\bar{\sigma}=\mathbf 0}
$$

而：

$$
\mathbf H^T\mathbf M\bar{\sigma}
$$

就是 fiber stresses 经过 section integration、beam Gauss integration、element transformation 和 global assembly 后的 free-DOF internal force vector。

所以：

$$
\boxed{\bar{\mathbf F}_{\mathrm{int},F}=\mathbf 0}
$$

即每一个 stress spatial mode 本身就在自由结构自由度上满足离散弱平衡。

---

# 21. self-equilibrated 的严格含义

不能简单说：

> “stress mode 的所有内力都为零。”

正确的是：

$$
\boxed{\bar{\mathbf F}_{\mathrm{int},F}=0}
$$

即 free DOFs 上 residual 为零。

对于固定塔底 constrained DOFs，可以存在：

$$
\boxed{\bar{\mathbf R}_C\neq0}
$$

也就是说 eigenstrain correction 可以产生 base reactions。

所以更准确的表述是：

> **stress correction 在无附加 external correction load 时满足离散弱平衡；自由 DOFs 上无不平衡节点力，而约束 DOFs 可以承担 reaction forces。**

---

# 22. weak equilibrium 不等于每个截面 N=M=0

即使：

$$
\mathbf H^T\mathbf M\bar{\sigma}=0
$$

局部仍然可以有：

$$
\bar N_{eg}\neq0
$$

以及：

$$
\bar M_{eg}\neq0
$$

真正成立的是：

$$
\boxed{\delta\mathbf U_F^T\bar{\mathbf F}_{\mathrm{int},F}=0\qquad\forall\delta\mathbf U_F}
$$

因此应称为：

$$
\boxed{\text{weakly equilibrated stress field}}
$$

而不是“zero section force field”。

---

# 23. projector 形式下的平衡证明

由于：

$$
\mathcal E_{\mathrm{tower}}=\mathbf H(\mathbf H^T\mathbf M\mathbf C_0\mathbf H)^{-1}\mathbf H^T\mathbf M\mathbf C_0
$$

所以：

$$
\mathbf H^T\mathbf M\mathbf C_0\mathcal E_{\mathrm{tower}}=\mathbf H^T\mathbf M\mathbf C_0
$$

因此：

$$
\mathbf H^T\mathbf M\mathbf C_0(\mathcal E_{\mathrm{tower}}-\mathbf I)=\mathbf 0
$$

对任意 $\bar{\varepsilon}^p$：

$$
\mathbf H^T\mathbf M\mathbf C_0(\mathcal E_{\mathrm{tower}}-\mathbf I)\bar{\varepsilon}^p=\mathbf 0
$$

又因为：

$$
\bar{\sigma}=\mathbf C_0(\mathcal E_{\mathrm{tower}}-\mathbf I)\bar{\varepsilon}^p
$$

所以严格得到：

$$
\boxed{\mathbf H^T\mathbf M\bar{\sigma}=\mathbf 0}
$$

因此 equilibrium property 是 operator definition 的数学结果，而不是经验观察。

---

# 24. 与现有 1D LATIN-PGD 代码的对应

当前一维杆 `latin/equilibrium_operator.py` 已经实现：

$$
\varepsilon=\mathcal E q
$$

以及：

$$
\sigma=C(\varepsilon-q)
$$

其 source strain $q$ 可以对应 plastic correction 或 damage residual correction。

当前 `latin/pgd_basis.py` 中一个 `PGDMode1D` 保存：

```text
spatial_plastic_strain
spatial_stress
temporal_amplitude
temporal_rate
```

并按：

$$
\Delta\dot\varepsilon^p_j=\dot\lambda_j(t)\bar\varepsilon^p_j
$$

$$
\Delta\sigma'_j=\lambda_j(t)\bar\sigma_j
$$

进行重构。

因此 tower-PGD 不需要改变 PGD mode 的基本概念，只需要把空间场从：

```text
n_elements
```

扩展为：

```text
n_elements × n_gauss × n_fibers
```

或 flatten 成：

```text
n_fiber_material_points
```

。

---

# 25. 与现有 tower FOM 的对应

当前已有模块：

```text
fem/fiber_section.py
fem/viscoplastic_fiber_section.py
fem/beam_column_2d.py
fem/viscoplastic_beam_column_2d.py
fem/viscoplastic_tower_system_2d.py
```

已经提供：

```text
fiber strain / stress
        ↓
section resultant [N, M]
        ↓
B^T [N, M]
        ↓
element internal force
        ↓
coordinate transformation
        ↓
global assembly
```

所以构造 $\mathcal E_{\mathrm{tower}}$ 不需要重写 beam finite element mechanics。

后续只需增加一个独立的 reference-elastic tower equilibrium projection operator。

---

# 26. 建议的 tower equilibrium operator 数据流

建议 future operator 按下列数据流组织：

```text
bar_eps_p[e, g, f]
        ↓
equivalent plastic section resultants
        ↓
bar_r_p[e, g, :] = [N_p, M_p]
        ↓
equivalent element plastic load
        ↓
bar_f_p_element[e, :]
        ↓
global equivalent plastic load
        ↓
bar_F_p
        ↓
reference elastic solve
        ↓
bar_U
        ↓
beam generalized strains
        ↓
bar_eta[e, g, :]
        ↓
compatible fiber strains
        ↓
bar_eps[e, g, f]
        ↓
stress spatial mode
        ↓
bar_sigma[e, g, f]
```

其中：

$$
\bar\sigma_{egf}=E_0(\bar\varepsilon_{egf}-\bar\varepsilon^p_{egf})
$$

---

# 27. 一个 tower PGD mode 的最终定义

主空间未知量：

$$
\boxed{\bar{\varepsilon}^p_j\in\mathbb R^{N_{\mathrm{fiberGP}}}}
$$

位移空间模态：

$$
\boxed{\bar{\mathbf U}_{F,j}=(\mathbf K^0_{FF})^{-1}\mathbf H^T\mathbf M\mathbf C_0\bar{\varepsilon}^p_j}
$$

compatible strain spatial mode：

$$
\boxed{\bar{\varepsilon}_j=\mathcal E_{\mathrm{tower}}\bar{\varepsilon}^p_j}
$$

stress spatial mode：

$$
\boxed{\bar{\sigma}_j=\mathbf C_0(\mathcal E_{\mathrm{tower}}-\mathbf I)\bar{\varepsilon}^p_j}
$$

temporal quantities：

$$
\lambda_j(t)
$$

以及：

$$
\dot\lambda_j(t)
$$

所以一个 mode 最小可理解为：

$$
\boxed{\{\bar{\varepsilon}^p_j,\bar{\sigma}_j,\lambda_j(t),\dot\lambda_j(t)\}}
$$

其中 $\bar{\sigma}_j$ 是派生量，不是独立未知量。

---

# 28. Eq. (47)–(53) 的 tower mapping 总表

| 原论文对象 | Tower counterpart | 数值位置 |
|---|---|---|
| $x$ | $(s,y)$ | FE 后为 $(e,g,f)$ |
| $\bar\varepsilon^p(x)$ | fiber axial plastic-strain mode | fiber material points |
| $\bar u(x)$ | tower displacement spatial mode | beam nodal DOFs |
| $\bar\varepsilon(x)$ | compatible fiber strain mode | fiber material points |
| $\mathcal E$ | $\mathcal E_{\mathrm{tower}}$ | reference elastic equilibrium projection |
| $\mathbb C$ | $\mathbf C_0$ | reference elastic fiber modulus operator |
| $\mathbb C(\mathcal E-I)$ | $\mathbf C_0(\mathcal E_{\mathrm{tower}}-I)$ | tower stress correction operator |
| $\bar\sigma(x)$ | fiber stress spatial mode | fiber material points |
| weak equilibrium | $\mathbf H^T\mathbf M\bar{\sigma}=0$ | free structural DOFs |
| $\lambda(t)$ | temporal amplitude | LATIN time grid |
| $\dot\lambda(t)$ | temporal rate | LATIN time grid |

---

# 29. 一维杆、二维连续体与塔筒的统一认识

一维杆：

```text
x → bar material point
```

二维连续体：

```text
x → continuum Gauss point
```

fiber beam-column tower：

```text
x → beam × section × fiber material point
```

真正保持不变的是：

$$
\boxed{\text{local material field}\longrightarrow\text{global equilibrium projection}\longrightarrow\text{equilibrated stress correction}}
$$

所以塔筒迁移本质上是：

$$
\boxed{\text{replace spatial discretisation/operator, not PGD dimensionality}}
$$

---

# 30. 本阶段关键理论结论

1. 继续采用原论文 x-t PGD 在数学上是可行的。
2. tower 中最自然的 PGD 主空间未知量是 fiber axial plastic-strain correction field。
3. compatible strain 不应作为独立 spatial PGD unknown。
4. stress spatial mode 不应作为独立 spatial PGD unknown。
5. 第一版 equilibrium operator 应使用固定 reference elastic operator。
6. Eq. (50) 可严格离散为 tower reference elastic equilibrium problem。
7. Eq. (51) 可严格定义为 $\mathcal E_{\mathrm{tower}}$。
8. Eq. (52) 可严格定义为 $\mathbf C_0(\mathcal E_{\mathrm{tower}}-\mathbf I)$。
9. 每个 stress spatial mode 自动满足 free DOFs 上的弱平衡。
10. 固定塔底可以存在 reaction forces。
11. 原论文 Eq. (53) 的 separated representation 在 tower 中保持不变。
12. Eq. (47)–(53) 的 tower spatial mapping 已经形成闭合数学链条。

---

# 31. 当前尚未解决的问题

下一阶段仍需解决：

1. fixed spatial basis 下 temporal amplitudes 如何更新；
2. search direction $H_\sigma$ 如何进入 tower temporal reduced system；
3. fiber weights 如何进入 Eq. (58)–(59)；
4. multiple spatial modes 的 coupled temporal equations 如何构造；
5. new mode enrichment 时 space-time alternating solver 如何写；
6. damage residual correction 如何复用 tower equilibrium operator；
7. tower Gram-Schmidt 应使用什么 weighted inner product；
8. spatial mode normalization 如何定义；
9. taper、非均匀 fiber area 与变截面如何进入 norm；
10. 如何构造 tower operator 的 numerical verification tests。

---

# 32. 下一阶段的明确入口

下一阶段正式进入原论文 Eq. (58)–(59)。

假设已经有 $m$ 个 spatial modes：

$$
\bar{\varepsilon}^p_1,\ldots,\bar{\varepsilon}^p_m
$$

及对应 stress modes：

$$
\bar{\sigma}_1,\ldots,\bar{\sigma}_m
$$

下一步需要求：

$$
\lambda_1(t),\ldots,\lambda_m(t)
$$

以及：

$$
\dot\lambda_1(t),\ldots,\dot\lambda_m(t)
$$

核心问题为：

> **在 tower fiber material-point space 中固定 spatial PGD basis 后，如何从 LATIN global-stage search direction 推导出原论文 Eq. (58)–(59) 对应的 reduced temporal system。**

这一步将标志着研究从“tower spatial operator 构造”正式进入“tower LATIN-PGD global-stage reduced solver”。

---

# 33. 建议后续代码模块边界

建议未来新增独立模块：

```text
latin/tower_equilibrium_operator.py
```

职责仅包括：

```text
input:
    fiber source strain field

reference spatial solve:
    fiber → section → element → global
    K0_FF * U_bar_F = F_source_F

output:
    compatible fiber strain
    equilibrated fiber stress
    nodal displacement mode
    base reactions
```

这样可以与现有：

```text
latin/equilibrium_operator.py
```

形成清晰的一一对应，并方便进行：

- unit tests；
- elastic reference comparison；
- Eq. (50) residual verification；
- stress equilibrium verification；
- PGD basis reuse。

---

# 34. 参考依据

## 原论文

Bhattacharyya, M., Fau, A., Nackenhorst, U., Néron, D., & Ladevèze, P. (2018).  
*A LATIN-based model reduction approach for the simulation of cycling damage*.  
Computational Mechanics, 62, 725–743.  
DOI: 10.1007/s00466-017-1523-z.

本阶段主要对应：

- Eq. (47)：plastic-strain separated mode；
- Eq. (48)–(49)：displacement/plastic temporal relation；
- Eq. (50)：spatial equilibrium problem；
- Eq. (51)：equilibrium operator；
- Eq. (52)：stress correction spatial operator；
- Eq. (53)：final separated representation。

## 当前项目关键代码

```text
material/viscoplastic_damage_1d.py

fem/fiber_section.py
fem/viscoplastic_fiber_section.py
fem/beam_column_2d.py
fem/viscoplastic_beam_column_2d.py
fem/viscoplastic_tower_system_2d.py

latin/equilibrium_operator.py
latin/pgd_basis.py
```

## 当前分支

```text
feature/offshore-wind-turbine-tower-fatigue
```

---

# 35. 一句话阶段定位

> **本阶段已经把原论文 Eq. (47)–(53) 从 bar/continuum spatial problem 严格迁移为 fiber beam-column tower spatial problem，并建立了从 fiber plastic-strain PGD mode 到 compatible strain mode、equilibrated stress mode 以及 global weak equilibrium 的完整离散数学链条；下一步可在此基础上进入 Eq. (58)–(59) temporal update。**
