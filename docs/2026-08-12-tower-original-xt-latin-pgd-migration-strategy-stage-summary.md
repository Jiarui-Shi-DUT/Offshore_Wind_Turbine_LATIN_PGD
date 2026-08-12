# 海上风机塔筒 LATIN-PGD：原论文 x-t 形式优先迁移策略阶段总结

**日期：2026-08-12**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**当前分支：`feature/offshore-wind-turbine-tower-fatigue`**  
**阶段定位：原论文 `x-t` LATIN-PGD 向纤维梁柱海上风机塔筒的直接迁移策略确认**

---

# 1. 本阶段讨论的背景

前期工作已经完成以下研究基础：

1. 成功复现 Bhattacharyya et al. (2018) 论文中一维三材料杆的 LATIN-PGD 循环损伤算例；
2. 在 GitHub 中建立并持续维护 LATIN-PGD 的模块化代码框架；
3. 建立海上风机钢塔筒二维纤维梁柱有限元模型；
4. 建立基于一维粘塑性-损伤本构的纤维材料积分点；
5. 完成塔筒弹性、非线性、反向加载、Bauschinger 效应、非对称循环、ratcheting、多周期损伤等逐级验证；
6. 建立稳定可运行的 100 周期塔筒全阶模型（FOM），并冻结 snapshot；
7. 对 100 周期塔筒响应开展 cycle-phase-space 张量化和 SVD/HOSVD 风格低秩诊断。

此前曾考虑进一步将原论文的

$$
q(x,t)
$$

扩展为

$$
q(n,\tau,x),
$$

其中 $n$ 表示慢循环坐标，$\tau$ 表示循环内快相位坐标。

经过本阶段重新讨论后，研究路线作出重要调整：

> **当前不应立即从原论文的 $x-t$ PGD 跨越到 $n-\tau-x$ 三变量 PGD。现阶段的首要目标应当是尽可能保持原论文 LATIN-PGD 数学结构不变，将其直接应用到由纤维梁柱单元组装而成的海上风机塔筒中。**

因此，当前研究主线重新明确为：

$$
\boxed{ \text{原论文 }x-t\text{ LATIN-PGD} \longrightarrow \text{纤维梁柱海上风机塔筒 }x-t\text{ LATIN-PGD} }
$$

而不是立即采用：

$$
\boxed{ x-t \longrightarrow n-\tau-x }
$$

---

# 2. 原论文方法的核心结构

参考论文：

> Bhattacharyya, M., Fau, A., Nackenhorst, U., Néron, D., & Ladevèze, P. (2018).  
> *A LATIN-based model reduction approach for the simulation of cycling damage*.  
> Computational Mechanics, 62, 725–743.

原论文的核心不是单纯对有限元计算结果进行后处理降阶，而是：

$$
\boxed{ \text{LATIN nonlinear solution framework} + \text{PGD model reduction embedded in the global stage} }
$$

LATIN 方法将问题划分为两个流形：

$$
A
$$

和

$$
\Gamma,
$$

并在整个时空域上交替迭代：

$$
s_i\in A \rightarrow \hat s_{i+1/2}\in\Gamma \rightarrow s_{i+1}\in A.
$$

其中：

- **Local stage**：处理局部、非线性的材料演化关系；
- **Global stage**：处理结构平衡、可容许条件及可线性处理的状态关系；
- **Search directions**：连接 local/global 两个阶段；
- **PGD**：嵌入 global stage，对主要修正场进行空间-时间分离；
- **LATIN indicator $\xi$**：衡量 local/global solution 的距离；
- **Saturation indicator $\zeta$**：判断已有 PGD basis 是否饱和、是否需要 enrichment。

---

# 3. 原论文 PGD 的基本形式仍应保留

原论文采用：

$$
\boxed{x-t}
$$

两变量分离。

例如，对塑性应变修正量：

$$
\Delta\varepsilon^p_{i+1}(x,t) = \lambda^p(t)\, \bar\varepsilon^p(x),
$$

对应塑性应变率：

$$
\Delta\dot\varepsilon^p_{i+1}(x,t) = \dot\lambda^p(t)\, \bar\varepsilon^p(x).
$$

在已有 $m$ 个 PGD 模态后：

$$
\Delta\dot\varepsilon^p_{i+1}(x,t) = \sum_{j=1}^{m} \Delta\dot\lambda_j(t) \bar\varepsilon^p_j(x).
$$

与塑性部分对应的应力修正表示为：

$$
\Delta\sigma'_{i+1}(x,t) = \sum_{j=1}^{m} \Delta\lambda_j(t) \mathbb{C} \bar\varepsilon^p_j(x).
$$

因此，现阶段应优先保持：

$$
\boxed{ \lambda_j(t)\,X_j(x) }
$$

这一原始 PGD 结构，而不是立刻改写为：

$$
N_j(n)T_j(\tau)X_j(x).
$$

---

# 4. 原论文已经证明 LATIN-PGD 不依赖于一维杆结构

原论文第 5 节依次给出：

1. 一维三材料杆；
2. 二维 L 形结构；
3. 二维开孔板。

这说明原论文的 LATIN-PGD 框架并不依赖于：

$$
\text{1D bar}
$$

这一特定结构类型。

在二维连续体算例中，方法的主要框架并没有改变：

$$
\boxed{ \text{Local stage} + \text{Global stage} + \text{Search directions} + \text{x-t PGD} + \text{Adaptive basis enrichment} }
$$

真正发生改变的是：

- 空间有限元离散；
- Gauss point 数量；
- 空间算子；
- search direction 的维度与计算量；
- global equilibrium problem 的结构。

因此，这为当前将方法进一步应用到纤维梁柱海上风机塔筒提供了直接的理论依据。

---

# 5. 当前塔筒问题不需要改变 PGD 的数学维数

对于一维杆：

$$
x
$$

表示杆轴线上的空间位置。

对于二维开孔板：

$$
x
$$

表示二维连续体中的材料位置。

对于当前纤维梁柱塔筒，同样可以继续使用：

$$
x,
$$

但需要把它理解为一个**广义空间坐标**。

在数值实现中：

$$
x \longleftrightarrow (e,g,f),
$$

其中：

- $e$：beam element index；
- $g$：beam integration Gauss point；
- $f$：section fiber index。

因此，塔筒中的 PGD 空间模态仍然可以写成：

$$
\bar\varepsilon^p_j(x),
$$

而程序中的离散存储可对应为：

$$
\bar\varepsilon^p_j(e,g,f).
$$

这意味着：

> **当前并不是把 PGD 从两变量扩展到三变量，而只是把原论文中 $x$ 所代表的空间离散，从 bar/continuum Gauss points 替换为 beam–section–fiber material points。**

这是本阶段最重要的理论定位之一。

---

# 6. 三类问题之间的对应关系

当前可以建立如下对应：

| 问题 | 结构空间离散 | 材料非线性位置 | Global equilibrium |
|---|---|---|---|
| 一维三材料杆 | 1D bar elements | element / material point | 杆轴向平衡 |
| 二维开孔板 | 2D continuum elements | continuum Gauss points | 二维连续体平衡 |
| 海上风机塔筒 | beam-column elements + fiber section | fiber material points | 梁柱结构整体平衡 |

因此，塔筒扩展的核心并不是：

$$
\text{新的 PGD 理论}
$$

而是：

$$
\boxed{ \text{continuum/bar spatial operator} \rightarrow \text{fiber beam-column spatial operator} }
$$

---

# 7. 纤维梁柱塔筒与 LATIN local/global 分裂具有天然对应关系

## 7.1 Local stage

塔筒每一个材料积分点可以表示为：

$$
(e,g,f).
$$

每个纤维材料点采用当前已经建立的一维粘塑性-损伤模型，其内部状态包括：

$$
\mathbf{s}_{\rm mat} = [\varepsilon_p,\alpha,\bar r,D].
$$

在 Local stage 中，可继续局部处理：

$$
\dot\varepsilon_p, \qquad \dot\alpha, \qquad \dot{\bar r}, \qquad \dot D,
$$

以及与其对应的：

$$
\sigma, \qquad \beta, \qquad R, \qquad Y.
$$

因此：

$$
\boxed{ \text{continuum Gauss point local stage} \rightarrow \text{fiber material point local stage} }
$$

这一迁移是自然的。

---

## 7.2 Global stage

塔筒 Global stage 的空间平衡路径为：

$$
\sigma_f \rightarrow \text{section force} \rightarrow \text{element internal force} \rightarrow \text{tower global equilibrium}.
$$

即：

$$
\sigma_f \rightarrow (N,M) \rightarrow \mathbf f_{\rm int}^{e} \rightarrow \mathbf F_{\rm int}.
$$

现有代码中已经建立：

```text
fem/viscoplastic_fiber_section.py
fem/viscoplastic_beam_column_2d.py
fem/viscoplastic_tower_system_2d.py
```

这些模块已经具备：

- fiber stress integration；
- section force assembly；
- beam element internal force；
- global tower assembly；
- Newton equilibrium；
- trial / commit / revert material-state management。

因此，当前塔筒 FOM 已经提供了构造 LATIN Global stage 所需的重要空间算子基础。

---

# 8. 塔筒 PGD 应首先围绕 fiber plastic-strain correction 建立

原论文并不是把所有变量都直接 PGD 化。

尤其需要注意：

$$
D, \quad \alpha, \quad \bar r
$$

等内部变量在原论文中仍然在 Gauss point 上局部更新，而不是直接作为 PGD separated variables。

因此，在塔筒第一版 LATIN-PGD 中，应继续遵循原论文的思路，优先对：

$$
\boxed{ \Delta\dot\varepsilon^p(x,t) }
$$

进行 PGD 分解。

塔筒离散形式可写为：

$$
\Delta\dot\varepsilon^p(e,g,f,t) = \sum_{j=1}^{m} \Delta\dot\lambda_j(t) \bar\varepsilon^p_j(e,g,f).
$$

而不是一开始就直接构造：

$$
D(x,t) = \sum_j D_j^t(t)D_j^x(x).
$$

因此，当前推荐保持：

$$
\boxed{ \text{PGD 主对象：plastic correction field} }
$$

同时保留：

$$
\boxed{ \text{hardening / damage：local material-point update} }
$$

这一原论文结构。

---

# 9. 当前 100 周期 SVD/HOSVD 工作仍然有重要价值

此前已完成：

$$
q(n,\tau,x)
$$

形式的张量化与 SVD/HOSVD 风格诊断。

本阶段调整路线以后，这部分工作的定位应重新表述为：

$$
\boxed{ \text{塔筒循环响应具有低秩可压缩性的离线证据} }
$$

而不是：

$$
\boxed{ \text{必须立即采用 }n-\tau-x\text{ PGD 的理由} }
$$

此前结果已经表明：

- displacement 高度低秩；
- stress 高度低秩；
- plastic strain 具有明显低秩结构；
- damage 具有明显低秩结构；
- 不可逆增量的复杂度高于 raw field；
- cycle、phase、space 方向均存在明显可压缩性。

因此，这些结果现在可以用于支撑：

> **原论文 $x-t$ PGD 在塔筒循环疲劳问题中具有应用潜力。**

也就是说，现有 SVD/HOSVD 工作应当保留，并作为后续方法迁移的可行性验证基础。

---

# 10. 为什么现阶段不应立即采用 n-τ-x PGD

如果现在直接从：

$$
x-t
$$

跳到：

$$
n-\tau-x,
$$

会同时引入多个新的理论问题：

1. $t\rightarrow(n,\tau)$ 后时间导数如何定义；
2. cycle index $n$ 是离散变量还是连续慢时间；
3. cycle boundary 的内部变量如何连续传递；
4. Add a pair 如何改写为 Add a triplet；
5. 原论文 time-update strategy 如何改写；
6. 原论文 Eq. (59) minimisation 如何变为三变量 alternating solver；
7. $\xi$ 的时空 norm 如何在三变量域中定义；
8. Gram-Schmidt 与 basis management 如何扩展；
9. 原论文已有数值验证基础与新方法之间将产生较大距离。

这样会使“塔筒结构扩展”和“时间多尺度方法创新”两个问题耦合在一起，难以判断数值问题究竟来源于：

- 塔筒空间离散；
- LATIN 结构迁移；
- 还是新的三变量 PGD。

因此，本阶段决定采用更加稳妥的渐进路线：

$$
\boxed{ \text{先验证原始 }x-t\text{ LATIN-PGD 可以用于塔筒} }
$$

再根据高周循环下的实际计算瓶颈决定是否进一步发展：

$$
\boxed{ x-t \rightarrow n-\tau-x }
$$

---

# 11. 当前收敛与 enrichment 策略

当前不重新设计原论文的基本 LATIN/PGD 判据。

## 11.1 LATIN convergence indicator

继续保留原论文：

$$
\boxed{ \xi = \frac{ \left\| \hat{\mathbf{s}}^p_{i+1/2} - \mathbf{s}^p_{i+1} \right\| }{ \left\| \hat{\mathbf{s}}^p_{i+1/2} \right\| + \left\| \mathbf{s}^p_{i+1} \right\| } }
$$

其物理意义为：

> local-stage solution 与 global-stage solution 之间的相对距离。

---

## 11.2 PGD saturation / enrichment indicator

继续采用原论文 Eq. (60)：

$$
\boxed{ \zeta = \frac{ \xi_i-\xi_{i+1} }{ \xi_i+\xi_{i+1} } }
$$

其中：

- $\zeta$ 较大：已有 basis 仍然能够有效降低 LATIN indicator；
- $\zeta$ 较小：已有 basis 趋于饱和，需要增加新的 PGD pair。

因此：

$$
\boxed{ \xi \rightarrow \zeta \rightarrow \text{Update / Enrich} }
$$

仍应作为塔筒版本的核心自适应机制。

---

# 12. 结合一维杆复现经验增强数值鲁棒性

在一维三材料杆复现过程中已经认识到，仅机械依赖单一 indicator 可能遇到：

- indicator stagnation；
- basis 持续增加但改进有限；
- damage/internal-variable 已趋稳但整体 indicator 波动；
- search direction 病态；
- 数值假收敛或不必要 enrichment。

因此，塔筒版本推荐采用：

$$
\boxed{ \text{论文主收敛准则} + \text{辅助数值稳定性检查} }
$$

需要特别强调：

> 这些辅助判据不是重新定义 LATIN 理论，而是工程实现中的 robustness control。

因此后续可继续保留我们在一维杆复现中形成的：

- LATIN relative indicator；
- saturation / enrichment control；
- stagnation detection；
- absolute tolerance；
- reduced residual check；
- damage/internal-variable stability check；
- search-direction regularisation；
- basis small-mode rejection。

---

# 13. 当前研究目标的正式定义

现阶段研究目标建议正式表述为：

> **在保持 Bhattacharyya et al. (2018) 原始 $x-t$ LATIN-PGD formulation 基本不变的前提下，将其从杆/连续体有限元问题扩展至由纤维梁柱单元组装的海上风机钢塔筒，并重点研究基于纤维材料积分点的 PGD 空间表示、LATIN global equilibrium operator、材料局部更新以及适用于塔筒循环疲劳问题的收敛和自适应增广策略。**

其核心数学迁移为：

$$
\boxed{ \text{原论文 continuum/bar spatial problem} \rightarrow \text{fiber beam-column tower spatial problem} }
$$

而保持：

$$
\boxed{ x-t\text{ PGD} }
$$

暂时不变。

---

# 14. 当前方法迁移的总体层级

可以将当前研究路线概括为：

```text
Bhattacharyya et al. (2018)
original x-t LATIN-PGD
        ↓
1D three-material bar reproduction
        ↓
convergence / enrichment / robustness experience
        ↓
2D fiber beam-column tower FOM
        ↓
fiber-level viscoplastic-damage material points
        ↓
100-cycle nonlinear reference solution
        ↓
offline low-rank feasibility evidence
        ↓
direct x-t LATIN-PGD migration to tower
```

当前的目标不是：

```text
invent a new PGD first
```

而是：

```text
make the original LATIN-PGD work on the tower first
```

---

# 15. 当前最重要的理论问题

在正式进入代码实现以前，首先需要回答：

$$
\boxed{ \text{原论文中的 } \bar\varepsilon^p(x) \text{ 在纤维梁柱塔筒中应如何严格定义？} }
$$

原论文 Eq. (47)：

$$
\Delta\varepsilon^p_{i+1} = \lambda^p(t) \bar\varepsilon^p(x)
$$

迁移到塔筒后，应明确：

1. $x$ 是否直接定义为 fiber material-point space；
2. $\bar\varepsilon^p(x)$ 的离散维数；
3. fiber spatial mode 与 section resultant 之间的映射；
4. fiber plastic-strain mode 如何进入 beam compatibility；
5. 如何由 fiber mode 构造对应的 stress correction；
6. 如何建立塔筒版 equilibrium operator；
7. 原论文 Eq. (50)–(53) 在纤维梁柱离散下如何重写。

---

# 16. 下一阶段工作原则

下一阶段遵循以下原则：

> **一次只解决一个数学问题，不直接同时修改理论与代码。**

因此下一步暂时：

- 不引入 $n-\tau-x$；
- 不重新设计 enrichment indicator；
- 不重新设计 damage PGD；
- 不直接编写新的三变量求解器；
- 不改变原论文 local/global 基本结构。

下一步正式从原论文 Eq. (47) 开始：

$$
\Delta\varepsilon^p_{i+1}(x,t) = \lambda^p(t) \bar\varepsilon^p(x),
$$

首先明确：

$$
\boxed{ \bar\varepsilon^p(x) }
$$

在纤维梁柱海上风机塔筒中的物理意义、离散定义及其与梁柱结构空间算子的对应关系。

---

# 17. 本阶段最终结论

1. 当前首要目标是**将原论文 LATIN-PGD 直接应用于海上风机塔筒**，而不是立即构造新的三变量 PGD。
2. 原论文已经从一维杆扩展到二维连续体，说明 LATIN-PGD 本身不受一维结构限制。
3. 海上风机塔筒可视为进一步更换空间离散形式，而不必立即改变 PGD 的 $x-t$ 分解结构。
4. 塔筒中的广义空间坐标可表示为：
$$
x\leftrightarrow(e,g,f).
$$
5. 原论文中的 continuum Gauss point 可自然对应到塔筒 fiber material point。
6. Local stage 可以继续在每个 fiber material point 上进行本构更新。
7. Global stage 的核心变化是建立：
$$
\sigma_f \rightarrow \text{section} \rightarrow \text{beam element} \rightarrow \text{tower equilibrium}.
$$
8. PGD 第一阶段仍应主要作用于 plastic-strain correction，而不是直接 PGD 分解 damage。
9. 已完成的 SVD/HOSVD 工作继续保留，但定位为 **offline low-rank feasibility evidence**。
10. LATIN convergence indicator $\xi$ 和 saturation indicator $\zeta$ 应继续继承原论文。
11. 一维杆复现得到的辅助收敛与鲁棒性判据可以用于增强塔筒算法，但不改变原 LATIN 理论。
12. 下一步应从原论文 Eq. (47) 出发，首先严格定义塔筒中的：
$$
\bar\varepsilon^p(x).
$$

---

# 18. 下一步

下一步正式开展：

> **原论文 Eq. (47) 在纤维梁柱海上风机塔筒中的空间模式定义。**

核心问题：

$$
\boxed{ \bar\varepsilon^p(x) \quad \longrightarrow \quad \bar\varepsilon^p(e,g,f) }
$$

并进一步明确它如何通过：

$$
\text{fiber} \rightarrow \text{section} \rightarrow \text{element} \rightarrow \text{global tower}
$$

进入 LATIN-PGD 的 Global stage。

---

# 参考依据

## 原论文

Bhattacharyya, M., Fau, A., Nackenhorst, U., Néron, D., & Ladevèze, P. (2018).  
**A LATIN-based model reduction approach for the simulation of cycling damage.**  
*Computational Mechanics, 62*, 725–743.  
DOI: 10.1007/s00466-017-1523-z.

## 当前项目关键代码

```text
material/viscoplastic_damage_1d.py

fem/viscoplastic_fiber_section.py
fem/viscoplastic_beam_column_2d.py
fem/viscoplastic_tower_system_2d.py

latin/local_stage.py
latin/global_stage.py
latin/search_directions.py
latin/pgd_basis.py
latin/pgd_time_update.py
latin/pgd_enrichment.py
latin/pgd_global_stage.py
latin/pgd_saturation.py
latin/pgd_solver.py
```

## 当前分支

```text
feature/offshore-wind-turbine-tower-fatigue
```

