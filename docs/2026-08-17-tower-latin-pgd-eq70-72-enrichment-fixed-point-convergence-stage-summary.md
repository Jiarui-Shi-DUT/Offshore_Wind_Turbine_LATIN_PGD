# Tower LATIN-PGD Eq. (70)–(72) enrichment fixed-point loop、normalization 与 convergence criterion 阶段总结

**日期：2026-08-17**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**当前研究路线：Bhattacharyya et al. 原论文 $x-t$ LATIN-PGD → 2D fiber beam-column offshore wind turbine tower**  
**阶段范围：在 Eq. (72) backward-Euler temporal update、$t_0$ initial treatment、BE temporal contractions 与 denominator conditioning 已完成的基础上，闭合 new separated pair 在 Eq. (70)–(72) 之间的 alternating fixed-point 主循环；明确 spatial/temporal half-step 的数据依赖；区分 fixed-point 内 scale normalization 与 fixed-point 后 Gram–Schmidt orthogonalisation；确定 fixed-point 应比较完整 separated mechanical correction，而非单独 spatial mode 或 temporal function；评估 current 1D `fixed_point_change` 的可迁移部分与不可直接迁移部分；给出 tower v1 的 primary fixed-point indicator、secondary residual-control diagnostic 及 tolerance 的数学与物理含义。**  
**上一阶段衔接：`2026-08-17-tower-latin-pgd-eq72-initial-time-be-contractions-conditioning-stage-summary.md`**  
**下一阶段：闭合 fixed-point 收敛后 Gram–Schmidt、existing temporal functions correction、new temporal function rescaling、insignificant-mode rejection、new-mode acceptance 与 existing-basis temporal re-optimisation 的完整先后顺序。**

---

# 1. 本阶段定位

上一阶段已经完成 Eq. (72) single-new-mode temporal half-step 的离散化，并明确了以下基础事实。

对于 $n \ge 1$，采用 right-endpoint backward Euler：

$$ \dot{\lambda}_n = \frac{\lambda_n-\lambda_{n-1}}{\Delta t_n} $$

定义：

$$ \vec{g}_n = \frac{\vec{p}}{\Delta t_n} - D_{H,n}\vec{s} $$

以及：

$$ \vec{b}_n = \frac{\vec{p}}{\Delta t_n}\lambda_{n-1} - \vec{\Delta}_n $$

则 one-step residual 为：

$$ \vec{r}_n = \vec{g}_n\lambda_n - \vec{b}_n $$

weighted objective 为：

$$ J_n = \frac{1}{2}\vec{r}_n^T M D_{H,n}^{-1}\vec{r}_n $$

single-mode scalar temporal update 为：

$$ \lambda_n = \frac{\vec{g}_n^T M D_{H,n}^{-1}\vec{b}_n}{\vec{g}_n^T M D_{H,n}^{-1}\vec{g}_n} $$

对于 virgin-state tower v1，初始 temporal amplitude 采用：

$$ \lambda_0 = 0 $$

但 initial rate 独立求解：

$$ \dot{\lambda}_0 = -\frac{\vec{p}^T M D_{H,0}^{-1}\vec{\Delta}_0}{\vec{p}^T M D_{H,0}^{-1}\vec{p}} $$

同时，Eq. (70)–(71) 中的关键 temporal contraction 已统一为 BE-consistent right-endpoint form：

$$ a_h = \sum_{n=1}^{N}\Delta t_n\dot{\lambda}_n\lambda_n $$

在 $\lambda_0=0$ 时：

$$ a_h = \frac{1}{2}\lambda_N^2 + \frac{1}{2}\sum_{n=1}^{N}(\lambda_n-\lambda_{n-1})^2 $$

因此对非平凡 discrete temporal mode：

$$ a_h > 0 $$

上一阶段还确认，Eq. (72) denominator：

$$ A_n = \vec{g}_n^T M D_{H,n}^{-1}\vec{g}_n $$

代表 current spatial candidate 在该 time step 对 residual 的 temporal controllability，而不是普通 stiffness singularity。tower v1 不采用 denominator clipping，而采用 detect → diagnose → solve original equation when admissible → verify → reject/restart if genuinely degenerate 的 safeguard 原则。

在这些基础上，本阶段需要回答的核心问题不再是某一个 temporal scalar equation，而是：

> Eq. (72) 求出 temporal function 后，如何返回 Eq. (70)–(71)，并把 new mode enrichment 组成一个完整、尺度一致、可判收敛的 alternating fixed-point algorithm？

---

# 2. 本阶段资料边界与理论地位

本阶段结论必须严格区分四类来源。

| 类别 | 本阶段含义 |
|---|---|
| 原论文明确内容 | Bhattacharyya et al. Eq. (61)–(72) 给出 new separated pair、fixed temporal 的 spatial Galerkin problem、fixed spatial 的 temporal residual minimisation，并在 Eq. (72) 后说明 spatial basis Gram–Schmidt 以及相应 temporal functions correction。 |
| 由原论文推导 | 将 Eq. (70)–(72) 写成 tower material-point 离散形式；分析 rank-one scaling freedom；推导 post-orthogonalisation temporal coefficient transformation；构造 scale-invariant complete-pair fixed-point indicator。 |
| current 1D implementation | `latin/pgd_enrichment.py` 中 deterministic spatial seed、inner loop、current `fixed_point_change`、inner orthogonalisation、final temporal recomputation 与 post-loop scalar line search。 |
| tower v1 engineering choice | fixed-point sweep 的精确执行顺序、每轮 convergence check 的位置、fixed-point 内只做 scale normalization、primary graph-norm indicator、secondary residual-control diagnostic、provisional tolerance 等。 |

需要特别强调：

> 原论文明确提出 fixed-point algorithm，但正文没有给出 inner fixed-point convergence formula。因此，本阶段给出的 $\chi_{\rm fp}$ 属于 tower v1 numerical specification，而不是原论文显式公式。

---

# 3. new separated pair 的基本未知量与物理含义

在一个 enrichment 中，新增 rank-one correction 写为：

$$ \Delta\dot{\vec{\varepsilon}}^p(t) = \dot{\lambda}(t)\vec{p} $$

以及：

$$ \Delta\vec{\sigma}'(t) = \lambda(t)\vec{s} $$

其中：

$$ \vec{p} = \vec{\bar\varepsilon}^p $$

为 material-point plastic-strain spatial mode，

$$ \vec{s} = \vec{\bar\sigma} $$

为与该 spatial source strain 对应的 equilibrated stress mode。

对于 tower fiber material-point discretisation：

$$ q=(e,g,f) $$

并将全部 material points flatten 为长度 $N_q$ 的空间向量。

reference equilibrium operator 已在前序阶段建立，因此：

$$ \vec{s} = C_0(\mathcal E_{\rm tower}-I)\vec{p} $$

其中：

$$ \mathcal E_{\rm tower} = H(H^TMC_0H)^{-1}H^TMC_0 $$

因此 $\vec{s}$ 不应作为 independent spatial unknown 存储或迭代，而是由 $\vec{p}$ 通过 reference equilibrium mapping 派生。

new pair 的 canonical data 为：

$$ \mathcal P = \{\vec{p},\vec{s},\lambda(t),\dot\lambda(t)\} $$

---

# 4. enrichment 内部保持固定的外层数据

一个 new-mode inner fixed point 是在同一次 LATIN global-stage enrichment 内进行的。

因此在 inner loop 中，下列数据必须视为固定：

$$ \vec{\bar\Delta}(t) $$

$$ D_H(t)=\operatorname{diag}(H_\sigma(t,q)) $$

$$ M $$

$$ H $$

$$ C_0 $$

以及当前已接受的 existing PGD basis。

这里的 shifted defect 来自 fixed-basis update 后的 Eq. (62) 结构：

$$ \vec{\bar\Delta} = D_H(\hat{\vec\sigma}-\vec\sigma^{\rm up})-(\hat{\dot{\vec\varepsilon}}^p-\dot{\vec\varepsilon}^{p,\rm up}) $$

在一个 new-pair alternating fixed point 内，不应重新执行 local stage，不应重新更新 $H_\sigma$，也不应重新定义 $\bar\Delta$。

这保证 inner fixed point 解的是同一个 frozen enrichment subproblem，而不是在迭代过程中不断改变目标方程。

---

# 5. 为什么不能以 zero temporal function 直接进入 Eq. (70)–(71)

如果初始化采用：

$$ \lambda^{(0)}(t)=0 $$

并直接进入 fixed-temporal spatial half-step，则 BE-consistent contractions 中：

$$ \sum_{n=1}^{N}\Delta t_n H_{\sigma,nq}(\lambda_n^{(0)})^2=0 $$

并且：

$$ a_h^{(0)}=\sum_{n=1}^{N}\Delta t_n\dot\lambda_n^{(0)}\lambda_n^{(0)}=0 $$

因此 Eq. (69)–(71) 中 effective spatial operator 退化。

同时：

$$ \bar\delta_q^{(0)}=\sum_{n=1}^{N}\Delta t_n\bar\Delta_{nq}\lambda_n^{(0)}=0 $$

Eq. (71) 还包含：

$$ \frac{1}{a_h^{(0)}} $$

因此 zero temporal initialization 无法生成第一轮有效 spatial update。

结论是：

**结论：new-mode fixed point 必须从 non-zero spatial seed 开始，而不是从 zero temporal mode 开始。**

---

# 6. tower v1 的初始化顺序

本阶段建议冻结以下初始化结构：

$$ \vec{p}^{(0)} \rightarrow \vec{s}^{(0)} \rightarrow \{\lambda^{(0)},\dot\lambda^{(0)}\} $$

首先构造一个 deterministic non-zero spatial seed：

$$ \vec{p}_{\rm seed}\ne0 $$

然后进行 scale normalization，得到：

$$ \vec{p}^{(0)} $$

再通过 reference equilibrium mapping 得到：

$$ \vec{s}^{(0)} = C_0(\mathcal E_{\rm tower}-I)\vec{p}^{(0)} $$

随后第一次执行 Eq. (72) temporal solve：

$$ (\vec{p}^{(0)},\vec{s}^{(0)}) \rightarrow \lambda^{(0)},\dot\lambda^{(0)} $$

至此才得到第一个完整 mutually consistent pair：

$$ \mathcal P^{(0)} = \{\vec{p}^{(0)},\vec{s}^{(0)},\lambda^{(0)},\dot\lambda^{(0)}\} $$

只有形成完整 pair 后，才进入正式 spatial ↔ temporal alternating fixed point。

---

# 7. deterministic spatial seed 的当前工程参考

current 1D `latin/pgd_enrichment.py` 的 `_initial_spatial_function()` 采用 weighted residual 最大时间截面构造 seed。

其思路为：

1. 对每一个 time row 计算 residual energy；
2. 找出 residual energy 最大的 time step；
3. 取该时刻 residual 的负方向作为 spatial seed；
4. 再做 spatial normalization。

对应的 1D weighted row energy 结构为：

$$ e_n = \sum_q \frac{V_q}{H_{\sigma,nq}}\bar\Delta_{nq}^2 $$

seed time index 为：

$$ n_* = \operatorname*{arg\,max}_n e_n $$

候选方向可写为：

$$ \vec{p}_{\rm seed}\propto-\vec{\bar\Delta}_{n_*} $$

这一策略具有 residual-driven 和 deterministic 两个优点，但其 exact tower version 尚未在本阶段最终冻结。

因此目前应表述为：

> current 1D seed construction 是 tower v1 的优先初始化候选方案，而不是 Bhattacharyya et al. 原论文明确规定。

---

# 8. 完整 fixed-point sweep 的定义

假设第 $k$ 轮开始时已经拥有一个完整 pair：

$$ \mathcal P^{(k)}=\{\vec{p}^{(k)},\vec{s}^{(k)},\lambda^{(k)},\dot\lambda^{(k)}\} $$

则 tower v1 将一轮完整 fixed-point sweep 定义为：

**完整 sweep：`Spatial half-step → Temporal half-step → Complete-pair convergence check`。**

这一定义非常重要。

convergence check 不应放在 spatial half-step 之后，因为此时新的 spatial mode 尚未对应新的 temporal function；只有重新执行 Eq. (72) 后，才能形成新的 mutually consistent pair。

---

# 9. Spatial half-step：固定 temporal function

第 $k$ 轮 spatial half-step 固定：

$$ \lambda^{(k)}(t) $$

以及：

$$ \dot\lambda^{(k)}(t) $$

首先计算 BE-consistent contraction：

$$ a_h^{(k)}=\sum_{n=1}^{N}\Delta t_n\dot\lambda_n^{(k)}\lambda_n^{(k)} $$

对于每一个 tower material point $q$，定义：

$$ \mathcal A_q^{(k)}=\sum_{n=1}^{N}\Delta t_n H_{\sigma,nq}(\lambda_n^{(k)})^2 $$

以及：

$$ \bar\delta_q^{(k)}=\sum_{n=1}^{N}\Delta t_n\bar\Delta_{nq}\lambda_n^{(k)} $$

于是 Eq. (69) 的 effective operator 为：

$$ (W_q^{(k)})^{-1}=\mathcal A_q^{(k)}+\frac{a_h^{(k)}}{E_{0,q}} $$

若 tower v1 first stage 使用 uniform scalar reference steel modulus，则：

$$ E_{0,q}=E_0 $$

定义：

$$ D_W^{(k)}=\operatorname{diag}(W_q^{(k)}) $$

---

# 10. Eq. (70) 的 tower spatial FE solve

固定 temporal function 后，求解：

$$ H^TMD_W^{(k)}H\vec{\bar{\tilde U}}^{(k+1)}=-H^TMD_W^{(k)}\vec{\bar\delta}^{(k)} $$

得到 compatible strain mode：

$$ \vec{\bar{\tilde\varepsilon}}^{(k+1)}=H\vec{\bar{\tilde U}}^{(k+1)} $$

恢复 temporary equilibrated stress mode：

$$ \vec{\bar\sigma}_{\rm raw}^{(k+1)}=D_W^{(k)}\left(\vec{\bar{\tilde\varepsilon}}^{(k+1)}+\vec{\bar\delta}^{(k)}\right) $$

根据 Eq. (71)，得到 raw plastic-strain spatial mode：

$$ \vec{p}_{\rm raw}^{(k+1)}=\frac{1}{a_h^{(k)}}\vec{\bar{\tilde\varepsilon}}^{(k+1)}-C_0^{-1}\vec{\bar\sigma}_{\rm raw}^{(k+1)} $$

这里的 $W^{(k)}$ 只是在 fixed-temporal spatial half-step 中出现的 effective enrichment operator。

它不是最终 accepted PGD basis 的 constitutive law。

一旦 Eq. (71) 恢复出 canonical plastic mode，accepted mode 的 stress spatial field 仍应通过 reference equilibrium mapping：

$$ \vec{s}=C_0(\mathcal E_{\rm tower}-I)\vec{p} $$

建立。

---

# 11. fixed-point 内的 scale normalization

rank-one representation 存在 reciprocal scaling freedom。

对于任意非零常数 $c$：

$$ \vec{p}^*=c\vec{p} $$

以及：

$$ \lambda^*(t)=\frac{\lambda(t)}{c} $$

则：

$$ \lambda^*(t)\vec{p}^*=\lambda(t)\vec{p} $$

同理：

$$ \dot\lambda^*(t)\vec{p}^*=\dot\lambda(t)\vec{p} $$

而由于 stress spatial mode 对 plastic mode 的 mapping 是线性的：

$$ \vec{s}^*=c\vec{s} $$

所以：

$$ \lambda^*(t)\vec{s}^*=\lambda(t)\vec{s} $$

因此 spatial factor 与 temporal factor 本身没有唯一尺度。

若不固定这个 scaling gauge，可能出现 spatial vector 越来越小而 temporal amplitude 越来越大，或者反过来，虽然 separated physical field 几乎不变，却使 numerical fixed-point comparison 失去意义。

---

# 12. tower v1 的 spatial normalization metric

建议在每次 Eq. (71) 得到 raw spatial mode 后，采用 material-point volume metric：

$$ \|\vec{p}\|_M=\sqrt{\vec{p}^TM\vec{p}} $$

定义：

$$ c_{k+1}=\sqrt{(\vec{p}_{\rm raw}^{(k+1)})^TM\vec{p}_{\rm raw}^{(k+1)}} $$

若：

$$ c_{k+1} $$

低于 minimum spatial norm，则应视为 current candidate spatial direction linearly dependent、numerically negligible 或 degenerate，而不是继续强行归一化。

若 $c_{k+1}$ 合法，则：

$$ \vec{p}^{(k+1)}=\frac{\vec{p}_{\rm raw}^{(k+1)}}{c_{k+1}} $$

并重新通过 canonical equilibrium mapping 得到：

$$ \vec{s}^{(k+1)}=C_0(\mathcal E_{\rm tower}-I)\vec{p}^{(k+1)} $$

从而固定：

$$ \|\vec{p}^{(k+1)}\|_M=1 $$

选择 $M$-norm 的原因是：

- $M$ 是 tower fiber material-point 的自然 spatial integration metric；
- 不依赖当前 time function；
- 不依赖当前 $H_\sigma(t)$；
- 与 element → Gauss point → fiber → material point 的离散结构一致；
- 可避免将 transient search-direction scaling 混入 spatial gauge definition。

---

# 13. normalization 后是否需要立即 reciprocal rescale 旧 temporal function

本阶段结论为：在 tower v1 的完整 sweep 定义下，通常不需要。

原因是 normalization 发生在 spatial half-step 末尾，随后立即重新执行 Eq. (72)：

$$ (\vec{p}^{(k+1)},\vec{s}^{(k+1)}) \rightarrow \lambda^{(k+1)},\dot\lambda^{(k+1)} $$

因此旧的：

$$ \lambda^{(k)} $$

不会与新的：

$$ \vec{p}^{(k+1)} $$

共同作为 complete-pair convergence comparison 的对象。

换言之，本轮流程不是：

不采用 `p_raw^(k+1) → normalize → 保持旧 λ^(k) 作为新 pair` 这一流程。

而是：

采用 `p_raw^(k+1) → normalize → s^(k+1) → Eq. (72) re-solve`。

如果某个中间 diagnostic 需要严格保持 normalization 前后的 instantaneous rank-one product 不变，则可暂时 reciprocal rescale temporal factor；但这不是 tower v1 主 fixed-point loop 所必需的操作。

---

# 14. Temporal half-step：固定新的 spatial pair

完成 normalization 后，固定：

$$ \vec{p}^{(k+1)} $$

以及：

$$ \vec{s}^{(k+1)} $$

重新执行 Eq. (72)。

对于 $t_0$：

$$ \lambda_0^{(k+1)}=0 $$

并独立求：

$$ \dot\lambda_0^{(k+1)}=-\frac{(\vec{p}^{(k+1)})^TMD_{H,0}^{-1}\vec{\bar\Delta}_0}{(\vec{p}^{(k+1)})^TMD_{H,0}^{-1}\vec{p}^{(k+1)}} $$

对于 $n\ge1$：

$$ \vec{g}_n^{(k+1)}=\frac{\vec{p}^{(k+1)}}{\Delta t_n}-D_{H,n}\vec{s}^{(k+1)} $$

$$ \vec{b}_n^{(k+1)}=\frac{\vec{p}^{(k+1)}}{\Delta t_n}\lambda_{n-1}^{(k+1)}-\vec{\bar\Delta}_n $$

然后：

$$ \lambda_n^{(k+1)}=\frac{(\vec{g}_n^{(k+1)})^TMD_{H,n}^{-1}\vec{b}_n^{(k+1)}}{(\vec{g}_n^{(k+1)})^TMD_{H,n}^{-1}\vec{g}_n^{(k+1)}} $$

并恢复：

$$ \dot\lambda_n^{(k+1)}=\frac{\lambda_n^{(k+1)}-\lambda_{n-1}^{(k+1)}}{\Delta t_n} $$

至此得到新的完整 pair：

$$ \mathcal P^{(k+1)}=\{\vec{p}^{(k+1)},\vec{s}^{(k+1)},\lambda^{(k+1)},\dot\lambda^{(k+1)}\} $$

只有在这一时刻，才进入 fixed-point convergence check。

---

# 15. 为什么 convergence check 必须放在 Eq. (72) 之后

如果在 spatial half-step 之后就比较，则 candidate 实际为：

$$ \{\vec{p}^{(k+1)},\vec{s}^{(k+1)},\lambda^{(k)},\dot\lambda^{(k)}\} $$

这一组合的 temporal function 并不是在新 spatial mode 上求得。

因此它不是 spatial and temporal mutually consistent fixed-point iterate。

真正的 alternating fixed-point map 应写为：

$$ \mathcal P^{(k)} \xrightarrow{\mathrm{Eq.70-71}} \{\vec{p}^{(k+1)},\vec{s}^{(k+1)}\} \xrightarrow{\mathrm{Eq.72}} \mathcal P^{(k+1)} $$

所以 convergence 应比较：

$$ \mathcal P^{(k+1)} $$

与：

$$ \mathcal P^{(k)} $$

而不是比较 half-updated intermediate state。

---

# 16. current 1D `pgd_enrichment.py` 的实际 inner-loop 顺序

current 1D enrichment 的重要逻辑是：

```text
current spatial function
    -> equilibrium stress
    -> temporal update
    -> spatial update
    -> equilibrium stress for updated spatial function
    -> build candidate mode using updated spatial + preceding temporal function
    -> fixed_point_change
```

因此 current 1D 在 inner loop 中计算 `fixed_point_change` 时，使用的是：

也就是说，current 1D 的该次比较对应 `new spatial + temporal function solved from preceding spatial`。

只有退出 inner loop 后，代码才针对 final spatial vector 再执行一次 temporal update，使最终 returned pair mutually consistent。

这一结构在现有 1D 三材料杆复现中已经获得了实际验证，因此不能简单称为“错误”。

但 tower v1 当前正在重新按照 Eq. (70)–(72) 建立 paper-derived spatial operator，因此可以采用更清晰的算法语义：

**tower v1 规定：每一次 convergence check 都在完成 Eq. (72) temporal re-solve 以后进行。**

这样 fixed-point history 中的每一个记录值都对应一个完整 pair。

---

# 17. fixed-point 内 normalization 与 fixed-point 后 Gram–Schmidt 必须区分

这是本阶段最重要的结构性结论之一。

fixed-point 内的 normalization 解决的是 rank-one factor scaling non-uniqueness：

$$ \vec{p}\leftrightarrow\lambda $$

其目的只是固定 scaling gauge。

fixed-point 后的 Gram–Schmidt 解决的是 new spatial direction 与 existing PGD basis 之间的 linear dependence 和 basis conditioning。

两者的数学目的不同，因此不应混在同一步完成。

建议 tower v1 采用：

**inner fixed point：只做 scale normalization。**

以及：

**fixed-point convergence 后：再对 existing basis 做 Gram–Schmidt。**

原论文在 Eq. (72) 后讨论 spatial basis orthonormalisation，并同时调整 temporal functions，这与“先求收敛 raw pair，再做 basis transformation”的结构一致。

---

# 18. current 1D 与建议 tower v1 在 orthogonalisation 位置上的差异

current 1D `_solve_spatial_function()` 在每一次 spatial update 后就调用 `_orthogonalize_and_normalise()`，即：

current 1D 顺序为 `inner spatial solve → orthogonalise against existing basis → normalize`。

建议 tower v1 则为：

建议 tower v1 的 inner 顺序为 `inner spatial solve → scale normalize only`。

当 fixed point 收敛后才：

随后在 fixed point 收敛后执行 `Gram–Schmidt against existing basis`。

这种调整的目的不是否定 current 1D，而是使 tower v1 的 Eq. (70)–(72) alternating subproblem 与 post-enrichment basis maintenance 分离得更清楚。

---

# 19. post-fixed-point Gram–Schmidt 的一般投影形式

设 existing plastic spatial basis 为：

$$ P=[\vec{p}_1,\ldots,\vec{p}_m] $$

raw converged new spatial mode 为：

$$ \vec{p}_* $$

若 existing basis 尚未严格满足 $M$-orthonormality，则首先定义 Gram matrix：

$$ G=P^TMP $$

projection coefficients 为：

$$ \vec{a}=G^{-1}P^TM\vec{p}_* $$

数值实现时不应显式求 $G^{-1}$，而应通过 stable linear solve 或 least-squares 求：

$$ G\vec{a}=P^TM\vec{p}_* $$

去除已有 basis 分量：

$$ \vec{p}_\perp=\vec{p}_*-P\vec{a} $$

定义：

$$ c=\|\vec{p}_\perp\|_M $$

若：

$$ c $$

过小，则 new spatial direction 基本落在 existing basis span 内，应视为 linearly dependent 或 insignificant candidate。

若 $c$ 合法，则新 basis mode 为：

$$ \vec{p}_{m+1}=\frac{\vec{p}_\perp}{c} $$

当 existing basis 已严格满足：

$$ P^TMP=I $$

则简化为：

$$ \vec{a}=P^TM\vec{p}_* $$

---

# 20. Gram–Schmidt 后 temporal functions 为什么必须同步变换

由：

$$ \vec{p}_*=P\vec{a}+c\vec{p}_{m+1} $$

原 new plastic correction 为：

$$ \lambda_*\vec{p}_* $$

因此：

$$ \lambda_*\vec{p}_*=P(\vec{a}\lambda_*)+\vec{p}_{m+1}(c\lambda_*) $$

为了在 basis transformation 后保持完整 separated field 完全不变，existing temporal amplitudes 必须更新：

$$ \lambda_j^{\rm new}=\lambda_j^{\rm old}+a_j\lambda_* $$

new mode temporal amplitude 为：

$$ \lambda_{m+1}=c\lambda_* $$

对于 rate 同理：

$$ \dot\lambda_j^{\rm new}=\dot\lambda_j^{\rm old}+a_j\dot\lambda_* $$

以及：

$$ \dot\lambda_{m+1}=c\dot\lambda_* $$

由于 stress mode mapping：

$$ \mathcal L = C_0(\mathcal E_{\rm tower}-I) $$

是线性的，

$$ \vec{s}_*=\mathcal L\vec{p}_*=S\vec{a}+c\vec{s}_{m+1} $$

因此相同的 temporal coefficient transformation 也保持 stress correction：

$$ \lambda_*\vec{s}_* $$

不变。

这给出了原论文所述“spatial basis orthonormalisation 后 former time functions 也需要修改”的具体 tower algebraic interpretation。

本阶段只闭合这一 basis transformation identity；其与 subsequent existing-basis temporal re-optimisation 的先后顺序仍留到下一阶段完整确定。

---

# 21. 为什么 fixed-point 不应单独比较 spatial mode

一种直观 criterion 是：

$$ \frac{\|\vec{p}^{(k+1)}-\vec{p}^{(k)}\|}{\|\vec{p}^{(k+1)}\|} $$

但它不适合作为 primary fixed-point criterion。

原因包括：

- spatial factor 本身存在 reciprocal scaling ambiguity；
- normalization convention 会直接影响该指标；
- sign flip 可能造成 apparent large change；
- spatial mode 接近不代表 temporal factor 已经同步稳定；
- fixed point 的真正对象是 space-time rank-one correction，而不是单独空间向量。

即使人为固定：

$$ \|\vec{p}\|_M=1 $$

单独比较 spatial mode 仍不能保证：

$$ \lambda(t) $$

已与之达到 mutual consistency。

---

# 22. 为什么 fixed-point 不应单独比较 temporal function

同理，若采用：

$$ \frac{\|\lambda^{(k+1)}-\lambda^{(k)}\|}{\|\lambda^{(k+1)}\|} $$

也存在问题。

首先 temporal factor 本身同样受 reciprocal scaling 影响。

其次 temporal amplitude 稳定并不保证 spatial mode 已经稳定。

第三，Eq. (72) sequential BE temporal update 是在 fixed spatial pair 下求解，因此其变化量只反映 temporal half-step，不代表完整 alternating map 已经收敛。

因此 fixed-point convergence 应直接作用在 separated physical correction 上。

---

# 23. complete separated mechanical correction

第 $k$ 次完整 pair 对应的 plastic-rate correction field 为：

$$ \Delta\dot{\vec\varepsilon}^{p,(k)}_n=\dot\lambda_n^{(k)}\vec{p}^{(k)} $$

stress correction field 为：

$$ \Delta\vec\sigma'^{(k)}_n=\lambda_n^{(k)}\vec{s}^{(k)} $$

定义 complete mechanical pair：

$$ z^{(k)}=\left(\Delta\dot{\vec\varepsilon}^{p,(k)},\Delta\vec\sigma'^{(k)}\right) $$

这一对象具有非常关键的性质。

对于 reciprocal scaling：

$$ \vec{p}^*=c\vec{p} $$

$$ \vec{s}^*=c\vec{s} $$

$$ \lambda^*=\lambda/c $$

$$ \dot\lambda^*=\dot\lambda/c $$

有：

$$ \Delta\dot{\vec\varepsilon}^{p,*}=\Delta\dot{\vec\varepsilon}^{p} $$

以及：

$$ \Delta\vec\sigma'^*=\Delta\vec\sigma' $$

因此 $z$ 对 PGD factor scaling 完全不变。

这正是 fixed-point primary convergence variable 应具有的性质。

---

# 24. tower v1 primary fixed-point graph norm

结合 LATIN search-direction metric，对 complete mechanical pair 定义：

$$ \|z\|_{\rm fp,h}^2=\sum_{n=1}^{N}\Delta t_n\left[(\Delta\dot{\vec\varepsilon}^{p}_n)^TMD_{H,n}^{-1}\Delta\dot{\vec\varepsilon}^{p}_n+(\Delta\vec\sigma'_n)^TMD_{H,n}\Delta\vec\sigma'_n\right] $$

该 norm 具有以下性质。

第一，只要：

$$ M>0 $$

以及：

$$ D_{H,n}>0 $$

就具有 positive-definite mechanical metric。

第二，plastic-rate correction 与 stress correction 分别计量，不会在 norm 内彼此 cancellation。

第三，其 weighting 与 LATIN search-direction 结构一致。

第四，使用与 Eq. (70)–(72) 当前离散策略一致的 right-endpoint time-slab summation：

$$ \sum_{n=1}^{N}\Delta t_n(\cdot)_n $$

而不在 inner fixed-point 中重新引入 trapezoidal time integration。

第五，它直接作用于 complete separated physical correction，因此对 reciprocal scaling 不敏感。

---

# 25. primary fixed-point indicator

定义 consecutive complete-pair difference：

$$ \delta z^{(k+1)}=z^{(k+1)}-z^{(k)} $$

建议 tower v1 使用 symmetric relative change：

$$ \chi_{\rm fp}^{(k+1)}=\frac{\|z^{(k+1)}-z^{(k)}\|_{\rm fp,h}}{\|z^{(k+1)}\|_{\rm fp,h}+\|z^{(k)}\|_{\rm fp,h}} $$

在 denominator 非零时，由 triangle inequality：

$$ 0\le\chi_{\rm fp}^{(k+1)}\le1 $$

fixed-point convergence 定义为：

$$ \boxed{\chi_{\rm fp}^{(k+1)}\le\varepsilon_{\rm fp}} $$

与仅用 current iterate norm 作 denominator 相比，symmetric form 不人为选择其中一轮作为绝对 reference，并在 amplitude 较大或较小时保持更平衡的 relative interpretation。

如果：

$$ \|z^{(k+1)}\|_{\rm fp,h}+\|z^{(k)}\|_{\rm fp,h} $$

本身已经低于 candidate significance threshold，则不应通过人为 denominator floor 把它解释为“fixed point 收敛”，而应进入 negligible-mode diagnostic。

---

# 26. $\chi_{\rm fp}$ 的数学含义

$\chi_{\rm fp}$ 衡量的是：

> 在相同 frozen shifted defect 与 search directions 下，连续两次完整 spatial–temporal alternating sweeps 所产生的 rank-one mechanical correction 在 LATIN metric 中还变化多少。

因此：

$$ \chi_{\rm fp}\rightarrow0 $$

意味着 alternating map 已接近 fixed point：

$$ \mathcal F(\mathcal P^*)\approx\mathcal P^* $$

这里 $\mathcal F$ 表示：

$$ \mathrm{Eq.70-71}\rightarrow\mathrm{normalization}\rightarrow\mathrm{Eq.72} $$

构成的一次完整 sweep。

需要注意，$\chi_{\rm fp}$ 不是对 exact FOM solution 的 error estimate，而只是 inner nonlinear alternating subproblem 的 iteration-change measure。

---

# 27. $\chi_{\rm fp}$ 的物理含义

在 material-point level，new mode 所做的物理修正主要由两部分构成：

$$ \Delta\dot{\vec\varepsilon}^p $$

代表新增的 plastic flow-rate pattern；

$$ \Delta\vec\sigma' $$

代表与该 spatial mode 对应的 equilibrated stress correction。

因此：

$$ \chi_{\rm fp}\ll1 $$

表示：

> 当前 new rank-one mode 所代表的 plastic-flow correction 与 equilibrated stress correction 在连续两轮 space-time alternating iteration 中已经基本不再改变。

这是一种 numerical fixed-point stability，而不是材料疲劳演化误差、结构响应误差或 total LATIN convergence error。

---

# 28. current 1D `fixed_point_change` 的核心思想

current 1D 代码不是单独比较 spatial mode，也不是单独比较 temporal amplitude。

它先构造 mode-controlled mechanical residual contribution：

$$ \vec{c}=\Delta\dot{\vec\varepsilon}^p-D_H\Delta\vec\sigma' $$

然后比较 consecutive correction：

$$ \vec{c}^{(k+1)}-\vec{c}^{(k)} $$

并采用 $H_\sigma^{-1}$ weighted space-time norm。

其 current relative structure大致为：

$$ \chi_{\rm 1D}=\frac{\|\vec{c}^{(k+1)}-\vec{c}^{(k)}\|_{H^{-1}}}{\|\vec{c}^{(k+1)}\|_{H^{-1}}} $$

这一设计的最重要优点是：

$$ \vec{c} $$

本身由完整 separated pair 构成，因此避免直接比较 $\vec{p}$ 或 $\lambda$ 的 reciprocal scaling ambiguity。

这一思想应继承到 tower。

---

# 29. 为什么 current 1D `fixed_point_change` 不建议原样迁移

本阶段结论是：

**结论：核心思想继承，但 exact formula 与 evaluation timing 不原样复制。**

原因一：current 1D spatial half-step 是 weighted least-squares spatial solve，而 tower v1 已决定使用 paper-derived Eq. (70)–(71) $W$-based FE solve。

原因二：current 1D inner-loop `fixed_point_change` 的 candidate pair 在 evaluation 时尚未针对最新 spatial mode 重做 temporal update；tower v1 改为 Eq. (72) 后比较 mutually consistent complete pairs。

原因三：current 1D weighted norm 使用 `np.trapz` 进行 time integration；tower Eq. (70)–(72) 当前已采用 right-endpoint BE-consistent contraction，因此 inner fixed-point metric 建议使用同一 time-slab summation。

原因四：current 1D 先把 plastic-rate correction 和 stress correction 合并成：

$$ \Delta\dot{\vec\varepsilon}^p-D_H\Delta\vec\sigma' $$

其中存在 cancellation 可能，而 tower v1 primary graph norm 将两部分分别计量。

---

# 30. residual-control combination 仍然具有重要诊断价值

虽然不再建议把：

$$ \vec{c}=\Delta\dot{\vec\varepsilon}^p-D_H\Delta\vec\sigma' $$

作为 tower v1 唯一 primary convergence variable，但它仍然具有直接的 Eq. (61) residual意义。

new-pair residual 可写为：

$$ \vec{r}=\vec{c}+\vec{\bar\Delta} $$

由于 $\vec{\bar\Delta}$ 在 inner fixed point 中固定：

$$ \vec{r}^{(k+1)}-\vec{r}^{(k)}=\vec{c}^{(k+1)}-\vec{c}^{(k)} $$

因此 $\vec{c}$ 的稳定程度直接表示 new pair 对 LATIN residual 的 control action 是否稳定。

---

# 31. secondary residual-control diagnostic

定义：

$$ \|\vec{c}\|_{H^{-1},h}^2=\sum_{n=1}^{N}\Delta t_n(\vec{c}_n)^TMD_{H,n}^{-1}\vec{c}_n $$

可进一步记录 symmetric relative diagnostic：

$$ \chi_c^{(k+1)}=\frac{\|\vec{c}^{(k+1)}-\vec{c}^{(k)}\|_{H^{-1},h}}{\|\vec{c}^{(k+1)}\|_{H^{-1},h}+\|\vec{c}^{(k)}\|_{H^{-1},h}} $$

于是两个指标的职责明确区分：

$\chi_{\rm fp}$ 回答 complete mechanical pair 是否稳定。

$\chi_c$ 回答 Eq. (61) residual-control action 是否稳定。

由于 $\chi_c$ 中先形成 plastic-rate 与 stress 的差，再做 norm，它可能受到 cancellation 影响，因此更适合作为 secondary diagnostic，而不是唯一停止判据。

---

# 32. fixed-point tolerance $\varepsilon_{\rm fp}$ 的准确含义

$\varepsilon_{\rm fp}$ 是 inner rank-one alternating solver 的 numerical convergence tolerance。

它不是：

- total LATIN convergence tolerance；
- PGD saturation tolerance；
- full-order model error tolerance；
- fatigue life error tolerance；
- material constitutive integration tolerance；
- Eq. (72) one-step scalar residual tolerance；
- new mode acceptance tolerance。

它只表示：

> consecutive complete spatial–temporal sweeps 所得到的 new separated mechanical correction，在 selected LATIN-weighted metric 下相对变化已经足够小。

若：

$$ \varepsilon_{\rm fp}=10^{-6} $$

其准确含义是：

> 连续两次完整 fixed-point sweeps 的 complete mechanical correction difference 相对于两轮 correction magnitude 的对称相对量已降至约 $10^{-6}$ 量级。

这不是“tower response 已达到 $10^{-6}$ 的物理精度”。

---

# 33. provisional tolerance 与 current 1D reference

current 1D `enrich_pgd_basis_once()` 默认：

$$ \varepsilon_{\rm fp}^{\rm 1D}=10^{-6} $$

并设置 maximum fixed-point iterations。

因此 tower v1 可以首先把：

$$ \varepsilon_{\rm fp}=10^{-6} $$

作为 reference starting value，以便 bar/tower numerical behaviour 具有初始可比性。

但目前不能把该数值写成 universal tolerance，更不能说它由原论文规定。

后续应通过以下数据进行 calibration：

- fixed-point history decay；
- residual reduction；
- accepted mode count；
- outer LATIN convergence；
- tower FOM comparison；
- runtime sensitivity；
- normalization scale；
- time-step size；
- $H_\sigma$ scaling；
- near-degenerate controllability events。

---

# 34. fixed-point convergence 与 new-mode acceptance 必须分离

即使：

$$ \chi_{\rm fp}\le\varepsilon_{\rm fp} $$

也不能直接说明 new mode 有价值。

可能存在一个 candidate 很快稳定，但它的 magnitude 几乎为零，或者对 residual 几乎没有 reduction。

因此必须区分：

**fixed-point convergence**

与：

**enrichment usefulness / mode acceptance**

前者回答：

> alternating iteration 对当前 rank-one candidate 是否稳定？

后者回答：

> 这个稳定 candidate 是否值得加入 PGD basis？

---

# 35. negligible-mode case 不能伪装成 fixed-point convergence

如果：

$$ \|z^{(k)}\|_{\rm fp,h}\approx0 $$

以及：

$$ \|z^{(k+1)}\|_{\rm fp,h}\approx0 $$

则 symmetric denominator 也接近零。

此时不能简单采用：

$$ \max(\|z^{(k+1)}\|+\|z^{(k)}\|,\varepsilon) $$

然后得到一个很小的 $\chi_{\rm fp}$，再宣布“converged”。

更合理的语义是：

**该情况应分类为：candidate pair numerically negligible。**

应进入 significance / acceptance diagnostic。

这与上一阶段 Eq. (72) denominator safeguard 的哲学一致：

> 不用任意 floor 修改原问题的数值含义，而是检测退化并明确分类。

---

# 36. Eq. (72) controllability diagnostics 在 fixed-point loop 中的位置

上一阶段已经定义：

$$ \eta_n=\frac{\|\vec{g}_n\|_{W_n}}{\|\vec{p}/\Delta t_n\|_{W_n}+\|D_{H,n}\vec{s}\|_{W_n}} $$

其中：

$$ \|\vec{x}\|_{W_n}=\sqrt{\vec{x}^TMD_{H,n}^{-1}\vec{x}} $$

$\eta_n$ 检测：

> temporal control column $\vec{g}_n$ 是否由两项强烈 cancellation 导致 near-degenerate。

同时：

$$ \rho_n=\frac{|\vec{g}_n^TW_n\vec{b}_n|}{\sqrt{(\vec{g}_n^TW_n\vec{g}_n)(\vec{b}_n^TW_n\vec{b}_n)}} $$

衡量 scalar temporal degree of freedom 对当前 residual 的 alignment / effectiveness。

这些指标不替代 $\chi_{\rm fp}$。

它们应在每次 temporal half-step 中作为 conditioning diagnostics 记录：

$\eta_n$ 与 $\rho_n$ 用于 temporal subproblem conditioning / effectiveness diagnostics。

而：

$\chi_{\rm fp}$ 用于 complete alternating fixed-point convergence。

---

# 37. fixed-point convergence 后的 residual reduction 仍需单独检查

设 enrichment 前 residual 为：

$$ \vec{r}_{\rm before}=\vec{\bar\Delta} $$

符号具体取决于代码采用 residual 或 forcing convention，但必须在实现中保持一致。

new mode correction 为：

$$ \vec{c}_*=\Delta\dot{\vec\varepsilon}^p_*-D_H\Delta\vec\sigma'_* $$

则 enrichment 后 residual 为：

$$ \vec{r}_{\rm after}=\vec{r}_{\rm before}+\vec{c}_* $$

至少应要求：

$$ \|\vec{r}_{\rm after}\|_{H^{-1},h}<\|\vec{r}_{\rm before}\|_{H^{-1},h} $$

才能说明该 candidate 对当前 residual 具有正向作用。

但是 exact acceptance threshold、是否保留 current 1D post-loop scalar line search，以及它与 paper-style temporal re-optimisation 的关系，本阶段尚未最终闭合，留待下一阶段处理。

---

# 38. current 1D post-loop scalar line search 的理论边界

current 1D 在 final spatial mode 上重新计算 temporal function 后，还对 complete separated correction 做一次 scalar line search。

其目的在于：

> current sequential backward-Euler temporal update 并不是 whole-time residual 的 coupled global minimizer，因此使用一个 global scalar coefficient保证 accepted correction 不增加 weighted residual norm。

这一 scalar line search 是 current 1D engineering safeguard，不是 Bhattacharyya et al. 原论文 Eq. (72) 后明确给出的步骤。

因此 tower v1 当前阶段不自动继承该步骤。

需要在下一阶段与以下操作一起重新排序和判断：

- fixed-point converged pair；
- Gram–Schmidt；
- existing temporal function correction；
- new temporal function rescaling；
- possible temporal re-optimisation；
- residual reduction test；
- insignificant-mode rejection。

---

# 39. outer LATIN convergence 与 inner fixed-point convergence 不能混淆

current project 已建立 outer LATIN convergence indicator：

$$ \xi $$

以及 basis adequacy / saturation-related indicator：

$$ \zeta $$

本阶段新增：

$$ \chi_{\rm fp} $$

三者作用层级不同。

$\xi$ 表示 local-stage / global-stage state distance。

$\zeta$ 表示 existing reduced basis 对当前 global correction 是否足够。

$\chi_{\rm fp}$ 表示当前一个 new rank-one pair 的 inner alternating convergence。

因此 tower v1 不应使用：

$$ \chi_{\rm fp}<\varepsilon $$

替代 outer LATIN stop，也不应使用：

$$ \xi<\varepsilon $$

来提前停止一个尚未形成稳定 space-time pair 的 inner enrichment iteration。

---

# 40. inner fixed-point time integration 与 outer LATIN norm 的关系

current 1D LATIN global convergence norm 与 residual norms 中部分使用 trapezoidal integration。

本阶段建议只对 Eq. (70)–(72) new-mode fixed-point 采用统一的 BE-consistent right-endpoint time-slab metric：

$$ \sum_{n=1}^{N}\Delta t_n(\cdot)_n $$

原因是该 inner problem 的 spatial contractions 和 temporal recurrence 已经按同一 discrete time interpretation 建立。

这并不意味着本阶段要求立即修改 current outer LATIN $\xi$ 的 integration convention。

准确策略是：

> 先保证同一个 Eq. (70)–(72) enrichment subproblem 内部离散自洽；outer LATIN norm 是否需要统一重构，只有在 tower benchmark 中出现明确 inconsistency 证据时再单独讨论。

---

# 41. sign ambiguity 与 complete-pair criterion

即使 spatial normalization 固定：

$$ \|\vec{p}\|_M=1 $$

仍然存在 sign ambiguity：

$$ \vec{p}^*=-\vec{p} $$

$$ \vec{s}^*=-\vec{s} $$

$$ \lambda^*=-\lambda $$

$$ \dot\lambda^*=-\dot\lambda $$

其 complete physical corrections仍保持不变：

$$ \dot\lambda^*\vec{p}^*=\dot\lambda\vec{p} $$

$$ \lambda^*\vec{s}^*=\lambda\vec{s} $$

因此 complete-pair criterion 天然消除了 sign flip 的影响。

这进一步说明：

> fixed-point stop 应基于 separated product，而不是 individual factors。

---

# 42. 推荐的 tower v1 fixed-point data dependency

一个 enrichment candidate 的 data dependency 可整理为：

```text
frozen outer-stage data
    shifted defect Delta_bar(t,q)
    H_sigma(t,q)
    M
    H
    C0
    existing basis

        |
        v

construct non-zero spatial seed p^(0)
        |
scale-normalise with M
        |
derive s^(0) = C0(E_tower - I)p^(0)
        |
Eq. (72) temporal solve
        |
complete pair P^(0)
        |
        v

fixed-point sweep k

    lambda^(k), lambda_dot^(k)
        |
        v
    BE temporal contractions
        |
        v
    Eq. (70) spatial FE solve
        |
        v
    Eq. (71) recover p_raw^(k+1)
        |
        v
    scale-normalise p^(k+1)
        |
        v
    derive s^(k+1)
        |
        v
    Eq. (72) temporal solve
        |
        v
    lambda^(k+1), lambda_dot^(k+1)
        |
        v
    build complete pair z^(k+1)
        |
        +--> chi_fp
        +--> chi_c
        +--> eta_n, rho_n
        +--> finite / degeneracy diagnostics
        |
        v
    converged?
        | no
        +--------> next sweep
        |
       yes
        v

post-fixed-point basis processing
    Gram-Schmidt
    temporal coefficient transformation
    significance / acceptance checks
    subsequent temporal re-optimisation strategy
```

---

# 43. 推荐的 tower v1 inner fixed-point pseudo-algorithm

```text
INPUT
    Delta_bar[n, q]
    H_sigma[n, q]
    time[n]
    material-point metric M
    tower strain operator H
    reference elasticity C0
    existing PGD basis

INITIALISE
    build deterministic non-zero p_raw^(0)
    scale-normalise p^(0) with ||p^(0)||_M = 1
    derive s^(0) from the reference equilibrium operator
    solve Eq. (72) for lambda^(0), lambda_dot^(0)
    build z^(0)

FOR k = 0, 1, 2, ...

    SPATIAL HALF-STEP
        hold lambda^(k), lambda_dot^(k) fixed
        evaluate BE-consistent temporal contractions
        build W^(k) and delta_bar^(k)
        solve Eq. (70)
        recover p_raw^(k+1) from Eq. (71)
        test spatial norm and finiteness
        scale-normalise p^(k+1)
        derive equilibrated s^(k+1)

    TEMPORAL HALF-STEP
        hold p^(k+1), s^(k+1) fixed
        enforce lambda_0^(k+1) = 0
        solve lambda_dot_0^(k+1) independently
        march Eq. (72) sequential BE for n >= 1
        record eta_n and rho_n

    COMPLETE-PAIR DIAGNOSTICS
        build z^(k+1)
        compute chi_fp^(k+1)
        compute chi_c^(k+1)
        check finite fields
        check near-negligible pair magnitude
        check temporal controllability diagnostics

    IF candidate is genuinely degenerate
        reject or restart candidate

    IF chi_fp^(k+1) <= fp_tolerance
        declare inner fixed point converged
        BREAK

    set z^(k) <- z^(k+1)

POST FIXED-POINT
    perform Gram-Schmidt against existing spatial basis
    transform existing and new temporal functions consistently
    perform new-mode significance and residual-reduction checks
    determine final acceptance
```

---

# 44. 本阶段相对于 current 1D 的继承与调整

| 内容 | current 1D | tower v1 本阶段决定 |
|---|---|---|
| spatial seed | weighted residual 最大时间截面 | 作为优先候选，尚待最终实现细节冻结 |
| first valid iterate | spatial seed → temporal solve | 继承 |
| spatial half-step | weighted least-squares spatial solve | 改为 paper-derived Eq. (70)–(71) FE solve |
| temporal half-step | sequential BE | 原则继承，并使用上一阶段 tower derivation |
| $t_0$ | amplitude zero，rate independent solve | 继承并已理论化 |
| fixed-point check 时机 | spatial update 后、final temporal recompute 前 | 改为每次 Eq. (72) 后 |
| inner spatial normalization | 有 | 保留 scaling normalization |
| inner Gram–Schmidt against existing basis | 每次 spatial solve 内执行 | tower v1 暂不在 inner fixed point 内执行 |
| primary fixed-point variable | combined correction $c$ | complete pair $z=(\Delta\dot\varepsilon^p,\Delta\sigma')$ |
| primary time integration | trapezoidal | BE-consistent right-endpoint sum |
| fixed-point relative form | current iterate denominator | symmetric relative change |
| residual-control correction $c$ | primary | secondary diagnostic |
| post-loop temporal recompute | 有 | tower complete sweep 已内置，不需作为补救步骤 |
| post-loop scalar line search | 有 | 暂不自动继承，下一阶段与 acceptance 一起决定 |
| default tolerance | $10^{-6}$ | 可作为 provisional reference，不视为 universal value |

---

# 45. 本阶段已关闭的问题

本阶段可以认为以下理论问题已获得明确答案。

## 45.1 一轮 fixed-point 的完整顺序

**固定顺序：`Spatial half-step → Temporal half-step → Complete-pair check`。**

## 45.2 convergence check 的位置

必须在 Eq. (72) temporal re-solve 后进行，以保证比较 mutually consistent complete pairs。

## 45.3 fixed-point 内 normalization 的目的

只固定 PGD reciprocal scaling gauge，不承担 basis orthogonality maintenance。

## 45.4 Gram–Schmidt 的位置

建议放在 inner fixed point 收敛后，再与 existing basis 做 orthogonalisation。

## 45.5 fixed-point comparison variable

不比较单独 $\vec{p}$，不比较单独 $\lambda$，而比较：

$$ z=(\Delta\dot{\vec\varepsilon}^p,\Delta\vec\sigma') $$

## 45.6 primary fixed-point metric

采用 search-direction-weighted complete mechanical graph norm。

## 45.7 current 1D `fixed_point_change` 的迁移结论

核心思想继承：比较 physical separated correction。

exact formula 不原样复制：tower 改为 complete-pair graph norm + Eq. (72) 后评价 + BE-consistent time weighting。

## 45.8 tolerance 的含义

$\varepsilon_{\rm fp}$ 是 inner rank-one alternating numerical convergence tolerance，不是 FOM error、fatigue error、outer LATIN error 或 mode acceptance tolerance。

---

# 46. 本阶段仍未关闭的问题

为了避免把本阶段结论延伸过度，下列内容仍需下一阶段专门闭合。

1. fixed-point converged raw new pair 做 Gram–Schmidt 后，existing temporal functions correction 与 new temporal function rescaling 的 exact implementation order；
2. Gram–Schmidt transformation 后是否立即对全部 temporal coefficients 做一次 Eq. (58)–(59) style re-optimisation；
3. paper 中“modified new time function insignificant”在 tower discrete setting 下的具体 norm 与 rejection threshold；
4. new mode acceptance 中 residual reduction、mode norm、$\eta_n$、$\rho_n$、linear dependence 与 outer saturation criterion 的 exact combination；
5. current 1D post-loop scalar optimal scaling 是否仍有必要；
6. fixed-point maximum iterations、restart strategy 与 seed replacement strategy；
7. $\varepsilon_{\rm fp}=10^{-6}$ 是否适用于 tower，仍需 benchmark calibration；
8. fixed-point 内若出现 $\eta_n\ll1$ 但 residual reduction 尚可，candidate 应继续、restart 还是 reject 的具体 threshold policy；
9. accepted new mode 加入 basis 后，是否立即重新计算 complete existing-basis temporal functions，再回到 Eq. (60) saturation check。

这些问题共同构成下一阶段：

> post-fixed-point basis insertion and acceptance sequence。

---

# 47. 对后续代码设计的直接约束

虽然当前阶段仍不进入 tower LATIN-PGD 正式编码，但未来代码接口至少应尊重以下结构。

new-mode enrichment routine 不应只返回一个 spatial vector 或一个 temporal array，而应返回一个完整 candidate object，至少包含：

```text
spatial_plastic_mode
spatial_stress_mode
temporal_amplitude
temporal_rate
fixed_point_history
residual_control_history
controllability_history
alignment_history
fixed_point_iterations
fixed_point_converged
candidate_degenerate
candidate_negligible
```

inner spatial solver 与 temporal solver 应拆分为独立函数，以便严格实现：

```text
spatial_half_step(lambda, lambda_dot, frozen_data)
temporal_half_step(p, s, frozen_data)
```

fixed-point controller 负责：

```text
normalization
data handoff
complete-pair diagnostics
convergence decision
restart/reject logic
```

而 Gram–Schmidt、basis transformation 与 final acceptance 应位于 inner controller 之外的 post-processing 层。

这将使代码结构直接映射本阶段理论层级，避免把 inner nonlinear solve、basis maintenance 与 mode acceptance 混在一个函数中。

---

# 48. 本阶段最核心的算法认识

本阶段最重要的认识不是某一个 tolerance 数值，而是 fixed-point object 的重新定义。

new-mode enrichment 的真正 fixed-point unknown 不是：

$$ \vec{p} $$

也不是：

$$ \lambda(t) $$

而是二者共同构造的 physical separated mechanical correction：

$$ \boxed{z=\left(\dot\lambda\vec{p},\lambda\vec{s}\right)} $$

spatial factor 与 temporal factor 只是这个 rank-one space-time correction 的一种非唯一分解。

因此：

- normalization 只负责选择一种 factor scaling representation；
- Eq. (70)–(71) 更新 spatial factor；
- Eq. (72) 更新 temporal factor；
- convergence 应评价 complete separated correction 是否稳定；
- Gram–Schmidt 是 fixed point 之后的 basis representation transformation；
- mode acceptance 则评价这个稳定 correction 是否值得加入 reduced basis。

这五个层次必须严格区分。

---

# 49. 推荐冻结的 tower v1 theoretical specification

截至本阶段结束，建议将 tower v1 new-mode inner fixed-point 规范冻结为：

核心 fixed-point map 为：

$$ \mathcal P^{(k)} \xrightarrow{\mathrm{Eq.70-71}} \vec{p}^{(k+1)} \xrightarrow{\mathrm{normalise}} \vec{s}^{(k+1)} \xrightarrow{\mathrm{Eq.72}} \mathcal P^{(k+1)} $$

随后用 $\chi_{\rm fp}$ 作 convergence decision。

其中 normalization 使用：

$$ \|\vec{p}\|_M=1 $$

primary fixed-point norm 使用：

$$ \|z\|_{\rm fp,h}^2=\sum_{n=1}^{N}\Delta t_n\left[(\Delta\dot{\vec\varepsilon}^{p}_n)^TMD_{H,n}^{-1}\Delta\dot{\vec\varepsilon}^{p}_n+(\Delta\vec\sigma'_n)^TMD_{H,n}\Delta\vec\sigma'_n\right] $$

primary indicator 使用：

$$ \boxed{\chi_{\rm fp}^{(k+1)}=\frac{\|z^{(k+1)}-z^{(k)}\|_{\rm fp,h}}{\|z^{(k+1)}\|_{\rm fp,h}+\|z^{(k)}\|_{\rm fp,h}}} $$

secondary residual-control diagnostic 使用：

$$ \vec{c}=\Delta\dot{\vec\varepsilon}^p-D_H\Delta\vec\sigma' $$

以及：

$$ \chi_c^{(k+1)}=\frac{\|\vec{c}^{(k+1)}-\vec{c}^{(k)}\|_{H^{-1},h}}{\|\vec{c}^{(k+1)}\|_{H^{-1},h}+\|\vec{c}^{(k)}\|_{H^{-1},h}} $$

provisional fixed-point stop 为：

$$ \chi_{\rm fp}\le\varepsilon_{\rm fp} $$

current 1D reference value 可暂取：

$$ \varepsilon_{\rm fp}=10^{-6} $$

但该 tolerance 仍需 tower benchmark calibration。

---

# 50. 与下一阶段的准确接口

本阶段结束后，Eq. (61)–(72) 已经从：

即从 **new-pair residual definition**

推进到：

推进到 **new-pair alternating fixed-point convergence**。

下一阶段不应重新讨论 Eq. (72) temporal scalar recurrence，也不应重新讨论 $t_0$、$a_h$ 或 denominator conditioning，除非新问题直接暴露矛盾。

下一阶段的准确起点是：

> 已有一个 fixed-point converged raw new pair $\mathcal P_*$。现在需要把它安全、等价且数值稳定地插入 existing PGD basis。

下一阶段建议连续闭合以下链条：

下一阶段目标链条为：`raw converged pair → Gram–Schmidt → temporal coefficient transformation → possible re-optimisation → significance test → residual reduction test → accept/reject → return to saturation check`。

完成该链条后，Eq. (61)–(72) 从“决定需要 enrichment”到“new pair 正式进入 basis”的全过程才算真正闭合。

---

# 51. 本阶段引用与代码依据

本阶段讨论主要基于以下项目资料与原论文内容：

- Bhattacharyya et al. (2018), LATIN-PGD cycling damage paper，重点 Eq. (61)–(72) 与 Eq. (76)–(77)；
- `docs/2026-08-17-tower-latin-pgd-eq70-71-spatial-solve-plastic-mode-stage-summary.md`；
- `docs/2026-08-17-tower-latin-pgd-eq72-temporal-minimisation-dg0-stage-summary.md`；
- `docs/2026-08-17-tower-latin-pgd-eq72-backward-euler-scalar-update-stage-summary.md`；
- `docs/2026-08-17-tower-latin-pgd-eq72-initial-time-be-contractions-conditioning-stage-summary.md`；
- `docs/2026-08-16-tower-latin-pgd-eq61-64-enrichment-stage-summary.md`；
- `docs/2026-08-16-tower-latin-pgd-eq60-saturation-and-latin-norm-stage-summary.md`；
- `latin/pgd_enrichment.py`；
- `latin/pgd_time_update.py`；
- `latin/equilibrium_operator.py`；
- `latin/search_directions.py`；
- `latin/pgd_basis.py`；
- `docs/section_5_1_convergence_criteria.md`。

---

# 52. 阶段结论

本阶段完成的核心工作，是把此前相对独立的 Eq. (70)–(71) spatial solve 与 Eq. (72) temporal solve 连接成一个完整 new-mode alternating fixed point，并给出与 PGD scaling freedom、LATIN metric 以及 current 1D numerical experience相协调的 convergence definition。

最终结构为：

最终主循环可概括为：`non-zero spatial seed → Eq. (72) initial temporal solve → [Eq. (70)–(71) spatial → scale normalization → Eq. (72) temporal → χ_fp] repeat`。

fixed point 收敛后才进入：

fixed point 收敛后再进入 **Gram–Schmidt and basis insertion**。

primary convergence object 为：

$$ \boxed{z=(\Delta\dot{\vec\varepsilon}^p,\Delta\vec\sigma')} $$

而不是单独的 spatial factor 或 temporal factor。

这使 fixed-point criterion 对 reciprocal scaling 和 sign ambiguity 保持不变，并且把以下四个层级清晰分开：

必须保持四个层级互不混淆：**factor normalization ≠ inner fixed-point convergence ≠ basis orthogonalisation ≠ mode acceptance**。

至此，tower LATIN-PGD new-mode enrichment 的 inner fixed-point 主体已经具备可转化为算法规格的理论闭合度。下一阶段只处理 fixed-point 之后如何把该 raw pair 正式、安全地纳入 existing reduced basis。
