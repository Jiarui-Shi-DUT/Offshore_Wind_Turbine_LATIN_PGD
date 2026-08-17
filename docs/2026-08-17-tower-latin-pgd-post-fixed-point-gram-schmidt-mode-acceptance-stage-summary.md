# Tower LATIN-PGD post-fixed-point Gram–Schmidt、temporal correction 与 new-mode acceptance 阶段总结

**日期：2026-08-17**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**当前研究路线：Bhattacharyya et al. 原论文 $x-t$ LATIN-PGD → 2D fiber beam-column offshore wind turbine tower**  
**阶段范围：在 Eq. (70)–(72) enrichment alternating fixed-point loop、scale normalization 与 fixed-point convergence criterion 已闭合的基础上，继续闭合 fixed-point 收敛后的 weighted Gram–Schmidt orthogonalisation、existing temporal functions 的精确坐标变换、new temporal function rescaling、modified-new-time-function significance、candidate rollback、enlarged-basis all-mode temporal re-optimisation、full residual benefit check，以及 tower v1 单次 enrichment event 接受一个 new pair 的流程。**  
**上一阶段衔接：`2026-08-17-tower-latin-pgd-eq70-72-enrichment-fixed-point-convergence-stage-summary.md`**  
**下一阶段：闭合 accepted new pair 之后 global stage 的 finishing sequence，包括 hardening Eq. (73)–(74)、damage Eq. (75)、damage-history inheritance、relaxation 以及 Eq. (76) outer LATIN convergence 之间的准确数据依赖。**

---

# 1. 本阶段定位

上一阶段已经把 new separated pair 在 Eq. (70)–(72) 之间的 alternating fixed-point 主循环闭合。

对于第 $k$ 次 inner fixed-point sweep，已经确定采用：

```text
complete pair at iteration k
    ↓
fixed temporal function
    ↓
Eq. (70)–(71) spatial half-step
    ↓
spatial scale normalization
    ↓
derive equilibrated stress mode
    ↓
Eq. (72) temporal half-step
    ↓
build mutually consistent complete pair
    ↓
fixed-point convergence check
```

new pair 的完整机械 correction 写为：

$$ \Delta\dot{\vec{\varepsilon}}^{p,(k)}_n = \dot{\lambda}^{(k)}_n\vec{p}^{(k)} $$

以及：

$$ \Delta\vec{\sigma}'^{(k)}_n = \lambda^{(k)}_n\vec{s}^{(k)} $$

其中：

$$ \vec{s}^{(k)} = C_0(\mathcal E_{\rm tower}-I)\vec{p}^{(k)} $$

上一阶段进一步确定，inner fixed-point convergence 不应单独比较 spatial factor 或 temporal factor，而应比较完整 separated mechanical pair。

primary fixed-point norm 建议采用：

$$ \|z\|_{\rm fp,h}^2 = \sum_{n=1}^{N}\Delta t_n\left[(\Delta\dot{\vec{\varepsilon}}^p_n)^TM D_{H,n}^{-1}\Delta\dot{\vec{\varepsilon}}^p_n+(\Delta\vec{\sigma}'_n)^TM D_{H,n}\Delta\vec{\sigma}'_n\right] $$

并定义 symmetric relative change：

$$ \chi_{\rm fp}^{(k+1)} = \frac{\|z^{(k+1)}-z^{(k)}\|_{\rm fp,h}}{\|z^{(k+1)}\|_{\rm fp,h}+\|z^{(k)}\|_{\rm fp,h}} $$

当：

$$ \chi_{\rm fp}^{(k+1)} \le \varepsilon_{\rm fp} $$

时，才认为 current rank-one alternating solver 已经收敛。

因此本阶段的起点不是重新讨论 Eq. (70)–(72)，而是：

> fixed point 已经得到一个 mutually consistent raw new pair 后，如何把它转换成真正可加入已有 PGD basis 的 accepted mode？

---

# 2. 本阶段需要闭合的问题

fixed-point 收敛后，仍然存在一系列尚未闭合的问题。

第一，原论文要求 spatial basis Gram–Schmidt orthonormalisation，但 tower material-point discretisation 下应采用什么 inner product？

第二，new spatial function 被正交化之后，为什么已有 temporal functions 必须修改？具体修改公式是什么？

第三，new temporal function 应如何 rescale，才能保证 physical correction 不变？

第四，原论文所说 modified new time function norm 很小时可以拒绝 pair，在 tower 离散框架中应如何解释？

第五，Gram–Schmidt 之后是否需要重新优化所有 temporal coefficients？

第六，joint temporal re-optimisation 应使用 inner enrichment 的 shifted defect，还是返回 Eq. (58)–(59) 的 full forcing？

第七，candidate 最终是否应该通过 full mechanical residual reduction 再做一次 acceptance check？

第八，若 candidate 被拒绝，应回滚 new mode 还是回滚整个 basis state？

第九，current 1D implementation 中的 post-loop scalar optimal scaling 是否继续迁移到 tower v1？

第十，一个 tower LATIN global-stage enrichment event 应允许连续加入多个 modes，还是先按照原论文流程一次只接受一个 pair？

本阶段围绕上述问题给出统一的 post-fixed-point algorithm specification。

---

# 3. 本阶段资料边界与结论类别

本阶段继续严格区分四类来源。

| 类别 | 本阶段含义 |
|---|---|
| 原论文明确内容 | Eq. (72) 后说明 spatial basis 采用 Gram–Schmidt orthonormalisation；former time functions 随之更新；new time function 也被修改；若 modified new time function norm insignificant，则 pair 可被拒绝。Fig. 2 给出 Add a pair → orthonormalise → correct time functions → test new time function → accept/reject 的流程。 |
| 由原论文与线性 basis transformation 推导 | weighted projection、orthogonal residual mode、existing/new temporal functions 的精确坐标变换、stress correction invariance、plastic correction invariance。 |
| current 1D implementation | `PGDBasis1D` 的 mode storage；accepted mode append 后调用 `update_pgd_time_functions()` 联合更新时间函数；单个 global-stage call 内可通过 while loop 连续加入多个 modes；`pgd_enrichment.py` 在最终 pair 后存在额外 scalar line search。 |
| tower v1 engineering choice | material-point $M$-weighted modified Gram–Schmidt；basis-health diagnostics；modified new temporal function 的 BE-consistent relative significance；exact transformation 后 field-invariance check；accepted candidate 后 Eq. (58)–(59) all-mode temporal re-optimisation；full residual benefit gate；完整 rollback；baseline 不继承 post-loop scalar line search；一次 enrichment event 最多接受一个 pair。 |

需要强调：

> 本阶段不是修改 Bhattacharyya et al. 的 LATIN-PGD 理论，而是把原论文 Eq. (72) 后只用文字描述的 basis-management 过程明确成可离散、可测试、可回滚的 tower algorithm specification。

---

# 4. fixed-point 收敛后的 raw new pair

设已有 $m$ 个 accepted PGD spatial modes：

$$ P_m = [\vec p_1,\ldots,\vec p_m] $$

对应的 equilibrated stress modes：

$$ S_m = [\vec s_1,\ldots,\vec s_m] $$

其中：

$$ \vec s_j = C_0(\mathcal E_{\rm tower}-I)\vec p_j $$

已有 temporal amplitudes 写为：

$$ \vec\lambda_m(t) = [\lambda_1(t),\ldots,\lambda_m(t)]^T $$

temporal rates 写为：

$$ \dot{\vec\lambda}_m(t) = [\dot\lambda_1(t),\ldots,\dot\lambda_m(t)]^T $$

inner fixed point 收敛后得到 raw candidate：

$$ \mathcal P_* = \{\vec p_*,\vec s_*,\lambda_*(t),\dot\lambda_*(t)\} $$

且：

$$ \vec s_* = C_0(\mathcal E_{\rm tower}-I)\vec p_* $$

上一阶段已经建议 fixed-point 内对 spatial mode 做 scale normalization，因此进入本阶段时通常有：

$$ \|\vec p_*\|_M = 1 $$

但：

$$ \vec p_* $$

此时并不保证与已有 basis 正交。

所以 raw candidate 仍不能直接作为第 $m+1$ 个 persistent PGD mode。

---

# 5. tower spatial inner product 的选择

tower 的 PGD spatial field 位于 flatten 后的 fiber material-point space。

material-point index 为：

$$ q=(e,g,f) $$

每个 material point 的 integration weight 为：

$$ v_q = A_{egf}w_gJ_e $$

因此空间积分 metric 已经定义为：

$$ M = \operatorname{diag}(v_q) $$

对于任意两个 material-point spatial vectors $\vec u$ 与 $\vec v$，最自然的离散空间 inner product 是：

$$ \langle\vec u,\vec v\rangle_M = \vec u^TM\vec v $$

对应 norm：

$$ \|\vec u\|_M = \sqrt{\vec u^TM\vec u} $$

因此 tower PGD spatial orthogonalisation 应满足：

$$ \vec p_i^TM\vec p_j = \delta_{ij} $$

即：

$$ P_m^TMP_m = I $$

这一选择与前序阶段的 spatial scale normalization 使用同一个 material-point integration metric。

---

# 6. 为什么不能使用普通 Euclidean Gram–Schmidt

普通 Euclidean inner product：

$$ \vec u^T\vec v $$

会把所有 material points 当成等权离散样本。

但 tower 中不同 element、Gauss point、fiber 对结构积分的实际贡献由：

$$ v_q=A_{egf}w_gJ_e $$

决定。

如果直接使用普通 Euclidean norm，则：

- 不同 beam element 长度造成的 integration contribution 被忽略；
- Gauss quadrature 权重被忽略；
- fiber area 权重被忽略；
- basis orthogonality 与实际结构空间积分不一致。

因此 tower v1 应冻结：

$$ \boxed{\text{spatial orthogonalisation metric}=M} $$

而不是 identity metric。

---

# 7. weighted Gram–Schmidt 的基本 projection

假设已有 basis 已经满足：

$$ P_m^TMP_m=I $$

raw candidate 在已有 basis 上的 projection coefficients 为：

$$ \vec a=P_m^TM\vec p_* $$

即：

$$ a_j=\vec p_j^TM\vec p_* $$

已有 basis 能够表示的部分为：

$$ P_m\vec a $$

因此真正新的 orthogonal spatial component 为：

$$ \vec p_\perp=\vec p_*-P_m\vec a $$

由 projection construction 可得：

$$ P_m^TM\vec p_\perp=0 $$

这说明：

$$ \vec p_\perp $$

才是 raw candidate 中真正不属于已有 spatial subspace 的新信息。

---

# 8. 当已有 basis 不是精确正交时的更一般公式

理论上：

$$ P_m^TMP_m=I $$

但浮点计算中不能假设这一关系永远严格成立。

定义 weighted Gram matrix：

$$ G_m=P_m^TMP_m $$

则更一般的 projection coefficients 应通过：

$$ G_m\vec a=P_m^TM\vec p_* $$

求得。

若：

$$ G_m\approx I $$

则退化回：

$$ \vec a\approx P_m^TM\vec p_* $$

因此未来实现时应至少保留：

$$ \|G_m-I\| $$

作为 basis-health diagnostic。

---

# 9. tower v1 建议采用 weighted Modified Gram–Schmidt

原论文只明确要求 Gram–Schmidt orthonormalisation，没有规定 classical Gram–Schmidt 还是 modified Gram–Schmidt。

tower v1 建议采用：

$$ \boxed{\text{weighted Modified Gram--Schmidt}} $$

必要时允许一次 re-orthogonalisation。

原因是 PGD enrichment 可能逐渐产生接近已有 subspace 的 candidate，classical Gram–Schmidt 在 mode 数增长时更容易积累 loss of orthogonality。

对于第 $j$ 个已有 mode，可以逐次执行：

$$ a_j=\vec p_j^TM\vec p_{\rm work} $$

然后：

$$ \vec p_{\rm work}\leftarrow\vec p_{\rm work}-a_j\vec p_j $$

所有已有 modes 处理完成后：

$$ \vec p_\perp=\vec p_{\rm work} $$

若 basis-health diagnostic 表明正交性不足，可再执行一次相同 projection。

这属于 numerical stability choice，不改变原论文的 PGD mathematical space。

---

# 10. new orthogonal spatial mode 的 normalization

定义：

$$ c=\|\vec p_\perp\|_M $$

若：

$$ c>0 $$

则 normalized new spatial mode 为：

$$ \vec p_{m+1}=\frac{\vec p_\perp}{c} $$

因此 raw candidate 可以严格写成：

$$ \vec p_*=P_m\vec a+c\vec p_{m+1} $$

这条 decomposition 是本阶段最重要的关系之一。

它说明 Gram–Schmidt 不是简单把 old projection 删除，而是把 raw candidate 改写成：

- 已有 basis 上的坐标贡献；
- 一个新的 orthogonal spatial direction。

因此若希望 physical field 完全不变，就必须把被投影到 old basis 上的部分重新转移到 existing temporal coefficients 中。

---

# 11. 当 raw spatial mode 已 normalization 时的 identity

若 fixed-point 收敛后已有：

$$ \|\vec p_*\|_M=1 $$

且 existing basis 严格 orthonormal，则：

$$ \|\vec p_*\|_M^2=\|P_m\vec a\|_M^2+\|\vec p_\perp\|_M^2 $$

因此：

$$ 1=\vec a^T\vec a+c^2 $$

即：

$$ c^2=1-\vec a^T\vec a $$

这一关系适合作为 consistency diagnostic。

但是程序计算 $c$ 时不建议使用：

$$ c=\sqrt{1-\vec a^T\vec a} $$

因为浮点误差可能导致：

$$ 1-\vec a^T\vec a<0 $$

的微小数值异常。

实际应直接计算：

$$ c=\sqrt{\vec p_\perp^TM\vec p_\perp} $$

再用：

$$ 1-\vec a^T\vec a-c^2 $$

作为 orthogonality consistency residual。

---

# 12. spatial novelty diagnostic

定义尺度无关 spatial novelty ratio：

$$ \gamma_{\rm sp}=\frac{\|\vec p_\perp\|_M}{\|\vec p_*\|_M} $$

在 fixed-point 已 normalization 时：

$$ \gamma_{\rm sp}=c $$

若：

$$ \gamma_{\rm sp}\approx1 $$

表示 raw candidate 与已有 spatial subspace 基本正交，提供明显新 spatial direction。

若：

$$ \gamma_{\rm sp}\ll1 $$

表示：

$$ \vec p_*\approx P_m\vec a $$

即 raw candidate 几乎完全位于已有 spatial subspace 中。

因此：

$$ \boxed{\gamma_{\rm sp}=\text{candidate spatial independence diagnostic}} $$

它反映的是 candidate 是否带来新的 spatial information。

---

# 13. $\gamma_{\rm sp}$ 的作用边界

$\gamma_{\rm sp}$ 不应直接等价于最终 mode usefulness。

因为即使 spatial direction 很新：

$$ \gamma_{\rm sp}\not\ll1 $$

其 temporal amplitude 也可能几乎为零。

反过来，一个 candidate 即使与旧 basis 接近，只要其 raw temporal amplitude 很大，Gram–Schmidt 后真正新增的 orthogonal component 是否仍有意义，也需要通过 modified temporal function 判断。

因此 tower v1 中：

$$ \boxed{\gamma_{\rm sp}=\text{novelty / linear-dependence diagnostic}} $$

而不是唯一 acceptance criterion。

---

# 14. Gram–Schmidt 前的完整 plastic correction

在 raw candidate 加入前，existing basis 的 accumulated plastic correction 为：

$$ \Delta\vec{\varepsilon}^p_m(t)=P_m\vec\lambda_m(t) $$

加入 raw candidate 后：

$$ \Delta\vec{\varepsilon}^p_{\rm raw}(t)=P_m\vec\lambda_m(t)+\lambda_*(t)\vec p_* $$

将：

$$ \vec p_*=P_m\vec a+c\vec p_{m+1} $$

代入：

$$ \Delta\vec{\varepsilon}^p_{\rm raw}=P_m\vec\lambda_m+P_m\vec a\,\lambda_*+c\lambda_*\vec p_{m+1} $$

整理得：

$$ \Delta\vec{\varepsilon}^p_{\rm raw}=P_m[\vec\lambda_m+\vec a\lambda_*]+\vec p_{m+1}[c\lambda_*] $$

因此 Gram–Schmidt 后，existing temporal functions 与 new temporal function 的更新公式不是经验选择，而是由 exact basis-coordinate transformation 强制决定。

---

# 15. existing temporal amplitudes 的精确修正

对于：

$$ j=1,\ldots,m $$

former temporal amplitudes 必须更新为：

$$ \lambda_j^+(t)=\lambda_j(t)+a_j\lambda_*(t) $$

vector form 为：

$$ \vec\lambda_m^+(t)=\vec\lambda_m(t)+\vec a\,\lambda_*(t) $$

new orthogonal mode 的 temporal amplitude 为：

$$ \lambda_{m+1}^+(t)=c\lambda_*(t) $$

于是：

$$ P_m\vec\lambda_m+\lambda_*\vec p_*=P_m\vec\lambda_m^++\lambda_{m+1}^+\vec p_{m+1} $$

这说明：

> Gram–Schmidt 后 former time functions 的 modification 本质上是 basis-coordinate change，而不是重新求解物理问题。

---

# 16. 为什么 new temporal function 必须乘 $c$

new spatial mode 被 normalization：

$$ \vec p_{m+1}=\frac{\vec p_\perp}{c} $$

如果 temporal amplitude 不做 reciprocal compensation，则新 orthogonal correction 会被人为放大：

$$ \lambda_*\vec p_{m+1}=\frac{\lambda_*}{c}\vec p_\perp $$

这与 raw candidate 中真正的 orthogonal component：

$$ \lambda_*\vec p_\perp $$

不同。

因此必须有：

$$ \lambda_{m+1}^+=c\lambda_* $$

才能保证：

$$ \lambda_{m+1}^+\vec p_{m+1}=\lambda_*\vec p_\perp $$

所以 $c$ 不是 arbitrary rescaling coefficient，而是 Gram–Schmidt normalization 与 temporal coordinate transformation 之间的严格 coupling coefficient。

---

# 17. temporal rate 必须同步做相同变换

new pair 的 plastic-rate correction 为：

$$ \Delta\dot{\vec{\varepsilon}}^p(t)=\sum_j\dot\lambda_j(t)\vec p_j $$

对 time derivative 进行同样 coordinate transformation：

$$ \dot\lambda_j^+(t)=\dot\lambda_j(t)+a_j\dot\lambda_*(t) $$

以及：

$$ \dot\lambda_{m+1}^+(t)=c\dot\lambda_*(t) $$

于是：

$$ \sum_{j=1}^{m}\dot\lambda_j\vec p_j+\dot\lambda_*\vec p_*=\sum_{j=1}^{m}\dot\lambda_j^+\vec p_j+\dot\lambda_{m+1}^+\vec p_{m+1} $$

因此 accumulated plastic correction 与 plastic-rate correction 均保持不变。

---

# 18. temporal amplitude 与 temporal rate 不应在 transformation 后独立重建

在 tower v1 中，Gram–Schmidt 的 exact coordinate transformation 应同时作用于：

$$ \lambda $$

与：

$$ \dot\lambda $$

而不应先修改 $\lambda$，再用另一套 arbitrary finite-difference rule 重新构造所有 $\dot\lambda$。

原因是当前 temporal solver 已经明确区分：

$$ t_0:\quad \lambda_0=0,\quad \dot\lambda_0\ \text{independently solved} $$

与：

$$ n\ge1:\quad \dot\lambda_n=\frac{\lambda_n-\lambda_{n-1}}{\Delta t_n} $$

尤其 $t_0$ 的 rate 并非由 amplitude finite difference 得到。

因此 exact basis transformation 时，已有：

$$ \dot\lambda_* $$

应直接同步线性变换，保证 current physical plastic-rate correction 严格不变。

后续若执行 all-mode temporal re-optimisation，再由 temporal solver 重新生成一套新的 mutually consistent amplitudes 与 rates。

---

# 19. stress spatial modes 的 transformation

定义 reference linear stress operator：

$$ \mathcal L=C_0(\mathcal E_{\rm tower}-I) $$

因此：

$$ \vec s_j=\mathcal L\vec p_j $$

以及：

$$ \vec s_*=\mathcal L\vec p_* $$

由于：

$$ \vec p_*=P_m\vec a+c\vec p_{m+1} $$

且 $\mathcal L$ 线性，所以：

$$ \vec s_*=S_m\vec a+c\vec s_{m+1} $$

其中：

$$ \vec s_{m+1}=\mathcal L\vec p_{m+1} $$

因此 same temporal coordinate transformation 自动保持 stress correction。

---

# 20. stress correction invariance

Gram–Schmidt 前：

$$ \Delta\vec\sigma'_{\rm raw}(t)=S_m\vec\lambda_m(t)+\lambda_*(t)\vec s_* $$

代入：

$$ \vec s_*=S_m\vec a+c\vec s_{m+1} $$

得到：

$$ \Delta\vec\sigma'_{\rm raw}=S_m[\vec\lambda_m+\vec a\lambda_*]+c\lambda_*\vec s_{m+1} $$

因此：

$$ \Delta\vec\sigma'_{\rm raw}=S_m\vec\lambda_m^++\lambda_{m+1}^+\vec s_{m+1} $$

也就是说：

$$ \boxed{\Delta\vec\sigma'_{\rm GS}=\Delta\vec\sigma'_{\rm raw}} $$

理论上 exact coordinate transformation 不应改变任何 stress correction。

---

# 21. Gram–Schmidt 本质上不是新的物理求解

这一阶段最重要的概念之一是：

$$ \boxed{\text{Gram--Schmidt + temporal correction}=\text{representation change}} $$

而不是：

$$ \boxed{\text{new physical minimisation}} $$

因此完成 exact transformation 后，以下三个 field 都必须保持不变：

$$ \Delta\vec\varepsilon^p_{\rm GS}=\Delta\vec\varepsilon^p_{\rm raw} $$

$$ \Delta\dot{\vec\varepsilon}^p_{\rm GS}=\Delta\dot{\vec\varepsilon}^p_{\rm raw} $$

$$ \Delta\vec\sigma'_{\rm GS}=\Delta\vec\sigma'_{\rm raw} $$

相应地，current mechanical residual 也必须保持不变。

---

# 22. field-invariance test

未来实现时建议显式记录：

$$ e_{\varepsilon^p}^{\rm GS}=\frac{\|\Delta\vec\varepsilon^p_{\rm GS}-\Delta\vec\varepsilon^p_{\rm raw}\|}{\max(\|\Delta\vec\varepsilon^p_{\rm raw}\|,\epsilon)} $$

$$ e_{\dot\varepsilon^p}^{\rm GS}=\frac{\|\Delta\dot{\vec\varepsilon}^p_{\rm GS}-\Delta\dot{\vec\varepsilon}^p_{\rm raw}\|}{\max(\|\Delta\dot{\vec\varepsilon}^p_{\rm raw}\|,\epsilon)} $$

$$ e_{\sigma}^{\rm GS}=\frac{\|\Delta\vec\sigma'_{\rm GS}-\Delta\vec\sigma'_{\rm raw}\|}{\max(\|\Delta\vec\sigma'_{\rm raw}\|,\epsilon)} $$

这些 quantity 应仅处于 floating-point roundoff / linear-solver tolerance level。

如果 Gram–Schmidt 后 residual 明显发生变化，则优先说明：

- projection coefficient 计算错误；
- temporal transformation 符号错误；
- stress mode 没有按 same spatial transformation 构造；
- rate transformation 不一致；
- existing basis 非正交而 projection 仍错误假设 $G_m=I$。

而不应把这种变化解释为“orthogonalisation 提升或降低了 PGD solution”。

---

# 23. basis-health diagnostics

除了 field invariance，未来应记录：

$$ e_{\rm ortho}=\|P_{m+1}^TMP_{m+1}-I\| $$

并可进一步检查：

$$ e_{{\rm stress},j}=\|\vec s_j-C_0(\mathcal E_{\rm tower}-I)\vec p_j\| $$

对于新 mode：

$$ e_{{\rm stress},m+1}=\|\vec s_{m+1}-C_0(\mathcal E_{\rm tower}-I)\vec p_{m+1}\| $$

这样可以把：

- spatial basis orthogonality；
- reference equilibrium mapping；
- temporal coordinate transformation；

三个不同层级的数值错误分离诊断。

---

# 24. 原论文 modified new time function significance 的离散解释

原论文在 Eq. (72) 后说明：

> Gram–Schmidt 后 new time function 被修改；如果 modified new time function norm 很小，则该 space-time pair 可以拒绝。

原论文没有明确：

- 使用什么 time norm；
- 是否 weighted；
- tolerance 是多少。

因此 tower v1 需要给出与 current BE discretisation 一致的离散解释。

---

# 25. BE-consistent temporal norm

当前 Eq. (70)–(72) 已经采用 right-endpoint backward-Euler time-slab interpretation。

因此 modified temporal amplitude 的 time norm 建议定义为：

$$ \|\lambda_j\|_{T,h}^2=\sum_{n=1}^{N}\Delta t_n\lambda_{j,n}^2 $$

这里不包含 $t_0$ 的独立 nodal weight，而是把 time integral 按每个 slab 的 right endpoint 近似。

这样 temporal significance metric 与 Eq. (70)–(72) 的离散 contraction 保持同一套 time convention。

---

# 26. temporal norm 与 orthogonal plastic correction 的关系

Gram–Schmidt 后：

$$ P_{m+1}^TMP_{m+1}=I $$

因此某个 normalized spatial mode：

$$ \vec p_j $$

满足：

$$ \|\vec p_j\|_M=1 $$

对应 accumulated plastic correction：

$$ \lambda_j(t)\vec p_j $$

在单个时间点有：

$$ \|\lambda_j(t)\vec p_j\|_M^2=\lambda_j(t)^2 $$

对时间进行 BE-consistent integration：

$$ \sum_{n=1}^{N}\Delta t_n\|\lambda_{j,n}\vec p_j\|_M^2=\|\lambda_j\|_{T,h}^2 $$

因此：

$$ \boxed{\|\lambda_j\|_{T,h}=\text{该 normalized spatial mode 的 space-time accumulated-plastic magnitude}} $$

这为原论文 modified new time function norm 提供了直接的 tower material-point interpretation。

---

# 27. new modified time function

Gram–Schmidt 后：

$$ \lambda_{m+1}^+=c\lambda_* $$

因此：

$$ \|\lambda_{m+1}^+\|_{T,h}=c\|\lambda_*\|_{T,h} $$

若：

$$ c\ll1 $$

则即使 raw $\lambda_*$ 不小，真正 orthogonal new direction 的 contribution 仍然可能很小。

这正是为什么不能在 Gram–Schmidt 前仅仅根据 raw temporal amplitude 判断 candidate significance。

---

# 28. relative modified-time-function significance

绝对 threshold：

$$ \|\lambda_{m+1}^+\|_{T,h}<\varepsilon $$

会依赖：

- loading magnitude；
- physical units；
- spatial normalization；
- current global-stage correction scale。

因此 tower v1 建议采用 relative significance：

$$ \gamma_\lambda=\frac{\|\lambda_{m+1}^+\|_{T,h}}{\sqrt{\sum_{j=1}^{m+1}\|\lambda_j^+\|_{T,h}^2}} $$

由于 spatial basis 已 M-orthonormal，分母对应当前全部 orthogonal plastic modes 的 combined space-time magnitude。

因此：

$$ 0\le\gamma_\lambda\le1 $$

并且：

$$ \boxed{\gamma_\lambda=\text{new orthogonal mode 在当前 accumulated plastic correction 中的相对份额}} $$

---

# 29. paper-faithful significance gate

可以把原论文的 qualitative statement 离散化为：

$$ \gamma_\lambda\le\varepsilon_{\rm mode}\quad\Rightarrow\quad\text{candidate insignificant} $$

这里：

$$ \varepsilon_{\rm mode} $$

是 tower v1 的 numerical tolerance，而不是原论文给定常数。

本阶段不应贸然指定：

$$ \varepsilon_{\rm mode}=10^{-x} $$

因为需要后续通过：

- validated 1D three-material bar；
- first tower benchmark；
- accepted/rejected candidate distributions；
- basis novelty distributions；
- residual reduction distributions；

共同标定。

---

# 30. $\gamma_{\rm sp}$ 与 $\gamma_\lambda$ 的区别

spatial novelty：

$$ \gamma_{\rm sp}=\frac{\|\vec p_\perp\|_M}{\|\vec p_*\|_M} $$

回答：

> raw candidate 的 spatial direction 与已有 spatial subspace 有多独立？

modified temporal significance：

$$ \gamma_\lambda $$

回答：

> 真正 orthogonalised new mode 在 current space-time plastic correction 中占多少？

所以：

$$ \boxed{\gamma_{\rm sp}\neq\gamma_\lambda} $$

前者主要是 linear dependence / basis novelty diagnostic。

后者更直接对应原论文所说 modified new time function significance。

---

# 31. 两类典型 candidate

第一类 candidate：

$$ \gamma_{\rm sp}\ll1 $$

说明 raw mode 几乎位于已有 basis 中。

这类 candidate 可能导致：

$$ c\approx0 $$

进而 normalized new spatial direction 对 roundoff 极其敏感。

应视为 spatially degenerate / linearly dependent candidate。

第二类 candidate：

$$ \gamma_{\rm sp}\not\ll1 $$

但：

$$ \gamma_\lambda\ll1 $$

说明 spatial direction 虽然独立，但当前 correction 对它几乎没有 temporal demand。

这类 candidate 在物理 reduced correction 中贡献很小，也没有必要加入 persistent basis。

---

# 32. exact transformation 与 temporal re-optimisation 必须严格分离

完成 Gram–Schmidt 后有两类完全不同的 operation。

第一类：

$$ \boxed{\text{exact coordinate transformation}} $$

包括：

$$ \vec p_*\rightarrow\{\vec a,c,\vec p_{m+1}\} $$

以及：

$$ \{\vec\lambda_m,\lambda_*\}\rightarrow\{\vec\lambda_m^+,\lambda_{m+1}^+\} $$

它不应改变 physical fields。

第二类：

$$ \boxed{\text{all-mode temporal re-optimisation}} $$

它会在 enlarged spatial basis 上重新最小化 current reduced residual，因此会真正改变 temporal solution。

这两类操作不能混成一步。

---

# 33. 为什么 exact transformation 后才允许 re-optimisation

若 Gram–Schmidt 与 temporal re-solve 同时执行，一旦 residual 改变，就无法判断原因来自：

- basis transformation error；
- temporal solver error；
- enlarged basis 的真实优化效果。

因此 tower v1 的 testable sequencing 应为：

```text
raw converged pair
    ↓
Gram-Schmidt
    ↓
exact temporal coordinate transformation
    ↓
field-invariance test
    ↓
significance test
    ↓
tentative enlarged basis
    ↓
all-mode temporal re-optimisation
```

这样每一个步骤的数学职责明确、diagnostic 可单独定位。

---

# 34. 为什么 enlarged basis 后应重新优化所有 temporal coefficients

inner new-mode fixed point 求的是：

$$ \text{existing reduced solution}+\text{one residual rank-one correction} $$

在寻找 raw candidate 时，已有 modes 的 temporal functions 暂时作为 current baseline。

但是 candidate accepted 后 spatial approximation space 已从：

$$ \operatorname{span}\{\vec p_1,\ldots,\vec p_m\} $$

扩展为：

$$ \operatorname{span}\{\vec p_1,\ldots,\vec p_m,\vec p_{m+1}\} $$

在这个 enlarged spatial subspace 中，原有：

$$ \lambda_1,\ldots,\lambda_m $$

一般不再是 jointly optimal temporal coordinates。

因此 tower v1 建议 accepted candidate 在正式进入 global-stage solution 前执行一次：

$$ \boxed{\text{Eq. (58)–(59) all-mode temporal re-optimisation}} $$

即固定 enlarged spatial basis，联合更新全部 temporal functions。

---

# 35. current 1D implementation 的对应做法

current 1D `pgd_global_stage.py` 在一个 enrichment 被接受并 append 后，会调用：

```text
update_pgd_time_functions(...)
```

重新更新 working basis 中全部 temporal coefficients。

因此 current 1D 已经采用：

$$ \boxed{\text{accepted new spatial mode}\rightarrow\text{joint temporal re-optimisation}} $$

这一思想。

tower v1 可以继承这个 architecture-level decision，但 temporal residual、spatial modes 与 discrete metrics 要使用当前已经重新推导的 tower formulations。

---

# 36. original paper 与 all-mode re-optimisation 的关系

原论文明确说 spatial orthonormalisation 后：

- former time functions are updated；
- new time function is modified。

但正文并未给出足够细节来唯一判断这种 update 是否仅指：

$$ \lambda_j^+=\lambda_j+a_j\lambda_* $$

这种 exact basis transformation，

还是还包含随后一次 Eq. (58)–(59) reduced temporal re-solve。

因此需要保持准确表述：

> exact temporal coordinate transformation 是保持 physical correction 不变所必需的数学步骤；其后进行 all-mode Eq. (58)–(59) re-optimisation，是 tower v1 推荐的 numerical choice，并与 current validated 1D implementation 的 architecture 一致。

---

# 37. re-optimisation 必须使用 full global-stage forcing

inner enrichment Eq. (61)–(72) 使用的是 current existing-basis update 后剩余的 shifted defect。

前序阶段定义：

$$ \bar\Delta=H_\sigma(\hat\sigma-\sigma^{\rm up})-(\hat{\dot\varepsilon}^p-\dot\varepsilon^{p,\rm up}) $$

new rank-one mode 的任务是针对：

$$ \bar\Delta $$

补充当前 basis 无法表示的 remaining defect。

但是 candidate 一旦进入 enlarged basis，joint temporal re-optimisation 的目标重新变成：

$$ \boxed{\text{Eq. (58)–(59) fixed-basis reduced global problem}} $$

因此必须使用 full plastic global-stage forcing：

$$ f $$

而不是继续使用：

$$ \bar\Delta $$

所以：

$$ \boxed{\bar\Delta\rightarrow\text{new-mode generation}} $$

而：

$$ \boxed{f\rightarrow\text{enlarged-basis all-mode temporal update}} $$

这一数据边界应在未来代码接口中明确分开。

---

# 38. re-optimisation 前后的 basis state

在 exact transformation 后，定义 tentative orthonormal basis：

$$ \tilde P_{m+1}=[P_m,\vec p_{m+1}] $$

corresponding stress basis：

$$ \tilde S_{m+1}=[S_m,\vec s_{m+1}] $$

exact transformed temporal coordinates 为：

$$ \tilde{\vec\lambda}_{m+1}^+ $$

与：

$$ \dot{\tilde{\vec\lambda}}_{m+1}^+ $$

此时 physical fields 与 raw candidate 加入前完全相同。

然后 all-mode temporal re-optimisation 产生新的：

$$ \vec\lambda_{m+1}^{\rm opt} $$

以及：

$$ \dot{\vec\lambda}_{m+1}^{\rm opt} $$

这一步才真正改变 reduced solution。

---

# 39. re-optimisation 后的 full residual

设 enrichment 前 existing basis 经 Eq. (58)–(59) temporal update 后的 full mechanical residual 为：

$$ R_m(t) $$

candidate 进入 enlarged basis 并完成 all-mode re-optimisation 后：

$$ R_{m+1}(t) $$

tower v1 建议使用与 Eq. (72) / enrichment diagnostics 一致的 discrete weighted norm：

$$ \|R\|_{r,h}^2=\sum_{n=1}^{N}\Delta t_n R_n^TM D_{H,n}^{-1}R_n $$

定义 residual ratio：

$$ q_{\rm res}=\frac{\|R_{m+1}\|_{r,h}}{\|R_m\|_{r,h}} $$

以及 residual reduction fraction：

$$ \Delta_{\rm res}=1-q_{\rm res} $$

---

# 40. residual benefit 的意义

若：

$$ \Delta_{\rm res}>0 $$

说明 enlarged basis 经 temporal re-optimisation 后，current full reduced mechanical residual 实际降低。

若：

$$ \Delta_{\rm res}\approx0 $$

说明 candidate 即使通过了 inner fixed point 和 orthogonalisation，对当前 global-stage reduced solution 的帮助仍然可以忽略。

若：

$$ \Delta_{\rm res}<0 $$

说明 re-optimised result 的 whole-time residual 反而增大，应进入 diagnostic / rejection 路径。

---

# 41. 为什么 sequential BE 下仍需要显式 full residual check

如果 temporal problem 被写成一个 whole-time coupled global least-squares，并且 enlarged space 严格包含 old space，则理论 minimizer 不应比 old-space solution 更差。

但是 current tower v1 temporal discretisation采用 sequential backward-Euler updates：

$$ \lambda_n\ \text{在给定}\ \lambda_{n-1}\ \text{条件下逐步求解} $$

它不是一个对全时间未知量同时求解的 global coupled minimisation。

因此不能仅凭“basis 扩大”就无条件保证：

$$ \|R_{m+1}\|_{r,h}\le\|R_m\|_{r,h} $$

对于当前 sequential implementation，post-reoptimisation full residual reconstruction 仍然是必要的 acceptance safeguard。

---

# 42. residual benefit threshold

最终可定义：

$$ \Delta_{\rm res}>\varepsilon_{\rm accept} $$

作为 actual reduced-solution benefit gate。

其中：

$$ \varepsilon_{\rm accept} $$

不等于：

$$ \varepsilon_{\rm fp} $$

也不等于：

$$ \varepsilon_{\rm mode} $$

三者控制完全不同的问题。

---

# 43. fixed-point convergence、mode significance 与 residual benefit 的三层结构

tower v1 建议 candidate screening 分为三层。

第一层：inner solver validity。

要求：

$$ \chi_{\rm fp}\le\varepsilon_{\rm fp} $$

并且：

- spatial / temporal variables finite；
- Eq. (72) controllability 没有真实退化；
- $\eta_n$、$\rho_n$ diagnostics 未显示不可接受异常；
- no exploding mode correction。

第二层：basis novelty / significance。

检查：

$$ \gamma_{\rm sp} $$

以及：

$$ \gamma_\lambda>\varepsilon_{\rm mode} $$

第三层：actual full reduced-solution benefit。

检查：

$$ \Delta_{\rm res}>\varepsilon_{\rm accept} $$

因此：

$$ \boxed{\text{fixed-point converged}\not\Rightarrow\text{mode accepted}} $$

---

# 44. 三个 tolerance 的不同含义

fixed-point tolerance：

$$ \varepsilon_{\rm fp} $$

表示：

> current raw rank-one separated pair 在 consecutive spatial-temporal sweeps 之间是否稳定。

mode-significance tolerance：

$$ \varepsilon_{\rm mode} $$

表示：

> Gram–Schmidt 后真正新增的 orthogonal mode 在 current space-time correction 中是否足够显著。

acceptance tolerance：

$$ \varepsilon_{\rm accept} $$

表示：

> enlarged basis 是否在 full reduced global-stage residual 上带来可识别改善。

三者必须分别记录，不应合并成一个通用“PGD tolerance”。

---

# 45. $\eta_n$ 与 $\rho_n$ 在 post-fixed-point 阶段的地位

前序阶段已经定义 temporal controllability ratio：

$$ \eta_n=\frac{\|\vec g_n\|_{W_n}}{\|\vec p/\Delta t_n\|_{W_n}+\|D_{H,n}\vec s\|_{W_n}} $$

以及 residual-alignment effectiveness：

$$ \rho_n=\frac{|\vec g_n^TW_n\vec b_n|}{\sqrt{(\vec g_n^TW_n\vec g_n)(\vec b_n^TW_n\vec b_n)}} $$

在本阶段不建议立即把：

$$ \eta_n $$

和：

$$ \rho_n $$

变成 hard final mode-acceptance thresholds。

更合适的定位是：

$$ \boxed{\eta_n,\rho_n=\text{conditioning / local effectiveness diagnostics}} $$

真正的最终 usefulness 由：

$$ \Delta_{\rm res} $$

直接衡量。

这样可以避免 tower v1 在尚未经过 benchmark 标定前叠加过多未经验证的 local thresholds。

---

# 46. current 1D post-loop scalar optimal scaling

current 1D `pgd_enrichment.py` 在 final temporal recomputation 后，还会对 complete separated correction 进行一次 scalar line search。

其目的，是针对 sequential BE temporal update 并非 whole-time global minimiser这一事实，通过一个额外 scalar：

$$ \alpha $$

调整 new mode contribution，使 whole space-time weighted residual 不增加。

这是 current 1D engineering safeguard。

原论文 Eq. (72) 后并未明确提出这一 scalar line search。

---

# 47. tower v1 baseline 是否继承 scalar line search

本阶段建议：

$$ \boxed{\text{tower v1 baseline 不继承 current 1D post-loop scalar line search}} $$

理由有三点。

第一，当前目标是首先保持 Eq. (58)–(72) algorithm correspondence 清楚。

第二，tower v1 已经增加：

$$ \text{all-mode temporal re-optimisation} $$

以及：

$$ \text{full residual benefit check} $$

因此 residual 是否改善可以直接检测，而无需额外修改 paper-derived temporal solution。

第三，scalar line search 会在 Eq. (72) 或 Eq. (58)–(59) 求得 temporal solution 后再次人工改变整个 temporal amplitude，增加理论解释层。

因此 tower v1 baseline 采用：

```text
solve original temporal equations
    ↓
reconstruct full residual
    ↓
accept or reject
```

而不是：

```text
solve original temporal equations
    ↓
modify solution by extra scalar
    ↓
accept
```

但 current 1D scalar line search 应保留为 future robustness / sensitivity variant，而不是永久删除这一经验。

---

# 48. candidate rejection 前为什么必须保存 basis snapshot

Gram–Schmidt 后：

$$ \lambda_j^+=\lambda_j+a_j\lambda_* $$

会修改已有 modes 的 temporal functions。

后续 all-mode re-optimisation 还会进一步修改：

$$ \lambda_1,\ldots,\lambda_m $$

因此 candidate rejection 时不能简单：

```text
delete mode m+1
```

因为 old modes 的 temporal coordinates 已经发生变化。

所以在任何 tentative post-fixed-point operation 前，应保存：

$$ \mathcal B_m^{\rm before} $$

包括：

- all existing spatial plastic modes；
- all existing stress modes；
- all existing temporal amplitudes；
- all existing temporal rates；
- basis dimension；
- relevant diagnostics / metadata。

---

# 49. rollback 的准确语义

若以下任何情况发生：

- spatial candidate degenerate；
- modified new time function insignificant；
- all-mode temporal solve fails；
- non-finite coefficients；
- full residual benefit insufficient；
- basis-health check fails；

则执行：

$$ \boxed{\mathcal B\leftarrow\mathcal B_m^{\rm before}} $$

即回滚 entire tentative basis state。

不能只执行：

$$ P_{m+1}\rightarrow P_m $$

却保留 modified former temporal functions。

否则 physical reduced solution 不再等于 candidate 进入前的 baseline。

---

# 50. rejection 与 restart 的区别

对于明显 spatial linear dependence：

$$ \gamma_{\rm sp}\approx0 $$

可以直接视为当前 seed / fixed-point candidate 不提供新 direction。

此时可选择：

$$ \boxed{\text{reject current candidate}} $$

或者 future version：

$$ \boxed{\text{restart enrichment with another deterministic seed}} $$

本阶段只冻结：

> 无论 reject 还是 restart，在开始下一 candidate 前都必须恢复 $\mathcal B_m^{\rm before}$。

是否在同一个 global-stage event 内尝试 multiple seeds，可留到 implementation/robustness 阶段。

---

# 51. current 1D global-stage multiple-mode enrichment

current 1D `pgd_global_stage.py` 的基本结构是：

```text
update existing temporal functions

while relative_residual > reduced_tolerance
      and modes_added < max_new_modes:

    enrich one mode

    if accepted:
        append mode

    jointly update all temporal functions

    recompute residual

    possibly continue adding another mode
```

因此一个 current 1D global-stage call 可以连续加入多个 spatial modes。

这是已经验证 1D reproduction 中采用的工程化策略。

---

# 52. 原论文 Fig. 2 对 Add-a-pair 的流程提示

原论文流程图把：

```text
PGD Update
    ↓
saturation decision
    ↓
Add a pair
    ↓
orthonormalise spatial basis
    ↓
correct time functions
    ↓
test new time function
    ↓
build global state / continue LATIN
```

画成一个明确的 pair-addition branch。

这强烈支持一种更 paper-faithful 的 first tower interpretation：

> 一个 Add-a-pair event 首先完成一个 new space-time pair 的生成、正交化、time-function correction 与 accept/reject，然后重新形成完整 LATIN global candidate。

需要注意：

> 这是基于论文流程结构得到的 algorithm interpretation，而不是正文中一句明确写出的“每个 global stage 严禁加入第二个 mode”。

---

# 53. tower v1 建议一次 enrichment event 最多接受一个 pair

第一版 tower LATIN-PGD 的首要目标是：

$$ \boxed{\text{paper correspondence + numerical traceability}} $$

而不是立即达到最大 enrichment speed。

因此建议 tower v1 冻结：

$$ \boxed{\text{one accepted new pair per LATIN enrichment event}} $$

也就是一个 accepted pair 完成后：

- 不立即再次进入 Eq. (61) 产生第二个 mode；
- 先完成当前 global stage；
- 形成 relaxed global candidate；
- 重新计算 outer LATIN convergence / saturation information；
- 下一次是否继续 enrichment 再由 outer control 决定。

---

# 54. 这一选择与 Eq. (60) saturation 的关系

Eq. (60) saturation indicator：

$$ \zeta_i=\frac{\xi_i-\xi_{i+1}}{\xi_i+\xi_{i+1}} $$

依赖 successive LATIN global candidates。

如果一个尚未完成的 global-stage call 内连续加入多个 modes，则：

- 哪一个 intermediate reduced state 对应 $\xi_{i+1}$；
- basis improvement 与 LATIN nonlinear improvement 如何区分；
- mode-level residual decrease 与 outer saturation 如何对应；

都会更难解释。

因此 first tower implementation 采用 one-pair-per-event，有利于保持：

$$ \boxed{\text{LATIN outer iteration}\leftrightarrow\text{PGD enrichment decision}} $$

之间的可追踪关系。

---

# 55. accepted pair 后的 immediate action

若 candidate 通过：

- fixed-point convergence；
- spatial validity；
- modified-time significance；
- all-mode temporal solve；
- residual benefit；

则正式：

$$ m\leftarrow m+1 $$

并将：

$$ \mathcal B_{m+1} $$

设为 persistent working basis。

此时 new pair 的 plastic PGD branch 可以结束。

下一步不是立即添加第二个 mode，而是返回 current LATIN global-stage finishing sequence。

---

# 56. accepted pair 后为什么还不能立刻计算 outer $\xi$

PGD plastic correction 只是 current global stage 的一个组成部分。

原论文 global stage 还包括：

- hardening variables；
- damage variables；
- energy-release variable；
- damage-related structural correction。

因此 accepted new pair 后：

$$ \boxed{\text{plastic branch complete}} $$

并不等于：

$$ \boxed{\text{entire LATIN global state complete}} $$

还必须完成 Eq. (73)–(75) 相关变量与 damage branch，再形成完整 current global state。

之后才适合：

- global relaxation；
- Eq. (76) LATIN indicator；
- Eq. (60) saturation logic。

这些内容留给下一阶段正式闭合。

---

# 57. post-fixed-point 完整 data dependency

本阶段推荐的数据依赖可以写成：

```text
INPUT
    existing accepted basis B_m
    current full forcing f
    current shifted defect Delta_bar
    H_sigma
    M
    reference tower operator
    fixed-point converged raw pair
        p_*
        s_*
        lambda_*
        lambda_dot_*

SNAPSHOT
    save B_m_before

SPATIAL ORTHOGONALISATION
    weighted projection of p_*
    obtain a
    obtain p_perp
    compute gamma_sp

    if spatially degenerate:
        rollback
        reject/restart

NORMALISATION
    c = ||p_perp||_M
    p_(m+1) = p_perp / c
    s_(m+1) = L p_(m+1)

EXACT TEMPORAL COORDINATE TRANSFORMATION
    lambda_j^+ = lambda_j + a_j lambda_*
    lambda_(m+1)^+ = c lambda_*

    lambda_dot_j^+ = lambda_dot_j + a_j lambda_dot_*
    lambda_dot_(m+1)^+ = c lambda_dot_*

FIELD INVARIANCE TEST
    plastic strain unchanged
    plastic strain rate unchanged
    stress correction unchanged
    residual unchanged

MODIFIED-TIME SIGNIFICANCE
    compute gamma_lambda

    if insignificant:
        rollback
        reject

TENTATIVE ENLARGED BASIS
    form B_(m+1)

ALL-MODE TEMPORAL RE-OPTIMISATION
    use full forcing f
    solve Eq. (58)-(59)
    update all temporal amplitudes/rates

FULL RESIDUAL RECONSTRUCTION
    reconstruct R_(m+1)
    compare with R_m

ACTUAL BENEFIT CHECK
    compute Delta_res

    if insufficient or unstable:
        rollback
        reject

ACCEPT
    persist B_(m+1)
    m <- m+1

RETURN TO GLOBAL-STAGE FINISHING
    hardening
    damage
    energy release
    relaxation
    outer LATIN indicator
```

---

# 58. source separation in the algorithm chain

需要再次明确每一步的理论地位。

## 58.1 原论文明确

原论文明确：

```text
new pair fixed point
    ↓
Gram-Schmidt orthonormalisation
    ↓
former time functions updated
    ↓
new time function modified
    ↓
insignificant modified new time function may be rejected
```

## 58.2 数学上由 basis transformation 必然得到

以下关系由 linear coordinate transformation 必然得到：

$$ \vec p_*=P_m\vec a+c\vec p_{m+1} $$

$$ \vec\lambda_m^+=\vec\lambda_m+\vec a\lambda_* $$

$$ \lambda_{m+1}^+=c\lambda_* $$

以及 identical rate transformation。

它们不是 arbitrary engineering heuristic。

## 58.3 tower v1 numerical specification

以下属于 tower v1：

- $M$-weighted modified Gram–Schmidt；
- $\gamma_{\rm sp}$；
- BE-consistent $\|\lambda\|_{T,h}$；
- relative $\gamma_\lambda$；
- exact field-invariance tests；
- all-mode Eq. (58)–(59) re-optimisation；
- full residual benefit $\Delta_{\rm res}$；
- complete rollback；
- one accepted pair per enrichment event；
- baseline 不使用额外 scalar line search。

---

# 59. 与 current 1D implementation 的相同点

tower v1 与 current 1D 保留以下 architecture-level 经验。

第一，PGD mode 仍由：

- spatial plastic-strain function；
- associated equilibrated stress function；
- temporal amplitude；
- temporal rate；

组成。

第二，new mode accepted 后，existing enlarged spatial basis 可固定，再 joint update all temporal functions。

第三，mode acceptance 不应仅依赖 fixed-point convergence。

第四，reduced residual reconstruction 必须作为 independent diagnostic。

这些原则已经在 1D reproduction 中证明具有实际数值价值。

---

# 60. 与 current 1D implementation 的主要差异

tower v1 不直接复制 current 1D 的以下细节。

第一，current 1D 在 inner spatial solve 中已经对 candidate 做 existing-basis orthogonalisation；tower v1 将 Gram–Schmidt 放到 fixed point 收敛之后。

第二，current 1D spatial half-step 是 weighted least-squares；tower v1 采用 paper-derived Eq. (70)–(71) $W$-based FE solve。

第三，current 1D fixed-point change 使用 residual-combination field 且在 half-updated pair 上评价；tower v1 采用 complete mutually consistent pair graph norm。

第四，current 1D post-loop 有 scalar optimal scaling；tower v1 baseline 暂不继承。

第五，current 1D 一个 global-stage call 可加入多个 modes；tower v1 baseline 先采用 one-pair-per-enrichment-event。

第六，tower v1 的 post-fixed-point temporal / residual norm 尽量与 right-endpoint BE contractions 保持一致，而不是混合 trapezoidal integration。

---

# 61. 为什么本阶段暂不赋 $\varepsilon_{\rm mode}$ 与 $\varepsilon_{\rm accept}$ 的具体值

原论文没有给出这两个 tower-discrete thresholds。

若现在直接指定：

$$ \varepsilon_{\rm mode}=10^{-6} $$

或：

$$ \varepsilon_{\rm accept}=10^{-8} $$

缺少足够依据。

后续应通过 diagnostic distributions 来标定，包括：

- accepted candidate 的 $\gamma_{\rm sp}$；
- rejected candidate 的 $\gamma_{\rm sp}$；
- accepted candidate 的 $\gamma_\lambda$；
- residual reduction $\Delta_{\rm res}$；
- final LATIN convergence；
- mode count；
- comparison with FOM；
- sensitivity to time-step size；
- sensitivity to spatial normalization；
- sensitivity to search directions。

因此本阶段只冻结 metric 与 decision architecture，不冻结最终 numerical threshold。

---

# 62. provisional fixed-point tolerance 的地位

上一阶段 current 1D：

$$ \varepsilon_{\rm fp}=10^{-6} $$

可以作为 first tower implementation 的 provisional reference。

但其准确含义是：

> inner rank-one alternating pair 的 relative stability tolerance。

不能与本阶段新增的：

$$ \varepsilon_{\rm mode} $$

或：

$$ \varepsilon_{\rm accept} $$

混淆。

---

# 63. future unit test 1：weighted orthogonality

构造任意已有 $M$-orthonormal basis：

$$ P_m $$

加入一个 raw candidate：

$$ \vec p_* $$

执行 weighted Gram–Schmidt 后检查：

$$ P_{m+1}^TMP_{m+1}\approx I $$

并检查：

$$ \vec p_j^TM\vec p_{m+1}\approx0 $$

对所有：

$$ j=1,\ldots,m $$

成立。

---

# 64. future unit test 2：exact accumulated-plastic invariance

Gram–Schmidt 前 reconstruct：

$$ \Delta\vec\varepsilon^p_{\rm raw}(t) $$

Gram–Schmidt 与 temporal coordinate transformation 后 reconstruct：

$$ \Delta\vec\varepsilon^p_{\rm GS}(t) $$

检查：

$$ \Delta\vec\varepsilon^p_{\rm GS}\approx\Delta\vec\varepsilon^p_{\rm raw} $$

到 floating-point tolerance。

---

# 65. future unit test 3：exact plastic-rate invariance

同样检查：

$$ \Delta\dot{\vec\varepsilon}^p_{\rm GS}\approx\Delta\dot{\vec\varepsilon}^p_{\rm raw} $$

尤其 $t_0$：

$$ \dot\lambda_0 $$

必须通过 direct coordinate transformation 保持一致。

---

# 66. future unit test 4：stress correction invariance

检查：

$$ \Delta\vec\sigma'_{\rm GS}\approx\Delta\vec\sigma'_{\rm raw} $$

同时检查 new stress spatial mode：

$$ \vec s_{m+1}=C_0(\mathcal E_{\rm tower}-I)\vec p_{m+1} $$

并满足 tower weak equilibrium。

---

# 67. future unit test 5：mechanical residual invariance before re-optimisation

在 full forcing 与 search direction 不变时：

$$ R_{\rm GS}=R_{\rm raw} $$

应成立。

如果在 all-mode re-optimisation 之前 residual 已改变，则 implementation 必有问题。

---

# 68. future unit test 6：nearly dependent candidate rejection

构造：

$$ \vec p_*=\vec p_1+\epsilon\vec r $$

其中 $\epsilon$ 很小。

应得到：

$$ \gamma_{\rm sp}\ll1 $$

并触发 spatial-dependence safeguard。

同时 rejection 后应验证：

$$ \mathcal B_{\rm after}=\mathcal B_{\rm before} $$

---

# 69. future unit test 7：insignificant modified-time rejection

构造 spatially independent candidate：

$$ \gamma_{\rm sp}=O(1) $$

但 temporal amplitude 很小：

$$ \|\lambda_*\|_{T,h}\ll1 $$

应得到：

$$ \gamma_\lambda\ll1 $$

并触发 insignificant-mode rejection。

---

# 70. future unit test 8：complete rollback

candidate tentative append 后主动制造 temporal solve failure 或 residual acceptance failure。

检查 rejection 后：

- mode count unchanged；
- existing spatial modes unchanged；
- existing stress modes unchanged；
- existing amplitudes unchanged；
- existing rates unchanged；
- reconstructed correction unchanged。

这将验证 rollback 是 transactional，而不是 partial delete。

---

# 71. future unit test 9：all-mode re-optimisation uses full forcing

通过一个 manufactured test 区分：

$$ f $$

与：

$$ \bar\Delta $$

确保 post-acceptance Eq. (58)–(59) solver 的 forcing interface 使用 full global-stage plastic forcing，而不是 enrichment residual defect。

这项测试可以防止未来代码中两类 residual 的语义混用。

---

# 72. future unit test 10：one-pair-per-event

tower v1 第一版应检查：

- 一个 enrichment event 最多使 basis dimension 增加 1；
- accepted pair 后函数返回 global-stage finishing flow；
- 不在同一个 event 内自动进入第二次 Eq. (61) candidate construction。

将来若添加 accelerated multiple-mode option，应作为显式 optional strategy，而不是悄然改变 baseline。

---

# 73. 本阶段冻结的 post-fixed-point algorithm

本阶段建议正式冻结以下 baseline。

```text
Given:
    existing accepted basis B_m
    fixed-point converged raw pair P_*
    full forcing f
    shifted defect Delta_bar
    current H_sigma and material-point metric M

A. Save complete basis snapshot B_m_before.

B. Weighted modified Gram-Schmidt:
       p_* -> projection coefficients a
       p_perp = p_* - P_m a

C. Spatial novelty diagnostic:
       gamma_sp = ||p_perp||_M / ||p_*||_M

       if degenerate:
           rollback
           reject/restart

D. Normalize:
       c = ||p_perp||_M
       p_(m+1) = p_perp / c
       s_(m+1) = L p_(m+1)

E. Exact temporal coordinate transformation:
       lambda_m^+ = lambda_m + a lambda_*
       lambda_(m+1)^+ = c lambda_*

       same transformation for temporal rates

F. Exact invariance checks:
       accumulated plastic correction unchanged
       plastic-rate correction unchanged
       stress correction unchanged
       mechanical residual unchanged

G. Modified-new-time significance:
       compute gamma_lambda

       if insignificant:
           rollback
           reject

H. Tentatively form enlarged basis B_(m+1).

I. Eq. (58)-(59) all-mode temporal re-optimisation:
       fixed enlarged spatial basis
       use full forcing f
       solve all m+1 temporal functions

J. Reconstruct whole mechanical residual R_(m+1).

K. Residual benefit:
       Delta_res = 1 - ||R_(m+1)|| / ||R_m||

       if insufficient or unstable:
           rollback
           reject

L. Accept:
       persist B_(m+1)
       m <- m+1

M. End this enrichment event.

N. Return to global-stage finishing:
       hardening
       damage
       energy release
       relaxation
       outer LATIN indicator
```

---

# 74. 本阶段核心结论

本阶段最重要的结论可以浓缩为以下十点。

第一，fixed-point converged raw pair 不能直接 append，必须先对 spatial mode 做 material-point $M$-weighted Gram–Schmidt。

第二，Gram–Schmidt 后 existing temporal functions 的修改不是经验操作，而是保持 physical field 不变的 exact basis-coordinate transformation。

第三，new orthogonal spatial mode normalization 系数：

$$ c=\|\vec p_\perp\|_M $$

必须同步乘到 new temporal amplitude 与 rate 上。

第四，same linear transformation 自动保持 associated equilibrated stress correction。

第五，Gram–Schmidt + temporal correction 本身不应改变 accumulated plastic strain、plastic rate、stress correction 或 mechanical residual。

第六，原论文 modified new time function insignificant 的说法，在 tower 中可通过 BE-consistent relative temporal significance $\gamma_\lambda$ 离散化。

第七，$\gamma_{\rm sp}$ 衡量 spatial novelty，$\gamma_\lambda$ 衡量真正 orthogonal new mode 的 space-time significance，两者不能混用。

第八，exact basis transformation 后，建议在 enlarged spatial basis 上重新执行 Eq. (58)–(59) all-mode temporal update，并使用 full global-stage plastic forcing $f$。

第九，candidate 最终还应经过 full residual benefit check，并且 rejection 必须完整 rollback 到 candidate 进入前的 basis snapshot。

第十，tower v1 baseline 建议一次 LATIN enrichment event 最多接受一个 pair，并暂不继承 current 1D post-loop scalar optimal line search。

---

# 75. Eq. (61)–(72) Add-a-pair 链条的当前完整状态

经过前序多个阶段，目前 Add-a-pair branch 已可以完整写为：

$$ \boxed{\text{Eq. (60): existing PGD basis judged insufficient}} $$

$$ \downarrow $$

$$ \boxed{\text{Eq. (61)–(64): shifted defect + rank-one enrichment ansatz}} $$

$$ \downarrow $$

$$ \boxed{\text{Eq. (65)–(71): fixed temporal function → spatial Galerkin solve}} $$

$$ \downarrow $$

$$ \boxed{\text{Eq. (72): fixed spatial mode → temporal residual minimisation}} $$

$$ \downarrow $$

$$ \boxed{\text{Eq. (70)–(72): alternating fixed-point loop}} $$

$$ \downarrow $$

$$ \boxed{\text{complete-pair fixed-point convergence}} $$

$$ \downarrow $$

$$ \boxed{\text{$M$-weighted Gram--Schmidt}} $$

$$ \downarrow $$

$$ \boxed{\text{exact temporal coordinate transformation}} $$

$$ \downarrow $$

$$ \boxed{\text{modified new temporal function significance}} $$

$$ \downarrow $$

$$ \boxed{\text{enlarged-basis Eq. (58)–(59) temporal re-optimisation}} $$

$$ \downarrow $$

$$ \boxed{\text{full mechanical residual benefit check}} $$

$$ \downarrow $$

$$ \boxed{\text{accept persistent mode or complete rollback}} $$

因此到本阶段结束：

> Eq. (61)–(72) 从“产生一个 residual-driven new pair”到“该 pair 正式进入 persistent PGD basis”的核心算法链已经基本闭合。

剩余问题主要属于：

- threshold calibration；
- implementation diagnostics；
- unit tests；
- solver robustness variants；

而不再是 Add-a-pair 主结构缺失。

---

# 76. 本阶段之后不应立即做的事情

虽然 Add-a-pair 主结构已经闭合，但当前仍不建议立即开始 tower LATIN-PGD 完整代码实现。

原因是整个 LATIN global stage 仍有一个关键 finishing block 尚未完全闭合：

- hardening Eq. (73)–(74) 的准确 update order；
- damage Eq. (75) 的准确 update order；
- $D$ history 到底复制 local integrated damage 还是复制 rate 后 reintegrate；
- damage-dependent structural correction 与 plastic PGD correction 的最终 recombination；
- final energy-release variable；
- relaxation 到底作用于哪些 global variables；
- relaxation 前后 Eq. (76) indicator 应使用哪个 state；
- outer $\xi$、Eq. (60) saturation 与 one-pair enrichment event 如何在程序控制层连接。

这些问题直接决定 global-state update interface，因此应在写 tower solver 代码前继续闭合。

---

# 77. 下一阶段建议

下一阶段建议一次性处理以下紧密相关内容：

1. accepted new pair 后 plastic PGD branch 的结束状态；
2. Eq. (73)–(74) hardening variables 的 global update；
3. Eq. (75) damage / energy-release variables 的 global update；
4. current paper structure 与 1D implementation 对 damage history 的差异；
5. tower v1 的 damage-history inheritance choice；
6. plastic correction、damage correction、hardening correction 的 final global-state assembly；
7. global relaxation 的位置；
8. Eq. (76) LATIN indicator 计算所使用的 relaxed/unrelaxed states；
9. Eq. (60) saturation 与下一次 enrichment event 的 outer control；
10. complete local → global → relaxation → convergence data-flow specification。

这一阶段完成后，tower LATIN-PGD 理论层面就接近从 Eq. (47) 到 Eq. (77) 的完整 algorithm specification。

---

# 78. 当前准确研究停点

截至本阶段结束，当前准确停点为：

$$ \boxed{\text{new pair accepted into persistent PGD basis}} $$

下一步应从：

$$ \boxed{\text{Eq. (73)–(75) global-stage finishing}} $$

开始，而不是继续添加第二个 PGD mode，也不是开始写 tower enrichment code。

这可以保证后续代码实现时：

- plastic reduced branch；
- hardening local-global branch；
- damage full-order branch；
- relaxation；
- outer LATIN convergence；

五个层级一次性具有清晰的数据边界和更新顺序。
