# Tower LATIN-PGD Eq. (72) backward-Euler 离散、标量 temporal update 与 whole-time minimisation 边界阶段总结

**日期：2026-08-17**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**当前研究路线：Bhattacharyya et al. 原论文 $x-t$ LATIN-PGD → 2D fiber beam-column offshore wind turbine tower**  
**阶段范围：自上一份 Eq. (72) temporal minimisation 总结之后，完成 single-new-mode backward-Euler residual 离散、$H_\sigma^{-1}$ weighted scalar least-squares、tower $(e,g,f)$ 展开、whole-time minimisation 与 sequential time marching 的严格区分、现有 1D `pgd_time_update.py` 的数学对应及后续初始时间点问题定位**  
**上一阶段衔接：`2026-08-17-tower-latin-pgd-eq72-temporal-minimization-discretization-stage-summary.md`**  
**下一阶段：闭合第一个时间点 $t_0$ 的 $\lambda_0$、$\dot{\lambda}_0$ 与 initial temporal treatment，再决定 tower temporal solver 的完整离散接口**

---

# 1. 本阶段定位

上一阶段已经完成原论文 Eq. (72) 的连续层面解释，并建立了以下共识：

- Eq. (72) 来自 Eq. (61) 的 remaining LATIN descent search-direction mechanical residual；
- fixed spatial mode 后，Eq. (72) 只更新 new pair 的 temporal function $\lambda(t)$；
- minimisation metric 为 $H_\sigma^{-1}$；
- tower continuum spatial coordinate $x$ 映射为 fiber material-point index $q=(e,g,f)$；
- 原论文明确使用 DG0 处理 temporal minimisation，但 2018 原文没有给出足够完整的 DG0 algebra；
- 当前 tower migration 第一版暂时继承已经在 1D reproduction 中验证的 backward-Euler temporal implementation；
- 不应把当前 backward-Euler implementation 声称为 exact original DG0。

上一阶段尚未闭合的问题是：

$$ \lambda(t) \rightarrow \lambda_n $$

即：如何把 Eq. (72) 在 single-new-mode tower 情形下真正离散成可计算的 temporal update。

本阶段完成的主要工作即围绕这一问题展开。

---

# 2. 本阶段核心结论

本阶段得到 single-new-mode backward-Euler residual：

$$ \vec{r}_n = \left( \frac{\vec{p}}{\Delta t_n} - D_{H,n}\vec{s} \right)\lambda_n - \frac{\vec{p}}{\Delta t_n}\lambda_{n-1} + \vec{\Delta}_n $$

定义：

$$ \vec{g}_n = \frac{\vec{p}}{\Delta t_n} - D_{H,n}\vec{s} $$

以及：

$$ \vec{b}_n = \frac{\vec{p}}{\Delta t_n}\lambda_{n-1} - \vec{\Delta}_n $$

则：

$$ \vec{r}_n = \vec{g}_n\lambda_n - \vec{b}_n $$

当前 time step 的 weighted objective 为：

$$ J_n = \frac{1}{2}\vec{r}_n^T M D_{H,n}^{-1}\vec{r}_n $$

最优条件：

$$ \frac{dJ_n}{d\lambda_n} = \vec{g}_n^T M D_{H,n}^{-1}\vec{r}_n = 0 $$

得到 single-new-mode scalar normal equation：

$$ \left( \vec{g}_n^T M D_{H,n}^{-1}\vec{g}_n \right)\lambda_n = \vec{g}_n^T M D_{H,n}^{-1}\left( \frac{\vec{p}}{\Delta t_n}\lambda_{n-1} - \vec{\Delta}_n \right) $$

若 denominator 非零，则：

$$ \lambda_n = \frac{\vec{g}_n^T M D_{H,n}^{-1}\left( \frac{\vec{p}}{\Delta t_n}\lambda_{n-1} - \vec{\Delta}_n \right)}{\vec{g}_n^T M D_{H,n}^{-1}\vec{g}_n} $$

但本阶段进一步确认：

**该公式是 sequential one-step weighted least-squares 的精确解，而不是 whole-time Eq. (72) 离散后的全局精确解。**

这是本阶段最重要的理论边界修正。

---

# 3. 记号统一

为避免与 Eq. (70)–(71) spatial half-step 中的临时量混淆，本阶段采用：

$$ \vec{p} = \vec{\bar{\varepsilon}}^p_{m+1} $$

表示 Eq. (71) 得到的新 plastic-strain spatial mode。

采用：

$$ \vec{s} = \vec{\bar{\sigma}}_{m+1} $$

表示对应的 equilibrated stress spatial mode。

采用 $D_H(t)$ 表示 tower fiber material-point space 中的 diagonal search-direction operator，其对角项为：

$$ D_H(t) = {\rm diag}\left( H_{\sigma,q}(t) \right) $$

采用：

$$ \vec{\Delta}(t) = \vec{\bar{\Delta}}_{i+1}(t) $$

表示 existing PGD update 后仍未消除的 shifted LATIN defect。

tower material-point quadrature metric 记为：

$$ M = {\rm diag}(v_q) $$

其中：

$$ v_q = A_{egf} w_g J_e $$

且：

$$ q = (e,g,f) $$

---

# 4. Eq. (72) continuous residual 的离散起点

上一阶段已经得到 tower vector form：

$$ \vec{r}(t) = \dot{\lambda}(t)\vec{p} - D_H(t)\lambda(t)\vec{s} + \vec{\Delta}(t) $$

这一 residual 直接来自原论文 Eq. (61)：

$$ \Delta\dot{\varepsilon}^p - H_\sigma \Delta\sigma' + \bar{\Delta} = 0 $$

其中 new rank-one pair：

$$ \Delta\dot{\varepsilon}^p = \dot{\lambda}\bar{\varepsilon}^p $$

以及：

$$ \Delta\sigma' = \lambda\bar{\sigma} $$

因此本阶段的 temporal discretisation 不重新定义 residual，只处理 $\dot{\lambda}(t)$ 的离散。

---

# 5. backward-Euler temporal derivative

设时间网格为：

$$ 0 = t_0 < t_1 < \cdots < t_n < \cdots < t_{N_t-1} = T $$

第 $n$ 个时间增量：

$$ \Delta t_n = t_n - t_{n-1} $$

继承当前 1D implementation 的 backward-Euler temporal derivative：

$$ \dot{\lambda}_n = \frac{\lambda_n - \lambda_{n-1}}{\Delta t_n} $$

这一步建立了当前 temporal amplitude $\lambda_n$ 与上一时刻 $\lambda_{n-1}$ 的因果联系。

---

# 6. backward Euler 代入 Eq. (72) residual

在 $t=t_n$：

$$ D_H(t_n) = D_{H,n} $$

$$ \vec{\Delta}(t_n) = \vec{\Delta}_n $$

因此：

$$ \vec{r}_n = \frac{\lambda_n-\lambda_{n-1}}{\Delta t_n}\vec{p} - D_{H,n}\lambda_n\vec{s} + \vec{\Delta}_n $$

这一式子是本阶段 temporal discretisation 的第一条核心公式。

它仍然只是 residual，并未执行 minimisation。

---

# 7. 将当前 unknown 单独整理

展开：

$$ \vec{r}_n = \frac{\vec{p}}{\Delta t_n}\lambda_n - \frac{\vec{p}}{\Delta t_n}\lambda_{n-1} - D_{H,n}\vec{s}\lambda_n + \vec{\Delta}_n $$

将所有含 $\lambda_n$ 的项组合：

$$ \vec{r}_n = \left( \frac{\vec{p}}{\Delta t_n} - D_{H,n}\vec{s} \right)\lambda_n - \frac{\vec{p}}{\Delta t_n}\lambda_{n-1} + \vec{\Delta}_n $$

定义：

$$ \vec{g}_n = \frac{\vec{p}}{\Delta t_n} - D_{H,n}\vec{s} $$

因此：

$$ \vec{r}_n = \vec{g}_n\lambda_n - \frac{\vec{p}}{\Delta t_n}\lambda_{n-1} + \vec{\Delta}_n $$

---

# 8. $\vec{g}_n$ 的数学意义

由 residual 对当前 temporal unknown 的导数：

$$ \frac{\partial\vec{r}_n}{\partial\lambda_n} = \vec{g}_n $$

因此 $\vec{g}_n$ 是当前 temporal amplitude $\lambda_n$ 对整个 material-point residual field 的 sensitivity direction。

改变 $\lambda_n$ 一个单位，tower fiber material-point residual 沿 $\vec{g}_n$ 方向变化。

---

# 9. $\vec{g}_n$ 的组成

有：

$$ \vec{g}_n = \frac{\vec{p}}{\Delta t_n} - D_{H,n}\vec{s} $$

第一项来源于 plastic-rate correction。

第二项来源于 stress contribution in the LATIN search direction。

因此 $\vec{g}_n$ 反映 temporal amplitude 对两个 residual contributions 的综合影响。

---

# 10. 定义已知 right-hand-side vector

定义：

$$ \vec{b}_n = \frac{\vec{p}}{\Delta t_n}\lambda_{n-1} - \vec{\Delta}_n $$

于是：

$$ \vec{r}_n = \vec{g}_n\lambda_n - \vec{b}_n $$

在 sequential one-step formulation 中：

- $\vec{g}_n$ 已知；
- $\vec{b}_n$ 已知；
- $\lambda_n$ 为唯一 scalar unknown。

因此当前问题退化为 single-variable weighted least-squares problem。

---

# 11. tower 的 $H_\sigma^{-1}$ spatial metric

原论文 Eq. (72) 使用 $H_\sigma^{-1}$ weighted residual norm。

在 tower material-point discretisation 中：

$$ \int_\Omega r H_\sigma^{-1} r \, d\Omega \rightarrow \vec{r}^T M D_H^{-1}\vec{r} $$

其中 $M$ 存储 material-point quadrature volume，$D_H^{-1}$ 存储每个 material point 的 $H_\sigma^{-1}$。

因此 time step $n$ 的 spatial weighting operator 为：

$$ W_{H,n} = M D_{H,n}^{-1} $$

这里 $W_{H,n}$ 只是 temporal least-squares weighting operator，不应与 Eq. (67)–(70) 中的 spatial effective operator $W$ 混淆。

---

# 12. 第 $n$ 个时间步的 weighted objective

定义：

$$ J_n = \frac{1}{2}\vec{r}_n^T M D_{H,n}^{-1}\vec{r}_n $$

代入：

$$ \vec{r}_n = \vec{g}_n\lambda_n - \vec{b}_n $$

得到：

$$ J_n = \frac{1}{2}\left( \vec{g}_n\lambda_n-\vec{b}_n \right)^T M D_{H,n}^{-1}\left( \vec{g}_n\lambda_n-\vec{b}_n \right) $$

该 objective 是原论文 Eq. (72) residual/metric 结构在当前 sequential one-step discretisation 下的对应形式。

---

# 13. 为什么不能要求所有 fiber residual 都严格为零

当前 $\lambda_n$ 只有一个 scalar unknown，而 tower 中存在大量 material points：

$$ q = (e,g,f) $$

因此一般无法同时满足：

$$ r_{n,q} = 0 $$

对所有 material points 都成立。

所以 Eq. (72) 采用的是 weighted best fit over the full material-point space，而不是 pointwise exact satisfaction。

---

# 14. 对 $\lambda_n$ 求导

由于：

$$ \frac{\partial\vec{r}_n}{\partial\lambda_n} = \vec{g}_n $$

因此：

$$ \frac{dJ_n}{d\lambda_n} = \vec{g}_n^T M D_{H,n}^{-1}\vec{r}_n $$

最优 temporal amplitude 满足：

$$ \frac{dJ_n}{d\lambda_n} = 0 $$

所以：

$$ \vec{g}_n^T M D_{H,n}^{-1}\vec{r}_n = 0 $$

---

# 15. weighted orthogonality condition

上式说明：

$$ \vec{g}_n^T M D_{H,n}^{-1}\vec{r}_n = 0 $$

即最优 residual 与当前唯一可调 residual direction $\vec{g}_n$ 在 $H_\sigma^{-1}$ weighted metric 下正交。

这是 weighted least-squares 的标准 orthogonality condition。

---

# 16. single-new-mode scalar normal equation

代入：

$$ \vec{r}_n = \vec{g}_n\lambda_n - \vec{b}_n $$

得到：

$$ \vec{g}_n^T M D_{H,n}^{-1}\left( \vec{g}_n\lambda_n-\vec{b}_n \right) = 0 $$

展开：

$$ \left( \vec{g}_n^T M D_{H,n}^{-1}\vec{g}_n \right)\lambda_n = \vec{g}_n^T M D_{H,n}^{-1}\vec{b}_n $$

恢复 $\vec{b}_n$：

$$ \left( \vec{g}_n^T M D_{H,n}^{-1}\vec{g}_n \right)\lambda_n = \vec{g}_n^T M D_{H,n}^{-1}\left( \frac{\vec{p}}{\Delta t_n}\lambda_{n-1}-\vec{\Delta}_n \right) $$

这是本阶段得到的 single-new-mode scalar normal equation。

---

# 17. explicit scalar update

若：

$$ \vec{g}_n^T M D_{H,n}^{-1}\vec{g}_n > 0 $$

则：

$$ \lambda_n = \frac{\vec{g}_n^T M D_{H,n}^{-1}\left( \frac{\vec{p}}{\Delta t_n}\lambda_{n-1}-\vec{\Delta}_n \right)}{\vec{g}_n^T M D_{H,n}^{-1}\vec{g}_n} $$

这一公式是 current sequential one-step weighted least-squares problem 的 exact scalar minimiser。

---

# 18. denominator 的意义

定义：

$$ A_n = \vec{g}_n^T M D_{H,n}^{-1}\vec{g}_n $$

则 $A_n$ 是 $\vec{g}_n$ 在当前 weighted metric 下的平方范数。

因此 $A_n$ 衡量当前 new spatial mode 在 time step $n$ 对 LATIN residual 的有效 weighted sensitivity magnitude。

若 $A_n$ 很小，则当前 temporal DOF 对 residual 的控制方向很弱，可能出现 conditioning 问题。

---

# 19. numerator 的意义

定义：

$$ c_n = \vec{g}_n^T M D_{H,n}^{-1}\left( \frac{\vec{p}}{\Delta t_n}\lambda_{n-1}-\vec{\Delta}_n \right) $$

则：

$$ A_n\lambda_n = c_n $$

$c_n$ 表示 current history/defect forcing 在 new-mode controllable residual direction 上的 weighted projection。

因此：

$$ \lambda_n = \frac{c_n}{A_n} $$

---

# 20. tower $(e,g,f)$ 形式的 $\vec{g}_n$

对于 material point：

$$ q = (e,g,f) $$

有：

$$ p_q = \bar{\varepsilon}^p_{egf} $$

$$ s_q = \bar{\sigma}_{egf} $$

$$ H_{\sigma,n,q} = H_{\sigma,egf}(t_n) $$

因此：

$$ g_{n,egf} = \frac{\bar{\varepsilon}^p_{egf}}{\Delta t_n} - H_{\sigma,egf}(t_n)\bar{\sigma}_{egf} $$

---

# 21. tower $(e,g,f)$ 形式的 residual

material point $(e,g,f)$ 上：

$$ r_{n,egf} = g_{n,egf}\lambda_n - \frac{\bar{\varepsilon}^p_{egf}}{\Delta t_n}\lambda_{n-1} + \bar{\Delta}_{egf}(t_n) $$

完整展开：

$$ r_{n,egf} = \left( \frac{\bar{\varepsilon}^p_{egf}}{\Delta t_n} - H_{\sigma,egf}(t_n)\bar{\sigma}_{egf} \right)\lambda_n - \frac{\bar{\varepsilon}^p_{egf}}{\Delta t_n}\lambda_{n-1} + \bar{\Delta}_{egf}(t_n) $$

这就是 Eq. (72) 从 continuum $x$ 到 tower fiber material-point space 的 direct discrete correspondence。

---

# 22. tower denominator 的显式求和

由：

$$ v_{egf} = A_{egf}w_gJ_e $$

可得：

$$ A_n = \sum_{e,g,f}\frac{A_{egf}w_gJ_e}{H_{\sigma,egf}(t_n)}g_{n,egf}^2 $$

进一步代入 $g_{n,egf}$：

$$ A_n = \sum_{e,g,f}\frac{A_{egf}w_gJ_e}{H_{\sigma,egf}(t_n)}\left( \frac{\bar{\varepsilon}^p_{egf}}{\Delta t_n} - H_{\sigma,egf}(t_n)\bar{\sigma}_{egf} \right)^2 $$

---

# 23. tower numerator 的显式求和

右端 scalar：

$$ c_n = \sum_{e,g,f}\frac{A_{egf}w_gJ_e}{H_{\sigma,egf}(t_n)}g_{n,egf}\left( \frac{\bar{\varepsilon}^p_{egf}}{\Delta t_n}\lambda_{n-1} - \bar{\Delta}_{egf}(t_n) \right) $$

因此 explicit scalar update 可以在 tower material-point loop 中通过 scalar accumulation 完成。

---

# 24. Eq. (72) temporal solve 不需要重新进行 tower FE structural solve

在 Eq. (70)–(71) spatial half-step 后：

- $\vec{p}$ 已经确定；
- $\vec{s}$ 已经确定；
- $D_{H,n}$ 已知；
- $\vec{\Delta}_n$ 已知。

Eq. (72) temporal half-step 不再求 nodal displacement，也不再构造新的 tower equilibrium equation。

所有 spatial information 被 contraction 为 temporal scalar coefficients。

因此 Eq. (72) 是 reduced temporal problem，而不是新的 tower FE equilibrium solve。

---

# 25. 与现有 1D `pgd_time_update.py` 的对应

当前 1D code 中 residual 为：

$$ r = P\dot{\lambda} - H_\sigma S\lambda - f $$

backward-Euler reduced matrix 为：

$$ G_n = \frac{P}{\Delta t_n} - H_{\sigma,n}S $$

single-mode 时：

$$ P \rightarrow \vec{p} $$

$$ S \rightarrow \vec{s} $$

因此：

$$ G_n \rightarrow \vec{g}_n $$

代码中的 right-hand side 为：

$$ f_n + \frac{P}{\Delta t_n}\lambda_{n-1} $$

而当前 Eq. (72) sign convention 为：

$$ f_n = -\vec{\Delta}_n $$

所以：

$$ f_n + \frac{\vec{p}}{\Delta t_n}\lambda_{n-1} = \frac{\vec{p}}{\Delta t_n}\lambda_{n-1} - \vec{\Delta}_n $$

与本阶段推导一致。

---

# 26. 当前 1D code 的 weighting

1D `pgd_time_update.py` 中每个 time step 使用：

$$ w_{n,e} = \frac{V_e}{H_{\sigma,n,e}} $$

tower 对应为：

$$ w_{n,egf} = \frac{A_{egf}w_gJ_e}{H_{\sigma,egf}(t_n)} $$

因此现有 1D weighted least-squares kernel 在数学结构上可以自然迁移到 tower material-point space。

---

# 27. 一个必须修正的表述

本阶段早期曾暂时采用以下直观说法：

> backward Euler 使 Eq. (72) 变为逐时间步 scalar minimisation。

经过进一步分析，这一表述不够严格。

更准确的说法是：

> 当前项目选择采用 causal sequential backward-Euler one-step weighted least-squares implementation。

原因是：若先离散整个 Eq. (72) whole-time objective，再对所有 $\lambda_n$ 联合最小化，则相邻时间节点仍然耦合。

---

# 28. whole-time discrete objective

若从连续 Eq. (72) 的 whole-time objective 出发，并采用 backward-Euler residual，则可写为：

$$ J_h = \frac{1}{2}\sum_{n=1}^{N_t-1}\omega_n\vec{r}_n^T W_n\vec{r}_n $$

其中：

$$ W_n = M D_{H,n}^{-1} $$

$\omega_n$ 为 temporal quadrature weight。

此时：

$$ \vec{r}_n = \vec{g}_n\lambda_n - \frac{\vec{p}}{\Delta t_n}\lambda_{n-1} + \vec{\Delta}_n $$

---

# 29. 为什么 whole-time minimisation 会产生时间节点耦合

对于内部节点 $\lambda_k$：

$$ \lambda_k \rightarrow \vec{r}_k $$

同时：

$$ \lambda_k \rightarrow \vec{r}_{k+1} $$

因为：

$$ \vec{r}_{k+1} = \vec{g}_{k+1}\lambda_{k+1} - \frac{\vec{p}}{\Delta t_{k+1}}\lambda_k + \vec{\Delta}_{k+1} $$

因此一个 temporal amplitude 不只影响当前 residual，还影响下一 time slab 的 backward-difference residual。

---

# 30. whole-time objective 对 $\lambda_k$ 的导数

对于内部时间节点：

$$ \frac{\partial\vec{r}_k}{\partial\lambda_k} = \vec{g}_k $$

以及：

$$ \frac{\partial\vec{r}_{k+1}}{\partial\lambda_k} = -\frac{\vec{p}}{\Delta t_{k+1}} $$

所以：

$$ \frac{\partial J_h}{\partial\lambda_k} = \omega_k\vec{g}_k^T W_k\vec{r}_k - \omega_{k+1}\left( \frac{\vec{p}}{\Delta t_{k+1}} \right)^T W_{k+1}\vec{r}_{k+1} $$

whole-time optimality condition 为：

$$ \omega_k\vec{g}_k^T W_k\vec{r}_k - \omega_{k+1}\left( \frac{\vec{p}}{\Delta t_{k+1}} \right)^T W_{k+1}\vec{r}_{k+1} = 0 $$

---

# 31. whole-time discrete system 的结构

由于 $\vec{r}_k$ 包含 $\lambda_{k-1}$ 和 $\lambda_k$，而 $\vec{r}_{k+1}$ 包含 $\lambda_k$ 和 $\lambda_{k+1}$，所以 whole-time optimality equation 对内部节点一般同时耦合：

$$ \lambda_{k-1},\lambda_k,\lambda_{k+1} $$

single-mode 情况下会形成类似 time-tridiagonal algebraic system。

因此：

**whole-time backward-Euler discretisation 不等于逐步 one-step scalar minimisation。**

---

# 32. whole-time coupling 与连续 Euler-Lagrange equation 的一致性

上一阶段从 continuous Eq. (72) 第一变分得到：

$$ \frac{d}{dt}\left( A\dot{\lambda}+d \right) - C\lambda + e = 0 $$

展开：

$$ A\ddot{\lambda} + \dot{A}\dot{\lambda} - C\lambda + \dot{d} + e = 0 $$

连续问题包含 $\ddot{\lambda}$。

因此时间离散后出现相邻 time nodes coupling 是自然的。

这一点说明，不能仅凭 backward Euler 就声称 whole-time minimisation 自动退化为单向 sequential scalar updates。

---

# 33. sequential one-step formulation 的准确数学地位

本阶段 explicit scalar formula：

$$ \lambda_n = \frac{\vec{g}_n^T W_n\left( \frac{\vec{p}}{\Delta t_n}\lambda_{n-1}-\vec{\Delta}_n \right)}{\vec{g}_n^T W_n\vec{g}_n} $$

应准确理解为：

> 在固定 $\lambda_{n-1}$ 后，仅对当前 $J_n$ 最小化得到的 exact one-step solution。

它不是所有 temporal amplitudes 联合最小化 whole-time objective 的 global solution。

---

# 34. sequential scheme 的因果结构

current implementation 的算法逻辑是：

$$ \lambda_{n-1} \rightarrow \lambda_n \rightarrow \lambda_{n+1} $$

即：

- past fixed；
- optimise present；
- propagate to future。

这种处理仍然保留 history dependence，但不会为了降低 future residual 再回头调整已经接受的 past amplitude。

---

# 35. sequential 与 whole-time minimisation 的核心区别

sequential route：

$$ \lambda_n = {\rm argmin}\ J_n $$

且按 $n=1,2,\ldots$ 顺序进行。

whole-time route：

$$ \lambda_0,\lambda_1,\ldots,\lambda_{N_t-1} = {\rm argmin}\ J_h $$

前者是 causal time marching。

后者是 global-in-time optimisation。

两者 residual structure 相同，但 optimisation scope 不同。

---

# 36. 为什么当前 tower migration 暂时保留 sequential route

本项目当前目标仍然是先验证 original $x-t$ LATIN-PGD spatial architecture 能否迁移到 offshore tower fiber space。

当前 1D reproduction 已经使用 sequential backward-Euler weighted least-squares temporal update 并完成验证。

因此 tower v1 保持该 temporal implementation，可以：

- 避免同时修改 spatial migration 与 temporal solver；
- 保留与现有 1D reference implementation 的直接可比性；
- 将算法差异控制在 tower spatial discretisation 本身。

因此当前研究决策仍为：

> tower v1 继承现有 1D sequential backward-Euler temporal kernel。

---

# 37. 不能对原论文 DG0 做过度声明

原论文明确说明 temporal minimisation 使用 zero-order discontinuous Galerkin method。

然而当前掌握的 2018 原论文正文并未给出完整的：

- DG0 trial/test space；
- jump terms；
- left/right traces；
- first slab；
- endpoint treatment；
- final algebraic system。

因此目前不能证明：

$$ {\rm current\ sequential\ backward\ Euler} = {\rm original\ exact\ DG0\ algebra} $$

后续代码或论文中不应声称 exact equivalence。

---

# 38. 推荐的严谨方法描述

当前 implementation 建议描述为：

> a sequential backward-Euler discretised $H_\sigma^{-1}$ weighted temporal residual minimisation consistent with the residual structure of Eq. (72)

中文可表述为：

> 基于原论文 Eq. (72) 的 $H_\sigma^{-1}$ 加权机械残差结构，采用与已验证 1D 实现一致的顺序 backward-Euler temporal least-squares 离散。

这一表述明确区分：

- original continuous LATIN-PGD structure；
- current project temporal discretisation choice。

---

# 39. 当前 `pgd_time_update.py` 的证据边界

现有 1D `latin/pgd_time_update.py` 中已经明确：

- residual 为 $P\dot{\lambda}-H_\sigma S\lambda-f$；
- temporal coefficients minimise the $H_\sigma^{-1}$ weighted mechanical residual；
- temporal derivative 使用 backward Euler；
- 每一个 time step 构造 weighted least-squares problem。

但 current docstring 将 backward Euler 描述为 DG0 的 discrete counterpart。

结合本阶段理论边界，后续进入代码整理阶段时，可以考虑收紧这一措辞。

本阶段暂不修改代码。

---

# 40. single-mode 与 multi-mode 的关系

本阶段集中于 new enrichment pair，所以每个 time step 只有一个 scalar unknown：

$$ \lambda_n $$

如果 fixed spatial basis 中已有 $m$ 个 modes，则：

$$ \vec{\lambda}_n \in R^m $$

对应 reduced matrix：

$$ G_n = \frac{P}{\Delta t_n} - D_{H,n}S $$

每个 time step 的 sequential weighted least-squares 变为一个 $m$ 维 reduced problem。

所以 single-new-mode scalar update 是 current multi-mode time update 的一维特例。

---

# 41. computational implication for tower

single-new-mode sequential Eq. (72) 不需要构造 dense material-point matrices。

实际数学操作可写为：

$$ A_n = \sum_q w_{n,q} g_{n,q}^2 $$

$$ c_n = \sum_q w_{n,q} g_{n,q} b_{n,q} $$

其中：

$$ w_{n,q} = \frac{v_q}{H_{\sigma,n,q}} $$

然后：

$$ \lambda_n = \frac{c_n}{A_n} $$

因此 spatial dimension $N_q$ 只进入 scalar contraction。

---

# 42. denominator positivity 与 conditioning

若：

$$ v_q > 0 $$

且：

$$ H_{\sigma,n,q} > 0 $$

则：

$$ A_n = \sum_q \frac{v_q}{H_{\sigma,n,q}} g_{n,q}^2 \ge 0 $$

若至少存在一个 material point 满足：

$$ g_{n,q} \ne 0 $$

则：

$$ A_n > 0 $$

因此 positivity of $H_\sigma$ 对 temporal weighted least-squares 的良定性同样重要。

这与当前 project 中 search-direction positivity/regularisation safeguards 一致。

---

# 43. 尚未解决的 denominator degeneracy 问题

如果：

$$ \vec{g}_n \approx \vec{0} $$

则：

$$ A_n \approx 0 $$

此时 scalar update 会出现 conditioning 问题。

本阶段不决定具体 safeguard。

后续候选策略包括：

- reject current temporal update；
- reuse previous amplitude；
- reinitialise current enrichment pair；
- trigger mode rejection；
- use tolerance-based pseudoinverse logic；
- interpret current new mode as having negligible residual-control authority。

在理论闭合前，不应简单做 denominator clipping。

---

# 44. 当前 Eq. (72) fixed-point 结构仍保持不变

本阶段 temporal discretisation 并未改变 enrichment alternating architecture：

$$ \lambda^{(k)} \rightarrow \bar{\varepsilon}^{p,(k+1)} \rightarrow \lambda^{(k+1)} \rightarrow \cdots $$

其中：

- fixed $\lambda^{(k)}$：Eq. (70)–(71) 求 spatial mode；
- fixed $\bar{\varepsilon}^{p,(k+1)}$：Eq. (72) temporal update；
- 重复直到 new separated pair fixed-point 收敛。

本阶段仅使 Eq. (72) temporal half-step 在 current project discretisation 下可计算。

---

# 45. 与 Eq. (70)–(71) 的职责边界

Eq. (70)–(71)：

> fixed temporal function → solve spatial FE problem → recover plastic spatial mode。

Eq. (72)：

> fixed spatial mode → solve temporal weighted residual minimisation。

因此 Eq. (70)–(71) 包含 tower structural compatibility/equilibrium，而 Eq. (72) 包含 reduced temporal optimisation。

不能将两者混成同一个 least-squares solve。

---

# 46. current tower v1 temporal strategy 的正式状态

截至本阶段，tower Eq. (72) temporal strategy 建议固定为：

1. 保留原论文 Eq. (72) residual；
2. 保留 $H_\sigma^{-1}$ metric；
3. 保留 fixed-spatial / updated-temporal architecture；
4. spatial integration 使用 tower fiber metric $M$；
5. temporal derivative 使用 current validated 1D backward Euler；
6. time marching 使用 current validated 1D sequential weighted least-squares structure；
7. 不声称 exact DG0 equivalence；
8. 若未来恢复 original DG0 algebra，再进行 term-by-term comparison。

---

# 47. 本阶段对理论忠实性的准确判断

当前路线不是完全重写 Eq. (72)。

被保留的是：

- LATIN search-direction residual structure；
- $H_\sigma^{-1}$ metric；
- new spatial mode fixed, temporal function updated；
- original $x-t$ separated PGD architecture。

被替换或尚未严格复原的是：

- original DG0 temporal algebra。

所以 current tower v1 应定位为：

> paper-faithful continuous LATIN-PGD structure + validated project temporal discretisation。

---

# 48. 为什么这一边界对后续论文写作重要

未来论文中如果直接写：

> the original DG0 scheme was employed

而当前代码实际是 sequential backward-Euler least squares，则方法描述不准确。

更稳妥的写法应明确：

- Eq. (72) residual/minimisation 来源于 Bhattacharyya et al.；
- temporal derivative/discrete marching 继承 validated 1D implementation；
- exact DG0 equivalence is not asserted unless later proven。

这可以避免审稿阶段对 algorithmic fidelity 的质疑。

---

# 49. 本阶段尚未解决的第一个核心问题：$t_0$

当前 sequential formula 适用于：

$$ n \ge 1 $$

因为它依赖：

$$ \lambda_{n-1} $$

所以第一个时间点：

$$ t_0 $$

无法直接使用同一公式。

现有 1D code 对 $t_0$ 采用特殊处理：

- amplitudes initially set to zero；
- first-point temporal rate independently solved by weighted least squares；
- subsequent amplitudes use backward-Euler marching。

因此下一阶段必须明确：

$$ \lambda_0 $$

以及：

$$ \dot{\lambda}_0 $$

在理论和 tower implementation 中分别是什么。

---

# 50. 为什么不能未经推导直接设 $\lambda_0=0$

虽然 current code 的 amplitudes array 初始为 zero，但这并不自动等价于理论上：

$$ \lambda(0) = 0 $$

上一阶段已经确认：

- 原论文没有明确给出 individual enrichment temporal mode 的 endpoint conditions；
- full correction initial admissibility 不等价于每一个 separated mode 都必须单独满足 zero initial amplitude；
- PGD mode scaling 和 orthogonalisation 也会影响 individual temporal amplitude 的解释。

所以必须区分：

- implementation initialisation；
- theoretical boundary condition。

---

# 51. 下一阶段必须回答的问题

下一阶段应严格围绕 $t_0$ 展开，依次确认：

- current 1D code 在 $t_0$ 实际解的是什么 least-squares problem；
- 为什么 first step 求的是 $\dot{\lambda}_0$ 而不是 $\lambda_0$；
- amplitudes[0] 保持 zero 是 algorithmic convention 还是 theoretical condition；
- Eq. (72) continuous residual 在 $t_0$ 应如何理解；
- tower spatial mode fixed 后，first-point rate solve 如何写成 $(e,g,f)$ contraction；
- 是否需要把 $\lambda_0=0$ 固定为 tower v1 convention；
- 如果固定，应明确这是 inherited discretisation convention，而不是 original paper endpoint condition；
- temporal rate reconstruction 如何与 later backward differences 保持一致。

---

# 52. 后续还需解决的数值问题

在 initial time treatment 之后，还需依次处理：

- denominator conditioning safeguard；
- single-new-mode temporal residual reduction check；
- temporal rate storage convention；
- Eq. (72) 更新后何时立即回到 Eq. (70)–(71)；
- fixed-point convergence quantity；
- spatial normalization 后 temporal rescaling；
- candidate mode acceptance/rejection；
- accepted new mode 后 all-mode temporal reoptimisation；
- existing 1D scalar line search 是否保留、删除或标记为 project stabilisation；
- exact original DG0 若未来恢复，如何与 current solver compare。

---

# 53. 当前不进入代码实现的原因

虽然 single-new-mode update 已经得到显式公式，但 temporal solver 尚未从 $t_0$ 到 $T$ 完全闭合。

尤其是 first time point treatment 仍未建立完整理论解释。

因此当前不应立即编写 tower Eq. (72) code。

继续坚持：

> theory first → discrete algorithm closed → code implementation。

---

# 54. 本阶段对 current 1D code 的结论

当前 `latin/pgd_time_update.py` 的核心 temporal structure 与本阶段推导高度一致：

- fixed spatial basis；
- backward-Euler temporal derivative；
- $H_\sigma^{-1}$ weighted residual；
- per-time-step weighted least-squares；
- previous temporal amplitude enters current RHS；
- residual reconstruction and weighted residual norm；
- conditioning diagnostics。

因此它仍然是 tower temporal implementation 的主要继承对象。

但 first-time-point special treatment 和 DG0 wording 需要在下一阶段进一步理论审查。

---

# 55. 本阶段最终结论

本阶段已经完成 Eq. (72) 从 continuous residual 到 current project sequential backward-Euler scalar temporal update 的主要代数闭合：

$$ \vec{r}_n = \vec{g}_n\lambda_n - \vec{b}_n $$

$$ \vec{g}_n = \frac{\vec{p}}{\Delta t_n} - D_{H,n}\vec{s} $$

$$ J_n = \frac{1}{2}\vec{r}_n^T M D_{H,n}^{-1}\vec{r}_n $$

$$ \vec{g}_n^T M D_{H,n}^{-1}\vec{r}_n = 0 $$

$$ \lambda_n = \frac{\vec{g}_n^T M D_{H,n}^{-1}\left( \frac{\vec{p}}{\Delta t_n}\lambda_{n-1}-\vec{\Delta}_n \right)}{\vec{g}_n^T M D_{H,n}^{-1}\vec{g}_n} $$

同时完成了一个关键理论边界修正：

**该 sequential scalar update 不是 whole-time Eq. (72) minimisation 的全局精确离散解，而是 current validated 1D causal discretisation 的 one-step exact minimiser。**

whole-time backward-Euler objective 会产生相邻 temporal amplitudes 的耦合，并形成 time-global algebraic system。

因此 tower v1 的正式方法定位应为：

> original Eq. (72) residual/metric + tower material-point contraction + inherited validated sequential backward-Euler temporal least squares。

而不是 exact original DG0 algebra。

---

# 56. 下一阶段

下一阶段只处理一个问题：

$$ t_0,\lambda_0,\dot{\lambda}_0 $$

即第一个时间点的 initial temporal treatment。

需要从 current `pgd_time_update.py` 的 first-point solve 出发，解释其数学来源、物理/算法意义与 tower $(e,g,f)$ 对应。

在这一问题闭合之前，不进入 tower temporal code implementation。
