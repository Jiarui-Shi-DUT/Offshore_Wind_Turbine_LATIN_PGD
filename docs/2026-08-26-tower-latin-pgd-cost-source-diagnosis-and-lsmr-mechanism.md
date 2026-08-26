# 海上风机塔筒 LATIN-PGD 计算成本来源诊断与 residual-LS / LSMR 机制阶段总结

**日期：** 2026-08-26  
**项目：** Offshore_Wind_Turbine_LATIN_PGD  
**分支：** `feature/offshore-wind-turbine-tower-fatigue`  
**阶段：** FOM-3D-5  
**当前状态：** 1、2、5、10 周期 wall-time scaling 已完成；源码级成本链审计与 residual-LS / LSMR 机制分析已完成；函数级 wall-time profiling 尚未开始  

---

## 1. 本文档的起点与范围

本文档承接上一份正式阶段总结：

`docs/2026-08-23-tower-latin-pgd-efficiency-scaling-1-10cycles.md`

上一份总结对应的 Git checkpoint 为：

`c404d49  docs: summarize LATIN wall-time scaling through 10 cycles`

上一阶段已经确认：

- FOM-3C 的 1、2、5、10 周期匹配精度验证已经关闭；
- FOM-3D 的 1、2、5、10 周期基础 wall-time scaling 已完成；
- 时间点数从 41 增长到 401 时，PGD rank 仅从 11 增长到 21；
- 当前实现的 LATIN-PGD wall time 却从约 23.73 s 增长到约 467.36 s；
- 1 周期和 2 周期时 LATIN 仍快于 FOM；
- 5 周期开始出现效率 crossover；
- 10 周期时 LATIN 已慢于 FOM；
- 当前实现尚未把较好的低秩表达能力转化为持续的 wall-time 优势；
- 在解释成本来源之前，不应直接推进 100-cycle whole-history LATIN 求解。

因此，上一份文档把下一阶段正式任务定义为：

> **FOM-3D-5：LATIN solver wall-time 成本来源诊断。**

本文档总结从该 checkpoint 之后完成的源码审计、算法推导和讨论，重点包括：

1. `solve_tower_latin_pgd()` 的完整 Trial A / Trial B 调用链；
2. outer LATIN iteration、Trial、PGD enrichment 三类计数之间的关系；
3. LATIN indicator `xi` 与 PGD saturation parameter `zeta` 的作用区别；
4. enrichment 的 fixed-point 结构及其接受/拒绝逻辑；
5. `fixed-point iterations` 与 `LSMR iterations` 的严格区别；
6. residual-LS 空间问题为什么是一个 whole-history 加权最小二乘问题；
7. 为什么数学上的大矩阵 `C` 存在，但当前实现没有显式组装它；
8. LSMR 为什么需要迭代，以及它怎样逐步得到空间函数 `p`；
9. 当前实现中两个最值得怀疑的主要耗时区域；
10. 下一阶段如何在不修改核心算法的前提下记录真实 fixed-point 与 LSMR 迭代历史。

本阶段仍然没有修改任何 `latin/` 核心算法。

---

## 2. 当前 benchmark 与数值 checkpoint 保持不变

### 2.1 空间离散

保持：

- 梁柱单元数：10；
- 每单元 Gauss 点数：2；
- 每 Gauss 点环向纤维数：16；
- 径向纤维层数：1；
- 总材料点数：`Nq = 320`；
- 总结构自由度：33；
- 自由自由度：30。

canonical material-point index 仍采用：

`q <-> (element, Gauss point, fiber)`

其排列顺序为：

`element-major -> Gauss-major -> fiber-major`

### 2.2 时间离散

| cycles | Nt |
|---:|---:|
| 1 | 41 |
| 2 | 81 |
| 5 | 201 |
| 10 | 401 |

### 2.3 LATIN 参数

当前正式 benchmark 继续使用：

- `spatial_strategy="residual_ls"`
- `tolerance = 1.0e-5`
- `saturation_enrichment_tolerance = 0.1`
- `saturation_stopping_tolerance = 1.0e-4`
- `fixed_point_tolerance = 1.0e-6`
- benchmark 中 `max_fixed_point_iterations = 200`
- `mode_significance_tolerance = 0.0`
- `acceptance_tolerance = 0.0`

residual-LS 空间求解器内部当前默认使用：

- `atol = 1.0e-12`
- `btol = 1.0e-12`
- `conlim = 1.0e12`
- `max_iterations = 2000`

### 2.4 已完成的 1–10 周期 scaling 数据

| cycles | Nt | rank | LATIN iterations | trials | FOM mean / s | LATIN mean / s | FOM/LATIN |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 41 | 11 | 18 | 29 | 37.142084 | 23.726349 | 1.565379 |
| 2 | 81 | 13 | 23 | 36 | 69.305289 | 59.022472 | 1.174250 |
| 5 | 201 | 17 | 33 | 50 | 161.272881 | 205.363929 | 0.785303 |
| 10 | 401 | 21 | 39 | 60 | 307.891095 | 467.361528 | 0.658699 |

1 周期到 10 周期的增长倍数为：

- `Nt`：约 `9.7805 倍`；
- rank：约 `1.9091 倍`；
- outer LATIN iterations：约 `2.1667 倍`；
- trials：约 `2.0690 倍`；
- FOM wall time：约 `8.2896 倍`；
- LATIN wall time：约 `19.698 倍`。

这一现象仍然是本阶段所有成本分析的起点。

---

## 3. FOM-3D-5 的核心问题重新定义

上一阶段的主要疑问是：

> 为什么时间域只增加约 9.78 倍、PGD rank 只增加约 1.91 倍，而当前 LATIN wall time 却增加约 19.70 倍？

经过本阶段讨论后，需要把这个问题进一步拆成：

1. 一个 outer LATIN iteration 实际调用了哪些 whole-history 操作？
2. Trial A 和 Trial B 分别意味着什么？
3. 一次 successful enrichment 内部到底做了多少次空间与时间求解？
4. residual-LS 空间求解中的一次 LSMR solve 又包含多少次迭代？
5. 每一次 LSMR iteration 是否都需要作用于完整时间历史？
6. fixed-basis temporal update 是否也随着 `Nt` 和 rank 同时变贵？
7. 最终 wall time 的主要来源究竟是 temporal update，还是 enrichment 内部的 fixed-point / LSMR，或者二者共同作用？

当前阶段的原则仍然是：

> **先定位成本，再考虑优化；在没有 profiling 数据以前不修改核心算法。**

---

## 4. `solve_tower_latin_pgd()` 的顶层调用链审计

源码：`latin/tower_latin_pgd_solver.py`

当前每一个 outer LATIN iteration 从一个 persistent baseline 开始：

`(baseline_state, baseline_basis, baseline_indicator)`

随后固定执行：

1. `solve_tower_local_stage(...)`
2. `compute_tower_descent_search_directions(...)`
3. `prepare_frozen_global_data(...)`
4. `update_tower_pgd_time_functions(...)`，形成 Trial A 的 fixed-basis PGD result
5. `build_unrelaxed_candidate(...)`
6. `evaluate_tower_trial(...)`，得到 Trial A
7. 检查 LATIN absolute convergence
8. 若未收敛，则计算 PGD saturation decision
9. 若当前 basis 仍足够有效，则接受 Trial A 并进入下一 outer LATIN iteration
10. 若当前 basis 已不足，则调用 `enrich_tower_pgd_basis_once(...)`
11. 若 enrichment 成功，则形成 Trial B
12. 对 Trial B 进行 candidate construction 与 evaluation
13. 只有合法的 Trial A 或 Trial B 才能 commit 为新的 persistent state/basis/indicator

因此，当前 solver 是严格 transactional 的：

> Trial A 和 Trial B 都是 provisional candidate；只有被接受以后才成为新的 persistent snapshot。

失败的 enrichment 不允许把 provisional Trial A 偷偷提升为 persistent state。

---

## 5. Trial A 与 Trial B 的算法含义

### 5.1 Trial A

Trial A 的核心特征是：

> **不增加新的空间 mode，只在当前已有 spatial basis 上重新求解全部 temporal coordinates。**

如果当前 basis 为：

$$ P=[p_1,p_2,\ldots,p_m] $$

那么 Trial A 保持 `p_1,...,p_m` 不变，只更新：

$$ \lambda_1(t),\lambda_2(t),\ldots,\lambda_m(t) $$

Trial A 回答的问题是：

> 当前这些空间形状是否已经足够，只需要重新调整它们随时间的组合权重？

### 5.2 Trial B

如果 Trial A 表明当前 basis 不够，则尝试增加一个新 pair：

$$ p_{m+1}\lambda_{m+1}(t) $$

成功 enrichment 后，basis 从 `m` 增加到 `m+1`。

但 Trial B 不是只求新 mode 的 `lambda_(m+1)`。

当前实现随后还会对扩大后的整个 basis 再执行一次：

`update_tower_pgd_time_functions(...)`

即重新优化：

$$ \lambda_1(t),\ldots,\lambda_m(t),\lambda_{m+1}(t) $$

因此 Trial B 比简单的“加一个 mode”更重。

---

## 6. outer iterations、modes added 与 trial count 的确定关系

当前所有 matched asymmetric successful runs 都从 empty basis 开始，并且没有 enrichment failure。

因此：

> 每一个 outer LATIN iteration 必定有一个 Trial A；每一次 successful enrichment 额外增加一个 Trial B，同时 rank 增加 1。

于是对于当前 successful runs：

$$ N_{\mathrm{trial}}=N_{\mathrm{outer}}+N_{\mathrm{mode\ added}} $$

实测：

- 1 cycle：`18 + 11 = 29`
- 2 cycles：`23 + 13 = 36`
- 5 cycles：`33 + 17 = 50`
- 10 cycles：`39 + 21 = 60`

特别是 10-cycle：

- 39 次 Local stage；
- 39 次 search-direction construction；
- 39 次 FrozenGlobalData construction；
- 39 次 Trial-A fixed-basis temporal update；
- 39 次 Trial-A candidate/evaluation；
- 21 次 successful enrichment；
- 21 次 Trial-B candidate/evaluation；
- 21 次 enrichment 后 enlarged-basis full temporal re-optimization。

因此至少可以确认：

> **10-cycle 路径中存在 39 + 21 = 60 次 full-basis temporal update。**

这里的“60 次”是源码调用链能够确认的事实，不是复杂度估计。

---

## 7. `xi` 与 `zeta`：两个不同层级的问题

本阶段讨论中一个重要澄清是：

> LATIN convergence indicator `xi` 与 PGD saturation parameter `zeta` 不能混为一谈。

### 7.1 `xi` 回答“离整体 LATIN 收敛还有多远”

当前 benchmark 使用 `tolerance = 1.0e-5`。

如果 Trial A 已经满足：

$$ \xi_A \le 10^{-5} $$

则 absolute LATIN convergence 优先级最高，Trial A 可以直接 commit 并终止求解。

### 7.2 `zeta` 回答“现有 basis 是否仍然有效”

当前 saturation parameter 为：

$$ \zeta=\frac{\xi_{\mathrm{previous}}-\xi_{\mathrm{current}}}{\xi_{\mathrm{previous}}+\xi_{\mathrm{current}}} $$

当前阈值：

- `zeta_enrich = 0.1`
- `zeta_stop = 1.0e-4`

概念上：

- `xi` 看整体收敛距离；
- `zeta` 看当前 fixed spatial basis 对降低 `xi` 是否仍然有效。

### 7.3 为什么 `xi` 仍较大时 Trial A 仍可能被接受

例如仅作说明：

- `xi_previous = 0.010`
- `xi_current = 0.005`

则：

$$ \zeta=\frac{0.010-0.005}{0.010+0.005}=0.333 $$

虽然 `xi_current` 仍远大于 `1e-5`，但 `zeta > 0.1` 表明已有 basis 仍然能够有效降低 LATIN indicator。

此时合理动作不是立刻 enrichment，而是接受 Trial A 并推进到下一 LATIN iteration。

因此：

> **“尚未整体收敛”和“当前 basis 已不足”是两个不同问题。**

---

## 8. enrichment 何时发生，以及成功后做什么

源码：`latin/tower_pgd_enrichment.py`

当 solver 判断当前 basis 不足并且 reduced residual 仍未解决时，才进入：

`enrich_tower_pgd_basis_once(...)`

当前 enrichment 的完整逻辑为：

1. 从 Trial-A mechanical residual 形成 shifted defect；
2. 从 residual energy 最大的时间位置构造 initial spatial seed；
3. 对 residual-LS 路径做 in-loop M-orthogonalisation；
4. 根据当前 `p` 求一次 temporal function；
5. 进入 fixed-point loop；
6. 根据当前 `lambda(t)` 求新的空间函数 `p`；
7. residual-LS 路径中，这一步调用 matrix-free LSMR；
8. 新 `p` 再做 in-loop M-orthogonalisation；
9. 根据新的 `p` 再求 `lambda(t)`；
10. 根据 complete-pair change `chi` 判断 fixed-point convergence；
11. fixed-point 收敛后做 M-weighted modified Gram-Schmidt；
12. 对已有 temporal coordinates 做 exact coordinate transformation；
13. 检查 field invariance、orthogonality、spatial novelty、temporal significance；
14. 将新 mode 加入 basis；
15. 对 enlarged basis 执行一次 full `update_tower_pgd_time_functions(...)`；
16. 比较 enrichment 前后的 full mechanical residual；
17. 只有 residual benefit 通过阈值才接受新 mode。

这说明：

> **一次 successful enrichment 远不只是“一次求 p + 一次求 lambda”。**

---

## 9. fixed-point iteration 到底是什么

对一个待新增的 PGD mode：

$$ \Delta\dot{\varepsilon}^{p}_{\mathrm{new}}(t)=\dot{\lambda}(t)p $$

空间函数 `p` 与时间函数 `lambda(t)` 相互依赖。

因此当前 enrichment 采用交替 fixed-point：

`seed p0 -> lambda0 -> p1 -> lambda1 -> p2 -> lambda2 -> ...`

一个 fixed-point iteration 可以通俗理解为：

> **根据当前时间函数重新求一次空间函数，再根据新的空间函数重新求一次时间函数。**

当前 convergence measure 为 complete-pair change `chi`，并使用：

`fixed_point_tolerance = 1.0e-6`

当：

$$ \chi \le 10^{-6} $$

时认为该新 pair 的空间和时间部分已经稳定。

源码中 `TowerEnrichmentResult` 已经保存：

- `fixed_point_history`
- `fixed_point_iterations`
- `fixed_point_converged`

因此单个 enrichment 的 fixed-point 次数实际上已经可观测。

---

## 10. enrichment acceptance 与 failure rollback

fixed-point 收敛并不自动意味着新 mode 一定被接受。

当前代码还要检查 enlarged-basis full residual benefit：

$$ \mathrm{benefit}=1-\frac{R_{\mathrm{after}}}{R_{\mathrm{before}}} $$

当前 benchmark：`acceptance_tolerance = 0.0`。

因此要求：

$$ \mathrm{benefit}>0 $$

也就是新 mode 至少必须真正降低 full residual。

如果：

$$ \mathrm{benefit}\le0 $$

则 enrichment 被拒绝。

### 10.1 enrichment failure 后 solver 做什么

当前 transactional solver 不会：

- 自动修改参数；
- 自动重试另一个 enrichment；
- 自动保留失败的 mode；
- 自动把 provisional Trial A 当成 persistent state。

相反，它返回最后可靠的 persistent baseline，并以 `ENRICHMENT_FAILED` 结束当前 solver call。

### 10.2 当前 1/2/5/10 周期是否发生过 enrichment failure

没有。

当前四组正式 matched asymmetric runs 均：

- `converged=True`；
- 最终 rank 分别为 11、13、17、21；
- 没有观测到 `ENRICHMENT_FAILED`。

因此：

> **当前 wall-time 增长不能解释为 enrichment 失败以后不断 retry。**


---

## 11. fixed-point 次数与 LSMR 迭代次数的严格区别

这是本阶段讨论中最重要的概念澄清之一。

### 11.1 fixed-point iteration 是外层

一次 fixed-point iteration 表示：

`当前 lambda -> 求 p -> 新 p -> 求新 lambda`

所以：

> `fixed_point_iterations` = 为了使一个新 PGD pair 稳定，空间与时间之间来回修正了多少轮。

### 11.2 LSMR iteration 是 fixed-point 内部的更深一层

在当前：

`spatial_strategy="residual_ls"`

下，每一次 fixed-point 中的“求 p”不是直接一次求完，而是调用：

`solve_tower_residual_ls_spatial(...)`

该函数内部又调用：

`scipy.sparse.linalg.lsmr(...)`

所以层级是：

`one enrichment -> many fixed-point iterations -> one LSMR solve per fixed-point -> many LSMR iterations per LSMR solve`

### 11.3 定量示意

下面仅为示意，不是当前实测：

如果一个 mode 需要 6 次 fixed-point，而每次 LSMR 平均 10 iterations，则约有：

`6 × 10 = 60`

次 LSMR inner iterations。

如果 fixed-point 仍是 6 次，但每次 LSMR 平均 50 iterations，则约有：

`6 × 50 = 300`

次 LSMR inner iterations。

因此相同的 fixed-point count 可能对应完全不同的实际空间求解成本。

---

## 12. 为什么 residual-LS 空间问题牵涉整个时间历史

源码：`latin/tower_residual_ls_spatial.py`

在一次 enrichment 的 fixed-point 中，当时间函数已经固定以后，需要求一个空间函数：

$$ p\in\mathbb{R}^{N_q} $$

当前：

$$ N_q=320 $$

这个 `p` 不是每个时间点分别求一个不同的空间形状。

PGD separation 要求同一个空间 mode 在整个时间域内共用：

$$ p\lambda(t) $$

因此同一个 `p` 必须同时对 `t_1,t_2,\ldots,t_{N_t}` 全部时间点尽可能合适。

10-cycle 时：

- `Nt = 401`
- `Nq = 320`

因此 space-time residual 一共有：

$$ 401\times320=128320 $$

个分量。

但未知的空间向量仍然只有：

$$ p=(p_1,\ldots,p_{320})^T $$

即 320 个未知量。

这正是一个高瘦的 overdetermined least-squares 问题。

---

## 13. residual 对 `p` 为什么是线性的

固定当前时间函数以后，在第 `n` 个时间点可以把 residual 写为：

$$ r_n(p)=\dot{\lambda}_n p-\lambda_n H_{\sigma,n}A_\sigma p+d_n $$

这里：

- `p` 是当前未知空间函数；
- `lambda_n` 已知；
- `lambda_dot_n` 已知；
- `H_sigma,n` 已知；
- defect `d_n` 已知；
- `A_sigma` 是 reference equilibrium stress action，对 `p` 是线性的。

因此可以定义：

$$ B_n=\dot{\lambda}_n I-\lambda_nH_{\sigma,n}A_\sigma $$

则：

$$ r_n(p)=B_np+d_n $$

所以在 fixed temporal function 下，空间问题对 `p` 是线性的。

需要特别区分：

> **“空间问题是线性的”并不意味着“数值上必须使用直接法”。**

当前只是选择用 iterative least-squares solver LSMR 求这个线性最小二乘问题。

---

## 14. 从每个时间点的 residual 到大系统 `Cp ≈ b`

为了理解矩阵结构，可以把每个时间点的关系堆起来。

假设教学例子只有两个空间未知量：

$$ p=(p_1,p_2)^T $$

某一个时间点可能产生：

$$ 2p_1+p_2\approx5 $$

$$ p_1+3p_2\approx7 $$

第二个时间点由于 `lambda`、`lambda_dot`、`H_sigma`、defect 改变，又可能产生：

$$ 2p_1+2p_2\approx8 $$

$$ 2p_1+4p_2\approx6 $$

虽然方程系数随时间变，但空间未知量仍是同一个 `p_1,p_2`。

因此可堆成：

$$ Cp\approx b $$

真实塔筒中：

- `C` 数学尺寸为 `(Nt × Nq) × Nq`；
- 10-cycle 为 `128320 × 320`；
- `p` 为 `320 × 1`；
- `b` 为 `128320 × 1`。

---

## 15. 为什么不能要求全部 residual 都严格为零

对于 overdetermined system：

$$ Cp\approx b $$

通常不存在一个 320 维 `p` 能让全部 128320 条条件都精确成立。

因此目标不是：

$$ Cp-b=0 $$

而是寻找一个折中空间模式，使所有 residual 的综合度量最小。

普通 least squares 写为：

$$ J(p)=\frac12\|Cp-b\|_2^2 $$

关键点是：

> 虽然 `J` 最终是一个标量，但它不是固定数，而是 320 个未知量的函数。

即：

$$ J=J(p_1,p_2,\ldots,p_{320}) $$

最优解满足：

$$ \nabla_pJ=0 $$

对于普通 least squares：

$$ C^T(Cp-b)=0 $$

即正规方程：

$$ C^TCp=C^Tb $$

需要注意：

> 最优解一般不要求 `Cp-b=0`，而要求 normal residual `C^T(Cp-b)=0`。

---

## 16. 当前 residual-LS 实际是加权最小二乘

当前代码不是简单最小化所有 residual 的普通平方和。

它使用：

- 时间积分权重；
- material-point metric 权重；
- `H_sigma^{-1}` 权重。

源码中：

`sqrt_weight = sqrt(time_weight * material_weight / H_sigma)`

因此可以把第 `n` 个时间点的权重矩阵概念上写成：

$$ W_n=w_nMH_{\sigma,n}^{-1} $$

于是当前空间 objective 可以理解为：

$$ J_h(p)=\frac12\sum_n r_n(p)^TW_nr_n(p) $$

代入：

$$ r_n(p)=B_np+d_n $$

得到：

$$ J_h(p)=\frac12\sum_n(B_np+d_n)^TW_n(B_np+d_n) $$

其一阶最优条件为：

$$ \sum_nB_n^TW_n(B_np+d_n)=0 $$

也就是：

$$ \left(\sum_nB_n^TW_nB_n\right)p=-\sum_nB_n^TW_nd_n $$

这说明从数学上完全可以得到一个 `Nq × Nq` 的正规系统。

对于当前塔筒：

$$ N_q\times N_q=320\times320 $$

---

## 17. `C` 到底“有没有”：数学存在与显式存储必须区分

本阶段讨论中另一个重要澄清是：

> **数学上的 `C` 当然存在。当前程序只是没有把它显式 materialize 成一个 `128320 × 320` 的 NumPy dense array。**

### 17.1 显式矩阵路线

如果真的组装：

`C.shape = (128320, 320)`

则元素数量约为：

$$ 128320\times320=41062400 $$

使用 float64 时，仅该 dense matrix 的原始数据约为：

`328 MB`

还没有包括其它工作数组和求解 workspace。

### 17.2 当前 matrix-free 路线

当前程序没有存储整个 `C`，而是只告诉 SciPy：

- 给一个 `p`，怎样算 `Cp`；
- 给一个 `y`，怎样算 `C^Ty`。

因此当前代码建立的是：

`scipy.sparse.linalg.LinearOperator`

其 shape 仍然是：

`(Nt*Nq, Nq)`

所以正确表述是：

> **`C` 数学上存在，程序也知道它如何作用于向量，但不显式保存它的全部矩阵元素。**

这就是当前 residual-LS 所谓的 matrix-free。

---

## 18. 当前 `Cp` 的具体计算

源码中的 `matvec(vector)` 对应：

$$ p\longmapsto Cp $$

首先通过 equilibrium operator 得到：

$$ A_\sigma p $$

然后对整个时间历史构造：

$$ \dot{\lambda}(t)p-\lambda(t)H_\sigma(t)A_\sigma p $$

最后乘以：

$$ \sqrt{w_tM H_\sigma^{-1}} $$

并 reshape 为长度：

$$ N_tN_q $$

的向量。

10-cycle 时该长度为 `128320`。

---

## 19. 当前 `C^Ty` 的具体计算

源码中的 `rmatvec(vector)` 对应：

$$ y\longmapsto C^Ty $$

它先把输入恢复为 `(Nt, Nq)`，然后：

1. 乘入 space-time square-root weights；
2. 对 `rate × z` 沿时间求和，形成 direct part；
3. 对 `lambda × H_sigma × z` 沿时间求和，形成 adjoint source；
4. 利用 equilibrium stress action 的 M-adjoint 关系计算 stress adjoint part；
5. 返回一个长度为 `Nq=320` 的向量。

因此当前 matrix-free implementation 已经完整提供了 LSMR 所需的 forward action 与 adjoint action。

---

## 20. `A_sigma` 也没有显式形成 dense material-point matrix

当前 equilibrium operator 同样采用 operator form。

概念上 reference structural projection 为：

$$ E=H(H^TMC_0H)^{-1}H^TMC_0 $$

对于 source strain `p`，得到：

$$ \varepsilon=Ep $$

$$ \sigma=C_0(\varepsilon-p) $$

因此概念上：

$$ A_\sigma=C_0(E-I) $$

但当前代码不会显式构造完整 `A_sigma`。

`apply_spatial(p)` 实际执行：

1. 计算 `H^T M C0 p`；
2. 使用已经 Cholesky factorized 的 reduced reference stiffness 解自由度空间方程；
3. 恢复 compatible strain；
4. 恢复 stress correction。

所以当前整个 residual-LS 设计思想是一致的：

> **优先定义 operator action，而不是显式形成 material-point dense operators。**

---

## 21. 为什么 LSMR 还要“迭代”求一个线性问题

需要明确区分两个概念。

### 21.1 线性问题

固定 `lambda(t)` 后：

$$ Cp\approx b $$

对 `p` 是线性的。

### 21.2 迭代求解器

“线性”只是问题的数学性质，并不要求必须使用 LU、QR、SVD 等直接法。

当前实现选择的是：

`scipy.sparse.linalg.lsmr`

LSMR 是针对大型 least-squares problem 的 iterative Krylov solver。

它的主要优势是：

> **只需要 `Cp` 和 `C^Ty`，不要求显式存储 `C`，也不要求显式形成 `C^TC`。**

因此 LSMR 的迭代并不是因为 fixed-`lambda` 空间问题是非线性的。

它迭代的原因是：

> **当前选择了 matrix-free iterative least-squares 路线。**

---

## 22. LSMR 怎样逐步得到 `p`

本阶段为了理解 LSMR，使用过只有两个未知量的教学例子。

设：

$$ Cp\approx b $$

LSMR 的直观过程不是逐个修改 `p_1,p_2,...`，而是逐步建立未知数空间中的搜索方向：

$$ v_1,v_2,v_3,\ldots $$

第 1 次迭代后，解被限制在：

$$ p^{(1)}=a_1v_1 $$

第 2 次迭代后：

$$ p^{(2)}=a_1v_1+a_2v_2 $$

第 `k` 次迭代后：

$$ p^{(k)}\in\mathrm{span}\{v_1,\ldots,v_k\} $$

这些方向由 Golub-Kahan bidiagonalization 通过交替的 `C v_k` 与 `C^T u_k` 逐步生成。

因此可以通俗理解为：

1. 从 `b` 的信息出发；
2. 通过 `C^T` 拉回未知数 `p` 的空间，找到一个最值得尝试的方向；
3. 通过 `C` 检查该方向在 residual space 中解释了什么；
4. 从剩余信息中产生新的方向；
5. 在不断扩大的 Krylov subspace 中更新最优 `p`；
6. 直到 residual 或 normal residual 达到停止条件。

当前 SciPy LSMR 返回值中，代码显式保存：

- `iterations`
- `residual_norm`
- `normal_residual_norm`
- `operator_norm_estimate`
- `condition_estimate`
- `istop`

并将 `istop in (1, 2, 4, 5)` 视为收敛。

---

## 23. fixed-point 与 LSMR 为什么是两层完全不同的迭代

### 23.1 fixed-point iteration 的原因

因为 `p` 与 `lambda(t)` 相互依赖，因此需要：

`lambda -> p -> lambda -> p -> ...`

### 23.2 LSMR iteration 的原因

在一次 fixed-point 中，`lambda(t)` 已经固定。

此时空间问题虽然是线性的，但采用 matrix-free LSMR 求解，所以内部还要：

`C/C^T actions -> Krylov search directions -> update p`

### 23.3 当前真实嵌套关系

`LATIN outer iteration`

`-> if enrichment`

`-> fixed-point iterations`

`-> each fixed-point calls one residual-LS spatial solve`

`-> each spatial solve calls one LSMR`

`-> each LSMR contains many LSMR iterations`

因此成本分析不能只看最终 PGD rank。

---

## 24. 为什么一次 LSMR iteration 会受到 `Nt` 的直接影响

当前 matrix-free operator 的 forward/adjoint action 都涉及完整 time-by-material history。

1-cycle：

$$ 41\times320=13120 $$

10-cycle：

$$ 401\times320=128320 $$

因此仅从 space-time field size 看：

$$ \frac{128320}{13120}\approx9.78 $$

也就是说，在其它条件完全不变时：

> **10-cycle 的一次 whole-history operator action 需要面对约为 1-cycle 9.78 倍的 time-material entries。**

LSMR 在每轮 bidiagonalization 中通常需要 forward/adjoint operator action，因此 `Nt` 增长会直接放大每一个 LSMR solve 的单位迭代成本。

更关键的是目前还不知道：

> 10-cycle 的 LSMR iteration count 是否也比 1-cycle 更高。

这正是后续 profiling 必须测量的量。

---

## 25. 为什么必须测每个 mode 的 fixed-point count

最终 rank 只告诉我们成功加入了多少个 mode。

例如 10-cycle：`rank = 21`。

但它不告诉我们每一个 mode 是容易找到还是很难找到。

示意：

- 若平均 3 个 fixed-point iterations / mode，则 21 个 mode 约有 63 次 spatial solve；
- 若平均 6 个，则约 126 次；
- 若平均 10 个，则约 210 次。

这些只是示意，不是当前实测。

因此：

> **rank growth 低，并不自动意味着 enrichment search cost 低。**

需要记录每个 enrichment 的 `fixed_point_iterations`，才能知道后期 mode 是否越来越难搜索。


---

## 26. 为什么还必须测每一次 spatial solve 的 LSMR iterations

仅知道 fixed-point count 仍然不够。

因为每一个 fixed-point 的 spatial solve 本身是 iterative LSMR。

例如同样都是 6 次 fixed-point：

- 若每次 LSMR 平均 10 iterations，则约 60 个 inner iterations；
- 若每次 LSMR 平均 50 iterations，则约 300 个 inner iterations。

因此真实空间搜索成本同时受到：

$$ J_{\mathrm{FP}} $$

与：

$$ J_{\mathrm{LSMR}} $$

控制。

可以把它记成：

> fixed-point count = “空间与时间来回修正多少轮”；

> LSMR count = “每一轮空间求解内部又迭代多少次”。

---

## 27. fixed-basis temporal update 的源码成本结构

源码：`latin/tower_pgd_time_update.py`

对于已有 basis：

$$ P=[p_1,\ldots,p_m] $$

$$ S=[s_1,\ldots,s_m] $$

其 mechanical residual 为：

$$ r(t)=P\dot{\lambda}(t)-H_\sigma(t)S\lambda(t)-f(t) $$

目标是在 `H_sigma^{-1}` material metric 下最小化 whole-history residual。

当前实现：

1. 构造 `spatial_plastic`，shape 为 `(Nq,m)`；
2. 构造 `spatial_stress`，shape 为 `(Nq,m)`；
3. 分配 `amplitudes` 与 `rates`，shape 为 `(Nt,m)`；
4. 对初始时间点做一次 weighted least squares；
5. 显式执行 `for step in range(1, basis.n_time)`；
6. 每一个时间步建立一个 `(Nq,m)` 的 `reduced_matrix`；
7. 每一个时间步调用 `np.linalg.lstsq(...)`；
8. 全部时间点结束后重建 full-history plastic correction、rate correction、stress correction；
9. 调用 equilibrium operator 的 `apply_history(...)`；
10. 构造 whole-history mechanical residual；
11. 计算加权 residual norm；
12. 返回一个 immutable `FixedBasisPGDResult`。

因此这一部分的成本不仅随 `Nt` 增长，也随当前 rank `m` 增长。

---

## 28. fixed-basis temporal update 的复杂度只能作为“推断”，不能当成 profiling 结果

对于 `Nq >> m` 的 dense least-squares，小型矩阵 `(Nq,m)` 的典型计算量可以粗略理解为与：

$$ O(N_qm^2) $$

同量级。

在 `Nt` 个时间点上，则可形成一个非常粗略的启发式：

$$ O(N_tN_qm^2) $$

这一表达式是源码结构推断，不是 wall-time 实测，也不是对当前全部 NumPy/LAPACK 实现的严格 Big-O benchmark。

为了直观说明 `Nt` 与 `m` 同时增长的效应，可以看一个非常粗略的代理量：

1-cycle 末端：

$$ 41\times11^2=4961 $$

10-cycle 末端：

$$ 401\times21^2=176841 $$

比值约 `35.6`。

但必须强调：

- 实际每次 temporal update 的 rank 在迭代过程中并不是一直等于最终 rank；
- 该数字不是 operation count；
- 更不是测得的 wall-time ratio；
- 它只说明 `Nt × m^2` 的组合增长值得警惕。

---

## 29. 一个有用但必须谨慎解释的粗粒度 wall-time 代理

本阶段还讨论了一个非常粗的 whole-history reuse 指标：

$$ N_{\mathrm{trial}}\times N_t $$

1-cycle：

`29 × 41 = 1189`

10-cycle：

`60 × 401 = 24060`

比值约 `20.2`。

而实测 LATIN wall-time 增长约 `19.7`。

两者数量级很接近。

但这个结果只能作为直觉提示：

> whole-history operations 被更多 Trial 反复调用，可能是 wall-time 增长的重要原因。

不能把 `trial × Nt` 解释成真实 FLOP count，也不能据此认定 Trial 是唯一瓶颈。

---

## 30. 当前最值得怀疑的两个主要耗时区域

经过源码调用链和算法机制分析，目前有两个优先级最高的 suspect。

### 30.1 嫌疑 A：反复的 full-basis temporal update

它的特点是：

- 每个 outer iteration 至少发生一次 Trial-A temporal update；
- 每次 successful enrichment 后还会发生一次 enlarged-basis temporal re-optimization；
- 10-cycle 因此至少发生 60 次 full temporal update；
- 每次 update 遍历整个 `Nt`；
- 每个时间点调用 dense `np.linalg.lstsq`；
- 成本还随 rank `m` 增长。

这是“确定大量发生”的成本。

### 30.2 嫌疑 B：enrichment 内部 fixed-point + residual-LS + LSMR

它的特点是：

- 10-cycle 有 21 次 successful enrichment；
- 每个 enrichment 需要若干 fixed-point iterations；
- 每个 fixed-point residual-LS 路径调用一次 LSMR；
- 每个 LSMR 又包含若干 inner iterations；
- 每个 LSMR iteration 的 operator action 与 whole-history `Nt × Nq` 有关。

这是“结构上可能非常昂贵，但目前缺少真实内部计数”的成本。

### 30.3 当前不能下的结论

目前不能写：

> `LSMR is the bottleneck.`

也不能写：

> `temporal update is the bottleneck.`

当前正确结论只能是：

> **这两个区域是目前源码层面最主要的优先诊断对象；谁真正占主导，必须用 profiling 数据回答。**

---

## 31. 为什么“低 rank”与“快”不能直接画等号

10-cycle：

- `Nt = 401`
- `rank = 21`

rank 只有 21 的确说明：

> whole-history response 可以被相对少量的 space-time separated modes 表达。

但是 solver 成本还包括：

- outer LATIN iterations；
- Trial A；
- Trial B；
- fixed-basis temporal re-optimization；
- mode search；
- fixed-point；
- LSMR；
- basis orthogonalisation；
- full residual checks；
- state construction/evaluation；
- Local/Global history operations。

因此：

> **低秩 representation 是表示效率；wall-time speedup 是算法与实现效率。二者相关，但不等价。**

---

## 32. 能否不使用 LSMR，直接形成 `320 × 320` 系统

从数学上：可以。

由：

$$ J_h(p)=\frac12\sum_n(B_np+d_n)^TW_n(B_np+d_n) $$

可以直接得到：

$$ K_{\mathrm{LS}}=\sum_nB_n^TW_nB_n $$

$$ f_{\mathrm{LS}}=-\sum_nB_n^TW_nd_n $$

然后求：

$$ K_{\mathrm{LS}}p=f_{\mathrm{LS}} $$

其中 `K_LS` 为 `320 × 320`。

一个 `320 × 320` float64 matrix 仅约 `0.82 MB`。

因此：

> **当前使用 LSMR 是实现选择，不是因为数学上无法构造正规系统。**

但目前不能直接因此判定“改成直接法一定更快”。

需要考虑：

1. 构造 `K_LS` 本身需要多少 operator work；
2. 是否能利用 equilibrium operator 的结构高效形成矩阵；
3. LSMR 实际需要多少 iterations；
4. normal equations 会带来条件数平方问题；
5. 直接形成正规方程可能降低数值稳定性；
6. 当前 `Nq=320` 较小，但未来更细 tower discretization 的 `Nq` 会增长。

因此这属于后续 potential optimization，不属于当前应立即修改的算法。

---

## 33. 正规方程的数值稳定性提醒

若显式把：

$$ Cp\approx b $$

转为：

$$ C^TCp=C^Tb $$

则经典关系为：

$$ \kappa(C^TC)=\kappa(C)^2 $$

因此即使 `C^TC` 尺寸只有 `320 × 320`，其 condition number 可能比原问题显著恶化。

这也是 LSMR/LSQR 类方法不显式形成 normal equations 的重要数值动机之一。

所以未来若比较 direct-normal-system 与 LSMR，不能只比较 wall time，还要比较：

- residual；
- normal residual；
- condition estimate；
- resulting mode；
- outer LATIN convergence path；
- 最终 FOM accuracy。

---

## 34. 当前 source-level cost model

在不声称精确 FLOP count 的前提下，可以把当前 LATIN 总成本概念化为：

$$ T_{\mathrm{LATIN}}\approx\sum_{i=1}^{N_{\mathrm{outer}}}\left(T_{\mathrm{Local}}+T_{\mathrm{search}}+T_{\mathrm{frozen}}+T_{\mathrm{time},m_i}+T_{\mathrm{TrialA}}\right)+\sum_{e=1}^{N_{\mathrm{enrich}}}\left(T_{\mathrm{FP},e}+T_{\mathrm{transform},e}+T_{\mathrm{time},m_e+1}+T_{\mathrm{TrialB},e}\right) $$

其中：

$$ T_{\mathrm{FP},e}\approx\sum_{j=1}^{J_{\mathrm{FP},e}}\left(T_{\mathrm{LSMR},e,j}+T_{\mathrm{single\ mode\ time},e,j}+T_{\mathrm{checks},e,j}\right) $$

而一个 LSMR solve 的成本又受到 `J_LSMR` 以及每次 operator action 的 whole-history size `Nt × Nq` 共同控制。

这一 cost model 的价值在于明确“需要测什么”，而不是提前给出谁最耗时的结论。

---

## 35. 本阶段哪些结论是实测，哪些只是源码推断

### 35.1 已实测

以下是正式运行结果：

- 1/2/5/10-cycle wall time；
- 1/2/5/10-cycle convergence；
- outer iterations；
- trial counts；
- final rank；
- final LATIN indicator；
- FOM/LATIN ratio；
- 1/2/5/10-cycle successful runs 没有 observed enrichment failure。

### 35.2 源码可确定

以下是调用链和实现事实：

- 每个 outer iteration 必有 Trial A；
- successful enrichment 增加 Trial B；
- 10-cycle 对应 39 Trial A + 21 Trial B；
- successful enrichment 最后会进行 enlarged-basis full temporal re-optimization；
- residual-LS fixed-point 每轮调用一次 LSMR spatial solver；
- `TowerEnrichmentResult` 保存单次 enrichment 的 fixed-point iterations；
- `TowerResidualLSSpatialResult` 保存单次 LSMR solve 的 `iterations`；
- 当前 top-level solver 只保存 `last_enrichment_result`，并没有保存全部 enrichment result history；
- `_raw_fixed_point` 当前也没有把每一轮 residual-LS result 的 LSMR `iterations` 向上传递成完整历史。

### 35.3 目前只是推断/示意

以下不能当成实测：

- `O(Nt*Nq*m^2)` 是 temporal update 的启发式复杂度描述；
- `trial × Nt` 是粗粒度 proxy；
- “平均每个 mode 6 次 fixed-point”只是教学示意；
- “平均每次 LSMR 10/20/50 次”只是教学示意；
- LSMR 与 temporal update 谁占主要 wall time 尚未测定。

这一分类后续必须保持。

---

## 36. 下一步为什么先测 fixed-point 与 LSMR count

当前 profiling 的第一目标不应该立即改算法，而应该回答：

> **enrichment 内部到底做了多少工作？**

需要的最小观测量为：

1. 每一次 enrichment 的 `fixed_point_iterations`；
2. 每一次 fixed-point spatial solve 的 LSMR `iterations`；
3. 每一次 LSMR solve 的 `istop`；
4. residual norm / normal residual norm；
5. condition estimate；
6. 最好同时记录各层 wall time，后续才能与 temporal update 比较。

有了这些数据以后，才能区分：

- 是 10-cycle 需要更多 fixed-point rounds；
- 还是每个 LSMR solve 本身需要更多 inner iterations；
- 还是二者都没有明显增加，主要成本其实来自 repeated full-basis temporal update；
- 或者多个成本同时增长。

---

## 37. 不修改核心算法的外部 recorder 思路

当前最安全的第一步是：

> **使用外部 diagnostic wrapper 记录结果，而不修改 `latin/` 生产代码。**

### 37.1 记录 enrichment fixed-point history

`enrich_tower_pgd_basis_once(...)` 本身已经返回 `TowerEnrichmentResult.fixed_point_iterations`。

因此外部 wrapper 可以：

1. 调用原始 `enrich_tower_pgd_basis_once(...)`；
2. 完整保留原函数行为；
3. 从返回值读取 `fixed_point_iterations`；
4. append 到外部 list；
5. 原样返回结果。

这样如果 10-cycle 一共有 21 次 successful enrichment，就可以得到 21 个真实 fixed-point counts。

### 37.2 记录每次 residual-LS 的 LSMR iteration

`solve_tower_residual_ls_spatial(...)` 已经返回 `TowerResidualLSSpatialResult.iterations`。

因此也可以用外部 wrapper：

1. 调用原始 residual-LS solver；
2. 记录 `iterations`、`istop`、residual 等；
3. 原样返回结果。

当前 `_raw_fixed_point` 对 residual-LS solver 使用 lazy import，因此原则上可以从外部替换对应 module attribute 来记录数据，而无需修改生产源码。

### 37.3 这一方案的性质

该方案只做 observation，不改变：

- spatial strategy；
- fixed-point tolerance；
- LSMR tolerance；
- PGD basis；
- LATIN convergence logic；
- Trial A/B transaction；
- material model；
- equilibrium operator。

因此适合作为下一阶段的第一项诊断。

---

## 38. profiling 不能直接从 10-cycle 开始，先用 1-cycle 做数值不扰动验证

下一步建议首先在 1-cycle 上运行 recorder。

原因是：

- wall time 短；
- 已有非常稳定的数值 checkpoint；
- 如果 wrapper 不小心改变调用路径，可以快速发现；
- 通过后再扩大到 2/5/10-cycle。

1-cycle 必须保持：

- `converged = True`
- `iterations = 18`
- `trials = 29`
- `rank = 11`
- `final xi = 7.918424536257e-06`

如果这些 invariants 发生变化，则 recorder 不能被视为纯 measurement instrumentation。

通过 1-cycle 后，再推进：

`2 -> 5 -> 10 cycles`

10-cycle 最终必须保持：

- `converged = True`
- `iterations = 39`
- `trials = 60`
- `rank = 21`
- `final xi = 8.941607234831e-06`

---

## 39. 后续真正需要的 profiling 表

### 39.1 enrichment-level table

建议字段：

| enrichment index | outer iteration | basis rank before | basis rank after | fixed-point iterations | accepted | residual benefit |
|---:|---:|---:|---:|---:|---|---:|

它回答：

> 后期 mode 是否越来越难通过 fixed-point 找到？

### 39.2 LSMR-level table

建议字段：

| enrichment index | fixed-point index | LSMR iterations | istop | residual norm | normal residual norm | condition estimate | wall time |
|---:|---:|---:|---:|---:|---:|---:|---:|

它回答：

> 同一个 enrichment 内，空间 least-squares solve 是否越来越难？不同 cycle count 下 LSMR 迭代数是否增长？

### 39.3 后续再增加 top-function timing

在上述内部计数确认以后，可以进一步统计：

- Local stage time；
- search-direction time；
- frozen-data time；
- Trial-A full temporal update time；
- enrichment total time；
- residual-LS spatial time；
- single-mode temporal solve time；
- enlarged-basis full temporal re-optimization time；
- Trial-B build/evaluation time。

这样才有条件回答：

> **当前约 467 s 的 10-cycle LATIN wall time，究竟有多少秒花在 temporal update，多少秒花在 LSMR，多少秒花在其它环节。**

---

## 40. 现阶段对“最耗时部分”的正式表述

截至本文档，不应写：

> `LSMR is the bottleneck.`

也不应写：

> `temporal update is the bottleneck.`

当前能够严谨写出的结论是：

> **源码审计表明，当前 whole-history tower LATIN-PGD 的两个最高优先级成本候选分别是 repeated full-basis temporal update，以及 enrichment 中嵌套的 fixed-point / matrix-free residual-LS / LSMR。二者都具有随时间历史增长被重复放大的结构。现阶段尚缺少函数级 wall-time 与完整 inner-iteration history，因此不能判断谁占主导。**

这就是下一步 profiling 的直接研究问题。

---

## 41. 与原论文理论的边界仍需保持

### 41.1 原论文直接支持的部分

原论文 Eq. (58)–(59) 明确要求在已有 spatial PGD basis 上更新 temporal functions，并在 `H_sigma^{-1}` metric 下最小化 mechanical residual。

原论文 Eq. (60) 定义 PGD saturation parameter。

原论文 Eq. (61)–(72) 给出 enrichment search direction / fixed-point 的理论框架。

### 41.2 当前 tower implementation 的工程选择

以下不能直接说成原论文原样算法：

- `spatial_strategy="residual_ls"`；
- matrix-free `LinearOperator + scipy LSMR`；
- transactional Trial A / Trial B；
- 当前 acceptance checks；
- 当前 runtime diagnostics。

尤其：

> **residual-LS 是当前 tower implementation 的工程扩展，不是 Eq. (65)–(71) 的直接代数恒等改写。**

### 41.3 仍未关闭的时间离散问题

原论文在 Eq. (59) 后明确说明 temporal multivariable differential problem 使用 discontinuous Galerkin order zero。

当前 tower fixed-basis update 使用已经从一维版本继承并验证过的 backward-Euler-like sequential temporal discretization。

二者严格等价性仍未数学证明。

因此 FOM-3D 效率诊断期间仍不应为了性能原因顺便修改时间离散理论。

---

## 42. 当前阶段没有改变的研究原则

1. 原论文是理论基准；
2. 一维三材料杆验证仍然是 tower 扩展的重要数值依据；
3. 当前阶段先诊断，不做未经验证的优化；
4. 不把低 rank 直接解释为 speedup；
5. 不把当前 wall-time 结果泛化成“LATIN-PGD 方法本身低效”；
6. 当前结果只说明“当前 whole-history transactional tower implementation”存在规模扩展效率问题；
7. 不在 FOM-3D-5 关闭前直接推进 100-cycle whole-history LATIN；
8. profiling instrumentation 必须首先证明不会改变数值路径。

---

## 43. 当前最重要的研究认识

经过本阶段推导和讨论，现在可以把此前较模糊的“LATIN 为什么慢”转化成明确的多层调用结构：

`outer LATIN`

`-> Trial A full temporal update`

`-> if enrichment`

`-> fixed-point`

`-> residual-LS spatial solve`

`-> LSMR inner iterations`

`-> whole-history C/C^T actions`

`-> single-mode temporal update`

`-> fixed-point convergence`

`-> enlarged-basis full temporal re-optimization`

`-> Trial B`

这说明当前 wall-time 不由一个单独的 rank 数决定，而由多个嵌套重复次数共同决定。

最核心的概念可浓缩为：

> **rank 告诉我们最终用了多少个 mode；fixed-point count 告诉我们每个新 mode 空间与时间来回修正了多少轮；LSMR count 告诉我们每一轮空间最小二乘内部又做了多少 Krylov iterations；Nt 决定每一次 whole-history operator action 的时间域规模。**

---

## 44. 下一次工作开始时应恢复的 checkpoint

下一次继续工作时，应从以下结论开始：

> **FOM-3C 精度验证已经关闭。**

> **FOM-3D 1/2/5/10-cycle 基础 wall-time scaling 已完成并正式总结。**

> **FOM-3D-5 的源码级 call-chain、temporal update、enrichment、fixed-point、residual-LS 与 LSMR 机制审计已经完成。**

当前已知：

- 10-cycle：`Nt=401`
- `Nq=320`
- rank `21`
- outer iterations `39`
- Trial A `39`
- successful enrichments / Trial B `21`
- total trials `60`
- `final xi = 8.941607234831e-06`
- LATIN mean wall time `467.361528 s`
- FOM mean wall time `307.891095 s`

当前最重要但尚未测出的量：

1. 每个 enrichment 的真实 `fixed_point_iterations`；
2. 每个 residual-LS spatial solve 的真实 `LSMR iterations`；
3. 这些 inner counts 随 cycle count 是否增长；
4. temporal update 与 enrichment/LSMR 各自在总 wall time 中的实际占比。

因此下一项正式操作应为：

> **建立外部 profiling recorder，不修改 `latin/` 核心代码；首先运行 1-cycle，记录全部 enrichment fixed-point history 与 LSMR iteration history，并验证数值 invariants 完全不变。**

在该 measurement 通过以前：

> **不修改 residual-LS 求解策略，不改成直接 normal-equation solver，也不直接运行 100-cycle LATIN。**

---

## 45. 当前阶段最终判定

截至本总结，可以正式记录：

> **FOM-3D-5 源码级成本来源审计与 residual-LS / LSMR 机制梳理：PASS。**

同时必须记录：

> **函数级 wall-time profiling 尚未完成，因此 FOM-3D-5 整体阶段尚未关闭。**

当前研究问题已经从：

> 为什么 10-cycle LATIN wall time 增长到约 467 s？

进一步具体化为：

> **10-cycle 的约 467 s 中，repeated full-basis temporal update 与 enrichment 内部 fixed-point / LSMR 分别占多少；这些成本是因为调用次数增加、单次 inner iteration 数增加、whole-history `Nt` 增长，还是三者共同作用？**

下一阶段不再依靠推测，而将通过外部 recorder 和函数级 timing 将这一问题转化为可量化数据。

