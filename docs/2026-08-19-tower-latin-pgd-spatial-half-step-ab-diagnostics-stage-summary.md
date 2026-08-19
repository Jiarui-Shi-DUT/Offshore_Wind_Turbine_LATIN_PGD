# Tower LATIN-PGD 第4模态：1D-style in-loop orthogonalisation 与 direct weighted-LS spatial half-step A/B/C 诊断阶段总结

**日期：2026-08-19**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**仓库：`Jiarui-Shi-DUT/Offshore_Wind_Turbine_LATIN_PGD`**  
**分支：`feature/offshore-wind-turbine-tower-fatigue`**  
**理论主线：Bhattacharyya et al. 原论文 $x-t$ LATIN-PGD → fiber beam-column offshore wind turbine tower**  
**本阶段性质：diagnostic-only；不修改 `latin/` 核心算法，不绕过 fixed-point convergence gate，不改变 transaction semantics。**

---

# 1. 本阶段起点

上一阶段已经针对 tower reversed benchmark 中第 4 个 PGD enrichment mode 的 `fixed_point_not_converged` 完成了以下诊断：

- iteration-cap 排查；
- complete-pair period-3 orbit 确认；
- constant spatial under-relaxation 排查；
- raw-map defect 诊断；
- whole-time coupled BE temporal minimisation 诊断；
- unconverged period-3 phase usefulness 诊断；
- top-10 residual-energy deterministic seed / basin-sensitivity 诊断。

此前结果共同表明：

```text
current tower fourth-mode raw fixed-point map
→ ordinary fixed point does not converge
→ complete pair approaches a robust period-3 orbit
→ three phases remain genuinely useful after diagnostic post-basis management
→ the failure is not explained by iteration cap, trivial basis saturation,
  constant damping, sequential-BE scope, or the tested residual-row seed family
```

在此基础上，本阶段暂时不优先恢复原论文 DG0，而采用一个更贴近当前研究目标的原则：

> 先把已经成功复现的一维三材料杆 LATIN-PGD numerical implementation 作为 tower v1 的 numerical reference，逐项检查 1D 与 tower enrichment fixed-point map 的实质差异。

这一选择不改变理论来源边界：

- continuous LATIN-PGD structure 仍以 Bhattacharyya et al. 原论文为理论基准；
- current tower temporal discretisation 仍明确属于已验证 1D implementation 所采用的 backward-Euler engineering choice；
- 不将 current BE 声称为 paper-exact DG0；
- DG0 fidelity verification 暂缓，不作为当前 tower v1 migration 的阻塞项。

---

# 2. 1D 与 tower enrichment 的第一项差异：fixed-point 内 orthogonalisation

current validated 1D `latin/pgd_enrichment.py` 在两个位置对 candidate spatial mode 做已有基底消除：

```text
residual-driven initial seed
→ orthogonalise against existing spatial PGD basis
→ volume-normalise
```

以及：

```text
each fixed-temporal spatial update
→ solve spatial weighted LS problem
→ orthogonalise against existing spatial PGD basis
→ volume-normalise
```

因此，在 1D enrichment fixed-point 内，candidate spatial mode 始终被限制在 existing basis 的加权正交补中。

而 current tower `latin/tower_pgd_enrichment.py` 的 raw fixed-point map 采用：

```text
seed
→ M-normalise only
→ temporal solve
→ tower Eq. (70)-(71) spatial solve
→ M-normalise only
→ temporal solve
→ complete-pair convergence check
```

existing-basis M-weighted Gram-Schmidt 被放在 raw fixed point 收敛之后。

因此，第一个 A/B 问题为：

> tower 第 4 mode 的 period-3 orbit 是否主要来自 fixed-point 内没有排除 existing rank-3 spatial subspace？

---

# 3. In-loop orthogonalisation diagnostic 的设计

新增 diagnostic script：

```text
examples/tower_latin_pgd_fourth_mode_inloop_orthogonalization_probe.py
```

共同冻结：

```text
benchmark
Trial-A state
rank-3 persistent basis
shifted defect
H_sigma
material-point metric M
tower equilibrium operator
tower Eq. (70)-(71) spatial half-step
backward-Euler Eq. (72) temporal half-step
complete-pair convergence criterion
fixed-point tolerance = 1e-6
max fixed-point iterations = 200
```

A/B 唯一区别：

```text
A) current tower raw map
   no in-loop basis orthogonalisation

B) same raw map
   + 1D-style M-orthogonalisation of initial seed
   + 1D-style M-orthogonalisation after every spatial half-step
```

B 中每次 orthogonalisation 后重新通过 tower reference equilibrium operator 构造 associated stress mode。

---

# 4. In-loop orthogonalisation A/B 数值结果

共同 baseline：

```text
termination = enrichment_failed
committed   = 7
attempted   = 8
rank        = 3
xi          = 8.242691499e-04
```

failing Trial-A relative reduced residual：

```text
9.825719402e-01
```

## 4.1 A：current tower raw map

```text
converged           = False
iterations          = 200
last chi            = 6.034870798845e-01
lag-3 distance      = 2.093710261832e-06
max basis overlap   = 3.387997185086e-01
```

tail：

```text
0.577410350
0.658244303
0.603487691
0.577410569
0.658244106
0.603487458
0.577410759
0.658243935
0.603487256
0.577410925
0.658243787
0.603487080
```

即：

$$ \chi^{(k)} \approx 0.57741,\ 0.65824,\ 0.60349,\ 0.57741,\ 0.65824,\ 0.60349,\ldots $$

并且：

$$ d\!\left(z^{k},z^{k-3}\right)\approx2.09\times10^{-6} $$

再次确认原 period-3 orbit。

## 4.2 B：1D-style in-loop M-orthogonalised map

```text
converged           = False
iterations          = 200
last chi            = 6.062615795343e-01
lag-3 distance      = 9.665845681520e-08
seed novelty        = 9.299125524165e-01
min novelty         = 9.299125524165e-01
final novelty       = 9.829716829867e-01
max basis overlap   = 5.551115123126e-16
```

tail：

```text
0.568634869
0.658373297
0.606261618
0.568634876
0.658373288
0.606261603
0.568634881
0.658373281
0.606261590
0.568634886
0.658373275
0.606261580
```

即：

$$ \chi^{(k)} \approx 0.56863,\ 0.65837,\ 0.60626,\ 0.56863,\ 0.65837,\ 0.60626,\ldots $$

并且：

$$ d\!\left(z^{k},z^{k-3}\right)\approx9.67\times10^{-8} $$

---

# 5. In-loop orthogonalisation 假设的结论

B 中：

$$ \max_j\left|\langle p_{\rm candidate},p_j\rangle_M\right|=O(10^{-16}) $$

说明 candidate spatial mode 已在数值精度内与 existing rank-3 basis 正交。

但 fixed point 仍不收敛，并且形成更清晰的 period-3 orbit。

因此：

$$ \boxed{\text{lack of in-loop basis orthogonalisation is not the principal cause of the fourth-mode period-3 failure}} $$

这一假设可以关闭。

更具体地说：

- original map 中确实存在一定 existing-basis overlap，最大约为 0.339；
- 1D-style orthogonalisation 成功将 overlap 降到 machine precision；
- 但 ordinary fixed-point convergence 没有恢复；
- lag-3 distance 反而从 $O(10^{-6})$ 降至 $O(10^{-8})$；
- 因此 period-3 attractor 对 existing-basis contamination 并不敏感。

这不意味着 in-loop orthogonalisation 在一般情况下没有数值价值，只意味着它不是当前 failing fourth-mode subproblem 的主因。

---

# 6. 第二项实质差异：spatial half-step formulation

1D 与 tower enrichment 的更本质差异在 spatial half-step。

## 6.1 Current validated 1D spatial half-step

对于固定 temporal pair：

$$ \lambda_n,\qquad \dot\lambda_n $$

定义 spatial stress linear map：

$$ s=A_\sigma p $$

则离散 mechanical residual 可写为：

$$ r_n(p)=\dot\lambda_n p-H_{\sigma,n}\lambda_nA_\sigma p+\Delta_n $$

即：

$$ r_n(p)=B_n p+\Delta_n $$

其中：

$$ B_n=\dot\lambda_n I-\lambda_nD_{H,n}A_\sigma $$

current 1D `pgd_enrichment.py` 直接在空间变量上最小化加权 space-time residual：

$$ J_h(p)=\frac12\sum_n w_n\,r_n(p)^TMD_{H,n}^{-1}r_n(p) $$

对应 normal equation：

$$ \left(\sum_n w_n B_n^TMD_{H,n}^{-1}B_n\right)p=-\sum_n w_n B_n^TMD_{H,n}^{-1}\Delta_n $$

随后再对 existing spatial basis 做 weighted orthogonalisation 与 normalization。

## 6.2 Current tower spatial half-step

tower current implementation 采用此前由 Bhattacharyya et al. Eq. (65)-(71) 推导出的 structural Galerkin route。

固定 temporal function 后构造：

$$ a_h=\sum_{n=1}^{N}\Delta t_n\dot\lambda_n\lambda_n $$

$$ A_q=\sum_{n=1}^{N}\Delta t_nH_{\sigma,nq}\lambda_n^2 $$

$$ \bar\delta_q=\sum_{n=1}^{N}\Delta t_n\Delta_{nq}\lambda_n $$

定义：

$$ W_q=\left(A_q+\frac{a_h}{C_{0,q}}\right)^{-1} $$

然后求：

$$ H^TMD_WH\,\bar{\tilde U}=-H^TMD_W\bar\delta $$

再恢复 compatible strain、stress mode 与 plastic spatial mode。

因此：

```text
1D spatial half-step
= direct weighted residual minimisation

tower spatial half-step
= Eq. (70)-(71) structural Galerkin elimination
```

---

# 7. Direct weighted-LS spatial A/B/C diagnostic 的设计

新增 diagnostic script：

```text
examples/tower_latin_pgd_fourth_mode_direct_ls_spatial_probe.py
```

比较：

```text
A) current tower Eq.(70)-(71)
   no in-loop orthogonalisation

B) tower Eq.(70)-(71)
   + 1D-style in-loop M-orthogonalisation

C) literal 1D-style direct weighted-LS spatial half-step
   + same in-loop M-orthogonalisation as B
```

B 与 C 共同保持：

```text
same residual-driven seed
same existing rank-3 basis
same in-loop M-orthogonalisation
same shifted defect
same H_sigma
same M
same backward-Euler temporal solve
same complete-pair chi
same fixed-point tolerance
same maximum iteration count
```

因此 B → C 隔离的主要变化是：

```text
tower Eq.(70)-(71) spatial half-step
        ↓
literal 1D-style direct weighted residual minimisation
```

为了在 coarse diagnostic benchmark 中忠实构造 1D-style normal equation，C 临时显式构造 material-point stress linear map：

$$ s=A_\sigma p $$

其中 material-point count 为：

```text
Nq = 320
```

dense mapping 仅用于 diagnostic，不代表 production tower solver 的推荐实现。

---

# 8. A/B/C 数值结果

共同条件：

```text
baseline termination             = enrichment_failed
committed LATIN iterations       = 7
attempted LATIN iterations       = 8
persistent rank                  = 3
persistent xi                    = 8.242691499e-04
failing Trial-A relative residual= 9.825719402e-01
max fixed-point iterations       = 200
fixed-point tolerance            = 1.0e-06
```

## 8.1 A：current tower Eq. (70)-(71)

```text
converged           = False
iterations          = 200
last chi            = 6.034870798845e-01
lag-3 distance      = 2.093710259844e-06
max basis overlap   = 3.387997185086e-01
final ||p||_M       = 1.000000000000e+00
```

保持原 period-3。

## 8.2 B：tower Eq. (70)-(71) + in-loop orthogonalisation

```text
converged           = False
iterations          = 200
last chi            = 6.062615795343e-01
lag-3 distance      = 9.665845681520e-08
max basis overlap   = 5.551115123126e-16
final ||p||_M       = 1.000000000000e+00
```

仍保持 period-3。

## 8.3 C：literal 1D-style direct weighted-LS spatial half-step

```text
converged           = True
iterations          = 27
last chi            = 7.552396390751e-07
lag-3 distance      = 3.702454207154e-06
seed novelty        = 9.299125524165e-01
min novelty         = 9.299125524165e-01
final novelty       = 9.615356082090e-01
max basis overlap   = 5.412337245048e-16
final ||p||_M       = 1.000000000000e+00
```

最后 12 个 complete-pair changes：

```text
8.5866e-05
5.5837e-05
3.6310e-05
2.3613e-05
1.5355e-05
9.9860e-06
6.4940e-06
4.2230e-06
2.7460e-06
1.7860e-06
1.1610e-06
7.5524e-07
```

因此：

$$ \chi_{\rm fp}^{(27)}=7.5524\times10^{-7}\lt10^{-6} $$

满足 ordinary fixed-point convergence gate。

而且 tail 呈稳定衰减，不是单次偶然跨阈值。

---

# 9. 本阶段最关键结论

在当前 failing fourth-mode frozen subproblem 中：

$$ \boxed{\text{changing only the spatial half-step from tower Eq. (70)-(71) to literal 1D-style direct weighted LS restores ordinary fixed-point convergence}} $$

这是目前所有诊断中，第一个被严格隔离且能够解除 period-3 failure 的算法差异。

因此当前排除链更新为：

```text
iteration cap
    ×

scalar-chi artefact
    ×

constant spatial under-relaxation
    ×

sequential-BE scope
    ×

basis saturation / mode insignificance
    ×

high-energy residual-row seed sensitivity
    ×

lack of in-loop orthogonalisation
    ×

spatial half-step formulation
    ✓ first isolated difference that removes period-3
```

---

# 10. 但不能直接得出“Eq. (70)-(71) 错了”

本结果只证明：

> 对当前 discrete tower enrichment problem，current Eq. (70)-(71) spatial map 与 current BE temporal map 组合形成 period-3，而 direct weighted-LS spatial map 与同一个 BE temporal map 可以收敛。

它尚不能区分：

1. tower Eq. (70)-(71) 代码实现存在 algebraic inconsistency；
2. Eq. (70)-(71) implementation 本身正确，但与 current project-BE Eq. (72) temporal discretisation 在离散层面不构成一致的 alternating minimisation；
3. 两种因素同时存在。

因此：

$$ \boxed{\text{do not modify core tower enrichment yet}} $$

下一阶段必须先做 derivation audit。

---

# 11. 一个值得检验的工作假设：common discrete objective consistency

Case C 的两个 half-steps 都可理解为对同一个 discrete residual objective 做交替最小化。

定义：

$$ J_h(p,\lambda)=\frac12\sum_n w_n\,r_n(p,\lambda)^TMD_{H,n}^{-1}r_n(p,\lambda) $$

其中：

$$ r_n(p,\lambda)=\dot\lambda_n p-D_{H,n}\lambda_nA_\sigma p+\Delta_n $$

固定 $\lambda$ 时，C 求：

$$ p^{k+1}=\operatorname*{arg\,min}_p J_h(p,\lambda^k) $$

固定 $p$ 时，current BE temporal solve 求：

$$ \lambda^{k+1}\approx\operatorname*{arg\,min}_{\lambda\ {\rm under\ current\ BE\ sequential\ treatment}} J_h(p^{k+1},\lambda) $$

因此 C 更接近一个离散 alternating least-squares / alternating minimisation map。

而 current tower B 使用：

```text
spatial half-step
= paper-derived Eq.(70)-(71) Galerkin route

temporal half-step
= project backward-Euler residual minimisation
```

因此一个新的、尚待证明的解释是：

> current tower mixed discrete map 的两个 half-steps 可能并非同一个 discrete objective 的 coordinate minimisers，因此不必具有单调下降或 contraction 性质，并可能出现 stable period-3 orbit。

这是数值结果支持的工作假设，不是已经完成的数学证明。

---

# 12. 下一阶段唯一优先任务：spatial-half-step consistency audit

下一阶段暂不增加新的数值 probe，也暂不修改 `latin/` core。

优先完成：

## 12.1 从 current BE residual 推导 direct-LS spatial normal equation

从：

$$ r_n=\dot\lambda_np-D_{H,n}\lambda_nA_\sigma p+\Delta_n $$

和：

$$ J_h=\frac12\sum_nw_n\,r_n^TMD_{H,n}^{-1}r_n $$

严格推导：

$$ \frac{\partial J_h}{\partial p}=0 $$

并确认 C 的代码是否逐项等于该离散 stationarity condition。

## 12.2 重新展开 tower Eq. (65)-(71)

从原论文 spatial Galerkin equation 出发，重新恢复：

- temporal contractions；
- $W_q$；
- $\bar\delta$；
- compatible mode；
- equilibrium constraint；
- Eq. (70) structural solve；
- Eq. (71) plastic mode recovery。

## 12.3 B 与 C 逐项比较

重点比较以下时间系数和 metric structure：

```text
dot(lambda)^2
lambda*dot(lambda)
lambda^2
H_sigma
H_sigma^{-1}
shifted defect coupling
time quadrature
initial-time contribution
```

目标不是要求 B 与 C 必须相同，而是确定：

> 它们从哪一步开始代表不同的 discrete stationarity condition。

## 12.4 再决定 tower v1 的路线

可能路线 A：

```text
paper-fidelity first
→ retain Eq.(70)-(71)
→ recover a temporally consistent discretisation
→ possibly return to exact DG0 audit
```

可能路线 B：

```text
migration-stability first
→ explicitly adopt validated 1D-style discrete alternating minimisation
→ generalise it to tower material points
→ label it as tower-v1 engineering discretisation
→ retain paper Eq.(70)-(71) as later fidelity comparison
```

在 derivation audit 完成前，不应提前选择。

---

# 13. 本阶段新增 diagnostic files

```text
examples/tower_latin_pgd_fourth_mode_inloop_orthogonalization_probe.py
examples/tower_latin_pgd_fourth_mode_direct_ls_spatial_probe.py
```

两者均为 diagnostic-only。

不得将其视为 persistent production algorithm 的正式修改。

---

# 14. 本阶段核心算法保持不变

本阶段没有修改：

```text
latin/tower_pgd_enrichment.py
latin/tower_pgd_time_update.py
latin/tower_latin_pgd_solver.py
latin/tower_global_finishing.py
latin/tower_equilibrium_operator.py
```

没有：

- 绕过 `fixed_point_not_converged`；
- 接受 unconverged period-3 phase；
- 改变 Trial-A / Trial-B transaction semantics；
- 改变 persistent basis；
- 改变 LATIN convergence indicator；
- 改变 saturation rule；
- 引入 Anderson / Aitken / cycle-aware acceptance；
- 引入 DG0；
- 引入 cycle-phase 或 multi-time-scale PGD。

---

# 15. 当前冻结结论

本阶段建议冻结以下结论。

**结论 1**

1D-style in-loop basis orthogonalisation 可以有效消除 candidate 对 existing rank-3 basis 的 overlap，但不能消除第四模态 period-3 fixed-point orbit。

**结论 2**

因此 existing-basis contamination 不是当前 fourth-mode fixed-point failure 的主要原因。

**结论 3**

将 tower Eq. (70)-(71) spatial half-step 替换为 literal 1D-style direct weighted residual minimisation，在保持同一 BE temporal half-step、同一 in-loop orthogonalisation 和同一 complete-pair convergence gate 的情况下，可以在 27 iterations 内恢复 ordinary fixed-point convergence。

**结论 4**

spatial-half-step formulation 是目前第一个被隔离出来、能够解除 period-3 failure 的算法差异。

**结论 5**

当前结果不能直接证明 tower Eq. (70)-(71) 理论或代码错误；下一阶段必须首先区分 implementation error 与 discrete-formulation inconsistency。

**结论 6**

在完成 spatial-half-step derivation audit 前，不修改 core tower enrichment。

---

# 16. Git checkpoint 建议

本阶段建议提交：

```text
docs/2026-08-19-tower-latin-pgd-spatial-half-step-ab-diagnostics-stage-summary.md
examples/tower_latin_pgd_fourth_mode_inloop_orthogonalization_probe.py
examples/tower_latin_pgd_fourth_mode_direct_ls_spatial_probe.py
```

建议 commit message：

```text
test: isolate tower fourth-mode spatial half-step mismatch
```

在 commit 前先完成：

```text
python -m py_compile examples/tower_latin_pgd_fourth_mode_inloop_orthogonalization_probe.py
python -m py_compile examples/tower_latin_pgd_fourth_mode_direct_ls_spatial_probe.py
```

然后再按项目 Git workflow 逐步保存。

---

# 17. 下一阶段标题建议

```text
Tower LATIN-PGD spatial-half-step consistency audit:
paper Eq. (65)-(71) Galerkin formulation
vs
validated 1D-style discrete residual minimisation
```

下一阶段只做理论与离散一致性审计，不立即改 production code。
