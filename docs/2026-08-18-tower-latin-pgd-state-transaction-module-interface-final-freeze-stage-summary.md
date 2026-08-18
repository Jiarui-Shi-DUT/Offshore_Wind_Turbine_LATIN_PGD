# 2026-08-18 Tower LATIN-PGD State Ownership、Transaction Semantics 与 Module Interface Final Freeze 阶段总结

- 项目：Offshore Wind Turbine and LATIN-PGD
- 仓库：`Jiarui-Shi-DUT/Offshore_Wind_Turbine_LATIN_PGD`
- 分支：`feature/offshore-wind-turbine-tower-fatigue`
- 日期：2026-08-18
- 阶段性质：code-before-final-interface-freeze
- 上一阶段文档：`docs/2026-08-18-tower-latin-state-field-by-field-contract-stage-summary.md`
- 当前阶段范围：在 14-entry state contract 已冻结后，闭合四层 state ownership、Trial A / Trial B transaction、mode rollback、LATIN iteration atomic commit、module I/O boundary、dependency direction 与 first-code-entry order
- 当前结论：代码前的理论、数据、transaction 与 module interface 已完成最终冻结
- 本次修订：区分 persistent `B_m`、Trial A fixed-basis temporal-update result `B_m^A`、Trial B enriched candidate `B_(m+1)^(B*)`；明确 mode rollback target 与 whole-iteration rollback target 不同
- 下一阶段：正式创建 `latin/tower_state.py` 与第一组 state/layout 单元测试
- PyCharm Preview 兼容性：所有 display-math 统一采用单物理行 `$$ ... $$`；数学环境不使用 raw `<`，小于关系统一使用 `\lt`

---

# 1. 本阶段的准确起点

上一阶段已经冻结 future tower state 的 canonical storage 与逐字段行为。

material-point coordinate：

$$ q\leftrightarrow(e,g,f). $$

除 `time` 外，全部 material histories 的 canonical shape：

$$ (N_t,N_q). $$

future tower state 的 13 个 material-point fields 为：

```text
primary
    plastic_strain_rate
    elastic_strain
    alpha_rate
    r_bar_rate
    damage_rate
    stress
    beta
    R_bar
    energy_release_rate

support histories
    plastic_strain
    alpha
    r_bar
    damage
```

加上：

```text
time
```

形成 14-entry conceptual contract。

本阶段不再修改这些字段本身的理论含义，而是回答：

> 谁拥有 accepted state？谁可以生成 trial？basis 在什么时候可以改变？Trial A 与 Trial B 是否可以 chain？enrichment 失败时到底 rollback 什么？最终哪些 module 接受什么 input、返回什么 output、绝不能修改什么？

---

# 2. 本阶段继续采用的四类来源

## 2.1 原论文明确内容

包括：

- LATIN local/global alternating structure；
- Eq. (58)–(59) fixed-basis PGD temporal update；
- Eq. (60) saturation parameter；
- Eq. (61)–(72) rank-one enrichment；
- post-enrichment Gram–Schmidt 与 time-function modification；
- Eq. (73)–(75) hardening/damage finishing；
- relaxation；
- Eq. (76)–(77) convergence indicator。

## 2.2 严格数学推导

包括：

- common-baseline Trial A / Trial B 对比的必要性；
- linear equilibrium / compatibility 在 convex relaxation 下保持；
- exact basis-coordinate transformation；
- mode acceptance 与 full LATIN state acceptance 的层级差异；
- atomic commit 的 state–basis consistency requirement。

## 2.3 current 1D implementation

关键参考：

```text
latin/state.py
latin/pgd_basis.py
latin/pgd_time_update.py
latin/pgd_enrichment.py
latin/pgd_global_stage.py
latin/iteration_control.py
latin/pgd_saturation.py
latin/pgd_solver.py
```

current 1D 已经提供：

- `LatinState` deep copy；
- `PGDBasis1D` 与 state 分离；
- existing basis temporal update；
- residual-driven enrichment；
- field-wise relaxation；
- Eq. (76)–(77) indicator；
- saturation logic；
- accepted-iteration stagnation control。

## 2.4 tower v1 engineering specification

本阶段新增并冻结：

- `s_i`、`B_m`、`\xi_i` 为 current LATIN iteration persistent immutable baseline；
- local snapshot、search directions、damage structural data 在 same iteration Trial A/B 间 frozen；
- Trial A 与 Trial B 必须从 same baseline 独立构造；
- candidate basis必须在独立 mode transaction 内工作；
- mode acceptance 与 LATIN iteration commit分离；
- normal iteration commit为 state + basis + indicator atomic commit；
- enrichment真正被要求但无法产生合法 mode时，不把 provisional Trial A 静默升级为 accepted state；
- reduced residual已经 solved 时，不制造 false enrichment failure；
- one accepted new pair per enrichment event；
- solver是唯一 persistent transaction owner；
- bottom modules均采用 pure builder / evaluator semantics。

---

# 3. 四层 LATIN state 再确认

上一 accepted relaxed global state：

$$ s_i. $$

current local half-step：

$$ \hat s_{i+1/2}. $$

complete unrelaxed global candidate：

$$ \breve s_{i+1}. $$

relaxed next global state：

$$ s_{i+1}. $$

正常单次 iteration 的状态链：

```text
accepted s_i
    ↓
local stage
    ↓
hat{s}_{i+1/2}
    ↓
global PGD + finishing
    ↓
breve{s}_{i+1}
    ↓
relaxation
    ↓
trial s_(i+1)
    ↓
commit decision
```

关键点：

> 只有经过 commit decision 的 relaxed state 才能成为下一轮 persistent `s_(i+1)`。

---

# 4. Persistent transaction baseline

第 $i$ 次 LATIN iteration 的 persistent accepted baseline 定义为：

$$ \mathcal T_i^{\rm base}=(s_i,\mathcal B_m,\xi_i). $$

其中：

- $s_i$：上一 accepted relaxed global state；
- $\mathcal B_m$：当前 accepted persistent PGD basis；
- $\xi_i$：上一 accepted LATIN indicator。

本阶段冻结：

> current iteration 未完成前，三者均为逻辑只读。

---

# 5. 为什么 basis 与 state 必须处于同一 persistent level

new mode accepted 后可能发生：

```text
append new spatial mode
modify old temporal functions
re-optimise all temporal functions
```

若 persistent basis被原地修改，而随后 Trial B失败，则：

```text
state
    still corresponds to old basis

basis
    already corresponds to new reduced solution
```

会造成 persistent pair不一致。

因此：

> accepted state 与 accepted basis必须作为一个 coupled persistent snapshot理解。

---

# 6. Algorithmic immutability 与 Python immutability 的区别

本阶段冻结的是：

> algorithmic immutability。

也就是任何 module都不允许把 accepted object当作 scratch buffer。

禁止语义：

```text
modify accepted state
try computation
undo on failure
```

采用：

```text
read accepted object
build new result
validate
commit or discard
```

至于 future code 是否使用：

`@dataclass(frozen=True)`

或者 NumPy read-only flags，是 implementation detail。

第一版建议优先 correctness，再根据 memory cost优化。

---

# 7. Local snapshot ownership

由：

$$ s_i $$

构造：

$$ \hat s_{i+1/2}. $$

其生命周期：

```text
iteration i starts
    ↓
build local snapshot
    ↓
freeze
    ↓
used by Trial A
    ↓
possibly reused by Trial B
    ↓
iteration ends
```

它不是 persistent accepted state。

---

# 8. Trial A 与 Trial B 使用同一 local state

冻结：

$$ \hat s_{i+1/2}^{A}=\hat s_{i+1/2}^{B}. $$

same iteration 中 enrichment发生后，不重新执行 local stage。

否则 Trial A 与 Trial B 不再是在：

> same nonlinear local problem

下比较不同 reduced basis。

---

# 9. Search directions 也必须 frozen

定义：

$$ \mathcal H_i=\{H_{\sigma,i},H_{\beta,i},H_{\bar R,i}\}. $$

same iteration 中：

$$ \mathcal H_i^{A}=\mathcal H_i^{B}. $$

禁止：

```text
Trial A
    ↓
new stress
    ↓
recompute search directions
    ↓
Trial B
```

否则 basis effect 与 search-direction change混在一起。

---

# 10. Damage structural data 也属于 frozen iteration data

damage residual strain：

$$ \Delta\varepsilon_i^R $$

damage compatible correction：

$$ \Delta\tilde\varepsilon_i $$

damage stress correction：

$$ \Delta\tilde\sigma_i $$

均由 baseline/local data决定。

same iteration中：

$$ \Delta\varepsilon_A^R=\Delta\varepsilon_B^R, $$

$$ \Delta\tilde\varepsilon_A=\Delta\tilde\varepsilon_B, $$

$$ \Delta\tilde\sigma_A=\Delta\tilde\sigma_B. $$

---

# 11. Frozen iteration data 的 conceptual package

本阶段定义 conceptual object：

$$ \mathcal F_i=(\hat s_{i+1/2},\mathcal H_i,\Delta\varepsilon_i^R,\Delta\tilde\varepsilon_i,\Delta\tilde\sigma_i,f_i). $$

其中：

$$ f_i $$

是 full plastic global-stage forcing。

这一 package在 same iteration Trial A/B间只读。

未来软件中对应 `FrozenGlobalData` + local state + search directions。

---

# 12. Fixed-basis Trial A

使用 accepted persistent basis：

$$ \mathcal B_m $$

作为 read-only input，执行 fixed-basis temporal update。

该步骤保持 spatial rank 与 spatial modes 不变，但更新 existing modes 的 temporal coordinates，因此返回新的 provisional fixed-basis representation：

$$ \mathcal B_m \longrightarrow \mathcal B_m^{A}. $$

其中 spatial quantities 满足：

$$ P_m^{A}=P_m,\qquad S_m^{A}=S_m, $$

而一般：

$$ \Lambda_m^{A}\neq\Lambda_m,\qquad \dot\Lambda_m^{A}\neq\dot\Lambda_m. $$

因此 Trial A 的 plastic correction 必须由：

$$ \mathcal B_m^{A} $$

重构，而不是继续由 old persistent：

$$ \mathcal B_m $$

重构。

由 $\mathcal B_m^{A}$ 得到：

$$ \Delta\varepsilon_A^p, $$

$$ \Delta\dot\varepsilon_A^p, $$

$$ \Delta\sigma'_A, $$

$$ \Delta\varepsilon'_A. $$

完成 global finishing：

$$ \breve s_{i+1}^{A}. $$

再 relaxation：

$$ s_{i+1}^{A}=(1-\mu)s_i+\mu\breve s_{i+1}^{A}. $$

因此 Trial A 实际上是一个 coupled provisional pair：

$$ \mathcal T_A^{\rm trial}=(s_A,\mathcal B_m^{A},\xi_A). $$

这里 $\mathcal B_m^{A}$ 仍未 persistent commit。

---

# 13. Trial A indicator

定义：

$$ \xi_A=\xi(\hat s_{i+1/2},s_{i+1}^{A}). $$

以及：

$$ \zeta_A=\frac{\xi_i-\xi_A}{\xi_i+\xi_A}. $$

reference：

$$ \xi_i $$

在 same iteration 内必须 frozen。

---

# 14. Trial A 不是 persistent state

即使 Trial A完成了：

- complete unrelaxed candidate；
- relaxation；
- Eq. (76)；
- Eq. (60)；

它仍是：

> provisional complete trial。

只有 decision logic明确允许 commit时才成为 accepted state。

---

# 15. 为什么 Trial A 不能先 commit再 enrichment

如果先：

$$ s_i\leftarrow s_A, $$

随后 enrichment：

Trial B将不得不从：

$$ s_A $$

出发。

这破坏：

$$ \text{same baseline Trial A / Trial B}. $$

因此：

> enrichment request出现时，Trial A只 hold，不 commit。

---

# 16. Mode transaction 的 baseline

必须区分两层 baseline。

LATIN iteration persistent rollback target 仍是：

$$ \mathcal B_m. $$

但是在 fixed-basis temporal update完成后，current reduced defect $R_A$ 对应的是 provisional fixed-basis solution：

$$ \mathcal B_m^{A}. $$

因此若 enrichment 被真正要求，one-mode transaction 的直接 input / snapshot 应为：

$$ \mathcal B_m^{A}. $$

工作对象：

$$ \mathcal B_{\rm work}. $$

所有：

- raw rank-one generation；
- Gram–Schmidt；
- temporal coordinate transformation；
- all-mode temporal re-optimisation；

均作用于 $\mathcal B_m^{A}$ 的 working copy。

因此本阶段冻结两个不同 rollback targets：

```text
mode-transaction rollback target
    = B_m^A

LATIN-iteration rollback target
    = persistent B_m
```

前者用于撤销 failed/rejected new-mode construction；后者用于撤销整个 failed LATIN iteration。

---

# 17. Mode transaction 的 accepted candidate

若 new pair通过：

- fixed-point validity；
- spatial novelty；
- modified-time significance；
- field invariance；
- all-mode temporal solve；
- full residual benefit；

则由：

$$ \mathcal B_m^{A} $$

得到 enlarged candidate：

$$ \mathcal B_{m+1}^{B*}. $$

星号表示：

> mode transaction内部 accepted，但还没有 persistent commit。

其中 existing modes 与 new mode 已经完成 enlarged-basis all-mode temporal re-optimisation。

---

# 18. Mode acceptance 与 LATIN iteration acceptance 不同

mode acceptance衡量：

> enlarged reduced basis 是否改善 current reduced plastic global problem。

iteration acceptance衡量：

> complete relaxed LATIN state是否可以成为下一 persistent nonlinear iterate。

因此：

$$ \text{mode accepted}\not\Rightarrow\text{iteration already committed}. $$

---

# 19. Trial B 从 same baseline重建

Trial B 的 **state baseline** 与 Trial A 完全相同，仍使用：

$$ s_i, $$

$$ \hat s_{i+1/2}, $$

$$ \mathcal H_i, $$

$$ \mathcal F_i. $$

但它使用 enriched candidate basis：

$$ \mathcal B_{m+1}^{B*}. $$

重新生成：

$$ \breve s_{i+1}^{B}. $$

禁止：

$$ \breve s_B=s_A+\text{new mode correction}. $$

也禁止把 $s_A$ 当作 Trial B 的新的 nonlinear state baseline。

因此：

> Trial A/B share the same LATIN state baseline $s_i$，但 reduced-basis working path 是 $\mathcal B_m\rightarrow\mathcal B_m^A\rightarrow\mathcal B_{m+1}^{B*}$。

---

# 20. Trial A/B common-baseline graph

正确关系：

```text
persistent accepted basis
        B_m
         │
         │ fixed-basis temporal update
         ▼
       B_m^A ─────────────── Trial A
         │
         │ residual R_A
         │ one-mode enrichment
         ▼
    B_(m+1)^(B*) ─────────── Trial B
```

与此同时，两个 complete LATIN trials 都从 same accepted state baseline 构造：

```text
accepted s_i ─────────────── Trial A
      │
      └───────────────────── Trial B
```

禁止的 state chaining：

```text
s_i
 ↓
Trial A
 ↓
Trial B
```

因此必须同时记住：

> common-baseline 指的是 LATIN state / local / search-direction / damage-global data；basis 则存在 fixed-basis temporal-update 的 provisional working evolution。

---

# 21. Same-iteration hardening invariance

same frozen local state/search directions下：

$$ \breve\alpha_A=\breve\alpha_B, $$

$$ \breve{\dot\alpha}_A=\breve{\dot\alpha}_B, $$

$$ \breve\beta_A=\breve\beta_B. $$

以及：

$$ \breve{\bar r}_A=\breve{\bar r}_B, $$

$$ \breve{\dot{\bar r}}_A=\breve{\dot{\bar r}}_B, $$

$$ \breve{\bar R}_A=\breve{\bar R}_B. $$

---

# 22. Same-iteration damage invariance

tower v1 已冻结：

$$ \breve D=\hat D, $$

$$ \breve{\dot D}=\hat{\dot D}. $$

因此：

$$ \breve D_A=\breve D_B, $$

$$ \breve{\dot D}_A=\breve{\dot D}_B. $$

---

# 23. Same-iteration Y 仍是 trial-dependent

mechanical PGD branch不同，一般：

$$ \breve\sigma_A\neq\breve\sigma_B. $$

因此：

$$ \breve Y_A=Y(\breve\sigma_A,\hat D), $$

$$ \breve Y_B=Y(\breve\sigma_B,\hat D). $$

Trial B必须重新计算 $Y$。

---

# 24. Atomic LATIN iteration commit

一次正常 commit定义：

$$ (s_i,\mathcal B_m,\xi_i)\longrightarrow(s_{i+1},\mathcal B_{\rm next},\xi_{i+1}). $$

三者在同一 logical commit point一起替换。

禁止：

```text
basis committed
state rejected
```

或：

```text
state committed
candidate basis discarded
```

---

# 25. Trial A 绝对收敛优先

若：

$$ \xi_A\le\varepsilon_{\rm LATIN}, $$

则：

$$ s_{i+1}=s_A, $$

$$ \mathcal B_{\rm next}=\mathcal B_m^{A}, $$

$$ \xi_{i+1}=\xi_A. $$

即 atomic commit：

$$ (s_i,\mathcal B_m,\xi_i)\longrightarrow(s_A,\mathcal B_m^{A},\xi_A). $$

并结束 nonlinear solve。

优先级：

> absolute LATIN convergence 高于 saturation/enrichment decision。

这里不能 commit old $\mathcal B_m$，因为 $s_A$ 是由 updated temporal coordinates $\mathcal B_m^A$ 构造的。

---

# 26. Trial A existing-basis improvement充分

若：

$$ \zeta_A\gt\zeta_{\rm enrich}, $$

则 fixed spatial basis的 temporal update已经带来足够改善。

commit：

$$ (s_i,\mathcal B_m,\xi_i)\rightarrow(s_A,\mathcal B_m^{A},\xi_A). $$

进入下一 LATIN iteration。

因此这里的“fixed basis”仅表示：

> spatial basis / rank unchanged。

并不表示 temporal coordinates保持为 old persistent values。

---

# 27. Trial A enrichment request

以下情况表示 basis可能不足。

改善不足：

$$ \zeta_{\rm stop}\lt\zeta_A\le\zeta_{\rm enrich}. $$

或者：

$$ \zeta_A\lt0. $$

negative saturation不是 convergence saturation，而是 current reduced approximation worsened。

---

# 28. Very-small saturation 不直接作为 nonlinear convergence

若：

$$ 0\le\zeta_A\le\zeta_{\rm stop}, $$

tower v1 不直接宣布 solver converged。

仍需区分：

- LATIN nonlinear error是否已经小；
- reduced plastic residual是否已经小。

---

# 29. Reduced residual escape condition

设 Trial A fixed-basis reduced residual为：

$$ r_{\rm red,A}. $$

若 saturation提示 enrichment，但：

$$ r_{\rm red,A}\le\varepsilon_{\rm red}, $$

则说明 current reduced plastic global problem已经充分求解。

剩余 gap主要属于：

> nonlinear LATIN alternation。

因此不强制制造一个不存在的 residual-driven mode。

此时 atomic commit Trial A：

$$ (s_i,\mathcal B_m,\xi_i)\rightarrow(s_A,\mathcal B_m^{A},\xi_A). $$

然后进入下一 LATIN iteration。

---

# 30. 真正打开 enrichment transaction 的条件

只有当：

- Trial A未绝对收敛；
- saturation / inadequacy logic请求 enrichment；
- 且：

$$ r_{\rm red,A}\gt\varepsilon_{\rm red}, $$

才真正打开 mode transaction。

---

# 31. Mode rejection 的定义

可能原因：

```text
fixed-point not converged
spatial degeneracy
modified-time significance too small
basis-health failure
temporal solve failure
non-finite mode
full residual benefit insufficient
```

此时先在 mode-transaction 层：

$$ \mathcal B_{\rm work}\rightarrow\text{discard}, $$

并恢复到 provisional fixed-basis snapshot：

$$ \mathcal B_m^{A}. $$

但由于这里讨论的是：

> enrichment 已被真正要求，且 $r_{\rm red,A}\gt\varepsilon_{\rm red}$，

所以 tower v1 不把 $\mathcal B_m^A$ 或 Trial A persistent commit。

整个 LATIN iteration 随后失败，persistent rollback target仍是：

$$ \mathcal B_m. $$

---

# 32. Enrichment真正被要求但 mode失败

tower v1冻结：

> 不把 provisional Trial A 静默升级为 accepted persistent state。

persistent baseline保持：

$$ s_i,\mathcal B_m,\xi_i. $$

返回：

```text
ENRICHMENT_FAILED
```

同时可保留 Trial A作为 diagnostic。

---

# 33. 与 current 1D enrichment failure 的差异

current 1D 在 enrichment失败时会返回 current candidate state，同时标记 `ENRICHMENT_FAILED`。

tower v1 更严格：

> failed iteration 不改变 persistent accepted baseline。

这一点属于：

> transaction-safety refinement。

不是原论文显式要求。

---

# 34. Trial B 正常形成后的 commit

mode accepted后形成合法 complete Trial B。

若：

- $\breve s_B$ finite；
- relaxation finite；
- state invariants valid；
- $\xi_B$ finite；

则：

$$ (s_i,\mathcal B_m,\xi_i)\rightarrow(s_B,\mathcal B_{m+1}^{B*},\xi_B). $$

正式：

$$ \mathcal B_{m+1}=\mathcal B_{m+1}^{B*}. $$

注意：

> persistent commit 的 enlarged basis 已包含 all-mode temporal re-optimisation后的 coordinates。

---

# 35. Trial B 的 xi 不要求一定小于 Trial A

mode acceptance已经使用：

> full reduced mechanical residual benefit。

不额外要求：

$$ \xi_B\lt\xi_A. $$

因为：

- reduced residual；
- complete LATIN indicator；

衡量不同层级。

因此不会再增加 second mode gate。

---

# 36. Trial B 未收敛不代表 mode rollback

若：

$$ \xi_B\gt\varepsilon_{\rm LATIN}, $$

但 Trial B本身合法，则仍然 atomic commit Trial B。

之后进入下一 LATIN iteration。

这符合：

> one accepted pair per enrichment event。

---

# 37. Trial B hard failure

若完整 Trial B出现：

```text
NaN / Inf
invalid damage
invalid hardening recursion
equilibrium operator failure
invalid relaxation
indicator non-finite
```

则：

$$ s_B\rightarrow\text{discard}, $$

$$ \mathcal B_{m+1}^{B*}\rightarrow\text{discard}. $$

mode-transaction working state可回到：

$$ \mathcal B_m^{A}, $$

用于 diagnostics。

但整个 failed LATIN iteration 不 commit Trial A，因此 persistent baseline仍为：

$$ s_i,\mathcal B_m,\xi_i. $$

---

# 38. Four rollback concepts

## 38.1 Mode rejection

丢弃 candidate enlarged basis，并在 mode-transaction 层恢复：

$$ \mathcal B_m^{A}. $$

这是 **local mode rollback**。

## 38.2 Trial B failure

丢弃：

$$ s_B,\mathcal B_{m+1}^{B*}. $$

provisional $\mathcal B_m^A$ 可仅保留为 diagnostics。

## 38.3 Iteration failure

不产生新的 persistent $s_{i+1}$，也不 persistent commit $\mathcal B_m^A$。

persistent snapshot保持：

$$ s_i,\mathcal B_m,\xi_i. $$

这是 **LATIN iteration rollback**。

## 38.4 Normal solver commit

同时替换：

$$ s,\mathcal B,\xi. $$

Trial A commit必须使用：

$$ (s_A,\mathcal B_m^A,\xi_A). $$

Trial B commit必须使用：

$$ (s_B,\mathcal B_{m+1}^{B*},\xi_B). $$

因此未来不使用一个含义模糊的万能 `rollback()`。

---

# 39. Stagnation history 的 transaction semantics

nonlinear stagnation仅基于：

> committed LATIN iterations。

provisional Trial A不更新 stagnation counter。

rejected mode不更新。

failed Trial B不更新。

只有 atomic commit后：

$$ \xi_i\rightarrow\xi_{i+1} $$

才更新 stagnation history。

---

# 40. 完整 state-machine 优先级

```text
1. build Trial A

2. check absolute LATIN convergence

3. if not converged:
       evaluate saturation / basis adequacy

4. if enrichment suggested:
       inspect reduced residual

5. only if unresolved reduced residual exists:
       attempt one mode

6. if mode rejected:
       fail iteration transaction

7. if mode accepted:
       build complete Trial B

8. if Trial B valid:
       atomic commit B

9. update accepted convergence/stagnation history
```

---

# 41. Persistent objects、operators、ephemeral results 三层划分

## 41.1 Persistent value objects

```text
MaterialPointLayout
MaterialPointMetric
LatinStateTower
PGDBasisTower
```

## 41.2 Immutable operators / discretisation data

```text
TowerEquilibriumOperator
material parameters
search directions
time grid conventions
```

## 41.3 Ephemeral iteration results

```text
local snapshot
FrozenGlobalData
FixedBasisPGDResult
TowerEnrichmentResult
TowerGlobalCandidate
TowerTrialEvaluation
diagnostics
```

---

# 42. `MaterialPointLayout` 最终职责

只负责：

> topology / indexing。

概念字段：

```text
n_elements
n_material_points
element_index[q]
gauss_index[q]
fiber_index[q]
mapping q <-> (e,g,f)
```

ordering：

```text
element-major
Gauss-major
fiber-major
```

---

# 43. `MaterialPointLayout` 不负责什么

不存：

```text
LATIN state
PGD basis
search directions
damage
stress
integration metric
solver diagnostics
```

这样 layout可以作为整个 simulation 共享的 immutable topology descriptor。

---

# 44. Material-point metric

定义：

$$ v_q=J_{eg}w_gA_f. $$

conceptual：

$$ M=\operatorname{diag}(v_q). $$

代码不构造 dense $N_q\times N_q$ matrix。

存储：

$$ \mathbf v\in\mathbb R^{N_q}. $$

---

# 45. Metric operations

对于 spatial vectors：

$$ \mathbf x^TM\mathbf y=\sum_qv_qx_qy_q. $$

对于：

$$ M\mathbf x, $$

实现为：

$$ \mathbf v\odot\mathbf x. $$

这避免上万个 material points下的 dense metric storage。

---

# 46. Layout 与 metric 分离

冻结：

```text
MaterialPointLayout
    topology

MaterialPointMetric
    integration weights
```

它们可以都定义在 tower基础模块中，但概念职责不混淆。

---

# 47. `LatinStateTower`

建议文件：

```text
latin/tower_state.py
```

拥有：

```text
time

plastic_strain_rate
elastic_strain
alpha_rate
r_bar_rate
damage_rate
stress
beta
R_bar
energy_release_rate

plastic_strain
alpha
r_bar
damage
```

---

# 48. State shape contract

$$ \operatorname{shape}(\text{time})=(N_t,). $$

其余：

$$ \operatorname{shape}(\text{field})=(N_t,N_q). $$

所有 fields：

- finite；
- common time grid；
- common $N_q$；
- damage admissible。

---

# 49. `LatinStateTower` 不存 topology

不存：

```text
element id
Gauss id
fiber id
fiber area
fiber y-coordinate
quadrature weight
```

这些由 layout / section discretisation外部持有。

---

# 50. `LatinStateTower` 不存 PGD data

不存：

```text
spatial modes
temporal functions
PGD residual
mode diagnostics
basis size
```

PGD representation与 LATIN material state严格分离。

---

# 51. State public immutability

public state API不提供 in-place mutator。

操作原则：

```text
input state
    ↓
return new state
```

第一版可以用 deep copy确保 transaction安全。

---

# 52. Writable aliasing 禁止

若 future性能优化共享 arrays，必须保证：

> shared array不可通过任何 public path修改。

否则 accepted baseline可能被 trial污染。

---

# 53. `PGDModeTower`

建议继续定义在：

```text
latin/pgd_basis.py
```

而不立即新建 tower-specific basis文件。

字段：

$$ p_j\in\mathbb R^{N_q}, $$

$$ s_j\in\mathbb R^{N_q}, $$

$$ \lambda_j\in\mathbb R^{N_t}, $$

$$ \dot\lambda_j\in\mathbb R^{N_t}. $$

以及：

```text
iteration_added
```

---

# 54. `PGDBasisTower`

是 mode collection。

但 persistent semantics采用：

> value-style basis。

不允许 external caller直接：

```text
basis.append(...)
basis.modes[j].temporal_amplitude[:] = ...
```

---

# 55. PGD reconstruction

设：

$$ P=[p_1,\ldots,p_m]\in\mathbb R^{N_q\times m}, $$

$$ \Lambda=[\lambda_1,\ldots,\lambda_m]\in\mathbb R^{N_t\times m}. $$

则：

$$ \Delta E^p=\Lambda P^T. $$

shape：

$$ (N_t,N_q). $$

---

# 56. Plastic-rate reconstruction

设：

$$ \dot\Lambda=[\dot\lambda_1,\ldots,\dot\lambda_m]. $$

则：

$$ \Delta\dot E^p=\dot\Lambda P^T. $$

---

# 57. Stress correction reconstruction

设：

$$ S=[s_1,\ldots,s_m]. $$

则：

$$ \Delta\Sigma'=\Lambda S^T. $$

---

# 58. `tower_equilibrium_operator.py`

唯一职责：

> reference structural projection。

理论 operator：

$$ \mathcal E_{\rm tower}=H(H^TMC_0H)^{-1}H^TMC_0. $$

但 public API不显式暴露 dense $\mathcal E_{\rm tower}$。

---

# 59. `TowerEquilibriumOperator`

初始化时绑定：

```text
tower discretisation
MaterialPointLayout
MaterialPointMetric
reference elastic stiffness C0
compatibility operator H
boundary conditions
free structural DOFs
```

构造后视为 immutable operator。

---

# 60. Spatial projection API

输入：

$$ r\in\mathbb R^{N_q}. $$

返回：

```text
compatible_strain_q
stress_q
displacement_free
```

conceptual：

```text
apply_spatial(source_q)
```

---

# 61. History projection API

输入：

$$ R\in\mathbb R^{N_t\times N_q}. $$

返回：

$$ E_{\rm comp}\in\mathbb R^{N_t\times N_q}, $$

$$ \Sigma\in\mathbb R^{N_t\times N_q}, $$

$$ U\in\mathbb R^{N_t\times N_{\rm dof,free}}. $$

conceptual：

```text
apply_history(source_tq)
```

---

# 62. `EquilibriumProjectionTower`

ephemeral result：

```text
compatible_strain
stress
displacement
```

它不包含：

```text
LatinStateTower
hardening
damage
PGD basis
xi
zeta
```

---

# 63. `tower_pgd_time_update.py`

唯一职责：

> fixed spatial basis 下 solve Eq. (58)–(59) temporal update。

不生成新 spatial mode。

---

# 64. Fixed-basis time-update inputs

```text
basis
time
full forcing f
H_sigma
MaterialPointMetric
TowerEquilibriumOperator
time discretisation convention
rcond / temporal solver controls
```

---

# 65. `FixedBasisPGDResult`

建议包含：

```text
basis
plastic_strain_correction
plastic_strain_rate_correction
plastic_projection
mechanical_residual
relative_residual
forcing_norm
reduced_converged
diagnostics
```

其中 result中的 basis：

> spatial basis不变，只更新 temporal coordinates。

---

# 66. Time-update module 禁止事项

不能：

```text
local constitutive integration
new-mode generation
Gram-Schmidt
mode append
LatinStateTower construction
relaxation
xi / zeta
commit
```

---

# 67. Empty basis semantics

若：

$$ m=0, $$

不报 numerical exception。

返回：

$$ \Delta\varepsilon^p=0, $$

$$ \Delta\dot\varepsilon^p=0, $$

$$ \Delta\sigma'=0. $$

机械 residual：

$$ R=-f. $$

这样 solver可以根据真实 residual决定是否 bootstrap first mode。

---

# 68. First-mode bootstrap

若 empty basis且：

$$ r_{\rm red}\gt\varepsilon_{\rm red}, $$

则 first enrichment是必需的。

不需要把 empty basis写成特殊 error。

---

# 69. `tower_pgd_enrichment.py`

唯一职责：

> one-mode transaction。

不拥有 persistent basis。

---

# 70. Enrichment inputs

```text
provisional fixed-basis B_m^A
current FixedBasisPGDResult_A
full forcing f
current shifted defect / residual R_A
H_sigma
MaterialPointMetric
TowerEquilibriumOperator
fixed-point tolerance
mode significance controls
acceptance tolerance
iteration id
```

这里不把 old persistent $\mathcal B_m$ 作为 enrichment working basis。

原因：

> current residual $R_A$ 本身对应已经完成 Eq. (58)–(59) temporal update 的 $\mathcal B_m^A$。

persistent $\mathcal B_m$ 由 outer solver保留，仅作为 whole-iteration rollback target。

---

# 71. Enrichment internal sequence

```text
raw rank-one fixed point
    ↓
weighted Modified Gram-Schmidt
    ↓
exact temporal coordinate transformation
    ↓
field-invariance test
    ↓
gamma_sp / gamma_lambda
    ↓
tentative enlarged basis
    ↓
all-mode temporal re-optimisation
    ↓
full residual reconstruction
    ↓
Delta_res acceptance
```

---

# 72. `TowerEnrichmentResult`

建议包含：

```text
accepted
failure_reason
candidate_fixed_basis_result
raw_mode_diagnostics
orthogonality_diagnostics
temporal_significance
residual_benefit
fixed_point_history
```

---

# 73. Enrichment rejection result

若 rejected：

```text
accepted = False
candidate_fixed_basis_result = None
```

输入：

$$ \mathcal B_m^{A} $$

必须完全不变。

也就是说 one-mode transaction可以失败，但不得污染 Trial A 的 provisional fixed-basis result。

---

# 74. Enrichment accepted result

若 accepted：

`candidate_fixed_basis_result.basis` 为：

$$ \mathcal B_{m+1}^{B*}. $$

它由 $\mathcal B_m^A$ enriched 并完成 enlarged all-mode temporal re-optimisation得到。

仍未 persistent commit。

---

# 75. Enrichment module 不读取 hardening/damage state

它只处理 reduced plastic mechanical problem。

因此不能：

```text
update alpha
update D
compute Y
relax state
compute xi
compute zeta
```

---

# 76. `tower_global_finishing.py`

负责：

1. prepare frozen global data；
2. build complete unrelaxed global candidate。

二者都采用 pure-function semantics。

---

# 77. `prepare_frozen_global_data(...)`

输入：

```text
s_i
hat{s}
search directions
materials
TowerEquilibriumOperator
MaterialPointLayout / Metric as needed
```

输出：

```text
FrozenGlobalData
```

---

# 78. `FrozenGlobalData`

至少包含：

```text
damage_residual_strain
damage_projection
full_plastic_forcing
damage_displacement_correction
```

必要时包含与 forcing reconstruction有关的已计算 quantities。

same iteration Trial A/B间 immutable。

---

# 79. 为什么 forcing 与 damage projection 应提前计算一次

same iteration中：

- baseline same；
- local same；
- search directions same；
- damage histories same。

因此 Trial A/B不需要重复求相同 full-order damage branch。

这同时保证严格 common-input comparison。

---

# 80. `build_unrelaxed_candidate(...)`

输入：

```text
baseline state s_i
local state hat{s}
search directions
FrozenGlobalData
FixedBasisPGDResult
materials
```

输出：

```text
TowerGlobalCandidate
```

---

# 81. Mechanical candidate assembly

$$ \breve\varepsilon^p=\varepsilon_i^p+\Delta\varepsilon^p. $$

$$ \breve{\dot\varepsilon}^p=\dot\varepsilon_i^p+\Delta\dot\varepsilon^p. $$

$$ \breve\sigma=\sigma_i+\Delta\sigma'+\Delta\tilde\sigma. $$

$$ \breve\varepsilon^e=\varepsilon_i^e+\Delta\varepsilon'-\Delta\varepsilon^p+\Delta\tilde\varepsilon. $$

---

# 82. Hardening candidate assembly

执行 Eq. (73)–(74) tower v1 BE finishing。

得到：

```text
breve alpha
breve alpha_rate
breve beta
breve r_bar
breve r_bar_rate
breve R_bar
```

---

# 83. Damage candidate assembly

tower v1：

$$ \breve D=\hat D. $$

$$ \breve{\dot D}=\hat{\dot D}. $$

不执行第二次 BE damage reintegration。

---

# 84. Energy-release candidate assembly

final stress完成后：

$$ \breve Y=Y(\breve\sigma,\breve D). $$

candidate上 nonlinear $Y$ relation精确成立。

---

# 85. `TowerGlobalCandidate`

建议包含：

```text
state
total_displacement_correction
plastic_displacement_correction
damage_displacement_correction
diagnostics
```

其中：

```text
state
```

就是：

$$ \breve s_{i+1}. $$

---

# 86. Global-finishing module 禁止事项

不能：

```text
modify basis
generate mode
Gram-Schmidt
relax state
compute xi
compute zeta
commit
```

---

# 87. `tower_iteration_control.py`

职责：

```text
relaxation
Eq. (76)
Eq. (77)
Eq. (60)
convergence/saturation data evaluation
```

不拥有 basis。

---

# 88. Relaxation API

输入：

$$ s_i,\breve s. $$

输出：

$$ s_{\rm trial}=(1-\mu)s_i+\mu\breve s. $$

所有 primary/support histories使用同一 $\mu$。

---

# 89. Relaxation 后不 reproject Y

禁止：

$$ Y_{\rm trial}\leftarrow Y(\sigma_{\rm trial},D_{\rm trial}). $$

relaxed state允许：

$$ Y_{\rm trial}\neq Y(\sigma_{\rm trial},D_{\rm trial}). $$

---

# 90. Eq. (77) reader fields

只读取：

$$ \sigma, $$

$$ \beta, $$

$$ \bar R, $$

$$ \dot\varepsilon^p, $$

$$ \varepsilon^e, $$

$$ \dot\alpha, $$

$$ \dot{\bar r}. $$

---

# 91. Tower Eq. (77) spatial integration

1D element-volume integration替换为 material-point metric：

$$ \sum_qv_q(\cdot). $$

time integration使用 tower v1统一的离散 convention。

---

# 92. `TowerTrialEvaluation`

建议包含：

```text
unrelaxed_state
relaxed_state
indicator
saturation
converged
finite
diagnostics
```

但：

> 不包含 basis commit semantics。

---

# 93. `tower_latin_pgd_solver.py`

唯一职责：

> transaction orchestration。

也是唯一 persistent owner。

---

# 94. Solver persistent variables

只由 solver持有并更新：

$$ s_i, $$

$$ \mathcal B_m, $$

$$ \xi_i. $$

---

# 95. Solver normal orchestration

```text
accepted baseline
    s_i, B_m, xi_i
        ↓
local stage
        ↓
search directions
        ↓
prepare FrozenGlobalData
        ↓
fixed-basis time update
        ↓
provisional B_m^A
        ↓
build Trial A candidate from s_i + B_m^A
        ↓
relax/evaluate Trial A
        ↓
decision
```

若 Trial A被接受，则 atomic commit：

```text
(s_i, B_m, xi_i)
        ↓
(s_A, B_m^A, xi_A)
```

---

# 96. Solver enrichment orchestration

```text
enrichment required
and r_red,A > eps_red
        ↓
open one-mode transaction on B_m^A
        ↓
TowerEnrichmentResult
    ├── rejected
    │      ↓
    │   mode rollback to B_m^A
    │      ↓
    │   fail current LATIN iteration
    │      ↓
    │   persistent remains (s_i, B_m, xi_i)
    │
    └── accepted
           ↓
       candidate B_(m+1)^(B*)
           ↓
       build Trial B again from SAME s_i
           ↓
       relax/evaluate Trial B
           ↓
       atomic commit
           ↓
       (s_B, B_(m+1)^(B*), xi_B)
```

---

# 97. Solver 唯一允许改变 persistent references

正常 commit block：

```text
accepted_state     = trial.relaxed_state
accepted_basis     = trial_basis
accepted_indicator = trial.indicator
```

三者必须属于同一 trial。

---

# 98. Local stage 与 search directions 的 first-code strategy

第一批不额外创建：

```text
tower_local_stage.py
tower_search_directions.py
tower_pgd_basis.py
```

优先复用 / 扩展：

```text
latin/local_stage.py
latin/search_directions.py
latin/pgd_basis.py
```

理由：

> 现有 architecture boundary已经验证；先避免重复实现。

后续若 1D/tower分支逻辑明显失控，再拆 module。

---

# 99. 正式新增的 seven-module set

```text
latin/tower_state.py

latin/tower_equilibrium_operator.py

latin/tower_pgd_time_update.py

latin/tower_pgd_enrichment.py

latin/tower_global_finishing.py

latin/tower_iteration_control.py

latin/tower_latin_pgd_solver.py
```

---

# 100. Module sole responsibility 表

| Module                          | Sole responsibility                                      |
|---------------------------------|----------------------------------------------------------|
| `tower_state.py`                | material-point layout + state value representation       |
| `tower_equilibrium_operator.py` | metric + reference structural projection                |
| `tower_pgd_time_update.py`      | fixed-basis temporal solve                               |
| `tower_pgd_enrichment.py`       | one-mode transaction                                     |
| `tower_global_finishing.py`     | global preparation + complete unrelaxed candidate        |
| `tower_iteration_control.py`    | relaxation + Eq. (76)–(77) + saturation evaluation      |
| `tower_latin_pgd_solver.py`     | transaction orchestration + commit / failure handling    |

---

# 101. Dependency direction

推荐：

```text
tower_state
    ↑
    │
pgd_basis       tower_equilibrium_operator
    ↑                    ↑
    │                    │
tower_pgd_time_update ───┘
    ↑
    │
tower_pgd_enrichment

tower_state ──────────────→ tower_global_finishing
tower_equilibrium_operator ─→ tower_global_finishing

tower_state ──────────────→ tower_iteration_control

all lower layers
        ↑
        │
tower_latin_pgd_solver
```

---

# 102. Circular dependency prohibition

底层 module绝不 import：

```text
tower_latin_pgd_solver.py
```

solver只能向下依赖。

不能让：

```text
global_finishing
    imports solver

iteration_control
    imports solver
```

---

# 103. Conceptual data containers 最终清单

```text
MaterialPointLayout
MaterialPointMetric
LatinStateTower
PGDModeTower
PGDBasisTower
EquilibriumProjectionTower
FrozenGlobalData
FixedBasisPGDResult
TowerEnrichmentResult
TowerGlobalCandidate
TowerTrialEvaluation
TowerLatinPGDResult
```

---

# 104. 不采用 mutable global context

不设计：

```text
TowerLatinContext
```

把：

```text
state
basis
layout
operator
directions
local state
residual
xi
zeta
```

全部塞进一个可修改对象。

原因：

> ownership会重新模糊。

---

# 105. Failure diagnostics 与 accepted state 分离

若 iteration失败，final solver result仍可保存：

```text
last_accepted_state
last_accepted_basis
last_accepted_indicator

last_trial_A
last_enrichment_result
failure_reason
last_reduced_residual
```

但：

> failed trial不得伪装成 accepted state。

---

# 106. Unit-test boundary 1：state/layout

检查：

```text
q mapping
shape validation
time validation
damage bounds
copy/value semantics
no accepted-state mutation
layout/state separation
```

---

# 107. Unit-test boundary 2：metric/equilibrium operator

检查：

$$ H^TM\Delta\sigma\approx0. $$

检查 compatible projection。

检查 spatial/history API一致性。

---

# 108. Unit-test boundary 3：fixed-basis time update

检查：

```text
empty basis
single mode
multi-mode
BE temporal consistency
residual reconstruction
no basis-size change
input basis immutability
```

---

# 109. Unit-test boundary 4：one-mode enrichment

检查：

```text
weighted orthogonality
exact field invariance
gamma_sp
gamma_lambda
all-mode re-optimisation
Delta_res
complete rollback on rejection
input basis immutability
one accepted mode only
```

---

# 110. Unit-test boundary 5：global finishing

检查：

```text
mechanical assembly
equilibrium
compatibility
Eq. (73)
Eq. (74)
direct D inheritance
direct D_rate inheritance
no second BE damage integration
candidate Y state law
Trial A/B hardening invariance
Trial A/B damage invariance
Trial-dependent Y
```

---

# 111. Unit-test boundary 6：iteration control

检查：

```text
field-wise relaxation
Eq. (77) seven-field reader
material-point metric integration
same mu on all stored histories
no relaxed-Y reprojection
xi finite
zeta uses same baseline xi_i
```

---

# 112. Unit-test boundary 7：transaction semantics

检查：

```text
s_i immutable
persistent B_m immutable
xi_i frozen

fixed-basis update returns B_m^A
B_m^A has same spatial rank/modes as B_m
B_m^A may have different temporal coordinates

Trial A and Trial B use same state baseline s_i
same local snapshot
same search directions
same FrozenGlobalData

enrichment works on B_m^A
mode rejection leaves B_m^A unmodified
mode rejection does not mutate persistent B_m

Trial B failure leaves persistent baseline untouched

commit A updates state + temporal-updated basis + indicator atomically
    -> (s_A, B_m^A, xi_A)

commit B updates state + enriched basis + indicator atomically
    -> (s_B, B_(m+1)^(B*), xi_B)

iteration failure restores persistent
    -> (s_i, B_m, xi_i)

stagnation updates only after commit
```

---

# 113. Minimal integration test

首个 end-to-end test不直接上完整 40 × 4 × 64 tower。

应采用 coarse tower：

```text
10 elements
2 Gauss points
16 fibers
```

先检查：

- storage；
- operator；
- local/global alternation；
- one pair；
- transaction；
- no mutation；
- finite outputs。

---

# 114. 正式 tower mesh 仍保持后续目标

candidate formal discretisation：

```text
40 elements
4 Gauss points
64 fibers
```

对应：

$$ N_q=40\times4\times64=10240. $$

因此本阶段关于：

- flattened `(Nt,Nq)`；
- vector metric；
- no dense `M`；

的选择具有直接必要性。

---

# 115. 当前 1D architecture 中应继承的部分

继承：

```text
state and basis separate
basis reconstruction
deep-copy safety direction
fixed-basis time update
residual-driven enrichment
joint temporal re-optimisation
field-wise relaxation
Eq. (76) / Eq. (77)
accepted-iteration convergence history
```

---

# 116. 当前 1D architecture 中不直接复制的部分

不直接复制：

```text
one global-stage function owns too many responsibilities
mutable working basis exposed too broadly
multiple accepted modes in one iteration
damage rate copy + BE damage reintegration
trapezoidal metric everywhere
failed enrichment returns provisional state as if it were accepted baseline
```

tower v1会在接口层拆开这些职责。

---

# 117. One-pair-per-event 再确认

tower v1：

$$ \text{max accepted new pairs per enrichment event}=1. $$

accepted pair后：

```text
complete Trial B
    ↓
commit
    ↓
next LATIN iteration
```

不在 same iteration继续生成第二个 mode。

---

# 118. 为什么 one-pair rule 有利于 transaction traceability

它保持：

```text
one saturation decision
    ↓
one pair attempt
    ↓
one complete enriched trial
    ↓
one commit
```

这样：

- mode diagnostics；
- $\xi$；
- $\zeta$；
- basis size；

之间关系清楚。

---

# 119. 正式进入 code stage 的门槛

截至本文档：

```text
canonical material-point storage
    frozen

14-entry state contract
    frozen

mechanical field contract
    frozen

hardening field contract
    frozen

damage-energy contract
    frozen

state ownership
    frozen

Trial A / Trial B semantics
    frozen

mode transaction
    frozen

atomic commit / rollback
    frozen

module I/O
    frozen

dependency direction
    frozen

unit-test boundaries
    frozen
```

因此：

> 代码前 specification 已完成。

---

# 120. 下一步不再继续扩展理论

下一步直接：

```text
create:
    latin/tower_state.py

create tests:
    state/layout focused unit tests
```

并采用：

> one implementation step at a time。

先确保：

- `MaterialPointLayout`；
- `MaterialPointMetric` 的 ownership边界；
- `LatinStateTower`；
- `PGDModeTower / PGDBasisTower` compatibility plan；

清楚并测试通过。

---

# 121. 第一批 code 不应包含的内容

在 `tower_state.py` 第一批实现中，不写：

```text
equilibrium operator
PGD time solver
enrichment
global finishing
relaxation
solver loop
```

先只实现 state/layout contract。

---

# 122. 第一批 tests 的最低目标

至少验证：

```text
layout flatten / unflatten
element-major ordering
state shape consistency
strictly increasing time
finite fields
damage bounds
deep-copy/value behavior
state does not own topology
state construction does not alias writable input arrays
```

---

# 123. Git workflow 继续保持标准流程

每个可独立验证的 code stage：

```text
git status --short
git add <file paths>
git status --short
git commit -m "<message>"
git status
git push origin feature/offshore-wind-turbine-tower-fatigue
git status
```

不要把多个理论/实现阶段混成一个大 commit。

---

# 124. Markdown/PyCharm Preview 规范继续继承

后续所有阶段总结：

- inline math：`$...$`
- display math：单物理行 `$$ ... $$`
- raw `<` 禁止进入数学环境
- less-than 使用 `\lt`
- greater-than 优先 `\gt`
- 长流程使用 fenced code block
- 不把复杂中文嵌入 LaTeX `array`
- preview 异常时优先最小复现与 binary split
- 一次只改变一个变量

---

# 125. 本阶段最终研究停点

截至 2026-08-18：

> Offshore-wind-turbine tower LATIN-PGD 的 code-before theory/data/transaction/module-interface specification 已完成 final freeze。

current exact stop：

```text
NEXT
    implement latin/tower_state.py

THEN
    add state/layout unit tests

ONLY AFTER PASS
    proceed to tower_equilibrium_operator.py
```

这一停点意味着：

> 后续不再通过继续增加理论层来推迟实现，而是进入受 unit tests 驱动的增量代码阶段。
