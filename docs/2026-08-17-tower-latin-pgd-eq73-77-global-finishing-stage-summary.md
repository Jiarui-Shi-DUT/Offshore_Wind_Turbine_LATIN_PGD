<!-- PyCharm-preview-safe regenerated copy: 2026-08-17 -->

# Tower LATIN-PGD global-stage finishing、damage inheritance、relaxation 与 outer control 阶段总结

**日期：2026-08-17**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**当前研究路线：Bhattacharyya et al. 原论文 $x-t$ LATIN-PGD → 2D fiber beam-column offshore wind turbine tower**  
**阶段范围：在 post-fixed-point Gram–Schmidt、temporal correction、new-mode significance、all-mode temporal re-optimisation 与 mode acceptance 已闭合的基础上，继续闭合 accepted new pair 之后 global stage 的 finishing sequence，包括 mechanical candidate assembly、Eq. (73)–(74) hardening variables、Eq. (75) damage treatment、damage-history inheritance、final energy-release-rate update、unrelaxed/relaxed state separation、support histories relaxation、Eq. (76) LATIN indicator、Eq. (60) saturation control，以及 one-pair-per-enrichment-event 下完整 outer iteration 的数据依赖。**  
**上一阶段衔接：`2026-08-17-tower-latin-pgd-post-fixed-point-gram-schmidt-mode-acceptance-stage-summary.md`**  
**下一阶段：对 Eq. (77) mechanical norm 与 future tower `LatinState` data structure 做 field-by-field contract，明确 primary LATIN fields、integrated support histories、derived fields 与 exact consistency checks；随后冻结 tower LATIN-PGD solver 模块接口，再进入代码实现。**

---

# 1. 本阶段定位

上一阶段已经将一个 new PGD pair 从 fixed-point converged raw candidate 推进到 persistent accepted basis。

完整 Add-a-pair branch 已经可以写成：

```text
Eq. (60): existing basis insufficient
    ↓
Eq. (61)–(64): shifted defect + rank-one ansatz
    ↓
Eq. (65)–(71): fixed temporal → spatial solve
    ↓
Eq. (72): fixed spatial → temporal minimisation
    ↓
Eq. (70)–(72): alternating fixed point
    ↓
complete-pair fixed-point convergence
    ↓
M-weighted Gram-Schmidt
    ↓
exact temporal coordinate transformation
    ↓
modified-new-time significance
    ↓
all-mode Eq. (58)–(59) temporal re-optimisation
    ↓
full residual benefit check
    ↓
accept persistent mode or rollback
```

因此本阶段的起点是：

$$ \boxed{\text{accepted PGD basis } \mathcal B_{m+1}} $$

以及该 enlarged basis 对 current global-stage plastic correction 的最终 reconstruction。

本阶段不再讨论 new mode 如何产生，而是回答：

> accepted new pair 之后，如何完成整个 LATIN global stage，并把 unrelaxed global candidate 变成真正参与 Eq. (76)–(77) convergence evaluation 的 relaxed global state？

---

# 2. 本阶段要闭合的主要问题

本阶段集中处理以下紧密耦合的问题。

第一，accepted PGD plastic correction 与 damage-dependent structural correction 如何组装为完整 unrelaxed mechanical candidate？

第二，Eq. (73) 中 normal hardening state law 在当前 tower material-point model 中具体是什么？

第三，Eq. (74) 如何离散成每个 fiber material point 上的 local linear ODE？

第四，hardening variables 在 $t_0$ 如何处理？

第五，Eq. (75) 的 $b^-=0$ 到底意味着只复制 local damage rate，还是 integrated damage history 也应直接继承？

第六，current 1D implementation 采用 local RK4 damage history + global backward-Euler re-integration，这与原论文连续理论之间有什么差异？

第七，tower v1 是否应该直接令 unrelaxed global damage history 等于 local integrated damage history？

第八，final energy-release rate $Y$ 应在什么时刻计算？

第九，unrelaxed candidate $\breve s_{i+1}$ 与 relaxed state $s_{i+1}$ 必须如何区分？

第十，原论文 formal state 中未显式存储的 integrated support histories 是否也应该 relaxation？

第十一，Eq. (76) 应使用 relaxation 前还是 relaxation 后的 global state？

第十二，Eq. (60) saturation 应如何与 fixed-basis trial、one-pair enrichment trial 和下一 LATIN iteration 连接？

第十三，one-pair-per-enrichment-event 在完整 outer loop 中的准确含义是什么？

---

# 3. 本阶段资料依据

本阶段主要依据以下资料。

原论文理论基准：

- Bhattacharyya et al. LATIN-PGD cycling damage paper；
- Eq. (47)–(57)：global mechanical correction decomposition；
- Eq. (58)–(72)：PGD update / enrichment；
- Eq. (73)–(75)：hardening / damage global-stage finishing；
- Eq. (76)–(77)：LATIN convergence indicator 与 norm；
- Fig. 2：overall LATIN-PGD algorithm flow；
- relaxation parameter $\mu=0.8$。

当前 repository 参考代码：

- `latin/global_stage.py`
- `latin/pgd_global_stage.py`
- `latin/local_stage.py`
- `latin/iteration_control.py`
- `latin/pgd_solver.py`
- `latin/state.py`
- `latin/search_directions.py`

当前已保存阶段文档：

- `2026-08-15-tower-latin-pgd-global-stage-temporal-damage-hardening-stage-summary.md`
- `2026-08-16-tower-latin-pgd-eq60-saturation-and-latin-norm-stage-summary.md`
- `2026-08-17-tower-latin-pgd-post-fixed-point-gram-schmidt-mode-acceptance-stage-summary.md`

---

# 4. 本阶段结论类别

继续严格区分四类来源。

本阶段继续严格区分以下四类来源。

**原论文明确内容：** Eq. (73) 给出 linear normal hardening state law；Eq. (74) 给出 hardening descent relation；Eq. (75) 在 $b^-=0$ 时给出 damage-rate inheritance；最终 energy release rate 由 Eq. (4d) 计算；global candidate 完成后应用 relaxation；Eq. (76) 比较 local half-step state 与 relaxed global state。

**由原论文连续关系推导：** 在 common initial damage 下，如果 global 与 local damage rate 一致，则 integrated damage histories 也一致；linear support-history relaxation 与 rate relaxation 在连续积分关系下相容。

**current 1D implementation：** local stage 使用 RK4 integration；global hardening 使用 backward Euler；global damage 复制 local damage rate 后重新用 backward Euler 积分；`relax_global_state()` 对 primary fields 与 integrated histories 全部做 convex blend；Eq. (76) 在 relaxation 后计算。

**tower v1 engineering choice：** unrelaxed global damage history 直接继承 local integrated damage；保留 local damage rate；不再对 local damage rate 做第二次 backward-Euler re-integration；support histories 与 primary state 使用相同 relaxation；relaxation 后不再强制将 energy-release rate重新投影到 nonlinear state law；一个 LATIN enrichment event 最多接受一个 new pair。

---

# 5. accepted new pair 后的准确数据状态

设 current LATIN iteration 为：

$$ i $$

上一 accepted relaxed global state 为：

$$ s_i $$

local stage 已得到：

$$ \hat s_{i+1/2} $$

search directions 已固定：

$$ H_\sigma,\quad H_\beta,\quad H_{\bar R} $$

damage-related structural source 已固定。

PGD branch 已完成 current enlarged basis：

$$ \mathcal B_{m+1} $$

并经 all-mode temporal re-optimisation 得到 final plastic correction。

因此现在已知：

$$ \Delta\varepsilon^p_{i+1}(t,q) $$

以及：

$$ \Delta\dot\varepsilon^p_{i+1}(t,q) $$

其中：

$$ q=(e,g,f) $$

表示 tower 中：

- beam element $e$；
- Gauss point $g$；
- section fiber $f$。

---

# 6. tower reference equilibrium operator 回顾

前序阶段已经建立 material-point source strain 到 globally compatible strain / equilibrated stress 的 reference operator。

定义：

$$ \mathcal E_{\rm tower}=H(H^TMC_0H)^{-1}H^TMC_0 $$

其中：

- $H$ 为 tower fiber kinematic operator；
- $M$ 为 material-point integration metric；
- $C_0=E_0I$ 为 reference elasticity。

对于任意 source strain：

$$ p $$

compatible strain 为：

$$ \varepsilon_{\rm comp}=\mathcal E_{\rm tower}p $$

equilibrated stress correction 为：

$$ \sigma_{\rm eq}=C_0(\mathcal E_{\rm tower}-I)p $$

并满足：

$$ H^TM\sigma_{\rm eq}=0 $$

因此 plastic PGD 与 damage residual branch 可以共享同一个 reference equilibrium operator。

---

# 7. plastic branch 的 final mechanical correction

accepted PGD basis reconstructs：

$$ \Delta\varepsilon^p_{i+1} $$

PGD temporal rates reconstruct：

$$ \Delta\dot\varepsilon^p_{i+1} $$

由 reference equilibrium operator 得到：

$$ \Delta\varepsilon'_{i+1}=\mathcal E_{\rm tower}\Delta\varepsilon^p_{i+1} $$

以及：

$$ \Delta\sigma'_{i+1}=C_0(\mathcal E_{\rm tower}-I)\Delta\varepsilon^p_{i+1} $$

这里 prime branch 表示由 plastic source strain 引起的 globally admissible mechanical correction。

---

# 8. damage-dependent structural branch

在原论文中，由 damage nonlinear elastic law 引起的 residual correction 不进入 PGD reduced plastic branch。

它通过 full reference-equilibrium projection 处理。

定义 damage residual strain source：

$$ \Delta\varepsilon^R_{i+1} $$

经过：

$$ \Delta\tilde\varepsilon_{i+1}=\mathcal E_{\rm tower}\Delta\varepsilon^R_{i+1} $$

以及：

$$ \Delta\tilde\sigma_{i+1}=C_0(\mathcal E_{\rm tower}-I)\Delta\varepsilon^R_{i+1} $$

得到 damage-related compatible strain 与 equilibrated stress correction。

因此 global mechanical correction 仍然保持：

$$ \boxed{\text{plastic PGD reduced branch}+\text{damage full projection branch}} $$

---

# 9. 为什么 damage branch 仍不能进入 PGD basis

即使当前 damage history 被直接继承 local stage，damage-dependent structural residual correction仍然必须存在。

原因是两个不同概念：

$$ D(t,q) $$

是 constitutive internal history。

而：

$$ \Delta\varepsilon^R\rightarrow\Delta\tilde\varepsilon,\Delta\tilde\sigma $$

是 nonlinear damaged elastic state law 相对于 reference elasticity 的 structural admissibility correction。

所以：

$$ \boxed{D_{\rm global}=\hat D} $$

并不意味着：

$$ \boxed{\Delta\tilde\sigma=0} $$

两者必须严格区分。

---

# 10. unrelaxed global candidate 的符号

本阶段统一定义：

$$ \breve s_{i+1} $$

为 current global stage 完成后、relaxation 前的 candidate。

真正用于下一 LATIN iteration 的 global state 为：

$$ s_{i+1} $$

并满足：

$$ s_{i+1}=\mu\breve s_{i+1}+(1-\mu)s_i $$

因此未来代码与文档中必须避免把：

```text
candidate_state
```

与：

```text
accepted relaxed global_state
```

混为一谈。

---

# 11. unrelaxed accumulated plastic strain

plastic accumulated history 更新为：

$$ \breve\varepsilon^p_{i+1}=\varepsilon^p_i+\Delta\varepsilon^p_{i+1} $$

这是 current PGD spatial-temporal correction 对 previous global accumulated plastic history 的增量更新。

注意：

$$ \varepsilon^p $$

不是原论文 Eq. (76) primary norm 中直接出现的 field，但数值实现需要显式保存，因为：

- local stage history integration 需要它；
- plastic-rate consistency test 需要它；
- material update 需要 initial condition。

---

# 12. unrelaxed plastic strain rate

rate field 更新为：

$$ \breve{\dot\varepsilon}^p_{i+1}=\dot\varepsilon^p_i+\Delta\dot\varepsilon^p_{i+1} $$

该 field 是 formal LATIN state 的 primary component。

因此它会直接进入 Eq. (77) mechanical norm。

---

# 13. unrelaxed stress assembly

stress correction由 plastic branch 与 damage branch 相加：

$$ \breve\sigma_{i+1}=\sigma_i+\Delta\sigma'_{i+1}+\Delta\tilde\sigma_{i+1} $$

这也是 current 1D `global_stage.py` 与 `pgd_global_stage.py` 的基本 assembly structure。

因此 tower v1 不需要重新定义 stress branch，只需要把 1D bar equilibrium operator 换成 tower reference equilibrium operator。

---

# 14. unrelaxed elastic strain assembly

total compatible strain correction 中 plastic source strain需要被扣除，因此：

$$ \breve\varepsilon^e_{i+1}=\varepsilon^e_i+\Delta\varepsilon'_{i+1}-\Delta\varepsilon^p_{i+1}+\Delta\tilde\varepsilon_{i+1} $$

这条式子保证：

$$ \Delta\varepsilon^{\rm total}=\Delta\varepsilon'+\Delta\tilde\varepsilon $$

而：

$$ \Delta\varepsilon^{\rm total}=\Delta\varepsilon^e+\Delta\varepsilon^p $$

所以：

$$ \Delta\varepsilon^e=\Delta\varepsilon'-\Delta\varepsilon^p+\Delta\tilde\varepsilon $$

与 previous global state 相加后得到上式。

---

# 15. mechanical assembly 后的职责边界

一旦：

$$ \breve\sigma_{i+1} $$

$$ \breve\varepsilon^e_{i+1} $$

$$ \breve{\dot\varepsilon}^p_{i+1} $$

$$ \breve\varepsilon^p_{i+1} $$

已知，global mechanical branch完成。

后续 hardening / damage update不再要求新的 tower global equilibrium solve。

也就是说：

$$ \boxed{\text{hardening Eq. (73)–(74) 是 pointwise local-in-space global update}} $$

以及：

$$ \boxed{\text{damage Eq. (75) 也是 pointwise local-in-space finishing update}} $$

---

# 16. Eq. (73) 的一般形式

原论文 normal formulation 写：

$$ Z_{i+1}=\mathcal LX_{i+1} $$

其中：

$$ X $$

是 transformed hardening internal variables，

$$ Z $$

是对应 thermodynamic forces。

对于 current viscoplastic damage material：

$$ X = \begin{bmatrix} \alpha \\ \bar r \end{bmatrix} $$

以及：

$$ Z = \begin{bmatrix} \beta \\ \bar R \end{bmatrix} $$

---

# 17. current material 的 normal hardening law

kinematic hardening branch：

$$ \beta=C\alpha $$

isotropic transformed branch：

$$ \bar R=R_\infty\bar r $$

所以：

$$ \mathcal L = \begin{bmatrix} C & 0 \\ 0 & R_\infty \end{bmatrix} $$

这是线性 normal state law。

---

# 18. 为什么 hardening 不需要 PGD

Eq. (73)–(74) 在每个 material point 上独立。

输入：

$$ \hat{\dot\alpha}(t,q) $$

$$ \hat\beta(t,q) $$

$$ H_\beta(t,q) $$

以及：

$$ \hat{\dot{\bar r}}(t,q) $$

$$ \hat{\bar R}(t,q) $$

$$ H_{\bar R}(t,q) $$

输出：

$$ \breve\alpha(t,q),\breve\beta(t,q),\breve{\dot\alpha}(t,q) $$

与：

$$ \breve{\bar r}(t,q),\breve{\bar R}(t,q),\breve{\dot{\bar r}}(t,q) $$

空间点之间没有 equilibrium coupling。

因此 hardening branch 直接保留 full material-point histories 更自然。

---

# 19. Eq. (74) 的一般 descent relation

原论文 hardening descent relation：

$$ -(\dot X_{i+1}-\hat{\dot X}_{i+1/2})=H_Z(Z_{i+1}-\hat Z_{i+1/2}) $$

等价于：

$$ \dot X_{i+1}+H_ZZ_{i+1}=\hat{\dot X}_{i+1/2}+H_Z\hat Z_{i+1/2} $$

代入：

$$ Z_{i+1}=\mathcal LX_{i+1} $$

得到：

$$ \boxed{\dot X_{i+1}+H_Z\mathcal LX_{i+1}=\hat{\dot X}_{i+1/2}+H_Z\hat Z_{i+1/2}} $$

这是 linear first-order ODE。

---

# 20. kinematic hardening branch 的 ODE

对：

$$ X=\alpha $$

以及：

$$ Z=\beta=C\alpha $$

得到：

$$ -(\dot\alpha-\hat{\dot\alpha})=H_\beta(\beta-\hat\beta) $$

整理：

$$ \boxed{\dot\alpha+H_\beta C\alpha=\hat{\dot\alpha}+H_\beta\hat\beta} $$

所有量均是：

$$ (t,q) $$

dependent。

---

# 21. kinematic hardening backward-Euler 离散

对于：

$$ n\ge1 $$

定义：

$$ \dot\alpha_n=\frac{\alpha_n-\alpha_{n-1}}{\Delta t_n} $$

代入：

$$ \frac{\alpha_n-\alpha_{n-1}}{\Delta t_n}+H_{\beta,n}C\alpha_n=\hat{\dot\alpha}_n+H_{\beta,n}\hat\beta_n $$

两边乘：

$$ \Delta t_n $$

得到：

$$ \alpha_n-\alpha_{n-1}+\Delta t_nH_{\beta,n}C\alpha_n=\Delta t_n\hat{\dot\alpha}_n+\Delta t_nH_{\beta,n}\hat\beta_n $$

所以：

$$ \boxed{\alpha_n=\frac{\alpha_{n-1}+\Delta t_n(\hat{\dot\alpha}_n+H_{\beta,n}\hat\beta_n)}{1+\Delta t_nH_{\beta,n}C}} $$

---

# 22. kinematic force 与 rate reconstruction

得到：

$$ \alpha_n $$

后：

$$ \boxed{\beta_n=C\alpha_n} $$

以及：

$$ \boxed{\dot\alpha_n=\frac{\alpha_n-\alpha_{n-1}}{\Delta t_n}} $$

因此：

- $\alpha$ 为 integrated support history；
- $\beta$ 为 primary force；
- $\dot\alpha$ 为 primary rate。

未来 state interface 应保留这三个字段或能够无歧义恢复它们。

---

# 23. isotropic transformed hardening branch

Eq. (74) 对：

$$ X=\bar r $$

以及：

$$ Z=\bar R=R_\infty\bar r $$

得到：

$$ -(\dot{\bar r}-\hat{\dot{\bar r}})=H_{\bar R}(\bar R-\hat{\bar R}) $$

整理：

$$ \boxed{\dot{\bar r}+H_{\bar R}R_\infty\bar r=\hat{\dot{\bar r}}+H_{\bar R}\hat{\bar R}} $$

---

# 24. $\bar r$ backward-Euler 离散

对于：

$$ n\ge1 $$

有：

$$ \dot{\bar r}_n=\frac{\bar r_n-\bar r_{n-1}}{\Delta t_n} $$

因此：

$$ \boxed{\bar r_n=\frac{\bar r_{n-1}+\Delta t_n(\hat{\dot{\bar r}}_n+H_{\bar R,n}\hat{\bar R}_n)}{1+\Delta t_nH_{\bar R,n}R_\infty}} $$

然后：

$$ \boxed{\bar R_n=R_\infty\bar r_n} $$

以及：

$$ \boxed{\dot{\bar r}_n=\frac{\bar r_n-\bar r_{n-1}}{\Delta t_n}} $$

---

# 25. hardening denominator positivity

kinematic denominator：

$$ 1+\Delta t_nH_{\beta,n}C $$

isotropic denominator：

$$ 1+\Delta t_nH_{\bar R,n}R_\infty $$

由于：

$$ \Delta t_n>0 $$

$$ H_{\beta,n}>0 $$

$$ H_{\bar R,n}>0 $$

$$ C>0 $$

$$ R_\infty>0 $$

因此：

$$ \boxed{1+\Delta t_nH_{\beta,n}C>1} $$

以及：

$$ \boxed{1+\Delta t_nH_{\bar R,n}R_\infty>1} $$

hardening BE solve不存在与 Eq. (72) 类似的 controllability cancellation denominator。

---

# 26. hardening $t_0$ 的 state treatment

initial integrated hardening variables由 initial condition 给定。

virgin tower：

$$ \alpha_0=0 $$

$$ \bar r_0=0 $$

更一般：

$$ \alpha_0=\alpha_{\rm init} $$

$$ \bar r_0=\bar r_{\rm init} $$

因此在 global stage：

$$ \boxed{\breve\alpha_0=\hat\alpha_0=\alpha_{\rm init}} $$

以及：

$$ \boxed{\breve{\bar r}_0=\hat{\bar r}_0=\bar r_{\rm init}} $$

---

# 27. hardening $t_0$ force treatment

Eq. (73) 直接给：

$$ \boxed{\breve\beta_0=C\breve\alpha_0} $$

以及：

$$ \boxed{\breve{\bar R}_0=R_\infty\breve{\bar r}_0} $$

如果 local state 与 global candidate 使用相同 initial history，通常：

$$ \breve\beta_0=\hat\beta_0 $$

$$ \breve{\bar R}_0=\hat{\bar R}_0 $$

---

# 28. hardening $t_0$ rate treatment

Eq. (74) 在 $t_0$ 给：

$$ \breve{\dot\alpha}_0=\hat{\dot\alpha}_0-H_{\beta,0}(\breve\beta_0-\hat\beta_0) $$

以及：

$$ \breve{\dot{\bar r}}_0=\hat{\dot{\bar r}}_0-H_{\bar R,0}(\breve{\bar R}_0-\hat{\bar R}_0) $$

若 initial force 也一致，则退化为：

$$ \breve{\dot\alpha}_0=\hat{\dot\alpha}_0 $$

以及：

$$ \breve{\dot{\bar r}}_0=\hat{\dot{\bar r}}_0 $$

---

# 29. current 1D hardening implementation

当前：

```text
latin/global_stage.py
```

中的 `_update_hardening_variables()` 已采用：

- $t_0$ integrated history 从 local state 继承；
- $t_0$ rate 由 descent relation计算；
- $n\ge1$ 使用 backward Euler；
- $\beta=C\alpha$；
- $\bar R=R_\infty\bar r$。

因此这部分可视为 mature baseline。

tower v1 主要只需要：

$$ \boxed{\text{element index}\rightarrow\text{fiber material-point index }q} $$

而不需要理论重构。

---

# 30. Eq. (75) damage branch 的原论文内容

原论文选择：

$$ b^-=0 $$

因此 damage descent relation退化为：

$$ \boxed{\dot D_{i+1}=\hat{\dot D}_{i+1/2}} $$

这表示 global stage 不对 local damage evolution rate施加额外 correction。

随后：

$$ Y_{i+1} $$

通过 nonlinear damage state relation重新计算。

---

# 31. $b^-=0$ 的准确含义

必须避免把：

$$ b^-=0 $$

错误理解为：

$$ D_{i+1}=D_i $$

它真正表示：

$$ \dot D_{\rm global}=\dot D_{\rm local} $$

也就是说：

> damage evolution 在 current global half-step 中不沿额外 descent direction 修正。

因此 damage仍然可以随时间持续增长。

---

# 32. 连续理论下 integrated damage history

如果：

$$ D_{\rm global}(0)=D_{\rm local}(0)=D_0 $$

以及：

$$ \dot D_{\rm global}(t)=\dot D_{\rm local}(t) $$

则：

$$ D_{\rm global}(t)=D_0+\int_0^t\dot D_{\rm global}(\tau)d\tau $$

以及：

$$ D_{\rm local}(t)=D_0+\int_0^t\dot D_{\rm local}(\tau)d\tau $$

因此：

$$ \boxed{D_{\rm global}(t)=D_{\rm local}(t)} $$

这是由 Eq. (75) 与 common initial condition直接推导出的连续结论。

---

# 33. damage-history inheritance 的两个 possible implementation

数值上存在两个可能路径。

路径 A：

```text
copy local damage rate
    ↓
re-integrate damage rate in global stage
    ↓
obtain candidate D
```

路径 B：

```text
copy local damage rate
and
copy local already-integrated D history
```

current 1D 使用路径 A。

本阶段建议 tower v1 使用路径 B。

---

# 34. current local-stage damage integration

current `latin/local_stage.py` 对：

$$ [\varepsilon^p,\alpha,\bar r,D] $$

作为一个 coupled local internal-state vector。

每个 time step 使用 classical RK4 integration。

因此 local damage history：

$$ \hat D_n $$

是 RK4 integration result。

其 nodal damage rate：

$$ \hat{\dot D}_n $$

是在 current state 处重新评价 local damage evolution law得到。

---

# 35. current global-stage damage implementation

current `latin/global_stage.py` 先做：

$$ \dot D^{\rm new}=\hat{\dot D} $$

然后重新使用：

$$ D_n=D_{n-1}+\Delta t_n\dot D_n $$

即 backward-Euler / right-endpoint style integration。

因此 current 1D 形成：

$$ \boxed{\text{local }D=\text{RK4 history}} $$

而：

$$ \boxed{\text{global candidate }D=\text{BE re-integration of local nodal rates}} $$

---

# 36. 为什么两种 discrete damage histories 不一定相等

RK4 integration 使用 step 内多个 intermediate evaluations：

$$ k_1,k_2,k_3,k_4 $$

而 BE re-integration只使用：

$$ \dot D_n $$

所以一般：

$$ \hat D_n^{\rm RK4}\neq D_n^{\rm BE} $$

即使：

$$ \dot D_n^{\rm global}=\hat{\dot D}_n $$

成立。

因此 current 1D 的：

$$ D_{\rm candidate} $$

与：

$$ \hat D $$

存在纯 numerical integration discrepancy。

---

# 37. 这一 discrepancy 的理论地位

原论文 Eq. (75) 本身没有要求：

$$ \text{local history solver}=\text{RK4} $$

也没有要求：

$$ \text{global history reconstruction}=\text{BE} $$

因此 current 1D 的 mismatch 属于 implementation choice，而不是原论文显式算法要求。

这点必须在 tower v1 中明确记录。

---

# 38. tower v1 damage-history choice

本阶段建议正式冻结：

$$ \boxed{\breve{\dot D}_{i+1}=\hat{\dot D}_{i+1/2}} $$

同时：

$$ \boxed{\breve D_{i+1}=\hat D_{i+1/2}} $$

这里：

$$ \breve{} $$

表示 unrelaxed global candidate。

也就是说：

> current global stage 不对 local damage rate 修正，也不重新积分 local damage history。

---

# 39. 选择 direct-copy integrated $D$ 的理由 1：continuous Eq. (75)

连续理论已经给出：

$$ \dot D_{\rm global}=\dot D_{\rm local} $$

且 common initial condition下：

$$ D_{\rm global}=D_{\rm local} $$

所以 direct-copy：

$$ \breve D=\hat D $$

是 Eq. (75) 最直接的 discrete interpretation。

---

# 40. 选择 direct-copy integrated $D$ 的理由 2：避免双重时间离散

local stage已经完成 nonlinear damage evolution integration。

如果 global stage再次对 sampled nodal rate积分，就相当于同一 current local damage evolution 被：

- RK4；
- BE；

分别离散两次。

direct-copy 可以避免这层不必要差异。

---

# 41. 选择 direct-copy integrated $D$ 的理由 3：保持 local internal-state consistency

local stage 的：

$$ \hat D_n $$

与同一 RK4 history 下的：

$$ \hat\varepsilon^p_n $$

$$ \hat\alpha_n $$

$$ \hat{\bar r}_n $$

共同来自同一个 coupled internal-state integration。

直接继承：

$$ \hat D $$

可以保留这一 nonlinear local integration 的 internal consistency。

---

# 42. 选择 direct-copy integrated $D$ 的理由 4：职责更加清楚

tower v1 可以明确：

$$ \boxed{\text{local stage负责 nonlinear damage history integration}} $$

而：

$$ \boxed{\text{global stage负责 structural admissibility correction，不重新积分 damage evolution}} $$

这比 current 1D 的 mixed responsibility 更清楚。

---

# 43. direct-copy $D$ 不改变 damage structural projection

再次强调：

$$ \breve D=\hat D $$

不意味着 damage branch没有 global structural correction。

因为 nonlinear damaged elasticity导致：

$$ \hat\varepsilon^e $$

与 reference linear relation之间仍存在 residual。

所以：

$$ \Delta\varepsilon^R $$

仍需经过 tower reference equilibrium projection。

---

# 44. final energy-release rate 的计算时机

原论文在 Eq. (75) 之后才计算：

$$ Y_{i+1} $$

因此必须先完成：

$$ \breve\sigma_{i+1} $$

与：

$$ \breve D_{i+1} $$

再计算：

$$ \boxed{\breve Y_{i+1}=Y(\breve\sigma_{i+1},\breve D_{i+1})} $$

不能提前使用：

$$ \sigma_i $$

或：

$$ \hat\sigma $$

代替 final candidate stress。

---

# 45. current 1D unilateral energy-release relation

当前一维 material-point relation：

拉伸：

$$ Y=\frac{\sigma^2}{2E(1-D)^2} $$

压缩：

$$ Y=\frac{h\sigma^2}{2E(1-hD)^2} $$

因此 tower fiber material point $q$：

若：

$$ \breve\sigma_q\ge0 $$

则：

$$ \boxed{\breve Y_q=\frac{\breve\sigma_q^2}{2E_q(1-\breve D_q)^2}} $$

若：

$$ \breve\sigma_q\lt 0 $$

则：

$$ \boxed{\breve Y_q=\frac{h_q\breve\sigma_q^2}{2E_q(1-h_q\breve D_q)^2}} $$

---

# 46. $Y$ 为什么不直接复制 local value

local stage 中：

$$ \hat Y $$

作为 thermodynamic force history由 ascent-direction choice保持 current global force。

但是 global mechanical candidate 的：

$$ \breve\sigma $$

发生了变化。

因此 global candidate必须重新通过 nonlinear state relation得到：

$$ \breve Y $$

而不能简单：

$$ \breve Y=\hat Y $$

---

# 47. unrelaxed hardening / damage / energy finishing 的准确顺序

建议冻结：

```text
final plastic PGD correction
    ↓
plastic + damage structural mechanical assembly
    ↓
final unrelaxed sigma, eps_e, eps_p, eps_p_dot
    ↓
Eq. (73)-(74) hardening update
    ↓
Eq. (75) damage rate/history inheritance
    ↓
final Y = Y(sigma, D)
    ↓
complete breve{s}_{i+1}
```

hardening 与 damage history本身不依赖 final $Y$。

但：

$$ Y $$

依赖 final：

$$ \sigma $$

与：

$$ D $$

所以必须放在 finishing block 的最后。

---

# 48. complete unrelaxed formal LATIN state

原论文 primary state 可写：

$$ s=\{\dot\varepsilon^p,\varepsilon^e,\dot X,\dot D,\sigma,Z,Y\} $$

对于 current model：

$$ \dot X=\{\dot\alpha,\dot{\bar r}\} $$

以及：

$$ Z=\{\beta,\bar R\} $$

所以：

$$ \boxed{\breve s_{i+1}=\{\breve{\dot\varepsilon}^p,\breve\varepsilon^e,\breve{\dot\alpha},\breve{\dot{\bar r}},\breve{\dot D},\breve\sigma,\breve\beta,\breve{\bar R},\breve Y\}} $$

---

# 49. numerical support histories

代码还需要保存：

$$ \varepsilon^p $$

$$ \alpha $$

$$ \bar r $$

$$ D $$

这些是 integrated internal/support histories。

因此 future tower state object 实际需要比 paper formal state 多保存一层数据。

但必须在语义上标注：

$$ \boxed{\text{primary LATIN fields}} $$

与：

$$ \boxed{\text{integrated support histories}} $$

不同。

---

# 50. relaxation 的原论文形式

完成：

$$ \breve s_{i+1} $$

后：

$$ \boxed{s_{i+1}=\mu\breve s_{i+1}+(1-\mu)s_i} $$

原论文使用：

$$ \boxed{\mu=0.8} $$

因此：

$$ 1-\mu=0.2 $$

---

# 51. primary fields 的 relaxation

对任意 primary field：

$$ q\in\{\dot\varepsilon^p,\varepsilon^e,\dot\alpha,\dot{\bar r},\dot D,\sigma,\beta,\bar R,Y\} $$

采用：

$$ \boxed{q_{i+1}=(1-\mu)q_i+\mu\breve q_{i+1}} $$

这正是 formal LATIN relaxation。

---

# 52. integrated support histories 是否 relaxation

原论文 formal state没有把：

$$ \varepsilon^p,\alpha,\bar r,D $$

显式写入 $s$。

因此原论文没有逐项给出它们的 relaxation formula。

但是代码为了下一 local stage必须保存这些 histories。

本阶段建议：

$$ \boxed{\text{integrated support histories 使用与 primary state 相同的 convex relaxation}} $$

---

# 53. accumulated plastic strain relaxation

定义：

$$ \boxed{\varepsilon^p_{i+1}=(1-\mu)\varepsilon^p_i+\mu\breve\varepsilon^p_{i+1}} $$

这与 current `relax_global_state()` 保持一致。

---

# 54. hardening histories relaxation

定义：

$$ \boxed{\alpha_{i+1}=(1-\mu)\alpha_i+\mu\breve\alpha_{i+1}} $$

以及：

$$ \boxed{\bar r_{i+1}=(1-\mu)\bar r_i+\mu\breve{\bar r}_{i+1}} $$

---

# 55. damage history relaxation

定义：

$$ \boxed{D_{i+1}=(1-\mu)D_i+\mu\breve D_{i+1}} $$

注意：

$$ \breve D_{i+1}=\hat D_{i+1/2} $$

是 unrelaxed candidate inheritance。

但 accepted relaxed state一般：

$$ D_{i+1}\neq\hat D_{i+1/2} $$

因为 LATIN relaxation 会把它与 previous global history进行 convex blend。

---

# 56. 为什么 support-history relaxation 与 rate relaxation相容

假设：

$$ X^{(a)}(t)=X_0+\int_0^t\dot X^{(a)}(\tau)d\tau $$

以及：

$$ X^{(b)}(t)=X_0+\int_0^t\dot X^{(b)}(\tau)d\tau $$

定义 relaxed rate：

$$ \dot X^{(r)}=(1-\mu)\dot X^{(a)}+\mu\dot X^{(b)} $$

则：

$$ X^{(r)}=X_0+\int_0^t\dot X^{(r)}d\tau $$

展开：

$$ X^{(r)}=X_0+(1-\mu)\int_0^t\dot X^{(a)}d\tau+\mu\int_0^t\dot X^{(b)}d\tau $$

所以：

$$ \boxed{X^{(r)}=(1-\mu)X^{(a)}+\mu X^{(b)}} $$

因此连续理论下，rate relaxation 与 integrated-history relaxation一致。

---

# 57. discrete implementation 的 caveat

current tower / 1D implementation存在：

- local histories可能由 RK4；
- global hardening histories由 BE；
- plastic histories由 PGD/BE reconstruction；
- relaxed rate/history都直接 convex blend。

因此离散层面：

$$ X_{n}-X_{n-1}=\Delta t_n\dot X_n $$

不一定对所有 relaxed support histories严格机器精度成立。

所以 future tower implementation应把：

$$ \boxed{\text{continuous consistency}} $$

与：

$$ \boxed{\text{discrete support-history consistency diagnostic}} $$

分开。

---

# 58. hardening linear state law 在 relaxation 后是否保持

如果：

$$ \beta_i=C\alpha_i $$

以及：

$$ \breve\beta_{i+1}=C\breve\alpha_{i+1} $$

则：

$$ \beta_{i+1}=(1-\mu)\beta_i+\mu\breve\beta_{i+1} $$

代入：

$$ \beta_{i+1}=C[(1-\mu)\alpha_i+\mu\breve\alpha_{i+1}] $$

所以：

$$ \boxed{\beta_{i+1}=C\alpha_{i+1}} $$

exactly。

---

# 59. isotropic linear state law 在 relaxation 后也保持

同理：

$$ \bar R_i=R_\infty\bar r_i $$

以及：

$$ \breve{\bar R}_{i+1}=R_\infty\breve{\bar r}_{i+1} $$

因此：

$$ \boxed{\bar R_{i+1}=R_\infty\bar r_{i+1}} $$

convex relaxation不会破坏 Eq. (73)。

---

# 60. damage bounds 在 relaxation 后保持

如果：

$$ 0\le D_i\lt D_{\max} $$

以及：

$$ 0\le\breve D_{i+1}\lt D_{\max} $$

且：

$$ 0\lt \mu\le1 $$

则 convex blend满足：

$$ \boxed{0\le D_{i+1}\lt D_{\max}} $$

因此 relaxation本身不会生成超出两端值域的新 damage。

---

# 61. damage irreversibility 与 whole-time convex combination

若 previous global history：

$$ D_i(t_n) $$

随时间非减，

而 unrelaxed candidate：

$$ \breve D_{i+1}(t_n) $$

也随时间非减，

则：

$$ D_{i+1}(t_n)=(1-\mu)D_i(t_n)+\mu\breve D_{i+1}(t_n) $$

仍随时间非减。

因此：

$$ \boxed{\text{convex relaxation保留 whole-time damage monotonicity}} $$

前提是两条输入 histories都已单调。

---

# 62. energy-release rate $Y$ 在 relaxation 后的特殊性

unrelaxed candidate满足：

$$ \breve Y=Y(\breve\sigma,\breve D) $$

previous accepted state一般来自前一次 relaxed LATIN state。

如果直接 formal relaxation：

$$ Y_{i+1}=(1-\mu)Y_i+\mu\breve Y_{i+1} $$

由于：

$$ Y(\sigma,D) $$

是 nonlinear function，一般：

$$ \boxed{Y_{i+1}\neq Y(\sigma_{i+1},D_{i+1})} $$

---

# 63. 这是否意味着 relaxation 后必须重算 $Y$

本阶段建议：

$$ \boxed{\text{不重算}} $$

即保持：

$$ Y_{i+1}=(1-\mu)Y_i+\mu\breve Y_{i+1} $$

而不再执行：

$$ Y_{i+1}\leftarrow Y(\sigma_{i+1},D_{i+1}) $$

---

# 64. 为什么 relaxation 后不强制重算 $Y$

原论文 formal relaxation作用于完整 state：

$$ s_{i+1}=\mu\breve s_{i+1}+(1-\mu)s_i $$

而：

$$ Y $$

属于 formal state。

如果 relaxation 后又强制：

$$ Y=Y(\sigma,D) $$

就增加了一次额外 local nonlinear state-law projection。

这会把原本 LATIN 的：

```text
global candidate
    ↓
relax
    ↓
next local projection
```

变成：

```text
global candidate
    ↓
relax
    ↓
partial nonlinear reprojection of Y only
    ↓
next local projection
```

理论上不够干净。

---

# 65. local/global manifold mismatch 是 LATIN 的正常状态

relaxed：

$$ s_{i+1} $$

属于下一 iteration 的 global iterate。

它并不要求在 relaxation 后仍严格满足所有 nonlinear local constitutive relations。

下一 local stage正是用来把：

$$ s_{i+1} $$

重新投影到：

$$ \Gamma $$

得到：

$$ \hat s_{i+3/2} $$

因此：

$$ Y_{i+1}\neq Y(\sigma_{i+1},D_{i+1}) $$

可以视为 LATIN alternating process 的正常 local/global mismatch。

---

# 66. current 1D relaxation implementation

当前 `latin/iteration_control.py` 的 `relax_global_state()` 对以下字段全部做同一 convex blend：

- `plastic_strain_rate`
- `elastic_strain`
- `alpha_rate`
- `r_bar_rate`
- `damage_rate`
- `stress`
- `beta`
- `R_bar`
- `energy_release_rate`
- `plastic_strain`
- `alpha`
- `r_bar`
- `damage`

因此 current 1D 已经采用：

$$ \boxed{\text{primary fields + support histories uniform relaxation}} $$

tower v1 可以继承这一 architecture-level choice。

---

# 67. Eq. (76) 使用哪个 global state

原论文顺序是：

```text
global stage
    ↓
breve{s}_{i+1}
    ↓
relaxation
    ↓
s_{i+1}
    ↓
Eq. (76)
```

所以：

$$ \boxed{\text{Eq. (76) 必须使用 relaxed }s_{i+1}} $$

而不是 unrelaxed：

$$ \breve s_{i+1} $$

---

# 68. Eq. (76) 的 relative LATIN indicator

定义：

$$ \boxed{\xi_{i+1}=\frac{\|\hat s^p_{i+1/2}-s^p_{i+1}\|}{\|\hat s^p_{i+1/2}\|+\|s^p_{i+1}\|}} $$

其中：

$$ s^p $$

表示 mechanical subset。

这个 indicator 测量：

> current local constitutive state 与 current relaxed globally admissible mechanical state 之间的 relative LATIN distance。

---

# 69. Eq. (77) mechanical subset

对于 current model：

$$ s^p=\{\dot\varepsilon^p,\varepsilon^e,\dot\alpha,\dot{\bar r},\sigma,\beta,\bar R\} $$

因此直接进入 norm 的 seven field groups为：

$$ \dot\varepsilon^p $$

$$ \varepsilon^e $$

$$ \dot\alpha $$

$$ \dot{\bar r} $$

$$ \sigma $$

$$ \beta $$

$$ \bar R $$

---

# 70. damage variables 为什么不直接进入 Eq. (77)

以下 fields：

$$ D $$

$$ \dot D $$

$$ Y $$

不直接出现在 Eq. (77) mechanical norm。

但是它们仍通过以下路径影响：

- local viscoplastic flow；
- nonlinear elastic strain；
- stress history；
- search directions；
- future local stage；
- damage structural residual。

因此：

$$ \boxed{\text{not in norm}\neq\text{irrelevant}} $$

---

# 71. Eq. (77) tower material-point norm 回顾

tower 离散 norm 应基于 material-point integration metric：

$$ M=\operatorname{diag}(v_q) $$

并使用 search-direction metrics。

continuous-style form可写：

$$ \|s^p\|^2=\int_{\Omega\times I}\left[\sigma^TH_\sigma\sigma+\beta^TH_\beta\beta+\bar R^TH_{\bar R}\bar R+\dot\varepsilon^{pT}H_\sigma^{-1}\dot\varepsilon^p+\varepsilon^{eT}C_0\varepsilon^e+\dot\alpha^TH_\beta^{-1}\dot\alpha+\dot{\bar r}^TH_{\bar R}^{-1}\dot{\bar r}\right]d\Omega dt $$

具体离散 quadrature 已在前序 Eq. (60)/Eq. (76) 阶段讨论。

---

# 72. Eq. (76) 计算前不能遗漏 finishing variables

因为 Eq. (77) 包含：

$$ \beta $$

$$ \bar R $$

$$ \dot\alpha $$

$$ \dot{\bar r} $$

所以不能在 plastic mechanical assembly刚完成后就直接计算：

$$ \xi $$

必须先完成 Eq. (73)–(74)。

这从数学上证明：

$$ \boxed{\text{hardening finishing不是附属 postprocessing，而是 Eq. (76) 前的必需步骤}} $$

---

# 73. $Y$ 不进入 Eq. (77) 是否意味着可以晚一点计算

虽然：

$$ Y $$

不直接进入 Eq. (77)，但完整：

$$ s_{i+1} $$

需要作为下一 local stage输入。

下一 local stage 的 damage evolution需要：

$$ Y_{i+1} $$

因此 current global iteration完成前必须更新：

$$ Y $$

不能把它延迟到下一 local stage临时再算。

---

# 74. Eq. (60) saturation indicator

定义：

$$ \boxed{\zeta_i=\frac{\xi_i-\xi_{i+1}}{\xi_i+\xi_{i+1}}} $$

其作用不是 nonlinear convergence本身，而是评价：

> 使用 current existing PGD basis 完成一次 LATIN global update 后，LATIN distance 是否取得足够 improvement。

---

# 75. baseline $\xi_i$ 的语义

在 current LATIN iteration $i$ 开始时：

$$ \xi_i $$

是 previous accepted relaxed global iterate的 indicator。

它在 current iteration 内必须视为固定 baseline。

也就是说，在同一个 current local state下做：

- fixed-basis trial；
- enriched-basis trial；

都应与相同：

$$ \xi_i $$

比较。

---

# 76. 为什么 trial indicator 不能相互当 baseline

如果 fixed-basis trial得到：

$$ \xi_{i+1}^{\rm up} $$

随后 enrichment trial得到：

$$ \xi_{i+1}^{\rm enr} $$

不能定义：

$$ \zeta=\frac{\xi_{i+1}^{\rm up}-\xi_{i+1}^{\rm enr}}{\xi_{i+1}^{\rm up}+\xi_{i+1}^{\rm enr}} $$

并把它称为原论文 Eq. (60)。

因为 Eq. (60) 的 baseline是 successive LATIN iteration indicator：

$$ \xi_i $$

而不是同一 iteration 内两个 trial candidates之间的 residual ratio。

---

# 77. fixed-basis Trial A

current iteration首先使用 existing accepted spatial basis：

$$ \mathcal B_m $$

保持 spatial modes固定。

执行：

$$ \boxed{\text{Eq. (58)–(59) all-mode temporal update}} $$

得到 current fixed-basis plastic correction。

---

# 78. Trial A 的 complete unrelaxed candidate

fixed-basis plastic correction必须经过完整 finishing：

```text
plastic PGD reconstruction
    ↓
damage structural projection
    ↓
mechanical assembly
    ↓
hardening Eq. (73)-(74)
    ↓
damage Eq. (75)
    ↓
Y update
    ↓
breve{s}_{i+1}^{up}
```

不能只用 reduced residual决定 Eq. (60)。

---

# 79. Trial A relaxation

定义：

$$ \boxed{s_{i+1}^{\rm up}=\mu\breve s_{i+1}^{\rm up}+(1-\mu)s_i} $$

然后计算：

$$ \boxed{\xi_{i+1}^{\rm up}} $$

---

# 80. Trial A saturation

定义：

$$ \boxed{\zeta_i^{\rm up}=\frac{\xi_i-\xi_{i+1}^{\rm up}}{\xi_i+\xi_{i+1}^{\rm up}}} $$

这一 quantity才用于判断：

$$ \mathcal B_m $$

是否需要 enrichment。

---

# 81. fixed-basis trial 的三种典型结果

第一种：

$$ \xi_{i+1}^{\rm up}\le\xi_{\rm tol} $$

说明 LATIN nonlinear solve 已满足 absolute convergence criterion。

第二种：

$$ \xi_{i+1}^{\rm up}>\xi_{\rm tol} $$

但：

$$ \zeta_i^{\rm up} $$

足够大。

说明 current basis仍能带来明显 LATIN improvement，可接受 fixed-basis global state并进入下一 LATIN iteration。

第三种：

$$ \xi_{i+1}^{\rm up}>\xi_{\rm tol} $$

且：

$$ \zeta_i^{\rm up} $$

太小。

说明 current basis可能已饱和，需要进入 Add-a-pair branch。

---

# 82. reduced residual 与 Eq. (60) saturation 的角色分工

Eq. (60)：

$$ \zeta $$

是 outer LATIN improvement indicator。

reduced mechanical residual：

$$ R_m $$

是 global plastic reduced solve的 internal diagnostic。

因此建议 enrichment触发继续使用“双信号”思路：

$$ \boxed{\text{outer saturation suggests basis deficiency}} $$

同时：

$$ \boxed{\text{reduced residual confirms there is unresolved plastic mechanical defect}} $$

避免只凭一个 quantity决定所有 enrichment。

---

# 83. current 1D 成熟经验

validated 1D implementation已经证明：

- 只用 Eq. (60) saturation作为 nonlinear stop并不稳健；
- $\xi$ absolute criterion需要保留；
- stagnation criterion需要保留；
- $\zeta$ 更适合作为 enrichment / basis adequacy indicator；
- reduced residual可辅助避免无意义 enrichment。

tower v1 不应丢弃这些经验。

---

# 84. enrichment event 开始时哪些数据必须冻结

若 fixed-basis Trial A判定 basis不足，进入 one-pair enrichment。

必须保持 current LATIN iteration的以下数据不变：

$$ s_i $$

$$ \hat s_{i+1/2} $$

$$ H_\sigma $$

$$ H_\beta $$

$$ H_{\bar R} $$

$$ \Delta\varepsilon^R $$

以及由 current local/global state定义的 plastic forcing。

不能重新调用 local stage。

---

# 85. 为什么 enrichment 内不能重做 local stage

Eq. (61)–(72) 的 new pair是针对：

$$ \text{current local/global pair} $$

的 remaining defect。

如果 inner enrichment期间重新计算 local stage，则：

$$ \hat s_{i+1/2} $$

改变，

从而：

$$ \bar\Delta $$

也改变。

这会把一个 rank-one fixed-point problem变成 moving-target problem。

因此：

$$ \boxed{\text{one enrichment event 内 local state 和 search directions frozen}} $$

---

# 86. one-pair enrichment branch

在 frozen current iteration data下执行：

```text
Eq. (61)-(64)
    ↓
Eq. (65)-(71)
    ↓
Eq. (72)
    ↓
alternating fixed point
    ↓
fixed-point convergence
    ↓
weighted Gram-Schmidt
    ↓
exact temporal coordinate transform
    ↓
gamma_sp
    ↓
gamma_lambda
    ↓
all-mode re-optimisation
    ↓
Delta_res
    ↓
accept / rollback one candidate pair
```

这部分已由前序阶段闭合。

---

# 87. accepted pair 后 basis state

如果 new mode通过：

- inner fixed-point validity；
- spatial novelty；
- modified temporal significance；
- full residual benefit；

则：

$$ \mathcal B_m\rightarrow\mathcal B_{m+1} $$

此时 basis-level acceptance完成。

但 current LATIN iteration还没有完成。

---

# 88. enriched candidate 为什么必须重新完整 assembly

new mode改变：

$$ \Delta\varepsilon^p $$

$$ \Delta\dot\varepsilon^p $$

$$ \Delta\sigma' $$

$$ \Delta\varepsilon' $$

因此 final:

$$ \breve\sigma $$

与：

$$ \breve\varepsilon^e $$

都会变化。

所以需要重新形成：

$$ \breve s_{i+1}^{\rm enr} $$

而不是只更新 reduced residual history。

---

# 89. same-iteration hardening branch 是否变化

Eq. (73)–(74) 的输入是：

$$ \hat{\dot X} $$

$$ \hat Z $$

$$ H_Z $$

这些在 same enrichment event 内冻结。

因此如果 hardening equations不显式依赖 final mechanical correction，则：

$$ \boxed{\breve X^{\rm enr}=\breve X^{\rm up}} $$

以及：

$$ \boxed{\breve Z^{\rm enr}=\breve Z^{\rm up}} $$

在数学上相同。

---

# 90. same-iteration damage history 是否变化

因为：

$$ \breve{\dot D}=\hat{\dot D} $$

以及 tower v1选择：

$$ \breve D=\hat D $$

local stage未重算，因此：

$$ \boxed{\breve D^{\rm enr}=\breve D^{\rm up}} $$

以及：

$$ \boxed{\breve{\dot D}^{\rm enr}=\breve{\dot D}^{\rm up}} $$

---

# 91. enriched candidate 中 $Y$ 必须重新计算

虽然：

$$ D $$

same，

但：

$$ \breve\sigma^{\rm enr}\neq\breve\sigma^{\rm up} $$

一般成立。

因此：

$$ \boxed{\breve Y^{\rm enr}=Y(\breve\sigma^{\rm enr},\breve D)} $$

必须重新计算。

这是 same-iteration enrichment finishing 中唯一明确依赖 new mechanical candidate 的 damage-related primary force。

---

# 92. enriched candidate relaxation

得到完整：

$$ \breve s_{i+1}^{\rm enr} $$

后：

$$ \boxed{s_{i+1}^{\rm enr}=\mu\breve s_{i+1}^{\rm enr}+(1-\mu)s_i} $$

必须使用 same：

$$ s_i $$

作为 relaxation baseline。

不能对：

$$ s_{i+1}^{\rm up} $$

再做第二次 incremental relaxation。

---

# 93. 为什么不能从 fixed-basis relaxed state继续 blend

错误做法：

$$ s_{i+1}^{\rm enr}=(1-\mu)s_{i+1}^{\rm up}+\mu\breve s_{i+1}^{\rm enr} $$

会把 same LATIN iteration的 previous trial state错误当作 outer baseline。

正确做法始终是：

$$ \boxed{s_{i+1}^{\rm trial}=(1-\mu)s_i+\mu\breve s_{i+1}^{\rm trial}} $$

无论 trial是 fixed-basis还是 enriched-basis。

---

# 94. enriched LATIN indicator

计算：

$$ \boxed{\xi_{i+1}^{\rm enr}=\frac{\|\hat s^p_{i+1/2}-s_{i+1}^{p,\rm enr}\|}{\|\hat s^p_{i+1/2}\|+\|s_{i+1}^{p,\rm enr}\|}} $$

local half-step：

$$ \hat s_{i+1/2} $$

保持 same。

---

# 95. enriched saturation quantity

若需要记录：

$$ \zeta_i^{\rm enr} $$

则仍使用 same baseline：

$$ \boxed{\zeta_i^{\rm enr}=\frac{\xi_i-\xi_{i+1}^{\rm enr}}{\xi_i+\xi_{i+1}^{\rm enr}}} $$

这可作为 enriched trial对 outer LATIN improvement 的诊断。

---

# 96. one-pair-per-enrichment-event 的准确含义

本阶段进一步明确：

$$ \boxed{\text{tower v1 一个 LATIN iteration 内最多接受一个 new PGD pair}} $$

也就是：

- fixed-basis trial可发生；
- 若 basis不足，可生成一个 new pair；
- new pair若接受，则形成 enriched candidate；
- current LATIN iteration结束；
- 下一 pair只能在下一 LATIN iteration重新判断后生成。

---

# 97. 为什么这个 baseline 更适合第一版 tower solver

第一，原论文 flow更容易逐步对应。

第二，每个 accepted LATIN iterate只有：

$$ 0 $$

或：

$$ 1 $$

个 spatial basis dimension increase。

第三，outer：

$$ \xi,\zeta $$

与 basis size变化之间更容易追踪。

第四，debug时可以明确定位：

- local nonlinear change；
- global fixed-basis change；
- one enrichment change；
- relaxation change。

第五，避免第一版 tower solver同时引入 multiple-mode acceleration。

---

# 98. current 1D 与 tower v1 在这一点上的差异

current 1D `pgd_global_stage.py` 和 `pgd_solver.py` 支持：

$$ \boxed{\text{same LATIN iteration 内 multiple enrichments}} $$

且 baseline：

$$ \xi_i $$

在 same iteration enrichment过程中保持不变。

tower v1选择：

$$ \boxed{\text{same iteration 内最多一个 accepted pair}} $$

属于 paper-traceable simplification，而不是对 current 1D成熟策略的否定。

---

# 99. negative $\zeta$ 的处理

如果：

$$ \xi_{i+1}>\xi_i $$

则：

$$ \zeta_i\lt 0 $$

它表示：

$$ \boxed{\text{current trial相对 previous accepted LATIN indicator没有 improvement}} $$

但不应立刻解释为：

$$ \boxed{\text{nonlinear solver必然失败}} $$

---

# 100. 为什么不强制 $\xi$ 单调

LATIN 是 alternating nonlinear method。

单个 iteration 中：

- local nonlinear projection；
- approximate reduced global projection；
- relaxation；

组合后未必保证每一步：

$$ \xi_{i+1}\lt \xi_i $$

严格成立。

因此 tower v1 不增加：

$$ \boxed{\xi_{i+1}\lt \xi_i} $$

作为 universal hard acceptance rule。

---

# 101. negative $\zeta$ 更适合做什么

negative或很小：

$$ \zeta $$

更适合作为：

- basis inadequacy signal；
- poor global improvement diagnostic；
- stagnation information；
- enrichment trigger input；

而不是 standalone nonlinear termination rule。

---

# 102. absolute LATIN convergence 仍然是主 nonlinear criterion

validated 1D经验支持：

$$ \boxed{\xi\le\xi_{\rm tol}} $$

作为 absolute nonlinear convergence criterion。

例如 current mature baseline：

$$ \xi_{\rm tol}=10^{-4} $$

这一数量级来自1D reproduction经验，不是原论文 universal constant。

tower v1初期可沿用为 provisional value，但最终仍需 benchmark calibration。

---

# 103. stagnation criterion 的必要性

current 1D还使用：

$$ \xi\le\xi_{\rm stag} $$

并且：

$$ |\xi_i-\xi_{i+1}|\le\varepsilon_{\rm stag} $$

持续若干 accepted iterations作为 practical stagnation stop。

这一经验在 tower v1 中仍值得保留。

因为 Eq. (60) saturation不能单独替代 nonlinear convergence。

---

# 104. basis saturation 与 nonlinear stagnation 必须区分

basis saturation：

$$ \zeta $$

主要回答：

> current PGD basis是否还能有效改善 global LATIN approximation？

nonlinear stagnation：

$$ |\xi_i-\xi_{i+1}| $$

主要回答：

> current accepted LATIN iterations是否已经几乎不再变化？

二者不能用同一个 threshold或 termination reason。

---

# 105. current iteration 的 state hierarchy

本阶段正式冻结四个核心 state layers：

第一：

$$ \boxed{s_i} $$

previous accepted relaxed global state。

第二：

$$ \boxed{\hat s_{i+1/2}} $$

current nonlinear local state。

第三：

$$ \boxed{\breve s_{i+1}} $$

current complete unrelaxed global candidate。

第四：

$$ \boxed{s_{i+1}} $$

current relaxed global state。

---

# 106. state hierarchy 的控制流

```text
accepted relaxed state
s_i
    ↓
local nonlinear projection
    ↓
hat{s}_{i+1/2}
    ↓
global correction + finishing
    ↓
breve{s}_{i+1}
    ↓
relaxation
    ↓
s_{i+1}
    ↓
Eq. (76)
    ↓
xi_{i+1}
```

这是 future tower solver最基本的 state machine。

---

# 107. same iteration fixed-basis / enriched trial hierarchy

fixed-basis trial：

$$ \breve s_{i+1}^{\rm up} $$

$$ s_{i+1}^{\rm up} $$

$$ \xi_{i+1}^{\rm up} $$

若 enrichment：

$$ \breve s_{i+1}^{\rm enr} $$

$$ s_{i+1}^{\rm enr} $$

$$ \xi_{i+1}^{\rm enr} $$

两套 trial共享：

$$ s_i $$

与：

$$ \hat s_{i+1/2} $$

不能把前一个 trial变成后一 trial的 outer baseline。

---

# 108. current fixed-basis trial 的 complete data flow

```text
INPUT
    s_i
    B_m
    xi_i

LOCAL
    s_i -> hat{s}_{i+1/2}

SEARCH DIRECTIONS
    hat{s}_{i+1/2}
        ->
    H_sigma, H_beta, H_Rbar

DAMAGE STRUCTURAL SOURCE
    s_i + hat{s}_{i+1/2}
        ->
    Delta eps_R
        ->
    full tower reference equilibrium projection

PGD FIXED-BASIS UPDATE
    B_m
        ->
    Eq. (58)-(59)
        ->
    updated temporal coefficients

PLASTIC RECONSTRUCTION
    Delta eps_p
    Delta eps_p_dot
    Delta eps_prime
    Delta sigma_prime

MECHANICAL ASSEMBLY
    + damage structural correction
        ->
    breve sigma
    breve eps_e
    breve eps_p
    breve eps_p_dot

HARDENING
    Eq. (73)-(74)

DAMAGE
    breve D_dot = hat D_dot
    breve D = hat D

ENERGY
    breve Y = Y(breve sigma, breve D)

COMPLETE
    breve{s}_{i+1}^{up}

RELAX
    s_{i+1}^{up}

INDICATOR
    xi_{i+1}^{up}

SATURATION
    zeta_i^{up}
```

---

# 109. current enrichment trial 的 complete data flow

```text
IF basis insufficient:

freeze
    s_i
    hat{s}_{i+1/2}
    H_sigma
    H_beta
    H_Rbar
    damage structural source
    full plastic forcing

ADD ONE PAIR
    Eq. (61)-(72)
        ->
    fixed point
        ->
    Gram-Schmidt
        ->
    mode significance
        ->
    enlarged-basis temporal re-optimisation
        ->
    mode accepted

RECONSTRUCT ENRICHED PLASTIC BRANCH

REBUILD COMPLETE
    breve{s}_{i+1}^{enr}

reuse
    same hardening candidate
    same D
    same D_dot

recompute
    breve Y from enriched stress

RELAX
    from same s_i

COMPUTE
    xi_{i+1}^{enr}
    zeta_i^{enr}

ACCEPT
    B_{m+1}
    s_{i+1}^{enr}

END current LATIN iteration
```

---

# 110. same-iteration hardening reuse 的 implementation consequence

理论上：

$$ \breve X^{\rm up}=\breve X^{\rm enr} $$

所以未来实现可以：

- 第一次 Trial A 时求 hardening；
- enrichment accepted后直接 reuse hardening fields。

但为了第一版代码可读性，也可以重新调用 deterministic hardening update，并用 unit test验证两者相同。

优先原则：

$$ \boxed{\text{correctness first, micro-optimisation later}} $$

---

# 111. same-iteration damage reuse 的 implementation consequence

tower v1：

$$ \breve D=\hat D $$

因此 fixed-basis与 enriched trial 的 unrelaxed damage history完全相同。

不需要重新积分。

只需重新计算：

$$ \breve Y $$

因为 stress changed。

---

# 112. damage-history choice 对 future code interface 的影响

future global-stage finishing函数可以接受：

```text
local_state
mechanical_candidate
search_directions
materials
```

然后：

- hardening从 local state + search directions求；
- candidate damage rate直接来自 local state；
- candidate damage直接来自 local state；
- $Y$ 用 mechanical candidate stress + copied damage求。

这样 damage ownership非常明确。

---

# 113. current 1D 是否需要立刻修改

本阶段不建议立即修改 current 1D。

原因：

- 当前 1D 是 validated reference；
- 其 convergence/history已经建立；
- 直接修改 damage reconstruction可能改变 baseline结果；
- 当前目标是 tower theory closure。

未来应做 controlled comparison：

```text
Variant A:
local RK4 D + global BE re-integration

Variant B:
local RK4 D + direct-copy D
```

比较：

- LATIN convergence；
- final damage；
- stress；
- mode count；
- runtime；
- FOM agreement。

---

# 114. tower v1 damage choice 的理论标签

未来文档与代码注释应写清：

$$ \boxed{\text{paper-consistent tower-v1 refinement}} $$

而不能写成：

$$ \boxed{\text{paper explicitly says copy integrated D}} $$

因为原论文显式写的是：

$$ \dot D_{i+1}=\hat{\dot D}_{i+1/2} $$

integrated $D$ direct-copy 是基于连续 relation与 common initial condition得到的离散策略选择。

---

# 115. future consistency test：damage candidate identity

tower v1 unrelaxed candidate应检查：

$$ \boxed{\breve D-\hat D=0} $$

到 machine precision。

同样：

$$ \boxed{\breve{\dot D}-\hat{\dot D}=0} $$

到 machine precision。

---

# 116. future consistency test：hardening state laws

每个 time-material point检查：

$$ \boxed{\breve\beta-C\breve\alpha=0} $$

以及：

$$ \boxed{\breve{\bar R}-R_\infty\breve{\bar r}=0} $$

---

# 117. future consistency test：hardening BE residual

对：

$$ n\ge1 $$

检查：

$$ r_{\alpha,n}=\frac{\alpha_n-\alpha_{n-1}}{\Delta t_n}+H_{\beta,n}C\alpha_n-\hat{\dot\alpha}_n-H_{\beta,n}\hat\beta_n $$

应满足：

$$ \boxed{|r_{\alpha,n}|\approx0} $$

---

# 118. future consistency test：isotropic hardening BE residual

定义：

$$ r_{\bar r,n}=\frac{\bar r_n-\bar r_{n-1}}{\Delta t_n}+H_{\bar R,n}R_\infty\bar r_n-\hat{\dot{\bar r}}_n-H_{\bar R,n}\hat{\bar R}_n $$

检查：

$$ \boxed{|r_{\bar r,n}|\approx0} $$

---

# 119. future consistency test：unrelaxed $Y$

每个 material point检查：

$$ \boxed{\breve Y-Y(\breve\sigma,\breve D)=0} $$

到 floating-point tolerance。

---

# 120. future consistency test：relaxed hardening state law

relaxation后继续检查：

$$ \boxed{\beta_{i+1}=C\alpha_{i+1}} $$

以及：

$$ \boxed{\bar R_{i+1}=R_\infty\bar r_{i+1}} $$

这两个应仍然 machine-precision成立。

---

# 121. future diagnostic：relaxed $Y$ manifold defect

由于 tower v1不在 relaxation后重算 $Y$，可记录：

$$ r_Y^{\rm relax}=Y_{i+1}-Y(\sigma_{i+1},D_{i+1}) $$

该 quantity不要求为零。

它可以作为：

$$ \boxed{\text{relaxed state 到 local damage-energy manifold 的 distance diagnostic}} $$

但不作为 immediate error。

---

# 122. future consistency test：relaxation formula

对所有 primary / support fields检查：

$$ q_{i+1}-(1-\mu)q_i-\mu\breve q_{i+1}=0 $$

这样可以确保 trial state没有 accidental cumulative blending。

---

# 123. future test：same baseline for multiple trials

在同一个 LATIN iteration内制造 fixed-basis与 enriched trial。

检查二者都使用：

$$ s_i $$

作为 relaxation baseline。

即：

$$ s_{i+1}^{\rm up}=(1-\mu)s_i+\mu\breve s_{i+1}^{\rm up} $$

以及：

$$ s_{i+1}^{\rm enr}=(1-\mu)s_i+\mu\breve s_{i+1}^{\rm enr} $$

---

# 124. future test：same local state for enrichment trial

记录：

$$ \hat s_{i+1/2}^{\rm before} $$

进入 enrichment后验证：

$$ \boxed{\hat s_{i+1/2}^{\rm after}=\hat s_{i+1/2}^{\rm before}} $$

search directions同样不应变化。

---

# 125. future test：hardening same across same-iteration trials

fixed-basis与 enriched trial中：

$$ \boxed{\breve\alpha^{\rm enr}=\breve\alpha^{\rm up}} $$

$$ \boxed{\breve\beta^{\rm enr}=\breve\beta^{\rm up}} $$

$$ \boxed{\breve{\bar r}^{\rm enr}=\breve{\bar r}^{\rm up}} $$

$$ \boxed{\breve{\bar R}^{\rm enr}=\breve{\bar R}^{\rm up}} $$

应成立。

---

# 126. future test：damage same across same-iteration trials

检查：

$$ \boxed{\breve D^{\rm enr}=\breve D^{\rm up}} $$

以及：

$$ \boxed{\breve{\dot D}^{\rm enr}=\breve{\dot D}^{\rm up}} $$

---

# 127. future test：$Y$ changes if enriched stress changes

若：

$$ \breve\sigma^{\rm enr}\neq\breve\sigma^{\rm up} $$

则一般：

$$ \boxed{\breve Y^{\rm enr}\neq\breve Y^{\rm up}} $$

这验证 $Y$ 在 each trial finishing中确实重新计算。

---

# 128. future test：Eq. (76) uses relaxed state

计算：

$$ \xi_{\rm relaxed} $$

与：

$$ \xi_{\rm unrelaxed} $$

未来测试应确保 solver history记录前者，而不是后者。

---

# 129. future test：Eq. (60) baseline fixed

在 one current LATIN iteration内：

$$ \zeta_i^{\rm up} $$

和：

$$ \zeta_i^{\rm enr} $$

都必须使用：

$$ \xi_i $$

作为 previous baseline。

---

# 130. future test：one-pair-per-iteration

tower v1 baseline应保证：

$$ m_{i+1}-m_i\in\{0,1\} $$

其中：

$$ m_i $$

为 accepted basis dimension。

不能在一个 LATIN iteration中：

$$ m_{i+1}-m_i>1 $$

除非未来显式开启 accelerated mode。

---

# 131. current 1D 与 tower v1 finishing 的相同点

保留：

- plastic / damage structural correction分离；
- hardening pointwise BE；
- final $Y$ 由 final unrelaxed stress与 damage计算；
- $\mu=0.8$ relaxation；
- Eq. (76) 在 relaxation后计算；
- mechanical norm不显式包含 damage variables；
- support histories显式存储；
- support histories与 primary fields一起 relaxation。

---

# 132. current 1D 与 tower v1 finishing 的不同点

主要差异：

第一，current 1D unrelaxed damage：

$$ D $$

由 local $\dot D$重新 BE integration。

tower v1：

$$ \boxed{\breve D=\hat D} $$

第二，current 1D global-stage framework支持 one call 内multiple enrichments。

tower v1 baseline：

$$ \boxed{\text{one accepted pair per LATIN iteration}} $$

第三，tower v1 time metrics尽量统一 BE/right-endpoint conventions。

第四，tower v1 future state/data contracts将更明确区分 primary state与 integrated histories。

---

# 133. 原论文 explicit 与 tower inference 的最终边界

原论文 explicit：

$$ \dot D_{i+1}=\hat{\dot D}_{i+1/2} $$

tower inference / choice：

$$ \breve D_{i+1}=\hat D_{i+1/2} $$

原论文 explicit：

$$ s_{i+1}=\mu\breve s_{i+1}+(1-\mu)s_i $$

tower implementation choice：

- integrated support histories也进行 same convex blend；
- relaxation后不再强制重算 nonlinear $Y$。

---

# 134. 为什么本阶段不把 $D$ 直接加入 Eq. (77)

即使 tower v1将 damage history处理得更严格，也不应擅自改变原论文 mechanical convergence norm。

因此：

$$ D,\dot D,Y $$

仍不直接进入 Eq. (77)。

如果未来发现 damage convergence需要额外监控，可增加：

$$ \boxed{\text{secondary damage diagnostic}} $$

但不能把它偷偷改写成“原论文 Eq. (77)”。

---

# 135. possible future damage diagnostic

未来可记录：

$$ \xi_D=\frac{\|\hat D-D\|}{\|\hat D\|+\|D\|} $$

或者：

$$ \xi_{\dot D}=\frac{\|\hat{\dot D}-\dot D\|}{\|\hat{\dot D}\|+\|\dot D\|} $$

但在 tower v1 current choice下，unrelaxed：

$$ \breve D=\hat D $$

而 relaxed：

$$ D_{i+1}\neq\hat D $$

因此 damage diagnostic主要反映 relaxation-induced local/global gap。

现阶段只建议保留为 future option。

---

# 136. current LATIN iteration 的完整 frozen/updated data map

固定：

$$ s_i $$

local stage后固定：

$$ \hat s_{i+1/2} $$

$$ H_\sigma,H_\beta,H_{\bar R} $$

global fixed-basis trial更新：

$$ \lambda_j,\dot\lambda_j $$

$$ \Delta\varepsilon^p $$

$$ \Delta\dot\varepsilon^p $$

$$ \breve\sigma,\breve\varepsilon^e $$

hardening：

$$ \breve\alpha,\breve\beta,\breve{\bar r},\breve{\bar R} $$

damage：

$$ \breve D=\hat D $$

$$ \breve{\dot D}=\hat{\dot D} $$

energy：

$$ \breve Y $$

relaxation：

$$ s_{i+1}^{\rm trial} $$

indicator：

$$ \xi_{i+1}^{\rm trial} $$

---

# 137. enriched trial 只改变哪些 finishing inputs

相对于 fixed-basis trial，accepted new PGD pair主要改变：

$$ \Delta\varepsilon^p $$

$$ \Delta\dot\varepsilon^p $$

$$ \Delta\sigma' $$

$$ \Delta\varepsilon' $$

所以改变：

$$ \breve\sigma $$

$$ \breve\varepsilon^e $$

$$ \breve\varepsilon^p $$

$$ \breve{\dot\varepsilon}^p $$

最终改变：

$$ \breve Y $$

以及 relaxed state与：

$$ \xi $$

---

# 138. enriched trial 不改变哪些 local/global inputs

same current iteration内不改变：

$$ \hat{\dot\alpha} $$

$$ \hat\beta $$

$$ \hat{\dot{\bar r}} $$

$$ \hat{\bar R} $$

$$ \hat{\dot D} $$

$$ \hat D $$

$$ H_\beta $$

$$ H_{\bar R} $$

所以 hardening candidate和 damage history可以保持 same。

---

# 139. Eq. (47)–(77) 主算法链的当前闭合状态

到本阶段结束，主链可写：

```text
elastic initialisation
    ↓
accepted global state s_i in A
    ↓
local stage
    ↓
hat{s}_{i+1/2} in Gamma
    ↓
search directions
    ↓
global damage residual projection
    ↓
PGD fixed-basis update Eq. (58)-(59)
    ↓
if needed:
    Eq. (61)-(72) add one pair
    ↓
mechanical assembly
    ↓
hardening Eq. (73)-(74)
    ↓
damage Eq. (75)
    ↓
Y from nonlinear state law
    ↓
breve{s}_{i+1}
    ↓
relaxation mu=0.8
    ↓
s_{i+1}
    ↓
Eq. (76)-(77) LATIN indicator
    ↓
Eq. (60) saturation / basis adequacy
    ↓
absolute / stagnation nonlinear control
    ↓
next LATIN iteration or stop
```

这意味着：

$$ \boxed{\text{Eq. (47)–(77) 的主结构现在已经基本贯通}} $$

---

# 140. 当前仍未完全冻结的最后理论接口

虽然主算法链贯通，但 future implementation前还有一个最后的 interface-level理论问题。

必须逐字段明确：

- 哪些字段属于 formal primary LATIN state；
- 哪些属于 integrated support histories；
- 哪些属于 derived fields；
- 哪些字段必须在 local stage后满足 exact constitutive relation；
- 哪些字段必须在 unrelaxed global stage后满足 exact global relation；
- 哪些 relation在 relaxation后仍保持；
- 哪些 relation在 relaxation后允许 temporary manifold mismatch；
- Eq. (77) exactly读取哪些字段；
- future `LatinStateTower` 的 shape 与 memory layout。

---

# 141. 下一阶段的建议目标

下一阶段建议专门建立：

$$ \boxed{\text{tower LATIN state field contract}} $$

包括：

1. `time`
2. `plastic_strain_rate`
3. `elastic_strain`
4. `alpha_rate`
5. `r_bar_rate`
6. `damage_rate`
7. `stress`
8. `beta`
9. `R_bar`
10. `energy_release_rate`
11. `plastic_strain`
12. `alpha`
13. `r_bar`
14. `damage`

并明确这些数组在 tower 中统一采用：

$$ (N_t,N_q) $$

还是保留：

$$ (N_t,N_e,N_g,N_f) $$

的结构化 storage。

---

# 142. future state contract 还要明确的 derived structural data

除了 material-point state，还需要区分：

- PGD basis data；
- beam nodal displacement correction；
- section strain/resultants；
- material-point compatible total strain；
- damage residual source；
- full mechanical residual；
- search directions；
- outer convergence histories。

这些不应全部塞进 `LatinStateTower`。

---

# 143. future module boundary

理论上已经可以预见的模块职责：

```text
tower_state.py
    primary/support LATIN histories

tower_equilibrium_operator.py
    H, M, C0
    source strain -> compatible strain / equilibrated stress

tower_pgd_time_update.py
    Eq. (58)-(59)

tower_pgd_enrichment.py
    Eq. (61)-(72)
    fixed point
    Gram-Schmidt
    acceptance

tower_global_finishing.py
    mechanical assembly
    Eq. (73)-(75)
    Y

tower_iteration_control.py
    relaxation
    Eq. (76)-(77)
    Eq. (60)
    stagnation

tower_latin_pgd_solver.py
    outer state machine
```

本阶段还不开始创建这些代码。

---

# 144. 本阶段最终冻结的 damage-history decision

这是本阶段最重要的新结论之一：

$$ \boxed{\breve{\dot D}_{i+1}=\hat{\dot D}_{i+1/2}} $$

以及：

$$ \boxed{\breve D_{i+1}=\hat D_{i+1/2}} $$

其中第二式属于：

$$ \boxed{\text{tower-v1 discrete interpretation / refinement}} $$

不是原论文显式写出的独立公式。

---

# 145. 本阶段最终冻结的 $Y$ decision

unrelaxed candidate：

$$ \boxed{\breve Y_{i+1}=Y(\breve\sigma_{i+1},\breve D_{i+1})} $$

relaxation：

$$ \boxed{Y_{i+1}=(1-\mu)Y_i+\mu\breve Y_{i+1}} $$

relaxation后：

$$ \boxed{\text{不额外强制 }Y_{i+1}=Y(\sigma_{i+1},D_{i+1})} $$

让下一 local stage处理 nonlinear manifold projection。

---

# 146. 本阶段最终冻结的 relaxation decision

formal primary fields：

$$ \boxed{q_{i+1}=(1-\mu)q_i+\mu\breve q_{i+1}} $$

integrated support histories：

$$ \boxed{x_{i+1}=(1-\mu)x_i+\mu\breve x_{i+1}} $$

provisional baseline：

$$ \boxed{\mu=0.8} $$

与原论文一致。

---

# 147. 本阶段最终冻结的 Eq. (76) decision

Eq. (76) 使用：

$$ \boxed{\hat s_{i+1/2}} $$

与：

$$ \boxed{s_{i+1}\text{ relaxed}} $$

而不是：

$$ \breve s_{i+1} $$

所以：

$$ \boxed{\text{relaxation precedes LATIN convergence evaluation}} $$

---

# 148. 本阶段最终冻结的 one-pair control decision

tower v1 baseline：

$$ \boxed{\text{one LATIN iteration最多接受一个 new PGD pair}} $$

如果 fixed-basis trial不足：

- add one pair；
- rebuild complete enriched global candidate；
- relax；
- compute indicator；
- accept current LATIN iterate；
- 下一 LATIN iteration再判断是否继续 enrichment。

---

# 149. 本阶段最终冻结的 same-iteration baseline decision

same LATIN iteration内所有 trial：

$$ \boxed{\text{relaxation baseline始终是 }s_i} $$

以及：

$$ \boxed{\text{Eq. (60) baseline始终是 }\xi_i} $$

不能把 fixed-basis trial变成 enriched trial的新 outer baseline。

---

# 150. 本阶段最终冻结的 current 1D compatibility boundary

可直接继承：

- hardening BE equations；
- relaxation architecture；
- Eq. (76) after relaxation；
- primary/support state storage concept。

不直接继承：

- global BE re-integration of local damage history；
- multiple accepted pairs per same LATIN iteration。

---

# 151. 本阶段最重要的 conceptual distinction

以后必须同时保持以下四个层级。

第一：

$$ \boxed{\text{constitutive local history}} $$

由 local stage nonlinear integration得到。

第二：

$$ \boxed{\text{globally admissible mechanical correction}} $$

由 reference equilibrium operator + PGD完成。

第三：

$$ \boxed{\text{unrelaxed global candidate}} $$

满足 current global-stage linear/admissibility equations。

第四：

$$ \boxed{\text{relaxed LATIN iterate}} $$

作为下一 alternating iteration起点，不要求仍处于 nonlinear local manifold。

---

# 152. 为什么这一 distinction 对 tower 特别重要

tower中 material-point数量远大于1D bar。

如果 future code不明确：

- local history；
- global candidate；
- relaxed iterate；
- PGD basis；

各自 ownership，极易出现：

- trial state原地覆盖 accepted state；
- enrichment失败后 rollback不完整；
- damage history被重复积分；
- same iteration重复 local stage；
- support histories与 primary rates不一致；
- $\xi$ 在错误 state上计算。

因此本阶段实际上完成的是 solver state machine 的理论奠基。

---

# 153. 本阶段建议的 future transaction semantics

每个 LATIN iteration建议采用 transaction-like结构：

```text
accepted inputs:
    s_i
    B_m
    xi_i

derive immutable local data:
    hat{s}
    H

build trial:
    breve{s}^{up}
    s^{up}
    xi^{up}

if enrichment:
    copy basis snapshot
    build candidate pair

    if candidate fails:
        rollback basis
        keep fixed-basis trial or terminate according to policy

    if candidate passes:
        build breve{s}^{enr}
        build s^{enr}
        compute xi^{enr}
        commit B_{m+1}, s_{i+1}

end iteration
```

这与前一阶段的 mode-level rollback可以自然衔接。

---

# 154. mode-level rollback 与 LATIN state rollback 的区别

mode-level rollback：

$$ \mathcal B_{m+1}^{\rm tentative}\rightarrow\mathcal B_m $$

发生在 candidate basis acceptance失败时。

LATIN state trial不是 accepted state，直到 iteration commit。

因此：

$$ s_i $$

在整个 current trial/enrichment过程中应保持 immutable。

只有 current iteration最终接受后：

$$ s_i\leftarrow s_{i+1} $$

---

# 155. why immutable previous global state matters

如果 current trial过程中原地修改：

$$ s_i $$

则：

- relaxation baseline被污染；
- damage residual source可能变化；
- fixed/enriched trial不可公平比较；
- rollback变复杂；
- Eq. (60) baseline语义失真。

因此 tower v1 state API应默认：

$$ \boxed{\text{copy / new-state construction rather than in-place mutation of accepted }s_i} $$

---

# 156. 本阶段对代码实现的直接要求

未来 tower global finishing代码至少应返回：

```text
unrelaxed_candidate_state
mechanical_correction diagnostics
hardening diagnostics
damage inheritance diagnostics
energy consistency diagnostics
```

iteration-control再单独返回：

```text
relaxed_state
xi
zeta
stagnation diagnostics
```

这样 global physics与 outer control不混在同一个函数中。

---

# 157. 本阶段不建议现在写代码的原因

虽然主链已经贯通，但 state field contract尚未最终冻结。

如果现在开始实现：

- array shape；
- storage ownership；
- primary/support distinction；
- trial/accepted state semantics；

很可能还要返工。

所以最合理顺序是：

```text
本阶段总结
    ↓
state field-by-field contract
    ↓
module interface specification
    ↓
then code
```

---

# 158. 本阶段核心公式总览

plastic candidate：

$$ \breve\varepsilon^p=\varepsilon_i^p+\Delta\varepsilon^p $$

$$ \breve{\dot\varepsilon}^p=\dot\varepsilon_i^p+\Delta\dot\varepsilon^p $$

$$ \breve\sigma=\sigma_i+\Delta\sigma'+\Delta\tilde\sigma $$

$$ \breve\varepsilon^e=\varepsilon_i^e+\Delta\varepsilon'-\Delta\varepsilon^p+\Delta\tilde\varepsilon $$

hardening：

$$ \beta=C\alpha $$

$$ \bar R=R_\infty\bar r $$

$$ \dot\alpha+H_\beta C\alpha=\hat{\dot\alpha}+H_\beta\hat\beta $$

$$ \dot{\bar r}+H_{\bar R}R_\infty\bar r=\hat{\dot{\bar r}}+H_{\bar R}\hat{\bar R} $$

damage：

$$ \breve{\dot D}=\hat{\dot D} $$

$$ \breve D=\hat D $$

energy：

$$ \breve Y=Y(\breve\sigma,\breve D) $$

relaxation：

$$ s_{i+1}=\mu\breve s_{i+1}+(1-\mu)s_i $$

LATIN indicator：

$$ \xi_{i+1}=\frac{\|\hat s^p_{i+1/2}-s^p_{i+1}\|}{\|\hat s^p_{i+1/2}\|+\|s^p_{i+1}\|} $$

saturation：

$$ \zeta_i=\frac{\xi_i-\xi_{i+1}}{\xi_i+\xi_{i+1}} $$

---

# 159. 本阶段核心数据依赖总览

$$ \boxed{\hat s,H\rightarrow\text{hardening candidate}} $$

$$ \boxed{\hat D,\hat{\dot D}\rightarrow\text{damage candidate}} $$

$$ \boxed{\text{PGD plastic branch}+\text{damage structural projection}\rightarrow\breve\sigma,\breve\varepsilon^e} $$

$$ \boxed{\breve\sigma+\breve D\rightarrow\breve Y} $$

$$ \boxed{s_i+\breve s\rightarrow s_{i+1}\text{ by relaxation}} $$

$$ \boxed{\hat s+s_{i+1}\rightarrow\xi_{i+1}} $$

$$ \boxed{\xi_i+\xi_{i+1}\rightarrow\zeta_i} $$

---

# 160. 本阶段准确研究停点

截至本阶段结束，当前理论停点为：

$$ \boxed{\text{complete relaxed tower LATIN iterate }s_{i+1}\text{ and outer indicators }\xi_{i+1},\zeta_i} $$

并且已经明确：

- accepted new PGD pair如何进入 complete global candidate；
- hardening如何 finishing；
- damage如何 inheritance；
- $Y$何时计算；
- relaxation如何作用；
- Eq. (76)使用哪个 state；
- Eq. (60)使用哪个 baseline；
- one-pair-per-iteration如何嵌入 outer loop。

---

# 161. 下一阶段准确任务

下一阶段不再继续推导 Eq. (73)–(76)。

下一阶段应专门回答：

> future tower `LatinState` 到底保存什么、每个 field 的 shape 是什么、谁负责更新、在哪个 stage 必须满足哪些 exact relation、哪些只是 diagnostic、Eq. (77) 从哪些字段读取。

重点包括：

- formal primary LATIN fields；
- integrated support histories；
- derived constitutive fields；
- structural correction fields；
- PGD basis是否属于 state；
- structured material-point indexing；
- flat indexing；
- memory and reconstruction；
- local/global/trial/accepted consistency contracts；
- unit-test contracts。

---

# 162. 本阶段结论压缩版

本阶段最重要的结论可以压缩为以下十二条。

第一，accepted PGD pair后必须先完成 complete unrelaxed global candidate，不能直接进入 Eq. (76)。

第二，mechanical candidate由 plastic PGD correction与 damage full projection共同组装。

第三，Eq. (73)–(74) hardening继续在每个 tower fiber material point上使用 linear BE solve，不需要 PGD。

第四，$t_0$ hardening state由 initial history固定，initial rate由 descent relation计算。

第五，原论文 Eq. (75) 显式给出 $\dot D_{\rm global}=\dot D_{\rm local}$。

第六，tower v1进一步采用 $\breve D=\hat D$，避免 current 1D 的 RK4-local / BE-global damage re-integration mismatch。

第七，direct-copy integrated $D$ 不取消 damage structural residual projection。

第八，unrelaxed $Y$ 必须由 final unrelaxed stress与 copied damage history重新计算。

第九，formal primary state与 integrated support histories均采用同一个 $\mu=0.8$ convex relaxation。

第十，relaxation后不额外把 $Y$ 强制 reprojection 到 nonlinear state law。

第十一，Eq. (76) 使用 relaxed global state，而不是 unrelaxed global candidate。

第十二，tower v1 一个 LATIN iteration最多接受一个 new PGD pair；same-iteration all trials共享相同 $s_i$ 与 $\xi_i$ baseline。

---

# 163. 最终 algorithm skeleton

```text
INITIAL
    accepted s_i
    accepted basis B_m
    accepted xi_i

LOCAL STAGE
    s_i
        ->
    hat{s}_{i+1/2}

SEARCH DIRECTIONS
    ->
    H_sigma, H_beta, H_Rbar

GLOBAL DAMAGE STRUCTURAL PROJECTION
    ->
    Delta eps_tilde
    Delta sigma_tilde

FIXED-BASIS PGD UPDATE
    Eq. (58)-(59)

PLASTIC RECONSTRUCTION
    ->
    Delta eps_p
    Delta eps_p_dot
    Delta eps_prime
    Delta sigma_prime

MECHANICAL ASSEMBLY
    ->
    breve sigma
    breve eps_e
    breve eps_p
    breve eps_p_dot

HARDENING
    Eq. (73)-(74)

DAMAGE
    breve D_dot = hat D_dot
    breve D = hat D

ENERGY
    breve Y = Y(breve sigma, breve D)

COMPLETE UNRELAXED TRIAL
    breve{s}_{i+1}^{up}

RELAX
    from s_i
    ->
    s_{i+1}^{up}

LATIN INDICATOR
    ->
    xi_{i+1}^{up}

SATURATION
    ->
    zeta_i^{up}

IF BASIS SUFFICIENT
    commit s_{i+1}^{up}

ELSE
    freeze current local/global data

    ADD ONE PAIR
        Eq. (61)-(72)
        fixed point
        Gram-Schmidt
        significance
        all-mode temporal re-optimisation
        residual benefit

    IF PAIR REJECTED
        rollback basis

    IF PAIR ACCEPTED
        B_m -> B_{m+1}

        reconstruct enriched mechanical correction

        reuse hardening candidate
        reuse D and D_dot
        recompute Y

        ->
        breve{s}_{i+1}^{enr}

        relax from same s_i

        ->
        s_{i+1}^{enr}

        compute xi_{i+1}^{enr}
        compute zeta_i^{enr}

        commit one enriched LATIN iterate

NEXT LATIN ITERATION
    new local stage
```

---

# 164. 本阶段结束语

经过本阶段，tower LATIN-PGD 已经不再只具备 PGD enrichment 的局部公式，而是已经形成从：

$$ \boxed{\text{local nonlinear projection}} $$

到：

$$ \boxed{\text{global reduced/full correction}} $$

再到：

$$ \boxed{\text{hardening / damage finishing}} $$

再到：

$$ \boxed{\text{relaxation}} $$

最后到：

$$ \boxed{\text{Eq. (76)–(77) convergence and Eq. (60) saturation}} $$

的完整 outer algorithm backbone。

下一步只剩 state/data contract 与 implementation interface 的最后冻结。

在这一点完成之前，仍不建议开始 tower LATIN-PGD 主求解器代码，以避免 state ownership 与 array semantics 在实现过程中反复返工。
