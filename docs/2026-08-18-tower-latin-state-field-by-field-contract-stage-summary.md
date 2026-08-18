# 2026-08-18 Tower LATIN State Field-by-Field Contract Stage Summary

- 项目：Offshore Wind Turbine and LATIN-PGD
- 仓库：`Jiarui-Shi-DUT/Offshore_Wind_Turbine_LATIN_PGD`
- 分支：`feature/offshore-wind-turbine-tower-fatigue`
- 日期：2026-08-18
- 阶段性质：theory–data-contract–interface freeze
- 当前状态：future tower `LatinStateTower` 的 13 个 material-point fields 已完成 field-by-field conceptual contract；尚未开始 tower LATIN-PGD solver 代码修改
- 上一阶段文档：`docs/2026-08-18-tower-latin-state-storage-field-classification-stage-summary.md`
- 本阶段范围：在 canonical material-point storage、immutable layout 与 14-entry 第一层分类已冻结的基础上，继续逐字段闭合 mechanical、hardening、damage-energy 三组 state contract
- 下一阶段：四层 state ownership、transaction semantics、trial/rollback/commit 与 module interfaces 的最终冻结；完成后进入代码编写、单元测试与调试
- PyCharm Preview 兼容性：本修正版除继续禁止数学环境中的 raw `<` 外，所有 display-math 块均压缩为单物理行 `$$ ... $$`，避免独立 `+` / `-` 行被 Markdown 误解析为列表标记

---

# 1. 本阶段的准确起点

上一阶段已经冻结：

```text
canonical material-point coordinate
    q

q <-> (element, Gauss point, fiber)

canonical space-time field shape
    (Nt, Nq)

MaterialPointLayout
    immutable

LatinStateTower
    separate from MaterialPointLayout

PGDBasis
    separate from LatinStateTower
```

并把 future tower state 的 14 个 candidate entries 分为：

```text
grid metadata
    time

formal primary LATIN fields
    plastic_strain_rate
    elastic_strain
    alpha_rate
    r_bar_rate
    damage_rate
    stress
    beta
    R_bar
    energy_release_rate

integrated support histories
    plastic_strain
    alpha
    r_bar
    damage
```

因此，本阶段不再讨论 storage 方案，而是回答：

> 这些 fields 在 accepted global state、local half-step、unrelaxed global candidate 和 relaxed global state 中分别由谁读取、谁更新、满足什么 exact relation、哪些 relation 在 relaxation 后继续成立、哪些允许暂时不成立，以及 Eq. (77) 如何读取这些 fields。

---

# 2. 本阶段继续采用的四层判断来源

后续所有结论继续严格区分：

## 2.1 原论文明确内容

包括：

- formal LATIN state；
- local/global alternating structure；
- ascent/descent search directions；
- Eq. (73)–(75) hardening/damage global finishing；
- relaxation；
- Eq. (76)–(77) convergence indicator 与 mechanical norm。

## 2.2 严格数学推导

包括：

- 线性 state relation 在相同 convex relaxation 下的保持性；
- 线性 time-difference operator 与 convex relaxation 的可交换性；
- 相同 damage rate 与相同 initial damage 在 continuous setting 下对应相同 integrated damage history；
- 线性 equilibrium / compatibility 在 convex relaxation 下的保持性。

## 2.3 current 1D implementation

包括：

- `latin/state.py`
- `latin/local_stage.py`
- `latin/global_stage.py`
- `latin/iteration_control.py`
- `material/viscoplastic_damage_1d.py`

用于继承已经验证的一维数值经验，但不能把实现选择反写成“论文明确规定”。

## 2.4 tower v1 engineering choice

本阶段新增并冻结的 tower v1 选择包括：

- global plastic correction 的 accumulated/rate consistency 继续采用当前 original x-t 路线下的 BE/right-endpoint convention；
- hardening global finishing 继续采用 BE；
- damage global candidate 直接继承 local RK4 damage history，不再进行第二次 BE reintegration；
- relaxed global state 不对 nonlinear energy-release relation做额外 reprojection；
- accepted baseline state在 current LATIN iteration 内 immutable；
- trial state采用 build-new-state / commit-later 语义。

---

# 3. 四个 LATIN state layer

本阶段所有 field contract 均围绕四个 state layer 定义。

上一 accepted relaxed global state：

$$ s_i $$

current local half-step state：

$$ \hat s_{i+1/2} $$

complete unrelaxed global candidate：

$$ \breve s_{i+1} $$

current relaxed global state：

$$ s_{i+1} $$

其基本关系：

```text
accepted s_i
    ↓
local projection
    ↓
hat{s}_{i+1/2}
    ↓
global PGD + global finishing
    ↓
breve{s}_{i+1}
    ↓
relaxation
    ↓
s_{i+1}
```

其中：

> `s_i` 一旦 current LATIN iteration 开始，就必须作为 immutable baseline。

---

# 4. Material-point canonical shape

除 `time` 外，本阶段所有 fields 统一采用：

$$ (N_t,N_q). $$

其中：

$$ q\leftrightarrow(e,g,f). $$

`time` 继续采用：

$$ (N_t,). $$

本阶段不改变上一阶段已经冻结的数据结构。

---

# 5. Mechanical field group

本阶段首先闭合：

$$ \varepsilon^p,\qquad \dot\varepsilon^p,\qquad \varepsilon^e,\qquad \sigma. $$

其分类为：

```text
plastic_strain
    integrated support history

plastic_strain_rate
    formal primary LATIN field

elastic_strain
    formal primary LATIN field

stress
    formal primary LATIN field
```

---

# 6. 必须区分 physical initial time 与 LATIN initial history

两个“初始”概念不能混淆。

physical initial time：

$$ t_0 $$

描述结构真实时间历史起点。

LATIN initial globally admissible history：

$$ s_0(t,q) $$

则是定义在完整 time-space grid 上的 initial global guess。

对于 virgin-state baseline，tower v1 继续采用：

$$ \varepsilon^p_0(t,q)=0, $$

$$ \dot\varepsilon^p_0(t,q)=0. $$

而：

$$ \sigma_0(t,q) $$

和：

$$ \varepsilon^e_0(t,q) $$

来自 tower elastic initialization。

这意味着 LATIN initialization 是完整时间域上的 globally admissible history，而不是只初始化单个 $t_0$ state。

---

# 7. `plastic_strain` 的 physical meaning

定义：

$$ \varepsilon^p(t,q) $$

为 material point $q$ 的 accumulated plastic strain history。

其职责：

- 为下一 local stage提供 integrated plastic initial/history information；
- 与 plastic-rate field维持 global discrete consistency；
- 承载 PGD plastic correction 的 accumulated form；
- 支持 trial / relaxation / rollback；
- 支持 material-state verification。

它不直接进入 Eq. (77)。

---

# 8. `plastic_strain` 在 accepted state 中的 ownership

上一 accepted state 保存：

$$ \varepsilon_i^p(t,q). $$

current LATIN iteration 开始后：

> `epsilon_p_i` 必须 immutable。

任何 fixed-basis trial、enriched trial、candidate finishing 都不能原地改写它。

---

# 9. Local stage 对 `plastic_strain` 的处理

local stage从 current accepted global state 的 physical initial internal state出发，在 fixed thermodynamic-force histories 下重新进行 constitutive history integration。

因此：

```text
s_i at t0
    ↓
local constitutive integration over complete physical time interval
    ↓
hat epsilon_p
```

current 1D local implementation与：

$$ \alpha,\quad \bar r,\quad D $$

一起对：

$$ \varepsilon^p $$

执行 RK4 integration。

因此 tower v1冻结：

> `hat plastic_strain` 是 local constitutive history，不是 global accepted history的简单复制。

---

# 10. Global candidate 对 `plastic_strain` 的处理

PGD mechanical branch求解：

$$ \Delta\varepsilon^p. $$

因此 unrelaxed global candidate：

$$ \breve\varepsilon^p_{i+1} = \varepsilon_i^p + \Delta\varepsilon^p_{i+1}. $$

必须明确：

$$ \breve\varepsilon^p \neq \hat\varepsilon^p $$

一般成立。

local history与global candidate plastic history之间的差异是 LATIN alternating manifold gap 的一部分。

---

# 11. `plastic_strain` relaxation

support history使用与 formal state相同的 relaxation：

$$ \varepsilon^p_{i+1} = (1-\mu)\varepsilon_i^p + \mu\breve\varepsilon^p_{i+1}. $$

代入 global correction：

$$ \varepsilon^p_{i+1} = \varepsilon_i^p + \mu\Delta\varepsilon^p_{i+1}. $$

因此，在 mechanical branch 上，relaxation 等价于对 global plastic correction 乘以 $\mu$。

---

# 12. `plastic_strain_rate` 的 physical meaning

定义：

$$ \dot\varepsilon^p(t,q) $$

为 material-point plastic strain rate。

它是：

> formal primary LATIN field。

它直接进入 Eq. (77)。

---

# 13. Local stage 中 `plastic_strain_rate`

local stage必须满足 current nonlinear viscoplastic flow law。

因此：

$$ \hat{\dot\varepsilon}^p = \mathcal F_p ( \hat\sigma, \hat\beta, \hat{\bar R}, \hat Y, \hat D ). $$

这里：

> `hat plastic_strain_rate` 由 local constitutive law拥有。

它不是 PGD reconstruction，也不是简单由 nodal plastic history backward difference 得到。

---

# 14. Global candidate 中 `plastic_strain_rate`

PGD branch同时得到：

$$ \Delta\varepsilon^p $$

与：

$$ \Delta\dot\varepsilon^p. $$

因此：

$$ \breve{\dot\varepsilon}^p_{i+1} = \dot\varepsilon_i^p + \Delta\dot\varepsilon^p_{i+1}. $$

不能令：

$$ \breve{\dot\varepsilon}^p = \hat{\dot\varepsilon}^p. $$

否则 global PGD correction 将被破坏。

---

# 15. Global plastic rate-history discrete consistency

tower v1冻结：

> accumulated plastic correction 与 plastic-rate correction 使用相同 time-difference convention。

对：

$$ n\ge1, $$

采用：

$$ \Delta\dot\varepsilon^p_n = \frac{ \Delta\varepsilon^p_n - \Delta\varepsilon^p_{n-1} }{ \Delta t_n }. $$

这里继续采用当前 original x-t 路线中已经形成的 BE/right-endpoint convention。

此处属于：

> tower v1 discrete consistency choice。

不是“原论文规定 Python 必须使用某个 array difference”。

---

# 16. Plastic rate-history consistency 与 relaxation

设 previous accepted global state满足：

$$ \dot\varepsilon_i^p = D_t\varepsilon_i^p, $$

candidate满足：

$$ \breve{\dot\varepsilon}^p = D_t\breve\varepsilon^p. $$

若 $D_t$ 为同一线性 difference operator，则：

$$ D_t \left[ (1-\mu)\varepsilon_i^p + \mu\breve\varepsilon^p \right] = (1-\mu)\dot\varepsilon_i^p + \mu\breve{\dot\varepsilon}^p. $$

因此相同 convex relaxation 不破坏 global plastic rate-history consistency。

---

# 17. `plastic_strain_rate` relaxation

定义：

$$ \dot\varepsilon^p_{i+1} = (1-\mu)\dot\varepsilon_i^p + \mu\breve{\dot\varepsilon}^p_{i+1}. $$

因此：

$$ \dot\varepsilon^p_{i+1} = \dot\varepsilon_i^p + \mu\Delta\dot\varepsilon^p. $$

---

# 18. `stress` 的 physical meaning

定义：

$$ \sigma(t,q) $$

为 material-point stress field。

它是：

> formal primary LATIN field。

它直接进入 Eq. (77)。

在 global stage 中，它也是 global equilibrium admissibility 的核心 field。

---

# 19. Local stage 中 stress ownership

current ascent-direction choice保持 thermodynamic-force histories fixed。

因此：

$$ \hat\sigma_{i+1/2} = \sigma_i. $$

local stage：

```text
reads sigma_i
does not solve structural equilibrium
copies sigma_i into hat sigma
```

---

# 20. Plastic branch stress correction

accepted PGD plastic correction：

$$ \Delta\varepsilon^p $$

经过 tower reference equilibrium operator：

$$ \Delta\varepsilon' = \mathcal E_{\rm tower} \Delta\varepsilon^p $$

以及：

$$ \Delta\sigma' = C_0 (\mathcal E_{\rm tower}-I) \Delta\varepsilon^p. $$

其中：

$$ \mathcal E_{\rm tower} = H(H^TMC_0H)^{-1}H^TMC_0. $$

---

# 21. Damage structural branch stress correction

damage-dependent residual source：

$$ \Delta\varepsilon^R $$

经过同一 reference equilibrium operator：

$$ \Delta\tilde\varepsilon = \mathcal E_{\rm tower} \Delta\varepsilon^R, $$

$$ \Delta\tilde\sigma = C_0 (\mathcal E_{\rm tower}-I) \Delta\varepsilon^R. $$

---

# 22. Unrelaxed global stress assembly

因此：

$$ \breve\sigma_{i+1} = \sigma_i + \Delta\sigma' + \Delta\tilde\sigma. $$

必须明确：

> plastic PGD branch 与 damage full structural projection branch 都进入 final stress candidate。

---

# 23. Global equilibrium preservation

reference equilibrium operator保证 correction满足：

$$ H^TM\Delta\sigma'=0, $$

以及：

$$ H^TM\Delta\tilde\sigma=0. $$

如果：

$$ \sigma_i $$

已经满足 current prescribed loading history的 global equilibrium，则：

$$ \breve\sigma_{i+1} $$

仍满足相同 equilibrium。

因此 stress candidate属于 global admissibility manifold的一部分。

---

# 24. Stress relaxation 后 equilibrium 仍精确保持

stress relaxation：

$$ \sigma_{i+1} = (1-\mu)\sigma_i + \mu\breve\sigma_{i+1}. $$

由于 equilibrium equation在线性结构空间中是线性的，而 previous state 与 candidate 都满足同一 external loading equilibrium，因此：

> relaxed global stress 继续精确满足 equilibrium。

这是 future unit test的重要 invariant。

---

# 25. `elastic_strain` 的 physical meaning

定义：

$$ \varepsilon^e(t,q) $$

为 material-point elastic strain。

它是：

> formal primary LATIN field。

它直接进入 Eq. (77)。

---

# 26. Local stage 中 `elastic_strain`

local stage在 fixed stress history下更新 damage history。

因此必须通过 nonlinear damaged elastic state law重新计算：

$$ \hat\varepsilon^e = \mathcal G_e ( \hat\sigma, \hat D ). $$

对于 current unilateral model：

- tension 与 compression 使用不同 degraded stiffness factor；
- 该 relation属于 local constitutive manifold的 exact relation。

---

# 27. Global candidate 中 `elastic_strain`

global stage不重新执行 nonlinear elastic projection。

而是通过 compatible structural corrections组装：

$$ \breve\varepsilon^e_{i+1} = \varepsilon_i^e + \Delta\varepsilon' - \Delta\varepsilon^p + \Delta\tilde\varepsilon. $$

---

# 28. 为什么必须减去 plastic correction

总应变满足：

$$ \varepsilon^{\rm total} = \varepsilon^e+\varepsilon^p. $$

new plastic history：

$$ \breve\varepsilon^p = \varepsilon_i^p+\Delta\varepsilon^p. $$

compatible total-strain correction为：

$$ \Delta\varepsilon' + \Delta\tilde\varepsilon. $$

因此：

$$ \breve\varepsilon^{\rm total} = \varepsilon_i^{\rm total} + \Delta\varepsilon' + \Delta\tilde\varepsilon. $$

要保持这一关系，必须有：

$$ \breve\varepsilon^e = \varepsilon_i^e + \Delta\varepsilon' + \Delta\tilde\varepsilon - \Delta\varepsilon^p. $$

所以负号来自 additive strain split，而不是 numerical trick。

---

# 29. Global stage 对 elastic strain 的 exact relation

global candidate要求的是：

> total-strain compatibility。

即：

$$ \breve\varepsilon^e + \breve\varepsilon^p = \varepsilon_i^e + \varepsilon_i^p + \Delta\varepsilon' + \Delta\tilde\varepsilon. $$

如果 previous total strain compatible，而两项 structural corrections来自 tower compatibility operator，则 new total strain继续 compatible。

---

# 30. Global candidate 不要求重新满足 local damaged elastic law

一般：

$$ \breve\varepsilon^e \neq \mathcal G_e ( \breve\sigma, \breve D ). $$

这是允许的。

global stage应满足：

> global admissibility。

而不是把 candidate再投影回 local constitutive manifold。

---

# 31. Elastic-strain relaxation 后 compatibility 仍保持

relaxation：

$$ \varepsilon^e_{i+1} = (1-\mu)\varepsilon_i^e + \mu\breve\varepsilon^e. $$

同时：

$$ \varepsilon^p_{i+1} = (1-\mu)\varepsilon_i^p + \mu\breve\varepsilon^p. $$

因此：

$$ \varepsilon^e_{i+1} + \varepsilon^p_{i+1} $$

是两个 compatible total-strain histories 的 convex combination。

由于 compatibility relation是线性的：

> relaxation后 total-strain compatibility精确保持。

但 nonlinear damaged elastic law一般仍允许不成立。

---

# 32. Mechanical group 的 LATIN manifold distinction

本阶段冻结：

```text
local state
    exact constitutive relations

global candidate
    exact global equilibrium
    exact global compatibility

relaxed global state
    equilibrium preserved
    compatibility preserved

not required after global stage
    nonlinear local flow law
    nonlinear damaged elastic law
```

这是 mechanical-field contract的核心。

---

# 33. Mechanical fields 与 Eq. (77)

Eq. (77)对本组读取：

```text
plastic_strain
    no direct participation

plastic_strain_rate
    direct participation

elastic_strain
    direct participation

stress
    direct participation
```

因此：

$$ \varepsilon^p $$

虽然不直接进入 Eq. (77)，但不可删除。

---

# 34. Fixed-basis Trial A 与 enriched Trial B

same current LATIN iteration 中：

$$ s_i $$

以及：

$$ \hat s_{i+1/2} $$

必须 frozen。

Trial A：

$$ \Delta\varepsilon_A^p, \qquad \Delta\dot\varepsilon_A^p. $$

Trial B：

$$ \Delta\varepsilon_B^p, \qquad \Delta\dot\varepsilon_B^p. $$

每个 trial都必须从同一个 immutable：

$$ s_i $$

重新构建 mechanical candidate。

禁止：

```text
Trial A
    ↓
mutate
    ↓
Trial B
```

正确语义：

```text
immutable s_i
    ├── Trial A
    └── Trial B
```

---

# 35. Mechanical group rollback semantics

candidate被拒绝时：

```text
candidate basis
    discard or rollback

candidate mechanical fields
    discard

accepted s_i
    untouched
```

因此 future implementation必须遵循：

> build-new-state / commit-later。

不采用：

> modify-accepted-state / undo-if-failed。

---

# 36. Mechanical group future invariants

未来至少检查：

- `hat sigma == sigma_i`；
- local flow-law consistency；
- local damaged elastic-law consistency；
- candidate plastic correction assembly；
- candidate plastic-rate correction assembly；
- candidate stress assembly；
- candidate total-strain compatibility；
- equilibrium preservation；
- relaxation-preserved equilibrium；
- relaxation-preserved compatibility；
- accepted baseline immutability；
- Trial A/B common-baseline semantics。

---

# 37. Hardening field group

第二组闭合：

$$ \alpha,\quad \dot\alpha,\quad \beta,\quad \bar r,\quad \dot{\bar r},\quad \bar R. $$

分为：

```text
kinematic hardening
    alpha
    alpha_rate
    beta

transformed isotropic hardening
    r_bar
    r_bar_rate
    R_bar
```

---

# 38. Hardening fields 的分类

```text
alpha
    integrated support history

alpha_rate
    formal primary LATIN field

beta
    formal primary thermodynamic-force field

r_bar
    integrated support history

r_bar_rate
    formal primary LATIN field

R_bar
    formal primary thermodynamic-force field
```

全部 shape：

$$ (N_t,N_q). $$

---

# 39. Global linear hardening state laws

current material model中：

$$ \beta=C\alpha, $$

以及：

$$ \bar R=R_\infty\bar r. $$

这些是 global hardening state-law relations。

---

# 40. 为什么 beta 与 R_bar 不能删除

即使：

$$ \beta $$

可由：

$$ \alpha $$

计算，

以及：

$$ \bar R $$

可由：

$$ \bar r $$

计算，它们仍然是 formal LATIN force coordinates。

它们直接参与：

- search-direction relation；
- Eq. (73)–(74) global hardening finishing；
- Eq. (77) mechanical norm。

所以它们必须作为真正存储的 primary fields。

---

# 41. Local stage 不要求满足 global linear hardening state laws

current ascent choice固定：

$$ \hat\beta=\beta_i, $$

以及：

$$ \hat{\bar R}=\bar R_i. $$

但 local stage重新积分：

$$ \hat\alpha, \qquad \hat{\bar r}. $$

因此一般：

$$ \hat\beta \neq C\hat\alpha, $$

以及：

$$ \hat{\bar R} \neq R_\infty\hat{\bar r}. $$

这不是错误。

---

# 42. Local kinematic hardening evolution

current local material law可概念写为：

$$ \hat{\dot\alpha} = \dot\lambda \left( n_p - \frac{a}{C}\hat\beta \right). $$

其中 $n_p$ 为 plastic flow direction。

因此：

> local `alpha_rate` 由 constitutive evolution law拥有。

---

# 43. Local transformed isotropic hardening evolution

current local law：

$$ \hat{\dot{\bar r}} = \dot\lambda \left( \sqrt{\gamma} - \frac{\gamma\hat{\bar R}} {2R_\infty} \right). $$

因此：

> local `r_bar_rate` 由 constitutive evolution law拥有。

---

# 44. Local integrated hardening histories

local stage沿完整 physical time interval积分：

$$ \hat\alpha $$

与：

$$ \hat{\bar r}. $$

current 1D implementation使用 RK4。

因此 local hardening contract是：

```text
beta_i, R_bar_i
    fixed as local forces

local evolution laws
    ↓
hat alpha_rate
hat r_bar_rate
    ↓
RK4
    ↓
hat alpha
hat r_bar
```

---

# 45. Eq. (74) kinematic branch

global hardening descent relation：

$$ \dot\alpha + H_\beta\beta = \hat{\dot\alpha} + H_\beta\hat\beta. $$

代入：

$$ \beta=C\alpha, $$

得到：

$$ \dot\alpha + H_\beta C\alpha = \hat{\dot\alpha} + H_\beta\hat\beta. $$

这是每个 material point上独立的 linear global-finishing ODE。

---

# 46. Eq. (74) transformed isotropic branch

同理：

$$ \dot{\bar r} + H_{\bar R}\bar R = \hat{\dot{\bar r}} + H_{\bar R}\hat{\bar R}. $$

利用：

$$ \bar R = R_\infty\bar r, $$

得到：

$$ \dot{\bar r} + H_{\bar R}R_\infty\bar r = \hat{\dot{\bar r}} + H_{\bar R}\hat{\bar R}. $$

---

# 47. Tower v1 hardening time discretisation

tower v1继续采用 current mature 1D global finishing中的 backward Euler。

对：

$$ n\ge1, $$

有：

$$ \dot\alpha_n = \frac{ \alpha_n-\alpha_{n-1} }{ \Delta t_n }. $$

---

# 48. Kinematic hardening BE formula

代入 global hardening equation：

$$ \breve\alpha_n = \frac{ \breve\alpha_{n-1} + \Delta t_n \left( \hat{\dot\alpha}_n + H_{\beta,n}\hat\beta_n \right) }{ 1+ \Delta t_nH_{\beta,n}C }. $$

随后：

$$ \breve{\dot\alpha}_n = \frac{ \breve\alpha_n-\breve\alpha_{n-1} }{ \Delta t_n }, $$

以及：

$$ \breve\beta_n = C\breve\alpha_n. $$

---

# 49. Isotropic hardening BE formula

同理：

$$ \breve{\bar r}_n = \frac{ \breve{\bar r}_{n-1} + \Delta t_n \left( \hat{\dot{\bar r}}_n + H_{\bar R,n}\hat{\bar R}_n \right) }{ 1+ \Delta t_nH_{\bar R,n}R_\infty }. $$

随后：

$$ \breve{\dot{\bar r}}_n = \frac{ \breve{\bar r}_n-\breve{\bar r}_{n-1} }{ \Delta t_n }, $$

以及：

$$ \breve{\bar R}_n = R_\infty\breve{\bar r}_n. $$

---

# 50. Hardening denominator positivity

若：

$$ H_\beta\gt0, $$

$$ H_{\bar R}\gt0, $$

并且：

$$ C\gt0, $$

$$ R_\infty\gt0, $$

则：

$$ 1+\Delta t_nH_{\beta,n}C\gt1, $$

以及：

$$ 1+\Delta t_nH_{\bar R,n}R_\infty\gt1. $$

因此 hardening scalar recursion denominator 不会接近零。

---

# 51. Hardening physical initial time

$t_0$ 不存在 previous physical time node，因此不能使用 backward difference。

tower v1冻结：

$$ \breve\alpha_0 = \hat\alpha_0, $$

$$ \breve{\bar r}_0 = \hat{\bar r}_0. $$

随后 enforcing global linear state laws：

$$ \breve\beta_0 = C\breve\alpha_0, $$

$$ \breve{\bar R}_0 = R_\infty\breve{\bar r}_0. $$

---

# 52. Hardening t0 rates

$t_0$ rate通过 descent relation直接得到：

$$ \breve{\dot\alpha}_0 = \hat{\dot\alpha}_0 - H_{\beta,0} ( \breve\beta_0-\hat\beta_0 ), $$

以及：

$$ \breve{\dot{\bar r}}_0 = \hat{\dot{\bar r}}_0 - H_{\bar R,0} ( \breve{\bar R}_0-\hat{\bar R}_0 ). $$

---

# 53. Virgin initial hardening state

current baseline从 virgin steel开始，因此：

$$ \alpha_0=0, $$

$$ \bar r_0=0, $$

进而：

$$ \beta_0=0, $$

$$ \bar R_0=0. $$

future contract保留一般性：

> 若未来采用 pre-hardened initial state，可允许非零 histories，但 accepted global initial state必须满足 linear state laws。

---

# 54. Unrelaxed global hardening exact relations

complete unrelaxed candidate必须满足：

$$ \breve\beta = C\breve\alpha, $$

$$ \breve{\bar R} = R_\infty\breve{\bar r}. $$

对：

$$ n\ge1, $$

还必须满足：

$$ \breve{\dot\alpha}_n = \frac{ \breve\alpha_n-\breve\alpha_{n-1} }{ \Delta t_n }, $$

$$ \breve{\dot{\bar r}}_n = \frac{ \breve{\bar r}_n-\breve{\bar r}_{n-1} }{ \Delta t_n }. $$

并满足 current Eq. (74) descent residual到 numerical precision。

---

# 55. Relaxation 后 linear hardening state law 的保持性

若 previous accepted state：

$$ \beta_i=C\alpha_i, $$

candidate：

$$ \breve\beta=C\breve\alpha, $$

并对两者使用同一 $\mu$：

$$ \alpha_{i+1} = (1-\mu)\alpha_i + \mu\breve\alpha, $$

$$ \beta_{i+1} = (1-\mu)\beta_i + \mu\breve\beta, $$

则：

$$ \beta_{i+1} = C\alpha_{i+1}. $$

因此 Eq. (73) kinematic hardening state law 在 relaxation 后精确保持。

---

# 56. Relaxation 后 isotropic hardening state law 的保持性

完全同理：

$$ \bar R_i = R_\infty\bar r_i, $$

$$ \breve{\bar R} = R_\infty\breve{\bar r} $$

可推出：

$$ \bar R_{i+1} = R_\infty\bar r_{i+1}. $$

因此 transformed isotropic hardening linear state law 也在 relaxation 后精确保持。

---

# 57. Hardening rate-history consistency 与 relaxation

若 previous state与candidate均满足同一线性 BE operator：

$$ \dot\alpha=D_t\alpha, $$

则：

$$ \dot\alpha_{i+1} = D_t\alpha_{i+1}. $$

同理：

$$ \dot{\bar r}_{i+1} = D_t\bar r_{i+1}. $$

因此 global accepted state可以保持：

```text
beta = C alpha
R_bar = R_inf r_bar
alpha_rate = D_t alpha
r_bar_rate = D_t r_bar
```

其中 rate/history relation针对 $n\ge1$。

---

# 58. Eq. (74) descent relation在 relaxation 后不要求继续成立

Eq. (74) 描述的是 current：

$$ \hat s_{i+1/2} $$

与：

$$ \breve s_{i+1} $$

之间的 descent relation。

relaxation后：

$$ s_{i+1} = (1-\mu)s_i + \mu\breve s_{i+1}. $$

previous：

$$ s_i $$

一般不满足“相对于 current local state”的同一 Eq. (74)。

因此：

> relaxed `s_(i+1)` 不要求 current Eq. (74) 继续精确成立。

---

# 59. Local hardening state 不应使用 BE consistency test

local stage使用 RK4。

因此不能强制：

$$ \hat{\dot\alpha}_n = \frac{ \hat\alpha_n-\hat\alpha_{n-1} }{ \Delta t_n } $$

作为 exact identity。

同理也不能强制：

$$ \hat{\dot{\bar r}}_n = \frac{ \hat{\bar r}_n-\hat{\bar r}_{n-1} }{ \Delta t_n }. $$

future tests必须区分：

```text
local state
    RK4 constitutive consistency

global candidate
    BE hardening consistency
```

---

# 60. Fixed-basis Trial A 与 enriched Trial B 的 hardening invariance

same current LATIN iteration中以下 quantities frozen：

$$ \hat{\dot\alpha}, \quad \hat\beta, \quad \hat{\dot{\bar r}}, \quad \hat{\bar R}, $$

以及：

$$ H_\beta, \quad H_{\bar R}. $$

PGD enrichment改变 mechanical plastic branch，但 Eq. (74) hardening finishing不依赖新的 plastic mode本身。

因此：

> Trial A 与 Trial B 的 unrelaxed hardening candidate相同。

---

# 61. Hardening fields 仍必须进入 complete trial state

即使 same-iteration可以复用 hardening candidate，每个 complete trial state仍必须包含：

$$ \dot\alpha, \quad \beta, \quad \dot{\bar r}, \quad \bar R, $$

因为 Eq. (77) 直接读取这些 fields。

---

# 62. Hardening group 与 Eq. (77)

```text
alpha
    no direct participation

alpha_rate
    direct participation

beta
    direct participation

r_bar
    no direct participation

r_bar_rate
    direct participation

R_bar
    direct participation
```

---

# 63. Hardening group ownership

accepted `s_i`：

```text
alpha_i
alpha_rate_i
beta_i
r_bar_i
r_bar_rate_i
R_bar_i
    immutable during current LATIN iteration
```

local stage：

```text
reads / fixes
    beta_i
    R_bar_i

computes
    hat alpha_rate
    hat r_bar_rate

integrates
    hat alpha
    hat r_bar

copies
    hat beta = beta_i
    hat R_bar = R_bar_i
```

global finishing：

```text
reads
    hat alpha_rate
    hat beta
    H_beta
    hat r_bar_rate
    hat R_bar
    H_R_bar

writes
    breve alpha
    breve alpha_rate
    breve beta
    breve r_bar
    breve r_bar_rate
    breve R_bar
```

---

# 64. Hardening group future invariants

future unit tests至少包括：

- local beta force-copy；
- local R_bar force-copy；
- local constitutive-rate correctness；
- Eq. (73) kinematic residual；
- Eq. (73) transformed isotropic residual；
- Eq. (74) kinematic descent residual；
- Eq. (74) isotropic descent residual；
- global BE consistency；
- relaxation-preserved linear state laws；
- Trial A/B hardening invariance；
- accepted baseline immutability。

---

# 65. Damage-energy field group

第三组闭合：

$$ D,\qquad \dot D,\qquad Y. $$

分类：

```text
damage
    integrated support history

damage_rate
    formal primary LATIN field

energy_release_rate
    formal primary thermodynamic-force field
```

全部 shape：

$$ (N_t,N_q). $$

---

# 66. Damage energy-release nonlinear state law

current unilateral model中，tension branch：

$$ Y = \frac{ \sigma^2 }{ 2E(1-D)^2 }. $$

compression branch：

$$ Y = \frac{ h\sigma^2 }{ 2E(1-hD)^2 }. $$

这是 nonlinear relation。

它与 hardening linear state laws的 relaxation性质不同。

---

# 67. Elastic initialization 中 damage fields

current virgin baseline采用：

$$ D_0(t,q)=0, $$

$$ \dot D_0(t,q)=0. $$

initial：

$$ Y_0(t,q) $$

由 elastic stress history 与 $D=0$ 计算。

注意：

> 这是 LATIN initial globally admissible guess，不意味着完整 time history 已经满足 nonlinear damage evolution。

---

# 68. Local stage 中 Y ownership

current ascent direction保持 thermodynamic-force history：

$$ \hat Y = Y_i. $$

因此 local stage不重新求 global energy force。

与：

$$ \hat\sigma=\sigma_i $$

同属 fixed-force history contract。

---

# 69. Local damage-rate evolution

current damage evolution law概念上：

$$ \hat{\dot D} = \mathcal F_D(\hat Y). $$

在当前材料模型中可写为：

$$ \hat{\dot D} = k_D \left[ \text{positive part of } ( \hat Y-Y_0 ) \right]^{n_D}. $$

因此：

> local `damage_rate` 由 damage evolution law拥有。

---

# 70. Local damage history

local stage使用 RK4对：

$$ D $$

沿完整 physical time interval积分。

所以：

```text
Y_i fixed
    ↓
local damage evolution law
    ↓
hat D_rate
    ↓
RK4
    ↓
hat D
```

---

# 71. Local state 一般不要求 Y = Y(sigma,D)

由于 local：

$$ \hat Y=Y_i $$

保持 fixed，

同时：

$$ \hat D $$

发生更新，

以及：

$$ \hat\sigma=\sigma_i, $$

一般：

$$ \hat Y \neq Y(\hat\sigma,\hat D). $$

这不是 local-stage error。

因此 future validator不得强制：

```text
hat Y == Y(hat sigma, hat D)
```

---

# 72. Local elastic relation 与 local Y relation必须区分

local mechanical field中：

$$ \hat\varepsilon^e = \mathcal G_e(\hat\sigma,\hat D) $$

是 exact constitutive relation。

但：

$$ \hat Y = Y(\hat\sigma,\hat D) $$

不是 current local-stage exact relation。

二者不能混淆。

---

# 73. Eq. (75) 的直接理论内容

current：

$$ b^-=0 $$

时，Eq. (75) 直接给出 damage-rate inheritance：

$$ \breve{\dot D} = \hat{\dot D}. $$

原论文直接约束的是：

> damage rate。

并没有逐字给出 tower discrete array上的：

$$ \breve D=\hat D. $$

---

# 74. Current 1D damage global finishing

current 1D implementation：

```text
local RK4
    ↓
hat D
hat D_rate
    ↓
global
    copy hat D_rate
    ↓
BE re-integrate
    ↓
new global D
```

因此一般：

$$ D_{\rm local}^{RK4} \neq D_{\rm global}^{BE}. $$

这是 current 1D implementation choice，而不是必须保留到 tower v1 的理论要求。

---

# 75. Continuous interpretation of Eq. (75)

若 continuous setting 中：

$$ \breve{\dot D}(t) = \hat{\dot D}(t), $$

且：

$$ \breve D(t_0) = \hat D(t_0), $$

则：

$$ \breve D(t) = D_0 + \int_{t_0}^{t} \breve{\dot D}(\tau)\,d\tau, $$

和：

$$ \hat D(t) = D_0 + \int_{t_0}^{t} \hat{\dot D}(\tau)\,d\tau $$

给出相同 continuous history。

因此：

$$ \breve D = \hat D $$

是由 common initial condition + identical continuous rate 推导出的自然结果。

---

# 76. Tower v1 damage inheritance frozen choice

tower v1正式冻结：

$$ \breve{\dot D}_{i+1} = \hat{\dot D}_{i+1/2}, $$

以及：

$$ \breve D_{i+1} = \hat D_{i+1/2}. $$

也就是：

```text
local RK4
    ↓
hat D_rate
hat D
    ↓
global finishing
    ↓
direct inheritance
```

不再执行第二次 BE damage reintegration。

---

# 77. Tower v1 damage change 的来源标签

必须明确：

$$ \breve{\dot D} = \hat{\dot D} $$

属于：

> paper explicit Eq. (75) consequence under current $b^-=0$ choice。

而：

$$ \breve D = \hat D $$

属于：

> strict continuous derivation + tower v1 discrete engineering choice。

不能把后一条写成“论文显式给出”。

---

# 78. 为什么 direct damage-history inheritance 更干净

主要收益：

- 避免 same LATIN half-step 内同时存在 RK4 local damage history 与 BE global damage history；
- 避免额外的 RK4-vs-BE discretisation defect；
- 更直接体现 $b^-=0$ 下 global damage branch不重新求解 damage evolution；
- 简化 trial invariance 与 regression tests。

---

# 79. Direct D inheritance 不等于 damage structural correction消失

即使：

$$ \breve D=\hat D, $$

仍然必须保留 damage-dependent structural residual projection：

$$ \Delta\varepsilon^R \rightarrow \Delta\tilde\varepsilon, \quad \Delta\tilde\sigma. $$

其中：

- $D$ 是 constitutive internal history；
- $\Delta\varepsilon^R$ branch 是 damaged nonlinear elasticity相对于 reference operator 的 structural admissibility correction。

因此不能推出：

$$ \Delta\tilde\sigma=0. $$

---

# 80. Global candidate 中 Y 的计算时机

必须先完成 final mechanical stress：

$$ \breve\sigma = \sigma_i + \Delta\sigma' + \Delta\tilde\sigma. $$

然后继承：

$$ \breve D=\hat D. $$

最后计算：

$$ \breve Y = Y(\breve\sigma,\breve D). $$

因此 finishing order：

```text
final plastic mechanical correction
    ↓
damage structural projection
    ↓
final breve sigma
    ↓
inherit breve D_rate
inherit breve D
    ↓
compute breve Y
```

---

# 81. Unrelaxed Y tension branch

material point $q$ 若：

$$ \breve\sigma_q\ge0, $$

则：

$$ \breve Y_q = \frac{ \breve\sigma_q^2 }{ 2E_q(1-\breve D_q)^2 }. $$

---

# 82. Unrelaxed Y compression branch

material point $q$ 若：

$$ \breve\sigma_q\lt0, $$

则：

$$ \breve Y_q = \frac{ h_q\breve\sigma_q^2 }{ 2E_q(1-h_q\breve D_q)^2 }. $$

这条 nonlinear state relation必须在 complete unrelaxed candidate上精确成立。

---

# 83. 为什么不能直接令 breve Y = hat Y

local：

$$ \hat Y=Y_i $$

是 fixed force history。

global mechanical finishing改变：

$$ \sigma_i \rightarrow \breve\sigma. $$

因此必须重新计算：

$$ \breve Y = Y(\breve\sigma,\breve D). $$

不能复制：

$$ \breve Y=\hat Y. $$

---

# 84. Local 与 global 对 Y 的 owner 不同

local stage：

> $Y$ 是 fixed thermodynamic-force input。

unrelaxed global finishing：

> $Y$ 是 nonlinear state-law output。

因此同一 field 在两个 half-stage 中有不同 owner semantics。

---

# 85. Damage history relaxation

定义：

$$ D_{i+1} = (1-\mu)D_i + \mu\breve D. $$

由于：

$$ \breve D=\hat D, $$

得到：

$$ D_{i+1} = (1-\mu)D_i + \mu\hat D. $$

---

# 86. Damage-rate relaxation

定义：

$$ \dot D_{i+1} = (1-\mu)\dot D_i + \mu\breve{\dot D}. $$

由于：

$$ \breve{\dot D} = \hat{\dot D}, $$

得到：

$$ \dot D_{i+1} = (1-\mu)\dot D_i + \mu\hat{\dot D}. $$

---

# 87. Damage group 不强制 BE identity

tower v1 candidate直接继承 local RK4 history 与 constitutive nodal rates。

因此不能要求：

$$ \breve{\dot D}_n = \frac{ \breve D_n-\breve D_{n-1} }{ \Delta t_n }. $$

否则等于重新要求 RK4 history满足 BE endpoint identity。

future damage validator应检查：

$$ \breve D=\hat D, $$

$$ \breve{\dot D}=\hat{\dot D}, $$

而不是 BE relation。

---

# 88. Relaxed damage state 也不强制 BE identity

同样，accepted relaxed：

$$ D_{i+1} $$

与：

$$ \dot D_{i+1} $$

是两个 LATIN states 的 convex blend。

在 discrete RK4 storage语义下，不强制：

$$ \dot D_{i+1,n} = \frac{ D_{i+1,n}-D_{i+1,n-1} }{ \Delta t_n }. $$

---

# 89. Energy-release-rate relaxation

formal LATIN relaxation：

$$ Y_{i+1} = (1-\mu)Y_i + \mu\breve Y. $$

tower v1继续使用这一 field-wise convex blend。

---

# 90. Relaxation 后不重新投影 Y

完成：

$$ \sigma_{i+1}, $$

$$ D_{i+1}, $$

$$ Y_{i+1} $$

的 field-wise relaxation后，不执行：

$$ Y_{i+1} \leftarrow Y( \sigma_{i+1}, D_{i+1} ). $$

这是当前 tower v1 已冻结的重要原则。

---

# 91. 为什么不能 relaxation 后 reproject Y

因为：

$$ Y=Y(\sigma,D) $$

是 nonlinear constitutive relation。

若 relaxation 后额外重新计算：

$$ Y $$

就会把：

```text
global candidate
    ↓
LATIN relaxation
    ↓
relaxed global state
```

改成：

```text
global candidate
    ↓
LATIN relaxation
    ↓
extra nonlinear constitutive projection
    ↓
another state
```

这会改变 LATIN relaxation定义的 state。

---

# 92. Allowed Y manifold mismatch

因此 relaxed global state一般允许：

$$ Y_{i+1} \neq Y( \sigma_{i+1}, D_{i+1} ). $$

这是：

> allowed temporary manifold mismatch。

不是 numerical inconsistency。

---

# 93. 为什么 Y 必须真正存储

如果 `energy_release_rate` 只是动态 property：

```text
Y = function(sigma, D)
```

那么任何读取都会自动把它重新投影到 nonlinear constitutive relation。

这样就无法表示 relaxed LATIN state 中允许的：

$$ Y_{i+1} \neq Y( \sigma_{i+1}, D_{i+1} ). $$

因此 `energy_release_rate` 必须是 stored formal primary field。

---

# 94. Trial A/B 中 D 与 D_rate 不变

same current LATIN iteration 中：

$$ \hat D, \qquad \hat{\dot D} $$

frozen。

tower v1又定义：

$$ \breve D=\hat D, $$

$$ \breve{\dot D}=\hat{\dot D}. $$

因此：

$$ \breve D_A = \breve D_B, $$

以及：

$$ \breve{\dot D}_A = \breve{\dot D}_B. $$

所以 same-iteration PGD enrichment不改变 unrelaxed damage history。

---

# 95. Trial A/B 中 Y 一般不同

Trial A与Trial B的 damage history相同，但 mechanical PGD correction不同，因此：

$$ \breve\sigma_A \neq \breve\sigma_B $$

一般成立。

所以：

$$ \breve Y_A = Y( \breve\sigma_A, \hat D ), $$

而：

$$ \breve Y_B = Y( \breve\sigma_B, \hat D ). $$

因此 Trial B 必须重新计算 $Y$。

---

# 96. Same-iteration damage-energy dependency graph

```text
frozen local damage data
    hat D
    hat D_rate
        │
        ├── Trial A
        │       ↓
        │   breve sigma_A
        │       ↓
        │   breve Y_A
        │
        └── Trial B
                ↓
            breve sigma_B
                ↓
            breve Y_B
```

因此：

```text
D
D_rate
    same-iteration invariant

Y
    trial-dependent output
```

---

# 97. Damage bounds

current material model要求：

$$ 0\le D\lt D_{\max}. $$

current default：

$$ D_{\max}=0.999. $$

tower v1继续要求：

- damage finite；
- damage位于 material admissible interval；
- local RK4负责 protection；
- unrelaxed candidate直接继承 local admissibility；
- convex relaxation不产生超出两端 convex hull 的新 overshoot。

---

# 98. Damage monotonicity

current damage law满足：

$$ \dot D\ge0. $$

因此 local integrated damage history应 non-decreasing。

若：

- previous accepted global history non-decreasing；
- local candidate history non-decreasing；

则相同 $\mu$ 的 convex blend仍然 non-decreasing。

因此 accepted tower damage history可将 monotonicity作为 current material model下的 secondary invariant。

但：

> 不能将该结论泛化成所有未来 damage models 的 universal LATIN rule。

---

# 99. Damage fields 与 Eq. (77)

current Eq. (77) mechanical norm不直接读取：

$$ D, \qquad \dot D, \qquad Y. $$

但它们通过：

- damaged elasticity；
- plastic flow；
- search-direction evaluation；
- stress；
- elastic strain；

间接影响 mechanical convergence。

未来若增加 damage-specific diagnostic，必须明确标记为：

> secondary implementation diagnostic。

不能把它写成 paper Eq. (77) 的一部分。

---

# 100. Damage-energy group ownership

accepted `s_i`：

```text
D_i
D_rate_i
Y_i
    immutable during current LATIN iteration
```

local stage：

```text
reads
    D_i at physical initial time
    Y_i over complete time history

fixes
    hat Y = Y_i

computes
    hat D_rate

integrates
    hat D by RK4
```

global finishing：

```text
reads
    hat D
    hat D_rate
    final breve sigma

writes
    breve D = hat D
    breve D_rate = hat D_rate
    breve Y = Y(breve sigma, breve D)
```

relaxation：

```text
D_(i+1)
    blend(D_i, breve D)

D_rate_(i+1)
    blend(D_rate_i, breve D_rate)

Y_(i+1)
    blend(Y_i, breve Y)

no nonlinear Y reprojection
```

---

# 101. Damage-energy future invariants

future unit tests至少包括：

- local Y force-copy；
- local damage-rate constitutive correctness；
- local RK4 damage integration；
- direct unrelaxed damage-rate inheritance；
- direct unrelaxed damage-history inheritance；
- no second BE damage reintegration；
- unrelaxed Y nonlinear state-law consistency；
- relaxed Y no-reprojection semantics；
- Trial A/B D invariance；
- Trial A/B D_rate invariance；
- Trial A/B Y recomputation；
- damage bounds；
- accepted baseline immutability。

---

# 102. 13 个 material-point fields 的完整分类再次确认

Mechanical group：

```text
plastic_strain
    support history

plastic_strain_rate
    primary

elastic_strain
    primary

stress
    primary
```

Hardening group：

```text
alpha
    support history

alpha_rate
    primary

beta
    primary

r_bar
    support history

r_bar_rate
    primary

R_bar
    primary
```

Damage-energy group：

```text
damage
    support history

damage_rate
    primary

energy_release_rate
    primary
```

---

# 103. Local-state exact relations 总览

current local half-step：

$$ \hat s_{i+1/2} $$

必须满足 / 使用：

```text
hat sigma = sigma_i

hat beta = beta_i

hat R_bar = R_bar_i

hat Y = Y_i

hat plastic_strain_rate
    from local viscoplastic evolution law

hat alpha_rate
    from local hardening evolution law

hat r_bar_rate
    from local hardening evolution law

hat damage_rate
    from local damage evolution law

hat plastic_strain
hat alpha
hat r_bar
hat damage
    from local RK4 integrated histories

hat elastic_strain
    from nonlinear damaged elastic state law
```

---

# 104. Local-state relations that are NOT required

local half-step一般不要求：

$$ \hat\beta = C\hat\alpha, $$

$$ \hat{\bar R} = R_\infty\hat{\bar r}, $$

$$ \hat Y = Y( \hat\sigma, \hat D ). $$

local RK4 state也不要求：

$$ \hat{\dot\alpha}=D_t^{BE}\hat\alpha, $$

$$ \hat{\dot{\bar r}}=D_t^{BE}\hat{\bar r}, $$

$$ \hat{\dot D}=D_t^{BE}\hat D. $$

这些“不要求”关系必须进入 future validator设计，避免错误测试。

---

# 105. Unrelaxed global candidate exact relations 总览

complete：

$$ \breve s_{i+1} $$

必须满足：

Mechanical：

$$ \breve\varepsilon^p = \varepsilon_i^p + \Delta\varepsilon^p, $$

$$ \breve{\dot\varepsilon}^p = \dot\varepsilon_i^p + \Delta\dot\varepsilon^p, $$

$$ \breve\sigma = \sigma_i + \Delta\sigma' + \Delta\tilde\sigma, $$

$$ \breve\varepsilon^e = \varepsilon_i^e + \Delta\varepsilon' - \Delta\varepsilon^p + \Delta\tilde\varepsilon. $$

Hardening：

$$ \breve\beta = C\breve\alpha, $$

$$ \breve{\bar R} = R_\infty\breve{\bar r}, $$

并满足 Eq. (74) 与 global BE recursion。

Damage：

$$ \breve{\dot D} = \hat{\dot D}, $$

$$ \breve D = \hat D. $$

Energy：

$$ \breve Y = Y( \breve\sigma, \breve D ). $$

---

# 106. Unrelaxed global candidate 的 global admissibility

complete candidate必须保持：

> stress equilibrium。

以及：

> total-strain compatibility。

但一般不要求：

> nonlinear plastic flow law。

也不要求：

> nonlinear damaged elastic law。

这正是 global manifold语义。

---

# 107. Relaxed global state exact relations 总览

定义：

$$ s_{i+1} = (1-\mu)s_i + \mu\breve s_{i+1}. $$

relaxation 后继续精确保持：

- global equilibrium；
- total-strain compatibility；
- linear hardening state law；
- global plastic rate-history linear consistency；
- global hardening BE rate-history consistency。

---

# 108. Relaxed global state 允许的 temporary mismatch

relaxation后不要求：

- current local plastic flow law；
- current local hardening evolution law；
- current Eq. (74) descent relation；
- nonlinear damaged elastic law；
- nonlinear energy-release relation。

特别是：

$$ Y_{i+1} \neq Y( \sigma_{i+1}, D_{i+1} ) $$

一般允许。

---

# 109. Eq. (77) 的 exact field reader

Eq. (77) mechanical norm直接读取：

$$ \sigma, $$

$$ \beta, $$

$$ \bar R, $$

$$ \dot\varepsilon^p, $$

$$ \varepsilon^e, $$

$$ \dot\alpha, $$

$$ \dot{\bar r}. $$

不直接读取：

$$ \varepsilon^p, $$

$$ \alpha, $$

$$ \bar r, $$

$$ D, $$

$$ \dot D, $$

$$ Y. $$

但后者仍通过 constitutive coupling影响前七个 mechanical fields。

---

# 110. Full state contract 与 Eq. (77) 的关系

判断 field 是否保存的依据不是：

> 是否进入 Eq. (77)。

而是：

> 它是否属于 formal LATIN coordinate、integrated numerical support history，或者下一阶段 local/global update所必须的 state information。

因此不进入 Eq. (77) 的 histories不能删除。

---

# 111. Trial A / Trial B 的 complete-field semantics

same current LATIN iteration：

```text
immutable
    s_i
    hat{s}
    H_sigma
    H_beta
    H_R_bar
    damage residual source

Trial A
    fixed-basis plastic correction
    mechanical finishing
    hardening finishing
    damage inheritance
    Y_A
    relaxation from s_i

Trial B
    enriched plastic correction
    mechanical finishing
    same hardening candidate
    same D and D_rate candidate
    recomputed Y_B
    relaxation from same s_i
```

---

# 112. Same-iteration quantities changed by enrichment

accepted new PGD pair主要改变：

$$ \Delta\varepsilon^p, $$

$$ \Delta\dot\varepsilon^p, $$

$$ \Delta\sigma', $$

$$ \Delta\varepsilon'. $$

进而改变：

$$ \breve\varepsilon^p, $$

$$ \breve{\dot\varepsilon}^p, $$

$$ \breve\sigma, $$

$$ \breve\varepsilon^e, $$

以及：

$$ \breve Y. $$

---

# 113. Same-iteration quantities unchanged by enrichment

在 frozen local state / search directions条件下，unrelaxed hardening candidate相同：

$$ \breve\alpha, \quad \breve{\dot\alpha}, \quad \breve\beta, $$

$$ \breve{\bar r}, \quad \breve{\dot{\bar r}}, \quad \breve{\bar R}. $$

damage candidate相同：

$$ \breve D, \qquad \breve{\dot D}. $$

但：

$$ \breve Y $$

必须随 trial stress重新计算。

---

# 114. State immutability principle

current LATIN iteration中：

$$ s_i $$

必须始终 immutable。

candidate creation、enrichment、mode rejection、full trial evaluation都不能对其进行 in-place mutation。

这是后续 transaction semantics 的核心前提。

---

# 115. Mode rollback 与 state rollback 必须区分

future implementation中至少存在两个不同 rollback层次。

Mode-level rollback：

```text
candidate PGD pair rejected
    ↓
basis restored
```

LATIN-state trial discard：

```text
complete trial not committed
    ↓
trial state discarded
```

二者不能混为同一个 `revert()`。

---

# 116. 本阶段已提前形成的 unit-test philosophy

future tests必须分层。

Local-manifold tests：

- fixed force histories；
- constitutive rates；
- RK4 integrated histories；
- local nonlinear elastic relation。

Global-unrelaxed tests：

- PGD mechanical correction assembly；
- equilibrium；
- compatibility；
- Eq. (73)；
- Eq. (74)；
- direct damage inheritance；
- unrelaxed $Y$ nonlinear relation。

Relaxed-global tests：

- exact convex blending；
- equilibrium preservation；
- compatibility preservation；
- linear hardening-state preservation；
- no post-relaxation Y reprojection。

Transaction tests：

- accepted baseline immutability；
- same-baseline Trial A/B；
- mode rollback不污染 state；
- trial discard不污染 accepted state。

---

# 117. 本阶段对 current 1D implementation 的继承内容

明确继承：

- 14-entry state naming philosophy；
- time-by-space state arrays；
- local full-history integration philosophy；
- fixed thermodynamic-force local-stage semantics；
- mechanical plastic + damage structural correction split；
- hardening BE global finishing；
- field-wise LATIN relaxation；
- Eq. (77) seven-field mechanical norm philosophy；
- accepted baseline 与 candidate state分离的方向。

---

# 118. 本阶段对 current 1D implementation 的明确修改

tower v1最明确的 state-contract修改是 damage history。

current 1D：

```text
copy local D_rate
    ↓
BE reintegrate global D
```

tower v1：

```text
copy local D_rate
copy local RK4 D
```

即：

$$ \breve{\dot D} = \hat{\dot D}, $$

$$ \breve D = \hat D. $$

这个改变必须在未来 regression tests 中锁定。

---

# 119. 本阶段仍然没有引入的新方法

本阶段没有引入：

- cycle-phase decomposition；
- multi-time-scale PGD；
- cycle jump；
- random wind-wave ROM；
- online SVD；
- new damage convergence norm；
- extra post-relaxation constitutive projection；
- new scalar line search；
- multi-pair enrichment per LATIN iteration。

当前路线仍然保持：

> original x-t LATIN-PGD first。

---

# 120. Future `LatinStateTower` conceptual field contract 已基本闭合

截至本阶段结束：

```text
time
    grid metadata

mechanical group
    plastic_strain
    plastic_strain_rate
    elastic_strain
    stress

hardening group
    alpha
    alpha_rate
    beta
    r_bar
    r_bar_rate
    R_bar

damage-energy group
    damage
    damage_rate
    energy_release_rate
```

均已明确：

- physical meaning；
- field class；
- canonical shape；
- local-stage semantics；
- global-stage semantics；
- relaxation semantics；
- exact relations；
- allowed mismatch；
- Eq. (77) participation；
- trial invariance / dependence；
- rollback implications；
- future unit-test invariants。

---

# 121. 还没有正式冻结的最后理论接口

虽然 field-by-field conceptual contract 已基本闭合，但正式写代码前还必须完成最后一个 interface-level阶段。

需要一次性冻结：

- 四层 state object ownership；
- immutable vs mutable arrays；
- `MaterialPointLayout` 的实际 interface；
- `LatinStateTower` 的 constructor / validation semantics；
- state cloning / readonly semantics；
- `PGDBasis` snapshot / candidate / commit interface；
- fixed-basis Trial A container；
- enriched Trial B container；
- mechanical correction result container；
- hardening finishing result；
- damage-energy finishing result；
- relaxation interface；
- Eq. (77) reader interface；
- mode-level rollback；
- LATIN iteration-level commit；
- solver orchestration order；
- module import dependencies；
- minimal unit-test boundary。

---

# 122. 下一阶段建议的 module boundary

上一阶段已提出，下一阶段应正式冻结以下候选 module responsibilities：

```text
tower_state.py

tower_equilibrium_operator.py

tower_pgd_time_update.py

tower_pgd_enrichment.py

tower_global_finishing.py

tower_iteration_control.py

tower_latin_pgd_solver.py
```

当前仍然只是 interface specification。

还没有开始创建这些文件。

---

# 123. 下一阶段准确目标

下一阶段不再继续推 material constitutive equations。

目标是：

> 四层 state ownership + transaction semantics + module interface final freeze。

核心要回答：

```text
who owns s_i?

who builds hat{s}?

who owns breve{s}?

who relaxes to s_(i+1)?

who owns basis snapshot?

when does basis commit?

when does state commit?

what exactly is rolled back?

which module may mutate which object?

which data are immutable common inputs?

what data may be cached across Trial A / Trial B?

what is the minimal object passed to Eq. (77)?

what are the minimal unit tests before integration?
```

---

# 124. 进入代码阶段的门槛

下一阶段完成后，应达到：

```text
state data contract
    frozen

field ownership
    frozen

trial semantics
    frozen

rollback semantics
    frozen

module I/O
    frozen

unit-test contracts
    frozen
```

一旦这些条件满足：

> 不再继续延长纯理论阶段，正式进入 tower LATIN-PGD 代码编写与调试。

---

# 125. 预计进入代码阶段后的实现顺序

当前建议顺序仍然是：

```text
1. tower_state.py
    ↓
2. tower equilibrium / layout bridge
    ↓
3. tower global finishing
    ↓
4. fixed-basis temporal update
    ↓
5. one-pair enrichment
    ↓
6. iteration control
    ↓
7. tower_latin_pgd_solver.py
    ↓
8. focused unit tests
    ↓
9. minimal tower integration test
    ↓
10. cyclic tower regression / diagnostics
```

具体文件创建仍需等下一阶段 interface freeze结束后再做。

---

# 126. Markdown / PyCharm Preview 规范

后续阶段总结继续继承：

`docs/2026-08-17-markdown-pycharm-math-preview-debugging-lessons.md`

规则：

1. 数学环境中不使用 raw less-than sign；
2. 小于关系使用 `\lt`；
3. 大于关系优先使用 `\gt`；
4. 独立数学公式使用 `$$ ... $$`；
5. 行内数学使用 `$...$`；
6. 不用复杂 LaTeX `array` 包裹中文；
7. 流程使用 Markdown code block；
8. 解释文字放在数学环境外；
9. Preview异常优先定位第一处失败公式；
10. 使用 minimal reproduction 与 binary split 排错；
11. 一次只改变一个变量。

本文档按该规则生成。

---

# 127. 本阶段最终研究停点

截至 2026-08-18 本阶段结束：

```text
canonical storage
    frozen

14-entry first-level classification
    frozen

mechanical field contract
    frozen

hardening field contract
    frozen

damage-energy field contract
    frozen
```

因此 current exact stop 为：

> 进入四层 state ownership、transaction semantics、trial/rollback/commit 与 module-interface final freeze。

该 final interface freeze结束后：

> 正式开始 tower LATIN-PGD 代码编写、单元测试与调试。
