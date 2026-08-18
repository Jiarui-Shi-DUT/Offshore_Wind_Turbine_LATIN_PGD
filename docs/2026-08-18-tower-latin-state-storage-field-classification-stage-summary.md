# 2026-08-18 Tower LATIN State Storage and Field Classification Stage Summary

- 项目：Offshore Wind Turbine and LATIN-PGD
- 仓库：`Jiarui-Shi-DUT/Offshore_Wind_Turbine_LATIN_PGD`
- 分支：`feature/offshore-wind-turbine-tower-fatigue`
- 日期：2026-08-18
- 阶段性质：理论—数据结构—接口冻结阶段
- 当前状态：本阶段结论已形成，尚未开始 tower LATIN-PGD solver 代码修改
- 上一阶段理论停点：Eq. (47)–(77) 主算法链基本贯通，下一任务为 future tower `LatinState` field contract 与 state/interface semantics
- 本文档目的：完整记录本轮关于 tower material-point canonical storage、material-point indexing、`LatinStateTower` 第一层字段分类及其理论边界的推导、讨论和冻结结论，为下一阶段逐字段 contract 奠定统一基础

---

# 1. 本阶段为什么必须先做数据结构冻结

截至上一阶段，tower LATIN-PGD 已经完成了从 Eq. (47) 到 Eq. (77) 的主理论链闭合，包括：

- tower reference equilibrium operator；
- Eq. (58)–(59) fixed-basis temporal update；
- Eq. (60) saturation；
- Eq. (61)–(72) rank-one enrichment；
- Eq. (70) tower spatial solve；
- Eq. (72) temporal update；
- enrichment fixed-point convergence；
- post-fixed-point weighted Gram-Schmidt；
- temporal coordinate transformation；
- mode significance / full residual acceptance；
- Eq. (73)–(74) hardening finishing；
- Eq. (75) damage treatment；
- complete unrelaxed global state；
- relaxation；
- Eq. (76)–(77) LATIN convergence indicator；
- fixed-basis / enriched trial 与 outer-control 的衔接。

因此，当前问题已经不再是继续推导新的 PGD 方程，而是：

> future tower `LatinState` 到底如何组织 material-point state、每个 field 保存什么、各 field 的 canonical shape 是什么，以及这些 state 如何与 tower FOM、PGD basis、equilibrium operator 和 outer LATIN iteration 对接。

如果这些接口没有先冻结就直接编码，后续极容易出现：

- 同一物理字段同时存在 flat 与 structured 两套 mutable representation；
- fixed-basis trial 和 enriched trial 对 state ownership 的理解不一致；
- accepted state 在 trial 构造过程中被原地污染；
- material-point indexing 与 FOM hierarchy 发生错位；
- PGD basis dimension、Eq. (77) metric dimension 与 state dimension 不统一；
- rollback 只能恢复部分数据；
- support histories 与 formal primary fields 混淆。

所以，本阶段首先冻结的是 solver 数据语义，而不是 solver 实现。

---

# 2. 本阶段遵循的资料与判断层级

后续所有判断继续严格区分四个层次。

## 2.1 原论文明确给出的内容

原论文明确给出：

- formal LATIN state 的组成；
- space-time LATIN alternating structure；
- Eq. (47)–(77) 的核心理论结构；
- PGD spatial / temporal separation；
- Eq. (58)–(59) fixed-basis temporal update；
- Eq. (60) saturation；
- Eq. (61)–(72) rank-one enrichment；
- Eq. (73)–(75) finishing；
- relaxation；
- Eq. (76)–(77) LATIN convergence indicator。

## 2.2 从原论文与 tower 离散严格推导得到的内容

当前已经接受的严格离散对应包括：

- tower 中每个 fiber material point 对应原论文 continuum material point；
- material-point hierarchy 为 `element -> Gauss point -> section fiber -> material point`；
- tower PGD spatial function 离散后应作用于全部 material points；
- tower material-point spatial coordinate 可统一记为 $q$；
- PGD spatial modes、material-point residuals、search directions 和 Eq. (77) spatial quantities应共享同一个 material-point dimension。

## 2.3 current 1D implementation

当前一维实现中：

- `LatinState` 的 time 之外所有 field 采用 `(N_t, N_e)`；
- 其中第二维代表 bar element / material point；
- local stage 使用 RK4；
- temporal update 使用当前实现的 BE-like 离散；
- global state 中显式保存 integrated histories；
- relaxation 同时作用于 formal fields 和 support histories。

这些属于已经验证的一维 Python 实现经验，不能自动写成“论文规定”。

## 2.4 tower v1 engineering choice

本阶段新增并冻结的内容主要属于 tower v1 engineering choice，包括：

- canonical material-point storage 采用 `(N_t, N_q)`；
- 使用 immutable `q <-> (e,g,f)` mapping；
- structured hierarchy 只作为 view / indexing convenience，不作为第二套 mutable state；
- `LatinStateTower`、`MaterialPointLayout`、`PGDBasis` 和 structural/equilibrium data 严格分离；
- 14 个 candidate entries 划分为 grid metadata、formal primary LATIN fields 和 integrated support histories 三类。

这些选择应在后续文档和代码注释中明确标记为 tower v1 design，而不是原论文原文。

---

# 3. Tower material-point hierarchy 的物理定义

当前 tower FOM 已经采用如下 material-point hierarchy：

```text
tower
└── element e
    └── Gauss point g
        └── section fiber f
            └── material point q
```

其中：

- $e$：beam-column element index；
- $g$：element 内 Gauss integration-point index；
- $f$：section fiber index；
- $q$：LATIN-PGD 中统一使用的 flat material-point index。

因此：

$$ q \leftrightarrow (e,g,f) $$

是一个一一对应关系。

本阶段不改变 tower FOM 的物理层级，只是为 LATIN-PGD 建立统一的 spatial coordinate。

---

# 4. 两种 candidate storage 方案

在进入 field contract 前，本阶段比较了两种主要 storage 方案。

## 4.1 方案 A：structured storage

所有 material-point histories 采用：

$$ (N_t,N_e,N_g,N_f) $$

优点：

- 与 tower FOM 对象 hierarchy 直观一致；
- post-processing 时容易按 element / Gauss point / fiber 访问；
- 人工检查局部响应方便。

缺点：

- PGD spatial mode 最终仍需要 flatten；
- Eq. (70) spatial solve 最终仍处理长度为 $N_q$ 的向量；
- Eq. (77) spatial metric 最终仍需要 material-point vector；
- weighted Gram-Schmidt 最终仍需要 material-point vector；
- equilibrium operator 的 source strain / compatible strain / equilibrated stress 最终仍需要统一 material-point dimension；
- repeated reshape / flatten 会增加实现复杂度；
- 若 structured array 与 flat array 同时可写，会产生双重数据 owner 风险。

## 4.2 方案 B：flat material-point storage

所有 material-point histories 采用：

$$ (N_t,N_q) $$

其中：

$$ N_q=\text{tower 中全部 fiber material points 的总数}. $$

element / Gauss / fiber hierarchy 通过单独的 immutable layout 保存。

优点：

- 直接对应原论文离散后的 spatial coordinate；
- PGD spatial mode天然为 $\mathbb R^{N_q}$；
- state、PGD、search directions、equilibrium operator 与 Eq. (77) 使用统一 spatial dimension；
- 避免 flat / structured 两套 mutable representation；
- trial / copy / rollback 只需处理一套 canonical arrays；
- 与 current 1D `(N_t,N_e)` implementation 在概念上保持同构。

本阶段最终选择方案 B。

---

# 5. Frozen decision 1：canonical material-point storage

tower v1 正式冻结：

> LATIN-PGD 内部统一采用 flat material-point coordinate $q$。所有 material-point space-time state fields 使用 `(N_t,N_q)` 作为 canonical storage。

因此 future tower state 中：

```text
stress.shape                  = (Nt, Nq)
elastic_strain.shape          = (Nt, Nq)
plastic_strain_rate.shape     = (Nt, Nq)
alpha_rate.shape              = (Nt, Nq)
r_bar_rate.shape              = (Nt, Nq)
damage_rate.shape             = (Nt, Nq)
beta.shape                    = (Nt, Nq)
R_bar.shape                   = (Nt, Nq)
energy_release_rate.shape     = (Nt, Nq)

plastic_strain.shape          = (Nt, Nq)
alpha.shape                   = (Nt, Nq)
r_bar.shape                   = (Nt, Nq)
damage.shape                  = (Nt, Nq)
```

而：

```text
time.shape = (Nt,)
```

注意：这里的 `(N_t,N_q)` 是 tower v1 canonical representation，不是原论文规定的 Python array shape。

---

# 6. 为什么 `(N_t,N_q)` 与原论文 x-t PGD 结构更一致

原论文 PGD spatial function 可写为：

$$ \bar{\varepsilon}^p(x). $$

在 tower 离散中，continuum material point $x$ 由 fiber material point $q$ 代表。

因此第 $j$ 个 tower plastic spatial mode 可写为：

$$ \bar{\boldsymbol\varepsilon}^{p}_j\in\mathbb R^{N_q}. $$

相应 temporal function 为：

$$ \lambda_j\in\mathbb R^{N_t}. $$

单个 separated mode 为：

$$ \Delta\varepsilon^p_j(t,q)=\lambda_j(t)\bar{\varepsilon}^{p}_j(q). $$

若当前 basis dimension 为 $m$，则：

$$ \Lambda\in\mathbb R^{N_t\times m} $$

以及：

$$ P\in\mathbb R^{N_q\times m}. $$

重构后的 plastic correction 为：

$$ \Delta E^p=\Lambda P^T $$

且：

$$ \Delta E^p\in\mathbb R^{N_t\times N_q}. $$

因此 `(N_t,N_q)` 并不是为了编程方便而任意选择，而是与当前坚持的 original $x-t$ PGD architecture 具有直接的离散同构关系。

---

# 7. 与 current 1D implementation 的关系

当前一维 `LatinState` 的基本形式是：

```text
time:
    (Nt,)

all state fields:
    (Nt, Ne)
```

在 1D bar 中，第二维的 element index 本质上同时承担 material-point spatial index 的角色。

因此 tower v1 的迁移可理解为：

```text
current 1D:
space index = element

tower v1:
space index = fiber material point q
```

也就是说，tower v1 不改变 LATIN-PGD 的 x-t data philosophy，而只是把 spatial discretisation 从 bar elements 扩展到 tower fiber material points。

---

# 8. Frozen decision 2：MaterialPointLayout

为了保留 tower FOM hierarchy，future solver 需要一个独立于 `LatinStateTower` 的 immutable material-point layout。

其理论职责是：

> 唯一管理 flat material-point index $q$ 与 FOM hierarchy `(element, Gauss point, fiber)` 之间的映射。

概念上：

```text
q = 0
    -> element 0
    -> gauss 0
    -> fiber 0

q = 1
    -> element 0
    -> gauss 0
    -> fiber 1

...

q = ...
    -> element e
    -> gauss g
    -> fiber f
```

正式冻结：

$$ q \leftrightarrow (e,g,f) $$

必须是可逆的一一映射。

---

# 9. Material-point ordering

tower v1 推荐并冻结的 canonical ordering 为：

```text
element-major
    ↓
Gauss-point-major
    ↓
fiber-major
```

即最内层先遍历 fiber，之后 Gauss point，再 element。

在所有 element 具有相同 $N_g$ 和 $N_f$ 的简单情形下，可以写出概念性索引：

$$ q=((eN_g)+g)N_f+f. $$

但后续 solver module 不应各自复制这一索引公式。

更稳健的接口原则是：

```text
(e, g, f) -> q
q -> (e, g, f)
```

均统一通过 material-point layout 提供。

这样未来即使：

- element 的 Gauss-point 数发生变化；
- section fiber 数发生变化；
- 某些 element 使用不同离散；

核心 LATIN-PGD module 仍然只认统一的 $q$。

---

# 10. Structured hierarchy 只作为 view

虽然 canonical storage 使用：

$$ (N_t,N_q), $$

但 post-processing 和 debugging 仍需要按：

```text
element
Gauss point
fiber
```

访问历史。

因此 future interface 可以提供 structured access，但其语义只能是：

> view / indexing convenience

而不能成为第二套独立可写的数据。

例如概念上：

```text
q = layout.index(e, g, f)
history = state.stress[:, q]
```

或者未来可能存在：

```text
state.at(e, g, f)
```

但 structured access 必须最终指向 canonical `(N_t,N_q)` arrays。

Frozen principle：

```text
single mutable owner:
    (Nt, Nq)

structured hierarchy:
    read/indexing view
```

---

# 11. 为什么必须坚持 single source of truth

当前 tower LATIN iteration 将至少出现：

```text
accepted state s_i
local state hat{s}
fixed-basis unrelaxed trial
fixed-basis relaxed trial
enriched unrelaxed trial
enriched relaxed trial
```

此外还存在 PGD basis candidate 与 rollback。

如果同时维护：

```text
flat state
```

和：

```text
structured mutable state
```

则很容易出现：

- 一个 trial 更新 flat data 但没有更新 structured data；
- rollback 只恢复其中一套；
- fixed-basis / enriched trial 比较时读取不同 representation；
- Eq. (77) 与 local stage 使用不同 physical state；
- FOM-to-LATIN mapping 出现 silent indexing error。

因此，本阶段明确冻结：

> tower LATIN state 只能有一个 canonical mutable representation。

---

# 12. 与 tower equilibrium operator 的统一空间维数

上一阶段已经建立 tower reference equilibrium operator 的结构：

$$ \mathcal E_{\rm tower}=H(H^TMC_0H)^{-1}H^TMC_0. $$

采用 $q$ 作为统一 spatial coordinate 后：

material-point source strain：

$$ p\in\mathbb R^{N_q}. $$

compatible material-point strain：

$$ \varepsilon_{\rm comp}\in\mathbb R^{N_q}. $$

equilibrated material-point stress：

$$ \sigma_{\rm eq}\in\mathbb R^{N_q}. $$

因此可以统一得到：

```text
LatinState spatial dimension
        =
PGD spatial-mode dimension
        =
search-direction spatial dimension
        =
M metric dimension
        =
equilibrium-operator material-point dimension
        =
Nq
```

这是 `(N_t,N_q)` canonical storage 的核心理论收益。

---

# 13. 与 material-point integration metric 的关系

对 material point：

$$ q=(e,g,f), $$

其空间积分权重概念上来自：

$$ v_q=J_e\,w_g\,A_f. $$

其中：

- $J_e$：element / Gauss mapping 的 Jacobian contribution；
- $w_g$：Gauss quadrature weight；
- $A_f$：section fiber area。

因此 material-point metric 可写成：

$$ M=\operatorname{diag}(v_1,\ldots,v_{N_q}). $$

随后 weighted Gram-Schmidt 可以统一写成：

$$ \langle p_i,p_j\rangle_M=p_i^TMp_j. $$

Eq. (77) 的 spatial integration 也应基于同一 material-point metric philosophy。

但本阶段尚未进一步冻结：

- $v_q$ 的最终单位；
- stress-to-force conversion；
- $M$ 是否显式构造成 dense / diagonal array；
- equilibrium operator 与 metric object 的具体代码 API。

这些属于后续 operator/interface contract，不在本阶段扩展。

---

# 14. Frozen decision 3：MaterialPointLayout 不属于 LatinState

`LatinStateTower` 描述：

> 当前 material-point LATIN solution state。

`MaterialPointLayout` 描述：

> 当前 tower discretisation 的 material-point topology。

因此 layout 不应作为每个 state copy 中重复保存、重复复制的 mutable data。

以下信息不应重复塞入每个 `LatinStateTower`：

```text
element_index
gauss_index
fiber_index
fiber_area
fiber_y
gauss_weight
element_jacobian
section geometry
```

这些属于 discretisation / layout / operator data。

四个主要 LATIN state layer：

$$ s_i $$

$$ \hat s_{i+1/2} $$

$$ \breve s_{i+1} $$

$$ s_{i+1} $$

应共享同一个 immutable material-point layout。

---

# 15. Frozen decision 4：PGD basis 不属于 LatinState

第 $j$ 个 PGD pair 包括：

$$ \bar{\varepsilon}^{p}_j\in\mathbb R^{N_q} $$

以及：

$$ \lambda_j\in\mathbb R^{N_t}. $$

这些属于 reduced representation：

```text
PGDBasis
```

而不是 physical LATIN state：

```text
LatinStateTower
```

原因是 same current LATIN iteration 中可能同时存在：

```text
fixed-basis Trial A
```

和：

```text
enriched Trial B
```

两者可对应不同 basis snapshot，但都基于相同：

$$ s_i $$

和：

$$ \hat s_{i+1/2}. $$

因此必须严格区分：

```text
physical LATIN state
```

与：

```text
reduced representation used to build the state
```

这也是 future rollback semantics 能否清晰实现的前提。

---

# 16. 进入 LatinState 字段分类前的理论边界

完成 canonical storage 后，下一问题是：

> 14 个 candidate entries 是否都属于同一种 state field？

答案是否定的。

本阶段进一步冻结：

- `time` 是 grid metadata；
- 9 个字段是 formal primary LATIN fields；
- 4 个字段是 integrated support histories；
- derived structural / PGD / diagnostic data 不属于这 14 个条目。

---

# 17. 原论文 formal LATIN state

原论文 formal state 可写为：

$$ s=\{\dot\varepsilon^p,\varepsilon^e,\dot X,\dot D,\sigma,Z,Y\}. $$

对于当前 mixed-hardening material：

$$ \dot X=\{\dot\alpha,\dot{\bar r}\} $$

以及：

$$ Z=\{\beta,\bar R\}. $$

展开后得到 9 个 formal primary fields：

$$
\dot\varepsilon^p,\;
\varepsilon^e,\;
\dot\alpha,\;
\dot{\bar r},\;
\dot D,\;
\sigma,\;
\beta,\;
\bar R,\;
Y.
$$

对应 future field names：

```text
plastic_strain_rate
elastic_strain
alpha_rate
r_bar_rate
damage_rate
stress
beta
R_bar
energy_release_rate
```

这些字段是 LATIN formal coordinates。

---

# 18. Frozen decision 5：Formal primary LATIN fields

future tower `LatinState` 中以下 9 个 material-point histories 被正式分类为：

> formal primary LATIN fields

```text
plastic_strain_rate
elastic_strain
alpha_rate
r_bar_rate
damage_rate
stress
beta
R_bar
energy_release_rate
```

每一个 field 的 canonical shape 为：

$$ (N_t,N_q). $$

注意：

> primary 的判断标准不是“这个变量能不能从别的变量算出来”，而是“它是否属于论文 formal LATIN state，并且是否作为 LATIN alternating/search-direction/convergence structure 的坐标直接参与算法”。

---

# 19. Integrated support histories

代码还必须保存以下 integrated variables：

$$ \varepsilon^p $$

$$ \alpha $$

$$ \bar r $$

$$ D $$

对应：

```text
plastic_strain
alpha
r_bar
damage
```

它们不属于原论文 formal state 的显式坐标，但对数值积分和下一 local stage 不可缺少。

因此本阶段冻结其语义为：

> integrated support histories

而不是：

> additional primary LATIN fields

---

# 20. Frozen decision 6：Integrated support histories

以下 4 个 material-point histories 归类为：

```text
plastic_strain
alpha
r_bar
damage
```

统一 shape：

$$ (N_t,N_q). $$

其主要职责包括：

- 为下一 local stage 提供 accumulated internal state；
- 保持 rate 与 history 的数值一致性；
- 支持 constitutive evaluation；
- 支持 trial / relaxation / rollback；
- 支持 diagnostics 与 consistency tests。

---

# 21. 为什么 beta 仍是 primary field

当前 state law 包含：

$$ \beta=C\alpha. $$

因此从纯代数角度：

$$ \alpha\rightarrow\beta $$

是可计算关系。

但不能因此删除 `beta`。

原因是：

- 原论文 formal state 使用 thermodynamic force set $Z$；
- 当前模型中 $Z$ 包含 $\beta$；
- search-direction relation直接作用于 $\beta$；
- Eq. (73) hardening finishing直接涉及 $\beta$；
- Eq. (77) mechanical norm直接读取 $\beta$。

因此：

```text
alpha
    = integrated support history

beta
    = formal primary LATIN field
```

即使在某些 admissible state 中满足：

$$ \beta=C\alpha, $$

它们仍然属于不同的 solver semantics。

---

# 22. 为什么 R_bar 仍是 primary field

同理：

$$ \bar R=R_\infty\bar r. $$

虽然 $\bar R$ 可由 $\bar r$ 得到，但：

- $\bar R$ 属于 formal thermodynamic force coordinate；
- search-direction relation直接使用 $H_{\bar R}$ 与 $\bar R$；
- Eq. (74) hardening finishing直接处理该 force；
- Eq. (77) 直接使用 $\bar R$。

因此：

```text
r_bar
    = integrated support history

R_bar
    = formal primary LATIN field
```

不能因为存在代数 relation 就合并两者。

---

# 23. 为什么 Y 不能降级为普通 derived property

current constitutive relation 中：

$$ Y=Y(\sigma,D). $$

例如一维材料点在拉压条件下分别有不同的 energy-release relation。

因此在 unrelaxed complete global candidate 上，应重新计算：

$$ \breve Y=Y(\breve\sigma,\breve D). $$

但这并不意味着 `energy_release_rate` 应从 state 中删除。

原因有两层。

第一，原论文 formal LATIN state 明确包含：

$$ Y. $$

第二，relaxation 后 tower v1 已经冻结：

$$ Y_{i+1}=(1-\mu)Y_i+\mu\breve Y_{i+1}. $$

而不执行额外的 nonlinear reprojection：

$$ Y_{i+1}\leftarrow Y(\sigma_{i+1},D_{i+1}). $$

因此 relaxed LATIN iterate 一般允许：

$$ Y_{i+1}\neq Y(\sigma_{i+1},D_{i+1}). $$

这个 mismatch 是 LATIN alternating state machine 的允许状态，而不是 immediate numerical error。

所以：

```text
energy_release_rate
```

必须是一个真正保存的 formal primary field，而不能仅是动态计算 property。

---

# 24. damage_rate 与 damage 的严格区分

原论文 Eq. (75) 直接规定的是 damage rate：

$$ \dot D. $$

因此：

```text
damage_rate
```

属于 formal primary LATIN field。

而：

```text
damage
```

是 integrated support history。

tower v1 当前已经进一步选择：

$$ \breve{\dot D}=\hat{\dot D} $$

以及：

$$ \breve D=\hat D. $$

但后一条是 tower v1 对 integrated history 的 discrete treatment，不改变 formal state classification。

所以必须始终保持：

```text
damage_rate
    = formal primary field

damage
    = integrated support history
```

---

# 25. time 的特殊地位

`time` 虽然计入 future state object 的 candidate entries，但它不是 material state variable。

它的准确语义应为：

> immutable common time-grid metadata

shape：

$$ (N_t,). $$

其他 13 个 material-point histories：

$$ (N_t,N_q). $$

因此 future state geometry 应理解为：

```text
time
    (Nt,)

formal primary fields
    9 x (Nt, Nq)

integrated support histories
    4 x (Nt, Nq)
```

而不是“14 个字段具有相同 shape”。

---

# 26. Frozen decision 7：14 个 candidate entries 的第一层分类

本阶段正式冻结如下。

## 26.1 Grid metadata

```text
time
```

数量：

```text
1
```

shape：

$$ (N_t,). $$

## 26.2 Formal primary LATIN fields

```text
plastic_strain_rate
elastic_strain
alpha_rate
r_bar_rate
damage_rate
stress
beta
R_bar
energy_release_rate
```

数量：

```text
9
```

shape：

$$ (N_t,N_q). $$

## 26.3 Integrated support histories

```text
plastic_strain
alpha
r_bar
damage
```

数量：

```text
4
```

shape：

$$ (N_t,N_q). $$

总计：

$$ 1+9+4=14. $$

---

# 27. 为什么本阶段不在 14 个 entries 内再划一组 derived fields

最初 interface-level问题中包含：

> 哪些字段属于 derived fields？

本阶段进一步明确：

> 真正的 derived structural / PGD / diagnostic quantities 不应该成为 `LatinStateTower` 的第四类 internal storage。

原因是 `LatinStateTower` 应保持语义单一：

> material-point LATIN solution state。

如果将所有方便访问的 derived data 都塞入 state，会使 state object 变成：

- material state；
- structural response；
- PGD basis；
- operator cache；
- residual container；
- convergence history；

的混合体。

这会破坏 trial / copy / rollback 的语义边界。

---

# 28. 不属于 LatinStateTower 的主要 derived / external quantities

以下量不属于本阶段冻结的 14 个 entries。

## 28.1 Structural kinematic / response quantities

```text
nodal displacement correction
section axial strain
section curvature
section axial force
section bending moment
compatible total material-point strain
```

## 28.2 Global-stage correction quantities

```text
damage residual strain
plastic forcing
plastic strain correction
plastic strain-rate correction
plastic stress correction
damage stress correction
compatible strain correction
```

## 28.3 PGD quantities

```text
PGD spatial modes
PGD temporal functions
PGD temporal rates
basis snapshots
candidate mode
mode significance
orthogonalisation diagnostics
```

## 28.4 Search-direction / metric quantities

```text
H_sigma
H_beta
H_R_bar
material-point integration metric M
reference equilibrium operator data
```

## 28.5 Residual / outer-control quantities

```text
full mechanical residual
reduced residual
LATIN indicator xi
saturation zeta
stagnation diagnostics
acceptance diagnostics
```

这些量分别应由独立 module / result / diagnostics object 管理。

---

# 29. LatinStateTower、MaterialPointLayout、PGDBasis 的语义边界

本阶段形成以下三层明确分工。

## 29.1 LatinStateTower

负责：

> 当前 material-point LATIN physical state。

包含：

- time grid reference / metadata；
- 9 formal primary fields；
- 4 integrated support histories。

## 29.2 MaterialPointLayout

负责：

> tower discretisation 到 flat material-point coordinate $q$ 的 immutable mapping。

包含：

- $q\leftrightarrow(e,g,f)$；
- 与 material-point topology 有关的固定信息；
- 后续可能关联 integration metadata。

## 29.3 PGDBasis

负责：

> reduced representation。

包含：

- accepted spatial modes；
- temporal functions；
- temporal rates；
- mode metadata；
- basis snapshots / candidate basis semantics。

三者不能互相替代。

---

# 30. 与 tower FOM trial / commit / revert 的关系

当前 tower FOM 已经具有成熟的 material-state transaction：

```text
committed material state
    ↓
Newton trial
    ↓
converged
    ↓
commit
```

或者：

```text
Newton trial fails
    ↓
revert to last commit
```

future LATIN state transaction 与此不同。

LATIN 层面至少存在：

```text
accepted relaxed state s_i
    ↓
local state hat{s}
    ↓
fixed-basis trial
    ↓
optional enriched trial
    ↓
iteration commit
```

因此：

> FOM material trial/commit semantics 与 LATIN outer trial/commit semantics 是两个不同层级。

本阶段的 `(N_t,N_q)` canonical state 与 field classification，就是为了下一阶段能够定义这两个 transaction layer 之间的明确接口，而不是把它们混成同一套 state flags。

---

# 31. 本阶段对 same-iteration trial semantics 的直接影响

此前已经冻结：

- current accepted $s_i$ 在整个 LATIN iteration 内保持 immutable；
- fixed-basis trial 与 enriched trial共享同一个 $s_i$ baseline；
- local state 与 search directions在 one enrichment event 内冻结；
- mode acceptance失败时 rollback basis；
- trial state不是 accepted state，直到 current iteration commit。

本阶段新增的 canonical storage 使这些规则更清晰：

```text
accepted state:
    one canonical (Nt, Nq) state

fixed-basis trial:
    another canonical (Nt, Nq) state

enriched trial:
    another canonical (Nt, Nq) state

all states:
    share same immutable MaterialPointLayout
```

因此 future implementation 不需要依赖 in-place mutation 来节省表示层级。

---

# 32. 与 Eq. (77) 的直接关系

Eq. (77) mechanical norm 使用的 formal fields包括：

$$
\sigma,\;
\beta,\;
\bar R,\;
\dot\varepsilon^p,\;
\varepsilon^e,\;
\dot\alpha,\;
\dot{\bar r}.
$$

因此这些 field 必须能够直接从 `LatinStateTower` 中读取。

而：

$$ D,\;\dot D,\;Y $$

不直接进入原论文 Eq. (77) mechanical norm。

但这不意味着它们不属于 state：

- $\dot D$ 和 $Y$ 仍属于 formal primary state；
- $D$ 属于 integrated support history；
- 它们通过 constitutive coupling 影响 $\sigma$、$\varepsilon^e$ 和 search directions。

所以后续不得通过“是否进入 Eq. (77)”来决定一个 field 是否应被 state 保存。

---

# 33. 当前仍未冻结的内容

本阶段完成的是 storage 和第一层 field classification，而不是完整 implementation contract。

以下内容尚未最终冻结：

- 4 个 support histories 的具体 initial-condition handling；
- 各 field 在 $s_i$、$\hat s$、$\breve s$、$s_{i+1}$ 中分别必须满足哪些 exact relation；
- local stage 对每个 field 的 read / write ownership；
- global finishing 对每个 field 的 read / write ownership；
- relaxation 后哪些 relation必须保持 machine precision；
- 哪些 manifold relations允许 temporary mismatch；
- `LatinStateTower` 的具体 class API；
- copy / clone / readonly view 的 Python 语义；
- memory order；
- dtype policy；
- structured accessor API；
- `MaterialPointLayout` 的实际 data members；
- $M$ 的最终离散单位与实现；
- tower FOM 与 LATIN state 的 conversion interface。

这些必须按下一阶段逐字段继续闭合。

---

# 34. 本阶段冻结结论总览

本阶段可压缩为以下核心结论。

第一：

> tower LATIN-PGD canonical spatial coordinate 采用 flat material-point index $q$。

第二：

> 每个 $q$ 唯一对应 `(element, Gauss point, section fiber)`。

第三：

> 所有 material-point state histories 统一采用 `(N_t,N_q)` canonical storage。

第四：

> `(e,g,f)` hierarchy 通过 immutable `MaterialPointLayout` 保存，不作为第二套 mutable state。

第五：

> `LatinStateTower` 与 `MaterialPointLayout` 分离。

第六：

> `PGDBasis` 与 `LatinStateTower` 分离。

第七：

> `time` 是 grid metadata，shape 为 `(N_t,)`。

第八：

> 9 个 formal primary LATIN fields 统一 shape 为 `(N_t,N_q)`。

第九：

> 4 个 integrated support histories 统一 shape 为 `(N_t,N_q)`。

第十：

> `beta` 与 `R_bar` 即使可由 support histories代数计算，仍属于 formal primary fields。

第十一：

> `energy_release_rate` 即使在 unrelaxed state 可由 $\sigma,D$ 计算，仍属于 formal primary LATIN field，并且 relaxation 后允许与 nonlinear manifold 暂时不一致。

第十二：

> `damage_rate` 是 formal primary field，`damage` 是 integrated support history。

第十三：

> structural corrections、PGD basis、search directions、residuals、Eq. (77) diagnostics 和 outer histories 不属于 `LatinStateTower`。

第十四：

> current stage 仍然不开始 tower LATIN-PGD solver 代码实现。

---

# 35. 本阶段 frozen state geometry

future tower state 的第一层 geometry 现在可以写成：

```text
MaterialPointLayout
    q <-> (element, gauss, fiber)
    immutable

LatinStateTower
    time
        (Nt,)

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

        each:
            (Nt, Nq)

    integrated support histories
        plastic_strain
        alpha
        r_bar
        damage

        each:
            (Nt, Nq)

PGDBasis
    separate from LatinStateTower

structural / equilibrium / diagnostics
    separate from LatinStateTower
```

---

# 36. 下一阶段准确任务

下一阶段不再讨论 storage 方案。

下一阶段应开始：

> field-by-field contract。

第一组建议集中处理：

$$
\varepsilon^p,\;
\dot{\varepsilon}^p,\;
\varepsilon^e,\;
\sigma.
$$

这是因为这四个字段直接连接：

```text
PGD plastic correction
    ↓
tower equilibrium projection
    ↓
mechanical assembly
    ↓
complete unrelaxed candidate
    ↓
relaxation
    ↓
Eq. (77)
```

下一阶段对这四个字段逐项冻结：

- physical meaning；
- canonical shape；
- initial value；
- local-stage read / write ownership；
- global-stage read / write ownership；
- fixed-basis trial behavior；
- enriched-trial behavior；
- relaxation behavior；
- exact relations；
- allowed manifold mismatch；
- Eq. (77) participation；
- rollback requirements。

在这一组闭合后，再进入 hardening fields。

---

# 37. 后续 Markdown 编写规范

后续所有阶段总结继续继承 `2026-08-17-markdown-pycharm-math-preview-debugging-lessons.md` 的规范。

特别注意：

1. 数学环境中禁止直接书写原始小于号；
2. 小于关系使用 `\lt`；
3. 大于关系优先使用 `\gt` 或 `\ge`；
4. 独立公式使用 `$$ ... $$`；
5. 行内公式使用 `$...$`；
6. 不用复杂 `array` 数学环境包裹中文说明；
7. 长流程优先使用 Markdown code block；
8. 数学公式只承载数学结构，解释写在公式外；
9. 若 PyCharm Preview 发生级联失效，优先寻找第一条异常公式；
10. 排错采用最小复现与 binary split；
11. 一次只改变一个变量。

本文档本身按上述规则编写。

---

# 38. 本阶段最终研究停点

截至本阶段结束，future tower LATIN state 已经完成：

```text
material-point coordinate choice
    ↓
canonical storage choice
    ↓
immutable layout separation
    ↓
PGD basis separation
    ↓
14-entry first-level classification
```

当前准确停点为：

> 开始逐字段冻结 `plastic_strain`、`plastic_strain_rate`、`elastic_strain` 和 `stress` 的 state contract。

在该 field-by-field contract 完成之前：

> 不开始 tower LATIN-PGD solver 代码修改。
