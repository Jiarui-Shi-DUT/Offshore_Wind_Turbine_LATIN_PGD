# 海上风机塔筒 LATIN-PGD：Eq. (61)–(64) 新模态富集起点、剩余缺陷与速率形式桥接阶段总结

**日期：2026-08-16**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**当前分支：`feature/offshore-wind-turbine-tower-fatigue`**  
**前置总结：`docs/2026-08-16-tower-latin-pgd-eq60-saturation-and-latin-norm-stage-summary.md`**  
**本阶段范围：原论文 Eq. (61)–(64)**  
**阶段目标：在不改变原论文 `x-t` LATIN-PGD 基本数学结构的前提下，厘清 fixed-basis temporal update 之后 new PGD mode enrichment 的起点、剩余 search-direction defect 的重新定义、rank-one space-time pair 的含义，以及 Eq. (64) 如何把材料点 search-direction equation 桥接到后续 Galerkin 空间问题。**  
**下一阶段：原论文 Eq. (65) kinematic admissibility in rate form，并进一步进入 Eq. (66)–(71) spatial Galerkin problem。**

---

# 1. 本阶段的定位

上一阶段已经完成：

1. 原论文 Eq. (58)–(59) fixed spatial basis temporal update；
2. Eq. (60) saturation criterion；
3. Eq. (76)–(77) LATIN relative indicator 与 mechanical norm；
4. 一维三材料杆中成熟 adaptive PGD / LATIN convergence control 的回查；
5. 确认塔筒阶段应继承一维杆已经调试成熟的两层控制逻辑，而不重新发明 convergence control。

因此，当前算法已经推进到原论文 Fig. 2 中：

```text
existing PGD basis
        ↓
Eq. (58)–(59)
fixed-basis temporal update
        ↓
Eq. (60)
basis adequacy / saturation decision
        ↓
existing basis insufficient
        ↓
PGD: Add a pair
        ↓
Eq. (61)–(72)
new space-time pair enrichment
```

本阶段正式进入：

$$ \boxed{\text{PGD: Add a pair}} $$
但本阶段只处理其起始部分：

$$ \boxed{\text{Eq. (61)–(64)}} $$
暂不进入 Eq. (65)–(72) 的完整空间 Galerkin 与 fixed-point 求解。

---

# 2. 本阶段坚持的总研究路线

当前研究路线继续保持：

$$ \boxed{ \text{original }x-t\text{ LATIN-PGD} \longrightarrow \text{fiber beam-column offshore wind turbine tower} } $$
现阶段不引入：

$$ n-\tau-x $$
三变量 PGD，也不立即引入：

- cycle jump；
- temporal homogenisation；
- multi-time-scale LATIN；
- 所有内部变量 simultaneous PGD；
- 基于 FOM snapshot 的 online SVD surrogate；
- 真实随机风浪荷载下的直接降阶；
- 超长寿命 $10^6\sim10^8$ 周期的跳跃策略。

当前第一目标仍然是：

> **尽可能忠实地把 Bhattacharyya et al. (2018) 的原始 LATIN-PGD，从 bar/continuum 空间离散迁移到 fiber beam-column tower 空间离散。**

因此，本阶段 Eq. (61)–(64) 的推导首先完全站在原论文框架内进行，只在最后讨论其对 tower 离散的直接含义。

---

# 3. 原论文 enrichment 阶段为什么从 Eq. (61) 开始

原论文 global-stage plastic branch 的基本 descent search-direction relation 已经在 Eq. (41) 写成：

$$ \boxed{ \Delta\dot{\varepsilon}^{p}_{i+1} - H_\sigma\Delta\sigma'_{i+1} + \bar{\Delta}_{i+1} =0 } $$
其中：

- $\Delta\dot{\varepsilon}^{p}_{i+1}$：plastic-strain-rate correction；
- $\Delta\sigma'_{i+1}$：plastic-dependent stress correction；
- $H_\sigma$：descent search-direction operator 的 plastic part；
- $\bar{\Delta}_{i+1}$：在给定 baseline state 上评估该 descent relation 后留下的 defect。

在 fixed spatial basis temporal update 之前，baseline 是上一个 accepted LATIN global state：

$$ s_i. $$
但 Eq. (58)–(59) 已经利用旧 PGD basis 得到了一个新的 update-stage candidate：

$$ \boxed{s_{i+1}^{\mathrm{up}}} $$
因此，若 Eq. (60) 判定：

$$ \boxed{\text{existing basis insufficient}} $$
则 enrichment 阶段不应该重新从 $s_i$ 开始逼近完整 local/global gap，而应该以：

$$ \boxed{s_{i+1}^{\mathrm{up}}} $$
作为新的 baseline，只补偿 update 后剩余的缺陷。

这就是 Eq. (61)–(62) 的真正算法背景。

---

# 4. Eq. (61)：new-mode enrichment 所求的不是“完整 correction”

原论文 Eq. (61)：

$$ \boxed{ \Delta\dot{\varepsilon}^{p}_{i+1} - H_\sigma\Delta\sigma'_{i+1} + \bar{\Delta}_{i+1} =0 } \tag{61} $$
形式上与前面的 Eq. (41) 几乎相同。

但在 enrichment subsection 中，$\Delta$ 的具体算法含义已经发生变化。

## 4.1 Eq. (41) 阶段

可理解为：

$$ s_i \longrightarrow s_{i+1} $$
所需的 global correction。

即：

$$ \Delta s = s_{i+1}-s_i. $$
---

## 4.2 Eq. (61) 阶段

fixed-basis temporal update 已经给出：

$$ s_{i+1}^{\mathrm{up}}. $$
因此 new-mode enrichment 实际上求的是：

$$ \boxed{ \delta s_{m+1} = s_{i+1}^{\mathrm{new}} - s_{i+1}^{\mathrm{up}} } $$
而不是：

$$ s_{i+1}^{\mathrm{new}}-s_i. $$
为了避免符号混淆，本阶段讨论中可临时写成：

$$ \delta\dot{\varepsilon}^{p}, \qquad \delta\sigma', $$
而原论文仍继续使用：

$$ \Delta\dot{\varepsilon}^{p}_{i+1}, \qquad \Delta\sigma'_{i+1}. $$
所以 Eq. (61) 最准确的理解是：

$$ \boxed{ \text{在 }s_{i+1}^{\mathrm{up}}\text{ 基础上，求一个新的 PGD enrichment correction} } $$
---

# 5. Update 与 enrichment 的严格关系

原论文 hybrid PGD strategy 不是：

$$ \boxed{\text{update 或 enrichment 二选一}} $$
而是：

$$ \boxed{\text{先 update，必要时再 enrichment}} $$
完整逻辑：

$$ \text{existing }m\text{ spatial modes} $$
$$ \downarrow $$
$$ \text{Eq. (58)–(59): update old temporal functions} $$
$$ \downarrow $$
$$ s_{i+1}^{\mathrm{up}} $$
$$ \downarrow $$
$$ \text{Eq. (60): evaluate basis adequacy} $$
若已有 basis 仍能提供足够 LATIN indicator improvement：

$$ \boxed{\text{accept current global approximation}} $$
若改善不足：

$$ \boxed{\text{add a new space-time pair}} $$
因此：

> **Eq. (61) 不是重新求一次 Eq. (41)，而是在 Eq. (59) 已经利用旧 basis 之后，对“旧 basis 解释不了的剩余部分”进行 rank-one enrichment。**

---

# 6. “Residual enrichment”的核心认识

设 local-stage target 与原 accepted global state 之间的初始 gap 为：

$$ \hat{s}_{i+1/2}-s_i. $$
经过 fixed-basis temporal update 后：

$$ s_i \longrightarrow s_{i+1}^{\mathrm{up}}. $$
已有 PGD basis 已经解释了部分 gap。

因此 new mode 应面对的是：

$$ \boxed{ \hat{s}_{i+1/2}-s_{i+1}^{\mathrm{up}} } $$
而不是再次面对：

$$ \hat{s}_{i+1/2}-s_i. $$
所以可以把 PGD enrichment 概括为：

$$ \boxed{\text{residual enrichment}} $$
即：

> **每一个新 PGD mode 只负责补偿当前 reduced basis 尚不能有效表示的剩余 LATIN search-direction defect。**

这一认识对 tower 实现非常重要，因为如果 new mode 仍针对完整 gap 构造，就会与已有 modes 产生大量重复表示，破坏 greedy enrichment 的意义。

---

# 7. Eq. (62)：为什么必须重新定义 $\bar{\Delta}_{i+1}$

## 7.1 原始 descent search direction

从原论文 global stage 的 descent search direction，取 plastic branch：

$$ \boxed{ \dot{\varepsilon}^{p}_{i+1} - \hat{\dot{\varepsilon}}^{p}_{i+1/2} - H_\sigma \left( \sigma_{i+1}-\hat{\sigma}_{i+1/2} \right) =0 } \tag{A} $$
---

## 7.2 以 $s_i$ 为 baseline 时

写：

$$ \dot{\varepsilon}^{p}_{i+1} = \dot{\varepsilon}^{p}_{i} + \Delta\dot{\varepsilon}^{p}_{i+1} $$
$$ \sigma_{i+1} = \sigma_i + \Delta\sigma_{i+1}. $$
代入 Eq. (A)：

$$ \Delta\dot{\varepsilon}^{p}_{i+1} - H_\sigma\Delta\sigma_{i+1} + \left[ H_\sigma (\hat{\sigma}_{i+1/2}-\sigma_i) - ( \hat{\dot{\varepsilon}}^{p}_{i+1/2} - \dot{\varepsilon}^{p}_{i} ) \right] =0. $$
于是得到原论文 Eq. (40)：

$$ \boxed{ \bar{\Delta}_{i+1}^{(40)} = H_\sigma (\hat{\sigma}_{i+1/2}-\sigma_i) - \left( \hat{\dot{\varepsilon}}^{p}_{i+1/2} - \dot{\varepsilon}^{p}_{i} \right) } $$
---

## 7.3 enrichment 阶段 baseline 已变成 $s_{i+1}^{\mathrm{up}}$

此时写：

$$ \boxed{ \dot{\varepsilon}^{p}_{i+1} = \dot{\varepsilon}^{p,\mathrm{up}}_{i+1} + \delta\dot{\varepsilon}^{p} } $$
$$ \boxed{ \sigma_{i+1} = \sigma^{\mathrm{up}}_{i+1} + \delta\sigma' } $$
代回 Eq. (A)：

$$ \dot{\varepsilon}^{p,\mathrm{up}}_{i+1} + \delta\dot{\varepsilon}^{p} - \hat{\dot{\varepsilon}}^{p}_{i+1/2} - H_\sigma \left[ \sigma^{\mathrm{up}}_{i+1} + \delta\sigma' - \hat{\sigma}_{i+1/2} \right] =0. $$
整理：

$$ \delta\dot{\varepsilon}^{p} - H_\sigma\delta\sigma' + \left[ H_\sigma \left( \hat{\sigma}_{i+1/2} - \sigma^{\mathrm{up}}_{i+1} \right) - \left( \hat{\dot{\varepsilon}}^{p}_{i+1/2} - \dot{\varepsilon}^{p,\mathrm{up}}_{i+1} \right) \right] =0. $$
于是自然定义：

$$ \boxed{ \bar{\Delta}_{i+1} = H_\sigma \left( \hat{\sigma}_{i+1/2} - \sigma^{\mathrm{up}}_{i+1} \right) - \left( \hat{\dot{\varepsilon}}^{p}_{i+1/2} - \dot{\varepsilon}^{p,\mathrm{up}}_{i+1} \right) } \tag{62} $$
---

# 8. Eq. (40) 与 Eq. (62) 的统一理解

两者实际上具有完全相同的通式：

$$ \boxed{ \bar{\Delta} = H_\sigma (\hat{\sigma}-\sigma_{\mathrm{baseline}}) - ( \hat{\dot{\varepsilon}}^p - \dot{\varepsilon}^p_{\mathrm{baseline}} ) } $$
区别仅在：

## Eq. (40)

$$ \boxed{ s_{\mathrm{baseline}}=s_i } $$
## Eq. (62)

$$ \boxed{ s_{\mathrm{baseline}}=s_{i+1}^{\mathrm{up}} } $$
因此：

$$ \boxed{ \text{Eq. (62) = Eq. (40) with shifted baseline} } $$
即：

> **Eq. (62) 并没有改变 LATIN descent search direction，本质上只是把 search-direction defect 的评估点从旧 accepted state 移到了 temporal update 后的 candidate state。**

---

# 9. $\bar{\Delta}$ 的准确物理/数值含义

$\bar{\Delta}$ 可以称为 residual，但必须避免与其他 residual 混淆。

更准确的名称是：

$$ \boxed{ \text{descent-search-direction defect} } $$
或者：

$$ \boxed{ \text{search-direction residual evaluated at the current baseline state} } $$
它不是以下三个量。

## 9.1 不是结构平衡残差

不是：

$$ F_{\mathrm{ext}}-F_{\mathrm{int}}. $$
因此它不应被解释成 Newton equilibrium residual。

---

## 9.2 不是 Eq. (59) reduced least-squares residual

Eq. (59) 中的 residual 是：

$$ r = \Delta\dot{\varepsilon}^p - H_\sigma\Delta\sigma' + \bar{\Delta} $$
在已有 reduced basis 上进行机械范数最小化的 residual。

$\bar{\Delta}$ 是该 residual 中的 known forcing / defect term，而不是 Eq. (59) residual 本身。

---

## 9.3 不是 LATIN indicator $\xi$

$$ \xi $$
衡量：

$$ \boxed{ \text{local-stage state 与 global-stage state 的整体 mechanical distance} } $$
而：

$$ \bar{\Delta}(x,t) $$
是定义在 material-point space-time domain 上的 search-direction defect field。

所以：

$$ \boxed{ \bar{\Delta} \neq r_{\mathrm{Eq.59}} \neq \xi } $$
---

# 10. Eq. (61) 的另一种有用写法

Eq. (61)：

$$ \delta\dot{\varepsilon}^{p} - H_\sigma\delta\sigma' + \bar{\Delta} =0 $$
可写成：

$$ \boxed{ \delta\dot{\varepsilon}^{p} - H_\sigma\delta\sigma' = -\bar{\Delta} } $$
因此 new mode 的任务不是“直接拟合 $\bar{\Delta}$”这个 scalar/field，而是：

> **生成一个 coupled plastic-rate / stress correction，使其在 LATIN search-direction operator 下产生 $-\bar{\Delta}$，从而抵消当前 baseline state 的 defect。**

这一点说明：

$$ \boxed{ \text{new PGD mode ≠ ordinary SVD mode of }\bar{\Delta} } $$
因为它必须同时满足：

- constitutive/search-direction relation；
- structural admissibility；
- stress equilibrium；
- PGD separated representation。

---

# 11. 一个重要极限检查

如果 temporal update 已经使：

$$ \sigma^{\mathrm{up}}_{i+1} = \hat{\sigma}_{i+1/2} $$
且：

$$ \dot{\varepsilon}^{p,\mathrm{up}}_{i+1} = \hat{\dot{\varepsilon}}^{p}_{i+1/2}, $$
则：

$$ \boxed{ \bar{\Delta}_{i+1}=0 } $$
此时 Eq. (61) 为：

$$ \delta\dot{\varepsilon}^{p} - H_\sigma\delta\sigma' =0. $$
其中：

$$ \delta\dot{\varepsilon}^{p}=0, \qquad \delta\sigma'=0 $$
自然满足。

所以：

$$ \boxed{ \bar{\Delta}=0 \Rightarrow \text{当前 update state 在该 descent relation 上无需 enrichment} } $$
这进一步验证了 Eq. (62) 的含义。

---

# 12. Eq. (63)：new rank-one space-time pair

在 Eq. (61)–(62) 得到 remaining defect 后，原论文规定新增 correction 采用一个 separable pair：

$$ \boxed{ \Delta\dot{\varepsilon}^{p}_{i+1}(x,t) = \Delta\dot{\lambda}_{m+1}(t) \bar{\varepsilon}^{p}_{m+1}(x) } \tag{63a} $$
以及：

$$ \boxed{ \Delta\sigma'_{i+1}(x,t) = \Delta\lambda_{m+1}(t) \bar{\mathbb C} \bar{\varepsilon}^{p}_{m+1}(x) } \tag{63b} $$
其中：

$$ \boxed{ \bar{\mathbb C} = \mathbb C(\mathcal E-I) } $$
由 Eq. (52) 已经定义。

---

# 13. 为什么叫 “Add a pair”

新增的不是一个单独 spatial vector，而是：

$$ \boxed{ \left\{ \bar{\varepsilon}^{p}_{m+1}(x), \; \Delta\lambda_{m+1}(t) \right\} } $$
即：

$$ \boxed{ \text{space function} \times \text{time function} } $$
组成一个新的 rank-one correction。

因此：

$$ \boxed{ \text{PGD: Add a pair} } $$
就是：

$$ \boxed{ \text{add one new rank-one space-time contribution} } $$
---

# 14. Eq. (63) 中真正独立的 spatial unknown 只有一个

虽然 Eq. (63) 同时出现 plastic-strain spatial mode 与 stress spatial mode，但它们不是两个独立未知量。

主空间未知量：

$$ \boxed{ \bar{\varepsilon}^{p}_{m+1}(x) } $$
而 stress spatial mode：

$$ \boxed{ \bar{\sigma}_{m+1}(x) = \bar{\mathbb C} \bar{\varepsilon}^{p}_{m+1}(x) } $$
由 equilibrium operator 派生。

逻辑仍然是：

$$ \boxed{ \bar{\varepsilon}^p \longrightarrow \mathcal E\bar{\varepsilon}^p \longrightarrow \mathbb C(\mathcal E-I)\bar{\varepsilon}^p } $$
因此：

- plastic-strain spatial mode：primary spatial unknown；
- compatible total-strain mode：derived quantity；
- stress spatial mode：derived quantity。

---

# 15. 为什么 plastic strain rate 用 $\dot{\lambda}$，stress 用 $\lambda$

这是 Eq. (63) 最容易混淆的问题之一。

原论文首先对 plastic strain correction 写：

$$ \boxed{ \Delta\varepsilon^p(x,t) = \lambda(t)\bar{\varepsilon}^p(x) } $$
对时间求导：

$$ \boxed{ \Delta\dot{\varepsilon}^p(x,t) = \dot{\lambda}(t)\bar{\varepsilon}^p(x) } $$
所以 plastic strain rate 前自然出现：

$$ \boxed{\dot{\lambda}} $$
---

另一方面，stress correction 来自：

$$ \Delta\sigma' = \mathbb C ( \Delta\varepsilon' - \Delta\varepsilon^p ). $$
而 compatible total-strain correction 具有同一个 time function：

$$ \Delta\varepsilon' = \lambda(t)\mathcal E\bar{\varepsilon}^p. $$
因此：

$$ \Delta\sigma' = \lambda(t) \mathbb C (\mathcal E-I) \bar{\varepsilon}^p. $$
即：

$$ \boxed{ \Delta\sigma' = \lambda(t)\bar{\sigma}(x) } $$
所以：

$$ \boxed{ \text{plastic strain rate}\rightarrow\dot{\lambda} } $$
而：

$$ \boxed{ \text{stress correction}\rightarrow\lambda } $$
二者来自同一个 underlying plastic-strain mode，而不是两个独立 PGD ansatz。

---

# 16. Eq. (63) 的统一关系链

$$ \boxed{ \Delta\varepsilon^p = \lambda(t)\bar{\varepsilon}^p(x) } $$
一方面：

$$ \frac{\partial}{\partial t} $$
得到：

$$ \boxed{ \Delta\dot{\varepsilon}^p = \dot{\lambda}(t)\bar{\varepsilon}^p(x) } $$
另一方面：

$$ \text{equilibrium + elasticity} $$
得到：

$$ \boxed{ \Delta\sigma' = \lambda(t)\bar{\mathbb C}\bar{\varepsilon}^p(x) } $$
所以 Eq. (63) 中的两个 correction 具有同一个空间母模态与相关联的时间函数结构。

---

# 17. 为什么 Eq. (63) 后必须采用 fixed-point algorithm

将 Eq. (63) 代入 Eq. (61)：

$$ \dot{\lambda}(t)\bar{\varepsilon}^p(x) - H_\sigma(x,t) \lambda(t) \bar{\mathbb C}\bar{\varepsilon}^p(x) + \bar{\Delta}(x,t) =0. $$
其中：

$$ \lambda(t) $$
未知，

同时：

$$ \bar{\varepsilon}^p(x) $$
也未知。

因此整个问题对这两个未知量联合来看具有 bilinear / nonlinear coupling。

不能像 Eq. (59) 那样在 spatial basis 已知的情况下直接只求 temporal coefficients。

所以原论文采用 alternating fixed-point strategy：

$$ \boxed{ \lambda^{(k)} \rightarrow \bar{\varepsilon}^{p,(k+1)} \rightarrow \lambda^{(k+1)} \rightarrow \bar{\varepsilon}^{p,(k+2)} \rightarrow\cdots } $$
其中：

- spatial function：通过 Galerkin problem 求；
- temporal function：通过 mechanical residual minimisation 求。

这就是 Eq. (64)–(72) 后续推导的总体结构。

---

# 18. Eq. (58) 与 Eq. (63) 的本质区别

## Eq. (58)：reuse existing basis

已知：

$$ \bar{\varepsilon}_1^p,\ldots,\bar{\varepsilon}_m^p. $$
未知仅为：

$$ \Delta\lambda_1(t),\ldots,\Delta\lambda_m(t). $$
所以：

$$ \boxed{ \text{spatial basis known} + \text{temporal coefficients unknown} } $$
---

## Eq. (63)：enrichment

新增：

$$ \bar{\varepsilon}_{m+1}^p(x) $$
与：

$$ \Delta\lambda_{m+1}(t) $$
均未知。

因此：

$$ \boxed{ \text{new spatial mode unknown} + \text{new temporal function unknown} } $$
这正是为什么：

- Eq. (59) 是 reduced temporal minimisation；
- Eq. (63)–(72) 必须采用 alternating fixed-point solve。

---

# 19. 一个新 mode 一般不能完全消除 remaining defect

对于复杂的：

$$ \bar{\Delta}(x,t), $$
一个 rank-one pair：

$$ \lambda_{m+1}(t)\bar{\varepsilon}_{m+1}^p(x) $$
一般无法一次完全解释。

因此：

$$ m \rightarrow m+1 \rightarrow m+2 \rightarrow\cdots $$
通过 greedy enrichment 逐步增加：

$$ \sum_{j=1}^{m} \lambda_j(t)\bar{\varepsilon}_j^p(x). $$
所以 PGD enrichment 的思想不是寻找一个“完美空间模态”，而是：

$$ \boxed{ \text{successive rank-one correction construction} } $$
即逐步构造能够表达当前 nonlinear global correction 的 reduced space-time representation。

---

# 20. Eq. (64)：从 strain partition 到 rate form

原论文 Eq. (64)：

$$ \boxed{ \Delta\dot{\varepsilon}' = \Delta\dot{\varepsilon}^{p} + \mathbb C^{-1}\Delta\dot{\sigma}' } \tag{64} $$
这一式本身没有引入新的物理假设。

---

## 20.1 从应变分解开始

对于 plastic-related branch：

$$ \boxed{ \Delta\varepsilon' = \Delta\varepsilon^e + \Delta\varepsilon^p } $$
而 linear elastic relation：

$$ \Delta\sigma' = \mathbb C \Delta\varepsilon^e. $$
因此：

$$ \Delta\varepsilon^e = \mathbb C^{-1}\Delta\sigma'. $$
代入：

$$ \boxed{ \Delta\varepsilon' = \Delta\varepsilon^p + \mathbb C^{-1}\Delta\sigma' } $$
---

## 20.2 对时间求导

在 global-stage reference linear problem 中：

$$ \mathbb C $$
为固定 Hooke operator，因此：

$$ \frac{d}{dt} \left( \mathbb C^{-1}\Delta\sigma' \right) = \mathbb C^{-1}\Delta\dot{\sigma}'. $$
于是：

$$ \boxed{ \Delta\dot{\varepsilon}' = \Delta\dot{\varepsilon}^{p} + \mathbb C^{-1}\Delta\dot{\sigma}' } $$
得到 Eq. (64)。

所以：

$$ \boxed{ \text{Eq. (64)} = \text{strain partition} + \text{linear elastic relation} + \text{time derivative} } $$
---

# 21. 为什么 Eq. (64) 必须写成 rate form

Eq. (61) 中直接出现：

$$ \Delta\dot{\varepsilon}^p. $$
如果继续保留：

$$ \Delta\varepsilon' = \Delta\varepsilon^p + \mathbb C^{-1}\Delta\sigma' $$
的非 rate 形式，就无法直接把 Eq. (61) 的 search-direction relation 插入后续 kinematic admissibility。

因此 Eq. (64) 的核心作用是：

$$ \boxed{ \text{把 total-strain relation 也转成与 Eq. (61) 相同的 rate level} } $$
使 local search-direction equation 与 global kinematic admissibility 可以直接连接。

---

# 22. Eq. (61) + Eq. (64) 的关键组合

由 Eq. (61)：

$$ \Delta\dot{\varepsilon}^p - H_\sigma\Delta\sigma' + \bar{\Delta} =0 $$
得到：

$$ \boxed{ \Delta\dot{\varepsilon}^p = H_\sigma\Delta\sigma' - \bar{\Delta} } $$
代入 Eq. (64)：

$$ \Delta\dot{\varepsilon}' = H_\sigma\Delta\sigma' - \bar{\Delta} + \mathbb C^{-1}\Delta\dot{\sigma}'. $$
整理：

$$ \boxed{ \Delta\dot{\varepsilon}' = H_\sigma\Delta\sigma' + \mathbb C^{-1}\Delta\dot{\sigma}' - \bar{\Delta} } \tag{B} $$
这是 Eq. (64) 最重要的结果。

左侧：

$$ \Delta\dot{\varepsilon}' $$
已经是结构 compatible total-strain-rate correction。

右侧：

- $H_\sigma\Delta\sigma'$：search-direction stress contribution；
- $\mathbb C^{-1}\Delta\dot{\sigma}'$：elastic stress-rate contribution；
- $-\bar{\Delta}$：known remaining search-direction defect。

因此 Eq. (64) 实际完成：

$$ \boxed{ \text{local constitutive/search-direction information} \longrightarrow \text{global kinematic-admissibility-ready relation} } $$
---

# 23. 将 Eq. (63) 代入 Eq. (64)

定义纯 spatial stress mode：

$$ \boxed{ \bar{\sigma}(x) = \bar{\mathbb C}\bar{\varepsilon}^p(x) } $$
并保持：

$$ \boxed{ \Delta\sigma'(x,t) = \lambda(t)\bar{\sigma}(x) } $$
则：

$$ \boxed{ \Delta\dot{\sigma}'(x,t) = \dot{\lambda}(t)\bar{\sigma}(x) } $$
代入上式：

$$ \boxed{ \Delta\dot{\varepsilon}' = H_\sigma \lambda \bar{\sigma} + \dot{\lambda} \mathbb C^{-1}\bar{\sigma} - \bar{\Delta} } \tag{C} $$
这就是后续 Eq. (65)–(66) 进入 Galerkin 积分前真正使用的关系。

---

# 24. 为什么 $\bar{\sigma}$ 必须是纯 spatial function

为了与 Eq. (52)、Eq. (53)、Eq. (63) 保持一致，本阶段统一采用：

$$ \boxed{ \bar{\sigma}(x) = \bar{\mathbb C} \bar{\varepsilon}^p(x) } $$
而：

$$ \boxed{ \Delta\sigma'(x,t) = \lambda(t)\bar{\sigma}(x) } $$
因此：

$$ \bar{\sigma} $$
本身不能再包含 $\lambda(t)$。

否则：

$$ \Delta\sigma' $$
会错误变成：

$$ \lambda^2 $$
形式，并与原论文 separated representation 冲突。

这也是后续 Eq. (66) 中出现：

$$ \langle H_\sigma\lambda^2\rangle $$
而不是更高次 $\lambda$ 的必要前提。

---

# 25. Eq. (64) 为什么是 Eq. (61) 到 Eq. (65) 的桥梁

到 Eq. (61) 为止，问题仍然主要写在：

$$ \boxed{ \text{material-point constitutive/search-direction space} } $$
Eq. (64) 之后：

$$ \Delta\dot{\varepsilon}' $$
被写成：

$$ H_\sigma\Delta\sigma' + \mathbb C^{-1}\Delta\dot{\sigma}' - \bar{\Delta}. $$
而：

$$ \Delta\dot{\varepsilon}' $$
是 compatible total-strain-rate correction。

因此它可以直接进入：

$$ \boxed{ \text{kinematic admissibility weak form} } $$
这就是下一步 Eq. (65)。

因此：

$$ \boxed{ \text{Eq. (64) 是 local search-direction relation 与 global Galerkin spatial problem 之间的桥梁} } $$
---

# 26. Eq. (61)–(64) 的完整逻辑链

本阶段可以压缩为：

$$ \boxed{ \text{Eq. (58)–(59): old basis temporal update} } $$
$$ \downarrow $$
$$ s_{i+1}^{\mathrm{up}} $$
$$ \downarrow $$
$$ \boxed{ \text{Eq. (60): existing basis insufficient} } $$
$$ \downarrow $$
重新计算 remaining descent-search-direction defect：

$$ \boxed{ \bar{\Delta}_{i+1} = H_\sigma ( \hat{\sigma}_{i+1/2} - \sigma_{i+1}^{\mathrm{up}} ) - ( \hat{\dot{\varepsilon}}^p_{i+1/2} - \dot{\varepsilon}^{p,\mathrm{up}}_{i+1} ) } $$
$$ \downarrow $$
Eq. (61)：

$$ \boxed{ \delta\dot{\varepsilon}^p - H_\sigma\delta\sigma' = -\bar{\Delta} } $$
$$ \downarrow $$
Eq. (63) rank-one ansatz：

$$ \boxed{ \delta\dot{\varepsilon}^p = \dot{\lambda}\bar{\varepsilon}^p } $$
$$ \boxed{ \delta\sigma' = \lambda\bar{\sigma} } $$
$$ \downarrow $$
Eq. (64)：

$$ \boxed{ \delta\dot{\varepsilon}' = \delta\dot{\varepsilon}^p + \mathbb C^{-1}\delta\dot{\sigma}' } $$
$$ \downarrow $$
组合后：

$$ \boxed{ \delta\dot{\varepsilon}' = H_\sigma\lambda\bar{\sigma} + \dot{\lambda}\mathbb C^{-1}\bar{\sigma} - \bar{\Delta} } $$
$$ \downarrow $$
下一阶段：

$$ \boxed{\text{Eq. (65): kinematic admissibility in rate form}} $$
---

> **PyCharm 公式兼容说明：** 当前预览器不支持 `\boldsymbol`。以下 tower material-point 离散向量统一用 `\vec{\cdot}` 表示；这只是显示记号调整，不改变前述数学定义与推导。

# 27. 本阶段对 tower migration 的直接含义

虽然本阶段尚未正式推导 tower Eq. (65)–(72)，但 Eq. (61)–(64) 已经可以基本原样迁移。

连续 spatial coordinate：

$$ x $$
在 tower 中离散为：

$$ (s,y)\rightarrow(e,g,f). $$
因此新的 plastic-strain spatial mode 可写为：

$$ \boxed{ \vec{\bar{\varepsilon}}^p_{m+1} = [ \bar{\varepsilon}^p_{egf} ] } $$
当前 coarse tower：

$$ N_e=10,\qquad N_g=2,\qquad N_f=16, $$
所以：

$$ \boxed{ N_q=N_eN_gN_f=320 } $$
一个 tower plastic spatial mode 为：

$$ \boxed{ \vec{\bar{\varepsilon}}^p_{m+1} \in\mathbb R^{320} } $$
---

# 28. Tower Eq. (63) 的直接对应

可直接写为：

$$ \boxed{ \delta\dot{\vec{\varepsilon}}^p(t) = \dot{\lambda}_{m+1}(t) \vec{\bar{\varepsilon}}^p_{m+1} } $$
对应 equilibrated stress mode：

$$ \boxed{ \vec{\bar{\sigma}}_{m+1} = C_0 ( \mathcal E_{\mathrm{tower}}-I ) \vec{\bar{\varepsilon}}^p_{m+1} } $$
因此：

$$ \boxed{ \delta\vec{\sigma}'(t) = \lambda_{m+1}(t) \vec{\bar{\sigma}}_{m+1} } $$
这里：

$$ \mathcal E_{\mathrm{tower}} = H ( H^TMC_0H )^{-1} H^TMC_0 $$
为前一阶段已经建立的 fixed reference tower equilibrium projection。

---

# 29. Tower Eq. (61)–(62) 应保持 material-point level

对于 flatten 后的 fiber material-point index：

$$ q=1,\ldots,N_q, $$
remaining defect 应定义在：

$$ \boxed{ (t,q) } $$
space-time grid 上：

$$ \boxed{ \vec{\bar{\Delta}}(t) = D_H(t) \left( \vec{\hat{\sigma}}(t) - \vec{\sigma}^{\mathrm{up}}(t) \right) - \left( \hat{\dot{\vec{\varepsilon}}}^p(t) - \dot{\vec{\varepsilon}}^{p,\mathrm{up}}(t) \right) } $$
其中：

$$ D_H(t) = \operatorname{diag} ( H_{\sigma,1}(t),\ldots,H_{\sigma,N_q}(t) ). $$
所以 tower enrichment 仍然遵循：

$$ \boxed{ \text{fiber material-point local information} \rightarrow \text{new PGD spatial mode} \rightarrow \text{tower-equilibrated stress mode} } $$
而不是直接从 nodal displacement residual 或 tower base moment residual 构造新 mode。

---

# 30. Tower Eq. (64) 的直接对应

tower material-point vector form：

$$ \boxed{ \delta\dot{\vec{\varepsilon}}' = \delta\dot{\vec{\varepsilon}}^p + C_0^{-1} \delta\dot{\vec{\sigma}}' } $$
代入 search direction：

$$ \boxed{ \delta\dot{\vec{\varepsilon}}' = D_H(t) \delta\vec{\sigma}' + C_0^{-1} \delta\dot{\vec{\sigma}}' - \vec{\bar{\Delta}} } $$
再代入：

$$ \delta\vec{\sigma}' = \lambda(t)\vec{\bar{\sigma}}, $$
得到：

$$ \boxed{ \delta\dot{\vec{\varepsilon}}' = \lambda(t)D_H(t)\vec{\bar{\sigma}} + \dot{\lambda}(t)C_0^{-1}\vec{\bar{\sigma}} - \vec{\bar{\Delta}}(t) } $$
这一形式就是下一阶段 tower Eq. (65)–(66) spatial weak form 的直接起点。

---

# 31. 当前阶段不应做的事情

为了保持与原论文一致，本阶段不应提前做以下改动：

1. 不把 new mode 改成 nodal displacement mode；
2. 不把 new mode 改成 curvature mode；
3. 不把 new mode 改成 section moment mode；
4. 不直接对 damage 做独立 enrichment；
5. 不用 FOM snapshot SVD mode 代替 LATIN-PGD enrichment mode；
6. 不跳过 Eq. (59) update 直接每次新增 mode；
7. 不使用完整 local/global gap 代替 Eq. (62) remaining defect；
8. 不把 $\bar{\Delta}$ 当作 Newton equilibrium residual；
9. 不把 Eq. (59) least-squares residual 与 Eq. (62) defect 混为同一个量；
10. 不把 $\xi$ 或 $\zeta$ 直接作为 material-point forcing；
11. 不在 $\bar{\sigma}$ 中再次乘入 time function；
12. 不在此阶段引入 $n-\tau-x$ PGD。

---

# 32. 本阶段已经解决的问题

经过 Eq. (61)–(64) 的逐式推导，当前已经明确：

1. Eq. (61) 为什么在 Eq. (60) 后重新出现；
2. Eq. (61) 中的 correction 为什么是 enrichment correction，而不是整个 LATIN correction；
3. fixed-basis temporal update 与 new-mode enrichment 的先后关系；
4. 为什么 enrichment 必须以 $s_{i+1}^{\mathrm{up}}$ 为新的 baseline；
5. Eq. (62) 如何严格从 descent search direction 推导；
6. Eq. (40) 与 Eq. (62) 的统一形式；
7. $\bar{\Delta}$ 的准确含义；
8. $\bar{\Delta}$ 与 Newton residual、Eq. (59) residual、LATIN indicator $\xi$ 的区别；
9. Eq. (63) 为什么代表一个 new rank-one space-time pair；
10. 为什么真正独立的 spatial unknown 仍只有 $\bar{\varepsilon}^p$；
11. 为什么 stress spatial mode 由 equilibrium operator 派生；
12. 为什么 plastic-strain-rate term 使用 $\dot{\lambda}$；
13. 为什么 stress correction 使用 $\lambda$；
14. 为什么 Eq. (63) 后需要 fixed-point alternating algorithm；
15. Eq. (58) reuse 与 Eq. (63) enrichment 的本质区别；
16. Eq. (64) 如何由 strain partition 与 linear elasticity 推出；
17. 为什么 Eq. (64) 必须写成 rate form；
18. Eq. (64) 如何把 local search-direction relation 桥接到 global kinematic admissibility；
19. Eq. (61)–(64) 为什么仍然可以基本原样迁移到 fiber beam-column tower；
20. tower enrichment 为什么仍应在 fiber material-point field 上构造，而不是在 nodal/section resultant space 上直接构造。

---

# 33. 本阶段尚未解决的问题

以下问题留到下一阶段：

1. 原论文 Eq. (65) 为什么用 statically admissible stress field $\sigma^*\in S_0$ 检验 kinematic admissibility；
2. 为什么 correction problem 在 Eq. (65) 中右端为 0；
3. Eq. (65) 如何代入：

$$ \Delta\dot{\varepsilon}' = H_\sigma\lambda\bar{\sigma} + \dot{\lambda}\mathbb C^{-1}\bar{\sigma} - \bar{\Delta}; $$
4. Eq. (66) 中：

$$ \langle H_\sigma\lambda^2\rangle, \qquad \langle\lambda\dot{\lambda}\rangle, \qquad \langle\bar{\Delta}\lambda\rangle $$
分别如何出现；
5. Eq. (67) 中 $W^{-1}$ 与 $\bar{\delta}$ 的物理和数值意义；
6. Eq. (68)–(70) 如何把 stress-admissibility problem 转成 classical displacement-like spatial solve；
7. Eq. (71) 如何恢复新的 plastic-strain spatial mode；
8. Eq. (72) 如何更新对应 temporal function；
9. Eq. (65)–(71) 在 fiber beam-column tower 上的完整离散形式；
10. 是否可直接复用现有 `equilibrium_operator.py` / tower reference stiffness 作为 Eq. (70) 空间求解骨架；
11. 新 tower mode 的 Gram–Schmidt orthonormalisation 与 acceptance / rejection 如何接入现有 PGD basis 管理。

---

# 34. 当前与代码实现的关系

本阶段仍属于：

$$ \boxed{\text{theory derivation / algorithm closure stage}} $$
当前不应立即修改 tower solver 主流程。

现有一维代码中的：

```text
latin/pgd_enrichment.py
latin/pgd_time_update.py
latin/pgd_basis.py
latin/pgd_global_stage.py
latin/pgd_solver.py
```

可作为后续 tower implementation 的成熟参考，但下一步仍应先完成原论文 Eq. (65)–(72) 的理论闭合。

尤其应避免：

> 在尚未完整推导 Eq. (65)–(71) 的情况下，直接把一维 `pgd_enrichment.py` 的矩阵形式机械复制到 tower。

原因是 tower 的真正变化集中在：

$$ \boxed{ \text{spatial admissibility / equilibrium operator} } $$
而这一部分正是 Eq. (65)–(71) 要解决的内容。

---

# 35. 本阶段最终结论

本阶段最核心的结论可以压缩为以下五条。

## 结论 1：new-mode enrichment 是 residual enrichment

$$ \boxed{ \text{new mode does not reconstruct the full LATIN correction} } $$
而是：

$$ \boxed{ \text{correct the remaining defect after fixed-basis temporal update} } $$
---

## 结论 2：Eq. (62) 只是把 Eq. (40) 的 baseline 从 $s_i$ 移到 $s_{i+1}^{up}$

$$ \boxed{ \bar{\Delta} = H_\sigma(\hat{\sigma}-\sigma_{\mathrm{baseline}}) - ( \hat{\dot{\varepsilon}}^p - \dot{\varepsilon}^p_{\mathrm{baseline}} ) } $$
---

## 结论 3：Eq. (63) 只增加一个 rank-one pair

$$ \boxed{ \delta\varepsilon^p(x,t) = \lambda(t)\bar{\varepsilon}^p(x) } $$
主 spatial unknown 是：

$$ \boxed{\bar{\varepsilon}^p(x)} $$
而 stress mode 由 equilibrium operator 派生。

---

## 结论 4：Eq. (64) 是 local → global 的桥梁

$$ \boxed{ \delta\dot{\varepsilon}' = \delta\dot{\varepsilon}^p + \mathbb C^{-1}\delta\dot{\sigma}' } $$
使 Eq. (61) 的 local search-direction relation 可以进入后续 global kinematic admissibility weak form。

---

## 结论 5：到 Eq. (64) 为止，tower 不需要改变原论文 `x-t` PGD 数学结构

只需将：

$$ x $$
的空间离散从 continuum/bar material points 换成：

$$ (e,g,f) $$
fiber material points。

因此仍可保持：

$$ \boxed{ \delta\dot{\vec{\varepsilon}}^p(t) = \dot{\lambda}(t) \vec{\bar{\varepsilon}}^p } $$
$$ \boxed{ \delta\vec{\sigma}'(t) = \lambda(t) C_0 ( \mathcal E_{\mathrm{tower}}-I ) \vec{\bar{\varepsilon}}^p } $$
---

# 36. 下一阶段唯一任务

下一阶段从原论文 Eq. (65) 开始：

$$ \boxed{ \int_{[0,T]\times\Omega} \Delta\dot{\varepsilon}' : \sigma^* \,d\Omega\,dt = 0, \qquad \forall \sigma^*\in S_0 } \tag{65} $$
首先只解决：

1. 为什么这是 kinematic admissibility；
2. 为什么 test field 是 statically admissible stress；
3. 为什么 correction problem 的右端为 0；
4. 为什么采用：

$$ \sigma^*(x,t) = \lambda(t)\bar{\sigma}^*(x) $$
后可以进入 Eq. (66)；
5. 这一 weak form 在 fiber beam-column tower 中应如何解释，但暂不立即写代码。

确认 Eq. (65) 后，再进入 Eq. (66)。

---

# 参考论文

Bhattacharyya, M., Fau, A., Nackenhorst, U., Néron, D., & Ladevèze, P. (2018).  
*A LATIN-based model reduction approach for the simulation of cycling damage*.  
**Computational Mechanics, 62**, 725–743.  
DOI: 10.1007/s00466-017-1523-z

本阶段重点对应原论文：

- Eq. (41)；
- Eq. (52)–(53)；
- Eq. (58)–(60)；
- Eq. (61)–(64)；
- Fig. 2 hybrid LATIN-PGD algorithm。