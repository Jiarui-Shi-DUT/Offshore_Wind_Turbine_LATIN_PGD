# 海上风机塔筒 LATIN-PGD：100 周期 FOM、低秩诊断与 PGD 增广判据阶段总结

**日期：2026-08-12**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**当前分支：`feature/offshore-wind-turbine-tower-fatigue`**

---

## 1. 本阶段工作的总体目标

本阶段是在已经成功复现一维三材料杆 LATIN-PGD 算法的基础上，将研究推进到**海上风机钢塔筒高周循环疲劳问题**。当前重点不是立即构造最终塔筒 LATIN-PGD 求解器，而是先完成以下基础工作：

1. 建立可稳定运行的 100 周期塔筒全阶模型（FOM）；
2. 将昂贵的 100 周期 FOM 结果冻结为可重复使用的 snapshot 数据；
3. 判断塔筒循环响应是否具有适合 PGD 的低秩可分离结构；
4. 分析低秩结构是否随疲劳演化发生变化；
5. 明确原论文已有的 PGD enrichment criterion 是否可以直接继承；
6. 明确后续真正需要创新和推导的部分。

目前已经形成较清晰的结论：

> **塔筒 100 周期 FOM 具有明显低秩结构；不可逆变量的低秩复杂度会随疲劳过程发生阶段性变化；但是在线 PGD enrichment criterion 不需要重新发明，应优先继承原论文基于 LATIN indicator 饱和度的 $\zeta$ 判据。真正需要研究的核心，是如何把原论文的 $x-t$ 两变量 PGD 扩展为适合大量循环问题的 $n-\tau-x$ 三变量分离表示。**

---

# 2. 已有研究基础

## 2.1 一维杆 LATIN-PGD 复现

前期已经完成原论文一维三材料杆算例复现。主要特点包括：

- 90 个一维杆单元；
- 20 个循环；
- 三种材料通过不同屈服应力区分；
- 粘塑性、混合硬化和损伤；
- LATIN local/global stage 交替；
- PGD 在 global stage 中用于塑性相关修正量；
- 已建立 search directions、内部变量更新、收敛判据和 PGD enrichment 控制逻辑；
- 已使用 Git/GitHub 进行版本管理。

最终损伤约为：

$$D_1\approx0.224,\qquad D_2\approx0.184,\qquad D_3\approx0.152.$$

因此，后续塔筒工作不是从零开始，而是在已有 LATIN-PGD 主框架上进行结构层级与时间尺度扩展。

---

# 3. 塔筒结构模型

## 3.1 几何

当前塔筒基于 NREL 5 MW 名义塔筒简化：

$$H=87.6\ \mathrm{m},$$

$$D_b=6.00\ \mathrm{m},\qquad D_t=3.87\ \mathrm{m},$$

$$t_b=0.027\ \mathrm{m},\qquad t_t=0.019\ \mathrm{m}.$$

直径和壁厚沿高度线性变化。

## 3.2 梁理论

当前采用二维 Euler-Bernoulli 梁。现阶段重点是纵向弯曲疲劳以及 LATIN-PGD 方法扩展，因此优先采用较清晰、计算代价较低的梁模型。后续如有需要，可再比较 Timoshenko 梁剪切效应。

## 3.3 纤维截面

正式较高分辨率候选方案：

$$40\ \text{beam elements} \times4\ \text{Gauss/element} \times64\ \text{fibers}.$$

当前 100 周期 FOM 为控制计算量采用：

$$10\times2\times16\times1=320$$

个纤维积分空间点，即：

- 10 beam elements；
- 2 Gauss points / element；
- 16 circumferential fibers；
- 1 radial layer。

---

# 4. 材料模型

每个纤维积分点采用一维粘塑性-损伤模型，状态变量为：

$$\mathbf{s}_{\rm mat} = [\varepsilon_p,\alpha,\bar r,D].$$

其中：

- $\varepsilon_p$：塑性应变；
- $\alpha$：运动硬化内部变量；
- $\bar r$：各向同性硬化相关内部变量；
- $D$：损伤。

当前材料参数主要用于算法验证和机制研究，尚未完成针对真实海上风机塔筒钢材的实验标定。因此现阶段结果不能直接解释为真实塔筒寿命预测。

---

# 5. 非对称循环加载

正式 100 周期分析采用：

$$F_{\max}=+1.0\ \mathrm{MN},$$

$$R_F=-0.5,$$

因此：

$$F_{\min}=-0.5\ \mathrm{MN}.$$

平均值和幅值：

$$F_{\rm mean}=0.25\ \mathrm{MN},$$

$$F_{\rm amp}=0.75\ \mathrm{MN}.$$

周期：

$$T=10.$$

加载路径：

$$F_{\rm mean} \rightarrow F_{\max} \rightarrow F_{\rm mean} \rightarrow F_{\min} \rightarrow F_{\rm mean}.$$

每周期：

$$40\ \text{increments/cycle}.$$

---

# 6. 100 周期 FOM 冻结参考

## 6.1 冻结原因

正式 100 周期 FOM 一次求解耗时约：

$$3039\ \mathrm{s} \approx50.65\ \mathrm{min}.$$

为避免每次低秩分析都重新计算 FOM，本阶段建立 snapshot 持久化机制，并冻结为：

```text
outputs/tower_100cycle_fom_reference_v1.npz
```

文件大小约：

$$14.543\ \mathrm{MiB}.$$

## 6.2 数据契约

保存字段：

```text
snapshot_format_version
cycle_numbers
phase_times
phase_fractions
phase_forces
analysis_times
nodal_displacements
fiber_strains
fiber_stresses
fiber_states
```

不保存可重新计算的派生量，例如：

- `fiber_plastic_strains`
- `fiber_alphas`
- `fiber_r_bars`
- `fiber_damages`
- cycle increment
- SVD
- rank

这样 snapshot 文件保持为“原始全阶参考数据”。

---

# 7. Snapshot 持久化与全局周期切片

在：

```text
examples/nonlinear_tower_snapshot_tensor.py
```

中已实现：

```text
save_tower_cycle_phase_snapshots(...)
load_tower_cycle_phase_snapshots(...)
select_tower_cycle_range(...)
```

`select_tower_cycle_range(...)` 按**真实全局周期编号**切片。例如 Stage II：

```text
cycles 21-46
```

切片后仍保持：

```text
cycle_numbers = [21, ..., 46]
```

而不是重新编号为 `[1, ..., 26]`。

当前周期编号允许：

```text
[1, ..., 100]
[1, ..., 20]
[21, ..., 46]
[47, ..., 100]
```

但不允许非正、间断或逆序编号。

---

# 8. 当前 Git 检查点

## 8.1 Snapshot 持久化

```text
2934a65
feat: persist 100-cycle tower FOM snapshots
```

## 8.2 全局周期切片

```text
a481a1f
feat: add global cycle stage selection
```

已推送到：

```text
origin/feature/offshore-wind-turbine-tower-fatigue
```

---

# 9. 100 周期 FOM 正式结果

冻结数据主要尺寸：

```text
u shape     = (100, 41, 33)
sigma shape = (100, 41, 10, 2, 16)
state shape = (100, 41, 10, 2, 16, 4)
```

即：

$$100\ \text{cycles} \times41\ \text{phase points} \times320\ \text{fiber points}.$$

最终：

$$\max|\varepsilon_p| = 1.363490242\times10^{-3},$$

$$\max D = 4.409564790\times10^{-2}.$$

临界纤维：

```text
element = 1
Gauss   = 1
fiber   = 5
```

临界高度约：

$$z\approx1.8512\ \mathrm{m},$$

纤维坐标：

$$y\approx2.9451\ \mathrm{m}.$$

最大 Newton 迭代次数：

$$4$$

最大自由自由度残差：

$$8.761\times10^{-3}\ \mathrm{N}.$$

说明冻结数据可作为后续离线低秩分析参考。

---

# 10. Cycle-phase 张量化

为利用循环载荷的重复结构，将完整时间响应重构为：

$$q(n,\tau,x),$$

其中：

- $n$：slow cycle coordinate；
- $\tau$：fast phase coordinate；
- $x$：space coordinate。

这意味着后续不再只考虑原论文：

$$q(x,t),$$

而进一步研究：

$$q(n,\tau,x).$$

这也是未来 cycle-phase-space PGD 的基础。

---

# 11. Cycle increment

定义：

$$\Delta q(n,\tau,x) = q(n,\tau,x)-q(n,0,x).$$

周期终点不人为归零，因此：

$$\Delta q(n,\tau_{\rm end},x)$$

仍保留：

- ratcheting；
- 塑性累积；
- 损伤累积；
- 不可逆漂移。

因此 cycle increment 特别适合观察疲劳不可逆演化。

---

# 12. 当前低秩诊断的性质

当前采用的是 mode-wise SVD / HOSVD-style diagnostic。

对：

$$q(n,\tau,x)$$

分别构造：

### Cycle unfolding

$$n\times(\tau x)$$

### Phase unfolding

$$\tau\times(nx)$$

### Space unfolding

$$x\times(n\tau)$$

并得到：

$$(r_n,r_\tau,r_x).$$

必须强调：

> **这些是 HOSVD multilinear ranks，不是 CP rank，也不是 PGD separation rank。**

因此 `(2,1,2)` 不能直接解释为 `PGD rank = 2`，只能说明三个方向具有很强的低维可压缩性。

---

# 13. 全 100 周期低秩结果

99.99% squared-Frobenius energy 下：

## Raw fields

```text
u       : (2, 2, 2)
sigma   : (3, 2, 4)
eps_p   : (3, 1, 3)
D       : (3, 1, 3)
```

## Cycle-increment fields

```text
Delta u       : (2, 2, 1)
Delta sigma   : (1, 1, 2)
Delta eps_p   : (6, 3, 7)
Delta D       : (6, 2, 6)
```

主要认识：

1. 位移和应力高度低秩；
2. 累积塑性和累积损伤仍较低秩；
3. 不可逆增量 `Delta eps_p`、`Delta D` 的 cycle 和 space 方向更复杂；
4. phase 方向始终较低维；
5. 塔筒问题总体具备 PGD 可分离基础。

---

# 14. Stage I / II / III 分段

根据 100 周期损伤增量演化，暂时划分：

$$\text{Stage I}:1-20,$$

$$\text{Stage II}:21-46,$$

$$\text{Stage III}:47-100.$$

此前发现 cycle 46 附近损伤增量达到局部最低值，之后损伤增长速率重新提高。当前机制描述为：

> **damage-rate acceleration with continuing net-drift stabilization**

即净漂移仍趋稳，但损伤速率在后期重新加速。

---

# 15. Stage-wise 低秩诊断

脚本：

```text
examples/nonlinear_tower_stagewise_low_rank_probe.py
```

## Raw

```text
field      Stage I      Stage II     Stage III
u          (2,2,2)      (2,2,2)      (2,2,2)
sigma      (3,2,3)      (2,2,3)      (2,2,3)
eps_p      (3,2,2)      (2,1,2)      (2,1,2)
D          (3,2,2)      (2,1,2)      (2,1,2)
```

## Cycle increment

```text
field          Stage I      Stage II     Stage III
Delta u        (2,2,1)      (1,1,1)      (1,1,1)
Delta sigma    (2,2,2)      (1,1,1)      (1,1,1)
Delta eps_p    (4,2,4)      (3,1,3)      (3,2,4)
Delta D        (4,2,4)      (2,1,2)      (3,1,3)
```

初步表现为：

$$\boxed{ \text{Stage I 较复杂} \rightarrow \text{Stage II 最低秩} \rightarrow \text{Stage III 不可逆变量复杂度回升} }$$

但三个阶段长度分别为 20、26、54 个周期，因此还需要控制窗口长度。

---

# 16. 等长度窗口控制

采用四个 20-cycle 窗口：

$$W_1=1-20, \quad W_2=27-46, \quad W_3=47-66, \quad W_4=81-100.$$

脚本：

```text
examples/nonlinear_tower_equal_window_low_rank_probe.py
```

## 16.1 Delta eps_p

99.99% rank：

```text
W1: (4,2,4)
W2: (2,1,3)
W3: (2,2,3)
W4: (2,2,3)
```

W2 -> W3：

Cycle：

$$s_2/s_1: 2.556\times10^{-2} \rightarrow 5.108\times10^{-2},$$

约增加 1.998 倍。

Phase：

$$6.609\times10^{-3} \rightarrow 2.164\times10^{-2},$$

约增加 3.274 倍。

Space：

$$2.524\times10^{-2} \rightarrow 5.118\times10^{-2},$$

约增加 2.028 倍。

说明等长度条件下 cycle 46 前后塑性增量场仍发生明显结构变化。

## 16.2 Delta D

99.99% rank：

```text
W1: (4,2,4)
W2: (2,1,2)
W3: (2,1,2)
W4: (2,1,2)
```

虽然整数 rank 不变，但：

Cycle：

$$2.274\times10^{-2} \rightarrow 3.169\times10^{-2},$$

增加约 39.3%。

Space：

$$2.272\times10^{-2} \rightarrow 3.177\times10^{-2},$$

增加约 39.8%。

因此整数 rank 对渐进结构变化不够敏感。

---

# 17. 为什么不能只看整数 rank

99.99% rank 是经过阈值离散后的整数。

例如：

```text
(2,1,2) -> (2,1,2)
```

可能看起来没有变化，但第二奇异值相对第一奇异值已经增加数倍。

因此：

- rank 适合描述总体可压缩阶数；
- `s2/s1` 更适合描述次级模态重要程度的连续变化。

---

# 18. 20-cycle 滑动窗口分析

脚本：

```text
examples/nonlinear_tower_sliding_window_low_rank_probe.py
```

窗口长度 20 周期，步长 5 周期：

```text
1-20
6-25
11-30
...
81-100
```

结果写入：

```text
outputs/tower_100cycle_sliding_window_low_rank_v1.csv
```

---

# 19. Delta eps_p 的滑动窗口规律

## 19.1 Cycle direction

关键 `s2/s1_n`：

```text
1-20      2.357e-01
6-25      5.615e-02
11-30     2.875e-02
16-35     2.176e-02
21-40     1.861e-02
26-45     2.302e-02
31-50     3.816e-02
36-55     5.179e-02
41-60     5.612e-02
46-65     5.232e-02
51-70     4.571e-02
56-75     3.946e-02
61-80     3.429e-02
66-85     3.037e-02
71-90     2.757e-02
76-95     2.535e-02
81-100    2.336e-02
```

规律：

1. 早期复杂度很高；
2. 随后快速降低；
3. `21-40` 附近达到低值；
4. `26-45` 后重新增强；
5. `41-60` 附近达到局部高值；
6. 后期 cycle-to-cycle complexity 再次下降。

因此模态结构变化不是 cycle 47 的瞬时开关，更像一个约 cycle 30-50 的过渡带。

## 19.2 Phase direction

后半程持续增加，例如：

```text
21-40     4.817e-03
31-50     9.125e-03
41-60     1.676e-02
51-70     2.491e-02
61-80     3.328e-02
71-90     4.215e-02
81-100    5.159e-02
```

这说明：

> 后期即使 cycle-to-cycle 变化重新趋于规则，循环内部塑性路径仍在持续重构。

## 19.3 Space direction

后期仍保持明显次级空间模态：

```text
81-100:
s2/s1_x = 5.294e-02
```

因此后期塑性空间分布并没有回到固定单一空间模式。

---

# 20. Delta D 的滑动窗口规律

Cycle direction：

```text
21-40     1.193e-02
26-45     2.037e-02
31-50     3.172e-02
36-55     3.766e-02
41-60     3.669e-02
46-65     3.259e-02
...
81-100    1.482e-02
```

说明：

- Stage II 中后期损伤增量高度规则；
- 约 `26-45` 开始重新复杂；
- `36-55` 附近达到局部高值；
- 后期再逐渐降低。

Phase direction 始终约：

$$4\times10^{-3} \sim5\times10^{-3}.$$

说明：

> 损伤增量的 fast-phase 结构总体非常稳定。

Space direction 与 cycle direction 基本同步，因此损伤场主要表现为：

$$\boxed{ \text{cycle-space modal reorganisation} }$$

而不是 fast-phase waveform 的持续变化。

---

# 21. 对 cycle 46 的重新解释

现在应区分两个概念。

## 21.1 标量损伤速率转折

cycle 46 仍可以作为：

> 损伤增量 / 损伤速率局部最低附近的标量转折点。

## 21.2 全场模态结构转折

滑动窗口显示：

$$\boxed{ \text{约 cycle 30-50 是更合理的低秩结构重组过渡带} }$$

由于窗口长度为 20 周期，本身存在平滑效应，因此当前不应把模态结构转折硬性定位到某一个单独周期。

---

# 22. SVD 分析的正确定位

当前 SVD/HOSVD 工作应定义为：

$$\boxed{ \text{offline separability / low-rank feasibility study} }$$

其任务是回答：

1. 塔筒 FOM 是否低秩；
2. 哪些变量最易压缩；
3. cycle / phase / space 三个方向分别有多复杂；
4. 低秩结构是否随疲劳演化变化；
5. 是否有必要允许 basis 自适应更新。

当前答案是：

> **有明显低秩结构，并且不可逆变量的低秩结构会随疲劳过程变化，所以 adaptive basis 是合理且必要的。**

但 SVD 不应直接替代原论文的在线 enrichment criterion。

---

# 23. 重要路线修正：不重新设计 enrichment indicator

此前曾考虑基于 SVD residual energy 定义：

$$\eta^{(r)} = \frac{\|A-A^{(r)}\|_F}{\|A\|_F}.$$

该指标作为离线诊断是合理的，但重新对照原论文后确认：

> **原论文已经存在用于判断是否新增 PGD pair 的 saturation / enrichment indicator。**

因此当前路线调整为：

$$\boxed{ \text{不重新发明新的 online enrichment criterion} }$$

而是优先继承文献的：

$$\boxed{\zeta}$$

判据。

---

# 24. 原论文 Hybrid PGD strategy

原论文 global stage 中采用：

> **先 Update existing basis，再在必要时 Add a pair。**

假设已有 $m$ 个 PGD pairs：

$$\Delta\dot\varepsilon_{i+1}^{p} = \sum_{j=1}^{m} \Delta\dot\lambda_j(t) \bar\varepsilon_j^p(x).$$

首先固定已有空间基：

$$\bar\varepsilon_j^p(x),$$

只更新时间函数：

$$\Delta\lambda_j(t).$$

这是因为更新时间函数比生成新的空间函数便宜。

如果仅 update 就能继续有效降低 LATIN error，则无需 enrich。

只有当前 basis 出现饱和时，才增加：

$$\lambda_{m+1}(t) \bar\varepsilon_{m+1}^{p}(x).$$

---

# 25. 文献 LATIN convergence indicator xi

原论文 Eq. (76)：

$$\xi = \frac{ \left\| \hat{\mathbf{s}}^p_{i+1/2} - \mathbf{s}^p_{i+1} \right\| }{ \left\| \hat{\mathbf{s}}^p_{i+1/2} \right\| + \left\| \mathbf{s}^p_{i+1} \right\| }.$$

其中：

- $\hat{\mathbf{s}}^p_{i+1/2}$：local stage；
- $\mathbf{s}^p_{i+1}$：global stage。

因此：

$$\boxed{ \xi = \text{local-stage solution 与 global-stage solution 的相对距离} }$$

当：

$$\xi\rightarrow0$$

说明 LATIN 两个 manifold 的距离减小，算法趋于收敛。

---

# 26. 文献 saturation / enrichment indicator zeta

原论文 Eq. (60)：

$$\boxed{ \zeta = \frac{ \xi_i-\xi_{i+1} }{ \xi_i+\xi_{i+1} } }$$

它衡量的不是“当前误差有多大”，而是：

> **使用当前 PGD basis 后，LATIN indicator 还能下降多少。**

## 当 zeta 较大

若：

$$\zeta>\zeta^{\rm tol},$$

则：

$$\xi_i\rightarrow\xi_{i+1}$$

下降明显，说明已有 basis 仍有效，因此：

$$\boxed{\text{不新增 PGD pair}}$$

## 当 zeta 较小

若：

$$\zeta\lt\zeta^{\rm tol},$$

则：

$$\xi_i\approx\xi_{i+1}.$$

表示只更新时间函数已经无法显著改善解，当前 basis 出现饱和，因此：

$$\boxed{ \text{Add a new PGD pair} }$$

即增加：

$$\lambda_{m+1}(t) \bar\varepsilon_{m+1}^{p}(x).$$

---

# 27. xi 与 zeta 的区别

## xi

回答：

> 当前 LATIN local/global solution 相差多少？

属于：

$$\boxed{ \text{LATIN convergence indicator} }$$

## zeta

回答：

> 当前已有 PGD basis 是否还具有继续降低 xi 的能力？

属于：

$$\boxed{ \text{PGD basis saturation / enrichment indicator} }$$

逻辑为：

$$\xi \rightarrow \zeta \rightarrow \text{是否新增 PGD pair}.$$

---

# 28. 新模态 rejection

原论文新增一个 PGD pair 后并不是无条件保留。

新生成：

$$\lambda_{m+1}(t) \bar\varepsilon_{m+1}^{p}(x)$$

之后：

1. 使用 Gram-Schmidt 与已有空间基正交化；
2. 修正所有时间函数；
3. 检查新时间函数的范数；
4. 若新 pair 贡献过小，则 rejection。

因此原论文已有完整的自适应 basis 管理链：

```text
Update existing basis
        ↓
zeta saturation check
        ↓
Add a new pair if needed
        ↓
orthonormalisation
        ↓
small-mode rejection
```

---

# 29. zeta_tol 不能机械照抄

原论文不同算例采用不同阈值。

一维杆：

$$\zeta^{\rm tol}=0.1.$$

二维 L-shaped structure：

$$\zeta^{\rm tol}=10^{-2}.$$

因此：

$$\boxed{ \zeta^{\rm tol}\text{ 不是普适常数} }$$

它可能受以下因素影响：

- 结构复杂度；
- 时间离散；
- 空间离散；
- 非线性程度；
- search directions；
- PGD basis；
- 收敛要求。

因此塔筒中可以优先继承 Eq. (60) 的形式，但 `zeta_tol` 需要重新标定和敏感性分析。

---

# 30. 原论文 PGD 形式

原论文采用 space-time separation：

$$q(x,t) \approx \sum_{j=1}^{m} \lambda_j(t)X_j(x).$$

对于塑性应变率：

$$\Delta\dot\varepsilon^p = \sum_{j=1}^{m} \dot\lambda_j(t) \bar\varepsilon_j^p(x).$$

因此原论文显式分离：

$$\boxed{x-t}$$

两个维度。

---

# 31. 高周疲劳塔筒真正需要的扩展

大量循环问题若仍把整个时间轴完全视为单一变量 $t$，不能充分利用周期重复结构。

当前塔筒已经建立：

$$q(n,\tau,x).$$

因此目标表示为：

$$\boxed{ q(n,\tau,x) \approx \sum_{k=1}^{r} N_k(n) T_k(\tau) X_k(x) }$$

其中：

- $N_k(n)$：slow-cycle function；
- $T_k(\tau)$：fast-phase function；
- $X_k(x)$：space function。

这才是后续真正的方法扩展。

---

# 32. 真正的研究创新点

现阶段不应把创新重点放在：

$$\boxed{ \text{重新设计 enrichment criterion} }$$

因为原论文已经有：

$$\zeta = \frac{\xi_i-\xi_{i+1}} {\xi_i+\xi_{i+1}}.$$

真正需要研究的是：

$$\boxed{ \lambda(t)X(x) \longrightarrow N(n)T(\tau)X(x) }$$

即：

> 将原论文完整时间轴上的 space-time PGD，扩展为适合大量循环的 slow-cycle / fast-phase / space 三变量 PGD。

---

# 33. zeta 在三变量 PGD 中的预期继承

后续总体思想暂定为：

## Step 1：Update existing basis

假设已有：

$$q_r = \sum_{k=1}^{r} N_k(n)T_k(\tau)X_k(x).$$

优先复用已有 basis，更新相应函数。

## Step 2：计算 LATIN indicator

构造与原论文一致物理意义的：

$$\xi.$$

## Step 3：计算 saturation indicator

$$\zeta = \frac{ \xi_i-\xi_{i+1} }{ \xi_i+\xi_{i+1} }.$$

## Step 4：判断

若：

$$\zeta>\zeta^{\rm tol},$$

已有 basis 仍有效。

若：

$$\zeta\lt\zeta^{\rm tol},$$

则新增：

$$N_{r+1}(n) T_{r+1}(\tau) X_{r+1}(x).$$

因此：

> **zeta 的思想可以保留，真正需要扩展的是“Add a pair”内部的数学求解。**

---

# 34. Damage 在原论文中并非直接 PGD 分解对象

这一点必须特别注意。

原论文第 4.3 节中，hardening 和 damage variables 在 global stage 中并不是直接采用 PGD approximation，而是在各 Gauss point 上局部更新。

因此当前塔筒对：

$$D(n,\tau,x)$$

和：

$$\Delta D(n,\tau,x)$$

进行 SVD 的意义主要是：

- 观察损伤全场结构；
- 识别疲劳状态转变；
- 研究内部变量复杂度；
- 判断 PGD basis 是否可能受到内部状态演化影响。

这并不意味着：

$$D$$

必须直接成为 PGD separated variable。

后续推导时必须保持这一点。

---

# 35. 当前新增分析脚本

本阶段已新增并运行：

```text
examples/nonlinear_tower_stagewise_low_rank_probe.py
examples/nonlinear_tower_equal_window_low_rank_probe.py
examples/nonlinear_tower_sliding_window_low_rank_probe.py
```

并生成：

```text
outputs/tower_100cycle_sliding_window_low_rank_v1.csv
```

当前尚未再次执行 `git status` 并形成新的 Git checkpoint，因此后续应先检查工作区，再提交这些分析脚本和阶段总结。

---

# 36. 下一阶段暂不应直接写代码

下一阶段首先应该完成：

$$\boxed{ \text{原论文 }x-t\text{ PGD} \rightarrow \text{塔筒 }n-\tau-x\text{ PGD} }$$

的逐式数学推导。

在三变量数学结构明确以前，不宜直接进入代码实现。

---

# 37. 下一阶段需要解决的理论问题

## 37.1 原始 PGD 如何重写

从：

$$\Delta\dot\varepsilon^p(x,t) = \sum_{j=1}^{m} \dot\lambda_j(t) \bar\varepsilon_j^p(x)$$

扩展到类似：

$$\Delta\dot\varepsilon^p(n,\tau,x) = \sum_{j=1}^{m} N_j(n) \dot T_j(\tau) \bar\varepsilon_j^p(x)$$

是否严格正确，需要推导。

## 37.2 时间导数

必须明确：

$$t\longleftrightarrow(n,\tau)$$

以后：

$$\frac{d}{dt}$$

如何处理。

特别需要明确：

- $n$ 是离散慢变量还是连续慢时间；
- $\tau$ 是局部时间还是归一化 phase；
- 周期边界处状态如何传递；
- 不可逆内部变量如何跨周期连续。

## 37.3 Global stage 三变量分离

原论文只需要交替求 space / time functions。

三变量后需考虑：

$$X(x),\quad N(n),\quad T(\tau)$$

的 alternating fixed-point / separated solver。

## 37.4 Update 应更新哪些方向

候选：

- 固定 $X$，联合更新 $N,T$；
- 固定 $X,T$，只更新 $N$；
- $N\leftrightarrow T\leftrightarrow X$ 全方向交替更新。

需要结合 SVD 结果、计算代价和原论文 hybrid strategy 决定。

## 37.5 Add a pair -> Add a triplet

原论文新增：

$$\lambda_{m+1}(t) \bar\varepsilon_{m+1}^p(x).$$

目标扩展：

$$N_{m+1}(n) T_{m+1}(\tau) X_{m+1}(x).$$

这将是后续核心数学问题。

## 37.6 塔筒版 xi

原论文 Eq. (76)-(77) 的 LATIN norm 需要映射到：

- 梁单元；
- 纤维截面；
- 一维材料积分点；
- tower global fields；
- 当前 search directions。

需要保证 $\xi$ 仍严格表示 local/global LATIN distance。

## 37.7 zeta 形式

当前首选保持：

$$\boxed{ \zeta = \frac{ \xi_i-\xi_{i+1} }{ \xi_i+\xi_{i+1} } }$$

真正需要重新确定的是：

$$\boxed{ \zeta^{\rm tol} }$$

而不是重新设计 $\zeta$。

---

# 38. 建议后续路线

## Phase A：理论推导

完成：

$$x-t \rightarrow n-\tau-x$$

三变量 PGD 数学形式。

## Phase B：与原论文逐式对应

重点对应：

```text
Eq. (47)
Eq. (53)
Eq. (58)
Eq. (59)
Eq. (60)
Eq. (61)-(72)
Eq. (76)-(77)
```

## Phase C：最小三变量 PGD prototype

先使用简化问题或冻结 FOM，验证：

$$N(n)T(\tau)X(x)$$

表示和 alternating solver。

## Phase D：Enrichment

继承：

$$\xi \rightarrow \zeta \rightarrow \text{Add mode}.$$

将：

```text
Add a space-time pair
```

扩展为：

```text
Add a cycle-phase-space triplet
```

## Phase E：与 FOM 对比

比较：

- displacement；
- stress；
- plastic strain；
- damage；
- critical fiber；
- cycle evolution；
- computational cost；
- PGD mode count；
- LATIN iterations；
- $\xi$；
- $\zeta$。

## Phase F：真正高周加速

三变量 PGD 稳定后，再研究：

- 更大循环数；
- cycle compression；
- time multiscale；
- cycle jumps；
- hyper-reduction；
- Reference Point Method；
- 真实塔筒钢材标定；
- 工程载荷谱。

---

# 39. 本阶段最终结论

1. **100 周期塔筒 FOM 已稳定运行并冻结为可重复使用的参考数据。**
2. **已建立 snapshot persistence 与全局 cycle slicing。**
3. **位移、应力和累积内部变量具有明显低秩性。**
4. **不可逆增量 Delta eps_p 和 Delta D 的 cycle-space 复杂度更高。**
5. **Stage II 附近存在低复杂度区，约 cycle 30-50 存在明显低秩结构重组过渡带。**
6. **Delta eps_p 与 Delta D 的低秩演化机制不同。**
7. **整数 rank 对渐进复杂度变化不够敏感，singular-value ratio 更适合离线观察。**
8. **SVD/HOSVD 的正确作用是 offline separability / low-rank feasibility analysis。**
9. **原论文已经有成熟的 xi-zeta adaptive enrichment logic，应优先继承。**
10. **真正需要创新和扩展的是从 x-t space-time PGD 到 n-tau-x cycle-phase-space PGD。**
11. **下一阶段应先完成三变量 PGD 数学推导，再进入代码实现。**

研究主线可概括为：

```text
1D LATIN-PGD reproduction
        ↓
tower nonlinear FOM
        ↓
100-cycle frozen reference
        ↓
cycle-phase-space low-rank feasibility
        ↓
x-t  →  n-tau-x PGD extension
        ↓
retain LATIN xi + retain enrichment zeta
        ↓
adaptive cycle-phase-space LATIN-PGD tower solver
```

---

# 40. 下一步

下一步建议正式开展：

> **原论文 x-t PGD 到海上风机塔筒 n-tau-x PGD 的逐式数学推导。**

重点从原论文 Eq. (47)、Eq. (53)、Eq. (58) 开始，逐步推导到：

$$N(n)T(\tau)X(x)$$

并明确：

- correction variables；
- time derivative；
- weak form；
- alternating solver；
- update；
- enrichment；
- $\xi$；
- $\zeta$；
- orthogonalisation；
- rejection。

三变量数学框架明确后，再进入下一轮代码实现。

---

# 参考文献

Bhattacharyya, M., Fau, A., Nackenhorst, U., Néron, D., & Ladevèze, P. (2018).  
**A LATIN-based model reduction approach for the simulation of cycling damage.**  
*Computational Mechanics, 62*, 725-743.  
DOI: 10.1007/s00466-017-1523-z.

---

# 附：本阶段关键文件

```text
examples/nonlinear_tower_snapshot_tensor.py
examples/nonlinear_tower_low_rank_diagnostics.py
examples/nonlinear_tower_100cycle_low_rank_probe.py

examples/nonlinear_tower_stagewise_low_rank_probe.py
examples/nonlinear_tower_equal_window_low_rank_probe.py
examples/nonlinear_tower_sliding_window_low_rank_probe.py

tests/test_nonlinear_tower_snapshot_tensor.py
tests/test_nonlinear_tower_low_rank_diagnostics.py

outputs/tower_100cycle_fom_reference_v1.npz
outputs/tower_100cycle_sliding_window_low_rank_v1.csv
```

建议仓库路径：

```text
docs/2026-08-12-tower-100cycle-low-rank-and-pgd-enrichment-stage-summary.md
```
