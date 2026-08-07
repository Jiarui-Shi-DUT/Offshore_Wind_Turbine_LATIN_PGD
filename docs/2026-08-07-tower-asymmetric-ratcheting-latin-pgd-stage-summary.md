# LATIN-PGD 海上风机塔筒非对称循环、棘轮效应与损伤机制阶段总结

> 日期：2026-08-07
> 项目：Offshore_Wind_Turbine_LATIN_PGD
> 分支：`feature/offshore-wind-turbine-tower-fatigue`
> 本文档起点：上一次阶段总结 `docs/2026-08-07-tower-multicycle-damage-mechanism-stage-summary.md`（提交 `224cd7e`）之后的工作
> 当前远程最新正式提交：`691cada feat: add asymmetric damage mechanism probe`

---

## 0. 本阶段工作的核心目的

上一阶段已经完成了 NREL 5 MW 塔筒非线性纤维梁柱模型、循环加载、多循环诊断以及 `coupled damage` 与 `damage-disabled` 的初步机制分离。但此前主要采用完全反向循环载荷，其平均荷载为零，不适合直接研究非零均值荷载下可能出现的单向塑性累积。

本阶段的核心任务因此转向：

1. 将正式循环基准从完全反向加载扩展为正负不对称循环；
2. 建立与非对称循环相匹配的初始化方式；
3. 定义能够区分“稳定循环塑性”和“棘轮累积”的逐圈诊断量；
4. 在同一网格、同一加载、同一塑性与硬化参数下，对比：
   - 完整粘塑性-损伤耦合模型；
   - 仅关闭损伤演化的模型（`k_damage = 0`）；
5. 通过 20 圈计算判断：
   - 非对称循环是否产生净塑性累积；
   - 该累积是否主要由塑性本身产生；
   - 损伤如何放大位移漂移和塑性漂移；
   - 当前响应是否已经进入稳定持续棘轮、循环稳定化或损伤加速阶段；
6. 最终回答一个对 LATIN-PGD 更根本的问题：

> **塔筒长循环响应究竟具有怎样的时间结构？它是渐近周期、持续漂移，还是“快周期 + 慢演化”的多时间尺度问题？**

这一步不是独立于 LATIN-PGD 的旁支，而是在为后续降阶方法选择正确的时间分解形式。

---

# 1. 正式循环基准由完全反向加载改为非对称循环

## 1.1 采用的正式载荷比

本阶段确定正式非对称循环载荷比：

$$ R_F=\frac{F_{\min}}{F_{\max}}=-0.5 $$

因此：

$$ F_{\min}=-0.5F_{\max} $$

对应平均力：

$$ F_{\mathrm{mean}} =\frac{F_{\max}+F_{\min}}{2} =0.25F_{\max} $$

循环幅值：

$$ F_a =\frac{F_{\max}-F_{\min}}{2} =0.75F_{\max} $$

周期载荷写成：

$$ F(t)=F_{\mathrm{mean}} +F_a\sin\left(\frac{2\pi t}{T}\right) $$

一个完整周期的关键检查点为：

$$ F_{\mathrm{mean}} \rightarrow F_{\max} \rightarrow F_{\mathrm{mean}} \rightarrow F_{\min} \rightarrow F_{\mathrm{mean}} $$

与完全反向加载相比，这种加载具有两个重要特征：

- 仍然经历正负反向加载；
- 平均荷载不为零，因此允许出现循环后同一参考荷载处的累积漂移。

这使其更适合研究非对称塑性、棘轮效应以及损伤-塑性耦合下的长期演化。

---

## 1.2 为什么不能直接从零荷载进入周期历史

由于周期函数在 $t=0$ 时满足：

$$ F(0)=F_{\mathrm{mean}} $$

如果模型仍从零荷载、零状态直接开始，那么第一个周期的初始状态与后续所有周期的起点状态并不一致，会把“从 0 加载到 $F_{\mathrm{mean}}$”的初始化过程混入第一圈。

这会污染逐圈比较，尤其是：

$$ \Delta u_n = u_{\mathrm{end},n}-u_{\mathrm{start},n} $$

和：

$$ \Delta\varepsilon_{p,n} = \varepsilon_{p,\mathrm{end},n} -\varepsilon_{p,\mathrm{start},n} $$

因此本阶段采用显式预加载：

$$ 0\rightarrow F_{\mathrm{mean}} $$

然后再进入正式周期：

$$ F_{\mathrm{mean}} \rightarrow F_{\max} \rightarrow F_{\mathrm{mean}} \rightarrow F_{\min} \rightarrow F_{\mathrm{mean}} $$

预加载规则为：

- 预加载时长：$T/4$；
- 预加载增量数：`increments_per_cycle / 4`；
- 时间步长与正式周期一致；
- 预加载不计入循环数；
- 不保留预加载中间状态；
- 仅把预加载终点保存为周期历史的第 0 个状态。

这样，第 $n$ 圈的起点和终点都严格对应同一个外力：

$$ F=F_{\mathrm{mean}} $$

因此逐圈漂移具有明确物理意义。

---

# 2. 本阶段代码演化与 Git 提交

从上一次阶段总结提交 `224cd7e` 之后，本阶段完成了以下 5 个关键提交。

---

## 2.1 提交 `d6267b5`：add asymmetric cyclic tower loading

新增非对称循环加载模型，核心内容包括：

- `AsymmetricCyclicTopForceHistory`
- `evaluate_asymmetric_cyclic_top_force`
- `create_asymmetric_cyclic_top_force_history`

默认：

```text
force_ratio = -0.5
```

限制：

$$ -1<R_F<0 $$

即当前类专门表示“正最大值 + 负最小值 + 正平均值”的符号反转非对称循环。

同时保留此前：

- `PulsatingTopForceHistory`
- `ReversedCyclicTopForceHistory`

以保证旧算例兼容。

该阶段完成后，加载相关测试全部通过，完整测试数达到 184 项。

---

## 2.2 提交 `7fea44a`：add asymmetric nonlinear tower response

新增：

```text
examples/nonlinear_tower_asymmetric_response.py
```

核心函数：

```text
run_nonlinear_asymmetric_analysis(...)
```

主要完成：

1. 非对称循环正式进入塔筒非线性求解器；
2. 加入显式 `0 -> Fmean` 预加载；
3. 保证预加载不被计为正式循环；
4. 正式周期从 `Fmean` 状态开始；
5. 保留完整材料内部变量状态；
6. 在整个周期历史中自动确定临界纤维位置。

同时对原有：

```text
examples/nonlinear_tower_reversed_response.py
```

做了泛化，使其响应结构可以同时容纳：

- `ReversedCyclicTopForceHistory`
- `AsymmetricCyclicTopForceHistory`

并保留兼容别名：

```text
NonlinearCyclicResponse = NonlinearReversedResponse
```

新增 6 项定向测试后，完整测试数达到 190 项。

---

## 2.3 提交 `b620126`：generalize multicycle diagnostics for asymmetric loading

这一提交非常重要，因为此前多循环诊断默认循环检查点为：

$$ 0\rightarrow +F_a\rightarrow0\rightarrow-F_a\rightarrow0 $$

对非对称加载不再成立。

本阶段把循环诊断统一为：

$$ F_{\mathrm{ref}} \rightarrow F_{\max} \rightarrow F_{\mathrm{ref}} \rightarrow F_{\min} \rightarrow F_{\mathrm{ref}} $$

其中：

- 完全反向循环：$F_{\mathrm{ref}}=0$；
- 非对称循环：$F_{\mathrm{ref}}=F_{\mathrm{mean}}$。

为避免破坏已有代码，保留旧字段名：

```text
first_zero
residual_displacement
```

同时增加更通用的语义：

```text
midpoint_after_positive_peak
cycle_end_displacement
cycle_end_displacements
```

需要特别注意：

> 在非对称循环中，`residual_displacement` 这个旧名字只具有历史兼容意义。其数值实际上表示 **同一 $F_{\mathrm{mean}}$ 下的循环末位移**，不能再称为“零荷载残余位移”。

完成后完整测试数达到 193 项。

---

## 2.4 提交 `df8cfb0`：add asymmetric multicycle ratcheting diagnostics

新增：

```text
examples/nonlinear_tower_asymmetric_multicycle_response.py
```

核心入口：

```text
run_asymmetric_multicycle_tower_analysis(...)
```

正式默认参数：

```text
maximum_force = 1.0e6 N
force_ratio = -0.5
period = 10.0
n_cycles = 5
increments_per_cycle = 40
```

新增三组与棘轮相关的诊断。

### 2.4.1 同力位移漂移

$$ \Delta u_n = u_{\mathrm{end},n} - u_{\mathrm{start},n} $$

起点与终点均为：

$$ F=F_{\mathrm{mean}} $$

因此该量表示每完成一圈后，结构在相同外荷载下发生了多少净位移迁移。

---

### 2.4.2 归一化位移漂移

$$ \eta_{u,n} = \frac{|\Delta u_n|} {\Delta u_{\mathrm{range},n}} $$

用于区分：

- 绝对漂移很大但循环振幅也很大；
- 漂移相对于当前循环幅值确实显著。

---

### 2.4.3 临界纤维塑性应变漂移

$$ \Delta\varepsilon_{p,n} = \varepsilon_{p,\mathrm{end},n} - \varepsilon_{p,\mathrm{start},n} $$

这是判断塑性内部变量是否每圈发生净迁移的核心指标。

该提交新增 10 项测试，完整测试数达到：

```text
203 tests
```

---

## 2.5 提交 `691cada`：add asymmetric damage mechanism probe

这是本阶段最关键的机制分离提交。

新增：

```text
examples/nonlinear_tower_asymmetric_damage_mechanism_probe.py
tests/test_nonlinear_tower_asymmetric_damage_mechanism_probe.py
```

建立严格成对分析：

### Case A：coupled

完整粘塑性-损伤耦合：

$$ k_{\mathrm{damage}}>0 $$

### Case B：damage-disabled

仅关闭损伤：

$$ k_{\mathrm{damage}}=0 $$

除此之外保持完全一致：

- 塔筒几何；
- 网格；
- Gauss 点；
- 纤维离散；
- 非对称循环；
- $R_F$；
- $F_{\max}$；
- 塑性参数；
- 硬化参数；
- 时间离散；
- Newton 设置。

为了保证比较的是同一物理位置：

> 先由 coupled 分析确定临界纤维，再把该纤维位置强制用于 damage-disabled 分析。

因此两组结果的差异可以主要归因于损伤耦合，而不是临界点自动选择不同导致的混淆。

新增诊断：

- coupled 每圈位移漂移；
- damage-disabled 每圈位移漂移；
- 两者差值；
- coupled 每圈塑性应变漂移；
- damage-disabled 每圈塑性应变漂移；
- 两者差值。

新增 9 项定向测试：

```text
Ran 9 tests in 6.299s
OK
```

随后完整回归：

```text
Ran 212 tests in 194.619s
OK
```

并已推送远程：

```text
691cada feat: add asymmetric damage mechanism probe
```

---

# 3. 为什么必须区分“循环塑性活动”和“棘轮效应”

塑性应变在一个循环中发生变化，并不自动意味着发生棘轮。

定义单圈塑性路径长度：

$$ L_{p,n} = \sum_k \left| \Delta\varepsilon_{p,k} \right| $$

以及单圈净塑性增量：

$$ \Delta\varepsilon_{p,n} = \varepsilon_{p,\mathrm{end},n} - \varepsilon_{p,\mathrm{start},n} $$

进一步定义：

$$ r_{p,n} = \frac{ \left| \Delta\varepsilon_{p,n} \right| }{ L_{p,n} } $$

其物理意义为：

### 当 $r_{p,n}\approx0$

说明一圈中塑性变量发生明显往返，但循环末基本回到原位置：

$$ L_{p,n}>0, \qquad |\Delta\varepsilon_{p,n}|\ll L_{p,n} $$

这更接近：

> 稳定循环塑性 / 往复塑性活动

而不是持续棘轮。

### 当 $r_{p,n}\approx1$

说明一圈内绝大多数塑性路径都转化成了同一方向的净累积：

$$ |\Delta\varepsilon_{p,n}| \approx L_{p,n} $$

这是非常强的单向塑性累积信号。

但即使如此，也不能只依据一两圈就宣布存在“持续稳定棘轮”。

真正需要判断的是：

$$ \Delta\varepsilon_{p,n} $$

随着 $n$ 增大后究竟：

1. 趋近于 0；
2. 趋近于非零常数；
3. 先衰减后因损伤重新增大。

---

# 4. 20 圈非对称循环机制探针

在机制框架通过 212 项完整测试后，本阶段运行了一次实际 20 圈计算。

## 4.1 数值配置

```text
NREL 5 MW tower
n_elements = 10
n_gauss = 2
n_circumferential = 16
n_radial = 1

Fmax = 1.0 MN
R_F = -0.5
Fmin = -0.5 MN
Fmean = 0.25 MN

period = 10.0
n_cycles = 20
increments_per_cycle = 40
max_iterations = 40
```

同时运行：

```text
coupled damage
damage-disabled (k_damage = 0)
```

临界位置最终为：

```text
critical location = (0, 0, 4)
critical height = 1.851206 m
critical y = 2.945089 m
```

即位于塔底附近的外缘纤维，符合塔筒弯曲主导下危险位置的基本物理预期。

---

# 5. 20 圈结果的核心结论

## 5.1 当前模型明确存在单向循环塑性累积

第 1 圈：

$$ r_{p,1}^{c}=1.00000 $$

$$ r_{p,1}^{0}=1.00000 $$

第 20 圈：

$$ r_{p,20}^{c}=0.997853 $$

$$ r_{p,20}^{0}=0.999994 $$

其中上标：

- $c$：coupled；
- $0$：damage-disabled。

也就是说，到第 20 圈时：

$$ |\Delta\varepsilon_p| \approx L_p $$

塑性应变在一圈内部几乎没有形成明显的正反往复闭合，而是主要向一个方向持续积累。

因此当前响应不是普通“闭合塑性滞回 + 零净漂移”的状态，而具有非常明确的：

> **ratcheting-like 单向累积特征。**

---

## 5.2 但不能立即称为“稳定持续棘轮”

coupled 位移漂移：

$$ \Delta u_1 = 6.14475\times10^{-2}\ \mathrm{m} $$

到第 20 圈：

$$ \Delta u_{20} = 1.22649\times10^{-2}\ \mathrm{m} $$

下降约：

$$ 80.0\% $$

damage-disabled：

$$ 6.12192\times10^{-2} \rightarrow 1.16193\times10^{-2}\ \mathrm{m} $$

下降约：

$$ 81.0\% $$

临界纤维塑性应变漂移同样显著衰减。

coupled：

$$ 1.38795\times10^{-4} \rightarrow 1.63326\times10^{-5} $$

下降约：

$$ 88.2\% $$

damage-disabled：

$$ 1.38284\times10^{-4} \rightarrow 1.51991\times10^{-5} $$

下降约：

$$ 89.0\% $$

因此当前最准确的结论是：

> **非对称循环下存在显著的单向塑性累积，但每圈净累积速率仍在明显下降。20 圈尚不足以判断其最终是趋近于零、趋近于稳定非零值，还是在损伤进一步发展后重新增大。**

所以目前应使用：

```text
ratcheting-like behavior
单向循环塑性累积
棘轮型漂移
```

而不宜过早使用：

```text
stable sustained ratcheting
稳定持续棘轮
```

---

# 6. 关闭损伤后仍存在明显漂移：棘轮并非由损伤制造

第 20 圈：

coupled：

$$ \Delta u_{20}^{c} = 1.22649\times10^{-2}\ \mathrm{m} $$

damage-disabled：

$$ \Delta u_{20}^{0} = 1.16193\times10^{-2}\ \mathrm{m} $$

并且：

$$ \Delta\varepsilon_{p,20}^{0} = 1.51991\times10^{-5} \neq0 $$

因此即使：

$$ k_{\mathrm{damage}}=0 $$

模型依然具有明显净位移漂移和净塑性应变漂移。

这说明当前机制链条首先是：

$$ \boxed{ \text{非对称循环} \rightarrow \text{塑性单向累积} } $$

而不是：

$$ \boxed{ \text{损伤} \rightarrow \text{制造棘轮} } $$

因此：

> **非零均值的反向循环 + 当前粘塑性硬化规律，本身已经足以产生 ratcheting-like 累积。损伤是在这个基础上进一步放大响应。**

这一机制分离对后续 LATIN-PGD 很重要，因为降阶方法需要同时表示：

- 基础塑性慢漂移；
- 损伤引起的进一步非平稳演化。

---

# 7. 损伤对棘轮的放大效应随循环增强

## 7.1 位移漂移放大

第 1 圈：

$$ \Delta u_1^c-\Delta u_1^0 = 2.28374\times10^{-4}\ \mathrm{m} $$

coupled 相对于 damage-disabled 仅增加约：

$$ 0.37\% $$

第 20 圈：

$$ \Delta u_{20}^c-\Delta u_{20}^0 = 6.45595\times10^{-4}\ \mathrm{m} $$

对应相对放大约：

$$ \boxed{5.56\%} $$

这说明虽然两组每圈漂移都在衰减，但损伤造成的“额外漂移比例”却逐渐变得更明显。

---

## 7.2 塑性应变漂移放大

第 20 圈：

$$ \Delta\varepsilon_{p,20}^{c} = 1.63326\times10^{-5} $$

$$ \Delta\varepsilon_{p,20}^{0} = 1.51991\times10^{-5} $$

因此损伤使第 20 圈净塑性应变漂移提高约：

$$ \boxed{7.46\%} $$

所以当前可以提出如下机制：

$$ \boxed{ \text{非对称循环} \rightarrow \text{塑性净累积} } $$

同时：

$$ \boxed{ D\uparrow \rightarrow K_{\mathrm{eff}}\downarrow \rightarrow \text{结构变形增加} \rightarrow \text{塑性累积被进一步放大} } $$

即：

> **损伤不是棘轮的起因，但会逐步放大棘轮型漂移。**

---

# 8. 当前损伤仍处于缓慢累积，而非加速失稳

最大损伤：

$$ D_1 = 1.96404\times10^{-3} $$

$$ D_{20} = 1.36327\times10^{-2} $$

因此第 20 圈最大损伤约为：

$$ 1.36\% $$

损伤始终单调增加。

但每圈损伤增量：

$$ \Delta D_1 = 1.96404\times10^{-3} $$

到：

$$ \Delta D_{20} = 4.52072\times10^{-4} $$

下降约：

$$ 77.0\% $$

所以当前并不是：

$$ \Delta D_n\uparrow $$

的损伤加速阶段，而是：

$$ \Delta D_n\downarrow $$

的缓慢累积阶段。

因此现阶段不能称为：

```text
damage instability
fatigue runaway
损伤失稳
疲劳加速破坏
```

更准确的表述是：

> **损伤持续增长，但每圈增量仍在衰减，系统尚未进入损伤加速阶段。**

---

# 9. 损伤已经造成可辨识的软化效应

第 20 圈位移范围：

$$ \Delta u_{\mathrm{range},20}^{c} = 1.70309\ \mathrm{m} $$

$$ \Delta u_{\mathrm{range},20}^{0} = 1.69990\ \mathrm{m} $$

coupled 增大约：

$$ +0.188\% $$

临界纤维应力范围：

$$ \Delta\sigma_{20}^{c} = 174.264\ \mathrm{MPa} $$

$$ \Delta\sigma_{20}^{0} = 174.869\ \mathrm{MPa} $$

coupled 降低约：

$$ -0.346\% $$

这形成一致的损伤软化信号：

$$ D\uparrow \Rightarrow K_{\mathrm{eff}}\downarrow $$

在相同外力条件下表现为：

$$ u\uparrow $$

以及局部应力重分配：

$$ \Delta\sigma_{\mathrm{critical}}\downarrow $$

损伤幅值目前还不大，所以这种分叉仍然较弱，但方向已经稳定可辨。

---

# 10. 外功结果支持损伤增加非线性耗散

第 20 圈：

$$ |W|_{20}^{c} = 1.20571\times10^4\ \mathrm{J} $$

$$ |W|_{20}^{0} = 1.13938\times10^4\ \mathrm{J} $$

coupled 高约：

$$ \boxed{5.82\%} $$

因此在同一循环数下：

$$ |W|^{c}>|W|^{0} $$

说明损伤耦合后产生了额外的非线性耗散与路径依赖。

需要注意：

- 两组每圈外功本身都随循环数下降；
- 这与每圈塑性增量逐渐减小是一致的；
- 但 coupled 相对于 damage-disabled 的额外耗散逐渐显现。

---

# 11. 数值收敛情况

20 圈计算中：

```text
max Newton iterations:
coupled = 4
damage-disabled = 4
```

最大残差：

```text
coupled:
8.761245e-03

damage-disabled:
8.270273e-03
```

因此当前观察到的位移漂移、塑性漂移和损伤演化并不是由明显的 Newton 不收敛导致的数值伪影。

这为后续进一步延长循环数提供了一个可靠的全阶参考基础。

---

# 12. 这一问题为什么对后续 LATIN-PGD 至关重要

这一部分是本阶段最重要的方法学认识。

我们真正需要搞清楚的并不是“塔筒有没有棘轮”本身，而是：

> **长循环全阶响应的时间结构究竟是什么。**

因为这直接决定 PGD 是否低秩、怎样分离时间变量，以及是否需要自适应 enrichment。

---

## 12.1 情况 A：每圈漂移最终趋于 0

假设更长循环后：

$$ \Delta u_n\rightarrow0 $$

$$ \Delta\varepsilon_{p,n}\rightarrow0 $$

则说明系统逐渐发生：

```text
shakedown
cyclic stabilization
循环稳定化
```

此时后期响应近似满足：

$$ \mathbf u(\mathbf x,t+T) \approx \mathbf u(\mathbf x,t) $$

如果损伤也演化得很慢，则整个长时间历史具有很强的重复结构。

PGD 表示：

$$ \mathbf u(\mathbf x,t) \approx \sum_{i=1}^{r} \mathbf U_i(\mathbf x)\lambda_i(t) $$

很可能只需要较小的秩：

$$ r\ll N_t $$

这属于 LATIN-PGD 最容易获得显著降阶收益的情况。

本质上，大量循环共享相近的空间响应，只需少数空间模态配合时间系数描述。

---

## 12.2 情况 B：每圈漂移趋于稳定非零值

若：

$$ \Delta u_n\rightarrow c_u\neq0 $$

$$ \Delta\varepsilon_{p,n} \rightarrow c_p\neq0 $$

则存在持续棘轮。

此时：

$$ u(t+T)\neq u(t) $$

即使每一圈的局部形状相似，整条响应仍在不断迁移。

这说明“单一时间轴上的简单周期 PGD”未必是最合理的表示。

更自然的结构是引入：

- 循环数 $n$：慢时间；
- 圈内相位 $\tau$：快时间。

写成：

$$ u(\mathbf x,n,\tau) = \bar u(\mathbf x,n) + \widetilde u(\mathbf x,n,\tau) $$

其中：

$$ \bar u(\mathbf x,n) $$

描述跨循环慢速漂移，

而：

$$ \widetilde u(\mathbf x,n,\tau) $$

描述单圈内的快速周期响应。

于是 LATIN-PGD 的目标就不再只是：

$$ \mathbf x\times t $$

分离，而可能演化成：

$$ \boxed{ \mathbf x \times n \times \tau } $$

或等价的慢—快时间结构。

---

## 12.3 情况 C：前期稳定化，后期因损伤重新加速

第三种情况对疲劳问题尤其重要：

$$ \Delta u_n\downarrow $$

前期塑性逐渐稳定；

但随着：

$$ D_n\uparrow $$

刚度退化逐渐增强，后期可能出现：

$$ \Delta u_n\uparrow $$

甚至：

$$ \Delta D_n\uparrow $$

此时长循环历史可以分为：

$$ \boxed{ \text{初始调整} \rightarrow \text{近稳定阶段} \rightarrow \text{损伤加速阶段} } $$

这对 PGD 有一个非常重要的含义：

> 早期和中期可能用很少的空间模态就能表示，但损伤局部化或加速后，原有空间基可能不再充分。

于是可能需要：

$$ r=r(N) $$

即 PGD 秩随疲劳演化增加。

这会自然引出：

```text
adaptive enrichment
自适应模态增广
```

而不是预先固定一个永远不变的低秩空间基。

---

# 13. 当前最值得发展的 LATIN-PGD 时间结构：快—慢分离

经过本阶段分析，后续真正值得发展的方向已经逐渐清楚。

循环疲劳问题天然包含两个时间尺度：

## 快时间

一圈内部：

$$ \tau\in[0,T] $$

描述：

- 正向加载；
- 峰值；
- 卸载；
- 反向加载；
- 再卸载。

它主要决定单圈滞回形状。

---

## 慢时间

循环数：

$$ N=1,2,\ldots $$

描述：

- 塑性净累积；
- 硬化变量演化；
- 损伤增长；
- 刚度退化；
- 临界位置变化；
- 可能的损伤局部化。

因此后续 LATIN-PGD 最有方法学意义的结构不是简单“把论文 1D 杆换成塔筒”，而是研究：

$$ \boxed{ \text{LATIN} + \text{PGD} + \text{循环快时间} + \text{疲劳慢时间} } $$

这比单纯做更复杂几何的算例更具有博士课题层面的扩展价值。

---

# 14. 为什么必须先识别长期行为，再设计 PGD

如果现在未经判断就直接假定：

$$ \text{每一圈近似重复} $$

然后构造 PGD，很可能出现一种危险情况：

单个循环内：

$$ F-u $$

曲线和全阶解非常接近；

但是每圈塑性漂移存在一个很小误差：

$$ \delta(\Delta\varepsilon_p) $$

经过大量循环后：

$$ \sum_{n=1}^{N} \delta(\Delta\varepsilon_{p,n}) $$

可能逐渐累积。

最终导致：

$$ \varepsilon_p(N) $$

错误，进而：

$$ D(N) $$

错误，再进一步使：

$$ u(N) $$

和疲劳寿命判断发生明显偏差。

所以未来 LATIN-PGD 验证不能只比较：

```text
某一圈滞回环
某一时刻位移
某一时刻应力
```

还必须比较慢变量演化。

至少包括：

$$ \boxed{ \Delta u_n } $$

$$ \boxed{ \Delta\varepsilon_{p,n} } $$

$$ \boxed{ D_{\max}(n) } $$

以及最终真正重要的空间场：

$$ \boxed{ D(\mathbf x,N) } $$

---

# 15. 对 LATIN-PGD 降阶质量评价指标的启示

后续全阶模型与 LATIN-PGD 的对比应至少分成四层。

## 第一层：单圈快时间响应

比较：

- 塔顶位移历史；
- 临界应力历史；
- 单圈滞回形状；
- 外功或耗散。

用于判断 PGD 是否正确捕捉单圈动力学/准静力非线性路径。

---

## 第二层：跨循环慢漂移

比较：

$$ \Delta u_n $$

$$ \Delta\varepsilon_{p,n} $$

以及：

$$ \varepsilon_{p,\mathrm{end}}(n) $$

用于判断 PGD 是否正确保持路径依赖和累积效应。

---

## 第三层：损伤演化

比较：

$$ D_{\max}(n) $$

$$ \Delta D_n $$

以及可能的累积耗散。

用于判断长期疲劳损伤是否被正确预测。

---

## 第四层：损伤空间场

不能只看单个临界纤维。

最终必须比较：

$$ D(\mathbf x,N) $$

或离散形式：

$$ D_{e,g,f}(N) $$

其中：

- $e$：梁单元；
- $g$：Gauss 点；
- $f$：截面纤维。

这是未来判断 LATIN-PGD 是否真正能够用于风机塔筒疲劳分析的关键指标。

---

# 16. 当前已经明确的科学认识

截至本阶段，可以认为以下几点已经比较明确。

### 16.1 非对称循环框架已经建立

正式基准：

$$ R_F=-0.5 $$

并具有一致的预加载与周期定义。

### 16.2 同力逐圈漂移定义已经建立

每圈起点和终点均位于：

$$ F=F_{\mathrm{mean}} $$

因此位移漂移和塑性漂移具有清晰可比性。

### 16.3 当前模型存在显著单向塑性累积

证据：

$$ r_p\approx1 $$

且：

$$ \Delta u_n>0 $$

$$ \Delta\varepsilon_{p,n}>0 $$

连续 20 圈保持同号。

### 16.4 棘轮型漂移不依赖损伤才能出现

即使：

$$ k_{\mathrm{damage}}=0 $$

仍然存在明显：

$$ \Delta u_n $$

和：

$$ \Delta\varepsilon_{p,n} $$

因此基础机制主要来自非对称循环下的塑性演化。

### 16.5 损伤会逐渐放大漂移

第 20 圈：

- 位移漂移放大约 5.56%；
- 塑性应变漂移放大约 7.46%。

### 16.6 当前损伤尚未失稳

虽然：

$$ D_n\uparrow $$

但：

$$ \Delta D_n\downarrow $$

尚未出现损伤增量加速。

### 16.7 当前仍无法判断长期极限状态

20 圈内：

$$ \Delta u_n\downarrow $$

$$ \Delta\varepsilon_{p,n}\downarrow $$

但尚未接近足够明确的渐近状态。

因此不能判断最终是：

- shakedown；
- 非零稳定棘轮；
- 损伤驱动的后期再加速。

---

# 17. 当前仍未解决的问题

以下问题应明确保留，避免后续重新进入项目时误以为已经解决。

## 17.1 20 圈不足以证明持续棘轮

当前只能说：

```text
ratcheting-like
单向循环塑性累积
```

不能严格说：

```text
stable sustained ratcheting
```

需要更长循环观察。

---

## 17.2 当前尚未研究长期渐近值

需要判断：

$$ \lim_{n\rightarrow\infty}\Delta u_n $$

和：

$$ \lim_{n\rightarrow\infty} \Delta\varepsilon_{p,n} $$

究竟是否为 0。

---

## 17.3 预加载敏感性尚未系统检查

当前采用：

$$ 0\rightarrow F_{\mathrm{mean}} $$

作为控制初始化。

其物理和数值逻辑合理，但后续正式论文阶段最好做一次初始化敏感性检查，以确认长期结论并非由预加载方式主导。

---

## 17.4 空间离散收敛尚未完成

当前 20 圈机制探针采用：

```text
10 elements
2 Gauss/element
16 circumferential fibers
1 radial layer
```

其目的是机制探索，并非最终空间收敛结果。

正式参考模型候选仍应向此前规划的更细离散推进，例如：

```text
40 elements
4 Gauss/element
32 circumferential x 2 radial fibers
```

但不应在机制尚未搞清楚之前直接把所有计算都升级到最细模型。

---

## 17.5 尚未提取全塔损伤场

目前重点仍是临界纤维诊断。

后续必须扩展到：

$$ D_{e,g,f}(N) $$

以便真正建立 LATIN-PGD 的空间降阶目标。

---

# 18. 下一阶段推荐路线

下一阶段不应马上进入 LATIN-PGD 编码，而应先把全阶参考问题的长期行为分类清楚。

推荐顺序如下。

---

## Step 1：延长当前完全相同算例的循环数

保持：

```text
10 elements
2 Gauss/element
16 circumferential fibers
R = -0.5
Fmax = 1.0 MN
40 increments/cycle
```

仅改变：

```text
n_cycles
```

优先观察：

```text
50 cycles
100 cycles
```

目的不是追求“高周”数量级，而是识别趋势。

---

## Step 2：重点观察后期漂移

关注：

$$ \Delta u_n $$

$$ \Delta\varepsilon_{p,n} $$

$$ \Delta D_n $$

特别是后 10～20 圈。

可以检查：

1. 是否保持同号；
2. 是否持续单调下降；
3. 是否趋于数值噪声；
4. 是否趋于稳定非零平台；
5. 是否因损伤发展重新增大。

---

## Step 3：形成长期行为分类

最终将全阶问题划分为：

### A. 渐近周期稳定

$$ \Delta u_n\to0 $$

### B. 稳定棘轮

$$ \Delta u_n\to c\neq0 $$

### C. 损伤驱动非平稳演化

$$ \Delta u_n $$

或：

$$ \Delta D_n $$

后期再次加速。

只有完成这一步，才应该正式确定 LATIN-PGD 的时间表示形式。

---

## Step 4：再设计 LATIN-PGD 时间分离

若 A：

> 优先利用强周期重复性，低秩 PGD 很可能直接有效。

若 B：

> 引入慢漂移 + 快周期的时间结构。

若 C：

> 需要慢—快时间分离 + 自适应 enrichment，且可能需要按损伤阶段更新空间基。

---

## Step 5：扩展到损伤空间场

在进入正式降阶前，应增加：

```text
element × Gauss × fiber × cycle
```

层面的损伤场输出。

后续 LATIN-PGD 的真正比较对象不应只是一条临界纤维历史，而应包括：

$$ D(\mathbf x,N) $$

以及可能的：

$$ \varepsilon_p(\mathbf x,N) $$

---

# 19. 对本课题方法学定位的最新认识

此前最直接的扩展思路是：

> 将论文中的 LATIN-PGD 一维杆算例推广到海上风机塔筒纤维梁柱模型。

经过本阶段工作，这个目标已经可以进一步提升。

更值得发展的研究问题是：

> **针对具有周期快响应和疲劳慢演化的海上风机塔筒非线性循环问题，研究 LATIN-PGD 如何利用循环重复性，同时准确表示塑性漂移、损伤累积和空间场演化。**

换句话说，研究重点不应只是“几何从杆变成塔筒”，而应转向：

$$ \boxed{ \text{多时间尺度循环疲劳问题中的 LATIN-PGD 降阶} } $$

其潜在创新点包括：

1. 快时间 $\tau$ 与慢循环变量 $N$ 的分离；
2. 塑性内部变量与损伤内部变量的慢演化表示；
3. 损伤发展过程中 PGD 秩的演化；
4. 损伤局部化时的自适应 enrichment；
5. 长循环累积误差控制；
6. 以疲劳损伤空间场而非单一位移响应作为降阶准确性评价对象。

这比单纯复现原论文框架并替换结构模型，更具有独立的方法学意义。

---

# 20. 当前仓库状态说明

截至本阶段正式提交：

```text
691cada feat: add asymmetric damage mechanism probe
```

已经推送到：

```text
origin/feature/offshore-wind-turbine-tower-fatigue
```

推送后：

```text
git status --short --branch
```

显示本地分支与远程同步，正式代码工作区干净。

随后使用的：

```text
run_asymmetric_20cycle_probe.py
```

属于数值探针脚本。

当前建议：

> 在确认长期分析输出格式之前，不急于把该临时脚本作为正式仓库接口提交。等 50/100 圈长期趋势指标确定后，再整理成正式的 long-cycle benchmark / report 脚本。

---

# 21. 阶段结束时必须记住的核心结论

以后重新进入该项目时，优先恢复以下 10 点：

1. 正式循环已经从完全反向加载改为：
   $$ R_F=-0.5 $$

2. 非对称循环前采用：
   $$ 0\rightarrow F_{\mathrm{mean}} $$
   显式预加载，预加载不计入循环数。

3. 每圈起点与终点都在：
   $$ F=F_{\mathrm{mean}} $$
   因此可以定义严格同力漂移。

4. 位移棘轮指标：
   $$ \Delta u_n = u_{\mathrm{end},n} - u_{\mathrm{start},n} $$

5. 塑性棘轮指标：
   $$ \Delta\varepsilon_{p,n} = \varepsilon_{p,\mathrm{end},n} - \varepsilon_{p,\mathrm{start},n} $$

6. 20 圈内：
   $$ |\Delta\varepsilon_p|/L_p\approx1 $$
   表明当前临界纤维具有非常明显的单向塑性累积。

7. `k_damage=0` 后仍存在明显漂移，因此：
   > 基础棘轮来自非对称塑性，而不是损伤制造。

8. 损伤会逐渐放大棘轮型漂移；第 20 圈约表现为：
   - 位移漂移 +5.56%；
   - 塑性应变漂移 +7.46%。

9. 但每圈漂移和每圈损伤增量都还在下降，所以：
   > 20 圈尚不能证明稳定持续棘轮，也没有进入损伤加速失稳。

10. 对 LATIN-PGD 最重要的问题已经变成：
    > **长期响应最终是渐近周期、稳定非零漂移，还是损伤驱动的慢变非平稳过程？**

这个答案决定后续 LATIN-PGD 是采用简单时间低秩分解，还是发展：

$$ \boxed{ \text{空间} \times \text{循环慢时间} \times \text{单圈快时间} } $$

的多时间尺度表示，并决定是否需要随损伤演化进行自适应 PGD enrichment。

---

# 附录 A：20 圈主要逐圈结果

## A.1 位移漂移、塑性漂移与损伤

<table>
  <thead>
    <tr>
      <th>Cycle</th>
      <th>Δu_c (m)</th>
      <th>Δu_0 (m)</th>
      <th>差值 (m)</th>
      <th>Δε_p,c</th>
      <th>Δε_p,0</th>
      <th>差值</th>
      <th>D_c,end</th>
      <th>ΔD_c</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>1</td>
      <td>6.14475e-02</td>
      <td>6.12192e-02</td>
      <td>2.28374e-04</td>
      <td>1.38795e-04</td>
      <td>1.38284e-04</td>
      <td>5.11312e-07</td>
      <td>1.96404e-03</td>
      <td>1.96404e-03</td>
    </tr>
    <tr>
      <td>2</td>
      <td>4.17216e-02</td>
      <td>4.13904e-02</td>
      <td>3.31272e-04</td>
      <td>7.43453e-05</td>
      <td>7.36014e-05</td>
      <td>7.43860e-07</td>
      <td>3.14193e-03</td>
      <td>1.17789e-03</td>
    </tr>
    <tr>
      <td>3</td>
      <td>3.32861e-02</td>
      <td>3.28901e-02</td>
      <td>3.95950e-04</td>
      <td>5.32967e-05</td>
      <td>5.24797e-05</td>
      <td>8.17050e-07</td>
      <td>4.05299e-03</td>
      <td>9.11056e-04</td>
    </tr>
    <tr>
      <td>4</td>
      <td>2.88328e-02</td>
      <td>2.83793e-02</td>
      <td>4.53460e-04</td>
      <td>4.33614e-05</td>
      <td>4.24616e-05</td>
      <td>8.99841e-07</td>
      <td>4.83729e-03</td>
      <td>7.84300e-04</td>
    </tr>
    <tr>
      <td>5</td>
      <td>2.60230e-02</td>
      <td>2.55222e-02</td>
      <td>5.00830e-04</td>
      <td>3.77323e-05</td>
      <td>3.67521e-05</td>
      <td>9.80171e-07</td>
      <td>5.55115e-03</td>
      <td>7.13866e-04</td>
    </tr>
    <tr>
      <td>6</td>
      <td>2.39811e-02</td>
      <td>2.34431e-02</td>
      <td>5.37996e-04</td>
      <td>3.40518e-05</td>
      <td>3.30048e-05</td>
      <td>1.04706e-06</td>
      <td>6.22038e-03</td>
      <td>6.69225e-04</td>
    </tr>
    <tr>
      <td>7</td>
      <td>2.23564e-02</td>
      <td>2.17896e-02</td>
      <td>5.66766e-04</td>
      <td>3.13561e-05</td>
      <td>3.02581e-05</td>
      <td>1.09798e-06</td>
      <td>6.85775e-03</td>
      <td>6.37367e-04</td>
    </tr>
    <tr>
      <td>8</td>
      <td>2.09931e-02</td>
      <td>2.04040e-02</td>
      <td>5.89034e-04</td>
      <td>2.92131e-05</td>
      <td>2.80785e-05</td>
      <td>1.13465e-06</td>
      <td>7.47011e-03</td>
      <td>6.12368e-04</td>
    </tr>
    <tr>
      <td>9</td>
      <td>1.98123e-02</td>
      <td>1.92060e-02</td>
      <td>6.06264e-04</td>
      <td>2.74151e-05</td>
      <td>2.62552e-05</td>
      <td>1.15985e-06</td>
      <td>8.06151e-03</td>
      <td>5.91397e-04</td>
    </tr>
    <tr>
      <td>10</td>
      <td>1.87690e-02</td>
      <td>1.81495e-02</td>
      <td>6.19518e-04</td>
      <td>2.58542e-05</td>
      <td>2.46780e-05</td>
      <td>1.17624e-06</td>
      <td>8.63456e-03</td>
      <td>5.73044e-04</td>
    </tr>
    <tr>
      <td>11</td>
      <td>1.78347e-02</td>
      <td>1.72051e-02</td>
      <td>6.29583e-04</td>
      <td>2.44702e-05</td>
      <td>2.32842e-05</td>
      <td>1.18601e-06</td>
      <td>9.19111e-03</td>
      <td>5.56551e-04</td>
    </tr>
    <tr>
      <td>12</td>
      <td>1.69899e-02</td>
      <td>1.63528e-02</td>
      <td>6.37075e-04</td>
      <td>2.32259e-05</td>
      <td>2.20351e-05</td>
      <td>1.19078e-06</td>
      <td>9.73259e-03</td>
      <td>5.41483e-04</td>
    </tr>
    <tr>
      <td>13</td>
      <td>1.62205e-02</td>
      <td>1.55781e-02</td>
      <td>6.42474e-04</td>
      <td>2.20968e-05</td>
      <td>2.09050e-05</td>
      <td>1.19180e-06</td>
      <td>1.02602e-02</td>
      <td>5.27570e-04</td>
    </tr>
    <tr>
      <td>14</td>
      <td>1.55158e-02</td>
      <td>1.48696e-02</td>
      <td>6.46169e-04</td>
      <td>2.10655e-05</td>
      <td>1.98755e-05</td>
      <td>1.19002e-06</td>
      <td>1.07748e-02</td>
      <td>5.14633e-04</td>
    </tr>
    <tr>
      <td>15</td>
      <td>1.48673e-02</td>
      <td>1.42188e-02</td>
      <td>6.48465e-04</td>
      <td>2.01186e-05</td>
      <td>1.89327e-05</td>
      <td>1.18598e-06</td>
      <td>1.12773e-02</td>
      <td>5.02543e-04</td>
    </tr>
    <tr>
      <td>16</td>
      <td>1.42682e-02</td>
      <td>1.36186e-02</td>
      <td>6.49603e-04</td>
      <td>1.92456e-05</td>
      <td>1.80657e-05</td>
      <td>1.17992e-06</td>
      <td>1.17685e-02</td>
      <td>4.91207e-04</td>
    </tr>
    <tr>
      <td>17</td>
      <td>1.37128e-02</td>
      <td>1.30630e-02</td>
      <td>6.49763e-04</td>
      <td>1.84376e-05</td>
      <td>1.72657e-05</td>
      <td>1.17186e-06</td>
      <td>1.22491e-02</td>
      <td>4.80546e-04</td>
    </tr>
    <tr>
      <td>18</td>
      <td>1.31964e-02</td>
      <td>1.25473e-02</td>
      <td>6.49074e-04</td>
      <td>1.76867e-05</td>
      <td>1.65253e-05</td>
      <td>1.16140e-06</td>
      <td>1.27196e-02</td>
      <td>4.70500e-04</td>
    </tr>
    <tr>
      <td>19</td>
      <td>1.27149e-02</td>
      <td>1.20672e-02</td>
      <td>6.47653e-04</td>
      <td>1.69868e-05</td>
      <td>1.58383e-05</td>
      <td>1.14856e-06</td>
      <td>1.31806e-02</td>
      <td>4.61021e-04</td>
    </tr>
    <tr>
      <td>20</td>
      <td>1.22649e-02</td>
      <td>1.16193e-02</td>
      <td>6.45595e-04</td>
      <td>1.63326e-05</td>
      <td>1.51991e-05</td>
      <td>1.13353e-06</td>
      <td>1.36327e-02</td>
      <td>4.52072e-04</td>
    </tr>
  </tbody>
</table>

---

## A.2 第 20 圈 paired comparison 摘要

<table>
  <thead>
    <tr>
      <th>指标</th>
      <th>Coupled</th>
      <th>Damage-disabled</th>
      <th>Coupled 相对变化</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>每圈位移漂移</td>
      <td>1.22649e-02 m</td>
      <td>1.16193e-02 m</td>
      <td>+5.56%</td>
    </tr>
    <tr>
      <td>每圈塑性应变漂移</td>
      <td>1.63326e-05</td>
      <td>1.51991e-05</td>
      <td>+7.46%</td>
    </tr>
    <tr>
      <td>位移范围</td>
      <td>1.70309 m</td>
      <td>1.69990 m</td>
      <td>+0.188%</td>
    </tr>
    <tr>
      <td>临界应力范围</td>
      <td>174.264 MPa</td>
      <td>174.869 MPa</td>
      <td>-0.346%</td>
    </tr>
    <tr>
      <td>外功绝对值</td>
      <td>1.20571e+04 J</td>
      <td>1.13938e+04 J</td>
      <td>+5.82%</td>
    </tr>
    <tr>
      <td>最大损伤</td>
      <td>1.36327e-02</td>
      <td>0</td>
      <td>—</td>
    </tr>
    <tr>
      <td>net/path</td>
      <td>0.997853</td>
      <td>0.999994</td>
      <td>均接近 1</td>
    </tr>
  </tbody>
</table>

---

# 附录 B：本阶段关键提交顺序

```text
224cd7e docs: summarize multicycle damage mechanism stage
    ↓
d6267b5 feat: add asymmetric cyclic tower loading
    ↓
7fea44a feat: add asymmetric nonlinear tower response
    ↓
b620126 feat: generalize multicycle diagnostics for asymmetric loading
    ↓
df8cfb0 feat: add asymmetric multicycle ratcheting diagnostics
    ↓
691cada feat: add asymmetric damage mechanism probe
```

完整测试演化：

```text
184
→ 190
→ 193
→ 203
→ 212 tests
```

最终：

```text
Ran 212 tests in 194.619s
OK
```
