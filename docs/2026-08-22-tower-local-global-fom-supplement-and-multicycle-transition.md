# Offshore Wind Turbine Tower LATIN-PGD
# Local / Global–FOM 补充诊断、时间步敏感性与多循环验证转段总结

日期：2026-08-22  
分支：`feature/offshore-wind-turbine-tower-fatigue`  
上一阶段总结文档：`docs/2026-08-21-tower-fom-temporal-consistency-stage-summary.md`  

> 本文档总结自上一份 FOM-2 时间离散一致性阶段总结之后所进行的全部关键推导、数值诊断、方法讨论与研究路线调整。
>
> 本阶段没有修改 `latin/` 核心算法。主要工作是进一步回答三个问题：
>
> 1. Local RK4 与 Global Backward Euler 得到的内部变量究竟相差多少；
> 2. Local 与 Global 中，谁更接近逐时间步 FOM；
> 3. 在 40 increments/cycle 下误差已约为 5% 量级后，后续研究是否应从单循环误差追踪转向多循环误差累积与计算效率评估。
>
> 本文档严格区分：
>
> - 已被当前代码和数值结果直接验证的事实；
> - 根据现有结果得到的数值解释；
> - 尚未解决的理论问题；
> - 下一阶段拟采用的研究策略。
>
> 为避免 PyCharm Markdown 显示问题，本文档使用普通 Unicode 和代码块，不使用 LaTeX 行内公式。

---

# 1. 上一阶段结束时的状态

上一阶段 FOM-2 已完成：

- FOM-2A：critical-point damage-chain diagnostic；
- FOM-2B：critical-point stress-error decomposition；
- FOM-2C：elastic-strain error cancellation；
- FOM-2D：matched FOM–LATIN temporal refinement。

当时已经得到一个很强的数值结论：

```text
40 → 80 → 160 increments/cycle
```

时，FOM–LATIN 的：

```text
total strain
elastic strain
stress
plastic strain
alpha
r_bar
damage
```

误差均随 `dt` 减半而近似减半，观测阶约：

```text
p ≈ 0.93 ~ 0.98
```

因此：

> 当前 one-cycle benchmark 上，FOM–LATIN discrepancy 主要表现为有限时间步导致的近一阶时间离散误差，而没有观察到一个在 `dt → 0` 后仍然保持不变的 persistent discrepancy。

但是上一阶段仍留下两个理论问题：

```text
Local RK4 / Global BE 是否应当统一？
```

以及：

```text
当前 Global Backward Euler
与原论文 Eq. (59) 所述 DG0 时间离散
到底是什么严格关系？
```

这两个问题在本阶段仍未被理论上解决。

---

# 2. 对当前 Local / Global 时间积分结构的重新梳理

本阶段首先重新讨论了：

> 为什么当前程序会采用 Local RK4 与 Global Backward Euler，而不是从一开始就采用完全统一的时间积分格式？

当前结构可以概括为：

```text
Local stage
=
局部材料本构演化积分

Global stage
=
整个时间域上的平衡 / 相容 / PGD global correction
```

Local stage 中需要顺序积分：

```text
plastic strain      eps_p
alpha
r_bar
damage              D
```

当前 tower local stage 与已经验证的 1D local stage 使用同一套 RK4 material-point integrator。

Global stage 中：

```text
eps_p
alpha
r_bar
```

的 integrated history 采用 Backward-Euler 型离散关系进行重构。

当前 tower v1 中：

```text
damage
```

是一个特殊情况：

> Global 并没有像 `eps_p / alpha / r_bar` 一样独立用 BE 完整重新积分 damage history，而是基本继承 Local damage history。

因此当前实际情况并不是简单的：

```text
所有变量：
Local RK4
vs
Global BE
```

而更准确地说是：

```text
eps_p      Local RK4  vs  Global BE-like history
alpha      Local RK4  vs  Global BE-like history
r_bar      Local RK4  vs  Global BE-like history

damage     Local RK4  vs  Global primarily inherits Local history
```

---

# 3. 为什么之前没有强制 Local / Global 使用同一种时间积分格式

早期开发时，Local 和 Global 被视为两个数学职责不同的子问题。

Local 关注：

```text
在给定 thermodynamic-force / search-direction 信息后
准确积分非线性材料演化
```

因此使用 RK4 有利于：

- 提高局部 constitutive integration 精度；
- 处理塑性、硬化和损伤等非线性演化；
- 复用前期已验证的 1D material-point 实现。

Global 关注：

```text
在整个时间域上满足：
equilibrium
+
compatibility
+
PGD separated representation
```

Global temporal problem 中仍然需要离散：

```text
dx/dt
```

因此早期实现选择了简单、稳定的一阶 Backward-Euler 型 temporal derivative。

所以当初的设计逻辑更接近：

```text
Local：
选择适合本构 ODE 的积分器

Global：
选择适合 PGD temporal equations 的时间离散
```

而不是：

```text
同一个 ODE
故意同时采用两种不同算法
```

但是，这种工程实现并没有自动保证：

> Local 与 Global 在有限 `dt` 下对 integrated history 完全离散一致。

---

# 4. 已经发现的 Local / Global history closure 特征

上一阶段已经得到：

在：

```text
one cycle
40 increments/cycle
tol = 1e-5
residual_ls
```

下，Local 与 Global 的 rate 已经非常接近：

```text
plastic strain rate
≈ 2.67505e-4 relative difference

alpha rate
≈ 6.31677e-5

r_bar rate
≈ 4.73829e-5

damage rate
≈ 3.59433e-5
```

但 integrated histories 明显不同：

```text
plastic strain history
≈ 8.47 %

alpha history
≈ 8.53 %

r_bar history
≈ 5.39 %

damage history
≈ 2.63e-5
```

这说明：

```text
Local 和 Global 对“当前应该怎么变化”
判断已经非常接近

但

对“如何把这些变化累计成 history”
采用了不同的离散规则
```

从而造成：

```text
rate close
but
history different
```

---

# 5. 本阶段新的核心问题：Local 与 FOM 到底有多接近

此前正式 FOM accuracy comparison 一直采用：

```text
result.state
```

即最终 accepted / committed Global state。

这在 LATIN 框架中是正确的，因为：

```text
Global state
```

满足结构整体：

```text
equilibrium
+
compatibility
```

而：

```text
Local state
```

只是局部 constitutively admissible state，并不是最终 globally admissible structural solution。

但是，为了判断：

> Global history reconstruction 是否是内部变量误差的主要来源，

必须补做：

```text
FOM
vs
Local
vs
Global
```

三方对比。

---

# 6. Local–FOM–Global 三方诊断

本阶段建立了 40-increment one-cycle 三方诊断：

```text
FOM
vs
result.local_state
vs
result.state
```

计算条件：

```text
loading                 = fully reversed
cycles                   = 1
increments/cycle        = 40
dt                       = 0.25 s
LATIN tolerance          = 1e-5
spatial strategy         = residual_ls
termination              = converged
committed iterations     = 20
PGD rank                 = 15
final xi                 = 8.006194456505e-06
```

FOM 缓存：

```text
outputs/tower_1cycle_reversed_fom_reference_v1.npz
```

字段：

```text
fiber_strains
fiber_stresses
fiber_states
```

---

# 7. Local–FOM–Global：全场 relative L2

结果：

| 场变量 | Local / FOM | Global / FOM | Local / Global | Local/FOM ÷ Global/FOM |
|---|---:|---:|---:|---:|
| total strain | 3.558346576e-03 | 3.575144269e-03 | 5.115414525e-03 | 0.995302 |
| elastic strain | 2.873920107e-03 | 2.873943696e-03 | 4.702505892e-06 | 0.999992 |
| stress | 2.871167752e-03 | 2.871192934e-03 | 4.703832963e-06 | 0.999991 |
| plastic strain | 3.934805197e-02 | 7.297194761e-02 | 8.253102463e-02 | 0.539222 |
| alpha | 4.164895720e-02 | 7.351065984e-02 | 8.305736395e-02 | 0.566570 |
| r_bar | 3.957700223e-02 | 3.463837253e-02 | 5.175123921e-02 | 1.142577 |
| damage | 2.973214467e-02 | 2.973622884e-02 | 1.074355119e-05 | 0.999863 |

换成百分数：

```text
total strain:
Local   ≈ 0.356 %
Global  ≈ 0.358 %

elastic strain:
Local   ≈ 0.2874 %
Global  ≈ 0.2874 %

stress:
Local   ≈ 0.2871 %
Global  ≈ 0.2871 %

plastic strain:
Local   ≈ 3.935 %
Global  ≈ 7.297 %

alpha:
Local   ≈ 4.165 %
Global  ≈ 7.351 %

r_bar:
Local   ≈ 3.958 %
Global  ≈ 3.464 %

damage:
Local   ≈ 2.973 %
Global  ≈ 2.974 %
```

---

# 8. 三方诊断得到的第一个重要结论：eps_p 和 alpha 的 Global history 明显更差

Plastic strain：

```text
Local / FOM
=
3.9348 %

Global / FOM
=
7.2972 %
```

即使用 Local history 后，误差相对 Global 降低约：

```text
46 %
```

Alpha：

```text
Local / FOM
=
4.1649 %

Global / FOM
=
7.3511 %
```

误差相对 Global 降低约：

```text
43 %
```

因此可以确认：

> 对于 `eps_p` 和 `alpha`，当前 Global temporal history reconstruction 确实引入了显著的额外 FOM discrepancy。

但是：

```text
Local eps_p / FOM
≈ 3.93 %

Local alpha / FOM
≈ 4.16 %
```

仍然不为零。

因此不能说：

```text
全部约 7 % 的 eps_p / alpha error
=
Global BE 单独造成
```

更准确的分层解释是：

```text
Local constitutive / search path 本身
已经存在约 4 % 的 FOM discrepancy

然后

Global temporal history reconstruction
进一步使 eps_p / alpha error
扩大到约 7 %
```

注意：

> relative L2 norm 不能简单线性相加，因此这里描述的是误差来源层次，而不是严格的百分数加法分解。

---

# 9. r_bar 的行为与 eps_p / alpha 不同

Full-field：

```text
r_bar Local / FOM
=
3.9577 %

r_bar Global / FOM
=
3.4638 %
```

所以：

> 对 `r_bar` 而言，Global 反而比 Local 更接近 FOM。

这说明不能把当前问题简单总结为：

```text
RK4 always better
BE always worse
```

因为 Global stage 并不是一个孤立的 BE material integrator。

它同时包含：

```text
search-direction correction
+
global equilibrium
+
compatibility
+
PGD temporal correction
+
history reconstruction
```

因此不同内部变量可能受到不同方向的补偿或偏移。

---

# 10. damage 再次确认不是 Local / Global history mismatch 的主要问题

Damage：

```text
Local / FOM
=
2.973214467e-02

Global / FOM
=
2.973622884e-02
```

Local / Global 自身差异：

```text
1.074355119e-05
```

几乎为零。

Final maximum damage：

```text
FOM
=
2.535861405195e-03

Local
=
2.434563263013e-03

Global
=
2.434540391614e-03
```

所以：

> 当前 damage 的约 3% full-field FOM discrepancy 并不是 Global history reconstruction 重新制造出来的。

这一结果与 FOM-2A 的 damage-chain diagnosis 一致：

```text
small stress-path discrepancy
↓
energy release rate Y discrepancy
↓
nonlinear damage-rate amplification
↓
accumulated damage discrepancy
```

---

# 11. Critical material point q = 12 的三方对比

FOM 临界材料点：

```text
(element, Gauss point, fiber)
=
(0, 0, 12)

canonical q
=
12
```

Critical-point relative L2：

| 场变量 | Local / FOM | Global / FOM | Local / Global |
|---|---:|---:|---:|
| total strain | 6.858214273e-03 | 5.171760092e-03 | 1.124790828e-02 |
| elastic strain | 4.881790166e-03 | 4.881231135e-03 | 1.644366368e-06 |
| stress | 4.871969081e-03 | 4.871414552e-03 | 1.649917335e-06 |
| plastic strain | 4.708037471e-02 | 6.839993121e-02 | 8.237690062e-02 |
| alpha | 5.015794136e-02 | 6.920085218e-02 | 8.292384535e-02 |
| r_bar | 4.620585138e-02 | 3.108198254e-02 | 4.954126123e-02 |
| damage | 4.434331194e-02 | 4.435226371e-02 | 9.460894791e-06 |

即：

```text
critical eps_p:
Local   ≈ 4.708 %
Global  ≈ 6.840 %

critical alpha:
Local   ≈ 5.016 %
Global  ≈ 6.920 %

critical r_bar:
Local   ≈ 4.621 %
Global  ≈ 3.108 %

critical damage:
Local   ≈ 4.434 %
Global  ≈ 4.435 %
```

因此 critical point 与 full field 给出了相同的总体趋势：

```text
eps_p / alpha:
Local 更好

r_bar:
Global 更好

damage:
Local ≈ Global
```

---

# 12. 为什么 Local 与 Global 的内部变量可以差很多，但 stress 几乎完全一样

这是本阶段最重要的机械解释之一。

Full-field：

```text
Local / Global plastic strain
≈ 8.25 %

Local / Global alpha
≈ 8.31 %

Local / Global r_bar
≈ 5.18 %
```

但是：

```text
Local / Global elastic strain
≈ 4.70e-06

Local / Global stress
≈ 4.70e-06
```

即机械响应几乎相同。

因为：

```text
elastic strain
=
total strain - plastic strain
```

Global stage 在改变 integrated plastic history 的同时，也会调整 total strain。

因此出现：

```text
plastic strain history 改变
+
total strain 同方向调整
↓
两者在 elastic strain 中高度抵消
↓
elastic strain 几乎不变
↓
stress 几乎不变
```

这进一步确认：

> 当前 Local / Global 的主要不闭合发生在 accumulated internal histories，而不是最终机械 stress state。

---

# 13. Fresh-Local 诊断：排除 result.local_state 滞后一轮的影响

进一步检查 solver transactional structure 后发现：

```text
result.local_state
```

是最后一个 accepted Global commit 之前那一轮 local projection 的 state。

严格来说，它不是：

```text
以最终 result.state
再重新做一次 local projection
```

因此为了避免“Local 好于 Global 只是因为滞后一轮”的疑问，本阶段又做了一次：

```text
final result.state
↓
solve_tower_local_stage(...)
↓
fresh_local_state
```

这一步：

- 不修改 `latin/`；
- 仅把最终 committed Global state 再投影一次到 local constitutive manifold；
- Local projection 内部继续采用 RK4。

---

# 14. 当前 tower Local stage 的实际数学行为

当前 `latin/tower_local_stage.py` 中：

Local projection 保持 Global thermodynamic-force histories 固定：

```text
stress_hat
=
stress_global

beta_hat
=
beta_global

R_bar_hat
=
R_bar_global

Y_hat
=
Y_global
```

然后每个 canonical material point `q` 独立沿时间顺序积分：

```text
[eps_p, alpha, r_bar, D]
```

积分器：

```text
RK4
```

因此 Fresh Local 是一个很干净的诊断：

```text
完全相同的 final Global force histories
↓
重新用 Local RK4 积分内部变量
```

特别地：

```text
Fresh Local stress
=
Final Global stress
```

这是算法定义直接决定的，而不是偶然的数值吻合。

---

# 15. Fresh-Local–FOM–Global 全场结果

结果：

| 场变量 | Returned Local / FOM | Fresh Local / FOM | Global / FOM | Returned / Fresh | Fresh / Global |
|---|---:|---:|---:|---:|---:|
| total strain | 3.558346576e-03 | 3.558340951e-03 | 3.575144269e-03 | 5.992187195e-06 | 5.115416136e-03 |
| elastic strain | 2.873920107e-03 | 2.873942972e-03 | 2.873943696e-03 | 4.708462816e-06 | 1.486090301e-08 |
| stress | 2.871167752e-03 | 2.871192934e-03 | 2.871192934e-03 | 4.703832963e-06 | 0 |
| plastic strain | 3.934805197e-02 | 3.934759997e-02 | 7.297194761e-02 | 3.224176402e-05 | 8.253100251e-02 |
| alpha | 4.164895720e-02 | 4.164850864e-02 | 7.351065984e-02 | 3.151295746e-05 | 8.305738617e-02 |
| r_bar | 3.957700223e-02 | 3.957682629e-02 | 3.463837253e-02 | 1.953640960e-05 | 5.175120176e-02 |
| damage | 2.973214467e-02 | 2.972277965e-02 | 2.973622884e-02 | 2.129594729e-05 | 2.627997080e-05 |

---

# 16. Fresh Local 诊断的关键结论

Returned Local 与 Fresh Local 的差异只有：

```text
plastic strain
≈ 3.22e-05

alpha
≈ 3.15e-05

r_bar
≈ 1.95e-05

damage
≈ 2.13e-05
```

因此：

> `result.local_state` 比最终 Global state 滞后一轮，对此前 Local-vs-FOM 判断没有实质影响。

Fresh Local 仍然得到：

```text
plastic strain / FOM
≈ 3.935 %

alpha / FOM
≈ 4.165 %
```

而 Global：

```text
plastic strain / FOM
≈ 7.297 %

alpha / FOM
≈ 7.351 %
```

所以：

> Local 对 `eps_p / alpha` 明显更接近 FOM 这一结论是稳健的，不是 transaction lag artefact。

---

# 17. Fresh Local critical-point 结果

Critical q = 12：

| 场变量 | Returned Local / FOM | Fresh Local / FOM | Global / FOM | Fresh / Global |
|---|---:|---:|---:|---:|
| total strain | 6.858214273e-03 | 6.858183842e-03 | 5.171760092e-03 | 1.124760327e-02 |
| elastic strain | 4.881790166e-03 | 4.881229146e-03 | 4.881231135e-03 | 2.583219234e-08 |
| stress | 4.871969081e-03 | 4.871414552e-03 | 4.871414552e-03 | 0 |
| plastic strain | 4.708037471e-02 | 4.708030190e-02 | 6.839993121e-02 | 8.237406704e-02 |
| alpha | 5.015794136e-02 | 5.016095662e-02 | 6.920085218e-02 | 8.292295871e-02 |
| r_bar | 4.620585138e-02 | 4.620916574e-02 | 3.108198254e-02 | 4.954403760e-02 |
| damage | 4.434331194e-02 | 4.433305094e-02 | 4.435226371e-02 | 2.031213695e-05 |

Critical-point 判断保持不变：

```text
eps_p:
Fresh Local ≈ 4.708 %
Global      ≈ 6.840 %

alpha:
Fresh Local ≈ 5.016 %
Global      ≈ 6.920 %

r_bar:
Fresh Local ≈ 4.621 %
Global      ≈ 3.108 %

damage:
Fresh Local ≈ 4.433 %
Global      ≈ 4.435 %
```

---

# 18. Fresh Local 与 Global stress 完全相同的意义

Fresh Local / Global：

```text
stress relative difference
=
0
```

这是当前 local projection 保持 Global stress history 固定的直接结果。

同时：

```text
Fresh Local / Global elastic strain
=
1.486090301e-08
```

几乎机器精度一致。

但：

```text
Fresh Local / Global plastic strain
=
8.253100251e-02
≈ 8.25 %
```

这意味着：

> 当前最终 LATIN 收敛点允许 Local 与 Global 用明显不同的 accumulated internal histories 表示几乎完全相同的机械 stress / elastic-strain state。

这是当前 Local RK4 / Global BE-like history inconsistency 最直接的数值表现。

---

# 19. Fresh Local total strain 为什么不能替代 Global total strain

Fresh Local 中：

```text
total strain
=
elastic strain + fresh-local plastic strain
```

但是 Fresh Local 只是 local constitutive projection。

它不重新求解：

```text
global structural compatibility
```

所以：

```text
Fresh Local total strain
```

不应被直接当作正式 structural solution。

因此必须区分：

```text
Fresh Local:
用于诊断 constitutive/history integration

Global:
用于表示 globally admissible structural response
```

这也是为什么即使 Local 某些内部变量更接近 FOM，也不能简单把整个 Local state 作为最终结构解。

---

# 20. 本阶段对“谁与 FOM 比较”的重新认识

此前正式 comparison 使用：

```text
Global result.state
```

这是严格的 LATIN structural solution。

本阶段经过三方诊断后，可以把结果分为两个层次。

## 20.1 Structural solution 层

下列量原则上应保持 Global 作为正式结构结果：

```text
tower displacement
total strain
elastic strain
stress
```

理由：

```text
Global
满足 equilibrium + compatibility

Local
不保证 global compatibility
```

## 20.2 Internal-variable accuracy 层

当前 one-cycle 40-step benchmark 上：

```text
eps_p:
Local 更接近 FOM

alpha:
Local 更接近 FOM

r_bar:
Global 更接近 FOM

damage:
Local ≈ Global
```

因此 Local state 对内部变量具有重要的 diagnostic / potentially preferable history value。

---

# 21. 用户提出的当前精度接受标准

经过上述结果讨论后，当前研究判断发生了一个重要转折。

用户提出：

> 如果暂时不纠结 Local 和 Global 谁是理论上唯一应该使用的 history，而是针对每个物理量采用当前更接近 FOM 的结果，那么 40 increments/cycle 下主要误差已经在约 5% 以内，这一精度目前可以接受。

按照当前 full-field best-available comparison：

| 场变量 | 当前较优结果 | relative L2 |
|---|---|---:|
| total strain | Local / Global 几乎相同 | ≈ 0.356% |
| elastic strain | Local / Global 几乎相同 | ≈ 0.287% |
| stress | Local / Global 几乎相同 | ≈ 0.287% |
| plastic strain | Local | ≈ 3.935% |
| alpha | Local | ≈ 4.165% |
| r_bar | Global | ≈ 3.464% |
| damage | Local / Global 几乎相同 | ≈ 2.972% |

所以：

> 在 full-field relative L2 意义下，当前 one-cycle、40 increments/cycle 的主要场变量误差都已经低于约 5%。

Critical q = 12：

```text
eps_p Local
≈ 4.708 %

alpha Local
≈ 5.016 %

r_bar Global
≈ 3.108 %

damage
≈ 4.43 %

stress
≈ 0.487 %
```

因此更严谨的表述应为：

> 40 increments/cycle 下，全场主要变量误差均约在 5% 以内，临界材料点内部变量误差也处于约 5% 量级。

不应表述为：

```text
所有任意位置、任意变量严格 < 5 %
```

因为 critical alpha 当前约为：

```text
5.016 %
```

---

# 22. 关于“哪个结果好就用哪个”的方法学边界

用户当前倾向采用：

```text
Local / Global
哪个更接近 FOM
就采用哪个作为该变量的精度代表
```

这一策略对当前研究阶段非常有价值，因为它可以回答：

> 当前 LATIN 框架能够提供的 best-available state accuracy 到底是多少？

但是需要明确一个论文方法学风险：

如果在未来每个循环都先看 FOM，然后再动态选择：

```text
第 10 循环选 Local
第 50 循环改选 Global
第 100 循环再次切换
```

那么这种选择属于：

```text
post-hoc selection
```

即：

> 只有知道 FOM 答案后才能决定使用哪个结果。

这会削弱方法作为独立预测工具的意义。

因此后续有两种严谨表达方式。

### 方案 A：固定选择规则

在进入 multicycle validation 前冻结：

```text
structural displacement / strain / stress
→ Global

eps_p
→ Local

alpha
→ Local

r_bar
→ Global

damage
→ Global 或 Local
```

其中 damage 两者几乎一致，可保持 Global 作为 final committed state。

之后：

> 无论 2 / 5 / 10 / 100 cycles 哪个更接近 FOM，都不再改选择规则。

### 方案 B：同时报告两者，并定义 best-available envelope

正式结构解仍报告：

```text
Global
```

同时对内部变量报告：

```text
Local error
Global error
min(Local error, Global error)
```

这里：

```text
min(Local error, Global error)
```

应被称为：

```text
best-available LATIN internal-state accuracy envelope
```

而不能直接称为唯一的 solver prediction。

当前用户的研究倾向更接近方案 B 的“哪个结果好就用哪个”评价思路。

后续进入论文定稿前，应再明确采用 A 还是 B。

---

# 23. 不同时间步长下 FOM 自身的变化

本阶段进一步提出：

> LATIN 与 FOM 的误差随 dt 变化之前，FOM 自己在 40 / 80 / 160 increments/cycle 下到底变化多少？

当前已经缓存：

```text
40 increments:
outputs/tower_1cycle_reversed_fom_reference_v1.npz

80 increments:
outputs/tower_1cycle_reversed_fom_reference_80inc_v1.npz

160 increments:
outputs/tower_1cycle_reversed_fom_reference_160inc_v1.npz
```

已有文档明确记录的 FOM final anchors：

| increments/cycle | dt | FOM final max \|eps_p\| | FOM final max D |
|---|---:|---:|---:|
| 40 | 0.25 s | 6.921951176224e-05 | 2.535861405195e-03 |
| 80 | 0.125 s | 7.028958289300e-05 | 2.567517301383e-03 |
| 160 | 0.0625 s | 7.055638197426e-05 | 2.575449350984e-03 |

---

# 24. FOM final max |eps_p| 的时间步敏感性

绝对差：

```text
40 → 80:
7.028958289300e-05
-
6.921951176224e-05
=
1.07007113076e-06

80 → 160:
7.055638197426e-05
-
7.028958289300e-05
=
2.6679908126e-07
```

以后一个更细网格为 reference 计算相对变化：

```text
40 → 80
≈ 1.522 %

80 → 160
≈ 0.378 %

40 → 160
≈ 1.895 %
```

所以：

> 40-step FOM 的 final max |eps_p| 相对 160-step FOM 仍有约 1.9% 的离散差异。

---

# 25. FOM final max D 的时间步敏感性

绝对差：

```text
40 → 80
≈ 3.1655896188e-05

80 → 160
≈ 7.932049601e-06
```

以后一个更细网格为 reference：

```text
40 → 80
≈ 1.233 %

80 → 160
≈ 0.308 %

40 → 160
≈ 1.537 %
```

因此：

> 40-step FOM 的 final max D 相对 160-step FOM 仍有约 1.5% 的离散差异。

---

# 26. FOM 自身 final-anchor refinement 的观测阶

对于 `max |eps_p|`：

```text
(40→80 difference)
/
(80→160 difference)
≈
4.0108
```

对应：

```text
p
=
log2(4.0108)
≈
2.004
```

对于 `max D`：

```text
difference ratio
≈
3.9909
```

对应：

```text
p
≈
1.997
```

因此这两个 final scalar anchors 表现出：

```text
approximately second-order refinement behaviour
```

但是必须限定：

> 这只是 `final max |eps_p|` 和 `final max D` 两个 scalar response indicators 的观测特征。

当前不能据此直接宣布：

```text
整个 FOM algorithm
=
严格二阶时间积分方法
```

因为 FOM 还包含：

- structural Newton solve；
- load-step interpolation；
- strain-driven material update；
- RK4 constitutive sub-integration；
- nonlinear state coupling。

---

# 27. FOM 时间步敏感性与 FOM–LATIN discrepancy 的相对量级

当前 40-step：

```text
FOM-only final max |eps_p|
40 vs 160
≈ 1.9 %
```

而：

```text
Global LATIN / FOM plastic-strain full-field error
≈ 7.3 %
```

Fresh Local：

```text
Fresh Local / FOM plastic-strain full-field error
≈ 3.93 %
```

这说明：

> FOM 本身的时间离散误差确实不能忽略，但仅靠 FOM reference 的 40-step 时间误差不足以解释原先约 7% 的 Global plastic-history discrepancy。

同样：

```text
FOM-only final damage
40 vs 160
≈ 1.5 %
```

而当前 full-field LATIN / FOM damage discrepancy：

```text
≈ 3 %
```

因此 FOM 自身的 finite-dt error 是 comparison error budget 的一部分，但不是全部。

---

# 28. 尚未完成的 FOM-only full-history refinement

目前对 FOM 自身只完成了已有结果的：

```text
final max |eps_p|
final max D
```

比较。

还没有正式计算：

```text
FOM40 vs FOM80
FOM80 vs FOM160
```

的全场 time-history relative L2：

```text
total strain
stress
eps_p
alpha
r_bar
damage
```

由于三个时间网格严格嵌套：

```text
40-grid nodes
⊂
80-grid nodes
⊂
160-grid nodes
```

理论上可以直接：

```text
80-step history 每隔 2 点取样
160-step history 每隔 4 点取样
```

在相同物理时刻比较，而不需要 interpolation。

这项诊断目前尚未执行。

考虑到当前研究重心已经准备转向 multicycle accuracy / efficiency，它不再作为进入下一阶段的强制 gate，但可以保留为后续补充验证。

---

# 29. 当前对 40 increments/cycle 的总体评价

综合现有数据：

```text
FOM itself:
40-step reference 仍有约 1–2% final-anchor temporal effect

LATIN best-available one-cycle full-field:
major fields ≲ 5%

critical material-point internal states:
≈ 5% level
```

因此当前研究判断是：

> 40 increments/cycle 并不是严格时间离散无关的网格，但在当前 tower benchmark 上已经可以作为 multicycle exploratory validation 的实用基准。

其理由不是：

```text
40 steps 已经完全收敛
```

而是：

```text
1. 当前 one-cycle main-field accuracy 已达到约 5% 量级；

2. 继续把 one-cycle 误差从约 4% 压到约 2%
   对当前研究主目标的边际价值开始降低；

3. LATIN-PGD 真正的研究价值不在单循环，
   而在大量循环下是否保持误差可控并显著降低计算成本。
```

---

# 30. 本阶段研究重心的正式转变

此前研究问题主要是：

```text
LATIN-PGD one-cycle
到底够不够准？
```

经过 I-0 至 I-6、FOM-1、FOM-2A 至 FOM-2D 以及本阶段 Local / Global supplemental diagnostics 后，当前问题已转为：

```text
随着循环次数增加：

1. 误差会不会累积？
2. 内部变量误差会不会放大？
3. damage error 会不会逐循环积累？
4. PGD rank 会如何增长？
5. LATIN outer iterations 会如何增长？
6. 计算时间能否显著低于 FOM？
7. 循环数越大，LATIN-PGD 的 speedup 是否越明显？
```

这比继续纠缠：

```text
one-cycle:
3.9 %
vs
2 %
```

更符合当前“LATIN-PGD 用于海上风机塔筒高周疲劳”的研究目标。

---

# 31. 下一阶段建议命名

下一阶段建议正式命名为：

```text
FOM-3
Multicycle accuracy and computational efficiency validation
```

即：

```text
多循环精度
+
计算效率
```

双主线验证。

---

# 32. FOM-3 的核心科学问题 A：误差是否随循环数累积

当前 one-cycle：

```text
best-available internal-state error
≈ 3–5 %
```

但这并不能保证：

```text
10 cycles
100 cycles
```

仍然是相同量级。

必须检查：

```text
cycle 1
cycle 2
cycle 5
cycle 10
...
cycle 100
```

时：

```text
stress
eps_p
alpha
r_bar
damage
residual displacement
```

的 discrepancy 是：

```text
保持稳定
```

还是：

```text
逐循环积累
```

甚至：

```text
发生误差失稳
```

---

# 33. FOM-3 的核心科学问题 B：LATIN-PGD 是否真正降低计算成本

LATIN-PGD 的目标不是单纯替代 FOM 得到近似相同结果，而是：

> 利用 x–t separated representation 和低秩 PGD，使大量循环问题的计算成本增长速度明显低于 conventional step-by-step nonlinear FOM。

因此后续必须记录：

```text
FOM wall-clock time

LATIN wall-clock time

speedup
=
FOM time / LATIN time
```

同时还应记录：

```text
FOM:
Newton iterations
load steps

LATIN:
outer iterations
trial evaluations
PGD rank
modes added
fixed-point iterations
final xi
```

才能判断：

> LATIN-PGD 的速度优势究竟来自哪里，以及这种优势是否随循环数增加而扩大。

---

# 34. 已有冻结 100-cycle FOM reference

当前已有正式冻结的 100-cycle FOM：

```text
outputs/tower_100cycle_fom_reference_v1.npz
```

该文件：

- 已完成完整 100-cycle nonlinear FOM；
- 位于 `outputs/`，被 Git ignore；
- 不应无必要重复计算；
- 不应删除或覆盖。

其 SHA256：

```text
b3230d577341e03e598463db4c89372e0ef8e21b51144ef59b2f3175b0b1f4e8
```

主要 shape：

```text
cycle_numbers
=
(100,)

phase_times
=
(41,)

nodal
=
(100, 41, 33)

strains / stresses
=
(100, 41, 10, 2, 16)

states
=
(100, 41, 10, 2, 16, 4)
```

Final anchors：

```text
final max |eps_p|
=
1.363490241501e-03

final max D
=
4.409564790481e-02
```

History max：

```text
max |eps_p|
=
1.366925100668e-03
```

---

# 35. 100-cycle frozen FOM 与当前 one-cycle benchmark 的荷载区别

必须特别注意：

当前 one-cycle accuracy benchmark 是：

```text
fully reversed

+1.0 MN
↔
-1.0 MN

R_F
=
-1
```

而冻结 100-cycle FOM 是：

```text
asymmetric cyclic

Fmax
=
+1.0 MN

Fmin
=
-0.5 MN

R_F
=
-0.5

mean force
=
+0.25 MN

force amplitude
=
0.75 MN
```

并且包含：

```text
explicit preload:
0
→
mean force
```

因此：

> 不能把当前 one-cycle fully reversed LATIN 结果直接与 frozen 100-cycle asymmetric FOM 的 cycle 1 比较。

进入 multicycle validation 后，必须构造与 frozen FOM 完全相同的：

```text
preload
+
asymmetric cyclic loading
```

LATIN problem。

---

# 36. FOM-3A 的建议第一步

在真正运行 multicycle LATIN 之前，先审计 frozen 100-cycle FOM 的 indexing。

需要明确：

```text
cycle_numbers
phase_times
preload
cycle-1 start
cycle-1 end
```

之间的物理对应。

重点确认：

```text
array index 0
```

到底表示：

```text
preload 后的 cycle-1 起点
```

还是：

```text
完整 cycle-1 history 的第一个 phase node
```

以及各 cycle 是否包含重复的：

```text
cycle end
=
next cycle start
```

节点。

---

# 37. 建议的 early multicycle validation sequence

在 100-cycle LATIN 一次性开跑之前，建议先比较：

```text
cycle 1
cycle 2
cycle 5
cycle 10
```

原因：

```text
1 cycle
→
检查荷载 / preload / state mapping 正确

2 cycles
→
首次检查跨循环 history propagation

5 cycles
→
检查误差是否出现累积趋势

10 cycles
→
检查短期多循环下的 rank / iteration / timing scaling
```

只有这些结果稳定后，再进入：

```text
20 / 50 / 100 cycles
```

会更加安全。

---

# 38. 多循环阶段建议固定记录的 accuracy metrics

建议每个指定 cycle 至少记录：

```text
tower-top displacement

total strain

elastic strain

stress

plastic strain

alpha

r_bar

damage
```

并同时给：

```text
full-field relative L2

critical-point relative L2

cycle-end scalar anchors

maximum over phase history
```

特别关注：

```text
max |eps_p|
max D
residual top displacement
critical stress
```

---

# 39. Local / Global 在多循环阶段的建议记录方式

为了不丢失本阶段发现的信息，multicycle 阶段建议同时保留：

```text
Global state

Local state
```

至少对内部变量：

```text
eps_p
alpha
r_bar
damage
```

分别记录：

```text
Local / FOM
Global / FOM
Local / Global
```

这样可以回答：

> Local / Global history gap 是否随着循环数增加而累积？

这是下一阶段非常有价值的附加结果。

---

# 40. “best result” 多循环评价策略的当前建议

用户当前明确提出：

> 对比 FOM 时，不必预先限定 Local 或 Global，哪个结果更好就用哪个，并以此判断当前方法能够达到的有效精度。

这一策略可继续用于：

```text
exploratory validation
```

但建议在结果表中始终同时保留：

```text
Local error
Global error
best error
```

即：

```text
best error
=
min(Local error, Global error)
```

这样不会掩盖：

```text
Local / Global disagreement
```

同时还能直观看出：

```text
当前 LATIN framework
在该变量上能够提供的 best-available accuracy
```

对于：

```text
displacement
total strain
stress
```

由于 Local 不具备完整 global admissibility，仍建议以 Global 为主要 structural result。

---

# 41. 多循环阶段建议记录的 computational-cost metrics

每个循环数：

```text
Ncycle
=
1
2
5
10
...
```

建议记录：

### FOM

```text
wall-clock time

number of load steps

total Newton iterations

maximum Newton iterations / step

average Newton iterations / step

peak memory
```

### LATIN-PGD

```text
wall-clock time

outer LATIN iterations

trial evaluations

accepted modes

final PGD rank

fixed-point iterations

final LATIN indicator

peak memory
```

最终计算：

```text
speedup
=
T_FOM / T_LATIN
```

以及必要时：

```text
memory reduction
```

---

# 42. 真正决定 LATIN-PGD 工程价值的结果组合

后续最有价值的结果不是单独：

```text
error = 3 %
```

也不是单独：

```text
speedup = 10x
```

而是二者组合。

例如若最终得到：

```text
100-cycle error
≈ 5 %

speedup
≈ 10–50x
```

则可以形成很强的结论：

> LATIN-PGD 能够在保持可接受疲劳历史精度的同时显著降低海上风机塔筒多循环非线性分析成本。

反之，如果：

```text
100-cycle error
→
20–30 %
```

即使 speedup 很高，也需要进一步研究：

```text
history correction
temporal refinement
windowing
restart
adaptive enrichment
```

等策略。

---

# 43. 尚未解决的理论问题必须继续保留

尽管当前准备进入 multicycle validation，但以下问题并没有被解决。

## 43.1 Local RK4 / Global BE 的理论一致性

当前已经数值证明：

```text
finite dt 下两者 history 不一致
```

并且：

```text
dt 减小时 gap 系统下降
```

但仍未证明：

```text
LATIN 理论是否允许 Local 与 Global
使用不同 temporal integration rule
```

以及：

```text
是否必须在统一离散时间空间中满足 search-direction consistency
```

---

# 44. 原论文 DG0 与当前 BE 的严格关系仍未知

已知原论文明确说明：

```text
DG0
```

被用于 Eq. (59) 的 time integration。

但当前 tower implementation 使用：

```text
Backward-Euler-type temporal discretisation
```

尚未完成：

```text
paper DG0 weak form
→
jump / trace treatment
→
discrete algebraic equations
→
current BE implementation
```

的严格推导。

因此当前不能声称：

```text
current BE
=
paper DG0 exactly
```

只能说：

> 当前 BE-like global temporal treatment 已在 one-cycle benchmark 上表现出稳定的一阶 FOM convergence，但与原论文 DG0 的严格等价性尚未建立。

---

# 45. 本阶段策略变化不等于理论问题已经关闭

本阶段决定把研究重心转向 multicycle，不意味着：

```text
Local / Global temporal inconsistency
=
solved
```

准确状态仍然是：

```text
误差来源诊断
=
已取得强数值证据

误差随 dt 收敛
=
已验证

Local / Global 统一时间积分
=
未解决

paper DG0 vs current BE
=
未解决

是否需要修改核心算法
=
暂不决定
```

当前选择是：

> 在 one-cycle accuracy 已达到约 5% 量级的情况下，先验证 multicycle error accumulation 与 computational benefit，再判断是否值得投入成本统一时间离散。

---

# 46. 为什么现在不立即把 Global BE 换成 RK4

当前没有充分依据直接执行：

```text
Global BE
→
Global RK4
```

原因：

```text
1. Global temporal problem 不是普通单材料点 ODE；

2. 它来自 PGD / LATIN global equations；

3. 原论文使用 DG0，而不是直接说明 RK4；

4. r_bar 当前 Global 反而比 Local 更接近 FOM；

5. 当前 mixed scheme 在 dt refinement 下已经表现出清晰收敛；

6. 未确认理论目标前修改 core algorithm
   可能破坏现有已验证的平衡、事务、PGD enrichment 和收敛链。
```

所以当前策略仍然是：

```text
diagnose first
validate multicycle value
then decide whether algorithmic unification is worth doing
```

---

# 47. 本阶段最终形成的误差因果图

当前 one-cycle benchmark 可以概括为：

```text
Local / Global temporal treatment differs
                ↓
finite-dt integrated history mismatch
                ↓
eps_p / alpha / r_bar histories differ
                ↓
global equilibrium / compatibility
adjust total strain
                ↓
large part of eps_p error cancels
inside elastic strain
                ↓
elastic strain discrepancy remains small
                ↓
stress discrepancy remains small
                ↓
Y discrepancy
                ↓
nonlinear Ddot amplification
                ↓
damage discrepancy
```

同时：

```text
Fresh Local stress
=
Global stress
```

而：

```text
Fresh Local eps_p
≠
Global eps_p
```

说明：

> 当前主要分歧集中在“历史变量表示”，而不是“机械应力状态”。

---

# 48. 当前 one-cycle accuracy 的推荐总结表述

后续文档或论文中建议使用：

> For the one-cycle tower benchmark with 40 increments per cycle, the globally admissible LATIN-PGD solution reproduces the displacement, elastic strain and stress fields with sub-percent discrepancies relative to the full-order step-by-step solution. The internal history variables exhibit larger finite-time-step discrepancies. When the local constitutive history is also examined, the best-available full-field errors of plastic strain, kinematic hardening, isotropic hardening and damage remain at approximately the 3–4% level, while the critical material-point discrepancies are approximately at the 5% level.

中文：

> 对于每循环 40 个时间增量的一循环塔筒算例，LATIN-PGD 的整体位移、弹性应变和应力与逐时间步全阶解的误差均低于 1%。内部历史变量受有限时间步和 Local / Global 时间离散差异影响更明显。综合考察 Local 与 Global history 后，塑性应变、运动硬化、各向同性硬化和损伤的全场最优误差约为 3% 至 4%，临界材料点的内部变量误差约为 5% 量级。

这比：

```text
所有误差都小于 5 %
```

更加准确。

---

# 49. 当前研究阶段的正式判断

截至本阶段，可以把项目状态更新为：

```text
Tower LATIN-PGD internal algorithm closure
=
PASS

One-cycle FOM contract / mapping
=
PASS

One-cycle global mechanical accuracy
=
STRONG PASS

One-cycle internal-state accuracy
=
approximately 3–5% best-available level

Temporal refinement behaviour
=
PASS

Local / Global history discrepancy diagnosis
=
STRONG DIAGNOSIS

Fresh-Local confirmation
=
PASS

Local / Global temporal integrator unification
=
OPEN

paper DG0 equivalence
=
OPEN

Multicycle accuracy
=
NOT YET VALIDATED

Computational speedup
=
NOT YET QUANTIFIED
```

---

# 50. 下一阶段进入条件

基于用户当前对约 5% one-cycle internal-state error 的接受度，已经具备进入：

```text
FOM-3
Multicycle accuracy and computational efficiency
```

的条件。

但进入前首先执行：

```text
FOM-3A
Frozen 100-cycle FOM indexing and preload audit
```

不直接修改 LATIN core。

---

# 51. FOM-3A 的具体目标

下一步只回答：

> 现有 `tower_100cycle_fom_reference_v1.npz` 中每一个 cycle / phase index 对应什么物理时刻？

需要明确：

```text
preload:
0 → +0.25 MN

cycle loading:
+1.0 MN ↔ -0.5 MN

cycle_numbers
phase_times
nodal
strains
stresses
states
```

之间的严格对应。

完成后才能安全提取：

```text
cycle 1
cycle 2
cycle 5
cycle 10
```

作为 multicycle LATIN comparison targets。

---

# 52. 本阶段没有完成或不应误写为完成的内容

以下内容必须明确标记为未完成：

```text
1. 没有修改 Local RK4。

2. 没有修改 Global BE。

3. 没有统一 Local / Global 时间积分。

4. 没有证明 paper DG0 ≡ current BE。

5. 没有完成 FOM40 / FOM80 / FOM160 的 full-history FOM-only relative-L2 comparison。

6. 没有进行 asymmetric multicycle LATIN-PGD accuracy comparison。

7. 没有得到 10 / 50 / 100-cycle LATIN error。

8. 没有得到 FOM vs LATIN speedup。

9. 没有证明 100-cycle accuracy 可接受。

10. 没有证明 high-cycle rank 会保持低增长。
```

---

# 53. 本阶段产生的临时诊断脚本

本阶段对话中使用了两个临时根目录诊断脚本：

```text
tower_local_fom_global_diagnostic.py

tower_fresh_local_fom_global_diagnostic.py
```

用途：

```text
第一个：
FOM vs returned Local vs final Global

第二个：
FOM vs returned Local vs fresh Local vs final Global
```

这两个脚本：

- 是诊断工具；
- 没有修改 `latin/`；
- 当前不应自动视为 production code；
- 是否保留到 `examples/` 或 `diagnostics/` 需在阶段整理时再决定；
- 如果仅用于一次性诊断，可以在 Git checkpoint 前删除，不必提交。

---

# 54. 本阶段最重要的三个新认识

## 54.1 Local history 不能忽略

此前：

```text
Global
=
唯一正式 FOM comparison object
```

现在更新为：

```text
Global
=
正式 globally admissible structural state

Local
=
非常有价值的 constitutive-history diagnostic state
```

尤其：

```text
eps_p
alpha
```

Local 明显比 Global 更接近 FOM。

---

## 54.2 当前误差已经从“能不能用”进入“多循环是否稳定”的阶段

One-cycle：

```text
macro response
< 1 %

best-available internal states
≈ 3–5 %
```

因此下一步最有研究价值的问题不再是：

```text
把 4 % 改成 2 %
```

而是：

```text
4 % 会不会在 100 cycles 变成 20 %？
```

---

## 54.3 LATIN-PGD 是否有价值最终必须由 accuracy × speedup 共同决定

最终研究判断必须同时回答：

```text
误差
+
速度
```

只有：

```text
long-cycle accuracy acceptable
+
computational cost substantially reduced
```

才能证明：

> LATIN-PGD 对海上风机塔筒高周循环疲劳分析具有实际数值价值。

---

# 55. 下一步工作

下一步建议严格按照：

```text
FOM-3A
Frozen asymmetric 100-cycle FOM indexing / preload audit
```

开始。

暂时：

```text
不修改 latin/
不重新跑 frozen 100-cycle FOM
不直接启动 100-cycle LATIN
```

先确认：

```text
cycle 1 / 2 / 5 / 10
```

在冻结 FOM 文件中的 exact indexing。

通过后，再构造 matched asymmetric LATIN-PGD benchmark，并同步开始记录：

```text
accuracy
+
wall-clock time
+
PGD rank
+
outer iterations
```

---

# 56. 当前 checkpoint 前建议

在进入 FOM-3A 之前，建议将本文档作为独立阶段总结提交。

建议文件名：

```text
docs/2026-08-22-tower-local-global-fom-supplement-and-multicycle-transition.md
```

建议 commit message：

```text
docs: summarize local-global FOM diagnostics and multicycle transition
```

建议 checkpoint 内容只包含：

```text
本阶段总结 Markdown
```

临时诊断脚本是否纳入 Git，应单独决定，不与总结文档自动混在同一个提交中。

---

# 57. 一句话阶段结论

> 在一循环 40 increments/cycle 的塔筒算例中，LATIN-PGD 已能以低于 1% 的误差复现主要机械响应；Local 与 Global 的内部历史变量虽存在明显有限时间步不一致，但综合使用当前更接近 FOM 的 history 后，全场内部变量误差约为 3% 至 4%，临界材料点约为 5% 量级。鉴于这一精度目前可接受，研究重点正式从继续压缩单循环误差转向验证多循环误差是否累积，以及 LATIN-PGD 相对于逐时间步 FOM 是否能够随着循环数增加体现出显著的计算效率优势。Local RK4 / Global BE 的理论统一以及原论文 DG0 与当前 BE 的严格关系仍作为开放问题保留，不应误写为已经解决。
