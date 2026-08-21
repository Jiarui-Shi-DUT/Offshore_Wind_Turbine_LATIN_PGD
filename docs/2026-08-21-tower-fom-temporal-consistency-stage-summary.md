# Offshore Wind Turbine Tower LATIN-PGD
# FOM-2 时间离散一致性诊断阶段总结

日期：2026-08-21  
分支：`feature/offshore-wind-turbine-tower-fatigue`  
前一阶段 checkpoint：`c5fd1b9`  
前一阶段文档：`docs/2026-08-20-tower-fom-one-cycle-validation-and-temporal-discretization-diagnosis.md`

---

## 1. 本阶段目的

FOM-1 已经证明：

- 塔筒整体位移、总应变和应力与 FOM 很接近；
- plastic strain、alpha、r_bar 和 damage 等内部变量仍存在几个百分点差异；
- LATIN 内部的 plastic strain、alpha、r_bar history gap 与 Local RK4 / Global Backward Euler 的时间离散不一致高度相关；
- 40 / 80 / 160 时间步细化时，LATIN Local / Global history gap 呈近似一阶下降。

但当时仍有两个问题没有完全闭环：

1. FOM 与 LATIN-PGD 之间实际的 damage 误差究竟从哪里产生；
2. FOM 与 LATIN-PGD 的真实误差是否也会随着 dt 减小而系统收敛。

因此本阶段 FOM-2 不修改任何 `latin/` 核心代码，而是围绕以下四个诊断展开：

- FOM-2A：critical-point damage-chain diagnostic；
- FOM-2B：critical-point stress-error decomposition；
- FOM-2C：elastic-strain error cancellation diagnostic；
- FOM-2D：matched FOM-LATIN temporal refinement。

本阶段核心目标不是继续追求更小的单次误差，而是回答：

> 当前 FOM-LATIN discrepancy 是可随时间网格细化消失的离散误差，还是存在不随 dt 消失的更深层算法不一致。

---

# 2. Benchmark 定义

本阶段继续采用与 I-5 和 FOM-1 相同的 fully reversed one-cycle tower benchmark。

## 2.1 结构离散

- tower elements：10；
- Gauss points / element：2；
- circumferential fibers / section：16；
- radial fibers：1；
- material points：320；
- structural DOFs：33；
- free DOFs：30。

Canonical material-point ordering：

```text
q
<->
(element, Gauss point, local fiber)
```

排序规则：

```text
element-major
→ Gauss-major
→ fiber-major
```

即：

```text
q = (element * n_gauss + gauss) * n_fibers + fiber
```

## 2.2 Loading

Fully reversed sinusoidal tower-top horizontal force：

```text
F(t) = Fa sin(2 pi t / T)
```

其中：

```text
Fa = 1.0 MN
T  = 10 s
```

一个完整循环：

```text
0
→ +1.0 MN
→ 0
→ -1.0 MN
→ 0
```

FOM-2D 中仅改变：

```text
increments_per_cycle
=
40
80
160
```

对应：

```text
dt
=
0.25
0.125
0.0625 s
```

除时间网格外，其余结构、材料、荷载和空间离散完全一致。

---

# 3. FOM 缓存策略

为避免重复执行 full-order nonlinear tower analysis，本阶段建立了 one-cycle FOM 缓存。

40-increment reference：

```text
outputs/tower_1cycle_reversed_fom_reference_v1.npz
```

80-increment reference：

```text
outputs/tower_1cycle_reversed_fom_reference_80inc_v1.npz
```

160-increment reference：

```text
outputs/tower_1cycle_reversed_fom_reference_160inc_v1.npz
```

这些文件位于已被 `.gitignore` 忽略的：

```text
outputs/
```

因此：

- 可以长期复用；
- 不污染 Git working tree；
- 不应提交到仓库；
- 后续 one-cycle refinement diagnostics 不需要再次执行相同 FOM。

FOM-2D 同时输出：

```text
outputs/tower_1cycle_fom_latin_temporal_refinement_v1.csv
```

---

# 4. 当前材料损伤链

当前材料模型中，damage evolution 可简化理解为：

```text
stress
↓
energy release rate Y
↓
positive part <Y - Y0>+
↓
damage rate Ddot
↓
integrated damage D
```

当前参数：

```text
Y0
=
2.388059701493e-02

k_damage
=
2.778

n_damage
=
2
```

因此：

```text
Ddot
=
k_damage * <Y - Y0>^2
```

这一平方关系意味着：

> 接近损伤激活区间时，小的 stress / Y discrepancy 可能被明显放大到 damage-rate discrepancy。

---

# 5. FOM-2A：critical material point damage-chain diagnostic

## 5.1 目的

FOM-1 中：

```text
damage relative L2 error
≈ 3 %
```

但 LATIN Local / Global damage 本身几乎一致，因此 damage discrepancy 不能简单归因于：

```text
Local RK4
vs
Global BE
```

FOM-2A 的目的，是沿完整 damage chain 逐层追踪：

```text
stress
→ Y
→ Ddot
→ D
```

并通过变量替换试验判断：

```text
Y error
```

究竟主要来自：

```text
stress difference
```

还是：

```text
already accumulated damage difference
```

## 5.2 Critical material point

FOM 临界材料点：

```text
(element, gauss, fiber)
=
(0, 0, 12)
```

Canonical index：

```text
q = 12
```

该点在本阶段后续 FOM-2B / FOM-2C 中继续固定使用。

---

# 6. FOM-2A solver state

LATIN-PGD：

```text
spatial strategy
=
residual_ls

tolerance
=
1e-5
```

结果：

```text
converged
=
True

termination
=
converged

iterations
=
20

PGD rank
=
15

final xi
=
8.006194456505e-06
```

说明：

> 后续 damage-chain discrepancy 不是明显的 outer LATIN non-convergence artefact。

---

# 7. FOM-2A：完整 damage-chain 误差

Critical point 全时程 relative L2：

```text
stress
=
4.8714145525e-03
≈ 0.487 %

Y
=
1.1533915127e-02
≈ 1.153 %

instantaneous Ddot
=
5.0728389153e-02
≈ 5.073 %

step-average Ddot
=
5.2312671171e-02
≈ 5.231 %

damage D
=
4.4352263713e-02
≈ 4.435 %
```

误差链：

```text
stress
0.487 %
↓
Y
1.153 %
↓
Ddot
5.073 %
↓
D
4.435 %
```

因此：

> damage error 并不是从 D 本身突然出现，而是上游 stress discrepancy 经过非线性 damage law 逐层放大。

---

# 8. FOM-2A：peak discrepancy

Peak values：

```text
stress:
FOM
=
1.1596359294e+02 MPa

LATIN
=
1.1515952009e+02 MPa

peak relative difference
=
6.9338386699e-03
≈ 0.693 %
```

```text
Y:
FOM
=
5.0291411759e-02

LATIN
=
4.9588379845e-02

peak relative difference
=
1.3979164429e-02
≈ 1.398 %
```

```text
instantaneous Ddot:
FOM
=
1.9377414943e-03

LATIN
=
1.8359526923e-03

peak relative difference
=
5.2529608448e-02
≈ 5.253 %
```

```text
damage:
FOM
=
2.5358614052e-03

LATIN
=
2.4345403916e-03

peak relative difference
=
3.9955264658e-02
≈ 3.996 %
```

---

# 9. 为什么 Y 的小误差会产生更大的 Ddot 误差

Critical point damage-rate peak 附近：

```text
FOM:
Y
=
0.05029141

LATIN:
Y
=
0.04958838
```

损伤阈值：

```text
Y0
=
0.02388060
```

因此：

```text
FOM:
Y - Y0
≈
0.02641081

LATIN:
Y - Y0
≈
0.02570778
```

损伤率满足：

```text
Ddot
∝
(Y - Y0)^2
```

因此：

```text
(0.02570778 / 0.02641081)^2
≈
0.9475
```

即 LATIN peak Ddot 预计约比 FOM 小：

```text
5.25 %
```

这与实际计算：

```text
5.25296 %
```

几乎完全一致。

结论：

> 当前约 5 % 的 damage-rate discrepancy 可以由 Y discrepancy 在平方损伤律中的非线性放大定量解释。

---

# 10. FOM-2A：damage activation timing

FOM：

```text
active nodes
=
10

first active node
=
25

last active node
=
34
```

LATIN：

```text
active nodes
=
10

first active node
=
25

last active node
=
34
```

因此：

> FOM 和 LATIN 对 damage activation window 的判断完全一致。

没有出现 FOM 与 LATIN 在损伤起始或终止时刻上的错位。

当前差异主要是：

> 在同一损伤激活时间区间内，LATIN stress / Y 略低，因此 Ddot 略低，最终 D 略低。

---

# 11. FOM-2A：Y-error decomposition

完整 LATIN Y error：

```text
1.153391512742e-02
```

使用：

```text
LATIN stress
+
FOM damage
```

得到：

```text
1.147186982018e-02
```

使用：

```text
FOM stress
+
LATIN damage
```

得到：

```text
1.609321635847e-04
```

即：

```text
full Y error
≈
stress-only Y error
```

而：

```text
damage-only Y error
≈
0
```

结论：

> Y discrepancy 几乎全部由 stress path discrepancy 引起，而不是由已经累计的 damage discrepancy 反馈引起。

---

# 12. FOM-2A：Ddot-error decomposition

完整 instantaneous Ddot error：

```text
5.072838915288e-02
```

仅保留 stress discrepancy：

```text
5.028777393968e-02
```

仅保留 damage discrepancy：

```text
6.804244140494e-04
```

因此：

```text
full Ddot error
≈
stress-only propagated Ddot error
```

结论：

> damage-rate discrepancy 主要是 stress error 通过 Y 和 damage law 传播形成的，而不是由 D 自反馈累积形成。

---

# 13. FOM-2A：LATIN stored-field consistency

LATIN stored energy-release rate 与根据最终 stress / damage 重新计算的 Y：

```text
relative L2
=
5.131158758230e-12
```

LATIN stored damage_rate 与 damage law 直接重算：

```text
relative L2
=
1.840081182861e-05
```

二者都非常小。

因此：

- `energy_release_rate` stored field 一致；
- `damage_rate` stored field 一致；
- 没有证据表明 damage-law implementation 本身存在显著代数错误。

---

# 14. FOM-2A 阶段结论

FOM-2A 将原先：

```text
damage ≈ 3 % discrepancy
原因未知
```

更新为：

```text
small stress-path discrepancy
↓
Y discrepancy
↓
nonlinear damage-law amplification
↓
Ddot discrepancy
↓
accumulated D discrepancy
```

因此：

> damage 不再是一个独立的未知问题。

下一步需要继续向上游追踪：

> stress discrepancy 从哪里产生。

---

# 15. FOM-2B：critical-point stress-error decomposition

## 15.1 目的

Stress relation 可简化理解为：

```text
stress
=
effective modulus
*
elastic strain
```

而：

```text
elastic strain
=
total strain
-
plastic strain
```

因此 critical-point stress discrepancy 可能来自：

- total strain；
- plastic strain；
- damage 对 effective modulus 的影响。

FOM-2B 通过 one-variable replacement 分别隔离三者影响。

---

# 16. FOM-2B：stress constitutive closure

Actual LATIN global stress relative to FOM：

```text
4.871414552459e-03
```

利用 LATIN total strain、plastic strain、damage 重新构造 stress relative to FOM：

```text
4.871414212671e-03
```

Actual LATIN stress 与 constitutive reconstruction：

```text
2.843532786300e-09
```

因此：

> LATIN global stress 与自身最终 constitutive state 高度一致。

可以基本排除 final global stress field assembly bug。

---

# 17. FOM-2B：one-variable replacement

仅替换 total strain 为 LATIN：

```text
stress error
=
5.365338895032e-03
≈ 0.537 %
```

仅替换 plastic strain 为 LATIN：

```text
stress error
=
9.905354031166e-03
≈ 0.991 %
```

仅替换 damage 为 LATIN：

```text
stress error
=
5.827568710325e-05
≈ 0.0058 %
```

完整 LATIN stress：

```text
4.871414552459e-03
≈ 0.487 %
```

因此：

> damage 对当前 stress discrepancy 的直接贡献几乎可以忽略。

Plastic strain 单独的潜在影响最大，但最终 actual stress error 反而小于 total-strain-only 和 plastic-strain-only 两个单独误差。

这说明：

> total strain 与 plastic strain discrepancy 之间存在明显误差抵消。

---

# 18. FOM-2B：underlying field differences

Critical point：

```text
total strain relative L2
=
5.171760092004e-03
≈ 0.517 %
```

```text
plastic strain relative L2
=
6.839993120920e-02
≈ 6.840 %
```

```text
damage relative L2
=
4.435226371292e-02
≈ 4.435 %
```

需要特别注意：

> plastic strain 的 6.84 % relative error 不能直接理解为 stress 也应该有 6.84 % error。

原因是：

- plastic strain 只是 total strain 的一部分；
- stress 由 elastic strain 控制；
- elastic strain 是 total strain 与 plastic strain 的差。

---

# 19. FOM-2C：elastic-strain error cancellation diagnostic

## 19.1 目的

FOM-2B 提示 delta total strain 与 delta plastic strain 可能高度同向。

由于：

```text
delta elastic strain
=
delta total strain
-
delta plastic strain
```

若二者同向，则会在形成 elastic strain 时部分抵消。

FOM-2C 直接验证这一假设。

---

# 20. FOM-2C：relative L2

Critical point：

```text
total strain
=
5.171760092004e-03
≈ 0.517 %
```

```text
plastic strain
=
6.839993120920e-02
≈ 6.840 %
```

```text
elastic strain
=
4.881231135390e-03
≈ 0.488 %
```

```text
stress
=
4.871414552459e-03
≈ 0.487 %
```

关键观察：

```text
elastic strain error
≈
stress error
```

---

# 21. FOM-2C：error-vector norms

```text
||delta total||
=
2.045819968594e-05
```

```text
||delta eps_p||
=
3.776798391140e-05
```

```text
||delta elastic||
=
1.861067332032e-05
```

Elastic-strain discrepancy 明显小于：

```text
||delta total||
+
||delta eps_p||
```

说明二者并非简单叠加。

---

# 22. FOM-2C：error-direction correlation

```text
cos(
    delta_total,
    delta_eps_p
)
=
9.697614528105e-01
```

即：

```text
≈ 0.970
```

非常接近 1。

说明：

> LATIN 相对于 FOM 的 total-strain error 与 plastic-strain error 在整个时程上高度同向。

---

# 23. FOM-2C：cancellation ratio

定义诊断比例：

```text
1
-
||delta elastic||
/
(
    ||delta total||
    +
    ||delta eps_p||
)
```

得到：

```text
6.803727778379e-01
```

即：

```text
≈ 68.0 %
```

因此可以通俗地理解为：

> 约 68 % 的 total-strain / plastic-strain combined discrepancy 在形成 elastic strain 时被相减抵消。

该比例是诊断量，不应理解成严格的能量正交分解比例。

---

# 24. FOM-2C：identity check

数值恒等式：

```text
delta elastic
=
delta total
-
delta eps_p
```

identity residual：

```text
2.159380034515e-19
```

接近机器精度。

因此误差抵消判断不是后处理近似，而是直接来自状态变量本身的精确关系。

---

# 25. FOM-2C：stress isolation

Actual LATIN stress error：

```text
4.871414552459e-03
```

采用：

```text
LATIN elastic strain
+
FOM damage
```

重构 stress：

```text
4.881153594081e-03
```

二者几乎一致。

因此：

> 当前 critical-point stress discrepancy 的直接来源是 elastic-strain discrepancy，而不是 damage discrepancy。

---

# 26. FOM-2C：selected phase points

例如：

```text
t/T = 0.25
F = +1.0 MN
```

有：

```text
delta total
=
-6.34129071e-06

delta eps_p
=
-9.66298202e-06

delta elastic
=
+3.32169131e-06
```

两者同号，相减后 residual error 明显减小。

又如：

```text
t/T = 0.725
F ≈ -0.9877 MN
```

有：

```text
delta total
=
+8.37802404e-06

delta eps_p
=
+1.44549874e-05

delta elastic
=
-6.07696341e-06
```

再次显示 total strain 与 plastic strain 会共同偏移，而 elastic strain 中出现明显抵消。

---

# 27. FOM-2C 阶段结论

可以形成如下误差传递关系：

```text
integrated plastic-history discrepancy
↓
global equilibrium / compatibility
adjust total strain
↓
delta total and delta eps_p
highly aligned
↓
elastic strain
=
total strain - plastic strain
↓
large cancellation
↓
small elastic-strain discrepancy
↓
small stress discrepancy
```

这解释了为什么：

```text
plastic strain FOM error
≈ 7 %
```

而：

```text
stress FOM error
≈ 0.3 ~ 0.5 %
```

并不矛盾。

---

# 28. FOM-2D：matched temporal refinement

## 28.1 核心目的

此前已经证明：

```text
LATIN Local
vs
LATIN Global
```

history discrepancy 会随 dt 减小而近似一阶下降。

但这还不能直接证明：

```text
FOM
vs
LATIN
```

实际 discrepancy 也会随 dt 减小。

FOM-2D 因此对 40、80、160 三个完全 matched 的时间网格分别执行：

```text
FOM(dt)
vs
LATIN-PGD(dt)
```

除 dt 外，所有空间离散、材料和荷载保持一致。

---

# 29. FOM-2D：LATIN solver results

## 40 increments

```text
dt
=
0.25

termination
=
converged

iterations
=
20

rank
=
15

xi
=
8.006194e-06
```

## 80 increments

```text
dt
=
0.125

termination
=
converged

iterations
=
20

rank
=
15

xi
=
9.972060e-06
```

## 160 increments

```text
dt
=
0.0625

termination
=
converged

iterations
=
21

rank
=
15

xi
=
9.761721e-06
```

三个时间网格均正常收敛。

因此：

> refinement trend 不是由 solver failure 或 rank change 造成。

尤其是三组 rank 都为 15，更有利于将变化主要归因于时间网格。

---

# 30. FOM-2D：matched full-field relative L2

当前采用：

```text
unweighted full-time
+
full 320-material-point
relative L2
```

作为 diagnostic metric。

这不是最终正式论文 metric，但足以判断 temporal convergence trend。

| increments | dt | total strain | elastic strain | stress | plastic strain | alpha | r_bar | damage |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 40 | 0.25000 | 3.5751443e-03 | 2.8739437e-03 | 2.8711929e-03 | 7.2971948e-02 | 7.3510660e-02 | 3.4638373e-02 | 2.9736229e-02 |
| 80 | 0.12500 | 1.8593526e-03 | 1.5026355e-03 | 1.5014464e-03 | 3.7533837e-02 | 3.7939414e-02 | 1.7868402e-02 | 1.5255037e-02 |
| 160 | 0.06250 | 9.5529001e-04 | 7.7245159e-04 | 7.7190582e-04 | 1.9230626e-02 | 1.9464937e-02 | 9.3633905e-03 | 7.7407568e-03 |

---

# 31. FOM-2D：百分比表达

## total strain

```text
40:
0.358 %

80:
0.186 %

160:
0.096 %
```

## elastic strain

```text
40:
0.287 %

80:
0.150 %

160:
0.077 %
```

## stress

```text
40:
0.287 %

80:
0.150 %

160:
0.077 %
```

## plastic strain

```text
40:
7.297 %

80:
3.753 %

160:
1.923 %
```

## alpha

```text
40:
7.351 %

80:
3.794 %

160:
1.946 %
```

## r_bar

```text
40:
3.464 %

80:
1.787 %

160:
0.936 %
```

## damage

```text
40:
2.974 %

80:
1.526 %

160:
0.774 %
```

所有主要变量均满足：

```text
error(40)
>
error(80)
>
error(160)
```

没有任何反常变量。

---

# 32. FOM-2D：observed error-reduction orders

定义：

```text
p
=
log2(
    error_coarse
    /
    error_fine
)
```

结果：

| field | p(40->80) | p(80->160) |
|---|---:|---:|
| total strain | 0.94320104 | 0.96078973 |
| elastic strain | 0.93553669 | 0.95997869 |
| stress | 0.93529727 | 0.95985627 |
| plastic strain | 0.95915018 | 0.96478605 |
| alpha | 0.95425606 | 0.96281977 |
| r_bar | 0.95496050 | 0.93230771 |
| damage | 0.96293605 | 0.97873913 |

全部：

```text
p ≈ 0.93 ~ 0.98
```

且全部 monotone：

```text
True
```

这是一组非常稳定的近一阶时间误差收敛结果。

---

# 33. FOM-2D：final plastic-strain anchors

## 40 increments

```text
FOM max |eps_p|
=
6.921951176224e-05

LATIN max |eps_p|
=
6.501006613767e-05
```

## 80 increments

```text
FOM max |eps_p|
=
7.028958289300e-05

LATIN max |eps_p|
=
6.769467404535e-05
```

## 160 increments

```text
FOM max |eps_p|
=
7.055638197426e-05

LATIN max |eps_p|
=
6.914156548167e-05
```

LATIN final plastic strain 随时间细化稳定向 FOM 靠近。

---

# 34. FOM-2D：final damage anchors

## 40 increments

```text
FOM max D
=
2.535861405195e-03

LATIN max D
=
2.434540391614e-03
```

## 80 increments

```text
FOM max D
=
2.567517301383e-03

LATIN max D
=
2.513842083813e-03
```

## 160 increments

```text
FOM max D
=
2.575449350984e-03

LATIN max D
=
2.547910217283e-03
```

同样表现为稳定靠拢，而不是不规则跳动。

---

# 35. FOM-2D 最重要的结论

FOM-2D 首次直接证明：

```text
FOM
vs
LATIN-PGD
```

实际 full-field discrepancy 本身会随着：

```text
dt
→
0
```

系统下降。

而且：

```text
error
≈
O(dt)
```

对所有主要变量都成立，包括：

- total strain；
- elastic strain；
- stress；
- plastic strain；
- alpha；
- r_bar；
- damage。

因此可以较有把握地判断：

> 当前 one-cycle benchmark 上 FOM-LATIN 的主要 discrepancy 属于有限时间步下的近一阶时间离散误差，而不是存在一个不随 dt 消失的 tower LATIN-PGD 根本性偏差。

---

# 36. FOM-2A～FOM-2D 的完整因果链

当前所有证据可以统一为：

```text
Local stage
uses RK4 history integration
        │
        ▼
Global stage
uses BE for eps_p / alpha / r_bar
        │
        ▼
finite-dt integrated-history discrepancy
        │
        ▼
eps_p / alpha / r_bar
show larger relative history error
        │
        ▼
global equilibrium + compatibility
adjust total strain
        │
        ▼
delta total strain
and
delta plastic strain
are highly aligned
        │
        ▼
elastic strain
=
total strain - plastic strain
        │
        ▼
large part of history discrepancy cancels
        │
        ▼
elastic-strain error remains small
        │
        ▼
stress error remains small
        │
        ▼
Y depends nonlinearly on stress
        │
        ▼
small stress error
produces larger Y error
        │
        ▼
Ddot ∝ <Y - Y0>^2
        │
        ▼
Ddot error is further amplified
        │
        ▼
finite accumulated damage error
        │
        ▼
dt is halved
        │
        ▼
all FOM-LATIN errors
approximately halve
```

这是当前 one-cycle tower validation 中最完整的误差解释。

---

# 37. 对此前几个疑问的正式更新

## 37.1 “plastic strain 误差 7 %，是不是 LATIN 算坏了？”

当前答案：不是。

更准确地说：

> integrated plastic history 在有限 dt 下受 Local RK4 / Global BE 离散差异影响较明显，但 Global equilibrium 会同步调整 total strain，使 elastic strain 和 stress 仍保持高精度。

并且：

```text
plastic strain FOM error
7.30 %
→
3.75 %
→
1.92 %
```

随 dt 近一阶下降。

因此不支持 persistent constitutive failure 这一解释。

## 37.2 “damage 误差 3 %，是不是 damage implementation 错了？”

当前没有证据支持这一判断。

证据：

- damage activation window 完全一致；
- stored Y 与重算 Y 一致；
- stored Ddot 与 damage law 重算一致；
- Y error 几乎全部来自 stress；
- Ddot error 几乎全部由 stress error 经损伤律放大；
- damage FOM-LATIN error 随 dt：

```text
2.974 %
→
1.526 %
→
0.774 %
```

近一阶下降。

因此：

> damage discrepancy 更合理地解释为上游 stress path finite-dt error 经 nonlinear damage law 放大，而不是 damage law coding bug。

## 37.3 “是不是 PGD rank 不够？”

目前不支持这一判断。

FOM-2D 三个网格均：

```text
PGD rank
=
15
```

同时所有 FOM-LATIN errors 随 dt 系统下降。

如果主要误差来自 rank insufficiency，则不会自然得到如此整齐的 O(dt) 收敛序列。

因此：

> 当前 principal error source 是 temporal discretisation，而不是 PGD rank。

## 37.4 “是不是 outer LATIN tolerance 不够严格？”

此前 tolerance sensitivity 已经证明：

```text
tol 1e-4
→
1e-5
→
1e-6
```

对最终状态影响远小于 FOM discrepancy。

FOM-2D 又显示 dt refinement 会系统改变误差。

因此：

> outer tolerance 不是主要误差来源。

---

# 38. 当前时间离散的理论边界

虽然 FOM-2D 结果非常强，但必须保持理论表述边界。

## 可以说

- 当前 tower residual-LS LATIN-PGD 在该 one-cycle benchmark 上表现出稳定的 temporal convergence；
- FOM-LATIN discrepancy 近似满足 O(dt)；
- mixed RK4 / BE implementation 在当前 benchmark 上没有表现出 persistent dt-independent inconsistency；
- 当前 finite-dt internal-variable discrepancy 可以主要解释为 temporal discretisation error。

## 不能说

- 已经严格证明当前 tower BE discretisation 与论文 DG0 完全等价；
- 已经数学证明该收敛阶对所有 tower fatigue problem 均成立；
- 已经完成 high-cycle fatigue accuracy validation；
- 已经证明 residual-LS 是论文 Eq. (65)-(71) 的等价重写；
- 已经证明 100-cycle LATIN-PGD 可以无条件直接使用当前时间步。

---

# 39. 与原论文 DG0 的关系

原论文对 PGD temporal problem 明确提到 DG0。

但当前 tower implementation 在若干 global / temporal updates 中使用的是 Backward Euler。

因此当前证据应解释为：

> project-level engineering discretisation 在 one-cycle tower benchmark 上具有清楚的一阶 temporal consistency。

而不是：

> 当前 Backward Euler 已经被证明与论文 DG0 完全相同。

如果后续论文理论部分需要严格声称 paper-consistent DG0，仍然需要单独做离散推导或实现核对。

---

# 40. 是否应该立即把 Global BE 改成 RK4

基于当前结果：

```text
不建议。
```

原因：

1. 当前 BE-based global formulation 已经表现出稳定一阶时间收敛；
2. 误差可通过 dt refinement 系统降低；
3. 直接改为 RK4 可能破坏：
   - global-stage discrete consistency；
   - 已验证的 1D implementation inheritance；
   - Eq. (58)-(59) / global correction structure；
   - 现有 convergence behaviour；
4. “数值上看起来更接近 FOM”并不等价于“理论上更正确”。

因此下一步不应是直接把 BE 换成 RK4。

---

# 41. 是否需要立即修改 LATIN convergence indicator

当前 scalar indicator xi 更关注：

- stress；
- thermodynamic force variables；
- rate-type fields；
- compatible/global fields。

Integrated histories 不是当前主 indicator 的直接组成部分。

FOM-2D 证明：

> 即使 finite-dt integrated histories 存在较大 relative difference，这些差异会随 dt 一阶收敛。

因此当前阶段也不建议直接将所有 integrated histories 强行加入主 LATIN norm。

更稳妥的后续方案是：

```text
main LATIN indicator
+
history-consistency diagnostic
```

即：

> 保留当前 theoretical / algorithmic convergence indicator，同时增加 fatigue-oriented integrated-history monitoring。

该项属于后续 robustness enhancement，而不是本阶段必须修改的 core solver requirement。

---

# 42. FOM-2 阶段总判定

## FOM-2A

```text
PASS
```

Damage-chain discrepancy 已被定位：

```text
stress
→
Y
→
Ddot
→
D
```

## FOM-2B

```text
PASS
```

Stress discrepancy 已被分解。

Damage 对 stress discrepancy 的直接贡献近乎可忽略。

## FOM-2C

```text
PASS
```

Total-strain error 与 plastic-strain error 高度同向：

```text
cos
≈
0.970
```

并在 elastic strain 中产生明显抵消。

## FOM-2D

```text
STRONG PASS
```

所有主要 FOM-LATIN full-field errors 均 monotone decrease，且 observed order：

```text
p ≈ 0.93 ~ 0.98
```

表现出非常稳定的近一阶 temporal convergence。

---

# 43. FOM-2 总结论

本阶段最重要的结论是：

> 当前 one-cycle tower benchmark 中 FOM 与 residual-LS LATIN-PGD 的主要 discrepancy 可以由有限时间步离散误差解释，并且该 discrepancy 随时间步细化近似一阶收敛。

同时：

> 当前没有发现需要立即修改 tower LATIN-PGD 核心算法的证据。

因此 FOM-2 的意义不是把误差强行调到零，而是：

```text
证明误差来源
+
证明误差随 dt 可控收敛
```

这比单纯获得一个小误差数字更重要。

---

# 44. 当前 one-cycle accuracy 状态

以 160 increments/cycle 为例：

```text
total strain
≈
0.096 %

elastic strain
≈
0.077 %

stress
≈
0.077 %

plastic strain
≈
1.92 %

alpha
≈
1.95 %

r_bar
≈
0.94 %

damage
≈
0.77 %
```

这表明：

- global mechanical response 已具有很高精度；
- internal history discrepancy 已显著缩小；
- damage 已降至 1 % 以下；
- plastic strain / alpha 仍比 macro fields 更敏感，但趋势完全稳定。

---

# 45. 100-cycle 之前的研究路线更新

FOM-2 之后，不再需要把“是否存在根本性 one-cycle temporal inconsistency”作为主要阻塞问题。

但进入 100-cycle LATIN-PGD 前，仍应明确：

1. 采用何种 time resolution；
2. 是否需要 cycle-level acceleration；
3. 100-cycle asymmetric frozen FOM 与当前 fully reversed benchmark 并非同一加载；
4. 如何定义 fatigue-oriented validation metrics；
5. 是否需要 integrated-history diagnostic；
6. 如何控制计算成本；
7. 是否继续沿 residual-LS route。

---

# 46. 100-cycle frozen FOM 的关系

现有 frozen 100-cycle FOM：

```text
outputs/tower_100cycle_fom_reference_v1.npz
```

对应 asymmetric loading：

```text
Fmax
=
+1.0 MN

Fmin
=
-0.5 MN

RF
=
-0.5
```

而当前 FOM-2 benchmark：

```text
fully reversed
+1.0 MN
↔
-1.0 MN
```

因此：

> FOM-2 的 convergence result 不能直接与现有 100-cycle asymmetric NPZ 做数值逐点比较。

进入 long-cycle stage 时，应：

```text
adapt LATIN-PGD benchmark
to the frozen asymmetric 100-cycle loading contract
```

而不是：

```text
discard or rerun the existing 100-cycle FOM
```

---

# 47. 下一阶段建议

本阶段完成后，建议下一阶段不是立即改 core equations，而是先做 Long-cycle entry design。

主要回答：

## A. 时间分辨率

基于 one-cycle refinement：

```text
40
80
160
```

综合 accuracy、LATIN rank、computational cost 和 long-cycle storage，选择 long-cycle baseline resolution。

## B. 与 frozen 100-cycle FOM loading 对齐

将 LATIN-PGD loading contract 从 fully reversed 扩展到：

```text
asymmetric
+1.0 MN
↔
-0.5 MN
```

并严格保留 explicit preload：

```text
0
→
Fmean
```

等 FOM contract。

## C. 先做短多循环 matched validation

在进入完整 100 cycles 前，优先做：

```text
2 cycles
5 cycles
10 cycles
```

或类似短多循环 benchmark。

目的：

> 检查单循环近一阶误差是否在多循环累计时保持可控。

## D. 再进入 100-cycle LATIN-PGD

只有在：

- loading contract 对齐；
- short multi-cycle internal histories 稳定；
- damage / plastic-history error 不异常积累；
- convergence / enrichment behaviour 正常；

之后，才正式进入：

```text
100-cycle LATIN-PGD
vs
frozen 100-cycle FOM
```

---

# 48. 当前推荐的 checkpoint 结论

本阶段适合独立 Git checkpoint。

建议文档作为：

```text
docs/2026-08-21-tower-fom-temporal-consistency-stage-summary.md
```

单独提交。

推荐 commit message：

```text
docs: summarize tower FOM temporal consistency diagnostics
```

该 checkpoint 应保持：

```text
documentation only
```

不同时混入 core solver 修改。

---

# 49. 最终一句话总结

FOM-2 最重要的科学结论可以压缩成：

> 当前海上风机塔筒 residual-LS LATIN-PGD 与 FOM 在 one-cycle benchmark 中的主要误差不是不收敛的算法偏差，而是随时间步细化近似一阶消失的 temporal discretisation error；其从 internal history 到 stress 再到 nonlinear damage 的误差传播机制已经被逐层定位并数值验证。

---

# 50. 阶段状态

```text
FOM-2A
damage-chain diagnosis
=
PASS

FOM-2B
stress-error decomposition
=
PASS

FOM-2C
elastic-strain cancellation diagnosis
=
PASS

FOM-2D
matched temporal refinement
=
STRONG PASS

Overall FOM-2 status
=
CLOSED / STRONG PASS
```

下一阶段：

```text
Long-cycle entry design
and
matched asymmetric multi-cycle LATIN-PGD validation
```
