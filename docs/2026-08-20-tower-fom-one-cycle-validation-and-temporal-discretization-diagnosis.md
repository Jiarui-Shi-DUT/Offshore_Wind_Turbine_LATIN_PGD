# 海上风机塔筒 LATIN-PGD 一循环 FOM 验证与时间离散诊断

日期：2026-08-20  
当前分支：`feature/offshore-wind-turbine-tower-fatigue`  
前一阶段闭环提交：`cd4355a354315ac6ec57729b1881c8bd9de506f0`

> 本文档记录塔筒 LATIN-PGD 在完成 I-0 至 I-6 内部数值验证之后，第一次进入与全阶非线性有限元模型 FOM 的直接精度对比阶段。
>
> 本阶段不修改 LATIN-PGD 核心算法，而是通过 matched one-cycle fully reversed benchmark，依次完成问题定义一致性、材料点映射、返回状态物理意义、FOM 误差、LATIN 容差敏感性、local/global history closure、Backward-Euler 与 RK4 离散差异以及时间步细化诊断。
>
> 为保证 PyCharm Markdown 稳定显示，本文档继续采用普通 Unicode 数学符号和代码块表示公式，不使用行内 LaTeX 语法。

---

# 1. 本阶段的核心问题

I-6 已经证明当前 residual-LS 路径能够：

- 保持塔筒平衡；
- 保持 PGD 模态平衡；
- 稳定完成空间富集；
- 保证 Trial-A / Trial-B 事务一致性；
- 完成完整 outer LATIN 收敛；
- 通过联合回归测试。

但是这些结果只能说明：

> 当前塔筒 LATIN-PGD 在自身算法框架内部已经形成稳定、可收敛的求解链条。

它还不能回答更重要的问题：

> LATIN-PGD 得到的塔筒非线性循环响应，与完整 Newton 全阶模型 FOM 到底有多接近？

因此，本阶段正式进入 FOM accuracy validation。

---

# 2. 本阶段 benchmark

为了避免直接使用已有 100-cycle asymmetric FOM 而造成荷载条件不一致，本阶段首先采用与 I-5 完全匹配的一循环 fully reversed benchmark。

## 2.1 结构离散

```text
塔筒单元数                         = 10
每单元 Gauss 点数                 = 2
每 Gauss 点环向 fiber 数          = 16
径向 fiber 数                     = 1
材料点总数 Nq                     = 320
结构总自由度 ndof                 = 33
自由自由度数                      = 30
```

材料点 canonical coordinate 为：

```text
q <-> (element, Gauss point, fiber)
```

并采用：

```text
element-major
→ Gauss-major
→ fiber-major
```

即：

```text
q = (element * n_gauss + gauss) * n_fibers + fiber
```

---

## 2.2 荷载

采用 fully reversed 正弦循环水平荷载：

```text
Fmax                    = +1.0 MN
Fmin                    = -1.0 MN
load ratio              = -1
period                  = 10 s
cycles                  = 1
increments / cycle      = 40
Nt                      = 41
dt                      = 0.25 s
```

完整路径：

```text
0
→ +1.0 MN
→ 0
→ -1.0 MN
→ 0
```

这一 benchmark 与当前 residual-LS outer LATIN-PGD I-5 benchmark 完全一致。

---

# 3. FOM-LATIN benchmark contract audit

在进行任何误差比较之前，首先验证 FOM 和 LATIN-PGD 是否真正求解同一个问题。

检查结果如下：

```text
material parameters identical              = True
time grid identical                        = True
load vectors identical                     = True
n_dof identical                            = True
n_material_points identical                = True
free DOFs identical                        = True
metric weights identical                   = True
reference modulus identical                = True
compatibility matrix identical             = True
reduced stiffness identical                = True
stress-to-force factor identical           = True
FOM physical initial strain zero           = True
FOM physical initial stress zero           = True
FOM physical internal state zero           = True
LATIN t0 plastic strain zero               = True
LATIN t0 damage zero                       = True
LATIN t0 stress zero                       = True
```

尺寸：

```text
Nt                          = 41
ndof                        = 33
Nq                          = 320
compatibility matrix        = (320, 30)
metric                      = (320,)
```

LATIN elastic initialization：

```text
maximum free equilibrium residual
= 2.361617263146e-06 N
```

结论：

```text
FOM-LATIN-PGD BENCHMARK CONTRACT AUDIT = PASS
```

因此，后续误差不是由材料参数、荷载、空间离散或材料点数量不一致造成的。

---

# 4. FOM-LATIN 材料点映射验证

第二步验证 FOM tensor 与 LATIN canonical q 是否能够直接比较。

检查：

```text
layout.flatten              = True
layout.unflatten            = True
FOM strain snapshot shape   = True
FOM state snapshot shape    = True
direct reshape to q         = True
q reshape round trip        = True
```

代表性映射：

```text
(0, 0, 0)  -> q =   0
(0, 0, 15) -> q =  15
(0, 1, 0)  -> q =  16
(1, 0, 0)  -> q =  32
(9, 1, 15) -> q = 319
```

因此可以直接使用：

```text
FOM:
(Nt, element, Gauss, fiber)

reshape

(Nt, Nq)
```

与 LATIN：

```text
(Nt, Nq)
```

逐材料点比较。

物理场映射为：

```text
FOM fiber_strains
<-> LATIN elastic_strain + plastic_strain

FOM fiber_stresses
<-> LATIN stress

FOM fiber_states[..., 0]
<-> LATIN plastic_strain

FOM fiber_states[..., 1]
<-> LATIN alpha

FOM fiber_states[..., 2]
<-> LATIN r_bar

FOM fiber_states[..., 3]
<-> LATIN damage
```

结论：

```text
FOM-LATIN-PGD FIELD-MAPPING AUDIT = PASS
```

---

# 5. LATIN 返回状态的物理意义

在进行 FOM comparison 前进一步确认：

```text
result.state
```

到底代表什么。

当前 transactional outer solver 中：

```text
accepted_state
accepted_basis
accepted_indicator
```

组成唯一 persistent snapshot。

最终返回：

```text
result.state
```

即最后正式 commit 的 relaxed global state。

本 benchmark 中：

```text
converged              = True
termination            = converged
committed iterations   = 12
final commit kind      = B
PGD rank               = 7
final LATIN indicator  = 9.424149148753e-05
```

最终返回状态与最后 committed relaxed state 的差异：

```text
0.000000000000e+00
```

---

# 6. 返回 global state 的可容许性

返回状态的外荷载平衡：

```text
maximum equilibrium residual
= 2.421520189770e-06 N

relative equilibrium residual
= 2.421520189770e-12
```

总应变：

```text
εtotal = εelastic + εplastic
```

其 compatibility residual：

```text
maximum compatibility residual
= 1.349831704744e-16

relative compatibility residual
= 1.474202370850e-13
```

重新计算最终 LATIN indicator：

```text
recomputed xi
= 9.424149148753e-05

stored xi
= 9.424149148753e-05
```

fresh local projection：

```text
xi
= 8.570337891448e-05
```

结论：

```text
RETURNED GLOBAL STATE AUDIT = PASS
```

因此：

> 后续 FOM 对比正式采用 `result.state` 作为 LATIN-PGD 解。

---

# 7. 一循环 FOM 结果

matched fully reversed FOM 采用完整 nonlinear Newton fiber-beam model。

主要结果：

```text
stress reversal                         = True
plastic-flow direction reversal         = True

final residual tower-top displacement
= -2.148196553228e-02 m

final max |eps_p|
= 6.921951176224e-05

final max damage
= 2.535861405195e-03

maximum Newton iterations
= 4

maximum free-DOF residual
≈ 9.98e-03 N
```

因此 FOM 本身数值收敛正常，同时正确出现应力反转与塑性流动方向反转。

---

# 8. FOM-1：第一次直接精度比较

采用：

```text
LATIN tolerance             = 1e-4
spatial strategy            = residual_ls
outer iterations            = 12
PGD rank                    = 7
final xi                    = 9.424149148753e-05
```

得到：

| 场变量 | relative L2 | relative max | max absolute error |
|---|---:|---:|---:|
| tower-top displacement | 3.276129007e-03 | 5.954668911e-03 | 6.954501183e-03 m |
| total strain | 3.574901030e-03 | 9.273589426e-03 | 8.447269966e-06 |
| stress | 2.873035723e-03 | 8.420707628e-03 | 9.764955117e-01 MPa |
| plastic strain | 7.298878773e-02 | 1.191522390e-01 | 1.489250011e-05 |
| alpha | 7.351178231e-02 | 1.236411042e-01 | 1.520170099e-05 |
| r_bar | 3.457941115e-02 | 3.463367476e-02 | 1.560944008e-05 |
| damage | 3.016969403e-02 | 4.296578785e-02 | 1.089552832e-04 |

主要特征非常明确：

```text
宏观结构响应误差
≈ 0.3 %

应力场误差
≈ 0.3 %

塑性与硬化 history
≈ 3–7 %

damage
≈ 3 %
```

因此：

> 当前 LATIN-PGD 已经能够很好地复现塔筒整体力学响应，但内部历史变量的精度仍明显低于位移、总应变和应力。

---

# 9. 塔顶位移关键相位比较

```text
t/T = 0.00
FOM    =  0.000000000e+00 m
LATIN  =  0.000000000e+00 m

t/T = 0.25
FOM    =  1.167907282e+00 m
LATIN  =  1.173401091e+00 m
diff   =  5.493809031e-03 m

t/T = 0.50
FOM    =  5.427117921e-02 m
LATIN  =  5.399591941e-02 m
diff   = -2.752598028e-04 m

t/T = 0.75
FOM    = -1.127748322e+00 m
LATIN  = -1.134187562e+00 m
diff   = -6.439240114e-03 m

t/T = 1.00
FOM    = -2.148196553e-02 m
LATIN  = -2.045157127e-02 m
diff   =  1.030394263e-03 m
```

峰值位移误差约为数毫米。

---

# 10. 最终状态比较

```text
FOM residual top ux
= -2.148196553228e-02 m

LATIN residual top ux
= -2.045157126879e-02 m
```

塑性应变：

```text
FOM final max |eps_p|
= 6.921951176224e-05

LATIN final max |eps_p|
= 6.504738665463e-05
```

损伤：

```text
FOM final max D
= 2.535861405195e-03

LATIN final max D
= 2.435685722670e-03
```

---

# 11. FOM 临界材料点

FOM peak-damage location：

```text
(element, Gauss, fiber)
= (0, 0, 12)

canonical q
= 12
```

最终值：

```text
stress:

FOM
= -5.312476318393 MPa

LATIN
= -4.965021814491 MPa
```

```text
eps_p:

FOM
= 6.835266117279e-05

LATIN
= 6.430877869873e-05
```

```text
damage:

FOM
= 2.535861405195e-03

LATIN
= 2.435685722670e-03
```

---

# 12. LATIN tolerance sensitivity

为了判断 3–7 % history error 是否只是 outer LATIN 收敛不充分导致，进一步采用：

```text
tolerance = 1e-4
tolerance = 1e-5
tolerance = 1e-6
```

结果：

```text
tol = 1e-4

termination = converged
iterations  = 12
rank        = 7
xi          = 9.424149148753e-05
```

```text
tol = 1e-5

termination = converged
iterations  = 20
rank        = 15
xi          = 8.006194456505e-06
```

```text
tol = 1e-6

termination = stagnated
iterations  = 25
rank        = 20
xi          = 2.678408345949e-06
```

相对于最严格 `1e-6` 解：

| tolerance | stress | eps_p | alpha | r_bar | damage |
|---|---:|---:|---:|---:|---:|
| 1e-4 | 8.96180511e-05 | 1.88317421e-03 | 2.32431600e-03 | 7.92309564e-04 | 9.05187643e-04 |
| 1e-5 | 5.73177247e-06 | 1.18949115e-04 | 1.55251597e-04 | 5.34578676e-05 | 5.26470714e-05 |
| 1e-6 | 0 | 0 | 0 | 0 | 0 |

最终：

```text
tol=1e-4

max |eps_p|
= 6.504738665463e-05

max D
= 2.435685722670e-03
```

```text
tol=1e-5

max |eps_p|
= 6.501006613767e-05

max D
= 2.434540391614e-03
```

```text
tol=1e-6

max |eps_p|
= 6.501516788022e-05

max D
= 2.434630416266e-03
```

结论：

> 继续收紧 outer LATIN tolerance 只会使最终内部变量发生约 1e-4 至 1e-3 相对量级的变化，远小于 FOM comparison 中 3–7 % 的误差。

因此：

```text
FOM history error
≠ outer LATIN convergence tolerance error
```

---

# 13. LATIN field-wise local/global closure

采用：

```text
tolerance = 1e-5
iterations = 20
rank = 15
xi = 8.006194456505e-06
```

对最终 returned global state 再进行 fresh local constitutive projection。

结果：

| field | local/global relative L2 |
|---|---:|
| elastic_strain | 1.486090298010e-08 |
| plastic_strain | 8.469947760328e-02 |
| plastic_strain_rate | 2.675053314149e-04 |
| alpha | 8.525719320267e-02 |
| alpha_rate | 6.316769995523e-05 |
| r_bar | 5.387336898273e-02 |
| r_bar_rate | 4.738292656417e-05 |
| damage | 2.627964300630e-05 |
| damage_rate | 3.594333942928e-05 |

这组结果揭示出非常重要的结构：

```text
rate fields:
基本闭合

integrated history:
eps_p   ≈ 8.47 %
alpha   ≈ 8.53 %
r_bar   ≈ 5.39 %

damage:
基本闭合
```

与此同时：

```text
max |beta - C alpha|
= 2.220446049250e-16

max |R_bar - R_inf r_bar|
= 1.734723475977e-18
```

说明 global hardening state law 本身严格闭合。

---

# 14. 当前 LATIN indicator 的解释边界

当前 scalar LATIN indicator 主要检查：

```text
stress
beta
R_bar
plastic_strain_rate
elastic_strain
alpha_rate
r_bar_rate
```

它并不直接检查：

```text
plastic_strain
alpha
r_bar
damage history
```

因此可能出现：

```text
LATIN xi
≈ 1e-5
```

同时：

```text
plastic_strain local/global history gap
≈ 8 %
```

这不是逻辑矛盾，而是当前 convergence norm 对 integrated histories 的覆盖范围有限。

因此：

> LATIN indicator 很小不等价于 FOM error 很小，也不等价于所有 integrated histories 已经 local/global 一致。

这一点在后续长期 fatigue accumulation 中必须特别注意。

---

# 15. Backward-Euler history closure

进一步检查 global 和 local history 是否满足：

```text
x_n - x_(n-1)
=
dt * xdot_n
```

## 15.1 GLOBAL

```text
plastic_strain BE residual L2
= 3.1005864858e-15

alpha BE residual L2
= 5.4707197946e-16

r_bar BE residual L2
= 1.0734780216e-15
```

因此：

```text
global eps_p
global alpha
global r_bar
```

严格满足 backward-Euler history relation。

damage：

```text
BE residual L2
= 1.9076666248e-01
```

说明 global damage 并不是 BE-integrated history。

---

## 15.2 LOCAL

```text
plastic_strain BE residual L2
= 1.9814891516e-01

alpha BE residual L2
= 1.9827348052e-01

r_bar BE residual L2
= 1.9816727701e-01

damage BE residual L2
= 1.9076620684e-01
```

这是合理的，因为 local constitutive stage 使用 RK4 顺序积分，而不是 BE nodal integration。

---

# 16. 用 local rate 重新构造 BE history

最关键的诊断是：

> 使用 fresh local nodal rate，但不使用 local RK4 history，而是人为按 backward Euler 重新积分，会得到什么？

结果：

| history | global vs local | global vs BE(local rate) | local vs BE(local rate) |
|---|---:|---:|---:|
| plastic_strain | 8.469947760328e-02 | 1.646361690413e-04 | 8.253115032980e-02 |
| alpha | 8.525719320267e-02 | 1.294231124620e-04 | 8.305684740377e-02 |
| r_bar | 5.387336898273e-02 | 4.134611842965e-05 | 5.175129985398e-02 |
| damage | 2.627964300630e-05 | 5.516576404640e-02 | 5.515825523947e-02 |

同时 rate differences：

```text
plastic_strain_rate
= 2.675053314149e-04

alpha_rate
= 6.316769995523e-05

r_bar_rate
= 4.738292656417e-05

damage_rate
= 3.594333942928e-05
```

因此，对于前三个变量：

```text
global vs local history
≈ 5–8 %
```

但：

```text
global vs BE(local rate)
≈ 0.004–0.016 %
```

这给出了非常强的数值证据：

> `eps_p`、`alpha`、`r_bar` 的 local/global history gap 主要来自不同时间积分定义，而不是 rate fields 没有收敛。

---

# 17. tower-v1 当前 history convention

根据实际代码和上述数值结果：

```text
LOCAL

eps_p    RK4 history
alpha    RK4 history
r_bar    RK4 history
damage   RK4 history
```

而：

```text
GLOBAL

eps_p    backward-Euler-consistent history
alpha    backward-Euler-consistent history
r_bar    backward-Euler-consistent history
damage   inherited local RK4 history
```

因此当前 tower-v1 混合了两类时间离散 history。

这一点对于一循环 benchmark 已经会产生几个百分点的 integrated-history gap。

对于未来几十至几百个循环的 fatigue accumulation，这种差异可能继续累积，因此不能简单忽略。

---

# 18. 与已验证一维实现的关系

已验证一维 LATIN-PGD global stage 中：

```text
eps_p
alpha
r_bar
damage
```

均采用 global-stage 离散时间关系进行更新。

当前 tower-v1 中前三个变量延续了 backward-Euler global history，但 damage 采用了：

```text
damage = local_state.damage
```

即直接继承 local RK4 history。

因此 tower-v1 damage finishing 是一个额外的工程化选择。

该选择目前：

- 内部 local/global damage closure 很好；
- 但与 FOM 仍存在约 3 % 的 damage error。

因此不能认为 damage 问题已经解决。

---

# 19. 时间步细化验证

为了判断前三个 history gap 是否真的是有限时间步误差，保持：

```text
结构不变
材料不变
空间离散不变
荷载周期不变
LATIN tolerance = 1e-5
```

只改变：

```text
increments/cycle
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

得到：

| increments | dt | termination | xi | iter | rank | eps_p gap | alpha gap | r_bar gap | D gap |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 40 | 0.25000 | converged | 8.006194e-06 | 20 | 15 | 8.46994776e-02 | 8.52571932e-02 | 5.38733690e-02 | 2.62796430e-05 |
| 80 | 0.12500 | converged | 9.972060e-06 | 20 | 15 | 4.06651408e-02 | 4.09432935e-02 | 2.28155312e-02 | 2.75873901e-05 |
| 160 | 0.06250 | converged | 9.761721e-06 | 21 | 15 | 2.01486480e-02 | 2.02854473e-02 | 1.05533633e-02 | 1.26191898e-05 |

时间步减半后的 gap ratio：

```text
40 -> 80

eps_p   = 0.480111
alpha   = 0.480233
r_bar   = 0.423503
```

```text
80 -> 160

eps_p   = 0.495477
alpha   = 0.495452
r_bar   = 0.462552
```

这表明：

```text
eps_p gap
8.47 %
→ 4.07 %
→ 2.01 %
```

```text
alpha gap
8.53 %
→ 4.09 %
→ 2.03 %
```

```text
r_bar gap
5.39 %
→ 2.28 %
→ 1.06 %
```

---

# 20. 时间离散收敛阶

对于 eps_p：

```text
p40->80
≈ log(8.4699 / 4.0665) / log(2)
≈ 1.06
```

```text
p80->160
≈ log(4.0665 / 2.0149) / log(2)
≈ 1.01
```

alpha 基本相同。

因此可以得到非常明确的数值判定：

> `eps_p` 和 `alpha` 的 local/global integrated-history gap 呈近似一阶时间收敛。

r_bar 同样呈明显时间收敛，阶数略高于 1 后逐渐接近 1。

这与 backward Euler 的一阶离散特征完全一致。

---

# 21. 当前已经可以确认的结论

## 21.1 benchmark 一致性已经排除

FOM 与 LATIN：

```text
材料
荷载
时间网格
有限元网格
材料点
积分权重
compatibility matrix
reference stiffness
```

均完全一致。

所以 FOM difference 不是 benchmark mismatch。

---

## 21.2 材料点映射已经排除

FOM tensor 与 LATIN canonical q 可直接一一对应。

所以误差不是 flatten / reshape / indexing 错误。

---

## 21.3 LATIN 返回 state 是正确 persistent global state

最终 `result.state`：

- 是最终正式 commit 的 relaxed global state；
- 满足外荷载平衡；
- 满足 compatibility；
- 与最终 stored indicator 一致。

所以误差不是因为比较了错误的 LATIN state。

---

## 21.4 outer LATIN tolerance 不是主要误差来源

从：

```text
1e-4
→ 1e-5
→ 1e-6
```

LATIN 解很快稳定。

因此 FOM 中 3–7 % 的 history error 无法由 outer iteration tolerance 解释。

---

## 21.5 rate fields 已基本收敛

当前最终状态：

```text
plastic_strain_rate
alpha_rate
r_bar_rate
damage_rate
```

local/global 差异已经约为：

```text
1e-4 至 1e-5
```

量级。

因此主要问题不是 rate field convergence。

---

## 21.6 eps_p、alpha、r_bar history gap 来源已经基本定位

通过：

1. global BE closure；
2. local RK4 closure；
3. local-rate -> BE-history reconstruction；
4. 40 / 80 / 160 时间步细化；

四组独立证据可以确认：

> 当前 `eps_p`、`alpha`、`r_bar` 的明显 local/global history gap，主要来自 local RK4 history 与 global backward-Euler history 的离散定义不统一。

并且这一 gap 随 dt 近似一阶下降。

---

# 22. 当前不能得出的结论

以下结论现阶段仍然不能声称。

## 22.1 不能说 FOM validation 已完全通过

目前：

```text
displacement
strain
stress
```

误差约 0.3 %，表现很好。

但：

```text
eps_p
alpha
r_bar
damage
```

仍有 3–7 % FOM difference。

对于 fatigue analysis，这些内部变量非常重要。

因此不能因为宏观响应接近，就宣布整个 FOM accuracy validation 完成。

---

## 22.2 不能说 residual-LS 已经获得 FOM 等价解

residual-LS 已经解决第四模态固定点失稳，并使 outer LATIN 收敛。

但是：

> 数值收敛到 LATIN 离散方程的解，不自动等价于与 FOM 完全一致。

两类结论必须分开。

---

## 22.3 不能把所有 FOM history error 都归因于 BE/RK4 mismatch

对于 eps_p、alpha、r_bar，当前证据非常强。

但是 damage：

```text
LATIN local/global gap
≈ 2.6e-05
```

同时：

```text
FOM vs LATIN damage error
≈ 3 %
```

因此 damage 的 FOM error 显然还存在其他来源。

---

## 22.4 不能直接把 global BE 全部替换成 RK4

当前 PGD temporal problem 本身采用 backward-Euler-like nodal update。

已验证的一维实现同样采用这一类 global discrete history。

因此不能为了使 local/global history 完全相同，就简单把 global history 改成 local RK4。

后续必须首先从统一离散理论出发决定修改方向。

---

# 23. 当前阶段的总体判定

可以把本阶段分成两个层次。

## 23.1 宏观 FOM accuracy

```text
tower-top displacement   ≈ 0.33 %
total strain             ≈ 0.36 %
stress                   ≈ 0.29 %
```

评价：

```text
STRONG
```

---

## 23.2 内部 history accuracy

```text
plastic strain           ≈ 7.3 %
alpha                    ≈ 7.35 %
r_bar                    ≈ 3.46 %
damage                   ≈ 3.02 %
```

评价：

```text
NOT YET CLOSED
```

---

## 23.3 时间离散诊断

```text
eps_p / alpha / r_bar
local-global history mismatch
```

已经表现出非常明确的：

```text
O(dt)
```

收敛趋势。

评价：

```text
ROOT CAUSE STRONGLY IDENTIFIED
```

更严谨地表述为：

> 对 `eps_p`、`alpha` 和 `r_bar`，当前主要 local/global history discrepancy 已被定位为 local RK4 与 global backward-Euler 时间离散 convention 不一致所导致的一阶时间离散误差。

---

# 24. 下一阶段建议

下一阶段不应立即进入 100-cycle LATIN-PGD 计算。

更合理的顺序是：

## FOM-2：时间离散一致性设计

需要重新回到：

```text
local stage
global stage
Eq. (58)-(59)
Eq. (73)-(75)
integrated histories
LATIN convergence indicator
```

逐项厘清：

1. 原论文中各 history / rate 的离散关系；
2. 当前 1D reproduction 实际采用的离散关系；
3. tower-v1 当前采用的离散关系；
4. 哪些部分必须保持与论文一致；
5. 哪些部分是当前工程实现；
6. local 与 global 是否必须采用同一种 nodal history convention；
7. damage 应继续继承 RK4 local history，还是应回到与 global rate 一致的离散更新；
8. integrated-history discrepancy 是否需要进入 LATIN convergence criterion。

在完成这些理论与数值选择之前，不建议直接修改核心代码。

---

# 25. 后续 FOM 验证路线

推荐顺序：

```text
FOM-1
matched one-cycle accuracy
+
temporal-discretization diagnosis
        ↓
当前阶段
```

```text
FOM-2
统一时间离散原则
        ↓
```

```text
FOM-3
重新进行 matched one-cycle FOM comparison
        ↓
```

重点确认：

```text
displacement
stress
eps_p
alpha
r_bar
damage
```

是否同时收敛。

之后再进入：

```text
FOM-4
time-grid convergence against FOM
```

以及：

```text
FOM-5
existing 100-cycle asymmetric frozen FOM
```

最后才进入：

```text
high-cycle / cycle-separable LATIN-PGD
```

---

# 26. 与已有 100-cycle FOM 的关系

项目中已经存在冻结的：

```text
outputs/tower_100cycle_fom_reference_v1.npz
```

其 benchmark 为：

```text
Fmax = +1.0 MN
Fmin = -0.5 MN
R_F  = -0.5
cycles = 100
increments/cycle = 40
```

该数据已经完成 integrity audit，因此：

> 后续绝不能重新计算或丢弃这一 100-cycle FOM reference。

但是它与本阶段 fully reversed：

```text
+1.0 MN
↔
-1.0 MN
```

不是同一个 load case。

因此当前一循环 fully reversed FOM validation 完成之前，不应把两套结果直接逐点比较。

---

# 27. 当前研究边界

到本阶段结束，可以安全陈述：

> 塔筒 residual-LS LATIN-PGD 已通过内部平衡、富集、事务与 outer convergence 验证，并在 matched fully reversed one-cycle benchmark 中以约 0.3 % 的 relative L2 error 复现 FOM 的宏观位移、总应变与应力响应。进一步诊断表明，plastic strain、alpha 与 r_bar 的主要 local/global history discrepancy 来源于 local RK4 与 global backward-Euler 时间离散 convention 不一致，该 discrepancy 随时间步细化呈近似一阶收敛。

同时必须附带：

> 内部变量的 FOM accuracy 尚未完全闭合，尤其 damage 的约 3 % FOM difference 尚不能由当前 BE/RK4 history mismatch 解释。因此本阶段不能视为完整 FOM validation 的最终结束点。

---

# 28. 本阶段最终状态

```text
Benchmark contract audit                       PASS
Material-point mapping audit                   PASS
Returned global-state audit                    PASS
One-cycle FOM macro-response comparison        STRONG
Outer-tolerance sensitivity diagnosis          CLOSED
Rate-field local/global closure                STRONG
BE/RK4 history mismatch identification         STRONG
40/80/160 temporal-refinement evidence         STRONG
Internal-history FOM accuracy                  NOT YET CLOSED
Damage FOM discrepancy                         OPEN
100-cycle LATIN-PGD comparison                 NOT STARTED
```

因此下一正式阶段应为：

```text
FOM-2
tower LATIN-PGD temporal discretisation consistency
```

而不是直接开始 100-cycle LATIN-PGD 求解。
