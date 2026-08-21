# 海上风机塔筒 LATIN-PGD 当前误差问题的通俗解释

日期：2026-08-20  
项目：Offshore Wind Turbine and LATIN-PGD  
当前阶段：FOM-1 结束，准备进入 FOM-2 时间离散一致性设计

> 本文档用尽量通俗的方式解释当前塔筒 LATIN-PGD 在与全阶模型 FOM 对比时遇到的问题，包括：哪里有误差、误差有多大、为什么会产生这些误差，以及为什么下一步不能直接进入 100-cycle LATIN-PGD。

---

# 1. 先用一句话说明当前问题

现在的 LATIN-PGD 已经能把塔筒：

- 整体怎么变形；
- 位移有多大；
- 应力有多大；

算得非常接近 FOM。

但是，对于材料内部：

- 累积了多少塑性变形；
- 发生了多少硬化；
- 累积了多少损伤；

还存在几个百分点的差别。

可以简单理解为：

> 外在表现已经很接近，但材料内部的“疲劳账本”还没有完全对上。

---

# 2. 当前到底哪里有误差

在 matched one-cycle fully reversed benchmark 中，LATIN-PGD 与 FOM 的主要误差如下。

## 2.1 宏观结构响应

```text
tower-top displacement   ≈ 0.33 %
total strain             ≈ 0.36 %
stress                   ≈ 0.29 %
```

这说明：

> 塔筒整体怎么弯、弯多少、应力多大，LATIN-PGD 已经和 FOM 非常接近。

这一部分目前表现很好。

---

## 2.2 材料内部历史变量

```text
plastic strain   ≈ 7.3 %
alpha            ≈ 7.35 %
r_bar            ≈ 3.46 %
damage           ≈ 3.02 %
```

这里的误差明显比位移和应力大。

也就是说：

> 当前真正需要继续解决的问题，不是塔筒整体响应，而是材料内部历史变量。

对于疲劳问题，这些变量非常重要，因此不能因为位移和应力已经很准，就直接认为整个方法已经完成验证。

---

# 3. 为什么“位移很准”，但“内部变量还差几个百分点”

可以把塔筒看成一个人。

位移和应力类似于：

> 这个人现在站在哪里、现在身上承受多大的力。

而 plastic strain、alpha、r_bar 和 damage 更像：

> 这个人过去到底累积受了多少“伤”。

前者主要和“当前状态”有关。

后者则依赖整个加载历史，需要把每一个时间步的变化不断累积起来。

因此：

> 即使每一时刻的变化速度已经很接近，只要两种方法采用不同的“累计方式”，最后得到的总历史仍然可能差几个百分点。

这正是当前出现的问题。

---

# 4. LATIN 内部实际上存在两套“记账方法”

LATIN 求解分为两个主要阶段：

```text
Local stage
+
Global stage
```

可以把它们简单理解成：

```text
Local stage：
每个材料点自己按照本构关系算“材料应该怎么发展”。

Global stage：
重新协调所有材料点，让整个塔筒满足结构平衡、兼容条件和 LATIN 搜索方向。
```

问题出现在：

> Local 和 Global 对历史变量采用了不同的时间积分方法。

---

# 5. Local stage 使用 RK4

当前 Local stage 对：

```text
plastic strain
alpha
r_bar
damage
```

采用 RK4 进行时间积分。

RK4 的特点是：

> 在一个时间步内部，不只是看结束点，而是会取多个中间信息来估计这一时间步的累积变化。

因此可以把它理解成一种比较精细的“记账方法”。

---

# 6. Global stage 对前三个变量使用 Backward Euler

当前 Global stage 对：

```text
plastic strain
alpha
r_bar
```

采用 Backward Euler，也就是后向欧拉方法。

它的基本思想可以简单理解成：

> 用这一时间步结束时的变化速度，去代表这一整个时间步的累积变化。

因此当前实际上是：

```text
Local：
RK4 history

Global：
Backward Euler history
```

虽然两边在算同一个物理量，但累计方法并不相同。

---

# 7. 用汽车行驶距离来理解 RK4 和 Backward Euler

假设汽车速度一直在变化，我们想计算 10 秒以后汽车一共走了多远。

速度对应：

```text
plastic_strain_rate
alpha_rate
r_bar_rate
```

而总路程对应：

```text
plastic_strain
alpha
r_bar
```

Backward Euler 可以粗略理解为：

> 这一秒走多远 = 这一秒结束时的速度 × 1 秒。

而 RK4 会考虑：

- 这一秒开始时的速度；
- 中间时刻的速度；
- 后面时刻的速度；

然后综合估算这一秒走了多远。

所以即使两种方法最后得到的“速度”几乎一样：

```text
rate difference ≈ 0.01 %
```

长期累积出来的“总路程”还是可能差几个百分点。

---

# 8. 我们已经验证：当前确实就是这个问题

前面做 local/global closure audit 时发现：

## 8.1 变化速率已经非常接近

```text
plastic_strain_rate difference   ≈ 0.027 %
alpha_rate difference            ≈ 0.006 %
r_bar_rate difference            ≈ 0.005 %
```

这些差异都已经很小。

说明：

> Local 和 Global 对“当前材料应该以多快速度发展”的判断已经非常接近。

---

## 8.2 但积分后的历史差异明显更大

```text
plastic_strain history gap   ≈ 8.47 %
alpha history gap            ≈ 8.53 %
r_bar history gap            ≈ 5.39 %
```

所以问题不是：

> 变化速度本身还没有收敛。

而是：

> 同样或非常接近的变化速度，经过不同时间积分方法以后，被累计成了不同的历史量。

---

# 9. 最直接的验证：用 Local rate 按 BE 重新积分

我们做了一个很关键的诊断。

做法是：

1. 取 Local stage 算出来的 rate；
2. 不再使用 Local 的 RK4 history；
3. 强制按 Backward Euler 重新积分；
4. 再与 Global history 比较。

结果如下。

## plastic strain

```text
Global vs Local
≈ 8.47 %

Global vs BE(Local rate)
≈ 0.016 %
```

## alpha

```text
Global vs Local
≈ 8.53 %

Global vs BE(Local rate)
≈ 0.013 %
```

## r_bar

```text
Global vs Local
≈ 5.39 %

Global vs BE(Local rate)
≈ 0.004 %
```

这说明：

> 一旦使用相同的积分方式，Local rate 几乎可以重构出 Global history。

因此可以比较明确地判断：

> plastic strain、alpha、r_bar 当前的主要 local/global history gap，不是因为 rate field 没有收敛，而是因为 RK4 与 Backward Euler 两种时间积分方式不一致。

---

# 10. 为什么还要做 40、80、160 时间步细化

如果误差真的是时间积分造成的，那么时间步越小，两种积分方法得到的结果应该越来越接近。

于是我们保持：

- 塔筒结构不变；
- 材料不变；
- 荷载不变；
- LATIN tolerance 不变；

只把每周期时间步数从：

```text
40
→ 80
→ 160
```

对应：

```text
dt = 0.25
→ 0.125
→ 0.0625 s
```

---

# 11. 时间步细化结果

## plastic strain

```text
40 steps     8.47 %
80 steps     4.07 %
160 steps    2.01 %
```

## alpha

```text
40 steps     8.53 %
80 steps     4.09 %
160 steps    2.03 %
```

## r_bar

```text
40 steps     5.39 %
80 steps     2.28 %
160 steps    1.06 %
```

可以非常清楚地看到：

> 时间步大约减半，误差也大约减半。

这就是非常典型的近似一阶时间离散误差。

因此现在已经有比较扎实的证据：

> plastic strain、alpha、r_bar 的主要 local/global history gap，来源于 Local RK4 与 Global Backward Euler 的时间离散不一致。

---

# 12. 为什么 damage 又不一样

damage 是当前最特殊的一个变量。

前三个变量在 Global stage 中都会重新按照 Global 的离散方式更新：

```text
plastic strain    BE
alpha             BE
r_bar             BE
```

但是 damage 当前不是这样。

当前 tower-v1 中：

```text
Global damage
≈ 直接继承 Local damage history
```

也就是说：

```text
Local damage    RK4
Global damage   inherited Local RK4 history
```

因此 LATIN 内部：

```text
Local vs Global damage gap
≈ 0.0026 %
```

几乎可以忽略。

---

# 13. 但 damage 和 FOM 仍然差约 3 %

这就是当前仍然没有解决的问题。

虽然：

```text
LATIN Local damage
≈
LATIN Global damage
```

但是：

```text
FOM damage
vs
LATIN damage
≈ 3 %
```

因此 damage 的误差不能用：

```text
RK4 vs Backward Euler
```

来解释。

也就是说，当前问题已经被拆成了两类。

---

# 14. 当前问题 A：plastic strain / alpha / r_bar

这一类问题已经基本定位。

```text
表现：
Local / Global history gap ≈ 5–8 %

主要原因：
Local RK4
vs
Global Backward Euler

证据：
1. rate fields 已经基本闭合
2. Local rate 用 BE 重积分后几乎等于 Global history
3. 40 / 80 / 160 时间细化后 gap 近似减半
```

因此：

> 这部分问题已经从“猜测”进入了“有明确数值证据支持”的阶段。

---

# 15. 当前问题 B：damage

damage 的情况不同：

```text
LATIN Local / Global damage
≈ 几乎一致
```

但：

```text
FOM vs LATIN damage
≈ 3 %
```

因此：

> damage 的约 3 % FOM 误差还存在其他原因，目前尚未完全定位。

这是下一阶段仍然需要继续追查的问题。

---

# 16. FOM 和 LATIN 本身的计算方式也不同

这一点非常重要。

## FOM

传统 FOM 的基本计算流程可以理解为：

```text
给定外荷载
↓
Newton 迭代
↓
得到结构应变
↓
每个 fiber 根据应变更新材料状态
↓
RK4 积分材料历史
↓
该荷载步收敛
↓
提交状态
```

可以简单称为：

```text
strain-driven constitutive integration
```

也就是：

> 结构先给材料一个应变历史，材料再根据应变更新状态。

---

## LATIN

LATIN 则不是按这个方式逐步 Newton。

它采用：

```text
Local stage
↕
Global stage
```

不断交替。

Local 负责：

```text
材料本构
```

Global 负责：

```text
平衡
兼容
搜索方向
PGD 修正
```

直到两边逐渐一致。

因此：

> LATIN 不是简单把 FOM 的 Newton 换成另一种求解器，而是采用了不同的 local/global 交替结构。

最终我们需要保证的是：

> 虽然两种算法走的路线不同，但在离散足够细、LATIN 充分收敛时，两者应该趋向同一个物理解。

当前宏观响应已经非常接近，但内部 history 还没有完全做到这一点。

---

# 17. PGD 是不是现在最大的嫌疑

目前来看，不是第一嫌疑。

前面已经验证：

```text
塔筒平衡算子正确
PGD 模态满足平衡
residual-LS enrichment 可以稳定工作
PGD rank 可以正常增长
outer LATIN 可以正常收敛
transaction semantics 正常
tolerance 收紧后结果已经基本稳定
rate fields local/global 已经基本闭合
```

如果主要问题是：

```text
PGD 模态太少
```

那么继续增加 rank 或收紧 LATIN tolerance 应该显著改变结果。

但实际发现：

```text
tol = 1e-4
→ 1e-5
→ 1e-6
```

最终内部状态变化已经非常小。

反而：

```text
dt 减半
```

会让 history gap 接近减半。

所以目前更明确的证据指向：

```text
时间离散
```

而不是：

```text
PGD rank 不够
```

---

# 18. 为什么 LATIN 显示已经收敛，但 history 还可以差 8 %

这也是最容易让人困惑的一点。

当前 LATIN convergence indicator xi 主要检查的是：

```text
stress
beta
R_bar
elastic_strain
plastic_strain_rate
alpha_rate
r_bar_rate
```

它更关注：

> 当前应力、热力学力和变化速率是否一致。

但是它没有直接强制检查：

```text
plastic_strain history
alpha history
r_bar history
damage history
```

因此可以出现：

```text
rate 已经很接近
↓
xi 很小
↓
LATIN 判定收敛
```

但因为：

```text
Local 用 RK4
Global 用 BE
```

所以：

```text
把这些 rate 积分成 history
↓
仍然可以差 5–8 %
```

因此必须明确：

> xi 很小，说明当前 LATIN 定义下 local/global 搜索已经收敛；它不代表所有累计历史变量都已经与 FOM 完全一致。

---

# 19. 当前问题最简单的逻辑图

```text
                     同一个塔筒循环问题
                            │
                 ┌──────────┴──────────┐
                 │                     │
                FOM                 LATIN-PGD
                 │                     │
          Newton + RK4            Local + Global
                                       │
                         ┌─────────────┴─────────────┐
                         │                           │
                       Local                       Global
                         │                           │
                    RK4 history                BE history
                                               eps_p / alpha / r_bar
                         │
                    damage RK4
                                                     │
                                    时间积分定义不同
                                                     │
                               eps_p / alpha / r_bar
                                   出现 5–8 % gap
                                                     │
                                   dt 减半
                                                     │
                                   gap 约减半
                                                     │
                                  主要原因已定位
```

另一条 damage 支线：

```text
LATIN Local damage
≈
LATIN Global damage
        │
        │
        但
        ↓
FOM vs LATIN damage
≈ 3 %
        │
        ↓
原因仍待进一步定位
```

---

# 20. 所以现在真正的问题是什么

当前已经不是：

```text
矩阵是不是写错
索引是不是错
平衡算子是不是错
PGD 是不是完全不收敛
```

这些基础问题已经基本排除。

现在真正的问题上升到了算法设计层面：

> LATIN 的 Local stage 和 Global stage 到底应该采用怎样一致的时间离散？

我们接下来需要决定：

```text
方案 A：
Local 继续 RK4
Global 继续 BE
并接受有限 dt 下的 history mismatch

方案 B：
Local 和 Global 统一采用同一种离散方式

方案 C：
重新回到原论文时间离散框架
从 DG0 / 时间弱式出发统一 Local 和 Global
```

但现在还不能因为：

```text
RK4 看起来更精细
```

就直接把 Global BE 改成 RK4。

因为：

> 让两个结果数值上更接近，不代表理论上就是正确的 LATIN-PGD 离散。

---

# 21. 为什么下一步不能直接进入 100-cycle LATIN-PGD

如果现在一循环中：

```text
plastic strain history gap
≈ 8 %
```

那么直接进入 100 个循环以后：

> 这种 history discrepancy 可能继续累积。

而我们的最终目标恰恰是：

```text
high-cycle fatigue
```

所以如果单循环 history convention 还没有厘清，就直接上 100-cycle，最终即使得到一条损伤曲线，也很难判断：

```text
它是物理疲劳累积
还是数值积分误差累积
```

因此必须先解决单循环时间离散一致性问题。

---

# 22. 当前最需要记住的四句话

## 第一句

> 当前 LATIN-PGD 已经可以非常准确地复现塔筒的位移、总应变和应力，误差约 0.3 %。

## 第二句

> 当前没有完全对上的，是 plastic strain、alpha、r_bar 和 damage 等材料内部历史变量，其 FOM difference 约为 3–7 %。

## 第三句

> 对 plastic strain、alpha 和 r_bar，我们已经基本查明主要原因：Local 使用 RK4，而 Global 使用 Backward Euler；时间步减半后，该差异也大约减半。

## 第四句

> damage 的约 3 % FOM difference 不能由上述 RK4/BE mismatch 解释，因此仍需单独排查。

---

# 23. 下一阶段真正要解决的问题

因此下一阶段 FOM-2 的核心不是：

```text
继续增加 PGD rank
```

也不是：

```text
继续收紧 LATIN tolerance
```

更不是：

```text
直接跑 100 cycles
```

而是：

```text
重新厘清和统一：
Local stage 的时间离散
Global stage 的时间离散
Eq. (58)-(59) 的时间更新
Eq. (73)-(75) 的历史变量更新
damage 的时间积分
LATIN convergence indicator 是否需要覆盖 integrated histories
```

最终目标是：

> 先建立一套理论上自洽、数值上可验证、尽可能接近原论文 LATIN-PGD 框架的统一时间离散方案，再重新进行 one-cycle FOM comparison。

---

# 24. 当前阶段的最终理解

可以把现在的状态概括为：

```text
结构整体响应：
已经很好
≈ 0.3 % FOM error

PGD / outer LATIN：
已经能够稳定工作

plastic strain / alpha / r_bar：
主要误差来源已经定位
→ Local RK4 vs Global BE

damage：
仍有约 3 % FOM difference
→ 原因尚未完全定位

下一步：
FOM-2 时间离散一致性设计

暂时不进入：
100-cycle LATIN-PGD
```

最终可以用一句话概括：

> 当前问题已经从“LATIN-PGD 能不能在塔筒上跑起来”，转变为“LATIN-PGD 的 Local 与 Global 时间离散怎样统一，才能让材料内部疲劳历史也与 FOM 一致”。

这正是进入高周疲劳计算之前必须解决的最后一类基础问题之一。
