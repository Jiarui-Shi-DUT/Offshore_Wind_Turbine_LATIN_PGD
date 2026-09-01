# OPT-5 阶段总结：塔筒 PGD 空间—时间固定点容差跨周期校准

**日期：2026-09-01**
**项目：Offshore Wind Turbine and LATIN-PGD**
**仓库：`Jiarui-Shi-DUT/Offshore_Wind_Turbine_LATIN_PGD`**
**当前分支：`perf/tower-local-stage-optimization`**
**OPT-5 代码提交：`8224819`，提交信息：`perf: calibrate tower PGD fixed-point tolerance`**

---

# 1. 阶段定位

OPT-5 是在 OPT-4 完成之后，对塔筒 LATIN-PGD 当前新性能热点开展的第一轮数值参数校准。

OPT-1 至 OPT-4 的共同特征是：

- 不修改 LATIN-PGD 数学算法；
- 不修改 PGD 富集停止精度；
- 不修改外层 LATIN 收敛判据；
- 主要针对 Python 实现中的高频标量调用、重复派生量、局部本构热路径以及材料点方向批量化开展实现级优化。

OPT-4 完成后，局部本构阶段已经不再是第一性能瓶颈。

新的函数级耗时分析与分层诊断表明：

> **当前第一 LATIN 性能热点已经迁移到 PGD 富集，主要累计成本来自空间—时间固定点迭代中重复触发的 residual-LS/LSMR 空间求解。**

因此本阶段不再继续优化 Local，而是进一步回答：

> **为什么后期 PGD 模态需要大量空间—时间固定点迭代，以及当前固定点停止精度是否存在过度求解。**

本阶段正式定义为：

> **OPT-5：塔筒 PGD 空间—时间固定点容差跨周期校准**

最终判定：

> **PASS**

需要特别强调：

> **OPT-5 与 OPT-1 至 OPT-4 的性质不同。OPT-5 修改的是数值算法内部停止容差，因此属于经过跨周期验证的数值参数校准，而不是严格意义上的“数学完全等价实现优化”。**

---

# 2. OPT-4 后为什么必须重新定位性能热点

OPT-4 正式 10 周期 benchmark 为：

```text
iterations         = 39
trial evaluations  = 60
PGD rank           = 21
modes added        = 21
final xi           = 8.941607234831e-06
LATIN solver       = 31.701352 s
LATIN total        = 32.700045 s
```

OPT-4 相对于 OPT-3 已经把 10 周期 LATIN total 从：

$$ 162.864288\ \mathrm{s} $$

降低到：

$$ 32.700045\ \mathrm{s} $$

因此 OPT-4 之前的热点占比已经失效。

重新进行 10 周期函数级耗时分析后，LATIN solver 的累计耗时结构发生明显迁移。

代表性 cProfile 结果为：

```text
solve_tower_latin_pgd                       38.236 s
PGD enrichment                              21.678 s
raw fixed point                             20.401 s
residual-LS spatial solves                  15.395 s
SciPy LSMR                                  14.897 s
tower local stage                            6.993 s
full time-function updates                   4.030 s
temporal reduced solves                      3.565 s
search directions                            3.347 s
```

以该次 cProfile 中 LATIN solver 的 `38.236 s` 为参考，主要模块约占：

```text
PGD enrichment      ≈ 56.7%
residual-LS         ≈ 40.3%
LSMR                ≈ 39.0%
local stage         ≈ 18.3%
time update         ≈ 10.5%
search directions   ≈  8.8%
```

这里的累计时间存在嵌套关系，因此不能把各项百分比直接相加。

可以确认的是：

> **OPT-4 后 Local 已从第一瓶颈下降为次要热点，PGD enrichment 成为新的第一模块级性能热点。**

---

# 3. 为什么不能仅根据 cProfile 就直接优化 LSMR

函数级耗时分析显示：

```text
SciPy LSMR cumulative time ≈ 14.897 s
```

这意味着 LSMR 当前确实占据很大的累计时间。

但是仅凭这一结果不能判断：

> “单次 LSMR 求解本身已经成为根本瓶颈。”

原因在于 LSMR 位于 PGD 空间—时间固定点内部。

如果固定点外层需要重复很多轮，则即使每一次 LSMR 都正常、稳定地收敛，其累计时间仍然可能很高。

因此真正需要区分的是：

```text
可能性 A：
每次 LSMR 本身越来越难
        |
        -> 单次 Krylov 迭代次数升高
        -> 应优先考虑线性求解器、预条件、算子实现

可能性 B：
每次 LSMR 正常
        |
        -> 但是固定点外层重复次数很多
        -> 应优先研究空间—时间固定点收敛
```

因此 OPT-5 前首先进行了固定点—LSMR 分层诊断。

---

# 4. 固定点—LSMR 分层诊断

使用不修改生产源码的外部诊断脚本，对 10 周期正式 LATIN-PGD 路径进行内存级包装记录。

诊断保持：

```text
cycles             = 10
Nt                 = 401
Nq                 = 320
termination        = converged
LATIN iterations   = 39
trial evaluations  = 60
PGD rank           = 21
modes added        = 21
final xi           = 8.941607234831e-06
```

即数值路径与 OPT-4 稳定 checkpoint 保持一致。

21 次成功 PGD 富集合计得到：

```text
fixed-point iterations total = 523
residual-LS solves total      = 523
LSMR solves total             = 523
```

因此实测确认：

$$ N_{\mathrm{FP}} = N_{\mathrm{residual-LS}} = N_{\mathrm{LSMR}} = 523 $$

也就是说：

> **当前每一轮空间—时间固定点恰好触发一次 residual-LS / LSMR 空间求解。**

---

# 5. LSMR 本身是否出现失效或明显恶化

523 次 LSMR 求解的统计为：

```text
count                = 523
sum iterations       = 12348
minimum iterations   = 18
mean iterations      = 23.610
median iterations    = 23
maximum iterations   = 56
```

全部求解：

```text
converged = True
istop     = 2
```

当前 residual-LS 源码中，`istop = 2` 被明确视为成功达到请求的最小二乘精度。

因此：

> **523 次 LSMR 均正常结束，没有出现未收敛、失败或停止码异常。**

更重要的是，后期模态的平均 LSMR 迭代次数并没有随 PGD rank 明显上升。

例如第 21 次富集：

```text
fixed-point iterations = 95
mean LSMR iterations   ≈ 23.48
```

也就是说，第 21 个模态之所以昂贵，不是因为某一次 LSMR 突然需要几百次迭代，而是因为：

> **正常的约 20 至 24 次 LSMR Krylov 迭代被固定点外层重复调用了 95 次。**

因此当前第一因果层问题应从：

```text
“LSMR 为什么慢”
```

进一步前移到：

```text
“PGD 空间—时间固定点为什么需要重复这么多轮”
```

---

# 6. 21 个成功富集的固定点次数

10 周期下，21 次成功富集的固定点迭代次数具有明显波动。

代表性统计为：

```text
number of enrichments = 21
fixed-point sum        = 523
minimum               = 4
mean                  = 24.905
median                = 21
maximum               = 95
```

后期若干模态尤其昂贵：

```text
enrichment 15   -> 45 fixed-point iterations
enrichment 18   -> 54 fixed-point iterations
enrichment 21   -> 95 fixed-point iterations
```

第 21 个模态单独占：

$$ \frac{95}{523}\approx 18.16\% $$

全部固定点迭代。

最后 7 个模态合计约占全部固定点工作的 53%。

因此高成本并不是均匀分布在所有模态上，而是明显向后期模态集中。

---

# 7. 为什么还不能直接引入 Anderson 或 Aitken 加速

观察到 45、54、95 次固定点迭代后，可以考虑：

- Aitken 加速；
- Anderson 加速；
- 松弛因子调整；
- 外推；
- 其他非线性固定点加速。

但是在正式引入这些方法之前，必须先判断当前固定点历史属于哪一种类型：

```text
A. 稳定慢收敛
B. 初期快、后期长尾
C. 持续振荡
D. 平台或停滞
```

不同形态对应不同策略。

因此下一步不是立即修改生产代码，而是展开最昂贵模态的完整固定点历史。

---

# 8. 第 15 个模态的固定点历史

第 15 次富集：

```text
LATIN iteration       = 18
rank before           = 14
fixed-point iterations= 45
initial chi           = 2.436649739345e-01
final chi             = 9.351773854878e-07
number of increases   = 0
```

全部 45 轮严格单调下降。

前期收缩较快，例如：

```text
chi_2 / chi_1 ≈ 0.5273
chi_3 / chi_2 ≈ 0.5215
chi_4 / chi_3 ≈ 0.5160
```

但是后期逐渐稳定到：

```text
chi_(k+1) / chi_k ≈ 0.8384
```

最后 20 轮几何平均收缩比：

$$ \rho_{15,\mathrm{tail}}\approx 0.83842 $$

因此第 15 模态属于：

> **完全单调，但后期线性慢收敛。**

---

# 9. 第 18 个模态的固定点历史

第 18 次富集：

```text
LATIN iteration        = 21
rank before            = 17
fixed-point iterations = 54
initial chi            = 1.166380346083e-01
final chi              = 8.703017449475e-07
number of increases    = 4
```

前 3 至 7 轮存在短暂反弹：

```text
chi_4 / chi_3 ≈ 1.2639
chi_5 / chi_4 ≈ 1.3231
chi_6 / chi_5 ≈ 1.2611
chi_7 / chi_6 ≈ 1.1286
```

但是第 8 轮以后进入稳定下降。

后期收缩比逐渐趋于：

```text
≈ 0.7746
```

最后 20 轮几何平均收缩比：

$$ \rho_{18,\mathrm{tail}}\approx 0.77460 $$

因此第 18 模态不是持续振荡。

更准确的描述为：

> **短暂初始非单调之后进入稳定线性收敛。**

---

# 10. 第 21 个模态的固定点历史

第 21 次富集是当前最昂贵的单个模态：

```text
LATIN iteration        = 24
rank before            = 20
fixed-point iterations = 95
initial chi            = 1.793256835592e-01
final chi              = 9.270810291572e-07
number of increases    = 4
```

第 6 至 9 轮存在短暂上升：

```text
chi_6 / chi_5 ≈ 1.2455
chi_7 / chi_6 ≈ 1.3596
chi_8 / chi_7 ≈ 1.3209
chi_9 / chi_8 ≈ 1.1437
```

但是之后连续下降，没有再次振荡。

后期收缩比逐渐稳定到：

```text
≈ 0.86986
```

最后 20 轮几何平均收缩比：

$$ \rho_{21,\mathrm{tail}}\approx 0.86987 $$

这意味着固定点误差后期每进行一次迭代，只缩小为原来的约 87%。

如果用线性收敛近似：

$$ \chi_{k+1}\approx \rho\chi_k $$

则降低一个十进数量级大约需要：

$$ N_{\mathrm{decade}}\approx\frac{\log(0.1)}{\log(0.86987)}\approx 16.5 $$

轮固定点迭代。

因此第 21 模态的 95 次固定点迭代并不来自数值发散，而是：

> **固定点映射的收缩因子已经较接近 1，从而产生稳定但缓慢的线性收敛。**

---

# 11. 固定点历史带来的第一项关键认识

第 15、18、21 个高成本模态没有表现出持续振荡。

因此当前不能简单判断：

> “需要更强阻尼。”

因为对于已经稳定单调下降的尾部，进一步减小更新幅度可能反而降低收敛速度。

同样，也没有证据显示固定点进入机器精度停滞。

`chi` 一直在持续下降并最终穿过 `1e-6`。

因此当前主要问题是：

$$ \boxed{\mathrm{slow\ linear\ contraction}} $$

即：

> **稳定的线性慢收敛。**

---

# 12. 为什么优先检查固定点容差是否过严

当前塔筒顶层 LATIN 收敛容差为：

$$ \varepsilon_{\mathrm{LATIN}} = 10^{-5} $$

OPT-4 时 PGD 富集内部固定点容差为：

$$ \varepsilon_{\mathrm{FP}} = 10^{-6} $$

因此内部固定点比外层 LATIN 收敛判据严格一个数量级。

这并不自动意味着 `1e-6` 错误。

但是在当前后期模态固定点收缩比接近 0.8 至 0.87 的情况下，从 `1e-5` 继续压到 `1e-6` 会额外消耗多轮 residual-LS / LSMR 空间求解。

因此在引入新的固定点加速算法之前，更低风险的问题是：

> **当前 `1e-6` 是否对塔筒 benchmark 存在内部过度求解。**

---

# 13. 固定点容差敏感性试验

在不修改生产源码的条件下，对 10 周期 benchmark 比较：

```text
1.0e-6   baseline
3.0e-6   intermediate
1.0e-5   candidate
```

结果为：

| 固定点容差 | LATIN 迭代 | Trial | PGD 阶数 | 模态数 | 固定点总次数 | 最终 xi |
|---:|---:|---:|---:|---:|---:|---:|
| `1e-6` | 39 | 60 | 21 | 21 | 523 | `8.941607234831e-06` |
| `3e-6` | 39 | 60 | 21 | 21 | 480 | `8.941608485939e-06` |
| `1e-5` | 39 | 60 | 21 | 21 | 431 | `8.941609245924e-06` |

因此：

```text
1e-6 -> 3e-6:
523 -> 480
saved = 43
reduction = 8.22%

1e-6 -> 1e-5:
523 -> 431
saved = 92
reduction = 17.59%
```

同时 `1e-5` 与 `1e-6` 的最终 LATIN 指标差：

$$ |\Delta \xi| = 2.011093169167\times 10^{-12} $$

离散求解路径保持一致。

---

# 14. 为什么“离散路径一致”仍然不够

仅仅得到：

```text
39 iterations
60 trial evaluations
rank 21
21 modes
same final xi
```

仍然不足以证明两个容差产生了近似相同的数值解。

因为不同 PGD 模态可能通过不同低秩分解组合得到相近的全局指标。

因此必须进一步比较：

- 最终接受 LATIN 状态；
- 最终 Local 状态；
- 全部 13 个材料场；
- PGD 重构修正场；
- LATIN 指标历史；
- Trial A/B 路径；
- PGD 基底阶数历史。

---

# 15. 10 周期最终接受状态全场等价性

比较：

```text
baseline  fixed_point_tolerance = 1e-6
candidate fixed_point_tolerance = 1e-5
```

最终接受 LATIN 状态包含 13 个材料场。

最差相对 $L_2$ 误差为：

```text
plastic_strain_rate
relative L2 = 3.134427e-08
```

最差相对最大幅值误差为：

```text
plastic_strain_rate
relative max-scale = 3.959759e-08
```

因此最终接受状态整体差异处于约：

$$ 10^{-8} $$

量级。

---

# 16. 10 周期关键物理场差异

## 16.1 应力

基线最大绝对应力约为：

$$ |\sigma|_{\max}\approx 110.8\ \mathrm{MPa} $$

两种固定点容差之间最大绝对应力差：

$$ \max|\Delta \sigma| = 9.068288\times10^{-7}\ \mathrm{MPa} $$

相对最大幅值误差：

$$ 8.184385\times10^{-9} $$

## 16.2 损伤

基线最大损伤：

$$ D_{\max}\approx 8.57511202\times10^{-3} $$

两种容差最大损伤绝对差：

$$ \max|\Delta D| = 1.576024\times10^{-10} $$

相对最大幅值误差：

$$ 1.837905\times10^{-8} $$

## 16.3 塑性应变

基线最大塑性应变幅值约：

$$ |\varepsilon^p|_{\max}\approx 4.99\times10^{-4} $$

最大绝对差约：

$$ 9.79\times10^{-12} $$

因此从最终物理状态看，`1e-5` 对当前 10 周期结果产生的扰动极小。

---

# 17. 最终 Local 状态全场等价性

最终 Local 状态 13 个场同样进行比较。

最差相对 $L_2$ 误差：

```text
damage_rate
1.444759e-08
```

最差相对最大幅值误差：

```text
R_bar
2.506794e-08
```

因此 Local 状态与最终接受状态具有相同数量级的一致性。

---

# 18. PGD 重构修正场为什么相对误差更大

PGD 重构修正场的比较结果为：

```text
plastic_strain_correction
relative L2 ≈ 7.675037e-06

plastic_strain_rate_correction
relative L2 ≈ 6.115049e-06

stress_correction
relative L2 ≈ 7.242956e-06
```

这些相对误差高于最终物理状态的 `1e-8` 量级。

但是必须结合修正场自身幅值解释。

例如应力修正场最大幅值仅约：

$$ 5.085519\times10^{-3}\ \mathrm{MPa} $$

两种容差之间最大绝对差只有：

$$ 2.776966\times10^{-8}\ \mathrm{MPa} $$

因此当前现象更合理的解释是：

> **低幅值 PGD 修正场对基底细微变化的相对误差更敏感，但最终重构后的物理场保持高度一致。**

不能因为 PGD 修正场的相对误差达到 `1e-5`，就直接判断最终塔筒响应发生了同量级物理误差。

---

# 19. 求解历史是否保持一致

10 周期两种容差下：

```text
trial_basis_size_history  exactly equal
modes_added_history       exactly equal
trial_kind_history        exactly equal
commit_kind_history       exactly equal
```

连续指标历史存在极小数值差异，例如：

```text
indicator_history relative L2
≈ 4.55e-06
```

但是离散 Trial A/B 事务路径完全一致。

因此可以确认：

> **`1e-5` 没有改变当前 10 周期 benchmark 的 LATIN 事务路径与 PGD rank 演化路径。**

---

# 20. 为什么必须进一步做 1 / 2 / 5 / 10 周期跨周期验证

如果只验证 10 周期，则仍存在一种可能：

> `1e-5` 只是对这一特定长时间网格恰好安全。

当前项目已经建立稳定的 1 / 2 / 5 / 10 周期 benchmark 家族：

```text
1 cycle   -> Nt = 41
2 cycles  -> Nt = 81
5 cycles  -> Nt = 201
10 cycles -> Nt = 401
```

因此正式校准必须跨越全部四个已验证时间尺度。

---

# 21. 1 周期跨周期校准结果

1 周期：

```text
same discrete path = True

LATIN iterations:
18 -> 18

trial evaluations:
29 -> 29

PGD rank:
11 -> 11

fixed-point sum:
173 -> 142

saved:
31 = 17.92%

final xi:
7.918424536257e-06
->
7.918467517050e-06
```

最终指标绝对差：

$$ 4.298\times10^{-11} $$

最终接受状态最差相对 $L_2$：

$$ 7.268\times10^{-9} $$

PGD 修正场最差相对 $L_2$：

$$ 3.527\times10^{-5} $$

---

# 22. 2 周期跨周期校准结果

2 周期：

```text
same discrete path = True

LATIN iterations:
23 -> 23

trial evaluations:
36 -> 36

PGD rank:
13 -> 13

fixed-point sum:
229 -> 185

saved:
44 = 19.21%
```

最终指标绝对差：

$$ 3.386\times10^{-12} $$

最终接受状态最差相对 $L_2$：

$$ 4.223\times10^{-8} $$

最终接受状态最差相对最大幅值误差：

$$ 6.871\times10^{-8} $$

PGD 修正场最差相对 $L_2$：

$$ 1.281\times10^{-4} $$

---

# 23. 5 周期跨周期校准结果

5 周期：

```text
same discrete path = True

LATIN iterations:
33 -> 33

trial evaluations:
50 -> 50

PGD rank:
17 -> 17

fixed-point sum:
391 -> 327

saved:
64 = 16.37%
```

最终指标绝对差：

$$ 5.014\times10^{-11} $$

最终接受状态最差相对 $L_2$：

$$ 3.153\times10^{-8} $$

PGD 修正场最差相对 $L_2$：

$$ 1.158\times10^{-4} $$

---

# 24. 10 周期跨周期校准结果

10 周期：

```text
same discrete path = True

LATIN iterations:
39 -> 39

trial evaluations:
60 -> 60

PGD rank:
21 -> 21

fixed-point sum:
523 -> 431

saved:
92 = 17.59%
```

最终指标绝对差：

$$ 2.011\times10^{-12} $$

最终接受状态最差相对 $L_2$：

$$ 3.134\times10^{-8} $$

PGD 修正场最差相对 $L_2$：

$$ 7.675\times10^{-6} $$

---

# 25. 跨周期总结果

四组 benchmark 汇总：

| 周期 | 路径一致 | LATIN 迭代 | PGD 阶数 | 固定点次数 `1e-6 -> 1e-5` | 固定点减少 | 最差最终状态相对 $L_2$ |
|---:|:---:|---:|---:|---:|---:|---:|
| 1 | True | `18 -> 18` | `11 -> 11` | `173 -> 142` | 17.92% | `7.268e-09` |
| 2 | True | `23 -> 23` | `13 -> 13` | `229 -> 185` | 19.21% | `4.223e-08` |
| 5 | True | `33 -> 33` | `17 -> 17` | `391 -> 327` | 16.37% | `3.153e-08` |
| 10 | True | `39 -> 39` | `21 -> 21` | `523 -> 431` | 17.59% | `3.134e-08` |

总固定点次数：

$$ 1316\longrightarrow1085 $$

减少：

$$ 231 $$

相对减少：

$$ 17.55\% $$

因为每一轮固定点对应一次 residual-LS / LSMR 空间求解，所以在这四组 benchmark 中也对应减少：

$$ 231 $$

次 residual-LS / LSMR 调用。

---

# 26. 跨周期验证最重要的结论

四个 benchmark 全部满足：

```text
same discrete path = True
```

最大最终接受状态相对 $L_2$ 误差：

$$ 4.223424\times10^{-8} $$

最大最终接受状态相对最大幅值误差：

$$ 6.871266\times10^{-8} $$

最大 PGD 修正场相对 $L_2$ 误差：

$$ 1.281216\times10^{-4} $$

因此可以确认：

> **在当前已验证的 1 / 2 / 5 / 10 周期 matched asymmetric tower benchmark 家族中，将固定点容差从 `1e-6` 校准至 `1e-5`，不会改变已验证的离散 LATIN-PGD 求解路径，并且最终物理状态扰动保持在约 `1e-8` 相对量级。**

这为生产参数校准提供了充分数值证据。

---

# 27. 为什么不是修改所有 PGD 模块的默认值

仓库中 `fixed_point_tolerance` 不只存在于塔筒求解器中。

还存在于：

```text
latin/pgd_enrichment.py
latin/pgd_global_stage.py
latin/pgd_solver.py
latin/tower_pgd_enrichment.py
```

以及多个示例与测试。

当前 1 / 2 / 5 / 10 周期验证严格针对的是：

> **塔筒 LATIN-PGD benchmark。**

因此不能把当前结果泛化为：

> **所有 PGD 问题都应该使用 `1e-5`。**

OPT-5 最终采用最小作用域设计。

---

# 28. OPT-5 正式代码修改

正式代码只修改：

```text
latin/tower_latin_pgd_solver.py
```

将塔筒顶层：

```python
fixed_point_tolerance: float = 1.0e-6
```

校准为：

```python
fixed_point_tolerance: float = 1.0e-5
```

而以下底层与通用模块默认值保持不变：

```text
latin/pgd_enrichment.py
latin/pgd_global_stage.py
latin/pgd_solver.py
latin/tower_pgd_enrichment.py
```

因此 OPT-5 的作用域是：

> **塔筒顶层求解器默认参数。**

不是：

> **全仓库统一降低固定点精度。**

---

# 29. 为什么保留底层 `tower_pgd_enrichment` 的保守默认值

`enrich_tower_pgd_basis_once()` 是比当前 benchmark 更底层的功能接口。

它可以被：

- 独立单元测试调用；
- integration test 调用；
- 其他实验脚本调用；
- 未来不同塔筒问题直接调用。

当前跨周期证据支持的是：

```text
solve_tower_latin_pgd()
```

这一完整塔筒算法链的参数校准。

因此保留底层富集函数的保守默认值可以避免把一个已在特定完整流程中验证的经验值错误泛化到其他调用场景。

这是 OPT-5 范围控制的重要部分。

---

# 30. 新增回归测试保护

正式测试文件：

```text
tests/test_tower_latin_pgd_solver.py
```

在现有 Trial-B 富集事务测试中捕获：

```python
enrich_tower_pgd_basis_once
```

的 mock 调用。

新增断言验证：

```python
enrich_mock.call_args.kwargs["fixed_point_tolerance"]
```

等于：

```text
1.0e-5
```

这条测试保护的不是固定点算法本身，而是：

> **塔筒顶层经过校准的默认参数必须被正确传递到实际 PGD enrichment。**

---

# 31. 针对性单元测试

运行：

```powershell
python -m pytest .\tests\test_tower_latin_pgd_solver.py
```

结果：

```text
8 passed in 0.24 s
```

因此塔筒顶层 Trial A/B 事务逻辑、回滚逻辑、富集调用以及新增参数传递断言均通过。

---

# 32. 全仓库回归测试

运行：

```powershell
python -m pytest
```

结果：

```text
314 passed, 10 warnings in 101.06 s
```

测试总数与 OPT-4 稳定 checkpoint 一致。

10 条 warning 均来自已有 Matplotlib / distutils 弃用提示，没有出现新的测试失败或回归。

因此仓库级回归：

> **PASS**

---

# 33. 正式 10 周期生产路径复验

完成代码修改后，重新使用未显式覆盖 `fixed_point_tolerance` 的正式效率 pilot：

```powershell
python .\tower_asymmetric_efficiency_scaling_pilot.py --cycles 10
```

这意味着本次运行实际使用的正是新的塔筒顶层默认值：

```text
fixed_point_tolerance = 1e-5
```

正式结果：

```text
termination_reason = converged
converged          = True
iterations         = 39
attempted          = 39
trial evaluations  = 60
PGD rank           = 21
modes added        = 21
final xi           = 8.941609245924e-06
```

因此生产代码路径与前面的外部容差敏感性诊断完全一致。

正式 10 周期数值路径：

> **PASS**

---

# 34. OPT-5 正式 10 周期 wall-time

本次生产 benchmark：

```text
FOM wall time          = 295.350239 s
LATIN setup            =   1.015139 s
LATIN solver           =  28.195112 s
LATIN total            =  29.210251 s
sample FOM/LATIN ratio =  10.111185
```

相对于 OPT-4 正式 benchmark：

```text
OPT-4 LATIN solver = 31.701352 s
OPT-5 LATIN solver = 28.195112 s
```

solver 时间降低：

$$ \frac{31.701352-28.195112}{31.701352}\times100\%\approx11.06\% $$

LATIN total：

```text
OPT-4 LATIN total = 32.700045 s
OPT-5 LATIN total = 29.210251 s
```

total 时间降低：

$$ \frac{32.700045-29.210251}{32.700045}\times100\%\approx10.67\% $$

---

# 35. 为什么固定点次数减少 17.55%，但 wall-time 只减少约 10%

固定点减少比例与总时间减少比例不同是合理的。

OPT-5 只减少空间—时间 fixed-point 内部的重复工作。

LATIN solver 中仍然包含：

- Local stage；
- search directions；
- Trial-A time update；
- Trial-B construction；
- complete trial evaluation；
- 全时域 PGD time-function reoptimization；
- 状态复制与事务逻辑；
- 其他 Python 与 NumPy 运算。

因此如果固定点部分减少约 17.5%，总求解时间只降低约 10%，符合整体成本组成。

---

# 36. 为什么 wall-time 仍然不能作为唯一验证依据

不同运行中的 FOM 与 LATIN wall-time 存在系统调度、缓存、CPU 频率等波动。

例如 OPT-4 正式 FOM：

```text
293.417296 s
```

OPT-5 正式 FOM：

```text
295.350239 s
```

两次 FOM 算法没有发生对应改变，但仍存在约秒级差异。

因此 OPT-5 最稳定的算法层证据不是：

> “某一次 wall-time 恰好减少多少秒。”

而是：

```text
1/2/5/10-cycle discrete path unchanged
1316 -> 1085 fixed-point iterations
231 residual-LS/LSMR calls removed
full accepted state error ~1e-8
```

wall-time 只作为进一步支持。

---

# 37. 当前 FOM/LATIN 比值应如何解释

本次实测：

$$ \frac{T_{\mathrm{FOM}}}{T_{\mathrm{LATIN}}}=10.111185 $$

这是真实的当前实现结果。

但是仍然只能解释为：

> **当前 FOM Python 实现与当前 LATIN-PGD Python 实现之间的实测时间比。**

不能解释为：

> **LATIN-PGD 方法理论上必然比 FOM 快约 10.1 倍。**

原因仍然包括：

- FOM 自身存在尚未统一优化的 Python 标量开销；
- 两套实现的代码路径、对象管理与向量化程度不同；
- 当前问题只有 10 个梁柱单元和 320 个材料点；
- 当前 benchmark 是特定不对称循环荷载与特定材料参数。

方法级公平效率比较仍需单独建立统一实现基线。

---

# 38. OPT-1 至 OPT-5 的 10 周期时间演化

10 周期 LATIN total 的主要性能演化为：

```text
original  ≈ 467.362 s
OPT-1     ≈ 223.161 s
OPT-2     ≈ 174.718 s
OPT-3     ≈ 162.864 s
OPT-4     =  32.700 s
OPT-5     =  29.210 s
```

从原始实现到 OPT-5：

$$ 467.362\ \mathrm{s}\longrightarrow29.210251\ \mathrm{s} $$

累计时间降低约：

$$ 93.75\% $$

累计时间比约：

$$ \frac{467.362}{29.210251}\approx16.00 $$

但是这一累计结果必须分成两类解释。

OPT-1 至 OPT-4：

> **主要是保持数学算法不变的实现级优化。**

OPT-5：

> **是经过严格数值验证的内部固定点停止容差校准。**

因此不能再把 OPT-1 至 OPT-5 整体描述为“完全不改变任何数值参数的纯实现加速”。

---

# 39. OPT-5 带来的第一项科研认识：性能热点会随实现优化发生迁移

原始实现中：

> Local stage 占据绝对主导。

OPT-4 之后：

> PGD enrichment 成为第一模块级热点。

这说明：

$$ \boxed{\mathrm{hotspot\ is\ implementation\ state\ dependent}} $$

也就是说：

> **性能瓶颈不是算法名称上的固定标签，而是会随着前一阶段优化而迁移。**

如果继续依据 OPT-1 前的 profiling 优化 Local，就会错过当前真正的主要成本来源。

---

# 40. OPT-5 带来的第二项科研认识：累计 LSMR 成本不等于单次 LSMR 困难

OPT-4 后 profiling 中 LSMR 累计耗时很高。

但是分层诊断发现：

```text
523 LSMR solves
all converged
all istop = 2
mean iterations ≈ 23.61
```

而固定点次数最高达到：

```text
95
```

因此：

> **当前 LSMR 高累计成本主要由固定点外层重复调用放大，而不是由后期单次 LSMR 迭代数恶化主导。**

这一区分对后续优化方向非常重要。

---

# 41. OPT-5 带来的第三项科研认识：后期 PGD 模态表现为线性慢收敛

最昂贵的几个模态后期固定点误差近似满足：

$$ \chi_{k+1}\approx\rho\chi_k $$

其中：

```text
mode 15 tail rho ≈ 0.8384
mode 18 tail rho ≈ 0.7746
mode 21 tail rho ≈ 0.8699
```

因此当前后期固定点具有清晰的线性收敛特征。

这意味着未来如果继续优化固定点，可优先研究：

- Aitken 型动态松弛；
- Anderson acceleration；
- 其他针对线性慢收敛的外推方式。

但是 OPT-5 本身没有引入这些方法。

---

# 42. OPT-5 带来的第四项科研认识：内部子问题容差不能脱离外层目标精度设置

当前外层：

$$ \varepsilon_{\mathrm{LATIN}}=10^{-5} $$

原内部固定点：

$$ \varepsilon_{\mathrm{FP}}=10^{-6} $$

跨周期验证表明，对于当前 benchmark，把内部容差校准为：

$$ \varepsilon_{\mathrm{FP}}=10^{-5} $$

仍然保持完整物理状态近似不变。

因此可以得到一个更一般但仍需谨慎表述的认识：

> **嵌套算法中的内部子问题不一定需要比外层目标精度严格一个固定数量级；其合理容差应通过对外层收敛路径和最终物理场的敏感性验证确定。**

不能把这一结论直接推广为所有 LATIN-PGD 问题都应该令两个容差完全相等。

---

# 43. 为什么 OPT-5 不是降低科学精度换速度

如果只看到：

```text
1e-6 -> 1e-5
```

很容易误解为：

> “通过降低精度换取速度。”

但本阶段实际做法是：

```text
先诊断固定点收敛形态
        |
再做 10-cycle tolerance sweep
        |
再做 13-field full-state equivalence
        |
再做 1/2/5/10-cycle cross-validation
        |
再修改生产默认值
        |
再跑 targeted tests
        |
再跑 314-test full regression
        |
再跑 production 10-cycle benchmark
```

因此更准确的表述是：

> **OPT-5 是通过数值敏感性与全场等价性证据识别并去除当前 benchmark 中不必要的内部过度求解。**

---

# 44. 当前已经证明的内容

截至 OPT-5，可以正式确认：

1. OPT-4 后第一 LATIN 模块级性能热点已经迁移到 PGD enrichment。
2. 10 周期 21 次成功富集合计执行 523 轮空间—时间固定点。
3. 每一轮固定点对应一次 residual-LS / LSMR 空间求解。
4. 523 次 LSMR 全部正常收敛，停止码均为 `istop = 2`。
5. LSMR 平均约 23.61 次 Krylov 迭代。
6. 后期昂贵模态的主要成本来自固定点重复次数，而不是单次 LSMR 明显恶化。
7. 第 15、18、21 个模态后期表现为稳定线性慢收敛。
8. 10 周期下 `1e-5` 固定点容差保持 `39 / 60 / 21` 求解路径。
9. 10 周期最终接受状态最差相对 $L_2$ 误差约为 `3.13e-8`。
10. 1 / 2 / 5 / 10 周期全部保持离散求解路径一致。
11. 四组 benchmark 最终接受状态最大相对 $L_2$ 误差约为 `4.22e-8`。
12. 四组 benchmark 固定点总次数从 1316 降到 1085。
13. 累计减少 231 轮固定点，即减少约 17.55%。
14. 正式塔筒顶层默认容差已从 `1e-6` 校准为 `1e-5`。
15. 通用 PGD 和底层 tower enrichment 默认值保持不变。
16. 塔筒 solver 针对性测试 `8 passed`。
17. 全仓库回归 `314 passed, 10 warnings`。
18. 修改后正式 10 周期生产路径仍为 `39 / 60 / rank 21`。
19. 修改后正式 10 周期 LATIN total 为 `29.210251 s`。
20. OPT-5 相对 OPT-4 正式 10 周期 total 减少约 10.67%。

---

# 45. 当前仍然没有证明的内容

以下结论仍然不能宣称：

1. 不能宣称所有塔筒问题的最佳固定点容差都是 `1e-5`。
2. 不能宣称所有 PGD 问题都应该令内部固定点容差等于外层 LATIN 容差。
3. 不能把当前 `1e-5` 直接推广到不同材料模型、不同荷载比、不同塔筒离散或不同损伤水平。
4. 不能宣称 `1e-6` 在数学上错误；只能说它在当前 benchmark 中表现出可验证的过度求解。
5. 不能宣称 PGD 修正模态本身逐模态完全一致；当前验证重点是完整重构场与求解路径。
6. 不能宣称 LSMR 已经无需进一步优化。
7. 不能宣称固定点加速算法已经无必要。
8. 不能宣称当前约 10.11 倍 FOM/LATIN 比值是方法理论加速比。
9. 不能宣称 10 周期以上更长时间尺度仍然保持完全相同的固定点节约比例。
10. 不能宣称当前 `1e-5` 是全局最优容差；本阶段只证明其为经过验证的安全校准点。

---

# 46. OPT-5 的正式验证证据链

当前证据链可以整理为：

```text
Level 1
OPT-4 post-profile hotspot migration
PGD enrichment becomes first module-level hotspot

Level 2
fixed-point / LSMR layered diagnostic
523 FP iterations
523 residual-LS / LSMR solves
all LSMR converged
mean LSMR iterations ≈ 23.61

Level 3
selected fixed-point history diagnostic
mode 15: 45 iterations
mode 18: 54 iterations
mode 21: 95 iterations
late-stage behavior = stable linear slow convergence

Level 4
10-cycle tolerance sensitivity
1e-6 -> 3e-6 -> 1e-5
same 39 / 60 / rank 21 path
523 -> 480 -> 431 fixed-point iterations

Level 5
10-cycle full-state equivalence
worst accepted-state relative L2 ≈ 3.13e-8

Level 6
1/2/5/10-cycle cross-validation
all same discrete path
aggregate FP: 1316 -> 1085
saved 17.55%
max accepted-state relative L2 ≈ 4.22e-8

Level 7
targeted solver tests
8 passed

Level 8
full repository regression
314 passed, 10 warnings

Level 9
formal 10-cycle production benchmark
39 iterations
60 trial evaluations
rank 21
final xi = 8.941609245924e-06
LATIN total = 29.210251 s
```

因此 OPT-5 的 PASS 建立在性能诊断、固定点收敛形态、容差敏感性、全场等价性、跨周期验证、单元测试、全仓库回归和正式生产 benchmark 的完整证据链上。

---

# 47. OPT-5 正式代码 checkpoint

代码提交：

```text
8224819
perf: calibrate tower PGD fixed-point tolerance
```

提交内容只包括：

```text
latin/tower_latin_pgd_solver.py
tests/test_tower_latin_pgd_solver.py
```

提交统计：

```text
2 files changed
6 insertions
2 deletions
```

未跟踪的：

- `.pstats` profiling 文件；
- 外部诊断脚本；
- 容差敏感性脚本；
- 全场等价性脚本；
- 跨周期验证脚本；
- 临时 patch 文件；

均没有进入正式代码提交。

---

# 48. 当前优化链

当前性能优化分支主要稳定节点：

```text
900ddf2  docs: summarize 10-cycle LATIN function-level profiling

c27af15  perf: optimize scalar damage clipping in local stage
55d578e  docs: summarize scalar clipping optimization

e2d15b9  perf: cache material derived constants
2f0c15a  docs: summarize derived constant caching optimization

2a88b68  perf: flatten RK4 local-stage hot path
f9507b5  docs: summarize RK4 hot-path optimization

423a9c4  perf: vectorize homogeneous tower local stage
91a94fa  docs: summarize homogeneous tower local-stage vectorization

8224819  perf: calibrate tower PGD fixed-point tolerance
```

OPT-5 文档提交将在本总结人工检查通过后单独完成。

---

# 49. OPT-5 与原论文的关系应如何表述

当前固定点容差 `1e-5` 是：

> **当前塔筒 implementation 经跨周期数值验证得到的工程数值参数选择。**

不能表述为：

> “Bhattacharyya 等原论文明确规定固定点容差应为 `1e-5`。”

本阶段没有新的证据证明原论文规定了这一数值。

因此在论文、汇报和文档中应严格区分：

```text
原论文理论框架
        !=
当前塔筒 implementation 的数值容差校准
```

OPT-5 属于后者。

---

# 50. 对导师汇报时的推荐表述

可以用以下逻辑解释：

> OPT-4 以后，我们重新做了函数级耗时分析，发现性能瓶颈已经从局部材料更新迁移到了 PGD 富集。进一步把富集拆开以后发现，10 周期中 21 次成功富集合计需要 523 轮空间—时间固定点，每轮都要做一次 residual-LS/LSMR 空间求解。523 次 LSMR 本身全部正常收敛，平均只需要约 24 次 Krylov 迭代，真正的问题是后期 PGD 模态固定点次数很多，第 21 个模态达到 95 次。展开固定点历史以后发现它不是持续振荡，而是稳定的线性慢收敛。由于当前外层 LATIN 容差是 `1e-5`，而内部固定点原来要求 `1e-6`，我们进一步验证内部是否存在过度求解。最终在 1、2、5、10 周期全部算例中，将内部固定点容差校准到 `1e-5` 后，LATIN 迭代次数、Trial 路径和 PGD 阶数全部不变，完整状态场相对误差仍只有约 `1e-8`，但固定点总次数减少约 17.5%。因此我们把塔筒顶层固定点容差正式校准为 `1e-5`，10 周期 LATIN 总时间进一步降到约 29.2 秒。

这一表述能够同时说明：

- 为什么优化对象发生变化；
- 为什么没有直接优化 LSMR；
- 为什么没有直接使用复杂固定点加速；
- 为什么 `1e-5` 具有数值依据；
- 为什么它仍然属于当前塔筒 implementation 的工程选择。

---

# 51. 下一阶段不应立即做什么

OPT-5 完成后，不建议马上：

- 直接把 LSMR tolerance 从 `1e-12` 大幅放宽；
- 直接引入 Anderson acceleration；
- 直接修改原有 Trial A/B 事务；
- 直接删除固定点循环；
- 直接把 `1e-5` 推广到所有 PGD 模块；
- 直接宣称当前求解器已经达到性能最优。

这些都需要新的独立证据。

---

# 52. 下一阶段推荐的科研问题

OPT-5 后，下一阶段最值得回答的问题有两个层次。

第一层：

> **在已经校准停止容差之后，PGD 后期模态的固定点线性收缩因子是否仍然是主要剩余成本。**

需要重新 profiling 当前 `1e-5` 正式代码。

第二层：

> **如果固定点仍然是首要热点，是否能够采用具有严格回滚与等价性保护的固定点加速方法。**

优先候选可以包括：

- Aitken 动态松弛；
- Anderson acceleration；
- 保守的自适应松弛。

但是这些属于新的算法修改，必须单独建立：

```text
baseline
    |
candidate acceleration
    |
fixed-point history comparison
    |
1/2/5/10-cycle path comparison
    |
full-state equivalence
    |
wall-time benchmark
```

的验证链。

---

# 53. 是否应该继续把下一阶段叫 OPT-6

在 OPT-5 文档和代码归档完成后，如果重新 profiling 仍然确认固定点为第一热点，则可以进入：

> **OPT-6：PGD 空间—时间固定点加速候选**

但在重新 profiling 前，不应预先假设 OPT-6 必然应该优化固定点。

因为 OPT-5 已经减少约 17.55% 的固定点工作，新的热点比例可能再次变化。

因此仍然坚持：

> **每轮优化后重新测量，再决定下一优化对象。**

---

# 54. 正式阶段结论

OPT-4 后的函数级耗时分析显示，塔筒 LATIN-PGD 的性能热点已经从 Local stage 迁移到 PGD enrichment。

进一步分层诊断确认：

```text
10 cycles
21 successful enrichments
523 fixed-point iterations
523 residual-LS / LSMR solves
all LSMR converged
mean LSMR iterations ≈ 23.61
```

昂贵后期模态的固定点历史表明，其主要特征不是持续振荡或数值停滞，而是稳定的线性慢收敛。

当前外层 LATIN 收敛容差为：

$$ 10^{-5} $$

原内部固定点容差为：

$$ 10^{-6} $$

经过 10 周期容差敏感性、完整 13 场等价性以及 1 / 2 / 5 / 10 周期跨周期验证后，确认将塔筒顶层固定点容差校准为：

$$ 10^{-5} $$

时：

```text
all 1/2/5/10-cycle discrete paths unchanged
aggregate fixed-point iterations:
1316 -> 1085
saved:
231
reduction:
17.55%

maximum accepted-state relative L2:
4.223424e-08
```

正式代码仅修改塔筒顶层默认值，不修改通用 PGD 与底层 enrichment 默认参数。

验证结果：

```text
targeted tower solver tests:
8 passed

full repository regression:
314 passed, 10 warnings
```

修改后正式 10 周期 benchmark：

```text
termination_reason = converged
iterations         = 39
trial evaluations  = 60
PGD rank           = 21
modes added        = 21
final xi           = 8.941609245924e-06

FOM                = 295.350239 s
LATIN setup        =   1.015139 s
LATIN solver       =  28.195112 s
LATIN total        =  29.210251 s
FOM/LATIN          =  10.111185
```

相对 OPT-4：

```text
LATIN solver:
31.701352 -> 28.195112 s
reduction ≈ 11.06%

LATIN total:
32.700045 -> 29.210251 s
reduction ≈ 10.67%
```

因此正式判定：

> **OPT-5：塔筒 PGD 空间—时间固定点容差跨周期校准 — PASS**

当前稳定代码 checkpoint：

```text
8224819
perf: calibrate tower PGD fixed-point tolerance
```

下一步应先完成本阶段文档归档，再重新进行 OPT-5 后的函数级性能分析，以决定是否需要进入新的固定点加速或其他性能优化阶段。
