# Tower LATIN-PGD 第 4 个 PGD mode 固定点周期振荡诊断阶段总结

**日期：2026-08-19**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**仓库：`Jiarui-Shi-DUT/Offshore_Wind_Turbine_LATIN_PGD`**  
**分支：`feature/offshore-wind-turbine-tower-fatigue`**  
**研究路线：Bhattacharyya et al. 原论文 $x-t$ LATIN-PGD → 2D fiber beam-column offshore wind turbine tower**  
**阶段范围：围绕 tower LATIN-PGD reversed benchmark 中第 4 个 PGD enrichment mode 出现的 `fixed_point_not_converged`，完成 fixed-point 周期轨道确认、iteration-cap 排查、spatial under-relaxation 排查、raw-map defect 诊断、whole-time BE temporal minimisation 诊断、period-3 三相候选的 post-Gram–Schmidt usefulness 诊断，以及 top-10 residual-energy deterministic seed / basin-sensitivity 诊断。**  
**本阶段不包含：正式修改 `latin/` 核心算法、不绕过 fixed-point convergence gate、不把未收敛 raw pair 接受为 persistent mode。**  
**下一阶段建议：seed / basin sensitivity 已完成并基本排除为主要原因。下一阶段优先回到 Bhattacharyya et al. 原论文 Eq. (59)/(72) 的 DG0 temporal discretisation，恢复其 jump / trace / element-wise temporal algebra，并在不改变 spatial problem、shifted defect 与 transaction semantics 的前提下检验 period-3 orbit 是否与当前 backward-Euler temporal discretisation 的偏离有关；在这一 paper-fidelity 路线闭合前，不优先引入 Anderson/Aitken 或 cycle-aware acceptance。**

---

# 1. 本阶段定位

在本阶段开始前，tower LATIN-PGD 的理论、状态、事务与模块接口已经完成冻结，Stage 1–7B 的核心模块与单元测试已经闭合，Stage 8A-1 完成了 tower elastic LATIN initialization，Stage 8A-2A 完成了真实 nonlinear activation probe。

Stage 8A-2A 已经确认：

- benchmark 不是 elastic-only；
- 第 1 个 PGD mode 可以正常生成；
- 第 1 个 enrichment 的 raw fixed point 可以收敛；
- tower LATIN-PGD transaction solver 可以完成 Trial A → enrichment → Trial B → atomic commit；
- plastic strain 与 damage 均被真实激活；
- 第 1 mode 可以降低 mechanical residual 与 LATIN indicator。

因此，本阶段并不是从“tower LATIN-PGD 是否能工作”开始，而是从更具体的问题开始：

> 当 outer LATIN iteration 推进到第 8 次 attempted iteration、已有 3 个 accepted PGD modes 后，为什么第 4 个 PGD enrichment mode 在 Eq. (70)–(72) alternating fixed-point loop 中无法收敛？

本阶段的目标不是立即给 fixed-point 加 accelerator，而是按“先诊断、后修改”的原则，逐层排除可能原因。

---

# 2. 本阶段资料边界与结论类别

本阶段继续严格区分四类来源。

| 类别 | 本阶段含义 |
| --- | --- |
| 原论文明确内容 | Bhattacharyya et al. 在 Eq. (61)–(72) 中采用 new separated pair enrichment；固定 temporal function 求 spatial problem，固定 spatial function 求 temporal minimisation，并通过 fixed-point algorithm 交替求解；fixed-point 后再进行 Gram–Schmidt orthonormalisation 与 time-function correction。 |
| 由原论文与当前离散结构推导 | complete-pair fixed-point indicator、lag-$k$ pair distance、relaxed-map 与 raw-map defect 的区别、periodic orbit 的判别。 |
| current tower implementation | residual-driven deterministic seed；Eq. (70)–(71) tower spatial solve；Eq. (72) sequential backward-Euler temporal solve；complete-pair graph-norm convergence；fixed-point convergence 通过后才进入 M-weighted Gram–Schmidt、temporal coordinate transformation、all-mode temporal re-optimisation 与 residual-benefit gate。 |
| 本阶段 diagnostic-only engineering experiments | 120-sweep pair-lag 诊断、constant spatial under-relaxation、$\omega$ sweep、raw-map defect、whole-time coupled BE temporal minimisation、对未收敛 period-3 raw pairs 的 post-basis-management usefulness 检查。 |

必须强调：

> 本阶段所有 relaxation、whole-time temporal solve 与 unconverged-pair usefulness test 都是诊断，不是已经接受的 tower LATIN-PGD 正式算法。

---

# 3. 与本阶段直接相关的既有理论冻结

## 3.1 new separated pair

对 Trial-A fixed-basis update 后仍未消除的 shifted defect，新 mode 写为：

$$ \Delta\dot{\vec{\varepsilon}}^p(t)=\dot{\lambda}(t)\vec p $$

以及：

$$ \Delta\vec{\sigma}'(t)=\lambda(t)\vec s $$

其中：

$$ \vec s=C_0(\mathcal E_{\rm tower}-I)\vec p $$

因此 raw new pair 的 canonical data 为：

$$ \mathcal P=\{\vec p,\vec s,\lambda(t),\dot\lambda(t)\} $$

## 3.2 tower v1 alternating fixed-point sweep

既有阶段已冻结一轮完整 inner sweep 为：

```text
complete pair at iteration k
    ↓
fix temporal function
    ↓
Eq. (70)–(71) spatial half-step
    ↓
M-norm scale normalization
    ↓
equilibrium stress reconstruction
    ↓
Eq. (72) temporal half-step
    ↓
new complete pair
    ↓
complete-pair convergence check
```

fixed-point 内不做 Gram–Schmidt；Gram–Schmidt 属于 fixed-point 收敛之后的 basis-management。

## 3.3 complete-pair convergence indicator

tower v1 采用 complete mechanical pair graph norm，而不是单独比较 $\vec p$ 或 $\lambda$。

定义：

$$ \|z\|_{\rm fp,h}^2=\sum_{n=1}^{N}\Delta t_n\left[(\Delta\dot{\vec{\varepsilon}}^p_n)^TM D_{H,n}^{-1}\Delta\dot{\vec{\varepsilon}}^p_n+(\Delta\vec{\sigma}'_n)^TM D_{H,n}\Delta\vec{\sigma}'_n\right] $$

symmetric relative pair change 为：

$$ \chi_{\rm fp}^{(k+1)}=\frac{\|z^{(k+1)}-z^{(k)}\|_{\rm fp,h}}{\|z^{(k+1)}\|_{\rm fp,h}+\|z^{(k)}\|_{\rm fp,h}} $$

provisional convergence tolerance 为：

$$ \varepsilon_{\rm fp}=10^{-6} $$

这个 $\chi_{\rm fp}$ 是 tower v1 numerical specification，不是原论文显式给出的公式。

## 3.4 fixed-point 后的正式 acceptance chain

既有 Stage 4 / post-fixed-point 冻结要求：

```text
raw fixed point converged
    ↓
M-weighted modified Gram–Schmidt
    ↓
exact temporal-coordinate transformation
    ↓
field-invariance checks
    ↓
spatial novelty / temporal significance
    ↓
tentative enlarged basis
    ↓
all-mode Eq. (58)–(59) temporal re-optimisation with full forcing
    ↓
full residual benefit
    ↓
accept / reject
```

所以正式实现中：

> `fixed_point_not_converged` 必须在进入 Gram–Schmidt 之前失败。

因此，本阶段后续对 period-3 raw pair 做 Gram–Schmidt 与 residual-benefit 检查，仅用于判断其“潜在 usefulness”，不能作为正式 acceptance。

---

# 4. Stage 8A reversed benchmark 与共同 baseline

本阶段所有诊断均围绕同一个 coarse but nonlinear tower benchmark：

```text
loading               : fully reversed sinusoidal top horizontal force
force amplitude       : 1.0 MN
period                : 10
cycles                : 1
increments / cycle    : 40
Nt                    : 41
elements              : 10
Gauss points/element  : 2
fibers/Gauss point    : 16
material points Nq    : 320
PGD format            : original x-t separated representation
```

elastic initialization：

```text
max |stress|                    = 1.165343815e+02 MPa
max free equilibrium residual   = 2.361617263e-06 N
```

当前 material reference parameters 中：

```text
E       = 134000 MPa
sigma_y = 80 MPa
```

因此 elastic peak stress 已超过 reference yield stress，后续 local stage 真实激活 viscoplasticity 与 damage 是合理的。

---

# 5. 多 outer-iteration convergence probe 中第 4 mode 问题的首次暴露

使用：

```text
max_iterations             = 20
max_fixed_point_iterations = 30
mode_significance_tolerance = 0
acceptance_tolerance        = 0
```

得到：

```text
termination_reason     = enrichment_failed
failure_reason         = fixed_point_not_converged
attempted iterations   = 8
committed iterations   = 7
trial evaluations      = 11
final basis modes      = 3
total modes added      = 3
final accepted xi      = 8.242691499e-04
```

Trial history：

```text
 trial   kind   rank        xi             zeta          reduced_residual
  1       A      0     1.927159484e-02   9.621855550e-01   1.000000000e+00
  2       B      1     1.608886240e-02   9.683317808e-01   8.127092612e-01
  3       A      1     1.528515774e-02   2.561688493e-02   9.911719342e-01
  4       B      2     1.343977647e-02   8.971242956e-02   1.608913901e-01
  5       A      2     5.717259968e-03   4.031164488e-01   5.105877147e-01
  6       A      2     3.728415875e-03   2.105560393e-01   8.832645905e-01
  7       A      2     3.331378110e-03   5.623928486e-02   1.007887944e+00
  8       B      3     3.077266559e-03   9.567729936e-02   2.260821376e-01
  9       A      3     1.195089283e-03   4.405478723e-01   6.050413282e-01
 10       A      3     8.242691499e-04   1.836326465e-01   8.944343629e-01
 11       A      3     7.633968263e-04   3.834076214e-02   9.825719402e-01
```

persistent commit history：

```text
 commit   kind   modes_added   accepted_xi
   1       B         1        1.608886240e-02
   2       B         1        1.343977647e-02
   3       A         0        5.717259968e-03
   4       A         0        3.728415875e-03
   5       B         1        3.077266559e-03
   6       A         0        1.195089283e-03
   7       A         0        8.242691499e-04
```

这组结果同时说明两件事。

第一，前 3 个 PGD enrichment modes 都是有实际作用的：

- mode 1：reduced residual $1.0000\rightarrow0.8127$；
- mode 2：$0.9912\rightarrow0.1609$；
- mode 3：$1.0079\rightarrow0.2261$。

第二，第 8 次 attempted outer iteration 的 Trial A：

```text
xi               = 7.633968263e-04
zeta             = 3.834076214e-02
reduced residual = 9.825719402e-01
```

虽然 Trial-A 的 $\xi$ 比 persistent $\xi_i=8.242691499\times10^{-4}$ 更小，但：

- $\zeta\lt0.1$；
- reduced residual 仍接近 1；
- 因而 saturation logic 要求继续 enrichment；
- 第 4 mode enrichment 随后发生 `fixed_point_not_converged`；
- 按 transaction freeze，整个当前 iteration 必须 rollback 到原 persistent baseline。

因此返回结果仍为：

```text
basis rank = 3
accepted xi = 8.242691499e-04
```

而不是提交 Trial A 的 $7.633968263\times10^{-4}$。

这实际验证了 strict transaction rollback 在真实 nonlinear tower case 中按设计工作。

---

# 6. 第 4 mode 在 30 fixed-point sweeps 下的异常特征

第一次失败时：

```text
fixed-point iterations = 30
fixed-point converged  = False
```

history：

```text
0.522627 0.590739 0.564442
0.659674 0.621418 0.576535
0.662147 0.618225 0.575457
0.662878 0.615822 0.574838
0.663034 0.613767 0.574520
...
0.661225 0.607721 0.575075
```

这里最初看到的现象是：

> $\chi_{\rm fp}$ 并没有单调下降，而是在大约 $0.57\sim0.66$ 的三个水平之间循环。

但仅凭 scalar $\chi_{\rm fp}$ 的三相变化，还不能证明 complete pair 本身存在 period-3 orbit，因为不同 pair 也可能恰好产生相似的 scalar distance。

所以第一步不能直接修改算法，而必须先区分：

```text
scalar indicator looks periodic
```

与：

```text
complete separated pair itself is periodic
```

---

# 7. iteration-cap 排查：30 → 120

为了排除“只是 30 次 inner iteration 太少”，保持所有方程不变，只把：

```text
max_fixed_point_iterations: 30 → 120
```

结果 outer solver history 完全不变：

```text
termination_reason   = enrichment_failed
committed iterations = 7
final basis modes    = 3
final accepted xi    = 8.242691499e-04
```

而 120-sweep tail 稳定趋向：

```text
0.658344 0.603606 0.577315
0.658330 0.603590 0.577328
0.658319 0.603576 0.577339
0.658309 0.603564 0.577348
0.658300 0.603554 0.577357
0.658293 0.603545 0.577364
0.658286 0.603537 0.577370
```

这说明：

> 30 次不够不是主要原因；迭代并非缓慢趋于 0，而是在形成稳定的多周期结构。

因此继续简单增加：

```text
120 → 300 → 1000
```

没有理论依据。

---

# 8. complete-pair lag diagnostic：确认真正的 period-3 orbit

为了判断 complete pair 是否真的三周期，定义 lag-$\ell$ distance：

$$ d_\ell^{(k)}=d\left(z^{(k)},z^{(k-\ell)}\right) $$

其中 $d$ 使用与 complete-pair convergence 一致的 graph-norm relative distance。

对最后 12 个 sweeps 得到：

```text
 sweep    lag-1 distance     lag-2 distance     lag-3 distance
 109      6.583089223e-01    6.034864660e-01    1.702011807e-04
 110      6.035644134e-01    5.772652195e-01    1.409671181e-04
 111      5.773482737e-01    6.582502141e-01    1.543888463e-04
 112      6.583002658e-01    6.034864196e-01    1.477882030e-04
 113      6.035541108e-01    5.772844543e-01    1.224165671e-04
 114      5.773565890e-01    6.582492788e-01    1.340919109e-04
 115      6.582927443e-01    6.034863729e-01    1.283458293e-04
 116      6.035451651e-01    5.773011635e-01    1.063215854e-04
 117      5.773638216e-01    6.582484575e-01    1.164768641e-04
 118      6.582862084e-01    6.034863275e-01    1.114756692e-04
 119      6.035373965e-01    5.773156800e-01    9.235361185e-05
 120      5.773701121e-01    6.582477372e-01    1.011860459e-04
```

因此：

$$ d_1=O(10^{-1}) $$

$$ d_2=O(10^{-1}) $$

但：

$$ d_3=O(10^{-4}) $$

并继续下降。

这给出了本阶段第一个强结论：

> **第 4 mode 的 raw Eq. (70)–(72) alternating iteration 不只是 scalar $\chi_{\rm fp}$ 呈三相变化，而是完整 separated mechanical pair 本身正在逼近 period-3 orbit。**

即：

$$ z^{(k)}\not\approx z^{(k-1)} $$

$$ z^{(k)}\not\approx z^{(k-2)} $$

但：

$$ z^{(k)}\approx z^{(k-3)} $$

---

# 9. spatial-only constant under-relaxation 诊断

确认 period-3 后，第一种最保守的 stabilization 尝试是只对 spatial half-step 做 under-relaxation。

设原始 spatial map 输出：

$$ \vec p_*^{(k+1)} $$

构造：

$$ \vec p_{\rm mix}^{(k+1)}=(1-\omega)\vec p^{(k)}+\omega\vec p_*^{(k+1)} $$

随后：

1. 做 sign alignment，避免 $\vec p$ 与 $-\vec p$ 的 gauge ambiguity；
2. 做 $M$-norm normalization；
3. 重新通过 equilibrium operator 计算 $\vec s$；
4. 再执行原 Eq. (72) temporal solve。

这个选择的目的，是保证每次 iterate 仍然是一个 rank-one separated pair。

必须强调：

> 直接对两个完整 rank-one corrections 做凸组合，一般会产生 rank-2 field，因此本阶段没有采用 complete-field convex mixing。

---

# 10. $\omega=0.5$：振荡幅值减小但形成约 6 相循环

取：

```text
omega = 0.5
max_fixed_point_iterations = 120
```

outer solver 仍然：

```text
termination_reason   = enrichment_failed
failure_reason       = fixed_point_not_converged
committed iterations = 7
final basis modes    = 3
final accepted xi    = 8.242606246e-04
```

与 unrelaxed baseline 几乎一致。

last fixed-point history：

```text
0.348180 0.363320 0.368210 0.361086 0.353841 0.352452
0.348181 0.363342 0.368200 0.361082 0.353833 0.352452
0.348182 0.363364 0.368190 0.361077 0.353826 0.352452
...
```

可以看到约 6 项重复。

因此：

> $\omega=0.5$ 把相邻 pair change 从约 $0.58\sim0.66$ 压到约 $0.35\sim0.37$，但没有形成 contraction；原 period-3 行为被转换成较小振幅、较长周期的 orbit。

---

# 11. $\omega=0.25$：进一步缩幅但形成约 13 相循环

取：

```text
omega = 0.25
max_fixed_point_iterations = 240
```

仍然：

```text
fixed-point converged  = False
fixed-point iterations = 240
basis rank             = 3
accepted xi            = 8.242494729e-04
```

tail：

```text
0.168213 0.177397 0.183189 0.185839 0.186788 0.184387
0.178870 0.172302 0.165197 0.158476 0.154183 0.154153 0.159176
0.168091 0.177300 0.183141 0.185817 0.186792 0.184441
0.178947 0.172387 0.165286 0.158549 0.154213 0.154123 0.159085
...
```

前后两组大约每 13 个 sweep 重复。

因此从：

```text
omega = 1.00 → period ≈ 3
omega = 0.50 → period ≈ 6
omega = 0.25 → period ≈ 13
```

可以看到一个非常明确的趋势：

> 更强的 damping 在缩小单步变化的同时拉长 orbit period，但并未显现真正 fixed-point convergence。

---

# 12. 多 $\omega$ sweep：relaxed step 越小不等于 fixed-point defect 越小

为了避免继续手工试：

```text
0.20, 0.10, 0.05, ...
```

对完全相同的第 4 mode failing Trial-A residual 开展：

```text
omega = 1.00, 0.75, 0.50, 0.35, 0.25, 0.20, 0.15, 0.10, 0.05
```

每个 $\omega$ 运行 400 fixed-point sweeps。

最初得到：

```text
omega   converged   fp_iters     last_chi        min_tail_chi
1.000     False       400      6.582428e-01     5.774120e-01
0.750     False       400      4.929815e-01     4.887826e-01
0.500     False       400      3.609195e-01     3.482261e-01
0.350     False       400      2.541611e-01     2.286154e-01
0.250     False       400      1.735754e-01     1.537749e-01
0.200     False       400      1.209854e-01     1.190932e-01
0.150     False       400      1.033554e-01     8.647513e-02
0.100     False       400      7.339993e-02     5.577954e-02
0.050     False       400      3.293132e-02     2.699633e-02
```

如果只看 `last_chi`，很容易得到一个错误印象：

> $\omega$ 越小，fixed point 越接近收敛。

但对于 relaxed iteration：

$$ z^{k+1}=z^k+\omega\left[F(z^k)-z^k\right] $$

relaxed step 本身天然满足：

$$ z^{k+1}-z^k=\omega\left[F(z^k)-z^k\right] $$

所以当 $\omega\rightarrow0$ 时，即使原始 map defect：

$$ F(z^k)-z^k $$

并没有下降，relaxed step 也会被人为压小。

因此 under-relaxation 后必须区分：

$$ \chi_{\rm relaxed}=d(z^{k+1},z^k) $$

与：

$$ \chi_{\rm raw}=d(F(z^k),z^k) $$

真正的 fixed point 要求：

$$ \chi_{\rm raw}\rightarrow0 $$

而不是仅仅：

$$ \chi_{\rm relaxed}\rightarrow0 $$

---

# 13. raw-map defect 诊断：正式排除 constant spatial-only under-relaxation

对每一个 relaxed iterate，额外计算原始 unrelaxed map 输出 $F(z^k)$，得到：

```text
omega   relaxed_last    raw_last        raw_min_tail     raw_mean_tail
1.000   6.5824280e-01   6.5824280e-01   5.7741202e-01   6.1323813e-01
0.750   4.9298154e-01   5.7874138e-01   5.7538345e-01   6.1726043e-01
0.500   3.6091952e-01   5.8308388e-01   5.8012701e-01   6.1776907e-01
0.350   2.5416112e-01   5.7775890e-01   5.7764116e-01   6.1833952e-01
0.250   1.7357537e-01   5.9227229e-01   5.7899007e-01   6.1978279e-01
0.200   1.2098538e-01   6.5471146e-01   5.7841870e-01   6.2098434e-01
0.150   1.0335537e-01   5.8524895e-01   5.7866142e-01   6.1749296e-01
0.100   7.3399931e-02   5.8064600e-01   5.7888543e-01   6.1841282e-01
0.050   3.2931317e-02   5.8913881e-01   5.7931697e-01   6.1435048e-01
```

最关键的是：

$$ \chi_{\rm relaxed}:0.658\rightarrow0.033 $$

但：

$$ \overline{\chi}_{\rm raw,tail}\approx0.613\sim0.621 $$

基本不变。

因此可以正式排除：

> **constant spatial-only under-relaxation 并没有解决第 4 mode 的 fixed-point instability；它只是按 $\omega$ 缩小 relaxed step。**

这也是为什么不能通过把 $\omega$ 继续降到 $10^{-2}$、$10^{-3}$ 来“制造”一个小 $\chi_{\rm relaxed}$ 并宣布收敛。

---

# 14. whole-time coupled BE temporal minimisation 诊断

在排除 simple damping 后，需要检查当前 Eq. (72) sequential BE temporal half-step 是否是 period orbit 的主要来源。

当前 tower v1 Eq. (72) 使用：

$$ \lambda_n=\frac{\vec g_n^TMD_{H,n}^{-1}\vec b_n}{\vec g_n^TMD_{H,n}^{-1}\vec g_n} $$

并按：

```text
lambda_{n-1} known
    ↓
solve lambda_n
    ↓
advance to n+1
```

逐时间步进行 causal least-squares。

为了只改变一个因素，本阶段构造 diagnostic whole-time BE solve：

- 保持相同 BE derivative；
- 保持 $\lambda_0=0$；
- 保持相同 shifted defect；
- 保持相同 Eq. (70)–(71) spatial solve；
- 只把 $\lambda_1,\ldots,\lambda_N$ 一次性作为 coupled unknowns；
- 对全部 time slabs 的 discrete BE residual 做一次 joint least-squares。

离散 residual：

$$ \vec r_n=\left(\frac{\vec p}{\Delta t_n}-D_{H,n}\vec s\right)\lambda_n-\frac{\vec p}{\Delta t_n}\lambda_{n-1}+\vec\Delta_n $$

然后求：

$$ \{\lambda_1,\ldots,\lambda_N\}=\operatorname*{arg\,min}\sum_{n=1}^{N}\Delta t_n\,\vec r_n^TMD_{H,n}^{-1}\vec r_n $$

必须强调：

> 这个 diagnostic 仍然是 BE-discretised whole-time minimisation，不等于恢复了 Bhattacharyya 原论文 exact DG0 temporal algebra。

---

# 15. whole-time BE 结果：sequential temporal marching 不是主要原因

结果：

```text
accepted                = False
failure_reason           = fixed_point_not_converged
fixed-point converged    = False
fixed-point iterations   = 200
last temporal condition  = 5.751257706e+01
```

fixed-point tail：

```text
0.556716 0.655343 0.591148
0.555916 0.653368 0.601165
0.556329 0.648622 0.610231
0.557004 0.641559 0.618318
0.557639 0.632280 0.625644
0.558243 0.620585 0.632396
0.559023 0.606319 0.638632
0.560461 0.590149 0.644271
0.563431 0.574597 0.649191
0.568940 0.563342 0.653158
...
```

它仍然在 $O(10^{-1})$ 量级发生周期性变化，没有趋于 $10^{-6}$。

同时 whole-time temporal least-squares 的 condition number 约：

$$ \kappa\approx57.5 $$

这个量级并没有显示严重 linear ill-conditioning。

因此可以得到本阶段第二个强排除结论：

> **把 current sequential BE temporal solve 改成 whole-time coupled BE least-squares，并不能消除第 4 mode 的 fixed-point cycle。**

所以：

> **sequential BE temporal marching 不是当前 period-3 instability 的主要原因。**

但不能据此声称：

> “原论文 DG0 已经被排除。”

因为该 diagnostic 没有恢复论文的 exact DG0 formulation。

---

# 16. 为什么不能把 `spatial novelty = 0` 等早期输出解释为第 4 mode 无效

所有 fixed-point failure 输出中都曾看到：

```text
spatial novelty       = 0
temporal significance = 0
orthogonality error   = 0
residual benefit      = 0
```

这些数字不能解释为：

```text
第 4 mode 和已有 basis 线性相关
```

或：

```text
第 4 mode 没有 residual benefit
```

原因是正式 enrichment 流程中：

```text
raw fixed point
    ↓
if not converged:
    return fixed_point_not_converged
```

只有 `converged=True` 后才会执行：

```text
Gram-Schmidt
gamma_sp
gamma_lambda
all-mode temporal reoptimization
residual benefit
```

所以 failure record 中的这些 0 只是“尚未计算”的 placeholder semantics。

---

# 17. period-3 orbit 三个 late phases 的 usefulness diagnostic

为判断第 4 mode 是：

```text
A. fixed point 不收敛，而且 candidate 本质无用
```

还是：

```text
B. fixed point 不收敛，但周期轨道中的 candidate 实际很有用
```

本阶段取 unrelaxed raw fixed-point 的最后三个 phases：

```text
sweep 118
sweep 119
sweep 120
```

虽然它们没有通过 fixed-point gate，但仅作为 diagnostic，分别执行：

```text
raw pair
    ↓
M-weighted modified Gram–Schmidt
    ↓
exact temporal-coordinate transformation
    ↓
field-invariance
    ↓
gamma_sp / gamma_lambda
    ↓
all-mode temporal re-optimisation with full forcing
    ↓
full mechanical residual benefit
```

得到：

```text
sweep    gamma_sp     orth_scale    gamma_lambda   ortho_error    resid_after     resid_benefit   relative_resid
 118    9.93134e-01   9.93134e-01   7.55466e-01   5.93031e-16   3.8349461e-03   9.8122483e-02   8.8615866e-01
 119    9.82551e-01   9.82551e-01   7.36538e-01   5.39853e-16   3.5520045e-03   1.6466285e-01   8.2078092e-01
 120    9.82616e-01   9.82616e-01   8.00417e-01   6.01531e-16   3.4602512e-03   1.8624078e-01   7.9956433e-01
```

---

# 18. usefulness diagnostic 的物理与数值解释

## 18.1 spatial novelty 很高

三个 phases：

$$ \gamma_{\rm sp}\approx0.983\sim0.993 $$

这说明 raw spatial candidate 在去除已有 3-mode spatial subspace 后，仍保留约 98%–99% 的 $M$-norm magnitude。

所以第 4 candidate 并不是：

```text
existing spatial basis 的近线性组合
```

相反，它在 spatial material-point space 中具有很强的新方向。

## 18.2 temporal significance 很高

三个 phases：

$$ \gamma_\lambda\approx0.736\sim0.800 $$

说明经过 exact coordinate transformation 后，新 temporal coordinate 并不 insignificant。

因此也不能解释为：

```text
Gram-Schmidt 后 new temporal function 几乎消失
```

## 18.3 basis orthogonality 非常健康

orthogonality error：

$$ O(10^{-16}) $$

说明 weighted MGS + exact transformation 本身没有暴露 basis-health 问题。

## 18.4 三个 phases 都能降低 full residual

full residual benefit：

$$ 9.81\% $$

$$ 16.47\% $$

$$ 18.62\% $$

其中 sweep 120 最好：

$$ {\rm relative\ residual}:0.98257194\rightarrow0.79956433 $$

所以第 4 mode 的 raw orbit 不是在几个“无用 candidate”之间循环。

本阶段第三个强结论是：

> **第 4 mode 的 period-3 orbit 中至少存在明显 useful candidates；failure 的本质不是 basis saturation，也不是 candidate insignificance，而是 raw alternating fixed-point map 无法选择并收敛到一个单一 mutually-consistent separated pair。**

---

# 19. 本阶段已经排除的解释

经过上述诊断，目前可以较有把握地排除以下解释。

## 19.1 “只是 fixed-point iteration 次数太少”

否。

30 → 120 后周期结构更清楚，并没有趋向 0。

## 19.2 “只是 scalar $\chi_{\rm fp}$ 看起来像三周期”

否。

complete-pair lag-3 distance 已降至 $O(10^{-4})$，而 lag-1 / lag-2 保持 $O(10^{-1})$。

## 19.3 “constant spatial under-relaxation 可以解决”

否。

当 $\omega$ 从 1.0 降到 0.05，relaxed step 大幅减小，但 raw-map defect 的 tail mean 始终约 $0.61\sim0.62$。

## 19.4 “当前 sequential BE temporal marching 导致周期振荡”

不是主要原因。

whole-time coupled BE temporal least-squares 仍然不收敛。

## 19.5 “第 4 spatial mode 已经被 existing 3-mode basis 表示”

否。

$$ \gamma_{\rm sp}\approx0.98\sim0.99 $$

说明 spatial novelty 很高。

## 19.6 “第 4 temporal coordinate insignificantly small”

否。

$$ \gamma_\lambda\approx0.74\sim0.80 $$

## 19.7 “第 4 mode 即使接受也几乎不降低 residual”

否。

三个 orbit phases 的 diagnostic residual benefit 为约 9.8%–18.6%。

---

# 20. 本阶段尚未排除的解释

虽然已经进行了多层诊断，但以下问题仍未闭合。

## 20.1 deterministic residual-row seed / basin sensitivity

当前 new-mode seed 优先取 shifted defect 中 residual energy 最大的 time row：

$$ n_*=\operatorname*{arg\,max}_n e_n $$

然后：

$$ \vec p_{\rm seed}\propto-\vec\Delta_{n_*} $$

这一 seed 是 current 1D / tower v1 engineering choice，不是原论文明确规定的唯一初始化。

因此仍有可能：

> 当前 argmax seed 落入 period-3 basin，而其他合理 deterministic seeds 能进入真实 fixed point basin。

这就是下一阶段最优先的诊断。

## 20.2 Eq. (70)–(72) 的当前 BE 离散与论文 DG0 的差异

whole-time BE 已经排除 sequentiality，但并未恢复 paper-exact DG0。

所以仍不能排除：

> 某些 BE-consistent contraction / temporal treatment 与论文 DG0 的差异，改变了 alternating map 的 dynamical properties。

但在 seed sensitivity 未闭合前，不应直接重写 DG0。

## 20.3 multiple fixed points / periodic attractor coexistence

当前 period-3 orbit 可能与一个或多个 fixed points 共存。

如果不同 seed 能收敛到 fixed point，则说明是 basin issue。

如果不同合理 seeds 都进入相同或相近 period-3 orbit，则说明 periodic attractor 更可能是该 discrete map 的鲁棒特征。

## 20.4 是否需要 fixed-point accelerator

Aitken、Anderson、quasi-Newton 或 cycle-aware methods 尚未测试。

目前不应优先进入这些方法，因为 seed dependence 这一更基础问题尚未排除。

---

# 21. 当前最强的阶段性数学判断

对当前第 4 mode，可以把 raw alternating solver 抽象为 fixed-point map：

$$ z^{k+1}=F(z^k) $$

当前观测到：

$$ d(z^{k+1},z^k)=O(10^{-1}) $$

$$ d(z^{k+2},z^k)=O(10^{-1}) $$

但：

$$ d(z^{k+3},z^k)\rightarrow O(10^{-4}) $$

因此最符合现有证据的描述是：

> current tower-v1 discrete Eq. (70)–(72) raw map 在该第 4 mode enrichment subproblem 上存在一个稳定或近稳定的 period-3 attracting orbit。

同时，orbit 上的三个 representatives 均具有：

$$ \gamma_{\rm sp}=O(1) $$

$$ \gamma_\lambda=O(1) $$

以及：

$$ {\rm residual\ benefit}>0 $$

所以不能把这个 orbit 解释为“围绕零 mode 的数值噪声”。

---

# 22. 为什么当前不能直接采用 sweep 120 的 candidate

sweep 120 的 candidate residual benefit 最大，约 18.62%。

直觉上可能会想到：

```text
既然它能降 residual，就直接接受它
```

但当前不能这样做。

原因有三层。

第一，formal algorithm freeze 要求：

```text
fixed-point converged
    ↓
basis management
```

而不是：

```text
periodic pair looks useful
    ↓
accept directly
```

第二，sweep 118、119、120 是三个不同 complete pairs。

如果直接选择 residual benefit 最大的 phase，相当于引入新的：

```text
cycle-phase selection rule
```

这不是原论文已有规则，也不是当前 tower v1 已冻结规则。

第三，如果未来真要设计 cycle-aware acceptance，就必须先回答：

- 如何稳定识别 cycle period？
- 为什么选最大 residual benefit phase？
- 是否会破坏 stationarity / separated-pair mutual consistency？
- 是否对 seed、mesh、load history 鲁棒？
- acceptance 后 outer LATIN convergence 是否仍稳定？
- 这种规则与原论文 fixed-point mathematical interpretation 如何对应？

这些都需要独立理论阶段，不能在本阶段暗中加入。

因此当前正确状态是：

> sweep 118–120 只能作为“第 4 mode 有潜在价值”的诊断证据，不能成为 persistent basis commit。

---

# 23. transaction semantics 在本阶段得到的再次验证

本阶段所有真实 solver runs 都继续遵守既有事务冻结。

persistent baseline：

$$ (s_i,B_m,\xi_i) $$

Trial A 建立 provisional fixed-basis：

$$ B_m\rightarrow B_m^A $$

若：

```text
Trial A saturation says enrichment required
```

则 enrichment 从：

$$ B_m^A $$

和：

$$ R_A $$

出发。

若 enrichment：

```text
fixed_point_not_converged
```

则：

- 不提交 Trial A；
- 不提交 unconverged new pair；
- 不保留 tentative basis changes；
- persistent state/basis/indicator 返回 iteration 开始前的 baseline。

因此第 8 attempted iteration failure 后：

```text
persistent rank = 3
persistent xi   = 8.242691499e-04
```

而不是 Trial-A：

```text
xi = 7.633968263e-04
```

这说明当前 strict rollback 逻辑在真正的 problematic enrichment case 中是有效的。

---

# 24. 与原论文流程的一致性与偏离边界

## 24.1 仍然一致的核心结构

当前研究仍保持：

- original $x-t$ PGD representation；
- fixed temporal → spatial problem；
- fixed spatial → temporal minimisation；
- alternating fixed-point enrichment；
- fixed-point 后 basis orthonormalisation；
- time-function correction；
- mode acceptance / rejection；
- existing-basis temporal re-optimisation；
- LATIN local/global outer alternation。

因此本阶段没有改变研究主线。

## 24.2 tower v1 numerical specification

以下属于 tower v1：

- BE temporal derivative；
- complete-pair graph-norm convergence；
- residual-row deterministic seed；
- material-point $M$ metric；
- strict transaction rollback；
- residual-benefit gate；
- one accepted pair maximum per enrichment event。

## 24.3 本阶段 diagnostic-only 偏离

以下没有进入正式算法：

- spatial under-relaxation；
- $\omega$ sweep；
- raw-map defect instrumentation；
- whole-time coupled BE temporal solve；
- 对 unconverged raw pairs 强制执行 post-Gram–Schmidt usefulness evaluation。

后续文档与论文叙述中必须继续保留这一来源层级，不能把 diagnostic experiment 写成“原论文算法”。

---

# 25. 当前代码层面的正式逻辑与诊断逻辑

正式模块：

```text
latin/tower_pgd_enrichment.py
```

仍然维持：

```text
seed
→ raw fixed point
→ require converged
→ post-fixed-point MGS
→ significance
→ all-mode temporal reoptimisation
→ residual benefit
→ candidate result
```

本阶段没有修改该文件。

所有新增逻辑均位于 `examples/` diagnostic probes 中，通过：

- monkey-patching private helper；
- local diagnostic reconstruction；
- read-only post-processing；

来观察算法行为。

所以本阶段的代码策略是：

> **diagnose around the core, do not mutate the core until the mechanism is identified.**

---

# 26. 本阶段新增/修改的 diagnostic files

当前工作区中围绕第 4 mode 诊断陆续出现了以下 example-level 文件。

```text
examples/tower_latin_pgd_reversed_convergence_probe.py
```

用途：

- multi-outer-iteration baseline；
- 30 / 120 fixed-point history；
- pair lag-1 / lag-2 / lag-3 instrumentation。

```text
examples/tower_latin_pgd_fixed_point_relaxation_probe.py
```

用途：

- spatial under-relaxation；
- $\omega=0.5$；
- 后改为 $\omega=0.25$。

```text
examples/tower_latin_pgd_fixed_point_relaxation_sweep.py
```

用途：

- $\omega$ sweep；
- raw-map defect 与 relaxed-step distinction。

```text
examples/tower_latin_pgd_whole_time_temporal_enrichment_probe.py
```

用途：

- whole-time coupled BE temporal minimisation diagnostic。

```text
examples/tower_latin_pgd_fourth_mode_cycle_usefulness_probe.py
```

用途：

- period-3 late phases；
- diagnostic MGS；
- $\gamma_{\rm sp}$；
- $\gamma_\lambda$；
- residual benefit。

```text
examples/tower_latin_pgd_fourth_mode_seed_sensitivity_probe.py
```

用途：

- 对同一个 failing fourth-mode Trial-A shifted defect；
- 按 residual-row energy 从高到低选取前 10 个 deterministic seeds；
- 每个 seed 独立运行 original unrelaxed Eq. (70)–(72) raw fixed-point map；
- 记录 convergence、last complete-pair change、lag-3 distance；
- 对最终 unconverged raw pair 仅作 diagnostic post-basis-management usefulness evaluation；
- 判断当前 period-3 orbit 是否主要由单一 argmax seed / basin-of-attraction sensitivity 引起。

---

# 27. Git 状态与阶段保存原则

本阶段开始前，已推送的稳定 baseline 包括：

```text
9e5bd0744bf18df273c8be6144c522883619ba1e
feat: add tower elastic LATIN initialization
```

以及：

```text
0813a6fa5d5c05525ecdf542650ddd14d21279b0
feat: add tower LATIN-PGD activation probe
```

本阶段围绕第 4 mode 的 convergence / relaxation / whole-time / usefulness diagnostics 尚未形成新的正式核心 commit。

因此建议：

> 在 seed-sensitivity 诊断开展前，先把本阶段总结文档加入 `docs/`，作为当前理论与数值诊断 checkpoint。

是否把所有 diagnostic examples 同时纳入同一个 commit，应在查看 `git status --short` 后再决定；不建议把下载到仓库根目录的 `.patch` 临时文件提交到 Git。

---

# 28. deterministic residual-row seed sensitivity：诊断设计

在前述诊断完成后，尚未排除的最基础因素是 current deterministic seed / basin-of-attraction sensitivity。

current tower enrichment seed 取 shifted defect 中 residual energy 最大的 time row：

$$ e_n=\sum_q\frac{v_q}{H_{\sigma,nq}}\Delta_{nq}^2 $$

current default seed time index 为：

$$ n_*=\operatorname*{arg\,max}_n e_n $$

然后构造：

$$ \vec p_{\rm seed}\propto-\vec\Delta_{n_*} $$

这一做法来自 current 1D / tower-v1 numerical choice，不是 Bhattacharyya et al. 原论文规定的唯一 new-pair initialization。

因此需要回答：

> 当前 period-3 orbit 是否只是因为 single argmax residual-row seed 恰好落入一个 periodic basin，而其他合理 deterministic seeds 可以进入真正的 fixed-point basin？

为了只改变 seed，不改变 enrichment subproblem，本诊断固定：

- persistent seven-commit baseline；
- failing Trial-A rank-3 basis；
- shifted defect；
- $H_\sigma$；
- material-point metric $M$；
- tower equilibrium operator；
- original unrelaxed Eq. (70)–(72) alternating map；
- current sequential-BE temporal half-step；
- fixed-point tolerance；
- post-fixed-point formal gate。

唯一变化是 initial spatial seed。

按 $e_n$ 从大到小，选择 top 10 residual-energy time rows：

```text
rank 1 ... rank 10
```

分别令：

$$ \vec p_{\rm seed}^{(j)}\propto-\vec\Delta_{n_j} $$

每个 seed 最多运行 200 raw fixed-point sweeps。

对于每个 seed，记录：

```text
time_idx
energy / max_energy
converged?
last_chi
lag3
gamma_sp
gamma_lambda
diagnostic residual benefit
candidate relative residual
```

其中：

- `last_chi` 是 complete-pair lag-1 distance；
- `lag3` 是 $d(z^k,z^{k-3})$；
- `gamma_sp`、`gamma_lambda` 与 residual benefit 对 unconverged final pair 仍只是 diagnostic，不能绕过 formal fixed-point gate。

---

# 29. top-10 deterministic seeds 的数值结果

共同 baseline：

```text
termination = enrichment_failed
committed   = 7
rank        = 3
xi          = 8.242691499e-04
```

同一个 failing Trial-A relative residual：

```text
9.825719402e-01
```

top 10 residual-energy rows 的结果为：

```text
 rank  time_idx  energy/max      conv       last_chi           lag3     gamma_sp    gamma_lam  resid_benefit   relative_resid
    1        26 1.00000e+00     False  6.0348708e-01  2.0937103e-06  9.82546e-01  7.36506e-01  1.6469793e-01    8.2074645e-01
    2        27 4.45159e-01     False  6.0348703e-01  2.0107669e-06  9.82546e-01  7.36506e-01  1.6469796e-01    8.2074642e-01
    3        33 4.27624e-01     False  5.7741129e-01  1.7573459e-06  9.82654e-01  8.00501e-01  1.8618887e-01    7.9961535e-01
    4        10 3.37113e-01     False  6.5824533e-01  6.4798094e-06  9.93119e-01  7.55369e-01  9.8132793e-02    8.8614854e-01
    5        30 3.29181e-01     False  6.5824341e-01  1.5626691e-06  9.93119e-01  7.55364e-01  9.8133282e-02    8.8614806e-01
    6        13 2.34966e-01     False  5.7741133e-01  1.6635101e-06  9.82654e-01  8.00501e-01  1.8618882e-01    7.9961540e-01
    7        32 2.28306e-01     False  5.7739904e-01  3.1213014e-05  9.82643e-01  8.00476e-01  1.8620430e-01    7.9960019e-01
    8        25 1.90008e-01     False  6.0348749e-01  2.8197257e-06  9.82546e-01  7.36506e-01  1.6469764e-01    8.2074673e-01
    9        31 1.81869e-01     False  6.0348850e-01  4.6311118e-06  9.82546e-01  7.36507e-01  1.6469694e-01    8.2074742e-01
   10        29 1.75200e-01     False  6.5824520e-01  6.1372585e-06  9.93119e-01  7.55368e-01  9.8132827e-02    8.8614851e-01
```

最直接的事实是：

```text
10 / 10 seeds:
converged = False
```

因此在这一 deterministic residual-row family 中，没有发现能够进入 ordinary fixed point 的 alternative initialization。

---

# 30. seed sensitivity 结果揭示的是“同一三周期吸引子”，而不是十个不同失败轨道

这组结果最重要的地方，不只是 10 个 seeds 都失败，而是不同 seeds 的 late state 几乎严格聚类到此前已经识别的三个 period-3 phases。

可以把结果分为三个 phase families。

## 30.1 Phase A

对应：

```text
time_idx = 26, 27, 25, 31
```

共同特征：

```text
last_chi       ≈ 0.603487
gamma_sp       ≈ 0.982546
gamma_lambda   ≈ 0.736506
resid_benefit  ≈ 0.164698
relative_resid ≈ 0.820746
```

## 30.2 Phase B

对应：

```text
time_idx = 33, 13, 32
```

共同特征：

```text
last_chi       ≈ 0.577411
gamma_sp       ≈ 0.98265
gamma_lambda   ≈ 0.80050
resid_benefit  ≈ 0.1862
relative_resid ≈ 0.7996
```

## 30.3 Phase C

对应：

```text
time_idx = 10, 30, 29
```

共同特征：

```text
last_chi       ≈ 0.658245
gamma_sp       ≈ 0.993119
gamma_lambda   ≈ 0.75537
resid_benefit  ≈ 0.09813
relative_resid ≈ 0.88615
```

这三个 families 与此前 unrelaxed sweep 118–120 的三个 late phases 一一对应。

因此改变 initial residual-row seed 后，主要改变的是：

> 最终在 period-3 orbit 的哪一个 phase 上被 max-iteration cutoff 截止。

而不是：

> 进入完全不同的 attractor。

这比单纯的 `10/10 nonconverged` 更强。

---

# 31. lag-3 distance 将 period-3 robustness 推进到 $O(10^{-6})$

top-10 seed test 中：

$$ d_1=O(10^{-1}) $$

但绝大多数 seed 的 lag-3 distance 已达到：

$$ d_3=O(10^{-6}) $$

最大一项约：

$$ 3.12\times10^{-5} $$

而最小项约：

$$ 1.56\times10^{-6} $$

因此对于这些 deterministic initializations，有：

$$ F^3(z)\approx z $$

但：

$$ F(z)\not\approx z $$

相比此前 120-sweep single-seed diagnosis 的：

$$ d_3=O(10^{-4}) $$

现在通过更长 sweep 与多 seed family，period-3 attractor 的数值证据进一步增强了约两个数量级。

所以当前更准确的描述是：

> **在 top-10 residual-energy deterministic seed family 内，当前 discrete Eq. (70)–(72) raw alternating map 鲁棒地收敛到同一个 period-3 attracting orbit，而不是 ordinary fixed point。**

---

# 32. seed / basin sensitivity 应如何表述：基本排除主要原因，但不能做无限泛化

本诊断支持：

$$ \boxed{\text{single argmax residual-row seed 不是当前 failure 的主要原因}} $$

因为：

- rank-1 seed；
- 9 个 alternative high-energy seeds；

共 10 个合理 deterministic initializations 均未收敛，并且进入相同 three-phase orbit。

因此不再建议继续机械扩大：

```text
top 10 → top 20 → all 41 rows
```

这类 residual-row restart sweep。

但不能宣称：

$$ \boxed{\text{所有可能的 initial conditions 都不可能收敛}} $$

因为本阶段没有穷举：

- arbitrary linear combinations of residual rows；
- random spatial seeds；
- existing-mode-orthogonal seeds；
- nonlinear optimized seeds；
- exact unstable fixed-point continuation。

所以科学上应写为：

> **在当前最自然且与 existing implementation 一致的 deterministic residual-row seed family 内，没有证据支持 basin sensitivity 是 period-3 failure 的主要来源。**

而不是：

> “数学上不存在其他 fixed-point basin”。

---

# 33. seed-sensitivity 诊断再次证明 period-3 orbit 中的 candidates 是 genuinely useful

top-10 seed test 仍然得到三类稳定的 post-basis-management diagnostic values。

spatial novelty：

$$ \gamma_{\rm sp}\approx0.9825\sim0.9931 $$

temporal significance：

$$ \gamma_\lambda\approx0.7365\sim0.8005 $$

residual benefit：

$$ 0.0981\sim0.1862 $$

candidate relative residual：

$$ 0.7996\sim0.8861 $$

所以无论 seed 最终落到 period-3 的哪一个 phase：

- new spatial direction 都明显不同于 existing rank-3 spatial subspace；
- new temporal function 都不 insignificant；
- enlarged basis all-mode temporal re-optimisation 都能降低 full residual。

因此 seed-sensitivity test 不仅排除了一个 initialization hypothesis，还独立重复验证了：

> **period-3 orbit 不是围绕一个 insignificant / linearly-dependent mode 的数值噪声。**

---

# 34. 当前已完成的排除链

截至本次 seed-sensitivity test，第 4 mode failure 的排除链可以写成：

```text
fixed-point iteration cap too small?
    ↓ NO
30 → 120 sweeps 仍形成稳定三周期

scalar chi only looks periodic?
    ↓ NO
complete-pair lag-3 distance → very small

constant spatial under-relaxation?
    ↓ NO
relaxed step shrinks, raw-map defect remains O(0.6)

sequential BE temporal marching?
    ↓ NOT MAIN CAUSE
whole-time coupled BE least-squares still cycles

basis saturation / insignificant fourth mode?
    ↓ NO
gamma_sp ≈ 0.98–0.99
gamma_lambda ≈ 0.74–0.80
residual benefit ≈ 9.8%–18.6%

single residual-row seed / basin sensitivity?
    ↓ NO EVIDENCE AS MAIN CAUSE
top-10 deterministic high-energy residual-row seeds
all return to same period-3 attractor
```

因此当前问题已经被收敛到比最初 `fixed_point_not_converged` 更具体的位置：

$$ \boxed{\text{current discrete Eq. (70)–(72) raw alternating map itself}} $$

但“current discrete”仍包含一个重要的 project-specific choice：

$$ \boxed{\text{backward-Euler temporal discretisation}} $$

而原论文明确采用的是 DG0 temporal discretisation。

---

# 35. seed sensitivity 闭合后，下一研究路线的重新排序

seed sensitivity 闭合后，后续候选方向需要重新排序。

## 35.1 第一优先级：paper-fidelity DG0 reconstruction

当前最值得检查的明确 source-layer discrepancy 是：

```text
paper:
Eq. (59)/(72) temporal minimisation
+ zero-order discontinuous Galerkin temporal discretisation

current tower v1:
same residual/minimisation structure
+ project-validated backward Euler discretisation
```

whole-time BE diagnostic 已经说明：

```text
sequential vs global-in-time BE scope
```

不是主要原因。

但它没有回答：

```text
BE vs original DG0 discrete temporal algebra
```

是否改变 fixed-point map 的动力学性质。

因此下一阶段应回到原论文与其引用的 temporal discretisation资料，逐步恢复：

- temporal finite-element trial/test space；
- DG0 slab-wise unknown；
- inter-slab jump term；
- left/right trace；
- weak temporal derivative；
- Eq. (59) 的 discrete stationarity；
- Eq. (72) single-new-mode temporal equation；
- 与现有 BE formula 的准确关系。

目标不是立即替换 core，而是先得到：

$$ \boxed{\text{paper-faithful DG0 discrete Eq. (72)}} $$

再构造 diagnostic implementation。

## 35.2 第二优先级：fixed-point accelerator

只有在 paper-faithful DG0 仍出现相同 period-3 attractor 时，再优先考虑：

```text
Aitken
Anderson
quasi-Newton
```

但 accelerator 的作用空间必须先解决 rank-one gauge：

$$ (\vec p,\lambda)\equiv(c\vec p,\lambda/c) $$

以及 sign ambiguity：

$$ (\vec p,\lambda)\equiv(-\vec p,-\lambda) $$

所以不能直接把 unprocessed pair coordinates 当普通 Euclidean vector 做 Anderson acceleration。

## 35.3 第三优先级：cycle-aware candidate handling

例如：

```text
detect period-3
→ inspect three phases
→ choose best residual-benefit candidate
```

虽然当前数据已经证明某些 phases 很 useful，但这种规则会改变原论文的 fixed-point acceptance logic。

因此仍是较后方案。

---

# 36. 为什么现在优先 DG0，而不是立即 Anderson

选择 paper-fidelity DG0 优先，不是因为已经证明 BE 是错误的。

相反，当前证据是：

> current BE implementation 在多个测试中是数值自洽的，而且 whole-time BE 也不能解除周期。

优先 DG0 的原因是研究路线本身：

1. 当前第一版 tower LATIN-PGD 的目标是尽可能忠实迁移 Bhattacharyya et al. 原方法；
2. seed sensitivity、simple damping 与 sequential-vs-whole-time BE 已经逐层排除；
3. DG0 是目前尚未恢复、且原论文明确写出的 temporal discretisation；
4. 在添加新的 nonlinear accelerator 之前，应先确认这一 paper-explicit component 是否会改变 enrichment fixed-point behavior。

所以正确的问题不是：

> “BE 已经失败，因此必须换 DG0。”

而是：

> “在增加新的工程算法之前，应先恢复 paper-explicit DG0，确认 period-3 是原方法本身的离散 fixed-point pathology，还是 current BE discretisation 的特定表现。”

---

# 37. 当前正式算法仍保持不变

尽管诊断越来越深入，本阶段仍然没有理由修改正式 `latin/` core。

正式 logic 保持：

```text
Trial A
    ↓
need enrichment
    ↓
raw fixed point
    ↓
if fixed_point_not_converged:
    reject current transaction
    rollback
```

尤其不能因为：

```text
Phase B residual benefit ≈ 18.6%
```

就直接接受该 raw pair。

seed-sensitivity 结果反而说明：

- three phases 是 attractor 的系统性组成；
- 不是某个特殊 seed 偶然得到的 candidate；
- 如果未来选择其中一相，就必须显式设计新的 cycle-selection theory。

所以 formal convergence gate 继续冻结。

---

# 38. 本阶段最终结论

截至 2026-08-19，第 4 个 PGD mode fixed-point pathology 已经得到以下阶段性结论：

1. current tower LATIN-PGD reversed nonlinear benchmark 中，前 3 个 PGD modes 均能正常生成并明显降低 reduced residual；
2. 第 8 次 attempted outer iteration 的 Trial A 需要继续 enrichment，因此触发第 4 mode；
3. 第 4 mode raw Eq. (70)–(72) fixed-point 在 30、120 乃至后续 200–400 sweeps 的各类诊断中都没有趋向 ordinary fixed point；
4. complete-pair lag analysis 证明该行为是真正的 period-3 orbit，而不是 scalar indicator 的假周期；
5. constant spatial-only under-relaxation 只能缩小 relaxed step，并不能降低 raw-map defect；
6. whole-time coupled BE temporal minimisation 仍然不能消除周期，因此 sequential BE time marching 不是主要原因；
7. period-3 三个 phases 的 spatial novelty 均约 0.98–0.99；
8. temporal significance 均约 0.74–0.80；
9. diagnostic full residual benefit 约为 9.8%–18.6%，所以第 4 mode 不是 basis saturation 或 insignificant enrichment；
10. top-10 residual-energy deterministic seeds 全部 `converged=False`；
11. 不同 seeds 并没有形成多个不同失败轨道，而是全部回到同一个 three-phase attractor；
12. top-10 test 中 lag-3 distance 已达到约 $10^{-6}$–$10^{-5}$，进一步加强 period-3 attracting-orbit 证据；
13. 因此 single argmax residual-row initialization / basin sensitivity 在这一合理 deterministic seed family 内已经基本排除为主要原因；
14. 但不能无限泛化为“所有可能 seeds 均不存在 fixed-point basin”；
15. 正式算法仍保持 `fixed_point_not_converged → reject / rollback`；
16. 当前下一项最重要的 paper-fidelity discrepancy 是 original DG0 temporal discretisation 尚未恢复；
17. 下一阶段应优先恢复 Eq. (59)/(72) 的 DG0 discrete algebra，再判断 period-3 是否仍存在；
18. 在 DG0 路线闭合前，不优先修改 core，不优先加入 Anderson/Aitken，也不采用 cycle-aware acceptance。

---

# 39. 更新后的阶段 checkpoint

本阶段现在应被冻结为：

> **Tower LATIN-PGD 第 4 mode 的 failure 已经从一般性的 `fixed_point_not_converged` 定位为一个对 top-10 high-energy deterministic residual-row seeds 均鲁棒的 complete-pair period-3 attracting orbit。增加 iteration cap、constant spatial under-relaxation、改变 sequential BE 为 whole-time coupled BE、以及更换合理 deterministic residual-row seed 均不能使 raw Eq. (70)–(72) map 收敛到 ordinary fixed point。与此同时，period-3 三个 phases 都具有很高 spatial novelty、temporal significance 和显著 residual benefit，因此 failure 不是 mode saturation 或 insignificance。formal fixed-point gate 仍不能绕过。下一阶段首先恢复 Bhattacharyya et al. 原论文 Eq. (59)/(72) 的 DG0 temporal discretisation，再决定是否需要 fixed-point accelerator 或 cycle-aware algorithm。**

---

# 40. 与既有阶段文档的衔接

本总结应与以下既有文档连续阅读：

```text
docs/2026-08-17-tower-latin-pgd-eq70-72-enrichment-fixed-point-convergence-stage-summary.md
```

用于理解：

- Eq. (70)–(72) alternating loop；
- deterministic seed 的来源层级；
- complete-pair convergence criterion；
- fixed-point 内 normalization。

以及：

```text
docs/2026-08-17-tower-latin-pgd-post-fixed-point-gram-schmidt-mode-acceptance-stage-summary.md
```

用于理解：

- fixed-point 后 weighted MGS；
- exact temporal-coordinate transformation；
- field invariance；
- $\gamma_{\rm sp}$ / $\gamma_\lambda$；
- all-mode temporal re-optimisation；
- full residual benefit；
- rollback 与 acceptance。

本文件当前已经覆盖：

```text
baseline fourth-mode failure
→ period-3 confirmation
→ iteration-cap exclusion
→ spatial under-relaxation exclusion
→ raw-map defect
→ whole-time BE temporal-scope diagnosis
→ three-phase usefulness diagnosis
→ top-10 deterministic seed / basin-sensitivity diagnosis
```

因此下一份独立理论阶段总结应从：

> **Bhattacharyya et al. Eq. (59)/(72) original DG0 temporal discretisation 的资料恢复与离散推导**

开始，而不再重复本文件已闭合的 fixed-point pathology 诊断链。
