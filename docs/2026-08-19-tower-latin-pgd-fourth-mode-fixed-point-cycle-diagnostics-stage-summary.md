# Tower LATIN-PGD 第 4 个 PGD mode 固定点周期振荡诊断阶段总结

**日期：2026-08-19**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**仓库：`Jiarui-Shi-DUT/Offshore_Wind_Turbine_LATIN_PGD`**  
**分支：`feature/offshore-wind-turbine-tower-fatigue`**  
**研究路线：Bhattacharyya et al. 原论文 $x-t$ LATIN-PGD → 2D fiber beam-column offshore wind turbine tower**  
**阶段范围：围绕 tower LATIN-PGD reversed benchmark 中第 4 个 PGD enrichment mode 出现的 `fixed_point_not_converged`，完成 fixed-point 周期轨道确认、iteration-cap 排查、spatial under-relaxation 排查、raw-map defect 诊断、whole-time BE temporal minimisation 诊断，以及 period-3 三相候选的 post-Gram–Schmidt usefulness 诊断。**  
**本阶段不包含：正式修改 `latin/` 核心算法、不绕过 fixed-point convergence gate、不把未收敛 raw pair 接受为 persistent mode。**  
**下一阶段建议：对相同第 4 mode shifted defect 开展 deterministic residual-row seed sensitivity / restart 诊断；在该诊断闭合前，不进入 Anderson/Aitken 或 cycle-aware acceptance 的正式设计。**

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

下一步计划文件：

```text
examples/tower_latin_pgd_fourth_mode_seed_sensitivity_probe.py
```

目前应在本阶段总结保存之后再继续应用/运行。

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

# 28. 下一阶段的唯一优先问题：seed sensitivity

在进入任何高级 accelerator 之前，下一阶段应先回答：

$$ \boxed{\text{period-3 orbit 是否依赖 current deterministic residual-row seed？}} $$

建议对完全相同的：

- persistent baseline；
- Trial-A state；
- $H_\sigma$；
- shifted defect $R_A$；
- existing rank-3 basis；
- Eq. (70)–(72) map；

仅改变 initial spatial seed。

优先采用 deterministic family：

```text
shifted-defect residual energy 最大的前 10 个 time rows
```

即：

$$ e_n=\sum_q\frac{v_q}{H_{\sigma,nq}}\Delta_{nq}^2 $$

按 $e_n$ 从大到小选取前 10 个 $n$，并令：

$$ \vec p_{\rm seed}^{(j)}\propto-\vec\Delta_{n_j} $$

每个 seed 都运行原始 unrelaxed fixed-point，不加 damping，不改 temporal solve。

应记录：

```text
seed time index
relative residual-row energy
fixed-point converged?
fixed-point iterations
last chi
lag-3 distance
gamma_sp
gamma_lambda
diagnostic residual benefit
candidate relative residual
```

判断规则：

```text
若至少一个 alternative seed 收敛
    ↓
说明存在 basin / initialization sensitivity
    ↓
优先研究 deterministic restart / multi-start
```

否则：

```text
若多个高能 residual seeds 均不收敛
且都进入相同/相近 period-3 orbit
    ↓
说明 periodic attractor 对合理 seed family 具有鲁棒性
    ↓
再进入 fixed-point stabilization / cycle handling 理论阶段
```

---

# 29. 如果 seed sensitivity 也失败，后续候选方向的优先级

只有在 seed family 诊断闭合后，才建议进入以下方向。

优先级 1：

```text
deterministic restart / multi-start
```

如果不同 seed 能进入不同 basin，这是最小改动方案。

优先级 2：

```text
Aitken / Anderson / quasi-Newton acceleration
```

前提是明确其作用对象：

- spatial factor？
- complete pair？
- gauge-fixed coordinates？
- raw map residual？

不能直接对任意 rank-one pair 做未经 gauge 处理的向量加速。

优先级 3：

```text
paper-exact DG0 temporal discretisation reconstruction
```

如果需要继续贴近原论文，应回到 paper Eq. (59)/(72) 的 DG0 formulation，重新推导 jump / trace / element-wise temporal algebra，而不是把当前 whole-time BE diagnostic 当成 paper DG0。

优先级 4：

```text
cycle-aware candidate handling
```

这是偏离最大的方案。

若未来考虑：

```text
detect period-p orbit
→ evaluate p orbit phases
→ choose candidate
```

必须单独建立理论与 transaction semantics，不能从本阶段 diagnostic 直接迁移。

---

# 30. 本阶段最终结论

截至 2026-08-19，第 4 个 PGD mode 问题可以总结为：

1. current tower LATIN-PGD 在 reversed 1-cycle nonlinear benchmark 中，前 3 个 PGD modes 均可正常生成并有效降低 residual；
2. 第 8 次 attempted outer iteration 需要第 4 mode enrichment；
3. 第 4 mode raw Eq. (70)–(72) fixed-point 在 30 和 120 sweeps 下均不收敛；
4. complete-pair lag analysis 证明它不是 scalar indicator 假周期，而是完整 separated pair 正在逼近 period-3 orbit；
5. constant spatial under-relaxation 只会缩小 relaxed step 并拉长 orbit period；
6. raw-map defect 在 $\omega=1.0\rightarrow0.05$ 范围内始终约为 $O(0.6)$，因此 constant under-relaxation 被排除；
7. whole-time coupled BE temporal minimisation 仍不能消除周期行为，因此 sequential BE marching 不是主要原因；
8. period-3 三个 late phases 均具有很高 spatial novelty 与 temporal significance；
9. 三个 phases 在 diagnostic post-basis-management 后都能降低 full residual，最佳约 18.6%；
10. 因此第 4 mode 不是 basis saturation、linear dependence 或 insignificant mode；
11. 当前最符合证据的解释是：**该 discrete raw alternating map 在当前 enrichment subproblem 上存在稳定/近稳定 period-3 attracting orbit，而 orbit 中包含 genuinely useful enrichment directions；**
12. 正式算法仍必须保持 `fixed_point_not_converged → reject/rollback`，不能直接提交某个周期 phase；
13. 下一步首先应排除 deterministic seed / basin-of-attraction sensitivity；
14. 在 seed sensitivity 闭合前，不应贸然加入 Anderson/Aitken、cycle-aware acceptance 或修改 core `latin/` 实现。

---

# 31. 阶段 checkpoint

本阶段应被冻结为以下一句话：

> **Tower LATIN-PGD 第 4 mode 的失败已经从“普通 fixed-point 不收敛”定位为“raw Eq. (70)–(72) complete-pair period-3 attracting orbit”；iteration cap、constant spatial under-relaxation、sequential BE temporal marching、basis saturation 与 mode insignificance 均已通过针对性诊断被排除或显著削弱。period-3 orbit 中存在可显著降低 full residual 的 useful candidates，但 formal fixed-point gate 仍不能绕过。下一阶段首先检验 seed/basin sensitivity，再决定是否需要 fixed-point accelerator 或更接近原论文 DG0 的重构。**

---

# 32. 与既有阶段文档的衔接

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

本文件的定位不是替代上述理论冻结，而是记录：

> **当这些理论与算法 specification 第一次在 nonlinear tower reversed benchmark 的第 4 enrichment mode 上出现真实 fixed-point pathology 时，我们如何逐层诊断、排除错误解释，并把下一步研究问题收敛到 seed/basin sensitivity。**

