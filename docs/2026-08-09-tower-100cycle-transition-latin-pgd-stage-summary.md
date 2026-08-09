# 海上风机塔筒 100 圈非对称循环转折分析与 LATIN-PGD 方法学阶段总结

> 日期：2026-08-09
> 项目：Offshore_Wind_Turbine_LATIN_PGD
> 分支：`feature/offshore-wind-turbine-tower-fatigue`
> 本文档性质：对 50 圈阶段总结之后新增的 100 圈非对称循环长期演化结果进行专项总结
> 关联前序文档：`docs/2026-08-08-tower-50cycle-transition-latin-pgd-summary.md`

---

# 0. 本阶段为什么需要单独形成总结

在前一阶段的 50 圈非对称循环分析中，已经观察到一个新的长期演化信号：

$$ \Delta D_n $$

在约第 46 圈附近达到局部最低值，随后到第 50 圈出现轻微回升。

当时的关键数据为：

$$ \Delta D_{46}=3.466680\times10^{-4} $$

$$ \Delta D_{50}=3.484715\times10^{-4} $$

从第 46 圈到第 50 圈的增幅仅约为：

$$ 0.52\% $$

因此，50 圈阶段只能谨慎地提出：

> 每循环损伤增量可能正在由持续衰减转向平台或回升，但现有数据不足以区分真实机制转折与局部波动。

这直接产生了本阶段 100 圈分析的核心问题：

> **第 46 圈附近出现的损伤增量最低点究竟是偶然波动，还是一个真实且持续的损伤演化转折点？**

同时，50 圈阶段还出现了第二个值得关注的现象：

$$ r_p=\frac{|\Delta\varepsilon_p|}{L_p} $$

即 `net/path` 已经从接近 1 明显下降到：

$$ r_{p,50}^{c}=0.767281 $$

这说明一个循环内部的塑性应变路径已经不再接近纯单向累积，而开始出现明显的反向塑性活动。

因此，100 圈分析不仅用于确认损伤增量的转折，还用于回答另一个对 LATIN-PGD 更关键的问题：

> **随着循环继续发展，循环内部的快尺度塑性路径是否仍可视为固定周期波形，还是会随慢尺度循环数持续演化？**

本阶段最终得到的答案是：

> **第 46 圈附近的损伤增量最低点是真实且持续的转折；此后每循环损伤增量稳定上升，但净位移漂移和净塑性应变漂移直到第 100 圈仍继续下降。与此同时，循环内部正反向塑性活动显著增强，使 `net/path` 在 coupled case 中下降至约 0.492。**

因此，本阶段标志着研究认识从：

```text
“是否存在长期转折？”
```

推进到：

```text
“损伤速率已经发生转折，但宏观净漂移尚未反转，
且快尺度循环内塑性路径正在持续演化。”
```

这一认识对后续 LATIN-PGD 的时间分离方式、阶段化基函数以及自适应增广策略具有直接方法学意义。

---

# 1. 100 圈全阶参考分析配置

本次 100 圈分析保持与此前 20 圈、50 圈分析相同的结构模型、材料模型、加载比例和时间离散，仅将循环数延长至 100 圈。

计算配置如下：

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
increments_per_cycle = 40
n_cycles = 100
```

非对称循环加载关系为：

$$ R_F=\frac{F_{\min}}{F_{\max}}=-0.5 $$

因此：

$$ F_{\min}=-0.5F_{\max} $$

$$ F_{\mathrm{mean}}=\frac{F_{\max}+F_{\min}}{2}=0.25F_{\max} $$

$$ F_a=\frac{F_{\max}-F_{\min}}{2}=0.75F_{\max} $$

周期加载可写为：

$$ F(t)=F_{\mathrm{mean}}+F_a\sin\left(\frac{2\pi t}{T}\right) $$

正式周期开始前仍采用显式预加载：

```text
0 → Fmean
```

随后周期历史从 `Fmean` 开始：

```text
Fmean → +Fmax → Fmean → Fmin → Fmean
```

其中：

```text
Fmax  = +1.0 MN
Fmin  = -0.5 MN
Fmean = +0.25 MN
```

预加载不计入循环数。

---

# 2. 机制对照设计

100 圈分析继续采用 paired comparison：

```text
Case A: coupled viscoplastic-damage
Case B: damage-disabled, k_damage = 0
```

两组分析保持：

- 相同塔筒几何；
- 相同有限元网格；
- 相同 Gauss 积分；
- 相同纤维离散；
- 相同加载历史；
- 相同塑性与硬化参数；
- 相同时间离散；
- 相同 Newton 求解参数；
- 使用同一个由 coupled case 确定的临界纤维位置。

因此，二者差异主要用于识别：

$$ \text{damage coupling contribution} $$

而不是由离散方式、临界点选取或加载差异造成。

---

# 3. 数值稳定性与临界位置

100 圈计算中临界位置仍保持为：

```text
critical location = (0, 0, 4)
critical height   = 1.851206 m
critical y        = 2.945089 m
```

该位置仍位于塔底附近的外缘纤维。

这意味着从前期到第 100 圈：

> 危险位置没有发生跳变，因此当前长期演化趋势不是由“临界纤维重新选择”导致的人为不连续。

最大 Newton 迭代次数：

```text
coupled        = 4
damage-disabled = 4
```

最大残差范数：

```text
coupled        = 8.761245e-03
damage-disabled = 9.711299e-03
```

因此，100 圈内没有出现明显的非线性求解失稳或迭代恶化。

这为后续对损伤转折和塑性路径演化进行物理解读提供了必要前提。

---

# 4. 100 圈关键结果总览

选定循环结果如下：

<table>
<thead>
<tr>
<th>Cycle</th>
<th>du_c (m)</th>
<th>du_0 (m)</th>
<th>dep_c</th>
<th>dep_0</th>
<th>D_end</th>
<th>dD</th>
<th>net/path_c</th>
<th>net/path_0</th>
</tr>
</thead>
<tbody>
<tr><td>1</td><td>6.144753e-02</td><td>6.121915e-02</td><td>1.387954e-04</td><td>1.382840e-04</td><td>1.964039e-03</td><td>1.964039e-03</td><td>1.000000</td><td>1.000000</td></tr>
<tr><td>10</td><td>1.876899e-02</td><td>1.814947e-02</td><td>2.585422e-05</td><td>2.467798e-05</td><td>8.634555e-03</td><td>5.730438e-04</td><td>1.000000</td><td>1.000000</td></tr>
<tr><td>20</td><td>1.226492e-02</td><td>1.161932e-02</td><td>1.633264e-05</td><td>1.519911e-05</td><td>1.363268e-02</td><td>4.520716e-04</td><td>0.997853</td><td>0.999994</td></tr>
<tr><td>30</td><td>8.982901e-03</td><td>8.378241e-03</td><td>1.156499e-05</td><td>1.062832e-05</td><td>1.775694e-02</td><td>3.855591e-04</td><td>0.956183</td><td>0.994194</td></tr>
<tr><td>40</td><td>7.034037e-03</td><td>6.446353e-03</td><td>8.838646e-06</td><td>7.953675e-06</td><td>2.139706e-02</td><td>3.512231e-04</td><td>0.856252</td><td>0.965085</td></tr>
<tr><td>50</td><td>5.940139e-03</td><td>5.249792e-03</td><td>7.752436e-06</td><td>6.396428e-06</td><td>2.487349e-02</td><td>3.484715e-04</td><td>0.767281</td><td>0.916519</td></tr>
<tr><td>60</td><td>5.336072e-03</td><td>4.537272e-03</td><td>7.365140e-06</td><td>5.652822e-06</td><td>2.841615e-02</td><td>3.598587e-04</td><td>0.702921</td><td>0.875313</td></tr>
<tr><td>70</td><td>4.956335e-03</td><td>4.076439e-03</td><td>7.119338e-06</td><td>5.252806e-06</td><td>3.209328e-02</td><td>3.743944e-04</td><td>0.644650</td><td>0.844198</td></tr>
<tr><td>80</td><td>4.689173e-03</td><td>3.744367e-03</td><td>6.925928e-06</td><td>4.960666e-06</td><td>3.592401e-02</td><td>3.903628e-04</td><td>0.589714</td><td>0.816471</td></tr>
<tr><td>90</td><td>4.489398e-03</td><td>3.483097e-03</td><td>6.772613e-06</td><td>4.711866e-06</td><td>3.992108e-02</td><td>4.075197e-04</td><td>0.538700</td><td>0.789642</td></tr>
<tr><td>100</td><td>4.334027e-03</td><td>3.267305e-03</td><td>6.654395e-06</td><td>4.489670e-06</td><td>4.409565e-02</td><td>4.257387e-04</td><td>0.492039</td><td>0.763286</td></tr>
</tbody>
</table>

其中：

- `du_c`：coupled case 的每循环位移漂移；
- `du_0`：damage-disabled case 的每循环位移漂移；
- `dep_c`：coupled case 的每循环净塑性应变漂移；
- `dep_0`：damage-disabled case 的每循环净塑性应变漂移；
- `D_end`：循环结束时的最大损伤；
- `dD`：每循环最大损伤增量；
- `net/path`：净塑性增量与总塑性路径长度之比。

---

# 5. 最重要的结论之一：第 46 圈确实是损伤增量转折点

程序自动检查得到：

```text
dD minimum at cycle 46: 3.466680e-04
cycle-100 / minimum = 1.228088
```

因此：

$$ \Delta D_{46}=3.466680\times10^{-4} $$

是前 100 圈范围内的最小值。

第 100 圈：

$$ \Delta D_{100}=4.257387\times10^{-4} $$

二者相比：

$$ \frac{\Delta D_{100}}{\Delta D_{46}}=1.228088 $$

即：

$$ \frac{\Delta D_{100}-\Delta D_{46}}{\Delta D_{46}}\approx22.81\% $$

因此，50 圈时观察到的轻微回升已经可以被确认不是偶然波动。

本阶段可以较为明确地写出：

> **每循环损伤增量在约第 46 圈达到最低值，此后进入持续增长阶段。**

需要注意，这里的结论是：

```text
damage-rate acceleration
```

即：

> 损伤增长率开始提高。

而不是：

```text
damage instability
```

或：

```text
global structural instability
```

因为当前：

$$ D_{100}=0.04409565 $$

仍然处于相对较低的损伤水平，且 Newton 求解仍保持稳定。

---

# 6. 第 46 圈之后的损伤增长不是短暂反弹，而是持续趋势

100 圈尾部趋势进一步证明这一点。

81–90 圈：

$$ \text{slope}(\Delta D)=1.721491\times10^{-6} $$

且：

$$ \Delta D_{81\rightarrow90}=+3.952\% $$

91–100 圈：

$$ \text{slope}(\Delta D)=1.827106\times10^{-6} $$

且：

$$ \Delta D_{91\rightarrow100}=+4.017\% $$

两段斜率均为正：

$$ \text{slope}(\Delta D)>0 $$

并且第二个窗口的增长斜率略高于前一个窗口：

$$ 1.721491\times10^{-6}\rightarrow1.827106\times10^{-6} $$

因此可以形成比 50 圈时更明确的判断：

> **系统已经进入持续的损伤增长率提升阶段。**

这说明塑性稳定化机制虽然仍存在，但已经无法继续使每循环损伤增量保持下降。

---

# 7. 第二个核心现象：损伤速率反转，但宏观净漂移尚未反转

这一点是本阶段最值得重视的物理现象。

第 100 圈程序检查得到：

```text
du_c  minimum at cycle 100
du_0  minimum at cycle 100
dep_c minimum at cycle 100
dep_0 minimum at cycle 100
```

因此截至第 100 圈：

$$ \Delta u_n^c\downarrow $$

$$ \Delta u_n^0\downarrow $$

$$ \Delta\varepsilon_{p,n}^c\downarrow $$

$$ \Delta\varepsilon_{p,n}^0\downarrow $$

但与此同时：

$$ \Delta D_n\uparrow $$

因此三个关键演化指标已经发生趋势分离：

$$ \boxed{\Delta u_n\downarrow,\quad \Delta\varepsilon_{p,n}\downarrow,\quad \Delta D_n\uparrow} $$

这一现象说明系统当前不能简单归类为：

```text
fully stabilized
```

也不能简单归类为：

```text
global ratcheting acceleration
```

更准确的描述为：

> **净塑性漂移仍在趋缓，但损伤增长率已经转入加速阶段。**

即：

$$ \boxed{\text{plastic stabilization}+\text{damage-rate acceleration}} $$

---

# 8. 为什么“净漂移继续下降”不等于已经 shakedown

严格的循环 shakedown 通常意味着：

> 循环经过足够调整后，塑性增量趋近于零，后续响应逐渐接近稳定的弹性循环或严格稳定循环。

但当前 coupled case 中：

$$ \Delta\varepsilon_{p,100}^{c}=6.654395\times10^{-6} $$

仍然非零。

更重要的是：

$$ r_{p,100}^{c}=0.492039 $$

远低于 1。

因此，一个循环内部仍存在显著的塑性活动。

当前看到的是：

```text
net plastic drift ↓
```

而不是：

```text
plastic activity → 0
```

二者不能等同。

所以目前更适合使用：

```text
toward stabilization
```

或：

```text
net-drift stabilization
```

而不宜直接写成：

```text
plastic shakedown
```

---

# 9. `net/path` 的长期下降揭示了循环内塑性路径的重构

定义：

$$ \Delta\varepsilon_p=\varepsilon_{p,\mathrm{end}}-\varepsilon_{p,\mathrm{start}} $$

表示一圈结束后留下的净塑性应变漂移。

定义塑性路径长度：

$$ L_p=\sum_i|\delta\varepsilon_{p,i}| $$

其中：

$$ \delta\varepsilon_{p,i} $$

为循环内部第 $i$ 个增量步的塑性应变变化。

则：

$$ r_p=\frac{|\Delta\varepsilon_p|}{L_p} $$

即程序中的：

```text
net/path
```

当：

$$ r_p\approx1 $$

意味着：

> 一圈内部的塑性活动几乎全部沿同一方向发展。

当：

$$ r_p<1 $$

意味着：

> 一圈内部已经同时存在正向与反向塑性增量，部分塑性路径在循环结束时相互抵消。

---

# 10. 从第 20 圈到第 100 圈：plastic path 发生了根本变化

coupled case：

$$ r_{p,20}^{c}=0.997853 $$

$$ r_{p,50}^{c}=0.767281 $$

$$ r_{p,100}^{c}=0.492039 $$

因此：

```text
cycle 20  : near-unidirectional plastic accumulation
cycle 50  : clear reverse plastic activity
cycle 100 : strong forward + reverse cyclic plastic activity
```

也就是说：

$$ \boxed{\text{近乎纯单向塑性累积}\rightarrow\text{显著循环内正反向塑性活动}} $$

这并不意味着净棘轮方向已经反转。

因为第 100 圈仍满足：

$$ \Delta\varepsilon_{p,100}^{c}>0 $$

所以总体上仍存在正方向净塑性累积。

真正发生变化的是：

> **单圈塑性应变路径本身已经由近单调演化转变为明显的正反向往复演化。**

---

# 11. 第 100 圈 coupled case 的塑性路径定量分析

第 100 圈：

$$ \Delta\varepsilon_p^{c}=6.654395\times10^{-6} $$

且：

$$ r_p^{c}=0.492039 $$

因此：

$$ L_p^c=\frac{|\Delta\varepsilon_p^c|}{r_p^c} $$

得到：

$$ L_p^c\approx1.3524\times10^{-5} $$

即：

> 第 100 圈 coupled case 中，临界纤维实际经历的总塑性路径约为 $1.35\times10^{-5}$，明显高于最终留下的净塑性增量 $6.65\times10^{-6}$。

这说明：

> 单独观察净塑性漂移会低估循环内部真实的塑性活动强度。

---

# 12. 正向与反向塑性路径的进一步分解

设：

$$ P^+ $$

为一圈内部所有正向塑性增量之和，

$$ P^- $$

为所有反向塑性增量绝对值之和。

则：

$$ L_p=P^++P^- $$

在净方向仍为正时：

$$ \Delta\varepsilon_p=P^+-P^- $$

因此：

$$ r_p=\frac{P^+-P^-}{P^++P^-} $$

反向塑性路径占总塑性路径的比例为：

$$ \frac{P^-}{P^++P^-}=\frac{1-r_p}{2} $$

第 100 圈 coupled case：

$$ r_p^c=0.492039 $$

因此反向塑性路径比例约为：

$$ \frac{1-0.492039}{2}\approx25.4\% $$

正向塑性路径比例约为：

$$ 74.6\% $$

因此可以直观理解为：

```text
cycle 100 coupled plastic path:

forward plastic activity ≈ 74.6%
reverse plastic activity ≈ 25.4%
```

这已经远离“近乎纯单向塑性累积”的早期状态。

---

# 13. damage-disabled case 同样出现反向塑性，但明显弱于 coupled case

第 100 圈 damage-disabled：

$$ r_{p,100}^{0}=0.763286 $$

其反向塑性路径比例约为：

$$ \frac{1-0.763286}{2}\approx11.8\% $$

因此：

```text
cycle 100 reverse plastic path fraction:

damage-disabled ≈ 11.8%
coupled         ≈ 25.4%
```

这一结果揭示了一个重要的因果关系：

> **反向塑性活动并不是损伤单独制造出来的。**

即使：

$$ D=0 $$

非对称循环、粘塑性演化、随动硬化和 Bauschinger 型效应也会使循环内部逐步出现反向塑性。

但是：

$$ 25.4\%>11.8\% $$

说明损伤显著放大了这种循环内塑性重构。

因此更合理的机制链为：

$$ \boxed{\text{非对称循环}+\text{塑性演化}\rightarrow\text{基础循环塑性重组}} $$

随后：

$$ \boxed{D\uparrow\rightarrow K_{\mathrm{eff}}\downarrow\rightarrow\text{循环内塑性活动进一步放大}} $$

---

# 14. 第 100 圈损伤放大效应显著增强

第 100 圈 paired comparison 得到：

```text
du amplification       = 32.648357%
dep amplification      = 48.215660%
plastic-path ratio c/0 = 2.299226
u-range change         = 0.619775%
stress-range change    = -1.097607%
work change            = 39.129716%
```

因此损伤作用已经不再是早期的小幅修正。

## 14.1 位移漂移放大

$$ \frac{\Delta u_{100}^{c}-\Delta u_{100}^{0}}{\Delta u_{100}^{0}}=32.65\% $$

## 14.2 净塑性漂移放大

$$ \frac{\Delta\varepsilon_{p,100}^{c}-\Delta\varepsilon_{p,100}^{0}}{\Delta\varepsilon_{p,100}^{0}}=48.22\% $$

## 14.3 循环塑性路径放大

$$ \frac{L_{p,100}^{c}}{L_{p,100}^{0}}=2.299 $$

即：

> coupled case 一圈内部实际经历的塑性路径已经达到 damage-disabled case 的约 2.30 倍。

这一指标比单独比较净塑性漂移更能揭示损伤对循环塑性活动的真实放大作用。

---

# 15. 第 50 圈与第 100 圈损伤放大效应比较

第 50 圈：

```text
du amplification   ≈ 13.15%
dep amplification  ≈ 21.20%
work amplification ≈ 14.79%
```

第 100 圈：

```text
du amplification   ≈ 32.65%
dep amplification  ≈ 48.22%
work amplification ≈ 39.13%
```

因此：

<table>
<thead>
<tr>
<th>指标</th>
<th>Cycle 50</th>
<th>Cycle 100</th>
</tr>
</thead>
<tbody>
<tr><td>位移漂移放大</td><td>13.15%</td><td>32.65%</td></tr>
<tr><td>净塑性漂移放大</td><td>21.20%</td><td>48.22%</td></tr>
<tr><td>循环外功放大</td><td>14.79%</td><td>39.13%</td></tr>
</tbody>
</table>

说明：

$$ \boxed{\text{damage coupling effect increases strongly with cycle number}} $$

也就是说，损伤的影响具有明显的后期累积放大特征。

---

# 16. 宏观响应变化仍然较小，但内部变量已经明显分叉

第 100 圈：

位移范围变化：

$$ +0.619775\% $$

应力范围变化：

$$ -1.097607\% $$

这两个宏观循环范围变化仍然较小。

但是：

$$ \Delta\varepsilon_p $$

已经放大约：

$$ 48.22\% $$

塑性路径长度已经放大约：

$$ 2.30\ \text{times} $$

循环外功增加约：

$$ 39.13\% $$

因此当前损伤演化首先表现为：

```text
internal-variable evolution
```

而不是：

```text
sudden global structural response amplification
```

这意味着：

> 宏观位移和应力循环范围尚未出现剧烈退化，但临界材料点内部的塑性与损伤状态已经明显脱离 damage-disabled 基准。

这一现象对于高周疲劳降阶尤其重要，因为：

> 真正决定后续退化的往往首先是内部变量场，而不是单纯的全局位移响应。

---

# 17. 100 圈累计漂移结果

截至第 100 圈：

```text
sum du_c  = 9.346771e-01 m
sum du_0  = 8.595774e-01 m

sum dep_c = 1.363490e-03
sum dep_0 = 1.217583e-03
```

累计位移漂移差异约为：

$$ \frac{0.9346771-0.8595774}{0.8595774}\approx8.74\% $$

累计塑性漂移差异约为：

$$ \frac{1.363490\times10^{-3}-1.217583\times10^{-3}}{1.217583\times10^{-3}}\approx11.98\% $$

需要注意：

> 第 100 圈瞬时单圈差异明显大于前 100 圈累计差异。

例如：

```text
cycle-100 du amplification  ≈ 32.65%
100-cycle cumulative du     ≈  8.74%

cycle-100 dep amplification ≈ 48.22%
100-cycle cumulative dep    ≈ 11.98%
```

这进一步证明：

$$ \boxed{\text{damage effect is strongly late-cycle amplified}} $$

损伤影响并不是从第一圈开始保持固定比例，而是随循环逐渐增强。

---

# 18. 81–100 圈：coupled 与 damage-disabled 的趋稳速度进一步分离

81–90 圈：

```text
du_c  change = -3.798%
du_0  change = -6.259%

dep_c change = -1.973%
dep_0 change = -4.508%

dD    change = +3.952%
```

91–100 圈：

```text
du_c  change = -3.088%
du_0  change = -5.562%

dep_c change = -1.552%
dep_0 change = -4.243%

dD    change = +4.017%
```

因此：

> 两组模型的净漂移都在下降，但 coupled case 的下降速度显著慢于 damage-disabled case。

这意味着：

$$ D\uparrow $$

正在持续抵消塑性硬化与循环调整带来的稳定化趋势。

因此 coupled case 可以理解为：

$$ \boxed{\text{stabilization tendency weakened by cumulative damage}} $$

---

# 19. 当前 0–100 圈可以划分为三个演化阶段

基于现有数据，可以形成一个初步、但具有明确物理依据的阶段划分。

---

## Stage I：早期快速调整阶段

大致范围：

$$ n=1\sim20 $$

主要特征：

$$ \Delta u_n\downarrow\downarrow $$

$$ \Delta\varepsilon_{p,n}\downarrow\downarrow $$

$$ \Delta D_n\downarrow $$

同时：

$$ r_p\approx1 $$

说明：

> 初始循环中存在明显的单向塑性累积，随后由于硬化和内部变量调整，每循环净漂移快速下降。

这一阶段的主导特征是：

```text
early ratcheting adjustment
+
near-unidirectional plastic evolution
```

---

## Stage II：稳定化—损伤竞争过渡阶段

大致范围：

$$ n=20\sim46 $$

主要特征：

$$ \Delta u_n\downarrow $$

$$ \Delta\varepsilon_{p,n}\downarrow $$

$$ \Delta D_n\downarrow\quad\text{but flattening} $$

同时：

$$ r_p<1 $$

并持续下降。

说明：

> 净漂移仍在趋缓，但循环内部反向塑性活动开始增强，损伤增量的下降趋势越来越弱。

此时两个机制开始明显竞争：

$$ \boxed{\text{plastic stabilization}\leftrightarrow\text{damage softening}} $$

---

## Stage III：损伤增长率加速阶段

大致范围：

$$ n>46 $$

主要特征：

$$ \Delta D_n\uparrow $$

但同时：

$$ \Delta u_n\downarrow $$

$$ \Delta\varepsilon_{p,n}\downarrow $$

以及：

$$ r_p\downarrow $$

因此该阶段不能称为：

```text
global ratcheting acceleration
```

而更准确地描述为：

```text
damage-rate acceleration
with continuing net-drift stabilization
```

即：

> 损伤增长率已经反转并持续提高，但净位移漂移和净塑性漂移尚未出现反转。

---

# 20. 当前系统还没有进入“整体损伤失稳”

尽管：

$$ \Delta D_n\uparrow $$

但以下现象均未发生：

- `du_c` 尚未出现最低点后回升；
- `dep_c` 尚未出现最低点后回升；
- 位移循环范围尚未急剧增加；
- Newton 迭代尚未恶化；
- 全局响应尚未出现跳跃或发散；
- $D$ 尚未接近上限。

因此目前应避免使用：

```text
damage instability
structural instability
failure acceleration
```

更合适的术语为：

```text
damage-rate acceleration
damage-growth transition
damage-softening-dominated transition tendency
```

---

# 21. 为什么本阶段对 LATIN-PGD 特别重要

此前一个较为简单的 slow-fast 思路可以写为：

$$ u(x,n,\tau)\approx\bar u(x,n)+\widetilde u(x,\tau) $$

其中：

- $x$：空间变量；
- $n$：循环编号，慢尺度；
- $\tau\in[0,T]$：单循环内部快尺度时间。

这种写法隐含一个重要假设：

> 单个循环内部的快尺度波形 $\widetilde u(x,\tau)$ 基本稳定，仅慢尺度平均状态随循环数发生漂移。

但 100 圈结果表明：

$$ r_p^c:1.000\rightarrow0.9979\rightarrow0.7673\rightarrow0.4920 $$

循环内部的塑性路径形状本身在持续发生变化。

因此更合理的认识应为：

$$ \boxed{\widetilde u=\widetilde u(x,n,\tau)} $$

同样，对于塑性应变：

$$ \boxed{\widetilde\varepsilon_p=\widetilde\varepsilon_p(x,n,\tau)} $$

对于损伤：

$$ \boxed{D=D(x,n,\tau)} $$

也就是说：

> **快尺度周期响应不能再被视为与慢尺度循环数完全无关的固定模板。**

---

# 22. 对 PGD 分离形式的直接影响

更一般的 PGD 表达可考虑：

$$ q(x,n,\tau)\approx\sum_{k=1}^{r}X_k(x)N_k(n)T_k(\tau) $$

其中：

- $X_k(x)$：空间模态；
- $N_k(n)$：慢循环尺度函数；
- $T_k(\tau)$：单循环快尺度函数；
- $r$：所需分离秩。

100 圈结果提出了两个必须实际验证的问题。

## 22.1 固定 rank 是否足够

由于第 1 圈、第 20 圈、第 50 圈和第 100 圈的循环内塑性路径明显不同：

> 一个固定且很低的 rank 是否能够同时表示三个不同演化阶段，需要通过 FOM 快照的低秩性分析验证。

不能预先假定：

$$ r=\text{constant and small} $$

---

## 22.2 是否需要阶段化基函数

当前 0–100 圈已经表现出至少三个不同阶段：

```text
Stage I   early adaptation
Stage II  stabilization-damage competition
Stage III damage-rate acceleration
```

因此后续很可能需要考虑：

```text
stagewise PGD basis
```

或：

```text
adaptive enrichment
```

即：

$$ r=r(n) $$

而不是从第 1 圈到高周阶段始终使用完全固定的模态集合。

---

# 23. 为什么 `net/path` 是 LATIN-PGD 降阶中的重要诊断量

如果只观察：

$$ \Delta\varepsilon_p $$

会看到它从早期到第 100 圈持续下降。

这可能误导我们认为：

> 循环塑性正在逐渐消失，快尺度响应将越来越简单。

但是：

$$ r_p $$

却从接近 1 持续下降到约 0.492。

这说明：

> 单圈内部的正反向塑性活动正在增强。

因此真正的快尺度复杂度并没有简单随净漂移下降而减弱。

这意味着后续评估 LATIN-PGD 的低秩性时，不能仅看：

```text
cycle-end drift
```

而必须保留：

```text
full fast-time trajectory within each cycle
```

否则可能错误高估时间可分离性。

---

# 24. 后续 FOM 快照应从“标量指标”升级到“全场状态”

截至本阶段，主要诊断量仍集中于临界纤维：

- 每圈位移漂移；
- 每圈净塑性漂移；
- 损伤增量；
- `net/path`；
- 位移范围；
- 应力范围；
- 循环外功。

这些指标足以完成机制识别。

但要真正进入 LATIN-PGD 降阶阶段，需要开始构造完整 FOM 快照。

推荐至少提取：

```text
element
× Gauss point
× fiber
× cycle n
× fast phase τ
```

对应的状态量至少包括：

$$ u $$

$$ \sigma $$

$$ \varepsilon_p $$

$$ D $$

以及与当前材料模型直接相关的硬化内部变量。

这样才能从：

```text
critical-point scalar diagnosis
```

转向：

```text
full-field low-rank structure assessment
```

---

# 25. 下一阶段真正需要回答的核心问题

下一阶段不应首先继续盲目增加 FOM 循环数，而应优先回答：

> **当前 100 圈全阶响应在空间—慢循环—快时间三个维度上到底具有怎样的低秩结构？**

具体应研究：

$$ \text{rank of }u $$

$$ \text{rank of }\sigma $$

$$ \text{rank of }\varepsilon_p $$

$$ \text{rank of }D $$

并比较：

```text
Stage I
Stage II
Stage III
```

之间所需秩是否发生变化。

这将直接决定后续 LATIN-PGD 应采用：

```text
fixed basis
```

还是：

```text
stagewise basis
```

或：

```text
adaptive enrichment
```

---

# 26. 建议的下一阶段技术路线

当前最自然的研究顺序为：

```text
Step 1
保留 100-cycle FOM benchmark

Step 2
增加全场内部变量 snapshot extraction

Step 3
按 cycle n × fast phase τ 重组时间维度

Step 4
构造位移、应力、塑性应变、损伤快照张量

Step 5
先做 SVD / low-rank diagnostic

Step 6
比较不同演化阶段的 singular-value decay

Step 7
判断固定 rank 是否足够

Step 8
再设计 LATIN-PGD slow-fast separated representation

Step 9
根据低秩分析结果决定是否需要 adaptive enrichment
```

这一顺序的意义在于：

> 先用 FOM 数据证明问题具有怎样的低秩结构，再决定 PGD 的具体分离形式，而不是先人为假设某种低秩形式，再强行把非线性疲劳响应塞进去。

---

# 27. 当前暂时不需要引入真实风浪载荷谱

本阶段仍采用规则非对称循环：

$$ R_F=-0.5 $$

目的不是模拟完整真实海上服役荷载，而是建立一个：

- 可控；
- 可重复；
- 可解释；
- 具有明显非对称循环效应；
- 能够激发塑性、棘轮与损伤竞争；
- 适合验证 LATIN-PGD 长循环降阶机制

的参考基准问题。

因此，当前阶段继续保持：

```text
regular asymmetric cyclic loading
```

是合理的。

真实风浪载荷谱应在 LATIN-PGD 对规则长循环问题的基本降阶框架验证完成后再进入。

---

# 28. 外功指标的使用注意事项

程序给出的循环外功变化为：

$$ +39.129716\% $$

当前仍应称为：

```text
cycle external work
```

或：

```text
path-related external work
```

而不宜直接称为：

```text
pure material dissipation
```

因为循环尚未完全稳定时：

$$ W_{\mathrm{ext}}=\oint F\,du $$

中除了塑性与损伤相关耗散外，还可能包含：

$$ \Delta U_{\mathrm{stored}} $$

即循环始末储存能变化。

因此正式论文中应避免将当前循环外功与纯耗散能完全等同。

---

# 29. 当前阶段可以形成的机制链

综合 0–100 圈结果，当前最合理的机制链可以写为：

$$ \boxed{\text{asymmetric cyclic loading}} $$

$$ \Downarrow $$

$$ \boxed{\text{persistent net plastic drift}} $$

$$ \Downarrow $$

$$ \boxed{\text{kinematic / viscoplastic cyclic adjustment}} $$

同时：

$$ \boxed{\text{reverse plastic activity gradually develops}} $$

进一步：

$$ \boxed{D\uparrow} $$

导致：

$$ \boxed{K_{\mathrm{eff}}\downarrow} $$

进而：

$$ \boxed{\text{cyclic plastic path amplified}} $$

再进一步：

$$ \boxed{\Delta D_n\text{ reaches a minimum near cycle 46}} $$

随后：

$$ \boxed{\Delta D_n\uparrow} $$

但当前仍保持：

$$ \boxed{\Delta u_n\downarrow,\quad \Delta\varepsilon_{p,n}\downarrow} $$

因此系统目前处于：

$$ \boxed{\text{net-drift stabilization}+\text{damage-rate acceleration}} $$

而不是整体失稳。

---

# 30. 本阶段对 LATIN-PGD 的最终方法学认识

100 圈全阶分析给出的最重要认识不是单纯“损伤增加了多少”，而是：

> **不同时间尺度上的响应演化并不同步。**

慢尺度上：

$$ D_n $$

持续累积，并在约第 46 圈后表现出增量加速。

快尺度上：

$$ \varepsilon_p(\tau) $$

的循环内部路径持续发生形状变化。

宏观净漂移上：

$$ \Delta u_n $$

和：

$$ \Delta\varepsilon_{p,n} $$

仍然下降。

因此后续 LATIN-PGD 不宜把时间维度简单视为：

```text
fixed periodic fast mode
+
single slowly drifting amplitude
```

更适合考虑：

$$ \boxed{\text{space }x\times\text{slow cycle }n\times\text{fast phase }\tau} $$

并允许：

```text
fast-cycle pattern evolves with n
```

这自然导向：

```text
slow-fast tensorization
+
stagewise basis
+
adaptive enrichment
```

---

# 31. 本阶段最终结论

本次 100 圈非对称循环分析可以形成以下正式结论。

第一，50 圈阶段观察到的第 46 圈附近损伤增量最低点已被 100 圈结果确认。

$$ \boxed{\Delta D_n\text{ reaches its minimum near cycle }46} $$

第二，第 46 圈之后：

$$ \boxed{\Delta D_n\uparrow} $$

并在 81–90、91–100 两个独立窗口中均表现出持续正斜率，说明系统已经进入稳定的损伤增长率提高阶段。

第三，尽管损伤增长率已经反转：

$$ \Delta u_n $$

和：

$$ \Delta\varepsilon_{p,n} $$

截至第 100 圈仍在下降。

因此当前不是宏观棘轮加速或结构失稳，而是：

$$ \boxed{\text{damage-rate acceleration with continuing net-drift stabilization}} $$

第四，`net/path` 从早期接近 1 下降至：

$$ r_{p,100}^{c}=0.492039 $$

证明循环内部塑性路径已经由近乎纯单向累积转变为显著的正反向塑性活动。

第五，第 100 圈 coupled case 的塑性路径长度约为 damage-disabled case 的：

$$ \boxed{2.299\ \text{times}} $$

说明损伤对循环内部塑性活动的放大远强于单纯从净塑性漂移所看到的差异。

第六，当前 FOM 结果已经说明：

> 快尺度循环波形本身会随着慢尺度循环数和损伤状态演化。

因此后续 LATIN-PGD 不能简单假定：

$$ \widetilde q=\widetilde q(x,\tau) $$

而应允许：

$$ \widetilde q=\widetilde q(x,n,\tau) $$

或采用等价的多尺度可分离表示。

第七，当前研究重点可以从：

```text
long-cycle scalar mechanism diagnosis
```

正式转向：

```text
full-field snapshot construction
+
low-rank diagnostics
+
slow-fast LATIN-PGD formulation
```

---

# 32. 当前阶段一句话总结

> **100 圈非对称循环结果确认，第 46 圈附近存在稳定的每循环损伤增量转折，此后损伤增长率持续提高，而净位移和净塑性漂移仍继续下降，说明系统已进入“塑性稳定化与损伤软化竞争、损伤速率增强但宏观净漂移尚未反转”的长期过渡阶段；与此同时，`net/path` 从接近 1 降至约 0.492，证明循环内部快尺度塑性路径正在显著演化，从而为后续 LATIN-PGD 采用 slow-fast 时间表示、阶段化基函数以及自适应模态增广提供了直接的全阶依据。**

---

# 33. 下一阶段入口

下一阶段建议正式定义为：

```text
FOM full-field snapshot construction
and low-rank structure diagnosis
```

首先不修改材料模型和加载方式，而是在当前 100 圈基准问题上扩展状态场输出。

目标是得到：

```text
element × Gauss × fiber × cycle × fast phase
```

形式的全阶数据，并据此研究：

- 位移场低秩性；
- 应力场低秩性；
- 塑性应变场低秩性；
- 损伤场低秩性；
- 不同循环阶段的秩变化；
- fixed basis 与 adaptive basis 的适用边界。

这将作为真正进入 LATIN-PGD 塔筒高周疲劳降阶实现之前的直接数值依据。
