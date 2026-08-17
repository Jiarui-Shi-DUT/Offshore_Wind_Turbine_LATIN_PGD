# Tower LATIN-PGD Eq. (72) 初始时间处理、BE temporal contraction 与 denominator conditioning 阶段总结

**日期：2026-08-17**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**当前研究路线：Bhattacharyya et al. 原论文 $x-t$ LATIN-PGD → 2D fiber beam-column offshore wind turbine tower**  
**阶段范围：自上一份 Eq. (72) backward-Euler temporal update 阶段总结之后，完成 $t_0$ 初始 temporal treatment、$\lambda_0$ 与 $\dot{\lambda}_0$ 的区分和推导、Eq. (70)–(71) temporal contraction 的 backward-Euler 一致性分析、fixed-point 第一轮初始化逻辑、Eq. (72) scalar denominator 的数学含义、尺度无关 controllability 指标以及不修改原方程的 safeguard 原则**  
**上一阶段衔接：`2026-08-17-tower-latin-pgd-eq72-backward-euler-scalar-update-stage-summary.md`**  
**下一阶段：明确 Eq. (72) temporal update 返回 Eq. (70)–(71) 后 enrichment fixed-point 的完整循环顺序，并确定 new separated pair 的 fixed-point convergence criterion**

---

# 1. 本阶段定位

上一阶段已经完成 single-new-mode Eq. (72) 的 sequential backward-Euler temporal update。

对于 $n \ge 1$，已有：

$$ \vec{r}_n = \vec{g}_n \lambda_n - \vec{b}_n $$

其中：

$$ \vec{g}_n = \frac{\vec{p}}{\Delta t_n} - D_{H,n}\vec{s} $$

以及：

$$ \vec{b}_n = \frac{\vec{p}}{\Delta t_n}\lambda_{n-1} - \vec{\Delta}_n $$

weighted objective：

$$ J_n = \frac{1}{2}\vec{r}_n^T M D_{H,n}^{-1}\vec{r}_n $$

single-new-mode scalar update：

$$ \lambda_n = \frac{\vec{g}_n^T M D_{H,n}^{-1}\left( \frac{\vec{p}}{\Delta t_n}\lambda_{n-1} - \vec{\Delta}_n \right)}{\vec{g}_n^T M D_{H,n}^{-1}\vec{g}_n} $$

上一阶段同时确认，该公式是 sequential one-step weighted least-squares 的精确解，而不是 whole-time Eq. (72) minimisation 的 global exact solution。

但上一阶段仍留下三个关键问题：

1. 第一个时间点 $t_0$ 没有上一节点，不能直接使用 backward Euler；
2. Eq. (70)–(71) 中的 temporal contraction $a=\langle \dot{\lambda}\lambda \rangle$ 如何与 current backward-Euler temporal solver 保持离散一致；
3. Eq. (72) scalar denominator 接近零时，应如何识别和处理，而不篡改原论文 residual minimisation 结构。

本阶段即围绕这三个问题展开。

---

# 2. 本阶段核心结论

本阶段最终形成以下三个主要结论。

第一，对于 virgin-state tower v1，可采用：

$$ \lambda_0 = 0 $$

但：

$$ \dot{\lambda}_0 \ne 0 $$

是完全允许的。

$t_0$ 的 initial rate 通过独立的 $H_\sigma^{-1}$ weighted least-squares 求得：

$$ \dot{\lambda}_0 = -\frac{\vec{p}^T M D_{H,0}^{-1}\vec{\Delta}_0}{\vec{p}^T M D_{H,0}^{-1}\vec{p}} $$

第二，若 tower v1 采用 right-endpoint backward-Euler time-slab interpretation，则 Eq. (70)–(71) 中：

$$ a = \langle \dot{\lambda}\lambda \rangle $$

应采用离散量：

$$ a_h = \sum_{n=1}^{N}\Delta t_n \dot{\lambda}_n \lambda_n $$

即：

$$ a_h = \sum_{n=1}^{N}\lambda_n(\lambda_n-\lambda_{n-1}) $$

在 $\lambda_0=0$ 时：

$$ a_h = \frac{1}{2}\lambda_N^2 + \frac{1}{2}\sum_{n=1}^{N}(\lambda_n-\lambda_{n-1})^2 $$

因此只要 temporal mode 非平凡，就有：

$$ a_h > 0 $$

第三，Eq. (72) denominator：

$$ A_n = \vec{g}_n^T M D_{H,n}^{-1}\vec{g}_n $$

接近零的本质不是普通“除零错误”，而是当前 spatial candidate 在该 time step 对 LATIN residual 的 temporal controllability 发生退化。

因此 tower v1 不应采用：

$$ A_n \leftarrow \max(A_n,A_{\min}) $$

而应采用尺度无关的 controllability detection，并在真正退化时 reject 或 restart candidate pair。

---

# 3. 原论文 initial state 对 $t_0$ 的约束

原论文 LATIN initialisation 从 elastic solution 开始。

在该 initial elastic state 中，plastic strain 和其他 internal variables 均从初始状态开始。

对于当前 virgin-state tower v1，可采用：

$$ \varepsilon^p(x,0) = 0 $$

new PGD pair 的 plastic-strain correction 为：

$$ \Delta\varepsilon^p(x,t) = \lambda(t)\bar{\varepsilon}^p(x) $$

因此为了保持 initial plastic correction 为零：

$$ \Delta\varepsilon^p(x,0) = 0 $$

即：

$$ \lambda_0 \bar{\varepsilon}^p(x) = 0 $$

若 current spatial mode 非零：

$$ \bar{\varepsilon}^p(x) \ne 0 $$

则自然得到：

$$ \lambda_0 = 0 $$

---

# 4. $\lambda_0=0$ 的准确理论地位

需要区分两种说法。

不能说：

> 原论文 Eq. (72) 明确规定 $\lambda(0)=0$。

因为原论文没有在 Eq. (72) 后直接给出这一 temporal boundary equation。

更准确的说法是：

> 对 virgin-state problem，zero initial temporal amplitude 来源于 elastic LATIN initialisation 所要求的 zero initial plastic correction，再结合 separated plastic correction ansatz 推出。

因此 tower v1 可正式继承：

$$ \lambda_0 = 0 $$

但应把它表述为 initial correction admissibility / inherited initialisation condition，而不是 original Eq. (72) explicit boundary condition。

---

# 5. 为什么 $\lambda_0=0$ 不要求 $\dot{\lambda}_0=0$

temporal amplitude 与 temporal rate 是两个不同量。

例如：

$$ \lambda(t)=ct $$

则：

$$ \lambda(0)=0 $$

但：

$$ \dot{\lambda}(0)=c $$

可以非零。

因此：

$$ \lambda_0=0 $$

只表示 initial accumulated plastic correction 为零。

它并不要求：

$$ \dot{\lambda}_0=0 $$

在 LATIN-PGD 中，这意味着：

> new PGD pair 在 initial point 不改变 accumulated plastic state，但可以具有非零 initial plastic-rate correction。

---

# 6. current 1D `pgd_time_update.py` 的实际 first-point structure

当前 1D implementation 中：

```text
amplitudes = np.zeros((n_time, n_modes))
rates = np.zeros((n_time, n_modes))
```

因此：

$$ \lambda_j(t_0)=0 $$

对所有 temporal modes 均成立。

但代码随后单独求：

```text
rates[0, :]
```

并没有对 $t_0$ 使用 backward difference。

只有：

```text
for step in range(1, n_time):
```

以后才使用：

$$ \dot{\lambda}_n = \frac{\lambda_n-\lambda_{n-1}}{\Delta t_n} $$

因此 current 1D temporal treatment 的准确结构是：

$$ t_0:\quad \lambda_0=0,\quad \dot{\lambda}_0\ {\rm independently\ solved} $$

以及：

$$ n\ge1:\quad \dot{\lambda}_n=\frac{\lambda_n-\lambda_{n-1}}{\Delta t_n} $$

---

# 7. $t_0$ residual 的推导

continuous tower residual 为：

$$ \vec{r}(t)=\dot{\lambda}(t)\vec{p}-D_H(t)\lambda(t)\vec{s}+\vec{\Delta}(t) $$

在 $t_0$：

$$ \lambda_0=0 $$

因此：

$$ D_{H,0}\lambda_0\vec{s}=0 $$

于是：

$$ \vec{r}_0=\dot{\lambda}_0\vec{p}+\vec{\Delta}_0 $$

这说明 first point 上 current temporal correction 只能通过：

$$ \dot{\lambda}_0\vec{p} $$

调节 residual。

new stress correction 在 $t_0$ 为零。

---

# 8. $t_0$ weighted objective

$t_0$ 的 weighting operator 为：

$$ W_{H,0}=M D_{H,0}^{-1} $$

定义 initial objective：

$$ J_0(\dot{\lambda}_0)=\frac{1}{2}\left(\dot{\lambda}_0\vec{p}+\vec{\Delta}_0\right)^T M D_{H,0}^{-1}\left(\dot{\lambda}_0\vec{p}+\vec{\Delta}_0\right) $$

此时 unknown 不再是 $\lambda_0$，因为：

$$ \lambda_0=0 $$

已经固定。

唯一 unknown 为：

$$ \dot{\lambda}_0 $$

---

# 9. 对 $\dot{\lambda}_0$ 求导

由：

$$ \frac{\partial\vec{r}_0}{\partial\dot{\lambda}_0}=\vec{p} $$

得到：

$$ \frac{dJ_0}{d\dot{\lambda}_0}=\vec{p}^T M D_{H,0}^{-1}\left(\dot{\lambda}_0\vec{p}+\vec{\Delta}_0\right) $$

optimality condition：

$$ \vec{p}^T M D_{H,0}^{-1}\left(\dot{\lambda}_0\vec{p}+\vec{\Delta}_0\right)=0 $$

所以：

$$ \dot{\lambda}_0=-\frac{\vec{p}^T M D_{H,0}^{-1}\vec{\Delta}_0}{\vec{p}^T M D_{H,0}^{-1}\vec{p}} $$

这就是 current 1D first-point weighted least-squares solve 在 single-mode 情况下的数学对应。

---

# 10. 与 current code sign convention 的对应

current 1D code residual 写为：

$$ r=P\dot{\lambda}-H_\sigma S\lambda-f $$

而 Eq. (72) residual 写为：

$$ r=\dot{\lambda}p-H_\sigma\lambda s+\Delta $$

所以：

$$ f=-\Delta $$

因此 initial rate 也可写为：

$$ \dot{\lambda}_0=\frac{\vec{p}^T M D_{H,0}^{-1}\vec{f}_0}{\vec{p}^T M D_{H,0}^{-1}\vec{p}} $$

这与 current `pgd_time_update.py` 中：

```text
matrix=spatial_plastic
right_hand_side=forcing_array[0, :]
weights=element_volumes / directions.H_sigma[0, :]
```

完全对应。

---

# 11. tower $(e,g,f)$ 形式的 initial-rate solve

tower material-point index：

$$ q=(e,g,f) $$

其中：

$$ p_q=\bar{\varepsilon}^p_{egf} $$

$$ H_{\sigma,0,q}=H_{\sigma,egf}(t_0) $$

$$ \Delta_{0,q}=\bar{\Delta}_{egf}(t_0) $$

以及：

$$ v_q=A_{egf}w_gJ_e $$

因此：

$$ \dot{\lambda}_0=-\frac{\sum_{e,g,f}\frac{A_{egf}w_gJ_e}{H_{\sigma,egf}(t_0)}\bar{\varepsilon}^p_{egf}\bar{\Delta}_{egf}(t_0)}{\sum_{e,g,f}\frac{A_{egf}w_gJ_e}{H_{\sigma,egf}(t_0)}\left(\bar{\varepsilon}^p_{egf}\right)^2} $$

这就是 tower v1 single-new-mode first-point rate solve。

---

# 12. initial-rate formula 的意义

分母：

$$ \sum_{e,g,f}\frac{v_{egf}}{H_{\sigma,egf}(t_0)}\left(\bar{\varepsilon}^p_{egf}\right)^2 $$

表示 initial $H_\sigma^{-1}$ metric 下 plastic spatial mode 的 weighted magnitude。

分子：

$$ -\sum_{e,g,f}\frac{v_{egf}}{H_{\sigma,egf}(t_0)}\bar{\varepsilon}^p_{egf}\bar{\Delta}_{egf}(t_0) $$

表示 initial shifted LATIN defect 在该 spatial mode 上的 weighted projection。

因此 $\dot{\lambda}_0$ 的作用是：

> 在保持 accumulated initial plastic correction 为零的前提下，为 new spatial mode 选择最合适的 initial plastic-rate amplitude，使 $t_0$ mechanical residual 最小。

---

# 13. initial stress correction 与 compatible correction

new stress correction：

$$ \Delta\sigma'(t)=\lambda(t)\bar{\sigma} $$

所以：

$$ \Delta\sigma'_0=0 $$

同样，若 compatible strain correction 采用：

$$ \Delta\varepsilon'(t)=\lambda(t)\bar{\varepsilon}' $$

则：

$$ \Delta\varepsilon'_0=0 $$

因此：

> new PGD pair 在 $t_0$ 不改变 accumulated initial mechanical correction state，只允许 nonzero temporal rate 控制 initial evolution direction。

---

# 14. $t_0$ treatment 的最终结构

current tower v1 temporal solver 可分成两部分。

初始节点：

$$ \lambda_0=0 $$

$$ \dot{\lambda}_0=-\frac{\vec{p}^T M D_{H,0}^{-1}\vec{\Delta}_0}{\vec{p}^T M D_{H,0}^{-1}\vec{p}} $$

后续节点：

$$ \dot{\lambda}_n=\frac{\lambda_n-\lambda_{n-1}}{\Delta t_n} $$

以及：

$$ \lambda_n=\frac{\vec{g}_n^T M D_{H,n}^{-1}\left(\frac{\vec{p}}{\Delta t_n}\lambda_{n-1}-\vec{\Delta}_n\right)}{\vec{g}_n^T M D_{H,n}^{-1}\vec{g}_n} $$

其中：

$$ \vec{g}_n=\frac{\vec{p}}{\Delta t_n}-D_{H,n}\vec{s} $$

这样，current sequential Eq. (72) temporal treatment 已经从 $t_0$ 闭合到 $T$。

---

# 15. 回到 Eq. (70)–(71) 中的 $a$

Eq. (67) 中：

$$ W^{-1}=\langle H_\sigma\lambda^2\rangle+\langle\dot{\lambda}\lambda\rangle C^{-1} $$

定义：

$$ a=\langle\dot{\lambda}\lambda\rangle $$

Eq. (71) 中：

$$ \bar{\varepsilon}^p=\frac{1}{a}\bar{\tilde{\varepsilon}}-C^{-1}\bar{\sigma} $$

因此 $a$ 同时影响：

- spatial effective operator $W$；
- Eq. (71) recovery 中的 $1/a$。

如果 $a$ 接近零，Eq. (71) 会发生明显退化。

---

# 16. continuous $a$ 的端点恒等式

连续定义：

$$ a=\int_0^T\dot{\lambda}(t)\lambda(t)\,dt $$

利用：

$$ \frac{d}{dt}\lambda^2=2\lambda\dot{\lambda} $$

可得：

$$ a=\frac{1}{2}\left[\lambda^2(T)-\lambda^2(0)\right] $$

由于：

$$ \lambda(0)=0 $$

所以：

$$ a=\frac{1}{2}\lambda^2(T) $$

这说明 $\lambda_0=0$ 本身并不会导致 $a=0$。

---

# 17. continuous formulation 对完整循环的潜在问题

若 temporal mode 在完整历史末端又回到：

$$ \lambda(T)=0 $$

则 continuous identity 给出：

$$ a=0 $$

即使中间：

$$ \lambda(t)\ne0 $$

且经历了显著循环演化。

因此对循环疲劳问题，不能未经离散分析就直接把：

$$ a=\frac{1}{2}\left[\lambda^2(T)-\lambda^2(0)\right] $$

作为 tower code 中 Eq. (71) 的 numerical coefficient。

---

# 18. 为什么 original DG0 信息仍然重要

原论文 temporal minimisation 使用 DG0。

但当前掌握的 2018 原文没有给出完整的：

- temporal trial/test space；
- jump terms；
- left/right traces；
- first slab treatment；
- endpoint treatment；
- final DG0 algebra。

因此：

$$ a_{\rm continuous} $$

与 original DG0 中真正使用的 discrete coefficient 是否完全相同，目前无法证明。

所以 tower v1 仍不能宣称 current backward Euler 与 original DG0 exact equivalent。

---

# 19. tower v1 的 BE-consistent temporal contraction

为了与 current sequential backward-Euler temporal solver 保持一致，可定义：

$$ a_h=\sum_{n=1}^{N}\Delta t_n\lambda_n\dot{\lambda}_n $$

其中：

$$ N=N_t-1 $$

代入：

$$ \dot{\lambda}_n=\frac{\lambda_n-\lambda_{n-1}}{\Delta t_n} $$

得到：

$$ a_h=\sum_{n=1}^{N}\lambda_n(\lambda_n-\lambda_{n-1}) $$

该定义采用 right-endpoint time-slab interpretation。

需要强调：

> 这是 current tower v1 为保持 backward-Euler 离散一致性而采用的 project discrete coefficient，不是目前已经证明的 original DG0 coefficient。

---

# 20. $a_h$ 的离散恒等式

利用：

$$ \lambda_n(\lambda_n-\lambda_{n-1})=\frac{1}{2}\left[\lambda_n^2-\lambda_{n-1}^2+(\lambda_n-\lambda_{n-1})^2\right] $$

求和得到：

$$ a_h=\frac{1}{2}\left(\lambda_N^2-\lambda_0^2\right)+\frac{1}{2}\sum_{n=1}^{N}(\lambda_n-\lambda_{n-1})^2 $$

由于：

$$ \lambda_0=0 $$

所以：

$$ a_h=\frac{1}{2}\lambda_N^2+\frac{1}{2}\sum_{n=1}^{N}(\lambda_n-\lambda_{n-1})^2 $$

因此：

$$ a_h\ge0 $$

且对非平凡 temporal mode：

$$ a_h>0 $$

---

# 21. 为什么 BE 离散后完整循环不必使 $a_h=0$

若：

$$ \lambda_0=0 $$

且：

$$ \lambda_N=0 $$

则：

$$ a_h=\frac{1}{2}\sum_{n=1}^{N}(\lambda_n-\lambda_{n-1})^2 $$

只要 temporal history 不是恒等于零：

$$ \exists n:\lambda_n\ne\lambda_{n-1} $$

就有：

$$ a_h>0 $$

因此 BE-consistent discrete coefficient 不会因为 initial 和 final amplitude 相同就自动退化。

---

# 22. BE 离散中的额外非负项

continuous endpoint identity 中只有：

$$ \frac{1}{2}\left(\lambda_N^2-\lambda_0^2\right) $$

而 backward-Euler discrete contraction 多出：

$$ \frac{1}{2}\sum_{n=1}^{N}(\lambda_n-\lambda_{n-1})^2 $$

这一项是非负的。

可将其理解为 current right-endpoint backward-Euler temporal discretisation 自然产生的 discrete jump / numerical dissipation contribution。

但目前不能把它直接称为 original DG0 jump term。

---

# 23. Eq. (67) 其他 temporal contractions 的离散一致性

若采用同一 right-endpoint time-slab interpretation，则对 material point $q$：

$$ A_{q,h}=\sum_{n=1}^{N}\Delta t_n H_{\sigma,q,n}\lambda_n^2 $$

同时：

$$ \bar{\delta}_{q,h}=\sum_{n=1}^{N}\Delta t_n\bar{\Delta}_{q,n}\lambda_n $$

于是：

$$ W_{q,h}^{-1}=A_{q,h}+\frac{a_h}{E_0} $$

这样 Eq. (67)、Eq. (70)–(71) 与 Eq. (72) 使用同一套 backward-Euler time-slab interpretation。

---

# 24. Eq. (71) 的 BE-consistent form

tower v1 中可写为：

$$ \vec{\bar{\varepsilon}}^p=\frac{1}{a_h}H\vec{\bar{\tilde{U}}}-C_0^{-1}\vec{\bar{\sigma}} $$

其中：

$$ a_h=\sum_{n=1}^{N}\lambda_n(\lambda_n-\lambda_{n-1}) $$

该表达保持 Eq. (71) 的 operator structure 不变，只明确 continuous temporal contraction 如何在 current project backward-Euler discretisation 中实现。

---

# 25. $\dot{\lambda}_0$ 是否进入 $a_h$

如果采用：

$$ a_h=\sum_{n=1}^{N}\Delta t_n\lambda_n\dot{\lambda}_n $$

求和从：

$$ n=1 $$

开始。

所以：

$$ \dot{\lambda}_0 $$

不直接进入 $a_h$。

同时：

$$ \lambda_0=0 $$

因此 $t_0$ 也不直接贡献：

$$ A_{q,h} $$

或：

$$ \bar{\delta}_{q,h} $$

所以 first-point rate solve 的主要职责是正确重建 $t_0$ mechanical residual，而不是为 Eq. (71) 人工制造 nonzero $a_h$。

---

# 26. fixed-point 第一轮为什么不能从 zero temporal function 开始 spatial solve

如果直接采用：

$$ \lambda^{(0)}(t)\equiv0 $$

进入 Eq. (70)–(71)，则：

$$ a_h=0 $$

同时：

$$ A_{q,h}=0 $$

以及：

$$ \bar{\delta}_{q,h}=0 $$

于是：

$$ W_{q,h}^{-1}=0 $$

Eq. (70) 无法形成有效 spatial stiffness-like operator。

Eq. (71) 同时出现：

$$ \frac{1}{a_h} $$

因此：

$$ \lambda^{(0)}(t)\equiv0 $$

不能直接作为 Eq. (70)–(71) spatial half-step 的启动 temporal function。

---

# 27. current 1D enrichment 实际如何避免这一问题

current `latin/pgd_enrichment.py` 的实际流程不是：

```text
zero temporal function
→ spatial solve
```

而是：

1. 从 current residual 建立 deterministic spatial seed；
2. 根据 spatial seed 计算 associated equilibrated stress mode；
3. 构造 temporary one-mode basis；
4. 先调用 `update_pgd_time_functions()` 求 temporal function；
5. 再调用 fixed-temporal spatial solver 更新 spatial function；
6. temporal / spatial 两个 half-steps 继续交替。

所以 current 1D flow 实际是：

$$ {\rm spatial\ seed}\rightarrow{\rm temporal\ solve}\rightarrow{\rm spatial\ update} $$

这自动避免 zero temporal function 直接进入 spatial solve 的退化。

---

# 28. tower v1 fixed-point 的初始化策略

tower v1 可继承同一 control-flow principle：

$$ {\rm remaining\ defect}\rightarrow{\rm initial\ plastic\ spatial\ seed} $$

然后：

$$ {\rm spatial\ seed}\rightarrow{\rm associated\ equilibrated\ stress\ seed} $$

随后：

$$ {\rm Eq.\ (72)\ temporal\ solve} $$

再进入：

$$ {\rm Eq.\ (70)-(71)\ spatial\ solve} $$

因此第一轮 fixed point 应理解为：

$$ \bar{\varepsilon}^{p,(0)}\rightarrow\lambda^{(1)}\rightarrow\bar{\varepsilon}^{p,(1)}\rightarrow\lambda^{(2)}\rightarrow\cdots $$

而不是：

$$ \lambda^{(0)}=0\rightarrow\bar{\varepsilon}^{p,(1)} $$

---

# 29. 关于 $a_h$ 的当前状态必须保持谨慎

本阶段得到：

$$ a_h=\sum_{n=1}^{N}\lambda_n(\lambda_n-\lambda_{n-1}) $$

是 current backward-Euler tower v1 的自然离散选择。

但仍需明确：

- 它不是原论文正文直接给出的公式；
- 它不是目前已证明的 original DG0 algebra；
- 它属于 paper residual structure 在 current project BE discretisation 下的一致实现；
- 如果未来恢复 exact DG0 algebra，应重新比较 $a_h$、$A_{q,h}$ 和 $\bar{\delta}_{q,h}$。

---

# 30. 进入 Eq. (72) denominator conditioning 问题

对于 $n\ge1$：

$$ \vec{r}_n=\vec{g}_n\lambda_n-\vec{b}_n $$

其中：

$$ \vec{g}_n=\frac{\vec{p}}{\Delta t_n}-D_{H,n}\vec{s} $$

以及：

$$ \vec{b}_n=\frac{\vec{p}}{\Delta t_n}\lambda_{n-1}-\vec{\Delta}_n $$

定义：

$$ W_n=M D_{H,n}^{-1} $$

normal-equation denominator：

$$ A_n=\vec{g}_n^T W_n\vec{g}_n $$

normal-equation numerator：

$$ c_n=\vec{g}_n^T W_n\vec{b}_n $$

因此：

$$ \lambda_n=\frac{c_n}{A_n} $$

---

# 31. $A_n$ 的准确数学含义

由于：

$$ W_n>0 $$

所以：

$$ A_n=\|\vec{g}_n\|_{W_n}^2 $$

因此：

$$ A_n\ge0 $$

且：

$$ A_n=0 $$

等价于：

$$ \vec{g}_n=\vec{0} $$

所以 denominator 为零不是一个任意 numerical accident，而表示 current temporal residual-control direction 本身消失。

---

# 32. $\vec{g}_n$ 的 controllability 含义

由：

$$ \vec{r}_n=\vec{g}_n\lambda_n-\vec{b}_n $$

有：

$$ \frac{\partial\vec{r}_n}{\partial\lambda_n}=\vec{g}_n $$

因此 $\vec{g}_n$ 表示改变 $\lambda_n$ 时 residual 在 material-point space 中的变化方向。

若：

$$ \vec{g}_n=\vec{0} $$

则：

$$ \frac{\partial\vec{r}_n}{\partial\lambda_n}=0 $$

此时改变 $\lambda_n$ 无法改变 current residual。

因此 $A_n\approx0$ 的本质是：

> current spatial candidate 在 time step $n$ 对 LATIN residual 的 temporal control authority 发生退化。

---

# 33. $\vec{g}_n$ 为什么可能很小

定义：

$$ \vec{u}_n=\frac{\vec{p}}{\Delta t_n} $$

以及：

$$ \vec{v}_n=D_{H,n}\vec{s} $$

则：

$$ \vec{g}_n=\vec{u}_n-\vec{v}_n $$

所以 $\vec{g}_n\approx0$ 可能来自：

$$ \frac{\vec{p}}{\Delta t_n}\approx D_{H,n}\vec{s} $$

即：

> plastic-rate sensitivity 与 stress/search-direction sensitivity 在 weighted material-point space 中发生近似抵消。

这不是普通 stiffness singularity，而是 reduced temporal controllability cancellation。

---

# 34. 若 $A_n=0$，objective 会发生什么

定义：

$$ B_n=\vec{b}_n^T W_n\vec{b}_n $$

则：

$$ J_n(\lambda_n)=\frac{1}{2}\left(A_n\lambda_n^2-2c_n\lambda_n+B_n\right) $$

若：

$$ A_n=0 $$

则：

$$ \vec{g}_n=0 $$

因此：

$$ c_n=0 $$

于是：

$$ J_n=\frac{1}{2}B_n $$

完全不依赖：

$$ \lambda_n $$

所以 current one-step minimisation 对 $\lambda_n$ 失去唯一辨识能力。

---

# 35. 为什么不能简单 clip denominator

如果采用：

$$ A_n^{\rm used}=\max(A_n,A_{\min}) $$

那么实际求解变成：

$$ A_n^{\rm used}\lambda_n=c_n $$

而不再是原来的：

$$ A_n\lambda_n=c_n $$

因此 denominator clipping 直接改变了 Eq. (72) 的 discrete normal equation。

这不是单纯 diagnostic，而是 numerical model modification。

---

# 36. PGD scaling invariance 进一步说明绝对 $A_n$ 阈值不可靠

PGD pair 存在 scaling freedom。

若：

$$ \vec{p}^*=c\vec{p} $$

则 associated stress mode 同时变为：

$$ \vec{s}^*=c\vec{s} $$

因此：

$$ \vec{g}_n^*=c\vec{g}_n $$

所以：

$$ A_n^*=c^2A_n $$

同时：

$$ c_n^*=c\,c_n $$

于是 temporal amplitude 相应变为：

$$ \lambda_n^*=\frac{1}{c}\lambda_n $$

但 physical separated correction：

$$ \lambda_n^*\vec{p}^*=\lambda_n\vec{p} $$

不变。

因此：

> $A_n$ 的绝对数值依赖 spatial-mode scaling，本身不能作为尺度无关的 degeneration criterion。

---

# 37. scale-invariant controllability ratio

定义：

$$ \eta_n=\frac{\|\vec{g}_n\|_{W_n}}{\|\vec{u}_n\|_{W_n}+\|\vec{v}_n\|_{W_n}} $$

其中：

$$ \vec{u}_n=\frac{\vec{p}}{\Delta t_n} $$

$$ \vec{v}_n=D_{H,n}\vec{s} $$

由于 numerator 和 denominator 在 spatial scaling 下同时乘 $|c|$，因此：

$$ \eta_n^*=\eta_n $$

所以 $\eta_n$ 是 PGD scaling invariant 的。

---

# 38. $\eta_n$ 的取值意义

由 triangle inequality：

$$ 0\le\eta_n\le1 $$

若：

$$ \eta_n\approx1 $$

说明两部分没有显著 cancellation。

若：

$$ \eta_n\ll1 $$

说明：

$$ \vec{u}_n\approx\vec{v}_n $$

即 current temporal residual-control direction 由两个较大 contribution 的 cancellation 形成。

因此 $\eta_n$ 可作为：

> temporal residual controllability ratio

或者：

> reduced-column cancellation indicator。

---

# 39. local residual-reduction capability

仅有 controllability 还不足以判断 candidate mode 是否有效。

定义：

$$ B_n=\vec{b}_n^T W_n\vec{b}_n $$

并定义：

$$ \rho_n=\frac{|c_n|}{\sqrt{A_nB_n}} $$

在：

$$ A_n>0,\quad B_n>0 $$

时，根据 Cauchy-Schwarz：

$$ 0\le\rho_n\le1 $$

---

# 40. $\rho_n^2$ 的准确意义

无 correction 时：

$$ J_n(0)=\frac{1}{2}B_n $$

最优 scalar amplitude 下：

$$ J_n^*=\frac{1}{2}\left(B_n-\frac{c_n^2}{A_n}\right) $$

因此：

$$ J_n(0)-J_n^*=\frac{c_n^2}{2A_n} $$

relative local residual reduction：

$$ \frac{J_n(0)-J_n^*}{J_n(0)}=\frac{c_n^2}{A_nB_n} $$

所以：

$$ \rho_n^2=\frac{J_n(0)-J_n^*}{J_n(0)} $$

因此 $\rho_n^2$ 表示 current single temporal DOF 在该 time step 理论上能够消除的 weighted residual fraction。

---

# 41. $\eta_n$ 与 $\rho_n$ 的职责不同

$\eta_n$ 回答：

> current residual-control direction 是否由于 cancellation 而退化？

这是 conditioning / identifiability 问题。

$\rho_n$ 回答：

> current residual-control direction 与需要修正的 residual 是否对齐？

这是 effectiveness / projection quality 问题。

因此不能用其中一个完全替代另一个。

---

# 42. 一个典型情况：conditioning 好但 mode 无效

可能出现：

$$ \eta_n\approx1 $$

但：

$$ \rho_n\approx0 $$

此时 current temporal control direction 数值上很清晰。

但它与：

$$ \vec{b}_n $$

几乎正交。

所以该 spatial mode 在这个 time step 对 residual 几乎没有 reduction capability。

---

# 43. 另一个典型情况：direction 很弱但 alignment 很好

也可能出现：

$$ \eta_n\ll1 $$

但：

$$ \rho_n\approx1 $$

这说明 $\vec{g}_n$ 很弱，但方向与 current residual 很对齐。

此时 exact scalar minimiser 可能需要很大的：

$$ \lambda_n $$

来生成有限 correction。

因为 PGD 本身具有 spatial-temporal scaling freedom，所以 large temporal amplitude 本身不应自动视为错误。

真正需要检查的是 correction 和 residual 是否保持 finite 和稳定。

---

# 44. current 1D `condition_history` 的 single-mode blind spot

current `_weighted_least_squares()` 使用 SVD singular values 计算 reduced matrix condition number：

$$ \kappa=\frac{\sigma_{\max}}{\sigma_{\min}} $$

对于 single-new-mode，weighted reduced matrix 只有一列。

只要该列不严格为 zero，就只有一个 nonzero singular value：

$$ \sigma_1 $$

所以：

$$ \kappa=\frac{\sigma_1}{\sigma_1}=1 $$

因此即使：

$$ \|\vec{g}_n\|_{W_n} $$

非常小，single-column condition number 仍可能为 1。

所以 current `condition_history` 对 multi-mode reduced solve 有意义，但不能识别 single-mode residual-control column 本身接近 zero 的退化。

---

# 45. tower v1 为什么需要额外 controllability diagnostic

由于 single-mode SVD condition number 无法识别 column magnitude / cancellation，所以 tower enrichment 需要额外记录：

$$ \eta_n $$

用于检查：

$$ \frac{\vec{p}}{\Delta t_n} $$

与：

$$ D_{H,n}\vec{s} $$

之间是否发生近 cancellation。

这不是替代 current `condition_history`。

未来 multi-mode temporal update 中，两个 diagnostics 可以同时保留：

- reduced matrix SVD condition number；
- per-mode / per-step controllability ratio。

---

# 46. exact degeneration case 1

若：

$$ \vec{g}_n=0 $$

且：

$$ \vec{b}_n\ne0 $$

则：

$$ \vec{r}_n=-\vec{b}_n $$

无论 $\lambda_n$ 取什么值，current residual 都不改变。

因此 current spatial candidate 在该 time step 没有 temporal control authority。

这种情况下不应人为设置 denominator。

更合理的是：

> 将 current candidate pair 标记为 temporally degenerate，并触发 candidate rejection 或 spatial-seed reinitialisation。

---

# 47. exact degeneration case 2

若：

$$ \vec{g}_n=0 $$

同时：

$$ \vec{b}_n=0 $$

则：

$$ \vec{r}_n=0 $$

而且：

$$ J_n=0 $$

此时任意：

$$ \lambda_n $$

都是 one-step minimiser。

这是 exact non-identifiability，而不是 residual error。

为了避免产生 artificial temporal jump，可选择：

$$ \lambda_n=\lambda_{n-1} $$

于是：

$$ \dot{\lambda}_n=0 $$

这一选择没有改变 objective minimiser，因为所有 $\lambda_n$ 本来都是 minimisers。

---

# 48. near-degenerate case

若：

$$ \eta_n\ll1 $$

但：

$$ A_n>0 $$

数学上 exact scalar minimiser：

$$ \lambda_n=\frac{c_n}{A_n} $$

仍然存在。

因此不能仅因为 $\eta_n$ 小就直接修改 denominator。

第一步仍应计算 exact scalar solution，然后执行 a posteriori checks。

---

# 49. near-degenerate case 的 a posteriori checks

至少应检查：

1. $\lambda_n$ 是否 finite；
2. $\dot{\lambda}_n$ 是否 finite；
3. resulting separated correction 是否 finite；
4. local weighted residual 是否相对于 zero correction 降低；
5. first-order stationarity 是否满足；
6. subsequent spatial half-step 是否保持 stable；
7. fixed-point correction 是否发生异常爆增。

只有在这些 checks 失败时，才需要 reject 或 restart candidate pair。

---

# 50. 为什么 large $\lambda_n$ 本身不能作为 rejection criterion

PGD pair 有 scaling invariance：

$$ \lambda(t)\bar{\varepsilon}^p(x) $$

若 spatial mode 被 normalization：

$$ \bar{\varepsilon}^{p,*}=c\bar{\varepsilon}^p $$

则 temporal amplitude 可以自动变成：

$$ \lambda^*=\frac{\lambda}{c} $$

所以 temporal amplitude 的绝对大小依赖 spatial normalization。

因此：

> large $\lambda_n$ 只能作为 warning signal，不能独立决定 candidate rejection。

需要结合 correction magnitude、residual reduction 和 controllability 一起判断。

---

# 51. tower v1 safeguard 的核心原则

当前可以确定：

> detect → diagnose → verify → reject/restart if necessary

而不是：

> detect → regularise denominator

即：

1. 计算 $\vec{g}_n$；
2. 计算 $A_n$；
3. 计算 $\eta_n$；
4. 正常情况下使用原始 scalar equation；
5. 计算 local residual reduction；
6. 检查 finite / stationarity / fixed-point stability；
7. 真正失去 controllability 时 reject 或 restart current pair；
8. 不修改 Eq. (72) denominator。

---

# 52. 目前不固定具体 $\eta$ threshold

当前不应立即写死：

```text
eta_tol = 1e-8
```

或其他任意值。

threshold 需要在 implementation/testing 阶段根据：

- float64 precision；
- tower fiber spatial normalization；
- $H_\sigma$ 的实际数量级；
- $\Delta t$；
- current 1D benchmark 的 $\eta_n$ 分布；
- tower benchmark 的 material-point distribution；

进行标定。

因此本阶段只固定 diagnostic definition，不固定 tolerance。

---

# 53. denominator clipping 为什么不符合当前研究路线

当前研究策略是：

> 尽可能保留 Bhattacharyya et al. 原论文的 LATIN-PGD residual/minimisation structure，只对 tower spatial discretisation 和必要 temporal discretisation 做清晰、可追踪的迁移。

denominator clipping：

$$ A_n\leftarrow\max(A_n,A_{\min}) $$

直接修改 discrete normal equation。

candidate rejection / reinitialisation 则只是拒绝一个无法形成可靠 reduced correction 的 candidate。

因此后者更符合当前 paper-faithful migration strategy。

---

# 54. 当前可正式采用的 Eq. (72) diagnostics

normal-equation quantities：

$$ A_n=\vec{g}_n^T W_n\vec{g}_n $$

$$ c_n=\vec{g}_n^T W_n\vec{b}_n $$

$$ B_n=\vec{b}_n^T W_n\vec{b}_n $$

controllability ratio：

$$ \eta_n=\frac{\|\vec{g}_n\|_{W_n}}{\|\vec{p}/\Delta t_n\|_{W_n}+\|D_{H,n}\vec{s}\|_{W_n}} $$

local projection quality：

$$ \rho_n=\frac{|c_n|}{\sqrt{A_nB_n}} $$

local residual-reduction fraction：

$$ \rho_n^2=\frac{J_n(0)-J_n^*}{J_n(0)} $$

---

# 55. 本阶段对 current 1D implementation 的结论

current 1D code 已经具备：

- first-point amplitude zero；
- first-point rate weighted least-squares；
- later backward-Euler sequential time marching；
- weighted least-squares solve；
- SVD-based reduced matrix condition diagnostics；
- full residual reconstruction；
- weighted residual norm；
- enrichment spatial-seed-first control flow。

这些都可以作为 tower v1 的重要 reference implementation。

但 tower v1 还需要补充：

- Eq. (70)–(71) explicit spatial W-FE solve；
- BE-consistent temporal contraction $a_h$ 等；
- single-mode controllability ratio $\eta_n$；
- exact / near-degenerate candidate handling；
- tower material-point $(e,g,f)$ integration。

---

# 56. 本阶段尚未进入代码的原因

虽然 initial-time treatment、BE contraction 和 conditioning logic 已基本闭合，但 enrichment fixed-point 的完整 convergence logic 尚未重新审查。

当前还需要明确：

- temporal update 后何时返回 Eq. (70)–(71)；
- spatial normalization 应发生在 fixed point 内还是外；
- normalization 后 temporal function 如何 rescale；
- fixed-point convergence 应比较 spatial mode、temporal mode，还是 full separated correction；
- current 1D `fixed_point_change` 是否应直接迁移到 tower；
- accepted pair 前是否保留 scalar line search；
- all-mode temporal reoptimisation 应发生在什么位置。

所以当前仍不进入 tower LATIN-PGD code implementation。

---

# 57. 本阶段最终结论一：initial temporal treatment

tower v1 对 virgin initial state 可采用：

$$ \lambda_0=0 $$

但：

$$ \dot{\lambda}_0 $$

通过 initial weighted least-squares 求得：

$$ \dot{\lambda}_0=-\frac{\vec{p}^T M D_{H,0}^{-1}\vec{\Delta}_0}{\vec{p}^T M D_{H,0}^{-1}\vec{p}} $$

这与 current validated 1D `pgd_time_update.py` 一致。

---

# 58. 本阶段最终结论二：BE temporal contraction

tower v1 可暂时采用：

$$ a_h=\sum_{n=1}^{N}\Delta t_n\lambda_n\dot{\lambda}_n $$

即：

$$ a_h=\sum_{n=1}^{N}\lambda_n(\lambda_n-\lambda_{n-1}) $$

并有：

$$ a_h=\frac{1}{2}\lambda_N^2+\frac{1}{2}\sum_{n=1}^{N}(\lambda_n-\lambda_{n-1})^2 $$

所以对非平凡 temporal mode：

$$ a_h>0 $$

这一结论有利于避免 Eq. (71) 的 $1/a_h$ 在完整循环端点相同情况下自动退化。

---

# 59. 本阶段最终结论三：fixed-point initialization

第一轮 enrichment 不应采用：

$$ \lambda^{(0)}(t)\equiv0\rightarrow{\rm Eq.\ (70)} $$

而应采用：

$$ {\rm spatial\ seed}\rightarrow{\rm Eq.\ (72)\ temporal\ solve}\rightarrow{\rm Eq.\ (70)-(71)\ spatial\ solve} $$

这与 current 1D enrichment control flow 一致。

---

# 60. 本阶段最终结论四：denominator conditioning

Eq. (72) denominator：

$$ A_n=\vec{g}_n^T W_n\vec{g}_n $$

接近零表示 current temporal residual-control direction 退化。

不能简单采用：

$$ A_n\leftarrow\max(A_n,A_{\min}) $$

而应监测尺度无关：

$$ \eta_n=\frac{\|\vec{g}_n\|_{W_n}}{\|\vec{p}/\Delta t_n\|_{W_n}+\|D_{H,n}\vec{s}\|_{W_n}} $$

并结合：

$$ \rho_n^2=\frac{c_n^2}{A_nB_n} $$

判断 controllability 与 local residual-reduction capability。

---

# 61. 本阶段最终结论五：safeguard philosophy

tower v1 safeguard 应优先采用：

> candidate detection + diagnostic + residual verification + reject/restart

而不是：

> denominator clipping / artificial stiffness regularisation。

这样可以最大限度保留 original Eq. (72) residual minimisation structure。

---

# 62. 下一阶段

下一阶段只处理一个问题：

> Eq. (72) temporal update 完成后，如何严格返回 Eq. (70)–(71)，形成完整 enrichment fixed-point loop，并定义 new separated pair 的 fixed-point convergence criterion。

需要重点回答：

- fixed-point iteration 的完整顺序；
- spatial 与 temporal half-step 的数据依赖；
- normalization 与 rescaling 的正确时机；
- 应比较 $\bar{\varepsilon}^p$、$\lambda$，还是 full correction；
- current 1D `fixed_point_change` 是否适合 tower；
- convergence tolerance 的物理和数值意义。

在这一问题闭合之前，不进入 tower LATIN-PGD code implementation。
