# 海上风机塔筒 100 周期 FOM、$n\times\tau$ 张量化与 SVD 低秩诊断阶段总结

**项目：** Offshore Wind Turbine and LATIN-PGD
**阶段主题：** 从 100 周期全阶非线性参考解到 LATIN-PGD 可分离性诊断
**日期：** 2026-08-09
**建议仓库路径：** `docs/2026-08-09-tower-100cycle-fom-tensorization-low-rank-stage-summary.md`

---

## 1. 本阶段工作的核心目的

前一阶段已经完成海上风机塔筒非对称循环加载下的黏塑性损伤响应分析，并通过 20、50、100 周期计算逐步识别出棘轮漂移、循环塑性路径变化和损伤速率转折等现象。本阶段进一步向 LATIN-PGD 降阶方法推进，目标不再只是“继续增加循环数”，而是回答一个更基础的问题：

> **当前塔筒高周循环非线性响应，是否本身具有足够强的低秩和可分离结构，从而值得进一步构造 LATIN-PGD 降阶求解框架？**

本阶段建立的分析链条为：

$$\text{100-cycle FOM}\;\longrightarrow\;t\mapsto(n,\tau)\;\longrightarrow\;q(n,\tau,x)\;\longrightarrow\;\text{mode-wise SVD}\;\longrightarrow\;\text{empirical low-rank diagnosis}$$

其中：

1. **100 周期 FOM**：提供不经过 LATIN-PGD 降阶的完整参考解；
2. **$n\times\tau$ 张量化**：把传统单一时间轴重新组织为“慢循环数 × 快周期相位”；
3. **SVD 低秩诊断**：判断位移、应力、塑性应变和损伤场在循环、周期相位和空间三个方向上的可压缩程度；
4. **LATIN-PGD 可行性判断**：据此判断后续是否值得采用少量分离模态表示完整高周循环响应。

本阶段核心结论为：

> **当前 100 周期塔筒 FOM 响应表现出明显的多方向低秩结构。位移和应力高度低秩；塑性应变与损伤的累计场同样高度低秩，而它们的单周期增量场具有更高但仍然有限的 cycle-space 复杂度。结果总体支持继续推进 LATIN-PGD。**

---

## 2. 什么是 100 周期 FOM

### 2.1 FOM 的含义

FOM 是 **Full-Order Model**，即全阶模型。

在当前项目中，FOM 指已经建立并验证的塔筒非线性有限元模型。它不使用 LATIN-PGD、POD、SVD 截断或其他降阶近似，而是沿完整时间历史逐增量求解。

对于每一个加载增量，程序都需要进行完整的非线性平衡迭代，并在所有梁单元、Gauss 点和纤维位置更新材料状态。形式上可概括为：

$$\mathbf{K}_{\mathrm{t}}(\mathbf{u},\mathbf{z})\,\Delta\mathbf{u}=\mathbf{R}$$

其中：

- $\mathbf{u}$：结构自由度；
- $\mathbf{z}$：材料内部变量集合；
- $\mathbf{K}_{\mathrm{t}}$：当前切线刚度；
- $\mathbf{R}$：当前非平衡残差。

当前黏塑性损伤材料模型中保存的主要内部变量包括：

$$\varepsilon_p,\quad \alpha,\quad \bar r,\quad D$$

分别对应塑性应变、运动硬化相关变量、各向同性硬化/累积相关变量和损伤变量。

因此，所谓“100 周期 FOM”并不是一种新的数值算法，而是：

> **用现有完整非线性有限元求解器，把 100 个循环逐时间增量全部算完，并将其作为后续所有降阶结果的参考标准。**

### 2.2 本次 100 周期 FOM 的参数

| 项目 | 参数 |
|---|---:|
| 梁单元数 | 10 |
| Gauss 点数/单元 | 2 |
| 环向纤维数/截面 | 16 |
| 径向纤维层数 | 1 |
| 最大水平力 | $F_{\max}=1.0\ \mathrm{MN}$ |
| 力比 | $R_F=F_{\min}/F_{\max}=-0.5$ |
| 最小水平力 | $F_{\min}=-0.5\ \mathrm{MN}$ |
| 周期 | $T=10$ |
| 循环数 | 100 |
| 每周期增量数 | 40 |
| 周期阶段保存点数 | 4001 |

加载不是完全反向循环，而是具有非零平均值的非对称循环。

$$R_F=\frac{F_{\min}}{F_{\max}}=-0.5$$

因此：

$$F_{\min}=-0.5F_{\max}$$

$$F_{\mathrm{mean}}=\frac{F_{\max}+F_{\min}}{2}=0.25F_{\max}$$

$$F_a=\frac{F_{\max}-F_{\min}}{2}=0.75F_{\max}$$

周期荷载为：

$$F(t)=F_{\mathrm{mean}}+F_a\sin\left(\frac{2\pi t}{T}\right)$$

正式周期加载前首先进行：

$$0\rightarrow F_{\mathrm{mean}}$$

的显式预加载。预加载不计入疲劳循环数，预加载终点作为周期历史的第 0 个保存状态。

### 2.3 为什么必须先有 FOM

LATIN-PGD 的目标是减少大量逐时间步非线性求解，但任何降阶方法都必须与可信的参考解进行比较。

以后如果要评价 LATIN-PGD 的位移、应力、塑性应变、损伤误差以及计算时间减少比例，就必须先有一个未经降阶的标准答案。

因此：

$$\boxed{\text{FOM}=\text{reference solution / benchmark solution}}$$

只有在 FOM 基础上，后续“加速多少”和“误差多少”才具有明确含义。

---

## 3. 本次 100 周期 FOM 的基本结果

程序输出为：

```text
Discretisation: 10 elements, 2 Gauss/element, 16 x 1 fibers/section
Loading: Fmax=1 MN, Fmin=-0.5 MN, R_F=-0.5, cycles=100, increments/cycle=40
Periodic stored points: 4001
Cycle-phase grid: 100 cycles x 41 phase points
Full nodal displacement tensor: (100, 41, 33)
Fiber stress tensor: (100, 41, 10, 2, 16)
Fiber state tensor: (100, 41, 10, 2, 16, 4)
Final max |eps_p|=1.363490242e-03
Final max D=4.409564790e-02
Critical fiber: element=1, Gauss=1, fiber=5
Critical height=1.851205821e+00 m
Critical y=2.945089465e+00 m
Max Newton iterations=4
Max free-DOF residual=8.761245319e-03 N
Elapsed FOM=3057.346 s
Elapsed tensorization=0.035 s
```

因此：

$$\max|\varepsilon_p|=1.363490242\times10^{-3}$$

$$D_{\max}=4.409564790\times10^{-2}$$

临界位置仍然位于塔底第一单元危险纤维区域：

```text
element = 1
Gauss   = 1
fiber   = 5
```

最大 Newton 迭代数为 4，最大自由 DOF 残差约为：

$$8.76\times10^{-3}\ \mathrm{N}$$

这说明保存完整场历史和张量化并没有破坏原有非线性求解的稳定性。

### 3.1 计算成本

本次 100 周期 FOM 耗时：

$$3057.346\ \mathrm{s}\approx50.96\ \mathrm{min}$$

而张量化仅耗时：

$$0.035\ \mathrm{s}$$

说明当前真正昂贵的是 4000 个周期增量上的完整非线性求解，而不是后续数据重排和低秩分析。这进一步说明了开展 LATIN-PGD 的实际必要性。

---

## 4. 什么是 $n\times\tau$ 张量化

### 4.1 传统时间表示

传统逐时间步有限元把整个循环过程看成一根长时间轴：

$$t_0,t_1,t_2,\ldots,t_{4000}$$

这种表示并不显式区分当前属于第几个循环、也不显式区分一个循环内部处于哪个相位，因此没有直接利用高周循环中“循环重复、状态缓慢演化”的结构特征。

### 4.2 慢时间与快时间

因此把单一时间坐标重新表示为：

$$t\mapsto(n,\tau)$$

其中：

- $n$：循环编号，表示**慢时间尺度**；
- $\tau$：单个循环内部的相位，表示**快时间尺度**。

当前：

$$n=1,2,\ldots,100$$

每个循环有 40 个增量并保留两个端点，因此：

$$N_\tau=41$$

一个完整循环可理解为：

```text
tau = 0      : Fmean
tau = T/4    : +Fmax
tau = T/2    : Fmean
tau = 3T/4   : Fmin = -0.5 Fmax
tau = T      : Fmean
```

张量化时有意为每个循环保留闭合区间 $0\leq\tau\leq T$，因此相邻循环的公共端点会分别作为前一循环的终点和后一循环的起点出现。

### 4.3 从 $q(t,x)$ 到 $q(n,\tau,x)$

原始历史为：

$$q(t,x)$$

张量化后成为：

$$q(n,\tau,x)$$

节点位移本次维度为：

$$100\times41\times33$$

纤维标量场的空间离散数为：

$$10\times2\times16=320$$

因此：

$$\sigma,\varepsilon_p,D\in\mathbb{R}^{100\times41\times320}$$

需要强调：

> **张量化没有改变任何 FOM 数值结果，也没有进行降阶。它只是重新组织完全相同的数据。**

---

## 5. 为什么 $n\times\tau$ 表示与 LATIN-PGD 有直接关系

LATIN-PGD 希望利用不同坐标方向之间的可分离性。

理想情况下，希望将某一场变量表示为：

$$q(n,\tau,x)\approx\sum_{k=1}^{r}N_k(n)T_k(\tau)X_k(x)$$

其中：

- $N_k(n)$：慢循环演化函数；
- $T_k(\tau)$：单周期快时间函数；
- $X_k(x)$：空间模式；
- $r$：所需分离模态数量。

因此，张量化的真正意义是把高周疲劳中原本隐藏在长时间轴中的两类时间规律显式暴露出来：

$$\boxed{\text{slow evolution in }n}$$

$$\boxed{\text{fast cyclic waveform in }\tau}$$

这是从传统逐时间步有限元向多时间尺度 LATIN-PGD 表示过渡的关键桥梁。

---

## 6. 什么是 SVD 低秩诊断

设矩阵：

$$\mathbf{A}$$

其奇异值分解为：

$$\mathbf{A}=\mathbf{U}\mathbf{\Sigma}\mathbf{V}^{T}$$

其中：

$$\mathbf{\Sigma}=\operatorname{diag}(\sigma_1,\sigma_2,\ldots)$$

且：

$$\sigma_1\geq\sigma_2\geq\cdots\geq0$$

如果奇异值快速衰减，说明大部分数据集中在前几个主方向上，即使原矩阵很大，也可以用少数模态近似，这就是**低秩结构**。

本阶段采用平方 Frobenius 能量：

$$E_r=\frac{\sum_{i=1}^{r}\sigma_i^2}{\sum_i\sigma_i^2}$$

并统计达到 99%、99.9%、99.99% 能量所需要的最小 rank。

相对 Frobenius 截断误差为：

$$e_r=\sqrt{1-E_r}$$

因此：

| 保留能量 | 对应相对 Frobenius 误差 |
|---|---:|
| 99% | 10% |
| 99.9% | 3.162% |
| 99.99% | 1% |

这里是离散快照矩阵的 Frobenius 范数误差，并不等同于结构力学意义上的能量范数误差。

---

## 7. 为什么要同时分析 cycle / phase / space 三个方向

张量化后：

$$q=q(n,\tau,x)$$

已经是三维数据对象。

三个方向代表三个不同的物理问题：

| 方向 | 坐标 | 要回答的问题 |
|---|---|---|
| cycle | $n$ | 不同循环之间需要多少种慢演化模式？ |
| phase | $\tau$ | 一个循环内部需要多少种快时间波形？ |
| space | $x$ | 整个塔筒需要多少种空间响应模式？ |

LATIN-PGD 最终希望写成：

$$q(n,\tau,x)\approx\sum_{k=1}^{r}N_k(n)T_k(\tau)X_k(x)$$

所以必须同时判断三个坐标方向是否都存在低维结构。

### 7.1 cycle mode

固定一个循环 $n$，$q(n,\tau,x)$ 表示该循环内所有相位、所有空间位置的完整状态。

把每个循环的 $(\tau,x)$ 信息摊平成一行：

$$\mathbf{Q}_n\in\mathbb{R}^{N_n\times(N_\tau N_x)}$$

对于纤维标量场：

$$\mathbf{Q}_n\in\mathbb{R}^{100\times13120}$$

cycle-mode SVD 问的是：

> **100 个完整循环是否真的需要 100 个独立状态，还是只需要少数慢循环演化模式？**

### 7.2 phase mode

固定一个 $\tau$，$q(n,\tau,x)$ 表示所有循环在同一周期相位上的状态。

把 $(n,x)$ 展开：

$$\mathbf{Q}_\tau\in\mathbb{R}^{N_\tau\times(N_nN_x)}$$

对于纤维标量场：

$$\mathbf{Q}_\tau\in\mathbb{R}^{41\times32000}$$

phase-mode SVD 问的是：

> **一个完整循环内部的响应是否只需要少数几个基本快时间波形？**

### 7.3 space mode

固定空间位置 $x$，$q(n,\tau,x)$ 表示该位置跨越全部循环和全部相位的完整历史。

把 $(n,\tau)$ 展开：

$$\mathbf{Q}_x\in\mathbb{R}^{N_x\times(N_nN_\tau)}$$

对于纤维变量：

$$\mathbf{Q}_x\in\mathbb{R}^{320\times4100}$$

space-mode SVD 问的是：

> **塔筒数百个纤维积分位置是否真的需要数百个独立空间模式？**

### 7.4 为什么不能只看 cycle rank

假设：

$$q(n,\tau,x)=N(n)G(\tau,x)$$

则 cycle rank 可以等于 1。

但是如果：

$$G(\tau,x)=\sum_{k=1}^{100}T_k(\tau)X_k(x)$$

那么完整问题仍然可能需要大量分离项。

因此：

$$\boxed{\text{low cycle rank}\not\Rightarrow\text{low PGD rank}}$$

phase 和 space 方向同理，所以必须同时分析三个 mode。

### 7.5 三个方向与 PGD 因子的对应关系

| SVD 方向 | 物理意义 | 对应 PGD 因子 |
|---|---|---|
| cycle mode | 慢循环演化 | $N_k(n)$ |
| phase mode | 周期内快速变化 | $T_k(\tau)$ |
| space mode | 塔筒空间模式 | $X_k(x)$ |

LATIN-PGD 真正希望看到的是：

$$\boxed{\text{cycle 低秩}+\text{phase 低秩}+\text{space 低秩}}$$


---

## 8. HOSVD-style mode rank 与 PGD rank 的区别

本阶段做的并不是直接求 PGD rank。

对于三维张量 $q(n,\tau,x)$，分别沿 cycle、phase、space 展开并进行矩阵 SVD，得到的是类似 Tucker/HOSVD 意义下的多线性秩：

$$\left(r_n,r_\tau,r_x\right)$$

它表示：

- cycle 子空间需要多少维；
- phase 子空间需要多少维；
- space 子空间需要多少维。

而真正的 PGD/CP 形式是：

$$q(n,\tau,x)\approx\sum_{k=1}^{r_{\mathrm{PGD}}}N_k(n)T_k(\tau)X_k(x)$$

因此 $r_{\mathrm{PGD}}$ 与 $r_n,r_\tau,r_x$ 不是同一个概念。

例如：

$$\left(r_n,r_\tau,r_x\right)=(6,3,7)$$

不能直接写成：

$$r_{\mathrm{PGD}}=7$$

当前 SVD 的正确作用是：

> **诊断可分离性，而不是直接给出未来 PGD 所需模态数。**

---

## 9. 为什么同时分析 raw 与 cycle_increment

### 9.1 raw 场

原始场为：

$$q(n,\tau,x)$$

它同时包含：

1. 第 $n$ 个循环开始前已经积累的历史状态；
2. 当前循环内部新增的变化。

对于塑性应变和损伤，raw 场可能主要受到累计基线控制。

### 9.2 cycle_increment 场

定义：

$$\Delta q(n,\tau,x)=q(n,\tau,x)-q(n,0,x)$$

其中 $q(n,0,x)$ 是当前循环开始时的状态。

这一处理去除了每个循环的起始基线，但**没有**强迫周期末增量为零。

因此：

$$\Delta q(n,T,x)$$

仍然保留真实的：

- 棘轮漂移；
- 净塑性累积；
- 损伤累积。

对于塑性应变：

$$\Delta\varepsilon_p(n,\tau,x)=\varepsilon_p(n,\tau,x)-\varepsilon_p(n,0,x)$$

对于损伤：

$$\Delta D(n,\tau,x)=D(n,\tau,x)-D(n,0,x)$$

这两个增量场比 raw 场更适合判断：

> **一个循环内部的非线性演化机制是否随着循环数发生变化。**

---

## 10. 100 周期低秩诊断总体结果

以下采用最严格的 99.99% 能量标准，对应约 1% 的最优矩阵相对 Frobenius 截断误差。

| 场变量 | 表示 | cycle rank | phase rank | space rank |
|---|---|---:|---:|---:|
| $u$ | raw | 2 | 2 | 2 |
| $u$ | cycle increment | 2 | 2 | 1 |
| $\sigma$ | raw | 3 | 2 | 4 |
| $\sigma$ | cycle increment | 1 | 1 | 2 |
| $\varepsilon_p$ | raw | 3 | 1 | 3 |
| $\varepsilon_p$ | cycle increment | 6 | 3 | 7 |
| $D$ | raw | 3 | 1 | 3 |
| $D$ | cycle increment | 6 | 2 | 6 |

原始主要尺寸为：

$$N_n=100,\quad N_\tau=41,\quad N_x=320$$

而最复杂的 $\Delta\varepsilon_p$ 在 99.99% 标准下仅得到：

$$\left(r_n,r_\tau,r_x\right)=(6,3,7)$$

因此从经验低秩角度看，压缩潜力非常明显。

---

## 11. 位移场 $u$ 的详细解读

### 11.1 raw displacement

结果：

```text
u | raw | tensor shape = (100, 41, 33)

cycle | 99.00%:2  99.90%:2  99.99%:2
phase | 99.00%:2  99.90%:2  99.99%:2
space | 99.00%:1  99.90%:2  99.99%:2
```

奇异值比：

$$\sigma_2^{(cycle)}/\sigma_1^{(cycle)}=1.036\times10^{-1}$$

$$\sigma_3^{(cycle)}/\sigma_1^{(cycle)}=1.458\times10^{-3}$$

$$\sigma_2^{(phase)}/\sigma_1^{(phase)}=1.082\times10^{-1}$$

$$\sigma_3^{(phase)}/\sigma_1^{(phase)}=1.746\times10^{-3}$$

$$\sigma_2^{(space)}/\sigma_1^{(space)}=3.326\times10^{-2}$$

$$\sigma_3^{(space)}/\sigma_1^{(space)}=2.701\times10^{-4}$$

前两个模态之后奇异值迅速衰减，因此：

$$\boxed{(r_n,r_\tau,r_x)=(2,2,2)}$$

说明全局位移场具有非常明确的低秩性。

### 11.2 displacement cycle increment

结果：

```text
u | cycle_increment | tensor shape = (100, 41, 33)

cycle | 99.00%:1  99.90%:1  99.99%:2
phase | 99.00%:1  99.90%:1  99.99%:2
space | 99.00%:1  99.90%:1  99.99%:1
```

其中：

$$\sigma_2^{(cycle)}/\sigma_1^{(cycle)}=1.237\times10^{-2}$$

$$\sigma_2^{(phase)}/\sigma_1^{(phase)}=1.247\times10^{-2}$$

$$\sigma_2^{(space)}/\sigma_1^{(space)}=2.116\times10^{-3}$$

尤其在 99.9% 能量下：

$$\boxed{(r_n,r_\tau,r_x)=(1,1,1)}$$

这表明去除每个循环起始漂移以后，单周期位移增量的形状在整个 100 周期内极其稳定。

从物理上可以近似理解为：

$$\Delta u(n,\tau,x)\approx N_1(n)T_1(\tau)X_1(x)$$

即：

- 空间形状基本一致；
- 周期内波形基本一致；
- 主要变化表现为随循环数缓慢变化的幅值。

这对 LATIN-PGD 极为有利。

---

## 12. 应力场 $\sigma$ 的详细解读

### 12.1 raw stress

结果：

```text
sigma | raw | tensor shape = (100, 41, 320)

cycle | 99.00%:1  99.90%:2  99.99%:3
phase | 99.00%:2  99.90%:2  99.99%:2
space | 99.00%:2  99.90%:3  99.99%:4
```

99.99% 能量下：

$$\boxed{(r_n,r_\tau,r_x)=(3,2,4)}$$

虽然比位移略复杂，但仍然远低于原始规模。特别是 phase rank 只有 2，说明一个周期内部的应力路径仍然高度规则。

### 12.2 stress cycle increment

结果：

```text
sigma | cycle_increment | tensor shape = (100, 41, 320)

cycle | 99.00%:1  99.90%:1  99.99%:1
phase | 99.00%:1  99.90%:1  99.99%:1
space | 99.00%:1  99.90%:1  99.99%:2
```

因此：

$$\boxed{(r_n,r_\tau,r_x)=(1,1,2)}$$

奇异值比：

$$\sigma_2^{(cycle)}/\sigma_1^{(cycle)}=8.327\times10^{-3}$$

$$\sigma_2^{(phase)}/\sigma_1^{(phase)}=9.575\times10^{-3}$$

$$\sigma_2^{(space)}/\sigma_1^{(space)}=9.147\times10^{-3}$$

这意味着：

> **去除循环起始基线以后，100 个循环中的应力增量场几乎都属于同一个主 cycle-phase 模式，只需要少量空间修正。**

因此，应力响应同样提供了非常强的低秩证据。

---

## 13. 塑性应变 $\varepsilon_p$ 的详细解读

### 13.1 raw plastic strain

结果：

```text
eps_p | raw | tensor shape = (100, 41, 320)

cycle | 99.00%:1  99.90%:2  99.99%:3
phase | 99.00%:1  99.90%:1  99.99%:1
space | 99.00%:1  99.90%:2  99.99%:3
```

因此：

$$\boxed{(r_n,r_\tau,r_x)=(3,1,3)}$$

表面上看，塑性应变甚至比应力更低秩。但是不能直接据此认为塑性演化非常简单。

原因在于 raw 塑性应变包含长期累计的塑性基线：

$$\varepsilon_p(n,0,x)$$

随着循环推进，这一累计分量会逐渐占据 raw 场的大部分平方范数能量。因此 raw SVD 很容易首先识别出一个随循环数缓慢增长的累计塑性场，从而掩盖真正的周期内塑性路径变化。

### 13.2 plastic strain cycle increment

去除循环起始基线后：

$$\Delta\varepsilon_p(n,\tau,x)=\varepsilon_p(n,\tau,x)-\varepsilon_p(n,0,x)$$

结果变为：

```text
eps_p | cycle_increment | tensor shape = (100, 41, 320)

cycle | 99.00%:2  99.90%:5  99.99%:6
phase | 99.00%:1  99.90%:1  99.99%:3
space | 99.00%:2  99.90%:5  99.99%:7
```

99.99% 能量下：

$$\boxed{(r_n,r_\tau,r_x)=(6,3,7)}$$

cycle 奇异值比：

$$\sigma_2/\sigma_1=2.514\times10^{-1}$$

$$\sigma_3/\sigma_1=7.063\times10^{-2}$$

$$\sigma_5/\sigma_1=3.259\times10^{-2}$$

space 方向几乎同样：

$$\sigma_2/\sigma_1=2.508\times10^{-1}$$

$$\sigma_3/\sigma_1=7.038\times10^{-2}$$

$$\sigma_5/\sigma_1=3.283\times10^{-2}$$

第二、第三甚至第五个 mode 已经不能简单视为数值噪声。

### 13.3 塑性应变结果的物理意义

这表明：

> **不同循环中的塑性增量并不是完全相同的单周期波形乘以一个缓慢变化的幅值。**

随着循环数增加，发生塑性的空间区域、正反向塑性分配以及局部塑性增量比例都在调整。

但是值得特别注意：

$$r_\tau=3$$

明显低于：

$$r_n=6,\qquad r_x=7$$

说明主要复杂性并不是来自 $\tau$ 方向。

更合理的理解是：

$$\Delta\varepsilon_p(n,\tau,x)\approx\sum_k A_k(n,x)T_k(\tau)$$

其中所需的 $T_k(\tau)$ 很少，但 $A_k(n,x)$ 随循环演化的空间结构需要更多模态。

因此真正主要的复杂性来自：

$$\boxed{n\leftrightarrow x}$$

即慢循环演化与空间塑性分布之间的耦合。

---

## 14. 损伤场 $D$ 的详细解读

### 14.1 raw damage

结果：

```text
D | raw | tensor shape = (100, 41, 320)

cycle | 99.00%:1  99.90%:2  99.99%:3
phase | 99.00%:1  99.90%:1  99.99%:1
space | 99.00%:1  99.90%:2  99.99%:3
```

因此：

$$\boxed{(r_n,r_\tau,r_x)=(3,1,3)}$$

与 raw 塑性应变几乎相同。这主要反映损伤作为不可逆累计变量所具有的强慢变化背景。

### 14.2 damage cycle increment

定义：

$$\Delta D(n,\tau,x)=D(n,\tau,x)-D(n,0,x)$$

结果：

```text
D | cycle_increment | tensor shape = (100, 41, 320)

cycle | 99.00%:2  99.90%:4  99.99%:6
phase | 99.00%:1  99.90%:1  99.99%:2
space | 99.00%:2  99.90%:4  99.99%:6
```

99.99% 能量下：

$$\boxed{(r_n,r_\tau,r_x)=(6,2,6)}$$

cycle 奇异值比：

$$\sigma_2/\sigma_1=1.844\times10^{-1}$$

$$\sigma_3/\sigma_1=7.963\times10^{-2}$$

$$\sigma_5/\sigma_1=2.423\times10^{-2}$$

space 奇异值比：

$$\sigma_2/\sigma_1=1.840\times10^{-1}$$

$$\sigma_3/\sigma_1=7.965\times10^{-2}$$

$$\sigma_5/\sigma_1=2.426\times10^{-2}$$

两组数值极为接近。

而 phase 方向：

$$\sigma_2/\sigma_1=1.373\times10^{-2}$$

且 99.99% 只需要 rank 2。

### 14.3 损伤结果的物理意义

这表明：

> **损伤单周期增量的快时间形状非常低维，而主要复杂性存在于“损伤如何随循环数在塔筒空间中重新分布和演化”。**

因此可近似理解为：

$$\Delta D(n,\tau,x)\approx T_1(\tau)A_1(n,x)+T_2(\tau)A_2(n,x)+\cdots$$

其中 $\tau$ 方向只需要少量模式。

未来 LATIN-PGD 的关键并不是解决一个极其复杂的周期内快时间问题，而更可能是：

> **如何处理损伤和塑性在慢循环方向与空间方向之间逐渐增强的耦合。**

---

## 15. 与前一阶段 100 周期疲劳机理结果的联系

此前 100 周期分析已经得到以下关键现象：

1. cycle 46 附近，单周期损伤增量 $dD$ 达到局部最小值；
2. cycle 46 以后，$dD$ 重新持续增加；
3. 但循环末位移漂移 $\Delta u_n$ 和净塑性漂移 $\Delta\varepsilon_{p,n}$ 仍总体继续下降；
4. 因此不能称为整体棘轮加速或损伤失稳；
5. 更准确的状态是 **damage-rate acceleration with continuing net-drift stabilization**。

此前关键数值还包括：

$$dD_{\min}\approx dD_{46}$$

以及 cycle 100 时临界纤维的 net/path ratio：

$$\frac{\Delta\varepsilon_p^{\mathrm{net}}}{L_p}\approx0.492$$

该比值从早期接近 1 逐渐下降到约 0.49，说明循环内部反向塑性路径已经显著增强。

同时，cycle 100 相较无损伤对照的瞬时放大已经较明显：

- 单周期位移漂移放大约 32.6%；
- 单周期净塑性增量放大约 48.2%；
- 塑性路径长度约为无损伤模型的 2.30 倍；
- 单周期外功路径量约放大 39.1%。

然而，cycle 100 时 $\Delta u_n$ 和 $\Delta\varepsilon_{p,n}$ 仍在下降，因此这仍然不是整体失稳。

### 15.1 SVD 结果如何补充这一物理发现

raw 内部变量的 rank 很低：

$$\varepsilon_p:\ (3,1,3)$$

$$D:\ (3,1,3)$$

说明从累计场整体看，塑性和损伤仍然表现为平滑、规则的慢演化。

但是 cycle increment 明显更复杂：

$$\Delta\varepsilon_p:\ (6,3,7)$$

$$\Delta D:\ (6,2,6)$$

这说明：

> **虽然全局累计场仍然很规则，但每一个循环内部新增的塑性和损伤机制已经发生了更明显的 cycle-space 重组。**

这与此前观察到的现象一致：

- 全局位移范围变化仍较小；
- 应力范围变化仍有限；
- 但内部塑性路径、反向塑性比例和损伤速率已经明显演化。

因此，SVD 不只是数据压缩工具，它还从“解空间复杂度”的角度提供了新的物理解释：

$$\boxed{\text{internal mechanism evolution precedes strong global deterioration}}$$

即：

> **内部演化复杂度的增加可能早于明显的全局结构响应退化。**

---

## 16. 三类变量呈现出的低秩层级

### 16.1 第一层：全局结构响应高度低秩

位移：

$$u_{\mathrm{raw}}:(2,2,2)$$

$$\Delta u:(2,2,1)$$

应力：

$$\sigma_{\mathrm{raw}}:(3,2,4)$$

$$\Delta\sigma:(1,1,2)$$

这些结果非常有利于降阶。

### 16.2 第二层：累计内部变量同样高度低秩

塑性应变：

$$\varepsilon_{p,\mathrm{raw}}:(3,1,3)$$

损伤：

$$D_{\mathrm{raw}}:(3,1,3)$$

说明不可逆累计状态本身具有很强的规律性。

### 16.3 第三层：单周期不可逆增量是主要复杂来源

塑性应变增量：

$$\Delta\varepsilon_p:(6,3,7)$$

损伤增量：

$$\Delta D:(6,2,6)$$

它们比 raw 场明显复杂，但与原始尺寸相比仍然非常低。

因此未来 LATIN-PGD 真正需要重点处理的是：

$$\boxed{\text{slow-cycle evolution}\times\text{spatial redistribution}}$$

而不是完全失控的快时间复杂性。

---

## 17. 为什么当前结果支持继续推进 LATIN-PGD

当前主要原始维数为：

$$100\times41\times320$$

对于一个纤维标量场，共有：

$$100\times41\times320=1\,312\,000$$

个离散快照值。

而最复杂的 $\Delta\varepsilon_p$ 在 99.99% 能量下只有：

$$\left(r_n,r_\tau,r_x\right)=(6,3,7)$$

虽然这不能直接转换为 PGD rank，但它说明这个超过一百万离散值的三维响应，在每个坐标方向上实际上都只占据远低于原始维数的子空间：

$$6\ll100$$

$$3\ll41$$

$$7\ll320$$

这是非常明确的 empirical separability evidence。

因此当前严谨的结论不是“已经证明 PGD 只需要几个 mode”，而是：

> **已经有充分数值证据表明当前塔筒循环响应具有强低维结构，因此继续构造 LATIN-PGD 分离表示具有明确依据。**

---

## 18. 当前结果对 LATIN-PGD 算法设计的直接启示

### 18.1 不宜假设固定单一快时间波形

虽然 phase rank 总体很低，但 $\Delta\varepsilon_p$ 在 99.99% 能量下仍需要 3 个 phase mode。

因此过于简单的：

$$q(n,\tau,x)\approx A(n,x)T_1(\tau)$$

可能不足以高精度描述完整塑性路径。

未来应允许少量快时间模态富集。

### 18.2 更需要关注 cycle-space 耦合

对于 $\Delta\varepsilon_p$：

$$r_n=6,\qquad r_x=7$$

对于 $\Delta D$：

$$r_n=6,\qquad r_x=6$$

说明慢循环方向和空间方向是主要复杂来源。

因此后续算法应重点研究：

- 多个 $N_k(n)$；
- 多个 $X_k(x)$；
- 循环演化过程中空间模式是否需要更新；
- 是否需要 adaptive enrichment；
- 是否需要 stage-wise basis。

### 18.3 固定全寿命 basis 未必最优

前期机理分析已经识别出近似三个阶段：

**Stage I：cycles 1–20**
早期快速调整阶段。

**Stage II：cycles 21–46**
净漂移继续稳定，损伤增量下降并接近局部最低点。

**Stage III：cycles 47–100**
损伤速率重新增加，但净漂移仍继续下降。

如果不同阶段对应不同有效 rank，那么用同一套固定 basis 覆盖全部疲劳历史可能不是最经济的做法。

这自然引出：

$$\boxed{\text{stage-wise PGD basis}}$$

或：

$$\boxed{\text{adaptive modal enrichment}}$$

的可能性。

---

## 19. 当前 SVD 诊断的边界与需要避免的过度解释

本阶段结果非常积极，但必须保持方法论上的严格性。

### 19.1 当前是离散快照 SVD，不是 LATIN-PGD 求解

当前只是对已经得到的 FOM 数据做后处理。

它回答的是：

> “FOM 解是否具有低秩结构？”

而不是：

> “LATIN-PGD 已经成功求出了这些解。”

真正的 LATIN-PGD 求解器仍需要后续构造。

### 19.2 mode rank 不是 CP/PGD rank

不能把 $(r_n,r_\tau,r_x)$ 中的最大值直接当成 $r_{\mathrm{PGD}}$。

后续仍需真正检验：

$$q(n,\tau,x)\approx\sum_k N_k(n)T_k(\tau)X_k(x)$$

需要多少个乘积分离项才能达到目标误差。

### 19.3 当前 SVD 是未加权离散 Euclidean/Frobenius 诊断

目前每个快照、Gauss 点和纤维点在 SVD 中按离散样本等权处理。

因此当前 rank 是：

> **离散快照空间中的经验 rank**

而不是严格的连续场力学能量范数 rank。

以后若用于更严格的论文结论，可进一步考虑：

- Gauss 权重；
- 纤维面积权重；
- 单元长度权重；
- $\tau$ 方向积分权重；
- 物理场量纲和尺度归一化。

### 19.4 节点位移场包含不同类型自由度

当前 `nodal_displacements` 包含完整节点自由度，其中平移和转角量纲不同。

因此当前 $u$ 的 SVD 主要用于判断整体数值可压缩性。

以后若需要严格物理空间范数，可：

- 单独分析水平位移；
- 单独分析转角；
- 或引入力学意义明确的加权内积。

### 19.5 phase 端点的离散权重问题

当前每个循环使用 41 个闭合相位点，其中 $\tau=0$ 和 $\tau=T$ 分别作为一个循环的起点和终点保留。

现阶段 SVD 使用未加权离散 Euclidean 范数，因此两个端点都按完整样本权重参与统计。

这对于初步经验低秩诊断是可以接受的，但如果后续希望构造更接近连续时间 $L^2$ 内积的诊断，可考虑 $\tau$ 方向的梯形积分权重，使周期端点获得半权重。

---

## 20. 下一阶段最合理的工作路线

根据当前结果，不建议立即盲目增加至 200 周期。

更有价值的是利用现有 100 周期数据研究：

> **低秩结构如何随疲劳阶段发生变化。**

建议阶段划分为：

$$\text{Stage I}:1\sim20$$

$$\text{Stage II}:21\sim46$$

$$\text{Stage III}:47\sim100$$

重点比较：

$$\Delta\varepsilon_p$$

和：

$$\Delta D$$

在三个阶段中的：

$$r_n,\quad r_\tau,\quad r_x$$

以及奇异值衰减曲线。

真正需要回答的是：

> **cycle 46 附近损伤速率发生转折以后，解空间的有效维数是否也出现系统变化？**

如果出现类似：

$$r^{\mathrm{I}}\lt r^{\mathrm{II}}\lt r^{\mathrm{III}}$$

则可以建立：

$$\boxed{\text{fatigue mechanism transition}\rightarrow\text{rank evolution}\rightarrow\text{adaptive PGD enrichment}}$$

这一条非常有价值的物理—数值逻辑链。

---

## 21. 在进行阶段性 SVD 前需要先解决的数据复用问题

本次 FOM 已耗时约 51 min，而张量化仅需 0.035 s。

如果每进行一次新的 SVD 分析、调整阈值或绘图都重新运行 FOM，会造成大量无意义重复计算。

因此下一步工程实现上应先增加：

> **`TowerCyclePhaseSnapshots` 的 `.npz` 保存与读取功能。**

目标是将一次完整 FOM 结果冻结为本地参考数据集。

以后：

- Stage I/II/III SVD；
- 奇异值曲线；
- rank 阈值比较；
- raw/increment 比较；
- 未来 CP/PGD 离线诊断；

都直接读取同一份 FOM 数据，而不再重新进行非线性有限元求解。

整体流程应变为：

$$\text{FOM solve once}\rightarrow\text{save snapshots}\rightarrow\text{repeated offline ROM diagnostics}$$

这样既节约时间，也保证所有后处理使用完全相同的参考解。

---

## 22. 本阶段结论

本阶段已经完成从传统逐时间步非线性塔筒分析向 LATIN-PGD 低秩结构分析的关键过渡。

### 结论 1：100 周期 FOM 提供了可信参考解

$$D_{\max}\approx0.0441$$

$$\max|\varepsilon_p|\approx1.3635\times10^{-3}$$

临界区域仍位于塔底第一单元危险纤维位置。

### 结论 2：$n\times\tau$ 张量化成功显式分离慢循环与快周期坐标

传统：

$$q(t,x)$$

被重新组织为：

$$q(n,\tau,x)$$

从而能够直接分析高周循环问题中特有的双时间尺度结构。

### 结论 3：位移与应力具有极强低秩性

$$\Delta u:(2,2,1)$$

$$\Delta\sigma:(1,1,2)$$

说明全局结构响应高度可压缩。

### 结论 4：累计塑性应变和损伤同样高度低秩

$$\varepsilon_{p,\mathrm{raw}}:(3,1,3)$$

$$D_{\mathrm{raw}}:(3,1,3)$$

说明不可逆状态总体呈现非常规则的低维慢演化。

### 结论 5：真正复杂的是单周期塑性与损伤增量的 cycle-space 演化

$$\Delta\varepsilon_p:(6,3,7)$$

$$\Delta D:(6,2,6)$$

说明主要困难集中在：

$$n\leftrightarrow x$$

而不是 $\tau$ 方向。

因此未来 LATIN-PGD 重点应放在慢循环—空间耦合和模态富集，而不是假定快时间波形完全失去规律性。

### 结论 6：当前结果明确支持继续推进 LATIN-PGD

最复杂的变量仍满足：

$$6\ll100,\qquad3\ll41,\qquad7\ll320$$

因此 FOM 解在三个坐标方向上都表现出显著低维结构。

当前最严谨的英文表述可写为：

> **The 100-cycle full-order tower response exhibits strong empirical multilinear compressibility in the slow-cycle, fast-phase, and spatial directions. The accumulated fields remain highly low-rank, whereas the cycle-wise plastic-strain and damage increments reveal a richer but still low-dimensional cycle-space evolution. These observations provide direct numerical motivation for subsequent LATIN-PGD separation and adaptive enrichment studies.**

---

## 23. 一句话回顾本阶段

> **FOM 是把完整参考答案算出来；$n\times\tau$ 张量化是把高周循环中的慢演化和快周期结构显式整理出来；cycle/phase/space 三方向 SVD 则分别检查慢时间、快时间和空间方向是否可压缩。100 周期结果表明三者均具有明显低秩性，真正需要后续 LATIN-PGD 重点处理的是塑性和损伤在慢循环—空间方向上的逐步演化与模态富集。**

---

## 24. 相关代码、提交与下一步

截至本阶段，直接相关的主要代码包括：

```text
examples/nonlinear_tower_asymmetric_response.py
examples/nonlinear_tower_snapshot_tensor.py
examples/nonlinear_tower_low_rank_diagnostics.py
examples/nonlinear_tower_100cycle_low_rank_probe.py
tests/test_nonlinear_tower_snapshot_tensor.py
tests/test_nonlinear_tower_low_rank_diagnostics.py
```

已完成的关键 Git 提交包括：

```text
d00d07c  feat: add cycle-phase tower snapshot tensorization
cb5e8b4  feat: add full nodal displacement snapshots
b852a7e  feat: add tower low-rank SVD diagnostics
```

其中 `b852a7e` 是本阶段进入实际 100 周期低秩诊断前的稳定基线提交。

下一阶段建议顺序：

1. 增加 `TowerCyclePhaseSnapshots` 的 `.npz` 保存/读取；
2. 冻结一份 100 周期 FOM 参考数据；
3. 分析 Stage I（1–20）、Stage II（21–46）、Stage III（47–100）；
4. 比较三个阶段的 cycle/phase/space rank 和奇异值衰减；
5. 判断是否需要 stage-wise basis 或 adaptive enrichment；
6. 在此基础上再进入真正的 PGD/CP 分离近似与 LATIN-PGD 求解器设计。

---

**阶段状态：** 100 周期 FOM → $n\times\tau$ 张量化 → HOSVD-style SVD 低秩诊断已完成。
**下一阶段：** 保存/读取 FOM 快照数据，并研究不同疲劳演化阶段中的低秩结构变化。
