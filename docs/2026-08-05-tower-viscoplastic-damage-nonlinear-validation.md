# 风机塔筒黏塑性损伤非线性分析与阶段验证

- 日期：2026-08-05
- 项目：Offshore Wind Turbine and LATIN-PGD
- 目标分支：`feature/offshore-wind-turbine-tower-fatigue`
- 文档目的：记录一维黏塑性损伤材料模型向 NREL 5 MW 海上风机塔筒纤维梁柱有限元模型扩展后的阶段性验证内容，包括模型组成、验证方法、关键参数、单调非线性响应、正平均值脉动循环响应、时间离散收敛性、自动化测试结果、当前限制及下一阶段工作。

---

## 1. 阶段目标

本阶段的核心任务是将现有一维循环黏塑性损伤材料模型逐级嵌入海上风机塔筒有限元体系，并建立后续 LATIN-PGD 降阶计算所需的全阶参考模型。

完整扩展层级为：

```text
一维黏塑性损伤材料点
└── 圆环纤维截面
    └── 二维 Euler-Bernoulli 纤维梁柱单元
        └── 锥形风机塔筒有限元模型
            └── 单调加载与循环加载历史
```

本阶段重点验证以下问题：

1. 每根纤维是否能够独立保存材料内部变量；
2. 截面轴力和弯矩是否能够由纤维应力正确积分；
3. 梁柱单元是否能够正确形成内力和一致切线；
4. 全局 Newton-Raphson 求解是否能够处理塑性和损伤演化；
5. 单调加载下能否识别弹性、初始非线性和明显非线性阶段；
6. 脉动循环中能否得到塑性应变、损伤、残余位移和滞回响应；
7. 时间步加密后结果是否稳定收敛；
8. 新增功能是否与已有测试体系兼容。

---

## 2. 当前塔筒模型

### 2.1 几何参数

采用 NREL 5 MW 风机名义塔筒几何：

| 参数 | 数值 |
|---|---:|
| 塔筒高度 | 87.6 m |
| 底部外径 | 6.0 m |
| 顶部外径 | 3.87 m |
| 底部壁厚 | 0.027 m |
| 顶部壁厚 | 0.019 m |

塔筒外径和壁厚均沿高度线性变化。

### 2.2 有限元假定

当前模型采用：

```text
Structural model:
Two-dimensional cantilever tower

Beam theory:
Euler-Bernoulli beam theory

Element type:
Two-node displacement-based beam-column element

Nodal degrees of freedom:
Horizontal displacement
Vertical displacement
Rotation

Boundary condition:
Fully fixed tower base

External loading:
Horizontal concentrated force at tower top

Section model:
Circular annular fiber section

Material model:
Independent uniaxial viscoplastic-damage model at every fiber

Global solution:
Newton-Raphson iteration

State management:
Trial / commit / revert / restart
```

### 2.3 候选正式离散

当前候选正式模型为：

| 项目 | 数值 |
|---|---:|
| 梁单元数 | 40 |
| 每单元 Gauss 点数 | 4 |
| 每截面环向纤维数 | 32 |
| 每截面径向纤维层数 | 2 |
| 每截面纤维总数 | 64 |

对应材料点总数为：

```text
40 elements
× 4 Gauss points per element
× 64 fibers per section
= 10240 uniaxial material points
```

### 2.4 当前非线性验证离散

为控制单周期非线性计算成本，本阶段主要采用：

| 项目 | 数值 |
|---|---:|
| 梁单元数 | 10 |
| 每单元 Gauss 点数 | 2 |
| 每截面环向纤维数 | 16 |
| 每截面径向纤维层数 | 1 |
| 每截面纤维总数 | 16 |

该离散仅用于算法和响应合理性验证，不能替代后续空间离散收敛分析。

---

## 3. 黏塑性损伤材料模型

### 3.1 默认材料参数

当前材料模型默认参数为：

| 参数 | 数值 |
|---|---:|
| 弹性模量 E | 134000 MPa |
| 泊松比 ν | 0.3 |
| 初始屈服应力 σ_y | 80 MPa |
| 各向同性硬化极限 R_∞ | 30 MPa |
| 各向同性硬化参数 γ | 2.0 |
| 随动硬化参数 C | 5500 MPa |
| 非线性随动硬化参数 a | 250 |
| Norton 参数 K | 1220 |
| Norton 指数 n | 2.5 |
| 损伤参数 k_D | 2.778 |
| 损伤指数 n_D | 2.0 |
| 压缩损伤恢复参数 h | 0.2 |
| 损伤上限 | 0.999 |

### 3.2 内部变量

每个纤维材料点独立保存：

$$\mathbf q_f = \left[ \varepsilon_f^p,\, \alpha_f,\, \bar r_f,\, D_f \right].$$

其中：

- $\varepsilon_f^p$：塑性应变；
- $\alpha_f$：随动硬化内部变量；
- $\bar r_f$：各向同性硬化内部变量；
- $D_f$：损伤变量。

### 3.3 塑性演化

黏塑性乘子为：

$$\dot{\lambda} = K^{-n} \left\langle f \right\rangle_+^n.$$

塑性应变率为：

$$\dot{\varepsilon}^p = \frac{\dot{\lambda}}{1-D} \operatorname{sign} \left( \frac{\sigma}{1-D}-\beta \right).$$

因此：

```text
f > 0:
plastic flow continues

f <= 0:
plastic flow stops
```

### 3.4 损伤演化

损伤率为：

$$\dot D = k_D \left\langle Y-Y_0 \right\rangle_+^{n_D}.$$

其中：

$$Y_0 = \frac{\sigma_y^2}{2E}.$$

因此：

```text
Y > Y0:
damage continues to evolve

Y <= Y0:
damage evolution stops
```

这意味着外荷载开始卸载，并不必然导致塑性和损伤立即停止。真正控制内部变量演化的是材料驱动力是否仍高于阈值。

---

## 4. 已完成代码模块

当前已完成并验证的主要模块包括：

| 文件 | 主要作用 |
|---|---|
| `fem/fiber_section.py` | 圆环纤维截面离散、几何积分及弹性截面响应 |
| `fem/beam_column_2d.py` | 二维 Euler-Bernoulli 梁柱单元 |
| `fem/tower_system_2d.py` | 线弹性塔筒组装与静力求解 |
| `fem/tower_response_2d.py` | 塔筒截面与纤维响应提取 |
| `fem/tower_loading.py` | 正平均值脉动塔顶荷载 |
| `fem/viscoplastic_fiber_section.py` | 黏塑性损伤纤维截面 |
| `fem/viscoplastic_beam_column_2d.py` | 黏塑性损伤梁柱单元 |
| `fem/viscoplastic_tower_system_2d.py` | 非线性塔筒组装和 Newton 求解 |
| `fem/viscoplastic_tower_history_2d.py` | 顺序脉动荷载历史求解 |
| `examples/nonlinear_tower_load_probe.py` | 单调非线性荷载水平探测 |
| `examples/nonlinear_tower_pulsating_response.py` | 单周期非线性脉动响应 |
| `tests/test_nonlinear_tower_examples.py` | 非线性示例自动化回归测试 |

---

## 5. 线弹性塔筒基准验证

### 5.1 验证设置

线弹性基准采用：

```text
Tower elements:
40

Gauss points per element:
4

Section fibers:
32 × 2 = 64

Elastic modulus:
210 GPa

Tower-top horizontal force:
1.0 MN
```

### 5.2 关键结果

| 指标 | 结果 |
|---|---:|
| 塔顶水平位移 | 0.7185989 m |
| 塔顶转角 | -0.01442629 rad |
| 塔底水平反力 | -1.0 MN |
| 塔底弯矩 | 8.76×10⁷ N·m |
| 临界积分点最大纤维应力 | 约 116.2381 MPa |

### 5.3 验证结论

线弹性阶段已确认：

1. 塔底水平反力与塔顶水平力满足整体平衡；
2. 塔底弯矩与 $F H$ 一致；
3. 位移和转角沿塔高分布符合悬臂塔筒规律；
4. 最大弯矩和最大纤维应力位于塔底附近；
5. 截面相对两侧纤维应力符号相反且幅值近似对称；
6. 塔筒可变截面几何已在各 Gauss 点处正确更新。

需要注意：

> 线弹性基准采用 $E=210$ GPa，而当前黏塑性损伤材料默认采用 $E=134$ GPa，因此两组塔顶位移不能直接比较。

---

## 6. 单调非线性荷载水平探测

### 6.1 验证目的

在开展正式循环计算之前，需要确定当前材料参数与塔筒模型组合下的非线性起始荷载区间。

采用脚本：

```text
examples/nonlinear_tower_load_probe.py
```

塔顶水平力依次施加为：

$$0.2,\ 0.4,\ 0.6,\ 0.8,\ 1.0,\ 1.2\ \text{MN}.$$

每一级荷载从上一收敛状态继续加载，并记录：

- 塔顶水平位移；
- 最大纤维应力；
- 最大塑性应变；
- 最大损伤；
- Newton 迭代次数。

### 6.2 验证结果

| 塔顶力/MN | 塔顶位移/m | 最大应力/MPa | 最大塑性应变 | 最大损伤 | Newton 次数 |
|---:|---:|---:|---:|---:|---:|
| 0.2 | 0.2274 | 23.30688 | 0 | 0 | 2 |
| 0.4 | 0.4548 | 46.61375 | 0 | 0 | 2 |
| 0.6 | 0.6823 | 69.92063 | 0 | 0 | 2 |
| 0.8 | 0.9101 | 93.03441 | 2.1144×10⁻⁶ | 3.6062×10⁻⁵ | 4 |
| 1.0 | 1.160 | 112.2234 | 6.2449×10⁻⁵ | 8.6878×10⁻⁴ | 5 |
| 1.2 | 1.533 | 125.0199 | 3.1276×10⁻⁴ | 4.1758×10⁻³ | 6 |

### 6.3 结果分析

0.2–0.6 MN 范围内：

```text
Maximum plastic strain:
0

Maximum damage:
0

Response type:
Linear elastic
```

0.8 MN 时首次检测到：

```text
Maximum plastic strain:
2.1144e-6

Maximum damage:
3.6062e-5
```

因此当前非线性起点位于：

$$0.6\ \text{MN} < F_{\mathrm{nl}} \le 0.8\ \text{MN}.$$

按照弹性阶段应力线性外推，达到 $\sigma_y=80$ MPa 时的塔顶力约为：

$$F_y \approx 0.6 \frac{80}{69.92063} = 0.687\ \text{MN}.$$

但当前模型属于率相关黏塑性模型，因此严格的可观测塑性起点还与加载持续时间及时间增量有关。

### 6.4 结论

- 0.6 MN 及以下可视为弹性范围；
- 0.8 MN 为初始非线性区；
- 1.0 MN 能产生清晰但仍可稳定收敛的塑性和损伤；
- 1.2 MN 非线性显著增强；
- 后续单周期验证采用 $F_{\max}=1.0$ MN 较为合适。

---

## 7. 正平均值脉动循环分析

### 7.1 荷载定义

脉动塔顶荷载定义为：

$$F(t) = 0.55F_{\max} + 0.45F_{\max} \sin \left( \frac{2\pi t}{T} \right).$$

采用：

| 参数 | 数值 |
|---|---:|
| 最大荷载 F_max | 1.0 MN |
| 最小荷载 F_min | 0.1 MN |
| 力比 R_F | 0.1 |
| 平均荷载 | 0.55 MN |
| 荷载幅值 | 0.45 MN |
| 周期 T | 10 |
| 循环数 | 1 |

该荷载始终为正，因此属于正平均值脉动加载，而不是零均值反向加载。

### 7.2 40 步/周期结果

| t/T | 塔顶力/MN | 塔顶位移/m | 最大应力/MPa | 最大塑性应变 | 最大损伤 | Newton 次数 |
|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.55 | 0.62540 | 64.09391 | 0 | 0 | 2 |
| 0.25 | 1.00 | 1.1817 | 109.8167 | 1.0540×10⁻⁴ | 1.4702×10⁻³ | 3 |
| 0.50 | 0.55 | 0.70140 | 56.17343 | 1.6573×10⁻⁴ | 2.3806×10⁻³ | 2 |
| 0.75 | 0.10 | 0.18951 | 16.58532 | 1.6573×10⁻⁴ | 2.3806×10⁻³ | 2 |
| 1.00 | 0.55 | 0.70140 | 56.17343 | 1.6573×10⁻⁴ | 2.3806×10⁻³ | 2 |

### 7.3 临界位置

临界纤维位置为：

| 项目 | 结果 |
|---|---:|
| 临界单元 | 第 1 单元 |
| 临界 Gauss 点 | 第 1 点 |
| 临界纤维 | 第 5 根 |
| 临界高度 | 1.851206 m |
| 临界纤维 y 坐标 | 2.945089 m |

该位置位于塔筒底部附近的截面外缘，与悬臂塔筒底部弯矩最大、截面外缘应力最大的规律一致。

---

## 8. 非线性循环响应分析

### 8.1 相同荷载下的残余位移

在初始平均荷载：

$$F=0.55\ \text{MN}$$

首次加载时：

$$u_{\mathrm{top}}^{(0)} = 0.62540\ \text{m}.$$

经历峰值荷载后，再次回到相同平均荷载时：

$$u_{\mathrm{top}}^{(1)} = 0.70140\ \text{m}.$$

位移偏移为：

$$\Delta u_{\mathrm{top}} = 0.70140-0.62540 = 0.07600\ \text{m}.$$

这说明塔筒经历峰值非线性后产生了不可恢复变形。

### 8.2 首周期开放滞回路径

当前仅计算第一个脉动周期，初始状态为未加载、无塑性和无损伤状态。

周期结束时虽然外荷载返回初始平均值，但内部变量已经发生变化，因此：

```text
Initial force:
0.55 MN

Final force:
0.55 MN

Initial internal state:
Undamaged and unplasticized

Final internal state:
Nonzero plastic strain and damage
```

所以首周期荷载—位移路径不闭合是合理的，不能将其理解为已经形成稳定滞回环。

### 8.3 塑性应变演化

最大塑性应变大约从：

$$t/T \approx 0.10$$

开始增加。

达到峰值荷载：

$$t/T = 0.25$$

后，塑性应变仍继续增长一段时间，并大约在：

$$t/T \approx 0.40$$

进入平台。

对应物理过程为：

```text
t/T ≈ 0.10:
Force increases to approximately 0.8 MN
Plasticity begins

t/T = 0.25:
Force reaches 1.0 MN
Maximum loading point

0.25 < t/T < approximately 0.40:
External force decreases
But the viscoplastic yield function remains positive at some fibers

t/T ≈ 0.40:
Force decreases to approximately 0.8 MN
Plastic flow stops

t/T > approximately 0.40:
Mainly elastic unloading and reloading
```

### 8.4 损伤演化

损伤发展趋势与塑性应变相近：

1. 荷载升至非线性区间后开始增长；
2. 峰值后的初始卸载阶段继续增长；
3. 荷载进一步下降后进入平台。

这不是因为整个卸载阶段持续损伤，而是因为卸载初期仍可能满足：

$$Y-Y_0>0.$$

只有当能量释放率降至阈值以下时：

$$Y-Y_0 \le 0,$$

损伤率才变为零。

### 8.5 临界纤维应力—应变路径

固定临界纤维的应力—应变响应表现为：

1. 初始加载阶段近似线性；
2. 高应力区出现切线刚度下降；
3. 峰值附近产生黏塑性流动和损伤；
4. 卸载后留下明显水平应变偏移；
5. 低应力区的卸载与后续再加载路径近似重合。

同一应力水平下的应变偏移约为：

$$1.6\times10^{-4},$$

与最终最大塑性应变：

$$1.6573\times10^{-4}$$

基本一致。

---

## 9. 时间离散收敛性

### 9.1 验证方法

保持以下条件不变：

```text
Tower elements:
10

Gauss points:
2 per element

Section fibers:
16 × 1

Maximum force:
1.0 MN

Force ratio:
0.1

Period:
10
```

分别采用：

```text
20 increments per cycle
40 increments per cycle
80 increments per cycle
```

### 9.2 结果对比

| 指标 | 20 步 | 40 步 | 80 步 |
|---|---:|---:|---:|
| 峰值塔顶位移/m | 1.1807 | 1.1817 | 1.1819 |
| 周期末塔顶位移/m | 0.69944 | 0.70140 | 0.70191 |
| 峰值最大应力/MPa | 109.9008 | 109.8167 | 109.7950 |
| 峰值最大塑性应变 | 1.0372×10⁻⁴ | 1.0540×10⁻⁴ | 1.0583×10⁻⁴ |
| 峰值最大损伤 | 1.4485×10⁻³ | 1.4702×10⁻³ | 1.4759×10⁻³ |
| 最终最大塑性应变 | 1.624876×10⁻⁴ | 1.657313×10⁻⁴ | 1.665741×10⁻⁴ |
| 最终最大损伤 | 2.337137×10⁻³ | 2.380638×10⁻³ | 2.392074×10⁻³ |

### 9.3 20 步与 40 步差异

| 指标 | 相对差异 |
|---|---:|
| 峰值塔顶位移 | 约 0.085% |
| 周期末塔顶位移 | 约 0.279% |
| 峰值最大应力 | 约 0.077% |
| 峰值最大塑性应变 | 约 1.59% |
| 峰值最大损伤 | 约 1.48% |
| 最终最大塑性应变 | 约 1.96% |
| 最终最大损伤 | 约 1.83% |

### 9.4 40 步与 80 步差异

| 指标 | 相对差异 |
|---|---:|
| 峰值塔顶位移 | 约 0.02% |
| 周期末塔顶位移 | 约 0.07% |
| 峰值最大应力 | 约 0.02% |
| 峰值最大塑性应变 | 约 0.41% |
| 峰值最大损伤 | 约 0.39% |
| 最终最大塑性应变 | 约 0.51% |
| 最终最大损伤 | 约 0.48% |

### 9.5 时间离散结论

当前单周期、中等空间离散模型中，40 步/周期能够在精度和计算效率之间取得较好平衡：

```text
Global displacement error:
Less than approximately 0.1%

Maximum stress error:
Approximately 0.02%

Plastic strain error:
Approximately 0.4% to 0.5%

Damage error:
Approximately 0.4% to 0.5%
```

因此建议后续常规非线性循环分析采用：

$$\boxed{ 40\ \text{increments per cycle} }$$

该结论只适用于当前 10 单元、2 Gauss 点和 16 根纤维模型。正式空间离散确定后仍需重新检查时间步收敛性。

---

## 10. Newton 收敛性

40 步/周期分析中：

| 指标 | 结果 |
|---|---:|
| 单步最大 Newton 迭代次数 | 4 |
| 最大自由自由度残差范数 | 8.684860×10⁻³ N |
| 外荷载量级 | 10⁶ N |

对应相对残差量级约为：

$$8.7\times10^{-9}.$$

80 步/周期分析中：

| 指标 | 结果 |
|---|---:|
| 单步最大 Newton 迭代次数 | 3 |
| 最大自由自由度残差范数 | 2.562738×10⁻⁴ N |

说明当前全局非线性求解器在单周期中具有稳定的收敛性能。

---

## 11. 自动化回归测试

### 11.1 新增测试

新增文件：

```text
tests/test_nonlinear_tower_examples.py
```

主要覆盖：

1. 单调探测荷载序列的 MN 到 N 转换；
2. 非法荷载序列检查；
3. 临界纤维优先按最终损伤选择；
4. 无损伤时按最终塑性应变选择；
5. 非线性循环历史数组形状；
6. 所有结果是否为有限值；
7. 循环中是否产生非零塑性应变；
8. 循环中是否产生非零损伤；
9. 相同荷载返回时是否存在残余位移；
10. 临界纤维是否位于塔底附近；
11. Newton 迭代是否满足收敛要求。

### 11.2 定向测试结果

```text
Ran 9 tests in 1.513s

OK
```

### 11.3 完整测试集结果

```text
Ran 110 tests in 161.578s

OK
```

这说明新增非线性单调探测、脉动循环示例和回归测试没有破坏已有材料点、纤维截面、梁柱单元、塔筒系统和历史求解功能。

---

## 12. 当前阶段已验证内容

截至本阶段，已经完成以下验证闭环：

### 12.1 材料点层级

- 黏塑性流动方程已实现；
- 各向同性硬化已实现；
- 随动硬化已实现；
- 单侧损伤影响已实现；
- 内部变量 RK4 时间积分已实现；
- 状态 trial / commit / revert 已实现。

### 12.2 截面层级

- 圆环截面能够离散为独立单轴纤维；
- 每根纤维具有独立内部变量；
- 轴力和弯矩能够由纤维应力积分；
- 截面切线能够通过数值扰动形成；
- 重复 trial 不会错误累积塑性和损伤。

### 12.3 单元层级

- 二维 Euler-Bernoulli 梁柱单元能够调用多个截面；
- Gauss 点处可采用局部外径和壁厚；
- 单元内力和切线可用于全局 Newton 迭代；
- 单元状态可统一提交和回滚。

### 12.4 结构层级

- 锥形塔筒线弹性响应满足整体平衡；
- 非线性单调加载能够识别塑性和损伤起点；
- 正平均值脉动加载能够产生残余位移；
- 塑性应变和损伤历史具有合理阈值特征；
- 临界位置稳定出现在塔底截面外缘；
- Newton 迭代稳定收敛；
- 时间步加密后结果稳定逼近。

---

## 13. 当前阶段结论

本阶段可以得到以下结论：

1. 一维黏塑性损伤材料模型已经成功嵌入二维纤维梁柱塔筒模型；
2. 纤维内部变量能够在材料、截面、单元和全局层级正确传递；
3. 全局 Newton-Raphson 求解能够处理塑性与损伤导致的刚度变化；
4. 单调加载下非线性起点位于约 0.6–0.8 MN；
5. $F_{\max}=1.0$ MN、$R_F=0.1$ 的脉动循环可以产生清晰但可控的非线性；
6. 首周期中出现塑性应变、损伤、残余位移和开放滞回路径；
7. 峰值后的卸载初期仍存在少量塑性和损伤演化，与率相关材料方程一致；
8. 40 步/周期可作为当前模型后续常规循环分析的推荐时间离散；
9. 当前全部 110 项测试通过，代码基础稳定。

因此，当前全阶塔筒模型已经具备继续开展反向循环、多循环累积和 LATIN-PGD 降阶研究的基础。

---

## 14. 当前限制

当前结果仍存在以下限制：

1. 当前材料参数来自已有一维参考模型，并未针对风机塔筒钢材标定；
2. 当前弹性模量 $E=134$ GPa，低于常规结构钢约 210 GPa；
3. 当前屈服应力 $\sigma_y=80$ MPa，低于常见塔筒钢材；
4. 当前非线性验证采用 10 单元、2 Gauss 点和 16 根纤维的较粗模型；
5. 当前仅验证一个正平均值脉动周期；
6. 当前塔顶力始终为正，固定纤维尚未经历应力反向；
7. 尚未验证包辛格效应和拉压交替塑性；
8. 尚未分析多循环棘轮效应；
9. 尚未分析损伤随循环数的长期累积；
10. 尚未与 OpenSees、商业有限元或试验结果进行交叉验证；
11. 尚未开展 LATIN-PGD 降阶计算；
12. 当前结果不能直接作为真实风机塔筒疲劳寿命预测结果。

现阶段结果应定义为：

> 黏塑性损伤塔筒有限元模型的数值实现验证与响应合理性验证。

---

## 15. 下一阶段工作

### 15.1 固化当前阶段代码

本阶段待统一纳入版本管理的文件包括：

```text
examples/nonlinear_tower_load_probe.py
examples/nonlinear_tower_pulsating_response.py
tests/test_nonlinear_tower_examples.py
docs/2026-08-05-tower-viscoplastic-damage-nonlinear-validation.md
```

### 15.2 建立零均值反向循环

建议采用：

$$F(t) = F_a \sin \left( \frac{2\pi t}{T} \right).$$

主要验证：

- 固定纤维应力是否反向；
- 拉压交替塑性是否正确；
- 随动硬化是否产生包辛格效应；
- 正负方向滞回路径是否合理；
- 材料和结构响应是否具有预期对称性。

### 15.3 开展多循环正平均值脉动分析

重点记录：

- 每周期最大和最小塔顶位移；
- 周期末残余位移；
- 塑性应变棘轮累积；
- 最大损伤随循环数变化；
- 临界纤维位置是否迁移；
- 滞回路径是否逐步稳定；
- Newton 迭代次数是否随损伤增长。

### 15.4 空间离散收敛

建议比较：

```text
Element counts:
10, 20, 40

Gauss points per element:
2, 4

Section fibers:
16, 32, 64, 256
```

重点比较：

- 塔顶位移；
- 塔底曲率；
- 最大纤维应力；
- 最大塑性应变；
- 最大损伤；
- 临界位置；
- 单周期计算时间；
- 内存占用。

### 15.5 材料参数标定

后续应采用适用于风机塔筒钢材的：

- 弹性模量；
- 屈服应力；
- 各向同性硬化参数；
- 随动硬化参数；
- 黏塑性参数；
- 损伤参数；
- 循环试验或文献标定数据。

### 15.6 与 LATIN-PGD 对接

全阶模型进一步稳定后，需要开展：

1. 时间函数与空间模式分离；
2. 全局阶段与局部阶段构造；
3. 材料内部变量的可分离表达；
4. 多循环时间域压缩；
5. 与逐时间步全阶解对比；
6. 误差、计算时间和内存占用评价。

---

## 16. 复现命令

### 16.1 单调非线性探测

```powershell
python -m examples.nonlinear_tower_load_probe
```

### 16.2 20 步单周期分析

```powershell
python -m examples.nonlinear_tower_pulsating_response --no-plot --increments 20
```

### 16.3 40 步单周期分析

```powershell
python -m examples.nonlinear_tower_pulsating_response --no-plot --increments 40
```

### 16.4 80 步单周期分析

```powershell
python -m examples.nonlinear_tower_pulsating_response --no-plot --increments 80
```

### 16.5 定向回归测试

```powershell
python -m unittest discover -s tests -p "test_nonlinear_tower_examples.py" -v
```

### 16.6 完整测试集

```powershell
python -m unittest discover -s tests -p "test_*.py" -v
```

---

## 17. 相关提交记录

当前相关提交包括：

```text
d755067 feat: add elastic tower section response example
7035f8c feat: add pulsating tower top force history
446f359 feat: add elastic tower pulsating response
5bd4b81 feat: add viscoplastic damage fiber section
9efdb52 feat: add viscoplastic damage beam element
4253397 feat: add nonlinear viscoplastic tower solver
707430f feat: add nonlinear pulsating tower history solver
```

本阶段新增的单调荷载探测示例、非线性脉动响应示例、回归测试和本总结文档尚待提交。
