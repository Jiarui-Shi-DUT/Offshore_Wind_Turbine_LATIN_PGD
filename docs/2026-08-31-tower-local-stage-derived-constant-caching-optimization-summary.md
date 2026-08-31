# 2026-08-31 塔筒 LATIN-PGD 局部本构固定材料派生量缓存优化阶段总结

## 1. 阶段定位

本阶段是在 **OPT-1：局部本构阶段标量损伤限幅优化与多周期等价性验证** 通过之后开展的第二轮实现级性能优化。

OPT-1 已经证明：当前塔筒 LATIN-PGD 在 10 周期算例中的低效率，存在显著的 Python 标量通用函数调用开销。仅将局部本构热循环中的标量 `np.clip` 改写为等价的标量条件限幅，就使 10 周期 LATIN 总时间由约 467.36 s 降至 223.16 s，并且没有改变 LATIN 外层迭代次数、PGD 模式数和最终收敛指标。

在此基础上，本阶段继续遵循同一原则：

> 不改变 LATIN-PGD 数学算法，不改变材料本构关系，不改变四阶龙格-库塔时间积分，不改变空间离散和收敛判据，只消除本构热循环中对固定材料派生量的高频重复求值。

本阶段名称：

**OPT-2：局部本构固定材料派生量缓存优化与多周期等价性验证**

阶段结论：

**PASS**

---

## 2. 当前性能优化分支与基准

性能优化分支：

`perf/tower-local-stage-optimization`

OPT-1 代码 checkpoint：

`c27af15`

OPT-1 阶段总结 checkpoint：

`55d578e`

本阶段 OPT-2 代码 checkpoint：

`e2d15b9`

提交信息：

`perf: cache material derived constants`

本阶段仍采用与前面一致的塔筒效率基准：

- 塔筒单元数：10
- 每单元高斯点数：2
- 每高斯点周向纤维数：16
- 径向层数：1
- 总材料点数：320
- 非对称循环荷载：
  - 最大荷载：$+1.0$ MN
  - 最小荷载：$-0.5$ MN
  - 应力比：$R=-0.5$
- 周期：10 s
- 每周期：40 个时间区间
- LATIN 收敛容限：$10^{-5}$
- 空间策略：`residual_ls`

1 周期：

$$
N_t=41
$$

10 周期：

$$
N_t=401
$$

---

## 3. OPT-1 后仍然存在的高频固定量重复计算

前一阶段函数级耗时分析已经确定，在 10 周期稳定求解路径中：

- LATIN 外层迭代：39 次
- 时间区间：400
- 材料点：320

因此每次局部阶段沿完整材料历史推进所形成的材料状态时间推进总次数为：

$$
39\times400\times320
=
4,992,000
$$

即约 **499 万次材料状态时间推进**。

每一步采用经典四阶龙格-库塔法，每个时间步需要 4 次本构速率计算：

$$
4\times4,992,000
=
19,968,000
$$

此外，每个正式时间点还会重新计算并保存一次材料速率：

$$
39\times401\times320
=
5,004,480
$$

因此 10 周期稳定路径中的局部本构速率函数总调用次数为：

$$
19,968,000+5,004,480
=
24,972,480
$$

即接近 **2500 万次本构速率评价**。

OPT-1 之后，标量损伤限幅的主要 NumPy 通用数组开销已经被消除，但源码检查表明，在这约 2497 万次 `local_rates_from_forces()` 调用中，仍然存在若干实际上与时间、材料状态和 LATIN 迭代无关的固定量被重复求值。

---

## 4. 源码确认：哪些量实际上是常数

### 4.1 `MaterialParameters` 的数据结构

当前材料参数类定义为：

```python
@dataclass(frozen=True)
class MaterialParameters:
    ...
```

`frozen=True` 表明一个 `MaterialParameters` 对象创建之后，其基础材料参数不会在正常计算过程中被修改。

因此，只依赖基础材料参数的派生量，在同一个材料参数对象生命周期中理论上也是固定的。

---

### 4.2 Norton 黏塑性系数

原实现：

```python
@property
def k_viscoplastic(self) -> float:
    return float(self.K ** (-self.n))
```

对应：

$$
k=K^{-n}
$$

其中 $K$ 和 $n$ 都是固定材料参数。

原来的普通 `@property` 不会缓存结果，因此每次访问：

```python
material.k_viscoplastic
```

都会重新执行一次幂运算。

在 10 周期局部本构速率路径中，该属性至少会被访问约：

$$
24,972,480
$$

次。

---

### 4.3 损伤阈值

原实现：

```python
@property
def Y0(self) -> float:
    return float(self.sigma_y**2 / (2.0 * self.E))
```

对应：

$$
Y_0=\frac{\sigma_y^2}{2E}
$$

其中 $\sigma_y$ 和 $E$ 都是固定材料参数。

同样，普通 `@property` 不缓存结果，因此该量在局部本构热循环中也被重复求值约：

$$
24,972,480
$$

次。

---

### 4.4 $\sqrt{\gamma}$

局部本构阶段原来有两处：

```python
np.sqrt(material.gamma)
```

其一用于从变换后的各向同性硬化力恢复物理硬化力，其二用于内部变量演化速率。

由于 $\gamma$ 是固定材料参数：

$$
\sqrt{\gamma}
$$

也是固定量。

在 10 周期局部本构路径中，两处合计理论调用规模约为：

$$
2\times24,972,480
=
49,944,960
$$

次。

---

### 4.5 浮点机器精度

原实现：

```python
np.finfo(float).eps
```

用于判断相对有效应力是否足够接近零。

该值只由浮点数据类型决定，与材料状态和时间完全无关。

因此：

```python
np.finfo(float).eps
```

没有必要在每一次本构速率计算中重新通过 NumPy 接口取得。

10 周期局部本构路径中的调用规模约为：

$$
24,972,480
$$

次。

---

## 5. 重复求值规模

仅考虑 `local_rates_from_forces()` 热循环中的四类固定量：

- `material.k_viscoplastic`：约 2497 万次访问
- `material.Y0`：约 2497 万次访问
- `np.sqrt(material.gamma)`：两处，约 4994 万次
- `np.finfo(float).eps`：约 2497 万次

因此对应的高频固定量求值/访问规模约为：

$$
24,972,480
+
24,972,480
+
49,944,960
+
24,972,480
=
124,862,400
$$

即超过 **1.24 亿次**。

需要强调：

> 这并不意味着 1.24 亿次操作都具有相同成本，也不能仅根据调用次数直接推导实际耗时比例。

但该数量级足以说明：在 Python 标量热循环中，把固定材料量反复通过属性函数、幂运算、开方和 NumPy 通用接口重新求值，是一个明确且可控的实现级优化对象。

---

## 6. 本阶段优化原则

本阶段保持以下内容完全不变：

- 黏塑性损伤本构方程
- Norton 黏塑性规律
- 屈服函数
- 塑性应变演化
- 运动硬化演化
- 各向同性硬化演化
- 损伤演化
- 单边弹性关系
- 四阶龙格-库塔积分
- 材料状态变量定义
- 时间步长
- 时间离散
- 材料点数量
- LATIN 搜索方向
- PGD 模式富集
- 时间函数更新
- 空间最小二乘求解策略
- 收敛容限
- 模式接受逻辑

本阶段只改变：

> 固定材料派生量和固定数值常量在 Python 中的求值位置与缓存方式。

因此这一阶段仍属于：

**实现级等价优化**

而不是：

**算法级修改**

---

## 7. 代码修改

本阶段修改两个文件：

```text
latin/local_stage.py
material/viscoplastic_damage_1d.py
```

---

### 7.1 `k_viscoplastic` 改为缓存属性

原实现：

```python
@property
def k_viscoplastic(self) -> float:
    return float(self.K ** (-self.n))
```

修改为：

```python
@cached_property
def k_viscoplastic(self) -> float:
    return float(self.K ** (-self.n))
```

外部调用方式不变：

```python
material.k_viscoplastic
```

区别仅在于：

- 第一次访问：计算 $K^{-n}$
- 后续访问：直接读取已经缓存的结果

---

### 7.2 `Y0` 改为缓存属性

原实现：

```python
@property
def Y0(self) -> float:
    return float(self.sigma_y**2 / (2.0 * self.E))
```

修改为：

```python
@cached_property
def Y0(self) -> float:
    return float(self.sigma_y**2 / (2.0 * self.E))
```

外部调用方式仍保持：

```python
material.Y0
```

因此现有调用代码不需要改变。

---

### 7.3 新增 `sqrt_gamma`

新增：

```python
@cached_property
def sqrt_gamma(self) -> float:
    return float(np.sqrt(self.gamma))
```

随后将相关位置：

```python
np.sqrt(material.gamma)
```

替换为：

```python
material.sqrt_gamma
```

在同一个 `MaterialParameters` 对象上，$\sqrt{\gamma}$ 只在首次访问时计算一次。

---

### 7.4 机器精度改为模块级常量

在 `latin/local_stage.py` 模块加载阶段计算：

```python
_FLOAT_EPS = float(np.finfo(float).eps)
```

本构热循环中原来的：

```python
if abs(relative_effective_stress) <= np.finfo(float).eps:
```

改为：

```python
if abs(relative_effective_stress) <= _FLOAT_EPS:
```

数值值本身不变。

---

## 8. 为什么采用 `cached_property`

当前 Python 版本：

```text
Python 3.8.20
```

支持标准库：

```python
from functools import cached_property
```

而 `MaterialParameters` 是：

```python
@dataclass(frozen=True)
```

这意味着基础材料参数不会在正常计算中被重新赋值。

因此，对于完全由固定基础材料参数决定的：

$$
K^{-n}
$$

$$
\frac{\sigma_y^2}{2E}
$$

$$
\sqrt{\gamma}
$$

采用首次访问计算、后续缓存复用的方式，与原公式在数学上保持一致。

此外，仓库级搜索确认 `k_viscoplastic` 和 `Y0` 的现有使用均为只读访问，没有代码尝试对这两个属性赋值，因此现有调用接口保持兼容。

---

## 9. 完整仓库回归测试

代码修改后执行：

```powershell
python -m pytest -q
```

结果：

```text
313 passed, 10 warnings in 105.83s
```

313 个测试全部通过。

10 个 warning 仍然来自 Matplotlib / setuptools 中 `distutils` 版本比较的弃用提示，与本阶段修改无关。

因此，从现有仓库测试覆盖范围看，本阶段修改没有引入数值或接口回归。

需要注意：

> 完整测试总耗时从前一次约 123 s 变为本次约 106 s，不能直接作为 OPT-2 的正式性能收益，因为完整测试受操作系统调度、文件缓存、CPU 状态和测试路径组成等因素影响。

正式性能判断仍以相同塔筒基准程序为准。

---

## 10. 1 周期塔筒等价性与性能验证

执行：

```powershell
python .\tower_asymmetric_efficiency_scaling_pilot.py --cycles 1
```

### 10.1 FOM

本次 FOM：

- 计算时间：

$$
34.605644\ \mathrm{s}
$$

- 最大牛顿迭代次数：4
- 最大塔顶位移绝对值：

$$
1.173015956970\ \mathrm{m}
$$

---

### 10.2 LATIN-PGD 数值路径

本次结果：

- 时间点：41
- 材料点：320
- 外层迭代：18
- attempted：18
- 试算次数：29
- PGD 模式数：11
- 新增模式数：11
- 最终收敛指标：

$$
\xi
=
7.918424536257\times10^{-6}
$$

与 OPT-1 的 1 周期稳定路径完全一致。

---

### 10.3 LATIN-PGD 时间

本次：

- 建模时间：

$$
0.122860\ \mathrm{s}
$$

- 求解器时间：

$$
9.310437\ \mathrm{s}
$$

- LATIN 总时间：

$$
T_{\mathrm{OPT2,1}}
=
9.433297\ \mathrm{s}
$$

同期：

$$
\frac{T_{\mathrm{FOM}}}{T_{\mathrm{LATIN}}}
=
3.668457
$$

---

### 10.4 相对于 OPT-1

OPT-1 的 1 周期 LATIN 总时间：

$$
13.692853\ \mathrm{s}
$$

OPT-2：

$$
9.433297\ \mathrm{s}
$$

减少：

$$
13.692853-9.433297
=
4.259556\ \mathrm{s}
$$

相对降幅：

$$
\frac{13.692853-9.433297}{13.692853}
\times100\%
\approx
31.1\%
$$

对应增量加速：

$$
\frac{13.692853}{9.433297}
\approx
1.45
$$

倍。

---

### 10.5 相对于最初未优化实现

最初 1 周期 LATIN 总时间约：

$$
23.726\ \mathrm{s}
$$

当前：

$$
9.433297\ \mathrm{s}
$$

累计降幅约：

$$
60.2\%
$$

累计加速约：

$$
2.52
$$

倍。

---

## 11. 10 周期塔筒等价性与性能验证

执行：

```powershell
python .\tower_asymmetric_efficiency_scaling_pilot.py --cycles 10
```

---

### 11.1 FOM

本次 FOM：

- 计算时间：

$$
294.307999\ \mathrm{s}
$$

- 最大牛顿迭代次数：4
- 最大塔顶位移绝对值：

$$
1.426595421862\ \mathrm{m}
$$

---

### 11.2 LATIN-PGD 数值路径

本次：

- 时间点：401
- 材料点：320
- 外层迭代：39
- attempted：39
- 试算次数：60
- PGD 模式数：21
- 新增模式数：21
- 最终收敛指标：

$$
\xi
=
8.941607234831\times10^{-6}
$$

该数值路径与 OPT-1 以及优化前的 10 周期稳定基准完全一致。

因此，本阶段修改没有改变：

- 外层收敛迭代数
- 试算次数
- PGD 模式数
- 模式富集数量
- 最终收敛指标

---

### 11.3 LATIN-PGD 时间

本次：

- 建模时间：

$$
1.011535\ \mathrm{s}
$$

- 求解器时间：

$$
173.706005\ \mathrm{s}
$$

- LATIN 总时间：

$$
T_{\mathrm{OPT2,10}}
=
174.717540\ \mathrm{s}
$$

同期 FOM/LATIN 比值：

$$
\frac{T_{\mathrm{FOM}}}{T_{\mathrm{LATIN}}}
=
1.684479
$$

因此当前 10 周期基准下，LATIN-PGD 明显快于 FOM。

---

## 12. OPT-2 相对于 OPT-1 的增量性能收益

OPT-1 的 10 周期 LATIN 总时间：

$$
223.160690\ \mathrm{s}
$$

OPT-2：

$$
174.717540\ \mathrm{s}
$$

减少：

$$
223.160690-174.717540
=
48.443150\ \mathrm{s}
$$

相对 OPT-1 降幅：

$$
\frac{223.160690-174.717540}{223.160690}
\times100\%
\approx
21.7\%
$$

增量加速约：

$$
\frac{223.160690}{174.717540}
\approx
1.28
$$

倍。

这说明，在 OPT-1 已经显著降低局部本构开销之后，仅进一步消除固定材料派生量的高频重复求值，10 周期仍然能够获得明显的额外性能收益。

---

## 13. 从原始实现到 OPT-2 的累计变化

### 13.1 1 周期

| 阶段 | FOM 时间 / s | LATIN 总时间 / s | FOM/LATIN |
|---|---:|---:|---:|
| 原始实现 | 约 37.142 | 约 23.726 | 约 1.565 |
| OPT-1 | 36.418586 | 13.692853 | 2.659679 |
| OPT-2 | 34.605644 | 9.433297 | 3.668457 |

原始 LATIN 到 OPT-2：

$$
23.726
\rightarrow
9.433297\ \mathrm{s}
$$

累计时间下降约：

$$
60.2\%
$$

累计加速约：

$$
2.52
$$

倍。

---

### 13.2 10 周期

| 阶段 | FOM 时间 / s | LATIN 总时间 / s | FOM/LATIN |
|---|---:|---:|---:|
| 原始实现 | 约 307.891 | 约 467.362 | 约 0.659 |
| OPT-1 | 304.466891 | 223.160690 | 1.364339 |
| OPT-2 | 294.307999 | 174.717540 | 1.684479 |

原始 LATIN 到 OPT-2：

$$
467.362
\rightarrow
174.717540\ \mathrm{s}
$$

累计减少：

$$
292.644460\ \mathrm{s}
$$

累计时间下降约：

$$
62.6\%
$$

累计加速约：

$$
2.67
$$

倍。

---

## 14. 如何理解不同测次中的 FOM 时间变化

需要严格区分代码优化收益与机器运行状态波动。

OPT-1 的 10 周期 FOM：

$$
304.466891\ \mathrm{s}
$$

OPT-2 测次中的 FOM：

$$
294.307999\ \mathrm{s}
$$

两次相差约：

$$
3.3\%
$$

说明当前桌面计算环境中的实际运行时间存在数个百分点的自然波动。

因此不能把 LATIN 从 223.16 s 到 174.72 s 的全部差值机械地解释为精确的、无误差的优化收益。

但是：

- FOM 变化约 3.3%
- LATIN 变化约 21.7%
- LATIN 数值路径完全一致
- 修改对象正位于此前确认的高频本构热循环中

因此可以合理确认：

> OPT-2 带来了真实且显著的性能改善，其幅度明显超过本次 FOM 所反映的机器状态波动。

---

## 15. 本阶段最重要的技术认识

### 15.1 低效率仍然主要受到实现细节影响

OPT-1 已经证明，标量 `np.clip` 的通用数组调用成本足以严重影响整个 10 周期 LATIN-PGD 的效率。

OPT-2 进一步证明，即使这些限幅开销已经消除，数千万次本构速率评价中对固定派生材料量的重复计算，仍然能够造成可观的额外成本。

因此当前结果持续支持：

> 塔筒 LATIN-PGD 此前观察到的性能劣化，至少有很大一部分来自 Python 标量热循环的实现方式，而不是 LATIN-PGD 数学结构自身必然具有的成本。

---

### 15.2 “公式很简单”不意味着高频实现成本可以忽略

例如：

$$
Y_0=\frac{\sigma_y^2}{2E}
$$

或者：

$$
k=K^{-n}
$$

单独计算一次几乎可以忽略。

但是当它们进入约 2497 万次本构速率评价的 Python 热循环之后，即使每次调用成本很小，也会累计成为可观的总成本。

这说明，对于长时程非线性降阶方法：

> 性能判断必须同时考虑算法复杂度和热点函数的实际实现方式。

---

### 15.3 当前 10 周期效率结论发生了根本变化

最初：

$$
T_{\mathrm{LATIN}}
\approx
467.36\ \mathrm{s}
>
T_{\mathrm{FOM}}
\approx
307.89\ \mathrm{s}
$$

当时看起来 LATIN-PGD 在 10 周期已经失去效率优势。

OPT-1 后：

$$
T_{\mathrm{LATIN}}
\approx
223.16\ \mathrm{s}
<
T_{\mathrm{FOM}}
\approx
304.47\ \mathrm{s}
$$

OPT-2 后：

$$
T_{\mathrm{LATIN}}
\approx
174.72\ \mathrm{s}
<
T_{\mathrm{FOM}}
\approx
294.31\ \mathrm{s}
$$

因此原先所谓“2 到 5 周期之间出现效率交叉”不能继续作为 LATIN-PGD 方法自身的性能结论。

更准确的表述应该是：

> 原始 Python 实现曾在 2 到 5 周期之间出现 FOM/LATIN 计算效率交叉，但经过两个完全不改变数学算法的局部本构实现级优化后，该交叉现象已不再成立。因此此前的交叉主要反映原始实现开销，而不能直接视为 LATIN-PGD 方法本身的固有效率边界。

---

## 16. 当前可以确认的事实

本阶段已经通过实际代码和数值测试确认：

1. `MaterialParameters` 为冻结数据类，基础材料参数在正常计算中保持不变。
2. 原来的 `k_viscoplastic` 和 `Y0` 是普通 `@property`，每次访问都会重新执行对应公式。
3. `k_viscoplastic`、`Y0` 和 $\sqrt{\gamma}$ 均只依赖固定材料参数。
4. `np.finfo(float).eps` 与材料状态和时间无关。
5. 当前 Python 3.8.20 支持 `functools.cached_property`。
6. 仓库现有调用对 `k_viscoplastic` 和 `Y0` 均为只读访问。
7. 修改后 313 个仓库测试全部通过。
8. 1 周期 LATIN 求解路径与 OPT-1 完全一致。
9. 10 周期 LATIN 求解路径与 OPT-1 完全一致。
10. OPT-2 的 1 周期 LATIN 总时间为 9.433297 s。
11. OPT-2 的 10 周期 LATIN 总时间为 174.717540 s。
12. OPT-2 相对 OPT-1 的 10 周期 LATIN 时间下降约 21.7%。
13. 从原始实现到 OPT-2，10 周期 LATIN 总时间累计下降约 62.6%。
14. 当前 10 周期 FOM/LATIN 时间比为 1.684479。

---

## 17. 当前仍然没有被证明的内容

以下内容不能因为 OPT-2 成功就直接宣称：

### 17.1 不能宣称当前已经达到最优实现

当前仍然是 Python 标量循环实现。

局部本构阶段还可能存在：

- Python 函数调用
- 标量 `np.sign`
- 多层内部函数调用
- 每步数组构造
- 每个材料点独立 Python 循环
- 四阶龙格-库塔阶段中重复的数据转换
- 其他固定量或中间量重复计算

因此当前实现仍有优化空间。

---

### 17.2 不能继续引用旧性能分析中的 95.6% 作为当前占比

此前“局部本构阶段约占 LATIN 求解器函数级统计时间 95.6%”是对 **OPT-1 之前的旧实现** 做出的函数级性能分析结果。

经过 OPT-1 和 OPT-2 后，热点分布已经发生变化。

因此不能把：

$$
95.6\%
$$

作为当前 OPT-2 实现的局部阶段耗时占比继续引用。

若后续需要重新建立当前热点分布，必须重新做性能分析或低干扰分模块计时。

---

### 17.3 不能宣称 LATIN-PGD 已经在长周期问题上普遍优于 FOM

当前确认的是：

- 当前塔筒模型
- 当前 320 个材料点
- 当前荷载形式
- 当前时间步长
- 当前误差容限
- 当前 10 周期范围
- 当前残差最小二乘空间策略
- 当前硬件与 Python 实现

下的实际结果。

尚未验证：

- 20 周期
- 50 周期
- 100 周期
- 更密纤维离散
- 更细塔筒有限元网格
- 更复杂外载
- 更复杂损伤演化
- 不同收敛容限

因此不能把当前 10 周期性能关系外推为普遍结论。

---

### 17.4 不能证明剩余主要瓶颈仍然一定是原来的 Local 95.6%

虽然理论上完整材料历史仍然需要被多次积分，局部阶段很可能继续是重要成本来源，但两轮优化已经大幅降低其单次计算成本。

因此当前剩余瓶颈究竟是：

- 局部本构积分
- PGD 富集
- 搜索方向
- 时间函数更新
- 空间最小二乘
- Python 控制流

中的哪一部分，需要新的测量才能确认。

---

## 18. 对下一阶段的建议

当前不建议马上进行算法结构的大修改。

更稳妥的路线是继续坚持：

> 先完成低风险、理论等价、可单独归因的实现级优化，再讨论是否需要改变 LATIN-PGD 算法结构。

下一阶段可以考虑两类路线。

### 路线 A：重新测量当前热点分布

经过 OPT-1 和 OPT-2 后，原来的函数级耗时比例已经失效。

可重新进行一次有限范围性能诊断，确认：

- Local 目前还占多少
- 四阶龙格-库塔目前还占多少
- PGD 富集占比是否开始上升
- 空间最小二乘求解是否开始成为次级瓶颈

这有助于避免继续凭旧热点数据盲目优化。

### 路线 B：继续优化局部本构 Python 热循环

如果重新确认 Local 仍为主瓶颈，可以继续检查：

- `_local_state_rate()` 层级函数调用
- `_interpolate()` 高频标量函数
- 每次 RK4 阶段的 NumPy 小数组构造
- `np.asarray(rates, dtype=np.float64)`
- `np.sign` 的标量调用
- 材料点空间维度上的向量化
- 独立材料点之间的并行化或编译加速

其中，材料点之间彼此独立，而单个材料点沿时间方向存在历史依赖，因此未来更适合：

> 保持时间方向顺序积分，同时在材料点方向进行向量化、编译或并行。

但这已经属于比 OPT-1、OPT-2 更大的实现结构调整，应当单独建立验证阶段。

---

## 19. Git 状态

本阶段代码已经提交并推送至：

`perf/tower-local-stage-optimization`

代码 checkpoint：

`e2d15b9`

提交信息：

`perf: cache material derived constants`

远端已经同步到该 checkpoint。

当前工作区中的以下诊断脚本和性能分析文件仍保持未跟踪状态，没有进入正式提交：

```text
tower_10cycle_profile.pstats
tower_asymmetric_10cycle_latin_diagnostic.py
tower_asymmetric_1cycle_efficiency_pilot.py
tower_asymmetric_1cycle_efficiency_pilot_v2.py
tower_asymmetric_1cycle_efficiency_pilot_v3.py
tower_asymmetric_1cycle_latin_diagnostic.py
tower_asymmetric_1cycle_latin_tight_diagnostic.py
tower_asymmetric_2cycle_latin_diagnostic.py
tower_asymmetric_2cycle_latin_diagnostic_fixed.py
tower_asymmetric_5cycle_latin_diagnostic.py
tower_asymmetric_efficiency_scaling_pilot.py
tower_fresh_local_fom_global_diagnostic.py
tower_local_fom_global_diagnostic.py
```

这些文件不应通过 `git add .` 意外进入阶段总结提交。

---

## 20. 阶段结论

本阶段完成了第二轮受控实现级性能优化。

通过将：

$$
k=K^{-n}
$$

$$
Y_0=\frac{\sigma_y^2}{2E}
$$

$$
\sqrt{\gamma}
$$

改为材料对象级缓存，并将机器浮点精度改为模块级一次性求值，在不改变任何 LATIN-PGD 数学流程、材料本构方程、四阶龙格-库塔积分和收敛路径的条件下：

- 313 个仓库测试全部通过；
- 1 周期路径保持完全一致；
- 10 周期路径保持完全一致；
- 1 周期 LATIN 总时间由 OPT-1 的 13.69 s 降至 9.43 s；
- 10 周期 LATIN 总时间由 OPT-1 的 223.16 s 降至 174.72 s；
- OPT-2 相对 OPT-1 的 10 周期时间进一步下降约 21.7%；
- 相对于最初实现，10 周期 LATIN 总时间累计下降约 62.6%；
- 当前 10 周期 FOM/LATIN 时间比达到 1.684。

因此：

**OPT-2：局部本构固定材料派生量缓存优化与多周期等价性验证 — PASS**

本阶段进一步强化了当前项目的一项核心认识：

> 在长时程循环损伤 LATIN-PGD 实现中，算法效率不能仅从理论公式和模式数量判断。数千万次局部本构评价会把原本看似微小的 Python 标量实现开销放大为主导性成本。只有先清除这些实现级开销，才能更公平地评价 LATIN-PGD 方法本身相对于完整有限元计算的真实效率。
