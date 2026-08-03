# NREL 5 MW 风机塔筒几何尺寸选型、简化措施及文献依据

- 日期：2026-08-03
- 项目：Offshore Wind Turbine and LATIN-PGD
- 目标分支：`feature/offshore-wind-turbine-tower-fatigue`
- 文档目的：记录海上风机塔筒第一阶段模型的几何尺寸选型、简化措施、采用理由及相关文献依据，便于后续建模、代码实现和论文写作时追溯。

---

## 1. 第一阶段模型定位

本项目拟将已在一维三材料杆算例中实现的 LATIN-PGD 方法，逐步扩展至海上风机钢制塔筒的风致疲劳分析。

从由易到难的开发原则出发，第一阶段暂不建立完整风机系统，而是采用：

> 基于 NREL 5 MW 参考塔筒几何的 tower-only 简化模型。

该模型只保留钢制锥形空心圆筒塔筒，不包含叶片、轮毂、机舱和基础等其他组成部分。

第一阶段模型的主要目的不是复现完整 NREL 5 MW 风机的动力特性，而是依次验证：

1. 变截面梁柱单元；
2. 圆环截面几何计算；
3. 纤维截面离散；
4. 单轴循环黏塑性损伤材料；
5. LATIN 局部阶段与全局平衡；
6. PGD 空间—时间分离表示。

## 2. 塔筒几何原型的选取

建议采用 NREL 5 MW 参考风机中被大量文献采用的锥形钢制塔筒几何。

该塔筒可理想化为：

- 竖直悬臂结构；
- 锥形空心圆筒；
- 外径沿高度线性减小；
- 壁厚沿高度线性减小；
- 塔底完全固结；
- 塔顶在第一阶段不设置集中质量。

### 2.1 建议采用的几何参数

| 参数 | 符号 | 数值 |
|---|---:|---:|
| 塔筒高度 | $H$ | $87.6\ \mathrm{m}$ |
| 塔底外径 | $D_b$ | $6.00\ \mathrm{m}$ |
| 塔顶外径 | $D_t$ | $3.87\ \mathrm{m}$ |
| 塔底壁厚 | $t_b$ | $0.027\ \mathrm{m}$ |
| 塔顶壁厚 | $t_t$ | $0.019\ \mathrm{m}$ |
| 截面类型 | — | 锥形空心圆环截面 |
| 塔底边界 | — | 完全固结 |
| 塔顶附加质量 | — | 第一阶段不考虑 |
| 外径变化 | — | 沿高度线性变化 |
| 壁厚变化 | — | 沿高度线性变化 |

本项目第一阶段将该组参数记为：

```text
Geometry G1: nominal NREL 5 MW tower geometry
H  = 87.6 m
Db = 6.00 m
Dt = 3.87 m
tb = 0.027 m
tt = 0.019 m
```

## 3. 外径和壁厚沿高度的变化

以塔底为坐标原点：

$$
z=0,
$$

塔顶为：

$$
z=H.
$$

### 3.1 外径变化

$$
D(z)=D_b+\frac{D_t-D_b}{H}z.
$$

代入参数：

$$
D(z)=6.00-0.024315z.
$$

### 3.2 壁厚变化

$$
t(z)=t_b+\frac{t_t-t_b}{H}z.
$$

代入参数：

$$
t(z)=0.027-9.1324\times10^{-5}z.
$$

### 3.3 内径

$$
D_i(z)=D(z)-2t(z).
$$

## 4. 截面几何性质

任意高度处的截面面积为：

$$
A(z)=\frac{\pi}{4}\left[D(z)^2-D_i(z)^2\right].
$$

关于任一水平主轴的截面二次矩为：

$$
I(z)=\frac{\pi}{64}\left[D(z)^4-D_i(z)^4\right].
$$

第一阶段纯弹性 Euler–Bernoulli 梁模型中可以直接采用：

$$
EA(z),\qquad EI(z).
$$

在后续纤维截面模型中，通过纤维应力积分获得截面轴力和弯矩：

$$
N=\sum_{f=1}^{n_f}\sigma_f A_f,
$$

$$
M=\sum_{f=1}^{n_f}\sigma_f y_f A_f.
$$

## 5. 第一阶段的简化措施

### 5.1 暂不考虑的结构组成

- 叶片；
- 轮毂；
- 机舱；
- 传动链和发电机；
- 塔顶集中质量；
- 单桩或其他海上支撑结构；
- 地基与土体。

### 5.2 暂不考虑的荷载和物理效应

- 叶轮气动推力的精细分布；
- 叶片旋转产生的周期荷载；
- 塔筒自重；
- 叶轮—机舱系统自重；
- 重力引起的恒定轴向压力；
- 波浪和海流荷载；
- 土—结构相互作用；
- 基础柔度；
- 结构质量和惯性力；
- 阻尼；
- 几何非线性；
- 局部壳屈曲；
- 焊缝和法兰等局部构造细节。

### 5.3 第一阶段保留的内容

- 塔筒真实量级的高度；
- 塔底和塔顶外径；
- 塔底和塔顶壁厚；
- 外径与壁厚沿高度的渐变；
- 钢材单轴本构；
- 塔底固结边界；
- 横向静力或准静力循环荷载；
- 塔筒弯曲引起的截面应变与应力分布。

## 6. 采用上述简化的理由

### 6.1 避免同时引入多类误差源

若第一版模型同时包含风机顶部质量、自重、动力惯性、阻尼、气动荷载、波浪荷载和基础柔度，则计算结果异常时难以判断问题来自：

- 梁单元推导；
- 变截面插值；
- 截面积分；
- 材料积分；
- LATIN 搜索方向；
- PGD 模态求解；
- 动力平衡；
- 荷载输入；
- 边界条件。

先建立 tower-only 准静力模型，可以显著提高错误定位效率。

### 6.2 保持与现有一维材料模型的衔接

钢制塔筒采用纤维梁柱方法后，每个截面纤维仍然主要承受沿塔筒轴线方向的单轴应变与应力。

因此，现有的一维循环黏塑性损伤模型可先作为纤维材料模型继续使用，而不需要立即发展三维钢材本构。

### 6.3 优先验证结构算法主链路

$$
\text{变截面弹性梁}
\rightarrow
\text{圆环纤维截面}
\rightarrow
\text{材料非线性}
\rightarrow
\text{准静力循环分析}
\rightarrow
\text{LATIN-PGD}
\rightarrow
\text{动力与随机风荷载}.
$$

## 7. NREL 壁厚参数的解释

NREL 5 MW 技术报告给出的基础几何为：

$$
D_b=6.00\ \mathrm{m},\quad
t_b=0.027\ \mathrm{m},
$$

$$
D_t=3.87\ \mathrm{m},\quad
t_t=0.019\ \mathrm{m}.
$$

报告同时说明，为了提高塔筒刚度并调整固有频率，最终分布属性的壁厚相对于上述基础数值进行了约 $30\%$ 的放大。

按比例换算：

$$
t_{b,\mathrm{eff}}=1.3\times0.027=0.0351\ \mathrm{m},
$$

$$
t_{t,\mathrm{eff}}=1.3\times0.019=0.0247\ \mathrm{m}.
$$

因此文献中可能出现两种处理：

| 处理方式 | 塔底壁厚 | 塔顶壁厚 | 主要用途 |
|---|---:|---:|---|
| 名义几何壁厚 | $27.0\ \mathrm{mm}$ | $19.0\ \mathrm{mm}$ | 直接建立几何模型，文献使用广泛 |
| 刚度调整等效壁厚 | 约 $35.1\ \mathrm{mm}$ | 约 $24.7\ \mathrm{mm}$ | 更接近 NREL 最终分布刚度 |

本项目第一阶段采用名义几何壁厚：

$$
t_b=27\ \mathrm{mm},\qquad t_t=19\ \mathrm{mm}.
$$

后续进行动力分析时，应重新讨论是否采用刚度调整后的等效壁厚，或直接采用 NREL 报告给出的分布质量与刚度数据。

## 8. 相关参考文献及其支撑内容

### [1] NREL 原始技术报告

Jonkman, J., Butterfield, S., Musial, W., and Scott, G.
**Definition of a 5-MW Reference Wind Turbine for Offshore System Development.**
NREL Technical Report, NREL/TP-500-38060, 2009.
DOI: `10.2172/947422`

支撑内容：

- 建立 NREL 5 MW 参考风机；
- 塔筒基础外径为 $6.00\rightarrow3.87\ \mathrm{m}$；
- 基础壁厚为 $0.027\rightarrow0.019\ \mathrm{m}$；
- 半径和壁厚沿高度线性变化；
- 说明最终分布属性中采用了约 $30\%$ 的壁厚放大。

### [2] Zuo et al. (2018)

Zuo, H., Bi, K., and Hao, H.
**Dynamic analyses of operating offshore wind turbines including soil-structure interaction.**
*Engineering Structures*, 157 (2018), 42–62.
DOI: `10.1016/j.engstruct.2017.12.001`

支撑内容：

- 采用 NREL 5 MW 海上风机；
- 塔顶和塔底外径分别为 $3.87\ \mathrm{m}$ 和 $6.00\ \mathrm{m}$；
- 对应壁厚分别为 $0.019\ \mathrm{m}$ 和 $0.027\ \mathrm{m}$；
- 外径和壁厚沿塔筒高度线性减小；
- 采用 ABAQUS 进行风—浪—结构动力分析。

### [3] Kazemi Esfeh and Kaynia (2020)

Kazemi Esfeh, P., and Kaynia, A. M.
**Earthquake response of monopiles and caissons for Offshore Wind Turbines founded in liquefiable soil.**
*Soil Dynamics and Earthquake Engineering*, 136 (2020), 106213.
DOI: `10.1016/j.soildyn.2020.106213`

支撑内容：

- 参数表明确给出塔筒高度 $87.6\ \mathrm{m}$；
- 塔顶外径和壁厚为 $3.87\ \mathrm{m}$、$0.019\ \mathrm{m}$；
- 塔底外径和壁厚为 $6.00\ \mathrm{m}$、$0.027\ \mathrm{m}$。

### [4] Wang et al. (2023)

Wang, P.-G., Lu, H.-Q., Wang, M., Nagarajaiah, S., and Du, X.-L.
**Experimental and numerical investigations on seismic responses of wind turbine structures with amplifying damping transfer system.**
*Soil Dynamics and Earthquake Engineering*, 175 (2023), 108277.
DOI: `10.1016/j.soildyn.2023.108277`

支撑内容：

- 明确描述 NREL 5 MW 锥形钢塔筒；
- 塔筒高度为 $87.6\ \mathrm{m}$；
- 外径和壁厚从塔底 $6.00\ \mathrm{m}$、$0.027\ \mathrm{m}$ 渐变至塔顶 $3.87\ \mathrm{m}$、$0.019\ \mathrm{m}$；
- 以该模型作为缩尺试验与数值模拟原型。

### [5] Sorge et al. (2024)

Sorge, E., Riascos, C., Caterino, N., Demartino, C., and Georgakis, C. T.
**Optimal design of a hinge-spring-friction device for enhancing wind induced structural response of onshore wind turbines.**
*Engineering Structures*, 314 (2024), 118305.
DOI: `10.1016/j.engstruct.2024.118305`

支撑内容：

- 采用 NREL 5 MW 风机作为算例；
- 塔筒高度为 $87.60\ \mathrm{m}$；
- 塔底外径和壁厚为 $6.00\ \mathrm{m}$、$0.027\ \mathrm{m}$；
- 塔顶外径和壁厚为 $3.87\ \mathrm{m}$、$0.019\ \mathrm{m}$；
- 研究风致响应、塔底弯矩及疲劳需求。

### [6] Duan et al. (2024)

Duan, L.-X., Wang, W.-D., Zheng, L., and Shi, Y.-L.
**Dynamic response analysis of monopile CFDST wind turbine tower system under wind-wave-seismic coupling action.**
*Thin-Walled Structures*, 202 (2024), 112089.
DOI: `10.1016/j.tws.2024.112089`

支撑内容：

- 采用 NREL 5 MW OC3 单桩风机塔筒作为原型；
- 塔筒总高度为 $87.6\ \mathrm{m}$；
- 塔底外径和壁厚为 $6.00\ \mathrm{m}$、$0.027\ \mathrm{m}$；
- 塔顶外径和壁厚为 $3.87\ \mathrm{m}$、$0.019\ \mathrm{m}$；
- 对传统钢制塔筒和混凝土充填双钢管塔筒进行对比。

## 9. 文献依据的综合判断

上述文献表明，以下几何参数已在多篇 Elsevier 期刊论文中重复采用：

$$
H=87.6\ \mathrm{m},
$$

$$
D_b=6.00\ \mathrm{m},\qquad D_t=3.87\ \mathrm{m},
$$

$$
t_b=0.027\ \mathrm{m},\qquad t_t=0.019\ \mathrm{m}.
$$

将其用于本项目第一阶段塔筒模型具有以下合理性：

- 几何来源可追溯；
- 已被多个海上和陆上风机研究采用；
- 适用于有限元、动力响应、地震响应和风致疲劳研究；
- 可作为统一基准，便于后续与文献结果对比；
- 结构尺寸具有真实工程量级；
- 几何变化规律简单，适合算法开发和逐步验证。

## 10. 当前阶段的最终建议

```text
Model name:
NREL-5MW tower-only simplified model

Geometry:
H  = 87.6 m
Db = 6.00 m
Dt = 3.87 m
tb = 0.027 m
tt = 0.019 m

Section:
Linearly tapered hollow circular steel section

Boundary:
Fully fixed tower base

Temporarily neglected:
Rotor
Blades
Hub
Nacelle
Top lumped mass
Tower self-weight
Gravity axial load
Foundation flexibility
Soil-structure interaction
Wave and current loads
Structural inertia
Damping
Geometric nonlinearity
Local shell buckling

Initial loading:
Static or quasi-static transverse cyclic load
```

该模型应明确表述为：

> 基于 NREL 5 MW 名义塔筒几何建立的 tower-only 结构算法验证模型。

不应将其表述为完整的 NREL 5 MW 海上风机模型。
