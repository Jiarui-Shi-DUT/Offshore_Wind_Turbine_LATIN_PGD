# PyCharm Markdown 数学公式预览异常：排查过程、根因与后续规避规范

**项目：** Offshore Wind Turbine and LATIN-PGD
**日期：** 2026-08-10
**问题类型：** PyCharm Markdown Preview / LaTeX 数学公式渲染异常
**建议仓库路径：** `docs/2026-08-10-pycharm-markdown-math-preview-debugging-lessons.md`

---

## 1. 问题背景

在查看以下阶段总结文档时：

```text
docs/2026-08-09-tower-100cycle-fom-tensorization-low-rank-stage-summary.md
```

PyCharm Markdown Preview 出现异常：

- 文档前部的数学公式突然不再正常渲染；
- `\text{...}`、`\longrightarrow`、`\boxed{...}`、`\varepsilon_p` 等 LaTeX 命令直接以源码形式显示；
- 异常并不是只发生在某一条公式附近，而是会导致文档中大量数学公式整体失去正常渲染。

最初容易误判为：

1. Markdown 文档过长；
2. `$$...$$` 数学块没有闭合；
3. 某条复杂 LaTeX 公式存在语法错误；
4. `\mathrm`、`\boxed`、`\text` 等命令不受 PyCharm 支持；
5. PyCharm 在约 1200 行后存在固定的 Preview 长度限制。

后续通过逐步缩小范围和最小化测试，最终证明上述判断均不是根因。

---

## 2. 最终确认的根因

最终定位到原文第 20 节中的公式：

```latex
$$r^{\mathrm{I}}<r^{\mathrm{II}}<r^{\mathrm{III}}$$
```

其中数学环境内直接使用了字面量：

```text
<
```

在当前 PyCharm Markdown Preview 环境中，这一字符组合会触发异常解析。

其表现具有迷惑性：

- 单独公式并不一定只在当前位置显示错误；
- 它可能导致整个 Markdown Preview 的数学渲染过程受到影响；
- 因此文档前面原本正确的公式也可能突然显示为 LaTeX 源码。

将其修改为 LaTeX 命令：

```latex
$$r^{\mathrm{I}}\lt r^{\mathrm{II}}\lt r^{\mathrm{III}}$$
```

后，完整约 1287 行 Markdown 文档恢复正常渲染。

因此本次问题的直接结论为：

> **在 PyCharm Markdown 数学环境中，应避免直接使用字面量 `<`，优先使用 LaTeX 命令 `\lt`。**

---

## 3. 本次排查过程得到的关键证据

### 3.1 初始二分定位

首先通过截取原文前若干行生成测试文件：

```text
前 966 行
前 1127 行
前 1207 行
前 1247 行
前 1227 行
前 1217 行
前 1212 行
前 1211 行
前 1208 行
```

观察到：

```text
前 1207 行：正常
前 1208 行：异常
```

这最初似乎指向第 1208 行。

但第 1208 行仅为普通中文：

```text
说明全局结构响应高度可压缩。
```

随后检查其 Unicode 字符，也没有发现：

- BOM；
- 零宽字符；
- 不间断空格；
- 控制字符；
- 其他异常编码。

因此第 1208 行并不是实际错误源。

---

### 3.2 用 `TEST` 替换第 1208 行

将第 1208 行替换为：

```text
TEST
```

结果仍然异常。

说明：

> **异常不是由第 1208 行文字内容造成的。**

---

### 3.3 检查数学块和代码块闭合情况

统计前 1208 行：

```text
math_delimiter_count = 282
parity = 0

code_fence_count = 22
parity = 0
```

说明：

- `$$...$$` 数学块整体成对；
- ``` 代码围栏整体成对。

进一步统计行内 `$...$`：

```text
inline_dollar_count = 158
parity = 0
```

因此也排除了：

- 块公式分隔符未闭合；
- 行内公式分隔符未闭合；
- Markdown 代码块未闭合。

---

### 3.4 排除第 1204、1206 行公式

当时重点怀疑：

```latex
$$\Delta u:(2,2,1)$$
```

和：

```latex
$$\Delta\sigma:(1,1,2)$$
```

先改成多行块公式，再直接替换为普通文字：

```text
FORMULA_A
FORMULA_B
FORMULA_C
```

Preview 仍然异常。

因此这些公式也不是根因。

---

### 3.5 排除“Markdown 总行数限制”

建立完全独立的 1300 行 Markdown control 文件：

- 只有一个简单数学公式；
- 后续约 1300 行均为普通文本。

结果：

```text
1300 control：正常
```

因此不存在简单的：

```text
超过约 1200 行 → PyCharm 数学 Preview 必然失效
```

这一固定阈值。

---

### 3.6 拆分文档继续定位

将原 1287 行总结拆为：

```text
Part 1：第 1–12 节
Part 2：第 13–24 节
```

结果：

```text
Part 1：正常
Part 2：异常
```

继续拆分：

```text
2A：第 13–18 节
2B：第 19–24 节
```

异常进一步定位到 2B。

随后：

```text
2B1：第 19–21 节
2B2：第 22–24 节
```

异常进一步定位到 2B1。

继续：

```text
2B1A：第 19–20 节
2B1B：第 21 节
```

异常位于 2B1A。

再拆为：

```text
Section 19
Section 20
```

最终确认：

```text
Section 20：异常
```

---

### 3.7 最小公式测试最终确认

第 20 节继续拆分后，最终只剩两条主要可疑公式。

#### 公式 1

```latex
$$r^{\mathrm{I}}<r^{\mathrm{II}}<r^{\mathrm{III}}$$
```

结果：

```text
异常
```

#### 公式 2

```latex
$$\boxed{\text{fatigue mechanism transition}\rightarrow\text{rank evolution}\rightarrow\text{adaptive PGD enrichment}}$$
```

结果：

```text
正常
```

因此问题明确锁定到公式 1。

进一步把：

```latex
r^{\mathrm{I}}
```

改成：

```latex
r^{I}
```

但保留 `<`：

```latex
$$r^{I}<r^{II}<r^{III}$$
```

结果仍然异常。

因此排除：

```latex
\mathrm
```

作为根因。

最后将 `<` 替换为：

```latex
\lt
```

即：

```latex
$$r^{\mathrm{I}}\lt r^{\mathrm{II}}\lt r^{\mathrm{III}}$$
```

结果正常。

至此完成根因确认。

---

## 4. 为什么这次问题很容易误判

本次最重要的经验是：

> **Markdown Preview 中“错误开始出现的位置”不一定等于“真正的错误源位置”。**

本次曾出现：

```text
前 1207 行正常
前 1208 行异常
```

但真正错误实际上位于更早的第 20 节公式结构中，而不是普通中文第 1208 行。

更准确地说：

- Preview 是对整个 Markdown 文档进行解析和渲染；
- 某些字符可能被 Markdown/HTML 层先解释；
- 一旦解析状态异常，后续甚至前面的数学渲染显示都可能受到影响；
- 所以不能仅根据“第一处肉眼看到异常的位置”直接判断错误行。

---

## 5. 后续 Markdown 数学公式书写规范

为了避免再次触发类似问题，本项目后续 `.md` 文档统一采用以下规则。

### 5.1 数学环境中避免直接使用 `<`

不推荐：

```latex
$$a<b<c$$
```

推荐：

```latex
$$a\lt b\lt c$$
```

---

### 5.2 同理，比较符号优先采用明确 LaTeX 命令

推荐：

```latex
\lt
\gt
\leq
\geq
\neq
\approx
```

而不是在复杂数学环境中依赖 Markdown/HTML 可能同时具有语法意义的字符组合。

例如：

```latex
$$r_1\lt r_2$$
```

```latex
$$D\geq D_{\mathrm{crit}}$$
```

---

### 5.3 `$$...$$` 必须严格成对

长公式推荐写成：

```latex
$$
q(n,\tau,x)
\approx
\sum_{k=1}^{r}N_k(n)T_k(\tau)X_k(x)
$$
```

短公式可以继续写成：

```latex
$$D_{\max}=0.0441$$
```

但无论采用哪种形式，都必须确保 `$$` 成对。

---

### 5.4 行内数学 `$...$` 同样必须成对

推荐：

```markdown
cycle rank 为 $r_n=6$。
```

避免遗漏结束 `$`。

---

### 5.5 Markdown 与 LaTeX 语法冲突时优先使用 LaTeX 命令

Markdown Preview 不是纯 LaTeX 编译器。

其处理过程通常同时涉及：

```text
Markdown
↓
HTML / Preview DOM
↓
数学公式渲染
```

因此具有 Markdown/HTML 特殊意义的字符，应尽可能通过 LaTeX 命令表达。

---

## 6. 后续遇到类似异常时的标准排查流程

以后若再次发生“整篇公式突然显示为源码”，不要直接修改大量公式。

建议严格按以下顺序排查。

### Step 1：先建立最小正常基线

生成文档前半部分，确认：

```text
某个截断位置之前 Preview 正常
```

### Step 2：二分定位异常区间

例如：

```text
1–1200 正常
1–1400 异常
```

则继续测试：

```text
1–1300
```

不断缩小。

### Step 3：不要把截断边界自动当成根因

若：

```text
前 N 行正常
前 N+1 行异常
```

必须进一步验证：

- 第 N+1 行换成 `TEST` 是否仍异常；
- 将 `TEST` 插到其他位置是否仍异常。

### Step 4：检查结构闭合

优先检查：

```text
$$
$
```

以及：

```text
```
```

是否成对。

### Step 5：拆分成独立章节

如果局部测试仍不能确定原因，则按：

```text
Part A / Part B
Section A / Section B
```

继续二分。

### Step 6：最终必须建立“单公式最小测试”

不能仅根据上下文猜测。

例如分别建立：

```text
formula1-test.md
formula2-test.md
```

每个文件只放一条公式。

只有这样才能确认：

> **到底是哪一条表达式能够独立复现错误。**

### Step 7：一次只改变一个语法因素

例如：

原公式：

```latex
$$r^{\mathrm{I}}<r^{\mathrm{II}}<r^{\mathrm{III}}$$
```

第一次只改：

```latex
\mathrm
```

第二次只改：

```text
<
```

通过单变量实验，才能可靠确认真正根因。

---

## 7. PowerShell 调试命令的额外教训

本次排查中还出现了一个与 Markdown 无关、但非常重要的 PowerShell 问题。

在 PowerShell 双引号字符串中直接写：

```text
$$
```

可能发生 PowerShell 自身的变量展开，从而破坏传递给：

```text
python -c "..."
```

的 Python 代码。

曾导致 Python 收到被意外替换的内容，并出现：

```text
SyntaxError
```

因此以后从 PowerShell 构造 Markdown 数学分隔符时，不建议直接在复杂双引号命令里写：

```python
"$$"
```

更稳妥的方式是：

```python
d = chr(36) * 2
```

其中：

```text
chr(36) = $
```

因此：

```python
d = chr(36) * 2
```

等价于：

```text
$$
```

同理，反斜杠也可以在特殊情况下使用：

```python
chr(92)
```

避免 Shell 与 Python 多层转义互相干扰。

---

## 8. 对以后调试策略的改进

本次早期排查存在一个值得纠正的问题：

> 在证据尚不足时，过早把异常归因于“文档过长”“某个数学块”“PyCharm Preview 长度阈值”等假设。

以后应优先采用：

```text
观察
→ 对照实验
→ 最小复现
→ 单变量修改
→ 确认根因
```

而不是：

```text
观察
→ 猜测原因
→ 大范围修改
```

特别是在科研代码和科研文档中，应尽量避免为了修复 Preview 而大面积改写本来正确的公式和正文。

---

## 9. 本次最终修复

原公式：

```latex
$$r^{\mathrm{I}}<r^{\mathrm{II}}<r^{\mathrm{III}}$$
```

修改为：

```latex
$$r^{\mathrm{I}}\lt r^{\mathrm{II}}\lt r^{\mathrm{III}}$$
```

除此之外，不需要修改文档中的：

```latex
\mathrm
\text
\boxed
\longrightarrow
\varepsilon
```

等正常 LaTeX 表达式。

修复后：

> **完整 1287 行阶段总结在 PyCharm Markdown Preview 中恢复正常。**

---

## 10. 项目级经验总结

本次问题最终形成以下项目级规范：

1. **PyCharm Preview 全局失效时，不要默认错误位于当前显示异常的位置。**
2. **先用二分和章节拆分确定能够独立复现问题的最小区域。**
3. **最终必须通过单公式最小复现验证根因。**
4. **数学环境中避免直接使用字面量 `<`，统一优先写成 `\lt`。**
5. **类似比较运算符优先采用明确 LaTeX 命令，降低 Markdown/HTML 与数学解析冲突。**
6. **不要因为 Preview 出错而大范围改动科研正文；先证明是哪一处语法真正触发问题。**
7. **PowerShell 中构造 `$$` 时注意 Shell 变量展开，复杂命令优先使用 `chr(36)*2`。**
8. **调试必须坚持“一个实验只改变一个变量”。**

---

## 11. 一句话结论

> **本次 PyCharm Markdown 数学预览异常最终不是文档长度、公式数量或 `\mathrm` 所致，而是数学公式中直接使用字面量 `<` 与 Preview 解析产生冲突；将其改为 LaTeX 命令 `\lt` 后，完整文档恢复正常。后续应通过最小复现和单变量实验定位类似问题，并统一采用更稳健的 LaTeX 写法。**
