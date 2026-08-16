# PyCharm Markdown + LaTeX 公式编写经验教训：避免公式断裂、列表误判与不兼容命令

**日期：2026-08-16**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**适用环境：PyCharm Markdown 编辑器 + 右侧 Preview**  
**背景：在编写 `2026-08-16-tower-latin-pgd-eq61-64-enrichment-stage-summary.md` 时，连续出现了公式块被 Markdown 误解析、LaTeX 命令不受支持等问题。本文档总结本次排错得到的稳定写法，作为后续所有项目 Markdown 文档的格式规范。**

---

# 1. 本文档的目的

本次阶段总结文档在内容推导上没有问题，但在 PyCharm Markdown Preview 中先后出现了以下格式错误：

1. `\[ ... \]` 行间公式无法稳定识别；
2. 多行 `$$ ... $$` 中的 `-`、`+` 被 Markdown 当成列表项目符号；
3. 公式虽然位于 `$$ ... $$` 中，但被拆成若干文本行；
4. `\boldsymbol` 在当前 PyCharm 数学渲染器中不受支持，被标红；
5. 同一文档中同时混用多种 LaTeX 分隔符，导致预览行为不一致；
6. 为了修复显示问题而引入新的格式说明，增加了文档噪声。

因此需要形成统一规范：

> **后续项目内所有 Markdown 数学公式，优先采用已经在当前 PyCharm 环境中验证可稳定渲染的写法，而不是只考虑标准 LaTeX 语法是否合法。**

---

# 2. 第一条核心经验：不要在项目 Markdown 中使用 `\[ ... \]`

本次最初采用：

```markdown
\[
\boxed{\text{PGD: Add a pair}}
\]
```

在当前 PyCharm Preview 中不能稳定识别。

因此后续统一改为：

```markdown
$$ \boxed{\text{PGD: Add a pair}} $$
```

结论：

\[
\boxed{\text{项目 Markdown 中不再使用 } \backslash[ \cdots \backslash]}
\]

统一采用：

```markdown
$$ ... $$
```

作为行间公式分隔符。

---

# 3. 第二条核心经验：行间公式必须尽量保持在一个“物理行”内

这是本次最重要的经验。

错误写法：

```markdown
$$
\boxed{
\Delta\dot{\varepsilon}^{p}_{i+1}
-
H_\sigma\Delta\sigma'_{i+1}
+
\bar{\Delta}_{i+1}
=0
}
$$
```

虽然这在标准 LaTeX 中完全合理，但当前 PyCharm Markdown Preview 会把公式内部的独立行：

```text
-
```

和：

```text
+
```

误当成 Markdown 列表符号。

最终预览结果会被拆成：

- 第一段公式；
- 一个项目符号；
- 第二段公式；
- 又一个项目符号；
- 第三段公式。

因此，后续所有重要公式统一写成：

```markdown
$$ \boxed{ \Delta\dot{\varepsilon}^{p}_{i+1} - H_\sigma\Delta\sigma'_{i+1} + \bar{\Delta}_{i+1} = 0 } \tag{61} $$
```

即：

\[
\boxed{\text{一个 display-math block 尽量只占 Markdown 文件中的一个物理行}}
\]

这里的“一个物理行”是指：

> 从开头的 `$$` 到结尾的 `$$`，在源文件中不要换行。

---

# 4. 为什么多行公式会被误判为 Markdown 列表

Markdown 本身会把行首：

```text
-
```

解释成无序列表。

例如：

```markdown
-
H_\sigma\Delta\sigma'
```

Markdown parser 并不知道这本来只是公式中的减号。

同理，某些 Preview renderer 对：

```text
+
```

也可能做特殊解释。

所以，即使外层已经写了：

```markdown
$$
...
$$
```

如果当前 Preview 对 display math block 的解析不是完全优先于 Markdown block parser，就可能发生误判。

因此后续不要依赖：

> “既然在 `$$` 里面，Markdown 就一定不会再解析内部行。”

而应直接通过格式规避。

---

# 5. 第三条核心经验：行内公式统一使用 `$ ... $`

行内变量采用：

```markdown
$H_\sigma$
```

而不是：

```markdown
\(H_\sigma\)
```

当前 PyCharm 中：

```markdown
$H_\sigma$
```

已经验证稳定。

因此统一规则：

```text
行内公式：$ ... $
行间公式：$$ ... $$
```

不再混用：

```text
\( ... \)
\[ ... \]
```

这样可以最大限度减少不同 renderer 对分隔符支持不一致的问题。

---

# 6. 第四条核心经验：不要在当前 PyCharm Preview 中使用 `\boldsymbol`

本次 tower material-point vector 最初写成：

```latex
\boldsymbol{\varepsilon}
```

例如：

```markdown
$$ \boxed{ \bar{\boldsymbol{\varepsilon}}^p_{m+1} \in \mathbb R^{320} } $$
```

当前 PyCharm Preview 中：

```latex
\boldsymbol
```

被直接标红，说明当前数学渲染器不支持或没有加载对应 AMS/KaTeX 扩展。

因此不能继续使用。

---

# 7. 当前项目推荐的向量写法：`\vec{}`

为了兼顾：

1. 数学上明确表示“离散向量”；
2. PyCharm Preview 能正常渲染；
3. 不引入额外 LaTeX 宏包依赖；

本次最终采用：

```latex
\vec{\varepsilon}
```

例如：

```markdown
$$ \boxed{ \vec{\bar{\varepsilon}}^p_{m+1} \in \mathbb R^{320} } $$
```

对应：

```markdown
$$ \boxed{ \delta\dot{\vec{\varepsilon}}^p(t) = \dot{\lambda}_{m+1}(t)\vec{\bar{\varepsilon}}^p_{m+1} } $$
```

因此项目内推荐：

\[
\boxed{
\text{material-point / nodal 离散向量优先用 }\vec{\cdot}
}
\]

而不再使用：

```latex
\boldsymbol
```

---

# 8. `\vec{}` 只改变显示记号，不改变数学含义

这一点必须明确。

例如：

```latex
\bar{\varepsilon}^p(x)
```

表示连续 spatial field。

离散后：

```latex
\vec{\bar{\varepsilon}}^p
```

表示：

\[
[
\bar{\varepsilon}^p_{egf}
]
\]

组成的离散 material-point vector。

因此：

```text
bar / hat / dot
```

仍然保留原来的数学含义；

`\vec{}` 只用于说明：

> 该量已经被离散并以向量形式存储。

---

# 9. 推荐的复合记号顺序

本次发现复合记号如果层级太多，容易出现可读性问题。

推荐：

```latex
\vec{\bar{\varepsilon}}
```

而不是：

```latex
\bar{\vec{\varepsilon}}
```

推荐：

```latex
\vec{\hat{\sigma}}
```

而不是：

```latex
\hat{\vec{\sigma}}
```

对于时间导数推荐：

```latex
\dot{\vec{\varepsilon}}
```

对于塑性应变率：

```latex
\dot{\vec{\varepsilon}}^p
```

当前项目建议统一使用：

```text
vector → bar/hat inside vector accent
time derivative → outside vector notation
```

即：

```latex
\vec{\bar{\varepsilon}}
\vec{\hat{\sigma}}
\dot{\vec{\varepsilon}}
```

---

# 10. 第五条核心经验：公式中的长文本要谨慎使用 `\text{}`

以下写法在当前环境中已经验证可用：

```latex
\boxed{\text{PGD: Add a pair}}
```

以及：

```latex
\mathcal E_{\mathrm{tower}}
```

因此：

```latex
\text{}
```

可以保留。

但是不要在一个公式中塞入过长的自然语言。

例如不推荐：

```markdown
$$ \boxed{\text{the current reduced basis cannot adequately describe the remaining LATIN correction and therefore a new pair must be added}} $$
```

推荐拆成正文说明，再保留一个短公式：

```markdown
当前 basis 已不足，需要进行 enrichment：

$$ \boxed{\text{Add a new PGD pair}} $$
```

原因：

- 公式框会变得过宽；
- Preview 容易出现横向滚动；
- 移动端和 GitHub 页面显示也不友好。

---

# 11. 第六条核心经验：避免公式过宽

本次 Eq. (62) tower vector form 较长：

```latex
\vec{\bar{\Delta}}(t)
=
D_H(t)
(
\vec{\hat{\sigma}}(t)
-
\vec{\sigma}^{\mathrm{up}}(t)
)
-
(
\hat{\dot{\vec{\varepsilon}}}^p(t)
-
\dot{\vec{\varepsilon}}^{p,\mathrm{up}}(t)
)
```

虽然单行写法可以避免 Markdown 误判，但可能导致 Preview 横向过宽。

处理原则：

1. **优先保证公式不会被 Markdown parser 拆开；**
2. 如果过宽，优先拆成“定义 + 子式”，而不是在一个 `$$` block 内手工换行。

例如：

```markdown
定义：

$$ \vec{\bar{\Delta}}(t) = D_H(t)\Delta\vec{\sigma}(t) - \Delta\dot{\vec{\varepsilon}}^p(t) $$

其中：

$$ \Delta\vec{\sigma}(t) = \vec{\hat{\sigma}}(t) - \vec{\sigma}^{\mathrm{up}}(t) $$

$$ \Delta\dot{\vec{\varepsilon}}^p(t) = \hat{\dot{\vec{\varepsilon}}}^p(t) - \dot{\vec{\varepsilon}}^{p,\mathrm{up}}(t) $$
```

而不推荐：

```markdown
$$
\begin{aligned}
...
\end{aligned}
$$
```

在当前 PyCharm 环境下，后者虽然 LaTeX 更漂亮，但兼容风险更高。

---

# 12. 当前项目中 `\tag{}` 可以继续使用

本次以下形式已经正常显示：

```markdown
$$ \boxed{ ... } \tag{61} $$
```

因此 Eq. 编号可以保留：

```latex
\tag{61}
```

但应注意：

> `\tag{}` 只用于明确对应原论文公式编号，不应滥用为文档内部自动编号系统。

---

# 13. 当前验证可正常使用的 LaTeX 命令

在本次 PyCharm Preview 中，以下命令已经正常显示：

```latex
\boxed{}
\text{}
\vec{}
\bar{}
\hat{}
\dot{}
\mathbb{}
\mathcal{}
\operatorname{}
\mathrm{}
\left(
\right)
\begin{bmatrix} ... \end{bmatrix}
\frac{}{}
\sum
\int
\partial
\in
\rightarrow
\longrightarrow
\approx
\neq
\le
\ge
```

当前明确发现不兼容：

```latex
\boldsymbol{}
```

后续如果遇到其他命令被标红，应优先判断：

> 是 PyCharm renderer 不支持，而不是公式推导错误。

---

# 14. Markdown 标题和公式之间必须留空行

推荐：

```markdown
# 4. Eq. (61)

原论文 Eq. (61)：

$$ ... $$

形式上与 Eq. (41) 类似。
```

不推荐：

```markdown
原论文 Eq. (61)：
$$ ... $$
形式上与 Eq. (41) 类似。
```

虽然某些 renderer 能解析，但为了兼容性和可读性，正文、公式、列表、标题之间统一保留空行。

---

# 15. 列表附近尤其要注意公式换行

错误风险最高的组合：

```markdown
其中：

$$
a
-
b
$$

- 第一项；
- 第二项。
```

因为公式内部的 `-` 和正文列表的 `-` 对 Markdown parser 来说非常相似。

推荐：

```markdown
其中：

$$ a-b $$

- 第一项；
- 第二项。
```

即：

\[
\boxed{
\text{列表上下文中的 display math 更应坚持单物理行}
}
\]

---

# 16. 不要为了“公式美观”牺牲 Markdown 稳定性

标准 LaTeX 文档中经常推荐：

```latex
\begin{aligned}
a &= b+c \\
  &= d+e
\end{aligned}
```

但是当前项目文档的第一目标不是排版成期刊 PDF，而是：

1. 在 PyCharm 中稳定编辑；
2. 在 Preview 中即时检查；
3. 在 GitHub 上长期维护；
4. 在代码推导过程中方便搜索与 diff；
5. 不因 renderer 差异导致公式变形。

所以当前项目应遵循：

\[
\boxed{
\text{Markdown robustness > LaTeX visual elegance}
}
\]

如果以后需要正式论文排版，可以再将 Markdown 推导迁移到 LaTeX。

---

# 17. 后续阶段总结文档推荐公式模板

## 17.1 行内公式

```markdown
其中 $H_\sigma$ 为 descent search-direction operator。
```

---

## 17.2 简单行间公式

```markdown
$$ \zeta_i = \frac{\xi_i-\xi_{i+1}}{\xi_i+\xi_{i+1}} $$
```

---

## 17.3 带公式编号

```markdown
$$ \Delta\dot{\varepsilon}^{p}_{i+1} - H_\sigma\Delta\sigma'_{i+1} + \bar{\Delta}_{i+1} = 0 \tag{61} $$
```

---

## 17.4 强调公式

```markdown
$$ \boxed{ \Delta\dot{\varepsilon}^{p}_{i+1} - H_\sigma\Delta\sigma'_{i+1} + \bar{\Delta}_{i+1} = 0 } $$
```

---

## 17.5 离散向量

```markdown
$$ \vec{\bar{\varepsilon}}^p_{m+1} \in \mathbb R^{N_q} $$
```

---

## 17.6 简短算法关系

```markdown
$$ \text{fixed-basis update} \rightarrow \text{saturation check} \rightarrow \text{enrichment} $$
```

---

# 18. 新建 Markdown 文档后的强制检查流程

以后每次生成阶段总结 md，不要直接认为文件正确。

至少做以下检查。

## Step 1：检查公式分隔符

搜索：

```text
\[
\]
\(
\)
```

原则上项目阶段总结中应为：

```text
0 occurrences
```

统一改成：

```text
$ ... $
$$ ... $$
```

---

## Step 2：检查多行 display math

搜索是否存在：

```markdown
$$
```

单独占一行。

如果存在，优先改成：

```markdown
$$ formula $$
```

一个物理行。

---

## Step 3：检查不兼容命令

至少搜索：

```text
\boldsymbol
```

当前环境中应为：

```text
0 occurrences
```

如果需要向量记号，改用：

```text
\vec{}
```

---

## Step 4：检查公式中的孤立 `-` / `+`

特别检查：

```text
-
+
```

是否独立出现在数学块的源代码行首。

如果存在，必须合并公式为单行。

---

## Step 5：检查 PyCharm Preview

至少抽查：

1. 文档开头第一组公式；
2. 最长的一条公式；
3. 含 `\boxed{}` 的公式；
4. 含向量记号的公式；
5. 含上下标和 dot/bar/hat 组合的公式；
6. 文档最后一节公式。

只有 Preview 正常，才算完成。

---

# 19. 建议的自动化文本检查规则

以后如果用 Python 自动生成 md，可以在保存前做简单检查。

这里要特别注意一个本次新增发现的问题：

> **PyCharm 会对带有 `python` 语言标记的 fenced code block 进行 Python 静态检查。示例代码中如果直接使用未定义变量，即使它只是“说明性代码”，Problems 面板仍会报 `Unresolved reference`。**

因此，文档中的 Python 示例不要写成只有若干独立 `assert` 的片段，而应尽量写成一个**变量、导入和上下文都完整的最小可运行示例**。

推荐：

```python
from pathlib import Path
import re

md_path = Path("docs/example.md")
md_text = md_path.read_text(encoding="utf-8")

assert r"\[" not in md_text
assert r"\]" not in md_text
assert r"\(" not in md_text
assert r"\)" not in md_text
assert r"\boldsymbol" not in md_text

multiline_math_blocks = re.findall(
    r"\$\$\s*\n.*?\n\$\$",
    md_text,
    flags=re.S,
)

assert not multiline_math_blocks
```

这样：

- `Path` 已显式导入；
- `re` 已显式导入；
- `md_text` 在使用前已经定义；
- `multiline_math_blocks` 在断言前已经定义；
- PyCharm 不会再因为示例代码上下文不完整而报告 `Unresolved reference`。

不推荐的写法可以作为“反例”展示，但**不要再使用 `python` fenced code block**，否则 PyCharm 仍会对反例进行静态检查并报告错误。

应改成普通文本代码块：

```text
assert r"\[" not in text
assert r"\]" not in text
```

这里的目的只是展示“未定义变量的错误示例”，并不希望 IDE 把它当成可执行 Python 代码检查。

因此形成一条新的统一规则：

> **正确、可运行的 Python 示例使用 ` ```python `；故意错误、仅用于说明反例的代码使用 ` ```text `，避免 PyCharm 对反例产生 `Unresolved reference`。**

这类自动检查不能证明所有公式都正确，但可以防止本次已经遇到的格式错误再次出现。

---

# 20. Git diff 友好原则

项目使用 Git 管理，因此 Markdown 数学写法还应考虑 diff 可读性。

推荐：

```markdown
$$ a=b+c $$
```

而不是：

```markdown
$$
a
=
b
+
c
$$
```

原因：

- 单行公式修改时 diff 更集中；
- 不会因为排版换行产生大量无意义 diff；
- merge conflict 更少；
- 更容易比较前后公式变化。

因此，本次为了 PyCharm 兼容采用单行公式，同时也更符合当前仓库的 Git 维护需求。

---

# 21. 文档生成时不要加入临时调试说明

本次中间版本曾加入：

> “本版本已将上一版中未被 Markdown 预览正确识别的……”

这种说明对于调试有用，但不适合作为最终科研阶段总结正文。

因此以后流程应是：

```text
debug version
      ↓
preview verification
      ↓
remove temporary debug note
      ↓
final repository version
```

最终文档只保留真正有长期价值的技术内容。

---

# 22. 文档文件名保持研究阶段语义，而不是调试语义

不建议最终仓库使用：

```text
xxx-fixed.md
xxx-final.md
xxx-pycharm-fixed.md
```

这些名称只适用于临时下载文件。

正式放入仓库时仍应使用：

```text
2026-08-16-tower-latin-pgd-eq61-64-enrichment-stage-summary.md
```

即：

\[
\boxed{
\text{filename describes research content, not debugging history}
}
\]

---

# 23. 本次错误的根本原因总结

本次问题并不是数学推导错误，而是三个不同层次被混在了一起：

## 层次 1：LaTeX 语法正确性

例如：

```latex
\boldsymbol{}
```

在完整 LaTeX 环境中是合理的。

---

## 层次 2：Markdown parser 行为

例如多行：

```text
-
```

可能被解释成列表。

---

## 层次 3：PyCharm Preview 数学 renderer 支持范围

例如：

```latex
\boldsymbol
```

当前 renderer 不支持。

因此以后不能只问：

> “这个 LaTeX 写法标准吗？”

还必须问：

> “这个写法在当前 PyCharm Markdown renderer 中稳定吗？”

---

# 24. 后续项目 Markdown 的统一规范

从本阶段开始，建议 Offshore Wind Turbine LATIN-PGD 项目统一采用以下规则：

1. 行内公式只用 `$ ... $`；
2. 行间公式只用 `$$ ... $$`；
3. 行间公式尽量在一个物理行完成；
4. 不使用 `\[ ... \]`；
5. 不使用 `\( ... \)`；
6. 不使用 `\boldsymbol`；
7. 离散向量使用 `\vec{}`；
8. 避免 display math 内独立行首 `-`、`+`；
9. 长公式优先拆成多个定义式，而不是多行 LaTeX block；
10. `\boxed{}`、`\text{}`、`\mathbb{}`、`\mathcal{}`、`\operatorname{}` 当前可以使用；
11. 公式与正文、标题、列表之间留空行；
12. 保存后必须检查 PyCharm Preview；
13. 最终提交仓库前删除调试说明；
14. 正式文件名不使用 `fixed`、`final`、`pycharm` 等调试后缀；
15. 代码自动生成 md 时加入基本格式 sanity checks；
16. 带语言标记的 fenced code block 必须保证示例代码自身变量与导入完整，避免 PyCharm 静态检查产生 `Unresolved reference`；
17. 故意错误的反例代码不要标记为 `python`，统一使用 `text` 代码块，避免 IDE 对反例执行静态检查。

---

# 25. 推荐的最终写作模板

后续理论阶段总结可以直接使用下面的模板：

```markdown
# 标题

**日期：YYYY-MM-DD**  
**项目：Offshore Wind Turbine and LATIN-PGD**  
**当前分支：`feature/offshore-wind-turbine-tower-fatigue`**

---

# 1. 阶段定位

正文中使用 $x$、$t$ 等行内公式。

关键公式：

$$ \boxed{ q(x,t)=\sum_{j=1}^{m}\lambda_j(t)X_j(x) } $$

---

# 2. 原论文公式

原论文 Eq. (XX)：

$$ A-B+C=0 \tag{XX} $$

其中：

- $A$：……
- $B$：……
- $C$：……

---

# 3. Tower 离散

连续形式：

$$ q(x,t) $$

离散后：

$$ \vec q(t)=[q_{egf}(t)] $$

---

# 4. 阶段结论

$$ \boxed{ \text{current conclusion} } $$

---

# 5. 下一阶段

下一阶段只处理 Eq. (XX+1)。
```

---

# 26. 最终结论

本次格式调试得到的最重要经验是：

\[
\boxed{
\text{科研 Markdown 文档必须同时满足数学正确、Markdown 稳定、IDE renderer 兼容三项要求}
}
\]

在当前 PyCharm 环境下，最稳定的策略是：

```text
inline math  → $ ... $
display math → $$ formula $$   （单物理行）
vector       → \vec{}
```

并避免：

```text
\[ ... \]
\( ... \)
多行 $$ block
\boldsymbol{}
```

后续所有 LATIN-PGD 阶段总结文档优先遵守这一规范，以避免再次出现同类格式问题。
