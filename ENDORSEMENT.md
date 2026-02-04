# arXiv endorsement（不线下找人版）+ Zenodo DOI（推荐流程）

你现在的情况（20 岁、无论文、暂时无法走 arXiv“自动 endorsement”）属于 **必须走 personal endorsement** 的典型场景。好消息是：**全程线上即可**，不需要线下遇到人。

> 重要背景：arXiv 2026-01-21 更新政策，新提交者不再仅靠“机构邮箱”自动通过。  
> - 说明页：`https://info.arxiv.org/help/endorsement.html`  
> - 政策公告：`https://blog.arxiv.org/2026/01/21/attention-authors-updated-endorsement-policy/`

---

## 0. 先做一件事：准备“固定版本链接”（强烈推荐）

**目的**：让 endorsers 一点开就能看到“冻结版本 + 可复现 + 可引用”，而不是一个随时变动的大仓库。

推荐选项（二选一）：

- **Zenodo DOI（最佳）**：GitHub Release 自动归档到 Zenodo → 生成 DOI
- **GitHub Release（最低门槛）**：即使不连 Zenodo，也能提供固定 tag + release 页面

你仓库里已经有一篇很适合做“首篇预印本”的窄稿：`papers/hubble_tension_environment_v0.md`（以及可复现脚本与输出）。

---

## 1. Zenodo DOI（GitHub 集成）最短流程

1. 登录 Zenodo，打开 GitHub integration（选择要归档的仓库）
2. 回到 GitHub：打一个 tag（例如 `hubble-v0.1`），创建 Release
3. Zenodo 会自动抓取该 Release 并生成 DOI（通常几分钟内）

建议在 Release 描述里放三样东西：

- **一句话主结论**：\(H_0(\delta)=H_{0,\mathrm{ref}}(1-\beta\delta)\)，当前 \(\beta\approx0.171\)（其中 \(H_{0,\mathrm{ref}}\) 为参考锚点，近似取 CMB/Planck 推断值）
- **复现命令**（见下文的“给 endorsers 的 3 行复现”）
- **关键图**：`math/hubble_env_scan/plot_v1.png`

> 备注：本仓库已添加 `CITATION.cff` 与 `LICENSE`，方便 GitHub/Zenodo 生成规范引用信息。

---

## 2. arXiv endorsement：全程线上、最小社交成本

arXiv 官方流程概括如下（不需要认识任何人）：

1. 在 arXiv **Start a new submission**，选择你要投的主分类（建议先只选 1 个）
2. 系统会给你发一封 **endorsement request email**（包含链接/代码）
3. 你去找 **确实有资格 endorse 的作者**，给他们发很短的邮件 + 你的 endorsement link/code

官方说明里有一个非常关键的按钮：

- 在任意一篇 arXiv 论文摘要页底部，点 **“Which of these authors are endorsers?”**  
  你就能看到哪些作者对该领域有 endorsement 权限。

> 官方说明：`https://info.arxiv.org/help/endorsement.html`

### 推荐主分类（给哈勃张力这篇）

优先建议：`astro-ph.CO`（cosmology）

原因：主题最贴合；且你在首投阶段**不要跨多个分类**（跨得越多，越难解释“我为什么属于这个领域”）。

你可以直接把这一页摘要发给 endorsers（建议配合 DOI/Release 链接）：

- `ONEPAGER_HUBBLE_TENSION.md`

---

## 3. 2/11 起的语言新规（你需要用英文稿）

arXiv 将要求**所有新提交必须包含完整英文版本**（可选择只交英文，不必中英双语）。  

说明页：`https://info.arxiv.org/help/faq/multilang.html`  
公告：`https://blog.arxiv.org/2026/01/13/non-english-paper-submission-guidelines/`

因此：给 endorsers 的版本也建议用 **英文稿**（通过率更高）。

---

## 3.1 arXiv 上传包（推荐 TeX 源码，避免 PDF-only 风险）

arXiv 虽然允许 PDF-only，但 **由 TeX/LaTeX 生成的 PDF 往往会被要求改交源文件**（官方 FAQ 有明确提示）。
因此建议你直接上传 TeX 源码包。

本仓库已准备好一个可直接上传的 zip（会把 `main.tex` + 关键图复制到 `fig/` 并打包）：

```bash
python papers/build_arxiv_hubble_tension_environment.py
```

输出：

- `papers/_arxiv_hubble_tension_environment/hubble_tension_environment_arxiv_v0.zip`

如果你想发给 endorsers 一个更易读的 PDF（而不是 zip），在该目录下编译即可：

```bash
cd papers/_arxiv_hubble_tension_environment
pdflatex -interaction=nonstopmode hubble_tension_environment_arxiv_v0.tex
```

会得到：`hubble_tension_environment_arxiv_v0.pdf`

---

## 4. 给 endorsers 的“3 行复现”（你可以直接贴在邮件里）

在仓库根目录：

```bash
pip install -r requirements.txt
python math/hubble_env_scan/run_hubble_env_scan_v1.py
python math/spectral_graph/run_relational_hubble_beta_poisson_3d_v0.py
```

期望输出（固定路径）：

- `math/hubble_env_scan/results_v1.json`
- `math/hubble_env_scan/plot_v1.png`
- `math/spectral_graph/relational_hubble_beta_poisson_3d_results_v0.json`
- `math/spectral_graph/relational_hubble_beta_poisson_3d_plot_v0.png`

---

## 5. 冷邮件模板（尽量短）

### 模板 A（最短版）

Subject: arXiv endorsement request (first submission to astro-ph.CO)

Dear Prof. [Name],

I’m preparing my first arXiv submission in astro-ph.CO and would like to request an endorsement.

Title: Environmental correction to the Hubble tension (v0)
One-line summary: H0 depends on observer overdensity via H0(δ)=H0_ref(1-βδ), with β≈0.171 (H0_ref anchored by CMB/Planck; reproducible pipeline + falsifiable predictions).
Fixed version: [Zenodo DOI or GitHub Release link]
Endorsement code/link: [from arXiv endorsement request email]

If you’re not comfortable endorsing, no worries at all.
Best regards,
[Your name]

### 模板 B（更容易被认真对待版）

Subject: arXiv endorsement request (astro-ph.CO) + reproducible code

Dear Prof. [Name],

I’m preparing my first arXiv submission in astro-ph.CO and would like to request an endorsement.

Manuscript (fixed): [Zenodo DOI or GitHub Release link]
Key figure: H0 vs δ with MC band (plot_v1.png)
Reproduce (3 commands):
  pip install -r requirements.txt
  python math/hubble_env_scan/run_hubble_env_scan_v1.py
  python math/spectral_graph/run_relational_hubble_beta_poisson_3d_v0.py

Endorsement code/link: [from arXiv endorsement request email]

Best regards,
[Your name]

---

## 6. 礼仪与成功率建议（很关键）

- **别群发**：一次找 3–5 位最相关的 endorsers 足够
- **只投一个主分类**：首投阶段尽量不要 cross-list
- **发“固定版本”**：Zenodo DOI / GitHub Release 比“仓库主分支”强很多
- **把“可证伪点”写在一行里**：比如“若在 δ>0 过密环境测得 H0 系统性高于 H0_ref (≥2σ)，模型失败”

