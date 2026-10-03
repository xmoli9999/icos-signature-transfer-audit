# Gate 2 预注册草案 —— composition rather than within-state suppression

Status: **DRAFT v0.9（假设与分析规则已定稿；待锁定独立队列后升为 FROZEN v1.0）**
Drafted: 2026-09-15，在 Gate 1 揭盲之后、在接触任何 Gate 2 候选队列数据之前。
Provenance: 本假设由 Gate 1（GSE290679）的**探索性**分解生成，不是 Gate 1 的预设假设。
Supersedes: 无。与 [[gate1_SAP.md]] 并列，不得回溯修改 Gate 1 任何内容。

---

## 0. 假设从何而来（如实记录）

Gate 1 的预设主检验为 null。其合并 memory CD4 均值的负号来自亚群构成，而非亚群内程序水平：

| | 值 |
|---|---|
| CD4 T CM 内 β | +0.00160（P = 0.404） |
| CD4 T Eff/EM 内 β | +0.00144（P = 0.440） |
| Eff/EM 占 memory CD4 的比例 | CI-Sep 0.420 vs CI-NS 0.607，Δ = −0.186，Wilcoxon P = 0.014 |

该比例读数在 Gate 1 SAP 中被预先规定为 **secondary/exploratory**，
因此以下假设是**探索性生成**的，必须在一个独立队列中预注册后检验才具有验证效力。

## 1. Gate 2 主假设

> Does human sepsis alter the composition of conventional memory CD4 states — particularly
> depletion of Eff/EM relative to CM — rather than suppressing the ICOS-associated
> transcriptional program within a fixed memory-CD4 state?

拆成两个层次，**两者都要检验、且都要报告**：

**H1（composition，primary）**：sepsis 患者 memory CD4 中 Eff/EM 占比低于 critical-illness control。
方向预设为**下降**；推断双侧。

**H2（within-state，co-primary 的证伪层）**：在 CM 与 Eff/EM **各自内部**分别计算冻结的
107-gene depth-calibrated 分数，预期**不再**出现 pooled-memory 那样的下降。
H2 是一个预期为 null 的检验，因此必须预先规定等价性边界（见 §5），
不得用"P > 0.05"直接宣称"无差异"。

支持性结果的形态是：**Eff/EM ↓，而 CM 与 Eff/EM 内 signature 不 ↓**。
结论句（若成立）：apparent loss of an ICOS-associated CD4 program at the compartment level
reflects redistribution of CD4 memory states rather than uniform transcriptional silencing
within those states.

## 2. 队列资格标准（**在查看任何候选队列的表达数据之前冻结**）

必须同时满足：
1. 人类外周血 scRNA-seq，patient-level 可识别，每组 ≥ 15 例；
2. 含明确的 sepsis 组与**危重非脓毒症对照**组（healthy-only 对照不满足，因为无法把
   "脓毒症特异"与"危重病共性"分开）；
3. 提供作者标注的 CD4 记忆亚型（CM / Eff-EM 可区分），或提供足以复现该标注的原始矩阵；
4. **未参与 P0 marker 生成**：`GSE279451` / `GSE279452` 因此被明确排除在独立验证之外，
   最多用于 contextual extension；
5. 采样时相可界定（早期 vs 晚期须能分层或至少能报告）。

候选队列的筛选只依据元数据（GEO 样本页、论文 Methods），**不得先跑打分再选队列**。
一旦锁定队列，本文件升为 FROZEN v1.0，此后任何改动记为 protocol deviation。

## 3. 统计单位与模型

- 统计单位是**患者**，每例一行；不做 cell-level 检验。
- H1 主模型：`empirical_logit(Eff/EM) ~ group + z(median nFeature)`，OLS + HC3 95% CI，
  其中 **empirical logit 现在冻结为** `log[(n_EffEM + 0.5) / (n_CM + 0.5)]`。
  采用 +0.5 连续性校正而非裸 logit，是为了在某例患者某一亚群细胞数为 0 时
  **不必在首次接触数据时临时决定 pseudocount**——那会留下一个不必要的分析自由度。
  该 pseudocount 现在写死，不得在看到新队列后调整；
  比例的 logit 变换在此处是 primary scale（与 Gate 1 对 proportion 的处理不同之处在于：
  这里比例本身就是预设的 primary endpoint，而非合并分数的附属描述）。
  每例 memory CD4 细胞数 < 100 的患者按 §4 规则处理。
- H2 模型：在 CM、Eff/EM 内分别 `mean ΔUCell ~ group + z(median nFeature)`，OLS + HC3。
- 深度校准沿用 Gate 1 Amendment 3 的完整规则（K=100、seed 20260915、
  detection frequency × mean normalized expression 双重匹配、signature 基因不入对照池），
  控制集在**新队列自身**的 primary memory CD4 内重新构建，不移植 GSE290679 的对照集。
- **107 个基因不变**：不增删、不重新加权、不重新选择。

## 4. 预设的质量与纳入规则

- 每例患者 memory CD4 < 100 细胞者排除；该阈值现在冻结，不得事后调整。
- 细胞类型标签覆盖度 < 90% 时，转用与 Gate 1 相同的 raw fallback 流程。
- 报告完整 attrition：候选样本 → 通过 QC → 进入分析的逐级人数。

## 5. 等价性边界（H2 用）

H2 预期为 null，故预先规定：若 CM 与 Eff/EM 内 β 的 95% CI 完全落在
±0.0058 ΔUCell 单位（= Gate 1 在 n=19/19、power 80% 下的最小可检出效应）之内，
判为 "within-state program preserved within the resolvable range"；
否则只报告点估计与 CI，**不得宣称等价**。

**该边界的性质必须如实标注**：±0.0058 是一个 *Gate-1-derived resolvable-range boundary*，
来自 Gate 1 的功效与方差结构，**不是**预先定义的最小生物学重要差异（MCID）。
因此结论只能写成"在本设计可分辨的范围内未见下降"，
不得写成"证明了生物学等价"或"程序未受影响"。

## 6. 判定规则

- H1 方向相反（Eff/EM 在 sepsis 中升高）→ 本条主线停止，不得改写假设方向。
- H1 成立但 H2 亦显示 within-state 下降 → 结论降级为
  "composition shift 与 within-state suppression 并存"，不得只报前者。
- H1 与 H2 均 null → 如实报告为第二次 negative validation，项目转为方法学/阴性结果报告。

## 7. 明确禁止

- 不再对 GSE290679 追加任何以寻找阳性结果为目的的分析。
- 不得用 GSE279451/279452 充当本假设的独立验证队列。
- 不得在看过新队列的结果后回头修改 §1–§6 的任何一条。

---

## 附录 A. 第一轮队列检索结果（2026-09-15，outcome-blind，仅元数据）

检索范围：GEO 及相关公开队列。**未接触任何表达矩阵**，仅读 GEO 样本页与论文 Methods。

| 队列 | Sepsis / critical non-sepsis | 判断 | 原因 |
|---|---:|---|---|
| GSE290679 | 19 / 19 | **不合格：已用于 Gate 1** | 设计完全符合（PBMC，48–72 h ICU 采样），但已是 Gate 1 队列 |
| GSE279452 | 成人/儿童 sepsis 与 ICU controls | **硬排除** | P0 marker 来源论文数据；其论文确实报道 central-memory CD4 亚群，正因诱人更不能充当独立验证 |
| GSE279451 | 32 / 5 | 不合格 | ICU 对照仅 5，且同属 P0 来源 |
| SCP548（Reyes） | ICU-SEP 8 / ICU-NoSEP 7 | 不合格 | 匹配的 ICU 对照设计仅 8 vs 7；并入 UTI non-sepsis 或 healthy 将违反冻结的 control 定义 |
| GSE167363 | 5 / 0 | 不合格 | 仅 2 healthy controls，样本量不足 |
| GSE175453 | 4 late sepsis / 0 | 不合格 | 仅 5 healthy controls，且为 post-sepsis day 14–21 |
| GSE151263 | 7（4 sepsis-only + 3 sepsis+ARDS）/ 0 | 不合格 | 主要是 ARDS 问题 |
| GSE217906 | 2 acute sepsis / 0 | 不合格 | 4 PICS + 2 acute sepsis + 3 healthy，设计不同 |
| GSE342074 | 3 / 0 | 不合格 | 2026 新数据，仅 3 sepsis + 3 healthy |

**结论：本轮未发现合格队列。Gate 2 保持 DRAFT v0.9，不锁队列。**

## 附录 B. 第二轮检索的范围决定（冻结）

下一轮**只放宽数据仓库来源，不放宽任何科学资格标准**：
继续以完全相同的 outcome-blind 规则检索 Single Cell Portal、dbGaP、EGA、
GSA-Human/HRA、cellxgene 等公开或受控访问仓库。扩大检索范围不属于修改假设。

**明确拒绝的两条降标准做法**（记录在案，以免日后回头）：
1. 不得把"≥15/组"下调以迁就 SCP548 的 8 vs 7；
2. 不得用 healthy controls 替代危重非脓毒症对照——那会使"脓毒症特异"与"危重病共性"无法分离。

若第二轮仍无合格队列，**宁可暂缓 Gate 2**，也不降低验证标准。

---

## 附录 C. 第二轮队列检索结果（2026-09-15，仅扩仓库，不改标准）

检索：Broad Single Cell Portal、cellxgene Discover、dbGaP、EGA、GSA-Human/HRA、
ArrayExpress/BioStudies、Zenodo，以及 2023–2026 年数据不在 GEO 的脓毒症单细胞论文。
仅读 accession 页与论文 Methods，**未接触任何表达矩阵**。

**结论：仍未发现合格队列。**

按接近程度排序的 near-miss：

| 队列 | Sepsis / crit. non-sepsis / healthy | 失败于 | 说明 |
|---|---:|---|---|
| Kwok 2023，EGAS00001006283 / EGAD00001010927（衍生数据 Zenodo 7723202） | 26 / **7**（心脏术后 D1）/ 6 | 标准 1、3 | 设计最接近：有真正的危重非脓毒症对照、采样时点明确（ICU D1/3/5）。但对照仅 7 例；且已发表标注中 CD4 记忆为单一 `LTB+IL7R+` 群，无 CM/EM 划分（Zenodo 的计数矩阵原则上可重新标注）。另需注意：择期心脏术后与内科危重并非同一类对照 |
| Reyes 2020，SCP548 | ICU-SEP ≈8 / ICU-NoSEP ≈7 / 19 | 标准 1、3 | 经典的 ICU-SEP vs ICU-NoSEP 设计，开放获取，但两臂各仅 7–8 例，且无 CD4 记忆标注 |
| COMBAT，EGAS00001005493 | sepsis 组 / **0** / 有 | 标准 2 | 对照为 COVID-19、流感与健康人，无非感染性危重臂 |
| HRA008452（= GSE279451/452 来源） | 281 (+164) / 0 / 仅健康 | 标准 2、4 | 已硬排除，且对照仅健康人 |
| HRA004458 | 39 / 0 / 15 | 标准 1、2 | 对照为"白细胞升高、无器官功能障碍"，非危重 |
| PRJEB96265 | 14 / 0 / 3 | 标准 1、2 | 有 CM CD4 标注，但无危重对照且样本量不足 |

阴性扫描：cellxgene Discover 的 collection 索引中无任何 sepsis/critical illness/ICU 相关集合；
ArrayExpress/BioStudies 的 sepsis 全表（100 项）均为芯片或 bulk RNA-seq，无人类血液单细胞脓毒症队列；
Zenodo 仅有上述 Kwok 衍生数据。两篇 2026 年综述列出的约 14 个脓毒症血液单细胞数据集中，
除 GSE290679 外**全部**使用健康志愿者或非危重对照。

**两处未穷尽的核查缺口（须在宣布最终结论前人工补查）**：
1. **dbGaP** 因 robots.txt 与 API 限制未能穷尽检索；
2. cellxgene 的阴性结论基于 collections 索引，而非渲染后的 disease facet 查询。

## 附录 D. 由第二轮结果引出的设计选项（尚未选定，记录以备决策）

三条路，**均需在看任何表达数据之前选定并冻结**：
1. **暂缓 Gate 2**，等待合格队列出现（保持标准不动，代价是时间）。
2. **预注册的个体患者数据合并设计**：将 Reyes（8/7）与 Kwok（26/7）合并为
   34 sepsis vs 14 critical non-sepsis，队列作为固定效应进入模型。
   注意这**仍未达到 ≥15/组**，且两者的对照性质不同（内科 ICU 非脓毒症 vs 择期心脏术后），
   异质性必须预先承认。**这是一次设计变更，不是标准放宽，须显式决策而非默认滑入。**
3. **自建验证队列**：在本单位前瞻收集 sepsis 与危重非脓毒症对照的 PBMC。
   周期最长，但唯一能同时满足全部五条标准。

选项 2 若采用，须在合并之前另行冻结：队列效应的处理、对照异质性的敏感性分析、
以及"对照 14 例 < 预设 15 例"这一偏离的明确记录。

---

## 附录 E. 检索缺口封闭（2026-09-15，仅元数据）

附录 C 遗留的两处未穷尽已补查完毕，**结论不变**。

**dbGaP**（经 NCBI E-utilities `db=gap`，网页 UI 受 robots 限制）：
以 sepsis / septic shock / septicemia / critically ill / critical illness / intensive care /
ICU / bacteremia / SIRS 及其与 single cell、transcriptome 的组合检索，
5,048 个 UID 归并为 204 个 phs 母研究，其中 88 个落在脓毒症/危重病词族内。
真正属于脓毒症或危重病的只有 4 项（phs000686、phs000631、phs000334、phs001441），
**全部是受控访问的胚系 DNA 病例对照遗传学研究，无 bulk RNA-seq、无单细胞、无血液转录组**。
索引灵敏度对照：`single cell` 命中 6,953、`single-cell RNA sequencing` 1,971、
`peripheral blood mononuclear` 278 —— 说明该索引确实能检出单细胞与 PBMC 措辞，
因此这是真阴性而非关键词漏检。

**cellxgene Discover**（经 curation API 的真实 facet 查询，而非 collections 索引）：
2,226 个 dataset、390 个 collection，按每个 dataset 的 `disease` 本体字段过滤，
**全库零个** dataset 带有 sepsis / septic shock / septicemia / SIRS / bacteremia /
critical illness 术语。血液/PBMC dataset 共 152 个、涵盖 65 个不同疾病术语
（normal 123、COVID-19 38、CMV 感染 24、类风湿 3……），其中无脓毒症、无 ICU、无 SIRS。

**一处 facet 假阴性，值得记录以备审稿人追问**：COMBAT（Ahern et al., *Cell* 2022；
cellxgene collection `8f126edf-…`，EGAS00001005493）标题中含 sepsis，
但其 curated `disease` 字段只列 COVID-19 / influenza / normal，纯 facet 查询会漏掉它。
按 STAR Methods 核对其分组人数：healthy 10、COVID-19 mild 12 / severe 20 / critical 18、
community COVID 12、**influenza acute in-patient 10**、**sepsis acute in-patient 15**。
它**仍不合格**：脓毒症臂恰好 15，但唯一干净的危重非脓毒症对照（需机械通气的流感）仅 10 例；
唯一 ≥15 的危重对照是 COVID-19 critical（18），而按 Sepsis-3 那本身就是脓毒症，不是非脓毒症对照。
其脓毒症臂亦被原文描述为"住院、涵盖 severe 与 critical"，并非一律 ICU。

**最终结论（可写入论文）**：在已检索的公开与受控访问资源
（GEO、Broad Single Cell Portal、cellxgene Discover、dbGaP、EGA、GSA-Human/HRA、
ArrayExpress/BioStudies、Zenodo）中，**不存在**同时满足全部五条资格标准的独立队列。
最接近者为 COMBAT（败于危重非脓毒症对照臂 n=10）与 Kwok 2023（败于对照臂 n=7）。

## 附录 F. 决策（2026-09-15）

**Gate 2 保持 DRAFT v0.9 / PARKED，五条资格标准一字不改。**
不因"找不到数据"而事后放宽；不把 Reyes + Kwok 升级为 Gate 2 的 confirmatory validation。
理由：该合并同时触发三层偏离——对照总数 14 < 已冻结的 15；两套对照临床性质不同；
且两队列都需要**重新构建** CM/EM 标注而非使用作者现成标注。这是在改变验证对象与可比性框架。

可用的小队列另行利用：新建 **Gate 2b — cross-cohort transportability analysis**
（见 `gate2b_preregistration_DRAFT.md`），定位为**预注册的支持性/探索性**分析，
证据等级明确低于 confirmatory，**不得用于宣称 Gate 2 已通过**。

真正的 confirmatory Gate 2 首选前瞻性自建队列：本命题（Eff/EM 构成 ↓ 而固定状态内程序不 ↓）
对 comparator 定义极其敏感，内科 ICU 非脓毒症与择期心脏术后 D1 的免疫背景并不等价。
