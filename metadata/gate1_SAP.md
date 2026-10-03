# Gate 1 · Statistical Analysis Plan（冻结版）

```text
Status: FROZEN
Version: 1.0
Frozen date: 2026-09-15
Authoritative analysis plan for: Gate 1 — GSE290679
Supersedes: analysis_results/Gate1_GSE290679_analysis_plan.md (reconciled 2026-09-15)
Amendments after first expression access must be documented as protocol deviations.
```

项目：Cross-species conservation and context-dependent remodeling of an ICOS-associated CD4 T-cell state
冻结日期：2026-09-15，在 GSE290679 的任何表达矩阵被读取之前

---

## 0. 独立性声明
- **P0 marker 来源**：Nat Immunol 10.1038/s41590-025-02345-x → GSE279448 / GSE279451 / GSE279452
- **Gate 1 primary validation**：Nat Immunol 10.1038/s41590-025-02390-6 → **GSE290679**
  （9 healthy / 19 ICU non-sepsis / 19 ICU sepsis，644,147 cells）
GSE290679 未参与 P0 的 marker 建立或富集检验，可称 independent external validation。
GSE279451 / GSE279452 来自 P0 的源论文，**只能作为 contextual extension，不得称第二个独立复制**。

## 1. Prespecified direction
The rat-derived ICOS-associated CD4 state decreased after CLP; therefore the prespecified human
validation direction is lower signature activity in ICU sepsis than in ICU non-sepsis
(`sepsis < ICU non-sepsis`). Statistical inference remains **two-sided**, with directional
concordance required for successful validation. 不得因已有方向改用单侧 P 值。

## 2. Primary cellular compartment
Author-annotated **memory CD4 T cells**, defined exclusively from the original GSE290679
annotations before examining any signature score. All annotated CD4 T cells constitute the
secondary population. No reclustering or signature-informed reannotation.
- mapping 冻结在 `metadata/memory_cd4_label_mapping.csv`（两列 `original_label,is_memory`），
  由 audit 输出填写，填完方可第一次读取 expression。
- **Contingency**：若作者注释中不存在可明确映射的 memory-CD4 类别，则退回
  all annotated CD4 作为 primary，**不得**用本研究 signature 或重新聚类临时制造一个 memory 隔室。

## 3. Primary readout
Patient-level **mean** cell-level UCell score within author-annotated memory CD4 T cells.
Median 作为 sensitivity。（mean/median 的选择在此固定，不留事后自由度。）

## 4. Prespecified decomposition（不是 co-primary）
State abundance（high-score cells 占 memory CD4 的比例）与 state intensity（high-score cells
内部的评分）分别报告，用于区分组成改变与状态内强度改变。
**连续评分永远是 primary；proportion 自始即为 secondary/exploratory**，除非存在一个独立于
GSE290679 primary contrast 的、预先定义的外部阈值。不使用 cell-level 多峰性检验来决定证据等级
（64 万细胞 pooled 且患者内高度相关，此类检验的 P 值由 n 主导）。

## 5. Statistical unit
The patient/sample, never the individual cell.

## 6. Primary contrast
ICU sepsis vs ICU non-sepsis。Healthy controls 提供生物学锚定，不参与主假设检验。

## 7. Signature coverage QC（pragmatic threshold，非生物学标准）
H = ortholog attrition 完成后冻结的 unique human ortholog 集合。
`coverage = GSE290679 中可用的冻结 human ortholog 数 / |H|`
- coverage ≥ **90%** → 直接使用该命名空间
- < 90% → 先检查同一 h5ad 内的 `/raw/var` 能否恢复；能则转 raw 命名空间
- raw 亦不足 → 保留实际可测的 frozen subset，**完整报告 attrition**；
  不重新挑基因、不补新基因、不重新训练权重
命名空间（var_names vs var 中的 symbol 列）由人在 audit 输出后决定，脚本不得自动切换。

## 8. Expression source
第二步必须显式指定读取哪一个矩阵。processed 的 `/X` 可能已 scale/regress，
不得默认用于 UCell；若 `/raw/X` 为整数 counts，pseudobulk 优先使用。
该选择在 audit 后冻结，写入本文件。

## 9. 执行顺序（盲验证）
audit（不读表达值）→ 冻结 memory_cd4_label_mapping.csv 与 expression source
→ 冻结命名空间与 coverage 判定 → **此时才第一次读取 expression** → 打分 → 患者层聚合 → 检验

---

## 10. 评分方法与稳健性（自旧 plan 合并，2026-09-15）

- **Primary score**：UCell。
- **Sensitivity scores**：AUCell，以及 patient-level pseudobulk signature score。
- **稳健性要求**：方向必须在 UCell、AUCell、pseudobulk 三者间一致；
  且不得由**单一 CD4 亚群**或**少数几名患者**驱动
  （需做 leave-one-patient-out 与逐亚群分解）。

## 11. 报告要求（自旧 plan 合并）

每个对比均报告：组间差值、标准化效应量、95% CI、Wilcoxon P，
并展示**逐患者散点**（不得只给箱线图与 P 值）。

## 12. Gate 1 判定规则（自旧 plan 合并）

Gate 1 通过的条件是：ICU sepsis vs ICU non-sepsis 的方向与大鼠发现一致，
且该效应得到预先规定的敏感性分析的实质支持。
- 方向一致但 P 值边缘 → supportive but underpowered，可继续，但须如实标注功效不足
- **方向相反 → 原主线叙事停止**，不得改写假设方向或转而强调其他对比

## 13. 可用的 GEO 输入（自旧 plan 合并）

GEO 提供 processed labeled filtered AnnData（约 28 GB）、raw labeled filtered AnnData（约 11 GB）
与 `GSE290679_RAW.tar`。优先使用 processed labeled 对象，前提是它同时包含作者原始 CD4 注释
与可用于打分的表达矩阵；若第 7 节 coverage 不达标，或第 8 节判定 `/X` 已被 scale/regress，
则转用 raw labeled 对象。该选择在 audit 后冻结并写回第 8 节。

---

## 附：与旧 plan 的冲突及裁定（2026-09-15 reconcile）

| 条目 | 旧 plan | 本 SAP | 裁定 |
|---|---|---|---|
| Primary population | 全部作者注释的 CD4-lineage T 细胞 | memory CD4 为 primary，all-CD4 为 secondary | **以 SAP 为准**（frozen signature 仍含 S100a4/Itgb1 等 memory 成分，all-CD4 易受 naive/memory 组成漂移混杂）；旧写法保留为 secondary |
| 患者层聚合 | median | **mean** 为 primary，median 为 sensitivity | **以 SAP 为准**（消除 mean/median 的事后自由度） |

旧文件已重命名为 `analysis_results/Gate1_GSE290679_analysis_plan_SUPERSEDED.md`，
顶部加 superseded 标记，不再用于任何分析决策。

---

# Amendment 1 —— expression source 冻结（2026-09-15）

**性质**：本修订在**任何按分组读取表达值之前**做出，依据的是盲审计（Step 1a / 1b：
只读元数据与不涉及分组的数值结构判定）。因此属于 pre-specification 的补充，
**不是 protocol deviation**。

## A1.1 命名空间（第 7 节的回退已触发）
- processed `/X`：644,147 × **5,000**（作者的 HVG 子集），107 个冻结 ortholog 仅命中
  34 个，**coverage 31.8% < 90%** → 按 SAP 第 7 节不可用。
  `layers/counts` 与 `layers/scVI_normalized` 形状相同，同样不可用。
  （`layers/scVI_normalized` 另有 scVI 模型插补成分，本研究一律不使用。）
- **`/raw/var` 36,601 基因，107/107 命中，coverage = 100.0%** → **冻结使用 raw 命名空间**。
  无需保留子集，无 attrition 需报告。

## A1.2 表达矩阵的性质（盲检结果）
`/raw/X` 非整数。取前 200,000 个非零值，其 `expm1` 为同一常数的整数倍
（0.4115 × 1,2,3,4,5…），即
`raw/X = log1p(counts × per-cell size factor)`
= 标准的 library-size 归一化 + log1p，**不是 scaled / regressed 数据**。

## A1.3 由此冻结的用法
- **Primary（UCell）**：直接使用 `/raw/X`。UCell 基于细胞内基因秩，
  library-size 归一化与 log1p 均为单调变换，**不改变秩**，
  因此在该矩阵上的 UCell 结果与在 counts 上完全一致。
- **Sensitivity（AUCell）**：同样使用 `/raw/X`（亦为秩/排序法）。
- **Sensitivity（patient-level pseudobulk）**：需 counts 量纲。
  先验证能否由 `expm1(raw/X)` 除以每个细胞的最小正值（或 obs 的
  `_scvi_raw_norm_scaling`）还原整数 counts；
  若能还原则用还原后的 counts 做 pseudobulk；
  若不能，则以 `expm1(raw/X)` 按样本求和后转 CPM 作为 pseudobulk 的近似，
  并在文中明确说明该近似及其理由。此判定同样须在按分组读取之前完成。

## A1.4 数据集结构（供后续引用）
- 47 samples = 47 patients（`orig.ident` 与 `SCARAB_ID` 均为 47 个唯一值，
  且 sample × Group 无跨组样本）→ 不存在同一患者多样本的伪重复问题
- `Group`：CI-Sep 19 / CI-NS 19 / NHC 9；primary contrast = **CI-Sep vs CI-NS**
- `Batch` 4 个批次，`Age`、`Biologic_Sex`、`ICU`(MICU/SICU/Outpatient)、
  `Source`(感染来源 7 类) 可用于敏感性与后续 contextual 分析
- 作者注释：`celltype_coarse` 9 类、`celltype_fine` 27 类

---

# Amendment 2 —— memory CD4 mapping 冻结（2026-09-15）

同样在按分组读取表达值之前做出，仅依据作者原始注释 `celltype_fine` / `celltype_coarse`，
未使用任何 signature 评分，未做重新聚类。明细见 `metadata/memory_cd4_label_mapping.csv`。

| 群体 | 定义 | 细胞数 |
|---|---|---|
| **Primary：memory CD4** | `CD4 T CM` + `CD4 T Eff/EM` | 49,426 + 71,562 = **120,988** |
| Sensitivity 1 | Primary + `CD4 CTL` | + 20,564 = 141,552 |
| **Secondary：all CD4** | `celltype_coarse == "CD4 T"`（作者划分已排除 Treg） | **222,053** |
| Sensitivity 2 | Secondary + `CD4 Treg` | + 13,941 = 235,994 |

理由：大鼠 discovery 的状态是 ICOS-high / Foxp3-low 的 **conventional 记忆 CD4**。
CD4 CTL 属抗原经历但走终末细胞毒分化，另一条路线，并入 primary 会稀释；
完全排除又可能漏信号，故置于 sensitivity。
Treg 在大鼠端即被排除于 discovery 与 comparator 之外，人端保持一致。

四个群体的评分在同一次读取中一并算出（并集 = 235,994 细胞），
避免因多次读取产生不一致；分组比较在后续独立脚本中进行。

---

# Amendment 3 —— blinded depth calibration of the primary signature score（2026-09-15）

Before any comparison of ICU sepsis versus ICU non-sepsis was performed, blinded
quality-control analyses showed substantial associations between raw single-cell
signature scores and transcript-detection depth (Spearman ρ = 0.527 for UCell versus
nFeature_RNA and ρ = 0.627 for AUCell versus nFeature_RNA), together with marked
between-sample differences in median detected genes (median nFeature_RNA per sample
732–2,088; median detected signature genes per sample 8–24). To reduce nonspecific
detectability effects while preserving the frozen 107-gene signature, the primary
cell-level score was amended to a background-calibrated UCell score. For each cell, the
observed signature UCell score will be centered against the mean score of 100 randomly
generated control gene sets matched to the frozen signature on gene detection frequency
and mean normalized expression, using a fixed random seed and a prespecified gene
universe. No clinical group labels will be used in constructing the matched control sets.
The patient-level primary endpoint remains the mean cell-level score within the
prespecified primary memory-CD4 compartment. The original uncorrected UCell score will be
retained as a sensitivity analysis. This amendment was finalized before inspection of the
primary ICU-sepsis versus ICU-non-sepsis contrast and is therefore treated as a
prespecified blinded QC amendment rather than a post-result protocol deviation.

**The frozen 107 genes themselves are unchanged; no genes are removed, reweighted, or replaced.**

## A3.1 定义
cell 层：`ΔUCell_i = UCell_signature,i − (1/K)·Σ_k UCell_control_k,i`
patient 层 primary endpoint：`mean_i∈memoryCD4(ΔUCell_i)`

## A3.2 冻结参数（运行前锁定）
- `K = 100`
- `random seed = 20260915`
- **control universe**：`/raw/var` 中、在 **primary memory CD4（120,988 细胞）** 内
  detection frequency > 0 的基因，**排除 107 个 signature 基因本身**。
  universe 的选取规则写死在脚本里，运行后不得人工增删。
- **matching rule**：对每个 signature 基因，在与其**同时匹配**以下两项的 bin 内抽取对照：
  ① gene detection frequency（在 primary memory CD4 内的检出细胞比例）
  ② mean normalized expression（同一细胞集内的平均 log1p 表达）
  分箱：detection frequency 取 20 个分位数箱 × mean expression 取 20 个分位数箱；
  若某 bin 内可用基因不足，则向相邻 bin 逐层扩展（扩展规则写死在脚本里）。
- **重复规则**：同一 replicate 内不放回抽样；**不同 replicate 之间允许重复**抽到同一基因。
- 两项匹配量均在 primary memory CD4 细胞内计算，**不使用任何 group label**。

## A3.3 术语（Methods 用词）
本校正**不是** UCell 或 AUCell 的原生步骤。其思路与 `scanpy.tl.score_genes` 的
control-gene correction 相近，但须明确称为本研究**预先指定的
depth-calibrated / background-corrected UCell score**，
不得写成 "standard UCell/AUCell procedure"。

## A3.4 深度相关性的表述
正文表述固定为：
"raw signature scores showed substantial dependence on cellular transcript detection,
creating a potential technical/biological-complexity confounding pathway."
**不得**写成 "raw UCell was biased by sequencing depth" ——
活化 T 细胞本身即可能具有更高的 RNA complexity，不能全部归因于技术误差。

## A3.5 三层证据
| | 读数 | 地位 |
|---|---|---|
| Primary | calibrated ΔUCell（patient mean, memory CD4） | 主结论 |
| Sensitivity 1 | raw UCell（未校正） | 方向一致性 |
| Sensitivity 2 | `patient raw UCell ~ group + median nFeature` | 另一种深度处理策略下的方向一致性 |
三者方向一致则可信度显著提高；不一致则须在正文中如实呈现并降级结论。

**Option C（限定 nFeature 窗口的 depth-overlap restricted analysis）不进入主流程**，
仅在需要时作为很后面的补充敏感性；理由：各患者被保留的细胞比例不同，
且转录复杂度本身可能与真实生物学状态相关，硬切窗口会连真实状态一起切掉。

---

# Blinded QC record —— Amendment 3 执行结果（2026-09-15，揭盲前记录）

脚本 `scripts/gate1_03_calibrate.py`，输出 `analysis_results/gate1_cell_scores_calibrated.csv.gz`、
`gate1_control_sets.csv.gz`（100×107 对照基因表，可逐一复现）、
`gate1_calibration_diagnostics.json`、`gate1_patient_level_blinded.csv`。
**本步骤未读取任何 group label。**

- signature UCell 与 Step 2 最大绝对差 **7.4e-09**（同一表达源、同一 107 基因，数值一致）。
- 对照抽样：universe 28,167 基因；10,600/10,700 次抽样在原 bin 内完成（r=0），
  100 次扩展到相邻 bin（r=1），无更远扩展；共用到 6,774 个不同对照基因。
- 匹配质量：detection frequency 中位数 signature 0.0915 vs control 0.0941；
  mean normalized expression 中位数 0.1331 vs 0.1385。
- `CLDN10` 在 primary memory CD4 内检出率为 0（107 基因不变，照常参与打分，
  对所有细胞贡献同一常数；须在 limitation 中如实说明）。
- cell 层深度相关性（Spearman vs nFeature_RNA）：raw UCell **0.527 → ΔUCell 0.202**。
- **patient 层深度相关性：raw 患者均值 +0.832 → ΔUCell 患者均值 +0.537。**
  即 cell 层校正有效，但患者层残余关联仍然可观（患者 median nFeature 跨度 732–2,088，2.9×）。
- 患者层分布：n=47，每例 primary memory CD4 细胞中位 2,486（498–7,285）；
  患者均值 ΔUCell +0.0032，SD 0.0055。
- 对照集 Monte-Carlo 噪声：每细胞 100 组 UCell 的 SD 中位数 0.0197，
  100 组均值的标准误约 0.0020。**该误差由同一套对照集产生，在细胞间共享、不随细胞数平均掉**，
  但对所有患者同向，故在组间对比中大部分抵消。

---

# Amendment 4 —— depth-adjusted primary patient-level model（2026-09-15，揭盲前冻结）

冻结时点：Amendment 3 执行完毕、盲态 QC 结果已知，**但尚未读取任何 group label、
尚未查看任何 group-specific signature-score 结果**。冻结依据是 QC 中可量化的
residual depth dependence（patient 层 ρ=+0.537），不是任何组间结果。

## A4.1 Primary model（正文原文）
Primary model: among ICU patients only, the patient-level mean background-calibrated
UCell score within the prespecified primary memory-CD4 compartment will be modeled as
`mean ΔUCell ~ group + median nFeature_RNA`, where `group` is ICU sepsis versus ICU
non-sepsis and `median nFeature_RNA` is calculated within the same prespecified
memory-CD4 compartment. The group coefficient is the primary estimand. Inference will be
two-sided, while successful validation additionally requires the prespecified direction
`sepsis < ICU non-sepsis`. This amendment was finalized before inspection of any
group-specific signature-score results.

## A4.2 统计实现（运行前写死）
- `median nFeature_RNA` 作为**连续**变量进入模型；**不分组、不搜索最佳阈值、不做样条**。
- 该协变量标准化为 z-score，**仅为系数尺度稳定**，不改变推断结论。
- 主模型为线性回归（OLS），报告 **HC3 robust 95% CI** 与双侧 P——
  真正的统计单位只有约 38 例 ICU 患者，必须用异方差稳健方差。
- group 编码：`ICU non-sepsis` 为参照，group 系数即 sepsis 相对 ICU non-sepsis 的差；
  预期方向为**负**。
- 统计单位是患者；每例患者一行；不做任何 cell-level 检验。

## A4.3 证据层级（固定，不得事后调整次序）
1. **Primary**：`mean ΔUCell ~ group + median nFeature`
2. **Sensitivity 1**：未调整的 patient mean ΔUCell 组间差
3. **Sensitivity 2**：raw UCell 的 depth-adjusted model
4. **Sensitivity 3**：AUCell 对应结果，主要看**方向一致性**

## A4.4 两条记入正文的细节
- `CLDN10` 在 primary memory CD4 中检出为 0：**不删除**。记录为"frozen signature 中该基因
  在 primary compartment 无检测，因此不贡献组间差异"，107-gene definition 保持不动。
- Amendment 3 与 Amendment 4 **均在任何 primary group comparison 之前完成**，
  两者都属于 prespecified blinded amendment，不是 post-result protocol deviation。

## A4.5 已记录的执行偏离（与建议不同之处）
建议中的 "protein-coding-like gene universe" **未实现**：本数据集 `/raw/var` 仅有 `_index`
一列，无 biotype 注释；加 protein-coding 过滤需引入外部基因列表，等于在冻结规则中
植入外部依赖。实际 universe = primary memory CD4 内 detection frequency > 0 且非 signature
的 **28,167** 个基因，纯表达驱动、写死在 `gate1_03_calibrate.py` 中，运行后未作任何人工增删。
本条在揭盲前记录。

## A4.6 拒绝的方案（记录在案）
"先看两组 median nFeature 是否平衡，再决定是否调整"——**不采用**。
一旦读取 group label，预注册纯度即下降；且是否预先调整协变量应由
"协变量与 outcome 的强相关"决定，而不是由它在两组间是否显著不平衡决定。

---

# Gate 1 结果记录（2026-09-15，首次且唯一一次揭盲）

脚本 `scripts/gate1_04_gate_test.py`；输出 `gate1_gate_test_results.json`、
`gate1_patient_level_unblinded.csv`、`gate1_primary_perpatient.png`。
Amendment 3 与 Amendment 4 均已在本次运行之前冻结。

## 主结果
- 覆盖度 100%（阈值 90%）；primary memory CD4 = 120,988 细胞；CI-Sep 19 例 / CI-NS 19 例。
- **PRIMARY：β = −0.00138，HC3 95% CI [−0.00531, +0.00255]，P = 0.493。**
- Sensitivity 1（未调整）：Δ = +0.00042，95% CI [−0.00366, +0.00450]，Hedges g = +0.07，Wilcoxon P = 0.977。
- Sensitivity 2（raw UCell，调整）：β = −0.00110，P = 0.610。
- Sensitivity 3（AUCell，调整）：β = −0.00105，P = 0.648；未调整 Δ = +0.00375，g = +0.34，P = 0.321。
- 四层方向**不一致**（primary/Sens2/Sens3 为负，Sens1 与 Sens3' 为正）。
- LOPO 38 次：β ∈ [−0.00223, −0.00054]，无变号——但这是在一个本身不显著的估计附近的稳定性。

## 判定
按 §12：方向未相反，故不触发"主线叙事停止"条款；但效应量与零无法区分
（未调整 d = +0.07），**不构成对大鼠发现的阳性外部支持**。
Gate 1 判定为 **null / 不通过**，不得改写为"方向一致但功效不足"以外的任何表述。

功效：n = 19/19、α = 0.05 双侧、power = 80% 下最小可检出 d ≈ 0.93
（|Δ| ≈ 0.0058 ΔUCell 单位；患者层合并 SD = 0.0062）。
本队列只能排除 |d| ≳ 0.93 的效应，更小的真实效应无法区分。此句须进入 limitation。

## 一个必须如实报告的组成效应（Simpson 型）
- 亚群 ΔUCell 水平：CD4 T CM −0.00677，CD4 T Eff/EM +0.00811。
- 亚群内部分解（非 co-primary）：CM β = +0.00160（P = 0.404），
  Eff/EM β = +0.00144（P = 0.440）——**两个亚群内均为正**。
- 组成：Eff/EM 占 memory CD4 的比例 CI-Sep 0.420 vs CI-NS 0.607，Δ = −0.186，
  Wilcoxon P = 0.014（**secondary/exploratory 读数**，比例永不作 primary）。
- 即：合并 compartment 均值的负号完全来自两组的 memory 亚群构成差异，
  而不是亚群内部的程序水平差异。预先指定的 primary 读数（合并均值）对构成敏感——
  这一点在盲态下冻结，故结论照此报告；不得事后改用构成不变的读数并称其为 primary。

## 不得做的事（记录在案）
- 不得事后更换 primary compartment、readout 或对比以求阳性。
- 不得把上述组成差异提升为本文主假设的验证——它是探索性发现，
  若要成为主线，须在另一个独立队列中预注册后检验。
- 其他预设群体（memory+CTL β = −0.00255，P = 0.232；all-CD4 β = −0.00135，P = 0.561）
  同样为 null，一并报告。

## Gate 1 正式结论（论文用词，冻结）

The prespecified cross-species validation was negative. In the primary memory-CD4
compartment, the depth-calibrated ICOS-associated signature did not differ detectably
between ICU sepsis and ICU non-sepsis patients after prespecified adjustment for
transcript complexity (β = −0.00138, 95% CI −0.00531 to 0.00255; P = 0.493).
Sensitivity analyses were concordantly null. Given the sample size, the study excluded
only large effects and remained underpowered for small-to-moderate effects.

**不得写成"方向一致但未显著"**——未调整差值为 +0.00042（d ≈ +0.07），整体证据即 null。

图注用词：主图标题为 "Gate 1 · GSE290679 — prespecified primary validation"，
判定（primary test was null）写在图注中，不写进标题。

**GSE290679 到此封存：不再追加任何以寻找阳性结果为目的的分析。**
