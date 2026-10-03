# 分析计划（冻结版）
项目：脓毒症中抗原经历 CD4 T 细胞程序的隔室依赖性重塑
底层数据：GSE285325（大鼠 CLP 胸导管淋巴液 scRNA-seq，Cell Reports 2025;44:115469）
冻结日期：2026-09-14　　状态：G0A–G0B 已执行，G0C/P0/P1 规则已锁、待输入

本文件的作用：把所有判定规则在看到结果之前固定下来。任何已执行步骤的阈值不得回改；
新增分析必须在本文件中先写规则、再运行。执行记录见 findings.md。

**术语边界（重要）**：本文件不是 preregistration。G0A–G0B 在本文件定稿前已执行，
因此对外表述只能是 "The remaining analyses (G0C, P0, P1) were prospectively specified and
time-stamped before execution."，不得写成 preregistered study。
外部存档时**必须保留第 5 节的已执行结果**——正是它标出了前瞻与回顾的边界，
删掉会让整份文件看起来像全程预注册，反而失去可信度。
存档方式：OSF Registration（不可变）或 Zenodo 版本 DOI（可证时间顺序）。

---

## 1. 状态的操作定义（已冻结）

- 仅用 Sham 三只动物定义状态，CLP 不参与任何 feature selection
- conventional CD4 = `Cd3e/Cd3d/Cd3g 任一 > 0` 且 `Cd8a = 0` 且 `Cd8b = 0` 且非 B（Cd79a/Ms4a1/Mzb1 全 0）
- 在此群内 leiden 1.0 重聚类；先以 Treg 模块分（Foxp3, Il2ra, Ctla4, Tnfrsf18, Ikzf2,
  Tnfrsf4, Nt5e, Itgae）识别 Treg 簇，再在非 Treg 中取 Icos 最高簇
- **主状态 = cluster 3**（Icos⁺ 61.6%、Foxp3⁺ 6.2%、Il2ra⁺ 7.4%，三动物 8.0/9.4/11.5%）
- **排除**：cluster 6（Treg）与 cluster 2（126 细胞的中间态）不进入 discovery 与 comparator；
  cluster 2 仅作 sensitivity（`cluster3 + cluster2` 是否改变结论）

## 2. signature 的选取规则（已冻结）

统计单位 = animal × state 的 pseudobulk（counts 求和 → CPM → log2），不以细胞为 n。
入选须同时满足：
1. 三只 Sham 每只 Δlog2 > 0.25
2. 三只中位数 ≥ 0.50
3. 在 cluster 3 中的检出率 ≥ 0.10
4. 属于冻结后的三物种 1:1 直系同源 universe（primary-frozen 分支）

技术黑名单：`^Mt-`、`^Rps`、`^Rpl`、`^Hba`、`^Hbb`、`Malat1`
主动排除：`Icos`、`Foxp3`、`Il2ra`、`Ctla4`、`Tnfrsf18`、`Ikzf2`
细胞周期基因保留（留作 sensitivity）

两个分支：
- **Exploratory（全基因）**：仅用于数据审计与规模估计，不参与最终 signature
- **Primary-frozen**：先锁 rat–mouse–human 1:1 universe，再在其中从头选

## 3. comparator（已冻结）

主分析 = memory/antigen-experience matched：
同一动物内 1:1 最近邻匹配，匹配变量为 (MemoryScore − NaiveScore)、log nCount、log nFeature，
caliper 0.5 SD，无放回。必须报告匹配前后的标准化均差与匹配成功率。
Memory 模块：S100a4, S100a6, S100a11, Anxa1, Ahnak, Cd44, Itgb1, Lgals1, Crip1, Vim, Emp3, Capg, Ccl5, Id2
Naive 模块：Sell, Ccr7, Lef1, Tcf7, Bach2, Satb1, Actn1, Igfbp4, Foxp1, Klf2, Dapl1, Rgs10

## 4. 判定阈值（已冻结，不得回改）

**G0B-2**
- PASS：≥50 基因 且 三只 held-out AUROC 均 ≥0.75
- GREY：30–49 基因 或 最低 AUROC 0.70–0.75
- FAIL：<30 基因 或 任一 AUROC <0.70

**ambient robustness**
- `retention = |S_uncorrected ∩ S_corrected| / |S_uncorrected|`（分母 = 155）
- Spearman ρ 在两版共同可评估的**全部 eligible genes** 上比较 Δlog2FC 中位数
- 稳健 = retention ≥0.60 且 ρ ≥0.70；任一不达标则以校正版为主，155 版降为 provisional
- effect sign concordance 仅描述，不设门槛
- 仅有 filtered 矩阵时，此步只能称 **ambient-sensitive gene exclusion sensitivity**，
  不得称 ambient correction；真正的 SoupX 需 raw_feature_bc_matrix

## 5. 已执行结果（详见 findings.md）

| Gate | 结果 |
|---|---|
| G0A | PASS |
| G0B-1 | PASS（717 基因，provisional；置换 FDR 0.0006；LOO AUROC 0.985/0.970/0.921；去 Sham1 重叠 88.4%） |
| G0B-2 | **STRONG PASS**（155 基因；LOO 0.928/0.872/0.889；深度敏感性 ρ=0.780；匹配成功率 49–62%） |
| ambient exclusion | PASS（retention 0.974、ρ 0.965、sign 0.983）——弱检验 |
| G0B-3 | **assay-informative failure**，封存不再调参 |
| G0C | **PASS**（2026-09-14 执行）：universe = rat–mouse–human 1:1 共 15,869；衰减 155→135→129→**107**；受限 universe 内原 cluster3 保留 90.1%（新簇 n=1,495，Icos⁺ 63.7%、Foxp3⁺ 7.4%、Il2ra⁺ 7.7%，三动物 8.3/9.7/11.1%）→ 状态被保留，表述不降级 |

**G0B-3 封存理由（已量化）**：从 155 基因中取检出率最接近 Icos(0.616) 的六个基因作探针，
用同一 detected/undetected 设计，只能找回程序中的 0.6–2.6%。Icos 的 0% 落在该区间内，
故该设计在此深度下无检出能力，阴性结果不可解释。不得以更换阈值、匹配方法或 cut-off 重试。

## 6. 待执行步骤的规则（先锁后跑）

### G0C：跨物种冻结
输入 `00_docs/ortholog_rat_mouse_human.csv`（第一列 rat symbol）。
先锁 universe 再选 signature，规则同第 2 节。必须输出 `ortholog_attrition.csv`：
候选数 → rat/mouse 1:1 → rat/mouse/human 1:1 → 最终冻结数。
**并须验证**：在受限 universe 内重做 Sham CD4 降维聚类，报告原 cluster3 的 confusion matrix、
细胞保留率、Icos/Treg 模块表型、三动物组成。若 cluster3 在受限 universe 中不再成群，
则表述降级为 "a reduced cross-species transcriptional projection of the original state"。

### P0：与人 NR4A2⁺ central-memory marker 的重叠

**背景 universe（冻结）**：
`P0 universe = G0B-2 中所有可进入 signature selection 的 eligible genes ∩ rat–mouse–human 1:1 ortholog universe`
即同时满足：在 GSE285325 中可评估、不在技术 blacklist、不属于主动排除基因、
通过预定义的表达/检出 eligibility、且有三物种 1:1 直系同源。
**155-gene 与 717-gene 两层必须使用同一个 background universe**，否则两层的富集不可比较。
说明（写入 limitation）：人类 marker 列表本身是在人类全基因 universe 上定义的，
限制到我们的 universe 会改变 K，这是不可避免的，需在文中说明。

**marker 列表长度敏感性（现在冻结，非事后添加）**：若 Supplementary Table 6 提供排序统计量，
额外报告 top 50 / 100 / 200 三个截断下的富集倍数与 P；主结论以完整列表为准，
三个截断用于说明结论不依赖列表长度。若该表无排序统计量，则仅用完整列表并说明原因。

必须三层同时报告，且每层给出超几何随机期望：
1. 155-gene residual（扣除记忆身份后）↔ marker 列表
2. 717-gene pre-matching 版 ↔ marker 列表
3. 两者相对随机基线的富集倍数
**判读**：有信息量的是第 1 层相对基线的富集。若第 2 层高而第 1 层落回随机水平，
说明共享的只是记忆身份，跨物种线索需重新评估。

### P1：GSE279451 CITE-seq 的 protein–transcript 耦合
**前置条件（未满足则不启动）**：直接核实 ADT panel 含 CD278/ICOS ——
查 GEO 的 ADT feature reference 文件或原文 Methods 抗体表。综述转述不作数。
设计：只取作者注释的成人 CD4 T 细胞，不重新整合；限定在 central-memory / antigen-experienced 隔室内。
主读数 = **每个 donor 一个 Spearman ρ**(ADT-ICOS, frozen program UCell score)，donor 为 n。
**ADT 标准化（冻结）**：primary = **CLR，且按 feature 跨细胞计算**（Seurat 的 margin = 2），
不是按细胞跨 feature（margin = 1）——本分析要比较同一蛋白在细胞间的相对水平，
两种 margin 会给出不同答案，故必须指定。
dsb 仅在确实能拿到 background droplets 与合适 isotype 对照时作为 sensitivity；
GEO 当前公开的是每样本处理后的 feature/count 矩阵，通常不含 background droplets。
报告 isotype 背景水平。

**主比较（冻结）**：`ICU non-sepsis (n=5) vs sepsis (n=32)` 为 primary；
healthy (n=3) 仅作 descriptive/reference。检验用 Wilcoxon 秩和（donor 为 n）。
注意 5 vs 32 极不平衡，对照侧功效很低——因此**可解释性主要由控制蛋白承担**，
而非由 ICOS 组间比较的 P 值承担。
**控制蛋白为主设计的一部分**：同一 donor、同一细胞集合内同时计算
- 稳定对照：CD3 或 CD4（各组之间 ρ 不应变化）
- 阳性对照：CD45RA 或 CD27（应随记忆分化稳定共变）
**判读规则**：只有当对照蛋白的 ρ 在 healthy / ICU non-sepsis / sepsis 三组间保持稳定、
而 ICOS 的 ρ 单独下降时，方可称 uncoupling；若所有蛋白 ρ 同降，判为技术梯度，不作生物学解释。

### P2：GSE290679 等人群泛化 —— 最后执行。

## 7. 新动物实验（设计已定，待 Figure 2E 数据估算样本量）

2×2 配对隔室设计：**Sham 6–8 只/组 + CLP 6–8 只/组**（不是总数），
每只同时采胸导管淋巴液与外周血，统一 24 h，同批 scRNA 流程；
同一样本分 aliquot 做 CD3/CD4/CD25/ICOS/FOXP3 流式。
主模型：`StateScore ~ Condition * Compartment + (1|Animal)`，主检验 H0: β_interaction = 0。
**不上 rat CITE-seq**（抗体克隆与背景需另做 pilot，风险不对等）；
protein–RNA linkage 由人类 CITE-seq（P1）承担。
样本量由 Figure 2E 的逐动物配对流式数据做 simulation-based power 后确定，
不用现有 scRNA 估（refA/refB 问题使可用样本实为每组 2 只）。

## 8. Figure 2E 交互分析（标度与模型先定）

结局 p = ICOS⁺CD4 / all CD4。

**数据可用性决策树（由数据保存形式决定，不看结果）**：
- 若拿到**原始分子分母 event counts**（ICOS⁺ 与 ICOS⁻ 的事件数）→ primary 为
  **过离散校正的二项 GLMM**：`cbind(ICOS+, ICOS-) ~ Condition*Compartment + (1|Animal) + (1|Observation)`，
  或等价的 beta-binomial。**不得使用朴素二项 GLMM**：流式每样本数万个 event，
  朴素二项会把技术重复当成生物学重复，任何微小差异都会得到 p<1e-10，
  必须以观测水平随机效应或 beta-binomial 吸收过离散。
- 若只有**逐动物百分比** → primary 为下面的 logit LMM。
- 绝对百分点的 difference-in-differences 一律作 sensitivity。

**Primary（仅有百分比时）：logit 标度** `logit(p_ij) = β0 + β1·Condition + β2·Compartment + β3·Condition×Compartment + b_i`
（b_i 为动物随机截距，前提是 lymph 与 blood 确为同一动物配对——须先核实）。
主检验 H0: β3 = 0。绝对百分点的 difference-in-differences 作 sensitivity，
ratio-of-ratios 作描述性补充。不得因某一标度不显著而改用另一标度。

**存档边界更新（2026-09-14）**：G0C 已执行，故对外表述改为
"P0 and P1 were prospectively specified and time-stamped before execution."
G0A–G0C 均属已执行，第 5 节如实列出。

## 9. 当前关键路径

1. ~~`00_docs/ortholog_rat_mouse_human.csv`~~ 已到位（Ensembl 116 / GRCr8，2026-09-14）；G0C 已完成
1b. `00_docs/human_NR4A2_CM_markers.csv`（阻塞 P0）
2. 邮件：raw Cell Ranger 输出（阻塞真正的 ambient correction）+ Figure 2E FCS/pairing ID
   （阻塞第 8 节与新实验的样本量）+ 统一重跑的整合对象（可解除 refA/refB 限制）
3. GSE279451 的 ADT panel 核实（阻塞 P1）
