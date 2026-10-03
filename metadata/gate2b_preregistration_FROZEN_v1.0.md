# Gate 2b 预注册草案 —— cross-cohort transportability analysis

Status: **FROZEN v1.0**
Frozen date: 2026-09-15（在读取 Reyes 与 Kwok 的任何表达矩阵之前）
此后任何改动记为 protocol deviation。
Drafted: 2026-09-15，在接触 Reyes 与 Kwok 两队列的任何表达数据之前。
证据等级：**预注册的支持性 / 探索性分析（supportive / exploratory）**。

## 0. 定位（首要条款，不得弱化）

本分析**不是** Gate 2，**不能**用来宣称 Gate 2 已通过，
**不得**在摘要或结论中被表述为对 [[gate2_preregistration_DRAFT.md]] 主假设的 confirmatory validation。
Gate 2 的五条资格标准保持不变、状态 PARKED；本文件是在"公开数据中不存在合格队列"
这一已记录事实之下，对现有小队列的一次透明、预注册、但证据等级较低的利用。

正文中的定位用词固定为：
"a prespecified cross-cohort transportability analysis in two small independent cohorts,
reported as supportive rather than confirmatory evidence."

## 1. 队列

| 队列 | Sepsis | Critical non-sepsis | 对照性质 | 来源 |
|---|---:|---:|---|---|
| Reyes 2020 | ICU-SEP ≈8 | ICU-NoSEP ≈7 | 内科 ICU 非脓毒症 | Broad SCP548（开放） |
| Kwok 2023 | 26 | 7 | **择期心脏术后 D1** | EGAS00001006283；衍生数据 Zenodo 7723202 |

两者均独立于 Gate 1（GSE290679）与 P0（GSE279451/452）。
合计 34 vs 14 —— **明确记录：对照 14 例低于 Gate 2 冻结的 ≥15/组**，
这正是本分析不能充当 confirmatory 的原因之一。

**对照异质性必须预先承认**：Kwok 的对照是择期心脏术后 D1，Reyes 是内科 ICU 非脓毒症。
两者虽同称 critical non-sepsis，免疫背景并不等价。该异质性在结论中必须显式陈述，
不得以"均为危重对照"一笔带过。

## 2. 两层分析（冻结）

1. **每个 cohort 内独立**估计 H1 与 H2 效应；
2. **pooled model**：`outcome ~ group + cohort + z(median nFeature)`，OLS + HC3；
3. **始终同时报告 cohort-specific effect**，不得只给 pooled β；
4. 预先检查 `group × cohort` 交互，但其用途**仅限 heterogeneity assessment**，
   **不得**据此重新筛选保留哪个 cohort；
5. **若两个 cohort 方向相反，禁止用 pooled β 把它们平均成"总体效应"**——
   此时只报告两个 cohort-specific 估计与该分歧本身。

H1、H2 的定义、empirical logit `log[(n_EffEM+0.5)/(n_CM+0.5)]`、深度校准规则
（K=100、seed 20260915、detection frequency × mean normalized expression 双重匹配、
signature 基因不入对照池、控制集在各队列自身内重建）、
以及 ±0.0058 的 Gate-1-derived resolvable-range boundary 及其非-MCID 属性，
**全部沿用 Gate 2 草案，不作任何改动**。107 个基因不变。

## 3. CM/EM 重新标注配方（关键自由度，须先冻结再看数据）

两个队列的已发表标注都不含 CM/EM 划分，因此必须重新标注。
**该配方必须完全独立于 107-gene signature，且不得依据组间差异调整 marker。**

- **Primary recipe**：以固定的外部参考做 reference mapping——
  Azimuth 人 PBMC CITE-seq 参考（Hao et al. 2021），取其 level-2 标签
  `CD4 TCM` 与 `CD4 TEM`。参考版本与映射参数在运行前写死，运行后不得调整。
- **Sensitivity recipe**：基于 canonical memory-state markers 的规则式标注
  （CCR7 / SELL 阳性且 IL7R 阳性 → CM；CCR7 / SELL 阴性且 GZMA / GZMK 阳性 → EM；
  naive 由 TCF7 / LEF1 高且 GZM 家族阴性界定并排除）。marker 清单现在冻结。
- **执行顺序**：标注**必须在计算任何 H1 composition 之前完成，且对 group label 盲态**。
  标注完成后固化存盘（逐细胞标签表 + 校验和），此后不得重做。
- **禁止**：不得为了得到"更好看"的 CM/EM 划分而重新聚类、调分辨率或改 marker；
  不得用 107-gene score 参与标注；不得在看过 H1 结果后更换 primary recipe。
- 两套 recipe 的一致性（Cohen's κ 或交叉表）作为描述性质量指标报告，
  **不作为选择 recipe 的依据**——primary 永远是 reference mapping。

## 4. 判定与报告

- 支持性结果的形态与 Gate 2 相同：Eff/EM ↓，而 CM 与 Eff/EM 内 signature 不 ↓。
- 无论结果如何，报告中必须同时出现：对照 14 < 15、对照异质性、标注为重新构建、
  以及本分析为 supportive 而非 confirmatory 这四项限制。
- 结果为阳性**不**改变 Gate 2 的 PARKED 状态；Gate 2 的通过只能由一个满足全部五条
  标准的队列（首选前瞻性自建）来实现。

## 4b. 采样时点与患者独立性（冻结）

两个队列都不是单次采样，因此在任何分析之前锁定：

- **Kwok**：primary 仅使用 **acute sepsis baseline** 与 **cardiac-surgery D1**。
  原研究含 26 sepsis、9 例 paired convalescent、6 healthy、7 cardiac-surgery controls，
  脓毒症侧有 D1/D3/D5 序列采样。
  - 每例患者**只取一个样本**：急性期最早可用的时点（有 D1 取 D1）。
  - **paired convalescent samples 不进入 primary Gate 2b**（可在描述性部分提及，不参与检验）。
  - 若某例无 D1，取其最早的急性期样本，并在 attrition 表中逐例记录该替代的发生次数。
  - **绝不把同一患者的不同时点当作独立个体。**
- **Reyes**：每例患者一个样本；若存在同一 donor 的多个分选组分
  （如 `P18` / `P18H` / `P18F` 之类后缀），须先确认其含义，
  同一 donor 的多个组分合并为一行，不得当作独立个体。
- healthy 组在两个队列中均**不参与** H1/H2 检验（与 Gate 1 对 NHC 的处理一致，仅作生物学锚定）。

## 4c. 每例最小细胞数阈值（Gate 2b 专用，冻结）

**规则：每例患者 memory CD4 ≥ 30 个细胞方可进入主分析。**

**推导**（外部标尺取自 Gate 1 / GSE290679，与 Reyes、Kwok 无关）：
Eff/EM 占 memory CD4 比例的患者间真实变异实测 **SD = 0.219**（组内合并，已扣除组间差异）；
Gate 1 实测组间差 **Δ = 0.186**。在此标尺下，三条判据各自给出：

| 判据 | 结果 |
|---|---|
| 二项抽样 SE ≤ 患者间 SD 的一半 | n ≥ 21 |
| 抽样方差 ≤ 总方差的 20% | n ≥ 21 |
| 单例 95% CI 半宽 < Gate 1 实测组间差 0.186 | n ≥ 28 |

取上界并向上取整为 **30**。留 30 而非 21，是为 §3 的 CM/EM 重新标注误分类留余量——
前两条判据假设标注准确，而在 Reyes 上标注由本研究自行完成，本身带噪声。

参考值：n=30 时抽样 SE = 0.091（占患者间 SD 的 0.42），抽样方差占总方差 14.8%，
95% CI 半宽 0.179。

**为何不沿用 Gate 2 的 <100**：n=100 时抽样方差仅占总方差 4.9%，
为把一个已经只占 5% 的噪声继续压低而剔除大半队列，是以功效换不需要的精度。
<100 对每例数千细胞的 GSE290679 无代价，对 Reyes 则是决定性的。
**Gate 2 自身的 <100 规则保持不动**；本条仅适用于 Gate 2b，并附上述推导。

**配套规则（同时冻结，防止阈值本身成为自由度）**：
- 低于 30 的 donor 进入 attrition 表并**逐例列出**，不进入主分析；
- 预先规定在 **n ≥ 20** 与 **n ≥ 50** 两个替代阈值下各重跑一次作为敏感性分析；
- 三个阈值下方向不一致时，必须如实呈现该不一致并降级结论，
  **不得事后择优选取阈值**。

**记录**：本阈值的必要性是在冻结前、仅凭 metadata 中的逐 donor 细胞计数发现的
（见附录 D）。细胞计数不是 outcome，查看它不破坏对 H1/H2 的盲态；
但"何时看了什么"已如实记录于附录 C、D。

## 5. 升 FROZEN 的前置条件（已全部满足）

1. 确认 Kwok 数据的实际可获取性（EGA 受控访问申请，或确认 Zenodo 7723202 的衍生计数矩阵
   是否足以支撑 §3 的重新标注）；
2. 确认 Reyes SCP548 的计数矩阵与患者级分组信息可用；
3. 两项确认**仅基于元数据与数据字典**，不得先行打分。
完成后本文件升 FROZEN v1.0，此后任何改动记为 protocol deviation。

---

## 附录 A. 可获取性核查结果（2026-09-15，仅元数据）

| 条目 | 结论 | 依据 |
|---|---|---|
| **Kwok 数据可用性** | **PASS** | Zenodo 7723202 **Version 1.1** 提供 whole-blood BD Rhapsody 的 RNA counts、barcodes、features 与 cell-level metadata（Version 1.0 缺 barcode/feature，1.1 已补齐）。重新做 reference mapping 所需的三件基本要素齐备，**无需**回到 FASTQ。EGAD00001010927 存在且有 DAC 可申请，但 Gate 2b 不以此为前提 |
| **Reyes 表达矩阵** | **PASS** | SCP548 公开列出 `scp_gex_matrix.csv.gz`（300.6 MB）与 `scp_gex_matrix_raw.csv.gz`（113.9 MB）；原文数据可用性声明即指向 SCP |
| **Reyes patient × 分组映射** | **PASS（含一处待人工确认）** | SCP 公共 REST API（无需登录）返回 study-scoped 分组注释：`donor_id`（65 个不同 ID）与 `Cohort`（Control、Leuk-UTI、URO、Int-URO、**ICU-NoSEP**、**ICU-SEP**、Bac-SEP）。cluster 端点按 cell 返回 `donor_id` 与 `Cohort` 的逐细胞映射，多细胞→单 donor 的 patient-level 结构成立。**待确认**：`scp_meta_updated.txt`（18.3 MB）本身需免费 Google 登录才能下载，其字面列名未经目视核对 |

**访问性质澄清**：SCP548 的下载门槛是**免费账号登录 + 条款接受，不是 DUOS**。
DUOS 控制的是 `matrix_clustered.h5ad`、`clinical_data.csv` 等对象；
Gate 2b 所需的处理后矩阵与 metadata 走 SCP 公开路径即可。
另注：SCP 研究页另附 `Clinical_Data_Reyes_et_al_NATURE_MED.xlsx`（170.8 KB，同一登录门槛）。

**人数（原始文献核实）**：Nature Medicine 2020 Fig. 3c 图注明确写出
ICU-SEP **8** 例、ICU-NoSEP **7** 例。其余队列人数未逐一核实。

## 附录 B. 升 FROZEN v1.0 的剩余动作 —— **已完成（2026-09-15）**

登录 Single Cell Portal 取得 `scp_meta_updated.txt` 并目视核对表头，
确认 patient 级 ID 列与 ICU-SEP / ICU-NoSEP 分组列存在。**结果见附录 C：PASS。**
实际列名为 `donor_id` 与 `Cohort`（与 SCP 注释名一致）；
文件 md5 `d2caac081bcf5cdc250f95ad0b605db2`，18,344,290 bytes。
全程未读取任何表达矩阵。

**本文件于 2026-09-15 升为 FROZEN v1.0**，此后任何改动记为 protocol deviation。

---

## 附录 C. Reyes metadata 核查结果（2026-09-15，仅 metadata，未接触表达矩阵）

文件：`metadata/scp_meta_updated.txt`，18,344,290 bytes，**md5 `d2caac081bcf5cdc250f95ad0b605db2`**
（经 SCP 网页 Save Link As 取得；SCP 元数据文件第二行为 `TYPE` 声明行，非数据行）。

**列名（实际核对）**：`NAME`（首列名含尾随空格）、`Cell_Type`、`Cell_State`、**`Cohort`**、
`biosample_id`、**`donor_id`**、`species(+__ontology_label)`、`disease(+__ontology_label)`、
`organ(+__ontology_label)`、`library_preparation_protocol(+__ontology_label)`、`sex`。

**判定：PASS。** patient 级 ID 列与分组列均在，且为逐细胞映射。

- 细胞总数 126,351；`donor_id` 唯一值 65。
- `Cohort` 七个取值与细胞数：Control 49,718（19 donor）、Leuk-UTI 22,491（10）、
  URO 16,565（10）、Int-URO 14,450（7）、**ICU-SEP 8,848（8 donor）**、
  **ICU-NoSEP 8,505（7 donor）**、Bac-SEP 5,774（4）。人数与 Nat Med 2020 Fig. 3c 图注一致。
- ICU-SEP donors：P633, P636, P640, P662, P669, P670, P671, P672。
- ICU-NoSEP donors：P634, P635, P639, P650, P657, P667, P668。
- **无任何 donor 跨多个 cohort**（此前担心的 H/F 后缀问题不影响 ICU 两臂：
  P17/P17H/P18/P18F/P20/P20H 等出现在 Control 与 Leuk-UTI 侧，ICU 两臂的 donor 为纯 P6xx 编号）。
- `disease` 列全为 `normal`，**不携带分组信息**；组别只在 `Cohort` 列 —— 写入数据字典，禁止误用。
- `biosample_id` 仅两值：CD45（106,545）、DC（19,806），即分选组分而非样本 ID。
- `Cell_Type` 为粗分类：Mono 58,557、T 32,341、DC 14,299、NK 9,390、B 7,970、Megakaryocyte 3,794。
  **无 CD4/CD8 划分，更无 CM/EM** —— 与 §3 的判断一致，必须重新标注。

## 附录 D. 由本次核查暴露的可行性问题（**冻结前必须决策**）

ICU 两臂内 `Cell_Type == "T"` 的**逐 donor** 细胞数：

| ICU-NoSEP | T 细胞 | | ICU-SEP | T 细胞 |
|---|---:|---|---|---:|
| P634 | 271 | | P633 | 67 |
| P635 | 303 | | P636 | 405 |
| P639 | 265 | | P640 | 129 |
| P650 | 180 | | P662 | 49 |
| P657 | 79 | | P669 | 191 |
| P667 | 529 | | P670 | 236 |
| P668 | 153 | | P671 | 1,082 |
| | | | P672 | 145 |
| 合计 | 1,780 | | 合计 | 2,304 |

这里的 T 是**全部 T 细胞**（含 CD8、含 naive）。memory CD4 只是其中一部分，
因此按 Gate 2 冻结的"每例 memory CD4 < 100 细胞者排除"规则，
**Reyes 相当一部分 donor 很可能被剔除**，个别 donor（P662 共 49 个 T 细胞）几乎必然出局。
memory CD4 的实际数目在完成 §3 的重新标注前无法得知，而重新标注需要表达矩阵。

**这是在冻结前、仅凭 metadata 发现的结构性问题，须显式决策，不得在跑出 attrition 后再改规则。**

---

## 附录 E. 事实性更正（2026-09-15，冻结后）

**性质**：数据来源记录号的更正，**不涉及任何分析规则**，因此不构成 protocol deviation；
按纪律仍公开记录。

附录 A 将 Kwok 衍生数据记为 Zenodo **7723202**。经核对，该 record 为 **Version 1.0**，
其 7 个文件中**不含** barcode 与 feature 列表（记录描述自身亦注明 v1.0 不完整）。
**Version 1.1 是另一个 record ID：7924238**（发布日 2023-05-17，共 15 个文件），
barcode 与 feature 齐备。Gate 2b 使用的是 **7924238**。

本研究所需文件（BD Rhapsody whole-blood）：

| 文件 | 字节 |
|---|---:|
| `Kwok-et-al_2023_BD-Rhapsody_whole-blood_RNA-counts.mtx` | 3,494,797,339 |
| `Kwok-et-al_2023_BD-Rhapsody_whole-blood_RNA-counts_features.tsv` | 177,721 |
| `Kwok-et-al_2023_BD-Rhapsody_whole-blood_RNA-counts_barcodes.tsv` | 9,270,174 |
| `Kwok-et-al_2023_BD-Rhapsody_whole-blood_cell-metadata.tsv` | 151,531,112 |

**附带记录一项已知但不采用的选项**：该记录同时提供 ADT（蛋白）矩阵
（`..._ADT-counts.mtx` 55,445,449 字节，features 仅 340 字节，即小型抗体панель）。
若该 panel 含 CD45RA / CD45RO / CD27 / CCR7，理论上比 RNA reference mapping
更适合界定 CM/EM。但 §3 的 primary recipe 已冻结为 Azimuth 参考映射，
**本次不改用 ADT**；此处记录该选项的存在与未采用的理由，
以备日后另行预注册时评估。

---

## 附录 F. Kwok 数据落地与结构审计（2026-09-15，仅 metadata 与基因名，未读取表达值）

**来源**：Zenodo record **7924238**（v1.1）。四个文件经浏览器下载后复制至
`Desktop/Kwok2023_whole_blood/source/`，字节数与 Zenodo API 记录**逐一精确吻合**，
且与下载副本 md5 一致：

| 文件 | 字节 | md5 |
|---|---:|---|
| `..._RNA-counts.mtx` | 3,494,797,339 | `719ef6a602474e78c2331e39be9bee02` |
| `..._RNA-counts_barcodes.tsv` | 9,270,174 | `9812e0d17d1a25c2b65006edfc068e86` |
| `..._RNA-counts_features.tsv` | 177,721 | `f927042500a6b4e9f3e4344b87815626` |
| `..._cell-metadata.tsv` | 151,531,112 | `3c71e9e7947bd8a5b962b4d4a7ba4664` |

**基因空间**：23,013 个基因（单列符号）。冻结的 107 基因命中 **106/107**，
缺 **`GDA`**。按冻结规则基因不增删：`GDA` 照常参与打分，在 Kwok 中对所有细胞贡献同一常数
（与 Gate 1 处理 `CLDN10` 的方式一致），写入 limitation。

**细胞与样本**：272,993 细胞；48 个 `sample_id`；6 个 batch。
`source` 四类：Sepsis 26 样本、**CS（心脏术后）7 样本**、HV 6 样本、Sepsis_conv 9 样本。
sepsis 侧 `diagnosis` 为 CAP 9、Bili 5、Uro 5、IAS 3、CNS 1、Bacteraemia 1、IE 1、NF 1。

**§4b 的事实性澄清**：本次发布的对象中**没有** D1/D3/D5 序列采样列——
每例患者仅一个急性期样本，另加可选的 `*_CONV` 恢复期样本。
因此 §4b 的"每例只取最早急性期样本"在此数据上等价于"排除 `source == Sepsis_conv`"，
无需替代时点规则，attrition 表中该项记为 0 例。**规则本身不变。**

**注释粒度**：`fine_annot` 只有 `Memory_CD4_T_cells` / `Naive_CD4_T_cells` / `Naive_CD8_T_cells`，
**无 CM/EM 划分** —— 与 §3 的判断一致，必须按冻结的 recipe 重新标注。

**组织学差异（须写入 limitation）**：Kwok 是**全血**，`broad_annot` 中
Mature_neutrophils 占 205,527/272,993（75%），CD4_T_cells 仅 10,748。
而 Reyes 与 Gate 1 的 GSE290679 均为 **PBMC**。
这是 cross-cohort transportability 分析中一项真实的可比性限制，不得略过不谈。

## 附录 G. 冻结阈值在 Kwok 上的实际代价（**不改阈值，如实记录**）

`Memory_CD4_T_cells` 逐样本细胞数（作者标注口径，非重新标注后的口径）：

- **Sepsis（26）**：177, 133, 116, 114, 91, 87, 78, 74, 70, 66, 65, 51, 45, 36, 33, 31, 28, 25, 20, 16, 13, 10, 10, 10, 7, 7
- **CS（7）**：65, 63, 37, 34, 22, 22, 17
- HV（6）：307, 200, 186, 181, 143, 131 ——（不参与 H1/H2，仅锚定）

按 §4c 冻结的 **≥30**：Sepsis **16/26** 通过，**CS 仅 4/7 通过**。
预先规定的两个敏感性阈值：≥20 → Sepsis 19/26、CS 6/7；≥50 → Sepsis 12/26、CS 2/7。

**注意**：以上基于作者的 `Memory_CD4_T_cells` 标注；§3 重新标注后按 CM+EM 计的口径
可能略有出入，最终 attrition 以重新标注后的计数为准，**阈值不因该结果调整**。

**后果必须提前说明**：在冻结阈值下，Kwok 臂的危重非脓毒症对照只剩约 4 例。
Gate 2b 本就定位为支持性/探索性、证据等级低于 confirmatory；
此处进一步确认其功效极为有限。该事实在任何结果产生之前记录，
不得在看到结果后用来解释或重新加权结论。

---

# Amendment 1 —— viability gate（2026-09-15，在任何 mapping 之前冻结）

冻结时点：Kwok 数据已落地并完成结构审计（附录 F、G），**但尚未执行任何 reference mapping、
未计算任何 107-gene score、未查看任何 H1/H2 相关量**。

## A1.1 分阶段执行与硬停点
Gate 2b 拆成两段，中间设**硬停点**：

**阶段 1（viability check，本次执行）**
1. 对 Kwok 完成冻结的 Azimuth CM/EM mapping，固化逐细胞标签并计算 checksum；
2. 对 Reyes 完成同样的 mapping，同样固化与 checksum；
3. 按已冻结规则统计每例 `n_CM + n_EffEM`，应用 `≥30`（主）与 `≥20` / `≥50`（敏感性），
   并应用下述 cohort viability rule。

**阶段 1 期间禁止**：计算 107-gene score（UCell/AUCell/ΔUCell）、查看 H1、
查看任何组间 Eff/EM 比例。attrition 计数按 arm 汇总是允许的，因为它是纳入统计而非 outcome。

**阶段 2（H1/H2 全流程）**：仅在通过下述 viability rule 后才执行。

## A1.2 Cohort viability rule（冻结）
以 §4c 的主阈值 ≥30 计算每个 cohort、每个 arm 的合格患者数：

- 某 cohort 任一 arm **< 5** → 该 cohort 判为 **not viable**，不进入 Gate 2b 主分析；
- **两个 cohort 均 not viable**，或仅剩一个 viable cohort（此时 "cross-cohort" 不成立）
  → **当场停止 Gate 2b**，不再执行 UCell/AUCell/H2。
  正文写法固定为：
  > available public cohorts were insufficient for a meaningful cross-cohort test of the
  > prespecified composition hypothesis.
- 两个 cohort 均 ≥ 5/arm → 继续阶段 2，即使人数仍小；
  此时满足 Gate 2b 自身冻结的"可做 supportive cross-cohort"最低门槛，
  但仍为 supportive 而非 confirmatory。

停表判定**只依据 ≥30 主阈值**；≥20 与 ≥50 仅作敏感性报告，
**不得用替代阈值把一个 not viable 的 cohort 救回主分析**。

## A1.3 Mapping 输入集（实现细节，冻结）
Azimuth 参考为 **PBMC**，而 Kwok 为**全血**（75% 中性粒细胞，参考中无对应类别），
若把全部细胞送入映射将产生无意义的强制指派。故冻结：

- **输入集 = 各队列作者标注的 T 细胞全体**，两个队列口径对称：
  - Kwok：`broad_annot ∈ {CD4_T_cells, CD8_T_cells, Cycling_TNK}`；
  - Reyes：`Cell_Type == "T"`（该数据无 CD4/CD8 划分，故送入全部 T 细胞）。
- CD4 TCM / CD4 TEM 的判定**一律来自 Azimuth level-2 标签**，
  不使用作者原有的 CD4/CD8 或 memory/naive 标注来预先筛选。
  作者标注仅用于界定"哪些细胞是 T 细胞"这一步。
- 该输入集规则在看到任何 mapping 结果之前写定，运行后不得调整。

## A1.4 组织学差异的固定表述
Methods / Discussion 必须包含：
> whole-blood versus PBMC preprocessing and cell-recovery differences may affect
> transportability of composition estimates.

## A1.5 signature 可得性
`GDA` 在 Kwok 中不存在，记录为 **106/107 availability**；
不补、不换、不重新训练，固定 signature 不变。

---

# Amendment 2 —— reference-based kNN label transfer（2026-09-15，在任何 mapping 之前冻结）

**定义**：保留同一个外部参考与同一套 level-2 标签体系，
仅把 Seurat/Azimuth 的 anchor-transfer 实现替换为**预先冻结的 reference-based kNN label transfer**。
这是最小、最透明的实现偏离。（选项 B 装 R/Seurat 在本环境技术上不可执行——
出口白名单仅放行 pypi，CRAN 不可达；选项 C 改用 marker 规则会直接改掉已冻结的 primary
annotation recipe，偏离更大。）

## A2.1 冻结的实现规范

| 项 | 冻结内容 |
|---|---|
| Reference | Hao 2021 human PBMC multimodal reference。**须先确认**浏览器下载的 h5ad 与 Azimuth PBMC reference 对应，且含固定的 level-2 标签；目标标签只认 `CD4 TCM` 与 `CD4 TEM` |
| Query input | 按 Amendment 1，Reyes 与 Kwok **对称**地只送入作者原始标注的全部 T cells；**不得**先用本研究 marker 或 107-gene signature 筛 CD4 |
| Features | 只用 reference 中预先定义/计算的 HVG 与 query 的交集；**HVG 由 reference 单独确定**，不得在 query 上重新挑 |
| Normalization | reference 与 query 使用同一套预先指定的 log-normalization；**不做 group-specific normalization** |
| Embedding | **PCA 仅在 reference 上拟合**；query 投影进该 reference PCA 空间 |
| Transfer | 在 reference PCA 空间做 kNN，**k = 30**（写死） |
| Vote | **简单多数投票**（不使用距离加权），自由度最少 |
| Confidence | 记录 top-label vote fraction，**但不得用置信度阈值删除细胞**；低置信度仅作描述性 QC |
| Primary labels | 固定使用 kNN transfer 得到的 level-2 `CD4 TCM` / `CD4 TEM` |
| Sensitivity | §3 已冻结的 canonical-marker rule 保持 sensitivity，**不升级、不改 marker** |
| Concordance | κ / 交叉表仅描述，**不得用于选择采用哪套结果** |

核心原则：**reference-only fit, query-only projection** ——
query 不参与定义特征空间或坐标系，只被投影进去并接受标签。

## A2.2 Methods 必写句（防止被指"写 Azimuth 却没跑 Azimuth"）

> This was not an Azimuth anchor-transfer analysis. We used the PBMC multimodal reference
> and level-2 label ontology employed by Azimuth, but implemented label transfer by
> reference-space k-nearest-neighbor projection because the Seurat/Azimuth software stack
> was unavailable in the locked analysis environment.

## A2.3 失败规则（防止从 A 悄悄滑向 C）

若出现下列任一情形，**停止阶段 1**：
1. 参考 h5ad 中找不到可靠的 level-2 `CD4 TCM` / `CD4 TEM` 标签；
2. reference 与 query 的共享 feature 数低到明显不足以支撑稳定投影。

此时**不得自动切换到 marker-rule primary**。A 实现不了就停下来报告，
任何替代方案须另行预注册。

## A2.4 阶段 1 的禁止事项（重申）
标签锁定之前：不计算任何 107-gene score，不查看 H1，不查看任何组间 Eff/EM 比例。
标签固化后计算 checksum，此后不得重做。

## A2.5 补冻结的实现常数（2026-09-15，在计算任何标签之前）

Amendment 2 未指定的三项，现在写死，运行后不得调整：

1. **HVG 数量 = 2,000**，scanpy `seurat` flavor（按均值分箱后对 dispersion 取箱内 z 分数），
   在 **reference 上单独计算**，query 不参与。
2. **PCA 主成分数 = 50**；标准化为用 reference 的逐基因均值/标准差做 z 分数并在 ±10 处截断，
   均值与标准差**只来自 reference**。
3. **特征交集按 query 分别取**：先在 reference 上定出 top-2000 HVG，
   再分别与 Kwok、Reyes 的基因名取交集；PCA 在 reference 上按该交集重新拟合。
   两个队列因此各有自己的 reference PCA 空间——这是允许的，因为标签迁移在各队列内部完成，
   两队列从不需要共享同一个嵌入；合并只发生在患者层。实际交集基因数逐队列报告。

**归一化（两侧一致）**：一律使用 `raw` 计数，CP10K + log1p。
（reference h5ad 的 `X` 已是作者的归一化结果，为保证 reference 与 query 走同一套预先指定的
归一化，改用 `raw/X` 的原始计数自行归一化；已核实 `raw/X` 为整数计数，`X` 非整数。）

---

# 阶段 1 结果记录 —— Kwok mapping 与 viability（2026-09-15）

**执行**：`Kwok2023_whole_blood/work/` 下 `stageA_ref_stats.py` → `stageB1_build.py`
→ `stageB2_pca.py` → `stageC_map_kwok.py`。全程未计算任何 107-gene score、
未查看 H1、未查看任何组间 Eff/EM 构成比。

## 参考与实现（均按 Amendment 2 / A2.5 冻结值执行）
- Reference：CELLxGENE `4078abd1-063d-4b35-aa04-324dddf8244e.h5ad`
  （collection `b0cf0afa-…`，Hao et al. 2021 *Cell*，"nygc multimodal pbmc"），
  161,764 细胞 × 20,264 基因；`obs["celltype.l2"]` 31 类，含 `CD4 TCM` 与 `CD4 TEM`
  → **A2.3 失败规则 1 不触发**。
- 归一化：reference 与 query 一律用 `raw` 计数、CP10K + log1p（已核实 `raw/X` 为整数计数）。
- HVG：reference 单独计算 top-2000（seurat flavor）；**与 Kwok 基因交集 1,558 个**
  → **A2.3 失败规则 2 不触发**（下限 500）。
- PCA：仅在 reference 上拟合，50 主成分，`random_state=20260915`，解释方差合计 0.359；
  query 只做投影。
- kNN：k=30，简单多数投票。
- Query 输入：Kwok 作者标注 `broad_annot ∈ {CD4_T_cells, CD8_T_cells, Cycling_TNK}`，
  共 **19,091** 个细胞（未用作者的 CD4/memory 标注预筛）。
- QC（描述性，未据此删除细胞）：vote_fraction 中位数 **0.90**，< 0.5 者占 **3.7%**。
- **逐细胞标签已固化**：`work/kwok_azimuth_l2_labels.csv`，
  md5 **`9a3fc5292432984e6bd30d938a520d65`**。此后不得重做。

## 每例 n_CM + n_TEM（纳入统计量；未拆分 CM/TEM，未看构成比）

- **Sepsis（26）**：180, 150, 115, 112, 89, 88, 84, 81, 80, 76, 72, 69, 54, 41, 29, 28, 28, 23, 20, 15, 14, 12, 12, 12, 11, 8
- **CS（7）**：76, 64, 36, 35, 26, 22, 16
- HV（6，不参与检验）：291, 207, 204, 202, 148, 144

按 §4c 主阈值 **≥30**：Sepsis **14/26**，**CS 4/7**。
敏感性阈值：≥20 → Sepsis 19/26、CS 6/7；≥50 → Sepsis 13/26、CS 2/7。

与作者原始 `Memory_CD4_T_cells` 标注的逐例计数高度一致
（CS 侧作者口径为 65, 63, 37, 34, 22, 22, 17，同样 4/7 通过 ≥30），
说明该判定不是标注方法造成的假象。

## Viability 判定（按 Amendment 1 A1.2，只依据 ≥30 主阈值）

**Kwok 的 critical non-sepsis（CS）臂 = 4 < 5 → Kwok 判为 not viable。**

因此无论 Reyes 结果如何，最多只剩一个 viable cohort，
"cross-cohort" 不成立 → **触发硬停点，Gate 2b 就此停止**。
不再执行 UCell / AUCell / H2，不再计算任何 107-gene score。

**正文写法（A1.2 预先固定）**：
> available public cohorts were insufficient for a meaningful cross-cohort test of the
> prespecified composition hypothesis.

**不得**改用 ≥20 把 CS 臂救回主分析——A1.2 已明确停表判定只依据 ≥30。

---

# 阶段 1 结果记录（续）—— Reyes mapping 与两队列完整 attrition（2026-09-16）

**执行**：与 Kwok 完全相同的冻结管线与参数。
Reference 侧针对 Reyes 基因空间重建：reference-only top-2000 HVG → **与 Reyes 基因交集 1,616 个特征**
（> 500 下限，A2.3 失败规则 2 不触发）；PCA 仅在 reference 上拟合，50 主成分，
`random_state=20260915`，解释方差 0.364；query 只做投影；kNN k=30 简单多数投票。
Query 输入：作者标注 `Cell_Type == "T"` 的全部 T 细胞 **32,341** 个
（未用作者标注预筛 CD4；该数据本就无 CD4/CD8 划分）。
归一化：`scp_gex_matrix_raw.csv.gz`（md5 `cd883fc026c01c694e766ed716cad233`，
113,908,822 字节）的原始计数，CP10K + log1p，库大小按每个细胞跨全部 22,858 个基因求和
（中位 1,766 UMI）。

- **逐细胞标签已固化**：`Reyes2020_SCP548/work/reyes_azimuth_l2_labels.csv`，
  md5 **`b84708ec4f584aca1ef49ec8ecfc4fbd`**。
- QC（描述性，未据此删除细胞）：vote_fraction 中位数 **0.57**，< 0.5 者占 **37.0%**。
  **明显低于 Kwok（0.90 / 3.7%）**，与 Reyes 细胞测序深度浅（中位 1,766 UMI）一致。
  该差异须在 limitation 中如实说明：浅文库使参考映射的置信度下降。

## 每例 n_CM + n_TEM（未拆分 CM/TEM，未看构成比）
- **ICU-SEP（8）**：186, 70, 69, 24, 24, 20, 16, 5
- **ICU-NoSEP（7）**：71, 68, 60, 54, 44, 33, 32

按主阈值 **≥30**：**ICU-SEP 3/8**，ICU-NoSEP 7/7。
敏感性：≥20 → 6/8、7/7；≥50 → 3/8、4/7。

**Reyes 的 sepsis 臂 = 3 < 5 → Reyes 判为 not viable。**

## Stage-1 完整 attrition（两队列）

| Cohort | Sepsis ≥30 | Critical non-sepsis ≥30 | Viability |
|---|---:|---:|---|
| Kwok | 14/26 | 4/7 | **Not viable**（对照臂 4 < 5） |
| Reyes | 3/8 | 7/7 | **Not viable**（脓毒症臂 3 < 5） |
| Gate 2b | — | — | **Stop** |

**两个 cohort 均 not viable** —— 这是 A1.2 中更强的那一分支。
前文"Kwok 触发硬停点"的记录原文保留不动；本节只是把预注册的 Stage 1 做完整。

**Stage-1 completion note.** Although the Kwok cohort alone was sufficient to establish that
a two-cohort analysis could no longer proceed, Reyes mapping was completed solely to fulfill
the prespecified Stage-1 attrition accounting. No CM/TEM composition comparison, 107-gene
score, UCell/AUCell analysis, or H1/H2 testing was performed.

# Gate 2b 正式关闭（2026-09-16）

正文结论句（A1.2 预先固定）：
> available public cohorts were insufficient for a meaningful cross-cohort test of the
> prespecified composition hypothesis.

Gate 2 保持 DRAFT v0.9 / PARKED，五条资格标准未作任何放宽。
confirmatory Gate 2 的路径仍为前瞻性自建队列。

## 正文措辞约束 —— Reyes QC（2026-09-16，仅措辞，不涉及任何分析）

limitation 中的固定写法：

> Reference-mapping confidence was substantially lower in Reyes than in Kwok
> (median vote fraction 0.57 vs 0.90; 37.0% vs 3.7% below 0.5), consistent with the
> shallow sequencing depth of the Reyes dataset.

**不得**写成"浅文库导致 mapping 置信度下降"这类确定因果的表述；
更重要的是，**该 QC 观察不得被表述为 Gate 2b 关闭的理由**。
Gate 2b 的正式关闭理由只有一条，即预注册的 viability rule：
在 ≥30 主阈值下两个 cohort 各有一臂 < 5。
QC 数字只作为额外的测量/可迁移性限制呈现，避免被读作把 post hoc QC 混入停表判定。
