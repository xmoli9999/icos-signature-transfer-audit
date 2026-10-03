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
- H1 主模型：`logit(Eff/EM fraction) ~ group + z(median nFeature)`，OLS + HC3 95% CI；
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

## 6. 判定规则

- H1 方向相反（Eff/EM 在 sepsis 中升高）→ 本条主线停止，不得改写假设方向。
- H1 成立但 H2 亦显示 within-state 下降 → 结论降级为
  "composition shift 与 within-state suppression 并存"，不得只报前者。
- H1 与 H2 均 null → 如实报告为第二次 negative validation，项目转为方法学/阴性结果报告。

## 7. 明确禁止

- 不再对 GSE290679 追加任何以寻找阳性结果为目的的分析。
- 不得用 GSE279451/279452 充当本假设的独立验证队列。
- 不得在看过新队列的结果后回头修改 §1–§6 的任何一条。
