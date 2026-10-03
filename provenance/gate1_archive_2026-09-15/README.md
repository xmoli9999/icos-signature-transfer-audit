# Gate 1 archive — a lymph-derived ICOS-associated CD4 T-cell program in human sepsis

Deposit date: 2026-09-15
Contact: Yan Li, ORCID 0009-0000-2236-5199, criticalcare@163.com
Department of Emergency Medicine, Shanghai Fourth People's Hospital,
School of Medicine, Tongji University, Shanghai 200434, China.

## Timestamp disclaimer (read first)

> **Timestamp obtained after Gate 1 unblinding; used to lock the completed Gate 1 analysis
> and all subsequent hypotheses, not as evidence of prospective registration of Gate 1.**

Gate 1 的分析计划（`metadata/gate1_SAP.md`）及其 Amendments 1–4 均在各自所用分析之前
以书面形式冻结，Amendment 3 与 Amendment 4 明确在任何 group comparison 之前完成；
但本存档本身是在揭盲之后建立的，因此**不构成 Gate 1 的前瞻性注册证明**。
它锁定的是：已完成的 Gate 1 分析全文、以及由其生成的后续假设（`gate2_preregistration_DRAFT.md`）。
Gate 2 的时间戳则是真正前瞻性的——该文件在接触任何 Gate 2 候选队列数据之前建立。

## 结果一句话

The prespecified cross-species validation was negative (β = −0.00138,
95% CI −0.00531 to 0.00255; P = 0.493); sensitivity analyses were concordantly null.

## 内容

- `metadata/gate1_SAP.md` — 权威分析计划（FROZEN v1.0）+ Amendments 1–4 + 盲态 QC 记录
  + Gate 1 结果记录 + 正式结论用词。
- `metadata/gate2_preregistration_DRAFT.md` — 由 Gate 1 探索性分解生成的新假设，
  在接触任何候选队列数据之前写定。
- `metadata/rat_107gene_frozen_signature.csv` — 冻结的 107 基因（大鼠→人直系同源）。
  **全程未增删、未重新加权、未重新选择。**
- `scripts/` — 从结构审计到揭盲检验的全部脚本，按 01 → 04 顺序即可复现。
- `analysis_results/` — 审计输出、校准诊断、100×107 对照基因表、
  盲态与揭盲的患者层表、主图、运行日志。
- `provenance/` — 上游大鼠端的冻结计划与 signature 推导产物。

## 复现方式

输入数据不在本存档内（体积原因），从 GEO 获取：
- 人类队列：**GSE290679**，`*_allsamples_processed_labeled_filtered.h5ad`
- 大鼠上游：**GSE285325**

然后依次运行（`ICOS_BASE` 指向本目录的上级项目目录）：

```
python3 scripts/gate1_01_audit_v2.py  <h5ad>     # 结构审计
python3 scripts/gate1_02_score.py     <h5ad>     # UCell / AUCell
ICOS_STAGE=1 python3 scripts/gate1_03_calibrate.py <h5ad>   # 对照集（seed 20260915）
python3 scripts/gate1_03_calibrate.py <h5ad>     # ΔUCell
python3 scripts/gate1_04_gate_test.py            # 揭盲检验
```

随机性只有一处：`gate1_03_calibrate.py` 的 `SEED = 20260915`。
对照集已随存档提供（`gate1_control_sets.csv.gz`），重跑应逐位一致；
脚本在检测到已存在的对照集时会直接载入，以保证与首次运行完全相同。

自检基线：signature UCell 与独立实现的最大绝对差 7.4e-09；
`MANIFEST.md` 给出全部文件的 MD5。
