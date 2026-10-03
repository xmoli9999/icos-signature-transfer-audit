#!/usr/bin/env python3
"""
Gate 1 · Step 1 —— 只读元数据审计 v2（审查修正版）

v2 相对 v1 的修正（均来自代码审查）：
  [必改1] human ortholog 列不再自动猜：只接受明确的列名，找不到就停止覆盖率审计并报错，
          绝不 fallback 到最后一列。
  [必改2] sample/patient 候选删除裸 "id"，改用显式名单；crosstab 前检查基数，
          并检查"一个 sample 是否只对应一个 group"。
  [补3]   同时报告 var.columns、var_names 样例与候选 gene-symbol 列，
          并分别给出各命名空间下的覆盖率；不自动切换命名空间。
  [补4]   /X、每个 /layers/*、/raw/X 分别报告 encoding / shape / dtype / nnz。
  [小]    稀疏矩阵内存估算计入 indices 的真实 dtype 与 indptr。
  [小]    不再宣称 obs 是"唯一"内存开销（unique/value_counts/crosstab 都有临时开销）。

安全性不变：结构信息走 h5py 属性；仅 obs/var 用 anndata backed 解码；
全文无 .X / .raw.X / .layers[...] 取值、无 .copy()、无 .to_memory()、无矩阵切片。
"""
import sys, os, json
import numpy as np
import pandas as pd

RESULTS = os.path.expanduser("~/Desktop/Context-dependent_ICOS_analysis_results/analysis_results")
SIGFILE = os.path.expanduser("~/Desktop/Context-dependent_ICOS_analysis_results/metadata/rat_107gene_frozen_signature.csv")

HUMAN_COL_NAMES = {"human_symbol","human symbol","human gene name","human_gene_name",
                   "human_gene_symbol","human","hsapiens_symbol"}
SAMPLE_COL_NAMES = ["sample_id","sample","patient_id","patient","donor_id","donor",
                    "subject_id","subject","orig.ident","orig_ident","library","specimen"]
GROUP_KEYS = ["group","condition","disease","diagnos","status","cohort","sepsis","icu","class"]
CELLTYPE_KEYS = ["celltype","cell_type","cell.type","annotation","anno","cluster",
                 "label","ident","subset","lineage","majority"]

def _attr(obj, key, default=""):
    v = obj.attrs.get(key, default)
    return v.decode() if isinstance(v, bytes) else str(v)

def describe_matrix(f, path):
    """只读 HDF5 属性与 dataset 形状/dtype，绝不读取数值。"""
    if path not in f: return None
    o = f[path]
    enc = _attr(o, "encoding-type")
    if isinstance(o, __import__("h5py").Group):        # 稀疏
        shp = o.attrs.get("shape", None)
        d_dt = o["data"].dtype; i_dt = o["indices"].dtype; p_dt = o["indptr"].dtype
        nnz = int(o["data"].shape[0]); nptr = int(o["indptr"].shape[0])
        est = (nnz*d_dt.itemsize + nnz*i_dt.itemsize + nptr*p_dt.itemsize)/1e9
        return dict(path=path, encoding=enc, shape=[int(x) for x in shp] if shp is not None else None,
                    dtype=str(d_dt), indices_dtype=str(i_dt), nnz=nnz, est_full_load_GB=round(est,2))
    else:                                              # 稠密
        est = float(np.prod(o.shape))*o.dtype.itemsize/1e9
        return dict(path=path, encoding=enc or "dense", shape=[int(x) for x in o.shape],
                    dtype=str(o.dtype), nnz=None, est_full_load_GB=round(est,2))

def main(path):
    if not os.path.exists(path):
        sys.exit(f"找不到文件：{path!r}\n用法：python3 gate1_01_audit_v2.py /path/to/file.h5ad")
    os.makedirs(RESULTS, exist_ok=True)
    print(f"文件：{path}")
    print(f"磁盘大小：{os.path.getsize(path)/1e9:.2f} GB\n")
    report = {"file": path, "file_GB": round(os.path.getsize(path)/1e9,2)}

    # ---------- 1. h5py 结构审计（不经过 anndata）----------
    import h5py
    with h5py.File(path, "r") as f:
        print("=== HDF5 顶层结构 ===")
        for k in f.keys():
            o = f[k]
            kind = "group" if isinstance(o, h5py.Group) else "dataset"
            print(f"  /{k:<12s} {kind:<8s} encoding-type={_attr(o,'encoding-type')}")
        report["top_level_keys"] = list(f.keys())

        mats = []
        m = describe_matrix(f, "X")
        if m: mats.append(m)
        if "layers" in f:
            for ln in f["layers"].keys():
                mm = describe_matrix(f, f"layers/{ln}")
                if mm: mats.append(mm)
        if "raw" in f and "X" in f["raw"]:
            mm = describe_matrix(f, "raw/X")
            if mm: mats.append(mm)
        report["matrices"] = mats
        print("\n=== 各矩阵（仅属性，未读取任何数值）===")
        for mm in mats:
            nnz = f"{mm['nnz']:,}" if mm['nnz'] is not None else "dense"
            print(f"  {mm['path']:<14s} enc={mm['encoding']:<18s} shape={mm['shape']} "
                  f"dtype={mm['dtype']:<8s} nnz={nnz:>14s}  全量载入约 {mm['est_full_load_GB']} GB")
        print("  ↑ 第二步必须显式指定从哪一个矩阵读取；processed 的 /X 可能已被 scale/regress，")
        print("    直接拿它算 UCell 是错的。若 /raw/X 存在且为整数 counts，pseudobulk 应优先用它。")

        report["raw_var_present"] = ("raw" in f and "var" in f["raw"])
        if report["raw_var_present"]:
            print("  /raw/var 存在 —— 若 processed 命名空间覆盖率不足，可回退到 raw 命名空间。")

    # ---------- 2. obs / var（backed 模式）----------
    import anndata as ad
    print("\n=== 以 backed='r' 打开（X 留在磁盘）===")
    a = ad.read_h5ad(path, backed="r")
    n_cells, n_genes = a.shape
    print(f"shape = {a.shape}")
    mem = a.obs.memory_usage(deep=True).sum()/1e6
    print(f"obs 表约 {mem:.0f} MB 进入内存（后续 unique/value_counts/crosstab 另有少量临时开销）")
    report.update(shape=[int(n_cells), int(n_genes)], obs_mem_MB=round(mem,1))

    print("\n=== obs 列 ===")
    obs_info=[]
    for c in a.obs.columns:
        s = a.obs[c]; nu = int(s.nunique(dropna=True))
        ex = list(pd.Series(s.dropna().unique()).astype(str)[:6])
        obs_info.append(dict(column=c, n_unique=nu, dtype=str(s.dtype)))
        print(f"  {c:<34s} n_unique={nu:<8d} {str(s.dtype):<12s} 例: {ex}")

    low = {c for c in a.obs.columns if a.obs[c].nunique(dropna=True) <= 500}
    sample_cand = [c for c in a.obs.columns
                   if c.lower() in SAMPLE_COL_NAMES or
                      any(c.lower().endswith(k) or c.lower()==k for k in SAMPLE_COL_NAMES)]
    sample_cand = [c for c in sample_cand if c in low]          # 排除 barcode 级高基数列
    group_cand  = [c for c in a.obs.columns if any(k in c.lower() for k in GROUP_KEYS) and c in low]
    ct_cand     = [c for c in a.obs.columns if any(k in c.lower() for k in CELLTYPE_KEYS) and c in low]
    print("\n=== 候选字段（已按基数 ≤500 过滤，裸 'id' 不参与匹配）===")
    print(f"  sample/patient: {sample_cand}")
    print(f"  group/condition: {group_cand}")
    print(f"  celltype: {ct_cand}")
    report["candidates"] = dict(sample=sample_cand, group=group_cand, celltype=ct_cand)

    print("\n=== 候选分组字段的完整取值 ===")
    for c in group_cand[:8]:
        print(f"\n  [{c}]"); print(a.obs[c].value_counts(dropna=False).to_string())

    print("\n=== 候选细胞类型字段的完整取值 ===")
    for c in ct_cand[:5]:
        vc = a.obs[c].value_counts(dropna=False)
        print(f"\n  [{c}]  共 {len(vc)} 类"); print(vc.to_string())

    print("\n=== sample × group 一致性检查 ===")
    for sid in sample_cand[:3]:
        ns = int(a.obs[sid].nunique(dropna=True))
        if ns > 500:
            print(f"  [{sid}] 基数 {ns} 过高，跳过"); continue
        for gid in group_cand[:3]:
            g = a.obs.groupby(sid, observed=True)[gid].nunique()
            bad = int((g > 1).sum())
            ct = pd.crosstab(a.obs[sid], a.obs[gid])
            print(f"  {sid} × {gid}: {ns} 个 sample，每组样本数 {(ct>0).sum(axis=0).to_dict()}"
                  f"；跨组样本 {bad} 个" + ("  ← 异常，需排查" if bad else ""))

    # ---------- 3. 冻结 signature 的覆盖率（多命名空间并报，不自动切换）----------
    print("\n=== 冻结 signature 的覆盖率 ===")
    print(f"  var_names 前 5 个: {list(pd.Index(a.var_names).astype(str)[:5])}")
    print(f"  var.columns: {list(a.var.columns)}")
    sym_cols = [c for c in a.var.columns
                if any(k in c.lower() for k in ["symbol","gene_name","gene name","genes","feature_name"])]
    print(f"  候选 gene-symbol 列: {sym_cols}（不自动切换，由人决定）")

    cov = {}
    if not os.path.exists(SIGFILE):
        print(f"  !! 未找到 {SIGFILE} —— 覆盖率审计跳过")
    else:
        sig = pd.read_csv(SIGFILE)
        hits = [c for c in sig.columns if c.strip().lower() in HUMAN_COL_NAMES]
        if len(hits) != 1:
            print(f"  !! 无法在 signature 文件中明确识别 human gene 列。")
            print(f"     现有列：{list(sig.columns)}")
            print(f"     可接受的列名：{sorted(HUMAN_COL_NAMES)}")
            print(f"     覆盖率审计已停止（拒绝猜测；请重命名该列后重跑）。")
            report["coverage_error"] = "human gene column not unambiguously identified"
        else:
            hcol = hits[0]
            H = {g for g in sig[hcol].dropna().astype(str) if g and g.lower() != "nan"}
            print(f"  冻结的 human ortholog 集合 H = {len(H)} 个（来自列 '{hcol}'）")
            spaces = {"var_names": set(pd.Index(a.var_names).astype(str))}
            for c in sym_cols:
                spaces[f"var['{c}']"] = set(a.var[c].dropna().astype(str))
            for name, universe in spaces.items():
                hit = H & universe
                cov[name] = dict(n_hit=len(hit), coverage=round(len(hit)/len(H),4),
                                 missing=sorted(H-universe))
                print(f"  {name:<28s} 命中 {len(hit):>3d}/{len(H)}  coverage = {100*len(hit)/len(H):.1f}%")
            best = max(cov, key=lambda k: cov[k]["coverage"])
            print(f"\n  最佳命名空间：{best}（{100*cov[best]['coverage']:.1f}%）")
            if cov[best]["coverage"] < 0.90:
                print("  ← 低于 SAP 预设的 90% 阈值：应先检查 /raw/var 能否恢复，")
                print("    若 raw 亦不足，则保留可测子集并完整报告 attrition，不得补基因或重训权重。")
            if cov[best]["missing"]:
                print(f"  缺失（相对最佳命名空间）：{', '.join(cov[best]['missing'][:40])}")
        report["coverage"] = cov

    pd.DataFrame(obs_info).to_csv(f"{RESULTS}/gate1_obs_columns.csv", index=False)
    json.dump(report, open(f"{RESULTS}/gate1_audit_summary.json","w"), ensure_ascii=False, indent=2)
    a.file.close()
    print(f"\n已写出：\n  {RESULTS}/gate1_obs_columns.csv\n  {RESULTS}/gate1_audit_summary.json")
    print("审计完成 —— 全程未读取任何表达数值。")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "")
