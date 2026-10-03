#!/usr/bin/env python3
"""
Gate 1 · Step 1 —— 只读元数据审计（16 GB 机器安全）

安全设计：
  * 结构信息（X 的编码、形状、nnz、layers、raw、obsm）全部用 h5py 读 HDF5 属性，
    不经过 anndata，不会触发任何矩阵加载。
  * 仅 obs / var 用 anndata 的 backed 模式读取（为了让 anndata 解码 categorical）。
    backed 模式下 .X 保持在磁盘上，本脚本从不引用 a.X / a.raw.X / a.layers[...]。
  * 全程没有 .copy()、没有 .to_memory()、没有对矩阵的切片。
  * 唯一进入内存的是 obs 表与 var 名称；脚本会先报告 obs 的预估内存占用。

用法：
  python3 gate1_01_audit.py /path/to/GSE290679_processed_labeled_filtered.h5ad
"""
import sys, os, json
import numpy as np
import pandas as pd

RESULTS = os.path.expanduser("~/Desktop/Context-dependent_ICOS_analysis_results/analysis_results")
SIGFILE = os.path.expanduser("~/Desktop/Context-dependent_ICOS_analysis_results/metadata/rat_107gene_frozen_signature.csv")

def main(path):
    if not os.path.exists(path):
        sys.exit(f"找不到文件：{path}")
    os.makedirs(RESULTS, exist_ok=True)
    print(f"文件：{path}")
    print(f"磁盘大小：{os.path.getsize(path)/1e9:.2f} GB\n")

    # ---------- 第一部分：纯 h5py 结构审计，完全不碰 anndata ----------
    import h5py
    struct = {}
    with h5py.File(path, "r") as f:
        print("=== HDF5 顶层结构 ===")
        for k in f.keys():
            obj = f[k]
            kind = "group" if isinstance(obj, h5py.Group) else "dataset"
            enc  = obj.attrs.get("encoding-type", b"")
            enc  = enc.decode() if isinstance(enc, bytes) else str(enc)
            print(f"  /{k:<12s} {kind:<8s} encoding-type={enc}")
        struct["top_level_keys"] = list(f.keys())

        if "X" in f:
            X = f["X"]
            enc = X.attrs.get("encoding-type", b"")
            enc = enc.decode() if isinstance(enc, bytes) else str(enc)
            shp = X.attrs.get("shape", None)
            struct["X_encoding"] = enc
            if isinstance(X, h5py.Group):           # 稀疏
                nnz = X["data"].shape[0]
                dt  = str(X["data"].dtype)
                struct.update(X_shape=list(map(int, shp)) if shp is not None else None,
                              X_nnz=int(nnz), X_dtype=dt)
                print(f"\n  X: {enc}  shape={shp}  nnz={nnz:,}  dtype={dt}")
                print(f"     → 全量载入约需 {nnz*(np.dtype(dt).itemsize+4)/1e9:.1f} GB（不要这么做）")
            else:                                    # 稠密
                struct.update(X_shape=list(X.shape), X_dtype=str(X.dtype))
                print(f"\n  X: dense  shape={X.shape}  dtype={X.dtype}")
                print(f"     → 全量载入约需 {np.prod(X.shape)*X.dtype.itemsize/1e9:.1f} GB（不要这么做）")

        struct["has_raw"]  = "raw" in f
        struct["layers"]   = list(f["layers"].keys()) if "layers" in f else []
        struct["obsm"]     = list(f["obsm"].keys())   if "obsm"   in f else []
        print(f"\n  raw 存在：{struct['has_raw']}")
        print(f"  layers：{struct['layers']}")
        print(f"  obsm：{struct['obsm']}")
        if struct["layers"]:
            print("  注意：存在 layers，后续任何读取都要显式指定用哪一层，不能依赖默认值。")

    # ---------- 第二部分：obs / var（backed 模式，矩阵留在磁盘） ----------
    import anndata as ad
    print("\n=== 以 backed='r' 打开（X 保持在磁盘上）===")
    a = ad.read_h5ad(path, backed="r")
    print(f"shape = {a.shape}")
    mem = a.obs.memory_usage(deep=True).sum()/1e6
    print(f"obs 表进入内存约 {mem:.0f} MB（这是本脚本唯一的内存开销）")

    print("\n=== obs 列 ===")
    for c in a.obs.columns:
        s = a.obs[c]
        nu = s.nunique(dropna=True)
        ex = list(pd.Series(s.dropna().unique()).astype(str)[:8])
        print(f"  {c:<34s} n_unique={nu:<8d} 例: {ex}")

    def guess(keys):
        return [c for c in a.obs.columns if any(k in c.lower() for k in keys)]
    cand = {
        "sample/patient": guess(["sample","patient","donor","subject","orig.ident","orig_ident","id"]),
        "group/condition": guess(["group","condition","disease","diagnos","status","cohort","sepsis","icu"]),
        "celltype":       guess(["cell","type","anno","cluster","label","ident","subset","lineage"]),
    }
    print("\n=== 候选字段 ===")
    for k, v in cand.items():
        print(f"  {k}: {v}")

    print("\n=== 候选分组字段的取值分布 ===")
    for c in cand["group/condition"][:6]:
        print(f"\n  [{c}]")
        print(a.obs[c].value_counts(dropna=False).head(15).to_string())

    print("\n=== 候选细胞类型字段的取值 ===")
    for c in cand["celltype"][:4]:
        vc = a.obs[c].value_counts(dropna=False)
        print(f"\n  [{c}]  共 {len(vc)} 类")
        print(vc.head(50).to_string())

    print("\n=== 样本 × 分组 交叉表（用于确认患者层单位）===")
    sid = cand["sample/patient"][0] if cand["sample/patient"] else None
    gid = cand["group/condition"][0] if cand["group/condition"] else None
    if sid and gid:
        ct = pd.crosstab(a.obs[sid], a.obs[gid])
        print(f"  {sid} × {gid}: {ct.shape[0]} 个样本 × {ct.shape[1]} 组")
        print(f"  每组样本数: {(ct>0).sum(axis=0).to_dict()}")

    # ---------- 第三部分：冻结 signature 的覆盖率 ----------
    print("\n=== 冻结的 107 基因在 var_names 中的命中 ===")
    if os.path.exists(SIGFILE):
        sig = pd.read_csv(SIGFILE)
        col = [c for c in sig.columns if "human" in c.lower()] or [sig.columns[-1]]
        hs  = {g for g in sig[col[0]].dropna().astype(str) if g}
        vn  = set(pd.Index(a.var_names).astype(str))
        hit, miss = sorted(hs & vn), sorted(hs - vn)
        print(f"  signature {len(hs)} 个 → 命中 {len(hit)}，缺失 {len(miss)}")
        if miss:
            print(f"  缺失：{', '.join(miss)}")
            print("  （若缺失较多，需在打分前决定是否按可用基因重标定，此决定要先写进 SAP）")
    else:
        hit, miss, hs = [], [], set()
        print(f"  未找到 {SIGFILE}")

    # ---------- 输出 ----------
    pd.DataFrame({
        "obs_column": list(a.obs.columns),
        "n_unique":   [a.obs[c].nunique(dropna=True) for c in a.obs.columns],
        "dtype":      [str(a.obs[c].dtype) for c in a.obs.columns],
    }).to_csv(f"{RESULTS}/gate1_obs_columns.csv", index=False)

    struct.update(shape=list(a.shape), obs_mem_MB=round(mem,1),
                  sig_total=len(hs), sig_hit=len(hit), sig_missing=miss,
                  candidates=cand)
    json.dump(struct, open(f"{RESULTS}/gate1_audit_summary.json","w"),
              ensure_ascii=False, indent=2)

    a.file.close()
    print(f"\n已写出：\n  {RESULTS}/gate1_obs_columns.csv\n  {RESULTS}/gate1_audit_summary.json")
    print("审计完成 —— 全程未加载任何表达矩阵。")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "")
