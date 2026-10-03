#!/usr/bin/env python3
"""
Gate 1 · Step 1b —— raw 命名空间审计（SAP 第 7/8 节触发）
processed /X 仅 5,000 HVG，coverage 31.8% < 90%，按 SAP 回退检查 /raw。

本步骤仍然是盲的：不读取任何 obs 分组标签，不按组切分，不计算任何分组差异。
只做两件事：
  1) 用 /raw/var 的基因名算 107 个冻结 ortholog 的覆盖率（纯元数据）
  2) 抽取 /raw/X 的前若干个非零值，判断它是整数 counts 还是已标准化
     —— 这是 SAP 第 8 节要求的 expression source 判定，与分组无关

用法: python3 gate1_01b_raw_audit.py /path/to/file.h5ad
"""
import sys, os, json
import numpy as np, pandas as pd, h5py

SIGFILE = os.path.expanduser("~/Desktop/Context-dependent_ICOS_analysis_results/metadata/rat_107gene_frozen_signature.csv")
RESULTS = os.path.expanduser("~/Desktop/Context-dependent_ICOS_analysis_results/analysis_results")

def names_from_dataframe(g):
    """从 h5ad 的 dataframe group 里取 _index，以及所有可能的符号列名。"""
    out = {}
    idx = g.attrs.get("_index", b"_index")
    idx = idx.decode() if isinstance(idx, bytes) else str(idx)
    if idx in g:
        out["_index"] = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in g[idx][:]])
    for k in g.keys():
        if k == idx: continue
        try:
            d = g[k]
            if isinstance(d, h5py.Dataset) and d.dtype.kind in "SOU" and d.shape[0] == out["_index"].shape[0]:
                out[k] = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in d[:]])
        except Exception:
            pass
    return out

def main(path):
    os.makedirs(RESULTS, exist_ok=True)
    rep = {}
    sig = pd.read_csv(SIGFILE)
    H = {g for g in sig["human_symbol"].dropna().astype(str) if g}
    print(f"冻结 signature H = {len(H)} 个\n")

    with h5py.File(path, "r") as f:
        if "raw" not in f: sys.exit("该文件没有 /raw")
        print("=== /raw/var 的命名空间 ===")
        spaces = names_from_dataframe(f["raw"]["var"])
        cov = {}
        for name, arr in spaces.items():
            hit = H & set(arr)
            cov[name] = dict(n_genes=int(arr.shape[0]), n_hit=len(hit),
                             coverage=round(len(hit)/len(H), 4), missing=sorted(H - set(arr)))
            print(f"  raw/var['{name}']  基因数={arr.shape[0]:>6d}  命中 {len(hit):>3d}/{len(H)}"
                  f"  coverage = {100*len(hit)/len(H):.1f}%")
        rep["raw_coverage"] = cov
        best = max(cov, key=lambda k: cov[k]["coverage"])
        print(f"\n  最佳：raw/var['{best}']  ({100*cov[best]['coverage']:.1f}%)")
        if cov[best]["missing"]:
            print(f"  仍缺失：{', '.join(cov[best]['missing'])}")
        else:
            print("  全部 107 个基因均可用。")

        print("\n=== /raw/X 是否为整数 counts（盲检，不涉及任何分组）===")
        X = f["raw"]["X"]
        n = min(200000, X["data"].shape[0])
        vals = X["data"][:n]
        is_int = np.allclose(vals, np.round(vals))
        print(f"  抽取前 {n:,} 个非零值")
        print(f"  最小 {vals.min():.4f}  最大 {vals.max():.4f}  均值 {vals.mean():.4f}")
        print(f"  全部为整数：{is_int}")
        print(f"  取值样例：{np.unique(vals[:2000])[:12]}")
        rep["raw_X_is_integer_counts"] = bool(is_int)
        rep["raw_X_sample_stats"] = dict(min=float(vals.min()), max=float(vals.max()), mean=float(vals.mean()))
        if is_int:
            print("\n  → /raw/X 是原始 counts。UCell 打分与 pseudobulk 均应使用该矩阵，")
            print("    UCell 前在细胞内做 CPM/log 标准化（UCell 基于秩，标准化不改变秩，仅为可读性）。")
        else:
            print("\n  → /raw/X 非整数，可能已标准化；需在 SAP 第 8 节记录并说明如何使用。")

    json.dump(rep, open(f"{RESULTS}/gate1_raw_audit.json", "w"), ensure_ascii=False, indent=2)
    print(f"\n已写出 {RESULTS}/gate1_raw_audit.json")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else sys.exit("需要 h5ad 路径"))
