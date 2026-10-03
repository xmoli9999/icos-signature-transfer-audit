#!/usr/bin/env python3
"""
Gate 1 · Step 2 —— 逐细胞打分（UCell + AUCell），16 GB 安全

冻结设定（见 gate1_SAP.md 与 Amendment 1/2）：
  表达源   = /raw/X（log1p(counts × size factor)，秩不变，UCell/AUCell 适用）
  命名空间 = /raw/var['_index']，107/107 覆盖
  细胞集   = celltype_coarse ∈ {CD4 T, CD4 Treg} 的并集（235,994），
             四个分析群体都是它的子集，一次读取算完
  UCell    : maxRank = 1500
  AUCell   : top 5% 基因为阈值（36,601 × 0.05 ≈ 1830）

内存策略：按全局行块（默认 10,000 行）流式读取 /raw/X 的 data/indices，
逐细胞算秩后即丢弃；峰值 < 1 GB。全程不构造完整矩阵、不做 copy。

本脚本【不做任何分组比较】，只输出逐细胞分数。分组检验在 Step 3。

用法: python3 gate1_02_score.py /path/to/file.h5ad [块大小]
"""
import sys, os, time
import numpy as np, pandas as pd, h5py

B = os.path.expanduser("~/Desktop/Context-dependent_ICOS_analysis_results")
SIGFILE = f"{B}/metadata/rat_107gene_frozen_signature.csv"
MAPFILE = f"{B}/metadata/memory_cd4_label_mapping.csv"
OUT     = f"{B}/analysis_results"
MAXRANK = 1500
AUC_TOP_FRAC = 0.05

def main(path, block=10000):
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()

    sig = pd.read_csv(SIGFILE)
    H = sorted({g for g in sig["human_symbol"].dropna().astype(str) if g})
    mp = pd.read_csv(MAPFILE)
    print(f"signature {len(H)} 个；mapping {len(mp)} 行\n")

    import anndata as ad
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs[["orig.ident","SCARAB_ID","Group","Batch","Age","Biologic_Sex",
                 "ICU","Source","celltype_coarse","celltype_fine",
                 "nCount_RNA","nFeature_RNA","percent.mt"]].copy()
    a.file.close(); del a
    print(f"obs 读入：{obs.shape}")

    coarse = obs["celltype_coarse"].astype(str)
    fine   = obs["celltype_fine"].astype(str)
    keep = coarse.isin(["CD4 T","CD4 Treg"]).values
    sel = np.where(keep)[0]
    print(f"待打分细胞（CD4 T ∪ CD4 Treg）= {len(sel):,}")
    print(fine[keep].value_counts().to_string())

    with h5py.File(path, "r") as f:
        gv = f["raw"]["var"]
        gidx_name = gv.attrs.get("_index", b"_index")
        gidx_name = gidx_name.decode() if isinstance(gidx_name, bytes) else str(gidx_name)
        genes = np.array([x.decode() if isinstance(x,bytes) else str(x) for x in gv[gidx_name][:]])
        n_genes = genes.shape[0]
        pos = {g:i for i,g in enumerate(genes)}
        sig_idx = np.array(sorted(pos[g] for g in H if g in pos))
        assert len(sig_idx) == len(H), f"覆盖不足：{len(sig_idx)}/{len(H)}"
        sig_set = set(sig_idx.tolist())
        AUC_T = int(round(n_genes * AUC_TOP_FRAC))
        print(f"\nraw 基因数 {n_genes:,}；signature 命中 {len(sig_idx)}；AUCell 阈值 rank ≤ {AUC_T}")

        X = f["raw"]["X"]
        indptr = X["indptr"][:]                    # 644,148 个 int，几 MB
        data_ds, ind_ds = X["data"], X["indices"]
        n_cells = indptr.shape[0]-1

        nsig = len(sig_idx)
        mn = nsig*(nsig+1)/2.0
        mx = nsig*MAXRANK - mn
        ucell = np.full(n_cells, np.nan, dtype=np.float32)
        aucell= np.full(n_cells, np.nan, dtype=np.float32)
        ndet  = np.zeros(n_cells, dtype=np.int32)   # signature 基因中被检出的个数

        selset = np.zeros(n_cells, dtype=bool); selset[sel] = True
        done = 0
        for s in range(0, n_cells, block):
            e = min(s+block, n_cells)
            if not selset[s:e].any(): continue
            lo, hi = indptr[s], indptr[e]
            dat = data_ds[lo:hi]; idc = ind_ds[lo:hi]
            for i in range(s, e):
                if not selset[i]: continue
                a0, a1 = indptr[i]-lo, indptr[i+1]-lo
                ix = idc[a0:a1]; vl = dat[a0:a1]
                k = ix.shape[0]
                if k == 0: continue
                order = np.argsort(-vl, kind="stable")
                r = np.empty(k, dtype=np.int32); r[order] = np.arange(1, k+1)
                inmask = np.fromiter((j in sig_set for j in ix), dtype=bool, count=k)
                hit_ranks = r[inmask]
                ndet[i] = hit_ranks.shape[0]
                # UCell：未检出者按 MAXRANK+1 计
                miss = nsig - hit_ranks.shape[0]
                tot = float(np.minimum(hit_ranks, MAXRANK+1).sum()) + miss*(MAXRANK+1)
                ucell[i] = 1.0 - (tot-mn)/mx
                # AUCell：top-T 恢复曲线下面积
                good = hit_ranks[hit_ranks <= AUC_T]
                aucell[i] = float((AUC_T - good + 1).sum())/(nsig*AUC_T) if good.size else 0.0
            done += int(selset[s:e].sum())
            if (s//block) % 10 == 0:
                el = time.time()-t0
                print(f"  {done:,}/{len(sel):,} 细胞  用时 {el/60:.1f} 分  "
                      f"预计剩余 {(el/max(done,1))*(len(sel)-done)/60:.1f} 分", flush=True)

    out = pd.DataFrame({
        "cell_index": sel,
        "sample": obs["orig.ident"].astype(str).values[sel],
        "patient": obs["SCARAB_ID"].astype(str).values[sel],
        "group": obs["Group"].astype(str).values[sel],
        "batch": obs["Batch"].astype(str).values[sel],
        "age": obs["Age"].values[sel],
        "sex": obs["Biologic_Sex"].astype(str).values[sel],
        "icu": obs["ICU"].astype(str).values[sel],
        "source": obs["Source"].astype(str).values[sel],
        "celltype_coarse": coarse.values[sel],
        "celltype_fine": fine.values[sel],
        "nCount_RNA": obs["nCount_RNA"].values[sel],
        "nFeature_RNA": obs["nFeature_RNA"].values[sel],
        "pct_mt": obs["percent.mt"].values[sel],
        "ucell": ucell[sel],
        "aucell": aucell[sel],
        "n_sig_detected": ndet[sel],
    })
    m = mp.set_index("original_label")
    for col, newname in [("is_memory_primary","pop_primary_memory"),
                         ("in_sens_ctl","pop_sens_with_CTL"),
                         ("in_secondary_allCD4","pop_secondary_allCD4"),
                         ("in_sens_allCD4_with_treg","pop_sens_allCD4_treg")]:
        out[newname] = out["celltype_fine"].map(m[col].astype(str).str.upper().eq("TRUE")).fillna(False)

    out.to_csv(f"{OUT}/gate1_cell_scores.csv.gz", index=False, compression="gzip")
    print(f"\n各群体细胞数：")
    for c in ["pop_primary_memory","pop_sens_with_CTL","pop_secondary_allCD4","pop_sens_allCD4_treg"]:
        print(f"  {c:<26s} {int(out[c].sum()):,}")
    print(f"\nsignature 基因平均检出 {out['n_sig_detected'].mean():.1f}/{len(H)}")
    print(f"总用时 {(time.time()-t0)/60:.1f} 分")
    print(f"已写出 {OUT}/gate1_cell_scores.csv.gz")
    print("\n本步骤未做任何分组比较。")

if __name__ == "__main__":
    if len(sys.argv) < 2: sys.exit("需要 h5ad 路径")
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 10000)
