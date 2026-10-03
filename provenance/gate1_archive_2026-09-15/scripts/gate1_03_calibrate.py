#!/usr/bin/env python3
"""
Gate 1 · Step 3 —— Amendment 3：expression/detection-matched empirical-null calibration

严格按 gate1_SAP.md Amendment 3 §A3.1–A3.5 实现，参数全部写死在下面的
FROZEN 区，运行后不得人工增删。

  ΔUCell_i = UCell_signature,i − (1/K)·Σ_k UCell_control_k,i
  K = 100, seed = 20260915
  control universe = /raw/var 中、在 primary memory CD4 内 detection frequency > 0
                     的基因，排除 107 个 signature 基因本身
  matching = detection frequency（20 分位箱） × mean normalized expression（20 分位箱）
             bin 内候选不足时按 Chebyshev 半径 r=1,2,3… 逐层向相邻 bin 扩展
  抽样 = 同一 replicate 内不放回；不同 replicate 之间允许重复

【本脚本不读取、不使用、不输出任何 group label。】
所有分组比较在 Step 4 揭盲时进行。

用法: python3 gate1_03_calibrate.py /path/to/file.h5ad [块大小]
"""
import sys, os, time, json
import numpy as np, pandas as pd, h5py

# ----------------------------- FROZEN -----------------------------
K_REPLICATES = 100
SEED         = 20260915
MAXRANK      = 1500
N_BINS_DET   = 20
N_BINS_EXPR  = 20
MAX_RADIUS   = 20          # Chebyshev 扩展上限；超出即报错，不静默降级
PRIMARY_FINE = ["CD4 T CM", "CD4 T Eff/EM"]          # primary memory CD4（SAP §2）
SCORE_COARSE = ["CD4 T", "CD4 Treg"]                 # 打分细胞集（与 Step 2 一致）
# ------------------------------------------------------------------

# 项目根目录：优先读环境变量 ICOS_BASE（本机 Desktop 在挂载点下，~ 不指向它）
B = os.environ.get("ICOS_BASE") or os.path.expanduser("~/Desktop/Context-dependent_ICOS_analysis_results")
SIGFILE = f"{B}/metadata/rat_107gene_frozen_signature.csv"
OUT     = f"{B}/analysis_results"


def read_obs_cols(f, names):
    """直接用 h5py 读 /obs 的指定列（含 categorical 解码）。
    不走 anndata backed 模式——该模式会尝试整体读入 /layers（scVI_normalized 12 GB）。"""
    g = f["obs"]
    out = {}
    for nm in names:
        o = g[nm]
        if isinstance(o, h5py.Group):          # categorical
            cats = o["categories"][:]
            cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats])
            codes = o["codes"][:]
            v = np.where(codes >= 0, cats[np.clip(codes, 0, None)], "nan")
        else:
            v = o[:]
            if v.dtype.kind in "SO":
                v = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v])
        out[nm] = v
    return pd.DataFrame(out)


def stream_cells(f, indptr, data_ds, ind_ds, selmask, block):
    """按全局行块流式产出 (i, indices, values)，峰值内存 < 1 GB。"""
    n = indptr.shape[0] - 1
    for s in range(0, n, block):
        e = min(s + block, n)
        if not selmask[s:e].any():
            continue
        lo, hi = indptr[s], indptr[e]
        dat = data_ds[lo:hi]; idc = ind_ds[lo:hi]
        for i in range(s, e):
            if not selmask[i]:
                continue
            a0, a1 = indptr[i] - lo, indptr[i + 1] - lo
            yield i, idc[a0:a1], dat[a0:a1]


def main(path, block=10000):
    t0 = time.time()
    rng = np.random.default_rng(SEED)

    sig = pd.read_csv(SIGFILE)
    H = sorted({g for g in sig["human_symbol"].dropna().astype(str) if g})
    print(f"frozen signature: {len(H)} 个 human symbol（不做任何增删/加权）")

    with h5py.File(path, "r") as _f:
        obs = read_obs_cols(_f, ["orig.ident", "SCARAB_ID", "celltype_coarse",
                                 "celltype_fine", "nCount_RNA", "nFeature_RNA"])
    coarse = obs["celltype_coarse"].astype(str).values
    fine   = obs["celltype_fine"].astype(str).values

    prim_mask  = np.isin(fine, PRIMARY_FINE)
    score_mask = np.isin(coarse, SCORE_COARSE)
    print(f"primary memory CD4（校正量的计算集）= {prim_mask.sum():,}")
    print(f"打分细胞集（CD4 T ∪ CD4 Treg）    = {score_mask.sum():,}")

    with h5py.File(path, "r") as f:
        gv = f["raw"]["var"]
        gi = gv.attrs.get("_index", b"_index")
        gi = gi.decode() if isinstance(gi, bytes) else str(gi)
        genes = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in gv[gi][:]])
        n_genes = genes.shape[0]
        pos = {g: i for i, g in enumerate(genes)}
        sig_idx = np.array(sorted(pos[g] for g in H if g in pos))
        assert len(sig_idx) == len(H), f"覆盖不足 {len(sig_idx)}/{len(H)}"
        print(f"raw 基因数 {n_genes:,}；signature 命中 {len(sig_idx)}/{len(H)}")

        X = f["raw"]["X"]
        indptr  = X["indptr"][:]
        data_ds, ind_ds = X["data"], X["indices"]
        n_cells = indptr.shape[0] - 1
        assert n_cells == obs.shape[0], "细胞数与 obs 不一致"

        CSET = f"{OUT}/gate1_control_sets_geneidx.npy"
        if os.path.exists(CSET):
            ctrl_sets = np.load(CSET)
            assert ctrl_sets.shape == (K_REPLICATES, sig_idx.size), "已存在的对照集形状不符"
            with open(f"{OUT}/gate1_calibration_diagnostics.json") as fh:
                diag = json.load(fh)
            print(f"\n[载入] 已有对照集 {CSET}，跳过 PASS 1 与采样（保证与首次运行完全一致）")
        else:
            # ---------- PASS 1：blinded，仅在 primary memory CD4 内统计基因属性 ----------
            print("\n[PASS 1] 统计 detection frequency 与 mean normalized expression …")
            det_n = np.zeros(n_genes, dtype=np.int64)
            sum_e = np.zeros(n_genes, dtype=np.float64)
            n_prim = int(prim_mask.sum()); done = 0
            # 按块聚合后一次 bincount，比逐细胞 np.add.at 快一个量级；结果完全等价
            for s0 in range(0, n_cells, block):
                e0 = min(s0 + block, n_cells)
                sub = np.flatnonzero(prim_mask[s0:e0])
                if sub.size == 0:
                    continue
                lo, hi = indptr[s0], indptr[e0]
                dat = data_ds[lo:hi]; idc = ind_ds[lo:hi]
                parts_i, parts_v = [], []
                for i in sub + s0:
                    a0, a1 = indptr[i] - lo, indptr[i + 1] - lo
                    parts_i.append(idc[a0:a1]); parts_v.append(dat[a0:a1])
                ii = np.concatenate(parts_i); vv = np.concatenate(parts_v)
                det_n += np.bincount(ii, minlength=n_genes)
                sum_e += np.bincount(ii, weights=vv.astype(np.float64), minlength=n_genes)
                done += sub.size
                if done % 20000 < block:
                    print(f"  {done:,}/{n_prim:,}  {(time.time()-t0)/60:.1f} 分", flush=True)

            det_f = det_n / float(n_prim)
            mean_e = sum_e / float(n_prim)
            print(f"PASS 1 完成，用时 {(time.time()-t0)/60:.1f} 分")

            # ---------- 构造匹配的 control universe ----------
            in_sig = np.zeros(n_genes, dtype=bool); in_sig[sig_idx] = True
            universe = np.flatnonzero((det_f > 0) & (~in_sig))
            print(f"\ncontrol universe = {universe.size:,} 个基因"
                  f"（det_f>0 且非 signature）")
            miss_sig = sig_idx[det_f[sig_idx] == 0]
            if miss_sig.size:
                print(f"  注意：{miss_sig.size} 个 signature 基因在 primary memory CD4 内检出为 0："
                      f"{', '.join(genes[miss_sig][:10])}")

            # 分位箱边界由 universe 定义，signature 基因用同一套边界归箱
            qd = np.quantile(det_f[universe],  np.linspace(0, 1, N_BINS_DET + 1)[1:-1])
            qe = np.quantile(mean_e[universe], np.linspace(0, 1, N_BINS_EXPR + 1)[1:-1])
            bin_d = np.searchsorted(qd, det_f,  side="right")
            bin_e = np.searchsorted(qe, mean_e, side="right")

            pool = {}
            for g in universe:
                pool.setdefault((bin_d[g], bin_e[g]), []).append(g)
            pool = {k: np.array(v) for k, v in pool.items()}

            # ---------- 抽取 K 组匹配对照 ----------
            print(f"[采样] K={K_REPLICATES}, seed={SEED} …")
            radius_log = []
            ctrl_sets = np.empty((K_REPLICATES, sig_idx.size), dtype=np.int32)
            for k in range(K_REPLICATES):
                used = set()
                for j, g in enumerate(sig_idx):
                    bd, be = bin_d[g], bin_e[g]
                    pick = None
                    for r in range(0, MAX_RADIUS + 1):
                        cand = []
                        for dd in range(-r, r + 1):
                            for ee in range(-r, r + 1):
                                if r > 0 and max(abs(dd), abs(ee)) != r:
                                    continue      # 只取该半径的"外环"，保证由近及远
                                p = pool.get((bd + dd, be + ee))
                                if p is not None:
                                    cand.append(p)
                        if cand:
                            cand = np.concatenate(cand)
                            cand = cand[~np.isin(cand, list(used))] if used else cand
                            if cand.size:
                                pick = int(rng.choice(cand))
                                radius_log.append(r)
                                break
                    if pick is None:
                        raise RuntimeError(f"replicate {k} gene {genes[g]} 在半径 {MAX_RADIUS} 内无候选")
                    used.add(pick)
                    ctrl_sets[k, j] = pick
            rl = np.array(radius_log)
            print(f"  扩展半径分布：r=0 {int((rl==0).sum()):,} 次，r≥1 {int((rl>0).sum()):,} 次，"
                  f"最大 {int(rl.max())}")

            # 匹配质量（blinded 诊断）
            mq = {
                "sig_det_f_median":   float(np.median(det_f[sig_idx])),
                "ctrl_det_f_median":  float(np.median(det_f[ctrl_sets.ravel()])),
                "sig_mean_e_median":  float(np.median(mean_e[sig_idx])),
                "ctrl_mean_e_median": float(np.median(mean_e[ctrl_sets.ravel()])),
                "n_unique_ctrl_genes": int(np.unique(ctrl_sets).size),
            }
            print(f"  匹配质量 det_f 中位数 sig {mq['sig_det_f_median']:.4f} vs "
                  f"ctrl {mq['ctrl_det_f_median']:.4f}；"
                  f"mean_expr 中位数 sig {mq['sig_mean_e_median']:.4f} vs "
                  f"ctrl {mq['ctrl_mean_e_median']:.4f}")
            print(f"  被用到的不同对照基因 {mq['n_unique_ctrl_genes']:,} 个")

            np.save(CSET, ctrl_sets)
            pd.DataFrame({"set_id": np.repeat(np.arange(K_REPLICATES), sig_idx.size),
                          "matched_to": np.tile(genes[sig_idx], K_REPLICATES),
                          "control_gene": genes[ctrl_sets.ravel()]}
                         ).to_csv(f"{OUT}/gate1_control_sets.csv.gz", index=False, compression="gzip")
            diag = dict(mq)
            diag.update({"K": K_REPLICATES, "seed": SEED, "maxrank": MAXRANK,
                         "n_bins_det": N_BINS_DET, "n_bins_expr": N_BINS_EXPR,
                         "n_primary_memory_cells": int(prim_mask.sum()),
                         "n_scored_cells": int(score_mask.sum()),
                         "universe_size": int(universe.size),
                         "radius_r0_frac": float((rl == 0).mean()),
                         "radius_max": int(rl.max())})
            with open(f"{OUT}/gate1_calibration_diagnostics.json", "w") as fh:
                json.dump(diag, fh, indent=2, ensure_ascii=False)
            print("[STAGE 1 完成] 对照集与诊断已落盘；再次运行本脚本即进入 PASS 2。")
            if os.environ.get("ICOS_STAGE") == "1":
                return
        # ---------- 基因 → set-id 隶属索引（signature = set 0，对照 = 1..K） ----------
        rows = [(int(g), 0) for g in sig_idx]
        for k in range(K_REPLICATES):
            rows += [(int(g), k + 1) for g in ctrl_sets[k]]
        rows.sort()
        memb_g = np.array([r[0] for r in rows], dtype=np.int64)
        memb_s = np.array([r[1] for r in rows], dtype=np.int32)
        memb_len = np.zeros(n_genes, dtype=np.int32)
        np.add.at(memb_len, memb_g, 1)
        memb_start = np.zeros(n_genes + 1, dtype=np.int64)
        np.cumsum(memb_len, out=memb_start[1:])
        has_memb = memb_len > 0
        n_sets = K_REPLICATES + 1

        # ---------- PASS 2：逐细胞一次排秩，累加到全部 101 个基因集 ----------
        print(f"\n[PASS 2] 对 {int(score_mask.sum()):,} 个细胞计算 UCell × {n_sets} 组 …")
        nsig = sig_idx.size
        mn = nsig * (nsig + 1) / 2.0
        mx = nsig * MAXRANK - mn
        CAP = MAXRANK + 1

        u_sig  = np.full(n_cells, np.nan, dtype=np.float32)
        c_mean = np.full(n_cells, np.nan, dtype=np.float32)
        c_sd   = np.full(n_cells, np.nan, dtype=np.float32)
        n_score = int(score_mask.sum()); done = 0; t1 = time.time()

        for i, ix, vl in stream_cells(f, indptr, data_ds, ind_ds, score_mask, block):
            k = ix.shape[0]
            if k == 0:
                continue
            order = np.argsort(-vl, kind="stable")
            r = np.empty(k, dtype=np.int32); r[order] = np.arange(1, k + 1)
            np.minimum(r, CAP, out=r)

            hit = np.flatnonzero(has_memb[ix])
            tot = np.zeros(n_sets, dtype=np.float64)
            cnt = np.zeros(n_sets, dtype=np.int32)
            if hit.size:
                gj = ix[hit]; rj = r[hit]
                lens = memb_len[gj]
                total = int(lens.sum())
                starts = memb_start[gj]
                off = np.arange(total) - np.repeat(np.cumsum(lens) - lens, lens)
                flat = np.repeat(starts, lens) + off
                sid = memb_s[flat]
                np.add.at(tot, sid, np.repeat(rj, lens))
                np.add.at(cnt, sid, 1)
            tot += (nsig - cnt) * float(CAP)
            uc = 1.0 - (tot - mn) / mx
            u_sig[i]  = uc[0]
            c_mean[i] = uc[1:].mean()
            c_sd[i]   = uc[1:].std(ddof=1)

            done += 1
            if done % 20000 == 0:
                el = time.time() - t1
                print(f"  {done:,}/{n_score:,}  用时 {el/60:.1f} 分  "
                      f"剩余约 {(el/done)*(n_score-done)/60:.1f} 分", flush=True)

    sel = np.flatnonzero(score_mask)
    out = pd.DataFrame({
        "cell_index": sel,
        "sample": obs["orig.ident"].astype(str).values[sel],
        "patient": obs["SCARAB_ID"].astype(str).values[sel],
        "celltype_fine": fine[sel],
        "nFeature_RNA": obs["nFeature_RNA"].values[sel],
        "ucell_sig": u_sig[sel],
        "ucell_ctrl_mean": c_mean[sel],
        "ucell_ctrl_sd": c_sd[sel],
        "delta_ucell": (u_sig - c_mean)[sel],
    })
    out.to_csv(f"{OUT}/gate1_cell_scores_calibrated.csv.gz", index=False, compression="gzip")
    print("对照基因表与 ΔUCell 已落盘。\n")

    # 与 Step 2 的一致性校验（同一 signature、同一表达源，应完全一致）
    prev = pd.read_csv(f"{OUT}/gate1_cell_scores.csv.gz", usecols=["cell_index", "ucell"])
    chk = out[["cell_index", "ucell_sig"]].merge(prev, on="cell_index", how="inner")
    dmax = float(np.nanmax(np.abs(chk["ucell_sig"] - chk["ucell"])))
    print(f"\n[一致性] 与 Step 2 的 signature UCell 最大绝对差 = {dmax:.3e}（应 ~1e-7）")

    # blinded 深度相关性（不涉及任何 group label）
    def spearman(x, y):
        x = np.asarray(x, float); y = np.asarray(y, float)
        rx = pd.Series(x).rank().values; ry = pd.Series(y).rank().values
        return float(np.corrcoef(rx, ry)[0, 1])
    m = out["pop_flag"] if "pop_flag" in out else np.isin(out["celltype_fine"], PRIMARY_FINE)
    sub = out[m]
    r_raw = spearman(sub["ucell_sig"], sub["nFeature_RNA"])
    r_del = spearman(sub["delta_ucell"], sub["nFeature_RNA"])
    print(f"[blinded QC] primary memory CD4 内 Spearman ρ vs nFeature_RNA："
          f"raw UCell {r_raw:.3f} → ΔUCell {r_del:.3f}")

    diag.update({
        "consistency_max_abs_diff_vs_step2": dmax,
        "spearman_nFeature_raw_ucell": float(r_raw),
        "spearman_nFeature_delta_ucell": float(r_del),
    })
    with open(f"{OUT}/gate1_calibration_diagnostics.json", "w") as fh:
        json.dump(diag, fh, indent=2, ensure_ascii=False)


    print(f"\n已写出：")
    print(f"  {OUT}/gate1_cell_scores_calibrated.csv.gz")
    print(f"  {OUT}/gate1_calibration_diagnostics.json")
    print(f"  {OUT}/gate1_control_sets.csv.gz（可复现的 100×107 对照基因表）")
    print(f"总用时 {(time.time()-t0)/60:.1f} 分")
    print("\n本步骤未读取、未使用、未输出任何 group label。")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("需要 h5ad 路径")
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 10000)
