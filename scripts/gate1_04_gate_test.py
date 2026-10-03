#!/usr/bin/env python3
"""
Gate 1 · Step 4 —— 第一次揭盲检验：CI-Sep vs CI-NS

严格按 gate1_SAP.md §6/§10/§11/§12 与 Amendment 3/4 执行。
所有规则在本脚本运行前已冻结，本脚本不做任何模型选择、阈值搜索或事后调整。

证据层级（A4.3，次序固定）：
  Primary       mean ΔUCell ~ group + z(median nFeature)      OLS + HC3
  Sensitivity 1 未调整的 patient mean ΔUCell 组间差
  Sensitivity 2 raw UCell 的 depth-adjusted model
  Sensitivity 3 AUCell（主要看方向一致性）
"""
import os, json
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy import stats

B = os.environ.get("ICOS_BASE") or os.path.expanduser("~/Desktop/Context-dependent_ICOS_analysis_results")
O = f"{B}/analysis_results"
PRIMARY_FINE = ["CD4 T CM", "CD4 T Eff/EM"]
REF, TEST = "CI-NS", "CI-Sep"          # 参照 = ICU non-sepsis；系数 = sepsis − non-sepsis
COVERAGE_MIN = 0.90                     # SAP §7


def fit(df, ycol, adjust):
    """OLS + HC3；返回 group 系数（TEST 相对 REF）。"""
    d = df.copy()
    d["g"] = (d.group == TEST).astype(float)
    X = [d["g"]]
    if adjust:
        z = (d.med_nfeat - d.med_nfeat.mean()) / d.med_nfeat.std(ddof=1)
        X.append(pd.Series(z.values, index=d.index, name="z_nfeat"))
    X = sm.add_constant(pd.concat(X, axis=1))
    m = sm.OLS(d[ycol].values, X).fit(cov_type="HC3")
    ci = m.conf_int().loc["g"] if hasattr(m.conf_int(), "loc") else m.conf_int()[1]
    ci = np.asarray(ci).ravel()
    return dict(beta=float(m.params["g"]), se=float(m.bse["g"]),
                ci_lo=float(ci[0]), ci_hi=float(ci[1]), p=float(m.pvalues["g"]), n=int(len(d)))


def unadjusted(df, ycol):
    a = df.loc[df.group == TEST, ycol].values
    b = df.loc[df.group == REF,  ycol].values
    diff = a.mean() - b.mean()
    sp = np.sqrt(((len(a)-1)*a.var(ddof=1) + (len(b)-1)*b.var(ddof=1)) / (len(a)+len(b)-2))
    d = diff / sp
    J = 1 - 3/(4*(len(a)+len(b)) - 9)            # Hedges 小样本校正
    t = stats.ttest_ind(a, b, equal_var=False)
    se = np.sqrt(a.var(ddof=1)/len(a) + b.var(ddof=1)/len(b))
    dfree = se**4 / (a.var(ddof=1)**2/(len(a)**2*(len(a)-1)) + b.var(ddof=1)**2/(len(b)**2*(len(b)-1)))
    tc = stats.t.ppf(0.975, dfree)
    return dict(n_test=len(a), n_ref=len(b), mean_test=float(a.mean()), mean_ref=float(b.mean()),
                diff=float(diff), ci_lo=float(diff-tc*se), ci_hi=float(diff+tc*se),
                hedges_g=float(d*J), p_welch=float(t.pvalue),
                p_wilcoxon=float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue))


def show(tag, r, direction_ok=None):
    if "beta" in r:
        print(f"  {tag:<44s} β = {r['beta']:+.5f}  95%CI [{r['ci_lo']:+.5f}, {r['ci_hi']:+.5f}]  P = {r['p']:.4f}")
    else:
        print(f"  {tag:<44s} Δ = {r['diff']:+.5f}  95%CI [{r['ci_lo']:+.5f}, {r['ci_hi']:+.5f}]  "
              f"g = {r['hedges_g']:+.2f}  P(Wilcoxon) = {r['p_wilcoxon']:.4f}")


def main():
    cal = pd.read_csv(f"{O}/gate1_cell_scores_calibrated.csv.gz",
                      usecols=["cell_index","patient","celltype_fine","nFeature_RNA","ucell_sig","delta_ucell"])
    raw = pd.read_csv(f"{O}/gate1_cell_scores.csv.gz",
                      usecols=["cell_index","group","icu","source","aucell",
                               "pop_primary_memory","pop_sens_with_CTL","pop_secondary_allCD4"])
    d = cal.merge(raw, on="cell_index", how="inner")
    assert len(d) == len(cal), "合并后行数不符"

    # ---- SAP §7 覆盖度检查 ----
    cov = d.pop_primary_memory.notna().mean()
    print(f"细胞总数 {len(d):,}；primary memory CD4 = {int(d.pop_primary_memory.sum()):,}")
    print(f"群体标签覆盖度 {cov*100:.1f}%（阈值 {COVERAGE_MIN*100:.0f}%）"
          f" → {'通过' if cov >= COVERAGE_MIN else '不通过，转 raw fallback'}")

    def patient_table(mask, name):
        s = d[mask]
        p = s.groupby(["patient","group"]).agg(
                n_cells=("delta_ucell","size"), mean_delta=("delta_ucell","mean"),
                mean_raw=("ucell_sig","mean"), mean_auc=("aucell","mean"),
                med_nfeat=("nFeature_RNA","median")).reset_index()
        print(f"\n[{name}] 患者 {len(p)} 例；" +
              "；".join(f"{g} {n}" for g, n in p.group.value_counts().items()))
        return p

    P = patient_table(d.pop_primary_memory.astype(bool), "primary memory CD4")
    P.to_csv(f"{O}/gate1_patient_level_unblinded.csv", index=False)
    two = P[P.group.isin([TEST, REF])].reset_index(drop=True)
    print(f"  纳入主对比：{len(two)} 例（{TEST} {int((two.group==TEST).sum())} / {REF} {int((two.group==REF).sum())}）")
    print(f"  每例细胞数 中位 {two.n_cells.median():.0f}（{two.n_cells.min()}–{two.n_cells.max()}）")
    print(f"  两组 median nFeature：{TEST} {two.loc[two.group==TEST,'med_nfeat'].median():.0f}"
          f" vs {REF} {two.loc[two.group==REF,'med_nfeat'].median():.0f}")

    res = {}
    print("\n================ 证据层级（次序在 A4.3 中冻结） ================")
    res["primary"] = fit(two, "mean_delta", adjust=True)
    show("PRIMARY   ΔUCell ~ group + z(nFeature)", res["primary"])
    res["sens1_unadjusted_delta"] = unadjusted(two, "mean_delta")
    show("SENS 1    ΔUCell 未调整组间差", res["sens1_unadjusted_delta"])
    res["sens2_raw_adjusted"] = fit(two, "mean_raw", adjust=True)
    show("SENS 2    raw UCell ~ group + z(nFeature)", res["sens2_raw_adjusted"])
    res["sens3_aucell_adjusted"] = fit(two, "mean_auc", adjust=True)
    show("SENS 3    AUCell ~ group + z(nFeature)", res["sens3_aucell_adjusted"])
    res["sens3_aucell_unadjusted"] = unadjusted(two, "mean_auc")
    show("SENS 3'   AUCell 未调整组间差", res["sens3_aucell_unadjusted"])

    signs = [np.sign(res["primary"]["beta"]), np.sign(res["sens1_unadjusted_delta"]["diff"]),
             np.sign(res["sens2_raw_adjusted"]["beta"]), np.sign(res["sens3_aucell_adjusted"]["beta"])]
    res["all_four_same_direction"] = bool(len(set(signs)) == 1)
    res["direction_matches_prespecified"] = bool(res["primary"]["beta"] < 0)
    print(f"\n  四层方向一致：{res['all_four_same_direction']}；"
          f"primary 方向符合预设（sepsis < non-sepsis）：{res['direction_matches_prespecified']}")

    # ---- LOPO 稳定性 ----
    bl = [fit(two.drop(i), "mean_delta", adjust=True)["beta"] for i in two.index]
    res["lopo"] = dict(min=float(np.min(bl)), max=float(np.max(bl)),
                       n_sign_flips=int(np.sum(np.sign(bl) != np.sign(res["primary"]["beta"]))))
    print(f"\n  LOPO（逐例剔除 {len(bl)} 次）：β 范围 {res['lopo']['min']:+.5f} … {res['lopo']['max']:+.5f}；"
          f"变号 {res['lopo']['n_sign_flips']} 次")

    # ---- 亚群分解（非 co-primary，A4/§4）----
    print("\n  亚群分解（非 co-primary，仅描述）：")
    for ct in PRIMARY_FINE:
        sub = patient_table((d.celltype_fine == ct), f"仅 {ct}")
        sub = sub[sub.group.isin([TEST, REF])]
        r = fit(sub, "mean_delta", adjust=True)
        res[f"decomp_{ct}"] = r
        show(f"  {ct}", r)

    # ---- 敏感性群体 ----
    print("\n  其他预设群体：")
    for flag, nm in [("pop_sens_with_CTL", "memory CD4 + CD4 CTL"),
                     ("pop_secondary_allCD4", "secondary all-CD4（不含 Treg）")]:
        sub = patient_table(d[flag].astype(bool), nm)
        sub = sub[sub.group.isin([TEST, REF])]
        r = fit(sub, "mean_delta", adjust=True)
        res[f"pop_{flag}"] = r
        show(f"  {nm}", r)

    # ---- 生物学锚定：NHC（不参与主假设检验，SAP §6）----
    if (P.group == "NHC").any():
        nhc = P[P.group == "NHC"]
        res["nhc_anchor_mean_delta"] = float(nhc.mean_delta.mean())
        print(f"\n  生物学锚定 NHC（n={len(nhc)}，不参与主假设检验）："
              f"mean ΔUCell {nhc.mean_delta.mean():+.5f}；"
              f"{TEST} {two.loc[two.group==TEST,'mean_delta'].mean():+.5f}；"
              f"{REF} {two.loc[two.group==REF,'mean_delta'].mean():+.5f}")

    with open(f"{O}/gate1_gate_test_results.json", "w") as fh:
        json.dump(res, fh, indent=2, ensure_ascii=False)
    print(f"\n已写出 {O}/gate1_gate_test_results.json 与 gate1_patient_level_unblinded.csv")


if __name__ == "__main__":
    main()
