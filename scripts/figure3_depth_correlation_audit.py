#!/usr/bin/env python3
"""Direct computation of every detection-depth correlation reported in Figure 3.

Each value is recomputed from the archived score files and asserted against the
value previously hard-coded in the plotting script. Analysis sets are stated
explicitly because they differ between panels:
  cell level   : all primary-memory cells, ALL 47 donors (healthy included)
  donor level  : all 47 donors (healthy included)
  within group : mean of the two ICU groups' within-group Spearman (19 + 19)
"""
import os, numpy as np, pandas as pd
from scipy.stats import spearmanr
B = os.environ.get("ICOS_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
O = f"{B}/analysis_results"
cal = pd.read_csv(f"{O}/gate1_cell_scores_calibrated.csv.gz",
                  usecols=["cell_index","patient","nFeature_RNA","ucell_sig","delta_ucell"])
raw = pd.read_csv(f"{O}/gate1_cell_scores.csv.gz", usecols=["cell_index","group","pop_primary_memory"])
m = cal.merge(raw, on="cell_index")
m = m[m.pop_primary_memory.astype(bool)]
P = pd.read_csv(f"{O}/gate1_patient_level_unblinded.csv")
ICU = ["CI-NS","CI-Sep"]

cell_raw  = spearmanr(m.nFeature_RNA, m.ucell_sig).statistic
cell_cal  = spearmanr(m.nFeature_RNA, m.delta_ucell).statistic
don_raw   = spearmanr(P.med_nfeat, P.mean_raw).statistic
don_cal   = spearmanr(P.med_nfeat, P.mean_delta).statistic
icu = P[P.group.isin(ICU)]
wg = lambda y: float(np.mean([spearmanr(s.med_nfeat, s[y]).statistic for _, s in icu.groupby("group")]))
wg_raw, wg_cal = wg("mean_raw"), wg("mean_delta")

EXPECT = {"cell_raw":0.527,"cell_calibrated":0.202,"donor_raw":0.832,"donor_calibrated":0.537,
          "within_group_raw":0.825,"within_group_calibrated":0.595}
got = {"cell_raw":cell_raw,"cell_calibrated":cell_cal,"donor_raw":don_raw,"donor_calibrated":don_cal,
       "within_group_raw":wg_raw,"within_group_calibrated":wg_cal}
for k,v in got.items():
    assert abs(round(v,3)-EXPECT[k]) < 5e-4, f"{k}: recomputed {v:.4f} does not match reported {EXPECT[k]}"

rows = [(k, got[k]) for k in EXPECT] + [
    ("n_cells_primary_memory_all_donors", float(len(m))),
    ("n_cells_primary_memory_icu_only", float((m.group.isin(ICU)).sum())),
    ("n_donors_all", float(len(P))), ("n_donors_icu", float(len(icu))),
    ("plot_subsample_n", 25000.0), ("plot_subsample_seed", 3.0)]
pd.DataFrame(rows, columns=["quantity","value"]).to_csv(f"{O}/figure3_depth_correlations.csv", index=False)
for k,v in rows: print(f"{k:38s} {v}")
