import gzip, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
BASE=Path('/Users/heavenlake/Desktop/Context-dependent_ICOS_analysis_results')
cell=pd.read_csv(gzip.open(BASE/'analysis_results/gate1_cell_scores_calibrated.csv.gz','rt'))
pat=pd.read_csv(BASE/'analysis_results/gate1_patient_level_unblinded.csv')
cell=cell.merge(pat[['patient','group','mean_raw','mean_delta','mean_auc','med_nfeat']],on='patient',how='left')
rng=np.random.default_rng(20260915); sub=cell.iloc[rng.choice(len(cell),size=min(25000,len(cell)),replace=False)]
fig,ax=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
col={'NHC':'#777777','CI-NS':'#4C78A8','CI-Sep':'#E45756'}
for g,d in sub.groupby('group'): ax[0,0].scatter(d.nFeature_RNA,d.ucell_sig,s=3,alpha=.18,color=col.get(g,'#777'),label=g)
ax[0,0].set(xlabel='Detected genes per cell (nFeature_RNA)',ylabel='Raw UCell score',title='A  Cell-level raw depth dependence'); ax[0,0].legend(frameon=False,markerscale=3,fontsize=8)
for g,d in sub.groupby('group'): ax[0,1].scatter(d.nFeature_RNA,d.delta_ucell,s=3,alpha=.18,color=col.get(g,'#777'),label=g)
ax[0,1].set(xlabel='Detected genes per cell (nFeature_RNA)',ylabel='Calibrated score (ΔUCell)',title='B  Cell-level calibrated dependence'); ax[0,1].text(.03,.95,'Spearman ρ = 0.527 → 0.202',transform=ax[0,1].transAxes,va='top',fontsize=9)
for g,d in pat.groupby('group'):
 ax[1,0].scatter(d.med_nfeat,d.mean_raw,s=28,color=col.get(g,'#777'),label=g)
 ax[1,0].scatter(d.med_nfeat,d.mean_delta,s=28,facecolors='none',edgecolors=col.get(g,'#777'))
ax[1,0].set(xlabel='Patient median nFeature_RNA',ylabel='Patient-level score',title='C  Donor-level depth leakage'); ax[1,0].text(.03,.95,'Raw ρ = 0.832\nCalibrated ρ = 0.537',transform=ax[1,0].transAxes,va='top',fontsize=9); ax[1,0].legend(frameon=False,fontsize=8)
diag=json.load(open(BASE/'analysis_results/gate1_calibration_diagnostics.json'))
labels=['Signature','Matched controls']; det=[diag['sig_det_f_median'],diag['ctrl_det_f_median']]; expr=[diag['sig_mean_e_median'],diag['ctrl_mean_e_median']]
x=np.arange(2); w=.34; ax[1,1].bar(x-w/2,det,w,label='Detection frequency',color='#59A14F'); ax[1,1].bar(x+w/2,expr,w,label='Mean expression',color='#F28E2B'); ax[1,1].set_xticks(x,labels); ax[1,1].set_ylabel('Median matched metric'); ax[1,1].set_title('D  Matched-control diagnostics'); ax[1,1].legend(frameon=False,fontsize=8)
for a in ax.ravel(): a.spines['top'].set_visible(False); a.spines['right'].set_visible(False)
fig.savefig(BASE/'figures/figure3_detection_depth_leakage.png',dpi=300,bbox_inches='tight')
