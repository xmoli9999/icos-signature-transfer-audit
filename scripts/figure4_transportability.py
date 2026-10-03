from pathlib import Path
import json, pandas as pd, numpy as np, matplotlib.pyplot as plt
BASE=Path('/Users/heavenlake/Desktop/Context-dependent_ICOS_analysis_results')
res=json.load(open(BASE/'analysis_results/gate1_gate_test_results.json'))
pat=pd.read_csv(BASE/'analysis_results/gate1_patient_level_unblinded.csv')
fig,ax=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
col={'NHC':'#777777','CI-NS':'#4C78A8','CI-Sep':'#E45756'}
for g,d in pat.groupby('group'): ax[0,0].scatter(np.full(len(d),{'NHC':0,'CI-NS':1,'CI-Sep':2}[g]),d.mean_delta,s=30,color=col[g],label=g)
ax[0,0].set_xticks([0,1,2],['Healthy','ICU non-sepsis','ICU sepsis']); ax[0,0].set_ylabel('Patient-level calibrated score'); ax[0,0].set_title('A  Primary external validation',loc='left',weight='bold'); ax[0,0].text(.03,.95,'β = −0.00138\n95% CI −0.00531 to 0.00255\nP = 0.493',transform=ax[0,0].transAxes,va='top',fontsize=9)
items=[('Primary UCell',res['primary']),('Unadjusted',res['sens1_unadjusted_delta']),('Raw-adjusted',res['sens2_raw_adjusted']),('AUCell',res['sens3_aucell_adjusted'])]
y=np.arange(len(items)); est=[x[1].get('beta',x[1].get('diff')) for x in items]; lo=[x[1].get('ci_lo') for x in items]; hi=[x[1].get('ci_hi') for x in items]
ax[0,1].axvline(0,color='0.3'); ax[0,1].errorbar(est,y,xerr=[np.array(est)-np.array(lo),np.array(hi)-np.array(est)],fmt='o',capsize=3,color='#4C78A8'); ax[0,1].set_yticks(y,[x[0] for x in items]); ax[0,1].set_xlabel('Effect estimate (95% CI)'); ax[0,1].set_title('B  Prespecified sensitivity forest',loc='left',weight='bold')
ax[1,0].axis('off'); ax[1,0].text(.02,.94,'C  Outcome-blind transportability audit',fontsize=12,weight='bold'); stages=[('8','repositories searched'),('2','closest candidate cohorts'),('0','meeting all 5 frozen criteria')]; ys=[.72,.48,.24]
for i,((n,s),yy) in enumerate(zip(stages,ys)):
 ax[1,0].text(.20,yy,n,ha='center',va='center',fontsize=18,weight='bold',color=['#4C78A8','#F28E2B','#E45756'][i]); ax[1,0].text(.30,yy,s,ha='left',va='center',fontsize=11)
 if i<2: ax[1,0].annotate('',xy=(.20,ys[i+1]+.07),xytext=(.20,yy-.07),arrowprops=dict(arrowstyle='->',lw=1.5,color='0.4'))
coh=['Kwok 2023','Reyes 2020']; sep=[14,3]; ctl=[4,7]; x=np.arange(2); w=.34; ax[1,1].bar(x-w/2,sep,w,label='Sepsis retained',color='#E45756'); ax[1,1].bar(x+w/2,ctl,w,label='Non-sepsis retained',color='#4C78A8'); ax[1,1].axhline(5,color='0.25',ls='--',label='Viability threshold (5)'); ax[1,1].set_xticks(x,coh); ax[1,1].set_ylabel('Patients passing ≥30 mapped memory-CD4 cells'); ax[1,1].set_title('D  Gate 2b viability stop',loc='left',weight='bold'); ax[1,1].legend(frameon=False,fontsize=8)
for a in ax.ravel(): a.spines['top'].set_visible(False); a.spines['right'].set_visible(False)
fig.savefig(BASE/'figures/figure4_transportability_audit.png',dpi=300,bbox_inches='tight')
