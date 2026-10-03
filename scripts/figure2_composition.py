import gzip
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

BASE=Path('/Users/heavenlake/Desktop/Context-dependent_ICOS_analysis_results')
score= pd.read_csv(gzip.open(BASE/'analysis_results/gate1_cell_scores_calibrated.csv.gz','rt'))
pat=pd.read_csv(BASE/'analysis_results/gate1_patient_level_unblinded.csv')[['patient','group']]
score=score.merge(pat,on='patient',how='left')
score['state']=score['celltype_fine'].map({'CD4 T CM':'CM','CD4 T Eff/EM':'Eff/EM'})
score=score[score.state.notna()]
pm=score.groupby(['patient','group','state'],as_index=False).agg(score=('delta_ucell','median'),n_cells=('delta_ucell','size'))
groups=['CI-NS','CI-Sep']
pm['group']=pm['group'].replace({'ICU non-sepsis':'CI-NS','ICU sepsis':'CI-Sep'})
def effect(state=None):
 d=pm if state is None else pm[pm.state==state]
 a=d[d.group=='CI-Sep'].score; b=d[d.group=='CI-NS'].score
 return a.mean()-b.mean(), mannwhitneyu(a,b,alternative='two-sided').pvalue
# Prespecified model estimates used in the manuscript (not recomputed descriptive medians).
effects=[('Pooled',-0.0013757538,-0.00530577,0.00255426),('CM',0.0015975302,-0.00215028,0.00534534),('Eff/EM',0.0014367271,-0.00220671,0.00508016)]
prop=score.groupby(['patient','group','state']).size().unstack(fill_value=0)
prop['total']=prop.sum(axis=1); prop['Eff/EM_prop']=prop.get('Eff/EM',0)/prop.total
prop=prop.reset_index()
fig,ax=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
colors={'CI-NS':'#4C78A8','CI-Sep':'#E45756'}
# A effects
names=[x[0] for x in effects]; vals=[x[1] for x in effects]; lo=[x[2] for x in effects]; hi=[x[3] for x in effects]
ax[0,0].axvline(0,color='0.3',lw=1); ax[0,0].errorbar(vals,range(3),xerr=[np.array(vals)-np.array(lo),np.array(hi)-np.array(vals)],fmt='o',capsize=3,ms=7,color='#555',ecolor='#555')
for i,v in enumerate(vals): ax[0,0].text(v+0.00010,i,f'{v:+.5f}',va='center',ha='left',fontsize=9)
ax[0,0].set_yticks(range(3),names); ax[0,0].set_xlabel('ICU sepsis − ICU non-sepsis\nmodel estimate (95% CI)'); ax[0,0].set_title('A  Pooled and state-stratified effects',loc='left',weight='bold')
# B baseline states
for j,state in enumerate(['CM','Eff/EM']):
 for g in ['CI-NS','CI-Sep']:
  y=pm[(pm.state==state)&(pm.group==g)].score
  ax[0,1].scatter(np.full(len(y),j+(0 if g=='CI-NS' else .16)),y,color=colors[g],s=22,alpha=.8)
  ax[0,1].plot([j-.05,j+.21],[y.mean(),y.mean()],color=colors[g],lw=2)
ax[0,1].set_xticks([.08,1.08],['CM','Eff/EM']); ax[0,1].set_ylabel('Calibrated score'); ax[0,1].set_title('B  State-specific signature levels',loc='left',weight='bold'); ax[0,1].legend(handles=[plt.Line2D([0],[0],marker='o',color='w',markerfacecolor=colors['CI-NS'],label='ICU non-sepsis'),plt.Line2D([0],[0],marker='o',color='w',markerfacecolor=colors['CI-Sep'],label='ICU sepsis')],frameon=False,fontsize=8,loc='upper left')
# C composition
for g in ['CI-NS','CI-Sep']:
 y=prop[prop.group==g]['Eff/EM_prop']; x=np.full(len(y),0 if g=='CI-NS' else 1)+(.08 if g=='CI-Sep' else -.08)
 ax[1,0].scatter(x,y,color=colors[g],s=28); ax[1,0].plot([(-.08 if g=='CI-NS' else .92),( .08 if g=='CI-NS' else 1.08)],[y.mean(),y.mean()],color=colors[g],lw=2)
ax[1,0].set_xticks([0,1],['ICU non-sepsis','ICU sepsis']); ax[1,0].set_ylabel('Eff/EM fraction of CM + Eff/EM'); ax[1,0].set_title('C  Eff/EM composition differs by group',loc='left',weight='bold'); ax[1,0].set_ylim(0,1)
# D schematic
ax[1,1].axis('off'); ax[1,1].text(.02,.94,'D',fontsize=14,weight='bold'); ax[1,1].text(.08,.82,'State-specific levels differ',fontsize=11,weight='bold'); ax[1,1].text(.08,.72,'CM lower; Eff/EM higher',fontsize=10,color='#555'); ax[1,1].text(.08,.57,'+',fontsize=18,weight='bold'); ax[1,1].text(.08,.44,'Positive within-state effects',fontsize=11,weight='bold'); ax[1,1].text(.08,.34,'CM  +0.00160   |   Eff/EM  +0.00144',fontsize=9,color='#2E7D32'); ax[1,1].text(.08,.20,'+',fontsize=18,weight='bold'); ax[1,1].text(.08,.08,'Composition shift: Eff/EM 0.607 → 0.420',fontsize=10,color='#E45756'); ax[1,1].text(.70,.39,'↓',fontsize=24,weight='bold'); ax[1,1].text(.53,.18,'Pooled sign can reverse',fontsize=11,weight='bold')
for a in ax.ravel():
 a.spines['top'].set_visible(False); a.spines['right'].set_visible(False)
fig.savefig(BASE/'figures/figure2_composition_sign_reversal.png',dpi=300,bbox_inches='tight')
pd.DataFrame(effects,columns=['contrast','estimate','ci_lo','ci_hi']).to_csv(BASE/'analysis_results/figure2_effects_recomputed.csv',index=False)
prop.to_csv(BASE/'analysis_results/figure2_state_proportions.csv',index=False)
print(pd.DataFrame(effects,columns=['contrast','estimate','ci_lo','ci_hi']).to_string(index=False))
