from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
BASE=Path('/Users/heavenlake/Desktop/Context-dependent_ICOS_analysis_results')
fig,ax=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
stages=['Candidates','Rat–mouse\n1:1','Rat–mouse–human\n1:1','Frozen transfer\nobject']; vals=[155,135,129,107]
ax[0,0].plot(range(4),vals,'o-',color='#4C78A8',lw=2,ms=8); ax[0,0].set_xticks(range(4),stages); ax[0,0].set_ylabel('Genes'); ax[0,0].set_title('A  Signature attrition',loc='left',weight='bold'); ax[0,0].set_ylim(0,175)
for i,v in enumerate(vals): ax[0,0].text(i,v+5,str(v),ha='center')
labels=['LOO\n(no Sham1)','LOO\n(no Sham2)','LOO\n(no Sham3)','Ambient\nexclusion','Ortholog\nrestriction']; overlap=[.884,.895,.912,.974,.832]; concord=[1,1,1,.983,1]; y=np.arange(len(labels))
ax[0,1].scatter(overlap,y,color='#59A14F',label='Overlap/retention'); ax[0,1].scatter(concord,y,marker='D',color='#F28E2B',label='Directional concordance'); ax[0,1].set_yticks(y,labels); ax[0,1].set_xlim(.75,1.03); ax[0,1].set_xlabel('Fraction'); ax[0,1].set_title('B  Prespecified robustness',loc='left',weight='bold'); ax[0,1].legend(frameon=False,fontsize=8,loc='lower left')
ax[1,0].axis('off'); ax[1,0].text(.03,.92,'C  Frozen transfer object',fontsize=12,weight='bold'); ax[1,0].text(.06,.72,'107 genes',fontsize=20,weight='bold',color='#4C78A8'); ax[1,0].text(.06,.55,'Icos excluded\nDirection fixed before human validation\nNo post hoc selection or reweighting',fontsize=11,linespacing=1.6); ax[1,0].text(.06,.20,'Sham1 exclusion: 88.4% overlap\n100% directional concordance',fontsize=10,color='#555')
ax[1,1].axis('off'); ax[1,1].text(.03,.92,'D  Assay-integrity boundary',fontsize=12,weight='bold'); ax[1,1].text(.08,.68,'G0B-3',fontsize=22,weight='bold',color='#E45756'); ax[1,1].text(.08,.50,'Assay-informative failure',fontsize=12,weight='bold'); ax[1,1].text(.08,.32,'Retained as a non-interpretable\ncondition; not reclassified as biology.',fontsize=10,color='#555')
for a in ax.ravel(): a.spines['top'].set_visible(False); a.spines['right'].set_visible(False)
fig.savefig(BASE/'figures/figure1_signature_provenance.png',dpi=300,bbox_inches='tight')
