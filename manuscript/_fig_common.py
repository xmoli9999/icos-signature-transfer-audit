import os,json,numpy as np,pandas as pd,matplotlib
import matplotlib.pyplot as plt
B=os.environ["ICOS_BASE"]; O=f"{B}/analysis_results"; FIG=f"{B}/figures"
NS,SEP,NHC="#2a78d6","#eb6834","#8a8a85"
INK,INK2,GRID="#0b0b0b","#52514e","#d9d8d4"
plt.rcParams.update({"font.size":7.5,"axes.linewidth":.6,"xtick.major.width":.6,
    "ytick.major.width":.6,"axes.edgecolor":INK2,"text.color":INK,
    "axes.labelcolor":INK,"xtick.color":INK2,"ytick.color":INK2,"figure.dpi":300,
    "savefig.dpi":400,"axes.spines.top":False,"axes.spines.right":False})
def tag(ax,s):
    ax.text(-0.26,1.10,s,transform=ax.transAxes,fontsize=9.5,fontweight="bold",va="bottom",ha="left")
def ttl(ax,s):
    ax.set_title(s,fontsize=8.5,loc="left",pad=5)
P=pd.read_csv(f"{O}/gate1_patient_level_unblinded.csv")
res=json.load(open(f"{O}/gate1_gate_test_results.json"))
cal=pd.read_csv(f"{O}/gate1_cell_scores_calibrated.csv.gz",
                usecols=["cell_index","patient","celltype_fine","nFeature_RNA","ucell_sig","delta_ucell"])
rawc=pd.read_csv(f"{O}/gate1_cell_scores.csv.gz",usecols=["cell_index","group","pop_primary_memory"])
d=cal.merge(rawc,on="cell_index"); m=d[d.pop_primary_memory.astype(bool)]
m2=m[m.group.isin(["CI-Sep","CI-NS"])]
rng=np.random.default_rng(7)
def jit(n,i,w=.085): return np.full(n,i)+rng.uniform(-w,w,n)

comp=m2.assign(e=(m2.celltype_fine=="CD4 T Eff/EM")).groupby(["patient","group"]).apply(
    lambda g: pd.Series({"f":g.e.mean(),"mCM":g.loc[~g.e,"delta_ucell"].mean(),
                         "mEM":g.loc[g.e,"delta_ucell"].mean()}),include_groups=False).reset_index()
def dots(ax,vals,i,c,lab=None):
    ax.scatter(jit(len(vals),i),vals,s=11,color=c,alpha=.85,linewidth=.5,edgecolor="white",zorder=3,label=lab)
    ax.hlines(np.mean(vals),i-.2,i+.2,color=INK,lw=1.4,zorder=4)

W180=7.09          # 180 mm double-column width, final print size
def save_fig(fig,name):
    for ext in ("pdf","svg"):
        fig.savefig(f"{FIG}/{name}.{ext}",bbox_inches="tight")
    fig.savefig(f"{FIG}/{name}.png",dpi=600,bbox_inches="tight")
