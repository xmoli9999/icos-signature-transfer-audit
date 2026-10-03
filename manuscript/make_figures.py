#!/usr/bin/env python3
import os,json,numpy as np,pandas as pd,matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
exec(open(os.environ["ICOS_BASE"]+"/manuscript/_fig_common.py").read())

# ---------- 数据 ----------
# ================= Figure 2 =================
fig,ax=plt.subplots(2,2,figsize=(W180,5.5)); ax=ax.ravel()
# 2A forest
rows=[("Pooled (adjusted)",res["primary"]),("CM stratum",res["decomp_CD4 T CM"]),
      ("Eff/EM stratum",res["decomp_CD4 T Eff/EM"])]
for j,(lab,r) in enumerate(rows):
    y=len(rows)-1-j
    ax[0].plot([r["ci_lo"],r["ci_hi"]],[y,y],color=INK2,lw=1.2,solid_capstyle="round")
    ax[0].scatter([r["beta"]],[y],s=22,color=INK,zorder=3)
ax[0].axvline(0,color=GRID,lw=.8,zorder=0)
ax[0].set_yticks(range(len(rows))); ax[0].set_yticklabels([r[0] for r in rows[::-1]])
ax[0].set_xlabel("β (sepsis − non-sepsis)"); ttl(ax[0],"Estimates by stratum"); tag(ax[0],"A")
# 2B state levels
for i,g in enumerate(["CI-NS","CI-Sep"]):
    s=comp[comp.group==g]
    dots(ax[1],s.mCM.values,i*2,NS if g=="CI-NS" else SEP)
    dots(ax[1],s.mEM.values,i*2+1,NS if g=="CI-NS" else SEP)
ax[1].axhline(0,color=GRID,lw=.8,zorder=0)
ax[1].set_xticks([0,1,2,3]); ax[1].set_xticklabels(["CM","Eff/EM","CM","Eff/EM"],fontsize=7.0)
ax[1].set_ylabel("ΔUCell"); ttl(ax[1],"State-specific level"); tag(ax[1],"B")
for xx,tt in [(0.25,"non-sepsis"),(0.75,"sepsis")]:
    ax[1].text(xx,-0.26,tt,transform=ax[1].transAxes,ha="center",fontsize=7.0,color=INK2)
# 2C Eff/EM fraction
for i,g in enumerate(["CI-NS","CI-Sep"]):
    dots(ax[2],comp.loc[comp.group==g,"f"].values,i,NS if g=="CI-NS" else SEP)
ax[2].set_xticks([0,1]); ax[2].set_xticklabels(["ICU\nnon-sepsis","ICU\nsepsis"],fontsize=7.0)
ax[2].set_ylabel("Eff/EM fraction"); ttl(ax[2],"State composition"); tag(ax[2],"C")
ax[2].text(.5,.04,"0.607 → 0.420",transform=ax[2].transAxes,ha="center",fontsize=7.0,color=INK2)
# 2D waterfall
_d=pd.read_csv(f"{O}/figure2_formal_decomposition.csv").set_index("component")["estimate"]
vals=[float(_d[k]) for k in ["within_state_at_NS_weights","composition_at_NS_levels",
      "interaction","between_patient_covariance"]]
_tot=float(_d["observed_patient_weighted_total"])
labs=["Within-state","Composition","Interaction","Covariance"]
base=0
for j,(v,l) in enumerate(zip(vals,labs)):
    ax[3].add_patch(Rectangle((j-.3,min(base,base+v)),.6,max(abs(v),4e-6),
        facecolor=SEP if v>0 else NS,edgecolor="white",lw=.8))
    if j: ax[3].plot([j-1+.3,j-.3],[base,base],color=GRID,lw=.7,zorder=0)
    base+=v
ax[3].plot([3+.3,4-.3],[base,base],color=GRID,lw=.7,zorder=0)
_cum=0
for j,v in enumerate(vals):
    if v>0: ax[3].text(j,_cum+v+0.00022,f"{v:+.5f}",ha="center",va="bottom",fontsize=7.0,color=INK2)
    else:   ax[3].text(j,_cum+v-0.00022,f"{v:+.5f}",ha="center",va="top",fontsize=7.0,color=INK2)
    _cum+=v
ax[3].text(4,_tot+0.00022,f"{_tot:+.5f}",ha="center",va="bottom",fontsize=7.0,color=INK2)
ax[3].add_patch(Rectangle((4-.3,0),.6,max(base,4e-6),facecolor="#3a3a38",edgecolor="white",lw=.8))
ax[3].axhline(0,color=GRID,lw=.8)
ax[3].set_ylim(-0.0022,0.0056)
ax[3].set_yticks([-0.002,0,0.002,0.004])
ax[3].set_xticks(range(5))
ax[3].set_xticklabels(labs+["Observed total"],fontsize=7.0,rotation=28,ha="right",rotation_mode="anchor")
ax[3].set_ylabel("ΔUCell contribution",labelpad=1); ttl(ax[3],"Decomposition"); tag(ax[3],"D")
ax[3].text(.5,-.46,f"within-state, composition, interaction and between-patient\ncovariance sum exactly to the observed {_tot:+.5f}",transform=ax[3].transAxes,ha="center",fontsize=7.0,color=INK2)
fig.subplots_adjust(wspace=.30,hspace=.55,bottom=.13,top=.92,left=.105,right=.985); save_fig(fig,"Figure2"); plt.close(fig)
print("Figure2 ok")
