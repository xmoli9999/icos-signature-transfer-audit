#!/usr/bin/env python3
import os,json,numpy as np,pandas as pd,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
exec(open(os.environ["ICOS_BASE"]+"/manuscript/_fig_common.py").read())

# ================= Figure 3 =================
_r=pd.read_csv(f"{O}/figure3_depth_correlations.csv").set_index("quantity")["value"]
_SUB,_SEED=int(_r["plot_subsample_n"]),int(_r["plot_subsample_seed"])
mm=m.sample(n=min(_SUB,len(m)),random_state=_SEED)   # same analysis set the correlations use
fig,ax=plt.subplots(2,2,figsize=(W180,5.4)); ax=ax.ravel()
for k,(col,ptitle,rho) in enumerate([("ucell_sig","Raw UCell (cell)",float(_r["cell_raw"])),("delta_ucell","ΔUCell (cell)",float(_r["cell_calibrated"]))]):
    ax[k].scatter(mm.nFeature_RNA,mm[col],s=1.2,color=INK2,alpha=.08,linewidth=0,rasterized=True)
    ax[k].set_xlabel("detected genes"); ax[k].set_ylabel("raw UCell" if k==0 else "ΔUCell")
    ttl(ax[k],ptitle); tag(ax[k],"AB"[k])
    ax[k].text(.97,.05,f"ρ = {rho:.3f}",transform=ax[k].transAxes,ha="right",fontsize=7.0,color=INK)
# C donor level（两种读数各自 z 标准化后同轴比较，避免量纲压缩）
z=lambda v:(v-np.mean(v))/np.std(v,ddof=1)
zraw=z(P.mean_raw.values); zdel=z(P.mean_delta.values)
for g,c in [("NHC",NHC),("CI-NS",NS),("CI-Sep",SEP)]:
    k=(P.group==g).values
    ax[2].scatter(P.med_nfeat[k],zraw[k],s=14,color=c,alpha=.9,edgecolor="white",linewidth=.5,marker="o")
    ax[2].scatter(P.med_nfeat[k],zdel[k],s=14,color=c,alpha=.9,edgecolor="white",linewidth=.5,marker="^")
ax[2].set_xlabel("donor median detected genes"); ax[2].set_ylabel("donor mean score (z)")
ttl(ax[2],"Donor level"); tag(ax[2],"C")
ax[2].text(.03,.99,f"Overall, 47 donors:   ○ raw ρ = {float(_r['donor_raw']):.3f}   △ ΔUCell ρ = {float(_r['donor_calibrated']):.3f}\nMean within ICU group:   raw ρ = {float(_r['within_group_raw']):.3f}   ΔUCell ρ = {float(_r['within_group_calibrated']):.3f}",
           transform=ax[2].transAxes,va="top",fontsize=7.0,color=INK,zorder=6,
           bbox=dict(boxstyle="square,pad=0.25",facecolor="white",edgecolor="none",alpha=0.85))
ax[2].text(.5,-.34,"circles raw, triangles calibrated; grey healthy,\nblue ICU non-sepsis, orange ICU sepsis",
           transform=ax[2].transAxes,ha="center",fontsize=7.0,color=INK2)
# D matched-control diagnostics
diag=json.load(open(f"{O}/gate1_calibration_diagnostics.json"))
lab=["detection\nfrequency","mean\nexpression"]
sigv=[diag["sig_det_f_median"],diag["sig_mean_e_median"]]
ctlv=[diag["ctrl_det_f_median"],diag["ctrl_mean_e_median"]]
x=np.arange(2); w=.34
ax[3].bar(x-w/2-0.01,sigv,w,color=SEP,edgecolor="white",lw=.8,label="signature")
ax[3].bar(x+w/2+0.01,ctlv,w,color=NS,edgecolor="white",lw=.8,label="matched controls")
ax[3].set_xticks(x); ax[3].set_xticklabels(lab,fontsize=7.0); ax[3].set_ylabel("median")
ttl(ax[3],"Control matching"); tag(ax[3],"D")
ax[3].legend(frameon=False,fontsize=7.0)
fig.subplots_adjust(wspace=.28,hspace=.52,bottom=.12,top=.92,left=.095,right=.985)
save_fig(fig,"Figure3"); plt.close(fig); print("Figure3 ok")

# ================= Figure 4 =================
fig,ax=plt.subplots(2,2,figsize=(W180,5.4)); ax=ax.ravel()
for i,(g,c) in enumerate([("NHC",NHC),("CI-NS",NS),("CI-Sep",SEP)]):
    v=P.loc[P.group==g,"mean_delta"].values
    ax[0].scatter(jit(len(v),i),v,s=12,color=c,alpha=.88,edgecolor="white",linewidth=.5,zorder=3)
    ax[0].hlines(v.mean(),i-.2,i+.2,color=INK,lw=1.4,zorder=4)
ax[0].axhline(0,color=GRID,lw=.8,zorder=0)
ax[0].set_xticks(range(3)); ax[0].set_xticklabels(["healthy\n9","ICU\nnon-sepsis\n19","ICU\nsepsis\n19"],fontsize=7.0)
ax[0].set_ylabel("patient mean ΔUCell"); ttl(ax[0],"Primary readout"); tag(ax[0],"A")
rows=[("Calibrated UCell\n(primary)",res["primary"]["beta"],res["primary"]["ci_lo"],res["primary"]["ci_hi"]),
      ("Unadjusted",res["sens1_unadjusted_delta"]["diff"],res["sens1_unadjusted_delta"]["ci_lo"],res["sens1_unadjusted_delta"]["ci_hi"]),
      ("Raw UCell, adjusted",res["sens2_raw_adjusted"]["beta"],res["sens2_raw_adjusted"]["ci_lo"],res["sens2_raw_adjusted"]["ci_hi"]),
      ("AUCell, adjusted",res["sens3_aucell_adjusted"]["beta"],res["sens3_aucell_adjusted"]["ci_lo"],res["sens3_aucell_adjusted"]["ci_hi"])]
for j,(lab_,b,lo,hi) in enumerate(rows):
    y=len(rows)-1-j
    ax[1].plot([lo,hi],[y,y],color=INK2,lw=1.2,solid_capstyle="round")
    ax[1].scatter([b],[y],s=20,color=INK,zorder=3)
ax[1].axvline(0,color=GRID,lw=.8,zorder=0)
ax[1].set_yticks(range(len(rows))); ax[1].set_yticklabels([r[0] for r in rows[::-1]],fontsize=7.0)
ax[1].set_xlabel("effect (sepsis − non-sepsis)"); ttl(ax[1],"Primary and sensitivities"); tag(ax[1],"B")
# C funnel
stages=["Repositories searched","Closest supportive cohorts","Meeting all five criteria"]
vals=[8,2,0]
ys=np.arange(len(vals))[::-1]
for y,(s_,v) in zip(ys,zip(stages,vals)):
    if v:
        ax[2].barh(y,v,height=.55,color=NS,edgecolor="white",lw=.8)
    else:
        # a count of zero is drawn as a zero-length endpoint, not as a bar of
        # visible length, so the visual encoding matches the value exactly
        ax[2].plot([0,0],[y-.275,y+.275],color=SEP,lw=1.8,solid_capstyle="butt",zorder=3)
    ax[2].text(max(v,0)+0.45,y,f"{s_}  ({v})",va="center",fontsize=7.0,color=INK)
ax[2].set_xlim(0,17); ax[2].set_ylim(-.7,2.7); ax[2].set_yticks([]); ax[2].set_xticks([0,5,10])
ax[2].set_xlabel("count"); ax[2].spines["left"].set_visible(False)
ttl(ax[2],"Eligibility audit"); tag(ax[2],"C")
# D viability
coh=["Kwok","Reyes"]; sep_=[14/26,3/8]; ctl=[4/7,7/7]; nsep=["14/26","3/8"]; nctl=["4/7","7/7"]
x=np.arange(2); w=.34
ax[3].bar(x-w/2-.01,sep_,w,color=SEP,edgecolor="white",lw=.8,label="sepsis arm")
ax[3].bar(x+w/2+.01,ctl,w,color=NS,edgecolor="white",lw=.8,label="critical non-sepsis")
for j in range(2):
    ax[3].text(x[j]-w/2-.01,sep_[j]+.03,nsep[j],ha="center",fontsize=7.0,color=INK)
    ax[3].text(x[j]+w/2+.01,ctl[j]+.03,nctl[j],ha="center",fontsize=7.0,color=INK)
ax[3].set_xticks(x); ax[3].set_xticklabels(coh,fontsize=7.0)
ax[3].set_ylabel("fraction ≥30 cells"); ax[3].set_ylim(0,1.42)
ttl(ax[3],"Gate 2b viability"); tag(ax[3],"D")
ax[3].legend(frameon=False,fontsize=7.0,loc="upper center",ncol=2,columnspacing=.8,handlelength=1.1,handletextpad=.3,bbox_to_anchor=(.5,1.02))
ax[3].text(.5,-.46,"both arms required ≥5 qualifying patients;\n4 (Kwok) and 3 (Reyes) observed → stop",
           transform=ax[3].transAxes,ha="center",fontsize=7.0,color=INK2)
fig.subplots_adjust(wspace=.30,hspace=.55,bottom=.12,top=.92,left=.10,right=.985)
save_fig(fig,"Figure4"); plt.close(fig); print("Figure4 ok")
