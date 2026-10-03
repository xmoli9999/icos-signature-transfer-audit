#!/usr/bin/env python3
import os,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Rectangle
exec(open(os.environ["ICOS_BASE"]+"/manuscript/_fig_common.py").read())

fig=plt.figure(figsize=(W180,2.9))
gs=fig.add_gridspec(1,2,width_ratios=[1.00,1.35],wspace=0.75,left=.055,right=.985,bottom=.20,top=.82)

# --- 1A attrition flow ---
a=fig.add_subplot(gs[0,0]); a.axis("off"); tag(a,"A")
steps=[("Memory-matched\ncandidates","155"),("One-to-one\northolog restriction","−48"),
       ("Frozen rat→human\nsignature","107")]
for j,(t,v) in enumerate(steps):
    y=2-j
    a.add_patch(FancyBboxPatch((0.0,y*0.315+0.045),0.88,0.20,boxstyle="round,pad=0.010",
        facecolor="#f2f1ee" if j<2 else NS,edgecolor="none"))
    a.text(0.045,y*0.315+0.145,t,fontsize=7.0,va="center",color="white" if j==2 else INK)
    a.text(0.855,y*0.315+0.145,v,fontsize=7.0,va="center",ha="right",fontweight="bold",
           color="white" if j==2 else INK)
    if j<2: a.annotate("",xy=(0.44,y*0.315+0.012),xytext=(0.44,y*0.315+0.042),
        arrowprops=dict(arrowstyle="-|>",color=INK2,lw=.9))
a.text(0.44,-0.17,"index gene $\\it{Icos}$ excluded\nbefore transfer",fontsize=7.0,ha="center",color=INK2)
a.set_xlim(0,1); a.set_ylim(-0.02,1.02); ttl(a,"Signature provenance")

# --- 1B robustness audit panel ---
b=fig.add_subplot(gs[0,1]); tag(b,"B")
rows=[("Ambient exclusion (retention)",None,0.974,0.60,True),
      ("Sign concordance (ambient)",None,0.983,None,True),
      ("Sham1-excluded overlap",None,0.884,None,True),
      ("Mean leave-one-rat-out\nheld-out AUROC",None,0.896,0.70,True),
      ("$\\it{Icos}$-detected stratification",None,np.nan,None,False)]
ys=np.arange(len(rows))[::-1]
for (lab,metric,v,thr,ok),y in zip(rows,ys):
    if np.isfinite(v):
        b.barh(y,v,height=.5,color=NS if ok else GRID,edgecolor="white",lw=.8)
        b.text(v+0.015,y,f"{v:.3f}",va="center",fontsize=7.0,color=INK)
        if thr: b.plot([thr,thr],[y-.3,y+.3],color=SEP,lw=1.1)
    else:
        b.add_patch(Rectangle((0,y-.25),1.0,.5,facecolor="#f2f1ee",edgecolor="none",hatch="////"))
        b.text(0.5,y,"assay-informative failure",ha="center",va="center",fontsize=7.0,color=INK2)
b.set_yticks(ys); b.set_yticklabels([r[0] for r in rows],fontsize=7.0)
b.set_xlim(0,1.12); b.set_xlabel("value"); ttl(b,"Prespecified robustness checks")
b.text(.99,-.32,"orange tick = prespecified threshold",transform=b.transAxes,ha="right",fontsize=7.0,color=INK2)

save_fig(fig,"Figure1"); plt.close(fig); print("Figure1 ok")
