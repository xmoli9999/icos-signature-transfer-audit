#!/usr/bin/env python3
import os,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Rectangle
exec(open(os.environ["ICOS_BASE"]+"/manuscript/_fig_common.py").read())

fig=plt.figure(figsize=(3.35,2.2))
# --- 1C placeholder for paired flow ---
c=fig.add_subplot(111); c.axis("off")
c.add_patch(FancyBboxPatch((0.02,0.10),0.96,0.78,boxstyle="round,pad=0.015",
    facecolor="none",edgecolor=GRID,lw=1.0,linestyle=(0,(3,3))))
c.text(0.5,0.66,"Paired lymph vs blood",ha="center",fontsize=7.0,color=INK)
c.text(0.5,0.46,"group × compartment interaction\n$outcome \\sim group * compartment + (1|animal)$",
       ha="center",fontsize=7.0,color=INK2)
c.text(0.5,0.26,"not estimable: animal-level pairing unavailable\n(model specification frozen 2026-09-16; no estimate was generated)",ha="center",
       fontsize=7.0,color=INK2,style="italic")
c.set_xlim(0,1); c.set_ylim(0,1); ttl(c,"Prespecified paired lymph-versus-blood analysis")
save_fig(fig,"FigureS1"); plt.close(fig); print("FigureS1 ok")
