#!/usr/bin/env python3
"""
Figure 1 (John #10): three panels — 2020-2025, 2020-2023 (pandemic), 2024-2025 (post) — of per-country
all-ages excess P-score under Fa, TTa, STTa (John's three preferred models), countries ordered by the mean
of the three in 2020-2025.  Reads output/model_option1_periods.csv.  Writes output/fig_option1_3panel.png .
"""
import os, csv, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUTD=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"output")
MODELS=["Fa","TTa","STTa","STTa+"]; COL={"Fa":"#7f8c8d","TTa":"#159090","STTa":"#c2185b","STTa+":"#6a3d9a"}   # STTa+ = slope-anchored STTa (purple, consistent w/ montages)
ORDERKEY=["Fa","TTa","STTa","STTa+"]   # order by all-ages 2020-25 P-score averaged over all four models (Michael, 2026-07-29)
PERIODS=[("2020-2025","2020–2025 (all)"),("2020-2023","2020–2023 (pandemic)"),("2024-2025","2024–2025 (post)")]
R=[r for r in csv.reader(open(os.path.join(OUTD,"model_option1_periods.csv")))]
hdr=R[0]; rows=[dict(zip(hdr,r)) for r in R[1:]]
def P(country,model,period):
    # pooled P% = 100*(sumO-sumE)/sumE, the same statistic every table in the paper
    # reports (Tables 2, 3, S2, S3, S5); the mean-annual P-score differs slightly.
    for r in rows:
        if r["country"]==country and r["group"]=="All" and r["model"]==model and r["period"]==period:
            return 100.0*(float(r["O"])-float(r["E"]))/float(r["E"])
    return np.nan
# name map
import sys; sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__))))
from excess_anchor_window import load
_,name=load()
countries=sorted({r["country"] for r in rows})
key={c:np.mean([P(c,m,"2020-2025") for m in ORDERKEY]) for c in countries}
order=sorted(countries,key=lambda c:key[c])           # ascending -> highest at top
yv={c:i for i,c in enumerate(order)}
fig,axes=plt.subplots(1,3,figsize=(13.5,11.0),sharey=True)
for ax,(pk,pt) in zip(axes,PERIODS):
    ax.axvline(0,color="#888",lw=.9,zorder=1)
    for c in order:
        for m in MODELS:
            ax.scatter(P(c,m,pk),yv[c],s=34,color=COL[m],edgecolor="white",linewidth=.5,zorder=3)
    ax.set_title(pt,fontsize=12,fontweight="bold")
    ax.grid(axis="x",alpha=.25); ax.set_axisbelow(True); ax.set_xlabel("excess P-score (%), per year",fontsize=10.5)
    for s in ["top","right"]: ax.spines[s].set_visible(False)
axes[0].set_yticks(range(len(order))); axes[0].set_yticklabels([name.get(c,c) for c in order],fontsize=7.5)
axes[0].set_ylim(-0.7,len(order)-0.3)
hp=[plt.Line2D([0],[0],marker="o",color="w",markerfacecolor=COL[m],markersize=9,label=m) for m in MODELS]
axes[2].legend(handles=hp,loc="lower right",fontsize=10,frameon=True,title="baseline")
# heading removed per John (C18) — descriptive text now lives in the figure legend
fig.tight_layout()
FIG=os.path.join(OUTD,"fig_option1_3panel_4model.png"); fig.savefig(FIG,dpi=150,bbox_inches="tight"); plt.close(fig)
print("wrote",FIG,f"({len(order)} populations)")
