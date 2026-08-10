#!/usr/bin/env python3
"""
Figure 3 (9-panel form, Michael 2026-07-29): wealth gradient as 3 indicators (rows: per-capita GDP,
Gini, % in poverty) x 3 periods (cols: 2020-2025 / 2020-2023 / 2024-2025).  Each panel scatters the
per-population all-ages excess P-score AVERAGED over the four baselines (Fa, TTa, STTa, STTa+) against
the indicator, coloured by vulnerability group, with Pearson r + OLS fit.  Periods as columns to parallel
Figure 2.  Reads output/model_option1_periods.csv + data/vuln_covariates_38_v3.csv.
Writes output/fig_option1_wealth_9panel_indperiod.png .
"""
import os, csv, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUTD=os.path.join(ROOT,"output"); DATA=os.path.join(ROOT,"data")
MODELS=["Fa","TTa","STTa","STTa+"]
IND=[("gdp_pc_2021_usd","per-capita GDP 2021 (USD, ×1000)",1e-3),
     ("gini","Gini (income inequality)",1.0),
     ("poverty_pct","population in poverty (%)",1.0)]
PERIODS=[("2020-2025","2020–2025 (all)"),("2020-2023","2020–2023 (pandemic)"),("2024-2025","2024–2025 (post)")]
COV={r["code"]:r for r in csv.DictReader(open(os.path.join(DATA,"vuln_covariates_38_v3.csv")))}
per=list(csv.DictReader(open(os.path.join(OUTD,"model_option1_periods.csv"))))
GCOL={"more":"#c0392b","less":"#2e7d32"}
def Pmodel(c,m,period):
    r=next((x for x in per if x["country"]==c and x["group"]=="All" and x["model"]==m and x["period"]==period),None)
    return 100.0*(float(r["O"])-float(r["E"]))/float(r["E"]) if r else None
def Pmean(c,period):                                        # 4-model mean of the pooled P%
    vals=[Pmodel(c,m,period) for m in MODELS]
    return np.mean(vals) if all(v is not None for v in vals) else None

fig,axes=plt.subplots(3,3,figsize=(13.5,13.2))
rtab={}
for i,(key,ilab,sc) in enumerate(IND):
    for j,(pk,pt) in enumerate(PERIODS):
        ax=axes[i][j]; xs=[];ys=[];cs=[]
        for c,cov in COV.items():
            if cov.get(key,"")=="" : continue
            y=Pmean(c,pk)
            if y is None: continue
            xs.append(float(cov[key])*sc); ys.append(y); cs.append(GCOL[cov["group"]])
        xs=np.array(xs); ys=np.array(ys)
        ax.axhline(0,color="#bbb",lw=.7)
        ax.scatter(xs,ys,c=cs,s=34,edgecolor="white",linewidth=.5,zorder=3)
        pr=stats.pearsonr(xs,ys)[0]; sl,ic,_,_,_=stats.linregress(xs,ys); rtab[(key,pk)]=pr
        xf=np.linspace(xs.min(),xs.max(),50); ax.plot(xf,ic+sl*xf,color="#333",lw=1.5,ls="--",zorder=2)
        ax.text(0.04,0.94,f"r = {pr:+.2f}",transform=ax.transAxes,va="top",fontsize=11,fontweight="bold",
                bbox=dict(boxstyle="round",fc="white",ec="#ccc",alpha=.85))
        if i==0: ax.set_title(pt,fontsize=12,fontweight="bold")
        if j==0: ax.set_ylabel("4-model-mean excess P (%)",fontsize=9.5,fontweight="bold")
        ax.set_xlabel(ilab,fontsize=10)   # indicator label below each row's horizontal axis (John C19)
        ax.grid(alpha=.2); ax.set_axisbelow(True)
        for s in ["top","right"]: ax.spines[s].set_visible(False)
hp=[plt.Line2D([0],[0],marker="o",color="w",markerfacecolor=GCOL[g],markersize=9,label=f"{g} vulnerable") for g in ("more","less")]
axes[0][2].legend(handles=hp,loc="upper right",fontsize=9,frameon=True)
# heading removed per John (C19) — descriptive text now lives in the figure legend
fig.tight_layout()
FIG=os.path.join(OUTD,"fig_option1_wealth_9panel_indperiod.png"); fig.savefig(FIG,dpi=150,bbox_inches="tight"); plt.close(fig)
print("wrote",FIG)
print("\nPearson r (4-model-mean excess vs indicator):")
print(f"  {'indicator':14}"+"".join(f"{pt:>16}" for _,pt in PERIODS))
for key,ilab,_ in IND:
    print(f"  {key:14}"+"".join(f"{rtab[(key,pk)]:>16.2f}" for pk,_ in PERIODS))
