#!/usr/bin/env python3
"""
STMF-stratified moving-slope trend-of-trends for John (email 7 Jul 2026), ESP2013-standardized.
Six groups: All-ages + the 5 STMF bands (0-14, 15-64, 65-74, 75-84, 85+). For each (location, group):
  ESP-standardized rate R(y) = 1000 * sum_{b in group} w_b^ESP*(D_b/P_b) / sum w_b^ESP
  moving 5-yr slopes (%/yr of window mean), end years 2004..2019 = 16 points
  + anchor slope 2019->2024 (2-point, skips pandemic) = +1 point
  trend-of-trends = OLS of the 16 (no anchor) and the 16+1 (with anchor).
Writes output/stmf_moving_slopes.json and output/fig_moving_slope_stmf.png .
"""
import os, sys, json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from excess_anchor_window import load, cell_for

OUTD=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"output")
W=int(sys.argv[1]) if len(sys.argv)>1 else 5           # moving-window length (years); default 5
SFX="" if W==5 else f"_w{W}"
FINE=["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64",
      "65-69","70-74","75-79","80-84","85-89","90+"]
GROUPS={
 "All ages":  FINE,
 "0-14":      ["0","1-4","5-9","10-14"],
 "15-64":     ["15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64"],
 "65-74":     ["65-69","70-74"],
 "75-84":     ["75-79","80-84"],
 "85+":       ["85-89","90+"],
}
ORDER=list(GROUPS)
ESP={"0":1000,"1-4":4000,"5-9":5500,"10-14":5500,"15-19":5500,"20-24":6000,"25-29":6000,
     "30-34":6500,"35-39":7000,"40-44":7000,"45-49":7000,"50-54":7000,"55-59":6500,"60-64":6000,
     "65-69":5500,"70-74":5000,"75-79":4000,"80-84":2500,"85-89":1500,"90+":1000}
ENDS=list(range(2000+(W-1),2020)); ANCHOR_END=2025    # window-ends; anchor spans 2019->2024->2025
data,name=load()
# HMD/parent code -> 2025-build child code carrying the new 2025 row (IRL has no 2025 source)
BUILD25={"TWN":"TWN_MOI","JPN":"JPN_ESTAT","HKG":"HKG_CSD","KOR":"KOR_KOSTAT","CHL":"CHL_DEIS",
 "USA":"USA_WONDER","DNK":"DNK","GBRTENW":"GBRTENW_ONS","DEUTNP":"DEU_EUROSTAT","ITA":"ITA_EUROSTAT",
 "ESP":"ESP_EUROSTAT","FRATNP":"FRA_EUROSTAT","POL":"POL_EUROSTAT","NLD":"NLD_EUROSTAT","AUT":"AUT_EUROSTAT",
 "BGR":"BGR_EUROSTAT","CZE":"CZE_EUROSTAT","HUN":"HUN_EUROSTAT","GRC":"GRC_EUROSTAT","SVN":"SVN_EUROSTAT",
 "BEL":"BEL_EUROSTAT","CHE":"CHE_EUROSTAT","EST":"EST_EUROSTAT","FIN":"FIN_EUROSTAT","HRV":"HRV_EUROSTAT",
 "LTU":"LTU_EUROSTAT","LUX":"LUX_EUROSTAT","LVA":"LVA_EUROSTAT","NOR":"NOR_EUROSTAT","PRT":"PRT_EUROSTAT",
 "SVK":"SVK_EUROSTAT","SWE":"SWE_EUROSTAT"}
def std_rate25(loc,bands):
    code=BUILD25.get(loc)
    if not code: return None
    try: return 1000*sum(ESP[b]*data[code][(2025,b)]["T"][0]/data[code][(2025,b)]["T"][1] for b in bands)/sum(ESP[b] for b in bands)
    except KeyError: return None
REMAP={"GBR":"GBRTENW","DEU":"DEUTNP","FRA":"FRATNP"}
# The analysis population set. Originally read from a legacy artifact (asedx_all.json); the repo
# ships the canonical 38-population list as data/canonical_populations.csv so this is self-contained.
# (Order is irrelevant downstream — the engine re-sorts — but this matches the paper's set exactly.)
_CANON=os.path.join(os.path.dirname(OUTD),"data","canonical_populations.csv")
try:
    labs=[REMAP.get(r[0],r[0]) for r in json.load(open(os.path.join(OUTD,"asedx_all.json")))]
except FileNotFoundError:
    import csv as _csv
    labs=[row["code"] for row in _csv.DictReader(open(_CANON))]
DISP={"KOR":"South Korea","USA":"United States","GBRTENW":"England & Wales","DEUTNP":"Germany","FRATNP":"France"}
def dnm(loc): return DISP.get(loc, "Luxembourg" if name.get(loc,loc)=="Luxemburg" else name.get(loc,loc))
def cD(loc,y,b):
    c=cell_for(data,loc,y,b,"T")
    if c is None: raise KeyError((loc,y,b))
    return c
def std_rate(loc,y,bands): return 1000*sum(ESP[b]*cD(loc,y,b)[0]/cD(loc,y,b)[1] for b in bands)/sum(ESP[b] for b in bands)
def slope_pct(yrs,r):
    x=np.array(yrs,float); r=np.array(r); xc=x-x.mean()
    b=(xc*(r-r.mean())).sum()/(xc**2).sum(); return 100*b/r.mean()
def ols(x,y):
    x=np.array(x,float); y=np.array(y); xc=x-x.mean(); b=(xc*(y-y.mean())).sum()/(xc**2).sum(); return float(b)
def moving(loc,bands):
    out={}
    for ye in ENDS:
        ys=list(range(ye-(W-1),ye+1))
        try: out[ye]=slope_pct(ys,[std_rate(loc,y,bands) for y in ys])
        except KeyError: pass
    try:                                              # anchor rate-slope through 2019, 2024 (+2025 if built)
        ay=[2019,2024]; ar=[std_rate(loc,2019,bands),std_rate(loc,2024,bands)]
        r25=std_rate25(loc,bands)
        if r25 is not None: ay.append(2025); ar.append(r25)
        out[ANCHOR_END]=slope_pct(ay,ar)
    except KeyError: pass
    return out
def analyse(bands):
    rows=[]; coh={ye:[] for ye in ENDS+[ANCHOR_END]}
    for loc in labs:
        mv=moving(loc,bands)
        pre=[(ye,mv[ye]) for ye in ENDS if ye in mv]
        if len(pre)<4: continue
        for ye,v in mv.items(): coh[ye].append(v)
        tot_pre=ols([p[0] for p in pre],[p[1] for p in pre])
        anch=pre+([(ANCHOR_END,mv[ANCHOR_END])] if ANCHOR_END in mv else [])
        tot_anch=ols([p[0] for p in anch],[p[1] for p in anch])
        rows.append(dict(loc=loc,nm=dnm(loc),mv=mv,tot_pre=tot_pre,tot_anch=tot_anch,npts=len(pre)))
    cohmean={ye:(float(np.mean(coh[ye])) if coh[ye] else None) for ye in coh}
    return rows,cohmean

RES={g:analyse(GROUPS[g]) for g in ORDER}

# ---- John's non-overlapping 5-year BLOCKS S1..S5 ----
BLOCKS=[(2000,2004),(2005,2009),(2010,2014),(2015,2019)]      # S1..S4; S5 = 2019->2024 anchor
BLOCK_LABELS=["2000-04","2005-09","2010-14","2015-19","2019-24"]
def block_slopes(loc,bands):
    S=[]
    for (a,b) in BLOCKS:
        ys=list(range(a,b+1))
        try: S.append(slope_pct(ys,[std_rate(loc,y,bands) for y in ys]))
        except KeyError: S.append(None)
    try: S.append(slope_pct([2019,2024],[std_rate(loc,2019,bands),std_rate(loc,2024,bands)]))
    except KeyError: S.append(None)
    return S
def block_analyse(bands):
    rows=[]; coh=[[] for _ in range(5)]
    for loc in labs:
        S=block_slopes(loc,bands)
        if sum(v is not None for v in S[:4])<3: continue
        for i in range(5):
            if S[i] is not None: coh[i].append(S[i])
        pre=[(i+1,S[i]) for i in range(4) if S[i] is not None]; tot4=ols([p[0] for p in pre],[p[1] for p in pre])
        allp=[(i+1,S[i]) for i in range(5) if S[i] is not None]; tot5=ols([p[0] for p in allp],[p[1] for p in allp])
        rows.append(dict(loc=loc,nm=dnm(loc),S=S,tot4=tot4,tot5=tot5))
    cohmean=[float(np.mean(coh[i])) if coh[i] else None for i in range(5)]
    return rows,cohmean
BRES={g:block_analyse(GROUPS[g]) for g in ORDER}

# ---- JSON ----
digest={"groups":ORDER,"end_years":ENDS,"anchor_end":ANCHOR_END,
        "data":{g:{"cohort_mean":{str(k):(round(v,4) if v is not None else None) for k,v in RES[g][1].items()},
                   "rows":[{"loc":r["loc"],"name":r["nm"],"npts":r["npts"],
                            "tot_no_anchor":round(r["tot_pre"],4),"tot_with_anchor":round(r["tot_anch"],4),
                            "slopes":{str(ye):round(v,4) for ye,v in r["mv"].items()}} for r in RES[g][0]]}
               for g in ORDER},
        "block_labels":BLOCK_LABELS,
        "blocks":{g:{"cohort_mean":[round(v,4) if v is not None else None for v in BRES[g][1]],
                     "rows":[{"loc":r["loc"],"name":r["nm"],
                              "S":[round(v,4) if v is not None else None for v in r["S"]],
                              "tot_no_anchor":round(r["tot4"],4),"tot_with_anchor":round(r["tot5"],4)} for r in BRES[g][0]]}
                  for g in ORDER}}
json.dump(digest,open(os.path.join(OUTD,f"stmf_moving_slopes{SFX}.json"),"w"),indent=1)

# ---- figure: 6 panels ----
HL={"USA":"#2ca02c","KOR":"#d62728","JPN":"#1f77b4","NLD":"#7570b3"}
fig,axes=plt.subplots(2,3,figsize=(17,9.4)); axes=axes.ravel()
for k,g in enumerate(ORDER):
    ax=axes[k]; rows,cohmean=RES[g]; byloc={r["loc"]:r for r in rows}
    ax.axvspan(2019.5,2024.5,color="#fbeaea",zorder=0); ax.axhline(0,color="black",lw=0.9,zorder=1)
    for r in rows:
        if r["loc"] in HL: continue
        xs=sorted(r["mv"]); ax.plot(xs,[r["mv"][x] for x in xs],color="#c6c6c6",lw=0.7,alpha=0.5,zorder=2)
    for loc,col in HL.items():
        if loc in byloc:
            r=byloc[loc]; xs=sorted(r["mv"]); ax.plot(xs,[r["mv"][x] for x in xs],color=col,lw=1.8,zorder=4)
    xe=[ye for ye in ENDS if cohmean[ye] is not None]; ym=[cohmean[ye] for ye in xe]
    ax.plot(xe,ym,color="black",lw=2.8,zorder=5)
    if cohmean[ANCHOR_END] is not None: ax.plot(ANCHOR_END,cohmean[ANCHOR_END],"o",color="black",ms=6,zorder=6)
    b=ols(xe,ym); a=np.mean(ym)-b*np.mean(xe); xr=np.array([2004,2019]); ax.plot(xr,a+b*xr,color="#e6a000",lw=2.2,ls="--",zorder=6)
    ndec=sum(1 for r in rows if r["tot_pre"]>0)
    ax.text(0.03,0.04,f"cohort ToT = {b:+.3f} %/yr per yr\ndecelerating {ndec}/{len(rows)}",transform=ax.transAxes,
            fontsize=8.5,va="bottom",bbox=dict(boxstyle="round",fc="white",ec="#ccc",alpha=0.9))
    ax.set_title(f"{g}",loc="left",fontsize=12,fontweight="bold")
    ax.set_xlim(2003.5,2025.4); ax.set_xticks(range(2004,2025,4)); ax.tick_params(labelsize=8); ax.grid(alpha=0.3); ax.set_axisbelow(True)
handles=[plt.Line2D([],[],color="black",lw=2.8,label="cohort mean"),
         plt.Line2D([],[],color="#e6a000",lw=2.2,ls="--",label="OLS fit (pre-pandemic)"),
         plt.Line2D([],[],marker="o",color="black",lw=0,label="anchor 2019/2024/2025")]+\
        [plt.Line2D([],[],color=c,lw=1.8,label=dnm(l)) for l,c in HL.items()]
fig.legend(handles=handles,loc="lower center",ncol=7,fontsize=9,frameon=False,bbox_to_anchor=(0.5,-0.01))
fig.supxlabel(f"window end year  ({W}-yr window; {2000+(W-1)} = 2000–{str(2000+(W-1))[2:]})",fontsize=11)
fig.supylabel(f"moving {W}-year slope of the ESP2013-standardised rate  (%/yr)",fontsize=11)
fig.suptitle(f"Moving {W}-year slope regressed over time (trend-of-trends), by STMF age stratum — ESP2013 age-standardised, per country\n"
             "the 15–64 stratum shows the working-age reversal (rising) that the All-ages summary blends away",
             fontsize=12,fontweight="bold",y=1.0)
fig.tight_layout(rect=[0.01,0.02,1,0.965])
fig.savefig(os.path.join(OUTD,f"fig_moving_slope_stmf{SFX}.png"),dpi=145,bbox_inches="tight")

print(f"wrote output/stmf_moving_slopes{SFX}.json and fig_moving_slope_stmf{SFX}.png  (W={W})")
for g in ORDER:
    rows,cohmean=RES[g]; b=ols([y for y in ENDS if cohmean[y] is not None],[cohmean[y] for y in ENDS if cohmean[y] is not None])
    md=np.median([r["mv"].get(2019,np.nan) for r in rows])
    print(f"  {g:9s} n={len(rows):2d}  cohort ToT(no anchor) {b:+.3f}  decel {sum(r['tot_pre']>0 for r in rows)}/{len(rows)}")
