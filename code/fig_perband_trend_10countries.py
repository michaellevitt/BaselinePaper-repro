#!/usr/bin/env python3
"""
Figure S2 (formerly S3), 10-country version — per-age-band pre-pandemic mortality trend (two-step) in
the 10 countries with the most deaths, per John's request (comments 93/156/942):
  * all 10 countries (USA, Japan, Germany, Italy, England/Wales, France, Spain, Poland, Korea, Canada);
  * a 4-band centred moving-average spline over the per-band point estimates, to cut band-to-band noise;
  * a SHARED y-axis range across countries (separately for slope and slope-of-slopes) for visual
    comparison;
  * NO I^2 in the panel titles (John, #169): I^2 scales with the precision of the estimates, so the
    elderly bands — most deaths, tightest CIs — score a high I^2 even where their absolute departure
    from the all-ages value is small. I^2 is still computed and printed to stdout for reference.

Same per-band construction as the 3-country Figure S3 (see fig_perband_trend_heterogeneity.py):
5-year centred moving log-slopes -> OLS on window centre -> slope-in-2019 (g19) and slope-of-slopes (q),
each with a 95% CI.

Reads : data/master_5x1_DPM_90plus.csv, output/slope_of_slopes_CI.csv
Writes: output/perband_trend_heterogeneity_10countries.csv
        output/fig_perband_trend_10countries.png
"""
import os, sys, csv, numpy as np
import scipy.stats as ss
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); OUTD=os.path.join(ROOT,"output")
sys.path.insert(0, HERE)
from excess_anchor_window import load, cell_for

FINE=["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
      "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
def astart(a): return 0 if a=="0" else (1 if a=="1-4" else int(a.split("-")[0].replace("+","")))
GE65=[a for a in FINE if astart(a)>=65]
RELIABLE={"BGR":2015,"CHL":2011,"HRV":2007,"EST":2005,"HUN":2006,"LVA":2007,"LTU":2006,"POL":2008,"SVK":2006}
TEN=[("USA","United States"),("JPN","Japan"),("DEUTNP","Germany"),("ITA","Italy"),
     ("GBRTENW","England & Wales"),("FRATNP","France"),("ESP","Spain"),("POL","Poland"),
     ("KOR","Korea"),("CAN","Canada")]
data,name=load()

def ols_ci(x,y,x_eval):
    x=np.asarray(x,float); y=np.asarray(y,float); n=len(x)
    if n<3: return None
    xbar=x.mean(); Sxx=((x-xbar)**2).sum()
    b=((x-xbar)*(y-y.mean())).sum()/Sxx; a=y.mean()-b*xbar
    f=a+b*x_eval; dof=n-2; s2=((y-(a+b*x))**2).sum()/dof; t=ss.t.ppf(0.975,dof)
    se_b=np.sqrt(s2/Sxx); se_f=np.sqrt(s2*(1.0/n+(x_eval-xbar)**2/Sxx))
    return dict(b=b,b_lo=b-t*se_b,b_hi=b+t*se_b,se_b=se_b,f=f,f_lo=f-t*se_f,f_hi=f+t*se_f,se_f=se_f)

def lnr(loc,a,y):
    c=cell_for(data,loc,y,a,"T")
    return None if (c is None or c[1]<=0 or c[0]<=0) else np.log(c[0]/c[1])
def mslope(loc,a,c):
    pts=[(y,lnr(loc,a,y)) for y in range(c-2,c+3)]; pts=[(y,v) for y,v in pts if v is not None]
    if len(pts)<3: return None
    x=np.array([p[0] for p in pts],float); yv=np.array([p[1] for p in pts]); xc=x-x.mean()
    return 100.0*(xc*(yv-yv.mean())).sum()/(xc**2).sum()
def perband(loc):
    rsd=RELIABLE.get(loc,2003); out={}
    for a in FINE:
        pre=[(c,mslope(loc,a,c)) for c in range(rsd+2,2018)]; pre=[(c,v) for c,v in pre if v is not None]
        out[a]=ols_ci([c for c,_ in pre],[v for _,v in pre],2019) if len(pre)>=3 else None
    return out
def het(est,se):
    est=np.asarray(est,float); w=1.0/np.asarray(se,float)**2
    m=np.isfinite(est)&np.isfinite(w); est,w=est[m],w[m]; k=len(est)
    if k<2: return np.nan
    tbar=(w*est).sum()/w.sum(); Q=(w*(est-tbar)**2).sum()
    return max(0.0,(Q-(k-1))/Q)*100.0
def mavg(vals, w=4):                                  # centred moving average over adjacent bands (window w)
    x=np.array([np.nan if v is None else v for v in vals],float); n=len(x); out=np.full(n,np.nan)
    half=w//2
    for i in range(n):
        lo=max(0,i-half); hi=min(n,i+half+ (w%2)); seg=x[lo:hi]; seg=seg[np.isfinite(seg)]
        if len(seg): out[i]=seg.mean()
    return out

# all-ages "paper" dashed value per country
paper={}
for r in csv.DictReader(open(os.path.join(OUTD,"slope_of_slopes_CI.csv"))):
    if r.get("country") and r.get("u_islope2019","") not in ("",None):
        paper[r["country"]]=(float(r["u_islope2019"]),float(r["u_sos"]))

res={l:perband(l) for l,_ in TEN}

# ---- CSV ----
with open(os.path.join(OUTD,"perband_trend_heterogeneity_10countries.csv"),"w",newline="") as f:
    w=csv.writer(f); w.writerow(["code","name","band","g19","g19_lo","g19_hi","se_g19","q","q_lo","q_hi","se_q"])
    for l,nm in TEN:
        for a in FINE:
            d=res[l][a]
            if d: w.writerow([l,nm,a,round(d["f"],4),round(d["f_lo"],4),round(d["f_hi"],4),round(d["se_f"],4),
                              round(d["b"],4),round(d["b_lo"],4),round(d["b_hi"],4),round(d["se_b"],4)])

# ---- shared y-limits from the 5th/95th percentile of CI endpoints (robust to a few huge young-band CIs) ----
def ylim(key,lo,hi):
    ends=[]
    for l,_ in TEN:
        for a in FINE:
            d=res[l][a]
            if d: ends+=[d[lo],d[hi]]
    ends=np.array(ends); lo_,hi_=np.nanpercentile(ends,2.5),np.nanpercentile(ends,97.5)
    m=1.08; c=(lo_+hi_)/2; h=(hi_-lo_)/2*m
    return c-h,c+h
YL={"f":ylim("f","f_lo","f_hi"),"b":ylim("b","b_lo","b_hi")}

# ---- figure: 5 rows x 4 cols; each country = adjacent (slope, slope-of-slopes) pair ----
plt.rcParams.update({"font.size":9})
fig,axes=plt.subplots(5,4,figsize=(17,19)); xpos=np.arange(len(FINE))
COL_LT="#3b6fb0"; COL_GE="#c0392b"
I2LOG=[]   # (country, parameter, I^2 all bands, I^2 65+) — reported below, not on the figure
PARAMS=[("slope in 2019","f","f_lo","f_hi","%/yr",0),("slope-of-slopes","b","b_lo","b_hi","%/yr²",1)]
for ci,(l,nm) in enumerate(TEN):
    d=res[l]; rrow=ci//2; cblk=ci%2
    for (lbl,key,klo,khi,unit,pidx) in PARAMS:
        ax=axes[rrow][cblk*2+pidx]
        pts=[d[a][key] if d[a] else None for a in FINE]
        for k,a in enumerate(FINE):
            dd=d[a]
            if not dd: continue
            c=COL_GE if a in GE65 else COL_LT
            ax.plot([k,k],[dd[klo],dd[khi]],color=c,lw=1.0,alpha=0.55,zorder=2)
            ax.plot(k,dd[key],"o",color=c,ms=3.4,zorder=3)
        ax.plot(xpos,mavg(pts,4),color="black",lw=1.6,zorder=4,label="4-band moving avg")   # spline
        pv=paper[l][pidx]; ax.axhline(pv,ls="--",color="0.35",lw=1.1,zorder=1)
        ax.axhline(0,color="0.7",lw=0.7,zorder=1); ax.axvline(13.5,ls=":",color="0.6",lw=0.9)
        ax.set_ylim(*YL[key])
        # I^2 is computed and printed for reference but deliberately NOT shown in the panel titles:
        # it scales with the precision of the underlying estimates, so the elderly bands (most deaths,
        # tightest CIs) score a high I^2 even where their absolute departure from the all-ages value is
        # small. The absolute spread around the dashed line is the informative comparison (John, #169).
        i2a=het([d[a][key] for a in FINE if d[a]],[d[a]["se_"+("f" if key=="f" else "b")] for a in FINE if d[a]])
        i2b=het([d[a][key] for a in GE65 if d[a]],[d[a]["se_"+("f" if key=="f" else "b")] for a in GE65 if d[a]])
        I2LOG.append((nm,lbl,i2a,i2b))
        ax.set_title(f"{nm} — {lbl}",fontsize=9.2,fontweight="bold")
        ax.set_xticks(xpos); ax.set_xticklabels(FINE,rotation=90,fontsize=5.6)
        if cblk*2+pidx==0 or cblk*2+pidx==2: ax.set_ylabel(f"({unit})",fontsize=8)
        if rrow==4: ax.set_xlabel("age band",fontsize=8)
h=[plt.Line2D([],[],marker="o",ls="",color=COL_LT,label="<65 age bands"),
   plt.Line2D([],[],marker="o",ls="",color=COL_GE,label="65+ age bands"),
   plt.Line2D([],[],color="black",lw=1.6,label="4-band moving average"),
   plt.Line2D([],[],ls="--",color="0.35",lw=1.1,label="all-ages value (applied to every band)")]
fig.legend(handles=h,loc="upper center",ncol=4,fontsize=10,bbox_to_anchor=(0.5,0.985),frameon=False)
fig.suptitle("Per–age–band pre-pandemic mortality trend (two-step) in the 10 countries with the most deaths\n"
             "points = per-band estimate · bars = 95% CI · black = 4-band moving average · shared vertical scale across countries · dotted = 65 cut-off",
             fontweight="bold",fontsize=12,y=0.999)
fig.tight_layout(rect=[0,0,1,0.965])
fig.savefig(os.path.join(OUTD,"fig_perband_trend_10countries.png"),dpi=140,bbox_inches="tight")
print("wrote fig_perband_trend_10countries.png + perband_trend_heterogeneity_10countries.csv")
print("\nI² (reference only; NOT shown on the figure — see note in code):")
for nm,lbl,a,b in I2LOG: print(f"   {nm:16} {lbl:16} all-bands {a:5.0f}%   65+ {b:5.0f}%")
print(f"shared y-lims: slope {YL['f'][0]:.1f}..{YL['f'][1]:.1f} %/yr ; slope-of-slopes {YL['b'][0]:.2f}..{YL['b'][1]:.2f} %/yr²")
