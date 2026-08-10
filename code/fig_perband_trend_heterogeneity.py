#!/usr/bin/env python3
"""
Figure S3 + Table S7 — per-age-band pre-pandemic mortality trend (two-step) in the three
largest populations (USA, Japan, Germany; each >80M).

For each country and each 5-year age band we repeat the paper's TWO-STEP trend construction on
that band's OWN log death rate (instead of the all-ages ESP-standardised rate):
  * 5-year moving log-slopes (%/yr), each indexed at its window CENTRE c, for centres
    c = reliable_start+2 .. 2017  (full 5-yr windows inside [reliable_start, 2019]);
  * OLS of those moving slopes on centre-year, giving
        slope-in-2019 (g19)     = fitted value at 2019   [%/yr]
        slope-of-slopes (q)      = OLS slope              [%/yr^2]
    each with a 95% CI and standard error.

This exposes how heterogeneous the per-band trends are relative to the single all-ages value the
main analysis applies to every band. We quantify it with a fixed-effect (inverse-variance) Q /
I-squared across the 20 bands and across the 65+ subset.

Reads : data/master_5x1_DPM_90plus.csv (via excess_anchor_window.load / cell_for)
        output/slope_of_slopes_CI.csv  (the all-ages "paper" g19/q per country, for the dashed line)
Writes: output/perband_trend_heterogeneity_bigpops.csv   (per country x band: g19,q + CIs + SEs)
        output/fig_perband_trend_heterogeneity_bigpops.png (Figure S3, 3x2 portrait)
        docs/Table_S7_perband_heterogeneity.xlsx           (Table S7: Q, I^2, p per country x param)
"""
import os, sys, csv, numpy as np
import scipy.stats as ss
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
OUTD=os.path.join(ROOT,"output"); DOCS=os.path.join(ROOT,"docs")
sys.path.insert(0, HERE)
from excess_anchor_window import load, cell_for

FINE=["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
      "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
def astart(a): return 0 if a=="0" else (1 if a=="1-4" else int(a.split("-")[0].replace("+","")))
GE65=[a for a in FINE if astart(a)>=65]
RELIABLE={"BGR":2015,"CHL":2011,"HRV":2007,"EST":2005,"HUN":2006,"LVA":2007,"LTU":2006,"POL":2008,"SVK":2006}
# the three largest populations, HMD/analysis code -> display name
LOCS=[("USA","United States"),("JPN","Japan"),("DEUTNP","Germany")]

data,name=load()

def ols_ci(x,y,x_eval):
    """OLS y~x; return slope b and fitted value f at x_eval, each with 95% CI + standard error."""
    x=np.asarray(x,float); y=np.asarray(y,float); n=len(x)
    if n<3: return None
    xbar=x.mean(); Sxx=((x-xbar)**2).sum()
    b=((x-xbar)*(y-y.mean())).sum()/Sxx; a=y.mean()-b*xbar
    f=a+b*x_eval; dof=n-2
    s2=((y-(a+b*x))**2).sum()/dof; t=ss.t.ppf(0.975,dof)
    se_b=np.sqrt(s2/Sxx); se_f=np.sqrt(s2*(1.0/n+(x_eval-xbar)**2/Sxx))
    return dict(b=b,b_lo=b-t*se_b,b_hi=b+t*se_b,se_b=se_b,
                f=f,f_lo=f-t*se_f,f_hi=f+t*se_f,se_f=se_f,n=n)

def lnr(loc,a,y):
    c=cell_for(data,loc,y,a,"T")
    if c is None or c[1]<=0 or c[0]<=0: return None
    return np.log(c[0]/c[1])

def mslope(loc,a,c):                       # centred 5-yr log-slope (%/yr) at window centre c
    pts=[(y,lnr(loc,a,y)) for y in range(c-2,c+3)]; pts=[(y,v) for y,v in pts if v is not None]
    if len(pts)<3: return None
    x=np.array([p[0] for p in pts],float); yv=np.array([p[1] for p in pts]); xc=x-x.mean()
    return 100.0*(xc*(yv-yv.mean())).sum()/(xc**2).sum()

def perband(loc):
    rsd=RELIABLE.get(loc,2003); out={}
    for a in FINE:
        pre=[(c,mslope(loc,a,c)) for c in range(rsd+2,2018)]
        pre=[(c,v) for c,v in pre if v is not None]
        if len(pre)<3: out[a]=None; continue
        out[a]=ols_ci([c for c,_ in pre],[v for _,v in pre],2019)
    return out

def het(est,se):                           # fixed-effect (inverse-variance) Q / I^2 / p across bands
    est=np.asarray(est,float); w=1.0/np.asarray(se,float)**2
    m=np.isfinite(est)&np.isfinite(w); est,w=est[m],w[m]; k=len(est)
    if k<2: return dict(Q=np.nan,I2=np.nan,p=np.nan,k=k)
    tbar=(w*est).sum()/w.sum(); Q=(w*(est-tbar)**2).sum(); df=k-1
    I2=max(0.0,(Q-df)/Q)*100.0; p=ss.chi2.sf(Q,df)
    return dict(Q=Q,I2=I2,p=p,k=k)

# ---- all-ages "paper" values (dashed line): u_islope2019 (g19) and u_sos (q) per country ----
paper={}
for r in csv.DictReader(open(os.path.join(OUTD,"slope_of_slopes_CI.csv"))):
    if r.get("country") in (l for l,_ in LOCS):
        paper[r["country"]]=(float(r["u_islope2019"]),float(r["u_sos"]))

# ---- compute + write CSV ----
res={l:perband(l) for l,_ in LOCS}
rows=[]
for l,nm in LOCS:
    for a in FINE:
        d=res[l][a]
        if d is None: continue
        rows.append([l,nm,a,round(d["f"],4),round(d["f_lo"],4),round(d["f_hi"],4),round(d["se_f"],4),
                     round(d["b"],4),round(d["b_lo"],4),round(d["b_hi"],4),round(d["se_b"],4)])
with open(os.path.join(OUTD,"perband_trend_heterogeneity_bigpops.csv"),"w",newline="") as f:
    w=csv.writer(f); w.writerow(["code","name","band","g19","g19_lo","g19_hi","se_g19","q","q_lo","q_hi","se_q"])
    w.writerows(rows)

# ---- heterogeneity table S7 (all-20 and 65+ subset, for g19 and q) ----
S7=[]
for l,nm in LOCS:
    d=res[l]; bands=[a for a in FINE if d[a] is not None]; b65=[a for a in bands if a in GE65]
    for lbl,param,key in (("slope in 2019","g19","f"),("slope-of-slopes","q","b")):
        allh=het([d[a][key] for a in bands],[d[a]["se_"+ ("f" if key=="f" else "b")] for a in bands])
        h65=het([d[a][key] for a in b65],[d[a]["se_"+ ("f" if key=="f" else "b")] for a in b65])
        S7.append([nm,lbl,allh["k"],round(allh["Q"],1),round(allh["I2"]),f"{allh['p']:.1e}",
                   h65["k"],round(h65["Q"],1),round(h65["I2"]),f"{h65['p']:.1e}"])

try:
    import openpyxl
    wb=openpyxl.Workbook(); ws=wb.active; ws.title="Table S7"
    ws.append(["Country","Parameter","k(all)","Q(all)","I2%(all)","p(all)","k(65+)","Q(65+)","I2%(65+)","p(65+)"])
    for r in S7: ws.append(r)
    os.makedirs(DOCS,exist_ok=True); wb.save(os.path.join(DOCS,"Table_S7_perband_heterogeneity.xlsx"))
except Exception as e:
    print("Table S7 xlsx skipped:",e)

# I^2 lookup for panel titles: (country,param)->(I2 all, I2 65+)
I2={(r[0],r[1]):(r[4],r[8]) for r in S7}

# ---- figure S3: 3 rows (countries) x 2 cols (params) ----
plt.rcParams.update({"font.size":10})
fig,axes=plt.subplots(3,2,figsize=(11,13)); xpos=np.arange(len(FINE))
COL_LT="#3b6fb0"; COL_GE="#c0392b"
for i,(l,nm) in enumerate(LOCS):
    d=res[l]
    for j,(lbl,key,unit,pidx) in enumerate((("slope in 2019","f","%/yr",0),("slope-of-slopes","b","%/yr²",1))):
        ax=axes[i][j]; pv=paper[l][pidx]
        for k,a in enumerate(FINE):
            dd=d[a]
            if dd is None: continue
            est=dd["f"] if key=="f" else dd["b"]; lo=dd[key+"_lo"]; hi=dd[key+"_hi"]
            c=COL_GE if a in GE65 else COL_LT
            ax.plot([k,k],[lo,hi],color=c,lw=1.4,alpha=0.75,zorder=2)
            ax.plot(k,est,"o",color=c,ms=5,zorder=3)
        ax.axhline(pv,ls="--",color="black",lw=1.4,label=f"all-ages (paper) = {pv:+.2f}")
        ax.axhline(0,color="0.6",lw=0.8,zorder=1)
        ax.axvline(13.5,ls=":",color="0.5",lw=1.0)              # 65 cut-off (between 60-64 and 65-69)
        i2a,i2b=I2[(nm,lbl)]
        ax.set_title(f"{nm} — {lbl}   (I² {i2a}% all·{i2b}% 65+)",fontweight="bold",fontsize=10.5)
        ax.set_xticks(xpos); ax.set_xticklabels(FINE,rotation=90,fontsize=7.5)
        ax.set_ylabel(f"{lbl}\n({unit})"); ax.legend(loc="upper left",fontsize=8.5,framealpha=0.9)
        if i==2: ax.set_xlabel("age band")
fig.suptitle("Per–age–band pre-pandemic mortality trend (two-step) in the three largest populations\n"
             "points = per-band estimate · bars = 95% CI · dashed = all-ages value applied to every band · dotted = 65 cut-off",
             fontweight="bold",fontsize=12,y=0.995)
h=[plt.Line2D([],[],marker="o",ls="",color=COL_LT,label="<65 age bands"),
   plt.Line2D([],[],marker="o",ls="",color=COL_GE,label="65+ age bands")]
fig.legend(handles=h,loc="upper center",ncol=2,fontsize=10,bbox_to_anchor=(0.5,0.95),frameon=False)
fig.tight_layout(rect=[0,0,1,0.945])
fig.savefig(os.path.join(OUTD,"fig_perband_trend_heterogeneity_bigpops.png"),dpi=150,bbox_inches="tight")
print("wrote perband_trend_heterogeneity_bigpops.csv (%d rows), fig_perband_trend_heterogeneity_bigpops.png, Table_S7"%len(rows))
for r in S7: print("  ",r[0],r[1],"I2 all=%s%% 65+=%s%%"%(r[4],r[8]))
