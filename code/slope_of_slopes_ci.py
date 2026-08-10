#!/usr/bin/env python3
"""
PHASE 1 diagnostic (John, July 2026): per-country slope-of-slopes regression with 95% CIs, to judge
how much between-country variation in the trend-of-trends is genuine vs chance — the pivot for the
ATT/STT decision.

For each of the 33 populations, on the ALL-AGES ESP2013-standardised death rate (the montage series in
output/stmf_moving_slopes.json), we fit the trend-of-trends line g(t)=q*t+p to the moving 5-year slopes
(regressed against window-END year, matching the montage), reliable pre-pandemic windows only:
  UNANCHORED:  fit to reliable pre-pandemic slopes (window-ends >= reliability, <= 2019).
               report  intercept-slope @2019 = g(2019)   and   slope-of-slopes = q     (each +/-95% CI)
  ANCHORED:    add ONE anchor point = the 3-point {2019,2024,2025} rate-slope (D3) at x=2025.
               report  intercept-slope @2025 = g(2025)   and   slope-of-slopes = q     (each +/-95% CI)
Then population-weighted (2019 total population) global averages of the four quantities = the would-be
ATT parameters, plus unweighted mean and spread.

CIs: OLS, t-based (scipy.stats.t), residual variance SSE/(n-2); slope SE=sqrt(s2/Sxx); fitted-value
SE=sqrt(s2*(1/n+(x0-xbar)^2/Sxx)).  Writes output/slope_of_slopes_CI.csv and docs/Slope_of_slopes_CI_v1.xlsx .
"""
import os, sys, csv, json, numpy as np
import scipy.stats as ss
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, HERE)
from excess_anchor_window import load, cell_for
OUTD=os.path.join(ROOT,"output"); DOCS=os.path.join(ROOT,"docs")
FINE=["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
      "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
RELIABLE={"BGR":2015,"CHL":2011,"HRV":2007,"EST":2005,"HUN":2006,"LVA":2007,"LTU":2006,"POL":2008,"SVK":2006}
MORE={"PRT","SVN","ESP","HUN","EST","GBRTENW","ITA","HRV","GRC","LVA","SVK","CZE","POL","LTU","CHL","USA","BGR"}
LESS={"SWE","DNK","NOR","KOR","LUX","CHE","FIN","BEL","DEUTNP","FRATNP","NLD","AUT"}
def vuln(loc): return "more" if loc in MORE else ("less" if loc in LESS else "unclassified")

# ---- OLS with 95% CIs ----
def ols_ci(x,y,x_eval):
    x=np.asarray(x,float); y=np.asarray(y,float); n=len(x)
    if n<2: return dict(b=np.nan,b_lo=np.nan,b_hi=np.nan,f=np.nan,f_lo=np.nan,f_hi=np.nan,n=n)
    xbar=x.mean(); Sxx=((x-xbar)**2).sum(); b=((x-xbar)*(y-y.mean())).sum()/Sxx; a=y.mean()-b*xbar
    f=a+b*x_eval; dof=n-2
    if dof<1:                                    # exact fit, no residual dof -> point est only
        return dict(b=b,b_lo=np.nan,b_hi=np.nan,f=f,f_lo=np.nan,f_hi=np.nan,n=n)
    s2=((y-(a+b*x))**2).sum()/dof; t=ss.t.ppf(0.975,dof)
    se_b=np.sqrt(s2/Sxx); se_f=np.sqrt(s2*(1.0/n+(x_eval-xbar)**2/Sxx))
    return dict(b=b,b_lo=b-t*se_b,b_hi=b+t*se_b,f=f,f_lo=f-t*se_f,f_hi=f+t*se_f,n=n)

# ---- ESP2013 log rate from the SAME data the engine uses (cell_for); CENTRE convention (2026-07-29) ----
# Each 5-yr moving slope is indexed at its window CENTRE (not end); the {2019,2024,2025} return-slope anchor
# is placed at its centroid (~2022.7).  Computed from cell_for so the trend and the excess baseline share one
# data source (fixes the CHL/IRL/AUS slope-JSON vs baseline mismatch).
data,name=load()
ESP={"0":1000,"1-4":4000,"5-9":5500,"10-14":5500,"15-19":5500,"20-24":6000,"25-29":6000,
     "30-34":6500,"35-39":7000,"40-44":7000,"45-49":7000,"50-54":7000,"55-59":6500,"60-64":6000,
     "65-69":5500,"70-74":5000,"75-79":4000,"80-84":2500,"85-89":1500,"90+":1000}; ESPsum=sum(ESP.values())
def lnR(loc,y):
    cs=[cell_for(data,loc,y,a,"T") for a in FINE]
    if any(c is None for c in cs): return None
    r=sum(ESP[a]*cs[i][0]/cs[i][1] for i,a in enumerate(FINE))/ESPsum
    return np.log(r) if r>0 else None
def mslope(loc,c):                     # centred 5-yr log-slope (%/yr) at window centre c
    pts=[(y,lnR(loc,y)) for y in range(c-2,c+3)]; pts=[(y,v) for y,v in pts if v is not None]
    if len(pts)<3: return None
    x=np.array([p[0] for p in pts],float); yv=np.array([p[1] for p in pts]); xc=x-x.mean()
    return 100.0*(xc*(yv-yv.mean())).sum()/(xc**2).sum()
def anchor(loc):                       # return slope through {2019,2024,2025} (log-slope %/yr) at its centroid x
    ys=[y for y in (2019,2024,2025) if lnR(loc,y) is not None]
    if len(ys)<2: return None,None
    x=np.array(ys,float); yv=np.array([lnR(loc,y) for y in ys]); xc=x-x.mean()
    return 100.0*(xc*(yv-yv.mean())).sum()/(xc**2).sum(), float(x.mean())
def pop2019(loc):
    cs=[cell_for(data,loc,2019,a,"T") for a in FINE]
    return None if any(c is None for c in cs) else sum(c[1] for c in cs)
LOCS=[(r["loc"],r["name"]) for r in json.load(open(os.path.join(OUTD,"stmf_moving_slopes.json")))["data"]["All ages"]["rows"]]

rows=[]
for loc,nm in LOCS:
    rsd=RELIABLE.get(loc,2003)
    pre=[(c,mslope(loc,c)) for c in range(rsd+2,2018)]              # centres up to 2017 = full 5-yr windows within [reliable,2019]
    pre=[(c,v) for c,v in pre if v is not None]
    xp=[c for c,_ in pre]; yp=[v for _,v in pre]
    U=ols_ci(xp,yp,2019)                                            # unanchored, evaluate at 2019
    asl,acx=anchor(loc)                                             # anchor slope + centroid x (~2022.7)
    if asl is not None and len(xp)>=1:
        A=ols_ci(xp+[acx],yp+[asl],2025)                            # anchored (anchor at centroid), evaluate at 2025 (engine back-computes to 2019)
    else:
        A=dict(b=np.nan,b_lo=np.nan,b_hi=np.nan,f=np.nan,f_lo=np.nan,f_hi=np.nan,n=len(xp))
    p19=pop2019(loc)
    rows.append(dict(loc=loc,name=nm,vuln=vuln(loc),pop=p19,npre=len(pre),
        u_is=U["f"],u_is_lo=U["f_lo"],u_is_hi=U["f_hi"],u_sos=U["b"],u_sos_lo=U["b_lo"],u_sos_hi=U["b_hi"],
        a_is=A["f"],a_is_lo=A["f_lo"],a_is_hi=A["f_hi"],a_sos=A["b"],a_sos_lo=A["b_lo"],a_sos_hi=A["b_hi"]))
rows.sort(key=lambda d:d["name"])

# ---- weighted / unweighted global averages of the four point estimates ----
def agg(field):
    v=np.array([r[field] for r in rows],float); w=np.array([r["pop"] if r["pop"] else 0.0 for r in rows])
    m=~np.isnan(v)
    wmean=np.average(v[m],weights=w[m]) if m.any() else np.nan
    umean=float(np.nanmean(v)); sd=float(np.nanstd(v,ddof=1))
    return wmean,umean,sd,float(np.nanmin(v)),float(np.nanmax(v))
FIELDS=[("u_is","intercept-slope @2019 (unanchored)"),("u_sos","slope-of-slopes (unanchored)"),
        ("a_is","intercept-slope @2025 (anchored)"),("a_sos","slope-of-slopes (anchored)")]
GLOB={f:agg(f) for f,_ in FIELDS}

# ---- CSV ----
cols=["country","name","vulnerable","pop2019_millions","n_pre",
      "u_islope2019","u_islope2019_lo","u_islope2019_hi","u_sos","u_sos_lo","u_sos_hi",
      "a_islope2025","a_islope2025_lo","a_islope2025_hi","a_sos","a_sos_lo","a_sos_hi"]
def rr(x,nd=3): return "" if (x is None or (isinstance(x,float) and np.isnan(x))) else round(float(x),nd)
with open(os.path.join(OUTD,"slope_of_slopes_CI.csv"),"w",newline="") as f:
    w=csv.writer(f); w.writerow(cols)
    for r in rows:
        w.writerow([r["loc"],r["name"],r["vuln"],rr(r["pop"]/1e6 if r["pop"] else None,3),r["npre"],
            rr(r["u_is"]),rr(r["u_is_lo"]),rr(r["u_is_hi"]),rr(r["u_sos"]),rr(r["u_sos_lo"]),rr(r["u_sos_hi"]),
            rr(r["a_is"]),rr(r["a_is_lo"]),rr(r["a_is_hi"]),rr(r["a_sos"]),rr(r["a_sos_lo"]),rr(r["a_sos_hi"])])
    w.writerow([])
    w.writerow(["WEIGHTED MEAN (2019 pop)","","","","",rr(GLOB["u_is"][0]),"","",rr(GLOB["u_sos"][0]),"","",
                rr(GLOB["a_is"][0]),"","",rr(GLOB["a_sos"][0]),"",""])
    w.writerow(["UNWEIGHTED MEAN","","","","",rr(GLOB["u_is"][1]),"","",rr(GLOB["u_sos"][1]),"","",
                rr(GLOB["a_is"][1]),"","",rr(GLOB["a_sos"][1]),"",""])
    w.writerow(["SD across countries","","","","",rr(GLOB["u_is"][2]),"","",rr(GLOB["u_sos"][2]),"","",
                rr(GLOB["a_is"][2]),"","",rr(GLOB["a_sos"][2]),"",""])

# ---- console summary ----
print("Per-country slope-of-slopes regression (all-ages, ESP2013-standardised). Units: %/yr and %/yr/yr.\n")
print("GLOBAL AVERAGES (the would-be ATT parameters):")
print(f"  {'quantity':38}{'wtd(2019 pop)':>14}{'unwtd':>9}{'SD':>8}{'min':>8}{'max':>8}")
for f,lab in FIELDS:
    wm,um,sd,mn,mx=GLOB[f]; print(f"  {lab:38}{wm:14.3f}{um:9.3f}{sd:8.3f}{mn:8.3f}{mx:8.3f}")
# how many countries' 95% CI excludes 0 -> a "genuine" (non-chance) signal
def ci_excl0(lo_f,hi_f):
    k=tot=0
    for r in rows:
        lo,hi=r[lo_f],r[hi_f]
        if np.isnan(lo) or np.isnan(hi): continue
        tot+=1; k+=(lo>0 or hi<0)
    return k,tot
print("\nFraction of countries whose 95% CI excludes 0 (rest are consistent with chance):")
for lab,lo,hi in [("slope-of-slopes (unanchored)","u_sos_lo","u_sos_hi"),
                  ("slope-of-slopes (anchored)","a_sos_lo","a_sos_hi"),
                  ("intercept-slope @2019","u_is_lo","u_is_hi"),
                  ("intercept-slope @2025","a_is_lo","a_is_hi")]:
    k,tot=ci_excl0(lo,hi); print(f"  {lab:32}: {k}/{tot}")

# ---- styled xlsx for John ----
try:
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
    from openpyxl.utils import get_column_letter
    wb=Workbook(); ws=wb.active; ws.title="slope-of-slopes CI"
    HDR=PatternFill("solid",fgColor="1F4E79"); WHT=Font(color="FFFFFF",bold=True)
    UBLK=PatternFill("solid",fgColor="DDEBF7"); ABLK=PatternFill("solid",fgColor="FCE4D6")
    MOREF=PatternFill("solid",fgColor="FDECEA"); SUMF=PatternFill("solid",fgColor="FFF2CC")
    thin=Side(style="thin",color="BFBFBF"); BORD=Border(left=thin,right=thin,top=thin,bottom=thin); CEN=Alignment(horizontal="center")
    ws.append(["","","","","","UNANCHORED (evaluated at 2019)","","","","","","ANCHORED (3-pt 2019/2024/2025, at 2025)","","","","",""])
    ws.append(["Country","Name","Vuln.","Pop 2019 (M)","n pre","islope@2019","95% lo","95% hi","SoS","95% lo","95% hi",
               "islope@2025","95% lo","95% hi","SoS","95% lo","95% hi"])
    for j in range(1,18):
        c=ws.cell(2,j); c.fill=HDR; c.font=WHT; c.alignment=CEN; c.border=BORD
        if 6<=j<=11: ws.cell(1,j).fill=UBLK
        elif j>=12: ws.cell(1,j).fill=ABLK
    ws.cell(1,6).font=Font(bold=True); ws.cell(1,12).font=Font(bold=True)
    ws.merge_cells(start_row=1,start_column=6,end_row=1,end_column=11); ws.merge_cells(start_row=1,start_column=12,end_row=1,end_column=17)
    ws.cell(1,6).alignment=CEN; ws.cell(1,12).alignment=CEN
    for r in rows:
        row=[r["loc"],r["name"],r["vuln"],(round(r["pop"]/1e6,2) if r["pop"] else None),r["npre"],
             rr(r["u_is"],3),rr(r["u_is_lo"],3),rr(r["u_is_hi"],3),rr(r["u_sos"],3),rr(r["u_sos_lo"],3),rr(r["u_sos_hi"],3),
             rr(r["a_is"],3),rr(r["a_is_lo"],3),rr(r["a_is_hi"],3),rr(r["a_sos"],3),rr(r["a_sos_lo"],3),rr(r["a_sos_hi"],3)]
        ws.append(row); rn=ws.max_row
        for j in range(1,18):
            c=ws.cell(rn,j); c.border=BORD
            if j>=3: c.alignment=CEN
            if r["vuln"]=="more" and j<=5: c.fill=MOREF
    ws.append([]);
    for label,key in [("WEIGHTED MEAN (2019 pop)",0),("UNWEIGHTED MEAN",1),("SD across countries",2)]:
        ws.append([label,"","","","",rr(GLOB["u_is"][key],3),"","",rr(GLOB["u_sos"][key],3),"","",
                   rr(GLOB["a_is"][key],3),"","",rr(GLOB["a_sos"][key],3),"",""]); rn=ws.max_row
        for j in range(1,18):
            c=ws.cell(rn,j); c.fill=SUMF; c.font=Font(bold=True); c.alignment=CEN; c.border=BORD
    ws.column_dimensions["A"].width=11; ws.column_dimensions["B"].width=16
    for j in range(3,18): ws.column_dimensions[get_column_letter(j)].width=10
    ws.freeze_panes="C3"
    n2=wb.create_sheet("Notes")
    for i,t in enumerate([
      "Per-country trend-of-trends (slope-of-slopes) regression — all-ages ESP2013-standardised rate.",
      "",
      "Moving 5-year slopes of the standardised rate (%/yr of window mean) are regressed against window-END year.",
      "  UNANCHORED: reliable pre-pandemic windows only (window-end >= data-reliable year, <= 2019); evaluated at 2019.",
      "  ANCHORED  : adds ONE anchor point = the 3-point {2019,2024,2025} rate-slope, placed at 2025; evaluated at 2025.",
      "islope = intercept-slope = fitted slope-of-slopes line g(t) at the reference year (%/yr).",
      "SoS = slope-of-slopes = the regression coefficient q (%/yr per year; positive = decelerating decline).",
      "95% CIs are OLS t-based (slope: sqrt(s2/Sxx); fitted value: sqrt(s2*(1/n+(x0-xbar)^2/Sxx))).",
      "Bulgaria has one reliable pre-pandemic window only, so its unanchored fit / CIs are undefined.",
      "The WEIGHTED MEAN row gives the population-weighted global parameters = the proposed ATT model inputs.",
    ],1): n2.cell(i,1,t).font=Font(bold=(i==1))
    n2.column_dimensions["A"].width=115
    XLSX=os.path.join(DOCS,"Slope_of_slopes_CI_v1.xlsx"); wb.save(XLSX); print(f"wrote {XLSX}")
except ImportError:
    pass
print("wrote output/slope_of_slopes_CI.csv")
