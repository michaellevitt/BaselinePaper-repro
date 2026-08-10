#!/usr/bin/env python3
"""
Table 2 (age-resolution) REBUILT on Option-1 (John #14: "update STTa").  Total 2020-2025 excess deaths for
each of the FOUR models (Fa, Ta, TTa, STTa; ATTa dropped per #7) at four baseline age resolutions —
CRUDE (1 band), 2 bands (<65/65+), 5 bands (STMF), 20 five-year bands (default).  Coarser baselines cannot
track the shifting age structure during 2020-2025, so their excess is biased (compare each to 20-band).

Option-1 estimator throughout: recency-weighted WLS log-quadratic (tau=6), FITTED-alpha anchor, clip mult
in [0.4,2.5].  Fa=flat 2017-19 mean; Ta=linear 2015-19; TTa=each coarse band's own log-quadratic trend;
STTa=each coarse band's fitted 2019 level carried by the population's PRECISION-SHRUNK age-standardised
trend (beta_shrunk,gamma_shrunk from option1_country_params.csv).  The 20-band column reproduces the
model_option1_periods headline by construction.
Writes docs/Age_resolution_excess_Option1_v1.xlsx + output/age_resolution_option1.csv .
"""
import os, sys, csv, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from excess_anchor_window import load, cell_for
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); OUTD=os.path.join(ROOT,"output"); DOCS=os.path.join(ROOT,"docs")
FINE=["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
GROUPINGS={"Crude (1 band)":[list(range(20))],
           "2 bands (<65/65+)":[list(range(0,14)),list(range(14,20))],
           "5 bands (STMF)":[[0,1,2,3],[4,5,6,7,8,9,10,11,12,13],[14,15],[16,17],[18,19]],
           "20 bands (default)":[[i] for i in range(20)]}
order=list(GROUPINGS); MODELS=["Fa","Ta","TTa","STTa","STTa+"]
W=5; SPAN=[2020,2021,2022,2023,2024,2025]; TAU=6.0; CLIP=(0.4,2.5)
RELIABLE={"BGR":2015,"CHL":2011,"HRV":2007,"EST":2005,"HUN":2006,"LVA":2007,"LTU":2006,"POL":2008,"SVK":2006}
data,name=load()
def DP(loc,y,i): return cell_for(data,loc,y,FINE[i],"T")
def has_year(loc,y): return all(DP(loc,y,i) is not None for i in range(20))
def has_base(loc):
    se=max(2003+(W-1),RELIABLE.get(loc,2003)+(W-1)); return all(DP(loc,y,i) is not None for y in range(se-(W-1),2020) for i in range(20))
PC=[r[0] for r in csv.reader(open(os.path.join(OUTD,"slope_of_slopes_CI.csv"))) if r and r[0]!="country" and not r[0].startswith(("WEIGHTED","UNWEIGHTED","SD"))]
countries=[l for l in sorted(PC) if has_base(l)]
# precision-shrunk standardised trend per country (natural-log beta,gamma)
def _f(x):
    try: return float(x)
    except (TypeError,ValueError): return None
shr={}
for r in csv.DictReader(open(os.path.join(OUTD,"option1_country_params.csv"))):
    shr[r["country"]]=dict(bo=_f(r["beta_i"]),go=_f(r["gamma_i"]),bs=_f(r["beta_shrunk"]),gs=_f(r["gamma_shrunk"]),
                           bas=_f(r["beta_anch_shrunk"]),gas=_f(r["gamma_anch_shrunk"]))
def baseyrs(l): return list(range(max(2003,RELIABLE.get(l,2003)),2020))
def band_series(loc,idxs):
    out={}
    for y in range(1990,2026):
        cs=[DP(loc,y,i) for i in idxs]
        if any(c is None for c in cs): continue
        out[y]=(sum(c[0] for c in cs),sum(c[1] for c in cs))
    return out
def clipmult(b,g,t): return min(max(np.exp(b*t+g*t*t),CLIP[0]),CLIP[1])

results={res:{m:0.0 for m in MODELS} for res in GROUPINGS}
for res,bands in GROUPINGS.items():
    for loc in countries:
        ys=baseyrs(loc); t=np.array(ys,float)-2019.0; om=np.ones(len(ys))   # uniform (recency dropped, 2026-07-27)
        span=[y for y in SPAN if has_year(loc,y)]
        brs=[band_series(loc,b) for b in bands]
        Fa=[];Ta=[];TT=[]
        for s in brs:
            Fa.append(sum(s[y][0] for y in (2017,2018,2019))/sum(s[y][1] for y in (2017,2018,2019)))
            Ta.append(np.polyfit(np.arange(2015,2020,dtype=float)-2012,[s[y][0]/s[y][1] for y in range(2015,2020)],1))
            r=np.array([s[y][0]/s[y][1] for y in ys]); TT.append(np.polyfit(t,np.log(np.maximum(r,1e-12)),2,w=np.sqrt(om)))  # [g,b,alpha]
        sp=shr[loc]; bs,gs=sp["bs"],sp["gs"]; bas,gas=sp["bas"],sp["gas"]
        bo,go=(sp["bo"],sp["go"]) if sp["bo"] is not None else (bs,gs)   # TTa own two-step (global fallback for BGR)
        for m in MODELS:
            O=E=0.0
            for y in span:
                for k,s in enumerate(brs):
                    P=s[y][1]; c=TT[k]; tt=y-2019.0                       # c[2] = band 2019 level; country trend applied to it
                    if m=="Fa": rt=Fa[k]
                    elif m=="Ta": rt=float(np.polyval(Ta[k],y-2012))
                    elif m=="TTa": rt=np.exp(c[2])*clipmult(bo,go,tt)     # own two-step (country)
                    elif m=="STTa+": rt=np.exp(c[2])*clipmult(bas,gas,tt) # slope-anchored shrunk trend
                    else: rt=np.exp(c[2])*clipmult(bs,gs,tt)             # STTa: band level, shrunk trend
                    E+=rt*P; O+=s[y][0]
            results[res][m]+=(O-E)

print(f"n countries: {len(countries)}\nTotal 2020-2025 excess (MILLIONS) by model x age resolution:")
print("  model  "+"".join(f"{r.split(' ')[0]:>11}" for r in order))
for m in MODELS: print(f"  {m:5} "+"".join(f"{results[r][m]/1e6:>11.2f}" for r in order))
with open(os.path.join(OUTD,"age_resolution_option1.csv"),"w",newline="") as f:
    w=csv.writer(f); w.writerow(["model"]+order)
    for m in MODELS: w.writerow([m]+[round(results[r][m]) for r in order])

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
wb=Workbook(); ws=wb.active; ws.title="Excess by age resolution"
HDR=PatternFill("solid",fgColor="1F4E79"); WHT=Font(color="FFFFFF",bold=True,size=10)
CEN=Alignment(horizontal="center"); thin=Side(style="thin",color="BFBFBF"); BORD=Border(left=thin,right=thin,top=thin,bottom=thin)
GOOD=PatternFill("solid",fgColor="C6E0B4")
ws["A1"]="Table 2. Total 2020–2025 excess deaths (MILLIONS) by baseline age resolution — Option-1 family (Fa, Ta, TTa, STTa), 38 populations"; ws["A1"].font=Font(bold=True,size=12)
ws["A2"]=("Each model refitted at each resolution (Option-1: recency-weighted WLS log-quadratic, tau=6, fitted-2019 anchor; STTa precision-shrunk). "
          "The 20-band column is the default; coarser baselines cannot track the shifting age structure during 2020–2025, so their excess is biased. "
          "% vs 20-band shown below each count.")
ws["A2"].font=Font(italic=True,size=9,color="555555"); ws.merge_cells("A2:E2"); ws.row_dimensions[2].height=54
for j,h in enumerate(["Model"]+order,1):
    c=ws.cell(4,j,h); c.fill=HDR; c.font=WHT; c.alignment=CEN; c.border=BORD; ws.column_dimensions[chr(64+j)].width=(9 if j==1 else 18)
r=5
for m in MODELS:
    ws.cell(r,1,m).font=Font(bold=True,size=10); ws.cell(r,1).border=BORD
    base=results["20 bands (default)"][m]
    for j,res in enumerate(order,2):
        v=results[res][m]; pct=100*(v-base)/base if base else 0
        c=ws.cell(r,j,round(v/1e6,2)); c.alignment=CEN; c.border=BORD; c.font=Font(size=10)
        if res=="20 bands (default)": c.fill=GOOD
        c2=ws.cell(r+1,j,("—" if res=="20 bands (default)" else f"{pct:+.0f}% vs 20-band")); c2.alignment=CEN; c2.border=BORD
        c2.font=Font(size=8,italic=True,color="777777")
    ws.cell(r+1,1,"").border=BORD; r+=2
OUT=os.path.join(DOCS,"Age_resolution_excess_Option1_v1.xlsx"); wb.save(OUT); print("wrote",OUT)
