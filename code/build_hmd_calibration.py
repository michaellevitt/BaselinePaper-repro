#!/usr/bin/env python3
"""Uniform HMD-basis calibration of the recent/national segments (Michael, 2026-07-23). For each location
with an HMD/national overlap year Y*, compute per-band ratio rho_a = HMD_rate_a(Y*) / national_rate_a(Y*)
(both sexes). The national segment is then put on the HMD (exposure) basis by dividing its POPULATION by
rho_a (deaths unchanged) — so the recent rates and the 2025 anchor sit on the same basis as the pre-2019
HMD baseline. No-overlap locations (GRC,SVN,JPN,KOR,TWN,USA) get rho=1 (documented). Ratios capped to
[0.80,1.25] to guard young-band noise. Writes output/hmd_calibration.csv ."""
import os, sys, csv, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import excess_anchor_window as E
data,name=E.load()
FINE=["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
CANON="AUS AUT BEL BGR CAN CHE CHL CZE DEUTNP DNK ESP EST FIN FRATNP GBRTENW GRC HKG HRV HUN IRL ISL ISR ITA JPN KOR LTU LUX LVA NLD NOR NZL_NP POL PRT SVK SVN SWE TWN USA".split()
def yrs(code): return sorted({y for (y,a) in data.get(code,{})})
def rate(code,y,a):
    c=data.get(code,{}).get((y,a));
    return (c["T"][0]/c["T"][1]) if c and "T" in c and c["T"][1]>0 else None
rows=[]; summ=[]
for c in CANON:
    p=set(yrs(c)); children=E.CHILD.get(c,[])
    # pick overlap year = max year in BOTH parent(HMD) and some child
    best=None
    for b in children:
        ov=sorted(p & set(yrs(b)))
        if ov and (best is None or ov[-1]>best[1]): best=(b,ov[-1])
    if best is None:
        for a in FINE: rows.append([c,a,"1.000","","no-overlap → rho=1"])
        summ.append((c,"—",1.0,1.0)); continue
    child,ystar=best; rr=[]
    for a in FINE:
        rh=rate(c,ystar,a); rn=rate(child,ystar,a)
        if rh and rn and rn>0:
            rho=rh/rn; rho=min(max(rho,0.80),1.25)
        else: rho=1.0
        rr.append(rho); rows.append([c,a,f"{rho:.4f}",str(ystar),child])
    # standardized-rate ratio (ESP) as the summary number
    summ.append((c,f"{child}@{ystar}",np.mean(rr),np.median(rr)))
with open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"output","hmd_calibration.csv"),"w",newline="") as f:
    w=csv.writer(f); w.writerow(["location","age","rho","overlap_year","child_source"]); w.writerows(rows)
print("Per-location calibration (rho = HMD_rate/national_rate at overlap; >1 => national rate was low, gets raised):")
print(f"  {'loc':9} {'source@yr':22} {'mean rho':>9} {'median rho':>11}")
for c,src,mn,md in sorted(summ,key=lambda x:-abs(x[2]-1)):
    flag=" <== large seam" if abs(md-1)>0.02 else ""
    print(f"  {c:9} {src:22} {mn:>9.3f} {md:>11.3f}{flag}")
