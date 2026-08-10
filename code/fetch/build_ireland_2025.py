#!/usr/bin/env python3
"""
Ireland 2025 fine-band deaths + population, appended to master_5x1_DPM_90plus.csv (2026-07-21).

Corrects the earlier "no age-stratified Irish 2025 source" gap. Ireland's national fine-band data DO
exist: CSO table VSA07 "Deaths by area, sex, and age-group" (PxStat JSON-stat) gives deaths OCCURRING by
5-year band (Under 1, 1-4, ... 80-84, 85+), both sexes / Male / Female, 2007-2025 with 2025 as a full year.
This is cleaner than the STMF broad-band route used for Israel -- it is real national fine-band data.

  DEATHS 2025 = VSA07 State-level deaths occurring, per 5-yr band and sex. VSA07 tops out at 85+, whereas
                the master series (IRL_EUROSTAT, used by the pipeline via the IRL_EUROSTAT->IRL remap) runs
                to 90+. We therefore split VSA07's 85+ deaths into 85-89 / 90+ using Ireland's OWN 2024
                death shares in those two bands, computed PER SEX (F skews to 90+, M to 85-89).
  POP    2025 = mid-year exposure extrapolated linearly from the master's existing mid-year series:
                P_2025(band,sex) = 2*P_2024 - P_2023. Same "extrapolated exposure" convention as every other
                2025 European addition (build_2025_eurostat_wk.py).

Written under location IRL_EUROSTAT (so the existing IRL_EUROSTAT->IRL pipeline remap picks it up), but with
source="CSO_VSA07" to record the true provenance honestly.

CAVEATS (surfaced, not buried): (1) 2025 is by occurrence and PROVISIONAL -- Ireland's ~3-month registration
window means late-2025 deaths may still rise slightly. (2) The 85-89/90+ split uses the 2024 age structure.
(3) 2025 population is extrapolated, not observed.

Appends 60 rows (IRL_EUROSTAT 2025 x 20 bands x {the row carries F,M,T together}). Backs up the master first.
"""
import os, csv, json, shutil, datetime
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
MASTER=os.path.join(ROOT,"data","master_5x1_DPM_90plus.csv")
VSA07=os.path.join(ROOT,"data","IRL_CSO","raw","VSA07.json")
LOC="IRL_EUROSTAT"; SRC="CSO_VSA07"; YEAR="2025"

FINE=[("0",0),("1-4",1),("5-9",5),("10-14",10),("15-19",15),("20-24",20),("25-29",25),("30-34",30),
      ("35-39",35),("40-44",40),("45-49",45),("50-54",50),("55-59",55),("60-64",60),("65-69",65),
      ("70-74",70),("75-79",75),("80-84",80),("85-89",85),("90+",90)]
# VSA07 age label -> our band key (85+ handled specially; All ages skipped)
V2B={"Under 1 year":"0","1 - 4 years":"1-4","5 - 9 years":"5-9","10 - 14 years":"10-14",
     "15 - 19 years":"15-19","20 - 24 years":"20-24","25 - 29 years":"25-29","30 - 34 years":"30-34",
     "35 - 39 years":"35-39","40 - 44 years":"40-44","45 - 49 years":"45-49","50 - 54 years":"50-54",
     "55 - 59 years":"55-59","60 - 64 years":"60-64","65 - 69 years":"65-69","70 - 74 years":"70-74",
     "75 - 79 years":"75-79","80 - 84 years":"80-84","85 years and over":"85+"}

# ---------- 1. read VSA07: State x 2025 deaths occurring, per band per sex ----------
d=json.load(open(VSA07)); ds=d.get("dataset",d)
order=ds["id"]; sizes=ds["size"]; dims=ds["dimension"]
val=ds["value"]
def idx(dim,label):
    cat=dims[dim]["category"]; ind=cat["index"]; lab=cat["label"]
    # index maps code->position; label maps code->text. invert label to find code for a text.
    code=next(c for c,t in lab.items() if t==label)
    return ind[code] if isinstance(ind,dict) else ind.index(code)
SZ={k:sizes[i] for i,k in enumerate(order)}
def get(stat,time,area,sex,agelabel):
    pos={"STATISTIC":idx("STATISTIC",stat),order[1]:idx(order[1],time),
         order[2]:idx(order[2],area),order[3]:idx(order[3],sex),order[4]:idx(order[4],agelabel)}
    fi=0                                    # row-major flat index over dim order
    for k in order: fi=fi*SZ[k]+pos[k]
    v=val[fi] if isinstance(val,list) else val.get(str(fi))
    return None if v in (None,"",".","..") else float(v)

SEXV={"F":"Female","M":"Male","T":"Both sexes"}
draw={}  # draw[(band,sex)] = deaths ; band includes "85+"
for sx in ("F","M","T"):
    for vlab,bkey in V2B.items():
        draw[(bkey,sx)]=get("Deaths Occurring",YEAR,"State",SEXV[sx],vlab)
# all-ages check value
allages={sx:get("Deaths Occurring",YEAR,"State",SEXV[sx],"All ages") for sx in ("F","M","T")}

# ---------- 2. read master: 2023 & 2024 pop (all bands) + 2024 85+ death split ----------
rows=list(csv.DictReader(open(MASTER))); FN=rows[0].keys()
def cell(yr,band):
    for r in rows:
        if r["location"]==LOC and r["year"]==yr and r["age"]==band: return r
    return None
def f(r,c):
    try: return float(r[c])
    except: return None
# 2024 death shares within 85+ per sex, for splitting VSA07's 85+
split={}
for sx in ("F","M","T"):
    a=cell("2024","85-89"); b=cell("2024","90+")
    da=f(a,f"D_{sx}"); db=f(b,f"D_{sx}")
    split[sx]=(da/(da+db), db/(da+db))   # (share_85-89, share_90+)
# location_name
LNAME=next(r["location_name"] for r in rows if r["location"]==LOC)

# ---------- 3. assemble 2025 fine-band D and P ----------
out={}
for band,astart in FINE:
    out[band]={}
    for sx in ("F","M","T"):
        # deaths
        if band=="85-89": D=draw[("85+",sx)]*split[sx][0]
        elif band=="90+": D=draw[("85+",sx)]*split[sx][1]
        else:             D=draw[(band,sx)]
        # population: 2*P2024 - P2023
        r24=cell("2024",band); r23=cell("2023",band)
        p24=f(r24,f"P_{sx}"); p23=f(r23,f"P_{sx}")
        P=2*p24-p23
        out[band][sx]=(D,P)

# ---------- 4. validations ----------
tot={sx:sum(out[b][sx][0] for b,_ in FINE) for sx in ("F","M","T")}
print("VSA07 State 2025 all-ages (deaths occurring):  F=%.0f  M=%.0f  T=%.0f"%(allages["F"],allages["M"],allages["T"]))
print("Rebuilt fine-band sum (after 85+ split):       F=%.0f  M=%.0f  T=%.0f"%(tot["F"],tot["M"],tot["T"]))
print("  -> split conserves total?  Tdiff=%.3f"%(tot["T"]-allages["T"]))
print("CSO published REGISTERED deaths 2025 = 35,587 (occurring here = %.0f; occurring<registered expected, provisional)"%tot["T"])
# 2024 cross-check: VSA07 State 2024 all-ages vs master IRL_EUROSTAT 2024 D_T
v24=get("Deaths Occurring","2024","State","Both sexes","All ages")
m24=sum(f(cell("2024",b),"D_T") for b,_ in FINE)
print("2024 cross-check: VSA07 State all-ages=%.0f  vs master IRL_EUROSTAT sum=%.0f  (diff %.0f)"%(v24,m24,v24-m24))
print("F+M vs T per-band consistency (should ~match T):  Tdiff=%.1f"%(sum(out[b]["F"][0]+out[b]["M"][0] for b,_ in FINE)-tot["T"]))

# ---------- 5. append ----------
bk=os.path.join(ROOT,"data","old",f"master_5x1_DPM_90plus_preIRL2025_{datetime.date(2026,7,21)}.csv")
os.makedirs(os.path.dirname(bk),exist_ok=True); shutil.copy2(MASTER,bk); print("backup ->",bk)
newrows=[]
for band,astart in FINE:
    (DF,PF),(DM,PM),(DT,PT)=out[band]["F"],out[band]["M"],out[band]["T"]
    row={c:"" for c in FN}
    row.update({"location":LOC,"location_name":LNAME,"source":SRC,"year":YEAR,"age":band,"age_start":astart,
                "D_F":f"{DF:.1f}","D_M":f"{DM:.1f}","D_T":f"{DT:.1f}",
                "P_F":f"{PF:.1f}","P_M":f"{PM:.1f}","P_T":f"{PT:.1f}",
                "M_F":f"{DF/PF:.8f}","M_M":f"{DM/PM:.8f}","M_T":f"{DT/PT:.8f}"})
    newrows.append(row)
with open(MASTER,"a",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(FN));
    for r in newrows: w.writerow(r)
print(f"appended {len(newrows)} rows: {LOC} {YEAR} (source {SRC})")
