#!/usr/bin/env python3
"""
Canada 2025 fine-band deaths + population, appended to master_5x1_DPM_90plus.csv (2026-07-21).

Michael's specified full-year completion, applied per STMF broad band (0-14,15-64,65-74,75-84,85+):
    D(1-52, b, 2025) = D(1-45, b, 2025) + D(46-52, b, 2024) * D(1-45, b, 2025) / D(1-45, b, 2024)
                     = D(1-45, b, 2025) * D(1-52, b, 2024) / D(1-45, b, 2024)
i.e. observed 2025 weeks 1-45 grossed up by Canada's OWN 2024 full-year-to-week-45 ratio.  STMF Canada 2025
runs to ISO week 45 in the 2026-06-29 pooled file (../21July2026_stmf.csv, downloaded 2026-07-21).

No HMD calibration: STMF Canada 2024 sits only 0.26% below HMD 2024 (per-band 0.98-1.00), so STMF-native is
already on the HMD baseline level -- unlike Australia (1.1% seam).  This follows Michael's formula verbatim.

  DEATHS 2025[a] = g_b * D(1-45, b, 2025) * (HMD_2024[a] / HMD_2024_broad[b])   for a in broad band b,
                   where g_b = D(1-52, b, 2024) / D(1-45, b, 2024)              (both from STMF)
  POP    2025[a] = 2*P_2024[a] - P_2023[a]   (extrapolated mid-year, same convention as the other 2025 adds)

Built TOTAL-only (Canada's 2023/2024 master rows are Total-only; the pipeline uses D_T/P_T).  STMF Sex="b".
Idempotent: strips any existing CAN 2025 rows before appending.  Backs up the master first.

CAVEATS: weeks 46-52 (~13%, Nov-Dec = Canadian winter, a HIGH-mortality tail) are estimated from 2024's
seasonal shape scaled to the 2025 level; population 2025 is extrapolated, not observed.
"""
import os, csv, shutil, datetime
from collections import defaultdict
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
MASTER=os.path.join(ROOT,"data","master_5x1_DPM_90plus.csv")
STMF=os.path.join(ROOT,"..","21July2026_stmf.csv")
LOC="CAN"; SRC="STMF_2425"; YEAR="2025"; LASTWK=45

FINE=[("0",0),("1-4",1),("5-9",5),("10-14",10),("15-19",15),("20-24",20),("25-29",25),("30-34",30),
      ("35-39",35),("40-44",40),("45-49",45),("50-54",50),("55-59",55),("60-64",60),("65-69",65),
      ("70-74",70),("75-79",75),("80-84",80),("85-89",85),("90+",90)]
BROAD={"D0_14":["0","1-4","5-9","10-14"],
       "D15_64":["15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64"],
       "D65_74":["65-69","70-74"],"D75_84":["75-79","80-84"],"D85p":["85-89","90+"]}
BBANDS=list(BROAD); band2broad={a:bb for bb,aa in BROAD.items() for a in aa}

# ---------- 1. STMF Canada broad-band sums (both sexes) ----------
rows=list(csv.reader(open(STMF))); hi=next(i for i,r in enumerate(rows[:6]) if r and r[0]=="CountryCode")
ix={c:i for i,c in enumerate(rows[hi])}
def stmf_sum(year,weeks):
    o=defaultdict(float)
    for r in rows[hi+1:]:
        if r[ix["CountryCode"]]!=LOC or int(r[ix["Year"]])!=year or r[ix["Sex"]]!="b": continue
        if int(r[ix["Week"]]) not in weeks: continue
        for b in BBANDS:
            v=r[ix[b]]
            if v not in ("","."): o[b]+=float(v)
    return o
W145=set(range(1,LASTWK+1)); W152=set(range(1,53))
d25_145=stmf_sum(2025,W145)         # D(1-45, b, 2025)
d24_145=stmf_sum(2024,W145)         # D(1-45, b, 2024)
d24_152=stmf_sum(2024,W152)         # D(1-52, b, 2024)
g={bb: d24_152[bb]/d24_145[bb] for bb in BBANDS}        # Michael's per-band gross-up ratio

# ---------- 2. HMD 2024 fine-band deaths (Total) + 2023/2024 pop (Total) ----------
M=list(csv.DictReader(open(MASTER)))
def cell(yr,band):
    return next((r for r in M if r["location"]==LOC and r["year"]==yr and r["age"]==band),None)
hmd24_fine={a:float(cell("2024",a)["D_T"]) for a,_ in FINE}
hmd24_broad={bb:sum(hmd24_fine[a] for a in BROAD[bb]) for bb in BBANDS}
LNAME=next(r["location_name"] for r in M if r["location"]==LOC)

# ---------- 3. build 2025 fine-band deaths (Total) & population ----------
out={}
for a,_ in FINE:
    bb=band2broad[a]
    Dbroad25=g[bb]*d25_145[bb]                       # Michael's completion, per broad band
    share=hmd24_fine[a]/hmd24_broad[bb]              # within-broad-band age structure from HMD 2024
    D=Dbroad25*share
    p24=float(cell("2024",a)["P_T"]); p23=float(cell("2023",a)["P_T"])
    out[a]=(D,2*p24-p23)

# ---------- 4. validations ----------
totT=sum(out[a][0] for a,_ in FINE)
print("Michael formula per-broad-band gross-up g_b = D(1-52,2024)/D(1-45,2024):")
for bb in BBANDS:
    print(f"   {bb:7} g={g[bb]:.4f}   D25(1-45)={d25_145[bb]:,.0f} -> full {g[bb]*d25_145[bb]:,.0f}   [D24: 1-45={d24_145[bb]:,.0f}, 1-52={d24_152[bb]:,.0f}]")
print(f"\nCAN 2025 rebuilt deaths (Total):  T={totT:,.0f}")
print(f"  raw STMF 2025 w1-45 = {sum(d25_145.values()):,.0f}; completed to {totT:,.0f}  (tail = {totT-sum(d25_145.values()):,.0f}, {100*(totT-sum(d25_145.values()))/totT:.1f}%)")
hmd24T=sum(hmd24_fine[a] for a,_ in FINE)
print(f"HMD CAN 2024 total = {hmd24T:,.0f}  ->  2025/2024 = {totT/hmd24T:.4f}")
print(f"STMF CAN 2024 full = {sum(d24_152.values()):,.0f}  (STMF/HMD 2024 = {sum(d24_152.values())/hmd24T:.4f})")
badP=[a for a,_ in FINE if out[a][1]<=0]; print("non-positive extrapolated pop:", badP or "none")

# ---------- 5. idempotent append (Total-only) ----------
FN=list(M[0].keys())
existing=[r for r in M if r["location"]==LOC and r["year"]==YEAR]
bk=os.path.join(ROOT,"data","old",f"master_5x1_DPM_90plus_preCAN2025_{datetime.date(2026,7,21)}.csv")
os.makedirs(os.path.dirname(bk),exist_ok=True); shutil.copy2(MASTER,bk); print("backup ->",bk)
if existing:
    keep=[r for r in M if not (r["location"]==LOC and r["year"]==YEAR)]
    with open(MASTER,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=FN); w.writeheader(); w.writerows(keep)
    print(f"stripped {len(existing)} existing {LOC} {YEAR} rows (idempotent re-run)")
new=[]
for a,astart in FINE:
    DT,PT=out[a]
    row={c:"" for c in FN}
    row.update({"location":LOC,"location_name":LNAME,"source":SRC,"year":YEAR,"age":a,"age_start":astart,
                "D_T":f"{DT:.1f}","P_T":f"{PT:.1f}","M_T":f"{DT/PT:.8f}"})
    new.append(row)
with open(MASTER,"a",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=FN)
    for r in new: w.writerow(r)
print(f"appended {len(new)} rows: {LOC} {YEAR} (source {SRC})")
