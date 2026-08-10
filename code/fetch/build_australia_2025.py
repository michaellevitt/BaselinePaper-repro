#!/usr/bin/env python3
"""
Australia 2025 fine-band deaths + population, appended to master_5x1_DPM_90plus.csv (2026-07-21).

Same STMF route used for Israel: STMF broad-band totals (0-14,15-64,65-74,75-84,85+) split into the 20 fine
bands using Australia's own latest-observed within-broad-band DEATH shares (HMD 2024).  Two Australia-specific
refinements over the bare Israel recipe, both surfaced (not buried):

  FULL YEAR.  The 2026-06-29 STMF pooled file (../21July2026_stmf.csv, downloaded 2026-07-21) carries AUS
  2025 for the FULL 52 ISO weeks -- no tail extrapolation is needed (this supersedes an earlier week-48
  estimate from the stale cached file).

  SOURCE SEAM.  STMF AUS 2024 totals run 1.1% below HMD AUS 2024 (per-band 0.91-0.99).  Since Australia's
  historical baseline (2015-2019) in the master is HMD, inserting STMF-native 2025 would manufacture a
  spurious ~1% deficit.  We therefore calibrate the full-year STMF 2025 to the HMD level with a single
  per-broad-band factor
     f_b = HMD_2024_full[b] / STMF_2024_full52[b]
applied to STMF_2025_full52[b].  (This is the Australia analogue of the negligible-seam Canada build, which
per Michael's formula uses no calibration because its seam is only 0.26%.)

  DEATHS 2025[a] = f_b * STMF_2025_full52[b] * (HMD_2024[a] / HMD_2024_broad[b])   for a in broad band b
  POP    2025[a] = 2*P_2024[a] - P_2023[a]   (extrapolated mid-year, same convention as the EUROSTAT 2025 adds)

Built TOTAL-only (both sexes), because Australia's existing 2023/2024 rows in the master are themselves
Total-only; the pipeline uses D_T/P_T throughout.  Uses STMF Sex="b".

Appends 20 rows (AUS 2025, source STMF_HMD_2425).  Idempotent: strips any existing AUS 2025 rows first.
Backs up the master first.

CAVEATS: the 0-14 band carries the largest STMF-vs-HMD gap (~10%) but is <1% of deaths; population 2025 is
extrapolated, not observed.
"""
import os, csv, shutil, datetime
from collections import defaultdict
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
MASTER=os.path.join(ROOT,"data","master_5x1_DPM_90plus.csv")
STMF=os.path.join(ROOT,"..","21July2026_stmf.csv")
LOC="AUS"; SRC="STMF_HMD_2425"; YEAR="2025"

FINE=[("0",0),("1-4",1),("5-9",5),("10-14",10),("15-19",15),("20-24",20),("25-29",25),("30-34",30),
      ("35-39",35),("40-44",40),("45-49",45),("50-54",50),("55-59",55),("60-64",60),("65-69",65),
      ("70-74",70),("75-79",75),("80-84",80),("85-89",85),("90+",90)]
BROAD={"D0_14":["0","1-4","5-9","10-14"],
       "D15_64":["15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64"],
       "D65_74":["65-69","70-74"],"D75_84":["75-79","80-84"],"D85p":["85-89","90+"]}
BBANDS=list(BROAD); band2broad={a:bb for bb,aa in BROAD.items() for a in aa}

# ---------- 1. STMF AUS broad-band sums, per sex ----------
rows=list(csv.reader(open(STMF))); hi=next(i for i,r in enumerate(rows[:5]) if r and r[0]=="CountryCode")
H=rows[hi]; ix={c:i for i,c in enumerate(H)}
def stmf_sum(year,sex,weeks=None):
    out=defaultdict(float)
    for r in rows[hi+1:]:
        if r[ix["CountryCode"]]!="AUS" or int(r[ix["Year"]])!=year or r[ix["Sex"]].lower()!=sex: continue
        if weeks and int(r[ix["Week"]]) not in weeks: continue
        for b in BBANDS:
            v=r[ix[b]]
            if v not in ("","."): out[b]+=float(v)
    return out
W152=set(range(1,53))
s24=stmf_sum(2024,"b",W152)          # STMF 2024 full 52 weeks, both sexes
s25=stmf_sum(2025,"b",W152)          # STMF 2025 full 52 weeks, both sexes

# ---------- 2. HMD 2024 fine-band deaths (Total) + 2023/2024 pop (Total) ----------
M=list(csv.DictReader(open(MASTER)))
def cell(yr,band):
    return next((r for r in M if r["location"]==LOC and r["year"]==yr and r["age"]==band),None)
hmd24_fine={a:float(cell("2024",a)["D_T"]) for a,_ in FINE}
hmd24_broad={bb:sum(hmd24_fine[a] for a in BROAD[bb]) for bb in BBANDS}
LNAME=next(r["location_name"] for r in M if r["location"]==LOC)

# ---------- 3. build 2025 fine-band deaths (Total) & population ----------
f_b={bb: hmd24_broad[bb]/s24[bb] for bb in BBANDS}   # HMD-calibration (full year, no tail)
out={}
for a,_ in FINE:
    bb=band2broad[a]
    Dbroad25=f_b[bb]*s25[bb]
    share=hmd24_fine[a]/hmd24_broad[bb]
    D=Dbroad25*share
    p24=float(cell("2024",a)["P_T"]); p23=float(cell("2023",a)["P_T"])
    P=2*p24-p23
    out[a]=(D,P)

# ---------- 4. validations ----------
totT=sum(out[a][0] for a,_ in FINE)
print("Per-broad-band HMD-calibration factor f_b = HMD2024full / STMF2024full52:")
for bb in BBANDS: print(f"   {bb:7} {f_b[bb]:.4f}   (STMF2024 full={s24[bb]:,.0f}  HMD2024 full={hmd24_broad[bb]:,.0f})")
print(f"\nAUS 2025 rebuilt deaths (Total):  T={totT:,.0f}")
raw25=sum(s25[b] for b in BBANDS)
print(f"  (raw STMF 2025 full52 = {raw25:,.0f}; HMD-calibration lifts to {totT:,.0f})")
hmd24T=sum(hmd24_fine[a] for a,_ in FINE)
print(f"HMD AUS 2024 total = {hmd24T:,.0f}  ->  2025/2024 = {totT/hmd24T:.4f}")
badP=[a for a,_ in FINE if out[a][1]<=0]
print("bands with non-positive extrapolated pop:", badP or "none")

# ---------- 5. idempotent append (Total-only, matching AUS 2023/2024 rows) ----------
FN=list(M[0].keys())
existing=[r for r in M if r["location"]==LOC and r["year"]==YEAR]
bk=os.path.join(ROOT,"data","old",f"master_5x1_DPM_90plus_preAUS2025full_{datetime.date(2026,7,21)}.csv")
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
