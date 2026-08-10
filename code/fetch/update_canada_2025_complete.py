#!/usr/bin/env python3
"""
Replace the extrapolated CAN 2025 (STMF weeks 1-45 grossed up) with COMPLETE 2025 data, now that
Statistics Canada table 13-10-0768 (provisional weekly deaths by age group and sex) carries all 52
weeks of 2025 (release ~Jul 2026).

Method (HMD level, StatCan growth, StatCan age resolution):
  D2025_fine[a] = HMD_2024_fine[a] * ( StatCan_2025_group(g) / StatCan_2024_group(g) )   for a in group g
where g is the StatCan 4-way age group (0-44, 45-64, 65-84, 85+).  This applies the OBSERVED complete-2025
year-over-year change per age group to the HMD 2024 fine-band deaths (so the level stays on the master's
HMD basis, removing the StatCan/HMD ~0.3% seam, while the shape/growth are the observed complete data).
Population 2025 is unchanged (already the extrapolated mid-year 2*P2024-P2023).

Idempotent: re-derives from HMD 2024 + StatCan each run.  Backs up the master first.
"""
import os, csv, shutil, datetime
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
MASTER=os.path.join(ROOT,"data","master_5x1_DPM_90plus.csv")
SC=os.path.join(ROOT,"scratchpad","stmf","13100768.csv")
LOC="CAN"; SRC="STATCAN768_2025complete"; YEAR="2025"
FINE=["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
      "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
def scgrp(a):
    i=FINE.index(a)
    return "0-44" if i<=9 else "45-64" if i<=13 else "65-84" if i<=17 else "85+"
GRPMAP={"Age at time of death, 0 to 44 years":"0-44","Age at time of death, 45 to 64 years":"45-64",
        "Age at time of death, 65 to 84 years":"65-84","Age at time of death, 85 years and over":"85+"}

# ---- StatCan full-year (52 wk) deaths by group, both sexes, 2024 & 2025 ----
sc={2024:{},2025:{}}
for row in csv.DictReader(open(SC,encoding="utf-8-sig")):
    if row["GEO"]!="Canada, place of occurrence" or row["Sex"]!="Both sexes": continue
    g=GRPMAP.get(row["Age at time of death"]); y=int(row["REF_DATE"][:4])
    if g is None or y not in (2024,2025): continue
    v=row["VALUE"]
    if v not in ("","..","..."): sc[y][g]=sc[y].get(g,0.0)+float(v)
growth={g: sc[2025][g]/sc[2024][g] for g in ["0-44","45-64","65-84","85+"]}

# ---- HMD 2024 fine deaths + current 2025 rows ----
M=list(csv.DictReader(open(MASTER)))
FN=list(M[0].keys())
def cell(y,a): return next((r for r in M if r["location"]==LOC and r["year"]==y and r["age"]==a),None)
hmd24={a:float(cell("2024",a)["D_T"]) for a in FINE}
new25={a: hmd24[a]*growth[scgrp(a)] for a in FINE}

print("StatCan group growth 2024->2025:", {g:round(growth[g],4) for g in growth})
old_tot=sum(float(cell("2025",a)["D_T"]) for a in FINE); new_tot=sum(new25.values())
print(f"CAN 2025 total  OLD={old_tot:,.0f}  NEW={new_tot:,.0f}  ({new_tot-old_tot:+,.0f}, {100*(new_tot/old_tot-1):+.2f}%)")

# ---- backup + write (update D_T, M_T; keep P_T) ----
bk=os.path.join(ROOT,"data","old",f"master_preCAN2025complete_{datetime.date(2026,8,7)}.csv")
os.makedirs(os.path.dirname(bk),exist_ok=True); shutil.copy2(MASTER,bk); print("backup ->",bk)
n=0
for r in M:
    if r["location"]==LOC and r["year"]==YEAR and r["age"] in FINE:
        PT=float(r["P_T"]); DT=new25[r["age"]]
        r["D_T"]=f"{DT:.1f}"; r["M_T"]=f"{DT/PT:.8f}"; r["source"]=SRC; n+=1
with open(MASTER,"w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=FN); w.writeheader(); w.writerows(M)
print(f"updated {n} CAN {YEAR} fine-band rows in master (source {SRC})")
