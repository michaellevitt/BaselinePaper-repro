#!/usr/bin/env python3
"""
Israel 2024 & 2025 fine-band deaths + population, appended to master_5x1_DPM_90plus.csv (2026-07-18).
Michael wants Israel in Figure 1 and Table S4, i.e. in the full 10-model analysis, which needs its
2024/2025 fine-band endpoints (the anchors).  Israel's national fine-band DEATHS stop at 2023, so:

  DEATHS 2024/2025  = STMF broad-band annual totals (0-14,15-64,65-74,75-84,85+), split into the 20 fine
                      bands using Israel's 2023 CBS within-broad-band death shares (the latest observed,
                      war-era age structure).  STMF agrees with the national 2023 series to ~1%.
  POP    2024       = CBS mid-year population (already published to 2024).
  POP    2025       = CBS 2024 extrapolated one year by each band's 2023->2024 growth ratio.

CAVEAT (reported to Michael, not silently buried): Israel's 2024-2025 mortality carries a large
non-pandemic component from the October-2023 war, concentrated in young adults, so its anchored
baselines absorb war deaths.  This is an approximation flagged in the paper's Data paragraph.

Appends 40 Total-only rows (ISR 2024, 2025 x 20 bands).  Backs up the master first.
"""
import os, csv, shutil
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
MASTER=os.path.join(ROOT,"data","master_5x1_DPM_90plus.csv")
CBS=os.path.join(ROOT,"data","ISR_CBS_2017_2024","ISRAEL_CBS_deaths_pop_by_5yr_band_2017_2024.csv")

BANDS=[("0",0),("1-4",1),("5-9",5),("10-14",10),("15-19",15),("20-24",20),("25-29",25),("30-34",30),
       ("35-39",35),("40-44",40),("45-49",45),("50-54",50),("55-59",55),("60-64",60),("65-69",65),
       ("70-74",70),("75-79",75),("80-84",80),("85-89",85),("90+",90)]
BROAD={"0_14":[0,1,5,10],"15_64":[15,20,25,30,35,40,45,50,55,60],"65_74":[65,70],"75_84":[75,80],"85p":[85,90]}
def broad_of(a):
    for b,starts in BROAD.items():
        if a in starts: return b
    raise ValueError(a)
STMF={2024:{"0_14":817,"15_64":7931,"65_74":8685,"75_84":13533,"85p":20387},
      2025:{"0_14":760,"15_64":7743,"65_74":8652,"75_84":13728,"85p":20128}}

# ---- read CBS: deaths 2023 (for within-band shares); pop 2023 & 2024 (collapse 90-94/95+ -> 90+) ----
raw={}
for r in csv.DictReader(open(CBS)):
    raw[(r["metric"],r["age_band"])]=r
def collapse(metric, year):
    """return {age_start: value}, folding 90-94 and 95+ into 90+."""
    out={}
    for lab,a in BANDS:
        if lab=="90+":
            v=0.0
            for sub in ("90-94","95+"):
                cell=raw.get((metric,sub),{}).get(str(year),"")
                v+=float(cell) if cell not in ("",None) else 0.0
            out[a]=v
        else:
            cell=raw.get((metric,lab),{}).get(str(year),"")
            out[a]=float(cell) if cell not in ("",None) else 0.0
    return out
d2023=collapse("deaths",2023)
p2023=collapse("pop_midyear",2023); p2024=collapse("pop_midyear",2024)

# within-broad-band death shares from CBS 2023
share={}
for b,starts in BROAD.items():
    tot=sum(d2023[a] for a in starts)
    for a in starts: share[a]=d2023[a]/tot

# ---- deaths 2024/2025 = STMF broad total * within-band share ----
D={2024:{},2025:{}}
for y in (2024,2025):
    for lab,a in BANDS: D[y][a]=STMF[y][broad_of(a)]*share[a]
# ---- pop: 2024 from CBS; 2025 = 2024 * (2024/2023) growth per band ----
P={2024:{a:p2024[a] for _,a in BANDS},
   2025:{a:(p2024[a]*(p2024[a]/p2023[a]) if p2023[a] else p2024[a]) for _,a in BANDS}}

# validation
for y in (2024,2025):
    ds=sum(D[y].values()); stmf_tot=sum(STMF[y].values())
    assert abs(ds-stmf_tot)<1.0, (y,ds,stmf_tot)
    print(f"  {y}: deaths sum {ds:,.0f} (STMF total {stmf_tot:,}); pop sum {sum(P[y].values()):,.0f}")

# ---- append ----
shutil.copy(MASTER, MASTER.replace(".csv",".preISR2425.bak.csv"))
with open(MASTER,"a",newline="") as f:
    w=csv.writer(f)
    for y in (2024,2025):
        for lab,a in BANDS:
            d=round(D[y][a]); p=round(P[y][a]); m=d/p*1000 if p else ""
            w.writerow(["ISR","Israel","STMF_CBS_2425",y,lab,a,"","",d,"","",p,"","",
                        (f"{m:.4f}" if m!="" else "")])
print(f"appended ISR 2024+2025 (40 rows) to {os.path.basename(MASTER)}; backup .preISR2425.bak.csv")
