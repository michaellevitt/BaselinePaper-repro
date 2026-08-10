#!/usr/bin/env python3
"""
Assemble the newly-sourced fine 5-year-band mortality for the 5 re-added populations into the master
20-band structure (0,1-4,5-9,...,85-89,90+; sex=Total), and APPEND the new years to
data/master_5x1_DPM_90plus.csv under the primary codes. New years:
  AUS 2022-2024 (ABS Data API); CAN 2024 (StatCan); ISL 2024-2025 (Eurostat, already in ISL_EUROSTAT);
  ISR 2017-2023 deaths / 2018-2024 pop (Israel CBS; 2017 pop interpolated from HMD-2016 & CBS-2018);
  NZL_NP 2022-2025 (Stats NZ, RR3-rounded, registration year).
Top bands 90-94/95-99/100+ (Canada, NZ) and 90-94/95+ (Israel) collapse to 90+. Combined-0-4 population
(Canada, NZ) is split into 0 and 1-4 by each country's most recent master ratio. Deaths that already give
0 and 1-4 separately are kept. Validates all-ages death sums against the agents' reported totals.
Writes rows to master and prints a validation table. Iceland copies ISL_EUROSTAT -> ISL.
"""
import os, csv, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); DATD=os.path.join(ROOT,"data")
FINE=["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
ASTART={b:(0 if b=="0" else 1 if b=="1-4" else int(b.split("-")[0].replace("+",""))) for b in FINE}
MASTER=os.path.join(DATD,"master_5x1_DPM_90plus.csv")

# ---- load master (for ratios + ISL_EUROSTAT + ISR 2016) ----
rows=list(csv.DictReader(open(MASTER))); cols=rows[0].keys()
def get(loc,y,a,col):
    for r in rows:
        if r["location"]==loc and r["year"]==str(y) and r["age"]==a: return r[col]
    return None
def ratio04(loc,y):  # split 0-4 pop into 0 and 1-4 by this loc/year
    p0=float(get(loc,y,"0","P_T")); p14=float(get(loc,y,"1-4","P_T")); return p0/(p0+p14)

# ---- gather per (loc,year): dict age->(D,P) ----
DATA={}   # (loc,year) -> {age:(D,P)}
def top90(d):  # collapse any 90-94/95-99/95+/100+ into 90+
    keep={}; nz=0
    for a,v in d.items():
        if a in ("90-94","95-99","100+","95+","90+"): nz+=v
        else: keep[a]=v
    keep["90+"]=nz; return keep

# Australia (JSON) — bands already 0..90+
AJ=json.load(open("/private/tmp/claude-501/-Users-levitt-Dropbox-win1-DB-NewProjects25-mortality-org-HMD-Excess-Death/e78eb5c5-d2a8-4782-876e-ecd463b8a89a/scratchpad/AUS_extract.json"))
for y in ["2022","2023","2024"]:
    dd=AJ["deaths_registrationYear_persons"][y]; pp=AJ["population_ERP_30June_persons"][y]
    DATA[("AUS",int(y))]={a:(dd[a],pp[a]) for a in FINE}

# Iceland — copy ISL_EUROSTAT rows (already 20 bands)
for y in (2024,2025):
    DATA[("ISL",y)]={a:(float(get("ISL_EUROSTAT",y,a,"D_T")),float(get("ISL_EUROSTAT",y,a,"P_T"))) for a in FINE}

# Canada 2024 (StatCan) — deaths 0 & 1-4 split already; pop 0-4 combined -> split; top 90-94/95-99/100+ -> 90+
CAN_D={"0":1661,"1-4":245,"5-9":190,"10-14":246,"15-19":794,"20-24":1462,"25-29":2091,"30-34":2961,"35-39":3375,"40-44":4224,"45-49":5026,"50-54":7179,"55-59":11220,"60-64":19238,"65-69":26652,"70-74":33519,"75-79":42582,"80-84":46119,"85-89":48531,"90-94":43165,"95-99":21006,"100+":5293}
CAN_P0_4=1874025; CAN_P={"5-9":2151780,"10-14":2232438,"15-19":2313964,"20-24":2741716,"25-29":3031487,"30-34":3145235,"35-39":2950432,"40-44":2796675,"45-49":2538433,"50-54":2444797,"55-59":2492843,"60-64":2712778,"65-69":2442866,"70-74":1973110,"75-79":1535283,"80-84":963055,"85-89":558478,"90-94":270684,"95-99":80331,"100+":11919}
r=ratio04("CAN",2023); CAN_P["0"]=CAN_P0_4*r; CAN_P["1-4"]=CAN_P0_4*(1-r)
CAN_Dc=top90(CAN_D); CAN_Pc=top90(CAN_P)
DATA[("CAN",2024)]={a:(CAN_Dc[a],CAN_Pc[a]) for a in FINE}

# New Zealand 2022-2025 (Stats NZ, RR3) — deaths 0 & 1-4 split, top ->90+; pop 0-4 combined -> split, 90+ ok
NZ_D={2022:{"0":204,"1-4":42,"5-9":21,"10-14":36,"15-19":108,"20-24":177,"25-29":216,"30-34":255,"35-39":288,"40-44":378,"45-49":609,"50-54":993,"55-59":1512,"60-64":2148,"65-69":2718,"70-74":3681,"75-79":4821,"80-84":5988,"85-89":6228,"90-94":5373,"95-99":2316,"100+":456},
2023:{"0":201,"1-4":51,"5-9":18,"10-14":36,"15-19":144,"20-24":186,"25-29":216,"30-34":267,"35-39":291,"40-44":423,"45-49":585,"50-54":936,"55-59":1455,"60-64":2142,"65-69":2718,"70-74":3606,"75-79":4905,"80-84":5997,"85-89":6036,"90-94":5142,"95-99":2145,"100+":381},
2024:{"0":339,"1-4":36,"5-9":18,"10-14":39,"15-19":123,"20-24":159,"25-29":198,"30-34":231,"35-39":309,"40-44":429,"45-49":567,"50-54":978,"55-59":1341,"60-64":2103,"65-69":2664,"70-74":3621,"75-79":4794,"80-84":6015,"85-89":6174,"90-94":4980,"95-99":2205,"100+":402},
2025:{"0":276,"1-4":33,"5-9":24,"10-14":39,"15-19":126,"20-24":156,"25-29":171,"30-34":252,"35-39":321,"40-44":426,"45-49":525,"50-54":900,"55-59":1320,"60-64":2022,"65-69":2736,"70-74":3534,"75-79":5154,"80-84":5886,"85-89":6210,"90-94":4785,"95-99":2184,"100+":402}}
NZ_P04={2022:302730,2023:302310,2024:300900,2025:297580}
NZ_Prest={2022:{"5-9":321110,"10-14":339750,"15-19":319910,"20-24":319980,"25-29":349530,"30-34":382900,"35-39":346550,"40-44":314900,"45-49":312350,"50-54":330060,"55-59":315110,"60-64":300050,"65-69":254100,"70-74":216410,"75-79":157630,"80-84":107870,"85-89":57180,"90+":33610},
2023:{"5-9":322920,"10-14":344740,"15-19":333950,"20-24":323830,"25-29":355040,"30-34":402430,"35-39":371090,"40-44":333850,"45-49":311650,"50-54":332190,"55-59":312150,"60-64":305950,"65-69":260010,"70-74":217520,"75-79":168950,"80-84":109300,"85-89":58850,"90+":33330},
2024:{"5-9":326870,"10-14":349130,"15-19":345410,"20-24":328590,"25-29":350810,"30-34":409040,"35-39":390290,"40-44":348690,"45-49":313880,"50-54":331010,"55-59":312170,"60-64":309380,"65-69":266640,"70-74":221640,"75-79":177320,"80-84":113360,"85-89":61260,"90+":33660},
2025:{"5-9":325440,"10-14":348240,"15-19":354350,"20-24":327760,"25-29":344080,"30-34":401820,"35-39":401290,"40-44":357420,"45-49":316820,"50-54":326640,"55-59":313260,"60-64":309870,"65-69":272530,"70-74":225340,"75-79":185530,"80-84":117040,"85-89":65300,"90+":34400}}
rN=ratio04("NZL_NP",2021)
for y in (2022,2023,2024,2025):
    P={"0":NZ_P04[y]*rN,"1-4":NZ_P04[y]*(1-rN)}; P.update(NZ_Prest[y])
    Dc=top90(NZ_D[y])
    DATA[("NZL_NP",y)]={a:(Dc[a],P[a]) for a in FINE}

# Israel 2017-2023 deaths / 2018-2024 pop; use 2018-2023 (both); interpolate 2017 pop (HMD 2016 & CBS 2018)
ic=list(csv.reader(open(os.path.join(DATD,"ISR_CBS_2017_2024","ISRAEL_CBS_deaths_pop_by_5yr_band_2017_2024.csv"))))
hdr=ic[0]; yrs={hdr[i]:i for i in range(2,len(hdr))}
ID={}; IP={}
for r in ic[1:]:
    metric,band=r[0],r[1]
    for y,i in yrs.items():
        if r[i]=="": continue
        (ID if metric=="deaths" else IP).setdefault(int(y),{})[band]=float(r[i])
def isr_collapse(d):  # 90-94 + 95+ -> 90+
    out={a:d.get(a,0.0) for a in FINE if a!="90+"}; out["90+"]=d.get("90-94",0)+d.get("95+",0); return out
# ISR 2017 pop = interpolate HMD-2016 (master) & CBS-2018
isr17p={}
for a in FINE:
    p16=get("ISR",2016,a,"P_T")
    if a=="90+":
        p18=IP[2018].get("90-94",0)+IP[2018].get("95+",0)
    else: p18=IP[2018].get(a,0)
    isr17p[a]=(float(p16)+p18)/2 if p16 is not None else p18
for y in range(2018,2024):
    Dc=isr_collapse(ID[y]); Pc=isr_collapse(IP[y]); DATA[("ISR",y)]={a:(Dc[a],Pc[a]) for a in FINE}
Dc17=isr_collapse(ID[2017]); DATA[("ISR",2017)]={a:(Dc17[a],isr17p[a]) for a in FINE}

# ---- validate all-ages death sums vs agent-reported totals ----
REPORT={("AUS",2022):190939,("AUS",2023):183131,("AUS",2024):187268,("CAN",2024):326779,
 ("ISR",2020):49006,("ISR",2021):50984,("ISR",2022):52055,("ISR",2023):49977}
print(f"{'loc':7}{'year':6}{'ΣD':>10}{'report':>10}{'ΣP':>13}")
for (loc,y),d in sorted(DATA.items()):
    sd=sum(v[0] for v in d.values()); sp=sum(v[1] for v in d.values()); rep=REPORT.get((loc,y),"")
    flag="" if rep=="" else ("  OK" if abs(sd-rep)<=20 else f"  ** off by {sd-rep:.0f}")
    print(f"{loc:7}{y:<6}{sd:10.0f}{str(rep):>10}{sp:13.0f}{flag}")

# ---- append new rows to master ----
NAMES={"AUS":"Australia","CAN":"Canada","ISL":"Iceland","ISR":"Israel","NZL_NP":"New Zealand"}
existing={(r["location"],r["year"],r["age"]) for r in rows}
new=[]
for (loc,y),d in DATA.items():
    for a in FINE:
        if (loc,str(y),a) in existing: continue
        D,P=d[a]; new.append({"location":loc,"location_name":NAMES[loc],"source":"NATL_2025_build","year":str(y),
            "age":a,"age_start":str(ASTART[a]),"D_F":"","D_M":"","D_T":f"{D:.0f}","P_F":"","P_M":"","P_T":f"{P:.2f}",
            "M_F":"","M_M":"","M_T":f"{1000*D/P:.4f}" if P else ""})
with open(MASTER,"a",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(cols)); w.writerows(new)
print(f"\nappended {len(new)} rows to master ({len(DATA)} location-years)")
