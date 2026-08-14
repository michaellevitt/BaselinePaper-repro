#!/usr/bin/env python3
"""Every table in the manuscript, with the manuscript's own numbering, in one workbook:
docs/PaperA_manuscript_tables.xlsx

  Table 1  Per-population pre-pandemic trend patterns (slope 2019, slope-of-slopes, w, STTa slopes)
  Table 2  Combined excess deaths + pooled P-score, by period x age group x model
  Table 3  Excess mortality, more- vs less-vulnerable populations, by period x model
  Table S1 Per-population excess P-score, all ages (3 periods x 4 models + mean of 4)
  Table S2 ... 65 and over
  Table S3 ... under 65
  Table S4 Excess deaths by single year and model, all ages, pooled
  Table S5 Effect of age adjustment (age-resolution) on the 2020-2025 total
  Table S6 Age-band-specific vs all-age slopes, 10 largest countries
  Table S7 Vulnerability indicators and their sources

Reads the engine outputs in output/ plus data/vuln_covariates_38_v3.csv.
Run the engine (model_option1_periods.py), build_option1_age_resolution.py and
build_table_S6_agespecific.py first — run_all.sh does this in order."""
import os, csv, numpy as np
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTD=os.path.join(ROOT,"output"); DATA=os.path.join(ROOT,"data"); DOCS=os.path.join(ROOT,"docs")
HDR=PatternFill("solid",fgColor="1F4E79"); WHT=Font(color="FFFFFF",bold=True)
SUM=PatternFill("solid",fgColor="FFF2CC")
thin=Side(style="thin",color="BFBFBF"); BORD=Border(left=thin,right=thin,top=thin,bottom=thin)
CEN=Alignment(horizontal="center",wrap_text=True)
MODELS=["Fa","TTa","STTa","STTa+"]
GROUP_LABEL={"All":"All ages","GE65":"65+ years","LT65":"<65 years"}
PERIODS=[("2020-2025","2020-2025"),("2020-2023","2020-2023"),("2024-2025","2024-2025")]

def _f(x):
    try: return float(x)
    except (TypeError,ValueError): return None
wb=Workbook(); wb.remove(wb.active)
def sheet(name,header,rows,widths=None,sumrows=0):
    ws=wb.create_sheet(name); ws.append(header)
    for j in range(1,len(header)+1):
        c=ws.cell(1,j); c.fill=HDR; c.font=WHT; c.alignment=CEN; c.border=BORD
    for r in rows:
        ws.append(r); rn=ws.max_row
        for j in range(1,len(header)+1):
            c=ws.cell(rn,j); c.border=BORD
            if j>1: c.alignment=Alignment(horizontal="center")
        if sumrows and rn>len(rows)+1-sumrows:
            for j in range(1,len(header)+1): ws.cell(rn,j).fill=SUM; ws.cell(rn,j).font=Font(bold=True)
    for i,wd in enumerate(widths or [],1): ws.column_dimensions[get_column_letter(i)].width=wd
    ws.freeze_panes="A2"; return ws

# ---------- load engine outputs ----------
par=list(csv.DictReader(open(os.path.join(OUTD,"option1_country_params.csv"))))
per=list(csv.DictReader(open(os.path.join(OUTD,"model_option1_periods.csv"))))
pool=list(csv.DictReader(open(os.path.join(OUTD,"pooled_option1_by_period.csv"))))
cov=list(csv.DictReader(open(os.path.join(DATA,"vuln_covariates_38_v3.csv"))))
NAME={r["country"]:r["name"] for r in par}
GRP={r["code"]:r["group"] for r in cov}
E={(r["country"],r["group"],r["model"],r["period"]):(float(r["O"]),float(r["E"])) for r in per}
POOL={(r["period"],r["group"],r["model"]):r for r in pool}
def P(c,g,m,p):
    O,Ee=E[(c,g,m,p)]; return 100*(O-Ee)/Ee

# ---------- Table 1 ----------
def ci(v,lo,hi): return f"{v:+.2f} ({lo:+.2f}, {hi:+.2f})" if None not in (v,lo,hi) else (f"{v:+.2f}" if v is not None else "")
rows=[]; IS=[];SOS=[];POP=[];Wt=[]
for r in par:
    is19=_f(r["islope2019"]); sos=_f(r["sos"]); pop=_f(r["pop2019"]); wv=_f(r["w"])
    rows.append([r["name"], round(pop/1e6,2) if pop else "", round(wv,3) if wv is not None else "",
        ci(is19,_f(r["islope2019_lo"]),_f(r["islope2019_hi"])), ci(sos,_f(r["sos_lo"]),_f(r["sos_hi"])),
        f"{100*_f(r['beta_shrunk']):+.2f}" if r["beta_shrunk"] else "",
        f"{200*_f(r['gamma_shrunk']):+.3f}" if r["gamma_shrunk"] else ""])
    if is19 is not None: IS.append(is19);SOS.append(sos);POP.append(pop);Wt.append(wv)
IS=np.array(IS);SOS=np.array(SOS);POP=np.array(POP)
rows.append(["Population-weighted mean","","",f"{np.average(IS,weights=POP):+.2f}",f"{np.average(SOS,weights=POP):+.3f}","",""])
rows.append(["Unweighted mean","",f"{np.mean(Wt):.3f}",f"{IS.mean():+.2f}",f"{SOS.mean():+.3f}","",""])
sheet("Table 1",["Population","Pop (M)","w","Slope in 2019 [95% CI] (%/yr)","Slope-of-slopes [95% CI] (%/yr²)","STTa slope 2019","STTa slope-of-slopes"],
      rows,[22,8,7,27,28,14,16],sumrows=2)

# ---------- Table 2 ----------
t2=[]
for pk,plabel in PERIODS:
    for g in ("All","GE65","LT65"):
        row=[plabel,GROUP_LABEL[g]]
        for m in MODELS:
            r=POOL[(pk,g,m)]
            row.append(f"{float(r['excess'])/1e6:+.2f}M ({float(r['pooled_Ppct']):+.1f}%)")
        t2.append(row)
sheet("Table 2",["Period","Age group"]+MODELS,t2,[13,12]+[18]*4)

# ---------- Table 3 ----------
MORE=[c for c in NAME if GRP.get(c)=="more"]; LESS=[c for c in NAME if GRP.get(c)=="less"]
def pooledP(codes,g,m,p):
    O=sum(E[(c,g,m,p)][0] for c in codes); Ee=sum(E[(c,g,m,p)][1] for c in codes)
    return 100*(O-Ee)/Ee
t3=[]
for pk,plabel in PERIODS:
    for m in MODELS:
        a=pooledP(MORE,"All",m,pk); b=pooledP(LESS,"All",m,pk)
        t3.append([plabel,m,f"{a:+.2f}",f"{b:+.2f}",f"{a-b:+.2f}"])
sheet("Table 3",["Period","Model",f"More vulnerable (n={len(MORE)})",f"Less vulnerable (n={len(LESS)})","Difference (pp)"],
      t3,[13,8,20,20,16])

# ---------- Tables S1 / S2 / S3 ----------
def si_table(group,name):
    order=sorted(NAME,key=lambda c: np.mean([P(c,group,m,"2020-2025") for m in MODELS]),reverse=True)
    rows=[]
    for c in order:
        row=[NAME[c]]
        for pk,_ in PERIODS:
            vals=[P(c,group,m,pk) for m in MODELS]
            row+= [round(v,1) for v in vals]+[round(float(np.mean(vals)),1)]
        rows.append(row)
    tot=["Pooled (38 populations)"]
    for pk,_ in PERIODS:
        vals=[]
        for m in MODELS:
            r=POOL[(pk,group,m)]; vals.append(float(r["pooled_Ppct"]))
        tot += [round(v,1) for v in vals]+[round(float(np.mean(vals)),1)]
    rows.append(tot)
    hdr=["Population"]+[f"{plabel} {m}" for pk,plabel in PERIODS for m in MODELS+["mean of 4"]]
    sheet(name,hdr,rows,[22]+[11]*15,sumrows=1)
si_table("All","Table S1"); si_table("GE65","Table S2"); si_table("LT65","Table S3")

# ---------- Table S4 ----------
YEARS=[str(y) for y in range(2020,2026)]
t4=[]
for m in MODELS:
    row=[m]
    for y in YEARS:
        r=POOL[(y,"All",m)]
        row.append(f"{float(r['excess'])/1e3:+.0f} ({float(r['pooled_Ppct']):+.1f})")
    t4.append(row)
sheet("Table S4",["Model"]+YEARS,t4,[10]+[15]*6)

# ---------- Table S5 ----------
ar=list(csv.reader(open(os.path.join(OUTD,"age_resolution_option1.csv"))))
hdr=ar[0]; arm={r[0]:r for r in ar[1:]}
sheet("Table S5",[hdr[0]+" (millions)"]+hdr[1:],
      [[m]+[f"{float(arm[m][i])/1e6:.2f}" for i in range(1,len(hdr))] for m in MODELS if m in arm],[22,14,16,15,17])

# ---------- Table S6 ----------
s6p=os.path.join(OUTD,"table_S6_agespecific_vs_allage.csv")
if os.path.exists(s6p):
    s6=list(csv.DictReader(open(s6p)))
    rows=[[r["name"] if r["code"]!="TOTAL" else "Total (10 countries)"]+
          [int(r[f"excess_{g}_{k}"]) for g in ("All","<65","65+") for k in ("agespecific","allage")] for r in s6]
    sheet("Table S6",["Population","All: age-specific","All: all-age","<65: age-specific","<65: all-age","65+: age-specific","65+: all-age"],
          rows,[22]+[15]*6,sumrows=1)

# ---------- Table S7 ----------
sheet("Table S7",["Code","Population","GDP per capita 2021 (US$)","Gini","Poverty (%)","Group","Trigger","Source"],
      [[r["code"],r["name"],r["gdp_pc_2021_usd"],r["gini"],r["poverty_pct"],r["group"],r.get("trigger",""),r.get("source","")]
       for r in sorted(cov,key=lambda r:r["name"])],[7,20,16,8,10,10,14,40])

os.makedirs(DOCS,exist_ok=True)
out=os.path.join(DOCS,"PaperA_manuscript_tables.xlsx"); wb.save(out)
print("wrote",out)
print("sheets:",[ws.title for ws in wb.worksheets])
print("\nTable 2, 2020-2025 all ages:", t2[0][2:])
print("Table 3, 2020-2025 Fa      :", t3[0])
print("Table S4, Fa               :", t4[0][1:])
