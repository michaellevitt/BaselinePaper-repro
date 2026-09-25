"""
Option-1 per-period engine, AGE-SPECIFIC TRENDS as the primary analysis (v2, 2026-09-19).

Identical to model_option1_periods.py except for where TTa / STTa / STTa+ get their trend.
The published engine used ONE (beta, gamma) per country, estimated from the age-standardised
rate, applied to all twenty bands.  Here each band carries its own two-step trend, per John's
18 Sep 2026 request: "we defeat ourselves when we say that age adjustment is so important but
then ignore it in the main analyses."

  Fa   : flat 2017-2019 mean rate per band            UNCHANGED (already per band)
  Ta   : linear-in-rate OLS 2015-2019 per band        UNCHANGED (already per band)
  TTa  : band's OWN two-step (islope@2019, sos)
  STTa : band's shrunk trend, shrinkage computed WITHIN band
  STTa+: band's shrunk ANCHORED trend
  ATTa : unchanged country-level global (retained for provenance, not in the paper)

Trend parameters are read from output/agespecific_band_params_v1.csv, written by
code/model_agespecific_trends_v1.py, which must run first.  This mirrors how the published
engine reads its country trends from slope_of_slopes_CI.csv.  That CSV stores 4 decimals, so
reproducing these numbers exactly requires the same rounding.

option1_country_params.csv keeps its country-level columns unchanged, so Figure 1, Figure S1 and
Table 1 still describe the age-standardised trend and its between-country heterogeneity, and
gains two columns, g19_bandwmean and sos_bandwmean: the baseline-deaths-weighted mean of the
band trends actually used, which is the age-specific model's country-level summary.

Writes the SAME filenames as the published engine, so every downstream figure and table picks
these up without modification.
"""
import os, sys, csv, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); OUTD=os.path.join(ROOT,"output")
sys.path.insert(0, HERE)
from excess_anchor_window import load, cell_for
FINE=["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
      "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
def astart(a): return 0 if a=="0" else (1 if a=="1-4" else int(a.split("-")[0].replace("+","")))
GROUPS={"All":FINE,"GE65":[a for a in FINE if astart(a)>=65],"LT65":[a for a in FINE if astart(a)<65]}
ESP={"0":1000,"1-4":4000,"5-9":5500,"10-14":5500,"15-19":5500,"20-24":6000,"25-29":6000,
     "30-34":6500,"35-39":7000,"40-44":7000,"45-49":7000,"50-54":7000,"55-59":6500,"60-64":6000,
     "65-69":5500,"70-74":5000,"75-79":4000,"80-84":2500,"85-89":1500,"90+":1000}; ESPsum=sum(ESP.values())
W=5; RELIABLE={"BGR":2015,"CHL":2011,"HRV":2007,"EST":2005,"HUN":2006,"LVA":2007,"LTU":2006,"POL":2008,"SVK":2006}
PERIODS={"2020-2025":[2020,2021,2022,2023,2024,2025],"2020-2023":[2020,2021,2022,2023],"2024-2025":[2024,2025]}
PERIODS.update({str(y):[y] for y in range(2020,2026)})   # single years, for the by-year SI table
TAU=float("inf"); CLIP=(0.4,2.5); MINYRS=6; MODELS=["Fa","Ta","TTa","ATTa","STTa","STTa+"]  # TAU=inf -> UNIFORM weighting (recency dropped per John 2026-07-26; was 6.0)
data,name=load()
def DP(l,y,a): return cell_for(data,l,y,a,"T")
def has_year(l,y): return all(DP(l,y,a) is not None for a in FINE)
def has_base(l):
    se=max(2003+(W-1),RELIABLE.get(l,2003)+(W-1)); return all(DP(l,y,a) is not None for y in range(se-(W-1),2020) for a in FINE)
PC=[r[0] for r in csv.reader(open(os.path.join(OUTD,"slope_of_slopes_CI.csv"))) if r and r[0]!="country" and not r[0].startswith(("WEIGHTED","UNWEIGHTED","SD"))]
countries=[l for l in sorted(PC) if has_base(l)]
def baseyrs(l): return list(range(max(2003,RELIABLE.get(l,2003)),2020))
def esp_series(l,ys): return np.array([sum(ESP[a]*DP(l,y,a)[0]/DP(l,y,a)[1] for a in FINE)/ESPsum for y in ys])
def esp_logvar(l,y):     # Poisson sampling variance of the ESP-standardised log rate at year y (~ 1/deaths); S cancels
    num=sum(ESP[a]**2*DP(l,y,a)[0]/DP(l,y,a)[1]**2 for a in FINE); den=sum(ESP[a]*DP(l,y,a)[0]/DP(l,y,a)[1] for a in FINE)
    return num/den**2
def base_deaths(l,ys): return sum(DP(l,y,a)[0] for y in ys for a in FINE)

# ---- per-country standardised trend (beta,gamma) from the TWO-STEP moving-slope construction (John, 2026-07-27) ----
# beta = islope@2019/100 and gamma = slope-of-slopes/200 come from the OLS of the 5-yr moving slopes of the
# ESP-standardised rate against window-END year (code/slope_of_slopes_ci.py -> slope_of_slopes_CI.csv; pop-weighted
# mean islope -0.79, sos 0.107).  The sos sampling variance (from that OLS's 95% CI) drives the shrinkage weight
# w = 1/sqrt(VR); STTa = w*own + (1-w)*global.  Applied at country level to every band (via the band's 2019 level).
def _f(x):
    try: return float(x)
    except (TypeError,ValueError): return None
TS={}
for r in csv.reader(open(os.path.join(OUTD,"slope_of_slopes_CI.csv"))):   # 5-10: unanch islope2019/sos(+CI); 11,14: anchored islope@2025, sos
    if not r or r[0]=="country" or r[0].startswith(("WEIGHTED","UNWEIGHTED","SD")): continue
    TS[r[0]]=dict(is19=_f(r[5]),is_lo=_f(r[6]),is_hi=_f(r[7]),sos=_f(r[8]),sos_lo=_f(r[9]),sos_hi=_f(r[10]),
                  a_is2025=_f(r[11]),a_sos=_f(r[14]))
Z=1.959964
cpar={}
for l in countries:
    pop=sum(DP(l,2019,a)[1] for a in FINE); D=base_deaths(l,baseyrs(l)); ts=TS.get(l,{})
    if ts.get("is19") is not None and ts.get("sos") is not None:
        b=ts["is19"]/100.0; g=ts["sos"]/200.0
        ai=ts.get("a_is2025"); asos=ts.get("a_sos")                       # SLOPE-anchored: adds the {2019,2024,2025} return slope
        ba=(ai-6.0*asos)/100.0 if (ai is not None and asos is not None) else np.nan   # anchored slope@2019 = slope@2025 - 6*sos
        ga=asos/200.0 if asos is not None else np.nan
        cpar[l]=dict(b=b,g=g,ba=ba,ga=ga,vg=1.0/D,is19=ts["is19"],is_lo=ts["is_lo"],is_hi=ts["is_hi"],   # vg=1/deaths -> w=sqrt(D/Dmax)
                     sos=ts["sos"],sos_lo=ts["sos_lo"],sos_hi=ts["sos_hi"],pop=pop,D=D)
    else:
        cpar[l]=dict(b=np.nan,g=np.nan,ba=np.nan,ga=np.nan,vg=np.inf,is19=np.nan,is_lo=np.nan,is_hi=np.nan,
                     sos=np.nan,sos_lo=np.nan,sos_hi=np.nan,pop=pop,D=D)
valid=[l for l in countries if np.isfinite(cpar[l]["b"])]
wp=np.array([cpar[l]["pop"] for l in valid])
bG=float(np.average([cpar[l]["b"] for l in valid],weights=wp)); gG=float(np.average([cpar[l]["g"] for l in valid],weights=wp))
va=[l for l in countries if np.isfinite(cpar[l]["ba"])]; wa=np.array([cpar[l]["pop"] for l in va])
baG=float(np.average([cpar[l]["ba"] for l in va],weights=wa)); gaG=float(np.average([cpar[l]["ga"] for l in va],weights=wa))
vmin=min(cpar[l]["vg"] for l in valid)                                    # smallest slope-of-slopes sampling variance
for l in countries:
    p=cpar[l]
    if np.isfinite(p["b"]):
        p["VR"]=p["vg"]/vmin; p["w"]=1.0/np.sqrt(p["VR"])                 # sqrt-deaths shrinkage (John #2)
        p["bs"]=p["w"]*p["b"]+(1-p["w"])*bG; p["gs"]=p["w"]*p["g"]+(1-p["w"])*gG
        if np.isfinite(p["ba"]): p["bas"]=p["w"]*p["ba"]+(1-p["w"])*baG; p["gas"]=p["w"]*p["ga"]+(1-p["w"])*gaG   # STTa+ = shrunk ANCHORED slope
        else: p["bas"],p["gas"]=baG,gaG
    else:
        p["VR"]=np.inf; p["w"]=0.0; p["bs"],p["gs"]=bG,gG; p["bas"],p["gas"]=baG,gaG
    p["wb"]=p["wg"]=p["w"]                                                # back-compat: single weight for both

# ---- AGE-SPECIFIC trends: read the per-band two-step parameters (v2) ----
BP={}
_bp_path=os.path.join(OUTD,"agespecific_band_params_v1.csv")
if not os.path.exists(_bp_path):
    sys.exit("run code/model_agespecific_trends_v1.py first: "+_bp_path+" is missing")
for r in csv.DictReader(open(_bp_path)):
    BP.setdefault(r["country"],{})[r["band"]]=dict(
        is19=_f(r["islope2019"]), sos=_f(r["sos"]),
        is19_s=_f(r["islope2019_shrunk"]), sos_s=_f(r["sos_shrunk"]),
        a_is_s=_f(r["anch_islope2025_shrunk"]), a_sos_s=_f(r["anch_sos_shrunk"]),
        D=float(r["base_deaths"]))
missing=[l for l in countries if l not in BP]
if missing: sys.exit(f"per-band trends missing for {missing}")
# population-level summary of the band trends actually used (deaths-weighted), for Table 1/Fig 1
for l in countries:
    wD=np.array([BP[l][a]["D"] for a in FINE],float)
    for src,dst in (("is19_s","g19_bw"),("sos_s","sos_bw")):
        v=np.array([BP[l][a][src] if BP[l][a][src] is not None else np.nan for a in FINE],float)
        m=np.isfinite(v)&(wD>0)
        cpar[l][dst]=float(np.average(v[m],weights=wD[m])) if m.any() else np.nan
_CLIPHIT=[0,0]

# ---- per-band baselines: cache band WLS logquad coeffs + Fa/Ta rates ----
band={}; Fa_rate={}; Ta_rate={}
for l in countries:
    ys=baseyrs(l); t=np.array(ys,float)-2019.0; om=np.exp(-(2019-np.array(ys,float))/TAU); band[l]={}
    Fa_rate[l]={}; Ta_rate[l]={}
    for a in FINE:
        r=np.array([DP(l,y,a)[0]/DP(l,y,a)[1] for y in ys])
        band[l][a]=np.polyfit(t,np.log(np.maximum(r,1e-9)),2,w=np.sqrt(om))     # [g,b,alpha]
        Fa_rate[l][a]=sum(DP(l,y,a)[0] for y in (2017,2018,2019))/sum(DP(l,y,a)[1] for y in (2017,2018,2019))
        Ta_rate[l][a]=np.polyfit(np.arange(2015,2020,dtype=float)-2012,[DP(l,y,a)[0]/DP(l,y,a)[1] for y in range(2015,2020)],1)

def rate(l,a,model,y):
    """v2: TTa/STTa/STTa+ take the BAND's own two-step trend, not the country's."""
    t=y-2019.0; c=band[l][a]; alpha=c[2]; p=BP[l][a]
    if model=="Fa": return Fa_rate[l][a]
    if model=="Ta": return float(np.polyval(Ta_rate[l][a],y-2012))
    if model=="ATTa": b,g=bG,gG                              # unchanged country-level global
    elif model=="STTa+":
        ai,asos=p["a_is_s"],p["a_sos_s"]                     # anchored slope@2025 -> back to 2019
        if ai is None or asos is None: b,g=cpar[l]["bas"],cpar[l]["gas"]
        else: b,g=(ai-6.0*asos)/100.0, asos/200.0
    elif model=="TTa":
        b,g=(p["is19"]/100.0, p["sos"]/200.0) if p["is19"] is not None else (p["is19_s"],p["sos_s"])
        if p["is19"] is None: b,g=(b/100.0 if b is not None else bG),(g/200.0 if g is not None else gG)
    else:
        b,g=(p["is19_s"]/100.0, p["sos_s"]/200.0) if p["is19_s"] is not None else (bG,gG)
    m=np.exp(b*t+g*t*t)
    _CLIPHIT[0]+=1
    if m<CLIP[0] or m>CLIP[1]: _CLIPHIT[1]+=1
    return np.exp(alpha)*min(max(m,CLIP[0]),CLIP[1])

# ---- compute per period/group/model ----
recs=[]
for l in countries:
    for pname,pyears in PERIODS.items():
        yrs=[y for y in pyears if has_year(l,y)]
        for g,bands in GROUPS.items():
            for m in MODELS:
                ps=[]; O=E=0.0
                for y in yrs:
                    o=sum(DP(l,y,a)[0] for a in bands); e=sum(rate(l,a,m,y)*DP(l,y,a)[1] for a in bands)
                    O+=o; E+=e; ps.append(100*(o-e)/e if e else np.nan)
                recs.append([l,g,m,pname,len(yrs),round(O),round(E),round(O-E),round(float(np.nanmean(ps)),3)])
with open(os.path.join(OUTD,"model_option1_periods.csv"),"w",newline="") as f:
    w=csv.writer(f); w.writerow(["country","group","model","period","n_years","O","E","excess","pscore_meanann"]); w.writerows(recs)

pool=[]
for pname in PERIODS:
    for g in GROUPS:
        for m in MODELS:
            rr=[x for x in recs if x[3]==pname and x[1]==g and x[2]==m]
            O=sum(x[5] for x in rr); E=sum(x[6] for x in rr)
            pool.append([pname,g,m,len(rr),round(O),round(E),round(O-E),round(100*(O-E)/E,3),round(float(np.mean([x[8] for x in rr])),3)])
with open(os.path.join(OUTD,"pooled_option1_by_period.csv"),"w",newline="") as f:
    w=csv.writer(f); w.writerow(["period","group","model","n_countries","O","E","excess","pooled_Ppct","mean_meanann"]); w.writerows(pool)

def _fm(x,nd=3): return f"{x:.{nd}f}" if (x is not None and np.isfinite(x)) else ""
with open(os.path.join(OUTD,"option1_country_params.csv"),"w",newline="") as f:
    w=csv.writer(f); w.writerow(["country","name","beta_i","gamma_i",
        "islope2019","islope2019_lo","islope2019_hi","sos","sos_lo","sos_hi",   # two-step table cols (intercept-slope, sos, each +/-95% CI)
        "beta_shrunk","gamma_shrunk","beta_anch_shrunk","gamma_anch_shrunk",    # STTa (unanchored) and STTa+ (slope-anchored) shrunk slopes
        "g19_pct","sos_pct","w","VR","pop2019","base_deaths","g19_bandwmean","sos_bandwmean"])
    for l in countries:
        p=cpar[l]; fin=np.isfinite(p["b"])
        w.writerow([l,name.get(l,l), _fm(p["b"],6),_fm(p["g"],7),
            _fm(p["is19"]),_fm(p["is_lo"]),_fm(p["is_hi"]),_fm(p["sos"]),_fm(p["sos_lo"]),_fm(p["sos_hi"]),
            _fm(p["bs"],6),_fm(p["gs"],7),_fm(p["bas"],6),_fm(p["gas"],7),
            _fm(p["is19"]),_fm(p["sos"],4), f"{p['w']:.4f}", f"{p['VR']:.1f}" if fin else "",
            round(p["pop"]), round(p["D"]), _fm(p.get("g19_bw"),4), _fm(p.get("sos_bw"),4)])

print(f"\n[v2 AGE-SPECIFIC] baseline multiplier hit the clip in {_CLIPHIT[1]:,} of {_CLIPHIT[0]:,} band-years")
print(f"wrote model_option1_periods.csv ({len(recs)} rows), pooled_option1_by_period.csv, option1_country_params.csv")
print(f"global trend: g19={100*bG:+.2f}%/yr  sos={200*gG:+.3f}%/yr^2\n")
print("POOLED all-ages excess by period x model (two-step trend, slope-anchored STTa+):")
print(f"  {'period':10}"+"".join(f"{m:>16}" for m in MODELS))
for pname in PERIODS:
    cells=[]
    for m in MODELS:
        row=next(x for x in pool if x[0]==pname and x[1]=="All" and x[2]==m)
        cells.append(f"{row[6]/1e6:>7.2f}M({row[7]:+.1f}%)")
    print(f"  {pname:10}"+"".join(f"{c:>16}" for c in cells))
