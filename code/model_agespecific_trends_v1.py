#!/usr/bin/env python3
"""
model_agespecific_trends_v1.py — TTa, STTa and STTa+ with AGE-SPECIFIC trends, all 38 populations.

WHY (John, 18 Sep 2026)
  "I am increasingly intrigued that using age-adjusted analyses for TTa (but also then STTa and
  STTa+) for all countries and not just the 10 largest would be the best approach. In a way, we
  defeat ourselves when we say that age adjustment is so important but then ignore it in the main
  analyses."

  He is right about the asymmetry.  In the published models the starting LEVEL is age-specific
  (alpha_a is band a's own fitted 2019 log rate) but the TREND is not: one (beta, gamma) per
  country, taken from the two-step moving-slope fit to the ESP-standardised rate, is applied to
  all twenty bands.  Figure S2 already shows per-band trends for the 10 largest populations.
  This script fits them for every band of every population and runs the whole family on them.

CONSTRUCTION — identical conventions to the published models, applied per band
  * 5-year CENTRED moving log-slope of the band's own log death rate, at window centres
    RELIABLE+2 .. 2017 (same as slope_of_slopes_ci.py and fig_perband_trend_10countries.py).
  * OLS of those slopes on window centre, evaluated at 2019: islope@2019 and slope-of-slopes.
    beta_a = islope@2019/100, gamma_a = sos/200.
  * Anchored variant adds the {2019, 2024, 2025} return slope of that band at its centroid
    (~2022.7), evaluates at 2025, and back-computes beta as (islope@2025 - 6*sos)/100.
  * Baseline rate = exp(alpha_a) * clip(exp(beta_a*t + gamma_a*t^2), 0.4, 2.5), t = y - 2019,
    with alpha_a unchanged from the engine (band's own fitted 2019 level).

SHRINKAGE — the direct analogue of the published STTa, computed WITHIN each band
  A single band carries far fewer deaths than a whole country, so its own two-step trend is much
  noisier: a band with a few dozen deaths cannot support a quadratic.  Published STTa shrinks
  country i toward the 38-population global trend with weight sqrt(D_i / D_max).  Here band a of
  country i is shrunk toward the population-weighted global trend FOR THAT BAND, with weight
  sqrt(D_ia / max_i D_ia).  So the smallest bands borrow most, exactly as the smallest countries
  do in the published model.  The alternative target, each country's own all-ages trend, is
  reported as a sensitivity (STTa_asC) because it is the other defensible choice and John may
  prefer it.

Writes
  output/model_agespecific_periods_v1.csv     long, models TTa_as / STTa_as / STTa+_as / STTa_asC
  output/pooled_agespecific_by_period_v1.csv  pooled by period x group x model
  output/agespecific_band_params_v1.csv       per country x band: islope, sos, anchored, w, deaths
  output/agespecific_vs_common_v1.csv         side-by-side with the published common-trend models
"""
import os, sys, csv, numpy as np
import scipy.stats as ss
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
OUTD = os.path.join(ROOT, "output")
sys.path.insert(0, HERE)
from excess_anchor_window import load, cell_for

FINE = ["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
        "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
def astart(a): return 0 if a == "0" else (1 if a == "1-4" else int(a.split("-")[0].replace("+","")))
GROUPS = {"All": FINE, "GE65": [a for a in FINE if astart(a) >= 65],
          "LT65": [a for a in FINE if astart(a) < 65]}
ESP = {"0":1000,"1-4":4000,"5-9":5500,"10-14":5500,"15-19":5500,"20-24":6000,"25-29":6000,
       "30-34":6500,"35-39":7000,"40-44":7000,"45-49":7000,"50-54":7000,"55-59":6500,
       "60-64":6000,"65-69":5500,"70-74":5000,"75-79":4000,"80-84":2500,"85-89":1500,"90+":1000}
W = 5
RELIABLE = {"BGR":2015,"CHL":2011,"HRV":2007,"EST":2005,"HUN":2006,"LVA":2007,"LTU":2006,
            "POL":2008,"SVK":2006}
PERIODS = {"2020-2025":[2020,2021,2022,2023,2024,2025], "2020-2023":[2020,2021,2022,2023],
           "2024-2025":[2024,2025]}
PERIODS.update({str(y): [y] for y in range(2020, 2026)})
CLIP = (0.4, 2.5)
MODELS = ["TTa_as", "STTa_as", "STTa+_as", "STTa_asC"]

data, name = load()
def DP(l, y, a): return cell_for(data, l, y, a, "T")
def has_year(l, y): return all(DP(l, y, a) is not None for a in FINE)
def has_base(l):
    se = max(2003 + (W-1), RELIABLE.get(l, 2003) + (W-1))
    return all(DP(l, y, a) is not None for y in range(se-(W-1), 2020) for a in FINE)
PC = [r[0] for r in csv.reader(open(os.path.join(OUTD, "slope_of_slopes_CI.csv")))
      if r and r[0] != "country" and not r[0].startswith(("WEIGHTED","UNWEIGHTED","SD"))]
countries = [l for l in sorted(PC) if has_base(l)]
def baseyrs(l): return list(range(max(2003, RELIABLE.get(l, 2003)), 2020))

def ols_eval(x, y, x_eval):
    x = np.asarray(x, float); y = np.asarray(y, float); n = len(x)
    if n < 3: return None
    xbar = x.mean(); Sxx = ((x - xbar)**2).sum()
    if Sxx <= 0: return None
    b = ((x - xbar)*(y - y.mean())).sum()/Sxx; a = y.mean() - b*xbar
    return dict(b=b, f=a + b*x_eval, n=n)

# ---------------------------------------------------------------- per-band two-step
LNR = {}
def lnr(l, a, y):
    k = (l, a, y)
    if k not in LNR:
        c = DP(l, y, a)
        LNR[k] = None if (c is None or c[1] <= 0 or c[0] <= 0) else np.log(c[0]/c[1])
    return LNR[k]

def mslope(l, a, c):
    pts = [(y, lnr(l, a, y)) for y in range(c-2, c+3)]
    pts = [(y, v) for y, v in pts if v is not None]
    if len(pts) < 3: return None
    x = np.array([p[0] for p in pts], float); yv = np.array([p[1] for p in pts])
    xc = x - x.mean()
    return 100.0*(xc*(yv - yv.mean())).sum()/(xc**2).sum()

def anchor(l, a):
    ys = [y for y in (2019, 2024, 2025) if lnr(l, a, y) is not None]
    if len(ys) < 2: return None, None
    x = np.array(ys, float); yv = np.array([lnr(l, a, y) for y in ys]); xc = x - x.mean()
    return 100.0*(xc*(yv - yv.mean())).sum()/(xc**2).sum(), float(x.mean())

print(f"fitting per-band two-step trends: {len(countries)} populations x {len(FINE)} bands")
BP = {}; nfail = 0
for l in countries:
    rsd = RELIABLE.get(l, 2003); BP[l] = {}
    for a in FINE:
        pre = [(c, mslope(l, a, c)) for c in range(rsd+2, 2018)]
        pre = [(c, v) for c, v in pre if v is not None]
        xp = [c for c, _ in pre]; yp = [v for _, v in pre]
        U = ols_eval(xp, yp, 2019)
        asl, acx = anchor(l, a)
        A = ols_eval(xp + [acx], yp + [asl], 2025) if (asl is not None and xp) else None
        D = sum(DP(l, y, a)[0] for y in baseyrs(l))
        if U is None: nfail += 1
        BP[l][a] = dict(is19=(U["f"] if U else np.nan), sos=(U["b"] if U else np.nan),
                        a_is=(A["f"] if A else np.nan), a_sos=(A["b"] if A else np.nan),
                        npre=len(xp), D=D)
print(f"  bands with no usable two-step fit: {nfail} of {len(countries)*len(FINE)}")

# ---------------------------------------------------------------- shrinkage within band
pop2019 = {l: {a: DP(l, 2019, a)[1] for a in FINE} for l in countries}
GB = {}
for a in FINE:
    v_is = np.array([BP[l][a]["is19"] for l in countries], float)
    v_so = np.array([BP[l][a]["sos"] for l in countries], float)
    v_ai = np.array([BP[l][a]["a_is"] for l in countries], float)
    v_as = np.array([BP[l][a]["a_sos"] for l in countries], float)
    wp = np.array([pop2019[l][a] for l in countries], float)
    def wm(v):
        m = np.isfinite(v) & (wp > 0)
        return float(np.average(v[m], weights=wp[m])) if m.any() else np.nan
    GB[a] = dict(is19=wm(v_is), sos=wm(v_so), a_is=wm(v_ai), a_sos=wm(v_as))
Dmax = {a: max(BP[l][a]["D"] for l in countries) for a in FINE}

# country-level common trend (the published STTa target), for the STTa_asC sensitivity
CTREND = {}
for r in csv.DictReader(open(os.path.join(OUTD, "option1_country_params.csv"))):
    try:
        CTREND[r["country"]] = (float(r["beta_shrunk"]), float(r["gamma_shrunk"]),
                                float(r["beta_anch_shrunk"]), float(r["gamma_anch_shrunk"]))
    except (ValueError, KeyError):
        pass

for l in countries:
    for a in FINE:
        p = BP[l][a]
        p["w"] = float(np.sqrt(p["D"]/Dmax[a])) if Dmax[a] > 0 else 0.0
        g = GB[a]
        def mix(own, glob):
            if not np.isfinite(own): return glob
            return p["w"]*own + (1 - p["w"])*glob
        p["is19_s"] = mix(p["is19"], g["is19"]); p["sos_s"] = mix(p["sos"], g["sos"])
        p["a_is_s"] = mix(p["a_is"], g["a_is"]); p["a_sos_s"] = mix(p["a_sos"], g["a_sos"])
        cb = CTREND.get(l)
        if cb:                                    # shrink toward this country's own all-ages trend
            p["is19_c"] = mix(p["is19"], cb[0]*100.0); p["sos_c"] = mix(p["sos"], cb[1]*200.0)
        else:
            p["is19_c"], p["sos_c"] = p["is19_s"], p["sos_s"]

# ---------------------------------------------------------------- alpha (unchanged from engine)
alpha = {}
for l in countries:
    ys = baseyrs(l); t = np.array(ys, float) - 2019.0; alpha[l] = {}
    for a in FINE:
        r = np.array([DP(l, y, a)[0]/DP(l, y, a)[1] for y in ys])
        alpha[l][a] = float(np.polyfit(t, np.log(np.maximum(r, 1e-9)), 2)[2])

CLIPHIT = [0, 0]
def rate(l, a, model, y):
    t = y - 2019.0; p = BP[l][a]
    if model == "TTa_as":
        b, g = p["is19"]/100.0, p["sos"]/200.0
        if not np.isfinite(b): b, g = GB[a]["is19"]/100.0, GB[a]["sos"]/200.0
    elif model == "STTa_as": b, g = p["is19_s"]/100.0, p["sos_s"]/200.0
    elif model == "STTa_asC": b, g = p["is19_c"]/100.0, p["sos_c"]/200.0
    else:
        ai, asos = p["a_is_s"], p["a_sos_s"]
        b, g = (ai - 6.0*asos)/100.0, asos/200.0
    m = np.exp(b*t + g*t*t)
    CLIPHIT[0] += 1
    if m < CLIP[0] or m > CLIP[1]: CLIPHIT[1] += 1
    return np.exp(alpha[l][a])*min(max(m, CLIP[0]), CLIP[1])

# ---------------------------------------------------------------- periods x groups
recs = []
for l in countries:
    for pname, pyears in PERIODS.items():
        yrs = [y for y in pyears if has_year(l, y)]
        for g, bands in GROUPS.items():
            for m in MODELS:
                ps = []; O = E = 0.0
                for y in yrs:
                    o = sum(DP(l, y, a)[0] for a in bands)
                    e = sum(rate(l, a, m, y)*DP(l, y, a)[1] for a in bands)
                    O += o; E += e; ps.append(100*(o-e)/e if e else np.nan)
                recs.append([l, g, m, pname, len(yrs), round(O), round(E), round(O-E),
                             round(float(np.nanmean(ps)), 3)])
with open(os.path.join(OUTD, "model_agespecific_periods_v1.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["country","group","model","period","n_years","O","E","excess","pscore_meanann"])
    w.writerows(recs)

pool = []
for pname in PERIODS:
    for g in GROUPS:
        for m in MODELS:
            rr = [x for x in recs if x[3] == pname and x[1] == g and x[2] == m]
            O = sum(x[5] for x in rr); E = sum(x[6] for x in rr)
            pool.append([pname, g, m, len(rr), round(O), round(E), round(O-E),
                         round(100*(O-E)/E, 6)])
with open(os.path.join(OUTD, "pooled_agespecific_by_period_v1.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["period","group","model","n_countries","O","E","excess","pooled_Ppct"])
    w.writerows(pool)

with open(os.path.join(OUTD, "agespecific_band_params_v1.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["country","name","band","n_centres","base_deaths","w",
                "islope2019","sos","islope2019_shrunk","sos_shrunk",
                "anch_islope2025","anch_sos","anch_islope2025_shrunk","anch_sos_shrunk"])
    for l in countries:
        for a in FINE:
            p = BP[l][a]
            def q(x, n=4): return "" if not np.isfinite(x) else round(float(x), n)
            w.writerow([l, name.get(l, l), a, p["npre"], round(p["D"]), round(p["w"], 4),
                        q(p["is19"]), q(p["sos"]), q(p["is19_s"]), q(p["sos_s"]),
                        q(p["a_is"]), q(p["a_sos"]), q(p["a_is_s"]), q(p["a_sos_s"])])

# ---------------------------------------------------------------- comparison with published
common = {}
for r in csv.DictReader(open(os.path.join(OUTD, "pooled_option1_by_period.csv"))):
    common[(r["period"], r["group"], r["model"])] = r
PAIR = [("TTa", "TTa_as"), ("STTa", "STTa_as"), ("STTa+", "STTa+_as")]
with open(os.path.join(OUTD, "agespecific_vs_common_v1.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["period","group","model_common","P_common","excess_common",
                "model_agespecific","P_agespecific","excess_agespecific","delta_P_pp"])
    for pname in ["2020-2025","2020-2023","2024-2025"]:
        for g in GROUPS:
            for mc, ma in PAIR:
                c = common.get((pname, g, mc))
                a = next(x for x in pool if x[0] == pname and x[1] == g and x[2] == ma)
                if not c: continue
                w.writerow([pname, g, mc, round(float(c["pooled_Ppct"]), 3), c["excess"],
                            ma, round(a[7], 3), a[6], round(a[7] - float(c["pooled_Ppct"]), 3)])

print(f"\nbaseline multiplier hit the [{CLIP[0]}, {CLIP[1]}] clip in "
      f"{CLIPHIT[1]:,} of {CLIPHIT[0]:,} band-years ({100*CLIPHIT[1]/CLIPHIT[0]:.2f}%)")
print("\nPOOLED ALL-AGES: age-specific trends vs the published common trend")
print(f"  {'period':12}{'model':9}{'common P%':>11}{'age-spec P%':>13}{'delta pp':>10}"
      f"{'common M':>11}{'age-spec M':>12}")
for pname in ["2020-2025","2020-2023","2024-2025"]:
    for mc, ma in PAIR:
        c = common[(pname, "All", mc)]
        a = next(x for x in pool if x[0] == pname and x[1] == "All" and x[2] == ma)
        print(f"  {pname:12}{mc:9}{float(c['pooled_Ppct']):>+11.2f}{a[7]:>+13.2f}"
              f"{a[7]-float(c['pooled_Ppct']):>+10.2f}{int(c['excess'])/1e6:>+11.2f}"
              f"{a[6]/1e6:>+12.2f}")
print("\n  by age group, 2020-2025:")
for g in GROUPS:
    for mc, ma in PAIR:
        c = common[("2020-2025", g, mc)]
        a = next(x for x in pool if x[0] == "2020-2025" and x[1] == g and x[2] == ma)
        print(f"    {g:6}{mc:8}{float(c['pooled_Ppct']):>+9.2f} ->{a[7]:>+9.2f}"
              f"   ({a[7]-float(c['pooled_Ppct']):>+.2f} pp)")
sc = next(x for x in pool if x[0] == "2020-2025" and x[1] == "All" and x[2] == "STTa_asC")
print(f"\n  sensitivity STTa_asC (shrunk toward each country's own all-ages trend): "
      f"{sc[7]:+.2f}%  ({sc[6]/1e6:+.2f}M)")
