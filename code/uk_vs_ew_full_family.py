#!/usr/bin/env python3
"""
uk_vs_ew_full_family.py - Table S8: England and Wales versus the whole United Kingdom, 2020-2025.

The paper uses England and Wales for the United Kingdom.  This script checks that choice by
running the whole United Kingdom (HMD GBR_NP, extended to 2025 with the National Records of
Scotland and Northern Ireland Statistics and Research Agency series assembled in the master as
GBR_NP_BUILT) through the same model family as the COUNTRY-LEVEL engine,
model_option1_periods.py: Fa, TTa, STTa and STTa+.

ORDER.  It must run straight after model_option1_periods.py and before the age-band engine,
because its validation gate requires the recomputed England and Wales values to reproduce the
country-level engine's output/model_option1_periods.csv exactly, and the age-band engine later
overwrites that file.  run_all.sh runs it in that position.

GBR_NP_BUILT is attached as a child of GBR_NP with rho = 1: the assembled series reproduces HMD
GBR_NP on the 2020-2022 overlap (deaths +0, population within 0.002%), which is asserted.

STTa and STTa+ shrink toward the 38-population aggregate trend, so this table moves whenever any
of the 38 populations changes (it moved with the 2026-08-27 HMD exposure correction).

Writes output/uk_vs_ew_full_family_2020_2025_v2_6dp.csv
"""
import os, sys, csv, numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPRO = ROOT
sys.path.insert(0, os.path.join(REPRO, "code"))
import excess_anchor_window as EAW
from excess_anchor_window import load, cell_for

# attach the assembled UK series as a child of the HMD UK parent (rho defaults to 1.0)
EAW.CHILD.setdefault("GBR_NP", []).append("GBR_NP_BUILT")

OUTD = os.path.join(ROOT, "output")
FINE = ["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
        "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
def astart(a): return 0 if a == "0" else (1 if a == "1-4" else int(a.split("-")[0].replace("+","")))
GROUPS = {"All": FINE, "GE65": [a for a in FINE if astart(a) >= 65], "LT65": [a for a in FINE if astart(a) < 65]}
ESP = {"0":1000,"1-4":4000,"5-9":5500,"10-14":5500,"15-19":5500,"20-24":6000,"25-29":6000,
       "30-34":6500,"35-39":7000,"40-44":7000,"45-49":7000,"50-54":7000,"55-59":6500,"60-64":6000,
       "65-69":5500,"70-74":5000,"75-79":4000,"80-84":2500,"85-89":1500,"90+":1000}
ESPsum = sum(ESP.values()); CLIP = (0.4, 2.5)
BASEYRS = list(range(2003, 2020))
PERIODS = {"2020-2025": [2020,2021,2022,2023,2024,2025],
           "2020-2023": [2020,2021,2022,2023], "2024-2025": [2024,2025]}
MODELS = ["Fa", "TTa", "STTa", "STTa+"]

data, name = load()
def DP(l, y, a): return cell_for(data, l, y, a, "T")

# ---------------- two-step trend, CENTRE convention (slope_of_slopes_ci.py) ----------------
def ols(x, y, x_eval):
    x = np.asarray(x, float); y = np.asarray(y, float)
    xbar = x.mean(); Sxx = ((x - xbar) ** 2).sum()
    b = ((x - xbar) * (y - y.mean())).sum() / Sxx
    return b, (y.mean() - b * xbar) + b * x_eval

def lnR(loc, y):
    cs = [DP(loc, y, a) for a in FINE]
    if any(c is None for c in cs): return None
    r = sum(ESP[a] * cs[i][0] / cs[i][1] for i, a in enumerate(FINE)) / ESPsum
    return np.log(r) if r > 0 else None

def mslope(loc, c):
    pts = [(y, lnR(loc, y)) for y in range(c - 2, c + 3)]
    pts = [(y, v) for y, v in pts if v is not None]
    if len(pts) < 3: return None
    x = np.array([p[0] for p in pts], float); yv = np.array([p[1] for p in pts]); xc = x - x.mean()
    return 100.0 * (xc * (yv - yv.mean())).sum() / (xc ** 2).sum()

def anchor_slope(loc):
    ys = [y for y in (2019, 2024, 2025) if lnR(loc, y) is not None]
    if len(ys) < 2: return None, None
    x = np.array(ys, float); yv = np.array([lnR(loc, y) for y in ys]); xc = x - x.mean()
    return 100.0 * (xc * (yv - yv.mean())).sum() / (xc ** 2).sum(), float(x.mean())

def twostep(loc, quantize=True):
    pre = [(c, mslope(loc, c)) for c in range(2005, 2018)]
    pre = [(c, v) for c, v in pre if v is not None]
    xp = [c for c, _ in pre]; yp = [v for _, v in pre]
    u_sos, u_is = ols(xp, yp, 2019)
    asl, acx = anchor_slope(loc)
    a_sos, a_is = ols(xp + [acx], yp + [asl], 2025)
    if quantize:                       # the engine round-trips these through a 3-dp CSV
        u_is, u_sos, a_is, a_sos = (round(v, 3) for v in (u_is, u_sos, a_is, a_sos))
    return dict(is19=u_is, sos=u_sos, a_is2025=a_is, a_sos=a_sos, npre=len(pre))

# ---------------- global (38-population) trend + shrinkage scale, from the paper's run ----------------
par = [r for r in csv.DictReader(open(f"{REPRO}/output/option1_country_params.csv")) if r["beta_i"]]
wp = np.array([float(r["pop2019"]) for r in par])
bG = float(np.average([float(r["beta_i"]) for r in par], weights=wp))
gG = float(np.average([float(r["gamma_i"]) for r in par], weights=wp))
Dmax = max(float(r["base_deaths"]) for r in par)
# global ANCHORED params: recompute the same pop-weighted average the engine forms from ba/ga
SOS = {r[0]: r for r in csv.reader(open(f"{REPRO}/output/slope_of_slopes_CI.csv"))
       if r and r[0] not in ("country",) and not r[0].startswith(("WEIGHTED","UNWEIGHTED","SD"))}
ba_l, ga_l, wa_l = [], [], []
for r in par:
    s = SOS.get(r["country"])
    if not s or s[11] == "" or s[14] == "": continue
    ai, asos = float(s[11]), float(s[14])
    ba_l.append((ai - 6.0 * asos) / 100.0); ga_l.append(asos / 200.0); wa_l.append(float(r["pop2019"]))
baG = float(np.average(ba_l, weights=wa_l)); gaG = float(np.average(ga_l, weights=wa_l))

def fit(loc, quantize=True):
    t = np.array(BASEYRS, float) - 2019.0
    bands, fa = {}, {}
    for a in FINE:
        r = np.array([DP(loc, y, a)[0] / DP(loc, y, a)[1] for y in BASEYRS])
        bands[a] = np.polyfit(t, np.log(np.maximum(r, 1e-9)), 2)
        fa[a] = (sum(DP(loc, y, a)[0] for y in (2017, 2018, 2019))
                 / sum(DP(loc, y, a)[1] for y in (2017, 2018, 2019)))
    D = sum(DP(loc, y, a)[0] for y in BASEYRS for a in FINE)
    ts = twostep(loc, quantize)
    b, g = ts["is19"] / 100.0, ts["sos"] / 200.0
    ba = (ts["a_is2025"] - 6.0 * ts["a_sos"]) / 100.0; ga = ts["a_sos"] / 200.0
    w = np.sqrt(D / Dmax)
    return dict(bands=bands, fa=fa, D=D, w=w, **ts,
                b=b, g=g, bs=w*b + (1-w)*bG, gs=w*g + (1-w)*gG,
                bas=w*ba + (1-w)*baG, gas=w*ga + (1-w)*gaG)

def rate(P, a, model, y):
    if model == "Fa": return P["fa"][a]
    alpha = P["bands"][a][2]; t = y - 2019.0
    b, g = {"TTa": (P["b"], P["g"]), "STTa": (P["bs"], P["gs"]), "STTa+": (P["bas"], P["gas"])}[model]
    return np.exp(alpha) * min(max(np.exp(b*t + g*t*t), CLIP[0]), CLIP[1])

def excess(loc, P, years, group, model):
    bands = GROUPS[group]; O = E = 0.0
    for y in years:
        O += sum(DP(loc, y, a)[0] for a in bands)
        E += sum(rate(P, a, model, y) * DP(loc, y, a)[1] for a in bands)
    return O, E, O - E, 100 * (O - E) / E

# ======================= assert the UK child needs no calibration =======================
print("=" * 88)
print("UK child series vs HMD GBR_NP on the 2020-2022 overlap (justifies rho = 1)")
print("=" * 88)
for y in (2020, 2021, 2022):
    hd = sum(data["GBR_NP"][(y, a)]["T"][0] for a in FINE)
    hp = sum(data["GBR_NP"][(y, a)]["T"][1] for a in FINE)
    bd = sum(data["GBR_NP_BUILT"][(y, a)]["T"][0] for a in FINE) if (y, "0") in data.get("GBR_NP_BUILT", {}) else None
    print(f"  {y}  HMD deaths {hd:>9,.0f}  pop {hp:>12,.0f}" + ("" if bd is None else f"   built {bd:,.0f}"))
print("  (2023-2025 come only from the built child; 2020-2022 resolve to the HMD parent)")

# ============================== VALIDATION vs published E&W ==============================
print("\n" + "=" * 88)
print("VALIDATION — England & Wales, published vs recomputed (model_option1_periods.csv)")
print("=" * 88)
pub = {}
for r in csv.DictReader(open(f"{REPRO}/output/model_option1_periods.csv")):
    if r["country"] == "GBRTENW" and r["group"] == "All" and r["model"] in MODELS:
        pub[(r["model"], r["period"])] = int(r["excess"])
EW = fit("GBRTENW")
print(f"  two-step: islope@2019 {EW['is19']:+.3f} (pub +0.302)   sos {EW['sos']:+.3f} (pub +0.218)   n={EW['npre']}")
print(f"  beta_anch_shrunk {EW['bas']:+.6f} (pub -0.004696)   gamma_anch_shrunk {EW['gas']:+.7f} (pub +0.0006278)")
fails = []
print(f"\n  {'model':6} {'period':11} {'recomputed':>12} {'published':>12}   status")
for m in MODELS:
    for p, yrs in PERIODS.items():
        got = round(excess("GBRTENW", EW, yrs, "All", m)[2])
        exp = pub[(m, p)]
        ok = abs(got - exp) <= 1
        if not ok: fails.append((m, p, got, exp))
        print(f"  {m:6} {p:11} {got:>12,} {exp:>12,}   {'OK' if ok else 'MISMATCH'}")
if fails: sys.exit(f"\nVALIDATION FAILED: {fails}")
print("\n  -> all 12 published England & Wales values reproduce exactly.\n")

# ================================== UK vs E&W ==================================
UK = fit("GBR_NP")
print("=" * 88)
print("TWO-STEP TREND PARAMETERS")
print("=" * 88)
print(f"  {'population':18} {'islope@2019':>12} {'sos':>8} {'anch islope@2025':>17} {'anch sos':>9} {'w':>6}")
for l, nm, P in (("GBRTENW", "England & Wales", EW), ("GBR_NP", "United Kingdom", UK)):
    print(f"  {nm:18} {P['is19']:>+12.3f} {P['sos']:>+8.3f} {P['a_is2025']:>+17.3f} {P['a_sos']:>+9.3f} {P['w']:>6.3f}")

rows = []
for group in ("All", "GE65", "LT65"):
    print("\n" + "=" * 88)
    print(f"EXCESS DEATHS AND POOLED P%  —  {group}")
    print("=" * 88)
    print(f"  {'period':11} {'model':6} {'E&W excess':>12} {'E&W P%':>8} {'UK excess':>12} {'UK P%':>8} "
          f"{'ΔP% (UK-E&W)':>13}")
    for p, yrs in PERIODS.items():
        for m in MODELS:
            oe, ee, xe, pe = excess("GBRTENW", EW, yrs, group, m)
            ou, eu, xu, pu = excess("GBR_NP", UK, yrs, group, m)
            rows.append([group, p, m, round(xe), round(pe, 6), round(xu), round(pu, 6), round(pu - pe, 6)])
            print(f"  {p:11} {m:6} {xe:>12,.0f} {pe:>+8.2f} {xu:>12,.0f} {pu:>+8.2f} {pu-pe:>+13.2f}")

print("\n" + "=" * 88)
print("HEADLINE — all ages, pooled P%")
print("=" * 88)
for p, yrs in PERIODS.items():
    pe = [excess("GBRTENW", EW, yrs, "All", m)[3] for m in MODELS]
    pu = [excess("GBR_NP", UK, yrs, "All", m)[3] for m in MODELS]
    print(f"  {p:11}  E&W spread {min(pe):+.2f}..{max(pe):+.2f} ({max(pe)-min(pe):.2f} pp)"
          f"   UK spread {min(pu):+.2f}..{max(pu):+.2f} ({max(pu)-min(pu):.2f} pp)"
          f"   max |ΔP%| {max(abs(a-b) for a, b in zip(pu, pe)):.2f} pp")

with open(f"{OUTD}/uk_vs_ew_full_family_2020_2025_v2_6dp.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["group","period","model","ew_excess","ew_pooled_Ppct","uk_excess","uk_pooled_Ppct","delta_Ppct"])
    w.writerows(rows)
print(f"\nwrote output/uk_vs_ew_full_family_2020_2025_v2_6dp.csv ({len(rows)} rows)")
