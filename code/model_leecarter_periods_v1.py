#!/usr/bin/env python3
"""
model_leecarter_periods_v1.py — a Lee-Carter baseline for the Paper A model family.

WHY
  Both PNAS Editorial Board members asked for an established demographic method, naming
  Lee-Carter.  Dowd et al. (medRxiv 2026, doi 10.64898/2026.02.25.26347112) use exactly this
  and reach the OPPOSITE conclusion for 2024 (31/34 countries still below their pre-pandemic
  life-expectancy trajectory, i.e. an excess), where Paper A finds death deficits in 33/38.
  Adding Lee-Carter puts that disagreement inside one framework, on one dataset, so the cause
  can be identified rather than argued.  Also requested for the group's EJE paper on Austria.

SPECIFICATION — deliberately matched to Dowd et al.
  Poisson Lee-Carter (Brouhns, Denuit & Vermunt 2002), the same model StMoMo fits:
      D_ay ~ Poisson(P_ay * exp(alpha_a + beta_a * kappa_y)),  sum_y kappa_y = 0,  sum_a beta_a = 1
  fitted by the standard uni-dimensional Newton updates, then forecast by a random walk with
  drift on kappa, drift = (kappa_last - kappa_first)/(n-1).  Expected deaths for a pandemic or
  post-pandemic year are P_ay * exp(alpha_a + beta_a * kappa_hat_y).

  DIFFERENCES FROM DOWD, both unavoidable here and both stated in the manuscript:
    * They ungroup to single years of age to 100+ with the Penalized Composite Link Model;
      this fits the paper's own twenty five-year bands to 90+.  Lee-Carter on five-year bands
      is standard, but the age grid is coarser.
    * They fit 2000-2019 for all countries.  The primary fit here uses the paper's own eligible
      window (2003, or the country's first reliable year, to 2019) so that Lee-Carter is
      comparable with Fa/TTa/STTa/STTa+ on identical input.  Pass --start 2010 for the short
      window they use as their sensitivity check.

  No clipping is applied.  The other trend models clip exp(beta t + gamma t^2) to [0.4, 2.5];
  a six-year Lee-Carter forecast never approaches those bounds, and imposing them would be an
  extra assumption the method does not carry.

OUTPUT (same long format as model_option1_periods.py, so tables can merge on country/group/period)
  output/model_leecarter_periods_v1.csv    country,group,model,period,n_years,O,E,excess,pscore_meanann
  output/pooled_leecarter_by_period_v1.csv period,group,model,n_countries,O,E,excess,pooled_Ppct,...
  output/leecarter_country_params_v1.csv   per country: drift, kappa range, fit window, deviance
"""
import os, sys, csv, argparse, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
OUTD = os.path.join(ROOT, "output")
sys.path.insert(0, HERE)
from excess_anchor_window import load, cell_for

FINE = ["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
        "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
def astart(a): return 0 if a == "0" else (1 if a == "1-4" else int(a.split("-")[0].replace("+","")))
GROUPS = {"All": FINE, "GE65": [a for a in FINE if astart(a) >= 65],
          "LT65": [a for a in FINE if astart(a) < 65]}
W = 5
RELIABLE = {"BGR":2015,"CHL":2011,"HRV":2007,"EST":2005,"HUN":2006,"LVA":2007,"LTU":2006,
            "POL":2008,"SVK":2006}
PERIODS = {"2020-2025":[2020,2021,2022,2023,2024,2025], "2020-2023":[2020,2021,2022,2023],
           "2024-2025":[2024,2025]}
PERIODS.update({str(y): [y] for y in range(2020, 2026)})

ap = argparse.ArgumentParser()
ap.add_argument("--start", type=int, default=0,
                help="fixed first fitting year (e.g. 2010); default 0 = the paper's own window")
ap.add_argument("--tag", default="v1", help="output filename tag")
args = ap.parse_args()

data, name = load()
def DP(l, y, a): return cell_for(data, l, y, a, "T")
def has_year(l, y): return all(DP(l, y, a) is not None for a in FINE)
def has_base(l):
    se = max(2003 + (W-1), RELIABLE.get(l, 2003) + (W-1))
    return all(DP(l, y, a) is not None for y in range(se-(W-1), 2020) for a in FINE)

PC = [r[0] for r in csv.reader(open(os.path.join(OUTD, "slope_of_slopes_CI.csv")))
      if r and r[0] != "country" and not r[0].startswith(("WEIGHTED","UNWEIGHTED","SD"))]
countries = [l for l in sorted(PC) if has_base(l)]

def baseyrs(l):
    lo = max(2003, RELIABLE.get(l, 2003))
    if args.start:
        lo = max(lo, args.start)
        while lo > 2003 and not all(DP(l, y, a) is not None for y in range(lo, 2020) for a in FINE):
            lo -= 1
    return list(range(lo, 2020))

# ---------------------------------------------------------------- Poisson Lee-Carter
def fit_lee_carter(Dm, Pm, iters=3000, tol=1e-10):
    """Dm, Pm: (n_age, n_year) deaths and exposures.  Returns alpha, beta, kappa, deviance, n_iter."""
    A, Y = Dm.shape
    with np.errstate(divide="ignore", invalid="ignore"):
        alpha = np.log(np.maximum(Dm.sum(1), 0.5) / np.maximum(Pm.sum(1), 1e-9))
    beta = np.full(A, 1.0 / A)
    kappa = np.linspace(1.0, -1.0, Y)           # non-degenerate start; sign fixed by the constraint below
    prev = np.inf
    for it in range(iters):
        # kappa | alpha, beta
        mu = Pm * np.exp(alpha[:, None] + np.outer(beta, kappa))
        num = ((Dm - mu) * beta[:, None]).sum(0); den = (mu * (beta**2)[:, None]).sum(0)
        kappa = kappa + np.where(den > 0, num / np.maximum(den, 1e-12), 0.0)
        kappa -= kappa.mean()                                            # sum_y kappa = 0
        # alpha | beta, kappa
        mu = Pm * np.exp(alpha[:, None] + np.outer(beta, kappa))
        alpha = alpha + (Dm - mu).sum(1) / np.maximum(mu.sum(1), 1e-12)
        # beta | alpha, kappa
        mu = Pm * np.exp(alpha[:, None] + np.outer(beta, kappa))
        num = ((Dm - mu) * kappa[None, :]).sum(1); den = (mu * (kappa**2)[None, :]).sum(1)
        beta = beta + np.where(den > 0, num / np.maximum(den, 1e-12), 0.0)
        s = beta.sum()
        if abs(s) < 1e-12: return None
        beta, kappa = beta / s, kappa * s                                # sum_a beta = 1
        kappa -= kappa.mean(); alpha = alpha + beta * 0.0                # keep kappa centred
        mu = Pm * np.exp(alpha[:, None] + np.outer(beta, kappa))
        with np.errstate(divide="ignore", invalid="ignore"):
            dev = 2 * np.nansum(np.where(Dm > 0, Dm * np.log(Dm / np.maximum(mu, 1e-12)), 0.0) - (Dm - mu))
        if abs(prev - dev) < tol * max(1.0, abs(dev)): break
        prev = dev
    return alpha, beta, kappa, float(dev), it + 1

WIN = f"fixed from {args.start}" if args.start else "the paper's own"
print(f"Poisson Lee-Carter, {len(countries)} populations, window = {WIN}")
LC = {}
for l in countries:
    ys = baseyrs(l)
    Dm = np.array([[DP(l, y, a)[0] for y in ys] for a in FINE], float)
    Pm = np.array([[DP(l, y, a)[1] for y in ys] for a in FINE], float)
    res = fit_lee_carter(Dm, Pm)
    if res is None:
        print(f"  {l}: Lee-Carter failed to identify beta — skipped"); continue
    alpha, beta, kappa, dev, nit = res
    drift = (kappa[-1] - kappa[0]) / (len(kappa) - 1)
    LC[l] = dict(alpha=alpha, beta=beta, kappa=kappa, drift=drift, years=ys,
                 dev=dev, nit=nit, k2019=kappa[-1])

def lc_rate(l, a, y):
    """Forecast (or fitted) central death rate for band a in year y."""
    p = LC[l]; i = FINE.index(a)
    k = p["k2019"] + p["drift"] * (y - p["years"][-1]) if y > p["years"][-1] \
        else p["kappa"][p["years"].index(y)]
    return float(np.exp(p["alpha"][i] + p["beta"][i] * k))

# ---------------------------------------------------------------- periods x groups
recs = []
for l in LC:
    for pname, pyears in PERIODS.items():
        yrs = [y for y in pyears if has_year(l, y)]
        for g, bands in GROUPS.items():
            ps = []; O = E = 0.0
            for y in yrs:
                o = sum(DP(l, y, a)[0] for a in bands)
                e = sum(lc_rate(l, a, y) * DP(l, y, a)[1] for a in bands)
                O += o; E += e; ps.append(100*(o-e)/e if e else np.nan)
            recs.append([l, g, "LCa", pname, len(yrs), round(O), round(E), round(O-E),
                         round(float(np.nanmean(ps)), 3)])
with open(os.path.join(OUTD, f"model_leecarter_periods_{args.tag}.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["country","group","model","period","n_years","O","E","excess","pscore_meanann"])
    w.writerows(recs)

pool = []
for pname in PERIODS:
    for g in GROUPS:
        rr = [x for x in recs if x[3] == pname and x[1] == g]
        O = sum(x[5] for x in rr); E = sum(x[6] for x in rr)
        pool.append([pname, g, "LCa", len(rr), round(O), round(E), round(O-E),
                     round(100*(O-E)/E, 3), round(float(np.mean([x[8] for x in rr])), 3)])
with open(os.path.join(OUTD, f"pooled_leecarter_by_period_{args.tag}.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["period","group","model","n_countries","O","E","excess","pooled_Ppct","mean_meanann"])
    w.writerows(pool)

with open(os.path.join(OUTD, f"leecarter_country_params_{args.tag}.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["country","name","fit_from","fit_to","n_years","kappa_drift","kappa_first",
                "kappa_2019","deviance","iterations"])
    for l in sorted(LC):
        p = LC[l]
        w.writerow([l, name.get(l, l), p["years"][0], p["years"][-1], len(p["years"]),
                    f"{p['drift']:.5f}", f"{p['kappa'][0]:.4f}", f"{p['k2019']:.4f}",
                    f"{p['dev']:.1f}", p["nit"]])

print(f"\nwrote model_leecarter_periods_{args.tag}.csv ({len(recs)} rows) + pooled + params")
print("\nPOOLED all-ages, Lee-Carter:")
print(f"  {'period':12}{'O':>12}{'E':>12}{'excess':>12}{'pooled P%':>11}")
for pname in ["2020-2025","2020-2023","2024-2025"] + [str(y) for y in range(2020, 2026)]:
    r = next(x for x in pool if x[0] == pname and x[1] == "All")
    print(f"  {pname:12}{r[4]:>12,}{r[5]:>12,}{r[6]:>+12,}{r[7]:>+11.2f}")
n_def = sum(1 for l in LC if next(x for x in recs if x[0]==l and x[1]=="All" and x[3]=="2024-2025")[7] < 0)
print(f"\n  countries with a Lee-Carter death DEFICIT in 2024-2025: {n_def}/{len(LC)}")
