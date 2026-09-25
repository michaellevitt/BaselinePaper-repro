#!/usr/bin/env python3
"""
backtest_leecarter_2014_2019_v1.py — validate the Lee-Carter baseline on a quiet period.

PURPOSE
  Two things at once.  (1) Confirm the Poisson Lee-Carter implementation in
  model_leecarter_periods_v1.py behaves like Lee-Carter rather than like a coding error.
  (2) Measure how much a constant-drift extrapolation mis-states expected deaths when the
  rate of mortality improvement is itself decelerating -- the paper's central empirical claim.

DESIGN
  Fit each population on [first eligible year .. 2013] and forecast 2014-2019.  There is no
  pandemic in that window, so a well-calibrated baseline should return a pooled P-score near
  zero.  Fa (flat mean of the last three training years, 2011-2013) is carried as the reference
  because it makes no trend assumption at all.

  A positive P-score here means the baseline expected FEWER deaths than actually occurred,
  i.e. it over-predicted mortality improvement.

NOTE ON PROVENANCE
  John deprecated back-testing in July 2026 as a way of choosing between models.  This is not
  used for that.  It is used only to check an implementation and to quantify one bias, in
  response to the PNAS Board's suggestion that Lee-Carter would do better than the paper's
  quadratic trends.
"""
import os, sys, csv, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
OUTD = os.path.join(ROOT, "output")
sys.path.insert(0, HERE)
from excess_anchor_window import load, cell_for
from model_leecarter_periods_v1 import fit_lee_carter          # reuse the exact same fitter

FINE = ["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
        "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
RELIABLE = {"BGR":2015,"CHL":2011,"HRV":2007,"EST":2005,"HUN":2006,"LVA":2007,"LTU":2006,
            "POL":2008,"SVK":2006}
TRAIN_END, TEST = 2013, list(range(2014, 2020))

data, name = load()
def DP(l, y, a): return cell_for(data, l, y, a, "T")
def ok(l, ys): return all(DP(l, y, a) is not None for y in ys for a in FINE)

PC = [r[0] for r in csv.reader(open(os.path.join(OUTD, "slope_of_slopes_CI.csv")))
      if r and r[0] != "country" and not r[0].startswith(("WEIGHTED","UNWEIGHTED","SD"))]

rows = []
for l in sorted(PC):
    lo = max(2003, RELIABLE.get(l, 2003))
    train = list(range(lo, TRAIN_END + 1))
    if len(train) < 8 or not ok(l, train) or not ok(l, TEST):
        continue                                     # needs a usable training run and a full hold-out
    Dm = np.array([[DP(l, y, a)[0] for y in train] for a in FINE], float)
    Pm = np.array([[DP(l, y, a)[1] for y in train] for a in FINE], float)
    res = fit_lee_carter(Dm, Pm)
    if res is None: continue
    alpha, beta, kappa, dev, nit = res
    drift = (kappa[-1] - kappa[0]) / (len(kappa) - 1)

    O = E_lc = E_fa = 0.0
    for y in TEST:
        k = kappa[-1] + drift * (y - TRAIN_END)
        for i, a in enumerate(FINE):
            d, p = DP(l, y, a)
            O += d
            E_lc += float(np.exp(alpha[i] + beta[i] * k)) * p
            fa = (sum(DP(l, yy, a)[0] for yy in (2011, 2012, 2013))
                  / sum(DP(l, yy, a)[1] for yy in (2011, 2012, 2013)))
            E_fa += fa * p
    rows.append([l, name.get(l, l), train[0], TRAIN_END, len(train), round(O),
                 round(E_lc), round(100*(O-E_lc)/E_lc, 2),
                 round(E_fa), round(100*(O-E_fa)/E_fa, 2), f"{drift:.5f}"])

hdr = ["country","name","train_from","train_to","n_train","O_2014_2019",
       "E_LeeCarter","P_LeeCarter_pct","E_Fa","P_Fa_pct","kappa_drift"]
with open(os.path.join(OUTD, "backtest_leecarter_2014_2019_v1.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(hdr); w.writerows(rows)

O = sum(r[5] for r in rows); Elc = sum(r[6] for r in rows); Efa = sum(r[8] for r in rows)
plc = [r[7] for r in rows]; pfa = [r[9] for r in rows]
print(f"BACK-TEST on 2014-2019 (no pandemic): {len(rows)} populations, "
      f"trained to {TRAIN_END}\n")
print(f"  {'':22}{'pooled P%':>11}{'median P%':>11}{'>0':>6}{'<0':>6}{'worst +':>10}")
print(f"  {'Lee-Carter (drift)':22}{100*(O-Elc)/Elc:>+11.2f}{np.median(plc):>+11.2f}"
      f"{sum(1 for x in plc if x > 0):>6}{sum(1 for x in plc if x < 0):>6}{max(plc):>+10.2f}")
print(f"  {'Fa (flat 2011-2013)':22}{100*(O-Efa)/Efa:>+11.2f}{np.median(pfa):>+11.2f}"
      f"{sum(1 for x in pfa if x > 0):>6}{sum(1 for x in pfa if x < 0):>6}{max(pfa):>+10.2f}")
print(f"\n  Lee-Carter expected {Elc:,.0f} deaths; {O:,.0f} occurred -> it under-predicted by "
       f"{O-Elc:+,.0f} ({100*(O-Elc)/Elc:+.2f}%) in a period with no pandemic.")
print(f"  Fa expected {Efa:,.0f} -> {O-Efa:+,.0f} ({100*(O-Efa)/Efa:+.2f}%).")
print("\n  Largest Lee-Carter over-predictions of improvement (top 8 by P%):")
for r in sorted(rows, key=lambda x: -x[7])[:8]:
    print(f"    {r[1][:26]:28}{r[7]:>+8.2f}%   (Fa {r[9]:>+6.2f}%,  drift {r[10]})")
