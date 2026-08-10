#!/usr/bin/env python3
"""
Cumulative 2020-2023 excess deaths under flat (M1) and trended (M4) baselines,
each plain and 2024-anchored (M1*, M4*), for two reference windows (2015-19, 2017-19).

M1 = pooled flat per-capita rate  (ΣD_base / ΣP_base) per band
M4 = linear trend of per-capita rate on year per band  (mirrors all_methods_outputs.py)
Both fit PER 5-YEAR AGE BAND then summed; country series = HMD parent merged with
its Eurostat/ONS build child (prefer HMD, fill missing years).

Input : data/master_5x1_DPM_90plus.csv
Output: output/excess_anchor_window_2020_2023.csv  (columns incl. `model` in M1/M1*/M4/M4*)
"""
import csv, os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MASTER = os.path.join(ROOT, "data", "master_5x1_DPM_90plus.csv")
OUT    = os.path.join(ROOT, "output", "excess_anchor_window_2020_2023.csv")

SEXES = ["F", "M", "T"]
AGES = ["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
        "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
PRED = [2020, 2021, 2022, 2023]
WINDOWS = {"2015-19": [2015,2016,2017,2018,2019], "2017-19": [2017,2018,2019]}

PARENT = {"GBRTENW_ONS":"GBRTENW","DEU_EUROSTAT":"DEUTNP","ITA_EUROSTAT":"ITA",
          "ESP_EUROSTAT":"ESP","FRA_EUROSTAT":"FRATNP","POL_EUROSTAT":"POL",
          "NLD_EUROSTAT":"NLD","AUT_EUROSTAT":"AUT","IRL_EUROSTAT":"IRL",
          "BGR_EUROSTAT":"BGR","CZE_EUROSTAT":"CZE","HUN_EUROSTAT":"HUN",
          "GRC_EUROSTAT":"GRC","SVN_EUROSTAT":"SVN",
          # child series carrying 2025 (second anchor) — synced with compute_TT_lines.py so the
          # headline pipeline ingests all available 2025 (was silently dropping ~18 populations)
          "BEL_EUROSTAT":"BEL","CHE_EUROSTAT":"CHE","EST_EUROSTAT":"EST","FIN_EUROSTAT":"FIN",
          "HRV_EUROSTAT":"HRV","JPN_ESTAT":"JPN","HKG_CSD":"HKG","LTU_EUROSTAT":"LTU",
          "LUX_EUROSTAT":"LUX","LVA_EUROSTAT":"LVA","PRT_EUROSTAT":"PRT","SVK_EUROSTAT":"SVK",
          "SWE_EUROSTAT":"SWE","TWN_MOI":"TWN","ISL_EUROSTAT":"ISL",
          "USA_WONDER":"USA","CHL_DEIS":"CHL","KOR_KOSTAT":"KOR","NOR_EUROSTAT":"NOR"}
CHILD = defaultdict(list)
for b, p in PARENT.items():
    CHILD[p].append(b)
BUILDS = set(PARENT)

# HMD-basis calibration of national/child segments (build_hmd_calibration.py): CALIB[loc][band]=rho.
# Applied in cell_for: a child cell's population is divided by rho so its rate = rho * (D/P) is on the
# HMD (exposure) basis, consistent with the pre-2019 HMD baseline. Absent file -> no calibration.
CALIB = defaultdict(dict)
try:
    _cf = os.path.join(ROOT, "output", "hmd_calibration.csv")
    with open(_cf) as _f:
        for _r in csv.DictReader(_f):
            try: CALIB[_r["location"]][_r["age"]] = float(_r["rho"])
            except (KeyError, ValueError): pass
except FileNotFoundError:
    pass

def load():
    data = defaultdict(dict); name = {}
    with open(MASTER) as f:
        for r in csv.DictReader(f):
            loc = r["location"]
            try: y = int(r["year"])
            except (TypeError, ValueError): continue
            a = r["age"]; name.setdefault(loc, r["location_name"])
            cell = {}
            for s in SEXES:                       # keep each sex that is present; Total-only rows are OK
                D, P = r[f"D_{s}"], r[f"P_{s}"]
                if D == "" or P == "" or float(P) == 0: continue
                cell[s] = (float(D), float(P))
            if "T" in cell: data[loc][(y, a)] = cell   # analysis uses sex=Total throughout
    return data, name

def cell_for(data, loc, y, a, sex):
    c = data.get(loc, {}).get((y, a))
    if c: return c[sex]                            # HMD parent — already the calibration basis
    for b in CHILD.get(loc, []):
        c = data.get(b, {}).get((y, a))
        if c:                                      # national child — put on HMD basis: rate *= rho <=> P /= rho
            rho = CALIB.get(loc, {}).get(a, 1.0)
            if rho != 1.0:
                D, P = c[sex]
                return (D, P / rho)
            return c[sex]
    return None

def merged_series(data, loc, sex):
    D, P = {}, {}
    for y in range(2015, 2025):
        for a in AGES:
            cell = cell_for(data, loc, y, a, sex)
            if cell: D[(y, a)], P[(y, a)] = cell
    return D, P

def have(D, years):
    return all((y, a) in D for y in years for a in AGES)

def cum_flat(D, P, base):                       # M1: pooled rate
    rate = {a: sum(D[(y,a)] for y in base) / sum(P[(y,a)] for y in base) for a in AGES}
    O = E = 0.0
    for py in PRED:
        for a in AGES:
            E += rate[a] * P[(py,a)]; O += D[(py,a)]
    return O, E

def cum_trend(D, P, base):                      # M4: OLS of rate on year
    yc = sum(base)/len(base); sxx = sum((y-yc)**2 for y in base)
    a_int, b_slope = {}, {}
    for a in AGES:
        mr = sum(D[(y,a)]/P[(y,a)] for y in base)/len(base)
        sxy = sum((y-yc)*(D[(y,a)]/P[(y,a)]) for y in base)
        b = sxy/sxx; b_slope[a] = b; a_int[a] = mr - b*yc
    O = E = 0.0
    for py in PRED:
        for a in AGES:
            E += (a_int[a] + b_slope[a]*py) * P[(py,a)]; O += D[(py,a)]
    return O, E

FAMILY = {"M1": cum_flat, "M4": cum_trend}

def main():
    data, name = load()
    countries = sorted(l for l in data if l not in BUILDS)
    rows = []
    for loc in countries:
        for sex in SEXES:
            D, P = merged_series(data, loc, sex)
            if not (have(D, PRED) and all((y,a) in P for y in PRED for a in AGES)):
                continue
            for win, base in WINDOWS.items():
                if not have(D, base): continue
                has24 = have(D, [2024])
                for fam, fn in FAMILY.items():
                    O, E = fn(D, P, base)
                    rows.append([loc, name[loc], sex, win, fam, O, E])
                    if has24:
                        O2, E2 = fn(D, P, base + [2024])
                        rows.append([loc, name[loc], sex, win, fam+"*", O2, E2])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["country","country_name","sex","window","model",
                    "cum_observed","cum_expected","cum_excess","pscore_pct"])
        for loc, nm, sex, win, model, O, E in rows:
            X = O - E; ps = 100*X/E if E else float("nan")
            w.writerow([loc, nm, sex, win, model, round(O), round(E), round(X), round(ps,3)])
    print(f"Wrote {OUT}  ({len(rows)} rows)")

if __name__ == "__main__":
    main()
