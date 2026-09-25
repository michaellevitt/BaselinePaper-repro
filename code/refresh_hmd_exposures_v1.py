#!/usr/bin/env python3
"""Refresh the master's HMD-parent cells from the 2026-08-27 HMD release.

WHY.  Comparing the master against that release shows the DEATHS agree exactly for all 38
populations, but the 90+ POPULATIONS disagree for exactly two: France and the United States.
The gap is zero before about 2015 and then grows monotonically, which is the signature of HMD
having revised its old-age exposures after the master was built:

    90+ rate overstated in the master
    year      2015   2017   2019   2021   2023   2024
    USA      +0.27  +0.65  +1.38  +2.78  +5.01  +6.28   %
    France   -0.05  +0.19  +0.71  +1.71  +3.56      -   %

This matters here more than a constant offset would.  A constant proportional error in
population cancels in an excess calculation, because it scales the fitted baseline and the
observed years alike.  This error is small inside the fitting window and large in the
projection years, so it does not cancel: it inflates the observed 90+ rate relative to its own
baseline, and 90+ is the band that carries the most deaths.

It is also the root cause of the France calibration anomaly: the stored rho for France 90+ is
1.0149, outside the range of the ratio computed at every overlap year (0.967 to 0.997),
because it was fitted against the stale master rate rather than the HMD rate.

WHAT THIS DOES.  Writes a NEW master with FRATNP and USA parent cells taken from the release
(deaths and exposures, all years, all twenty bands, 90-94 through 110+ collapsed to 90+),
adding any year the release has and the master lacks.  The original is never modified.

Reads : data/master_5x1_DPM_90plus.csv, ~/Downloads/hmd_countries_20260827/
Writes: data/master_5x1_DPM_90plus_hmd20260827.csv
"""
import csv
import os
import re
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "data", "master_5x1_DPM_90plus.csv")
DST = os.path.join(ROOT, "data", "master_5x1_DPM_90plus_hmd20260827.csv")
HMD = os.path.expanduser("~/Downloads/hmd_countries_20260827")
FIX = ["FRATNP", "USA"]
COLLAPSE = {"90-94", "95-99", "100-104", "105-109", "110+"}
FINE = ["0", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44",
        "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90+"]
ASTART = {a: (0 if a == "0" else int(a.split("-")[0].replace("+", ""))) for a in FINE}


def read_hmd(path):
    """year -> band -> (female, male, total), with the oldest bands collapsed to 90+."""
    out = defaultdict(lambda: defaultdict(lambda: [0.0, 0.0, 0.0]))
    with open(path, errors="replace") as f:
        for line in f:
            p = line.split()
            if len(p) < 5 or not re.fullmatch(r"\d{4}", p[0]):
                continue
            a = "90+" if p[1] in COLLAPSE else p[1]
            if a not in FINE:
                continue
            y = int(p[0])
            for k, col in enumerate(p[2:5]):
                try:
                    out[y][a][k] += float(col)
                except ValueError:
                    pass
    return out


new = {}
for code in FIX:
    d = read_hmd(os.path.join(HMD, code, "STATS", "Deaths_5x1.txt"))
    e = read_hmd(os.path.join(HMD, code, "STATS", "Exposures_5x1.txt"))
    for y in sorted(set(d) & set(e)):
        for a in FINE:
            if a in d[y] and a in e[y]:
                new[(code, y, a)] = (d[y][a], e[y][a])
    print(f"{code}: {len(d)} death-years, {len(e)} exposure-years, "
          f"{max(set(d) & set(e))} latest")

rows = []
seen = set()
changed = defaultdict(int)
dchanged = defaultdict(int)
with open(SRC) as f:
    rd = csv.DictReader(f)
    cols = rd.fieldnames
    for r in rd:
        key = (r["location"], int(r["year"]) if r["year"].isdigit() else None, r["age"])
        if key in new:
            seen.add(key)
            (df, dm, dt), (pf, pm, pt) = new[key]
            if r["P_T"] and abs(float(r["P_T"]) - pt) > 0.5:
                changed[key[0]] += 1
            if r["D_T"] and abs(float(r["D_T"]) - dt) > 0.5:
                dchanged[key[0]] += 1
            r["D_F"], r["D_M"], r["D_T"] = f"{df:.0f}", f"{dm:.0f}", f"{dt:.0f}"
            r["P_F"], r["P_M"], r["P_T"] = f"{pf:.2f}", f"{pm:.2f}", f"{pt:.2f}"
            r["source"] = "HMD_20260827"
        rows.append(r)

# years the release has that the master lacks (France 2024)
NAME = {r["location"]: r["location_name"] for r in rows if r["location"] in FIX}
added = 0
for key in sorted(set(new) - seen):
    code, y, a = key
    (df, dm, dt), (pf, pm, pt) = new[key]
    rows.append({"location": code, "location_name": NAME.get(code, code),
                 "source": "HMD_20260827", "year": str(y), "age": a,
                 "age_start": str(ASTART[a]),
                 "D_F": f"{df:.0f}", "D_M": f"{dm:.0f}", "D_T": f"{dt:.0f}",
                 "P_F": f"{pf:.2f}", "P_M": f"{pm:.2f}", "P_T": f"{pt:.2f}"})
    added += 1

rows.sort(key=lambda r: (r["location"], int(r["year"]), ASTART.get(r["age"], 999)))
with open(DST, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)

print(f"\npopulation cells changed: {dict(changed)}")
print(f"death cells changed:      {dict(dchanged)}  (expected empty: deaths already matched)")
print(f"cells added (year absent from the master): {added}")
print(f"wrote {DST}")
if not os.path.exists(SRC):
    sys.exit("source master missing")
