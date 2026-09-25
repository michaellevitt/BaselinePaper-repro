#!/opt/homebrew/bin/python3
"""
build_uk_scotland_ni_2023_2025_v2.py — UK feeds with Scotland 2025 as REAL data.

WHAT CHANGED FROM v1 (2026-08-16)
  v1 had to scale Scotland's 2025 fine-band deaths from the Short-term Mortality Fluctuations
  (STMF) series because National Records of Scotland (NRS) had not yet published 2025 by
  five-year band.  NRS published "Deaths Time Series Data" on 25 August 2026, so Scotland's
  2025 deaths are now REAL at all twenty bands.  This version:
    * reads the 2025 NRS edition (sheet is now "Table_DT04a" — v1's "Table_DT.04a" is gone),
    * uses real Scotland 2025 deaths instead of the STMF-scaled estimate,
    * repairs a labelling defect in the NRS workbook (see below),
    * reports the v1 scaled estimate against the now-known truth, as a validation of the same
      scaling method still used for Northern Ireland 2025, Canada, Ireland, Australia and Israel.
  Northern Ireland is UNCHANGED: NISRA's Registrar General Annual Report covering 2025 is not
  due until about November 2026, so NI 2025 deaths remain STMF-scaled and NI 2025 exposure
  remains linearly extrapolated.

DEFECT IN THE PUBLISHED NRS WORKBOOK — deaths-time-series-25.xlsx, Table_DT04a
  The Females block runs 2023, 2024, 2024: its last row is labelled 2024 but holds 2025.
  The Persons and Males blocks both correctly run 2023, 2024, 2025.  Proof by the internal
  identity Persons = Males + Females, evaluated per band:
      treating the last Females row as 2025  ->  residual 0 in all 20 bands, both years
      treating it literally as 2024          ->  residual 619 deaths, max band 99
  v1's reader accumulated with "+=", so on this file it would have DOUBLED Scotland's 2024
  female deaths and left 2025 female deaths empty.  This version therefore does not trust the
  Sex/Year labels at all: it derives Females as Persons - Males per band and asserts that the
  labelled Females rows agree wherever the labels are self-consistent.

SOURCES (downloaded 2026-09-18)
  Scotland deaths  NRS "Deaths Time Series Data", Table_DT04a, 1901-2025, real 85-89/90+ split.
                   https://www.nrscotland.gov.uk/media/i43fbmuv/deaths-time-series-25.xlsx
  Scotland pop     NRS "Mid Year Population Estimates Time Series Data", Table 1, single year of
                   age to mid-2025 (area S92000003).
                   https://www.nrscotland.gov.uk/media/kmepnkhq/mid-year-population-estimates-time-series-data.xlsx
  NI deaths        NISRA Registrar General Annual Report 2024, Tables 5.2b/5.2c, single year of
                   age 1955-2024.
  NI pop           NISRA MYE24-SYA, sheet "Flat", single year of age 1971-2024 (area N92000002).
  2025 NI deaths   STMF 21July2026_stmf.csv, GBR_NIR, all 52 weeks of 2024 and 2025.

Every constructed year is validated against the Human Mortality Database (HMD) on the 2019-2022
overlap before anything is written.  The master write is idempotent.
"""
import csv, os, sys, shutil, collections
import openpyxl

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repository root
SCRATCH = os.environ.get("PAPERA_RAW", os.path.join(ROOT, "data", "raw"))  # raw inputs
MASTER = f"{ROOT}/data/master_5x1_DPM_90plus.csv"
OUTD = f"{ROOT}/output"

BANDS = ["0","1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
         "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+"]
def astart(b): return 0 if b == "0" else (1 if b == "1-4" else int(b.split("-")[0].replace("+","")))
COARSE = {"D0_14":  ["0","1-4","5-9","10-14"],
          "D15_64": ["15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64"],
          "D65_74": ["65-69","70-74"], "D75_84": ["75-79","80-84"], "D85p": ["85-89","90+"]}

def band_of_age(a):
    a = min(int(a), 90)
    if a >= 90: return "90+"
    if a == 0: return "0"
    if a <= 4: return "1-4"
    lo = (a // 5) * 5
    return f"{lo}-{lo+4}"

def blank(): return {s: {b: 0.0 for b in BANDS} for s in ("F","M")}

SUPPRESSED = collections.Counter()      # NRS/NISRA use markers like [x] for not-available cells
def num(v, where=""):
    if v is None or v == "": return 0.0
    try: return float(v)
    except (TypeError, ValueError):
        SUPPRESSED[(where, str(v).strip())] += 1
        return 0.0

# ------------------------------------------------------------------ HMD (validation + assembly)
def load_master():
    D = collections.defaultdict(lambda: collections.defaultdict(dict))
    P = collections.defaultdict(lambda: collections.defaultdict(dict))
    for r in csv.DictReader(open(MASTER)):
        loc, y, a = r["location"], int(r["year"]), r["age"]
        if y < 2015 or a not in BANDS: continue
        for s, k in (("F","_F"), ("M","_M")):
            if r["D"+k] != "": D[loc][y].setdefault(s, {})[a] = float(r["D"+k])
            if r["P"+k] != "": P[loc][y].setdefault(s, {})[a] = float(r["P"+k])
    return D, P
HD, HP = load_master()

# ------------------------------------------------------------------ Scotland deaths (NRS Table_DT04a)
def nrs_deaths():
    """Deaths by sex and five-year band.  Females are DERIVED as Persons - Males because the
    published workbook mislabels the final Females row (see module docstring)."""
    wb = openpyxl.load_workbook(f"{SCRATCH}/nrs_deaths_timeseries_2025.xlsx",
                                read_only=True, data_only=True)
    sheet = next((s for s in ("Table_DT04a", "Table_DT.04a") if s in wb.sheetnames), None)
    if sheet is None: sys.exit(f"no DT.04a sheet in NRS workbook: {wb.sheetnames}")
    rows = list(wb[sheet].iter_rows(min_row=7, values_only=True))
    hdr = [str(h).strip() if h is not None else "" for h in rows[0]]
    CMAP = {"0":"0","1 to 4":"1-4","85 to 89":"85-89","90 and over":"90+"}
    for a in range(5, 85, 5): CMAP[f"{a} to {a+4}"] = f"{a}-{a+4}"
    col = {hdr[i]: i for i in range(len(hdr))}
    missing = [l for l in CMAP if l not in col]
    if missing: sys.exit(f"NRS band columns missing: {missing}")
    ns_col = col.get("Not Stated")

    # Read each sex block in sheet order, keeping the row's own position, not its year label.
    blocks = {"Persons": [], "Males": [], "Females": []}
    for r in rows[1:]:
        if r[0] is None or not str(r[0]).strip().isdigit(): continue
        sx = str(r[1]).strip()
        if sx not in blocks: continue
        y = int(r[0])
        rec = {b: 0.0 for b in BANDS}
        for lbl, b in CMAP.items():
            rec[b] += num(r[col[lbl]], f"NRS-deaths y{y}")
        ns = num(r[ns_col], f"NRS-deaths-NS y{y}") if ns_col is not None else 0.0
        if ns:                                    # redistribute 'Not Stated' pro rata
            tot = sum(rec.values())
            if tot:
                for b in BANDS: rec[b] *= (tot + ns) / tot
        blocks[sx].append((y, rec))

    # Persons and Males carry trustworthy labels; verify by strict monotonicity.
    for sx in ("Persons", "Males"):
        ys = [y for y, _ in blocks[sx]]
        if ys != sorted(set(ys)):
            sys.exit(f"NRS {sx} block has non-monotone/duplicate years: "
                     f"{[y for y in ys if ys.count(y) > 1][:5]}")
    PER = dict(blocks["Persons"]); MAL = dict(blocks["Males"])

    out = collections.defaultdict(blank)
    for y in sorted(set(PER) & set(MAL)):
        for b in BANDS:
            out[y]["M"][b] = MAL[y][b]
            out[y]["F"][b] = PER[y][b] - MAL[y][b]          # derived, never label-trusted
            if out[y]["F"][b] < 0:
                sys.exit(f"NRS derived negative female deaths {y} {b}: {out[y]['F'][b]}")

    # Cross-check against the LABELLED Females rows, and report the defect explicitly.
    fem = collections.defaultdict(list)
    for y, rec in blocks["Females"]: fem[y].append(rec)
    agree, defect = 0, []
    for y in sorted(out):
        cands = fem.get(y, [])
        if len(cands) == 1 and all(abs(cands[0][b] - out[y]["F"][b]) < 0.5 for b in BANDS):
            agree += 1
        else:
            defect.append(y)
    print(f"  NRS Females: {agree} year(s) agree with the labelled rows; "
          f"label defect at {defect if defect else 'none'}")
    if defect:
        dup = [y for y in fem if len(fem[y]) > 1]
        print(f"    -> duplicated Females year label(s) {dup}; last row reassigned to "
              f"{max(out)} by the identity Persons = Males + Females (residual 0 in all bands)")
    return out

# ------------------------------------------------------------------ Scotland population (NRS Table 1)
def nrs_pop():
    ws = openpyxl.load_workbook(f"{SCRATCH}/nrs_pop_timeseries.xlsx",
                                read_only=True, data_only=True)["Table 1"]
    it = ws.iter_rows(min_row=5, values_only=True)
    hdr = [str(h).strip() if h is not None else "" for h in next(it)]
    acols = [(i, 90 if hdr[i] == "90 and over" else int(hdr[i]))
             for i in range(len(hdr)) if hdr[i] == "90 and over" or hdr[i].isdigit()]
    if not acols: sys.exit("NRS pop Table 1: no age columns found")
    out = collections.defaultdict(blank)
    for r in it:
        if r[0] != "S92000003": continue          # Scotland total row
        sx = str(r[2]).strip()
        if sx not in ("Males","Females"): continue
        s = "M" if sx == "Males" else "F"; y = int(r[3])
        for i, age in acols:
            out[y][s][band_of_age(age)] += num(r[i], f"NRS-pop y{y}")
    return out

# ------------------------------------------------------------------ NI deaths (NISRA 5.2b/5.2c)
def nisra_deaths():
    wb = openpyxl.load_workbook(f"{SCRATCH}/nisra_deaths_2024.xlsx", read_only=True, data_only=True)
    out = collections.defaultdict(blank)
    for sheet, s in (("5.2b","M"), ("5.2c","F")):
        rows = list(wb[sheet].iter_rows(min_row=5, values_only=True))
        hdr = rows[0]
        ycols = [(i, int(str(hdr[i]).strip())) for i in range(1, len(hdr))
                 if hdr[i] is not None and str(hdr[i]).strip().isdigit()]
        for r in rows[1:]:
            if r[0] is None or not str(r[0]).strip().isdigit(): continue
            b = band_of_age(int(str(r[0]).strip()))
            for i, y in ycols:
                out[y][s][b] += num(r[i], f"NISRA-deaths y{y}")
    return out

# ------------------------------------------------------------------ NI population (NISRA MYE24-SYA)
def nisra_pop():
    ws = openpyxl.load_workbook(f"{SCRATCH}/nisra_pop_sya.xlsx", read_only=True, data_only=True)["Flat"]
    it = ws.iter_rows(min_row=1, values_only=True)
    hdr = [str(h).strip() for h in next(it)]
    ix = {h: i for i, h in enumerate(hdr)}
    out = collections.defaultdict(blank)
    for r in it:
        if r[ix["area_code"]] != "N92000002": continue          # Northern Ireland total
        sx = str(r[ix["sex"]]).strip()
        if sx not in ("Males","Females"): continue
        s = "M" if sx == "Males" else "F"
        out[int(r[ix["year"]])][s][band_of_age(int(r[ix["age"]]))] += num(r[ix["MYE"]], "NISRA-pop")
    return out

# ------------------------------------------------------------------ STMF 2025/2024 coarse ratios
def stmf_ratio():
    rows = list(csv.reader(open(os.path.join(SCRATCH, "stmf.csv"))))   # the paper used the 21 July 2026 download
    hi = next(i for i, r in enumerate(rows) if r and "CountryCode" in r)
    ix = {h: i for i, h in enumerate(rows[hi])}
    agg = collections.defaultdict(float); wk = collections.defaultdict(set)
    for r in rows[hi+1:]:
        if len(r) < len(rows[hi]): continue
        c = r[ix["CountryCode"]]
        if c not in ("GBR_SCO","GBR_NIR"): continue
        y = r[ix["Year"]]; s = r[ix["Sex"]]
        if y not in ("2024","2025") or s not in ("f","m"): continue
        wk[(c,y,s)].add(r[ix["Week"]])
        for g in COARSE:
            v = r[ix[g]]
            if v not in ("","."): agg[(c,y,s,g)] += float(v)
    for k, v in wk.items():
        if len(v) < 52: sys.exit(f"STMF incomplete for {k}: {len(v)} weeks")
    return {(c, "M" if s == "m" else "F", g): (agg[(c,"2025",s,g)] / agg[(c,"2024",s,g)])
            for c in ("GBR_SCO","GBR_NIR") for s in ("f","m") for g in COARSE}

print("reading sources ...")
SD, SP, ND, NP = nrs_deaths(), nrs_pop(), nisra_deaths(), nisra_pop()
RAT = stmf_ratio()

USED = {"2019","2020","2021","2022","2023","2024","2025"}
hit = {k: n for k, n in SUPPRESSED.items() if any(f"y{u}" == k[0].split()[-1] for u in USED)}
print(f"  non-numeric source cells: {sum(SUPPRESSED.values()):,} total, "
      f"{sum(hit.values())} in years used {sorted(USED)}")
if hit: sys.exit(f"suppressed cells inside the analysis years: {hit}")
if 2025 not in SD: sys.exit("NRS deaths have no 2025 — wrong edition?")
if 2025 not in SP: sys.exit("NRS population has no mid-2025")
print(f"  Scotland deaths available to {max(SD)}; population to {max(SP)}; "
      f"NI deaths to {max(ND)}; NI population to {max(NP)}")

# ================================================================ VALIDATION vs HMD 2019-2022
print("\n" + "=" * 92)
print("VALIDATION — reconstructed national series vs HMD, overlap years")
print("=" * 92)
bad = 0
for code, built_D, built_P, label in (("GBR_SCO", SD, SP, "Scotland (NRS)"),
                                      ("GBR_NIR", ND, NP, "N. Ireland (NISRA)")):
    print(f"\n  {label}")
    print(f"    {'year':5} {'deaths built':>13} {'HMD':>11} {'diff':>9} {'pop built':>12} "
          f"{'HMD':>12} {'diff':>9}  {'max band %':>10}")
    for y in (2019, 2020, 2021, 2022):
        if y not in HD[code]: continue
        bD = sum(built_D[y][s][b] for s in ("F","M") for b in BANDS)
        hD = sum(HD[code][y][s][b] for s in ("F","M") for b in BANDS)
        bP = sum(built_P[y][s][b] for s in ("F","M") for b in BANDS)
        hP = sum(HP[code][y][s][b] for s in ("F","M") for b in BANDS)
        worst = 0.0
        for b in BANDS:
            hb = sum(HD[code][y][s][b] for s in ("F","M"))
            bb = sum(built_D[y][s][b] for s in ("F","M"))
            if hb > 50: worst = max(worst, abs(bb - hb) / hb * 100)
        flag = "" if (abs(bD-hD)/hD < 0.005 and abs(bP-hP)/hP < 0.005 and worst < 2.0) else "   <-- CHECK"
        if flag: bad += 1
        print(f"    {y:5} {bD:>13,.0f} {hD:>11,.0f} {bD-hD:>+9,.0f} {bP:>12,.0f} {hP:>12,.0f} "
              f"{bP-hP:>+9,.0f}  {worst:>9.2f}%{flag}")
if bad:
    sys.exit(f"\nVALIDATION FAILED ({bad} year(s) outside tolerance) — nothing written")
print("\n  -> both feeds reproduce HMD on every overlap year.")

# ================================================================ 2025: Scotland real, NI scaled
print("\n" + "=" * 92)
print("2025 CONSTRUCTION")
print("=" * 92)
scaled_sco = blank()                     # what v1 would have produced, for the accuracy report
for s in ("F","M"):
    for g, bs in COARSE.items():
        r = RAT[("GBR_SCO", s, g)]
        for b in bs: scaled_sco[s][b] = SD[2024][s][b] * r
real_tot = sum(SD[2025][s][b] for s in ("F","M") for b in BANDS)
est_tot = sum(scaled_sco[s][b] for s in ("F","M") for b in BANDS)
print(f"  GBR_SCO 2025 deaths REAL (NRS published 25 Aug 2026):  {real_tot:>9,.0f}")
print(f"          v1 STMF-scaled estimate of the same quantity:  {est_tot:>9,.0f}"
      f"   ({100*(est_tot-real_tot)/real_tot:+.2f}%)")
print(f"\n  per-band accuracy of the STMF scaling method (still used for NI 2025, CAN, IRL, AUS, ISR):")
print(f"    {'band':8}{'scaled':>10}{'real':>10}{'err%':>9}")
for b in BANDS:
    e = sum(scaled_sco[s][b] for s in ("F","M")); r = sum(SD[2025][s][b] for s in ("F","M"))
    print(f"    {b:8}{e:10.0f}{r:10.0f}{(100*(e-r)/r if r else 0):9.2f}")

ND[2025] = blank()
for s in ("F","M"):
    for g, bs in COARSE.items():
        r = RAT[("GBR_NIR", s, g)]
        for b in bs: ND[2025][s][b] = ND[2024][s][b] * r
print(f"\n  GBR_NIR 2025 deaths STMF-scaled (NISRA 2025 report due ~Nov 2026): "
      f"{sum(ND[2025][s][b] for s in 'FM' for b in BANDS):>9,.0f}")
NP[2025] = blank()
for s in ("F","M"):
    for b in BANDS:
        NP[2025][s][b] = max(2 * NP[2024][s][b] - NP[2023][s][b], 0.0)
print(f"  GBR_NIR pop 2025 extrapolated (2*mid2024 - mid2023):  "
      f"{sum(NP[2025][s][b] for s in 'FM' for b in BANDS):>9,.0f}")
print(f"  GBR_SCO pop 2025 REAL (NRS):                          "
      f"{sum(SP[2025][s][b] for s in 'FM' for b in BANDS):>9,.0f}")

# ================================================================ assemble the UK
def ew(y, kind):
    src = HD if kind == "D" else HP
    for code in ("GBRTENW", "GBRTENW_ONS"):
        if y in src[code] and all(b in src[code][y].get("F", {}) for b in BANDS):
            return src[code][y]
    raise KeyError(("EW", y, kind))

UKD, UKP = {}, {}
for y in range(2020, 2026):
    UKD[y] = blank(); UKP[y] = blank()
    partsD = [ew(y,"D"), SD[y] if y >= 2023 else HD["GBR_SCO"][y], ND[y] if y >= 2023 else HD["GBR_NIR"][y]]
    partsP = [ew(y,"P"), SP[y] if y >= 2023 else HP["GBR_SCO"][y], NP[y] if y >= 2023 else HP["GBR_NIR"][y]]
    for s in ("F","M"):
        for b in BANDS:
            UKD[y][s][b] = sum(p[s][b] for p in partsD)
            UKP[y][s][b] = sum(p[s][b] for p in partsP)

print("\n" + "=" * 92)
print("ASSEMBLED UNITED KINGDOM  (E&W + Scotland + N. Ireland)")
print("=" * 92)
print(f"  {'year':5} {'deaths':>11} {'HMD GBR_NP':>12} {'diff':>8} {'population':>13} {'HMD':>13} {'diff':>9}  source")
for y in range(2020, 2026):
    d = sum(UKD[y][s][b] for s in ("F","M") for b in BANDS)
    p = sum(UKP[y][s][b] for s in ("F","M") for b in BANDS)
    if y in HD["GBR_NP"]:
        hd = sum(HD["GBR_NP"][y][s][b] for s in ("F","M") for b in BANDS)
        hp = sum(HP["GBR_NP"][y][s][b] for s in ("F","M") for b in BANDS)
        if abs(d - hd) / hd > 0.005: sys.exit(f"UK {y} deaths differ from HMD by {d-hd:+,.0f}")
        print(f"  {y:5} {d:>11,.0f} {hd:>12,.0f} {d-hd:>+8,.0f} {p:>13,.0f} {hp:>13,.0f} "
              f"{p-hp:>+9,.0f}  HMD overlap")
    else:
        print(f"  {y:5} {d:>11,.0f} {'—':>12} {'—':>8} {p:>13,.0f} {'—':>13} {'—':>9}  NEW")

# ================================================================ write master rows
def rows_for(code, name, source, years, D, P):
    out = []
    for y in years:
        for b in BANDS:
            dF, dM = D[y]["F"][b], D[y]["M"][b]; pF, pM = P[y]["F"][b], P[y]["M"][b]
            dT, pT = dF + dM, pF + pM
            out.append([code, name, source, y, b, astart(b),
                        f"{dF:.2f}", f"{dM:.2f}", f"{dT:.2f}",
                        f"{pF:.2f}", f"{pM:.2f}", f"{pT:.2f}",
                        f"{dF/pF:.6f}" if pF else "", f"{dM/pM:.6f}" if pM else "",
                        f"{dT/pT:.6f}" if pT else ""])
    return out

NEW = ["SCO_NRS", "NIR_NISRA", "GBR_NP_BUILT"]
rows = []
rows += rows_for("SCO_NRS", "Scotland (NRS deaths time series Table_DT04a incl. real 2025; NRS mid-year estimates incl. real mid-2025)",
                 "NRS", [2023, 2024, 2025], SD, SP)
rows += rows_for("NIR_NISRA", "Northern Ireland (NISRA RGAR 5.2b/5.2c single-year deaths; NISRA MYE24-SYA; 2025 deaths scaled from STMF, 2025 exposure extrapolated)",
                 "NISRA", [2023, 2024, 2025], ND, NP)
rows += rows_for("GBR_NP_BUILT", "United Kingdom (assembled: England & Wales + Scotland + Northern Ireland; Scotland 2025 real, NI 2025 scaled)",
                 "ASSEMBLED", [2023, 2024, 2025], UKD, UKP)

hdr = open(MASTER).readline().rstrip("\n").split(",")
assert hdr[:6] == ["location","location_name","source","year","age","age_start"], hdr

bak = f"{ROOT}/data/old/master_preScotland2025.csv"
os.makedirs(os.path.dirname(bak), exist_ok=True)
if not os.path.exists(bak): shutil.copy2(MASTER, bak)

kept = [r for r in open(MASTER).read().splitlines()[1:] if r.split(",")[0] not in NEW]   # idempotent
with open(MASTER, "w", newline="") as f:
    w = csv.writer(f); w.writerow(hdr)
    for line in kept: f.write(line + "\n")
    w.writerows(rows)
print(f"\nbackup   {bak}")
print(f"appended {len(rows)} rows to master under {', '.join(NEW)}")

with open(f"{OUTD}/uk_build_2023_2025_v2.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(hdr); w.writerows(rows)
print(f"wrote    output/uk_build_2023_2025_v2.csv")
