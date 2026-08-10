#!/opt/homebrew/bin/python3
"""
build_2025_asia.py — extend the 5x1 masters to 2025 for TWN / HKG / JPN
from national sources (HMD stops at 2024 for all three; verified 2026-07-12).

New location codes (kept distinct from HMD codes, like the *_EUROSTAT builds):
  TWN_MOI   Taiwan  — MOI household-registration deaths by single age x sex
                      (occurrence-based, data.gov.tw 146625); population =
                      mean of end-Dec-2024 / end-Dec-2025 registry stocks
                      (dataset 77132, village files aggregated). Year 2025.
  HKG_CSD   Hong Kong — C&SD table 115-01022 registered deaths by sex x 5-yr
                      band (top 85+); population = mid-year single-age stocks
                      from table 110-01002 (thousands). Years 2020-2025
                      (overlap years kept for validation vs HMD).
                      2025 deaths provisional.
  JPN_ESTAT Japan   — MHLW vital statistics 2025 annual preliminary (gaisu,
                      Japanese nationals only) 5-yr bands to 100+; exposure =
                      mean of Jan-1-2025 (HMD Population.txt) and Jan-1-2026
                      (e-Stat monthly estimate). Year 2025.

Outputs: data/<CODE>/STATS/{Deaths,Exposures,Mx}_5x1.txt (HMD 5x1 layout,
native top band) + rows appended to master_5x1_DPM_wide.csv (native bands)
and master_5x1_DPM_90plus.csv (harmonized 20-band 0..90+ grid; HKG 85+ is
split into 85-89/90+ using same-year HMD shares, 2024 shares for 2025).
Masters are backed up to data/old/ first. Raw inputs copied to data/<CODE>/raw/.
"""

import csv, json, os, shutil, sys
from collections import defaultdict

DATA = "/Users/levitt/Dropbox/win1_DB/NewProjects25/mortality.org/HMD_Excess_Death/data"
SCRATCH = ("/private/tmp/claude-501/-Users-levitt-Dropbox-win1-DB-NewProjects25-"
           "mortality-org-HMD-Excess-Death/4618e575-676e-4bef-abdf-f884d5099d91/scratchpad")
BUILT = "12 Jul 2026"
BACKUP_TAG = "pre2025asia_2026-07-12"

BANDS20 = (["0", "1-4"] + [f"{a}-{a+4}" for a in range(5, 90, 5)] + ["90+"])
BANDS22 = (["0", "1-4"] + [f"{a}-{a+4}" for a in range(5, 100, 5)] + ["100+"])
BANDS19 = (["0", "1-4"] + [f"{a}-{a+4}" for a in range(5, 85, 5)] + ["85+"])

def age_start(band):
    return int(band.split("-")[0].rstrip("+"))

def band_of_age(a, top):
    """Map single age to band with top band 'top' (85/90/100)."""
    if a >= top: return f"{top}+"
    if a == 0: return "0"
    if a <= 4: return "1-4"
    lo = (a // 5) * 5
    return f"{lo}-{lo+4}"

def single_ages_to_bands(d, top):
    out = defaultdict(float)
    for a, v in d.items():
        out[band_of_age(a, top)] += v
    return dict(out)

# ---------------------------------------------------------------- HMD readers
def read_hmd_5x1(path):
    """-> {(year, ageband): (F, M, T)} ; '.' -> 0"""
    out = {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 5 and p[0][:4].isdigit():
                yr = int(p[0][:4])
                vals = [0.0 if v == "." else float(v) for v in p[2:5]]
                out[(yr, p[1])] = tuple(vals)
    return out

def read_hmd_pop_single(path, year):
    """Population.txt -> ({age:F}, {age:M}) for Jan-1 of `year`. 110+ -> 110."""
    F, M = {}, {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 5 and p[0].rstrip("+-").isdigit() and int(p[0].rstrip("+-")[:4]) == year:
                a = 110 if p[1] == "110+" else int(p[1])
                F[a] = F.get(a, 0.0) + (0.0 if p[2] == "." else float(p[2]))
                M[a] = M.get(a, 0.0) + (0.0 if p[3] == "." else float(p[3]))
    return F, M

# ---------------------------------------------------------------- output
def write_hmd_style(code, name_line, kind, rows, top_note):
    """rows: list of (year, band, F, M, T)"""
    os.makedirs(f"{DATA}/{code}/STATS", exist_ok=True)
    path = f"{DATA}/{code}/STATS/{kind}_5x1.txt"
    label = {"Deaths": "Deaths (period 5x1)",
             "Exposures": "Exposure to risk (period 5x1)",
             "Mx": "Death rates (period 5x1)"}[kind]
    fmt = "%.6f" if kind == "Mx" else "%.2f"
    with open(path, "w") as f:
        f.write(f"{name_line}, {label},\tBuilt: {BUILT}; {top_note}\n\n")
        f.write("  Year          Age             Female            Male           Total\n")
        for yr, band, F, M, T in rows:
            f.write(f"  {yr}       {band:<10s}{fmt % F:>14s}{fmt % M:>16s}{fmt % T:>16s}\n")
    return path

def dpm_rows(code, name, source, per_year_bands, band_order):
    """per_year_bands: {year: {band: dict(D_F..P_T)}} -> master rows"""
    rows = []
    for yr in sorted(per_year_bands):
        for band in band_order:
            r = per_year_bands[yr][band]
            m = {s: (r[f"D_{s}"] / r[f"P_{s}"] if r[f"P_{s}"] > 0 else "")
                 for s in ("F", "M", "T")}
            rows.append([code, name, source, yr, band, age_start(band)]
                        + [f"{r[f'D_{s}']:.2f}" for s in ("F", "M", "T")]
                        + [f"{r[f'P_{s}']:.2f}" for s in ("F", "M", "T")]
                        + [(f"{m[s]:.6f}" if m[s] != "" else "") for s in ("F", "M", "T")])
    return rows

# ================================================================ TAIWAN
def build_twn():
    def parse_deaths(path):
        F, M = defaultdict(float), defaultdict(float)
        with open(path, encoding="utf-8-sig") as f:
            rdr = csv.DictReader(f)
            for row in rdr:
                for a in range(0, 100):
                    F[a] += float(row[f"{a}歲_女"] or 0)
                    M[a] += float(row[f"{a}歲_男"] or 0)
                F[100] += float(row["100歲以上_女"] or 0)
                M[100] += float(row["100歲以上_男"] or 0)
        return F, M

    def parse_pop(path):
        F, M = defaultdict(float), defaultdict(float)
        with open(path, encoding="utf-8-sig") as f:
            rdr = csv.DictReader(f)
            for row in rdr:
                for a in range(0, 100):
                    F[a] += float(row[f"{a}歲-女"] or 0)
                    M[a] += float(row[f"{a}歲-男"] or 0)
                F[100] += float(row["100歲以上-女"] or 0)
                M[100] += float(row["100歲以上-男"] or 0)
        return F, M

    dF25, dM25 = parse_deaths(f"{SCRATCH}/twn_114.csv")
    dF24, dM24 = parse_deaths(f"{SCRATCH}/twn_113.csv")
    pF24, pM24 = parse_pop(f"{SCRATCH}/twn_11312.csv")   # end-Dec-2024 ~ Jan-1-2025
    pF25, pM25 = parse_pop(f"{SCRATCH}/twn_11412.csv")   # end-Dec-2025 ~ Jan-1-2026

    # --- validation: 2024 occurrence deaths vs HMD TWN 2024 (bands to 100+)
    hmd = read_hmd_5x1(f"{DATA}/TWN/STATS/Deaths_5x1.txt")
    val = []
    b24F = single_ages_to_bands(dF24, 100); b24M = single_ages_to_bands(dM24, 100)
    for band in BANDS22:
        h = [0.0, 0.0]
        for hb in ([band] if band != "100+" else ["100-104", "105-109", "110+"]):
            if (2024, hb) in hmd:
                h[0] += hmd[(2024, hb)][0]; h[1] += hmd[(2024, hb)][1]
        val.append((band, b24F.get(band, 0) + b24M.get(band, 0), h[0] + h[1]))
    tot_moi = sum(dF24.values()) + sum(dM24.values())
    tot_hmd = sum(v[0] + v[1] for (y, b), v in hmd.items() if y == 2024 and "-" in b or (y, b) == (2024, "0") or (y, b) == (2024, "110+"))
    # registry pop vs HMD Jan-1-2025 population
    hpF, hpM = read_hmd_pop_single(f"{DATA}/TWN/STATS/Population.txt", 2025)
    reg_tot = sum(pF24.values()) + sum(pM24.values())
    hmd_tot = sum(hpF.values()) + sum(hpM.values())

    per_year = {}
    exF = {a: (pF24[a] + pF25[a]) / 2 for a in pF24}
    exM = {a: (pM24[a] + pM25[a]) / 2 for a in pM24}
    for top, bands in ((100, BANDS22), (90, BANDS20)):
        bdF = single_ages_to_bands(dF25, top); bdM = single_ages_to_bands(dM25, top)
        bpF = single_ages_to_bands(exF, top);  bpM = single_ages_to_bands(exM, top)
        per_year[top] = {2025: {b: dict(D_F=bdF.get(b, 0), D_M=bdM.get(b, 0),
                                        D_T=bdF.get(b, 0) + bdM.get(b, 0),
                                        P_F=bpF.get(b, 0), P_M=bpM.get(b, 0),
                                        P_T=bpF.get(b, 0) + bpM.get(b, 0)) for b in bands}}

    name = "Taiwan (MOI registry deaths by occurrence; exposure=mean end-2024/end-2025 registry)"
    rows = [(2025, b, per_year[100][2025][b]["D_F"], per_year[100][2025][b]["D_M"],
             per_year[100][2025][b]["D_T"]) for b in BANDS22]
    write_hmd_style("TWN_MOI", name, "Deaths", rows, "to 100+; occurrence-based, 2025")
    rows = [(2025, b, per_year[100][2025][b]["P_F"], per_year[100][2025][b]["P_M"],
             per_year[100][2025][b]["P_T"]) for b in BANDS22]
    write_hmd_style("TWN_MOI", name, "Exposures", rows, "to 100+; registered population")
    rows = [(2025, b,
             per_year[100][2025][b]["D_F"] / max(per_year[100][2025][b]["P_F"], 1e-9),
             per_year[100][2025][b]["D_M"] / max(per_year[100][2025][b]["P_M"], 1e-9),
             per_year[100][2025][b]["D_T"] / max(per_year[100][2025][b]["P_T"], 1e-9)) for b in BANDS22]
    write_hmd_style("TWN_MOI", name, "Mx", rows, "to 100+")

    wide = dpm_rows("TWN_MOI", name, "MOI-Taiwan", per_year[100], BANDS22)
    p90 = dpm_rows("TWN_MOI", name, "MOI-Taiwan", per_year[90], BANDS20)
    report = {
        "moi_2025_total_deaths": sum(dF25.values()) + sum(dM25.values()),
        "moi_2024_total_deaths": tot_moi, "hmd_2024_total_deaths": tot_hmd,
        "registry_pop_end2024": reg_tot, "hmd_pop_jan1_2025": hmd_tot,
        "band_check_2024_moi_vs_hmd": val,
    }
    return wide, p90, report

# ================================================================ HONG KONG
def build_hkg():
    deaths = json.load(open(f"{SCRATCH}/hkg_deaths.json"))["dataSet"]
    pop = json.load(open(f"{SCRATCH}/hkg_pop.json"))["dataSet"]
    A_MAP = {"1-": "0"}
    years = range(2020, 2026)

    # deaths[year][sex][band]; sex '' = Total (includes unknown sex)
    D = {y: {s: defaultdict(float) for s in ("F", "M", "T")} for y in years}
    unk = {y: {s: 0.0 for s in ("F", "M", "T")} for y in years}
    for r in deaths:
        y = int(r["period"])
        if y not in D: continue
        s = {"F": "F", "M": "M", "": "T"}[r["SEX"]]
        a = r["AGE"]
        if a == "": continue                      # all-ages total row
        v = float(r["figure"])
        if a == "Unknown": unk[y][s] += v
        else:
            band = A_MAP.get(a, a.replace("85_and_over", "85+"))
            D[y][s][band] += v
    for y in years:                               # spread unknown-age pro rata
        for s in ("F", "M", "T"):
            tot = sum(D[y][s].values())
            if tot > 0 and unk[y][s] > 0:
                for b in D[y][s]: D[y][s][b] *= (1 + unk[y][s] / tot)

    # population: mid-year (period YYYY06), single age, thousands
    P = {y: {"F": defaultdict(float), "M": defaultdict(float)} for y in years}
    for r in pop:
        if r["SEX"] == "" or r["AGE"] == "": continue
        per = str(r["period"])
        if not per.endswith("06"): continue
        y = int(per[:4])
        if y not in P: continue
        a = r["AGE"]
        a = 85 if a == "85_and_over" else (0 if a == "zero" else int(a))
        P[y][r["SEX"]][a] += float(r["figure"]) * 1000.0

    per_year85 = {}
    for y in years:
        bpF = single_ages_to_bands(P[y]["F"], 85)
        bpM = single_ages_to_bands(P[y]["M"], 85)
        per_year85[y] = {b: dict(D_F=D[y]["F"].get(b, 0), D_M=D[y]["M"].get(b, 0),
                                 D_T=D[y]["T"].get(b, 0),
                                 P_F=bpF.get(b, 0), P_M=bpM.get(b, 0),
                                 P_T=bpF.get(b, 0) + bpM.get(b, 0)) for b in BANDS19}

    # 90+ grid: split 85+ using HMD same-year shares (2024 shares for 2025)
    hD = read_hmd_5x1(f"{DATA}/HKG/STATS/Deaths_5x1.txt")
    hE = read_hmd_5x1(f"{DATA}/HKG/STATS/Exposures_5x1.txt")
    def shares(tbl, y, si):
        y_ref = min(y, 2024)
        a8589 = tbl.get((y_ref, "85-89"), (0, 0, 0))[si]
        a90 = sum(tbl.get((y_ref, b), (0, 0, 0))[si]
                  for b in ("90-94", "95-99", "100-104", "105-109", "110+"))
        tot = a8589 + a90
        return (a8589 / tot, a90 / tot) if tot > 0 else (0.5, 0.5)
    per_year90 = {}
    for y in years:
        g = {}
        for b in BANDS20[:-2]:      # up to 80-84 identical
            g[b] = dict(per_year85[y][b])
        top = per_year85[y]["85+"]
        for newb in ("85-89", "90+"):
            idx = 0 if newb == "85-89" else 1
            g[newb] = {}
            for k, tbl in (("D", hD), ("P", hE)):
                for si, s in ((0, "F"), (1, "M"), (2, "T")):
                    g[newb][f"{k}_{s}"] = top[f"{k}_{s}"] * shares(tbl, y, si)[idx]
        per_year90[y] = g

    # validation vs HMD 2020-2024: total deaths + Mx 85+ (T)
    val = []
    for y in range(2020, 2025):
        hmd_tot = sum(v[2] for (yy, b), v in hD.items() if yy == y and b not in ("Total",))
        csd_tot = sum(per_year85[y][b]["D_T"] for b in BANDS19)
        h85 = sum(hD.get((y, b), (0, 0, 0))[2] for b in ("85-89", "90-94", "95-99", "100-104", "105-109", "110+"))
        e85 = sum(hE.get((y, b), (0, 0, 0))[2] for b in ("85-89", "90-94", "95-99", "100-104", "105-109", "110+"))
        c = per_year85[y]["85+"]
        val.append((y, csd_tot, hmd_tot,
                    c["D_T"] / max(c["P_T"], 1e-9), h85 / max(e85, 1e-9)))

    name = "Hong Kong (C&SD registered deaths 115-01022; exposure=mid-year pop 110-01002)"
    for kind, keys in (("Deaths", ("D_F", "D_M", "D_T")),
                       ("Exposures", ("P_F", "P_M", "P_T"))):
        rows = [(y, b, per_year85[y][b][keys[0]], per_year85[y][b][keys[1]],
                 per_year85[y][b][keys[2]]) for y in years for b in BANDS19]
        write_hmd_style("HKG_CSD", name, kind, rows, "to 85+; 2025 provisional")
    rows = [(y, b, per_year85[y][b]["D_F"] / max(per_year85[y][b]["P_F"], 1e-9),
             per_year85[y][b]["D_M"] / max(per_year85[y][b]["P_M"], 1e-9),
             per_year85[y][b]["D_T"] / max(per_year85[y][b]["P_T"], 1e-9))
            for y in years for b in BANDS19]
    write_hmd_style("HKG_CSD", name, "Mx", rows, "to 85+; 2025 provisional")

    wide = dpm_rows("HKG_CSD", name, "CSD-HK", per_year85, BANDS19)
    name90 = name + " [85+ split by HMD shares]"
    p90 = dpm_rows("HKG_CSD", name90, "CSD-HK", per_year90, BANDS20)
    report = {"validation_year_csdD_hmdD_csdMx85_hmdMx85": val,
              "deaths_2025_total": sum(per_year85[2025][b]["D_T"] for b in BANDS19)}
    return wide, p90, report

# ================================================================ JAPAN
def build_jpn(jan2026_bands=None):
    """jan2026_bands: {('F'|'M'): {band(<=85+ or finer): persons}} from e-Stat.
       If None, skip Japan (population for Jan-2026 not yet fetched)."""
    if jan2026_bands is None:
        return None, None, {"skipped": "awaiting Jan-1-2026 population"}

    raw = open(f"{SCRATCH}/jpn_deaths_2025.csv", encoding="cp932").read().splitlines()
    rows = list(csv.reader(raw))
    hdr_band, hdr_sex = rows[2], rows[3]
    # forward-fill band names
    cur = ""
    bands_by_col = []
    for i, c in enumerate(hdr_band):
        if c.strip(): cur = c.strip()
        bands_by_col.append(cur)
    nat = next(r for r in rows if any("全" in c and "国" in c for c in r[:3]))

    JMAP = {f"{a:02d}-{a+4:02d}歳": (f"{a}-{a+4}" if a else "00-04") for a in range(0, 100, 5)}
    JMAP["100歳以上"] = "100+"
    dF, dM = defaultdict(float), defaultdict(float)
    unkF = unkM = 0.0
    d0F = d0M = 0.0
    s14F, s14M = 0.0, 0.0
    for i, c in enumerate(nat):
        band = bands_by_col[i]; sex = hdr_sex[i].strip() if i < len(hdr_sex) else ""
        v = float(c.replace(",", "")) if c.strip().replace(",", "").isdigit() else None
        if v is None: continue
        if band == "0歳" and sex == "男": d0M = v
        elif band == "0歳" and sex == "女": d0F = v
        elif band in (f"{k}歳" for k in (1, 2, 3, 4)):
            if sex == "男": s14M += v
            elif sex == "女": s14F += v
        elif band in JMAP:
            b = JMAP[band]
            if sex == "男": dM[b] += v
            elif sex == "女": dF[b] += v
        elif "不詳" in band:
            if sex == "男": unkM += v
            elif sex == "女": unkF += v
    # native bands: 0, 1-4, 5-9..100+
    dF["0"], dM["0"] = d0F, d0M
    dF["1-4"], dM["1-4"] = s14F, s14M
    dF.pop("00-04", None); dM.pop("00-04", None)
    for d, u in ((dF, unkF), (dM, unkM)):        # unknown age pro rata
        tot = sum(d.values())
        if u > 0 and tot > 0:
            for b in d: d[b] *= (1 + u / tot)

    # exposure: Japanese-nationals population (consistent with gaisu deaths,
    # which cover Japanese nationals only). Jan-1-2025 and Jan-1-2026 e-Stat
    # 5-yr blocks; the 0-4 block is split 0/1-4 using the HMD Jan-2025
    # (total-population) within-block structure.
    hpF, hpM = read_hmd_pop_single(f"{DATA}/JPN/STATS/Population.txt", 2025)
    hmd25F = single_ages_to_bands(hpF, 100); hmd25M = single_ages_to_bands(hpM, 100)
    jan2025_bands = json.load(open(f"{SCRATCH}/jpn_pop_jan2025_bands.json"))
    # expand Jan-2026 to the 22-band grid using Jan-2025 within-block structure
    def expand(j26, p25):
        def rng(band):
            if band.endswith("+"): return (int(band[:-1]), 110)
            if "-" in band: lo, hi = band.split("-"); return (int(lo), int(hi))
            return (int(band), int(band))
        out = {}
        for blk, v in j26.items():
            if blk in BANDS22: out[blk] = out.get(blk, 0) + v
            else:  # coarse block (e.g. '0-4' -> '0'+'1-4'): split by Jan-2025 structure
                lo, hi = rng(blk)
                sub = [b for b in BANDS22 if rng(b)[0] >= lo and rng(b)[1] <= hi]
                tot = sum(p25[b] for b in sub)
                for b in sub: out[b] = out.get(b, 0) + v * (p25[b] / tot if tot else 0)
        return out
    p25F = expand(jan2025_bands["F"], hmd25F); p25M = expand(jan2025_bands["M"], hmd25M)
    p26F = expand(jan2026_bands["F"], hmd25F); p26M = expand(jan2026_bands["M"], hmd25M)
    exF = {b: (p25F[b] + p26F[b]) / 2 for b in BANDS22}
    exM = {b: (p25M[b] + p26M[b]) / 2 for b in BANDS22}

    per_year = {100: {2025: {}}, 90: {2025: {}}}
    for b in BANDS22:
        per_year[100][2025][b] = dict(D_F=dF.get(b, 0), D_M=dM.get(b, 0),
                                      D_T=dF.get(b, 0) + dM.get(b, 0),
                                      P_F=exF[b], P_M=exM[b], P_T=exF[b] + exM[b])
    for b in BANDS20:
        if b != "90+":
            per_year[90][2025][b] = dict(per_year[100][2025][b])
        else:
            agg = {k: sum(per_year[100][2025][bb][k] for bb in ("90-94", "95-99", "100+"))
                   for k in ("D_F", "D_M", "D_T", "P_F", "P_M", "P_T")}
            per_year[90][2025][b] = agg

    name = "Japan (MHLW vital stats 2025 prelim gaisu; nationals deaths / nationals exposure, mean Jan-2025/Jan-2026 e-Stat)"
    rows = [(2025, b, per_year[100][2025][b]["D_F"], per_year[100][2025][b]["D_M"],
             per_year[100][2025][b]["D_T"]) for b in BANDS22]
    write_hmd_style("JPN_ESTAT", name, "Deaths", rows, "to 100+; PRELIMINARY, nationals only")
    rows = [(2025, b, per_year[100][2025][b]["P_F"], per_year[100][2025][b]["P_M"],
             per_year[100][2025][b]["P_T"]) for b in BANDS22]
    write_hmd_style("JPN_ESTAT", name, "Exposures", rows, "to 100+; Japanese nationals (HMD uses total residents, ~2.4% larger)")
    rows = [(2025, b,
             per_year[100][2025][b]["D_F"] / max(per_year[100][2025][b]["P_F"], 1e-9),
             per_year[100][2025][b]["D_M"] / max(per_year[100][2025][b]["P_M"], 1e-9),
             per_year[100][2025][b]["D_T"] / max(per_year[100][2025][b]["P_T"], 1e-9)) for b in BANDS22]
    write_hmd_style("JPN_ESTAT", name, "Mx", rows, "to 100+; nationals-consistent rates, PRELIMINARY")

    wide = dpm_rows("JPN_ESTAT", name, "e-Stat", per_year[100], BANDS22)
    p90 = dpm_rows("JPN_ESTAT", name, "e-Stat", per_year[90], BANDS20)
    report = {"deaths_2025_total": sum(dF.values()) + sum(dM.values()),
              "unknown_age_redistributed": unkF + unkM}
    return wide, p90, report

# ================================================================ masters
def append_masters(new_wide, new_90, codes):
    os.makedirs(f"{DATA}/old", exist_ok=True)
    for fn in ("master_5x1_DPM_wide.csv", "master_5x1_DPM_90plus.csv"):
        bak = f"{DATA}/old/{fn.replace('.csv', '')}_{BACKUP_TAG}.csv"
        if not os.path.exists(bak):
            shutil.copy2(f"{DATA}/{fn}", bak)
    for fn, rows in (("master_5x1_DPM_wide.csv", new_wide),
                     ("master_5x1_DPM_90plus.csv", new_90)):
        path = f"{DATA}/{fn}"
        with open(path) as f:
            lines = [l for l in f if l.split(",", 1)[0] not in codes]
        with open(path, "w", newline="") as f:
            f.writelines(lines)
            w = csv.writer(f)
            w.writerows(rows)

if __name__ == "__main__":
    jan26 = None
    jp_path = f"{SCRATCH}/jpn_pop_jan2026_bands.json"
    if os.path.exists(jp_path):
        jan26 = json.load(open(jp_path))
    all_wide, all_90, reports = [], [], {}
    w, n, r = build_twn(); all_wide += w; all_90 += n; reports["TWN_MOI"] = r
    w, n, r = build_hkg(); all_wide += w; all_90 += n; reports["HKG_CSD"] = r
    w, n, r = build_jpn(jan26)
    if w: all_wide += w; all_90 += n
    reports["JPN_ESTAT"] = r
    codes = {"TWN_MOI", "HKG_CSD"} | ({"JPN_ESTAT"} if w else set())
    append_masters(all_wide, all_90, codes)
    print(json.dumps(reports, indent=1, default=str))
