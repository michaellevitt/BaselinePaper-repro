#!/opt/homebrew/bin/python3
"""
build_2025_eurostat_wk.py — extend the 5x1 masters to 2025 for 23 European
countries using Eurostat, since annual demo_magec still ends at 2024
(verified 2026-07-12) while weekly demo_r_mwk_05 has all 52 weeks of 2025.

2025 deaths  = sum of ISO-2025 weekly deaths by 5-yr band (Y_LT5..Y_GE90) x sex
               (demo_r_mwk_05, country level). UNK ages redistributed pro rata;
               the 0-4 block is split 0 vs 1-4 using the country's 2024
               demo_magec single-age structure. Top band 90+.
2025 exposure = EXTRAPOLATED mid-2025: 1.5*P(Jan-2025) - 0.5*P(Jan-2024) per
               band x sex from demo_pjan single ages (Jan-2026 not yet
               published; replace when it is).

Existing *_EUROSTAT codes (12: AUT BGR CZE DEU ESP FRA GRC HUN ITA NLD POL SVN)
get a 2025 row appended (masters + STATS txt). IRL is excluded — Irish weekly
data has no age detail.

New codes (11: BEL CHE EST FIN HRV LTU LUX LVA PRT SVK SWE _EUROSTAT) get full
2020-2025 series: 2020-2024 from demo_magec + demo_pjan (mean consecutive
Jan-1 exposure — same method as the June 2026 builds), 2025 from weekly.
These overlap HMD 2020-2024 for validation.

QA: per country, ISO-2024 weekly sum is compared to annual demo_magec 2024
(week/calendar-year mismatch + provisionality, expect within ~1%), and new
codes are compared to HMD (D_T totals and M_T by band).
"""

import csv, json, os, shutil, time, urllib.request
from collections import defaultdict

import os as _os
# Repo-relative: this script lives in code/fetch/, so the data dir is ../../data .
# (Provenance script — see code/fetch/README.md; raw inputs are not shipped with the repo.)
DATA = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))), "data")
API = "https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data"
BUILT = "12 Jul 2026"
BACKUP_TAG = "pre2025eu_2026-07-12"

EXISTING = {"AT": "AUT_EUROSTAT", "BG": "BGR_EUROSTAT", "CZ": "CZE_EUROSTAT",
            "DE": "DEU_EUROSTAT", "ES": "ESP_EUROSTAT", "FR": "FRA_EUROSTAT",
            "EL": "GRC_EUROSTAT", "HU": "HUN_EUROSTAT", "IT": "ITA_EUROSTAT",
            "NL": "NLD_EUROSTAT", "PL": "POL_EUROSTAT", "SI": "SVN_EUROSTAT"}
NEW = {"BE": ("BEL_EUROSTAT", "Belgium"), "CH": ("CHE_EUROSTAT", "Switzerland"),
       "EE": ("EST_EUROSTAT", "Estonia"), "FI": ("FIN_EUROSTAT", "Finland"),
       "HR": ("HRV_EUROSTAT", "Croatia"), "LT": ("LTU_EUROSTAT", "Lithuania"),
       "LU": ("LUX_EUROSTAT", "Luxembourg"), "LV": ("LVA_EUROSTAT", "Latvia"),
       "PT": ("PRT_EUROSTAT", "Portugal"), "SK": ("SVK_EUROSTAT", "Slovakia"),
       "SE": ("SWE_EUROSTAT", "Sweden")}
HMD_OF_NEW = {"BEL_EUROSTAT": "BEL", "CHE_EUROSTAT": "CHE", "EST_EUROSTAT": "EST",
              "FIN_EUROSTAT": "FIN", "HRV_EUROSTAT": "HRV", "LTU_EUROSTAT": "LTU",
              "LUX_EUROSTAT": "LUX", "LVA_EUROSTAT": "LVA", "PRT_EUROSTAT": "PRT",
              "SVK_EUROSTAT": "SVK", "SWE_EUROSTAT": "SWE"}

BANDS20 = (["0", "1-4"] + [f"{a}-{a+4}" for a in range(5, 90, 5)] + ["90+"])
BANDS22 = (["0", "1-4"] + [f"{a}-{a+4}" for a in range(5, 100, 5)] + ["100+"])
WK_BAND = {f"Y{a}-{a+4}": f"{a}-{a+4}" for a in range(5, 90, 5)}
WK_BAND["Y_GE90"] = "90+"

def age_start(band): return int(band.split("-")[0].rstrip("+"))

def get_json(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return json.load(r)
        except Exception:
            if i == tries - 1: raise
            time.sleep(3)

def js_records(d):
    """JSON-stat v2 -> yields ({dim: code}, value)."""
    dims = d["id"]; sizes = d["size"]
    cats = [list(d["dimension"][dim]["category"]["index"].keys()) for dim in dims]
    strides = [1] * len(sizes)
    for i in range(len(sizes) - 2, -1, -1):
        strides[i] = strides[i + 1] * sizes[i + 1]
    for k, v in d["value"].items():
        k = int(k); rec = {}
        for i, dim in enumerate(dims):
            rec[dim] = cats[i][(k // strides[i]) % sizes[i]]
        yield rec, v

def pjan_single(geo, years):
    """demo_pjan -> {year: {sex: {age_int: pop}}} single ages, Y_OPEN->100."""
    t = "&".join(f"time={y}" for y in years)
    d = get_json(f"{API}/demo_pjan?format=JSON&lang=en&geo={geo}&{t}")
    out = defaultdict(lambda: defaultdict(dict))
    for rec, v in js_records(d):
        a = rec["age"]
        if a in ("TOTAL", "UNK"): continue
        age = 0 if a == "Y_LT1" else (100 if a == "Y_OPEN" else int(a[1:]))
        out[int(rec["time"])][rec["sex"]][age] = \
            out[int(rec["time"])][rec["sex"]].get(age, 0) + v
    return out

def magec_single(geo, years):
    """demo_magec deaths -> {year: {sex: {age_int: deaths}}} (UNK kept as -1)."""
    t = "&".join(f"time={y}" for y in years)
    d = get_json(f"{API}/demo_magec?format=JSON&lang=en&geo={geo}&{t}")
    out = defaultdict(lambda: defaultdict(dict))
    for rec, v in js_records(d):
        a = rec["age"]
        if a == "TOTAL": continue
        age = -1 if a == "UNK" else (0 if a == "Y_LT1" else (100 if a == "Y_OPEN" else int(a[1:])))
        out[int(rec["time"])][rec["sex"]][age] = \
            out[int(rec["time"])][rec["sex"]].get(age, 0) + v
    return out

def weekly_bands(geo, year):
    """demo_r_mwk_05 ISO-year sum -> ({sex: {band: D}}, {sex: total_row}, nweeks).
       '0-4' kept as one block; UNK redistributed pro rata."""
    d = get_json(f"{API}/demo_r_mwk_05?format=JSON&lang=en&geo={geo}"
                 f"&sinceTimePeriod={year}-W01&untilTimePeriod={year}-W53")
    acc = {s: defaultdict(float) for s in ("F", "M", "T")}
    tot = {s: 0.0 for s in ("F", "M", "T")}
    unk = {s: 0.0 for s in ("F", "M", "T")}
    weeks = set()
    for rec, v in js_records(d):
        s = rec["sex"]; a = rec["age"]; weeks.add(rec["time"])
        if a == "TOTAL": tot[s] += v
        elif a == "UNK": unk[s] += v
        elif a == "Y_LT5": acc[s]["0-4"] += v
        else: acc[s][WK_BAND[a]] += v
    for s in acc:
        known = sum(acc[s].values())
        if known > 0 and unk[s] > 0:
            for b in acc[s]: acc[s][b] *= (1 + unk[s] / known)
    return acc, tot, len(weeks)

def bands_from_single(d, top):
    out = defaultdict(float)
    for a, v in d.items():
        if a < 0: continue
        if a >= top: out[f"{top}+"] += v
        elif a == 0: out["0"] += v
        elif a <= 4: out["1-4"] += v
        else:
            lo = (a // 5) * 5
            out[f"{lo}-{lo+4}"] += v
    return out

def read_hmd_5x1(path):
    out = {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 5 and p[0][:4].isdigit():
                out[(int(p[0][:4]), p[1])] = tuple(0.0 if v == "." else float(v) for v in p[2:5])
    return out

def write_hmd_style(code, name_line, kind, rows, top_note):
    os.makedirs(f"{DATA}/{code}/STATS", exist_ok=True)
    label = {"Deaths": "Deaths (period 5x1)",
             "Exposures": "Exposure to risk (period 5x1)",
             "Mx": "Death rates (period 5x1)"}[kind]
    fmt = "%.6f" if kind == "Mx" else "%.2f"
    with open(f"{DATA}/{code}/STATS/{kind}_5x1.txt", "w") as f:
        f.write(f"{name_line}, {label},\tBuilt: {BUILT}; {top_note}\n\n")
        f.write("  Year          Age             Female            Male           Total\n")
        for yr, band, F, M, T in rows:
            f.write(f"  {yr}       {band:<10s}{fmt % F:>14s}{fmt % M:>16s}{fmt % T:>16s}\n")

def append_txt_2025(code, rows):
    """Append 2025 rows to an existing STATS txt trio (drop old 2025 first)."""
    for kind, cols in (("Deaths", ("D_F", "D_M", "D_T")),
                       ("Exposures", ("P_F", "P_M", "P_T")), ("Mx", None)):
        path = f"{DATA}/{code}/STATS/{kind}_5x1.txt"
        with open(path) as f:
            lines = [l for l in f if not l.lstrip().startswith("2025")]
        fmt = "%.6f" if kind == "Mx" else "%.2f"
        with open(path, "w") as f:
            f.writelines(lines)
            for band in BANDS20:
                r = rows[band]
                if kind == "Mx":
                    F, M, T = (r["D_F"] / max(r["P_F"], 1e-9), r["D_M"] / max(r["P_M"], 1e-9),
                               r["D_T"] / max(r["P_T"], 1e-9))
                else:
                    F, M, T = r[cols[0]], r[cols[1]], r[cols[2]]
                f.write(f"  2025       {band:<10s}{fmt % F:>14s}{fmt % M:>16s}{fmt % T:>16s}\n")

def master_rows(code, name, source, year, bands_dict, band_order):
    rows = []
    for band in band_order:
        r = bands_dict[band]
        row = [code, name, source, year, band, age_start(band)]
        row += [f"{r[f'D_{s}']:.2f}" for s in ("F", "M", "T")]
        row += [f"{r[f'P_{s}']:.2f}" for s in ("F", "M", "T")]
        row += [(f"{r[f'D_{s}'] / r[f'P_{s}']:.6f}" if r[f"P_{s}"] > 0 else "") for s in ("F", "M", "T")]
        rows.append(row)
    return rows

def year_of_line(line):
    parts = line.split(",")
    for p in parts[3:6]:
        if p.isdigit() and len(p) == 4: return int(p)
    return None

def update_masters(rows_wide, rows_90, drop_all_codes, drop_2025_codes):
    os.makedirs(f"{DATA}/old", exist_ok=True)
    for fn, rows in (("master_5x1_DPM_wide.csv", rows_wide),
                     ("master_5x1_DPM_90plus.csv", rows_90)):
        bak = f"{DATA}/old/{fn.replace('.csv', '')}_{BACKUP_TAG}.csv"
        if not os.path.exists(bak):
            shutil.copy2(f"{DATA}/{fn}", bak)
        path = f"{DATA}/{fn}"
        with open(path) as f:
            keep = []
            for l in f:
                loc = l.split(",", 1)[0]
                if loc in drop_all_codes: continue
                if loc in drop_2025_codes and year_of_line(l) == 2025: continue
                keep.append(l)
        with open(path, "w", newline="") as f:
            f.writelines(keep)
            csv.writer(f).writerows(rows)

# ---------------------------------------------------------------- build one
def build_2025(iso2, p24, p25, magec24):
    """-> {band: {D_*,P_*}} on BANDS20 grid for 2025, + QA dict."""
    wk, tot, nweeks = weekly_bands(iso2, 2025)
    assert nweeks >= 52, f"{iso2}: only {nweeks} weeks in 2025"
    out = {}
    # infant share of 0-4 deaths from magec 2024 per sex
    def infant_share(sex):
        d = magec24[2024][sex]
        blk = sum(d.get(a, 0) for a in range(0, 5))
        return d.get(0, 0) / blk if blk > 0 else 0.8
    exp = {}
    for s in ("F", "M"):
        b24 = bands_from_single(p24[s], 90); b25 = bands_from_single(p25[s], 90)
        exp[s] = {b: max(1.5 * b25.get(b, 0) - 0.5 * b24.get(b, 0), 0.0) for b in BANDS20}
    for b in BANDS20:
        r = {}
        for s in ("F", "M"):
            if b == "0": D = wk[s]["0-4"] * infant_share(s)
            elif b == "1-4": D = wk[s]["0-4"] * (1 - infant_share(s))
            else: D = wk[s].get(b, 0.0)
            r[f"D_{s}"] = D; r[f"P_{s}"] = exp[s][b]
        r["D_T"] = r["D_F"] + r["D_M"]; r["P_T"] = r["P_F"] + r["P_M"]
        out[b] = r
    # QA: weekly 2024 vs magec 2024
    wk24, tot24, _ = weekly_bands(iso2, 2024)
    m24 = sum(v for a, v in magec24[2024]["T"].items() if a >= 0) + magec24[2024]["T"].get(-1, 0)
    qa = {"wk2025_T": tot["T"], "sum_bands_2025_T": out and sum(out[b]["D_T"] for b in BANDS20),
          "wk2024_T": tot24["T"], "magec2024_T": m24,
          "wk_vs_magec_2024_pct": (tot24["T"] / m24 - 1) * 100 if m24 else None}
    return out, qa

def build_new_code_history(iso2, code, name):
    """2020-2024 from magec+pjan (native to 100+), like the June builds."""
    magec = magec_single(iso2, range(2020, 2025))
    pjan = pjan_single(iso2, range(2020, 2026))
    hist = {}
    for y in range(2020, 2025):
        bands = {}
        for s in ("F", "M"):
            d = dict(magec[y][s]); unkv = d.pop(-1, 0)
            bd = bands_from_single(d, 100)
            known = sum(bd.values())
            if unkv > 0 and known > 0:
                for b in bd: bd[b] *= (1 + unkv / known)
            e24 = bands_from_single(pjan[y][s], 100)
            e25 = bands_from_single(pjan[y + 1][s], 100)
            for b in BANDS22:
                bands.setdefault(b, {})[f"D_{s}"] = bd.get(b, 0.0)
                bands[b][f"P_{s}"] = (e24.get(b, 0.0) + e25.get(b, 0.0)) / 2
        for b in BANDS22:
            bands[b]["D_T"] = bands[b]["D_F"] + bands[b]["D_M"]
            bands[b]["P_T"] = bands[b]["P_F"] + bands[b]["P_M"]
        hist[y] = bands
    return hist, magec, pjan

def collapse90(bands22_dict):
    out = {}
    for b in BANDS20[:-1]: out[b] = dict(bands22_dict[b])
    agg = {k: sum(bands22_dict[bb][k] for bb in ("90-94", "95-99", "100+"))
           for k in ("D_F", "D_M", "D_T", "P_F", "P_M", "P_T")}
    out["90+"] = agg
    return out

if __name__ == "__main__":
    all_wide, all_90 = [], []
    qa_all, val_all = {}, {}

    # ---- existing 12 codes: add 2025
    for iso2, code in EXISTING.items():
        pj = pjan_single(iso2, (2024, 2025))
        mg = magec_single(iso2, (2024,))
        b25, qa = build_2025(iso2, pj[2024], pj[2025], mg)
        qa_all[code] = qa
        name = f"{code[:3]} (Eurostat demo_magec/demo_pjan)"  # placeholder, fixed below
        # reuse the location_name already in the wide master
        with open(f"{DATA}/master_5x1_DPM_wide.csv") as f:
            for l in f:
                if l.startswith(code + ","):
                    name = l.split(",")[1]; break
        src = "Eurostat-wk"
        all_wide += master_rows(code, name, src, 2025, b25, BANDS20)
        all_90 += master_rows(code, name, src, 2025, b25, BANDS20)
        append_txt_2025(code, b25)
        print(f"{code}: 2025 D_T={sum(b25[b]['D_T'] for b in BANDS20):,.0f} "
              f"wk24-vs-magec24 {qa['wk_vs_magec_2024_pct']:+.2f}%", flush=True)

    # ---- new 11 codes: 2020-2024 history + 2025
    for iso2, (code, cname) in NEW.items():
        hist, magec, pjan = build_new_code_history(iso2, code, cname)
        b25, qa = build_2025(iso2, pjan[2024], pjan[2025], magec)
        qa_all[code] = qa
        name = f"{cname} (Eurostat demo_magec/demo_pjan; 2025 wk)"
        src = "Eurostat"
        rows_txt = []
        for y in range(2020, 2025):
            for b in BANDS22:
                r = hist[y][b]
                rows_txt.append((y, b, r["D_F"], r["D_M"], r["D_T"]))
            all_wide += master_rows(code, name, src, y, hist[y], BANDS22)
            all_90 += master_rows(code, name, src, y, collapse90(hist[y]), BANDS20)
        all_wide += master_rows(code, name, "Eurostat-wk", 2025, b25, BANDS20)
        all_90 += master_rows(code, name, "Eurostat-wk", 2025, b25, BANDS20)
        # STATS txt files (native: 100+ for 2020-24, 90+ for 2025)
        for kind, keys in (("Deaths", ("D_F", "D_M", "D_T")),
                           ("Exposures", ("P_F", "P_M", "P_T"))):
            rows = [(y, b, hist[y][b][keys[0]], hist[y][b][keys[1]], hist[y][b][keys[2]])
                    for y in range(2020, 2025) for b in BANDS22]
            rows += [(2025, b, b25[b][keys[0]], b25[b][keys[1]], b25[b][keys[2]])
                     for b in BANDS20]
            write_hmd_style(code, name, kind, rows,
                            "to 100+ (2020-24) / 90+ (2025, weekly-based)")
        rows = [(y, b, hist[y][b]["D_F"] / max(hist[y][b]["P_F"], 1e-9),
                 hist[y][b]["D_M"] / max(hist[y][b]["P_M"], 1e-9),
                 hist[y][b]["D_T"] / max(hist[y][b]["P_T"], 1e-9))
                for y in range(2020, 2025) for b in BANDS22]
        rows += [(2025, b, b25[b]["D_F"] / max(b25[b]["P_F"], 1e-9),
                  b25[b]["D_M"] / max(b25[b]["P_M"], 1e-9),
                  b25[b]["D_T"] / max(b25[b]["P_T"], 1e-9)) for b in BANDS20]
        write_hmd_style(code, name, "Mx", rows, "to 100+ (2020-24) / 90+ (2025)")
        # validation vs HMD
        hmd_code = HMD_OF_NEW[code]
        hD = read_hmd_5x1(f"{DATA}/{hmd_code}/STATS/Deaths_5x1.txt")
        hM = read_hmd_5x1(f"{DATA}/{hmd_code}/STATS/Mx_5x1.txt")
        val = []
        for y in range(2020, 2025):
            hmd_tot = sum(v[2] for (yy, b), v in hD.items() if yy == y)
            if hmd_tot == 0: continue
            eu_tot = sum(hist[y][b]["D_T"] for b in BANDS22)
            r7579 = (hist[y]["75-79"]["D_T"] / max(hist[y]["75-79"]["P_T"], 1e-9)) / hM[(y, "75-79")][2] if (y, "75-79") in hM and hM[(y, "75-79")][2] > 0 else None
            val.append((y, round(eu_tot / hmd_tot, 4), round(r7579, 4) if r7579 else None))
        val_all[code] = val
        print(f"{code}: 2020-24 built, 2025 D_T={sum(b25[b]['D_T'] for b in BANDS20):,.0f} "
              f"| D vs HMD {val} | wk24-vs-magec24 {qa['wk_vs_magec_2024_pct']:+.2f}%", flush=True)
        time.sleep(1)

    update_masters(all_wide, all_90,
                   drop_all_codes=set(c for c, _ in NEW.values()),
                   drop_2025_codes=set(EXISTING.values()))
    print(json.dumps({"qa": qa_all, "hmd_validation": val_all}, indent=1, default=str))
