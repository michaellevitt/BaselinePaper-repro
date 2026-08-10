#!/opt/homebrew/bin/python3
"""
build_2025_others.py — extend the 5x1 masters to 2025 for the remaining
locations, from national sources fetched 2026-07-12. Companion to
build_2025_asia.py (TWN/HKG/JPN) and build_2025_eurostat_wk.py (Europe).

Run with country args: build_2025_others.py CHL [KOR USA EW DEU ...]
Each country function returns (wide_rows, p90_rows, report) and registers
which existing master rows to drop (idempotent re-runs).

CHL_DEIS — deaths: DEIS preliminary microdata (individual records, FECHA_DEF
  2024-2025; EDAD_TIPO 2/3/4 -> age 0, unknown redistributed pro rata; sex
  'Indeterminado' allocated pro rata F/M within band). Exposure: INE base-2017
  projections, June-30 single age x sex (exactly mid-year). Years 2024 (HMD
  overlap validation) + 2025. Native bands to 100+.
"""

import csv, json, os, sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_2025_eurostat_wk as B

SCRATCH = ("/private/tmp/claude-501/-Users-levitt-Dropbox-win1-DB-NewProjects25-"
           "mortality-org-HMD-Excess-Death/4618e575-676e-4bef-abdf-f884d5099d91/scratchpad")

def band_of_age(a, top):
    if a >= top: return f"{top}+"
    if a == 0: return "0"
    if a <= 4: return "1-4"
    lo = (a // 5) * 5
    return f"{lo}-{lo+4}"

def finish_bands(D, P, bands):
    out = {}
    for b in bands:
        r = {}
        for s in ("F", "M"):
            r[f"D_{s}"] = D[s].get(b, 0.0); r[f"P_{s}"] = P[s].get(b, 0.0)
        r["D_T"] = r["D_F"] + r["D_M"]; r["P_T"] = r["P_F"] + r["P_M"]
        out[b] = r
    return out

def collapse90(b22):
    out = {b: dict(b22[b]) for b in B.BANDS20[:-1]}
    out["90+"] = {k: sum(b22[bb][k] for bb in ("90-94", "95-99", "100+"))
                  for k in ("D_F", "D_M", "D_T", "P_F", "P_M", "P_T")}
    return out

def emit(code, name, src, per_year, native_bands, top_note):
    """per_year: {year: bands_dict on native_bands}. Returns wide, p90 rows."""
    wide, p90 = [], []
    for y in sorted(per_year):
        wide += B.master_rows(code, name, src, y, per_year[y], native_bands)
        g = collapse90(per_year[y]) if native_bands == B.BANDS22 else per_year[y]
        p90 += B.master_rows(code, name, src, y, g, B.BANDS20)
    for kind, keys in (("Deaths", ("D_F", "D_M", "D_T")),
                       ("Exposures", ("P_F", "P_M", "P_T"))):
        rows = [(y, b, per_year[y][b][keys[0]], per_year[y][b][keys[1]],
                 per_year[y][b][keys[2]]) for y in sorted(per_year) for b in native_bands]
        B.write_hmd_style(code, name, kind, rows, top_note)
    rows = [(y, b, per_year[y][b]["D_F"] / max(per_year[y][b]["P_F"], 1e-9),
             per_year[y][b]["D_M"] / max(per_year[y][b]["P_M"], 1e-9),
             per_year[y][b]["D_T"] / max(per_year[y][b]["P_T"], 1e-9))
            for y in sorted(per_year) for b in native_bands]
    B.write_hmd_style(code, name, "Mx", rows, top_note)
    return wide, p90

# ================================================================ CHILE
def build_chl():
    # ---- deaths from microdata
    D = {y: {"F": defaultdict(float), "M": defaultdict(float)} for y in (2024, 2025)}
    unk = {y: {"F": 0.0, "M": 0.0} for y in (2024, 2025)}
    indet = {y: defaultdict(float) for y in (2024, 2025)}   # sex-indeterminate by band
    path = f"{SCRATCH}/chl_deis_micro/DEFUNCIONES_FUENTE_DEIS_2024_2026_07072026.csv"
    with open(path, encoding="latin-1") as f:
        rdr = csv.reader(f, delimiter=";")
        hdr = next(rdr)
        i_fecha = hdr.index("FECHA_DEF"); i_sex = hdr.index("SEXO_NOMBRE")
        i_tipo = hdr.index("EDAD_TIPO"); i_cant = hdr.index("EDAD_CANT")
        for row in rdr:
            y = int(row[i_fecha][:4])
            if y not in D: continue
            tipo = row[i_tipo].strip(); cant = row[i_cant].strip()
            sex = row[i_sex].strip()
            s = {"Mujer": "F", "Hombre": "M"}.get(sex)
            if tipo == "1" and cant.isdigit():
                age = min(int(cant), 105)          # cap implausible >105
                band = band_of_age(age, 100)
                if s: D[y][s][band] += 1
                else: indet[y][band] += 1
            elif tipo in ("2", "3", "4"):
                if s: D[y][s]["0"] += 1
                else: indet[y]["0"] += 1
            else:                                   # unknown age
                if s: unk[y][s] += 1
    for y in (2024, 2025):
        for b, v in indet[y].items():               # split indeterminate-sex pro rata
            fm = D[y]["F"].get(b, 0) + D[y]["M"].get(b, 0)
            for s in ("F", "M"):
                D[y][s][b] += v * (D[y][s].get(b, 0) / fm if fm else 0.5)
        for s in ("F", "M"):                        # unknown age pro rata
            tot = sum(D[y][s].values())
            if unk[y][s] > 0 and tot > 0:
                for b in D[y][s]: D[y][s][b] *= (1 + unk[y][s] / tot)

    # ---- INE June-30 population, single age blocks (AMBOS/HOMBRES/MUJERES)
    P = {y: {"F": defaultdict(float), "M": defaultdict(float)} for y in (2024, 2025)}
    block = None; col_of_year = {}
    with open(f"{SCRATCH}/chl_ine_proyecciones_1992-2050_base2017_pais.csv",
              encoding="latin-1") as f:
        for row in csv.reader(f, delimiter=";"):
            if not row: continue
            c0 = row[0].strip()
            if "AMBOS SEXOS" in c0: block = None; continue
            if c0 == "HOMBRES": block = "M"; continue
            if c0 == "MUJERES": block = "F"; continue
            if c0 == "EDAD":
                col_of_year = {int(v): i for i, v in enumerate(row) if v.strip().isdigit()}
                continue
            if block and (c0.isdigit() or c0 == "100+"):
                age = 100 if c0 == "100+" else int(c0)
                for y in (2024, 2025):
                    v = row[col_of_year[y]].replace(".", "").strip()
                    if v: P[y][block][band_of_age(age, 100)] += float(v)

    per_year = {y: finish_bands(D[y], P[y], B.BANDS22) for y in (2024, 2025)}

    # ---- validation vs HMD 2024
    hD = B.read_hmd_5x1(f"{B.DATA}/CHL/STATS/Deaths_5x1.txt")
    hM = B.read_hmd_5x1(f"{B.DATA}/CHL/STATS/Mx_5x1.txt")
    hmd24 = sum(v[2] for (yy, b), v in hD.items() if yy == 2024)
    val = {"deis_2024_D_T": round(sum(per_year[2024][b]["D_T"] for b in B.BANDS22)),
           "hmd_2024_D_T": round(hmd24),
           "mx7579_ratio_2024": round((per_year[2024]["75-79"]["D_T"] /
                                       per_year[2024]["75-79"]["P_T"]) /
                                      hM[(2024, "75-79")][2], 4) if (2024, "75-79") in hM else None,
           "deis_2025_D_T": round(sum(per_year[2025][b]["D_T"] for b in B.BANDS22))}

    name = "Chile (DEIS preliminary microdata; exposure=INE base-2017 June-30 projections)"
    wide, p90 = emit("CHL_DEIS", name, "DEIS-INE", per_year, B.BANDS22,
                     "to 100+; 2024-2025 preliminary registry")
    return wide, p90, val, {"CHL_DEIS"}

# ================================================================ GERMANY
def build_deu():
    """2025 deaths from Destatis GENESIS 12613-0003 (year x sex x single age,
       provisional 'e'); append to existing DEU_EUROSTAT (2020-2024). Exposure =
       extrapolated mid-2025 from Eurostat demo_pjan (1.5*Jan2025 - 0.5*Jan2024)."""
    d = json.load(open(f"{SCRATCH}/deu_12613-0003_data.json"))["data"][0]
    sizes = d["size"]; val = d["value"]
    strides = [1] * 6
    for i in range(4, -1, -1): strides[i] = strides[i + 1] * sizes[i + 1]
    GES = {"M": 0, "F": 1}; JAHR2025 = 69
    def cell(s, alt): return val[GES[s] * strides[3] + JAHR2025 * strides[4] + alt * strides[5]] or 0.0
    D = {"F": defaultdict(float), "M": defaultdict(float)}
    for s in ("F", "M"):
        for a in range(0, 100):
            D[s][band_of_age(a, 100)] += cell(s, a)
        D[s]["100+"] += cell(s, 100)                 # ALT100UM
        unk = cell(s, 101)                            # ALTNN
        tot = sum(D[s].values())
        if unk > 0 and tot > 0:
            for b in D[s]: D[s][b] *= (1 + unk / tot)

    pj = B.pjan_single("DE", (2024, 2025))
    P = {"F": defaultdict(float), "M": defaultdict(float)}
    for s in ("F", "M"):
        b24 = _bands100(pj[2024][s]); b25 = _bands100(pj[2025][s])
        for b in B.BANDS22:
            P[s][b] = max(1.5 * b25.get(b, 0) - 0.5 * b24.get(b, 0), 0.0)
    per_year = {2025: finish_bands(D, P, B.BANDS22)}

    hD = B.read_hmd_5x1(f"{B.DATA}/DEUTNP/STATS/Deaths_5x1.txt")
    hmd_last = max((y for (y, _) in hD), default=0)
    val_rep = {"deu_2025_D_T": round(sum(per_year[2025][b]["D_T"] for b in B.BANDS22)),
               "destatis_2025_total_raw": 1006588,
               "hmd_DEUTNP_last_year": hmd_last}

    name = ("Germany (Eurostat demo_magec/demo_pjan; 2025 Destatis GENESIS 12613-0003 prov.)")
    # keep 2020-2024 rows already in masters/txt; append only 2025
    wide = B.master_rows("DEU_EUROSTAT", name, "Destatis", 2025, per_year[2025], B.BANDS22)
    p90 = B.master_rows("DEU_EUROSTAT", name, "Destatis", 2025,
                        collapse90(per_year[2025]), B.BANDS20)
    _append_txt_native("DEU_EUROSTAT", 2025, per_year[2025], B.BANDS22,
                       "Deaths to 100+; 2025 Destatis provisional")
    return wide, p90, val_rep, set()  # drop handled via drop_2025

def _bands100(single):
    out = defaultdict(float)
    for a, v in single.items():
        if a < 0: continue
        out[band_of_age(a, 100)] += v
    return out

def _append_txt_native(code, year, bands_dict, band_order, note):
    for kind, keys in (("Deaths", ("D_F", "D_M", "D_T")),
                       ("Exposures", ("P_F", "P_M", "P_T")), ("Mx", None)):
        path = f"{B.DATA}/{code}/STATS/{kind}_5x1.txt"
        if not os.path.exists(path): continue
        with open(path) as f:
            lines = [l for l in f if not l.lstrip().startswith(str(year))]
        fmt = "%.6f" if kind == "Mx" else "%.2f"
        with open(path, "w") as f:
            f.writelines(lines)
            for b in band_order:
                r = bands_dict[b]
                if kind == "Mx":
                    vals = (r["D_F"] / max(r["P_F"], 1e-9), r["D_M"] / max(r["P_M"], 1e-9),
                            r["D_T"] / max(r["P_T"], 1e-9))
                else:
                    vals = (r[keys[0]], r[keys[1]], r[keys[2]])
                f.write(f"  {year}       {b:<10s}{fmt % vals[0]:>14s}{fmt % vals[1]:>16s}{fmt % vals[2]:>16s}\n")

# ================================================================ ENGLAND & WALES
def build_ew():
    """2025 deaths from ONS weekly age-sex CSV (England+Wales, Registrations,
       full 52 weeks); append to existing GBRTENW_ONS (2020-2024, top 90+).
       Exposure = extrapolated mid-2025 = 2*mid2024 - mid2023 (Nomis NM_2002_1)."""
    import csv as _csv
    AMAP = {"0-1": "0", "1-4": "1-4"}
    for a in range(5, 90, 5): AMAP[f"{a}-{a+4}"] = f"{a}-{a+4}"
    for a in ("90-94", "95-99", "100+"): AMAP[a] = "90+"
    D = {"F": defaultdict(float), "M": defaultdict(float)}
    with open(f"{SCRATCH}/ew_deaths_weekly_agesex_2025_v45.csv") as f:
        for r in _csv.DictReader(f):
            if r["Geography"] not in ("England", "Wales"): continue
            if r["RegistrationOrOccurrence"] != "Registrations": continue
            if r["Sex"] not in ("Female", "Male"): continue
            if r["AgeGroups"] in ("All ages",): continue
            s = "F" if r["Sex"] == "Female" else "M"
            D[s][AMAP[r["AgeGroups"]]] += float(r["v4_0"])

    def read_pop(path):
        P = {"F": defaultdict(float), "M": defaultdict(float)}
        with open(path) as f:
            for r in _csv.DictReader(f):
                s = "F" if r["GENDER_NAME"] == "Female" else "M"
                code = int(r["C_AGE_CODE"])          # 101=age0 ... 191=90+
                age = code - 101
                P[s][band_of_age(min(age, 90), 90)] += float(r["OBS_VALUE"])
        return P
    p24 = read_pop(f"{SCRATCH}/ew_pop_mid2024_nomis_NM2002.csv")
    p23 = read_pop(f"{SCRATCH}/ew_pop_mid2023_nomis.csv")
    P = {"F": {}, "M": {}}
    for s in ("F", "M"):
        for b in B.BANDS20:
            P[s][b] = max(2 * p24[s].get(b, 0) - p23[s].get(b, 0), 0.0)
    per_year = {2025: finish_bands(D, P, B.BANDS20)}

    hD = B.read_hmd_5x1(f"{B.DATA}/GBRTENW/STATS/Deaths_5x1.txt")
    val_rep = {"ew_2025_D_T": round(sum(per_year[2025][b]["D_T"] for b in B.BANDS20)),
               "ew_2025_pop_T": round(sum(P["F"][b] + P["M"][b] for b in B.BANDS20))}

    name = "England and Wales (ONS weekly age-sex deaths 2025 Registrations; exposure=extrapolated mid-2025)"
    wide = B.master_rows("GBRTENW_ONS", name, "ONS-weekly", 2025, per_year[2025], B.BANDS20)
    p90 = B.master_rows("GBRTENW_ONS", name, "ONS-weekly", 2025, per_year[2025], B.BANDS20)
    _append_txt_native("GBRTENW_ONS", 2025, per_year[2025], B.BANDS20,
                       "to 90+; 2025 ONS weekly Registrations, extrapolated exposure")
    return wide, p90, val_rep, set()

# ================================================================ KOREA
# 2025 deaths from KOSTAT provisional release (출생·사망통계 잠정, table "성·연령별
# 사망", unit=thousands, 10-yr bands, top 90+). Values transcribed from the HWPX
# (kor_births_deaths_2025_provisional.hwpx, Table 50, 2025p columns). Finer 5-yr x
# sex counts are not public until the Sept-2026 detailed release / KOSIS (API key).
# The 10-yr totals are split to 5-yr bands using 2024 HMD Korea death shares.
KOR_2025_10YR_K = {   # band(10yr): (Female_k, Male_k)  in thousands
    "0":      (0.3, 0.3),   "1-9":   (0.1, 0.2),  "10-19": (0.3, 0.4),
    "20-29":  (0.8, 1.4),   "30-39": (1.5, 2.9),  "40-49": (3.5, 6.8),
    "50-59":  (7.3, 17.6),  "60-69": (13.3, 35.3),"70-79": (24.8, 47.8),
    "80-89":  (69.6, 63.1), "90+":   (48.1, 18.0),
}
KOR_10TO5 = {"0": ["0"], "1-9": ["1-4", "5-9"], "10-19": ["10-14", "15-19"],
             "20-29": ["20-24", "25-29"], "30-39": ["30-34", "35-39"],
             "40-49": ["40-44", "45-49"], "50-59": ["50-54", "55-59"],
             "60-69": ["60-64", "65-69"], "70-79": ["70-74", "75-79"],
             "80-89": ["80-84", "85-89"], "90+": ["90+"]}

def build_kor():
    import csv as _csv
    hD = B.read_hmd_5x1(f"{B.DATA}/KOR/STATS/Deaths_5x1.txt")
    def share(subbands, sex_idx, y=2024):
        w = {b: hD.get((y, b), (0, 0, 0))[sex_idx] for b in subbands}
        tot = sum(w.values())
        return {b: (w[b] / tot if tot else 1.0 / len(subbands)) for b in subbands}
    D = {"F": defaultdict(float), "M": defaultdict(float)}
    for tenb, (fk, mk) in KOR_2025_10YR_K.items():
        subs = KOR_10TO5[tenb]
        for s, val_k, si in (("F", fk, 0), ("M", mk, 1)):
            sh = share(subs, si)
            for b in subs:
                D[s][b] += val_k * 1000.0 * sh[b]

    # population: MOIS resident registration, mean(Jan-2025, Jan-2026), national row
    def read_mois(path, ym):
        with open(path, encoding="euc-kr") as f:
            rows = list(_csv.reader(f))
        hdr = rows[0]; nat = rows[1]                 # 전국 (national) is first data row
        out = {"F": defaultdict(float), "M": defaultdict(float)}
        for i, col in enumerate(hdr):
            for klab, s in ((f"{ym}_남_", "M"), (f"{ym}_여_", "F")):
                if col.startswith(klab) and col.endswith("세") or col.endswith("세 이상"):
                    if not col.startswith(klab): continue
                    agestr = col[len(klab):]
                    if agestr == "100세 이상": age = 100
                    elif agestr.endswith("세"): age = int(agestr[:-1])
                    else: continue
                    v = float(nat[i].replace(",", "")) if nat[i].replace(",", "").strip().isdigit() else 0
                    out[s][band_of_age(age, 100)] += v
        return out
    p25 = read_mois(f"{SCRATCH}/kor_mois_pop_2025-01_single_age.csv", "2025년01월")
    p26 = read_mois(f"{SCRATCH}/kor_mois_pop_2026-01_single_age.csv", "2026년01월")
    def to90(bands100):  # collapse 90-94/95-99/100+ into 90+
        out = defaultdict(float)
        for b, v in bands100.items():
            out["90+" if B.age_start(b) >= 90 else b] += v
        return out
    P = {"F": {}, "M": {}}
    for s in ("F", "M"):
        m = to90({b: (p25[s].get(b, 0) + p26[s].get(b, 0)) / 2 for b in B.BANDS22})
        for b in B.BANDS20:
            P[s][b] = m.get(b, 0)

    per_year = {2025: finish_bands(D, P, B.BANDS20)}
    hmd24 = sum(v[2] for (y, b), v in hD.items() if y == 2024)
    hMx = B.read_hmd_5x1(f"{B.DATA}/KOR/STATS/Mx_5x1.txt")
    val_rep = {"kor_2025_D_T": round(sum(per_year[2025][b]["D_T"] for b in B.BANDS20)),
               "kostat_2025_total_raw": 363400, "hmd_2024_D_T": round(hmd24)}

    name = ("South Korea (KOSTAT 2025 provisional deaths, 10-yr bands->5-yr via 2024 HMD "
            "shares; exposure=MOIS resident reg mean Jan-2025/Jan-2026)")
    wide, p90 = emit("KOR_KOSTAT", name, "KOSTAT-MOIS", per_year, B.BANDS20,
                     "to 90+; 2025 PROVISIONAL, coarse (10-yr native, 0.1k resolution)")
    return wide, p90, val_rep, {"KOR_KOSTAT"}

# ================================================================ USA
def build_usa():
    """2025 deaths from CDC WONDER Provisional Mortality D176 (Five-Year Age Groups
       x Sex, provisional; retrieved via browser 2026-07-12 -> usa_wonder_deaths_2025.json),
       full detail to 100+, NO suppressed cells. Exposure = Census Vintage-2025 July-1
       resident population by single age x sex (usa_nc-est2025-agesex-res.csv)."""
    import csv as _csv
    dj = json.load(open(f"{SCRATCH}/usa_wonder_deaths_2025.json"))["deaths"]
    D = {"F": defaultdict(float), "M": defaultdict(float)}
    ns = {"F": dj["NotStated"][0], "M": dj["NotStated"][1]}
    for k, (f, m) in dj.items():
        if k == "NotStated": continue
        band = "0" if k == "<1" else k
        D["F"][band] += f; D["M"][band] += m
    for s, si in (("F", 0), ("M", 1)):          # distribute Age Not Stated pro rata
        tot = sum(D[s].values())
        if ns[s] and tot:
            for b in D[s]: D[s][b] *= (1 + ns[s] / tot)

    P = {"F": defaultdict(float), "M": defaultdict(float)}
    with open(f"{SCRATCH}/usa_nc-est2025-agesex-res.csv") as f:
        for r in _csv.DictReader(f):
            if r["SEX"] == "0" or r["AGE"] == "999": continue
            s = "M" if r["SEX"] == "1" else "F"
            age = int(r["AGE"])                  # 0..100 (100 = 100+)
            P[s][band_of_age(age, 100)] += float(r["POPESTIMATE2025"])

    per_year = {2025: finish_bands(D, P, B.BANDS22)}
    hD = B.read_hmd_5x1(f"{B.DATA}/USA/STATS/Deaths_5x1.txt")
    hMx = B.read_hmd_5x1(f"{B.DATA}/USA/STATS/Mx_5x1.txt")
    hmd24 = sum(v[2] for (y, b), v in hD.items() if y == 2024)
    val_rep = {"usa_2025_D_T": round(sum(per_year[2025][b]["D_T"] for b in B.BANDS22)),
               "wonder_total_raw": 3095562, "hmd_2024_D_T": round(hmd24),
               "pop_2025_T": round(sum(P["F"][b] + P["M"][b] for b in B.BANDS22)),
               "mx7579_ratio_vs_hmd24": round((per_year[2025]["75-79"]["D_T"] /
                     per_year[2025]["75-79"]["P_T"]) / hMx[(2024, "75-79")][2], 4)
                     if (2024, "75-79") in hMx else None}

    name = "USA (CDC WONDER Provisional Mortality D176, 2025 prov.; exposure=Census Vintage-2025 July-1)"
    wide, p90 = emit("USA_WONDER", name, "WONDER-Census", per_year, B.BANDS22,
                     "to 100+; 2025 PROVISIONAL (deaths occurring, as of Jul 2026)")
    return wide, p90, val_rep, {"USA_WONDER"}

BUILDERS = {"CHL": build_chl, "DEU": build_deu, "EW": build_ew, "KOR": build_kor, "USA": build_usa}

DROP2025 = {"DEU": "DEU_EUROSTAT", "EW": "GBRTENW_ONS"}

if __name__ == "__main__":
    wide, p90, drop = [], [], set()
    drop2025 = set(); reports = {}
    for arg in sys.argv[1:]:
        w, n, r, d = BUILDERS[arg]()
        wide += w; p90 += n; drop |= d; reports[arg] = r
        if arg in DROP2025: drop2025.add(DROP2025[arg])
    B.update_masters(wide, p90, drop_all_codes=drop, drop_2025_codes=drop2025)
    print(json.dumps(reports, indent=1, default=str))
