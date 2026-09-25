# `code/fetch/` — raw → master builders (provenance)

These scripts document **how `data/master_5x1_DPM_90plus.csv` was assembled** from raw national and
HMD sources. They are included for transparency and audit.

> **They are not part of the one-command reproduction.** The reproducible pipeline starts from the
> vendored master CSV (`../../data/master_5x1_DPM_90plus.csv`) and is driven by `../../run_all.sh`.
> The builders here read raw inputs that are **not** shipped in this repo (HMD is not
> redistributable; several national files are login-gated). To re-run them, obtain the raw sources
> per [`../../docs/data_sources.md`](../../docs/data_sources.md) and put them in `data/raw/` (git-ignored),
> or point the `PAPERA_RAW` environment variable at the folder that holds them.

## What each builder adds

| Script | Adds to master | Raw inputs (see docs/data_sources.md) |
|---|---|---|
| `build_2025_eurostat_wk.py` | 2025 for existing `*_EUROSTAT` series + 13 new full 2020–2025 `_EUROSTAT` series (BEL, CHE, EST, FIN, HRV, LTU, LUX, LVA, NOR, ISL, PRT, SVK, SWE). Defines the shared row/band helpers the other two `build_2025_*` scripts import. | Eurostat `demo_r_mwk_05` (weekly deaths), `demo_pjan`, `demo_magec` |
| `build_2025_others.py` | 2025 for `CHL_DEIS`, `DEU_EUROSTAT`, `GBRTENW_ONS`, `KOR_KOSTAT`, `USA_WONDER` | Chile DEIS microdata + INE; Destatis GENESIS 12613-0003; UK ONS weekly + Nomis mid-year pop; KOSTAT/MOIS; CDC WONDER + US Census |
| `build_2025_asia.py` | `TWN_MOI` (2025), `HKG_CSD` (2020–2025), `JPN_ESTAT` (2025) | Taiwan MOI; Japan e-Stat; Hong Kong C&SD |
| `assemble_new_locations.py` | AUS 2022–24, CAN 2024, ISL 2024–25, ISR 2017–23, NZL_NP 2022–25 | ABS Data API, Stats NZ, StatCan, Israel CBS (some values inlined) |
| `build_israel_2024_2025.py` | ISR 2024 & 2025 | Israel CBS 5-year-band CSV |
| `build_australia_2025.py` | AUS 2025 | HMD STMF pooled file |
| `build_canada_2025.py` | CAN 2025 (through ISO week 45, grossed up) | HMD STMF pooled file |
| `update_canada_2025_complete.py` | Replaces CAN 2025 with full-year data | Statistics Canada table 13-10-0768 |
| `build_ireland_2025.py` | IRL 2025 | CSO Ireland VSA07 |
| `build_uk_scotland_ni_2023_2025.py` | `SCO_NRS`, `NIR_NISRA`, `GBR_NP_BUILT` 2023–2025 (Table S8 only) | NRS deaths time series and mid-year estimates; NISRA deaths and mid-year estimates; HMD STMF pooled file |
| `refresh_hmd_exposures.py` | Replaces the `FRATNP` and `USA` HMD rows with the 27 August 2026 release. Writes `master_5x1_DPM_90plus_hmd20260827.csv`, which was then installed as the master. | HMD "all countries" bundle, unzipped |

## Order (if rebuilding from raw)

```
(HMD + Eurostat + ONS backbone, pre-existing base master)
  → build_2025_eurostat_wk.py
  → build_2025_others.py          # imports build_2025_eurostat_wk
  → build_2025_asia.py
  → assemble_new_locations.py
  → build_israel_2024_2025.py
  → build_australia_2025.py
  → build_canada_2025.py
  → update_canada_2025_complete.py
  → build_ireland_2025.py
  → build_uk_scotland_ni_2023_2025.py
  → refresh_hmd_exposures.py      # then copy its output over the master
```

Each builder is idempotent (it strips its own prior rows before re-appending) and writes a dated
backup of the master before modifying it.

> **Note.** The from-scratch builder that first assembled the HMD/Eurostat/ONS backbone master
> (2020–2024, all populations) predates this collection and is not included; the backbone is embodied
> in the vendored master. The scripts here cover the incremental extensions (chiefly the 2025 rows and
> the added non-HMD populations).
