# `data/` — inputs to the analysis

This folder holds the **processed analytic inputs** the pipeline reads. It does **not** contain the
raw national/HMD source files (those are large and, for HMD, not redistributable — see
[`../docs/data_sources.md`](../docs/data_sources.md)).

## Files

| File | Role | Regenerable? |
|---|---|---|
| `master_5x1_DPM_90plus.csv` | **The master analytic dataset.** Deaths, exposures and central death rates by 5-year age band, calendar year and population, for the 38 study populations plus their national "child" source series. Everything downstream is computed from this file. | Built incrementally from raw sources by the scripts in `../code/fetch/` (see below). Vendored here so the analysis runs without the raw data. |
| `canonical_populations.csv` | The 38 study populations (analysis code → display name), in the paper's order. Defines the analysis set. | Static. |
| `vuln_covariates_38_v3.csv` | Per-population vulnerability group (more/less vulnerable) and economic covariates (GDP per capita 2021, Gini, poverty rate). Used by Table 3, Table S7 and Figures 3 and S3. | Static (curated from World Bank / OECD; see `docs/data_sources.md`). |

## Master schema

CSV header (identical for every row):

```
location,location_name,source,year,age,age_start,D_F,D_M,D_T,P_F,P_M,P_T,M_F,M_M,M_T
```

| Column | Meaning |
|---|---|
| `location` | population / source code (e.g. `USA`, `DEUTNP`, `JPN`; national child series carry a suffix such as `USA_WONDER`, `DEU_EUROSTAT`, `CHL_DEIS`) |
| `location_name` | human-readable name |
| `source` | provenance tag for that row: `HMD`; `HMD_20260827` (France and United States rows refreshed from that HMD release); `Eurostat`; a national office such as `NRS` or `NISRA`; `ASSEMBLED` (the whole UK built from its parts); or a build tag such as `STATCAN768_2025complete` |
| `year` | calendar year |
| `age` | 5-year band label: `0`, `1-4`, `5-9`, …, `85-89`, `90+` |
| `age_start` | integer lower edge of the band |
| `D_F,D_M,D_T` | deaths — female, male, total |
| `P_F,P_M,P_T` | mid-year exposure / population — female, male, total |
| `M_F,M_M,M_T` | central death rate D/P per person, six decimals — female, male, total. Not read by the analysis, which computes every rate from D and P. |

The analysis uses **Total** (`_T`) throughout. Many recently-appended national rows populate only
the `_T` columns (F/M blank); the loader (`../code/excess_anchor_window.py`) keeps whichever sexes are
present and analyses Total.

The `_90plus` naming means all ages ≥90 are collapsed into a single `90+` band, giving a uniform
20-band grid across all populations.

## How the master is built (provenance)

The master is a **stateful, incrementally-appended** artifact:

1. An HMD + Eurostat + ONS backbone (1900s–2024) forms the base.
2. The scripts in [`../code/fetch/`](../code/fetch) then append recent years and additional
   populations (notably the 2025 rows), each writing a dated backup and re-deriving idempotently.
3. `build_uk_scotland_ni_2023_2025.py` adds Scotland (`SCO_NRS`), Northern Ireland (`NIR_NISRA`)
   and the whole United Kingdom (`GBR_NP_BUILT`) for 2023–2025. They are used only by Table S8.
4. `refresh_hmd_exposures.py` replaced the France (`FRATNP`) and United States (`USA`) HMD rows
   with the release of 27 August 2026, which revised their population estimates.

Those scripts are provided **for provenance and audit**, not as a one-command rebuild: several read
raw inputs from national portals, some login-gated (e.g. HMD and CDC WONDER). Each looks for its raw
inputs in `data/raw/` (git-ignored), or in the folder named by the `PAPERA_RAW` environment variable. See [`../code/fetch/README.md`](../code/fetch/README.md)
and [`../docs/data_sources.md`](../docs/data_sources.md) for exactly what each one consumes and how to
obtain it. The reproducible boundary of this repository is the vendored `master_5x1_DPM_90plus.csv`;
everything from the master onward regenerates with a single `./run_all.sh`.
