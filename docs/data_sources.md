# Data sources

The analysis covers 38 populations, 20 five-year age bands, calendar years ~2003–2025. All death and
population counts trace to official statistical offices. This page documents every source and how to
obtain it. The processed result is `data/master_5x1_DPM_90plus.csv`, which is vendored so the analysis
runs without re-downloading anything.

> **Redistribution.** The Human Mortality Database prohibits redistribution of its files, so no raw
> HMD data is included here. Other providers have their own terms. The single derived, aggregated
> `master_5x1_DPM_90plus.csv` is provided for reproducibility; reuse must respect the upstream
> providers' terms (in particular the HMD user agreement).

## Backbone (all populations, through ~2024)

| Source | Coverage | Access |
|---|---|---|
| **Human Mortality Database** (mortality.org) | Deaths and exposures by 5-year age band, historical → ~2020–2024 depending on country, for all HMD members | **Free after registration; login-gated; not redistributable.** Bulk `Deaths_5x1` / `Exposures_5x1` per country, or the zipped country bundle. |
| **HMD Short-Term Mortality Fluctuations (STMF)** | Weekly deaths by broad age band | Free download from the STMF page on mortality.org (pooled `stmf.csv`). Used to extend a few 2025 series. |
| **Eurostat** `demo_r_mwk_05` (weekly deaths), `demo_pjan` (population on 1 Jan), `demo_magec` | European populations, recent years incl. 2025 | Open, no login. Eurostat API / bulk download. |
| **UK ONS** + **Nomis** | England & Wales weekly deaths by age/sex; mid-year population (`NM_161_1` / `NM_2002_1`) | Open. |

## 2025 and additional-population sources

| Population(s) | Source | Access / notes |
|---|---|---|
| 13 European (`*_EUROSTAT`) + 2025 rows | Eurostat `demo_r_mwk_05` / `demo_pjan` / `demo_magec` | Open API. |
| USA (`USA_WONDER`) | **CDC WONDER** provisional deaths + US Census population | **Gated:** WONDER sits behind an Akamai TLS check and a click-through data-use agreement; retrieved via a real browser and hand-saved to `data/USA_WONDER/raw/`. (Ben Marten's independent USA 2025 extract matched WONDER exactly, as a cross-check.) |
| Germany 2025 (`DEU_EUROSTAT`) | **Destatis GENESIS** table `12613-0003` | Open portal; JSON export. |
| Chile (`CHL_DEIS`) | **DEIS** deaths microdata + **INE** population projections | Open portal. |
| South Korea (`KOR_KOSTAT`) | **KOSTAT / MOIS** | Partly transcribed from HWPX; finer detail needs a KOSIS API key or the later detailed release. |
| Japan (`JPN_ESTAT`) | **e-Stat** (statInfId 000040458974) | Open portal (cp932 CSV). |
| Taiwan (`TWN_MOI`) | **Ministry of the Interior** (data.gov.tw 146625) | Open; raw copied to `data/TWN_MOI/raw/`. |
| Hong Kong (`HKG_CSD`) | **Census & Statistics Department** | Open tables. |
| Australia (`AUS`) | **ABS Data API** + HMD STMF | Open API. |
| Canada (`CAN`) | **Statistics Canada** table `13-10-0768` (provisional weekly deaths by age/sex) | Open. Full-year 2025 replaces the earlier week-1–45 gross-up. |
| New Zealand (`NZL_NP`) | **Stats NZ** | Open; some values inlined. |
| Ireland (`IRL`) | **CSO** `VSA07` | Open; JSON. |
| Israel (`ISR`) | **CBS** deaths & population by 5-year band | Open CSV. |

## Covariates (static, curated)

| File | Contents | Source |
|---|---|---|
| `data/vuln_covariates_38_v3.csv` | vulnerability group + GDP per capita 2021, Gini, poverty rate | World Bank Open Data (`data.worldbank.org`), OECD |
| `data/methods_comparison.json` | prior baseline-method studies for Table 1 | curated from the cited literature |

## Reproducing the master from raw

See [`../code/fetch/README.md`](../code/fetch/README.md). In brief: obtain the HMD backbone and the
national/Eurostat sources above, place them where each `code/fetch/build_*` script expects them (repoint
the input-path constants), then run the builders in the documented order. Each is idempotent and backs
up the master before writing. Because HMD is login-gated and not redistributable, and CDC WONDER is
gated behind an interactive agreement, a fully unattended raw→master rebuild is not possible from this
repository alone — which is why the processed master is vendored.
