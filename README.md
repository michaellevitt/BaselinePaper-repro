# Baseline choice and pandemic-era excess mortality — reproducibility repository

Reproduces every figure and table in *Paper A* (Ioannidis & Levitt et al.), which shows that the choice
of **baseline model** — age resolution and pre-pandemic temporal trend — dominates estimates of
COVID-era excess mortality across **38 populations, 20 five-year age bands, 2020–2025**.

From a single vendored dataset, `./run_all.sh` regenerates all 6 manuscript figures and a workbook
with every computed table (Tables 1–3 and S1–S8, numbered as in the paper).

## Quick start

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
./run_all.sh
```

Outputs land in `output/` (figures, intermediate CSVs) and `docs/` (table workbooks). They are
git-ignored because they are fully reproducible; `output/EXPECTED_OUTPUTS.md` lists exactly what
should appear and gives the headline numbers to check against.

On macOS with a Homebrew Python, you may need `PYTHON=/opt/homebrew/bin/python3 ./run_all.sh`.

## What the analysis does

For each population it builds a family of pre-pandemic **baseline** mortality models and projects each
into 2020–2025 to get expected deaths; **excess = observed − expected**, reported as counts and as a
P-score (100·(O−E)/E). The baselines differ only in how they treat the pre-pandemic trend:

| Model | Baseline |
|---|---|
| **Fa** | flat (2017–2019 mean rate per age band) |
| **TTa** | two-step log-quadratic trend (moving 5-year slopes → slope-of-slopes) |
| **STTa** | TTa with the trend shrunk toward the 38-population trend (√deaths weight) |
| **STTa⁺** | STTa fitted with a {2019, 2024, 2025} return-slope anchor |

Each trend model is run two ways: **age-band trend modeling**, a separate trend fitted within each
five-year age band (the paper's headline numbers), and **country-level trend modeling**, one trend
per country applied to every band. Table 2 gives both.

Everything is fit **per 5-year age band** and summed, which is the paper's central methodological point:
coarser age bands bias the excess. Pooled all-ages excess over 2020–2025 ranges from
**0.92M (flat)** to **2.14–2.49M (trend models)**, a spread driven by baseline choice.

## Pipeline

```
data/master_5x1_DPM_90plus.csv   (vendored analytic dataset: deaths/exposures by band, year, population)
   │
   ├─ build_hmd_calibration.py, stmf_moving_slopes.py, slope_of_slopes_ci.py   (trend parameters)
   │
   ├─ 1. model_option1_periods.py              country-level engine
   │        └─ uk_vs_ew_full_family.py         Table S8 (validates against step 1)
   ├─ 2. model_agespecific_trends_v1.py        per-band trends; reads step 1 as its country-level arm
   ├─ 3. model_option1_periods_v2_agespecific.py   age-band engine: the headline numbers,
   │                                             overwrites step 1's output files
   ├─ model_leecarter_periods_v1.py            Lee-Carter contrast
   │
   ├─ fig_*.py ────────────▶ 6 manuscript figures      (output/*.png)
   └─ build_*tables*.py ───▶ Tables 1–3 and S1–S8       (docs/*.xlsx)
```

**The order of steps 1–3 matters.** The two engines write the same filenames, so the age-band engine
must run last, and the country-level engine must run first, because step 2 reads its output. Run out of
order, every script still succeeds but the country-level comparison in Table 2 is wrong. `run_all.sh`
explains this and runs them correctly; see `output/EXPECTED_OUTPUTS.md` for the artifact map.

## Layout

```
├── run_all.sh                 # one-command reproduction (master → all figures + tables)
├── requirements.txt
├── data/
│   ├── master_5x1_DPM_90plus.csv   # the analytic dataset (see data/README.md for schema)
│   ├── canonical_populations.csv   # the 38 study populations
│   ├── vuln_covariates_38_v3.csv   # vulnerability groups + economic covariates
│   └── README.md                   # schema + provenance
├── code/
│   ├── excess_anchor_window.py     # master loader (deaths/exposures, HMD-basis calibration)
│   ├── stmf_moving_slopes.py       # moving-slope trend-of-trends
│   ├── slope_of_slopes_ci.py       # two-step trend parameters + 95% CIs
│   ├── model_option1_periods.py    # country-level engine (runs first; see Pipeline)
│   ├── model_agespecific_trends_v1.py            # per-band trends
│   ├── model_option1_periods_v2_agespecific.py   # age-band engine (headline numbers)
│   ├── model_leecarter_periods_v1.py             # Lee-Carter contrast
│   ├── uk_vs_ew_full_family.py     # Table S8
│   ├── fig_*.py                    # the 6 manuscript figures
│   ├── build_*.py                  # the table workbooks (build_manuscript_tables.py: Tables 1–3, S1–S8)
│   └── fetch/                      # raw → master builders (provenance; see fetch/README.md)
├── docs/
│   ├── data_sources.md             # every raw source + how to obtain the gated ones
│   └── *.xlsx                      # generated table workbooks (git-ignored)
└── output/
    ├── EXPECTED_OUTPUTS.md          # what run_all.sh should produce (+ headline numbers)
    └── *.png / *.csv                # generated figures + intermediates (git-ignored)
```

## Data & reproducibility boundary

The reproducible boundary is the **processed master dataset** (`data/master_5x1_DPM_90plus.csv`).
Its France and United States rows were refreshed from the Human Mortality Database release of
27 August 2026, which revised population estimates (most at ages 80+, and nearly every age band in
the latest years); `code/fetch/refresh_hmd_exposures.py` documents that step.
Everything from the master onward regenerates with one command. The raw national/HMD source files are
**not** included: the Human Mortality Database prohibits redistribution, and several national sources
are login-gated (e.g. CDC WONDER). The builders that assembled the master from raw sources are provided
under `code/fetch/` for provenance, and every source is documented in `docs/data_sources.md`.

## License

Analysis **code** is MIT-licensed. The **data** is governed by the original providers' terms — see
`LICENSE` and `docs/data_sources.md`. In particular, redistribution of HMD data is not permitted.
