# Baseline choice and pandemic-era excess mortality — reproducibility repository

Reproduces every figure and table in *Paper A* (Ioannidis & Levitt et al.), which shows that the choice
of **baseline model** — age resolution and pre-pandemic temporal trend — dominates estimates of
COVID-era excess mortality across **38 populations, 20 five-year age bands, 2020–2025**.

From a single vendored dataset, `./run_all.sh` regenerates all 6 manuscript figures, the 8-sheet table
workbook (Tables 1–3 and S1–S5), and the supplementary heterogeneity tables.

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
| **TTa** | two-step log-quadratic trend (moving 5-year slopes → slope-of-slopes), per band |
| **STTa** | TTa with the country trend shrunk toward the global trend (√deaths weight) |
| **STTa⁺** | STTa fitted with a {2019, 2024, 2025} return-slope anchor |

Everything is fit **per 5-year age band** and summed, which is the paper's central methodological point:
coarser age bands bias the excess. Pooled all-ages excess ranges from **1.16M (flat)** to
**2.4–2.6M (trended/anchored)** over 2020–2025 — a spread driven entirely by baseline choice.

## Pipeline

```
data/master_5x1_DPM_90plus.csv   (vendored analytic dataset — deaths/exposures by band·year·population)
        │
        ├─ build_hmd_calibration.py ──▶ output/hmd_calibration.csv
        ├─ stmf_moving_slopes.py ─────▶ output/stmf_moving_slopes.json
        │        └─ slope_of_slopes_ci.py ─▶ output/slope_of_slopes_CI.csv   (two-step trend + 95% CIs)
        ▼
   model_option1_periods.py   (THE engine)
        └─▶ output/{model_option1_periods, pooled_option1_by_period, option1_country_params}.csv
                 │
                 ├─ fig_*.py ─────────────▶ 6 manuscript figures  (output/*.png)
                 └─ build_*tables*.py ─────▶ Tables 1–3 + S1–S7    (docs/*.xlsx)
```

See `output/EXPECTED_OUTPUTS.md` for the artifact-by-artifact map.

## Layout

```
├── run_all.sh                 # one-command reproduction (master → all figures + tables)
├── requirements.txt
├── data/
│   ├── master_5x1_DPM_90plus.csv   # the analytic dataset (see data/README.md for schema)
│   ├── canonical_populations.csv   # the 38 study populations
│   ├── vuln_covariates_38_v3.csv   # vulnerability groups + economic covariates
│   ├── methods_comparison.json     # curated methods landscape (Table 1)
│   └── README.md                   # schema + provenance
├── code/
│   ├── excess_anchor_window.py     # master loader (deaths/exposures, HMD-basis calibration)
│   ├── stmf_moving_slopes.py       # moving-slope trend-of-trends
│   ├── slope_of_slopes_ci.py       # two-step trend parameters + 95% CIs
│   ├── model_option1_periods.py    # the excess-mortality engine
│   ├── fig_*.py                    # the 6 manuscript figures
│   ├── build_*.py                  # the table workbooks
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
Everything from the master onward regenerates with one command. The raw national/HMD source files are
**not** included: the Human Mortality Database prohibits redistribution, and several national sources
are login-gated (e.g. CDC WONDER). The builders that assembled the master from raw sources are provided
under `code/fetch/` for provenance, and every source is documented in `docs/data_sources.md`.

## License

Analysis **code** is MIT-licensed. The **data** is governed by the original providers' terms — see
`LICENSE` and `docs/data_sources.md`. In particular, redistribution of HMD data is not permitted.
