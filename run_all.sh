#!/usr/bin/env bash
# Reproduce every Paper A figure and table from the vendored master dataset.
#
#   ./run_all.sh
#
# Override the interpreter if needed:  PYTHON=/opt/homebrew/bin/python3 ./run_all.sh
# Requires the packages in requirements.txt.  Runs in dependency order; stops on first error.
set -euo pipefail
cd "$(dirname "$0")"
PYTHON="${PYTHON:-python3}"
run() { echo; echo "==> $*"; "$PYTHON" "$@"; }

echo "Paper A — reproducing all outputs with: $("$PYTHON" --version 2>&1)"

# ---- 1. upstream intermediates (all derived from data/master_5x1_DPM_90plus.csv) ----
run code/build_hmd_calibration.py      # -> output/hmd_calibration.csv  (national→HMD basis ratios)
run code/stmf_moving_slopes.py         # -> output/stmf_moving_slopes.json (moving-slope trend-of-trends)
run code/slope_of_slopes_ci.py         # -> output/slope_of_slopes_CI.csv (two-step trend params + 95% CIs)

# ---- 2. the engine: excess by model × period × age group × population ----
run code/model_option1_periods.py      # -> output/{model_option1_periods,pooled_option1_by_period,option1_country_params}.csv

# ---- 3. a derived table input ----
run code/build_option1_age_resolution.py   # -> output/age_resolution_option1.csv  (feeds Table S1)

# ---- 4. figures ----
run code/fig_slopeofslopes_montage.py           # Figure 1
run code/fig_option1_3panel.py                  # Figure 2
run code/fig_option1_wealth_9panel_indperiod.py # Figure 3
run code/fig_option1_logquad_montage.py         # Figure S1
run code/fig_perband_trend_10countries.py       # Figure S2  (per-age-band trends, 10 countries)
run code/fig_option1_vulnerability_deaths.py    # Figure S3  (vulnerability, excess deaths)

# ---- 5. tables ----
run code/build_slope_heterogeneity_table.py    # between-country Q / I^2 heterogeneity (SI text)
run code/build_table_S6_agespecific.py         # Table S6 (age-band-specific vs all-age TTa excess)
run code/build_manuscript_tables.py            # ALL manuscript tables, manuscript numbering
                                               #   -> docs/PaperA_manuscript_tables.xlsx
run code/build_paperA_all_tables.py            # working workbook (own sheet naming), kept for continuity

# Superseded, retained for provenance (the 3-country precursor of Figure S2); not part of the
# current manuscript. Uncomment to regenerate:
# run code/fig_perband_trend_heterogeneity.py

echo
echo "DONE.  See output/EXPECTED_OUTPUTS.md for the full list of figures/tables produced."
