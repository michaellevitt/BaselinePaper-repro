#!/usr/bin/env bash
# Reproduce every Paper A figure and table from the vendored master dataset.
#
#   ./run_all.sh
#
# Override the interpreter if needed:  PYTHON=/opt/homebrew/bin/python3 ./run_all.sh
# Requires the packages in requirements.txt.  Runs in dependency order; stops on first error.
#
# ---------------------------------------------------------------------------------------
# ORDER MATTERS, and not in the obvious way.
#
# The two engines write the SAME output filenames.  model_option1_periods.py fits one trend
# per country; model_option1_periods_v2_agespecific.py fits a separate trend within each age
# band.  The paper reports the age-band results, so the v2 engine must run LAST.
#
# But the country-level engine must still run, and must run BEFORE
# model_agespecific_trends_v1.py, because that script reads pooled_option1_by_period.csv to
# get the country-level arm of the age-band-versus-country-level comparison in Table 2.  Run
# it afterwards, or not at all, and that comparison silently reports whatever happened to be
# in the file: either stale numbers, or the age-band results compared against themselves, in
# which case every difference comes out as exactly zero.  Both have happened.
#
# So the sequence is: country-level engine -> age-band trends -> age-band engine.
# Table S8 (uk_vs_ew_full_family.py) uses the country-level model family, and validates itself
# against the country-level engine's output, so it also runs between the first two.
#
# fig_slopeofslopes_montage.py, fig_option1_logquad_montage.py and
# fig_perband_trend_10countries.py are the superseded forms of Figures 1, S1 and S2 and are
# not run.
# ---------------------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")"
PYTHON="${PYTHON:-python3}"
run() { echo; echo "==> $*"; "$PYTHON" "$@"; }

echo "Paper A — reproducing all outputs with: $("$PYTHON" --version 2>&1)"

# ---- 1. upstream intermediates (all derived from data/master_5x1_DPM_90plus.csv) ----
run code/build_hmd_calibration.py      # -> output/hmd_calibration.csv  (national→HMD basis ratios)
run code/stmf_moving_slopes.py         # -> output/stmf_moving_slopes.json (moving-slope trend-of-trends)
run code/slope_of_slopes_ci.py         # -> output/slope_of_slopes_CI.csv (two-step trend params + 95% CIs)

# ---- 2. country-level engine: writes the country-level arm that step 3 reads ----
# Its results are overwritten in step 4; what survives is the copy step 3 takes of them.
run code/model_option1_periods.py      # -> output/pooled_option1_by_period.csv (country-level)
run code/uk_vs_ew_full_family.py       # Table S8; its gate needs THIS engine's England and Wales output

# ---- 3. per-age-band trends, and the comparison of the two modelings ----
# Reads the country-level pooled file from step 2 as its "common" arm.
run code/model_agespecific_trends_v1.py   # -> output/{agespecific_band_params_v1,agespecific_vs_common_v1}.csv

# ---- 4. the headline engine: excess by model × period × age group × population ----
# Fits a separate trend within each age band and overwrites step 2's files.  Must run LAST
# of the three.
run code/model_option1_periods_v2_agespecific.py   # -> output/{model_option1_periods,pooled_option1_by_period,option1_country_params}.csv

# ---- 5. the Lee-Carter contrast reported in Results ----
run code/model_leecarter_periods_v1.py    # -> output/{model_leecarter_periods_v1,pooled_leecarter_by_period_v1}.csv

# ---- 6. a derived table input ----
run code/build_option1_age_resolution.py  # -> output/age_resolution_option1.csv  (feeds Table S5)

# ---- 7. figures ----
run code/fig_slopeofslopes_montage_v2.py        # Figure 1   (adds the age-band trend line)
run code/fig_option1_3panel.py                  # Figure 2
run code/fig_option1_wealth_9panel_indperiod.py # Figure 3
run code/fig_option1_logquad_montage_v2.py      # Figure S1  (adds the age-band baseline)
run code/fig_perband_trend_38countries.py       # Figure S2  (per-age-band trends, all 38)
run code/fig_option1_vulnerability_deaths.py    # Figure S3  (vulnerability, excess deaths)

# ---- 8. tables ----
run code/build_slope_heterogeneity_table.py    # between-country Q / I^2 heterogeneity (SI text)
run code/build_table_S6_agespecific.py         # Table S6 (age-band vs country-level TTa excess)
run code/build_manuscript_tables.py            # Tables 1-3 and S1-S8, numbered as in the paper
                                               #   -> docs/PaperA_manuscript_tables.xlsx
run code/build_paperA_all_tables.py            # working workbook (own sheet naming), kept for continuity

# Not run, and not part of the current manuscript:
#   code/backtest_leecarter_2014_2019_v1.py     back-test, cut from the paper on 21 Sep 2026
#   code/fig_perband_trend_heterogeneity.py     3-country precursor of Figure S2
#   code/refresh_hmd_exposures_v1.py            one-off correction of the vendored master from
#                                               an HMD release; reads outside the repository, so
#                                               it sits outside the reproducible boundary (see
#                                               docs/data_sources.md)

echo
echo "DONE.  See output/EXPECTED_OUTPUTS.md for the full list of figures/tables produced."
