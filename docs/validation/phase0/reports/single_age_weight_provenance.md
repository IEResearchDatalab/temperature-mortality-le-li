# Phase 0 — Single-age weight provenance

## Verdict
BLOCKED

## File traced
- `results/le_li_input/lloyd_fig_s1_baseline_city_year_single_age.csv`
- Producer: `R_pipeline/11_build_lloyd_fig_s1_baseline.R`
- Input source: `temp_results_baseline/*.rds` filtered to `gcm == "ERA5"` and `sim == 0`

## What the producer does
1. Reads 65+ grouped AN baseline city files from `temp_results_baseline/*.rds`.
2. Disaggregates baseline population and baseline deaths with `pclm_disaggregate_nonnegative()`.
3. Disaggregates grouped AN by range with `pclm_disaggregate_signed()`.
4. Writes a single-age table containing `pop`, `death_baseline`, and signed AN ranges.

## Provenance conclusion
- The file is **fixed-baseline** and **ERA5 historical**; it is not SSP-specific, and it is not a direct authoritative future demographic-weight source. It is year-indexed output built from fixed-baseline historical inputs.
- Its construction uses PCLM on grouped historical inputs, so it is best treated as a validated baseline shape/reference file.
- That evidence supports a future-method choice in the **B** family: fit PCLM only to nonnegative projected all-cause deaths, using the validated baseline single-age profile as the reference shape or offset.

## Status
BLOCKED
