# Project Status -- Last Updated: 2026-07-25

This file is the single source of truth for current project state. Use it before reading archived session notes or historical debugging reports.

## ✅ COMPLETED & VERIFIED
- [x] Basis function fix completed: the project now uses `bs` with degree `2` instead of the earlier spline setup. Evidence: `R_pipeline/01_initialize.R`, `R_pipeline/03_attribution.R`, `R_pipeline/03_attribution_baseline.R`.
- [x] 2023 attributable-number formula implemented for baseline replication. Evidence: `R_pipeline/03_attribution_baseline.R` and exact match in `results/masselot_validation_city/CITY_VALIDATION_SUMMARY.md`.
- [x] 4-range temperature disaggregation implemented: `ExtrCold`, `ModCold`, `ModHeat`, `ExtrHeat`. Evidence: `R_pipeline/03_attribution.R` and `R_pipeline/03_attribution_baseline.R`.
- [x] ERA5 baseline replication completed for all 854 cities. Evidence: `temp_results_baseline/` contains 854 `.rds` files.
- [x] City-level baseline validation passes exactly. Evidence: `results/masselot_validation_city/CITY_VALIDATION_SUMMARY.md` reports `0.00%` error and `1.0000` correlation for cold and heat across 4,270 city-age combinations.
- [x] Script 06 local code defects fixed. Evidence: `R_pipeline/06_le_li_extraction.R` no longer has the earlier `data.table` scoping bug or simulation inflation logic.
- [x] Script 07 local code defects fixed at the access-pattern level. Evidence: `R_pipeline/07_le_li_decomposition.R` uses `pclm(... )$fitted` instead of the earlier invalid slot access.

## 🚧 IN PROGRESS
- [ ] No active coding task is in progress. Baseline replication is complete.

## ❌ BLOCKED / NOT STARTED
- [ ] Script 04 and Script 05 have not been re-validated against the new baseline workflow.
  - **Blocker:** `R_pipeline/04_aggregate_results.R` and `R_pipeline/05_generate_tables.R` still assume `temp_results/`, while the verified baseline outputs live in `temp_results_baseline/`.
- [ ] GCM projection rerun has not been regenerated with a current `temp_results/` directory.
  - **Blocker:** the live worktree has no current `temp_results/`; only backups and the validated `temp_results_baseline/` exist.
- [ ] Script 05 table generation remains unverified in the current workflow.
  - **Blocker:** no current table outputs in `results/tables/`, and prior runs are not trustworthy after the baseline-method changes.
- [ ] Script 07 PCLM decomposition remains unverified end-to-end.
  - **Blocker:** upstream `results/le_li_input/le_li_input_ans.csv` has not been regenerated from a current workflow.
- [ ] Total all-cause mortality data has not been acquired.
  - **Blocker:** required for counterfactual life-table construction and LE/LI analysis.
- [ ] Demographic methods document has not been drafted.
  - **Blocker:** PI requested methodology design approval before further LE/LI coding.

## 📊 CURRENT VALIDATION STATE

| Component | Status | Error | Evidence |
|-----------|--------|-------|----------|
| Baseline attribution (ERA5, `03_attribution_baseline.R`) | ✅ Complete | `0.00%` cold / `0.00%` heat | `results/masselot_validation_city/CITY_VALIDATION_SUMMARY.md` |
| 4-range additivity | ✅ Complete | Exact to floating-point tolerance | `R_pipeline/09_validate_city_level.R` and validation outputs |
| GCM projection path (`03_attribution.R`) | ⚠️ Code updated, not rerun | Not currently validated from fresh outputs | No current `temp_results/` directory |
| Aggregation (`04_aggregate_results.R`) | ⚠️ Historical output present, not current | Unknown | `results/attribution_aggregated/final_attribution_results.csv` should be treated as stale until rerun |
| Tables (`05_generate_tables.R`) | ❌ Not current | No current output | `results/tables/` is empty |
| LE/LI input (`06_le_li_extraction.R`) | ❌ Not current | No current output | `results/le_li_input/` is empty |
| LE/LI decomposition (`07_le_li_decomposition.R`) | ❌ Not current | No current output | `results/le_li_decomposition/` is empty |

## 🎯 NEXT SESSION PRIORITY
Test `R_pipeline/04_aggregate_results.R` and `R_pipeline/05_generate_tables.R` against `temp_results_baseline/` with the lightweight harness before attempting any full GCM rerun.

## 📁 KEY FILE LOCATIONS & STATUS
- `temp_results_baseline/` -- ✅ Current. Verified ERA5 baseline results for all 854 cities.
- `temp_results/` -- ❌ Missing in current worktree. Must be regenerated for the GCM projection path.
- `temp_results_backup_before_bs_fix/` -- ⚠️ Historical backup. Do not use as current output.
- `temp_results_backup_before_fix/` -- ⚠️ Historical backup. Do not use as current output.
- `temp_results_backup_wrong_formula_2024-07-24/` -- ⚠️ Historical backup created before the final formula fix.
- `results/masselot_validation_city/` -- ✅ Current baseline validation outputs.
- `results/attribution_aggregated/final_attribution_results.csv` -- ⚠️ Present but not confirmed current after the workflow split.
- `results/tables/` -- ❌ Empty.
- `results/le_li_input/` -- ❌ Empty.
- `results/le_li_decomposition/` -- ❌ Empty.
- `references/archive_sessions/` -- ✅ Archive location for superseded diaries, TODOs, and PI correction notes.

## 🚨 CRITICAL FACTS TO REMEMBER
- Baseline replication uses ERA5 observations, not GCM historical series.
- Baseline validation target is exact reproduction of `references/2025-masselot-zenodo/results/cityage.csv`.
- GCM projections and ERA5 baseline are separate tracks and should not be validated against the same expectation.
- The 4-range split is implemented and validated; do not re-open that problem unless a new regression appears.
- Read this file first. Do not start a new session by re-reading archived diaries unless this file is insufficient.
- Treat historical reports in `results/` and archived notes in `references/archive_sessions/` as context only, not as live status.

## 🗂️ DOCUMENTATION RULE
- Current status: `PROJECT_STATUS.md`
- Next-task planning: `SESSION_GOAL_YYYY-MM-DD.md`
- Reusable workflow guidance: `.copilot-instructions`
- Fast checks: `test_pipeline_component.R`
- Historical debugging notes: `references/archive_sessions/`