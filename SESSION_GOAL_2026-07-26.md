# Session Goal -- 2026-07-26

## PRIMARY OBJECTIVE
✅ Test Script 04 aggregation and Script 05 table logic against the validated baseline outputs in `temp_results_baseline/`.

## PRE-RUN VALIDATION
- [ ] Verify `temp_results_baseline/` contains 854 `.rds` files.
- [ ] Run `Rscript test_pipeline_component.R script04`.
- [ ] Run `Rscript test_pipeline_component.R script05`.
- [ ] Only attempt a full script run if the sample checks pass.

## SUCCESS CRITERIA
- [ ] Script 04 aggregation logic works on baseline-shaped inputs without structural errors.
- [ ] A baseline aggregation plan is defined clearly enough to produce `final_attribution_results_baseline.csv` or an equivalent validated output.
- [ ] Script 05 completes at least one sample batch without crashing.
- [ ] Any Script 05 failure is documented with the exact error and whether memory pressure is implicated.

## IF BLOCKED
- [ ] If Script 04 fails, compare its assumptions against `temp_results_baseline/` schema before editing downstream files.
- [ ] If Script 05 fails, retry only with a smaller sample or sequential processing; do not jump to a full run.
- [ ] Update `PROJECT_STATUS.md` before ending the session.

## FILES TO MODIFY
- [ ] `R_pipeline/04_aggregate_results.R` only if baseline input assumptions are wrong.
- [ ] `R_pipeline/05_generate_tables.R` only if the sample test exposes a concrete defect.

## PARKING LOT
- [ ] Script 06 and Script 07 remain out of scope until Script 05 is behaving.
- [ ] Future GCM projection reruns remain out of scope for this session.
- [ ] Total mortality data acquisition remains a separate blocker.