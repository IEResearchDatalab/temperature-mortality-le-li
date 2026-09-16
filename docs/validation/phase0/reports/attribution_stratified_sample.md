# Phase 0 — Attribution stratified sample

## Verdict
PASS

## Sample design
- Deterministic sample cells: 20
- First city in each sorted region, with all available age groups for those cities
- Same historical ERA5 / central-coefficient / fixed-death validation logic as the AT001C fixture
- The sample is used only for absolute-difference screening across regions and ages.

## Current diagnosis
- Reference total AN quantiles: min 2.230, median 66.127, max 388.608
- Total abs-diff quantiles: min 0, 25% 3.020e-14, median 7.105e-14, 75% 4.192e-13, 90% 7.077e-13, max 1.535e-12
- Total rel-diff quantiles: min 0, 25% 7.764e-16, median 1.456e-15, 75% 2.807e-15, 90% 3.231e-15, max 3.949e-15
- Mean absolute total difference on the sample: 2.903e-13
- Mean relative total difference on the sample: 1.705e-15
- Sign agreement on total ANs: 20 / 20

## Interpretation
- The code path is internally unit-consistent: raw daily sums are scaled to annual values by the year-specific day count.
- The earlier mismatch came from the averaging operator across years: equal-year means versus year-days weighted means.
- Once year-days weighting is applied, the sample matches the Masselot fixture to floating-point tolerance across all 20 cells.
- The result is a computational/reproductive correction, not a scientifically material discrepancy.
- The earlier `1e-6` screening threshold was a post hoc placeholder; the corrected sample comparison uses floating-point tolerance instead.

## Status
PASS
