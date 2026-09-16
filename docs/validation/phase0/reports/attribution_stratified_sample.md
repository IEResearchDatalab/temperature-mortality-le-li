# Phase 0 — Attribution stratified sample

## Verdict
FAIL

## Sample design
- Deterministic sample cells: 20
- First city in each sorted region, with all available age groups for those cities
- Same historical ERA5 / central-coefficient / fixed-death validation logic as the AT001C fixture
- The sample is used only for absolute-difference screening across regions and ages.

## Current diagnosis
- Max absolute total difference on the sample: 6.003e-03
- Mean absolute total difference on the sample: 1.874e-03

## Interpretation
- The code path is internally unit-consistent: raw daily sums are scaled to annual values by the year-specific day count.
- Leap-year annualization is not the source of the failure; using a fixed 365-day divisor makes the AT001C mismatch larger.
- Coefficient-order, basis, and centering mistakes were separately tested and produce much larger errors than the current fixture mismatch.
- The remaining discrepancy is therefore a genuine reference mismatch against the Masselot cityage fixture, not a tolerance artifact.

## Status
FAIL
