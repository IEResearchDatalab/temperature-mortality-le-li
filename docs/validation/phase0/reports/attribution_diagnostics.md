# Phase 0 — Attribution diagnostics

## Verdict
FAIL

## Exact AT001C/20-44 fixture
- Reproduced total: 8.915557
- Masselot total: 8.915087
- Absolute difference: 0.000470
- Relative difference: 5.277e-05
- 1990 raw/annualized ratio: 365.0
- The full-period raw sum / mean annual total is a different scope and is not used as the unit check.

## Contract artifacts
- `attribution_diagnostics_checks.csv` and `attribution_diagnostics_failures.csv` are written alongside this report.

## Counterfactuals tested
Index: <variant>
        variant reproduced     abs_diff     rel_diff
         <char>      <num>        <num>        <num>
1:      current   8.915557 4.704293e-04 5.276778e-05
2:     fixed365   8.920786 5.699142e-03 6.392694e-04
3:    no_center  88.022311 7.910722e+01 8.873410e+00
4: reverse_beta  17.703602 8.788515e+00 9.858025e-01

## Diagnosis
- The current code path is unit-consistent; the raw daily sum is reduced by the year-specific day count and the 365× scale is expected.
- Replacing the actual-day divisor with fixed 365 makes the match worse, so leap-year handling is not the source of the failure.
- Removing centering or reversing coefficient order produces much larger errors, so the basis/centering/coefficient order logic is not the source.
- The mismatch is therefore a genuine reference discrepancy against the Masselot fixture, not a tolerance artifact.

## Status
FAIL
