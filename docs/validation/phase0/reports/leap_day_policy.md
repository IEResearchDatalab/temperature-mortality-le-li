# Phase 0 — Leap-day policy

## Verdict
PASS

- Chosen rule: divide annualized daily sums by the actual calendar days in each year (365 or 366).
- Alternative compared: fixed 365-day divisor.
- City/age fixture: AT001C / 20-44

## Sensitivity
    abs_diff           rel_diff        
 Min.   :0.000000   Min.   :0.0000000  
 1st Qu.:0.000000   1st Qu.:0.0000000  
 Median :0.000000   Median :0.0000000  
 Mean   :0.005229   Mean   :0.0006393  
 3rd Qu.:0.000000   3rd Qu.:0.0000000  
 Max.   :0.029146   Max.   :0.0027397  

## Interpretation
- The fixed-365 rule changes only leap-year annual totals and introduces at most a sub-percent shift on this fixture.
- The actual-days rule is the auditable deterministic mapping already used by the validated historical path.
