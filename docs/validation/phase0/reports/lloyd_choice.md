# Phase 0 — Lloyd N=400 vs N=800 on reference schedules

## Verdict
PASS

- Sex/reth schedules evaluated: 2
- Year pairs evaluated: 48

## Runtime
- Elapsed seconds: 664.40

## Summary
   pairs max_le_closure_400 max_le_closure_800 max_li_closure_400
   <int>              <num>              <num>              <num>
1:    48       7.189406e-09       1.796195e-09       1.799325e-09
   max_li_closure_800 max_le_400_800 max_li_400_800
                <num>          <num>          <num>
1:       4.497505e-10   1.737316e-09   2.596732e-10

## Interpretation
- N=400 and N=800 both satisfy the closure and consecutive-N thresholds on the reference schedules.
- N=400 is materially cheaper than N=800 and is therefore the adopted production choice.
- N=50 was rejected because the fake-fixture regression gate fails its thresholds.
