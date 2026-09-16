# Phase 0 — Summary

## Gate classifications
- Population growth scaling: PASS
- ASSR-to-annual-deaths methodology: PASS
- Attribution units: FAIL
- Clamp validation: PASS
- Signed grouped-AN decision: PASS
- Lloyd convergence: PASS
- Coefficient draws: PASS
- Leap-day policy: PASS
- Negative grouped-AN census: BLOCKED
- PCLM zero/NA behavior: BLOCKED
- IPF summary-column behavior: BLOCKED

## Evidence summary
Population growth scaling preserves city shares and sample coverage to numerical precision, and the intended growth-ratio formula matches the prototype exactly.
The ASSR bridge reproduces Spain projected deaths to floating-point tolerance using the local five-year survival-ratio bridge and the audited city calibration path.
Attribution units show that the raw daily AF×deaths accumulation is about 365× the annualized result, but the Masselot comparison exceeds the explicit 1e-6 threshold, so the attribution-units gate fails on this fixture.
Clamp comparison shows unclamped AF is closer to the Masselot city-age totals for total and cold attributable numbers, while heat is a tie at displayed precision.
Signed grouped-AN Method A is numerically identical to Method B under the current linear split/recombine rule; Method C changes the estimand on negative cells. The remaining grouped-AN blocker is the broader negative-cell census, not the Method A/B choice itself.
Lloyd N=400 and N=800 both satisfy the closure and consecutive-N thresholds on the fake fixture and on the reference schedules; N=50 was rejected because it fails the regression thresholds.
The published coefficient-draw file contains 1000 draws and is tractable for the audited pipeline.
Leap-day annualization is fixed to actual calendar days (365 or 366), with only small sensitivity versus a fixed-365 rule on the audited fixture.
The clean-room validation pack reproduced the current validation pack exactly at the checksum level for all 43 generated artifacts.
The remaining unresolved blockers are the broader negative grouped-AN census, PCLM zero/NA behavior, and IPF summary-column behavior.

## Recommendation
Full implementation should **not** begin yet. The validated gates now support the five-script specification and the production choices for ASSR, coefficient draws, leap days, signed grouped AN, and Lloyd N=400, but the spec still has unresolved blockers and one failed reproduction gate that require a separate scientific decision.

## Status
BLOCKED
