# Phase 0 — Signed grouped-AN prototype

## Verdict
PASS

## Fixture coverage
- large_negative: UK001C 1992 75-84 ModHeat
- small_negative: UK018C 1997 65-74 ModHeat
- positive_control: FR001C 2008 85+ ModCold

## Key identities
- Method A vs B max LE delta: 0.000e+00
- Method A vs B max LI delta: 0.000e+00
- Method C changes LE on negative cells by up to 3.562e-04
- Method C changes LI on negative cells by up to 1.780e-05

## Interpretation
- Method A and Method B are numerically identical under the current linear split/recombine rule.
- Method A preserves grouped totals and sign while carrying the downstream single-age mortality, rest-deaths, LE, and LI prototype through completion.
- Method C is rejected because it changes the estimand on negative cells.
- Uniform or population weights remain sensitivity alternatives only; they are not adopted because they alter LE/LI materially on the representative cells.
