# Phase 0 — Specification blockers

## Confirmed blockers or unresolved items
1. **Negative grouped-AN census remains blocked.** A global census found negative grouped annual AN cells in the historical baseline. That remains a separate safety blocker until a downstream policy is finalized.
2. **Attribution units fail the explicit reproduction threshold.** The historical Masselot reproduction is above the `1e-6` absolute-difference gate.
3. **PCLM zero/NA behavior unresolved.** The current pipeline utilities handle zero and signed counts with custom guards, but the specification still does not define the approved behavior for all-NA grouped cells and signed counts through PCLM in production.
3. **IPF summary-column behavior unresolved.** The single-age disaggregation step requires a jointly constrained redistribution; the exact summary columns to retain and the failure policy if a cell does not converge are not finalized in the implemented path.
4. **Legacy Lloyd external-output files remain absent.** The repository still does not contain the saved `decomposition_LE_65plus.rds` / `decomposition_LI_65plus.rds` regression files. This is no longer a Phase 0 blocker because the approved gate now uses the validated N=400/N=800 fixture evidence, but the absence should remain documented.

## Resolved or supported items in this phase
- Historical ERA5 temperatures are present in `data/prep_data.RData` as `obs_data$tmean_obs`.
- The published ERA5 CSV matches the RData object exactly.
- Demographic growth-rescaling preserves city shares and sample coverage to numerical precision.
- The ASSR-to-annual-deaths bridge reproduces Spain projected deaths to floating-point tolerance.
- The unclamped AF gate is closer to the Masselot city-age fixture than clamped AF.
- Leap-day handling is now defined as actual-calendar-day annualization.
- Coefficient-draw handling is now defined as use of the published 1000 draws.
- Signed grouped-AN Method A is algebraically equivalent to B under the current linear rule; C is rejected.
- Lloyd N=400 and N=800 both satisfy the closure and consecutive-N thresholds on the fake fixture and reference schedules.

## Status
BLOCKED
