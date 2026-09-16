# Phase 0 — IPF summary-column behavior

## Verdict
PASS

## Evidence base
- `docs/pipeline_rebuild_spec.md` defines the five production scripts.
- `R_pipeline/17_build_lloyd_fig_s1_future_constrained.R` is the only repository file that defines and uses `ipf_with_margins()`.
- The approved production design selects Method A for signed grouped-AN allocation and does not use the future-constrained IPF branch.

## Code trace
1. `R_pipeline/17_build_lloyd_fig_s1_future_constrained.R` contains the only IPF implementation in the repository.
2. No production script among `00_demography.R` through `04_le_li_decomposition.R` calls IPF.
3. The single-age baseline producer `R_pipeline/11_build_lloyd_fig_s1_baseline.R` uses PCLM helpers only; it does not call IPF.
4. Therefore the IPF summary-column issue is not part of the approved five-script production estimand.

## Interpretation
- IPF remains relevant only for the rejected/non-production constrained prototype.
- Under the approved five-script production scope, IPF is unnecessary and the summary-column behavior issue does not block pilot readiness.
- This does not change the estimand because the selected method avoids the IPF branch entirely.

## Status
PASS
