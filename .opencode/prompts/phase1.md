# Phase 1: Madrid SSP3 Diagnostic Implementation

## Scope

Implement the five-script production pipeline for **Madrid, SSP3-7.0, one GCM, central coefficients, both climate variants** (with-CC and without-CC). This is a diagnostic pilot — results are non-final until all checks pass.

Do NOT expand scope beyond this target. No other cities, no coefficient draws, no multi-GCM runs.

## Branch

Continue on `agent/pipeline-rebuild` from commit `c8e4408`. The working tree must be clean before starting.

## Implementation order — strict sequential

Complete each script before starting the next. Commit after each script passes its own checks.

### Script 1: `00_demography.R`

Produce evolving grouped and single-age population and all-cause deaths for Madrid, SSP3-7.0.

Checks:
- Anchored to source data (no invented values)
- Complete city × year × age keys with no duplicates
- All population and death values nonneg and finite
- Grouped age totals = sum of constituent single ages

### Script 2: `01_attribution.R`

Produce with-CC and without-CC grouped signed attributable numbers (AN) for one GCM, central coefficients.

Checks:
- One annualization rule, documented and consistently applied
- Year-specific day counts (no constant 365)
- Coverage (city × year × age-group × temp-range) matches the demographic domain exactly

### Script 3: `02_single_age.R`

Allocate grouped AN to single ages using nonneg demographic weights via corrected Method A.

**Critical fixes — these are the main blockers from Phase 0:**

1. **Age-band mapping must be exact:**
   - Source group 65–74 → allocate to single ages 65, 66, …, 74 only
   - Source group 75–84 → allocate to single ages 75, 76, …, 84 only
   - Source group 85+ → allocate to single ages 85, 86, …, 100 only
   - No cross-band leakage. Ever.

2. **PCLM input contract (fail-fast):**
   - Grouped demographic inputs: finite, complete, nonneg
   - All-zero grouped deaths → return all-zero single-age vector (no PCLM call)
   - Any NA/NaN/Inf/negative demographic value → fail with city, year, scenario, age-group context
   - Zero deaths with nonzero AN → fail
   - Convergence check on PCLM output
   - Exact grouped-total reconstruction within tolerance

3. **Signed AN allocation order:**
   - Construct nonneg demographic weights from PCLM output FIRST
   - Apply signed AN × weights AFTER
   - Never pass signed values into nonneg PCLM

Checks:
- Grouped AN totals preserved within 1e-9 after allocation
- No single-age AN allocated outside its source age band
- All weight vectors sum to 1.0 within 1e-12
- No nonzero AN allocated with zero weight

### Script 4: `03_master_table.R`

Assemble the complete analysis table joining demographics, AN, and single-age outputs.

Checks:
- Primary keys are unique and form a complete grid
- Mortality conservation: total deaths ≥ attributable deaths where rest mortality must be nonneg
- No NA/NaN/Inf in any required column
- With-CC and without-CC branches use identical demographics

### Script 5: `04_le_li_decomposition.R`

Run Horiuchi decomposition with N=400.

Checks:
- Closure: sum of age-specific LE contributions = total LE difference (within 1e-6)
- Closure: sum of age-specific LI contributions = total LI difference (within 1e-6)
- Contribution signs are directionally plausible

## Output contract for each script

Every script produces exactly:

1. Its primary output file(s) under `results/`
2. One checks table: `results/checks/NN_scriptname_checks.csv` — one row per invariant, columns: `check_name`, `status`, `value`, `threshold`
3. One failures table: `results/checks/NN_scriptname_failures.csv` — empty if all pass; otherwise: key columns + `failing_check` + `observed_value`
4. One diagnostic figure: `results/figures/NN_scriptname_diagnostic.png`

Commit message format: `phase1: NN_scriptname passes checks`

## Hard stop criteria

**Halt immediately** if any of these occur. Report the first failing invariant, the affected keys, and observed vs expected values. Do NOT start a general audit or Phase 0 re-investigation.

- Duplicate or missing primary keys
- NA / NaN / Inf in any required field
- Negative population or all-cause deaths
- Nonzero AN with zero allocation weight
- Weight sums outside [1 − 1e-12, 1 + 1e-12]
- Grouped reconstruction error > 1e-9
- Negative rest mortality
- Decomposition closure failure beyond 1e-6
- With-CC and without-CC branches using different demographics

On failure: report the invariant, the affected key values, the observed value, and the expected bound. Then stop. One targeted fix, one rerun of affected scripts. If the same failure repeats, stop and report — do not loop.

## Forbidden actions

- Do NOT re-read, revise, or regenerate Phase 0 validation reports, dashboards, or gate labels
- Do NOT call the scientific reviewer until all 5 scripts complete and pass their own checks
- Do NOT create one-off validation or investigation scripts separate from production code
- Do NOT modify `.opencode/` configuration files during the session
- Do NOT experiment with different thresholds, aggregation rules, or annualization methods — use the documented choices from Phase 0
- Do NOT expand scope to other cities, GCMs, or coefficient draws
- Do NOT spend more than 2 consecutive steps diagnosing any single check failure — report and stop instead
- Do NOT re-read large context files (validation packs, prior reports) — work from the checks tables this session produces

## Known repository context

| Item | Path | Note |
|---|---|---|
| Validation runner | `validation/run_validation.R` | Do not modify |
| Historical single-age file | `results/le_li_input/lloyd_fig_s1_baseline_city_year_single_age.csv` | Historical ERA5 only — NOT usable as future demographic weights |
| Phase 0 evidence | `docs/validation/phase0/` | Frozen — do not modify |
| Agent-output scratch | `agent-output/` | Gitignored |

## After all 5 scripts pass their checks

1. Regenerate the validation pack once from the exact current commit
2. Request one GPT-5.6-sol scientific review
3. If review passes → commit, push, tag as `phase1-madrid-diagnostic`
4. If review finds a defect → fix the specific defect, rerun only the affected script(s) and downstream, re-review once. Maximum one correction cycle. If it fails again, stop and escalate.
