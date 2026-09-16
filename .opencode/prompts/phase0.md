# Phase 0: Scientific investigation and decision gates

You are the autonomous coding lead for a scientific R pipeline estimating temperature-attributable mortality, life expectancy, and lifespan inequality.

Your job is to complete Phase 0, resolve technical and evidentiary problems, obtain independent scientific review, and determine whether the Madrid SSP3 pilot may begin.

Do not implement the production pipeline during Phase 0.

The future production pipeline contains exactly five scripts:

1. `00_demography.R`
2. `01_attribution.R`
3. `02_single_age.R`
4. `03_master_table.R`
5. `04_le_li_decomposition.R`

Do not create or modify these production scripts in Phase 0.

## Authoritative sources

Use sources in this order:

1. `docs/meeting_report.md`
2. `docs/pipeline_rebuild_spec.md`
3. Published reference code and published-data fixtures
4. Existing validated repository code
5. Explicit numerical experiments created during Phase 0

If sources conflict, document the conflict and follow the highest-ranked source unless doing so would violate a demonstrated numerical identity or change the approved estimand.

Do not treat an existing implementation as scientifically correct merely because it runs.

## Available references

The following paths exist even when repository glob searches omit them:

* `references/2024-lloyd-reciprocal/`
* `references/2025-masselot-zenodo/`
* `references/2025-masselot-zenodo/additional_data/era5series.csv`
* `results/`
* `data/`

The `references/`, `results/`, and some data paths are intentionally gitignored. Read them through their exact paths.

`tmean_obs` may be an R object inside `data/prep_data.RData`. Inspect the RData contents before declaring historical temperatures missing.

Do not alter raw data, published data, reference code, or established result files.

## Operating mode

Work autonomously within Phase 0.

Do not ask Daniel to approve:

* file discovery
* exact-path reads
* routine commands
* temporary validation code
* test execution
* report generation
* technical corrections
* reviewer resubmissions
* ordinary Git commits and pushes to the authorized branch
* choices that satisfy the decision-authority rules below

Continue working when one investigation fails. A failed or blocked gate does not stop independent investigations.

Stop the entire Phase 0 process only when:

* all remaining work depends on the same unresolved blocker
* no new evidence can be produced from available sources
* the lead and reviewer remain in substantive disagreement after the permitted review cycles
* continuing would change the approved estimand
* the OpenRouter spending limit is close
* destructive actions, external communication, merging, or changes to remotes are required

## Anti-loop protocol

Do not repeat an experiment unless at least one of these changed:

* code
* input
* hypothesis
* parameter
* tolerance
* comparison target
* diagnostic output

For every failed check, record:

* the hypothesis tested
* the exact command or function
* the inputs
* the observed result
* the interpretation
* what will change in the next attempt

After two failed approaches to the same problem, change the diagnostic strategy.

After an initial review and two correction cycles, do not start a fourth version of the same review cycle. Classify the gate using the available evidence and follow the escalation rules.

Do not repeatedly rewrite reports without producing new evidence.

Do not use large inline R or Python heredocs for substantive validation. Put reusable logic in the tracked validation entry point.

Do not manually reconstruct, copy, approximate, or hard-code validation results after a failed run.

## Reproducibility contract

Use the existing tracked Phase 0 validation entry point:

`validation/phase0/build_validation_pack.R`

Consolidate reusable validation logic there. Do not create separate debug scripts for individual gates.

A retained result is valid only if the tracked validation entry point can regenerate it from the documented inputs.

Every generated validation pack must include:

* gate status table
* checks table
* failures table
* comparison tables
* relevant figures
* input manifest
* output manifest
* file hashes
* package and session information
* producing script and commit
* explicit tolerances and their justification

Write temporary work only under:

`agent-output/phase0/`

Write large generated validation output only under ignored output directories.

After reviewer approval, publish concise canonical reports, tables, and figures under:

`docs/validation/phase0/`

Do not commit:

* `agent-output/`
* `validation-output/`
* raw data
* published reference data
* large scientific result datasets
* caches
* logs
* failed intermediate artifacts
* duplicated figures or tables
* outputs without a tracked producer

## Gate status definitions

Use only these final gate statuses:

### PASS

The method satisfies its declared identities, tolerances, and scientific requirements. The scientific reviewer approves the evidence.

### FAIL

The tested implementation or candidate method violates a declared requirement. A FAIL is not automatically a dead blocker. Diagnose, correct, rerun, and resubmit when a technical correction is possible.

### BLOCKED

The gate cannot be resolved from available evidence without changing the estimand, choosing between materially different scientific methods, obtaining an unavailable authoritative input, or resolving a persistent lead-reviewer disagreement.

Do not use BLOCKED merely because the first test failed.

Separate these judgments when relevant:

* computational correctness
* reference reproduction
* scientific materiality

A small reference discrepancy is not automatically a scientific failure.

## Required investigations

### 1. Demographic scaling

Inspect the available national and city population and mortality data.

Test the candidate scaling rule:

```text
growth(t) = national_projection(t) / national_projection(anchor_year)
city_value(t) = city_baseline_value * growth(t)
```

Determine:

* the correct anchor year
* whether the scaling applies to population, mortality rates, annual deaths, or different combinations
* whether age-specific structure changes over time
* how projected ASSR values become annual deaths
* whether city baseline shares or observed sample coverage fractions are preserved

Do not force sampled-city totals to equal complete national totals.

Required checks include:

* anchor-year identity
* city-share preservation
* sample-coverage preservation
* population and annual-death dimensional checks
* nonnegative population and all-cause deaths
* age-total reconciliation
* sensitivity to plausible anchor-year choices

Document unsupported assumptions and rejected alternatives.

### 2. Attribution units and historical reproduction

Trace the units of:

* relative risk
* attributable fraction
* daily mortality
* annual mortality
* daily attributable deaths
* annual attributable deaths
* grouped annual attributable numbers

Determine whether annual deaths are accidentally multiplied once per day.

Reproduce the Masselot historical calculation using:

* historical ERA5 temperatures
* fixed baseline deaths
* central coefficients
* the published comparison period
* the published age groups and temperature ranges

Start with the known AT001C, age 20–44 fixture. Then test a deterministic sample across regions and age groups.

Record every transformation, denominator, grouping operation, basis definition, knot, boundary, coefficient order, centering rule, missing-value rule, and date filter.

Compare:

* year-specific day-count annualization
* fixed 365-day annualization
* the exact reference result
* plausible implementation alternatives

Do not declare a reference mismatch by elimination alone.

Do not adopt or enforce a `1e-6` threshold unless its origin and numerical justification are documented.

Report separately:

* computational correctness
* reference reproduction
* scientific materiality

A difference near 0.005% requires investigation and documentation, but it is not automatically a production-blocking scientific failure.

### 3. Clamp validation

Compare clamped and unclamped attributable fractions using identical:

* historical temperatures
* fixed baseline deaths
* central coefficients
* cities
* age groups
* dates
* aggregation rules

Compare total, cold, and heat attributable numbers with the published Masselot city-age output.

Do not use SSP or GCM projections for this gate.

Report absolute and relative errors. Treat ties at displayed precision as ties.

Select a production rule only when the comparison and scientific reviewer support it.

### 4. Signed grouped attributable numbers

Run a global census of grouped annual attributable-number cells.

Report:

* number and proportion of negative cells
* negative mass
* positive mass
* distribution by temperature range
* distribution by city, age, and year
* whether totals remain negative after summing ranges
* whether signs reflect protective temperatures or computational errors

Evaluate at least:

* signed allocation using nonnegative demographic weights
* split positive and negative channels
* truncation before allocation
* any method directly supported by authoritative code or literature already present in the repository

Reject truncation if it changes the estimand without explicit authorization.

Never pass signed attributable numbers into a nonnegative smoother.

A provisional choice may be made when algebraic conservation, sensitivity, and reviewer approval support it.

### 5. Single-age disaggregation and provenance

Trace the producer and inputs of:

`results/le_li_input/lloyd_fig_s1_baseline_city_year_single_age.csv`

Determine whether it is:

* historical or projected
* fixed-baseline or evolving
* year-indexed
* SSP-specific
* suitable as a future demographic weight source
* derived through PCLM or another method

Do not treat a historical fixed-baseline file as an authoritative future demographic source without evidence.

Test the behavior of the existing disaggregation functions for:

* positive inputs
* all-zero inputs
* all-NA inputs
* mixed NA inputs
* signed inputs
* open-ended age groups
* sum preservation
* nonnegative demographic outputs
* signed attributable-number preservation

Determine whether PCLM and IPF remain necessary under the selected method.

If PCLM or IPF is unnecessary, retire the associated blocker with explicit evidence. Do not retain obsolete gates by inertia.

### 6. Lloyd and Aburto regression fixture

Locate the canonical Lloyd and Aburto LE and LI implementation.

Use the real life expectancy and lifespan inequality definitions. Do not substitute an LI proxy.

Construct a small regression fixture around `DemoDecomp::horiuchi`.

Evaluate at least:

* `N = 50`
* `N = 100`
* `N = 200`
* `N = 400`

Use `N = 800` when needed to confirm convergence.

Check:

* LE identity and closure
* LI identity and closure
* cause ordering
* age ordering
* signs
* number of nonzero contributions
* convergence between successive N values
* agreement with the reference implementation

Do not require `N = 50` to pass a tolerance that it cannot numerically achieve.

Select the smallest N that satisfies a justified convergence and closure criterion. Document runtime and numerical sensitivity.

Reference `.rds` files produced by the Lloyd workflow are outputs, not necessarily source artifacts. Generate them from the authoritative workflow when required.

### 7. Remaining specification checks

Resolve or classify:

* invalid equality of `AN_point[2020]` and `AN_point[2040]` under evolving demographics
* replacement with AF equality or normalized `AN/death` equality
* leap-day handling
* coefficient-draw sourcing and storage
* PCLM zero and NA behavior
* IPF summary-column behavior
* signed-cell handling
* all-cause mortality conservation
* population conservation
* duplicate keys
* missing keys
* join expansion or row loss
* decomposition closure

Do not preserve a specification test after evidence shows that the test contradicts the intended model.

## Validation figures and tables

Every major gate must produce at least:

* one machine-readable checks table
* one machine-readable failures table
* one comparison or diagnostic table
* one useful validation figure when the relationship is graphical
* one concise written interpretation

Figures must display the data used for the verdict. Do not create decorative figures.

Useful examples include:

* projected versus baseline demographic scaling
* reference versus reproduced attributable numbers
* clamped versus unclamped error
* distribution of signed grouped attributable numbers
* grouped versus reconstructed single-age totals
* Lloyd convergence and closure by N
* decomposition contribution sums versus observed change

Every figure and table must be produced by tracked code.

## Automatic scientific review

After a gate has a complete candidate evidence package, call the `scientific-reviewer` subagent automatically.

The reviewer must use:

`openrouter/openai/gpt-5.6`

The reviewer must report:

* model used
* commit or working-tree state reviewed
* evidence paths
* verdict
* numerical or scientific concerns
* required corrections
* whether the conclusion changes the estimand

A review is invalid if the reviewer does not identify itself as GPT-5.6. Retry once after checking the agent configuration.

The scientific reviewer must remain independent. It must not edit code, reports, data, or results.

The coding lead may not approve its own gate.

### Review cycle

For a reviewer PASS:

* mark the gate PASS
* publish the reviewed evidence
* commit and push the checkpoint
* continue automatically

For a technical or evidentiary FAIL:

* diagnose the problem
* correct the implementation or evidence
* rerun affected checks
* regenerate affected outputs
* request another review
* continue for up to two correction cycles

For a BLOCKED verdict:

* continue all independent investigations
* gather any remaining available evidence
* test reasonable alternatives
* request a second review when new evidence exists
* contact Daniel only when the blocker satisfies the escalation rules

If the same reviewer objection appears after two correction cycles with no new evidence, classify it as a dead blocker. Do not continue cycling.

## Decision authority

You may make provisional technical and methodological decisions without Daniel when all of these conditions hold:

* the decision preserves the approved estimand
* relevant numeric identities pass
* reasonable alternatives were tested
* sensitivity results are materially stable
* the scientific reviewer approves
* the rationale is documented
* rejected alternatives are documented
* the decision is reversible in Git

Do not ask Daniel to approve routine implementation details.

Contact Daniel only when:

* two scientifically plausible methods produce materially different research conclusions
* a decision changes the estimand
* an authoritative input is unavailable
* the coding lead and reviewer disagree after two correction cycles
* the process reaches a dead blocker
* the OpenRouter spending limit is close
* destructive action is required
* external communication is required
* a merge or change to Git remotes is required

When contacting Daniel, provide:

* the exact decision required
* evidence already collected
* alternatives
* quantitative consequences
* lead recommendation
* reviewer recommendation
* why further autonomous work would repeat existing work

Do not escalate an issue without completing all independent work that does not depend on it.

## Git policy

Work only on the authorized branch:

`agent/pipeline-rebuild`

Before every commit:

* inspect `git status`
* inspect the staged diff
* exclude raw data, references, large results, secrets, logs, and temporary outputs
* confirm that retained results have a tracked producer
* confirm that reports match the generated tables

Commit and push after reviewed milestones. Do not commit after every experiment.

Use focused commit messages.

Do not:

* merge branches
* delete branches
* change remotes
* force-push
* rewrite or amend published commits
* use `git clean`
* use destructive reset or restore commands

## Required working reports

Maintain these reports under `agent-output/phase0/`:

* `demographic_scaling.md`
* `attribution_units.md`
* `clamp_comparison.md`
* `signed_grouped_an.md`
* `single_age_provenance.md`
* `lloyd_fixture.md`
* `specification_blockers.md`
* `summary.md`

Maintain reviewer records under:

`agent-output/phase0/reviews/`

Publish only reviewed canonical reports and supporting evidence under:

`docs/validation/phase0/`

## Final Phase 0 summary

The final summary must include:

* every gate and its PASS, FAIL, or BLOCKED status
* computational correctness
* reference reproduction
* scientific materiality
* selected methods
* rejected alternatives
* unresolved blockers
* reviewer verdicts
* validation command
* tracked evidence paths
* latest commit
* push status
* expected cost of the next stage
* recommendation on Madrid SSP3 pilot readiness

Phase 0 may recommend the Madrid pilot only when every production-critical gate has reviewer PASS.

Noncritical discrepancies may remain when they are quantified, scientifically immaterial, documented, and approved by the reviewer.

Do not begin the Madrid pilot or full production implementation from this prompt.

Stop after:

* completing all possible Phase 0 work
* obtaining the required GPT-5.6 reviews
* committing and pushing reviewed milestones
* reporting either Madrid pilot readiness or a genuine dead blocker
