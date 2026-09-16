You are investigating a scientific R pipeline for temperature-attributable
mortality, life expectancy, and lifespan inequality.

Authoritative sources, in order:

1. docs/meeting_report.md
2. docs/pipeline_rebuild_spec.md
3. Existing validated code and published-data fixtures

Execute Phase 0 only. Do not implement Scripts 00 through 05 and do not
perform the full pipeline rebuild.

Required investigations:

1. Demographic scaling
   - Inspect the available national and city population and mortality data.
   - Determine the correct anchor-year scaling formula.
   - The intended candidate is:
     growth(t) = national_projection(t) / national_projection(anchor_year)
     city_value(t) = city_baseline_value * growth(t)
   - Do not force sampled-city totals to equal complete national totals.
   - Test preservation of each city's baseline share or sample coverage
     fraction.
   - Report any unsupported scientific assumptions.

2. Attribution-unit validation
   - Trace the units of AF, annual deaths, and annual attributable deaths.
   - Determine whether the current calculation would accidentally multiply
     annual deaths once per day.
   - Reproduce one known Masselot city-age annual attributable-number result
     using historical ERA5, fixed baseline deaths, and central coefficients.
   - Record every transformation and unit.

3. Clamp validation
   - Compare clamped and unclamped AF using the identical historical inputs,
     fixed baseline deaths, central coefficients, and comparison window.
   - Compare both results with the existing Masselot city-age output.
   - Do not use future SSP or GCM data for this gate.

4. Lloyd regression fixture
   - Locate the existing Lloyd/Aburto LE and LI implementation.
   - Construct a small regression fixture around DemoDecomp::horiuchi using
     N = 50.
   - Confirm that LE and LI identities and signs match the reference code.
   - Do not substitute an LI proxy.

5. Specification checks
   - Flag the invalid test that expects AN_point[2020] to equal
     AN_point[2040] when demographics evolve.
   - Replace it conceptually with an AF equality test or normalized
     AN/death equality test.
   - Identify unresolved PCLM zero/NA, leap-day, coefficient-draw, and IPF
     summary-column behavior.

`tmean_obs` may be an R object stored inside `data/prep_data.RData`, not a
filesystem directory named `data/tmean_obs`. Inspect the RData contents before
classifying historical temperatures as missing.

The published historical-temperature source is also available at:
references/2025-masselot-zenodo/additional_data/era5series.csv

The references/ and results/ directories are intentionally gitignored.
Repository glob results may omit them, but exact-path reads are available.

Rules:

- Do not make scientific decisions silently.
- Do not alter raw or published data.
- Do not modify production pipeline scripts during Phase 0.
- Write temporary validation code and all reports only under
  agent-output/phase0/.
- Do not commit, push, delete branches, or change remotes.
- Do not install system packages.
- Use existing project dependencies where possible.
- If required data or code is missing, document the exact missing path and
  stop that investigation.
- Prefer small reproducible checks over running the entire dataset.

Required outputs:

- agent-output/phase0/demographic_scaling.md
- agent-output/phase0/attribution_units.md
- agent-output/phase0/clamp_comparison.md
- agent-output/phase0/lloyd_fixture.md
- agent-output/phase0/specification_blockers.md
- agent-output/phase0/summary.md

The summary must classify every gate as PASS, FAIL, or BLOCKED and recommend
whether full implementation may begin.

Operate autonomously within the authorized Phase 0 scope.

After completing each investigation, delegate an independent review to the
scientific-reviewer subagent.

If the reviewer finds a technical or evidentiary problem:
- correct it,
- rerun the relevant test,
- request another review,
- repeat for up to two correction cycles.

Do not ask Daniel to approve routine commands, file discovery, test execution,
report writing, or corrections.

Stop and request Daniel's decision only when:
- two scientifically plausible methods require a research choice,
- an authoritative input is unavailable,
- changing an approved scientific definition is necessary,
- the OpenRouter spending limit is close,
- destructive or external actions are required.

A reviewer PASS is required before marking a gate complete.

Decision authority

You are authorized to make provisional technical and methodological decisions
without waiting for Daniel when all of the following hold:

- the decision preserves the approved estimand,
- numeric identities and validation gates pass,
- reasonable alternatives were tested,
- sensitivity results are materially stable,
- the scientific reviewer approves,
- the rationale and rejected alternatives are documented,
- the decision remains reviewable and reversible in Git.

Escalate only when:

- lead and reviewer disagree after two correction cycles,
- methods produce materially different research conclusions,
- the decision changes the estimand,
- authoritative evidence is unavailable,
- cost or token limits are close,
- external communication, destructive action, or merging is required.

Do not ask Daniel to approve routine implementation details.

Stop after producing these outputs.