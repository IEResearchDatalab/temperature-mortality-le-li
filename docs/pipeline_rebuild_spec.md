# Pipeline Rebuild Specification

**Status:** Phase 0 audited draft
**Scope:** Five production scripts only; validation code lives outside the production tree.

## Goals
Build a traceable city-level pipeline from raw data to temperature-attributable mortality, life expectancy, lifespan inequality, and decomposition outputs.

## Production scripts

### 00_demography.R
Purpose: build city-level annual SSP-specific demography and annual deaths.

Inputs:
- `data/wittgenstein_pop.csv`
- `data/wittgenstein_assr.csv`
- `data/city_results.csv`
- `results/projdata/projdata_prototype.csv` as the audited historical bridge reference

Core rules:
- Apply the validated demographic scaling rule separately from the ASSR bridge.
- Growth scaling formula:

```text
growth(t) = national_projection(t) / national_projection(anchor_year)
city_value(t) = city_baseline_value * growth(t)
```

- Preserve each city baseline share or sample coverage fraction; do **not** force sampled-city totals to equal full national totals.
- Derive annual deaths from the audited ASSR bridge only after the annual-death formula has been validated.
- Annual-death bridge formula:

```text
death = 1000 * pop * (1 - assr)
```

- Aggregate over sex to country totals, map Wittgenstein age labels into repository age groups (`20-44`, `45-64`, `65-74`, `75-84`, `85+`), then divide five-year death totals by 5 after aggregation.
- Do not silently impute missing demography.

Outputs:
- `pipeline_v2/00_demography/*.csv` or partitioned tables at city-year-age-SSP level
- `checks.csv`, `failures.csv`, `summary.md`, `overview.png`

Hard checks:
- positive population
- positive deaths
- population growth scaling preserves city shares and sample coverage to numerical precision
- ASSR bridge reproduces Spain projected deaths to floating-point tolerance
- internal reconciliation of annual demography against the source bridge
- hard failure if any city-country-age-SSP cell lacks required support

### 01_attribution.R
Purpose: compute temperature-attributable deaths for each city/year/SSP/GCM/simulation/age group and collapse simulation draws in the same script.

Inputs:
- `pipeline_v2/00_demography/*`
- `data/tmeanproj.gz.parquet`
- historical ERA5 observations from `data/prep_data.RData$obs_data$tmean_obs`
- `data/coefs.csv`
- `data/coef_simu.csv`
- `data/city_results.csv`

Core rules:
- Use the validated unclamped AF rule for historical evidence.
- Preserve the Masselot temperature basis, centering, and temperature-range definitions.
- Collapse simulation draws inside this script; do not create a separate collapsing script.
- Use the published 1000 coefficient draws by default unless runtime or memory evidence shows they are impractical.
- For historical annualization, divide by the actual calendar days in each year (365 or 366), not a fixed 365-day denominator.

Outputs:
- attributable-number tables by city-year-SSP-GCM-sim-age-range
- simulation summary tables (mean, sd, quantiles, point estimate)
- `checks.csv`, `failures.csv`, `summary.md`, `overview.png`

Hard checks:
- historical reproduction against Masselot fixture cells
- no accidental annual-deaths-per-day multiplication
- no unapproved clamping in the validated historical path
- coefficient-draw storage and usage follow the published 1000-draw convention
- leap-day handling matches the approved actual-calendar-day rule

### 02_single_age.R
Purpose: disaggregate grouped ANs to single ages 65:100.

Inputs:
- grouped AN outputs from `01_attribution.R`
- single-age baseline death schedule from the audited validation pack

Core rule:
- Use the validated signed grouped-AN allocation rule: **Method A**, i.e. nonnegative baseline-death weights applied to signed grouped AN totals.
- Method B is algebraically equivalent under the current linear split/recombine rule; keep it as a bookkeeping check only.
- Method C (clamp-to-zero) is rejected because it changes the estimand.

Hard checks:
- grouped-total reconstruction
- sign preservation for signed cells
- nonnegative age-specific all-cause mortality and rest mortality
- stop on any hard failure

Outputs:
- single-age AN tables by city-year-SSP-GCM-age-range
- `checks.csv`, `failures.csv`, `summary.md`, `overview.png`

### 03_master_table.R
Purpose: combine single-age ANs with demography to create the city-level master table.

Inputs:
- `02_single_age.R` outputs
- `00_demography.R` outputs

Outputs:
- city-level master table with causes, populations, deaths, and mortality rates
- `checks.csv`, `failures.csv`, `summary.md`, `overview.png`

Hard checks:
- rest identity
- nonnegative mortality rates
- exact accounting identities for variant-specific totals

### 04_le_li_decomposition.R
Purpose: build life tables, compute LE65/LI, and run Horiuchi decomposition.

Inputs:
- `03_master_table.R` outputs
- canonical LE/LI wrappers from Lloyd reference code
- `DemoDecomp::horiuchi`

Hard checks:
- canonical wrapper identity on the fake fixture
- closure error threshold 1e-7
- consecutive-N convergence threshold 1e-8
- production choice: N=400, validated against N=800 on the fake fixture and on representative reference schedules
- stop if any hard check fails

Outputs:
- LE/LI levels
- temporal decomposition tables
- climate-increment decomposition tables
- `checks.csv`, `failures.csv`, `summary.md`, `overview.png`

## Mandatory validation contract
Every production stage must write:
- `checks.csv`
- `failures.csv`
- `summary.md`
- `overview.png`

A hard failure in `checks.csv` must stop execution immediately and prevent later stages from running.

## Validation pack
Validation and research prototypes live outside the production tree under `validation/phase0/`. They may read saved evidence and generated tables, but they must not be imported by the production scripts.

## Phase 0 outcomes that inform this draft
- population growth scaling: PASS
- ASSR-to-annual-deaths methodology: PASS
- attribution units: FAIL
- clamp comparison: PASS for the historical fixture, with a tie on heat at displayed precision
- signed grouped-AN allocation: PASS for Method A; Method B is equivalent; Method C is rejected
- Lloyd convergence: PASS with N=400 adopted; N=50 rejected
- coefficient draws: PASS; use the published 1000 draws
- leap-day policy: PASS; use actual calendar days
- negative grouped-AN census: BLOCKED
- unresolved PCLM zero/NA behavior: BLOCKED
- unresolved IPF summary-column behavior: BLOCKED
