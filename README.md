# Temperature-Mortality → Life Expectancy & Lifespan Inequality

A refactored pipeline that computes **temperature-attributable mortality** from CMIP6 climate projections under multiple SSPs and GCMs, and decomposes the impact of changing temperatures on **remaining life expectancy at age 65 (`e65`)** and **lifespan inequality (LI) above age 65** across 854 European cities.

## Project structure

Two separate but linked analytical branches:

- **Branch A — Attributable burden / absolute risk**: climate-only attributable numbers (ANs) under fixed vs rising longevity, answering how much of the change in temperature-related mortality risk is due to demographic ageing.
- **Branch B — Life-table and decomposition**: full mortality schedules → period life tables → `e65` and LI → Aburto/Horiuchi decomposition by age and temperature range over time.

## Method overview

The pipeline follows a standard attributable-risk framework (Gasparrini & Leone 2014, Masselot et al. 2023) using:

- **City-specific ERFs** from the MCC study, represented as B-spline basis coefficients (b1–b5) per age group, estimated from a two-stage meta-analysis of 854 European cities.
- **Daily mean temperature** from 21 CMIP6 GCMs under three SSP scenarios (SSP1-2.6, SSP3-7.0, SSP5-8.5), downscaled and bias-corrected to city level (ISIMIP3b).
- **EUROSTAT demographic data**: historical age-specific death rates (`demo_mlifetable`) and projected population (`proj_19np`, EUROPOP2019).
- **Period life tables** from age 65 onward, using the Chiang method for `qx` conversion.
- **Lifespan inequality**: SD of attained age at death above 65:

$$ \text{SD} = \sqrt{\frac{\sum d_x (x + a_x - \bar{A})^2}{l_{65}}}, \quad \bar{A} = \frac{\sum d_x (x + a_x)}{l_{65}} $$

- **Decomposition**: Horiuchi linear integral decomposition of annual changes in `e65` and LI by age and five temperature-defined causes (residual non-temperature, extreme cold, moderate cold, moderate heat, extreme heat).

## Pipeline: script-by-script

All scripts are in `scripts/`, numbered by execution order. Validation diagnostics for each step are in `validation/`.

### Step 00 — Setup (`00_Packages_Parameters.R`)

Defines global constants only — no modelling assumptions:

```r
ssp_keep <- c("hist", "1", "2", "3")          # SSP1-2.6, SSP3-7.0, SSP5-8.5
agebreaks <- c(20, 45, 65, 75, 85)            # Masselot age groups
varfun <- "bs"; vardegree <- 2; varper <- c(10, 75, 90)  # B-spline basis
cold_extreme_pct <- 0.025; heat_extreme_pct <- 0.975     # Range thresholds
```

Controls analysis scope via `city_subset`: set to `"ES001C"` for Madrid-only demo, or `NULL` for all 854 cities.

### Step 01 — Raw input preparation (`01_PrepData.R`)

Loads and validates raw inputs only. No future mortality schedules are constructed here.

**Inputs:**
- `data/coefs.csv`: ERF B-spline coefficients (b1–b5) per city × age group.
- `data/coef_simu.csv`: 1000 simulated coefficient vectors per city × age group for uncertainty.
- `data/city_results.csv`: City metadata, baseline demography (population, deaths by age group), MMT, MMP.
- `data/tmeanproj.gz.parquet`: Daily mean temperature for 854 cities × 21 GCMs × 3 SSPs × 1981–2099 (3.2 GB).
- Eurostat API: `proj_19np` (population projections) and `demo_mlifetable` (historical death rates).

**Outputs (RDS to `results/production/`):**
- `coefs_all.rds`, `coef_simu.rds`, `city_res.rds`
- `pop_eu_raw_filtered.rds`, `lt_eu_raw_filtered.rds`

**Validation (01):** ERF curves (RR vs temperature by age group), city map, baseline age distribution, historical mx trends.

### Step 02 — Climate-only AN computation (`02_ComputeAN.R`)

Computes **attributable numbers under fixed baseline deaths** — the epidemiological attribution layer only.

**Method:**

For each city × GCM × SSP × age group:

1. Build prediction basis matching the ERF estimation (B-spline, degree 2, knots at 10/75/90 percentiles of historical temperatures).

2. Find **MMT** (minimum mortality temperature) by searching for the minimum log-RR within the 25th–99th percentile range of the **empirical** temperature distribution (not the artificial temperature grid):

```r
mmt_lower <- quantile(hist_temps, 0.25, na.rm = TRUE)
mmt_upper <- quantile(hist_temps, 0.99, na.rm = TRUE)
mmt <- temp_seq[which(log_rr_seq == min(log_rr_seq[temp_seq >= mmt_lower & temp_seq <= mmt_upper]))]
```

3. Compute daily attributable fraction (Gasparrini & Leone 2014):

$$ \text{AF}_t = 1 - \exp(-\text{ basis}(T_t) \cdot \beta + \text{ basis}(\text{MMT}) \cdot \beta) $$

Negative AFs truncated to zero (range approach).

4. Compute daily AN = AF × (annual_baseline_deaths / 365), aggregate annually.

5. **Temperature range classification** using historical percentiles:

| Range | Definition |
|-------|-----------|
| Extreme cold | T < p2.5 |
| Moderate cold | p2.5 ≤ T < MMT |
| Moderate heat | MMT ≤ T ≤ p97.5 |
| Extreme heat | T > p97.5 |

6. **Annual uncertainty**: For each of the 1000 simulated coefficient draws, recompute annual ANs per temperature range (not aggregated across all years). CIs from 2.5th and 97.5th percentiles of the per-year simulation distribution.

**Key fix from original**: AF truncation now consistent between point estimate and simulations. Uncertainty computed per year, not as a single CI pasted across all years.

**Output:** `ans_annual_all_cities.csv` (city × year × GCM × SSP × age_group × temp_range × an_est/an_low/an_hi).

**Validation (02):** AN time series by range, age profile, MMT check vs reference, temperature distribution with thresholds.

### Step 03 — Mortality backbone (`03_MortalityBackbone.R`)

Builds **future country × year × age`mx` schedules** (non-temperature mortality) using an explicit, documented method.

**Method:**

For each country and single age (20–100):

1. Extract historical death rates from Eurostat `demo_mlifetable` (1990–2019).
2. Fit log-linear trend:

$$ \ln(mx_{a,t}) = \alpha_a + \beta_a t + \varepsilon $$

where $\beta_a$ is the age-specific log-linear slope (annual improvement rate).

3. Set baseline year = most recent historical year with data.
4. Project forward:

$$ mx_{a,y} = mx_{a, \text{baseline}} \cdot \exp[\beta_a (y - y_{\text{baseline}})] $$

5. Handle open age (100+) via forward fill.

**Why not the old method**: The previous `load_eurostat_mortality()` inferred mortality improvements from population cohort ratios ($P_{a+1, t+1} / P_{a, t}$), which confounds mortality change with migration and compositional shifts. The log-linear trend on observed death rates is transparent and avoids this confound.

**Output:** `mortality_backbone.rds` (country × year × age × mx_proj).

**Validation (03):** Historical + projected mx curves, improvement rates by age, mx heatmap.

### Step 04 — City mortality schedules (`04_CityMortalitySchedules.R`)

Translates the national mortality backbone into city-level all-cause mortality schedules.

**Method:**

1. Disaggregate city baseline population and deaths from the 5 Masselot age groups to single ages (20–100) using PCLM (Penalized Composite Link Model, Rizzi et al. 2015).

2. Compute city baseline mx:

$$ mx_{a}^{\text{city}} = \frac{\text{deaths}_a}{\text{population}_a} $$

3. Compute national improvement ratio:

$$ r_{a,y} = \frac{mx_{a,y}^{\text{national}}}{mx_{a, y_{\text{baseline}}}^{\text{national}}} $$

4. Apply to city:

$$ mx_{a,y}^{\text{city}} = mx_{a}^{\text{city}} \times r_{a,y} $$

This approach preserves city-specific baseline mortality differences while applying national demographic trends — consistent with the principle that mortality improvement is primarily a national-level process in Eurostat projections.

**Output:** `city_mortality_schedules.rds` (city × year × age × population × mx × deaths).

**Validation (04):** City vs national mx, projected deaths over time, population age structure.

### Step 05 — AN disaggregation (`05_Disaggregate.R`)

Converts grouped ANs (5 age groups) to single ages (20–100) using PCLM, preserving grouped totals.

**Method:**

For each city × year × GCM × SSP × temperature range, apply PCLM with age breaks at 20, 45, 65, 75, 85 and a last interval of 16 years (85–100).

**Key improvements from original:**
- Age groups explicitly ordered before PCLM.
- Grouped totals checked before and after disaggregation (warning if >5% mismatch).
- Missing age groups trigger a warning instead of silent skip.
- Uncertainty bounds disaggregated using PCLM-derived weights.

**Output:** `ans_single_age_all_cities.csv` (city × year × GCM × SSP × age × temp_range × AN).

**Validation (05):** Single-age AN profile, total preservation histogram.

### Step 06 — Cause-specific bookkeeping (`06_AnalysisDataset.R`)

The critical identity-enforcement step. Combines projected all-cause mortality with temperature-attributable ANs.

**Identity:**

$$ \text{deaths} = \text{residual} + \text{extreme\_cold} + \text{moderate\_cold} + \text{moderate\_heat} + \text{extreme\_heat} $$

**Rules:**
- All-cause deaths from city mortality schedules (Step 04).
- Temperature-attributable deaths from single-age ANs (Step 05).
- If residual < −0.001 for any row, the script **stops with an error** (no silent truncation).
- Small negatives within tolerance are clipped to zero.

**Output:** `analysis_dataset_all_cities.csv` (city × GCM × SSP × year × age × population × mx × deaths × residual + 4 temp range columns).

**Validation (06):** Stacked death decomposition area plot, residual age profile, identity check line.

### Step 07 — Life tables for ages 65+ (`07_LifeTables.R`)

Builds period life tables from age 65 onward. All calculations are **explicit in the script** — no hidden package helpers.

**Life-table functions:**

$$ q_x = \frac{n_x m_x}{1 + n_x (1 - a_x) m_x}, \quad q_{\omega} = 1 $$

$$ l_{x+n} = l_x (1 - q_x), \quad d_x = l_x q_x $$

$$ L_x = n_x (l_x - (1 - a_x) d_x), \quad L_{\omega} = \frac{l_{\omega}}{m_{\omega}} $$

$$ e_x = \frac{T_x}{l_x}, \quad T_x = \sum_{t \geq x} L_t $$

**LI definition (SD of attained age at death above 65):**

$$ \text{LI} = \sqrt{\frac{\sum d_x (x + a_x - \bar{A})^2}{l_{65}}}, \quad \bar{A} = \frac{\sum d_x (x + a_x)}{l_{65}} $$

**Output:** `lifespan_inequality_all_cities.csv` (city × GCM × SSP × year × e65 × sd).

**Validation (07):** e65 trajectory by GCM/SSP, SD trajectory, mx vs age profile, life-table column faceting.

### Step 08 — Decomposition (`08_Decomposition.R`)

Decomposes annual changes in `e65` and LI by age (65+) and temperature cause using the **Horiuchi linear integral method**.

**Method:**

For each consecutive year-pair (t, t+1) within a city-GCM-SSP:

1. Build midpoint mortality schedule (average of t and t+1 cause-specific mx).

2. Compute baseline `e65` and SD from midpoint total mx.

3. For each age $x \geq 65$ and each cause $c$:

   - Perturb total mx at age $x$ by a small $\varepsilon$:
   
   $$ \text{mx}_{\text{pert}}(x) = \text{mx}_{\text{mid}}(x) + \varepsilon $$
   
   - Compute perturbed `e65` and SD.
   - Numerical sensitivity:
   
   $$ \frac{\partial I}{\partial \mu_c(x)} \approx \frac{I_{\text{pert}} - I_{\text{mid}}}{\varepsilon} $$
   
   - Contribution:
   
   $$ \Delta I_c(x) = \frac{\partial I}{\partial \mu_c(x)} \times \left( \mu_{c, t+1}(x) - \mu_{c, t}(x) \right) $$

4. Causes sum to total change at each age by construction (since $\sum_c \Delta \mu_c(x) = \Delta \mu_{\text{total}}(x)$).

5. Annual contributions summed into reporting periods (2030, 2050, 2090).

**What this methodology avoids:**
- Decomposing the within-year climate vs no-climate gap (wrong estimand).
- Decomposing from age 20 (wrong target).
- Using AN shares to allocate age effects to temperature ranges.

**Outputs:**
- `decomposition_annual_all_cities.csv`: annual contributions.
- `decomposition_period_all_cities.csv`: period-summed contributions.

**Validation (08):** Age-cause decomposition bar plots for e65 and SD, cause-summary punch card, additivity check.

### Step 09 — Aggregation (`09_Aggregation.R`)

Aggregates city-level outputs to country and Europe level.

**Method:**
- **Deaths/ANs**: summed across cities within country.
- **Life tables**: population-weighted average of city-level e65 and SD.

**Outputs:**
- `country_analysis.csv`, `europe_analysis.csv`
- `country_lifespan_inequality.csv`, `europe_lifespan_inequality.csv`
- `country_decomposition.csv`, `europe_decomposition.csv`

### Step 10 — Absolute risk analysis (`10_AbsoluteRisk.R`)

Branch A analysis comparing ANs under fixed vs rising longevity.

**Method:**
- Fixed longevity ANs from Step 02 (fixed baseline deaths).
- Rising longevity ANs = same attributable fraction (AF = AN / baseline_deaths) applied to projected all-cause deaths from city schedules.
- Summary: percent change in total ANs due to demographic change.

**Outputs:**
- `absolute_risk_comparison.csv`: AF and ANs under both scenarios.
- `absolute_risk_summary.csv`: percent change by city-GCM-SSP-year.

## Validation scripts

Located in `validation/`, numbered to match pipeline steps. Each script:
1. Loads the corresponding step's outputs.
2. Prints diagnostic summary statistics.
3. Generates labelled plots with detailed captions explaining the method and interpretation.

Plots are saved to `results/validation_plots/`.

## How to run

### Madrid demo mode (single city for validation)

```r
# In scripts/00_Packages_Parameters.R, set:
city_subset <- "ES001C"   # Madrid

# Then from the project root:
Rscript scripts/00_RunAll.R

# Or via the shell wrapper:
bash run_pipeline.sh
```

### Full production mode (all 854 cities)

```r
# Set in scripts/00_Packages_Parameters.R:
city_subset <- NULL

Rscript scripts/00_RunAll.R
```

### Step-by-step validation (recommended for development)

Run each pipeline step, then its validation script, inspecting plots before proceeding:

```r
source("scripts/00_Packages_Parameters.R")
source("scripts/01_PrepData.R")
source("validation/01_validate_prep.R")
# Review plots in results/validation_plots/
source("scripts/02_ComputeAN.R")
source("validation/02_validate_an.R")
# ... etc
```

## Data sources

| File | Source | Description |
|------|--------|-------------|
| `data/coefs.csv` | Masselot et al. 2023 (Zenodo) | B-spline coefficients b1–b5 per city × age group |
| `data/vcov.csv` | Same Zenodo | Lower-triangle 5×5 variance–covariance per city × age |
| `data/city_results.csv` | Same Zenodo | City metadata, baseline deaths/pop, MMT, RR, AF |
| `data/tmeanproj.gz.parquet` | ISIMIP3b CMIP6 (3.2 GB) | Daily mean temperature, 854 cities × 21 GCMs × 3 SSPs |
| `data/coef_simu.csv` | Same Zenodo (470 MB) | 1000 simulated coefficient sets per city × age |
| `proj_19np` (Eurostat API) | EUROPOP2019 | Projected population by age, sex, country, year |
| `demo_mlifetable` (Eurostat API) | Eurostat | Historical death rates by age, sex, country |

## References

- Gasparrini A, Leone M. "Attributable risk from distributed lag models." *BMC Medical Research Methodology* 14:55, 2014. [DOI: 10.1186/1471-2288-14-55](https://doi.org/10.1186/1471-2288-14-55)
- Masselot P et al. "Excess mortality attributed to heat and cold: a health impact assessment study in 854 cities in Europe." *The Lancet Planetary Health* 7(4):e271–e281, 2023. [DOI: 10.1016/S2542-5196(23)00023-2](https://doi.org/10.1016/S2542-5196(23)00023-2)
- Horiuchi S, Wilmoth JR, Pletcher SD. "A decomposition method based on a model of continuous change." *Demography* 45(4):785–801, 2008. [DOI: 10.1353/dem.0.0033](https://doi.org/10.1353/dem.0.0033)
- Rizzi S, Gampe J, van der Gaag N. "An estimator for the pairwise score function." *Demographic Research* 32:625–656, 2015. [DOI: 10.4054/DemRes.2015.32.21](https://doi.org/10.4054/DemRes.2015.32.21)
- Preston SH, Heuveline P, Guillot M. *Demography: Measuring and Modeling Population Processes*. Blackwell, 2001.
- Chiang CL. *The Life Table and Its Applications*. Krieger, 1984.