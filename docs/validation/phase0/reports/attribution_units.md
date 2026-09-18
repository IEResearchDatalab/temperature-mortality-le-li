# Phase 0 — Attribution units

## Verdict
PASS

## Target fixture
- City: AT001C
- Age group: 20-44
- Temperature source: historical ERA5 in `data/prep_data.RData` (`obs_data$tmean_obs`)
- Coefficients: central coefficients from `data/coefs.csv` (sim = 0)
- Reference output: `references/2025-masselot-zenodo/results/cityage.csv`

## Unit trace
1. `tmean_obs` — daily temperature in °C.
2. `onebasis(..., fun = 'bs', degree = 2, knots = quantile(ERA5, c(10,75,90)%))` — spline basis, dimensionless.
3. `scale(basis, center = basis_at_MMT)` — centered basis, dimensionless.
4. `log_rr = centered_basis %*% beta` — log relative risk, dimensionless.
5. `af = 1 - exp(-log_rr)` — attributable fraction, dimensionless.
6. `death_annual` — annual baseline deaths, units deaths/year.
7. `af * death_annual` — daily attributable-count contribution on a deaths/year scale before annualization.
8. `sum(af * death_annual) / days_in_year` — annualized per-year contribution.
9. `sum(yearly annualized contributions * days_in_year) / sum(days_in_year)` — year-days weighted full-period annual mean, which matches the Masselot fixture.

## Numeric reproduction
- Raw daily `sum(af * death_annual)` for 1990 in the fixture: 2627.867
- Annualized 1990 total after dividing by days in year: 7.200
- Raw/annualized ratio: 365.0
- Simple unweighted mean annual total across all years: 8.91555730241802
- Year-days weighted annual mean: 8.91508687309519
- Negative daily AF values encountered: 513 (minimum AF = -0.268650)

### Masselot comparison for the full period mean
- Masselot `excess_total_est`: 8.91508687309519
- Reproduced year-days weighted annual mean total: 8.91508687309519
- Absolute difference: 3.553e-15
- Relative difference: 3.985e-16
- Masselot `excess_cold_est`: 4.4116633784664
- Reproduced weighted mean cold total: 4.4116633784664
- Masselot `excess_heat_est`: 4.50342349462879
- Reproduced weighted mean heat total: 4.50342349462879

### Source of the earlier discrepancy
- The earlier `0.000470` difference came from averaging annualized yearly totals with equal year weights instead of year-days weights.
- Once the year-days weighting is used, the fixture reproduces to floating-point tolerance.

## Full reproduced annual means by range
      range mean_annual_total
     <char>             <num>
1: ExtrCold         0.3376098
2: ExtrHeat         1.7086710
3:  ModCold         4.1080256
4:  ModHeat         2.9089232

## Interpretation
- The baseline code does **not** multiply annual deaths once per day in its final output; it explicitly annualizes by dividing by the number of days in each year.
- The raw daily accumulation is roughly 365× larger than the annualized result, which is the expected pre-normalization scale.
- The previous mismatch was computational: it came from an equal-year mean rather than a year-days weighted mean.
- The corrected year-days weighted reproduction matches the Masselot fixture to floating-point tolerance, so the attribution-units gate passes.
- The earlier `1e-6` gate was a post hoc engineering threshold; the corrected comparison is now at floating-point tolerance and does not rely on that threshold.

## Status
PASS
