# Phase 0 — Attribution units

## Verdict
FAIL

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
8. `sum(af * death_annual) / days_in_year` — final annual attributable number, units deaths/year.

## Numeric reproduction
- Raw daily `sum(af * death_annual)` for 1990 in the fixture: 2627.867
- Annualized 1990 total after dividing by days in year: 7.200
- Raw/annualized ratio: 365.0
- Mean annual total across all years: 8.915557
- Negative daily AF values encountered: 513 (minimum AF = -0.268650)

### Masselot comparison for the full period mean
- Masselot `excess_total_est`: 8.915087
- Reproduced annual mean total from historical ERA5 and central coefficients: 8.915557
- Absolute difference: 0.000470
- Masselot `excess_cold_est`: 4.411663
- Reproduced mean cold total: 4.411874
- Masselot `excess_heat_est`: 4.503423
- Reproduced mean heat total: 4.503683

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
- The reproduced city-age annual attributable number does not meet the explicit 1e-6 absolute-difference threshold used in the validation pack, so this gate is classified as FAIL.

## Status
FAIL
