# Non-EU coefficient-transfer precursor

`non_eu_transfer_benchmark.R` implements a leave-one-country-out benchmark of
simple ERF coefficient-transfer rules. It advances the non-EU workstream
without claiming that the European validation exercise is a Türkiye result.

## Question

If every city in one European country is treated as unobserved, how accurately
can simple rules reconstruct its central age-specific exposure-response
function from cities in other countries?

The country, rather than an individual city, is the held-out unit. This blocks
the strongest leakage path: a nearby city from the same national mortality and
data system cannot donate its coefficient vector.

## Candidate rules

The script compares four deliberately transparent rules for each age group:

1. **Country-excluded regional mean.** Mean coefficient vector among cities in
   the target city's European region after removing every city in its country.
2. **Nearest ERA5 climate analogue.** Coefficients from the closest city in
   standardized ERA5 mean, standard deviation and 5th, 50th and 95th
   temperature percentiles, excluding the target country.
3. **Nearest geographic city.** Coefficients from the nearest city by
   great-circle distance, excluding the target country.
4. **Country-excluded global mean.** Mean coefficient vector across all cities
   outside the target country.

Donor selection never uses target coefficients, MMT, attributable mortality or
another ERF outcome. ERA5 summaries and coordinates would be observable for a
new geography before its temperature-mortality ERF is estimated.

## Evaluation

The five transferred coefficients are evaluated on each target city's local
temperature-percentile spline basis. This basis exactly reuses the production
settings: quadratic B-splines, knots at the 10th, 75th and 90th percentiles,
and boundaries at the observed ERA5 range.

The primary metrics are:

- root mean squared error of the centered log-relative-risk curve over the 1st
  through 99th temperature percentiles;
- absolute error in the minimum-mortality temperature, recomputed over the
  25th through 99th percentiles.

Coefficient RMSE, coefficient RMSE standardized by between-city variation,
curve correlation and cold/heat tail error are retained as diagnostics. The
script reports method ranks separately for ages 65–74, 75–84 and 85+ and
calculates pairwise rank correlations.

## Run

The large ERA5 file is ignored by Git. On Vangelis, point `ERA5_FILE` to the
authoritative shared copy and keep outputs in the auxiliary results layer:

```bash
ERA5_FILE=/home/SHARED/temperature-mortality-le-li/data/era5series.gz.parquet \
OUT_DIR=results/europe/auxiliary/non_eu_transfer_benchmark_20261008 \
Rscript analysis/non_eu_transfer_benchmark.R
```

The script records server, repository, branch, commit, command, input paths,
input MD5 hashes, outputs and validation in the generated `README.md` and
`input_manifest.csv`.

## Acceptance checks

The run stops unless all of the following hold:

- 854 cities, all countries and all three older age groups are covered;
- every city-age pair has one prediction from each transfer rule;
- all coefficient vectors have five finite dimensions;
- all country-excluded regional pools contain donors;
- nearest donors never belong to the target country or equal the target city;
- ERA5 features exist for every city;
- all primary numerical metrics are finite.

## Interpretation limits

This is an internal transportability benchmark within the 30-country European
study domain. It is not a validation in Türkiye and cannot establish external
validity. It transfers central estimates only, without coefficient covariance.
The climate analogue omits socioeconomic conditions, healthcare, adaptation,
air conditioning and mortality-system quality. European region labels may not
map naturally beyond Europe. A method that ranks first here still requires
external data checks, uncertainty inflation and prospective validation before
operational use in a non-EU VIG market.
