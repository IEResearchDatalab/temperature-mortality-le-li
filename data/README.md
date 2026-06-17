# Data sources

This directory contains small source files needed to run the pipeline.
Large files (`tmeanproj.gz.parquet` ~3.2 GB) must be downloaded separately.

## Included files (committed to git)

| File | Source | Description |
|---|---|---|
| `coefs.csv` | Masselot et al. (2023) | City-age-group spline coefficients for ERFs |
| `city_results.csv` | Masselot et al. (2023) | City metadata, baseline population and deaths by age group |
| `mortality_projections/ES_mortality_male.csv` | EUROPOP2019 | Spain national single-age mx projections (2022–2100) |
| `codebook.md` | — | Full data documentation |

## Files that must be obtained separately

| File | Approx. size | Where to get it |
|---|---|---|
| `tmeanproj.gz.parquet` | 3.2 GB | Contact the project authors or see `00_download_data.R` in the original project for download instructions. Contains daily temperature projections for 854 European cities (CMIP6, 21 GCMs, SSP1-3, 1990–2099). |

## Derived data

The notebooks generate these files automatically in `results/`:

- `{city}_ans_annual.csv` — annual attributable numbers by age group and temperature range
- `{city}_ans_daily.csv` — daily attribution detail
- `{city}_ans_single_age.csv` — single-age ANs (PCLM disaggregation)
- `{city}_analysis_dataset.csv` — combined dataset for life-table analysis
- `{city}_lifetable_inputs.csv` — period life-table inputs by year and age
- `{city}_lifespan_inequality.csv` — lifespan inequality time series
- `{city}_le_decomposition.csv` — LE decomposition by age and temperature range
- `{city}_li_decomposition.csv` — LI decomposition summary
- `{city}_pclm_diagnostics.csv` — PCLM diagnostics
- `{city}_population_single_age.csv` — single-age population distribution