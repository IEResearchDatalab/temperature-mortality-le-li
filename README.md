# Temperature-Mortality → Life Expectancy & Lifespan Inequality

This pipeline estimates **temperature-attributable deaths** from climate projections and decomposes their effect on **remaining life expectancy at 65 (LE65)** and **lifespan inequality among people aged 65+ (LI65+, SD of age at death)**. The decomposition is by single year of age and by temperature range (extreme cold, moderate cold, moderate heat, extreme heat).

The design follows two reference studies:

1. **Masselot et al. (2025, Nat Med)** for the temperature-attributable deaths (ANs): same data, same code logic. The only change is that heat and cold are each split into moderate and extreme.
2. **Lloyd et al. (2024, Environ Int)**, using the Aburto et al. (2022) code, for the life tables and the Horiuchi decomposition of LE65 and LI65+.

Current stage: **single-city validation (Madrid, ES001C, SSP3-7.0, one GCM, central ERF coefficients)**. After that, the pipeline extends to all 854 cities.

## Pipeline (`R_pipeline/`)

Run the scripts in order from the repository root:

```bash
Rscript R_pipeline/00a_prep_temperature.R
Rscript R_pipeline/00_demography.R
Rscript R_pipeline/01_attribution.R
Rscript R_pipeline/02_single_age.R
Rscript R_pipeline/03_master_table.R
Rscript R_pipeline/04_le_li_decomposition.R
```

| Script | Step |
|---|---|
| `00a_prep_temperature.R` | Builds `data/prep_data.RData`: the observed ERA5-Land series and per-city thresholds (MMT, 2.5th/97.5th percentiles of 1990–2019) |
| `00_demography.R` | Wittgenstein Centre population and deaths (SSP-specific), scaled to the city, disaggregated to single ages 65–100+ |
| `01_attribution.R` | Daily ANs by age group (65–74, 75–84, 85+) and temperature range, with and without climate change |
| `02_single_age.R` | Grouped ANs allocated to single ages 65–100+ |
| `03_master_table.R` | The "dataset for analysis" (Lloyd et al. 2024, Fig S1): population and deaths by cause (4 ranges + rest) by single age |
| `04_le_li_decomposition.R` | Period life tables, LE65, LI65+, and Horiuchi decomposition by age × cause |

Each script stops with an error if any of its invariant checks fails. Check results are written to `results/checks/`, outputs to `results/phase1_madrid/`, and diagnostic figures to `results/figures/`.

## Data

Every input is read from `data/`. Large or third-party files are **not committed**: download them and place them in `data/` yourself (they are listed in `.gitignore`).

### Committed (in `data/`)

| File | Source | Description |
|---|---|---|
| `coefs.csv` | Masselot et al. 2023, Zenodo [10.5281/zenodo.10288665](https://doi.org/10.5281/zenodo.10288665) | B-spline ERF coefficients (b1–b5) per city × age group |
| `vcov.csv` | Same record | Variance–covariance of the coefficients |
| `city_results.csv` | Same record (`results/cityage.csv`) | City metadata, baseline population and deaths, MMT, historical excess deaths (used for calibration and as a validation fixture) |

### Not committed: download into `data/`

All of these come from the Masselot et al. (2025) data archive, Zenodo [10.5281/zenodo.14004322](https://doi.org/10.5281/zenodo.14004322) (`data.zip`). Unzip it and copy the files below into `data/`. The archive's `00_download_data.R` and `codebook.md` document how each file was produced.

| File | Size | Used by | Description |
|---|---|---|---|
| `tmeanproj.gz.parquet` | 3.2 GB | 01 | Daily mean temperature, 854 cities × 21 GCMs × {hist, SSP1-3}, 1990–2099 |
| `era5series.gz.parquet` | 31 MB | 01 | Observed ERA5-Land daily mean temperature per city, 1990–2019 (the series used to estimate the ERFs) |
| `wittgenstein_pop.csv` | 7 MB | 00 | Wittgenstein Centre population by country × sex × 5-year age group × SSP |
| `wittgenstein_assr.csv` | 5 MB | 00 | Wittgenstein Centre age-specific survival ratios, same breakdown |
| `coef_simu.csv` | 470 MB | (uncertainty, not yet used) | Monte Carlo draws of the ERF coefficients (Masselot et al. 2023 record, 10.5281/zenodo.10288665) |

`data/prep_data.RData` is a derived file built by `R_pipeline/00a_prep_temperature.R` from `era5series.gz.parquet` and `city_results.csv`. It is also not committed.

## Legacy code

`notebook/`, `scripts/` and `R/` hold an earlier implementation. It used EUROPOP2019/Eurostat demography downloaded on the fly with the `eurostat` package, and it was used for the July–August Europe runs. That implementation is superseded by `R_pipeline/` and is kept only for reference.

## References

- Gasparrini A, Leone M. Attributable risk from distributed lag models. *BMC Med Res Methodol* 14:55, 2014. doi:10.1186/1471-2288-14-55
- Masselot P et al. Excess mortality attributed to heat and cold: a health impact assessment study in 854 cities in Europe. *Lancet Planet Health* 7:e271–e281, 2023. doi:10.1016/S2542-5196(23)00023-2
- Masselot P et al. Estimating future heat-related and cold-related mortality under climate change, demographic and adaptation scenarios in 854 European cities. *Nat Med* 31:1294–1302, 2025.
- Lloyd SJ et al. The reciprocal relation between rising longevity and temperature-related mortality risk in older people, Spain 1980–2018. *Environ Int* 193:109050, 2024. doi:10.1016/j.envint.2024.109050
- Aburto JM et al. Significant impacts of the COVID-19 pandemic on race/ethnic differences in US mortality. *PNAS* 119:e2205813119, 2022.
- Horiuchi S, Wilmoth JR, Pletcher SD. A decomposition method based on a model of continuous change. *Demography* 45:785–801, 2008.
- Rizzi S, Gampe J, Eilers PHC. Efficient estimation of smooth distributions from coarsely grouped data. *Am J Epidemiol* 182:138–147, 2015.
- Pascariu MD et al. ungroup: An R package for efficient estimation of smooth distributions from coarsely binned data. *JOSS* 3:937, 2018.
