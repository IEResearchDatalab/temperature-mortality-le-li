# Temperature-Mortality → Life Expectancy & Lifespan Inequality

A pipeline that computes **temperature-attributable mortality** from climate projections and decomposes its impact on **life expectancy (LE)** and **lifespan inequality (LI)** by age and temperature range.

## Pipeline

| Notebook | Step |
|---|---|
| `notebook/LI1_AN.Rmd` | Attributable numbers (ANs) by age group and temperature range |
| `notebook/LI2_disaggregate.Rmd` | PCLM disaggregation of ANs from wide groups to single ages |
| `notebook/LI3_analysis.Rmd` | Combine ANs with population and all-cause mortality; period life tables |
| `notebook/LI4_decomposition.Rmd` | Decompose ΔLE and ΔLI by age and temperature range |

## Quick start

```r
# 1. Install dependencies
install.packages(c("data.table", "arrow", "dlnm", "splines",
                   "ggplot2", "scales", "ungroup", "MASS", "eurostat"))

# 2. Render all notebooks (data is downloaded on the fly via Eurostat API)
notebooks <- c("LI1_AN", "LI2_disaggregate", "LI3_analysis", "LI4_decomposition")
for (nb in notebooks) {
  rmarkdown::render(file.path("notebook", paste0(nb, ".Rmd")),
                    output_dir = "results/demo",
                    knit_root_dir = ".")
}
```

## Configuration

Each notebook has a config block at the top. Change `city`, `city_code`, `ssp_label`, and `demo_gcm` to run for a different European city or scenario. All outputs are written to `results/demo/` with filenames derived from the city name.

## Production pipeline (multi-city, multi-GCM, multi-SSP)

After validating with the notebooks, run the full pipeline for all 854 cities:

```bash
Rscript scripts/00_RunAll.R
```

This uses a cascading source chain (Masselot-style) and parallel `foreach` loops over cities. See `scripts/` for the individual numbered scripts.

## Method

Implements the standard attributable-risk framework (Gasparrini & Leone 2014, Masselot et al. 2023) using city-specific exposure-response functions from the MCC study, daily CMIP6 temperature projections, and EUROPOP2019 mortality projections. Life tables follow standard demographic methods. Lifespan inequality is measured via the standard deviation of age at death. Decomposition uses a stepwise replacement approach attributable by age and cause.

## Data sources

### Files committed to git (in `data/`)

| File | Source | Description |
|---|---|---|
| `coefs.csv` | [Masselot et al. 2023 — Zenodo](https://doi.org/10.5281/zenodo.8320789) | B-spline coefficients (b1–b5) per city × age group for reconstructing ERFs |
| `vcov.csv` | Same Zenodo record | Lower-triangle 5×5 variance–covariance per city × age group |
| `city_results.csv` | Same Zenodo record | City metadata, population, deaths, MMT, MMP, RR, historical attributable fractions |

### Files that must be downloaded separately (in `data/`, ignored by git)

| File | Size | Source | Description |
|---|---|---|---|
| `tmeanproj.gz.parquet` | 3.2 GB | ISIMIP3b CMIP6 (see Masselot et al. 2023 data notice) | Daily mean temperature for 854 cities, 21 GCMs, 3 SSPs, 1990–2099 |
| `coef_simu.csv` | 470 MB | Same Zenodo record as `coefs.csv` | 1000 simulated coefficient vectors per city × age group for empirical CIs |

### Data downloaded on the fly via Eurostat API (no file needed)

| Dataset | API code | Source | Description |
|---|---|---|---|
| Population projections | `proj_19np` | [Eurostat EUROPOP2019](https://ec.europa.eu/eurostat/web/population-demography/population-projections/database) | Single-age population by sex, country, year (2019–2100). Used to derive mortality improvement trends. |
| Life tables | `demo_mlifetable` | [Eurostat](https://ec.europa.eu/eurostat/data/database) | Historical age-specific death rates (DEATHRATE) by sex and country (1960–2024). Used as baseline mx. |

The notebooks call `load_eurostat_mortality(country_code, sex)` in `R/load_data.R`, which downloads both datasets via `get_eurostat()` and combines them into a projected mx time series for any European country and gender.

### Repository structure

```
.
├── data/                          # Input data (see table above)
│   ├── coefs.csv
│   ├── vcov.csv
│   ├── city_results.csv
│   ├── tmeanproj.gz.parquet       (download)
│   └── coef_simu.csv              (download)
├── notebook/                      # Validation notebooks (R Markdown)
│   ├── LI1_AN.Rmd
│   ├── LI2_disaggregate.Rmd
│   ├── LI3_analysis.Rmd
│   └── LI4_decomposition.Rmd
├── scripts/                       # Production pipeline (R scripts, Masselot-style)
│   ├── 00_Packages_Parameters.R
│   ├── 01_PrepData.R
│   ├── 02_ComputeAN.R
│   ├── 03_Disaggregate.R
│   ├── 04_AnalysisDataset.R
│   ├── 05_LifeTables.R
│   ├── 06_Decomposition.R
│   └── 00_RunAll.R               # Cascading master
├── R/                             # Shared helper functions
│   ├── rr_basis.R
│   ├── impact.R
│   ├── simulation.R
│   ├── load_data.R
│   ├── load_coefficients.R
│   ├── period_lifetable.R
│   ├── cohort_lifetable.R
│   ├── epv.R
│   ├── isimip3.R
│   └── utils.R
├── references/                    # Original reference code (Gasparrini, Masselot)
├── results/                       # Outputs (gitignored, except .gitkeep)
└── README.md
```

## References

- Gasparrini A, Leone M. "Attributable risk from distributed lag models." *BMC Medical Research Methodology* 14:55, 2014. [DOI: 10.1186/1471-2288-14-55](https://doi.org/10.1186/1471-2288-14-55)
- Masselot P et al. "Excess mortality attributed to heat and cold: a health impact assessment study in 854 cities in Europe." *The Lancet Planetary Health* 7(4):e271–e281, 2023. [DOI: 10.1016/S2542-5196(23)00023-2](https://doi.org/10.1016/S2542-5196(23)00023-2)
- Rizzi S, Gampe J, van der Gaag N. "An estimator for the pairwise score function." *Demographic Research* 32:625–656, 2015. [DOI: 10.4054/DemRes.2015.32.21](https://doi.org/10.4054/DemRes.2015.32.21)