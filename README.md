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
# Install dependencies
install.packages(c("data.table", "arrow", "dlnm", "splines",
                   "ggplot2", "scales", "ungroup"))

# Ensure tmeanproj.gz.parquet is in data/ (see data/README.md)
# Then render all notebooks
notebooks <- c("LI1_AN", "LI2_disaggregate", "LI3_analysis", "LI4_decomposition")
for (nb in notebooks) {
  rmarkdown::render(file.path("notebook", paste0(nb, ".Rmd")),
                    output_dir = "results/demo")
}
```

## Configuration

Each notebook has a config block at the top. Change `city`, `city_code`, `ssp_label`, and `demo_gcm` to run for a different European city or scenario. All outputs are written to `results/demo/` with filenames derived from the city name.

## Data

See `data/README.md` for data sources and requirements.

## Method

Implements the standard attributable-risk framework (Gasparrini & Leone 2014, Masselot et al. 2023) using city-specific exposure-response functions from the MCC study, daily CMIP6 temperature projections, and EUROPOP2019 mortality projections. Life tables follow standard demographic methods. Lifespan inequality is measured via the standard deviation of age at death. Decomposition uses a stepwise replacement approach attributable by age and cause.