# Temperature-Mortality to Life Expectancy & Lifespan Inequality

This repository reproduces and extends the health impact projection methodology from **Masselot et al. (2025)** to quantify the burden of heat and cold-related mortality across 854 European cities under various climate change scenarios.

## Simulation Pipeline

The work is organized into a robust **8-stage pipeline** located in [R_pipeline/](R_pipeline/), designed for high-performance parallel processing and incremental execution.

1.  **[01_initialize.R](R_pipeline/01_initialize.R)**: Global environment setup. Defines simulation parameters (500 Monte Carlo iterations), core usage (8 cores), and paths to climate/health datasets.
2.  **[02_prep_data.R](R_pipeline/02_prep_data.R)**: Data preparation. Merges city metadata with historical temperature baselines and calculates age-specific thresholds (P2.5, P97.5) and Minimum Mortality Temperatures (MMT).
3.  **[03_attribution.R](R_pipeline/03_attribution.R)**: The core simulation engine. Performs parallelized impact projections for each city. It implements:
    *   **ISIMIP3 Bias Correction**: Detrended monthly quantile mapping for GCM daily temperatures.
    *   **Health Impact Modeling**: Distributed Lag Non-linear Models (DLNM) centered at the city-specific MMT.
    *   **Iterative Resilience**: Skips already-processed cities to allow resuming long-running background jobs.
4.  **[04_aggregate_results.R](R_pipeline/04_aggregate_results.R)**: Aggregation and validation. Summarizes the ~85GB of city-level Monte Carlo simulations into decade-level summaries with 95% uncertainty intervals.
5.  **[05_generate_tables.R](R_pipeline/05_generate_tables.R)**: Country-level aggregation. Generates summary tables in Masselot 2023 Table S6 format and prepares input data for LE/LI analysis.
6.  **[06_le_li_extraction.R](R_pipeline/06_le_li_extraction.R)**: Extract attributable numbers for age groups 65+ comparing baseline (2010s) vs. future (2090s) periods.
7.  **[07_le_li_decomposition.R](R_pipeline/07_le_li_decomposition.R)**: Disaggregate age-grouped ANs to single years using PCLM (Rizzi method) and perform Aburto/Lloyd LE and LI decomposition. ⚠️ *Currently failing due to pclm data type issue.*
8.  **[08_validate_masselot.R](R_pipeline/08_validate_masselot.R)**: **✓ Reproducibility validation** against Masselot et al. 2025 baseline (2000-2014). Confirms 4-range disaggregation sums correctly (r > 0.92). See [results/masselot_validation/VALIDATION_SUMMARY.md](results/masselot_validation/VALIDATION_SUMMARY.md) for detailed results.

### Validation Status

**✓ REPRODUCIBILITY CONFIRMED** — Our 4-range temperature disaggregation successfully reproduces Masselot baseline with excellent spatial correlation:
- **Cold mortality**: r = 0.9745
- **Heat mortality**: r = 0.9241

See [results/masselot_validation/VALIDATION_SUMMARY.md](results/masselot_validation/VALIDATION_SUMMARY.md) for complete validation results.

## Methodology

*   **Temperature Ranges**: Impacts are disaggregated into four mutually exclusive ranges to allow for additive analysis:
    *   **Extreme Cold**: Temperatures below the 2.5th percentile.
    *   **Moderate Cold**: Temperatures between the 2.5th percentile and the MMT.
    *   **Moderate Heat**: Temperatures between the MMT and the 97.5th percentile.
    *   **Extreme Heat**: Temperatures above the 97.5th percentile.
*   **Uncertainty**: Uses 500 Monte Carlo simulations per city/age group to propagate uncertainty from the exposure-response relationship.
*   **Climate Drivers**: Forced by 4 GCMs (GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0) and 3 SSP scenarios (126, 245, 585).

## Data sources

### Files committed to git (in `data/`)

| File | Source | Description |
|---|---|---|
| `coefs.csv` | [Masselot et al. 2023 — Zenodo](https://doi.org/10.5281/zenodo.8320789) | B-spline coefficients (b1–b5) per city × age group for reconstructing ERFs |
| `vcov.csv` | Same Zenodo record | Lower-triangle 5×5 variance–covariance per city × age group |
| `city_results.csv` | Same Zenodo record | City metadata, population, deaths, MMT, MMP, RR, historical attributable fractions |

### Generated results (in `results/`, organized by script)

| Folder | Description |
|---|---|
| `results/attribution_aggregated/` | Final summarized output: decade-level attributable deaths (`an`) with mean and 95% eCI. |
| `results/tables/` | Summary tables in Masselot 2023 (Table S6) and Lloyd 2024 (Table S1) formats. |
| `results/le_li_input/` | Country-level ANs for 65+ age groups prepared for LE/LI decomposition. |
| `results/masselot_validation/` | Reproducibility validation: comparison CSV, scatter plots, and summary report. |

See [results/README.md](results/README.md) for detailed output descriptions.

### Files that must be downloaded separately (in `data/`, ignored by git)

| File | Size | Source | Description |
|---|---|---|---|
| `tmeanproj.gz.parquet` | ~3.2 GB | ISIMIP3b CMIP6 | Daily mean temperature parquet (compressed) for 854 cities. |
| `coef_simu.csv` | 470 MB | Same Zenodo record | 1000 simulated coefficient vectors per city × age group. |

## Usage

The pipeline should be run sequentially from the root of the repository:

```bash
# 1. Initialize environment and check requirements
Rscript R_pipeline/01_initialize.R

# 2. Prepare city data and thresholds
Rscript R_pipeline/02_prep_data.R

# 3. Run attribution simulations (highly recommended to run in background)
nohup Rscript R_pipeline/03_attribution.R > simulation.log 2>&1 &

# 4. Aggregate results and perform validation
Rscript R_pipeline/04_aggregate_results.R

# 5. Generate summary tables
Rscript R_pipeline/05_generate_tables.R

# 6-7. Life Expectancy/Lifespan Inequality analysis
Rscript R_pipeline/06_le_li_extraction.R
Rscript R_pipeline/07_le_li_decomposition.R  # ⚠️ Currently failing

# 8. Validate against Masselot baseline
Rscript R_pipeline/08_validate_masselot.R
```

**Output files** are organized in the `results/` directory:
- `results/attribution_aggregated/` — Decade-level attributable deaths
- `results/tables/` — Summary tables (Masselot S6, Lloyd S1 formats)
- `results/le_li_input/` — Input for LE/LI decomposition
- `results/le_li_decomposition/` — Life expectancy decomposition results
- `results/masselot_validation/` — Reproducibility validation (plots & CSV)

See [results/README.md](results/README.md) for detailed output descriptions.

## References

- Gasparrini A, Leone M. "Attributable risk from distributed lag models." *BMC Medical Research Methodology* 14:55, 2014. [DOI: 10.1186/1471-2288-14-55](https://doi.org/10.1186/1471-2288-14-55)
- Masselot P et al. "Excess mortality attributed to heat and cold: a health impact assessment study in 854 cities in Europe." *The Lancet Planetary Health* 7(4):e271–e281, 2023. [DOI: 10.1016/S2542-5196(23)00023-2](https://doi.org/10.1016/S2542-5196(23)00023-2)
- Rizzi S, Gampe J, van der Gaag N. "An estimator for the pairwise score function." *Demographic Research* 32:625–656, 2015. [DOI: 10.4054/DemRes.2015.32.21](https://doi.org/10.4054/DemRes.2015.32.21)