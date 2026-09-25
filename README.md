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
Rscript R_pipeline/05_figures.R
```

| Script | Step |
|---|---|
| `00a_prep_temperature.R` | Builds `data/prep_data.RData`: the observed ERA5-Land series and per-city thresholds (MMT, 2.5th/97.5th percentiles of 1990–2019) |
| `00_demography.R` | Wittgenstein Centre population and deaths (SSP-specific), scaled to the city, disaggregated to single ages 65–100+ |
| `01_attribution.R` | Daily ANs by age group (65–74, 75–84, 85+) and temperature range, with and without climate change |
| `02_single_age.R` | Grouped ANs allocated to single ages 65–100+ |
| `03_master_table.R` | The "dataset for analysis" (Lloyd et al. 2024, Fig S1): population and deaths by cause (4 ranges + rest) by single age |
| `04_le_li_decomposition.R` | Period life tables and LE65/LI65+ levels (`04_le_li_levels.csv`). Horiuchi decomposition by age × cause of the year-on-year change within each branch (`04_le/li_decomposition.csv`), and of the with − without CC difference per 5-year period (`04_between_branch_decomposition.csv`). Steps are cached, so a run can be resumed |
| `05_figures.R` | Summary figures: (1) LE65/LI65+ trajectories with vs without CC and their gap (dual axis); (2) contributions by temperature range, 5-year age band and 20-year block (Lloyd 2024 Figs 3–4 layout); (3) age profile of the climate-change effect, ~2050 vs ~2090. Also `05_summary.csv` with headline numbers, including the CC effect as a % of the LE65 gain |

Each script stops with an error if any of its invariant checks fails. Check results are written to `results/checks/`, outputs to `results/phase1_madrid/`, and diagnostic figures to `results/figures/`.

## Method choices and their source

Every choice below follows a reference implementation. If something deviates, it is flagged here.

| Step | Choice | Source |
|---|---|---|
| Demography | Wittgenstein Centre population and survival ratios (SSP-specific), both sexes summed; annual deaths = pop × (1 − ASSR) / 5, constant within each 5-year period | Masselot 2025 `02_prep_data.R` |
| City calibration | Age-group-specific factor = EUcityTRM city baseline (`city_results.csv`) ÷ **mean national Wittgenstein value over 2000–2014**, fixed over time | Masselot 2025 `02_prep_data.R` |
| Single ages (population, deaths) | PCLM (`ungroup::pclm`, BIC λ, person-scale counts) on the national 5-year bands 65–69 … 95–99, 100+. The open group is spread over 100–110 (`nlast = 11`) and collapsed back into 100+. Checked against observed Eurostat single-age data (LE65 error 0.007 y) | Rizzi et al. 2015; Simon's methods draft 2.5; meeting 10 Sep §2.4 |
| Single-age ANs | Each group's attributable fraction applied to the PCLM single-age deaths, so AN ≤ deaths at every age | Simon's methods draft 2.5, option (b) (confirmed by Daniel 23 Sep) |
| Life tables | Period life tables 65–100+, piecewise-constant hazard; LE65; LI65+ = SD of age at death conditional on reaching 65 | Lloyd 2024 `Code_1.R`, `Code_2.R` (Aburto 2022) |
| Decomposition | Horiuchi (`DemoDecomp`, N = 400) by single age × cause (4 ranges + rest): (i) consecutive years within each branch, summable over periods and ages; (ii) with vs without CC on 5-year-period mean rates, where rest contributes 0 by construction | Lloyd 2024; Simon's methods draft 2.7 |
| ERF basis | `bs`, degree 2; knots at the 10/75/90th percentiles and boundaries at the range of the city's **full ERA5-Land series 1990–2019** (the series the ERFs were estimated on) | Masselot 2025 `03_attribution.R` (`tper`) |
| Temperature projections (with climate change) | ISIMIP3BASD trend-preserving bias correction against ERA5 2000–2014, applied **by month × calibration period** (2015–29, 2030–39, …, 2090–99) | Masselot 2025 `03_attribution.R`, `functions/isimip3.R` |
| Without-climate-change counterfactual | Masselot's `demo` series: each 5-year block of the calibrated GCM series is re-mapped with ISIMIP3 onto the calibrated 2010–2014 distribution. Day-to-day weather is kept and the warming trend removed. Option `counterfactual = "era5_cycle"` (observed 2000–2019 repeated) is kept for sensitivity analysis | Masselot 2025 `03_attribution.R`; Simon's methods draft 2.4.3; Simon 23 Sep ("do whatever Masselot did") |
| Attributable number | Daily AN = (1 − 1/RR) × daily deaths, with **RR clamped at ≥ 1**; the ERF is held constant (no adaptation) | Masselot 2025 `03_attribution.R`; Simon's methods draft 2.4.1 |
| MMT | Recomputed with the Masselot 2025 rule: the argmin of the ERF over percentiles 25–99 of the full ERA5 series. The Masselot 2023 value is kept as `mmt_2023`, and the fixture uses it | Masselot 2025 `03_attribution.R` |
| Temperature ranges | Split at the MMT first, then at the 2.5th/97.5th percentiles of ERA5 1990–2019, fixed over time. MMT > p97.5 means all heat is extreme | Lloyd 2024 `09_0_Attr_Number.R`; Simon's methods draft 2.4.2 (**open point:** the draft says 2000–2014) |
| Calendar | 29 February removed from all daily series; 365-day years | Masselot 2025 `01_pkg_params.R` (`dayvec`) |

Checks built into the scripts: `04` verifies each step's closure, the whole-period closure (Simon's 11 Sep test: the life-table ΔLE65/ΔLI65+ equals the sum of all contributions), between-branch closure, and zero rest contribution. Validation fixture: `01_attribution.R` reproduces the Masselot et al. (2023) published historical heat and cold excess deaths (`city_results.csv`) for the city's 65+ age groups, and stops if the relative error exceeds 0.1%.

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
