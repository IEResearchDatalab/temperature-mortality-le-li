# Temperature-Mortality → Life Expectancy & Lifespan Inequality

This pipeline estimates **temperature-attributable deaths** from climate projections and decomposes their effect on **remaining life expectancy at 65 (LE65)** and **lifespan inequality among people aged 65+ (LI65+, SD of age at death)**. The decomposition is by single year of age and by temperature range (extreme cold, moderate cold, moderate heat, extreme heat).

The design follows two reference studies:

1. **Masselot et al. (2025, Nat Med)** for the temperature-attributable deaths (ANs): same data, same code logic. The only change is that heat and cold are each split into moderate and extreme.
2. **Lloyd et al. (2024, Environ Int)**, using the Aburto et al. (2022) code, for the life tables and the Horiuchi decomposition of LE65 and LI65+.

Current stage: **SSP3-7.0 is complete for all 854 cities × 19 GCMs and all 35 pooled geographies** (30 countries, four regions, Europe), using central ERF coefficients. SSP2-4.5 and SSP1-2.6 are incomplete in production; see *Limitations and warnings*. The pipeline was validated on Madrid (ES001C) first; that single-city stage is tagged `phase1-madrid-diagnostic`.

## Limitations and warnings

> **Read before using any result.** Status on 9 October 2026.

- **Central estimates only.** All results use the central ERF coefficients and the mean over 19 GCMs. No uncertainty intervals exist yet (coefficient draws, GCM spread and PCLM sensitivity are planned in issues #4–#12).
- **SSP1-2.6 is incomplete in production.** The Wittgenstein (WCDE v2) survival ratios are published to three decimals. Late in the century many old-age values are exactly `1.000`, which implies zero deaths in whole broad age groups (65–74 and, in several countries, 75–84). This affects 658 of the 854 cities, in 18 countries: AT, BE, CH, DE, DK, EL, ES, FI, FR, IE, IT, LU, NL, NO, PT, SE, SI and UK. Their demography stops at Part 00, and no pooled SSP1 geography exists. The unrounded WCDE v2 source has not been found in the public releases. Runs that replace exact ones with `0.9995` (all cities) or `0.99999` (Madrid pilot) are exploratory sensitivities only, **not a production fix**. The full `0.9995` run exceeds the 5% mortality-shape threshold in 562 cities.
- **SSP2-4.5 production excludes the 12 Swiss cities.** Their reconstructed single-age mortality falls by about 11% between consecutive ages, failing the 5% mortality-shape check, so no pooled SSP2 geography exists in production. A run with that check bypassed covers all 854 cities. It is an inclusion-versus-exclusion sensitivity and **does not validate the Swiss reconstruction**.
- **The 5% mortality-shape threshold is an internal project safeguard,** not an established demographic standard. Neither Masselot's code nor the PCLM method uses it. It flags any decline of more than 5% in single-age mortality between consecutive ages within a broad age group.
- **Sensitivity runs are labelled and isolated.** They carry `WARNING` rows in their Part 00 check tables and in `00_sensitivity_warnings.csv`, record the switches in `run_manifest.txt`, and are written below `results/europe/sensitivities/<name>/`, never into the canonical `results/europe/ssp<k>`, `collected` or `geographies` outputs. See "Explicit sensitivity runs" below.

## Pipeline (`pipeline/`)

The scientific pipeline has exactly five ordered scripts. This deliberate structure stays close to Masselot's code and makes each transformation easy to review. `pipeline/00_pkg_params.R` contains configuration and is sourced by each step; it is not a sixth step.

| Script | Scientific step | Established output |
|---|---|---|
| `00_prep_data.R` | Prepare the common ERA5 series and thresholds, and prepare SSP-specific city demography at grouped and single ages | `data/prep_data.RData`, `00_demography_grouped.csv`, `00_demography_single_age.csv` |
| `01_attribution.R` | Calculate grouped temperature-attributable deaths for one GCM, or the ensemble mean when `GCM=ENSEMBLE` | `01_attribution_grouped.csv` |
| `02_single_age.R` | Allocate grouped attributable deaths to single ages 65–100+ | `02_single_age_an.csv` |
| `03_master_table.R` | Construct the cause-specific analysis dataset for a city; collect Object 1; or sum city counts to a country, region or Europe | `03_master_table.csv`, `object1_dataset.parquet` |
| `04_le_li_decomposition.R` | Construct LE65/LI65+ and Horiuchi decompositions from any Part 03 table; collect Objects 2 and 3 | `04_le_li_levels.csv`, decomposition CSVs, `object2_contributions.parquet`, `object3_levels.parquet` |

The same Parts 03 and 04 implement every geographical level. For pooled geographies, Part 03 sums population and deaths by cause across cities before Part 04 constructs mortality rates, life tables and decompositions. City ERFs, LE/LI values and decomposition contributions are never averaged.

Each step writes invariant checks next to its scientific output and stops on failure. The main pipeline does not generate manuscript figures or tables.

## Batch runs: all cities × 19 GCMs × SSPs

```bash
./run_pipeline.sh 3 all 32                         # SSP3-7.0, all levels, 32 cores
./run_pipeline.sh 3 my_cities.txt 8                # partial city run; pooled levels are skipped
./run_pipeline.sh --config run_config.csv           # enabled SSPs, sequential and resumable
nohup analysis/run_city_figures.sh 3 32 > figures_ssp3.log 2>&1 &
nohup analysis/run_geography_results.sh 3 4 > geography_results_ssp3.log 2>&1 &
Rscript analysis/moderate_heat_erf_diagnostic.R                # ERF/threshold diagnostic requested on 25 Sep
CITY_ID=ES001C Rscript analysis/temperature_pattern_diagnostic.R # 19-GCM temperature pattern diagnostic
CITY_ID=ES001C Rscript analysis/demography_temperature_gap_diagnostic.R # exact two-factor LE-gap diagnostic
```

`run_pipeline.sh` is orchestration only. It repeats the five scripts across cities and GCMs, calculates the ensemble in Part 01, collects Object 1 in Part 03, and then applies Parts 03–04 to countries, regions and Europe. Settings come from environment variables (`CITY_ID`, `GEO_LEVEL`, `GEO_ID`, `SSP`, `GCM`, `OUT_DIR`, `DEMOG_DIR`, …) that override `00_pkg_params.R`.

`run_config.csv` is the tracked production run plan. Each enabled row is one
atomic SSP run; the runner processes those rows sequentially and parallelises
cities and GCMs within the active SSP. Scientific scripts therefore never loop
over scenarios or share scenario state. A failed city/ensemble stage stops
before canonical objects or pooled geographies are updated, and a failed SSP
stops the config before later SSPs begin. The tracked plan enables only
SSP3, at the validated `N_HORIUCHI=400`. SSP1 and SSP2 are disabled because
their production runs are known to stop at Part 00 (see *Limitations and
warnings*); enabling them would stop the config before SSP3. Each scenario
output records the Git revision and resolved city/GCM domains in
`run_manifest.txt`, `run_cities.txt` and `run_gcms.txt`.

Rerunning into an existing scenario root resumes it: finished jobs are skipped.
The runner first checks `run_manifest.txt`. It stops without deleting anything
if `N_HORIUCHI`, `DECOMP_ANNUAL`, the sensitivity switches or the `pipeline/`
code differ from the recorded run, or if `pipeline/` has uncommitted changes.
Otherwise it keeps the original manifest and appends a `resumed_at` record.
A root with finished outputs but no manifest (canonical SSP3 predates
manifests) is reused only with `ACCEPT_UNRECORDED_OUTPUTS=1`, and the new
manifest then says `prior_outputs=unrecorded_provenance`.

### Explicit sensitivity runs

Two opt-in Part 00 settings exist only for the documented SSP1/SSP2 input
sensitivities. `SKIP_MX_PLAUSIBILITY_CHECK=1` continues after the 5% within-age-
group mortality-shape guard, while preserving the observed failure as a
`WARNING`. `ASSR_ONE_REPLACEMENT=<value>` (the exploratory runs used `0.9995`
and `0.99999`) replaces source ratios reported as exactly one and records every
changed row. Neither setting is a production fix. Defaults leave the source and
guard unchanged. The runner refuses either setting unless `BATCH_ROOT_OVERRIDE`
resolves (after symlinks, `..` and trailing slashes) outside `results/europe`
or strictly below `results/europe/sensitivities/`, so sensitivity outputs
cannot overwrite or mix with the canonical outputs.

Other behaviour:

- Output goes to `results/europe/ssp<k>/<city>/{demography,<GCM>,ENSEMBLE}/`.
- Finished jobs leave a `.done` marker, so a run can be restarted and resumes. Failures are listed in `failed.txt`.
- Result scripts under `analysis/` generate figures and manuscript summaries only after the pipeline is complete.
- Horiuchi steps (`N_HORIUCHI`): `run_config.csv` uses N = 400; the runner's command-line default is N = 50; both run only the 5-year-period decompositions (`DECOMP_ANNUAL=0`). A direct `Rscript` run uses the `00_pkg_params.R` defaults, N = 400 with annual decompositions. On Madrid, N = 50 matches N = 400 to 6 decimals, but at N = 50 29 SSP1/SSP2 city ensembles missed the 1e-6 closure tolerance. The validated outputs therefore use N = 400: the SSP3 pooled run, and the SSP2 and SSP1 sensitivity runs (`run_manifest.txt`). Valletta and pooled Malta in the SSP1 sensitivity were rerun at N = 800. Canonical SSP1/SSP2 ran at N = 50.

Parts 03 and 04 write the three data objects agreed on 10 Sep:

| File | Contents |
|---|---|
| `object1_dataset.parquet` | Deaths by cause and population by city × SSP × scenario × year × single age |
| `object2_contributions.parquet` | Horiuchi contributions by city × SSP × age × cause, with vs without CC and between periods |
| `object3_levels.parquet` | LE65 and LI65+ by city × SSP × scenario × year, for the ensemble and each GCM |

Pooled outputs are under `results/europe/geographies/ssp<k>/`. Figures, smoothing, rankings, headline summaries and paper tables are deliberately outside `pipeline/` in `analysis/`.

Cost measured on 2 cores: about 30 s per city × GCM job and 1–2 min per city ensemble. One SSP for all 854 cities is therefore about 140 CPU-hours (roughly 4–5 h on 32 cores), and about 15 GB of disk.

## Method choices and their source

Every choice below follows a reference implementation. If something deviates, it is flagged here.

| Step | Choice | Source |
|---|---|---|
| Demography | Wittgenstein Centre population and survival ratios (SSP-specific), both sexes summed; annual deaths = pop × (1 − ASSR) / 5, constant within each 5-year period | Masselot 2025 `02_prep_data.R` |
| City calibration | Age-group-specific factor = EUcityTRM city baseline (`city_results.csv`) ÷ **mean national Wittgenstein value over 2000–2014**, fixed over time | Masselot 2025 `02_prep_data.R` |
| Single ages (population, deaths) | PCLM (`ungroup::pclm`, BIC λ, person-scale counts) on the national 5-year bands 65–69 … 95–99, 100+. The open group is spread over 100–110 (`nlast = 11`) and collapsed back into 100+. Checked against observed Eurostat single-age data (LE65 error 0.007 y) | Rizzi et al. 2015; Simon's methods draft 2.5; meeting 10 Sep §2.4 |
| Single-age ANs | Each group's attributable fraction applied to the PCLM single-age deaths, so AN ≤ deaths at every age | Simon's methods draft 2.5, option (b) (confirmed by Daniel 23 Sep) |
| Life tables | Period life tables 65–100+, piecewise-constant hazard; LE65; LI65+ = SD of age at death conditional on reaching 65 | Lloyd 2024 `Code_1.R`, `Code_2.R` (Aburto 2022) |
| Decomposition | Horiuchi (`DemoDecomp`, N = `N_HORIUCHI`; validated outputs use 400, see *Batch runs*) by single age × cause (4 ranges + rest): (i) consecutive years within each branch, summable over periods and ages, plus consecutive 5-year-period means (used for block figures, so single-year weather at block endpoints does not drive the result; Lloyd 2024 also decomposed multi-year average ANs); (ii) with vs without CC on 5-year-period mean rates, where rest contributes 0 by construction | Lloyd 2024; Simon's methods draft 2.7 |
| ERF basis | `bs`, degree 2; knots at the 10/75/90th percentiles and boundaries at the range of the city's **full ERA5-Land series 1990–2019** (the series the ERFs were estimated on) | Masselot 2025 `03_attribution.R` (`tper`) |
| Temperature projections (with climate change) | ISIMIP3BASD trend-preserving bias correction against ERA5 2000–2014, applied **by month × calibration period** (2015–29, 2030–39, …, 2090–99) | Masselot 2025 `03_attribution.R`, `functions/isimip3.R` |
| Without-climate-change counterfactual | Masselot's `demo` series: each 5-year block of the calibrated GCM series is re-mapped with ISIMIP3 onto the calibrated 2010–2014 distribution. Day-to-day weather is kept and the warming trend removed. Option `counterfactual = "era5_cycle"` (observed 2000–2019 repeated) is kept for sensitivity analysis | Masselot 2025 `03_attribution.R`; Simon's methods draft 2.4.3; Simon 23 Sep ("do whatever Masselot did") |
| Attributable number | Daily AN = (1 − 1/RR) × daily deaths, with **RR clamped at ≥ 1**; the ERF is held constant (no adaptation) | Masselot 2025 `03_attribution.R`; Simon's methods draft 2.4.1 |
| MMT | Recomputed with the Masselot 2025 rule: the argmin of the ERF over percentiles 25–99 of the full ERA5 series. The Masselot 2023 value is kept as `mmt_2023`, and the fixture uses it | Masselot 2025 `03_attribution.R` |
| Temperature ranges | Split at the MMT first, then at the 2.5th/97.5th percentiles of ERA5 1990–2019, fixed over time. MMT > p97.5 means all heat is extreme | Lloyd 2024 `09_0_Attr_Number.R`; Simon's methods draft 2.4.2 (**open point:** the draft says 2000–2014) |
| Calendar | 29 February removed from all daily series; 365-day years | Masselot 2025 `01_pkg_params.R` (`dayvec`) |
| GCM ensemble | 19 of the 21 GCMs in `tmeanproj` (CMCC_CM2_SR5 and TaiESM1 excluded; `01` refuses them). IITM_ESM SSP3 2099 is filled with its 2098 series. The pilot uses GFDL_ESM4 only | Masselot 2025 `01_pkg_params.R`, `03_attribution.R` |

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

`data/prep_data.RData` is a derived file built by the temperature mode of `pipeline/00_prep_data.R` from `era5series.gz.parquet` and `city_results.csv`. It is also not committed.

## Environment

Validated on 2026-09-23 with R 4.4.3 and these packages:

| Package | Version |
|---|---|
| data.table | 1.18.6.1 |
| arrow | 25.0.0 |
| dplyr | 1.2.1 |
| dlnm | 2.4.10 |
| ungroup | 1.4.4 |
| DemoDecomp | 1.14.1 |
| ggplot2 | 4.0.3 |
| patchwork | 1.3.2 |

```r
install.packages(c("data.table", "arrow", "dplyr", "dlnm", "ungroup", "DemoDecomp", "ggplot2", "patchwork"))
```

Runtime for one city, one GCM and one SSP on 2 cores: steps 00–03 take about 2 minutes, and 04 about 30 minutes with the script defaults (N = 400, `DECOMP_ANNUAL=1`). Step 04 caches every Horiuchi step, so an interrupted run resumes.

## Legacy code

An earlier implementation (`notebook/`, `scripts/`, `R/`), which used EUROPOP2019/Eurostat demography and produced the July–August Europe runs, has been removed; it is superseded by `pipeline/`. It is kept in git history under the tag `legacy-v1` (`git checkout legacy-v1 -- scripts` restores it).

## References

- Gasparrini A, Leone M. Attributable risk from distributed lag models. *BMC Med Res Methodol* 14:55, 2014. doi:10.1186/1471-2288-14-55
- Masselot P et al. Excess mortality attributed to heat and cold: a health impact assessment study in 854 cities in Europe. *Lancet Planet Health* 7:e271–e281, 2023. doi:10.1016/S2542-5196(23)00023-2
- Masselot P et al. Estimating future heat-related and cold-related mortality under climate change, demographic and adaptation scenarios in 854 European cities. *Nat Med* 31:1294–1302, 2025.
- Lloyd SJ et al. The reciprocal relation between rising longevity and temperature-related mortality risk in older people, Spain 1980–2018. *Environ Int* 193:109050, 2024. doi:10.1016/j.envint.2024.109050
- Aburto JM et al. Significant impacts of the COVID-19 pandemic on race/ethnic differences in US mortality. *PNAS* 119:e2205813119, 2022.
- Horiuchi S, Wilmoth JR, Pletcher SD. A decomposition method based on a model of continuous change. *Demography* 45:785–801, 2008.
- Rizzi S, Gampe J, Eilers PHC. Efficient estimation of smooth distributions from coarsely grouped data. *Am J Epidemiol* 182:138–147, 2015.
- Pascariu MD et al. ungroup: An R package for efficient estimation of smooth distributions from coarsely binned data. *JOSS* 3:937, 2018.
