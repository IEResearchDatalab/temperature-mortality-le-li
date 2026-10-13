# Simon v3 figures and table: audited specification

Audit date: 2026-10-13. Source material: `references/simon_v3/`, Simon's v3 skeleton, the 13 October decisions, and canonical SSP3 Objects 1–3. All outputs are central estimates without uncertainty intervals.

## Available data and conventions

Canonical SSP3 Object 1 has 24,595,200 rows for 854 cities, 2020–2099, two scenarios (`with_cc`, `without_cc`), ages 65–100+, and five current causes (`ExtrCold`, `ModCold`, `ModHeat`, `ExtrHeat`, `rest`). Object 2 has 7,071,120 age/cause decomposition rows and Object 3 has 2,732,800 LE65/LI65 level rows for 19 GCMs plus `ENSEMBLE`. The 35 pooled geographies are Europe, four regions, and 30 countries. The new run adds `NormExtrHeat`, making six causes.

`rest` is derivable in Part 03 as all-cause deaths minus without-climate-change temperature-attributable deaths at each age/year. It is held identical between branches; the with-CC total is altered only by the difference between with- and without-CC temperature deaths. Part 04 decomposes the resulting six cause-specific mortality rates.

## Figure 1 — absolute risk by age and temperature range

Layout: three SSP rows by two period columns (2055–59 and 2095–99). Within every age group, the upper bar is with climate change and the lower pale bar is without climate change. Cold extends left and heat right. Five plotted ranges are ExtrCold, ModCold, ModHeat, NormExtrHeat, and ExtrHeat.

Status: derivable after the rerun from Object 1. Aggregate single ages to 65–69, 70–74, 75–79, 80–84, 85–89, 90–94, and 95+. For each SSP/scenario/period/geography/age group/cause, compute `100000 * mean annual attributable deaths / mean annual population`; equivalently, sum deaths over the five years, divide by the sum of the five annual populations, and multiply by 100,000. Units are attributable deaths per 100,000 person-years. NormExtrHeat is structurally zero without climate change.

## Table 1 — SSP3 end-century geography summary

Rows: Europe, four regions, and 30 countries. Period: 2095–99 unless stated otherwise.

| Column | Status and formula | Unit |
|---|---|---|
| Geography | available from pooled geography IDs/labels | text |
| Temperature-attributable deaths | derivable from Object 1, five temperature causes, with CC | deaths per 100,000 aged 65+ person-years |
| Climate-change-attributable share | derivable as `(temperature deaths with CC - temperature deaths without CC) / temperature deaths with CC` | percent |
| ModHeat share of heat deaths | derivable as `ModHeat / (ModHeat + NormExtrHeat + ExtrHeat)` with CC | percent |
| NormExtrHeat share of heat deaths | derivable analogously | percent |
| Heat-only LE65 erosion ratio | derivable by recomputing the requested heat-only life-table comparison from existing cause-specific deaths and populations, or marked unavailable if it would require a new estimand | ratio |
| All-temperature LE65 erosion ratio | derivable from Object 3 gains: `(gain_without_cc - gain_with_cc) / gain_without_cc`, 2020–24 to 2095–99 | ratio |

Ratios are displayed as `n/a` when the without-CC gain is below one month; affected rows go in the footnote.

## Figure 2 — LE65 contributions

Layout: SSP1 top, SSP3 middle, SSP2 bottom. Panel a is a stacked time series of within-branch changes between consecutive 5-year-period means, converted to months per decade; the stack sums to a black net diamond. Panels b and c show contributions by 5-year age group for 2020–24 to 2055–59 and 2055–59 to 2095–99. Units are months of LE65 per decade. Status: available/derivable from Object 2 rows with `type=within` and `scenario=with_cc`; longer intervals are sums of consecutive 5-year-period decompositions. Panel d is optional and is omitted unless a transparent index can be made solely from existing Object 1/2/3 columns. A Southern Europe version uses the pooled regional outputs.

## Figure 3 — LI65+ contributions and LE/LI plane

Left: six panels (three SSPs by two periods), bars by 5-year age group, cold = ExtrCold + ModCold and heat = ModHeat + NormExtrHeat + ExtrHeat. Units are days of LI65+ per decade. The crossover age is the first age/group at which the cold contribution changes sign; ages above it are shaded. The all-age net is printed in each panel. Status: derivable from Object 2 within-branch contributions.

Right: cumulative contributions since 2020–24, x = months of LE65 and y = days of LI65+. Each SSP has one cold and one heat path, a dot at 2055–59, and an arrowhead at 2095–99. The grey social-progress direction is the normalized direction of the `rest` contribution only. Status: derivable by cumulative sums of Object 2 within-branch contributions. Quadrant labels follow the mock.

## Labels and captions

Every SSP1 output says “survival ratios patched (WCDE rounding)”; every SSP2 output says “including 12 Swiss cities”. Captions state the static ERA5-Land 2000–2014 baseline, the preceding 15-year moving P97.5 window, 5-year-period-mean decompositions, units, and “central estimate, no uncertainty intervals”. Annual decompositions are optional (`DECOMP_ANNUAL=1`) and are not used in the paper figures.

## Follow-up issue state at audit

1. Canonical SSP3 has no contemporaneous `run_manifest.txt`; provenance remains retrospective/incomplete.
3. Canonical SSP1 `failed.txt` was missing at audit and was restored byte-for-byte from the documented sensitivity seed (MD5 `5a598393ebbcf41ac38fca80dbecc57d`).
5. The shared `data/prep_data.RData` had no fingerprint validation and was reused by path alone.
7. `analysis/run_city_figures.sh` hard-coded `results/europe/ssp<k>` and ignored `BATCH_ROOT_OVERRIDE`.
