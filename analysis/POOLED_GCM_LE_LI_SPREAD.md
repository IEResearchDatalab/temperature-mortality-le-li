# Pooled per-GCM LE65/LI65 spread

This analysis quantifies climate-model spread in pooled European and regional
LE65 and LI65 under SSP3-7.0. It is an analysis-layer extension and does not
modify the five-script production pipeline.

## Method

For every one of 854 cities and 19 GCMs, the script reads the production
`01_attribution_grouped.csv`. It reconstructs single-age attributable deaths
with the same city/year/age-band death weights used by pipeline Part 02. It then
sums population and death counts for Europe and the four Masselot regions before
forming mortality rates or life tables. LE65 and LI65 use the unchanged
Lloyd/Aburto functions in pipeline Part 04.

The without-climate-change temperature deaths define the non-temperature
`rest` schedule, as in pipeline Part 03. The same `rest` schedule is combined
with each GCM's with- and without-climate-change attributable deaths.

Outputs include annual and five-year-period levels, with-minus-without climate
effects, change from the 2020-2024 climate effect, and empirical summaries
across the 19 GCMs. These are climate-model spreads, not confidence intervals.

## Run

```bash
cd /home/SHARED/temperature-mortality-le-li-vig-phase2-20261007
NCORES=8 Rscript analysis/pooled_gcm_le_li_spread.R \
  /home/SHARED/temperature-mortality-le-li \
  results/vig_phase2_20261007/pooled_gcm_le_li_spread \
  3
```

Generated outputs remain under the gitignored `results/` tree. The output
directory contains an input inventory, run manifest and validation table.
It also contains meeting-ready figures for the end-century LE65 spread, the
European LE65 trajectory and the end-century LI65 spread.

## Required validation

- complete 854-city × 19-GCM input coverage;
- complete grouped-AN domains and finite values;
- death weights sum to one and reproduce grouped deaths;
- reconstructed single-age ANs reproduce every grouped AN;
- pooled populations are positive and `rest` mortality is nonnegative;
- annual and period output domains are complete and finite;
- averaging pooled cause counts over GCMs reproduces the existing production
  ENSEMBLE pooled counts and LE65/LI65 levels;
- all empirical quantiles are ordered and based on 19 GCMs.

The exact validated values and code revision are written by the script to the
output manifest and `validation_checks.csv`.

## Validated SSP3 result (7 October 2026)

The full run read 16,226 grouped-AN files (854 cities × 19 GCMs) and produced
15,200 annual and 3,040 five-year-period level rows. All 16 checks passed.
Maximum errors were `3.638e-12` for grouped-AN reconstruction, `7.458e-11`
for production ENSEMBLE pooled-count reconstruction, `8.882e-14` years for
ENSEMBLE LE65 and `1.510e-14` SD-years for ENSEMBLE LI65.

For Europe, the change in the LE65 climate effect from 2020-2024 to 2095-2099
was negative in all 19 GCMs: mean `-0.2054` years, with empirical 2.5th and
97.5th percentiles of `-0.4197` and `-0.0244` years. Southern Europe was also
negative in all 19 GCMs (mean `-0.4214` years). Northern Europe had a near-zero
mean (`+0.00385` years) and a spread crossing zero.

## Companion attributable-death spread

`analysis/gcm_spread_attributable_deaths.R` is the previously validated
central-ERF mortality diagnostic, now tracked on the VIG branch. Run it with:

```bash
NCORES=8 Rscript analysis/gcm_spread_attributable_deaths.R \
  . \
  /home/SHARED/temperature-mortality-le-li/results/europe/ssp3 \
  results/vig_phase2_20261007/gcm_spread_attributable_deaths
```

Its run manifest, input inventory and four coverage/finite-value checks are
written beside the CSVs and two figures. Generated outputs remain gitignored.
