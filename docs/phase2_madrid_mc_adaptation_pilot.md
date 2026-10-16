# Madrid ERF Monte Carlo and adaptation pilot

## Purpose

This is an isolated Phase 2 pilot for Madrid (`ES001C`) under SSP3-7.0 and
GFDL-ESM4. It advances two commitments without modifying the five-step
production pipeline:

1. propagate the Masselot exposure-response coefficient draws into LE65 and
   LI65; and
2. demonstrate the published heat-only excess-relative-risk attenuation rule
   at 0%, 10%, 50% and 90%.

The pilot is an experiment conditional on one city, one SSP and one GCM. Its
empirical intervals quantify ERF coefficient uncertainty only. They do not
include climate-model spread, demographic uncertainty or PCLM structural
uncertainty.

## Reused production methods

- Temperature calibration and the without-climate-change `demo`
  counterfactual follow `pipeline/01_attribution.R`.
- Daily attributable deaths use `AN = (1 - 1/RR) x daily deaths`, with
  `RR >= 1` and leap days removed.
- The simulation coefficients come from Masselot's `coef_simu.csv`. All 1,000
  available draws are used; the original article commonly reports 500.
- Adaptation follows Masselot's published code:
  `RR_adapted = 1 + (RR - 1) x (1 - attenuation)` on heat days only.
- Grouped attributable deaths are allocated to single ages using the same
  annual death weights as `pipeline/02_single_age.R`.
- LE65 and LI65 use the unchanged Lloyd/Aburto functions in
  `pipeline/04_le_li_decomposition.R`.

The headline estimand is the change in the with-minus-without climate effect
between 2020-2024 and 2095-2099. For LE65 it is equivalent to the climate-driven
change in the projected LE65 gain.

## Command

```bash
cd /home/SHARED/temperature-mortality-le-li-vig-phase2-20261007
CITY_ID=ES001C SSP=3 GCM=GFDL_ESM4 \
SOURCE_ROOT=/home/SHARED/temperature-mortality-le-li \
CANONICAL_ROOT=/home/SHARED/temperature-mortality-le-li/results/europe \
OUTPUT_DIR=results/europe/auxiliary/vig_phase2_20261007_madrid_mc_adaptation \
N_DRAWS=0 DRAW_CHUNK_SIZE=100 \
Rscript analysis/madrid_mc_adaptation_pilot.R
```

## Headline results

The values below are empirical ERF coefficient-draw results for one GCM.

| Heat excess-RR attenuation | Mean change in LE65 climate effect | Empirical 2.5%-97.5% | Draws below zero |
|---:|---:|---:|---:|
| 0% | -4.57 months | -5.68 to -3.20 | 100.0% |
| 10% | -4.33 months | -5.45 to -2.96 | 100.0% |
| 50% | -3.00 months | -4.14 to -1.72 | 100.0% |
| 90% | -0.28 months | -1.04 to +0.39 | 77.3% |

The 0%, 10% and 50% settings remain negative in all 1,000 coefficient draws.
At 90% attenuation, the empirical interval crosses zero. This is a useful
engineering demonstration of the adaptation machinery, not evidence that 90%
attenuation is achievable or a preferred VIG assumption.

For LI65, the mean change is -0.081 SD-months without attenuation, with an
empirical interval of -0.170 to -0.005 SD-months. The LI65 result is much
smaller and becomes inconclusive under stronger attenuation; it should remain
secondary to LE65.

## Validation

All ten checks in `validation_checks.csv` pass:

- 3,000/3,000 Madrid age-group/draw rows are present (1,000 draws x 3 age
  groups), with no duplicate keys or non-finite coefficients.
- The mean draw coefficients reproduce the central coefficients closely:
  maximum absolute difference 0.004156; minimum within-coefficient correlation
  across the three age groups 1.000000.
- All daily-to-annual attributable deaths are finite and non-negative.
- Single-age death-allocation weights sum to one within `1.11e-16`.
- The independent central/no-adaptation calculation reproduces the existing
  GFDL production LE65 and LI65 levels with a maximum absolute error of
  `6.04e-14`.
- All life-table outputs are finite and within broad plausibility bounds; the
  minimum climate-adjusted single-age death count is 57.248.
- LE65 improves monotonically with greater heat-risk attenuation for the
  central estimate and all 1,000 simulation draws.
- Every central estimand lies inside its corresponding empirical 95% interval.

The spline package emits the same warning as the production pipeline when
future temperatures fall beyond historical boundary knots. This is expected
tail extrapolation, not a failed check; it is precisely why ERF uncertainty
widens under late-century heat.

## Outputs

- `pilot_draw_estimands.csv` — one row per coefficient draw and attenuation
  setting, including baseline and late-century LE65/LI65 contrasts.
- `pilot_summary.csv` — mean, standard deviation, empirical limits, sign shares
  and central values.
- `central_annual_levels.csv` — annual central levels for the ten comparison
  years.
- `coefficient_validation.csv` — central-versus-draw-mean coefficient audit.
- `production_reproduction.csv` — exact comparison with existing production
  levels.
- `validation_checks.csv` — acceptance checks.
- `01_le65_mc_adaptation.png` and `02_li65_mc_adaptation.png` — meeting-ready
  figures.
- `provenance.txt` — server, repository, branch, commit, command, inputs,
  outputs, interpretation and limitations.

## Limitations and safe wording

Safe: “In a Madrid/GFDL-ESM4 pilot, the late-century climate penalty to the
LE65 gain remains negative across all 1,000 ERF draws under 0%-50% heat-risk
attenuation; the interval crosses zero at 90% attenuation.”

Do not call the interval a combined confidence interval, generalize it to
Europe, or present the attenuation values as forecasts. A production Europe
analysis must cross ERF draws with all 19 GCMs, pool counts before constructing
life tables, and preserve draw identifiers across cities.
