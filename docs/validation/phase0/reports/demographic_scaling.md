# Phase 0 — Demographic scaling

## Verdict
PASS

## What was checked
- Source table: `results/projdata/projdata_prototype.csv` (230580 rows, 854 cities)
- Anchor year used for the numeric check: 2015
- Comparison year used for the numeric check: 2095

## Correct scaling rule
The observed city projections in the published prototype behave as a simple country-agegroup growth rescaling:

```text
growth(t) = national_projection(t) / national_projection(anchor_year)
city_value(t) = city_value(anchor_year) * growth(t)
```

In the prototype, the same country-agegroup growth factor is applied to every city in that country-agegroup cell.

## Numeric results
- Max absolute change in city share within a country-agegroup-SSP cell from 2015 to 2095: 1.665e-15
- Max absolute change in sample coverage fraction `sum(city values)/national projection`: 1.998e-15
- Max absolute discrepancy between city growth and country growth: 8.171e-14
- Max absolute discrepancy between country growth and `wittpop` growth: 5.684e-14
- Max span in coverage fraction across years within a country-agegroup-SSP cell: 4.330e-15

Representative cells with the largest observed share differences (all numerical noise):
   URAU_CODE CNTR_CODE agegroup    ssp   share_diff coverage_diff city_growth
      <char>    <char>   <char> <char>        <num>         <num>       <num>
1:    DK001C        DK    45-64      2 1.665335e-15  5.551115e-16   1.0605053
2:    CZ001C        CZ    65-74      3 1.387779e-15  7.216450e-16   0.8322096
3:    HU001C        HU    75-84      1 1.332268e-15  8.881784e-16   1.5561939
4:    SI001C        SI    65-74      3 1.332268e-15  4.718448e-16   0.9781726
5:    SI002C        SI    65-74      3 1.276756e-15  4.718448e-16   0.9781726
   country_growth witt_growth
            <num>       <num>
1:      1.0605053   1.0605053
2:      0.8322096   0.8322096
3:      1.5561939   1.5561939
4:      0.9781726   0.9781726
5:      0.9781726   0.9781726

## Interpretation
- City baseline shares are preserved over time to numerical precision.
- The sample coverage fraction is also preserved; the method does **not** force sampled-city totals to equal full national totals.
- This is consistent with the intended growth-ratio rescaling and with the existing prototype outputs.

## Unsupported scientific assumptions
1. Within each country-agegroup-SSP cell, all cities share the same demographic growth factor.
2. The city age structure within the baseline sample is held fixed apart from the multiplicative growth factor.
3. Any uncovered national population outside the Masselot urban sample is intentionally left out; the approach preserves sample coverage but does not reconstruct full national totals from the city sample.

## Status
PASS
