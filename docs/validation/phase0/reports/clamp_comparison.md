# Phase 0 — Clamp comparison

## Verdict
PASS

## Fixture
- City: AT001C
- Age group: 20-44
- Historical observations: `data/prep_data.RData` (`obs_data$tmean_obs`)
- Baseline deaths: fixed Masselot `cityage.csv$death` for the selected city-age row
- Coefficients: central `data/coefs.csv` row for the same city-age
- Aggregation: daily AF contributions summed by year and temperature range, then annualized by days in year

## Numeric comparison against Masselot
   metric unclamped  clamped masselot abs_diff_unclamped abs_diff_clamped
   <char>     <num>    <num>    <num>              <num>            <num>
1:  total  8.915557 8.986107 8.915087       0.0004704293     0.0710198769
2:   cold  4.411874 4.482424 4.411663       0.0002110726     0.0707605202
3:   heat  4.503683 4.503683 4.503423       0.0002593567     0.0002593567
   rel_diff_unclamped rel_diff_clamped    better
                <num>            <num>    <char>
1:       5.276778e-05      0.007966257 unclamped
2:       4.784423e-05      0.016039420 unclamped
3:       5.759100e-05      0.000057591 unclamped

## Behavior
- Negative daily AF count: 513 (min AF = -0.268650)
- Raw daily sum before annualization: 2627.867
- Annualized 1990 total: 7.200
- Raw/annualized ratio: 365.0
- Mean annual total unclamped: 8.915557
- Mean annual total clamped: 8.986107

## Interpretation
- The unclamped calculation is closer to the published Masselot values for total and cold attributable numbers; heat is a tie at the displayed precision.
- Clamping removes negative daily AF values and shifts the annual totals upward; in this fixture it worsens agreement with the reference.
- Therefore the historical validation gate supports the unclamped AF behavior for this exact Masselot city-age fixture.

## Status
PASS
