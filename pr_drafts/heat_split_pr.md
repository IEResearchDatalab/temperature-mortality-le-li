# Draft PR: three-way heat split with anchored 15-year threshold

## Summary

This branch implements the 15-year moving extreme-heat split, the six-cause downstream contract, prep-cache fingerprinting, and alternate-root figure handling. It remains a draft while the full SSP rerun and paper outputs are produced and checked.

## Threshold decision

The moving threshold uses an observed-baseline-plus-model-change anchor: `ERA5 P97.5(2000–2014) + [calibrated GCM P97.5(t-15:t-1) - calibrated GCM P97.5(2000–2014)]`. This follows Masselot's calibration principle, preserves the calibrated GCM trajectory, removes residual model-specific baseline offsets, and makes the 2015 threshold equal the static ERA5 threshold by construction.

## Madrid checks

- R and shell syntax pass.
- Masselot 2023 fixture: maximum relative error `4.56e-05` (tolerance `0.001`).
- Complete 15-year windows; no-CC NormExtrHeat AN exactly zero; two anchored-threshold clamps in 85 evaluated years.
- Old-versus-new heat/cold totals agree within `1.10e-11`.
- Five-year and annual N=400 decompositions pass all six-cause closure checks; annual runtime was 19.8 minutes.
- Canonical outputs remain untouched; the full rerun writes only below `results/europe/heatsplit_20261013/`.

## Remaining work

Complete SSP3, SSP2 and SSP1; verify pooled outputs; produce the paper figures, table, provisional PDF and methods text; then replace this draft status with the final results summary.
