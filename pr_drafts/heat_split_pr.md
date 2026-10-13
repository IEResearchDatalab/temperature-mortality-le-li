# Draft PR: WIP — three-way heat split blocked by 2015 threshold contradiction

## Summary

This branch contains the audited implementation scaffold for the 15-year moving extreme-heat split, the six-cause downstream contract, prep-cache fingerprinting, and alternate-root figure handling. It is intentionally marked WIP and must not be merged.

## Blocking test

Madrid/GFDL-ESM4/SSP3 shows that Masselot's calibrated historical GCM 2000–2014 P97.5 (28.2899418625876°C) does not equal the raw ERA5-Land static 2000–2014 P97.5 (28.201°C). The 0.0889418625876°C difference violates the requested numerical-equality check. Using ERA5 for the moving window or adding an anchor would contradict another explicit requirement unless the method is clarified.

## Completed checks

- R and shell syntax pass.
- Masselot 2023 fixture: maximum relative error 4.56e-05 (tolerance 0.001).
- Complete 15-year windows; no-CC NormExtrHeat AN exactly zero; one moving-threshold clamp in 85 evaluated years.
- Canonical outputs untouched; no full batch started.

## Review request

Please decide whether the moving historical window remains the calibrated GCM series without exact equality, whether 2015 uses ERA5, or whether a documented calibration anchor is added. See `agent_log/20261013_final_report.md`.
