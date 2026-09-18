# Validation contract

Every stage, whether production or validation, must write:

1. `checks.csv`
2. `failures.csv`
3. `summary.md`
4. `overview.png`

Rules:
- `checks.csv` lists every hard check, its threshold, and the observed value.
- `failures.csv` lists each failed hard check row, with enough keys to reproduce it.
- `summary.md` states the sample, units, thresholds, and verdict.
- `overview.png` must be a compact visual summary that can be inspected without reading code.
- Any hard failure stops execution immediately.

This contract is required for all five production scripts and all validation stages.
