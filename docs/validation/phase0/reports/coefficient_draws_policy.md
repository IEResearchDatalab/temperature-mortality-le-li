# Phase 0 — Coefficient-draw policy

## Verdict
PASS

- Published coefficient-draw file: data/coef_simu.csv
- Rows: 4270000
- Unique draws: 1000
- Max draw index: 1000
- File size (bytes): 491247391

## Decision
- Use the published 1000 coefficient draws by default.
- They are large but still tractable in the audited pipeline and preserve the intended uncertainty representation.
- The deterministic central-coefficient path remains available for historical validation gates.
