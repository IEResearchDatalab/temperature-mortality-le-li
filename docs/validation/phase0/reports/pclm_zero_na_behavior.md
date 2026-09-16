# Phase 0 — PCLM zero/NA behavior

## Verdict
BLOCKED

## Current helper behavior
- Positive inputs return finite nonnegative single-age weights and preserve the group sum.
- All-zero inputs return a zero vector without fitting PCLM.
- All-NA inputs stop with `Input contains only NA values.`
- Mixed NA inputs stop with the underlying PCLM error `'y' contains NA values`.
- Signed inputs are handled only by the signed wrapper, which splits positive and negative parts before PCLM and then recombines.

## Desired invariants from the Phase 0 request
- grouped deaths > 0 and finite: derive nonnegative weights
- grouped deaths = 0 and grouped AN = 0: emit zero single-age values without fitting PCLM
- grouped deaths = 0 and grouped AN != 0: hard failure because inputs are inconsistent
- missing grouped deaths: propagate NA, record the complete cell key, and apply the existing coverage stop rule
- missing AN with available deaths: preserve NA and record it
- PCLM nonconvergence or invalid weights: hard failure
- weights must be finite, nonnegative, and sum to one within 1e-12
- grouped AN reconstruction must be within 1e-10

## Conclusion
- The current helper already satisfies the positive and zero cases, but it does **not yet** satisfy the explicit NA-propagation behavior requested for production.
- Because of that, PCLM zero/NA behavior remains a real blocker.

## Status
BLOCKED
