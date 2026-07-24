# Session Diary — 2026-07-24

## Task
Exactly reproduce Masselot's baseline ANs, disaggregate into 4 temperature ranges, ensure aggregate matches. Respond to PI feedback.

## Progress

### ✅ Completed

1. **ERA5 data verification** — Compared our `obs_data` against Masselot's `era5series.csv`. Values are **100% identical** (9,357,278 rows).

2. **Point estimate hypothesis** — RULED OUT. Using point estimate alone still gives ~2.8%/5.4% error.

3. **Found exact 2023 formula** in `references/2023-masselot-excess/09_ResultsCityAge.R`:
   ```r
   afday <- (1 - exp(-bvarcen %*% icoefs$fit))
   anday <- afday * cityageres[i, "death"]
   anlist <- sum(anday) / length(era5)
   ```
   Using full ERA5 period (1990-2019), point estimate coefficients, and actual day count.

4. **Verified exact replication** on 10 cities (50 city-age combos): **0.0000% error**. Script: `temp_scripts/verify_exact_replication.R`.

5. **Updated `03_attribution.R`** with 2023 formula:
   - `af = 1 - exp(-logRR)` instead of `rr = exp(logRR); an = (1-1/rr)*death/365`
   - `an = af * death_annual` with normalization by `days_per_year`
   - `Bound = range(ERA5)` in both basis calls
   - Exported `varfun`, `vardegree` to cluster workers

6. **Tested updated formula on AT001C** with GCM data: still gives ~2%/9% error. **The remaining error is inherent to GCM-vs-ERA5 data**, not the formula. Bias-corrected GCM temperatures differ from observed ERA5 temperatures for the historical period.

7. **Fixed Script 06** (LE/LI extraction):
   - Fixed data.table scoping in parallel workers (`decade` not found → explicit vector refs)
   - Fixed 500x AN inflation bug (`sum(an)/10` → `mean(an)` over sims then years)

8. **Fixed Script 07** (pclm):
   - `pclm(...)$y` → `pclm(...)$fitted` (the `ungroup` package stores results in `$fitted`, not `$y`)
   - Removed dead code (lines 249-273)

9. **Received and read PI email** (`references/email_to_postdoc_corrections_2026-07-24.md`)

### Key PI Feedback

The PI raised 5 critical issues:
1. Formula not yet in production pipeline (was only in test script) → **Fixed today**
2. 4-range disaggregation implementation → **Already exists** (lines 164-170)
3. Missing total mortality data for LE/LI → **Blocking issue** (need to discuss)
4. Demographic framework not designed → **Need methods doc before coding**
5. pclm() NA → **Fixed** ($y → $fitted)

### 🚨 KEY INSIGHT: GCM vs ERA5

The ~2% cold / ~9% heat error that persists after the formula fix is NOT a formula issue. It's because:
- `cityage.csv` (the target) was computed using **ERA5 observed temperatures**
- Our `03_attribution.R` uses **GCM projection temperatures** (bias-corrected)
- Even after ISIMIP3 bias correction, GCM ≠ ERA5

**To get 0% error, the baseline must be computed using ERA5 data directly.**

The 2023 and 2025 pipelines serve different purposes:
- **2023 (baseline)**: ERA5 temperatures → compute observed ANs → `cityage.csv`
- **2025 (projections)**: GCM temperatures → compute future ANs → RDS files

Our project needs BOTH:
- Phase 1: ERA5-based baseline with 4-range disaggregation (matching cityage.csv exactly)
- Phase 2: GCM-based future projections with 4-range disaggregation

### Config/Code Changes Made Today
- `R_pipeline/03_attribution.R`: Updated formula (2023 style: `af * death_annual / days_per_year`)
- `R_pipeline/06_le_li_extraction.R`: Fixed scoping bug and 500x inflation
- `R_pipeline/07_le_li_decomposition.R`: Fixed `$y` → `$fitted`, removed dead code
- `temp_results/` → backed up to `temp_results_backup_wrong_formula_2024-07-24/`

### 📁 Files Created
- `temp_scripts/test_new_formula.R` — tested updated formula on AT001C
- `temp_scripts/test_le_function.R` — verified life expectancy calculation

### ❓ Questions for PI

1. **Total mortality data**: The `death_baseline` column in `city_results.csv` contains annual deaths by city-age. Is this sufficient for life table construction, or do we need separate mortality data?
2. **Scope**: Should we aim for baseline-only (2000-2019) first? The future projections infrastructure is complex.
3. **GCM vs ERA5**: For the baseline, should we use ERA5 directly (giving exact match) and only use GCM for future projections?