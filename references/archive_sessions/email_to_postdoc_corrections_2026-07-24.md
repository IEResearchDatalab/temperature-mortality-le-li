# Email: Critical Corrections and Path Forward

**To:** Post-doc Researcher  
**From:** PI (Demography Lead)  
**Date:** 2026-07-24  
**Subject:** Research Status Review – Critical Issues and Immediate Action Plan

---

Dear [Post-doc],

I've reviewed your session diaries (July 2, 23, 24), the project documentation, and the Masselot 2023 paper. First, let me acknowledge the **excellent detective work** you've done identifying the basis function mismatch (ns → bs) and discovering the formula discrepancy between the 2023 and 2025 code. That level of methodological rigor is exactly what we need.

However, after careful review, I've identified **critical gaps** that need immediate attention before we can proceed. I'm writing to clarify the current status and provide a corrected action plan.

---

## 🚨 **CRITICAL ISSUE #1: The Pipeline Still Uses the Wrong Formula**

Your July 24 diary states:

> "Verified exact replication across 10 random cities (50 city-age combos): **100% match**"

However, this validation was done in a **temporary test script** (`verify_exact_replication.R`), not in the production pipeline. 

### Current Status:
- ❌ `R_pipeline/03_attribution.R` **still uses the 2025 formula** (wrong)
- ❌ The 854 RDS files in `temp_results/` were generated with the **wrong formula** (July 23 run)
- ❌ All downstream results (aggregation, validation) are based on **incorrect data**

### What Needs to Happen:
1. **Update `03_attribution.R`** to use the correct 2023 formula:
   ```r
   # CORRECT (2023):
   afday <- (1 - exp(-bvarcen %*% icoefs$fit))
   anday <- afday * cityageres[i, "death"]
   anlist <- sum(anday) / length(era5)  # divide by TOTAL DAYS including leap years
   
   # NOT:
   # an <- (1 - 1/rr) * daily_deaths
   # an_annual <- sum(an) / n_years
   ```

2. **Back up the current (wrong) `temp_results/`** folder
3. **Re-run Script 03 for all 854 cities** with the corrected code
4. **Re-run Scripts 04 (aggregation) and 09 (validation)** with the new results
5. **Verify**: Validation should show ~0% error (not 1.7% cold / 6.8% heat)

**Action:** Do NOT proceed with any other work until this is complete. All current results are invalid.

---

## 🚨 **CRITICAL ISSUE #2: The 4-Range Disaggregation Has Not Been Implemented**

Your task (per `task_understanding.md`) is to disaggregate temperature-attributed mortality into:

| Range | Definition |
|-------|------------|
| **Extreme Cold** | T < P₂.₅ |
| **Moderate Cold** | P₂.₅ ≤ T < MMT |
| **Moderate Heat** | MMT ≤ T ≤ P₉₇.₅ |
| **Extreme Heat** | T > P₉₇.₅ |

### Current Status:
- ❌ The current code only computes **total cold** and **total heat**
- ❌ No logic exists to classify temperatures into 4 ranges
- ❌ No code to compute separate AN values for each range

### What Needs to Happen:
Modify `03_attribution.R` to:

1. **Classify each temperature observation** into one of 4 ranges:
   ```r
   temp_range <- case_when(
     tmean < p2.5 ~ "ExtrCold",
     tmean < mmt ~ "ModCold",
     tmean <= p97.5 ~ "ModHeat",
     tmean > p97.5 ~ "ExtrHeat"
   )
   ```

2. **Compute 4 separate AN values** (not just cold/heat):
   - Loop over the 4 ranges separately
   - Apply the attribution formula to each subset
   - Store results as: `an_extrcold`, `an_modcold`, `an_modheat`, `an_extrheat`

3. **Validate additivity**:
   ```r
   # These should match within floating-point tolerance:
   an_extrcold + an_modcold ≈ excess_cold_est
   an_modheat + an_extrheat ≈ excess_heat_est
   ```

**This is the core innovation of your project.** It cannot be skipped.

---

## 🚨 **CRITICAL ISSUE #3: Missing Total Mortality Data**

For life expectancy and lifespan inequality analysis, you need to compute:

$$\text{LE}_{\text{counterfactual}} = f(\text{total mortality} - \text{temperature-attributed mortality})$$

### Current Status:
- ✅ You have temperature-attributed deaths (once Script 03 is fixed)
- ❌ You do **NOT** have total all-cause mortality counts by city, age, year

### What Needs to Happen:
1. **Identify the data source**:
   - Do we have daily mortality counts for 854 cities (2000-2019)?
   - Or annual deaths by 5-year age groups (minimum requirement)?
   - Source: Eurostat NUTS3? MCC network? National statistical offices?

2. **Acquire and validate**:
   - Total deaths by city, age, year
   - Must cover the same 854 cities and age ranges as the temperature data
   - Validate: Do city-level deaths match Eurostat vital statistics?

3. **Compute counterfactual mortality**:
   ```r
   deaths_counterfactual <- deaths_observed - deaths_temp_attributed
   ```

**Action:** We need to discuss data availability before you proceed further. This may be a project-blocking issue.

---

## 🚨 **CRITICAL ISSUE #4: Demographic Framework Not Designed**

The paper skeleton mentions:
- Aburto/Lloyd methodology for LE/LI decomposition
- Testing the "perverse equalizer" hypothesis
- Geographic aggregation (city → regional → European)

### Current Status:
- ❌ No code implements life table construction
- ❌ No code implements Aburto decomposition
- ❌ No framework for computing lifespan inequality (Gini? SD? e-dagger?)
- ❌ No design for testing the "perverse equalizer" hypothesis

### What Needs to Happen:
Before writing code, we need a **methods document** that specifies:

1. **Life table construction**:
   - Period or cohort? (I recommend period for baseline, cohort for projections)
   - Radix (100,000?)
   - Age range (0-100+ or 65-100+?)
   - How to handle open age interval (85+ or 100+?)

2. **Cause-removed life tables**:
   - Method for removing temperature-attributed deaths
   - Do we assume independence, or use Chiang's competing risk method?

3. **Lifespan inequality metric**:
   - Which measure? (Gini coefficient, SD of age at death, e-dagger, IQR?)
   - Computed from $l_x$, $d_x$, or directly from observed deaths?

4. **Decomposition framework**:
   - Decompose ΔLE and ΔLI by:
     - Age (single years: 65, 66, 67, ...)
     - Temperature range (ExtrCold, ModCold, ModHeat, ExtrHeat)
     - Geography (city, region, country)
   - Aburto method for continuous age? Or stepwise discrete?

**Action:** Draft a 2-3 page "Demographic Methods" document for my review **before coding**. Include equations and pseudocode.

---

## 🚨 **CRITICAL ISSUE #5: pclm() Returning NA**

Script 07 is blocked because `pclm()` returns NA when trying to disaggregate 5-year age groups into single years.

### Likely Causes:
1. **Negative AN values** in some city-age-scenario combinations
   - The RR clamping fix (July 2) was supposed to prevent this
   - Re-check after applying the correct 2023 formula
   
2. **Zero or missing death counts** in some age groups
   - Small cities might have no deaths in youngest age groups (20-24)
   
3. **Improper pclm specification**:
   - Are you passing the correct offset (exposure)?
   - Are age intervals correctly specified?

4. **Trying to disaggregate negative numbers**:
   - pclm assumes non-negative counts
   - If some ANs are negative (e.g., cold "saves" lives in moderate range?), pclm will fail

### What Needs to Happen:
1. **After fixing Script 03**, check for negative ANs:
   ```r
   summary(an_extrcold)  # should be all >= 0
   summary(an_modcold)   # might have negatives if RR < 1 in this range
   ```

2. **If negatives exist**, decide on handling:
   - Set negatives to zero? (loses information)
   - Disaggregate positive and negative separately, then recombine?
   - Use absolute values, then reapply sign?

3. **Test pclm on simple case first**:
   - Try disaggregating **observed total deaths** (not ANs) for one city
   - If this works, the issue is with the AN data
   - If this fails, the issue is with pclm setup

**Action:** Do NOT attempt to fix pclm until Script 03 is corrected and you've checked for negative values.

---

## ⚠️ **Secondary Issues**

### **6. Confusion Between Two Masselot Papers**

You're working with:
- **Masselot et al. (2023)**: 854 cities, baseline period 2000-2019 (the PDF I provided)
- **Masselot et al. (2025)**: Future projections, ISIMIP3, SSPs, adaptation scenarios (referenced in Zenodo)

Your `references/2025-masselot-temp-related/` folder contains the **future projection code**, which uses a different (simpler) formula suitable for projecting to new climate scenarios.

**Our project needs BOTH**:
- **Phase 1**: Replicate 2023 baseline with 4-range disaggregation ← **you are here**
- **Phase 2**: Apply 2023 ERFs to 2025 future projections (SSPs, adaptation)

The 2025 formula is NOT wrong—it's just for a different purpose. But for the baseline replication, you must use the 2023 approach.

### **7. Future Projections Infrastructure Missing**

The paper skeleton states the goal is **2020-2100 projections** with:
- 3 SSPs (1-2.6, 2-4.5, 3-7.0)
- 4 adaptation scenarios (0%, 10%, 50%, 90%)
- Multiple GCMs

**Current status**: None of this exists. Specifically missing:

- ISIMIP3 bias-corrected temperature projections for 854 cities
- SSP population projections (age-specific, by city)
- Code to apply baseline ERFs to future temperature
- Adaptation parameterization framework

**We need to discuss**: Is this in scope for the current paper, or a follow-up study?

### **8. Script 05 Silent Crash**

You noted Script 05 (table generation) crashes without error output.

**Likely cause**: Out of memory (OOM) during Monte Carlo simulation with 500 iterations × 854 cities × 5 age groups.

**Suggested fix**:
- Reduce `n_cores` to 4 (from 16)
- Process one batch at a time
- Add memory monitoring: `pryr::mem_used()` before and after each step
- Consider reducing to 100 simulations for debugging, then scale back up

**Priority**: Low (fix after Scripts 03, 04, 06, 07 are working)

---

## 📋 **Corrected Action Plan**

### **Phase 1: Fix the Baseline (URGENT – This Week)**

| Priority | Task | Est. Time | Blocking? |
|----------|------|-----------|-----------|
| 🔴 **P0** | Update `03_attribution.R` with correct 2023 formula | 2 hours | YES |
| 🔴 **P0** | Re-run Script 03 for all 854 cities | 1 hour | YES |
| 🔴 **P0** | Re-run Scripts 04 (aggregation) and 09 (validation) | 30 min | YES |
| 🔴 **P0** | Verify: Validation error < 0.1% (not 1.7%/6.8%) | 30 min | YES |
| 🔴 **P0** | Implement 4-range disaggregation in Script 03 | 4 hours | YES |
| 🔴 **P0** | Re-run Script 03 with 4-range logic | 1 hour | YES |
| 🔴 **P0** | Validate additivity: sum(4 ranges) = total cold/heat | 1 hour | YES |

**Deliverable**: 854 RDS files with **correct** baseline ANs for 4 temperature ranges.

### **Phase 2: Acquire Missing Data (URGENT – Next Week)**

| Priority | Task | Est. Time | Blocking? |
|----------|------|-----------|-----------|
| 🔴 **P0** | Identify source for total all-cause mortality (854 cities, 2000-2019) | 1 day | YES |
| 🔴 **P0** | Acquire and validate total mortality data | 3 days | YES |
| 🟡 **P1** | Check for negative AN values after formula fix | 1 hour | MAYBE |
| 🟡 **P1** | Debug pclm() NA issue | 4 hours | YES |
| 🟡 **P1** | Test pclm on observed deaths (not ANs) | 2 hours | - |

**Deliverable**: Total mortality data ingested and validated; pclm working for at least one test city.

### **Phase 3: Design Demographic Framework (Next 2 Weeks)**

| Priority | Task | Est. Time | Blocking? |
|----------|------|-----------|-----------|
| 🟡 **P1** | Draft "Demographic Methods" document (2-3 pages) | 1 week | YES |
| 🟡 **P1** | PI review and feedback | 2 days | YES |
| 🟡 **P1** | Revise methods document | 2 days | YES |
| 🟢 **P2** | Implement life table construction | 1 week | - |
| 🟢 **P2** | Implement Aburto decomposition | 1 week | - |
| 🟢 **P2** | Compute LE and LI for baseline period | 3 days | - |

**Deliverable**: Working demographic analysis pipeline for 2000-2019 baseline.

### **Phase 4: Future Projections (Future Work – TBD)**

| Priority | Task | Est. Time | Blocking? |
|----------|------|-----------|-----------|
| 🟢 **P2** | Acquire ISIMIP3 temperature projections | 1 week | YES |
| 🟢 **P2** | Acquire SSP population projections | 1 week | YES |
| 🟢 **P2** | Design adaptation parameterization | 1 week | YES |
| 🟢 **P2** | Apply ERFs to future temperature | 3 days | - |
| 🟢 **P2** | Compute future LE and LI (2020-2100) | 1 week | - |

**Note**: We need to discuss whether future projections are in scope for this paper or a follow-up study.

---

## 🎯 **Immediate Next Steps (Starting Tomorrow)**

### **Step 1: Update Script 03 (Morning)**

1. Open `R_pipeline/03_attribution.R`
2. Locate the attribution calculation loop (around line 150-200?)
3. Replace the current formula with the 2023 approach:
   ```r
   # For each city-age-scenario-sim:
   # 1. Build basis matrix for observed temperatures
   bvar <- onebasis(tmean, fun="bs", degree=2, knots=knots)
   
   # 2. Center at MMT
   bvarcen <- bvar - bvar_mmt
   
   # 3. Compute attributable fraction per day
   af_daily <- 1 - exp(-bvarcen %*% coefs)
   
   # 4. Compute attributable deaths per day
   an_daily <- af_daily * deaths_annual  # NOT deaths_annual/365
   
   # 5. Sum over all days and divide by total days (including leap years)
   an_annual <- sum(an_daily, na.rm=TRUE) / n_days_total
   ```

4. **Do NOT implement 4-range disaggregation yet** (wait for validation)

### **Step 2: Re-run and Validate (Afternoon)**

5. Back up current `temp_results/` → `temp_results_backup_2026-07-24_wrong_formula/`
6. Run `Rscript R_pipeline/03_attribution.R` (854 cities, ~60 min with 16 cores)
7. Run `Rscript R_pipeline/04_aggregate_results.R`
8. Run `Rscript R_pipeline/09_validate_city_level.R`
9. **Expected results**:
   - Cold correlation: > 0.999 (currently 0.9999)
   - Cold median error: < 0.5% (currently 1.70%)
   - Heat correlation: > 0.999 (currently 0.9930)
   - Heat median error: < 1.0% (currently 6.81%)

10. If validation passes, email me the results and proceed to 4-range implementation
11. If validation fails, send me the log files and **DO NOT proceed**

### **Step 3: 4-Range Disaggregation (Day 2)**

Only proceed if Step 2 validation passes.

12. Modify Script 03 to classify temperatures into 4 ranges
13. Compute 4 separate AN values
14. Re-run for all 854 cities
15. Validate additivity
16. Email me the validation report

---

## 📧 **Communication Protocol**

From now on, please:

1. **Daily updates** (end of day):
   - What you completed
   - What's blocking you
   - What you'll work on tomorrow

2. **Before major decisions**, email me for approval:
   - "Should I set negative ANs to zero, or...?"
   - "Should I use period or cohort life tables?"
   - "Should I acquire population projections from source X or Y?"

3. **Validation results**:
   - Send me validation reports after each major pipeline run
   - Include plots, summary statistics, and any anomalies

4. **Blockers**:
   - If you're blocked for >2 hours, email me immediately
   - Don't spend a full day debugging without checking in

---

## 🎓 **Learning Opportunity**

This is a complex project bridging climate science, epidemiology, and demography. The issues I've identified are **normal** in multi-disciplinary research—they're not a reflection of your abilities.

Key lessons:
- **Test scripts ≠ production code**: Your verification script worked, but you didn't update the pipeline
- **Validation is continuous**: Re-validate after every change, not just at the end
- **Check additivity**: For disaggregation, always verify the parts sum to the whole
- **Document assumptions**: Life tables require many choices—document them before coding

I'm confident you'll get this on track. The detective work on the formula discrepancy was excellent, and once we fix the implementation, the rest should flow.

---

## ❓ **Questions for You**

Before you start, please respond to these:

1. **Total mortality data**: Do you know where to get all-cause deaths (854 cities, 2000-2019)? Is it already in the `data/` folder?

2. **Scope**: Should we aim for a "baseline only" paper (2000-2019 with 4-range disaggregation) first, then add future projections in a follow-up? Or is the full 2020-2100 analysis required?

3. **Timeline**: Given the issues identified, how long do you think Phase 1 (corrected baseline) will take? Be realistic.

4. **Blockers**: Is anything else blocking you that I haven't identified?

---

Let's schedule a **30-minute Zoom call tomorrow at 2 PM** to discuss any questions and clarify the path forward.

Best regards,

**[PI Name]**  
Professor of Demography  
[Institution]

---

**Attachments:**
- This email (save for reference)
- Original Masselot 2023 paper (re-read the Methods section)

**Action by EOD today:**
- Reply with answers to the 4 questions above
- Block your calendar for focused work on Script 03 tomorrow