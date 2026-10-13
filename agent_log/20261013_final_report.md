# Final report — stopped at the Madrid heat-window gate

Date: 2026-10-13

Branch: `heat-split-15y-window`

Base: `21d30438ac939dff17a7f233404ee8f4cc414c90`

## Outcome

The run stopped in Phase 2 before any full SSP batch, exactly as required by the prompt's stop conditions. Canonical SSP roots, collected objects and pooled geographies were not modified. The only intentional change outside the branch/new scratch root was restoration of the documented canonical SSP1 `failed.txt`; its MD5 matches the preserved source (`5a598393ebbcf41ac38fca80dbecc57d`).

The blocker is an irreconcilable numerical conflict between two explicit requirements. The moving 2015 P97.5 must use Masselot's calibrated historical GCM `full` series for 2000–2014, while the static P97.5 must use raw ERA5-Land 2000–2014, and the two must agree to numerical precision. On Madrid/GFDL-ESM4/SSP3, the direct implementation gives:

- calibrated-GCM moving P97.5: 28.2899418625876°C;
- static ERA5-Land P97.5: 28.2010000000000°C;
- absolute difference: 0.0889418625876°C (required tolerance: 1e-8°C).

Masselot's ISIMIP3BASD code trend-preserves and distribution-maps the historical GCM series but does not force its empirical P97.5 to equal the ERA5 empirical P97.5. Replacing the window with ERA5 would violate the explicit “do not substitute ERA5” instruction. Shifting/anchoring the moving series or special-casing 2015 would add a method not specified by Masselot or Simon. No such adjustment was improvised.

## Work completed before the stop

- Audited the skeleton, mocks, specified vault files, all five pipeline scripts, Masselot's archived code and Lloyd's attribution script.
- Added `references/simon_figures_spec.md`; copied the three mocks and skeleton text to `references/simon_v3/`.
- Audited canonical SSP3: 854 cities, 19 GCMs plus ensemble, 2020–2099, and all 35 pooled geographies.
- Implemented issue 5's batch-root-specific prep cache fingerprint guard and issue 7's city-figure root override.
- Added the 2000–2014 static threshold parameter, 15-year window parameter, `NormExtrHeat`, six-cause downstream contract, no-CC empty-band rule, clamp diagnostics, and default 5-year-period decomposition setting.
- Updated `run_config.csv` for SSP3 → SSP2 → SSP1, with N=400, SSP2 Swiss inclusion bypass, and SSP1 0.9995 WCDE rounding patch.
- Added `PIPELINE_RULES.md` with the six-cause data contract. Shell and R syntax checks pass.

## Tests and resources

The Madrid run passed temperature preparation, demography, complete domains, non-negative AN, the Masselot 2023 fixture (maximum relative error `4.56e-05`, limit `0.001`), counterfactual drift, complete 15-year windows, the no-CC zero NormExtrHeat check, and clamp counting (1/85 years). It failed only `moving_2015_equals_static`. The gate command took about 21 seconds wall time and wrote 41 MB under `results/europe/heatsplit_20261013_test/`. Free disk remained about 2.4 TB. Parts 02–04, annual decomposition, full reruns, figures, tables, PDF and methods text were not run because the gate failure requires stopping.

## Defaults and departures

Central ERF estimates only; 28 batch cores planned; N_HORIUCHI=400; DECOMP_ANNUAL=0; SSP1 ASSR 1.000→0.9995 with warnings; SSP2 Swiss mortality-shape bypass with warnings. The only scientific departure implemented is the requested three-range heat split, clearly marked in Part 01. No departure was added to reconcile the failed equality.

## Next decision and commands

Daniel/Simon must choose one rule: accept calibrated-GCM and static-ERA5 P97.5 as close but not identical; use ERA5 for the first window; or define an explicit anchoring adjustment. After that decision, amend the test and resume.

```bash
cd /home/SHARED/temperature-mortality-le-li
git push -u origin heat-split-15y-window
# Then open a draft PR using pr_drafts/heat_split_pr.md.
```

No promotion command is provided because there is no completed `results/europe/heatsplit_20261013/` root to promote.
