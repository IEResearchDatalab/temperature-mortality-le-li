# Cleanup and merge report

Date: 20 October 2026

## Job 1: cleanup

- Ran `git fetch --all --prune`. PR #14 was closed without merging. Its head branch was deleted remotely and locally after the backup verified. The comment says “Closing: this work will live in a separate repository” but omits the requested final period.
- Created and verified `/home/SHARED/backups/phase2-2026-10-07.bundle`. SHA-256: `670fed11a6543b8285a19e92b24ca36b37adccde988feee465c8233553d7a963`.
- Replaced the external prototype-path comment in [analysis/object1_overview.Rmd](../analysis/object1_overview.Rmd) with neutral text. Commit: `fdb96792d94512afe6466804c93f283b112eaa07`.
- The case-insensitive whole-word search of tracked text files returned no matches after the edit. A raw binary search also matched an isolated, unreadable byte sequence in `references/simon_v3/fig2_le65_contributions_mock.png`; it is not readable text.
- The commit-subject and author search over remaining refs returned no matches. The old comment remains in history from `670e6b2d27265cd5a80f4393e1ed1e51db0b62cf`; no history was rewritten.
- At the initial PR audit, PR #13 was the only other open PR. Its title had no restricted-word match, its description was empty, and it had no comments. It was later merged by Daniel.
- No other branches were deleted. `codex/simon-followup` is an ancestor of main. `pipeline/five-script-rebuild` was 19 commits behind `origin/main` at the initial audit, with tip `b61e29c` dated 25 September, rather than the prompt's stated 28 commits and 24 September. After the merge it is 33 commits behind. Both remain local deletion candidates.

## Job 2: validation

Daniel merged PR #13 before the validation was complete, so these checks ran against the merged `main`; they were not pre-merge gates.

- `git status --short --branch` was clean. `git diff --ignore-cr-at-eol --stat` was empty. `core.autocrlf` was unset; no line-ending settings or tracked files were changed for this check.
- R parsing passed with `Rscript -e 'invisible(lapply(list.files(c("pipeline","analysis"), pattern="[.]R$", full.names=TRUE), parse))'`. Every shell script passed `bash -n`.
- Ran Parts 00–04 for Madrid (`ES001C`, SSP3, `GFDL_ESM4`) in isolated output root `results/europe/heatsplit_20261020_madrid_gate`, with `N_HORIUCHI=400` and `DECOMP_ANNUAL=1`. All 60 generated check rows passed.
- The Masselot 2023 fixture maximum relative error was `4.56e-05` against tolerance `0.001`. All 85 moving windows contain 15 years and 5,475 days. The 2015 anchored threshold equals the static threshold exactly. The no-climate-change `NormExtrHeat` attributable number is zero. Two thresholds were clamped.
- Six-cause and annual decomposition checks passed. Maximum errors were `2.770e-09` for LE steps, `1.415e-09` for LI steps, `3.947e-09` for whole-period LE, `4.859e-09` for whole-period LI, `1.606e-09` between branches, and `5.183e-10` within periods. The `rest` contribution was zero. All were below `1e-06`.
- Audited manifests, run logs, output files, and checks under `results/europe/heatsplit_20261013`:
  - SSP3: 16,226 of 16,226 city-GCM jobs, 854 of 854 ensembles, and 35 of 35 pooled geographies; final log reports zero failures.
  - SSP2: 16,226 of 16,226 city-GCM jobs, 854 of 854 ensembles, and 35 of 35 pooled geographies; final log reports zero failures. All 889 ensemble/geography decompositions passed both within-period and between-branch checks; maximum error was `9.592e-07`, below `1e-06`.
  - SSP1: 16,226 of 16,226 city-GCM jobs, 854 of 854 ensembles, and 35 of 35 pooled geographies; the final precision-repair run reports zero failures. The manifest records the labelled `1` to `0.9995` survival-ratio sensitivity, only for SSP1. At least one source row changed in 771 cities, with 26,183 rows changed in total. The 5% mortality-shape check reports warnings in 562 cities; the run bypasses that guard by design and retains the warning. The other 292 cities pass the shape check.
- For SSP2, 12 Swiss cities have a within-age-band mortality decline above 5%; the maximum is 11.09%. The run keeps these cities with a warning. The other 842 cities pass this shape check.
- The SSP1 run manifest records `git_worktree=dirty`; SSP2 and SSP3 manifests record clean worktrees. The pipeline sources are unchanged between the result revisions and merged `main`. The separate attached 24 September run log with 23 failures refers to the older `results/europe/ssp3` root, not the audited heat-split root.
- `git diff origin/main^1 origin/main --stat -- results data` was empty. The merge changed no canonical result or source CSV.

## Job 3: merge

- PR #13 was merged by Daniel on 20 October 2026 at 10:50:34 UTC, before Jobs 1 and 2 were complete. Merge commit: `dc98029df216dce392d3191dd4943c0836d5d082`. Local `main` is fast-forwarded to `origin/main`. The merge contains the heat-split changes and the Rmd cleanup commit.
- The tracked-text search on `main` returned no matches. The commit-subject and author search also returned no matches. The historical source comment remains at the commit listed in Job 1.
- Created and pushed tag `heat-split-merged-20261020` to `dc98029df216dce392d3191dd4943c0836d5d082`.
- `git log --oneline -5` on main:

  ```text
  dc98029 Merge pull request #13 from IEResearchDatalab/heat-split-15y-window
  fdb9679 Remove reference to external prototype path
  1d440a0 Merge branch 'main' into heat-split-15y-window
  ac0137a Wrap SSP1 figure caption notes
  99a5809 Move SSP1 calibration note to figure captions
  ```

- No pull requests are open. The only remote branch is `main`. Local branches retained are `agent/pipeline-phase0`, `agent/pipeline-rebuild`, `alt_pipeline`, `backup/five-script-rebuild-orig`, `codex/simon-followup`, `full-refactor`, `heat-split-15y-window`, `main`, and `pipeline/five-script-rebuild`. No branch other than the closed PR #14 head was deleted.
- PR #13's description remains empty. It had already been merged, and the available browser session was signed out, so I could not update its description after the merge.

### Open questions for Daniel

- Decide whether to edit PR #13's description after the merge.
- Decide whether to leave the historical source comment in place; history was not rewritten.
- Decide whether to delete the two local branch candidates listed in Job 1.
- Decide whether the missing final period on PR #14's closing comment needs correction.
