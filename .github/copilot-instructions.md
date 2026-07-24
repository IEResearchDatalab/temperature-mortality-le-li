# Research Agent Protocol

## Session Start Checklist
1. Read `references/PROJECT_STATUS.md` first.
2. Read the active `references/SESSION_GOAL_YYYY-MM-DD.md` file.
3. Confirm that the session has one clear, falsifiable objective.
4. If no session-goal file exists, create one before changing code.
5. Do not start by re-reading archived diaries unless `references/PROJECT_STATUS.md` is missing required detail.

## Testing Protocol
- Test one city first whenever a component touches city-level outputs.
- Test a five-city sample next for workflow changes.
- Only run the full 854-city workload after the sample test passes.
- Use `Rscript temp_scripts/test_pipeline_component.R ...` before full runs whenever possible.
- Treat a passing narrow test as a gate for the next wider run.

## Stop Rules
- Never run the full dataset without a narrow sample check first.
- Never move to Script N+1 while Script N is failing validation.
- Never treat backup folders as live outputs.
- Never work on multiple pipeline stages at once unless the earlier stage is already validated.
- Never use archived session docs as the primary state source.

## Always Rules
- Keep one active objective per session.
- State a falsifiable hypothesis before the first edit.
- Run the narrowest available validation immediately after the first substantive edit.
- Update `references/PROJECT_STATUS.md` before ending the session.
- Create or refresh the next `references/SESSION_GOAL_YYYY-MM-DD.md` before stopping.
- Document blockers as soon as they are confirmed.

## Key Project Facts
- Baseline replication and GCM projections are separate validation tracks.
- `R_pipeline/03_attribution_baseline.R` is the exact ERA5 baseline path.
- `R_pipeline/03_attribution.R` is the GCM projection path.
- The validated baseline target is `references/2025-masselot-zenodo/results/cityage.csv`.
- Current downstream blockers are Script 04/05 baseline integration, total mortality data, and the demographic methods design.

## Session End Checklist
1. Update `references/PROJECT_STATUS.md`.
2. Move completed work out of the in-progress section.
3. Record any new blocker with an exact cause.
4. Create the next session-goal file.
5. Archive throwaway test artifacts if they are no longer useful.
6. Only then consider a commit.

## When to Contact the PI
- Validation misses its expected threshold.
- A design ambiguity changes the methodological interpretation.
- A required external dataset is missing.
- The requested scope needs to change.

## Mental Model: State Machine
- `Baseline validated` -> test baseline aggregation -> test baseline tables -> re-enter LE/LI preparation.
- `Projection path not rerun` -> regenerate `temp_results/` only after the baseline-side workflow is under control.
- `LE/LI blocked` -> do not code past the blocker without mortality data and an approved methods note.

## Good Session Pattern
- Read `references/PROJECT_STATUS.md`.
- Read the session goal.
- Run a one-city or sample test.
- Fix one concrete defect.
- Re-run the same focused test.
- Widen only after the focused test passes.
- Update the project state before stopping.

## Anti-Pattern
- Start from archived notes.
- Edit multiple scripts before testing any of them.
- Launch a full run before a sample run.
- Continue downstream after a failed validation.
- Stop without updating the project state.