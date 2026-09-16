# Phase 0 — ASSR-to-annual-deaths methodology

## Verdict
PASS

- Local evidence: `results/projdata/build_projdata_prototype.R` and `results/projdata/projdata_prototype.csv`.
- Interval mapping: five-year period survival ratio in `assr`, annualized by dividing five-year deaths by 5 after aggregation.
- Sex aggregation: sum over sex to country totals before age-group collapse.
- Age transitions: repository groups `20-44`, `45-64`, `65-74`, `75-84`, `85+`.
- 85+ handling: age labels with lower bound 85 are retained as the terminal open-ended group.

## Spain reproduction check
- Spain rows compared: 24300
- Max |pop difference|: 4.889e-09
- Max |death difference|: 5.093e-11
- Max |wittdeath difference|: 1.164e-10

## Interpretation
- The Spain subset reproduces exactly against the saved prototype to floating-point tolerance.
- The audited historical bridge therefore supports the current annual-death method and the chosen deterministic aggregation order.
- The method is PASS for the Phase 0 gate.
