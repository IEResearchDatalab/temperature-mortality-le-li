# Research pipeline rules

These rules govern the code for the temperature-mortality, life-expectancy and lifespan-inequality paper with Simon Lloyd. Simon is the leading researcher on the project; his methodological objective and reviewing needs are the guiding compass for all code decisions.

## Core structure

1. The scientific pipeline has exactly five ordered scripts.
2. Each script performs one clear, named scientific step and produces a stable, documented output that is consumed by the following step.
3. The five scripts are the authoritative implementation. Scientific functions and transformations should be written in those scripts rather than split across helper files or parallel implementations.
4. A configuration file may hold parameters, paths, geography selections and run settings. Configuration selects an analysis; it must not contain a second implementation of an analysis.
5. Shell runners may orchestrate repeated execution, but they must not contain scientific transformations.

## Scientific scope

6. The same pipeline must work at every required geographical level: city, country, region and Europe. Geography may be selected through configuration, or all levels may be produced together, but separate analytical implementations are not allowed.
7. Aggregation follows Simon's rule: estimate attributable deaths at city level, sum population and deaths by cause across cities, and only then construct pooled mortality rates, life tables and decompositions. Never average city ERFs, life expectancies, lifespan inequalities or decomposition contributions to obtain pooled results.
8. Every stage must preserve a clear data contract: required inputs, primary keys, units, geographical level and output fields must be evident from the script and documentation.
9. Validation checks belong next to the scientific step they validate. They may write check tables, but they must not create a parallel analytical pipeline.

## Style and reviewability

10. Keep the structure, naming and coding style close to Masselot's reference code wherever the scientific method is the same. This is important for transparent comparison and an efficient review.
11. Avoid redundant code. There must be one implementation of each scientific calculation, especially the temperature attribution, life table and Horiuchi decomposition.
12. Prefer friendly, explicit and readable code over compact or highly optimized code. A reviewer should be able to follow the analysis from top to bottom without reconstructing hidden abstractions.
13. Use short sections, descriptive object names and visible intermediate calculations. Do not introduce abstractions unless they remove genuine duplication without obscuring the method.
14. Changes that depart from Masselot, Lloyd or Simon's agreed specification must be clearly identified and justified in the code and methods documentation.

## Results and presentation

15. Figures, manuscript tables, rankings, smoothing choices, colour palettes and headline summaries are not pipeline stages. They belong in clearly separated auxiliary analysis/result scripts that read completed pipeline outputs.
16. The main pipeline should produce scientific data products, not paper presentation. Diagnostic plots are optional checks and must not be required outputs of a pipeline stage.
17. Presentation scripts must never recompute or alter the underlying scientific estimands.

## Heat-split data contract (13 October 2026)

- Part 00 writes a batch-root-specific, fingerprinted temperature cache. Static P2.5/P97.5 thresholds use ERA5-Land 2000–2014; ERF knots, boundaries and MMT estimation continue to use the full ERA5-Land 1990–2019 series.
- Part 01 keys are branch × year × age group × temperature range. The five temperature ranges are `ExtrCold`, `ModCold`, `ModHeat`, `NormExtrHeat`, and `ExtrHeat`. `NormExtrHeat` is the band from the static P97.5 to the clamped preceding-15-year P97.5; it is empty in the without-CC branch. Units are annual attributable deaths after summing daily AN and dividing by 365.
- Part 02 preserves those five temperature-range keys while allocating grouped attributable deaths to single ages 65–100+.
- Part 03 adds `rest`, so Object 1 and all master-table cause keys have exactly six levels: `ExtrCold`, `ModCold`, `ModHeat`, `NormExtrHeat`, `ExtrHeat`, `rest`. Population is persons and deaths are annual counts.
- Part 04 consumes the six-cause grid. Object 2 contributions are years of LE65 or LI65+ unless presentation code converts them to months/days; the standard within-branch decomposition compares consecutive 5-year-period mean mortality schedules. Object 3 stores LE65 and LI65+ levels in years. The cause-major, age-minor vector order is fixed by `cause_levels`.

## Decision test for future changes

Before adding a file or stage, ask:

- Does this produce a new, necessary scientific data product?
- Can the same result be obtained through the existing five stages using configuration?
- Does it duplicate a calculation that already exists?
- Is it actually a figure, table, diagnostic or report rather than pipeline work?
- Will its structure make Simon's review easier and keep comparison with Masselot transparent?

If it is only a parameter change, put it in configuration. If it is presentation or diagnosis, put it outside the pipeline. If it is a necessary scientific transformation, incorporate it into the relevant one of the five stages and revise the documented output contract rather than adding a sixth stage.
