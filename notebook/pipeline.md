## Proposed notebook series

I recommend **four demo notebooks**, all in **R Markdown**.

### Notebook 1

**Projected temperature-attributable mortality for one city, one GCM, one SSP, using wide age groups**

Purpose:
- implement the standard **AN** workflow
- stay as close as possible to **Masselot 2023** and **Gasparrini & Leone 2014**
- stop **before** single-age disaggregation and life tables

### Notebook 2

**Disaggregate annual ANs from wide age groups to single years of age using `pclm`**

Purpose:
- convert age-group ANs into single-age ANs
- build the age-specific mortality-impact quantities needed downstream

### Notebook 3

**Construct the “dataset for analysis” and derive period life-table inputs**

Purpose:
- combine:
    - all-cause mortality
    - population
    - single-age ANs by temperature range
- construct the table structure Simon referred to from Lloyd’s supplementary material

### Notebook 4

**Decompose changes in LE and LI by age and temperature range for one city**

Purpose:
- apply the Lloyd-style decomposition
- compute:
    - LE
    - LI
    - contributions by age and temperature range

Only after these four are stable would I expand to:
- multiple GCMs
- multiple cities
- country / Europe summaries
- adaptation scenarios
