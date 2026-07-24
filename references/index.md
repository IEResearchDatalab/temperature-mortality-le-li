# References Index

This folder contains the scientific papers, code repositories, and datasets that serve as the foundation for the current project.

## Project Documents
- [../PROJECT_STATUS.md](../PROJECT_STATUS.md): Canonical live status for the project. Read this first.
- [../SESSION_GOAL_TEMPLATE.md](../SESSION_GOAL_TEMPLATE.md): Template for defining one clear objective per session.
- [paper_skeleton.md](paper_skeleton.md): The structural outline and methodological sections for the current study.
- [task_understanding.md](task_understanding.md): Internal documentation on project requirements and implementation goals.

## Archived Session Documents
- [archive_sessions/](archive_sessions/): Historical diaries, TODOs, and PI correction notes. Use only for background context after reading `PROJECT_STATUS.md`.

## Core Reference Literatures & Codebases

### 1. Masselot et al. (2025) - Main Reproducibility Target
- **Folder**: [2025-masselot-temp-related/](2025-masselot-temp-related/)
- **Description**: The primary reference for the 4-part simulation pipeline (Prep, Attribution, Aggregation, Plotting).
- **Key Functions**: ISIMIP3 bias correction and health impact projections.
- **Data Source**: [2025-masselot-zenodo/](2025-masselot-zenodo/) contains the multi-city coefficients (`coefs.csv`), variance-covariance matrices (`vcov.csv`), and simulation coefficients (`coef_simu.csv`).

### 2. Lloyd et al. (2024) - Life Expectancy (LE) & Life Years Lost (LI)
- **Folder**: [2024-lloyd-reciprocal/](2024-lloyd-reciprocal/)
- **Description**: Sample code and fake data for life expectancy decomposition using the Aburto methodology.
- **Key Scripts**: `Code_1.R` contains the life table and decomposition logic.
- **Related**: [2024_lloyd_spain/](2024_lloyd_spain/) contains additional scripts for attributable number calculations.

### 3. Masselot et al. (2023) - Excess Mortality & Tables
- **Folder**: [2023-masselot-excess/](2023-masselot-excess/)
- **Description**: Methodology for age-stratified health impacts in European cities.
- **Reference Tables**: Source of the "Table S6" format used in our aggregation script.

### 4. Gasparrini et al. (2014) - DLNM Foundations
- **Folder**: [2014_gasparrini_BMCmrm_Rcodedata-master/](2014_gasparrini_BMCmrm_Rcodedata-master/)
- **Description**: The original implementation of temperature-attributable mortality using Distributed Lag Non-linear Models (DLNM).
- **Key Script**: `attrdl.R` contains the logic for calculating attributable numbers from DLNM objects.

---
*Last Updated: 2026-07-25*
