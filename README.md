# Meinungsbild

Subnational public opinion estimates for Germany using multilevel regression and poststratification (MRP).

## Overview

Meinungsbild estimates public opinion on 43 policy issues at three geographic levels:

- **Bundeslaender** (16 federal states)
- **Wahlkreise** (299 electoral districts)
- **Kreise** (400 counties)

Estimates are derived from ~118,000 survey respondents across five survey programs, combined with Zensus 2022 population data for poststratification.

## Data sources

| Source | N | Years | Geographic IDs |
|--------|---|-------|----------------|
| GLES Tracking (ZA6832) | 52,336 | 2009--2023 | Bundesland, Wahlkreis |
| GLES Cross-Section 2025 (ZA10100) | 7,337 | 2025 | Bundesland, Wahlkreis |
| GLES RCS 2025 (ZA10101) | 8,561 | 2025 | Bundesland, Wahlkreis |
| GLES Cumulation (ZA6835) | 21,040 | 2009--2021 | Bundesland |
| ALLBUS (ZA8974) | 29,112 | 2023--2024 | Bundesland |

Raw survey data are included for collaborators. These are restricted-access files from [GESIS](https://www.gesis.org/) and must not be redistributed publicly.

## Repository structure

```
meinungsbild/
├── code/                        # R pipeline scripts (run in order)
│   ├── 01_harmonize_all.R       # Pool and harmonize 5 survey datasets
│   ├── 02_build_poststrat_frame.R
│   ├── 02b_download_zensus.R    # Download Zensus 2022 cross-tabs
│   ├── 03_load_covariates.R     # Election results + INKAR covariates
│   ├── 03b_build_adjacency.R    # County adjacency matrix
│   ├── 04b_fit_all_lme4.R      # Fit MRP models (lme4) — primary
│   ├── 04_fit_mrp.R            # Experimental: brms/Stan models
│   ├── 05_poststratify.R       # Poststratification step
│   ├── 07_export_estimates.R    # Export to JSON/CSV
│   ├── 08_check_pipeline.R     # Validation checks (135 tests)
│   └── 09_validate_vote_shares.R
├── data/
│   ├── raw/                     # Survey microdata (GLES, ALLBUS, Zensus)
│   ├── harmonized/              # Pooled individual-level data
│   ├── covariates/              # Geographic covariates (.rds)
│   ├── estimates/               # Model estimates (.rds)
│   ├── poststrat/               # Poststratification frames (.rds)
│   ├── issue_concordance.csv    # Issue definitions + binary coding rules
│   └── variable_inventory.csv   # Variable mapping across surveys
├── docs/
│   └── methodology_notes.md     # Full model specification
├── harmonization/               # Variable-specific harmonization notes
├── output/
│   ├── checks/                  # Validation results (correlations, RMSE)
│   └── tables/                  # CSV estimates for download
└── web/                         # Next.js interactive map
    ├── public/data/             # GeoJSON boundaries + JSON estimates
    └── src/                     # React/TypeScript components
```

## Setup

### R pipeline

Scripts reference external GERDA data (election results, shapefiles) via `gerda_root`. Update paths in `03_load_covariates.R` and `03b_build_adjacency.R` if running outside the `german_election_data` directory.

```bash
Rscript code/01_harmonize_all.R
Rscript code/02_build_poststrat_frame.R
Rscript code/03_load_covariates.R
Rscript code/04b_fit_all_lme4.R
Rscript code/07_export_estimates.R
Rscript code/08_check_pipeline.R
```

### Website

```bash
cd web
npm install
npm run dev
```

Open [http://localhost:3000](http://localhost:3000).

## Methodology

MRP with `lme4::glmer()`. See [`docs/methodology_notes.md`](docs/methodology_notes.md) for the full model specification, design choices, and validation results.

**Validation:** Median correlation r = 0.899 and median RMSE = 5.5pp against direct Bundesland survey estimates across 43 issues.

## Key files

| File | Description |
|------|-------------|
| `data/issue_concordance.csv` | Defines all 43 issues: variable names per survey, binary coding rules, response scales |
| `code/04b_fit_all_lme4.R` | Primary estimation script — fits one glmer per issue, poststratifies, exports |
| `code/01_harmonize_all.R` | Harmonization — pools 5 surveys, applies binary coding from concordance |
| `output/checks/validation_bundesland.csv` | Bundesland-level validation (r, RMSE per issue) |

## Output formats

- **JSON** (`web/public/data/`): Consumed by the interactive map
- **CSV** (`output/tables/`): Downloadable estimates at all three geographic levels

## Note on paths

The R scripts were originally written to run from within the `german_election_data` repository, where `here::here()` resolved to the GERDA root. Scripts that only use meinungsbild data (`04b_fit_all_lme4.R`, `07_export_estimates.R`, etc.) use `mb_root <- file.path(here::here(), "meinungsbild")` which will need updating to `mb_root <- here::here()` for standalone use. Scripts `03_load_covariates.R` and `03b_build_adjacency.R` also read from GERDA directories and need an explicit path to that data.
