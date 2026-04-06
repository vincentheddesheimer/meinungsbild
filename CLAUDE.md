# Meinungsbild — Claude Code Instructions

## Project overview

Subnational German public opinion estimates using MRP (multilevel regression and poststratification). Estimates 43 binary policy issues at three geographic levels (Bundeslaender, Wahlkreise, Kreise) from ~118,000 survey respondents.

## Architecture

### R pipeline (`code/`)

Scripts run sequentially (`01_` through `08_`). The primary estimation script is `04b_fit_all_lme4.R` (lme4/frequentist). `04_fit_mrp.R` is an experimental brms/Stan alternative.

**Path convention:** Scripts were written for the `german_election_data` repo where `here::here()` resolved to the GERDA root and `mb_root <- file.path(here::here(), "meinungsbild")`. These paths need updating for standalone use — `mb_root` should become `here::here()`.

**External dependency:** `03_load_covariates.R` and `03b_build_adjacency.R` read election results and shapefiles from GERDA (`german_election_data/data/`). These require access to that directory.

### Issue concordance (`data/issue_concordance.csv`)

**Source of truth** for all variable coding. Each row maps an `issue_id` to a specific `variable` in a specific `dataset`, with a `binary_rule` (e.g., `y=1 if <= 2`) and `response_type` (e.g., `scale_1_5`).

Critical coding details:
- **ALLBUS** uses 1=yes/2=no for `binary_yn` items (not 0/1). Rules must be `<= 1` for "yes", `>= 2` for "no".
- **GLES tracking** (e0113* battery) uses **reversed** scale: 1=strongly disagree, 5=strongly agree. This is the **opposite** of GLES cross-sections (1=strongly agree, 5=strongly disagree).
- **Fear scales** (scale_1_7): Both tracking and cross-section use 1=no fear, 7=great fear.
- The `apply_binary_rule()` function in `01_harmonize_all.R` matches by both `issue_id` AND `variable`, so different datasets correctly get different rules.

### Web app (`web/`)

Next.js 15 + MapLibre GL JS choropleth map. TypeScript/React.

- `web/src/app/page.tsx` — main page, loads all data, manages state
- `web/src/components/Map.tsx` — MapLibre GL map with three geo levels
- `web/src/lib/i18n.tsx` — DE/EN translations (categories, labels, directions)
- `web/src/lib/types.ts` — TypeScript interfaces for estimates and geo levels
- `web/src/lib/colors.ts` — choropleth color scale and MapLibre expressions
- `web/public/data/` — GeoJSON boundaries + JSON estimates (served statically)

**GeoJSON property names:** `state_code` (string) for Bundeslaender, `wkr_nr` (integer) for Wahlkreise, `county_code` (string) for Kreise. The `choroplethExpression()` function uses `numericKeys=true` for Wahlkreise since MapLibre's `match` expression requires type-consistent keys.

## Common tasks

### Re-running the pipeline

After changing `issue_concordance.csv` or harmonization code:
```bash
Rscript code/01_harmonize_all.R
Rscript code/04b_fit_all_lme4.R
Rscript code/07_export_estimates.R
Rscript code/08_check_pipeline.R
```

### Adding a new issue

1. Add rows to `data/issue_concordance.csv` (one per dataset that measures it)
2. Run the full pipeline above
3. The website picks up new issues automatically from `issues.json`

### Adding translations

All UI strings are in `web/src/lib/i18n.tsx`. Categories need a `cat:` prefixed key in both `de` and `en` dictionaries.

## Data sensitivity

- `data/raw/` contains restricted GESIS survey microdata — do not redistribute publicly
- `data/harmonized/` contains individual-level data derived from restricted sources
- `output/tables/` and `web/public/data/` contain only aggregated estimates (safe to share)

## Validation

`output/checks/validation_bundesland.csv` has per-issue Bundesland correlations. Current median r = 0.899, median RMSE = 5.5pp. The pipeline check script (`08_check_pipeline.R`) runs 132 automated checks.
