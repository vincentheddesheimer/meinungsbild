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

### Berlin Bezirk-level poststratification

Berlin is a single county (AGS 11000) containing 12 Wahlkreise. Without sub-county demographics, all 12 WKRs would get identical poststratification weights. To fix this, `04b_fit_all_lme4.R` loads Bezirk-level demographics from `data/poststrat/poststrat_wkr_berlin.rds` and substitutes them for Berlin WKRs 75–86.

**Data sources:**
- **Age × sex**: EWR (Einwohnerregisterstatistik) 2020 at LOR Planungsraum level, aggregated to Bezirk
- **Education (Schulabschluss)**: Zensus 2022 table `2000S-3041` at Bezirk level (age × sex × Schulabschluss)
- **University split**: Zensus 2022 table `2000S-4028` (Berufsabschluss) to distinguish Abitur-only from university degree holders

**Bezirk→WKR mapping** (`data/poststrat/berlin_wkr_bezirk_mapping.csv`): 1:1 approximation using each WKR's primary Bezirk. WKR 78 (Spandau–Charlottenburg Nord) and WKR 83 (Friedrichshain-Kreuzberg–Prenzlauer Berg Ost) span Bezirk boundaries; they use Spandau (BEZ 05) and Friedrichshain-Kreuzberg (BEZ 02) respectively.

**Age group alignment**: Zensus uses 10-year groups (20–29, 30–39, etc.) while our model uses 18–29, 30–44, 45–59, 60–74, 75+. The 40–49 group is assigned to 30–44 and 70–79 to 60–74 as approximations.

### Web app (`web/`)

Next.js 15 + MapLibre GL JS choropleth map. TypeScript/React.

- `web/src/app/page.tsx` — main page, loads all data, manages state
- `web/src/components/Map.tsx` — MapLibre GL map with three geo levels
- `web/src/lib/i18n.tsx` — DE/EN translations (categories, labels, directions)
- `web/src/lib/types.ts` — TypeScript interfaces for estimates and geo levels
- `web/src/lib/colors.ts` — choropleth color scale and MapLibre expressions
- `web/public/data/` — GeoJSON boundaries + JSON estimates (served statically)

**GeoJSON property names:** `state_code` (string) for Bundeslaender, `wkr_nr` (integer) for Wahlkreise, `county_code` (string) for Kreise. The `choroplethExpression()` function uses `numericKeys=true` for Wahlkreise since MapLibre's `match` expression requires type-consistent keys.

**Deployment:** Data files are copied to `awiedem.github.io/assets/data/meinungsbild/` for the live website at german-elections.com. The JS loads from `/assets/data/meinungsbild/` (not raw.githubusercontent.com, which doesn't work with Git LFS).

## Common tasks

### Re-running the pipeline

After changing `issue_concordance.csv` or harmonization code:
```bash
Rscript code/01_harmonize_all.R
Rscript code/04b_fit_all_lme4.R    # ~45-90 min
Rscript code/07_export_estimates.R
Rscript code/08_check_pipeline.R
```

### Updating the live website

After re-exporting estimates, copy data to the website repo and push:
```bash
cp web/public/data/*.json web/public/data/*.geojson \
   /path/to/awiedem.github.io/assets/data/meinungsbild/
cd /path/to/awiedem.github.io && git add assets/data/meinungsbild/ && git commit && git push
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
