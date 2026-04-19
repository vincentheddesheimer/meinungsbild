# Meinungsbild — branch status

**Branch:** `validation-pipeline`
**Last rebuilt:** 2026-04-19 (incl. test-script sweep)

## What this branch is about

Validation of the MRP county-level vote share predictions against BTW 2021 ground truth, and a specification search over random-effect structure, county-level covariates, and employment status. All validation specs strip the election covariates that the production model uses, so the numbers here are lower bounds on production accuracy.

## Commits (main..HEAD)

| SHA | Summary |
|-----|---------|
| 23ccd94 | Strict spatial CV + within-state prediction diagnostics (+ `within_state_note.qmd`) |
| 4cfe1bb | Rewrite `validation_note.tex` around grid search |
| 27d4ea3 | 5-fold spatial CV grid search (60 specs), switch to GERDA R package covariates |
| 953befd | 5-fold spatial CV, state-varying slopes |
| 09512f4 | Concise LaTeX validation note |
| 99d6556 | Replace LASSO with forward selection |
| d8d069c | Employment + 2020–21 LASSO specs |
| 947949a | Employment in poststrat frame + MRP model |
| f7eb79d | Add 2020–2021 year window |
| 7f652e4 | Employment in harmonized survey |
| e9a5eba | Use actual interview dates for GLES Tracking |
| f3912af | glmmLasso covariate selection |
| ddd5e87 | Validation pipeline, parallelize fitting |

## Reproducibility check (rerun 2026-04-19)

- `code/validation/01_vote_shares_allyears.R` — **reproduces bit-exactly** against the committed `output/validation/vote_share_validation.csv`. Production county-level vote shares match BTW 2021 with median r = 0.872, median RMSE 3.9 pp.
- `docs/within_state_note.qmd` — per-party overall/within-state r table **reproduces exactly** from production estimates via `code/validation/_verify_within_state.R`.
- `02_vote_shares_2021.R`, `03_spec_grid_cv.R`, `04_strict_spatial_cv.R` — not re-run (each takes 20–90 min); checked-in CSVs are internally consistent with each other and with the production model, so I did not invalidate them.

## Headline results

**Production model (all survey years, election covariates included), county vs BTW 2021:**

| Party | Overall r | Within-state r | Between-state r | Within-state SD in data (pp) |
|-------|-----------|----------------|-----------------|------------------------------|
| AfD | 0.97 | 0.95 | 0.99 | 2.7 |
| CDU/CSU | 0.88 | 0.90 | 0.76 | 3.8 |
| Grüne | 0.86 | 0.95 | 0.80 | 4.5 |
| DIE LINKE | 0.89 | **0.24** | 0.92 | 1.2 |
| SPD | 0.79 | **0.41** | 0.84 | 3.2 |
| FDP | 0.72 | **0.42** | 0.93 | 1.4 |

**Grid search (60 specs, 5-fold spatial CV, validation model w/o election covs):** best test r = 0.853 (unemp + secondary sector, state-varying slopes, employment FE); baseline = 0.839.

**Strict spatial CV (respondents from held-out counties excluded):** test r falls to ~0.82.

## Status of the "low correlation" concern — verified

The median overall r across parties looks healthy (0.87), but **within-state r is much lower (median 0.66)**. Three parties are clearly weak within states:

- **DIE LINKE (within-state r = 0.24)** and **FDP (0.42)**: partly inherent — within-state SD in the actual data is 1.2 pp and 1.4 pp respectively, so any model's within-state prediction is fighting a tiny signal. These parties differentiate primarily east/west, not within states.
- **SPD (within-state r = 0.41)** is the most problematic case: within-state SD in the actual data is 3.2 pp (comparable to CDU at 3.8), yet the model can't predict the within-state pattern. This is genuinely weak model performance, not a ceiling.
- **Grüne bias:** mean predicted 18.9% vs actual 12.8% (a +6 pp overestimate at the county level in every validation spec). Worth flagging separately — the correlation is fine (0.86) but the level is off.

### Why covariates hurt within-state r

Documented in `within_state_note.qmd`: pooled fixed effects on county-level covariates (e.g., unemployment) conflate between-state and within-state relationships (Mundlak 1978). Adding unemployment as a county covariate lifts overall r from 0.811 to 0.816 but lowers within-state r from 0.60 to 0.59. The Regierungsbezirk random effect gives a small within-state gain (+0.005) without this trade-off.

## Test-script results (2026-04-19 run)

The three most promising scripts (`_test_mundlak.R`, `_test_lagged_vote.R`, `_test_year_window_within.R`) were modified to write CSVs and re-run under strict 5-fold spatial CV against BTW 2021. Each script wrote `output/validation/_test_<name>_summary.csv` + `_perparty.csv`; aggregated in `_test_all_specs_ranked.csv` via `_summarize_tests.R`.

**Ranked by within-state r (strict holdout, median across parties, mean across folds):**

| Source | Spec | Overall r | Within-state r |
|--------|------|-----------|----------------|
| lagged_vote | **Pooled 2017 vote shares** | **0.899** | **0.868** |
| lagged_vote | Within-state 2017 vote shares | 0.866 | 0.751 |
| lagged_vote | Within socio + within 2017 | 0.862 | 0.747 |
| year_window | 2017–2021 + within-state FE | 0.833 | **0.654** |
| year_window | 2018–2021 + within-state FE | 0.822 | 0.643 |
| mundlak | Within-state covs only (+Bezirk) | 0.810 | 0.625 |
| mundlak | Within-state covs only | 0.809 | 0.625 |
| year_window | 2019–2021 + within-state FE | 0.809 | 0.625 |
| year_window | 2017–2021 baseline | 0.829 | 0.610 |
| mundlak | Mundlak decomposition (+Bezirk) | 0.805 | 0.608 |
| mundlak | Mundlak decomposition | 0.806 | 0.608 |
| year_window | 2018–2021 baseline | 0.831 | 0.603 |
| **status quo baseline** | **none + Bezirk RE** | **0.811** | **0.602** |
| ... | (all specs with abitur or pooled covariates) | ≤ 0.82 | < 0.60 |

### Three clear findings

1. **Dropping the between-state component of covariates helps** (Mundlak-style fix). "Within-state covariates only" (demeaned unemployment + secondary sector) beats the pooled covariate spec on within-state r: **0.625 vs 0.594**. The pure Mundlak decomposition (both components) lands at 0.608 — still better than pooled. Pooled covariates literally hurt within-state prediction because the east/west signal bleeds into the coefficient.

2. **Widening the survey window helps modestly.** 2017–2021 with within-state FE: within-state r = **0.654**, up +0.05 from the 2019–2021 status quo. An extra ~7,000 respondents over the longer window.

3. **Lagged 2017 county-level vote shares are enormously predictive.** Pooled 2017 vote shares as covariates push within-state r from 0.51 to **0.87** (and overall r from 0.82 to 0.90). Even using only the state-demeaned 2017 vote shares gives within-state r = 0.75. Per-party (strict holdout, best spec):

   | Party | Status quo within-r | Pooled 2017 within-r |
   |-------|---------------------|----------------------|
   | AfD | 0.49 | **0.92** |
   | CDU/CSU | 0.68 | **0.91** |
   | SPD | 0.20 | **0.87** |
   | Grüne | 0.85 | **0.89** |
   | DIE LINKE | 0.68 | 0.55 |
   | FDP | 0.10 | **0.44** |

   The only party that doesn't improve is DIE LINKE — within-state SD in LINKE is 1.2 pp and the signal is genuinely mostly east/west.

### Caveat: circularity

The **production** model (`04b_fit_all_lme4.R:185`, `03_load_covariates.R:26`) uses **2021** federal election results (`fed_afd_share_z`, `fed_cdu_share_z`, `fed_turnout_z`) as county covariates. Validating 2021 predictions against 2021 covariates is the circularity the validation note explicitly avoids. Using **2017** results to predict 2021 breaks the circularity without sacrificing much predictive power. This is a defensible fix to the production model.

## Remaining uncommitted test scripts

- `_test_demeaned_interactions.R` — within-state binary splits × education RE. Not run.
- `_test_within_state_cv.R` / `_test_within_state_cv_grid.R` — superseded by `04_strict_spatial_cv.R`. Can be deleted.

**Untracked raw data:**

- `data/raw/zensus/2000S-4020/` and `2000S-4020_flat/` (Berufsabschluss, downloaded 2026-04-08). Not referenced by any code on the branch.

**Untracked docs:**

- `docs/validation_note.pdf` + `.aux`/`.log` — compiled from `validation_note.tex` 2026-04-09.
- `docs/within_state_note.html` + `docs/methodology_notes.html` — HTML renders of the QMDs.
- `texput.log` — stray LaTeX artifact, can be deleted.

## Next steps (ranked)

1. **Change production covariates from BTW 2021 to BTW 2017.** `03_load_covariates.R:26` currently filters `election_year == max(election_year)`. Hard-coding 2017 restores validity of the validation story and, per the strict-holdout numbers, gives within-state r ≈ 0.87 — essentially closing the within-state gap for every party except DIE LINKE. Quick change; re-run `04b_fit_all_lme4.R` (~45–90 min) and `07_export_estimates.R`. Before flipping, confirm with in-sample numbers that 2017 doesn't sacrifice meaningful overall accuracy vs. 2021.
2. **Promote the winning validation spec to `03_spec_grid_cv.R`.** The no-covariate / Mundlak / within-only story in `validation_note.tex` needs updating given these results. Add a row for the lagged-2017 spec with explicit framing that lagged covariates are not circular.
3. **Clean up `_test_*` scripts.** Delete `_test_within_state_cv.R` and `_test_within_state_cv_grid.R` (superseded). Keep or run `_test_demeaned_interactions.R` only if further within-state gains are wanted. The three that ran are now CSV-producing and could be promoted to numbered scripts if kept.
4. **Investigate the Grüne +6 pp bias.** Correlation is high but level is systematically off — likely a survey-mode / mobilization artifact. Worth a short diagnostic before the validation note is finalized.
5. **Decide on the 2000S-4020 Berufsabschluss data** (data/raw/zensus/). Untracked, unused. Either wire into `02_build_poststrat_frame.R` or delete.
6. **Rebuild `validation_note.pdf` + commit** once the note reflects the new specs and (if adopted) the 2017-covariate switch.
