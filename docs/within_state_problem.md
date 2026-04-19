# The within-state correlation problem

## What's going on

County-level MRP vote-share predictions correlate with BTW 2021 actuals at median r = 0.87 overall. Once state-level variation is stripped out (demean both predictions and actuals by state), the median drops to 0.66. Three parties collapse much more sharply.

| Party | Overall r | Within-state r | Within-state SD in actual data |
|-------|-----------|----------------|---------------------------------|
| AfD | 0.97 | 0.95 | 2.7 pp |
| CDU/CSU | 0.88 | 0.90 | 3.8 pp |
| Grüne | 0.86 | 0.95 | 4.5 pp |
| DIE LINKE | 0.89 | 0.24 | 1.2 pp |
| SPD | 0.79 | 0.41 | 3.2 pp |
| FDP | 0.72 | 0.42 | 1.4 pp |

## Why it matters

The headline overall r is flattered by 16 state random effects. For 400 counties, a useful model should differentiate counties within a state, not just sort counties across states. Median within-state r = 0.66 says the model does this only modestly; for DIE LINKE, SPD, and FDP it is barely better than guessing.

## Why the three weak parties are weak

- **DIE LINKE (0.24)** and **FDP (0.42)**: within-state SD in actual votes is 1.2 pp and 1.4 pp. Tiny signal; any model fights noise. Both parties differentiate east/west but are close to uniform within a state. Low within-state r here is partly inherent.
- **SPD (0.41)**: this one is genuine model weakness. Within-state SD = 3.2 pp, comparable to CDU (3.8 pp, r = 0.90). The signal exists; the model is not capturing it.

## Why covariates hurt (Mundlak problem)

Pooled fixed effects on county covariates like unemployment conflate between-state and within-state relationships. East German states have both high unemployment and distinctive voting patterns, so the pooled coefficient borrows east/west sorting. For a held-out county the model then pushes within-state predictions toward "how east-German does this county look" rather than local deviation. Adding unemployment lifts overall r (0.811 to 0.816) but drops within-state r (0.60 to 0.59).

## Fixes that work (strict 5-fold holdout, 2026-04-19)

| Spec | Overall r | Within-state r |
|------|-----------|----------------|
| Status quo (no covariates, Bezirk RE) | 0.811 | 0.602 |
| Mundlak / within-only covariates | 0.810 | 0.625 |
| 2017–2021 window + within-state FE | 0.833 | 0.654 |
| Within-state 2017 vote shares | 0.866 | 0.751 |
| **Pooled 2017 vote shares** | **0.899** | **0.868** |

Lagged 2017 county vote shares almost eliminate the problem. Per party (strict holdout): SPD 0.20 to 0.87, CDU 0.68 to 0.91, AfD 0.49 to 0.92, Grüne 0.85 to 0.89. DIE LINKE does not improve (inherent variance issue).

## The production-model implication

The production model (`03_load_covariates.R:26`, `04b_fit_all_lme4.R:185`) uses **BTW 2021** results as covariates. Validating 2021 predictions against 2021 covariates is the circularity the validation note was explicitly designed to avoid, which is why the note strips all election covariates and reports the much weaker within-state numbers above.

Switching production covariates to **BTW 2017** (lagged, not circular) recovers within-state r of 0.87 for most parties. This is the single biggest available improvement, and it also restores the validation story's validity.

## Bottom line

Two layers:

1. **Inherent** low within-state r for DIE LINKE and partly FDP. Their signal is mostly between-state.
2. **Fixable** low within-state r for SPD, CDU, AfD, Grüne. Current covariate specification pushes them toward east/west sorting.

Moving from BTW 2021 to BTW 2017 as the production covariate block is the most defensible fix and closes most of the within-state gap. Dropping the between-state component of socioeconomic covariates (Mundlak-style) is a smaller but independently useful change. Nothing else tested comes close.
