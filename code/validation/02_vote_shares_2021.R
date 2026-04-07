# ==============================================================================
# 02_vote_shares_2021.R
# Fit MRP models for vote_* issues using ONLY 2021 survey data,
# WITHOUT election covariates, then validate against BTW 2021 actual results.
# Runs two specs: (A) demographics + pop density, (B) + INKAR covariates.
# ==============================================================================

set.seed(20260407)

library(tidyverse)
library(lme4)
library(future.apply)

mb_root    <- here::here()
gerda_root <- Sys.getenv("GERDA_ROOT",
                         unset = normalizePath(file.path(here::here(), "..", "german_election_data"),
                                               mustWork = FALSE))
output_dir <- file.path(mb_root, "output", "validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plan(multisession, workers = max(1L, parallelly::availableCores() - 1L))
options(future.globals.maxSize = 1024 * 1024^2)

# ---- 1. Load data ------------------------------------------------------------

survey    <- readRDS(file.path(mb_root, "data", "harmonized", "survey_pooled.rds"))
kreis_cov <- readRDS(file.path(mb_root, "data", "covariates", "kreis_covariates.rds"))
poststrat <- readRDS(file.path(mb_root, "data", "poststrat", "poststrat_kreis.rds"))

vote_issues <- c("vote_cdu", "vote_spd", "vote_gruene",
                  "vote_fdp", "vote_afd", "vote_linke")

survey_2021 <- survey |>
  filter(year == 2021, issue_id %in% vote_issues)

message("Survey 2021 vote data: ", nrow(survey_2021), " obs, ",
        n_distinct(survey_2021$respondent_id), " respondents")

# ---- 2. Prepare covariates ---------------------------------------------------

kreis_cov_std <- kreis_cov |>
  mutate(
    log_pop_density_z  = scale(log(cty_pop_density + 1))[, 1],
    foreigner_share_z  = scale(Ausländeranteil_inkar)[, 1],
    median_income_z    = scale(Medianeinkommen_inkar)[, 1],
    rent_z             = scale(Mietpreise_inkar)[, 1],
    refugee_share_z    = scale(Schutzsuchende_an_Bevölkerung_inkar)[, 1]
  )

inkar_vars <- c("foreigner_share_z", "median_income_z", "rent_z", "refugee_share_z")

# ---- 3. Build poststrat prediction data --------------------------------------

pred_kreis <- poststrat |>
  left_join(kreis_cov_std |> select(county_code, log_pop_density_z, all_of(inkar_vars)),
            by = "county_code") |>
  mutate(male = as.integer(male), wkr_nr = NA_character_)

# ---- 4. Fit + poststratify function ------------------------------------------

fit_and_poststratify_kreis <- function(issue, survey_data, use_inkar = FALSE) {
  d <- survey_data |>
    filter(issue_id == !!issue) |>
    drop_na(y, age_cat, male, educ_label, state_code)

  if (nrow(d) < 300) {
    message("  Skipping ", issue, ": only ", nrow(d), " obs")
    return(NULL)
  }

  cov_cols <- c("county_code", "log_pop_density_z")
  if (use_inkar) cov_cols <- c(cov_cols, inkar_vars)

  d <- d |>
    left_join(kreis_cov_std |> select(all_of(cov_cols)), by = "county_code") |>
    mutate(
      across(any_of(c("log_pop_density_z", inkar_vars)), ~ replace_na(.x, 0)),
      age_cat       = factor(age_cat, levels = c("18-29","30-44","45-59","60-74","75+")),
      educ_label    = factor(educ_label, levels = c("no_degree","hauptschule","realschule",
                                                     "abitur","university")),
      state_code    = factor(state_code),
      survey_source = factor(survey_source),
      legperiod     = factor(legperiod),
      county_code   = ifelse(is.na(county_code), "missing", county_code),
      wkr_nr        = as.character(wkr_nr)
    ) |>
    filter(!is.na(wkr_nr)) |>
    droplevels()

  if (nrow(d) < 300) {
    message("  Skipping ", issue, " after filtering: only ", nrow(d), " obs")
    return(NULL)
  }

  n_kreise  <- n_distinct(d$county_code[d$county_code != "missing"])
  n_wkr     <- n_distinct(d$wkr_nr)
  n_sources <- nlevels(d$survey_source)
  n_legper  <- nlevels(d$legperiod)
  n_states  <- nlevels(d$state_code)

  if (use_inkar) {
    fe <- "y ~ male + log_pop_density_z + foreigner_share_z + median_income_z + rent_z + refugee_share_z"
  } else {
    fe <- "y ~ male + log_pop_density_z"
  }

  re <- "(1 | age_cat) + (1 | educ_label) + (1 | educ_label:age_cat) + (1 | male:age_cat) + (1 | male:educ_label)"
  if (n_sources > 1) re <- paste0(re, " + (1 | survey_source)")
  if (n_legper > 1)  re <- paste0(re, " + (1 | legperiod)")
  if (n_states > 1)  re <- paste0(re, " + (1 | state_code)")
  if (n_kreise > 10) re <- paste0(re, " + (1 | county_code)")
  if (n_wkr > 10)    re <- paste0(re, " + (1 | wkr_nr)")

  formula_str <- paste(fe, "+", re)
  message("  N=", nrow(d), " | Kreise=", n_kreise, " | inkar=", use_inkar)

  fit <- tryCatch(
    glmer(as.formula(formula_str), data = d, family = binomial(link = "logit"),
          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5),
                                 calc.derivs = FALSE),
          nAGQ = 0),
    error = function(e) {
      message("  First attempt failed: ", e$message)
      formula_fallback <- gsub(" \\+ \\(1 \\| wkr_nr\\)", "", formula_str)
      message("  Retrying without (1 | wkr_nr)...")
      tryCatch(
        glmer(as.formula(formula_fallback), data = d, family = binomial(link = "logit"),
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5),
                                     calc.derivs = FALSE),
              nAGQ = 0),
        error = function(e2) { message("  ERROR: ", e2$message); NULL }
      )
    }
  )

  if (is.null(fit)) return(NULL)

  pk <- pred_kreis |>
    mutate(
      survey_source = factor(levels(d$survey_source)[1], levels = levels(d$survey_source)),
      legperiod     = factor(tail(sort(unique(as.character(d$legperiod))), 1),
                             levels = levels(d$legperiod)),
      county_code   = county_code,
      wkr_nr        = NA_character_
    )

  pk$.pred <- predict(fit, newdata = pk, type = "response", allow.new.levels = TRUE)

  est_kreis <- pk |>
    group_by(county_code) |>
    summarise(estimate = weighted.mean(.pred, N, na.rm = TRUE),
              pop = sum(N), .groups = "drop")

  list(kreis = est_kreis)
}

# ---- 5. Load GERDA ground truth (BTW 2021) -----------------------------------

fed_cty <- readRDS(file.path(gerda_root, "data", "federal_elections",
                              "county_level", "final", "federal_cty_harm.rds"))

btw21_long <- fed_cty |>
  filter(election_year == 2021) |>
  select(county_code, cdu_csu, spd, gruene, fdp, afd, linke_pds) |>
  pivot_longer(-county_code, names_to = "party_col", values_to = "actual") |>
  mutate(issue_id = case_when(
    party_col == "cdu_csu"   ~ "vote_cdu",
    party_col == "spd"       ~ "vote_spd",
    party_col == "gruene"    ~ "vote_gruene",
    party_col == "fdp"       ~ "vote_fdp",
    party_col == "afd"       ~ "vote_afd",
    party_col == "linke_pds" ~ "vote_linke"
  )) |>
  select(county_code, issue_id, actual) |>
  mutate(county_code = as.character(county_code))

party_labels <- c(vote_cdu = "CDU/CSU", vote_spd = "SPD", vote_gruene = "Grüne",
                  vote_fdp = "FDP", vote_afd = "AfD", vote_linke = "DIE LINKE")

# ---- 6. Helper: run a spec and validate --------------------------------------

run_spec <- function(spec_name, use_inkar) {
  message("\n=== ", spec_name, " ===\n")

  results <- future_lapply(vote_issues, function(issue) {
    fit_and_poststratify_kreis(issue, survey_2021, use_inkar = use_inkar)
  }, future.seed = TRUE)
  names(results) <- vote_issues

  successful <- names(compact(results))
  message(length(successful), " / ", length(vote_issues), " fitted")

  mrp_kreis <- map_dfr(successful, ~ mutate(results[[.x]]$kreis, issue_id = .x)) |>
    mutate(county_code = as.character(county_code))

  val <- inner_join(mrp_kreis, btw21_long, by = c("county_code", "issue_id"))
  if (n_distinct(val$county_code) < 100) {
    val2 <- inner_join(
      mrp_kreis |> mutate(county_code = str_pad(county_code, 5, pad = "0")),
      btw21_long |> mutate(county_code = str_pad(county_code, 5, pad = "0")),
      by = c("county_code", "issue_id"))
    if (nrow(val2) > nrow(val)) val <- val2
  }

  metrics <- val |>
    group_by(issue_id) |>
    summarise(
      n_counties  = n(),
      r           = cor(estimate, actual, use = "complete.obs"),
      rmse_pp     = sqrt(mean((estimate - actual)^2, na.rm = TRUE)) * 100,
      bias_pp     = mean(estimate - actual, na.rm = TRUE) * 100,
      mean_mrp    = mean(estimate, na.rm = TRUE) * 100,
      mean_actual = mean(actual, na.rm = TRUE) * 100,
      .groups     = "drop"
    ) |>
    mutate(party = party_labels[issue_id], spec = spec_name) |>
    arrange(desc(r)) |>
    select(spec, party, issue_id, n_counties, r, rmse_pp, bias_pp, mean_mrp, mean_actual)

  list(metrics = metrics, val = val)
}

# ---- 7. glmmLasso covariate selection spec -----------------------------------

library(glmmLasso)

select_and_fit_lasso <- function(issue, survey_data) {
  d <- survey_data |>
    filter(issue_id == !!issue) |>
    drop_na(y, age_cat, male, educ_label, state_code)

  if (nrow(d) < 300) {
    message("  Skipping ", issue, ": only ", nrow(d), " obs")
    return(NULL)
  }

  d <- d |>
    left_join(kreis_cov_std |> select(county_code, log_pop_density_z, all_of(inkar_vars)),
              by = "county_code") |>
    mutate(
      across(any_of(c("log_pop_density_z", inkar_vars)), ~ replace_na(.x, 0)),
      age_cat       = factor(age_cat, levels = c("18-29","30-44","45-59","60-74","75+")),
      educ_label    = factor(educ_label, levels = c("no_degree","hauptschule","realschule",
                                                     "abitur","university")),
      state_code    = factor(state_code),
      survey_source = factor(survey_source),
      legperiod     = factor(legperiod),
      county_code   = ifelse(is.na(county_code), "missing", county_code),
      wkr_nr        = as.character(wkr_nr)
    ) |>
    filter(!is.na(wkr_nr)) |>
    droplevels()

  if (nrow(d) < 300) return(NULL)

  # glmmLasso: L1-penalize fixed effects, keep state_code as random effect
  # (glmmLasso supports fewer RE than glmer, so we use state_code only for selection)
  candidate_covs <- c("log_pop_density_z", inkar_vars)

  # BIC grid search over lambda
  lambdas <- c(500, 200, 100, 50, 20, 10, 5, 2, 1)
  best_bic <- Inf
  best_lambda <- lambdas[1]

  for (lam in lambdas) {
    fit_l <- tryCatch(
      glmmLasso(
        fix = as.formula(paste("y ~ male +", paste(candidate_covs, collapse = " + "))),
        rnd = list(state_code = ~1),
        data = d,
        lambda = lam,
        family = binomial(link = "logit"),
        control = list(print.iter = FALSE)
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_l)) {
      bic_val <- fit_l$bic
      if (bic_val < best_bic) {
        best_bic <- bic_val
        best_lambda <- lam
      }
    }
  }

  # Refit at best lambda
  fit_best <- tryCatch(
    glmmLasso(
      fix = as.formula(paste("y ~ male +", paste(candidate_covs, collapse = " + "))),
      rnd = list(state_code = ~1),
      data = d,
      lambda = best_lambda,
      family = binomial(link = "logit"),
      control = list(print.iter = FALSE)
    ),
    error = function(e) { message("  glmmLasso failed: ", e$message); NULL }
  )

  if (is.null(fit_best)) return(NULL)

  # Extract selected covariates (non-zero fixed effects, excluding intercept and male)
  coefs <- fit_best$coefficients
  selected <- names(coefs)[coefs != 0 & !names(coefs) %in% c("(Intercept)", "male")]
  message("  lambda=", best_lambda, " | BIC=", round(best_bic),
          " | selected: ", paste(selected, collapse = ", "))

  # Refit with glmer using only selected covariates (for proper poststratification)
  if (length(selected) == 0) {
    fe <- "y ~ male"
  } else {
    fe <- paste("y ~ male +", paste(selected, collapse = " + "))
  }

  n_kreise  <- n_distinct(d$county_code[d$county_code != "missing"])
  n_wkr     <- n_distinct(d$wkr_nr)
  n_sources <- nlevels(d$survey_source)
  n_legper  <- nlevels(d$legperiod)
  n_states  <- nlevels(d$state_code)

  re <- "(1 | age_cat) + (1 | educ_label) + (1 | educ_label:age_cat) + (1 | male:age_cat) + (1 | male:educ_label)"
  if (n_sources > 1) re <- paste0(re, " + (1 | survey_source)")
  if (n_legper > 1)  re <- paste0(re, " + (1 | legperiod)")
  if (n_states > 1)  re <- paste0(re, " + (1 | state_code)")
  if (n_kreise > 10) re <- paste0(re, " + (1 | county_code)")
  if (n_wkr > 10)    re <- paste0(re, " + (1 | wkr_nr)")

  formula_str <- paste(fe, "+", re)

  fit <- tryCatch(
    glmer(as.formula(formula_str), data = d, family = binomial(link = "logit"),
          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5),
                                 calc.derivs = FALSE),
          nAGQ = 0),
    error = function(e) {
      formula_fallback <- gsub(" \\+ \\(1 \\| wkr_nr\\)", "", formula_str)
      tryCatch(
        glmer(as.formula(formula_fallback), data = d, family = binomial(link = "logit"),
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5),
                                     calc.derivs = FALSE),
              nAGQ = 0),
        error = function(e2) NULL
      )
    }
  )

  if (is.null(fit)) return(NULL)

  pk <- pred_kreis |>
    mutate(
      survey_source = factor(levels(d$survey_source)[1], levels = levels(d$survey_source)),
      legperiod     = factor(tail(sort(unique(as.character(d$legperiod))), 1),
                             levels = levels(d$legperiod)),
      county_code   = county_code,
      wkr_nr        = NA_character_
    )

  pk$.pred <- predict(fit, newdata = pk, type = "response", allow.new.levels = TRUE)

  est_kreis <- pk |>
    group_by(county_code) |>
    summarise(estimate = weighted.mean(.pred, N, na.rm = TRUE),
              pop = sum(N), .groups = "drop")

  list(kreis = est_kreis, selected = selected, lambda = best_lambda)
}

# ---- 8. Run all three specs --------------------------------------------------

res_base  <- run_spec("baseline",       use_inkar = FALSE)
res_inkar <- run_spec("baseline_inkar", use_inkar = TRUE)

# glmmLasso spec (sequential — glmmLasso not compatible with future_lapply)
message("\n=== lasso_select ===\n")
lasso_results <- lapply(vote_issues, function(issue) {
  message("  ", issue)
  select_and_fit_lasso(issue, survey_2021)
})
names(lasso_results) <- vote_issues

successful_lasso <- names(compact(lasso_results))
message(length(successful_lasso), " / ", length(vote_issues), " fitted")

# Report selected covariates per issue
message("\nSelected covariates per issue:")
for (issue in successful_lasso) {
  sel <- lasso_results[[issue]]$selected
  message("  ", issue, ": ", if (length(sel) == 0) "(none)" else paste(sel, collapse = ", "))
}

mrp_lasso <- map_dfr(successful_lasso, ~ mutate(lasso_results[[.x]]$kreis, issue_id = .x)) |>
  mutate(county_code = as.character(county_code))

val_lasso <- inner_join(mrp_lasso, btw21_long, by = c("county_code", "issue_id"))
if (n_distinct(val_lasso$county_code) < 100) {
  val2 <- inner_join(
    mrp_lasso |> mutate(county_code = str_pad(county_code, 5, pad = "0")),
    btw21_long |> mutate(county_code = str_pad(county_code, 5, pad = "0")),
    by = c("county_code", "issue_id"))
  if (nrow(val2) > nrow(val_lasso)) val_lasso <- val2
}

metrics_lasso <- val_lasso |>
  group_by(issue_id) |>
  summarise(
    n_counties  = n(),
    r           = cor(estimate, actual, use = "complete.obs"),
    rmse_pp     = sqrt(mean((estimate - actual)^2, na.rm = TRUE)) * 100,
    bias_pp     = mean(estimate - actual, na.rm = TRUE) * 100,
    mean_mrp    = mean(estimate, na.rm = TRUE) * 100,
    mean_actual = mean(actual, na.rm = TRUE) * 100,
    .groups     = "drop"
  ) |>
  mutate(party = party_labels[issue_id], spec = "lasso_select") |>
  arrange(desc(r)) |>
  select(spec, party, issue_id, n_counties, r, rmse_pp, bias_pp, mean_mrp, mean_actual)

res_lasso <- list(metrics = metrics_lasso, val = val_lasso)

# ---- 9. Print comparison -----------------------------------------------------

combined <- bind_rows(res_base$metrics, res_inkar$metrics, res_lasso$metrics)

message("\n", strrep("=", 70))
message("COMPARISON: baseline vs. +INKAR vs. LASSO-selected")
message(strrep("=", 70))
message(sprintf("\n%-12s  %s  %s  %s",
                "Party", "r(base)", "r(+INKAR)", "r(LASSO)"))
for (p in sort(unique(combined$party))) {
  r_b <- combined$r[combined$party == p & combined$spec == "baseline"]
  r_i <- combined$r[combined$party == p & combined$spec == "baseline_inkar"]
  r_l <- combined$r[combined$party == p & combined$spec == "lasso_select"]
  if (length(r_b) == 1 && length(r_i) == 1 && length(r_l) == 1) {
    message(sprintf("%-12s  %.3f    %.3f      %.3f", p, r_b, r_i, r_l))
  }
}

med_b <- median(res_base$metrics$r)
med_i <- median(res_inkar$metrics$r)
med_l <- median(res_lasso$metrics$r)
message(sprintf("\n%-12s  %.3f    %.3f      %.3f", "MEDIAN", med_b, med_i, med_l))

med_rmse_b <- median(res_base$metrics$rmse_pp)
med_rmse_i <- median(res_inkar$metrics$rmse_pp)
med_rmse_l <- median(res_lasso$metrics$rmse_pp)
message(sprintf("%-12s  %.1fpp    %.1fpp      %.1fpp", "MEDIAN RMSE", med_rmse_b, med_rmse_i, med_rmse_l))

# ---- 10. Save outputs --------------------------------------------------------

write_csv(combined, file.path(output_dir, "vote_share_validation_2021only.csv"))
message("\nSaved: ", file.path(output_dir, "vote_share_validation_2021only.csv"))

# Scatter plot: all three specs
spec_labels <- c(baseline = "Baseline (pop density only)",
                 baseline_inkar = "All INKAR covariates",
                 lasso_select = "LASSO-selected covariates")

scatter_data <- bind_rows(
  res_base$val |> mutate(spec = "Baseline (pop density only)"),
  res_inkar$val |> mutate(spec = "All INKAR covariates"),
  res_lasso$val |> mutate(spec = "LASSO-selected covariates")
) |>
  mutate(estimate_pct = estimate * 100, actual_pct = actual * 100,
         party = party_labels[issue_id],
         spec = factor(spec, levels = c("Baseline (pop density only)",
                                         "All INKAR covariates",
                                         "LASSO-selected covariates")))

label_data <- combined |>
  mutate(party = party_labels[issue_id],
         spec = factor(spec_labels[spec], levels = levels(scatter_data$spec)),
         label = paste0("r = ", round(r, 3), "\nRMSE = ", round(rmse_pp, 1), "pp"))

p <- ggplot(scatter_data, aes(x = actual_pct, y = estimate_pct)) +
  geom_point(alpha = 0.25, size = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "steelblue") +
  facet_grid(spec ~ party, scales = "free") +
  geom_text(data = label_data,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.3, size = 2.2, color = "grey30",
            inherit.aes = FALSE) +
  labs(x = "Actual vote share (BTW 2021, %)",
       y = "MRP predicted vote share (%)",
       title = "MRP vote share validation (2021 only): baseline vs. INKAR vs. LASSO") +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(output_dir, "vote_share_scatter_2021only.pdf"),
       p, width = 14, height = 9)
message("Saved: ", file.path(output_dir, "vote_share_scatter_2021only.pdf"))
