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
    foreigner_share_z  = scale(share_foreign)[, 1],
    median_income_z    = scale(median_income)[, 1],
    unemployment_z     = scale(unemployment_rate)[, 1],
    share_65plus_z     = scale(share_65plus)[, 1],
    gdp_per_capita_z   = scale(gdp_per_capita)[, 1],
    share_abitur_z     = scale(share_abitur)[, 1],
    purchasing_power_z = scale(purchasing_power)[, 1],
    share_secondary_z  = scale(share_secondary_sector)[, 1],
    urban              = as.integer(log_pop_density_z > 0)
  )

inkar_vars <- c("foreigner_share_z", "median_income_z", "unemployment_z",
                "share_65plus_z", "gdp_per_capita_z", "share_abitur_z",
                "purchasing_power_z", "share_secondary_z")

# ---- 3. Build poststrat prediction data --------------------------------------

# 100-cell frame (with employment) for specs that use it
pred_kreis_100 <- poststrat |>
  left_join(kreis_cov_std |> select(county_code, log_pop_density_z, all_of(inkar_vars), urban),
            by = "county_code") |>
  mutate(male = as.integer(male), employed = as.integer(employed), wkr_nr = NA_character_)

# 50-cell frame (without employment) for baseline specs
pred_kreis_50 <- poststrat |>
  group_by(county_code, age_cat, male, educ_label, state_code) |>
  summarise(N = sum(N), .groups = "drop") |>
  left_join(kreis_cov_std |> select(county_code, log_pop_density_z, all_of(inkar_vars), urban),
            by = "county_code") |>
  mutate(male = as.integer(male), wkr_nr = NA_character_)

# Default pred_kreis for existing specs (50-cell)
pred_kreis <- pred_kreis_50

# ---- 4. Fit + poststratify function ------------------------------------------

fit_and_poststratify_kreis <- function(issue, survey_data, use_inkar = FALSE,
                                       use_employed = FALSE, extra_covs = NULL,
                                       use_urban_ix = FALSE, drop_wkr = FALSE) {
  drop_vars <- c("y", "age_cat", "male", "educ_label", "state_code")
  if (use_employed) drop_vars <- c(drop_vars, "employed")

  d <- survey_data |>
    filter(issue_id == !!issue) |>
    drop_na(all_of(drop_vars))

  if (nrow(d) < 300) {
    message("  Skipping ", issue, ": only ", nrow(d), " obs")
    return(NULL)
  }

  all_covs <- unique(c("log_pop_density_z",
                        if (use_inkar) inkar_vars,
                        if (use_urban_ix) "urban",
                        extra_covs))
  cov_cols <- c("county_code", all_covs)

  d <- d |>
    left_join(kreis_cov_std |> select(all_of(cov_cols)), by = "county_code") |>
    mutate(
      across(any_of(all_covs), ~ replace_na(.x, 0)),
      age_cat       = factor(age_cat, levels = c("18-29","30-44","45-59","60-74","75+")),
      educ_label    = factor(educ_label, levels = c("no_degree","hauptschule","realschule",
                                                     "abitur","university")),
      state_code    = factor(state_code),
      survey_source = factor(survey_source),
      legperiod     = factor(legperiod),
      county_code   = ifelse(is.na(county_code), "missing", county_code),
      wkr_nr        = as.character(wkr_nr)
    )
  if (use_employed) d <- d |> mutate(employed = as.integer(employed))
  d <- d |>
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

  fe <- "y ~ male + log_pop_density_z"
  if (use_inkar) fe <- paste(fe, "+", paste(inkar_vars, collapse = " + "))
  if (!is.null(extra_covs)) fe <- paste(fe, "+", paste(extra_covs, collapse = " + "))
  if (use_employed) fe <- paste0(fe, " + employed")

  re <- "(1 | age_cat) + (1 | educ_label) + (1 | educ_label:age_cat) + (1 | male:age_cat) + (1 | male:educ_label)"
  if (use_employed) re <- paste0(re, " + (1 | employed:age_cat)")
  if (use_urban_ix) re <- paste0(re, " + (1 | educ_label:urban) + (1 | age_cat:urban)")
  if (n_sources > 1) re <- paste0(re, " + (1 | survey_source)")
  if (n_legper > 1)  re <- paste0(re, " + (1 | legperiod)")
  if (n_states > 1)  re <- paste0(re, " + (1 | state_code)")
  if (n_kreise > 10) re <- paste0(re, " + (1 | county_code)")
  if (n_wkr > 10 && !drop_wkr) re <- paste0(re, " + (1 | wkr_nr)")

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

  pred_df <- if (use_employed) pred_kreis_100 else pred_kreis_50
  pk <- pred_df |>
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

run_spec <- function(spec_name, use_inkar = FALSE, use_employed = FALSE,
                     extra_covs = NULL, use_urban_ix = FALSE, drop_wkr = FALSE,
                     survey_data = survey_2021) {
  message("\n=== ", spec_name, " ===\n")

  results <- future_lapply(vote_issues, function(issue) {
    fit_and_poststratify_kreis(issue, survey_data, use_inkar = use_inkar,
                               use_employed = use_employed, extra_covs = extra_covs,
                               use_urban_ix = use_urban_ix, drop_wkr = drop_wkr)
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

# ---- 7. Forward selection on end-to-end validation metric --------------------

# Core forward selection engine. If train_counties is NULL, uses all counties.
# Otherwise, selection criterion uses only train_counties; predictions for all 400.
forward_select_engine <- function(survey_data, candidate_covs, use_urban_ix,
                                  train_counties = NULL) {

  eval_covs <- function(covs) {
    results <- future_lapply(vote_issues, function(issue) {
      fit_and_poststratify_kreis(issue, survey_data, extra_covs = covs,
                                 use_urban_ix = use_urban_ix)
    }, future.seed = TRUE)
    names(results) <- vote_issues
    successful <- names(compact(results))
    if (length(successful) < 4) return(list(median_r = -Inf, results = results))

    mrp <- map_dfr(successful, ~ mutate(results[[.x]]$kreis, issue_id = .x)) |>
      mutate(county_code = as.character(county_code))
    val <- inner_join(mrp, btw21_long, by = c("county_code", "issue_id"))

    # Evaluate only on training counties if specified
    val_eval <- if (!is.null(train_counties)) {
      val |> filter(county_code %in% train_counties)
    } else val

    party_r <- val_eval |>
      group_by(issue_id) |>
      summarise(r = cor(estimate, actual, use = "complete.obs"), .groups = "drop")

    list(median_r = median(party_r$r), results = results)
  }

  baseline <- eval_covs(NULL)
  best_r <- baseline$median_r
  best_results <- baseline$results
  selected <- character(0)
  message("  Baseline median r: ", round(best_r, 4))

  remaining <- candidate_covs
  round_num <- 0

  repeat {
    round_num <- round_num + 1
    if (length(remaining) == 0) break

    message("\n  Round ", round_num, ": testing ", length(remaining), " candidates")
    round_results <- list()

    for (cov in remaining) {
      test_covs <- c(selected, cov)
      res <- eval_covs(test_covs)
      round_results[[cov]] <- res
      message("    + ", cov, ": median r = ", round(res$median_r, 4))
    }

    round_r <- sapply(round_results, function(x) x$median_r)
    best_candidate <- names(which.max(round_r))
    best_candidate_r <- max(round_r)

    if (best_candidate_r > best_r) {
      selected <- c(selected, best_candidate)
      best_r <- best_candidate_r
      best_results <- round_results[[best_candidate]]$results
      remaining <- setdiff(remaining, best_candidate)
      message("  => Selected: ", best_candidate, " (median r: ", round(best_r, 4), ")")
    } else {
      message("  => No improvement. Stopping.")
      break
    }
  }

  list(selected = selected, best_r = best_r, results = best_results)
}

# Non-CV forward selection (original behavior)
forward_select_spec <- function(spec_name, survey_data, candidate_covs = inkar_vars,
                                use_urban_ix = FALSE) {
  message("\n=== ", spec_name, " (forward selection) ===\n")

  fs <- forward_select_engine(survey_data, candidate_covs, use_urban_ix)

  message("\nFinal selected covariates: ",
          if (length(fs$selected) == 0) "(none)" else paste(fs$selected, collapse = ", "))
  message("Final median r: ", round(fs$best_r, 4))

  # Build validation metrics from all 400 counties
  successful <- names(compact(fs$results))
  mrp <- map_dfr(successful, ~ mutate(fs$results[[.x]]$kreis, issue_id = .x)) |>
    mutate(county_code = as.character(county_code))
  val <- inner_join(mrp, btw21_long, by = c("county_code", "issue_id"))

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

  list(metrics = metrics, val = val, selected = fs$selected)
}

# K-fold spatial CV forward selection: select covariates on K-1 folds of counties,
# evaluate on held-out fold. Reports honest out-of-sample performance.
forward_select_cv <- function(spec_name, survey_data, candidate_covs = inkar_vars,
                              use_urban_ix = FALSE, K = 5) {
  message("\n=== ", spec_name, " (", K, "-fold spatial CV) ===\n")

  # Split 400 ground-truth counties into K folds
  all_counties <- sort(unique(btw21_long$county_code))
  set.seed(20260409)
  fold_ids <- sample(rep(1:K, length.out = length(all_counties)))
  folds <- split(all_counties, fold_ids)

  fold_results <- list()

  for (k in seq_len(K)) {
    test_counties  <- folds[[k]]
    train_counties <- setdiff(all_counties, test_counties)
    message("\n--- Fold ", k, "/", K, " (train: ", length(train_counties),
            " counties, test: ", length(test_counties), ") ---")

    # Forward selection using only training counties for evaluation
    fs <- forward_select_engine(survey_data, candidate_covs, use_urban_ix,
                                train_counties = train_counties)

    message("  Fold ", k, " selected: ",
            if (length(fs$selected) == 0) "(none)" else paste(fs$selected, collapse = ", "))

    # Evaluate on held-out test counties
    successful <- names(compact(fs$results))
    mrp <- map_dfr(successful, ~ mutate(fs$results[[.x]]$kreis, issue_id = .x)) |>
      mutate(county_code = as.character(county_code))
    val_test <- inner_join(mrp, btw21_long, by = c("county_code", "issue_id")) |>
      filter(county_code %in% test_counties)

    test_r <- val_test |>
      group_by(issue_id) |>
      summarise(r = cor(estimate, actual, use = "complete.obs"), .groups = "drop")

    fold_results[[k]] <- list(
      selected = fs$selected,
      train_r  = fs$best_r,
      test_median_r = median(test_r$r),
      test_party_r = test_r,
      val_test = val_test
    )

    message("  Fold ", k, " train median r: ", round(fs$best_r, 4),
            " | test median r: ", round(median(test_r$r), 4))
  }

  # Aggregate across folds
  train_rs <- sapply(fold_results, function(x) x$train_r)
  test_rs  <- sapply(fold_results, function(x) x$test_median_r)
  all_selected <- lapply(fold_results, function(x) x$selected)

  message("\n=== CV Summary ===")
  message("Mean train median r: ", round(mean(train_rs), 4))
  message("Mean test median r:  ", round(mean(test_rs), 4),
          " (SD: ", round(sd(test_rs), 4), ")")
  message("Selected covariates per fold:")
  for (k in seq_len(K)) {
    sel <- all_selected[[k]]
    message("  Fold ", k, ": ", if (length(sel) == 0) "(none)" else paste(sel, collapse = ", "))
  }

  # Count how often each covariate is selected
  cov_counts <- table(unlist(all_selected))
  cov_counts <- sort(cov_counts, decreasing = TRUE)
  message("Covariate selection frequency (across ", K, " folds):")
  for (i in seq_along(cov_counts)) {
    message("  ", names(cov_counts)[i], ": ", cov_counts[i], "/", K)
  }

  # Build combined test-set validation (each county appears exactly once)
  val_all_test <- bind_rows(lapply(fold_results, function(x) x$val_test))

  metrics <- val_all_test |>
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

  list(metrics = metrics, val = val_all_test, selected = all_selected,
       cv_summary = tibble(fold = 1:K, train_r = train_rs, test_r = test_rs))
}

# ---- 8. Run all specs --------------------------------------------------------

# 2020-2021 window (more data, slight temporal mismatch)
survey_2020_2021 <- survey |>
  filter(year %in% 2020:2021, issue_id %in% vote_issues)
message("Survey 2020-2021 vote data: ", nrow(survey_2020_2021), " obs, ",
        n_distinct(survey_2020_2021$respondent_id), " respondents")

# 2019-2021 window
survey_2019_2021 <- survey |>
  filter(year %in% 2019:2021, issue_id %in% vote_issues)
message("Survey 2019-2021 vote data: ", nrow(survey_2019_2021), " obs, ",
        n_distinct(survey_2019_2021$respondent_id), " respondents")

res_3yr     <- run_spec("baseline_2019_2021", survey_data = survey_2019_2021)
res_2yr     <- run_spec("baseline_2020_2021", survey_data = survey_2020_2021)
res_2yr_emp <- run_spec("baseline_2020_2021_emp", use_employed = TRUE, survey_data = survey_2020_2021)
res_base    <- run_spec("baseline")
res_emp     <- run_spec("baseline_employed", use_employed = TRUE)
res_inkar   <- run_spec("baseline_inkar", use_inkar = TRUE)

# Forward selection specs
res_fwd     <- forward_select_spec("forward_select",              survey_2021)
res_fwd_2yr <- forward_select_spec("forward_select_2020_2021",    survey_2020_2021)
res_fwd_3yr <- forward_select_spec("forward_select_2019_2021",    survey_2019_2021)

# 5-fold spatial CV forward selection (best year window only)
res_fwd_3yr_cv <- forward_select_cv("forward_select_cv_2019_2021", survey_2019_2021)

# ---- 9. Print comparison -----------------------------------------------------

combined <- bind_rows(res_3yr$metrics,
                      res_2yr$metrics, res_2yr_emp$metrics,
                      res_base$metrics, res_emp$metrics, res_inkar$metrics,
                      res_fwd$metrics, res_fwd_2yr$metrics, res_fwd_3yr$metrics,
                      res_fwd_3yr_cv$metrics)

all_specs <- c("baseline_2019_2021",
               "forward_select_2019_2021", "forward_select_cv_2019_2021",
               "baseline_2020_2021", "baseline_2020_2021_emp", "forward_select_2020_2021",
               "baseline", "baseline_employed", "baseline_inkar", "forward_select")

message("\n", strrep("=", 80))
message("COMPARISON: all specs (sorted by median r)")
message(strrep("=", 80))

spec_summary <- combined |>
  group_by(spec) |>
  summarise(median_r = median(r), median_rmse = median(rmse_pp), .groups = "drop") |>
  arrange(desc(median_r))

message(sprintf("\n%-35s  %-8s  %-8s", "Spec", "Med r", "Med RMSE"))
for (i in seq_len(nrow(spec_summary))) {
  message(sprintf("%-35s  %.3f     %.1fpp",
                  spec_summary$spec[i], spec_summary$median_r[i], spec_summary$median_rmse[i]))
}

message("\nPer-party correlations (top 5 specs):")
top_specs <- spec_summary$spec[1:min(5, nrow(spec_summary))]
message(sprintf("%-12s  %s", "Party", paste(sprintf("%-25s", top_specs), collapse = "  ")))
for (p in sort(unique(combined$party))) {
  vals <- sapply(top_specs, function(s) {
    v <- combined$r[combined$party == p & combined$spec == s]
    if (length(v) == 1) sprintf("%.3f", v) else "  NA "
  })
  message(sprintf("%-12s  %s", p, paste(sprintf("%-25s", vals), collapse = "  ")))
}

# Report forward-selected covariates
fwd_specs <- list(
  "2021" = res_fwd, "2020-2021" = res_fwd_2yr, "2019-2021" = res_fwd_3yr
)
message("\nForward-selected covariates:")
for (nm in names(fwd_specs)) {
  sel <- fwd_specs[[nm]]$selected
  message("  ", nm, ": ", if (length(sel) == 0) "(none)" else paste(sel, collapse = ", "))
}
message("\nCV forward selection (2019-2021, 5-fold):")
message("  Mean test median r: ", round(mean(res_fwd_3yr_cv$cv_summary$test_r), 4),
        " (train: ", round(mean(res_fwd_3yr_cv$cv_summary$train_r), 4), ")")
for (k in seq_along(res_fwd_3yr_cv$selected)) {
  sel <- res_fwd_3yr_cv$selected[[k]]
  message("  Fold ", k, ": ", if (length(sel) == 0) "(none)" else paste(sel, collapse = ", "))
}

# ---- 10. Save outputs --------------------------------------------------------

write_csv(combined, file.path(output_dir, "vote_share_validation_2021only.csv"))
message("\nSaved: ", file.path(output_dir, "vote_share_validation_2021only.csv"))

# Scatter plot: all specs
spec_labels <- c(baseline_2019_2021 = "Baseline (2019-2021)",
                 forward_select_2019_2021 = "Fwd-selected (2019-2021)",
                 forward_select_cv_2019_2021 = "Fwd-selected CV (2019-2021)",
                 baseline_2020_2021 = "Baseline (2020-2021)",
                 baseline_2020_2021_emp = "+ Employment (2020-2021)",
                 forward_select_2020_2021 = "Fwd-selected (2020-2021)",
                 baseline = "Baseline (2021 only)",
                 baseline_employed = "+ Employment (2021)",
                 baseline_inkar = "All INKAR covariates",
                 forward_select = "Fwd-selected (2021)")

spec_order <- names(spec_labels)
scatter_data <- combined |>
  select(spec, issue_id) |>
  distinct() |>
  left_join(
    bind_rows(
      res_3yr$val          |> mutate(spec = "baseline_2019_2021"),
      res_fwd_3yr$val      |> mutate(spec = "forward_select_2019_2021"),
      res_fwd_3yr_cv$val   |> mutate(spec = "forward_select_cv_2019_2021"),
      res_2yr$val      |> mutate(spec = "baseline_2020_2021"),
      res_2yr_emp$val  |> mutate(spec = "baseline_2020_2021_emp"),
      res_fwd_2yr$val  |> mutate(spec = "forward_select_2020_2021"),
      res_base$val     |> mutate(spec = "baseline"),
      res_emp$val      |> mutate(spec = "baseline_employed"),
      res_inkar$val    |> mutate(spec = "baseline_inkar"),
      res_fwd$val      |> mutate(spec = "forward_select")
    ),
    by = c("spec", "issue_id")
  ) |>
  mutate(estimate_pct = estimate * 100, actual_pct = actual * 100,
         party = party_labels[issue_id],
         spec = factor(spec_labels[spec], levels = unname(spec_labels)))

label_data <- combined |>
  mutate(party = party_labels[issue_id],
         spec = factor(spec_labels[spec], levels = unname(spec_labels)),
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
       title = "MRP vote share validation: year window, covariates, and selection compared") +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(output_dir, "vote_share_scatter_2021only.pdf"),
       p, width = 14, height = 23)
message("Saved: ", file.path(output_dir, "vote_share_scatter_2021only.pdf"))
