# ==============================================================================
# 03_spec_grid_cv.R
# Systematic grid search over MRP specification space with 5-fold spatial CV.
# Evaluates 60 specs (state-varying slopes × county covariates × employment) on
# 2019-2021 survey data, validated against BTW 2021 county-level vote shares.
# ==============================================================================

set.seed(20260409)

library(tidyverse)
library(glmmTMB)
library(future.apply)
library(cli)

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

vote_issues <- c("vote_cdu", "vote_spd", "vote_gruene", "vote_fdp", "vote_afd", "vote_linke")
survey_data <- survey |> filter(year %in% 2019:2021, issue_id %in% vote_issues)
message("Survey 2019-2021: ", nrow(survey_data), " obs, ",
        n_distinct(survey_data$respondent_id), " respondents")

kreis_cov_std <- kreis_cov |>
  mutate(
    log_pop_density_z  = scale(log(cty_pop_density + 1))[, 1],
    unemployment_z     = scale(unemployment_rate)[, 1],
    share_abitur_z     = scale(share_abitur)[, 1],
    foreigner_share_z  = scale(share_foreign)[, 1],
    share_secondary_z  = scale(share_secondary_sector)[, 1]
  )

cov_select_cols <- c("county_code", "log_pop_density_z", "unemployment_z",
                     "share_abitur_z", "foreigner_share_z", "share_secondary_z")

# 50-cell prediction frame (baseline)
pred_kreis_50 <- poststrat |>
  group_by(county_code, age_cat, male, educ_label, state_code) |>
  summarise(N = sum(N), .groups = "drop") |>
  left_join(kreis_cov_std |> select(all_of(cov_select_cols)), by = "county_code") |>
  mutate(male = as.integer(male), wkr_nr = NA_character_)

# 100-cell prediction frame (employment specs)
pred_kreis_100 <- poststrat |>
  left_join(kreis_cov_std |> select(all_of(cov_select_cols)), by = "county_code") |>
  mutate(male = as.integer(male), employed = as.integer(employed), wkr_nr = NA_character_)

# Ground truth
fed_cty <- readRDS(file.path(gerda_root, "data", "federal_elections",
                              "county_level", "final", "federal_cty_harm.rds"))
btw21_long <- fed_cty |>
  filter(election_year == 2021) |>
  select(county_code, cdu_csu, spd, gruene, fdp, afd, linke_pds) |>
  pivot_longer(-county_code, names_to = "party_col", values_to = "actual") |>
  mutate(issue_id = case_when(
    party_col == "cdu_csu" ~ "vote_cdu", party_col == "spd" ~ "vote_spd",
    party_col == "gruene" ~ "vote_gruene", party_col == "fdp" ~ "vote_fdp",
    party_col == "afd" ~ "vote_afd", party_col == "linke_pds" ~ "vote_linke"
  )) |> select(county_code, issue_id, actual) |> mutate(county_code = as.character(county_code))

# ---- 2. Define specification grid --------------------------------------------

spec_grid <- expand.grid(
  state_slopes = c("none", "age", "educ", "both"),
  covariates   = c("none", "unemp", "unemp_abitur", "unemp_foreigner", "unemp_secondary"),
  employment   = c("none", "fe_only", "fe_and_re"),
  stringsAsFactors = FALSE
) |> as_tibble() |>
  mutate(spec_id = row_number())

message("\nSpecification grid: ", nrow(spec_grid), " specs")
print(spec_grid)

# ---- 3. Define CV folds ------------------------------------------------------

K <- 5
all_counties <- sort(unique(btw21_long$county_code))
fold_ids <- sample(rep(1:K, length.out = length(all_counties)))
folds <- split(all_counties, fold_ids)

# ---- 4. Fit function ---------------------------------------------------------

fit_issue <- function(issue, state_slopes, covariates, employment) {
  use_emp <- employment != "none"
  drop_vars <- c("y", "age_cat", "male", "educ_label", "state_code")
  if (use_emp) drop_vars <- c(drop_vars, "employed")

  d <- survey_data |>
    filter(issue_id == !!issue) |>
    drop_na(all_of(drop_vars))

  if (nrow(d) < 300) return(NULL)

  # Determine which county-level covariates to include
  cov_cols <- c("county_code", "log_pop_density_z")
  extra_fe <- character(0)
  if (grepl("unemp", covariates)) {
    cov_cols <- c(cov_cols, "unemployment_z")
    extra_fe <- c(extra_fe, "unemployment_z")
  }
  if (covariates == "unemp_abitur") {
    cov_cols <- c(cov_cols, "share_abitur_z")
    extra_fe <- c(extra_fe, "share_abitur_z")
  }
  if (covariates == "unemp_foreigner") {
    cov_cols <- c(cov_cols, "foreigner_share_z")
    extra_fe <- c(extra_fe, "foreigner_share_z")
  }
  if (covariates == "unemp_secondary") {
    cov_cols <- c(cov_cols, "share_secondary_z")
    extra_fe <- c(extra_fe, "share_secondary_z")
  }

  d <- d |>
    left_join(kreis_cov_std |> select(all_of(cov_cols)), by = "county_code") |>
    mutate(
      across(any_of(c("log_pop_density_z", "unemployment_z", "share_abitur_z",
                       "foreigner_share_z", "share_secondary_z")),
             ~ replace_na(.x, 0)),
      age_cat       = factor(age_cat, levels = c("18-29", "30-44", "45-59", "60-74", "75+")),
      educ_label    = factor(educ_label, levels = c("no_degree", "hauptschule", "realschule",
                                                     "abitur", "university")),
      state_code    = factor(state_code),
      survey_source = factor(survey_source),
      legperiod     = factor(legperiod),
      county_code   = ifelse(is.na(county_code), "missing", county_code),
      wkr_nr        = as.character(wkr_nr)
    )
  if (use_emp) d <- d |> mutate(employed = as.integer(employed))
  d <- d |> filter(!is.na(wkr_nr)) |> droplevels()

  if (nrow(d) < 300) return(NULL)

  # Build formula
  fe <- "y ~ male + log_pop_density_z"
  if (length(extra_fe) > 0) fe <- paste(fe, "+", paste(extra_fe, collapse = " + "))
  if (use_emp) fe <- paste0(fe, " + employed")

  re_parts <- c("(1|age_cat)", "(1|educ_label)", "(1|educ_label:age_cat)",
                "(1|male:age_cat)", "(1|male:educ_label)")
  if (employment == "fe_and_re") re_parts <- c(re_parts, "(1|employed:age_cat)")
  if (state_slopes %in% c("age", "both"))  re_parts <- c(re_parts, "(1|age_cat:state_code)")
  if (state_slopes %in% c("educ", "both")) re_parts <- c(re_parts, "(1|educ_label:state_code)")
  if (nlevels(d$survey_source) > 1) re_parts <- c(re_parts, "(1|survey_source)")
  if (nlevels(d$legperiod) > 1) re_parts <- c(re_parts, "(1|legperiod)")
  re_parts <- c(re_parts, "(1|state_code)", "(1|county_code)", "(1|wkr_nr)")

  formula_str <- paste(fe, "+", paste(re_parts, collapse = " + "))

  fit <- tryCatch(
    glmmTMB(as.formula(formula_str), data = d, family = binomial),
    error = function(e) {
      re_parts2 <- setdiff(re_parts, "(1|wkr_nr)")
      formula_str2 <- paste(fe, "+", paste(re_parts2, collapse = " + "))
      tryCatch(
        glmmTMB(as.formula(formula_str2), data = d, family = binomial),
        error = function(e2) NULL
      )
    }
  )

  if (is.null(fit)) return(NULL)

  # Poststratify
  pred_df <- if (use_emp) pred_kreis_100 else pred_kreis_50
  pk <- pred_df |> mutate(
    survey_source = factor(levels(d$survey_source)[1], levels = levels(d$survey_source)),
    legperiod     = factor(tail(sort(unique(as.character(d$legperiod))), 1),
                           levels = levels(d$legperiod)),
    wkr_nr        = NA_character_
  )
  pk$.pred <- predict(fit, newdata = pk, type = "response", allow.new.levels = TRUE)
  pk |> group_by(county_code) |>
    summarise(estimate = weighted.mean(.pred, N, na.rm = TRUE), .groups = "drop")
}

# ---- 5. Evaluate one spec on one fold ----------------------------------------

eval_spec_fold <- function(spec_row, train_counties, test_counties) {
  results <- future_lapply(vote_issues, function(issue) {
    fit_issue(issue,
              state_slopes = spec_row$state_slopes,
              covariates   = spec_row$covariates,
              employment   = spec_row$employment)
  }, future.seed = TRUE)
  names(results) <- vote_issues
  successful <- names(compact(results))
  if (length(successful) < 4) return(tibble(train_median_r = NA_real_, test_median_r = NA_real_))

  mrp <- map_dfr(successful, ~ mutate(results[[.x]], issue_id = .x)) |>
    mutate(county_code = as.character(county_code))
  val <- inner_join(mrp, btw21_long, by = c("county_code", "issue_id"))

  train_r <- val |> filter(county_code %in% train_counties) |>
    group_by(issue_id) |>
    summarise(r = cor(estimate, actual, use = "complete.obs"), .groups = "drop")

  test_r <- val |> filter(county_code %in% test_counties) |>
    group_by(issue_id) |>
    summarise(r = cor(estimate, actual, use = "complete.obs"), .groups = "drop")

  tibble(train_median_r = median(train_r$r), test_median_r = median(test_r$r))
}

# ---- 6. Run grid search with progress bar ------------------------------------

n_total <- nrow(spec_grid) * K
message("\nRunning ", nrow(spec_grid), " specs × ", K, " folds = ", n_total, " evaluations\n")

results_list <- list()
pb <- cli_progress_bar("Grid search CV", total = n_total, clear = FALSE,
                       format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}")

t_start <- Sys.time()

for (i in seq_len(nrow(spec_grid))) {
  spec_row <- spec_grid[i, ]
  spec_label <- paste0(
    ifelse(spec_row$age_x_state, "age×state", "no_age×state"), " | ",
    spec_row$covariates, " | ",
    ifelse(spec_row$employment, "employed", "no_empl")
  )

  for (k in seq_len(K)) {
    test_counties  <- folds[[k]]
    train_counties <- setdiff(all_counties, test_counties)

    res <- eval_spec_fold(spec_row, train_counties, test_counties)

    results_list[[length(results_list) + 1]] <- tibble(
      spec_id      = spec_row$spec_id,
      state_slopes = spec_row$state_slopes,
      covariates   = spec_row$covariates,
      employment   = spec_row$employment,
      fold         = k,
      train_median_r = res$train_median_r,
      test_median_r  = res$test_median_r
    )

    cli_progress_update(id = pb)
  }
}

cli_progress_done(id = pb)
elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1)
message("\nCompleted in ", elapsed, " minutes")

# ---- 7. Aggregate and report -------------------------------------------------

fold_results <- bind_rows(results_list)
write_csv(fold_results, file.path(output_dir, "spec_grid_cv_results.csv"))

summary_results <- fold_results |>
  group_by(spec_id, state_slopes, covariates, employment) |>
  summarise(
    mean_train_r = mean(train_median_r, na.rm = TRUE),
    mean_test_r  = mean(test_median_r, na.rm = TRUE),
    sd_test_r    = sd(test_median_r, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(mean_test_r))

write_csv(summary_results, file.path(output_dir, "spec_grid_cv_summary.csv"))

message("\n", strrep("=", 90))
message("SPECIFICATION GRID SEARCH: 5-fold spatial CV (2019-2021)")
message(strrep("=", 90))
message(sprintf("\n%-5s  %-14s  %-18s  %-12s  %-10s  %-10s  %-8s",
                "Rank", "State slopes", "Covariates", "Employment", "Train r", "Test r", "SD"))

for (i in seq_len(nrow(summary_results))) {
  r <- summary_results[i, ]
  message(sprintf("%-5d  %-14s  %-18s  %-12s  %-10.4f  %-10.4f  %-8.4f",
                  i, r$state_slopes, r$covariates, r$employment,
                  r$mean_train_r, r$mean_test_r, r$sd_test_r))
}

# Best spec
best <- summary_results[1, ]
message("\nBest spec: state_slopes=", best$state_slopes,
        ", covariates=", best$covariates,
        ", employment=", best$employment)
message("Mean test r: ", round(best$mean_test_r, 4),
        " (train: ", round(best$mean_train_r, 4),
        ", SD: ", round(best$sd_test_r, 4), ")")

# Baseline for comparison
baseline <- summary_results |>
  filter(state_slopes == "none", covariates == "none", employment == "none")
message("Baseline test r: ", round(baseline$mean_test_r, 4))
message("Improvement: ", round(best$mean_test_r - baseline$mean_test_r, 4))

message("\nSaved: ", file.path(output_dir, "spec_grid_cv_results.csv"))
message("Saved: ", file.path(output_dir, "spec_grid_cv_summary.csv"))
