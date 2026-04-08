## 04b_fit_all_lme4.R  — PRIMARY ESTIMATION SCRIPT
## Fit deep MRP models for all issues using lme4, then poststratify
## to Kreis (county), Bundesland (state), and Wahlkreis (electoral district).
##
## Model specification (per issue):
##   y ~ male + employed +
##       fed_afd_share_z + fed_cdu_share_z + fed_turnout_z + log_pop_density_z +
##       (1 | age_cat) + (1 | educ_label) +
##       (1 | educ_label:age_cat) + (1 | male:age_cat) + (1 | male:educ_label) +
##       (1 | employed:age_cat) +
##       (1 | survey_source) + (1 | legperiod) +
##       (1 | state_code) + (1 | county_code) + (1 | wkr_nr)
##
## References:
##   Ghitza & Gelman (2013) — deep demographic interactions
##   Goplerud (2024) — deep hierarchical models match ML ensembles
##
## Outputs: estimates_kreis.rds, estimates_bundesland.rds, estimates_wkr.rds

library(tidyverse)
library(lme4)
library(future.apply)

mb_root <- here::here()
dir.create(file.path(mb_root, "data", "estimates", "models"), showWarnings = FALSE, recursive = TRUE)

plan(multisession, workers = max(1L, parallelly::availableCores() - 1L))
options(future.globals.maxSize = 1024 * 1024^2)  # 1 GiB

# ---- 1. Load data -----------------------------------------------------------

survey    <- readRDS(file.path(mb_root, "data", "harmonized", "survey_pooled.rds"))
kreis_cov <- readRDS(file.path(mb_root, "data", "covariates", "kreis_covariates.rds"))
poststrat <- readRDS(file.path(mb_root, "data", "poststrat", "poststrat_kreis.rds"))
conc      <- read_csv(file.path(mb_root, "data", "issue_concordance.csv"),
                      show_col_types = FALSE)

# WKR crosswalk for WKR-level poststratification
wkr_cw <- read_csv(
  file.path(mb_root, "data", "raw", "gles", "crosswalk", "btw_21_to_krs21.csv"),
  show_col_types = FALSE
) |>
  transmute(
    wkr_nr = as.integer(WKR_NR),
    county_code = str_pad(as.character(AGS), 5, pad = "0"),
    pct_area = PERCENTAGE / 100
  )

# ---- 2. Prepare covariates --------------------------------------------------

kreis_cov_std <- kreis_cov |>
  mutate(
    fed_afd_share_z   = scale(fed_afd)[, 1],
    fed_cdu_share_z   = scale(fed_cdu_csu)[, 1],
    fed_turnout_z     = scale(fed_turnout)[, 1],
    log_pop_density_z = scale(log(cty_pop_density + 1))[, 1]
  )

# ---- 3. Build poststrat prediction data --------------------------------------

# Kreis-level
pred_kreis <- poststrat |>
  left_join(kreis_cov_std |> select(county_code, fed_afd_share_z, fed_cdu_share_z,
                                     fed_turnout_z, log_pop_density_z),
            by = "county_code") |>
  mutate(male = as.integer(male), employed = as.integer(employed), wkr_nr = NA_character_)

# WKR-level poststrat (non-Berlin: county crosswalk; Berlin: Bezirk-level demographics)
poststrat_wkr_base <- poststrat |>
  inner_join(wkr_cw, by = "county_code", relationship = "many-to-many") |>
  mutate(N_wkr = N * pct_area) |>
  group_by(wkr_nr, age_cat, male, educ_label, employed) |>
  summarise(N = sum(N_wkr, na.rm = TRUE), .groups = "drop")

# Replace Berlin WKRs (75-86) with Bezirk-level demographics from Zensus + EWR
berlin_wkr_file <- file.path(mb_root, "data", "poststrat", "poststrat_wkr_berlin.rds")
if (file.exists(berlin_wkr_file)) {
  poststrat_wkr_berlin <- readRDS(berlin_wkr_file) |>
    mutate(age_cat = factor(age_cat, levels = levels(poststrat_wkr_base$age_cat)),
           educ_label = factor(educ_label, levels = levels(poststrat_wkr_base$educ_label)))
  berlin_wkr_nrs <- unique(poststrat_wkr_berlin$wkr_nr)

  # Berlin Bezirk file may not yet have employed column — if missing, expand it
  if (!"employed" %in% names(poststrat_wkr_berlin)) {
    # Use Berlin county (11000) employment shares from main poststrat to split
    berlin_emp <- poststrat |>
      filter(state_code == "11") |>
      group_by(age_cat, male, educ_label, employed) |>
      summarise(N_emp = sum(N), .groups = "drop") |>
      group_by(age_cat, male, educ_label) |>
      mutate(emp_share = N_emp / sum(N_emp)) |>
      ungroup() |>
      select(age_cat, male, educ_label, employed, emp_share)

    poststrat_wkr_berlin <- poststrat_wkr_berlin |>
      left_join(berlin_emp, by = c("age_cat", "male", "educ_label"),
                relationship = "many-to-many") |>
      mutate(N = round(N * emp_share)) |>
      select(-emp_share)
  }

  poststrat_wkr <- bind_rows(
    poststrat_wkr_base |> filter(!(wkr_nr %in% berlin_wkr_nrs)),
    poststrat_wkr_berlin |> select(wkr_nr, age_cat, male, educ_label, employed, N)
  )
  message("  Berlin WKRs (", paste(range(berlin_wkr_nrs), collapse="-"),
          "): using Bezirk-level demographics from Zensus 2022")
} else {
  poststrat_wkr <- poststrat_wkr_base
  message("  NOTE: poststrat_wkr_berlin.rds not found — Berlin WKRs use county-level demographics")
}

# WKR-level covariates (pop-weighted average of constituent Kreise)
wkr_covs <- wkr_cw |>
  left_join(kreis_cov_std |> select(county_code, fed_afd_share_z, fed_cdu_share_z,
                                     fed_turnout_z, log_pop_density_z),
            by = "county_code") |>
  group_by(wkr_nr) |>
  summarise(across(c(fed_afd_share_z, fed_cdu_share_z, fed_turnout_z,
                     log_pop_density_z),
                   ~ weighted.mean(.x, pct_area, na.rm = TRUE)),
            .groups = "drop")

pred_wkr <- poststrat_wkr |>
  left_join(wkr_covs, by = "wkr_nr") |>
  mutate(
    male = as.integer(male),
    employed = as.integer(employed),
    county_code = NA_character_,
    state_code = NA_character_
  )

# ---- 4. Fit + poststratify function -----------------------------------------

fit_and_poststratify <- function(issue, survey_data, pred_kreis_df, pred_wkr_df, conc_df) {
  d <- survey_data |>
    filter(issue_id == !!issue) |>
    drop_na(y, age_cat, male, educ_label, employed, state_code)

  # Skip continuous outcomes (lr_self)
  if (any(conc_df$binary_rule[conc_df$issue_id == issue] == "continuous", na.rm = TRUE)) {
    message("  Skipping continuous issue: ", issue)
    return(NULL)
  }

  if (nrow(d) < 300) {
    message("  Skipping ", issue, ": only ", nrow(d), " obs")
    return(NULL)
  }

  # Add covariates
  d <- d |>
    left_join(kreis_cov_std |> select(county_code, fed_afd_share_z, fed_cdu_share_z,
                                       fed_turnout_z, log_pop_density_z),
              by = "county_code") |>
    mutate(
      across(c(fed_afd_share_z, fed_cdu_share_z, fed_turnout_z, log_pop_density_z),
             ~ replace_na(.x, 0)),
      employed   = as.integer(employed),
      age_cat    = factor(age_cat, levels = c("18-29","30-44","45-59","60-74","75+")),
      educ_label = factor(educ_label, levels = c("no_degree","hauptschule","realschule",
                                                  "abitur","university")),
      state_code = factor(state_code),
      survey_source = factor(survey_source),
      legperiod  = factor(legperiod),
      county_code = ifelse(is.na(county_code), "missing", county_code),
      wkr_nr = as.character(wkr_nr)
    ) |>
    filter(!is.na(wkr_nr)) |>
    droplevels()

  if (nrow(d) < 300) {
    message("  Skipping ", issue, " after filtering: only ", nrow(d), " obs")
    return(NULL)
  }

  # Determine if we have enough variation for each random effect
  n_kreise  <- n_distinct(d$county_code[d$county_code != "missing"])
  n_wkr     <- n_distinct(d$wkr_nr)
  n_sources <- nlevels(d$survey_source)
  n_legper  <- nlevels(d$legperiod)
  n_states  <- nlevels(d$state_code)

  # Build formula adaptively — only add RE when > 1 level observed
  fe <- "y ~ male + employed + fed_afd_share_z + fed_cdu_share_z + fed_turnout_z + log_pop_density_z"
  re <- "(1 | age_cat) + (1 | educ_label) + (1 | educ_label:age_cat) + (1 | male:age_cat) + (1 | male:educ_label) + (1 | employed:age_cat)"

  if (n_sources > 1) re <- paste0(re, " + (1 | survey_source)")
  if (n_legper > 1)  re <- paste0(re, " + (1 | legperiod)")

  if (n_states > 1)  re <- paste0(re, " + (1 | state_code)")
  if (n_kreise > 10) re <- paste0(re, " + (1 | county_code)")
  if (n_wkr > 10) re <- paste0(re, " + (1 | wkr_nr)")

  formula_str <- paste(fe, "+", re)

  message("  N=", nrow(d), " | Kreise=", n_kreise, " | WKR=", n_wkr)

  # Fit (with fallback: drop wkr_nr RE if first attempt fails)
  fit <- tryCatch(
    glmer(
      as.formula(formula_str),
      data = d,
      family = binomial(link = "logit"),
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5), calc.derivs = FALSE),
      nAGQ = 0
    ),
    error = function(e) {
      message("  First attempt failed: ", e$message)
      # Retry without wkr_nr random effect
      formula_fallback <- gsub(" \\+ \\(1 \\| wkr_nr\\)", "", formula_str)
      message("  Retrying without (1 | wkr_nr)...")
      tryCatch(
        glmer(
          as.formula(formula_fallback),
          data = d,
          family = binomial(link = "logit"),
          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5), calc.derivs = FALSE),
          nAGQ = 0
        ),
        error = function(e2) {
          message("  ERROR (retry): ", e2$message)
          NULL
        }
      )
    }
  )

  if (is.null(fit)) return(NULL)

  # --- Poststratify to Kreis ---
  pk <- pred_kreis_df |>
    mutate(
      survey_source = factor(levels(d$survey_source)[1], levels = levels(d$survey_source)),
      legperiod = factor(tail(sort(unique(as.character(d$legperiod))), 1),
                         levels = levels(d$legperiod)),
      county_code = county_code,
      wkr_nr = NA_character_
    )

  pk$.pred <- predict(fit, newdata = pk, type = "response", allow.new.levels = TRUE)

  est_kreis <- pk |>
    group_by(county_code) |>
    summarise(estimate = weighted.mean(.pred, N, na.rm = TRUE),
              pop = sum(N), .groups = "drop")

  est_land <- pk |>
    group_by(state_code) |>
    summarise(estimate = weighted.mean(.pred, N, na.rm = TRUE),
              pop = sum(N), .groups = "drop")

  # --- Poststratify to WKR ---
  pw <- pred_wkr_df |>
    mutate(
      survey_source = factor(levels(d$survey_source)[1], levels = levels(d$survey_source)),
      legperiod = factor(tail(sort(unique(as.character(d$legperiod))), 1),
                         levels = levels(d$legperiod)),
      wkr_nr = as.character(wkr_nr)
    )

  pw$.pred <- predict(fit, newdata = pw, type = "response", allow.new.levels = TRUE)

  est_wkr <- pw |>
    group_by(wkr_nr) |>
    summarise(estimate = weighted.mean(.pred, N, na.rm = TRUE),
              pop = sum(N), .groups = "drop")

  list(kreis = est_kreis, bundesland = est_land, wahlkreis = est_wkr, fit = fit)
}

# ---- 5. Run all issues -------------------------------------------------------

issue_ids <- survey |>
  distinct(issue_id) |>
  pull()

# Remove lr_self (continuous)
issue_ids <- setdiff(issue_ids, "lr_self")

message("\n=== Fitting ", length(issue_ids), " issues with lme4 (",
        parallelly::availableCores() - 1L, " workers) ===\n")

results <- future_lapply(issue_ids, function(issue) {
  fit_and_poststratify(issue, survey, pred_kreis, pred_wkr, conc)
}, future.seed = TRUE)
names(results) <- issue_ids

# ---- 6. Stack and save estimates ---------------------------------------------

successful <- names(compact(results))
message("\n", length(successful), " / ", length(issue_ids), " issues fitted successfully")

estimates_kreis <- map_dfr(successful, ~ mutate(results[[.x]]$kreis, issue_id = .x))
estimates_land  <- map_dfr(successful, ~ mutate(results[[.x]]$bundesland, issue_id = .x))
estimates_wkr   <- map_dfr(successful, ~ mutate(results[[.x]]$wahlkreis, issue_id = .x))

saveRDS(estimates_kreis, file.path(mb_root, "data", "estimates", "estimates_kreis.rds"))
saveRDS(estimates_land,  file.path(mb_root, "data", "estimates", "estimates_bundesland.rds"))
saveRDS(estimates_wkr,   file.path(mb_root, "data", "estimates", "estimates_wkr.rds"))

message("\nSaved estimates:")
message("  Kreis:      ", nrow(estimates_kreis), " rows (", n_distinct(estimates_kreis$issue_id), " issues)")
message("  Bundesland: ", nrow(estimates_land), " rows")
message("  Wahlkreis:  ", nrow(estimates_wkr), " rows")

# ---- 7. Quick validation summary --------------------------------------------

bl_names <- c("01"="SH","02"="HH","03"="NI","04"="HB","05"="NW","06"="HE",
              "07"="RP","08"="BW","09"="BY","10"="SL","11"="BE","12"="BB",
              "13"="MV","14"="SN","15"="ST","16"="TH")

direct <- survey |>
  filter(issue_id != "lr_self") |>
  group_by(issue_id, state_code) |>
  summarise(direct = mean(y, na.rm = TRUE), n = n(), .groups = "drop")

val <- estimates_land |>
  left_join(direct, by = c("issue_id", "state_code")) |>
  filter(!is.na(direct))

val_summary <- val |>
  group_by(issue_id) |>
  summarise(
    r = cor(estimate, direct),
    rmse = sqrt(mean((estimate - direct)^2)),
    n_states = n(),
    .groups = "drop"
  ) |>
  arrange(desc(r))

message("\n=== Bundesland-level validation ===")
print(as.data.frame(val_summary), row.names = FALSE)
message("\nMedian correlation: ", round(median(val_summary$r, na.rm = TRUE), 3))
message("Median RMSE: ", round(100 * median(val_summary$rmse, na.rm = TRUE), 1), " pp")
