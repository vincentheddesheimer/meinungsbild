## Test year window effect on within-state r
## Compare 2019-2021 vs 2018-2021 with best within-state spec (demeaned covs)

set.seed(20260409)
library(tidyverse)
library(glmmTMB)
library(future.apply)

mb_root <- here::here()
gerda_root <- normalizePath(file.path(here::here(), "..", "german_election_data"))
plan(multisession, workers = max(1L, parallelly::availableCores() - 1L))
options(future.globals.maxSize = 1024 * 1024^2)

survey    <- readRDS(file.path(mb_root, "data", "harmonized", "survey_pooled.rds"))
kreis_cov <- readRDS(file.path(mb_root, "data", "covariates", "kreis_covariates.rds"))
poststrat <- readRDS(file.path(mb_root, "data", "poststrat", "poststrat_kreis.rds"))

vote_issues <- c("vote_cdu", "vote_spd", "vote_gruene", "vote_fdp", "vote_afd", "vote_linke")

kreis_cov_std <- kreis_cov |>
  mutate(
    log_pop_density_z = scale(log(cty_pop_density + 1))[, 1],
    unemployment_z    = scale(unemployment_rate)[, 1],
    share_secondary_z = scale(share_secondary_sector)[, 1],
    state_code        = substr(county_code, 1, 2)
  ) |>
  group_by(state_code) |>
  mutate(
    unemployment_within    = unemployment_z - mean(unemployment_z, na.rm = TRUE),
    share_secondary_within = share_secondary_z - mean(share_secondary_z, na.rm = TRUE)
  ) |>
  ungroup()

pred_kreis <- poststrat |>
  group_by(county_code, age_cat, male, educ_label, state_code) |>
  summarise(N = sum(N), .groups = "drop") |>
  left_join(kreis_cov_std |> select(county_code, log_pop_density_z,
                                     unemployment_within, share_secondary_within),
            by = "county_code") |>
  mutate(male = as.integer(male), wkr_nr = NA_character_)

fed_cty <- readRDS(file.path(gerda_root, "data", "federal_elections",
                              "county_level", "final", "federal_cty_harm.rds"))
btw21 <- fed_cty |>
  filter(election_year == 2021) |>
  select(county_code, cdu_csu, spd, gruene, fdp, afd, linke_pds) |>
  pivot_longer(-county_code, names_to = "p", values_to = "actual") |>
  mutate(issue_id = case_when(
    p == "cdu_csu" ~ "vote_cdu", p == "spd" ~ "vote_spd",
    p == "gruene" ~ "vote_gruene", p == "fdp" ~ "vote_fdp",
    p == "afd" ~ "vote_afd", p == "linke_pds" ~ "vote_linke"
  ), county_code = as.character(county_code),
     state_code = substr(county_code, 1, 2)) |>
  select(county_code, state_code, issue_id, actual)

all_counties <- sort(unique(btw21$county_code))
fold_ids <- sample(rep(1:5, length.out = length(all_counties)))
folds <- split(all_counties, fold_ids)

fit_issue <- function(issue, survey_d, use_within_covs, holdout = NULL) {
  d <- survey_d |> filter(issue_id == !!issue) |>
    drop_na(y, age_cat, male, educ_label, state_code)
  if (!is.null(holdout)) d <- d |> filter(is.na(county_code) | !(county_code %in% holdout))
  if (nrow(d) < 300) return(NULL)

  cov_cols <- c("county_code", "log_pop_density_z")
  fe_extra <- character(0)
  if (use_within_covs) {
    cov_cols <- c(cov_cols, "unemployment_within", "share_secondary_within")
    fe_extra <- c("unemployment_within", "share_secondary_within")
  }

  d <- d |>
    left_join(kreis_cov_std |> select(all_of(cov_cols)), by = "county_code") |>
    mutate(
      across(any_of(c("log_pop_density_z", "unemployment_within", "share_secondary_within")),
             ~ replace_na(.x, 0)),
      age_cat = factor(age_cat, levels = c("18-29", "30-44", "45-59", "60-74", "75+")),
      educ_label = factor(educ_label, levels = c("no_degree", "hauptschule", "realschule",
                                                  "abitur", "university")),
      state_code = factor(state_code), survey_source = factor(survey_source),
      legperiod = factor(legperiod),
      county_code = ifelse(is.na(county_code), "missing", county_code),
      wkr_nr = as.character(wkr_nr)
    ) |> filter(!is.na(wkr_nr)) |> droplevels()
  if (nrow(d) < 300) return(NULL)

  fe <- "y ~ male + log_pop_density_z"
  if (length(fe_extra) > 0) fe <- paste(fe, "+", paste(fe_extra, collapse = " + "))

  re_parts <- c("(1|age_cat)", "(1|educ_label)", "(1|educ_label:age_cat)",
                "(1|male:age_cat)", "(1|male:educ_label)")
  if (nlevels(d$survey_source) > 1) re_parts <- c(re_parts, "(1|survey_source)")
  if (nlevels(d$legperiod) > 1) re_parts <- c(re_parts, "(1|legperiod)")
  re_parts <- c(re_parts, "(1|state_code)", "(1|county_code)", "(1|wkr_nr)")

  fit <- tryCatch(
    glmmTMB(as.formula(paste(fe, "+", paste(re_parts, collapse = " + "))),
            data = d, family = binomial),
    error = function(e) {
      re2 <- setdiff(re_parts, "(1|wkr_nr)")
      tryCatch(glmmTMB(as.formula(paste(fe, "+", paste(re2, collapse = " + "))),
                        data = d, family = binomial),
               error = function(e2) NULL)
    }
  )
  if (is.null(fit)) return(NULL)

  pk <- pred_kreis |> mutate(
    survey_source = factor(levels(d$survey_source)[1], levels = levels(d$survey_source)),
    legperiod = factor(tail(sort(unique(as.character(d$legperiod))), 1),
                       levels = levels(d$legperiod)),
    wkr_nr = NA_character_)
  pk$.pred <- predict(fit, newdata = pk, type = "response", allow.new.levels = TRUE)
  pk |> group_by(county_code) |>
    summarise(estimate = weighted.mean(.pred, N, na.rm = TRUE), .groups = "drop")
}

# ---- Run tests ----

year_windows <- list(
  "2019-2021" = 2019:2021,
  "2018-2021" = 2018:2021,
  "2017-2021" = 2017:2021
)

message("\n=== Year window x within-state covariates (strict holdout) ===\n")
message(sprintf("%-12s  %-20s  %-10s  %-12s  %-10s  %-10s",
                "Years", "Covariates", "N resp", "Overall r", "Within r", "SD within"))

all_party_rows   <- list()
all_summary_rows <- list()

for (yr_name in names(year_windows)) {
  yrs <- year_windows[[yr_name]]
  sv <- survey |> filter(year %in% yrs, issue_id %in% vote_issues)
  n_resp <- n_distinct(sv$respondent_id)

  for (use_within in c(FALSE, TRUE)) {
    fold_res <- list()
    for (k in 1:5) {
      test_c <- folds[[k]]
      res <- future_lapply(vote_issues, function(i) {
        fit_issue(i, sv, use_within, holdout = test_c)
      }, future.seed = TRUE)
      names(res) <- vote_issues
      ok <- names(compact(res))
      if (length(ok) < 4) next

      mrp <- map_dfr(ok, ~ mutate(res[[.x]], issue_id = .x)) |>
        mutate(county_code = as.character(county_code))
      val <- inner_join(mrp, btw21, by = c("county_code", "issue_id")) |>
        filter(county_code %in% test_c)

      overall <- val |> group_by(issue_id) |>
        summarise(r = cor(estimate, actual, use = "complete.obs"), .groups = "drop")
      within_r <- val |>
        group_by(issue_id, state_code) |>
        mutate(est_dm = estimate - mean(estimate), act_dm = actual - mean(actual)) |>
        ungroup() |> group_by(issue_id) |>
        summarise(r_w = cor(est_dm, act_dm, use = "complete.obs"), .groups = "drop")

      fold_res[[k]] <- tibble(median_r = median(overall$r),
                              median_r_within = median(within_r$r_w))

      all_party_rows[[length(all_party_rows) + 1]] <- overall |>
        left_join(within_r, by = "issue_id") |>
        mutate(window = yr_name, covariates = ifelse(use_within, "within-state FE", "none"),
               fold = k, n_resp = n_resp)
    }
    fr <- bind_rows(fold_res)
    cov_label <- ifelse(use_within, "within-state FE", "none")
    message(sprintf("%-12s  %-20s  %-10d  %-12.3f  %-10.3f  %-10.3f",
                    yr_name, cov_label, n_resp,
                    mean(fr$median_r), mean(fr$median_r_within), sd(fr$median_r_within)))
    all_summary_rows[[length(all_summary_rows) + 1]] <- tibble(
      window = yr_name, covariates = cov_label, n_resp = n_resp,
      mean_r_overall = mean(fr$median_r),
      mean_r_within  = mean(fr$median_r_within),
      sd_r_within    = sd(fr$median_r_within)
    )
  }
}

out_dir <- file.path(mb_root, "output", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(bind_rows(all_party_rows),   file.path(out_dir, "_test_year_window_perparty.csv"))
readr::write_csv(bind_rows(all_summary_rows), file.path(out_dir, "_test_year_window_summary.csv"))
message("\nSaved: output/validation/_test_year_window_perparty.csv / _summary.csv")
