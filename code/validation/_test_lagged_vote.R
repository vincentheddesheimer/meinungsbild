## Test lagged (2017) county-level vote shares as covariates
## Both pooled and within-state demeaned versions

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
survey_data <- survey |> filter(year %in% 2017:2021, issue_id %in% vote_issues)

# Load 2017 county-level results
fed_cty <- readRDS(file.path(gerda_root, "data", "federal_elections",
                              "county_level", "final", "federal_cty_harm.rds"))
btw17 <- fed_cty |>
  filter(election_year == 2017) |>
  transmute(county_code = as.character(county_code),
            afd_17_z   = scale(afd)[, 1],
            cdu_17_z   = scale(cdu_csu)[, 1],
            spd_17_z   = scale(spd)[, 1],
            turnout_17_z = scale(turnout)[, 1])

# Build covariates
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
  ungroup() |>
  left_join(btw17, by = "county_code")

# State-demean 2017 vote shares
state_means_17 <- kreis_cov_std |>
  group_by(state_code) |>
  summarise(across(c(afd_17_z, cdu_17_z, spd_17_z, turnout_17_z), mean, na.rm = TRUE,
                   .names = "{.col}_state_mean"), .groups = "drop")

kreis_cov_std <- kreis_cov_std |>
  left_join(state_means_17, by = "state_code") |>
  mutate(
    afd_17_within    = afd_17_z - afd_17_z_state_mean,
    cdu_17_within    = cdu_17_z - cdu_17_z_state_mean,
    spd_17_within    = spd_17_z - spd_17_z_state_mean,
    turnout_17_within = turnout_17_z - turnout_17_z_state_mean
  )

all_cov_cols <- c("county_code", "log_pop_density_z",
                  "unemployment_within", "share_secondary_within",
                  "afd_17_z", "cdu_17_z", "spd_17_z", "turnout_17_z",
                  "afd_17_within", "cdu_17_within", "spd_17_within", "turnout_17_within")

pred_kreis <- poststrat |>
  group_by(county_code, age_cat, male, educ_label, state_code) |>
  summarise(N = sum(N), .groups = "drop") |>
  left_join(kreis_cov_std |> select(all_of(all_cov_cols)), by = "county_code") |>
  mutate(male = as.integer(male), wkr_nr = NA_character_)

# Ground truth: 2021
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

party_labels <- c(vote_cdu = "CDU/CSU", vote_spd = "SPD", vote_gruene = "Grüne",
                  vote_fdp = "FDP", vote_afd = "AfD", vote_linke = "DIE LINKE")

all_counties <- sort(unique(btw21$county_code))
fold_ids <- sample(rep(1:5, length.out = length(all_counties)))
folds <- split(all_counties, fold_ids)

# ---- Specs -------------------------------------------------------------------

specs <- list(
  baseline = list(
    fe = character(0), label = "No covariates"
  ),
  within_socio = list(
    fe = c("unemployment_within", "share_secondary_within"),
    label = "Within-state socio (best so far)"
  ),
  pooled_2017 = list(
    fe = c("afd_17_z", "cdu_17_z", "spd_17_z", "turnout_17_z"),
    label = "Pooled 2017 vote shares"
  ),
  within_2017 = list(
    fe = c("afd_17_within", "cdu_17_within", "spd_17_within", "turnout_17_within"),
    label = "Within-state 2017 vote shares"
  ),
  within_socio_plus_2017 = list(
    fe = c("unemployment_within", "share_secondary_within",
           "afd_17_within", "cdu_17_within", "spd_17_within", "turnout_17_within"),
    label = "Within socio + within 2017"
  ),
  within_2017_afd_only = list(
    fe = c("afd_17_within", "turnout_17_within"),
    label = "Within 2017 AfD + turnout only"
  )
)

# ---- Fit function ------------------------------------------------------------

fit_issue <- function(issue, fe_extra, holdout) {
  d <- survey_data |> filter(issue_id == !!issue) |>
    drop_na(y, age_cat, male, educ_label, state_code) |>
    filter(is.na(county_code) | !(county_code %in% holdout))
  if (nrow(d) < 300) return(NULL)

  cov_cols <- unique(c("county_code", "log_pop_density_z", fe_extra))
  d <- d |>
    left_join(kreis_cov_std |> select(all_of(cov_cols)), by = "county_code") |>
    mutate(
      across(any_of(c("log_pop_density_z", "unemployment_within", "share_secondary_within",
                       "afd_17_z", "cdu_17_z", "spd_17_z", "turnout_17_z",
                       "afd_17_within", "cdu_17_within", "spd_17_within", "turnout_17_within")),
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
                        data = d, family = binomial), error = function(e2) NULL)
    })
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

# ---- Run CV ------------------------------------------------------------------

message("\n=== 2017 lagged vote shares: strict holdout CV (2017-2021 survey) ===\n")

all_party_rows   <- list()
all_summary_rows <- list()

for (spec_name in names(specs)) {
  s <- specs[[spec_name]]
  all_val <- list()

  for (k in 1:5) {
    test_c <- folds[[k]]
    res <- future_lapply(vote_issues, function(i) fit_issue(i, s$fe, holdout = test_c),
                         future.seed = TRUE)
    names(res) <- vote_issues
    for (iss in names(compact(res))) {
      mrp <- res[[iss]] |> mutate(county_code = as.character(county_code))
      val <- inner_join(mrp, btw21 |> filter(issue_id == iss), by = "county_code") |>
        filter(county_code %in% test_c)
      all_val[[length(all_val) + 1]] <- val |> mutate(fold = k)
    }
  }

  val_all <- bind_rows(all_val)
  val_dm <- val_all |>
    group_by(issue_id, state_code) |>
    mutate(est_dm = estimate - mean(estimate), act_dm = actual - mean(actual)) |>
    ungroup()

  # Per-party, per-fold r
  perfold <- val_all |>
    group_by(issue_id, fold) |>
    summarise(r = cor(estimate, actual, use = "complete.obs"), .groups = "drop")
  perfold_w <- val_dm |>
    group_by(issue_id, fold) |>
    summarise(r_w = cor(est_dm, act_dm, use = "complete.obs"), .groups = "drop")
  party_rows <- left_join(perfold, perfold_w, by = c("issue_id", "fold")) |>
    mutate(spec = spec_name, label = s$label)
  all_party_rows[[spec_name]] <- party_rows

  cat("\n--- ", s$label, " ---\n")
  cat(sprintf("%-14s  %8s  %10s\n", "Party", "r_total", "r_within"))
  for (iss in vote_issues) {
    sub <- val_all |> filter(issue_id == iss)
    sub_dm <- val_dm |> filter(issue_id == iss)
    cat(sprintf("%-14s  %8.3f  %10.3f\n", party_labels[iss],
                cor(sub$estimate, sub$actual, use = "complete.obs"),
                cor(sub_dm$est_dm, sub_dm$act_dm, use = "complete.obs")))
  }
  med_t <- val_all |> group_by(issue_id) |>
    summarise(r = cor(estimate, actual, use = "complete.obs"), .groups = "drop")
  med_w <- val_dm |> group_by(issue_id) |>
    summarise(r = cor(est_dm, act_dm, use = "complete.obs"), .groups = "drop")
  cat(sprintf("%-14s  %8.3f  %10.3f\n", "MEDIAN", median(med_t$r), median(med_w$r)))

  all_summary_rows[[spec_name]] <- tibble(
    spec = spec_name, label = s$label,
    median_r_overall = median(med_t$r),
    median_r_within  = median(med_w$r)
  )
}

out_dir <- file.path(mb_root, "output", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(bind_rows(all_party_rows),   file.path(out_dir, "_test_lagged_vote_perparty.csv"))
readr::write_csv(bind_rows(all_summary_rows), file.path(out_dir, "_test_lagged_vote_summary.csv"))
message("\nSaved: output/validation/_test_lagged_vote_perparty.csv / _summary.csv")
