## Test Mundlak device: decompose county covariates into within-state and
## between-state components. Compare against pooled covariates and no covariates.

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
survey_data <- survey |> filter(year %in% 2019:2021, issue_id %in% vote_issues)

# Build covariates with Mundlak decomposition
kreis_cov_std <- kreis_cov |>
  mutate(
    log_pop_density_z = scale(log(cty_pop_density + 1))[, 1],
    unemployment_z    = scale(unemployment_rate)[, 1],
    share_secondary_z = scale(share_secondary_sector)[, 1],
    share_abitur_z    = scale(share_abitur)[, 1],
    bez_code          = substr(county_code, 1, 3),
    state_code        = substr(county_code, 1, 2)
  )

# State means for Mundlak device
state_means <- kreis_cov_std |>
  group_by(state_code) |>
  summarise(
    unemployment_state_mean    = mean(unemployment_z, na.rm = TRUE),
    share_secondary_state_mean = mean(share_secondary_z, na.rm = TRUE),
    share_abitur_state_mean    = mean(share_abitur_z, na.rm = TRUE),
    .groups = "drop"
  )

kreis_cov_std <- kreis_cov_std |>
  left_join(state_means, by = "state_code") |>
  mutate(
    # Within-state (demeaned) components
    unemployment_within    = unemployment_z - unemployment_state_mean,
    share_secondary_within = share_secondary_z - share_secondary_state_mean,
    share_abitur_within    = share_abitur_z - share_abitur_state_mean
  )

all_cov_cols <- c("county_code", "log_pop_density_z",
                  "unemployment_z", "share_secondary_z", "share_abitur_z",
                  "unemployment_state_mean", "share_secondary_state_mean", "share_abitur_state_mean",
                  "unemployment_within", "share_secondary_within", "share_abitur_within",
                  "bez_code")

pred_kreis <- poststrat |>
  group_by(county_code, age_cat, male, educ_label, state_code) |>
  summarise(N = sum(N), .groups = "drop") |>
  left_join(kreis_cov_std |> select(all_of(all_cov_cols)), by = "county_code") |>
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

# Folds
all_counties <- sort(unique(btw21$county_code))
fold_ids <- sample(rep(1:5, length.out = length(all_counties)))
folds <- split(all_counties, fold_ids)

# ---- Specs -------------------------------------------------------------------

# Each spec defines: extra FE terms, whether to use bezirk RE
specs <- list(
  # 1. Baseline: no county covariates
  baseline = list(
    fe_extra = character(0), bezirk = FALSE, label = "No covariates"
  ),
  # 2. Baseline + bezirk
  baseline_bezirk = list(
    fe_extra = character(0), bezirk = TRUE, label = "No covariates + Bezirk RE"
  ),
  # 3. Pooled covariates (current approach)
  pooled = list(
    fe_extra = c("unemployment_z", "share_secondary_z"),
    bezirk = FALSE, label = "Pooled covariates"
  ),
  # 4. Mundlak: within-state + state-mean components
  mundlak = list(
    fe_extra = c("unemployment_within", "share_secondary_within",
                 "unemployment_state_mean", "share_secondary_state_mean"),
    bezirk = FALSE, label = "Mundlak decomposition"
  ),
  # 5. Within-state only (drop between-state component entirely)
  within_only = list(
    fe_extra = c("unemployment_within", "share_secondary_within"),
    bezirk = FALSE, label = "Within-state covariates only"
  ),
  # 6. Mundlak + bezirk
  mundlak_bezirk = list(
    fe_extra = c("unemployment_within", "share_secondary_within",
                 "unemployment_state_mean", "share_secondary_state_mean"),
    bezirk = TRUE, label = "Mundlak + Bezirk RE"
  ),
  # 7. Within-only + bezirk
  within_bezirk = list(
    fe_extra = c("unemployment_within", "share_secondary_within"),
    bezirk = TRUE, label = "Within-state + Bezirk RE"
  ),
  # 8. Pooled + abitur (within-state relevant)
  pooled_abitur = list(
    fe_extra = c("unemployment_z", "share_secondary_z", "share_abitur_z"),
    bezirk = FALSE, label = "Pooled + Abitur share"
  ),
  # 9. Mundlak with abitur
  mundlak_abitur = list(
    fe_extra = c("unemployment_within", "share_secondary_within", "share_abitur_within",
                 "unemployment_state_mean", "share_secondary_state_mean", "share_abitur_state_mean"),
    bezirk = FALSE, label = "Mundlak + Abitur share"
  ),
  # 10. Mundlak + abitur + bezirk
  mundlak_abitur_bezirk = list(
    fe_extra = c("unemployment_within", "share_secondary_within", "share_abitur_within",
                 "unemployment_state_mean", "share_secondary_state_mean", "share_abitur_state_mean"),
    bezirk = TRUE, label = "Mundlak + Abitur + Bezirk RE"
  )
)

# ---- Fit function ------------------------------------------------------------

fit_issue <- function(issue, fe_extra, use_bezirk, holdout = NULL) {
  d <- survey_data |> filter(issue_id == !!issue) |>
    drop_na(y, age_cat, male, educ_label, state_code)

  if (!is.null(holdout)) {
    d <- d |> filter(is.na(county_code) | !(county_code %in% holdout))
  }
  if (nrow(d) < 300) return(NULL)

  cov_cols <- unique(c("county_code", "log_pop_density_z", fe_extra))
  if (use_bezirk) cov_cols <- unique(c(cov_cols, "bez_code"))

  d <- d |>
    left_join(kreis_cov_std |> select(all_of(cov_cols)), by = "county_code") |>
    mutate(
      across(any_of(c("log_pop_density_z", "unemployment_z", "share_secondary_z",
                       "share_abitur_z", "unemployment_within", "share_secondary_within",
                       "share_abitur_within", "unemployment_state_mean",
                       "share_secondary_state_mean", "share_abitur_state_mean")),
             ~ replace_na(.x, 0)),
      age_cat = factor(age_cat, levels = c("18-29", "30-44", "45-59", "60-74", "75+")),
      educ_label = factor(educ_label, levels = c("no_degree", "hauptschule", "realschule",
                                                  "abitur", "university")),
      state_code = factor(state_code), survey_source = factor(survey_source),
      legperiod = factor(legperiod),
      county_code = ifelse(is.na(county_code), "missing", county_code),
      wkr_nr = as.character(wkr_nr)
    )
  if (use_bezirk) d <- d |> mutate(bez_code = ifelse(is.na(bez_code), "missing", bez_code))
  d <- d |> filter(!is.na(wkr_nr)) |> droplevels()
  if (nrow(d) < 300) return(NULL)

  fe <- "y ~ male + log_pop_density_z"
  if (length(fe_extra) > 0) fe <- paste(fe, "+", paste(fe_extra, collapse = " + "))

  re_parts <- c("(1|age_cat)", "(1|educ_label)", "(1|educ_label:age_cat)",
                "(1|male:age_cat)", "(1|male:educ_label)")
  if (nlevels(d$survey_source) > 1) re_parts <- c(re_parts, "(1|survey_source)")
  if (nlevels(d$legperiod) > 1) re_parts <- c(re_parts, "(1|legperiod)")
  re_parts <- c(re_parts, "(1|state_code)")
  if (use_bezirk) re_parts <- c(re_parts, "(1|bez_code)")
  re_parts <- c(re_parts, "(1|county_code)", "(1|wkr_nr)")

  formula_str <- paste(fe, "+", paste(re_parts, collapse = " + "))

  fit <- tryCatch(
    glmmTMB(as.formula(formula_str), data = d, family = binomial),
    error = function(e) {
      re_parts2 <- setdiff(re_parts, "(1|wkr_nr)")
      tryCatch(
        glmmTMB(as.formula(paste(fe, "+", paste(re_parts2, collapse = " + "))),
                data = d, family = binomial),
        error = function(e2) NULL
      )
    }
  )
  if (is.null(fit)) return(NULL)

  pk <- pred_kreis |> mutate(
    survey_source = factor(levels(d$survey_source)[1], levels = levels(d$survey_source)),
    legperiod = factor(tail(sort(unique(as.character(d$legperiod))), 1),
                       levels = levels(d$legperiod)),
    wkr_nr = NA_character_
  )
  if (use_bezirk && !"bez_code" %in% names(pk)) pk$bez_code <- substr(pk$county_code, 1, 3)
  pk$.pred <- predict(fit, newdata = pk, type = "response", allow.new.levels = TRUE)
  pk |> group_by(county_code) |>
    summarise(estimate = weighted.mean(.pred, N, na.rm = TRUE), .groups = "drop")
}

# ---- Run CV ------------------------------------------------------------------

message("\n=== 5-fold strict spatial CV: Mundlak device test ===\n")
message(sprintf("%-35s  %-10s  %-12s  %-10s  %-10s",
                "Spec", "Overall r", "Within-st r", "SD overall", "SD within"))

all_party_rows <- list()
all_summary_rows <- list()

for (spec_name in names(specs)) {
  s <- specs[[spec_name]]
  fold_res <- list()

  for (k in 1:5) {
    test_c <- folds[[k]]
    res <- future_lapply(vote_issues, function(i) {
      fit_issue(i, s$fe_extra, s$bezirk, holdout = test_c)
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
      ungroup() |>
      group_by(issue_id) |>
      summarise(r_w = cor(est_dm, act_dm, use = "complete.obs"), .groups = "drop")

    fold_res[[k]] <- tibble(median_r = median(overall$r),
                            median_r_within = median(within_r$r_w))

    all_party_rows[[length(all_party_rows) + 1]] <- overall |>
      left_join(within_r, by = "issue_id") |>
      mutate(spec = spec_name, label = s$label, fold = k)
  }

  fr <- bind_rows(fold_res)
  message(sprintf("%-35s  %-10.3f  %-12.3f  %-10.3f  %-10.3f",
                  s$label,
                  mean(fr$median_r), mean(fr$median_r_within),
                  sd(fr$median_r), sd(fr$median_r_within)))

  all_summary_rows[[spec_name]] <- tibble(
    spec = spec_name, label = s$label,
    mean_r_overall   = mean(fr$median_r),
    mean_r_within    = mean(fr$median_r_within),
    sd_r_overall     = sd(fr$median_r),
    sd_r_within      = sd(fr$median_r_within),
    n_folds          = nrow(fr)
  )
}

out_dir <- file.path(mb_root, "output", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(bind_rows(all_party_rows),  file.path(out_dir, "_test_mundlak_perparty.csv"))
readr::write_csv(bind_rows(all_summary_rows), file.path(out_dir, "_test_mundlak_summary.csv"))
message("\nSaved: output/validation/_test_mundlak_perparty.csv / _summary.csv")
