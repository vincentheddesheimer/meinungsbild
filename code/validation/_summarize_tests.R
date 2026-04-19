## Compare within-state r across test-script specs and the status-quo baselines.
## Reads: output/validation/_test_*_summary.csv and _test_*_perparty.csv.
## Prints: a ranked table of all specs by mean_r_within; per-party within-state r
##         for the best spec vs. baseline; and writes a consolidated CSV.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr)
})

out_dir <- file.path(here::here(), "output", "validation")

# ---- Load baseline (status quo) ----
# from within_state_cv_summary.csv (commit 23ccd94): best no-covariate baseline
# achieved mean_r_within ~= 0.602. We take that as the reference.
baseline_summary <- read_csv(file.path(out_dir, "within_state_cv_summary.csv"),
                              show_col_types = FALSE)
baseline_ref <- baseline_summary |>
  transmute(script = "04_strict_cv", spec = covs_name,
            label = paste0(covs_name, ifelse(geo == "bezirk", " + Bezirk RE", "")),
            mean_r_overall = mean_r,
            mean_r_within  = mean_r_within,
            sd_r_within    = sd_r_within)

# ---- Load test-script summaries ----
load_summary <- function(file, script_name) {
  if (!file.exists(file)) return(NULL)
  df <- read_csv(file, show_col_types = FALSE)
  # Harmonize columns: ensure mean_r_overall, mean_r_within are present
  df |> mutate(script = script_name)
}

mun <- load_summary(file.path(out_dir, "_test_mundlak_summary.csv"),     "mundlak")
lag <- load_summary(file.path(out_dir, "_test_lagged_vote_summary.csv"), "lagged_vote")
yr  <- load_summary(file.path(out_dir, "_test_year_window_summary.csv"), "year_window")

# Normalize column names across scripts
mun2 <- if (!is.null(mun)) mun |>
  transmute(script, spec, label,
            mean_r_overall, mean_r_within, sd_r_within) else NULL

# lagged_vote summary only stores median (not mean across folds). Compute from per-party
# if we want SD; for now use median as proxy for mean_r.
lag2 <- if (!is.null(lag)) lag |>
  transmute(script, spec, label,
            mean_r_overall = median_r_overall,
            mean_r_within  = median_r_within,
            sd_r_within    = NA_real_) else NULL

yr2 <- if (!is.null(yr)) yr |>
  transmute(script, spec = paste(window, covariates, sep = "_"),
            label = paste0(window, " / ", covariates),
            mean_r_overall, mean_r_within, sd_r_within) else NULL

all_sum <- bind_rows(baseline_ref, mun2, lag2, yr2) |>
  arrange(desc(mean_r_within))

cat("\n===== All specs ranked by mean_r_within (test-set, strict holdout) =====\n\n")
all_sum |> print(n = Inf)

write_csv(all_sum, file.path(out_dir, "_test_all_specs_ranked.csv"))

# ---- Per-party r for best Mundlak vs. baseline ----
if (file.exists(file.path(out_dir, "_test_mundlak_perparty.csv"))) {
  mpp <- read_csv(file.path(out_dir, "_test_mundlak_perparty.csv"),
                   show_col_types = FALSE)
  pp_avg <- mpp |>
    group_by(spec, label, issue_id) |>
    summarise(r_overall_mean = mean(r, na.rm = TRUE),
              r_within_mean  = mean(r_w, na.rm = TRUE),
              .groups = "drop")

  # Best overall and best within by within-state r
  best_spec_within <- pp_avg |>
    group_by(spec, label) |>
    summarise(median_r_within = median(r_within_mean), .groups = "drop") |>
    arrange(desc(median_r_within)) |>
    slice_head(n = 1) |> pull(spec)

  cat("\n===== Per-party r (Mundlak best within-state spec vs. baseline) =====\n\n")
  cmp <- pp_avg |>
    filter(spec %in% c("baseline", "baseline_bezirk", best_spec_within)) |>
    select(spec, issue_id, r_overall_mean, r_within_mean) |>
    pivot_wider(names_from = spec, values_from = c(r_overall_mean, r_within_mean),
                names_glue = "{spec}_{.value}")
  print(cmp, width = Inf)
}

# ---- Per-party r for lagged-vote specs ----
if (file.exists(file.path(out_dir, "_test_lagged_vote_perparty.csv"))) {
  lpp <- read_csv(file.path(out_dir, "_test_lagged_vote_perparty.csv"),
                   show_col_types = FALSE)
  cat("\n===== Per-party r (lagged 2017 votes, all specs) =====\n\n")
  lpp |>
    group_by(spec, label, issue_id) |>
    summarise(r_overall = mean(r, na.rm = TRUE),
              r_within  = mean(r_w, na.rm = TRUE),
              .groups = "drop") |>
    arrange(spec, issue_id) |>
    print(n = Inf)
}

cat("\nSaved: output/validation/_test_all_specs_ranked.csv\n")
