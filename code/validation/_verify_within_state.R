## Verify the within-state r table in docs/within_state_note.qmd
## Compares production county-level vote-share estimates against BTW 2021 actual.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
})

mb_root    <- here::here()
gerda_root <- normalizePath(file.path(mb_root, "..", "german_election_data"))

est <- readRDS(file.path(mb_root, "data", "estimates", "estimates_kreis.rds"))
fed <- readRDS(file.path(gerda_root, "data", "federal_elections",
                         "county_level", "final", "federal_cty_harm.rds"))

btw21 <- fed |>
  filter(election_year == 2021) |>
  transmute(county_code = as.character(county_code),
            vote_cdu    = cdu_csu,
            vote_spd    = spd,
            vote_gruene = gruene,
            vote_fdp    = fdp,
            vote_afd    = afd,
            vote_linke  = linke_pds) |>
  pivot_longer(-county_code, names_to = "issue_id", values_to = "actual")

vote_est <- est |>
  filter(issue_id %in% c("vote_afd","vote_cdu","vote_spd","vote_gruene","vote_fdp","vote_linke")) |>
  transmute(county_code = as.character(county_code),
            issue_id, estimate = estimate * 100)

m <- inner_join(vote_est, btw21, by = c("county_code","issue_id")) |>
  mutate(state_code = substr(county_code, 1, 2))

cat(sprintf("n obs per party = %d\n\n", nrow(m) / 6))

tab <- m |>
  group_by(issue_id) |>
  summarise(r_overall = cor(estimate, actual),
            r_between = cor(ave(estimate, state_code),
                            ave(actual,   state_code)),
            .groups = "drop")

within <- m |>
  group_by(issue_id, state_code) |>
  mutate(est_dm = estimate - mean(estimate),
         act_dm = actual   - mean(actual)) |>
  ungroup() |>
  group_by(issue_id) |>
  summarise(r_within = cor(est_dm, act_dm),
            .groups = "drop")

share_between <- m |>
  group_by(issue_id) |>
  summarise(var_state = var(ave(actual, state_code)),
            var_total = var(actual),
            .groups = "drop") |>
  mutate(share_between_pct = round(100 * var_state / var_total))

out <- tab |>
  left_join(within, by = "issue_id") |>
  left_join(share_between |> select(issue_id, share_between_pct), by = "issue_id") |>
  arrange(desc(r_overall))

print(out, n = 20)
cat("\nmedians:\n")
cat(sprintf("  r_overall = %.3f\n",  median(out$r_overall)))
cat(sprintf("  r_within  = %.3f\n",  median(out$r_within)))
cat(sprintf("  r_between = %.3f\n",  median(out$r_between)))
