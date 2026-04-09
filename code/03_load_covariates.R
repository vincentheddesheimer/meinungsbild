## 03_load_covariates.R
## Load geographic-level covariates for MRP from GERDA package and election data
##
## Outputs: data/covariates/kreis_covariates.rds
##          data/covariates/bundesland_covariates.rds

library(tidyverse)
library(gerda)

# ---- Paths ----------------------------------------------------------------
mb_root    <- here::here()
gerda_root <- Sys.getenv("GERDA_ROOT",
                         unset = normalizePath(file.path(here::here(), "..", "german_election_data"),
                                               mustWork = FALSE))

dir.create(file.path(mb_root, "data", "covariates"), showWarnings = FALSE, recursive = TRUE)

# ---- 1. Federal election results (county level) ---------------------------

fed_cty <- readRDS(file.path(
  gerda_root, "data", "federal_elections", "county_level", "final",
  "federal_cty_harm.rds"
))

fed_cty_latest <- fed_cty |>
  filter(election_year == max(election_year)) |>
  transmute(
    county_code = county_code,
    fed_year     = election_year,
    fed_turnout  = turnout,
    fed_cdu_csu  = cdu_csu,
    fed_spd      = spd,
    fed_gruene   = gruene,
    fed_fdp      = fdp,
    fed_linke    = linke_pds,
    fed_afd      = afd
  )

message("Election data: BTW ", fed_cty_latest$fed_year[1],
        " (", nrow(fed_cty_latest), " counties)")

# ---- 2. County-level covariates (population, area) -----------------------

cty_covars <- readRDS(file.path(
  gerda_root, "data", "covars_county", "final",
  "cty_area_pop_emp.rds"
))

cty_covars_latest <- cty_covars |>
  filter(year == max(year)) |>
  transmute(
    county_code    = county_code_21,
    county_name    = county_name_21,
    cty_population = population_cty,
    cty_area       = area_cty,
    cty_pop_density = pop_density_cty
  )

# ---- 3. INKAR covariates from gerda package ------------------------------

gerda_cov <- gerda_covariates()

# Use 2021 (full coverage); fall back to most recent available year per county
gerda_latest <- gerda_cov |>
  filter(year <= 2021) |>
  group_by(county_code) |>
  filter(year == max(year)) |>
  ungroup() |>
  select(-year)

message("GERDA covariates: ", nrow(gerda_latest), " counties, ",
        ncol(gerda_latest) - 1, " variables (year <= 2021)")

# Check coverage
na_counts <- gerda_latest |>
  summarise(across(-county_code, ~ sum(is.na(.x)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_na") |>
  filter(n_na > 0)

if (nrow(na_counts) > 0) {
  message("Variables with NAs:")
  for (i in seq_len(nrow(na_counts))) {
    message("  ", na_counts$variable[i], ": ", na_counts$n_na[i], " / ", nrow(gerda_latest))
  }
}

# ---- 4. Derive East/West and state identifiers ----------------------------

derive_state_vars <- function(df) {
  df |>
    mutate(
      state_code = str_sub(county_code, 1, 2),
      east = state_code %in% c("12", "13", "14", "15", "16"),
      berlin = state_code == "11"
    )
}

# ---- 5. Merge all county-level covariates ---------------------------------

kreis_covariates <- cty_covars_latest |>
  left_join(fed_cty_latest, by = "county_code") |>
  left_join(gerda_latest, by = "county_code") |>
  derive_state_vars()

message("\nMerged covariates: ", nrow(kreis_covariates), " counties, ",
        ncol(kreis_covariates), " columns")

# ---- 6. Aggregate to Bundesland level -------------------------------------

bundesland_covariates <- kreis_covariates |>
  group_by(state_code) |>
  summarise(
    n_kreise        = n(),
    bl_population   = sum(cty_population, na.rm = TRUE),
    bl_pop_density  = sum(cty_population, na.rm = TRUE) / sum(cty_area, na.rm = TRUE),
    bl_fed_turnout  = weighted.mean(fed_turnout, cty_population, na.rm = TRUE),
    bl_fed_cdu_csu  = weighted.mean(fed_cdu_csu, cty_population, na.rm = TRUE),
    bl_fed_spd      = weighted.mean(fed_spd, cty_population, na.rm = TRUE),
    bl_fed_gruene   = weighted.mean(fed_gruene, cty_population, na.rm = TRUE),
    bl_fed_fdp      = weighted.mean(fed_fdp, cty_population, na.rm = TRUE),
    bl_fed_linke    = weighted.mean(fed_linke, cty_population, na.rm = TRUE),
    bl_fed_afd      = weighted.mean(fed_afd, cty_population, na.rm = TRUE),
    east            = first(east),
    .groups = "drop"
  )

# ---- 7. Save ---------------------------------------------------------------

saveRDS(kreis_covariates,
        file.path(mb_root, "data", "covariates", "kreis_covariates.rds"))
saveRDS(bundesland_covariates,
        file.path(mb_root, "data", "covariates", "bundesland_covariates.rds"))

message("\nSaved ", nrow(kreis_covariates), " Kreis-level covariates (",
        ncol(kreis_covariates), " cols)")
message("Saved ", nrow(bundesland_covariates), " Bundesland-level covariates")
