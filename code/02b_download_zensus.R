## 02b_download_zensus.R
## Parse Zensus 2022 tables and build the poststratification frame.
##
## Input:  data/raw/zensus/2000S-4020_flat/2000S-4020_de_flat.csv
##         (Age × Gender × Erwerbsstatus × Schulabschluss at Kreis level,
##          flat format from ergebnisse.zensus2022.de)
## Output: data/poststrat/poststrat_kreis.rds      (100 cells/county: 5 age × 2 gender × 5 educ × 2 employed)
##         data/poststrat/poststrat_bundesland.rds

library(tidyverse)

set.seed(20260408)

mb_root <- here::here()

zensus_path <- file.path(mb_root, "data", "raw", "zensus",
                         "2000S-4020_flat", "2000S-4020_de_flat.csv")
stopifnot(file.exists(zensus_path))

# ---- 1. Read flat-format CSV ------------------------------------------------

d <- read_delim(zensus_path, delim = ";", show_col_types = FALSE)
message("Read ", nrow(d), " rows from Zensus 4020 flat file")

# ---- 2. Mappings ------------------------------------------------------------

# Age: 5-year groups → MRP age categories (15-19 included as proxy for 18-19)
age_mapping <- c(
  "ALT015B019" = "18-29", "ALT020B024" = "18-29", "ALT025B029" = "18-29",
  "ALT030B034" = "30-44", "ALT035B039" = "30-44", "ALT040B044" = "30-44",
  "ALT045B049" = "45-59", "ALT050B054" = "45-59", "ALT055B059" = "45-59",
  "ALT060B064" = "60-74", "ALT065B069" = "60-74", "ALT070B074" = "60-74",
  "ALT075B079" = "75+",   "ALT080B084" = "75+",   "ALT085B089" = "75+",
  "ALT090BXXX" = "75+"
)

# Education: Schulabschluss → MRP categories (POS merged into Realschule)
educ_mapping <- c(
  "ABSCH-HAUPT"   = "hauptschule",
  "ABSCH-POS"     = "realschule",
  "ABSCH-REAL-01" = "realschule",
  "ABSCH-FHRABI"  = "abitur",
  "ABSCH-SCHUL-X" = "no_degree"
)

gender_mapping <- c("GESM" = 1L, "GESW" = 0L)

# Employment: use ERWERBSTAETIG for employed count, Insgesamt (NA code) for total
# Not employed = total - employed

# ---- 3. Filter and aggregate ------------------------------------------------

# Keep only: adult age groups, specific gender, specific education,
# and either ERWERBSTAETIG or Insgesamt (total) for employment
cells <- d |>
  filter(
    `2_variable_attribute_code` %in% names(age_mapping),
    `3_variable_attribute_code` %in% names(gender_mapping),
    `5_variable_attribute_code` %in% names(educ_mapping),
    `4_variable_attribute_code` == "ERWERBSTAETIG" | is.na(`4_variable_attribute_code`)
  ) |>
  mutate(
    county_code = `1_variable_attribute_code`,
    age_cat     = age_mapping[`2_variable_attribute_code`],
    male        = gender_mapping[`3_variable_attribute_code`],
    educ_label  = educ_mapping[`5_variable_attribute_code`],
    emp_status  = if_else(is.na(`4_variable_attribute_code`), "total", "employed"),
    N           = suppressWarnings(as.numeric(value))
  ) |>
  mutate(N = replace_na(N, 0))

message("Filtered cells: ", nrow(cells))

# Aggregate 5-year age groups to MRP bins, POS+Realschule merged
cells_agg <- cells |>
  group_by(county_code, age_cat, male, educ_label, emp_status) |>
  summarise(N = sum(N), .groups = "drop") |>
  pivot_wider(names_from = emp_status, values_from = N, values_fill = 0) |>
  mutate(not_employed = pmax(total - employed, 0))  # pmax guards against rounding

message("Aggregated cells: ", nrow(cells_agg),
        " (expect ", n_distinct(cells_agg$county_code), " × 40 = ",
        n_distinct(cells_agg$county_code) * 40, ")")

# ---- 4. Pivot to long: one row per (county, age, gender, educ, employed) ----

ps_employed <- cells_agg |>
  select(county_code, age_cat, male, educ_label, N = employed) |>
  mutate(employed = 1L)

ps_not_employed <- cells_agg |>
  select(county_code, age_cat, male, educ_label, N = not_employed) |>
  mutate(employed = 0L)

poststrat_raw <- bind_rows(ps_employed, ps_not_employed)

# ---- 5. Split Abitur into abitur + university --------------------------------

# Among adults with Abitur, ~55% have a university degree (Zensus 2022 national estimate)
uni_share_of_abitur <- 0.55

ps_abitur <- poststrat_raw |> filter(educ_label == "abitur")
ps_other  <- poststrat_raw |> filter(educ_label != "abitur")

ps_abitur_only <- ps_abitur |>
  mutate(N = round(N * (1 - uni_share_of_abitur)))

ps_university <- ps_abitur |>
  mutate(educ_label = "university",
         N = round(N * uni_share_of_abitur))

poststrat_kreis <- bind_rows(ps_other, ps_abitur_only, ps_university)

# ---- 6. Final formatting ----------------------------------------------------

poststrat_kreis <- poststrat_kreis |>
  mutate(
    age_cat    = factor(age_cat, levels = c("18-29", "30-44", "45-59", "60-74", "75+")),
    educ_label = factor(educ_label, levels = c("no_degree", "hauptschule", "realschule",
                                                "abitur", "university")),
    state_code = str_sub(county_code, 1, 2)
  ) |>
  arrange(county_code, age_cat, male, educ_label, employed)

# ---- 7. Diagnostics ---------------------------------------------------------

n_counties <- n_distinct(poststrat_kreis$county_code)
n_cells <- nrow(poststrat_kreis)

message("\n=== Zensus 2022 Poststratification Frame ===")
message("Cells: ", n_cells, " (", n_counties, " Kreise × ", n_cells / n_counties, " cells)")
message("Total adult population: ", format(sum(poststrat_kreis$N), big.mark = ","))

message("\nAge distribution:")
poststrat_kreis |>
  group_by(age_cat) |>
  summarise(N = sum(N)) |>
  mutate(pct = round(100 * N / sum(N), 1)) |>
  print()

message("\nGender distribution:")
poststrat_kreis |>
  group_by(male) |>
  summarise(N = sum(N)) |>
  mutate(pct = round(100 * N / sum(N), 1)) |>
  print()

message("\nEducation distribution:")
poststrat_kreis |>
  group_by(educ_label) |>
  summarise(N = sum(N)) |>
  mutate(pct = round(100 * N / sum(N), 1)) |>
  print()

message("\nEmployment distribution:")
poststrat_kreis |>
  group_by(employed) |>
  summarise(N = sum(N)) |>
  mutate(pct = round(100 * N / sum(N), 1)) |>
  print()

message("\nPer-Kreis stats:")
kreis_stats <- poststrat_kreis |>
  group_by(county_code) |>
  summarise(n_cells = n(), pop = sum(N), min_cell = min(N), zero_cells = sum(N == 0))
message("  Pop range: ", format(min(kreis_stats$pop), big.mark = ","),
        " - ", format(max(kreis_stats$pop), big.mark = ","))
message("  Smallest cell: ", min(kreis_stats$min_cell))
message("  Counties with zero cells: ", sum(kreis_stats$zero_cells > 0),
        " (total zero cells: ", sum(kreis_stats$zero_cells), ")")

# ---- 8. Marginal check against old 50-cell frame ----------------------------

old_ps_path <- file.path(mb_root, "data", "poststrat", "poststrat_kreis.rds")
if (file.exists(old_ps_path)) {
  old_ps <- readRDS(old_ps_path)
  # Sum new frame over employment → should match old frame
  new_marginal <- poststrat_kreis |>
    group_by(county_code, age_cat, male, educ_label) |>
    summarise(N_new = sum(N), .groups = "drop") |>
    arrange(county_code, age_cat, male, educ_label)
  old_check <- old_ps |>
    arrange(county_code, age_cat, male, educ_label) |>
    select(county_code, age_cat, male, educ_label, N_old = N)

  comparison <- inner_join(new_marginal, old_check,
                           by = c("county_code", "age_cat", "male", "educ_label"))
  message("\n=== Marginal check vs old 50-cell frame ===")
  message("Matched cells: ", nrow(comparison), " / ", nrow(old_check))
  message("Total pop (new): ", format(sum(comparison$N_new), big.mark = ","))
  message("Total pop (old): ", format(sum(comparison$N_old), big.mark = ","))

  # They may differ slightly because the data source is different (table 4020 vs 3044)
  # but should be very close
  comparison <- comparison |> mutate(diff = N_new - N_old, pct_diff = 100 * diff / pmax(N_old, 1))
  message("Median cell difference: ", median(comparison$diff))
  message("Max |cell difference|: ", max(abs(comparison$diff)))
  message("Correlation: ", round(cor(comparison$N_new, comparison$N_old), 5))
}

# ---- 9. Aggregate to Bundesland level ----------------------------------------

poststrat_bundesland <- poststrat_kreis |>
  group_by(state_code, age_cat, male, educ_label, employed) |>
  summarise(N = sum(N), .groups = "drop")

# ---- 10. Save ----------------------------------------------------------------

dir.create(file.path(mb_root, "data", "poststrat"), showWarnings = FALSE, recursive = TRUE)

saveRDS(poststrat_kreis,
        file.path(mb_root, "data", "poststrat", "poststrat_kreis.rds"))
saveRDS(poststrat_bundesland,
        file.path(mb_root, "data", "poststrat", "poststrat_bundesland.rds"))

message("\nSaved poststrat frames (with employment) to data/poststrat/")
message("100 cells per county: 5 age × 2 gender × 5 educ × 2 employed")
