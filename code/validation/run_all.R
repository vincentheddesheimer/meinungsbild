# ==============================================================================
# run_all.R — Run all validation analyses
# Execute from repo root: Rscript code/validation/run_all.R
# ==============================================================================

mb_root <- here::here()

scripts <- c(
  "01_vote_shares_allyears.R",
  "02_vote_shares_2021.R"
  # Add future validation scripts here:
  # "03_eurobarometer.R",
  # "04_state_polls.R",
)

for (script in scripts) {
  path <- file.path(mb_root, "code", "validation", script)
  message("\n", strrep("=", 70))
  message("Running: ", script)
  message(strrep("=", 70), "\n")
  source(path, local = new.env(parent = globalenv()))
}

message("\n", strrep("=", 70))
message("All validation scripts complete.")
message("Outputs: output/validation/")
message(strrep("=", 70))
