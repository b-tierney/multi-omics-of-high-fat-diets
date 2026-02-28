# Consolidated entry-point for core publication analyses.
# Run from repository root so relative paths inside scripts resolve correctly.

core_scripts <- c(
  "scripts/generate_intermediate_files_mousediet.R",
  "scripts/cohort2_diet_time_sex_overall.R",
  "scripts/metabolite_analysis.R",
  "scripts/serum_metabolite_analysis.R",
  "scripts/lipid_analysis.R",
  "scripts/all_metabolomics_reversal_analysis.R",
  "scripts/clean_dirty_analysis.R",
  "scripts/diet_trajectory_clusters_final_c_choice.R"
)

for (f in core_scripts) {
  message("Running: ", f)
  source(f, chdir = TRUE)
}

message("Completed all core analysis scripts.")
