# Scripts Directory Guide

This folder contains the active analysis workflows used for publication outputs in this repository.

## Primary Scripts

- `all_metabolomics_reversal_analysis.R`: End-to-end reversal-focused metabolomics analysis workflow, including differential/regression summaries used downstream for app tables and figures.
- `clean_dirty_analysis.R`: Clean-versus-dirty taxon/pathway behavior analysis and comparative association summaries.
- `cohort2_diet_time_sex_overall.R`: Main cohort 2 taxonomic association workflow across diet, time, and sex; source logic for key diet-associated taxa filters.
- `diet_trajectory_clusters_final_c_choice.R`: Diet trajectory clustering workflow and cluster selection/finalization script.
- `generate_intermediate_files_mousediet.R`: Builds merged/intermediate RDS objects used across analysis and app workflows.
- `lipid_analysis.R`: Lipidomics-specific association analysis and output generation.
- `metabolite_analysis.R`: Core metabolomics association analysis workflow for fecal/serum/lipid features.
- `serum_metabolite_analysis.R`: Serum-focused metabolite association and summary analysis pipeline.
- `run_all_analyses.R`: Lightweight entry-point runner to execute core analysis scripts in sequence.

## Consolidated Execution

To run the main publication analyses in one pass:

```r
source("scripts/run_all_analyses.R")
```

## Notes

- Execution order in `run_all_analyses.R` is intentionally explicit and can be edited for partial workflows.
- Some scripts expect project-relative paths and precomputed intermediate objects.
