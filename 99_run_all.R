# 99_run_all.R
# Master runner. Execute from the root folder with: source("R/99_run_all.R")
# Long-running null simulations are intentionally NOT triggered automatically.
steps <- c("00_setup_and_validate.R","01_chapter1_descriptive_tests.R","02_chapter1_figures.R",
           "03_community_overlap_beta_diversity.R","04_distance_decay_and_multivariate.R",
           "05_network_site_similarity.R","06_network_bipartite_taxon_connectors.R",
           "08_export_gephi_files.R")
for (s in steps) {
  message("\n===== Running ",s," =====")
  source(file.path("R",s),local=new.env(parent=globalenv()))
}
writeLines(capture.output(sessionInfo()),"Outputs/Diagnostics/sessionInfo.txt")
message("Core workflow complete. Run R/07_network_nulls_modules_sensitivity.R separately for null-model analyses.")
