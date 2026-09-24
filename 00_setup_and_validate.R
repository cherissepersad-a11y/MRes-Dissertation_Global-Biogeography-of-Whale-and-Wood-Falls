# 00_setup_and_validate.R
# MRes organic-fall dissertation: reproducible setup and read-only validation
#
# Purpose
# -------
# This script establishes the computational environment, imports the finalised
# analysis-ready database and performs checks that should pass before any
# ecological analysis is run. The source CSV files are NEVER edited in place.
#
# Run this script first from the root of the reproducibility package.

required_packages <- c(
  "tidyverse", "vegan", "betapart", "geosphere", "ggplot2", "ggrepel",
  "patchwork", "igraph", "tidygraph", "ggraph", "scales", "broom"
)

missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  message("Install the following packages before continuing: ",
          paste(missing_packages, collapse = ", "))
  message("Run: install.packages(c(",
          paste(sprintf('"%s"', missing_packages), collapse = ", "), "))")
  stop("Required R packages are missing.")
}

invisible(lapply(required_packages, library, character.only = TRUE))
set.seed(20260917)

DATA_DIR <- "Database"
OUT_DIR  <- "Outputs"
dir.create(file.path(OUT_DIR, "Figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "Diagnostics"), recursive = TRUE, showWarnings = FALSE)

prefix <- "MRes_Dissertation_Database_FINAL_AnalysisReady_"
read_final_csv <- function(stem) {
  readr::read_csv(file.path(DATA_DIR, paste0(prefix, stem, ".csv")), show_col_types = FALSE)
}

# -------------------------------------------------------------------------
# Analytical site universe
#
# SITE0223 (Panama Basin; 5.344 N, 81.9365 W) retains the inherited
# Verbatim_Basin value for provenance, but its analytical geographic
# assignment is corrected to the Pacific. This correction affects derived
# geographic summaries only; the underlying biological occurrence evidence
# is unchanged.
# -------------------------------------------------------------------------

site_universe <- read_final_csv("Analysis_SiteUniverse") %>%
  mutate(
    Longitude = readr::parse_number(
      na_if(Longitude, "unknown")
    ),
    
    Latitude = readr::parse_number(
      na_if(Latitude, "unknown")
    ),
    
    # Correct previously verified Panama Basin geographic assignment.
    Assigned_Basin = if_else(
      Site_ID == "SITE0223",
      "Pacific",
      Assigned_Basin
    ),
    
    Ocean = if_else(
      Site_ID == "SITE0223",
      "North Pacific",
      Ocean
    ),
    
    Region = if_else(
      Site_ID == "SITE0223",
      "NE Pacific",
      Region
    ),
    
    Node_ID = paste(
      Site_ID,
      Fall_Type,
      sep = "__"
    ),
    
    Basin = coalesce(
      Assigned_Basin,
      Ocean,
      "Unresolved"
    )
  )

# Explicit validation of the previously verified SITE0223 correction.
stopifnot(
  site_universe %>%
    filter(Site_ID == "SITE0223") %>%
    summarise(
      ok = all(
        Assigned_Basin == "Pacific",
        Basin == "Pacific",
        Ocean == "North Pacific",
        Region == "NE Pacific"
      )
    ) %>%
    pull(ok)
)
events <- readr::read_csv(
  file.path(DATA_DIR, paste0(prefix, "Canonical_Events.csv")),
  col_types = cols(
    .default = col_guess(),
    Wood_Length_cm = col_character(),
    Wood_Width_cm = col_character(),
    Wood_Height_cm = col_character(),
    Wood_Mass = col_character(),
    Animal_Group = col_character()
  ),
  show_col_types = FALSE
)
sites <- read_final_csv("Canonical_Sites")
occurrences <- read_final_csv("Canonical_Occurrences")
visits <- read_final_csv("Canonical_Visits")
species_long <- read_final_csv("Species_OTU_Presence_Long")
genus_long <- read_final_csv("Genus_Presence_Long")
family_long <- read_final_csv("Family_Presence_Long")
species_matrix <- read_final_csv("Species_OTU_SitePresence_Matrix")
genus_matrix <- read_final_csv("Genus_SitePresence_Matrix")
family_matrix <- read_final_csv("Family_SitePresence_Matrix")

# The dissertation analytical universe contains Whale and Wood only. Other
# organic-fall records may exist in the master database but are not analysed.
stopifnot(all(site_universe$Fall_Type %in% c("Whale", "Wood")))
stopifnot(!anyDuplicated(site_universe$Node_ID))

# Presence matrices use five metadata columns followed by taxon columns.
# A zero means "taxon not recorded at this node". It must NOT be described as
# confirmed ecological absence in the dissertation or in downstream scripts.
validate_presence_matrix <- function(x, name) {
  taxon_data <- x[, -(1:5)]
  values <- unique(unlist(taxon_data, use.names = FALSE))
  values <- values[!is.na(values)]
  if (!all(values %in% c(0, 1))) stop(name, " contains non-binary values.")
  node_ids <- paste(x$Site_ID, x$Fall_Type, sep = "__")
  if (anyDuplicated(node_ids)) stop(name, " contains duplicated Site x Fall_Type nodes.")
  invisible(TRUE)
}
validate_presence_matrix(species_matrix, "Species/OTU matrix")
validate_presence_matrix(genus_matrix, "Genus matrix")
validate_presence_matrix(family_matrix, "Family matrix")

validation_summary <- tibble(
  Item = c("Analytical nodes", "Whale nodes", "Wood nodes", "Canonical events",
           "Canonical occurrences", "Species/OTU taxa", "Genera", "Families"),
  Value = c(nrow(site_universe), sum(site_universe$Fall_Type == "Whale"),
            sum(site_universe$Fall_Type == "Wood"), nrow(events), nrow(occurrences),
            ncol(species_matrix)-5, ncol(genus_matrix)-5, ncol(family_matrix)-5)
)
write_csv(validation_summary, file.path(OUT_DIR, "Diagnostics", "validation_summary.csv"))

message("Read-only validation completed. Review Outputs/Diagnostics/validation_summary.csv.")
