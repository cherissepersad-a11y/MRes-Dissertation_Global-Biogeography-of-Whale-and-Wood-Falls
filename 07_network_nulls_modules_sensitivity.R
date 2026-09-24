# 07_network_nulls_modules_sensitivity.R
#
# Network robustness, threshold sensitivity, and incidence-preserving
# null-model analysis for the MRes organic-fall dissertation.
#
# This script is separate from 99_run_all.R because the 999-permutation
# null models can require substantial computation time.
#
# USAGE
# -----
# First source this script:
#
#   source("R/07_network_nulls_modules_sensitivity.R")
#
# Then run any resolution independently:
#
#   run_null_analysis("Family")
#   run_null_analysis("Genus")
#   run_null_analysis("Species_OTU")
#
# Or reproduce all three sequentially:
#
#   run_all_null_analyses()
#
# Each resolution uses its own fixed random seed, so results are
# reproducible independently of run order.
#
# Null-model method:
#   - binary site-by-taxon incidence matrix
#   - empty sites removed
#   - vegan quasiswap incidence-preserving randomisation
#   - Jaccard site similarity
#   - positive-similarity site network
#   - weighted Louvain community detection
#   - 999 permutations
#   - upper-tail empirical P value with +1 correction


# =========================================================================
# 0. SETUP
# =========================================================================

source("R/00_setup_and_validate.R")
source("R/05_network_site_similarity.R")

NPERM <- 999

NULL_SEEDS <- c(
  Family = 2026092201,
  Genus = 2026092202,
  Species_OTU = 2026092203
)


# =========================================================================
# 1. BASIC NETWORK ANALYSIS
# =========================================================================

analyse_network <- function(edges, threshold = 0) {
  
  e <- edges %>%
    filter(
      is.finite(Weight),
      Weight > threshold
    )
  
  if (nrow(e) == 0) {
    return(NULL)
  }
  
  g <- igraph::graph_from_data_frame(
    e,
    directed = FALSE
  )
  
  if (igraph::ecount(g) == 0) {
    return(NULL)
  }
  
  cl <- igraph::cluster_louvain(
    g,
    weights = igraph::E(g)$Weight
  )
  
  tibble(
    Nodes = igraph::vcount(g),
    Edges = igraph::ecount(g),
    Density = igraph::edge_density(g),
    Mean_degree = mean(
      igraph::degree(g)
    ),
    Clustering = igraph::transitivity(
      g,
      type = "global"
    ),
    Components = igraph::components(g)$no,
    Modularity = igraph::modularity(cl),
    Modules = length(
      unique(
        igraph::membership(cl)
      )
    )
  )
}


# =========================================================================
# 2. THRESHOLD-SENSITIVITY ANALYSIS
# =========================================================================

run_threshold_sensitivity <- function(edges, seed) {
  
  qs <- c(
    0,
    0.25,
    0.50,
    0.75,
    0.90
  )
  
  bind_rows(
    lapply(
      seq_along(qs),
      function(i) {
        
        q <- qs[i]
        
        threshold <- if (q == 0) {
          0
        } else {
          unname(
            quantile(
              edges$Weight,
              q,
              na.rm = TRUE
            )
          )
        }
        
        # Independent deterministic seed for each threshold.
        set.seed(seed + i)
        
        z <- analyse_network(
          edges,
          threshold
        )
        
        if (is.null(z)) {
          return(tibble())
        }
        
        z %>%
          mutate(
            Quantile = q,
            Threshold = threshold
          )
      }
    )
  )
}


# =========================================================================
# 3. NULL-MODEL FUNCTION
# =========================================================================

null_modularity <- function(
    presence_matrix,
    nperm = NPERM,
    resolution = "Unknown"
) {
  
  m <- as.matrix(
    presence_matrix[, -(1:5)]
  )
  
  # Retain sites with at least one recorded taxon.
  m <- m[
    rowSums(m) > 0,
    ,
    drop = FALSE
  ]
  
  # Current vegan API.
  # "quasiswap" preserves row and column totals of the incidence matrix.
  null_model <- vegan::nullmodel(
    m,
    method = "quasiswap"
  )
  
  vals <- rep(
    NA_real_,
    nperm
  )
  
  cat(
    "\nStarting ",
    resolution,
    " null model: ",
    nperm,
    " permutations\n",
    sep = ""
  )
  
  for (i in seq_len(nperm)) {
    
    mr <- simulate(
      null_model,
      nsim = 1
    )
    
    # simulate.nullmodel may return a 3D array when nsim = 1.
    # Convert the single simulated matrix back to an ordinary matrix.
    if (length(dim(mr)) == 3) {
      mr <- mr[, , 1]
    }
    
    mr <- as.matrix(mr)
    
    d <- as.matrix(
      vegan::vegdist(
        mr,
        method = "jaccard",
        binary = TRUE
      )
    )
    
    s <- 1 - d
    diag(s) <- 0
    
    idx <- which(
      upper.tri(s) & s > 0,
      arr.ind = TRUE
    )
    
    if (nrow(idx) == 0) {
      
      vals[i] <- NA_real_
      
    } else {
      
      e <- tibble(
        Source = rownames(s)[idx[, 1]],
        Target = colnames(s)[idx[, 2]],
        Weight = s[idx]
      )
      
      a <- analyse_network(
        e,
        threshold = 0
      )
      
      vals[i] <- if (is.null(a)) {
        NA_real_
      } else {
        a$Modularity
      }
    }
    
    # Progress output so long analyses do not appear frozen.
    if (
      i %% 25 == 0 ||
      i == nperm
    ) {
      
      cat(
        resolution,
        ": ",
        i,
        " / ",
        nperm,
        " permutations complete\n",
        sep = ""
      )
      
      flush.console()
    }
  }
  
  vals
}


# =========================================================================
# 4. READ OBSERVED COMBINED-NETWORK EDGE TABLES
# =========================================================================

species_edges <- readr::read_csv(
  "Outputs/Tables/Species_OTU_Combined_site_similarity_edges.csv",
  show_col_types = FALSE
)

genus_edges <- readr::read_csv(
  "Outputs/Tables/Genus_Combined_site_similarity_edges.csv",
  show_col_types = FALSE
)

family_edges <- readr::read_csv(
  "Outputs/Tables/Family_Combined_site_similarity_edges.csv",
  show_col_types = FALSE
)


# =========================================================================
# 5. RESOLUTION LOOKUP
# =========================================================================

get_resolution_objects <- function(resolution) {
  
  valid_resolutions <- c(
    "Family",
    "Genus",
    "Species_OTU"
  )
  
  if (!resolution %in% valid_resolutions) {
    
    stop(
      paste0(
        "Unknown resolution '",
        resolution,
        "'. Use one of: ",
        paste(
          valid_resolutions,
          collapse = ", "
        ),
        "."
      )
    )
  }
  
  if (resolution == "Family") {
    
    return(
      list(
        matrix = family_matrix,
        edges = family_edges,
        seed = unname(
          NULL_SEEDS["Family"]
        )
      )
    )
  }
  
  if (resolution == "Genus") {
    
    return(
      list(
        matrix = genus_matrix,
        edges = genus_edges,
        seed = unname(
          NULL_SEEDS["Genus"]
        )
      )
    )
  }
  
  list(
    matrix = species_matrix,
    edges = species_edges,
    seed = unname(
      NULL_SEEDS["Species_OTU"]
    )
  )
}


# =========================================================================
# 6. OBSERVED NETWORKS
# =========================================================================

set.seed(
  NULL_SEEDS["Species_OTU"]
)

species_observed <- analyse_network(
  species_edges
)

set.seed(
  NULL_SEEDS["Genus"]
)

genus_observed <- analyse_network(
  genus_edges
)

set.seed(
  NULL_SEEDS["Family"]
)

family_observed <- analyse_network(
  family_edges
)


observed_network_summary <- bind_rows(
  Species_OTU = species_observed,
  Genus = genus_observed,
  Family = family_observed,
  .id = "Resolution"
)

readr::write_csv(
  observed_network_summary,
  "Outputs/Tables/network_observed_summary.csv"
)


# =========================================================================
# 7. FULL-NETWORK NODE / COMPONENT ACCOUNTING
# =========================================================================

nonempty_nodes <- tibble(
  Resolution = c(
    "Species_OTU",
    "Genus",
    "Family"
  ),
  
  Nonempty_nodes = c(
    
    sum(
      rowSums(
        as.matrix(
          species_matrix[, -(1:5)]
        )
      ) > 0
    ),
    
    sum(
      rowSums(
        as.matrix(
          genus_matrix[, -(1:5)]
        )
      ) > 0
    ),
    
    sum(
      rowSums(
        as.matrix(
          family_matrix[, -(1:5)]
        )
      ) > 0
    )
  )
)


network_accounting <- observed_network_summary %>%
  select(
    Resolution,
    Linked_nodes = Nodes,
    Linked_components = Components
  ) %>%
  left_join(
    nonempty_nodes,
    by = "Resolution"
  ) %>%
  mutate(
    Isolates =
      Nonempty_nodes -
      Linked_nodes,
    
    Full_components =
      Linked_components +
      Isolates
  ) %>%
  select(
    Resolution,
    Nonempty_nodes,
    Linked_nodes,
    Isolates,
    Linked_components,
    Full_components
  )


readr::write_csv(
  network_accounting,
  "Outputs/Tables/network_node_component_accounting.csv"
)


# =========================================================================
# 8. THRESHOLD SENSITIVITY
# =========================================================================

species_threshold <- run_threshold_sensitivity(
  species_edges,
  seed = NULL_SEEDS["Species_OTU"]
) %>%
  mutate(
    Resolution = "Species_OTU"
  )


genus_threshold <- run_threshold_sensitivity(
  genus_edges,
  seed = NULL_SEEDS["Genus"]
) %>%
  mutate(
    Resolution = "Genus"
  )


family_threshold <- run_threshold_sensitivity(
  family_edges,
  seed = NULL_SEEDS["Family"]
) %>%
  mutate(
    Resolution = "Family"
  )


threshold_sensitivity <- bind_rows(
  species_threshold,
  genus_threshold,
  family_threshold
) %>%
  select(
    Resolution,
    everything()
  )


readr::write_csv(
  threshold_sensitivity,
  "Outputs/Tables/network_threshold_sensitivity.csv"
)


# =========================================================================
# 9. RUN ONE NULL ANALYSIS
# =========================================================================

run_null_analysis <- function(
    resolution,
    nperm = NPERM
) {
  
  if (nperm != NPERM) {
    
    warning(
      paste0(
        "The dissertation analysis specifies ",
        NPERM,
        " permutations. ",
        "You requested ",
        nperm,
        "."
      )
    )
  }
  
  obj <- get_resolution_objects(
    resolution
  )
  
  presence_matrix <- obj$matrix
  edges <- obj$edges
  seed <- obj$seed
  
  
  cat(
    "\n============================================================\n",
    "NULL-MODEL ANALYSIS: ",
    resolution,
    "\n",
    "Permutations: ",
    nperm,
    "\n",
    "Seed: ",
    seed,
    "\n",
    "============================================================\n",
    sep = ""
  )
  
  
  # -----------------------------------------------------------------------
  # Observed modularity
  # -----------------------------------------------------------------------
  
  set.seed(seed)
  
  observed <- analyse_network(
    edges,
    threshold = 0
  )
  
  if (is.null(observed)) {
    
    stop(
      paste0(
        "Observed ",
        resolution,
        " network contains no positive-similarity edges."
      )
    )
  }
  
  
  # -----------------------------------------------------------------------
  # Null distribution
  # -----------------------------------------------------------------------
  
  # Reset the seed before the null simulation so it is independent
  # of the stochastic Louvain calculation for the observed network.
  
  set.seed(seed)
  
  null_values <- null_modularity(
    presence_matrix,
    nperm = nperm,
    resolution = resolution
  )
  
  
  # -----------------------------------------------------------------------
  # Null-model summary
  # -----------------------------------------------------------------------
  
  valid <- null_values[
    is.finite(null_values)
  ]
  
  if (length(valid) == 0) {
    
    stop(
      paste0(
        "No valid null modularity values were produced for ",
        resolution,
        "."
      )
    )
  }
  
  
  # Upper-tail empirical permutation P value with +1 correction.
  empirical_p <- (
    sum(
      valid >= observed$Modularity
    ) + 1
  ) / (
    length(valid) + 1
  )
  
  
  summary <- tibble(
    Resolution = resolution,
    Seed = seed,
    
    Observed_nodes = observed$Nodes,
    Observed_edges = observed$Edges,
    Observed_density = observed$Density,
    Observed_mean_degree = observed$Mean_degree,
    Observed_clustering = observed$Clustering,
    Observed_linked_components = observed$Components,
    Observed_modularity = observed$Modularity,
    Observed_modules = observed$Modules,
    
    Requested_permutations = nperm,
    Valid_permutations = length(valid),
    
    Null_mean = mean(
      valid
    ),
    
    Null_SD = sd(
      valid
    ),
    
    Null_min = min(
      valid
    ),
    
    Null_max = max(
      valid
    ),
    
    Null_q025 = unname(
      quantile(
        valid,
        0.025
      )
    ),
    
    Null_median = median(
      valid
    ),
    
    Null_q975 = unname(
      quantile(
        valid,
        0.975
      )
    ),
    
    Null_ge_observed = sum(
      valid >= observed$Modularity
    ),
    
    Empirical_p_upper = empirical_p
  )
  
  
  # -----------------------------------------------------------------------
  # Output filenames
  # -----------------------------------------------------------------------
  
  prefix <- switch(
    resolution,
    Family = "family",
    Genus = "genus",
    Species_OTU = "species"
  )
  
  
  rds_file <- file.path(
    "Outputs",
    "Diagnostics",
    paste0(
      prefix,
      "_modularity_null_",
      nperm,
      ".rds"
    )
  )
  
  
  csv_file <- file.path(
    "Outputs",
    "Tables",
    paste0(
      prefix,
      "_modularity_null_",
      nperm,
      ".csv"
    )
  )
  
  
  summary_file <- file.path(
    "Outputs",
    "Tables",
    paste0(
      prefix,
      "_modularity_null_summary.csv"
    )
  )
  
  
  # -----------------------------------------------------------------------
  # Save completed resolution immediately
  # -----------------------------------------------------------------------
  
  saveRDS(
    null_values,
    rds_file
  )
  
  
  readr::write_csv(
    tibble(
      Resolution = resolution,
      Permutation = seq_along(
        null_values
      ),
      Null_modularity = null_values
    ),
    csv_file
  )
  
  
  readr::write_csv(
    summary,
    summary_file
  )
  
  
  # -----------------------------------------------------------------------
  # Console report
  # -----------------------------------------------------------------------
  
  cat(
    "\n------------------------------------------------------------\n",
    resolution,
    " NULL MODEL COMPLETE\n",
    "------------------------------------------------------------\n",
    sep = ""
  )
  
  print(
    summary,
    width = Inf
  )
  
  cat(
    "\nSaved:\n",
    "  ",
    rds_file,
    "\n",
    "  ",
    csv_file,
    "\n",
    "  ",
    summary_file,
    "\n",
    sep = ""
  )
  
  
  invisible(
    list(
      observed = observed,
      null = null_values,
      summary = summary
    )
  )
}


# =========================================================================
# 10. COMBINE COMPLETED NULL-MODEL SUMMARIES
# =========================================================================

combine_null_summaries <- function() {
  
  summary_files <- c(
    Family =
      "Outputs/Tables/family_modularity_null_summary.csv",
    
    Genus =
      "Outputs/Tables/genus_modularity_null_summary.csv",
    
    Species_OTU =
      "Outputs/Tables/species_modularity_null_summary.csv"
  )
  
  
  existing <- summary_files[
    file.exists(
      summary_files
    )
  ]
  
  
  if (length(existing) == 0) {
    
    message(
      "No completed null-model summary files were found."
    )
    
    return(
      invisible(NULL)
    )
  }
  
  
  combined <- bind_rows(
    lapply(
      existing,
      readr::read_csv,
      show_col_types = FALSE
    )
  )
  
  
  readr::write_csv(
    combined,
    "Outputs/Tables/network_modularity_null_summary.csv"
  )
  
  
  cat(
    "\nCombined null-model summary updated:\n",
    "  Outputs/Tables/network_modularity_null_summary.csv\n",
    sep = ""
  )
  
  
  invisible(
    combined
  )
}


# =========================================================================
# 11. RUN ALL THREE NULL ANALYSES
# =========================================================================

run_all_null_analyses <- function(
    nperm = NPERM
) {
  
  cat(
    "\n============================================================\n",
    "RUNNING ALL NETWORK NULL MODELS\n",
    "Order: Family -> Genus -> Species/OTU\n",
    "Permutations per resolution: ",
    nperm,
    "\n",
    "============================================================\n",
    sep = ""
  )
  
  
  family_result <- run_null_analysis(
    "Family",
    nperm = nperm
  )
  
  combine_null_summaries()
  
  
  genus_result <- run_null_analysis(
    "Genus",
    nperm = nperm
  )
  
  combine_null_summaries()
  
  
  species_result <- run_null_analysis(
    "Species_OTU",
    nperm = nperm
  )
  
  combine_null_summaries()
  
  
  cat(
    "\n============================================================\n",
    "ALL NETWORK NULL MODELS COMPLETED SUCCESSFULLY\n",
    "============================================================\n"
  )
  
  
  invisible(
    list(
      Family = family_result,
      Genus = genus_result,
      Species_OTU = species_result
    )
  )
}


# =========================================================================
# 12. READY MESSAGE
# =========================================================================

cat(
  "\n============================================================\n",
  "SCRIPT 07 LOADED SUCCESSFULLY\n",
  "============================================================\n",
  "No null model has been started automatically.\n\n",
  "Run one analysis with:\n",
  '  run_null_analysis("Family")\n',
  '  run_null_analysis("Genus")\n',
  '  run_null_analysis("Species_OTU")\n\n',
  "Or run all three sequentially with:\n",
  "  run_all_null_analyses()\n",
  "============================================================\n",
  sep = ""
)