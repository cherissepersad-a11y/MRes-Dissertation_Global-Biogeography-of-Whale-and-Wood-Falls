# ============================================================================
# 06_network_site_similarity_figure.R
#
# Main-text Figure 6:
# Observed weighted site-similarity networks across taxonomic resolutions.
#
# IMPORTANT:
# - Uses the existing positive-Jaccard edge tables produced by script 05.
# - Reproduces the observed weighted Louvain analysis used in script 07.
# - Does NOT run null-model permutations.
# - Validates observed topology against network_observed_summary.csv.
# ============================================================================

source("R/00_setup_and_validate.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(igraph)
  library(ggplot2)
  library(ggraph)
  library(patchwork)
  library(scales)
})

dir.create("Outputs/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("Outputs/Tables", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Finalised inputs
# ---------------------------------------------------------------------------

edge_files <- c(
  Species_OTU = "Outputs/Tables/Species_OTU_Combined_site_similarity_edges.csv",
  Genus       = "Outputs/Tables/Genus_Combined_site_similarity_edges.csv",
  Family      = "Outputs/Tables/Family_Combined_site_similarity_edges.csv"
)

obs <- read_csv(
  "Outputs/Tables/network_observed_summary.csv",
  show_col_types = FALSE
)

# Deterministic seeds used for the observed network analyses.
seeds <- c(
  Species_OTU = 2026092203,
  Genus       = 2026092202,
  Family      = 2026092201
)

pretty_resolution <- c(
  Species_OTU = "Species/OTU",
  Genus       = "Genus",
  Family      = "Family"
)

# ---------------------------------------------------------------------------
# 2. Build exactly the observed network used for modularity analysis
# ---------------------------------------------------------------------------

build_observed_network <- function(resolution) {
  
  edges <- read_csv(
    edge_files[[resolution]],
    show_col_types = FALSE
  ) %>%
    filter(
      is.finite(Weight),
      Weight > 0
    )
  
  g <- graph_from_data_frame(
    edges,
    directed = FALSE
  )
  
  # Match the deterministic observed analysis.
  set.seed(seeds[[resolution]])
  
  cl <- cluster_louvain(
    g,
    weights = E(g)$Weight
  )
  
  V(g)$Module <- membership(cl)
  V(g)$Fall_Type <- ifelse(
    grepl("__Whale$", V(g)$name),
    "Whale",
    ifelse(
      grepl("__Wood$", V(g)$name),
      "Wood",
      "Unknown"
    )
  )
  
  # -------------------------------------------------------------------
  # Validate against the already finalised/verified observed statistics
  # --------------------------------------------------------------------
  
  expected <- obs %>%
    filter(Resolution == resolution)
  
  if (nrow(expected) != 1) {
    stop(
      paste(
        "Could not identify exactly one observed-summary row for",
        resolution
      )
    )
  }
  
  got_nodes <- vcount(g)
  got_edges <- ecount(g)
  got_components <- components(g)$no
  got_modules <- length(unique(membership(cl)))
  got_q <- modularity(cl)
  
  cat(
    "\n", pretty_resolution[[resolution]], " validation:\n",
    "  Nodes:      ", got_nodes, "\n",
    "  Edges:      ", got_edges, "\n",
    "  Components: ", got_components, "\n",
    "  Modules:    ", got_modules, "\n",
    "  Modularity: ", format(got_q, digits = 10), "\n",
    sep = ""
  )
  
  if (got_nodes != expected$Nodes) {
    stop(paste("Node validation failed for", resolution))
  }
  
  if (got_edges != expected$Edges) {
    stop(paste("Edge validation failed for", resolution))
  }
  
  if (got_components != expected$Components) {
    stop(paste("Component validation failed for", resolution))
  }
  
  if (got_modules != expected$Modules) {
    stop(paste("Module validation failed for", resolution))
  }
  
  if (!isTRUE(all.equal(
    got_q,
    expected$Modularity,
    tolerance = 1e-7
  ))) {
    stop(
      paste(
        "Modularity validation failed for",
        resolution,
        "- obtained",
        got_q,
        "expected",
        expected$Modularity
      )
    )
  }
  
  list(
    graph = g,
    community = cl,
    q = got_q,
    modules = got_modules
  )
}

networks <- lapply(
  names(edge_files),
  build_observed_network
)

names(networks) <- names(edge_files)

cat("\nAll observed networks reproduced successfully.\n")

# ---------------------------------------------------------------------------
# 3. Save observed module membership
#
# Useful for provenance and later interpretation.
# ---------------------------------------------------------------------------

module_table <- bind_rows(
  lapply(
    names(networks),
    function(resolution) {
      
      g <- networks[[resolution]]$graph
      
      tibble(
        Resolution = pretty_resolution[[resolution]],
        Node = V(g)$name,
        Fall_Type = V(g)$Fall_Type,
        Module = as.integer(V(g)$Module),
        Degree = degree(g),
        Strength = strength(
          g,
          weights = E(g)$Weight
        )
      )
    }
  )
)

write_csv(
  module_table,
  "Outputs/Tables/network_observed_module_membership.csv"
)

# ---------------------------------------------------------------------------
# 4. Plotting function
#
# IMPORTANT:
# The analytical network is NOT thresholded for display.
# Every positive-similarity edge remains present.
#
# Edge alpha/width only alter visual prominence.
# ---------------------------------------------------------------------------

make_network_panel <- function(
    resolution,
    panel_letter
) {
  
  g <- networks[[resolution]]$graph
  
  # Deterministic force-directed layout.
  set.seed(seeds[[resolution]])
  
  lay <- create_layout(
    g,
    layout = "fr",
    weights = E(g)$Weight,
    niter = 1500
  )
  
  # Convert module to factor so colour represents community identity,
  # not a continuous numerical quantity.
  lay$Module_plot <- factor(lay$Module)
  
  q_value <- networks[[resolution]]$q
  n_modules <- networks[[resolution]]$modules
  
  ggraph(lay) +
    
    geom_edge_link(
      aes(
        edge_alpha = Weight,
        edge_width = Weight
      ),
      colour = "grey55",
      show.legend = FALSE
    ) +
    
    scale_edge_alpha(
      range = c(0.03, 0.35)
    ) +
    
    scale_edge_width(
      range = c(0.05, 0.65)
    ) +
    
    geom_node_point(
      aes(
        colour = Module_plot,
        shape = Fall_Type
      ),
      size = 2.15,
      alpha = 0.92,
      stroke = 0.25
    ) +
    
    scale_shape_manual(
      values = c(
        Whale = 16,
        Wood = 17,
        Unknown = 15
      ),
      name = "Fall type"
    ) +
    
    guides(
      colour = "none",
      shape = guide_legend(
        override.aes = list(
          size = 4,
          alpha = 1
        )
      )
    ) +
    
    labs(
      title = paste0(
        panel_letter,
        "  ",
        pretty_resolution[[resolution]]
      ),
      subtitle = paste0(
        "Q = ",
        sprintf("%.3f", q_value),
        "   |   ",
        n_modules,
        " modules"
      )
    ) +
    
    theme_void(base_size = 11) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = 9.5,
        hjust = 0,
        margin = margin(b = 4)
      ),
      legend.position = "bottom",
      legend.title = element_text(
        face = "bold"
      ),
      plot.margin = margin(
        8, 8, 8, 8
      )
    )
}

# ---------------------------------------------------------------------------
# 5. Build Figure 6
# ---------------------------------------------------------------------------

p_species <- make_network_panel(
  "Species_OTU",
  "A"
)

p_genus <- make_network_panel(
  "Genus",
  "B"
)

p_family <- make_network_panel(
  "Family",
  "C"
)

figure6 <- (
  p_species |
    p_genus |
    p_family
) +
  plot_layout(
    guides = "collect",
    widths = c(1, 1, 1)
  ) +
  plot_annotation(
    title = "Observed site-similarity network structure across taxonomic resolutions",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0.5
      )
    )
  ) &
  theme(
    legend.position = "bottom"
  )

# ---------------------------------------------------------------------------
# 6. Export
# ---------------------------------------------------------------------------

ggsave(
  "Outputs/Figures/Figure_6_site_similarity_network.png",
  figure6,
  width = 15,
  height = 6.3,
  dpi = 450,
  bg = "white"
)

ggsave(
  "Outputs/Figures/Figure_6_site_similarity_network.pdf",
  figure6,
  width = 15,
  height = 6.3,
  bg = "white"
)

cat(
  "\nFigure 6 written to Outputs/Figures/.\n",
  "Observed module membership written to:\n",
  "  Outputs/Tables/network_observed_module_membership.csv\n",
  "No null-model permutations were rerun.\n",
  sep = ""
)
