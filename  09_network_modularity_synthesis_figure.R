# 09_network_modularity_synthesis_figure.R
#
# Chapter 2 synthesis:
# observed network modularity relative to incidence-preserving null models,
# plus descriptive changes in network topology across taxonomic resolution.
#
# IMPORTANT:
# This script READS completed null-model outputs.
# It does NOT rerun permutations.

source("R/00_setup_and_validate.R")

# ------------------------------------------------------------
# 1. Read completed observed and null-model results
# ------------------------------------------------------------

observed <- readr::read_csv(
  "Outputs/Tables/network_observed_summary.csv",
  show_col_types = FALSE
)

species_null <- readr::read_csv(
  "Outputs/Tables/species_modularity_null_999.csv",
  show_col_types = FALSE
)

genus_null <- readr::read_csv(
  "Outputs/Tables/genus_modularity_null_999.csv",
  show_col_types = FALSE
)

family_null <- readr::read_csv(
  "Outputs/Tables/family_modularity_null_999.csv",
  show_col_types = FALSE
)

species_summary <- readr::read_csv(
  "Outputs/Tables/species_modularity_null_summary.csv",
  show_col_types = FALSE
)

genus_summary <- readr::read_csv(
  "Outputs/Tables/genus_modularity_null_summary.csv",
  show_col_types = FALSE
)

family_summary <- readr::read_csv(
  "Outputs/Tables/family_modularity_null_summary.csv",
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 2. Inspect column names
# ------------------------------------------------------------
#
# The null-model files were generated earlier, so this section
# identifies the modularity column rather than assuming its name.

find_modularity_column <- function(x) {
  candidates <- names(x)[
    grepl("modular|^Q$|null_q", names(x), ignore.case = TRUE)
  ]
  
  if (length(candidates) == 0) {
    stop(
      "Could not identify modularity column. Columns are: ",
      paste(names(x), collapse = ", ")
    )
  }
  
  candidates[1]
}

species_q_col <- find_modularity_column(species_null)
genus_q_col   <- find_modularity_column(genus_null)
family_q_col  <- find_modularity_column(family_null)

message("Species null modularity column: ", species_q_col)
message("Genus null modularity column: ", genus_q_col)
message("Family null modularity column: ", family_q_col)

# ------------------------------------------------------------
# 3. Standardise null distributions
# ------------------------------------------------------------

null_data <- dplyr::bind_rows(
  species_null %>%
    dplyr::transmute(
      Resolution = "Species/OTU",
      Null_Q = .data[[species_q_col]]
    ),
  
  genus_null %>%
    dplyr::transmute(
      Resolution = "Genus",
      Null_Q = .data[[genus_q_col]]
    ),
  
  family_null %>%
    dplyr::transmute(
      Resolution = "Family",
      Null_Q = .data[[family_q_col]]
    )
) %>%
  dplyr::filter(is.finite(Null_Q)) %>%
  dplyr::mutate(
    Resolution = factor(
      Resolution,
      levels = c("Species/OTU", "Genus", "Family")
    )
  )

# ------------------------------------------------------------
# 4. Use the LOCKED observed values
# ------------------------------------------------------------
#
# These values were independently validated against
# network_observed_summary.csv after the 999-permutation runs.

observed_plot <- tibble::tribble(
  ~Resolution,   ~Observed_Q, ~Null_mean, ~Empirical_p,
  "Species/OTU",  0.7131087,   0.6216540,  0.001,
  "Genus",        0.4181168,   0.3766254,  0.064,
  "Family",       0.5287217,   0.4259086,  0.001
) %>%
  dplyr::mutate(
    Resolution = factor(
      Resolution,
      levels = c("Species/OTU", "Genus", "Family")
    ),
    annotation = paste0(
      "Observed Q = ", sprintf("%.3f", Observed_Q),
      "\np = ", sprintf("%.3f", Empirical_p)
    )
  )

# Confirm each null file contains all 999 valid permutations.
null_check <- null_data %>%
  dplyr::count(Resolution, name = "Valid_permutations")

print(null_check)

stopifnot(all(null_check$Valid_permutations == 999))

# ------------------------------------------------------------
# 5. Panels A-C: observed modularity against null distributions
# ------------------------------------------------------------

make_null_panel <- function(resolution_name) {
  
  nd <- null_data %>%
    dplyr::filter(Resolution == resolution_name)
  
  od <- observed_plot %>%
    dplyr::filter(Resolution == resolution_name)
  
  ggplot2::ggplot(nd, ggplot2::aes(x = Null_Q)) +
    
    ggplot2::geom_histogram(
      ggplot2::aes(y = after_stat(density)),
      bins = 30,
      fill = "grey75",
      colour = "white",
      linewidth = 0.25
    ) +
    
    ggplot2::geom_density(
      linewidth = 0.8,
      colour = "grey25"
    ) +
    
    ggplot2::geom_vline(
      data = od,
      ggplot2::aes(xintercept = Observed_Q),
      colour = "#D55E00",
      linewidth = 1.1
    ) +
    
    ggplot2::annotate(
      "text",
      x = od$Observed_Q,
      y = Inf,
      label = od$annotation,
      hjust = 1.05,
      vjust = 1.25,
      size = 3.5
    ) +
    
    ggplot2::labs(
      title = resolution_name,
      x = "Modularity (Q)",
      y = "Density"
    ) +
    
    ggplot2::theme_bw(base_size = 11) +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      )
    )
}

p7a <- make_null_panel("Species/OTU")
p7b <- make_null_panel("Genus")
p7c <- make_null_panel("Family")

# ------------------------------------------------------------
# 6. Panel D: descriptive network topology across resolutions
# ------------------------------------------------------------
#
# These are observed network properties. They are descriptive;
# this panel is NOT a statistical test among taxonomic resolutions.

topology <- tibble::tribble(
  ~Resolution,   ~Linked_nodes, ~Density,    ~Components,
  "Species/OTU", 747,           0.03375791,  29,
  "Genus",       269,           0.15846418,   5,
  "Family",      110,           0.26538782,   3
) %>%
  dplyr::mutate(
    Resolution = factor(
      Resolution,
      levels = c("Species/OTU", "Genus", "Family")
    )
  )

# Because nodes, density and components have very different scales,
# display them as a compact labelled summary rather than forcing them
# onto an inappropriate common quantitative axis.

topology_long <- topology %>%
  tidyr::pivot_longer(
    cols = c(Linked_nodes, Density, Components),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Metric = factor(
      Metric,
      levels = c("Linked_nodes", "Components", "Density"),
      labels = c("Linked nodes", "Linked components", "Network density")
    ),
    Label = dplyr::case_when(
      Metric == "Network density" ~ sprintf("%.3f", Value),
      TRUE ~ sprintf("%.0f", Value)
    )
  )

p7d <- ggplot2::ggplot(
  topology_long,
  ggplot2::aes(x = Resolution, y = Metric)
) +
  
  ggplot2::geom_tile(
    fill = "grey95",
    colour = "white",
    linewidth = 1
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(label = Label),
    size = 4
  ) +
  
  ggplot2::labs(
    title = "Observed network topology",
    x = NULL,
    y = NULL
  ) +
  
  ggplot2::theme_bw(base_size = 11) +
  
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5
    )
  )

# ------------------------------------------------------------
# 7. Assemble Figure 7
# ------------------------------------------------------------

fig7 <- (
  p7a | p7b | p7c
) / p7d +
  patchwork::plot_layout(
    heights = c(2.2, 1)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  )

# ------------------------------------------------------------
# 8. Export publication-resolution files
# ------------------------------------------------------------

ggplot2::ggsave(
  "Outputs/Figures/Figure_7_network_modularity_synthesis.png",
  fig7,
  width = 12,
  height = 7.5,
  dpi = 450
)

ggplot2::ggsave(
  "Outputs/Figures/Figure_7_network_modularity_synthesis.pdf",
  fig7,
  width = 12,
  height = 7.5
)

message(
  "Figure 7 written to Outputs/Figures/. ",
  "No null-model permutations were rerun."
)