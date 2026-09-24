# 03_community_overlap_beta_diversity.R

# Chapter 2: overlap, beta diversity and taxonomic-resolution comparison.

source("R/00_setup_and_validate.R")


# -------------------------------------------------------------------------
# Taxonomic overlap
# -------------------------------------------------------------------------

matrix_to_numeric <- function(x) {
  as.matrix(x[, -(1:5)])
}

pooled_overlap <- function(x, resolution) {
  m <- matrix_to_numeric(x)
  
  whale <- colSums(
    m[x$Fall_Type == "Whale", , drop = FALSE],
    na.rm = TRUE
  ) > 0
  
  wood <- colSums(
    m[x$Fall_Type == "Wood", , drop = FALSE],
    na.rm = TRUE
  ) > 0
  
  shared <- sum(whale & wood)
  union <- sum(whale | wood)
  
  tibble(
    Resolution = resolution,
    Whale_taxa = sum(whale),
    Wood_taxa = sum(wood),
    Shared_taxa = shared,
    Union_taxa = union,
    Jaccard_overlap = shared / union
  )
}

overlap <- bind_rows(
  pooled_overlap(species_matrix, "Species_OTU"),
  pooled_overlap(genus_matrix, "Genus"),
  pooled_overlap(family_matrix, "Family")
)

write_csv(
  overlap,
  "Outputs/Tables/taxonomic_overlap.csv"
)


# -------------------------------------------------------------------------
# Beta diversity
# -------------------------------------------------------------------------

beta_one <- function(x, resolution) {
  m <- matrix_to_numeric(x)
  
  keep <- rowSums(m) > 0
  m <- m[keep, , drop = FALSE]
  
  b <- betapart::beta.pair(
    m,
    index.family = "sorensen"
  )
  
  # Save pairwise components in long form.
  # beta.sim = turnover component
  # beta.sne = nestedness-resultant component
  # beta.sor = total dissimilarity
  
  as.data.frame(
    as.table(as.matrix(b$beta.sor))
  ) %>%
    rename(
      i = Var1,
      j = Var2,
      beta_sor = Freq
    ) %>%
    filter(
      as.integer(i) < as.integer(j)
    ) %>%
    mutate(
      beta_sim = as.matrix(b$beta.sim)[
        cbind(as.integer(i), as.integer(j))
      ],
      beta_sne = as.matrix(b$beta.sne)[
        cbind(as.integer(i), as.integer(j))
      ],
      Resolution = resolution
    )
}

beta <- bind_rows(
  beta_one(species_matrix, "Species_OTU"),
  beta_one(genus_matrix, "Genus"),
  beta_one(family_matrix, "Family")
)

write_csv(
  beta,
  "Outputs/Tables/beta_diversity_components.csv"
)


# -------------------------------------------------------------------------
# Chapter 2 taxonomic-resolution colour convention
#
# Use this mapping whenever colour represents taxonomic resolution:
# Species/OTU = light blue
# Genus       = medium blue
# Family      = dark blue
# -------------------------------------------------------------------------

resolution_palette <- c(
  "Species/OTU" = "#9ECAE1",
  "Genus"       = "#4292C6",
  "Family"      = "#08519C"
)


# -------------------------------------------------------------------------
# Figure 3: pooled whale-wood taxonomic overlap across resolutions
# -------------------------------------------------------------------------

overlap_plot <- overlap %>%
  mutate(
    Resolution = factor(
      Resolution,
      levels = c(
        "Species_OTU",
        "Genus",
        "Family"
      ),
      labels = c(
        "Species/OTU",
        "Genus",
        "Family"
      )
    )
  )

p3 <- ggplot(
  overlap_plot,
  aes(
    x = Resolution,
    y = Jaccard_overlap,
    fill = Resolution
  )
) +
  geom_col(
    width = 0.72
  ) +
  scale_fill_manual(
    values = resolution_palette,
    guide = "none"
  ) +
  geom_text(
    aes(
      label = scales::percent(
        Jaccard_overlap,
        accuracy = 0.1
      )
    ),
    vjust = -0.5,
    size = 4
  ) +
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    limits = c(0, 0.27),
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  theme_bw(
    base_size = 11
  ) +
  labs(
    x = NULL,
    y = "Whale-wood Jaccard overlap"
  )

ggsave(
  "Outputs/Figures/Figure_3_taxonomic_overlap.png",
  p3,
  width = 6,
  height = 4.5,
  dpi = 450
)

ggsave(
  "Outputs/Figures/Figure_3_taxonomic_overlap.pdf",
  p3,
  width = 6,
  height = 4.5
)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Figure 4: beta-diversity partition across taxonomic resolutions
#
# Means are used because the pairwise distributions are strongly concentrated
# at complete dissimilarity. Median beta.sor and beta.sim are therefore 1,
# and median beta.sne is 0 at all resolutions, which obscures the
# cross-resolution pattern.
# -------------------------------------------------------------------------

beta_summary <- beta %>%
  group_by(Resolution) %>%
  summarise(
    Turnover = mean(beta_sim, na.rm = TRUE),
    Nestedness = mean(beta_sne, na.rm = TRUE),
    Total = mean(beta_sor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Resolution = factor(
      Resolution,
      levels = c(
        "Species_OTU",
        "Genus",
        "Family"
      ),
      labels = c(
        "Species/OTU",
        "Genus",
        "Family"
      )
    )
  )

beta_plot <- beta_summary %>%
  select(
    Resolution,
    Turnover,
    Nestedness
  ) %>%
  pivot_longer(
    cols = c(Turnover, Nestedness),
    names_to = "Component",
    values_to = "Mean"
  ) %>%
  mutate(
    Component = factor(
      Component,
      levels = c("Nestedness", "Turnover")
    )
  )

beta_component_palette <- c(
  "Turnover" = "#4292C6",
  "Nestedness" = "#F28E8B"
)

p4 <- ggplot(
  beta_plot,
  aes(
    x = Resolution,
    y = Mean,
    fill = Component
  )
) +
  geom_col(
    width = 0.72
  ) +
  scale_fill_manual(
    values = beta_component_palette,
    labels = c(
      "Nestedness-resultant",
      "Turnover"
    ),
    name = "Beta-diversity component"
  ) +
  geom_text(
    data = beta_summary,
    aes(
      x = Resolution,
      y = Total,
      label = sprintf("%.3f", Total)
    ),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 4
  ) +
  scale_y_continuous(
    limits = c(0, 1.05),
    expand = expansion(
      mult = c(0, 0.01)
    )
  ) +
  theme_bw(
    base_size = 11
  ) +
  labs(
    x = NULL,
    y = "Mean pairwise Sørensen dissimilarity"
  )

ggsave(
  "Outputs/Figures/Figure_4_beta_partition.png",
  p4,
  width = 7,
  height = 4.5,
  dpi = 450
)

ggsave(
  "Outputs/Figures/Figure_4_beta_partition.pdf",
  p4,
  width = 7,
  height = 4.5
)