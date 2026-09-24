# ============================================================
# INTRODUCTORY FIGURE
# Whale- and wood-fall Species/OTU overlap
# UpSet-style representation
# ============================================================

library(ggplot2)
library(patchwork)

# ---- Derive values directly from frozen matrix --------------

taxa_cols <- names(species_matrix)[6:ncol(species_matrix)]

whale_present <- colSums(
  species_matrix[
    species_matrix$Fall_Type == "Whale",
    taxa_cols,
    drop = FALSE
  ],
  na.rm = TRUE
) > 0

wood_present <- colSums(
  species_matrix[
    species_matrix$Fall_Type == "Wood",
    taxa_cols,
    drop = FALSE
  ],
  na.rm = TRUE
) > 0

whale_only <- sum(whale_present & !wood_present)
shared     <- sum(whale_present & wood_present)
wood_only  <- sum(!whale_present & wood_present)

whale_total <- sum(whale_present)
wood_total  <- sum(wood_present)
union_total <- sum(whale_present | wood_present)

# Frozen-data validation
stopifnot(
  whale_only == 687,
  shared == 80,
  wood_only == 732,
  whale_total == 767,
  wood_total == 812,
  union_total == 1499
)

# ============================================================
# A. Intersection sizes
# ============================================================

intersection_df <- data.frame(
  intersection = factor(
    c("Whale only", "Shared", "Wood only"),
    levels = c("Whale only", "Shared", "Wood only")
  ),
  n = c(whale_only, shared, wood_only)
)

p_intersections <- ggplot(
  intersection_df,
  aes(x = intersection, y = n)
) +
  geom_col(
    width = 0.62,
    fill = "grey30"
  ) +
  geom_text(
    aes(label = n),
    vjust = -0.45,
    fontface = "bold",
    size = 5
  ) +
  scale_y_continuous(
    limits = c(0, 820),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = NULL,
    y = "Number of Species/OTUs"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(5, 15, 0, 5)
  )

# ============================================================
# B. Membership matrix
# ============================================================

membership_df <- data.frame(
  intersection = factor(
    rep(
      c("Whale only", "Shared", "Wood only"),
      each = 2
    ),
    levels = c("Whale only", "Shared", "Wood only")
  ),
  fall_type = factor(
    rep(c("Whale falls", "Wood falls"), 3),
    levels = c("Whale falls", "Wood falls")
  ),
  present = c(
    TRUE,  FALSE,   # Whale only
    TRUE,  TRUE,    # Shared
    FALSE, TRUE     # Wood only
  )
)

# Lines connecting the shared membership only
shared_line <- data.frame(
  intersection = factor(
    c("Shared", "Shared"),
    levels = c("Whale only", "Shared", "Wood only")
  ),
  fall_type = factor(
    c("Whale falls", "Wood falls"),
    levels = c("Whale falls", "Wood falls")
  )
)

p_matrix <- ggplot(
  membership_df,
  aes(x = intersection, y = fall_type)
) +
  geom_line(
    data = shared_line,
    aes(group = 1),
    linewidth = 1,
    colour = "grey30"
  ) +
  geom_point(
    aes(
      fill = present,
      colour = present
    ),
    shape = 21,
    size = 5,
    stroke = 0.8
  ) +
  scale_fill_manual(
    values = c(
      `TRUE` = "grey20",
      `FALSE` = "white"
    )
  ) +
  scale_colour_manual(
    values = c(
      `TRUE` = "grey20",
      `FALSE` = "grey75"
    )
  ) +
  scale_y_discrete(
    limits = rev(levels(membership_df$fall_type))
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "none",
    plot.margin = margin(0, 15, 5, 5)
  )

# ============================================================
# C. Set sizes
# ============================================================

set_df <- data.frame(
  fall_type = factor(
    c("Whale falls", "Wood falls"),
    levels = c("Whale falls", "Wood falls")
  ),
  n = c(whale_total, wood_total)
)

p_sets <- ggplot(
  set_df,
  aes(x = n, y = fall_type)
) +
  geom_col(
    width = 0.55,
    fill = "grey55"
  ) +
  geom_text(
    aes(label = n),
    hjust = -0.25,
    fontface = "bold",
    size = 4.3
  ) +
  scale_x_continuous(
    limits = c(0, 900),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_discrete(
    limits = rev(levels(set_df$fall_type))
  ) +
  labs(
    x = "Total Species/OTUs",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(0, 5, 5, 5)
  )

# ============================================================
# Assemble figure
# ============================================================

right_panel <-
  p_intersections /
  p_matrix +
  plot_layout(heights = c(3.3, 1))

p_upset <-
  p_sets | right_panel +
  plot_layout(widths = c(1.05, 2.8))

p_upset <- p_upset +
  plot_annotation(
    title = "Whale- and wood-fall Species/OTU overlap",
    subtitle = "1,499 Species/OTUs recorded across the global dataset",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 16
      ),
      plot.subtitle = element_text(
        size = 12
      )
    )
  )

p_upset

# ---- Export --------------------------------------------------

ggsave(
  "Outputs/Figures/Figure_Intro_whale_wood_upset.png",
  p_upset,
  width = 9,
  height = 6,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Outputs/Figures/Figure_Intro_whale_wood_upset.pdf",
  p_upset,
  width = 9,
  height = 6,
  bg = "white"
)