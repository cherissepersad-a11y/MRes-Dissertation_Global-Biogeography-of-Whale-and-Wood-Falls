# ============================================================
# INTRODUCTORY FIGURE
# Whale- and wood-fall Species/OTU overlap
# Simplified UpSet-style representation
# ============================================================

library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# 1. Derive values directly from finalised Species/OTU matrix
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 2. Validate against finalised analysis
# ------------------------------------------------------------

stopifnot(
  whale_only == 687,
  shared == 80,
  wood_only == 732,
  whale_total == 767,
  wood_total == 812,
  union_total == 1499,
  whale_only + shared + wood_only == union_total
)

cat(
  "\nValidated Species/OTU overlap:\n",
  "Whale only =", whale_only, "\n",
  "Shared     =", shared, "\n",
  "Wood only  =", wood_only, "\n",
  "Whale total =", whale_total, "\n",
  "Wood total  =", wood_total, "\n",
  "Union       =", union_total, "\n\n"
)

# ------------------------------------------------------------
# 3. Colour palette
# ------------------------------------------------------------

whale_col  <- "#4477AA"
shared_col <- "#8B7D6B"
wood_col   <- "#CC8844"

# ------------------------------------------------------------
# 4. Intersection-size data
# ------------------------------------------------------------

intersection_df <- data.frame(
  intersection = factor(
    c("Whale only", "Shared", "Wood only"),
    levels = c("Whale only", "Shared", "Wood only")
  ),
  n = c(
    whale_only,
    shared,
    wood_only
  )
)

# ------------------------------------------------------------
# 5. Main intersection-size bars
# ------------------------------------------------------------

p_intersections <- ggplot(
  intersection_df,
  aes(
    x = intersection,
    y = n,
    fill = intersection
  )
) +
  geom_col(
    width = 0.60
  ) +
  geom_text(
    aes(label = n),
    vjust = -0.45,
    fontface = "bold",
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "Whale only" = whale_col,
      "Shared" = shared_col,
      "Wood only" = wood_col
    )
  ) +
  scale_y_continuous(
    limits = c(0, 820),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = NULL,
    y = "Intersection size"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_text(
      size = 12,
      margin = margin(r = 10)
    ),
    axis.text.y = element_text(size = 10),
    legend.position = "none",
    plot.margin = margin(5, 15, 0, 5)
  )

# ------------------------------------------------------------
# 6. Membership matrix
# ------------------------------------------------------------

membership_df <- data.frame(
  intersection = factor(
    rep(
      c("Whale only", "Shared", "Wood only"),
      each = 2
    ),
    levels = c("Whale only", "Shared", "Wood only")
  ),
  fall_type = factor(
    rep(
      c("Whale falls", "Wood falls"),
      times = 3
    ),
    levels = c("Whale falls", "Wood falls")
  ),
  present = c(
    TRUE,  FALSE,   # Whale only
    TRUE,  TRUE,    # Shared
    FALSE, TRUE     # Wood only
  )
)

# Background/open circles for all matrix positions
p_matrix <- ggplot(
  membership_df,
  aes(
    x = intersection,
    y = fall_type
  )
) +
  
  # Open circles show absence
  geom_point(
    shape = 21,
    size = 5,
    stroke = 0.8,
    fill = "white",
    colour = "grey75"
  ) +
  
  # Connector for the shared intersection
  annotate(
    "segment",
    x = 2,
    xend = 2,
    y = 1,
    yend = 2,
    linewidth = 1,
    colour = shared_col
  ) +
  
  # Whale-present dots
  geom_point(
    data = subset(
      membership_df,
      present & fall_type == "Whale falls"
    ),
    shape = 21,
    size = 5,
    stroke = 0.8,
    fill = whale_col,
    colour = whale_col
  ) +
  
  # Wood-present dots
  geom_point(
    data = subset(
      membership_df,
      present & fall_type == "Wood falls"
    ),
    shape = 21,
    size = 5,
    stroke = 0.8,
    fill = wood_col,
    colour = wood_col
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
      size = 11,
      margin = margin(t = 8)
    ),
    
    axis.text.y = element_text(
      face = "bold",
      size = 11
    ),
    
    legend.position = "none",
    
    plot.margin = margin(0, 15, 5, 5)
  )

# ------------------------------------------------------------
# 7. Set-size data
# ------------------------------------------------------------

set_df <- data.frame(
  fall_type = factor(
    c("Whale falls", "Wood falls"),
    levels = c("Whale falls", "Wood falls")
  ),
  n = c(
    whale_total,
    wood_total
  )
)

# ------------------------------------------------------------
# 8. Set-size bars
# ------------------------------------------------------------

p_sets <- ggplot(
  set_df,
  aes(
    x = n,
    y = fall_type,
    fill = fall_type
  )
) +
  geom_col(
    width = 0.55
  ) +
  geom_text(
    aes(label = n),
    hjust = -0.25,
    fontface = "bold",
    size = 4.5
  ) +
  scale_fill_manual(
    values = c(
      "Whale falls" = whale_col,
      "Wood falls" = wood_col
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1000),
    breaks = c(0, 250, 500, 750),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_discrete(
    limits = rev(levels(set_df$fall_type))
  ) +
  labs(
    x = "Set size",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    
    axis.title.x = element_text(
      size = 12,
      margin = margin(t = 8)
    ),
    
    axis.text.x = element_text(size = 10),
    
    legend.position = "none",
    
    plot.margin = margin(0, 5, 5, 5)
  )

# ------------------------------------------------------------
# 9. Assemble right-hand UpSet panel
# ------------------------------------------------------------

right_panel <- (
  p_intersections /
    p_matrix
) +
  plot_layout(
    heights = c(3.3, 1)
  )

# ------------------------------------------------------------
# 10. Assemble complete figure
# ------------------------------------------------------------

p_upset <- (
  p_sets | right_panel
) +
  plot_layout(
    widths = c(1.05, 2.8)
  ) +
  plot_annotation(
    subtitle = "1,499 Species/OTUs in total",
    theme = theme(
      plot.subtitle = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 8)
      )
    )
  )

# Display
p_upset

# ------------------------------------------------------------
# 11. Export
# ------------------------------------------------------------

ggsave(
  filename = "Outputs/Figures/Figure_Intro_whale_wood_upset.png",
  plot = p_upset,
  width = 9,
  height = 6,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "Outputs/Figures/Figure_Intro_whale_wood_upset.pdf",
  plot = p_upset,
  width = 9,
  height = 6,
  bg = "white"
)

cat(
  "Figure exported to:\n",
  "Outputs/Figures/Figure_Intro_whale_wood_upset.png\n",
  "Outputs/Figures/Figure_Intro_whale_wood_upset.pdf\n"
)