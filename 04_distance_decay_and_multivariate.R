# 04_distance_decay_and_multivariate.R
# Chapter 2: geographic distance decay and optional multivariate tests.
source("R/00_setup_and_validate.R")

pairwise_distance_similarity <- function(x, resolution) {
  m <- as.matrix(x[, -(1:5)])
  
  longitude <- suppressWarnings(as.numeric(x$Longitude))
  latitude  <- suppressWarnings(as.numeric(x$Latitude))
  
  keep <- rowSums(m) > 0 &
    is.finite(longitude) &
    is.finite(latitude)
  
  m <- m[keep, , drop = FALSE]
  meta <- x[keep, 1:5] %>%
    mutate(
      Longitude = longitude[keep],
      Latitude = latitude[keep]
    )
  
  dis <- as.matrix(
    vegan::vegdist(m, method = "jaccard", binary = TRUE)
  )
  
  geo <- geosphere::distm(
    as.matrix(meta[, c("Longitude", "Latitude")]),
    fun = geosphere::distHaversine
  ) / 1000
  
  idx <- which(upper.tri(dis), arr.ind = TRUE)
  
  tibble(
    Resolution = resolution,
    Node_1 = paste(meta$Site_ID[idx[,1]], meta$Fall_Type[idx[,1]], sep = "__"),
    Node_2 = paste(meta$Site_ID[idx[,2]], meta$Fall_Type[idx[,2]], sep = "__"),
    Distance_km = geo[idx],
    Jaccard_similarity = 1 - dis[idx]
  )
}
dd <- bind_rows(pairwise_distance_similarity(species_matrix,"Species_OTU"),
                pairwise_distance_similarity(genus_matrix,"Genus"),
                pairwise_distance_similarity(family_matrix,"Family"))
write_csv(dd,"Outputs/Tables/distance_decay_pairs.csv")

# Pairwise observations are not statistically independent because every site
# occurs in many pairs. The smooth below is therefore descriptive. Do not use a
# simple ordinary least-squares p-value as if all pairwise rows were independent.
# Figure 5: geographic distance decay across taxonomic resolutions.
# Pairwise observations are not statistically independent because every site
# occurs in many pairs. GAM smoothers are therefore descriptive and should not
# be interpreted as independent inferential tests or mechanistic dispersal models.

dd_plot <- dd %>%
  mutate(
    Resolution = factor(
      Resolution,
      levels = c("Species_OTU", "Genus", "Family"),
      labels = c("Species/OTU", "Genus", "Family")
    )
  )

p5 <- ggplot(
  dd_plot,
  aes(
    x = Distance_km,
    y = Jaccard_similarity
  )
) +
  geom_point(
    alpha = 0.04,
    size = 0.4,
    colour = "grey25"
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 5),
    se = TRUE,
    colour = "#2C6BED",
    fill = "grey75",
    linewidth = 0.9
  ) +
  facet_wrap(
    ~Resolution,
    nrow = 1
  ) +
  scale_x_continuous(
    breaks = seq(0, 20000, by = 5000),
    labels = scales::label_number(big.mark = ",")
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25)
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(
      fill = "grey90",
      colour = "grey40"
    ),
    strip.text = element_text(size = 11),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Great-circle distance (km)",
    y = "Pairwise Jaccard similarity"
  )

ggsave(
  "Outputs/Figures/Figure_5_distance_decay.png",
  p5,
  width = 10,
  height = 4.5,
  dpi = 450
)

ggsave(
  "Outputs/Figures/Figure_5_distance_decay.pdf",
  p5,
  width = 10,
  height = 4.5
)

# PERMANOVA template. Run separately by resolution only after checking that the
# grouping variable has adequate replication. adonis2 is paired with betadisper
# because PERMANOVA can respond to differences in multivariate dispersion.
run_permanova <- function(x) {
  m <- as.matrix(x[, -(1:5)]); keep <- rowSums(m)>0
  m <- m[keep,,drop=FALSE]; meta <- x[keep,1:5]
  d <- vegan::vegdist(m,method="jaccard",binary=TRUE)
  perm <- vegan::adonis2(d ~ Fall_Type, data=meta, permutations=9999)
  bd <- vegan::betadisper(d, meta$Fall_Type)
  list(permanova=perm, dispersion=anova(bd), dispersion_permutation=vegan::permutest(bd,permutations=9999))
}
# Results are intentionally not auto-written here because they must be reviewed
# for group sizes and dispersion before being quoted in the thesis.
