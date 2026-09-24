# 02_chapter1_figures.R
# Chapter 1 publication figures: global maps, depth-latitude panels and research footprint.
source("R/00_setup_and_validate.R")

world <- ggplot2::map_data("world")
basin_palette <- c("Atlantic"="#D55E00", "Indian"="#7570B3", "Pacific"="#009E73",
                   "Southern"="#0072B2", "Arctic"="#CC79A7", "Unresolved"="#999999")

# Species/OTU richness is calculated at the same Site x Fall_Type analytical
# node used in the community matrices. It is a count of recorded taxa and is
# therefore sensitive to research effort; captions should state this caveat.
richness <- species_long %>% distinct(Site_ID, Fall_Type, Taxon) %>%
  count(Site_ID, Fall_Type, name="Species_OTU_richness") %>%
  left_join(site_universe, by=c("Site_ID","Fall_Type"))

# Event depths are summarised to Site x Fall_Type for the spatial-bathymetric
# figure so the plotted unit matches the map. Multiple events at one node are
# represented by their median midpoint depth rather than pseudoreplicated.
event_depth <- events %>% filter(Analysis_Include==1, Fall_Type %in% c("Whale","Wood")) %>%
  mutate(Depth_m=coalesce((Event_Depth_Min_m+Event_Depth_Max_m)/2,
                          Event_Depth_Min_m, Event_Depth_Max_m)) %>%
  group_by(Site_ID, Fall_Type) %>%
  summarise(Depth_m=median(Depth_m,na.rm=TRUE), .groups="drop") %>%
  left_join(site_universe, by=c("Site_ID","Fall_Type"))

map_panel <- function(ft) {
  ggplot() +
    geom_polygon(
      data = world,
      aes(long, lat, group = group),
      fill = "grey85",
      colour = "grey55",
      linewidth = .15
    ) +
    geom_point(
      data = filter(richness, Fall_Type == ft),
      aes(
        Longitude,
        Latitude,
        size = Species_OTU_richness,
        colour = Basin
      ),
      alpha = .75
    ) +
    scale_colour_manual(
      values = basin_palette,
      drop = FALSE
    ) +
    scale_size_area(max_size = 7) +
    coord_quickmap(
      xlim = c(-180, 180),
      ylim = c(-70, 85),
      expand = FALSE
    ) +
    labs(
      title = paste0(ft, "-fall sampling distribution"),
      x = NULL,
      y = NULL,
      colour = "Ocean basin",
      size = "Species/OTU richness"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank()
    )
}

depth_panel <- function(ft) {
  d <- filter(
    event_depth,
    Fall_Type == ft,
    is.finite(Depth_m)
  )
  
  ggplot(d, aes(Latitude, Depth_m, colour = Basin)) +
    geom_point(alpha = .72, size = 2) +
    
    # LOESS smoother is descriptive only. Its purpose is to show broad
    # depth-latitude structure; it is not an inferential test of a
    # mechanistic latitudinal depth relationship.
    geom_smooth(
      aes(group = 1),
      method = "loess",
      se = TRUE,
      colour = "black",
      linewidth = .8
    ) +
    scale_colour_manual(
      values = basin_palette,
      drop = FALSE
    ) +
    scale_y_reverse() +
    labs(
      title = paste0(ft, "-fall depth–latitude distribution"),
      x = "Latitude",
      y = "Depth (m)",
      colour = "Ocean basin"
    ) +
    theme_bw(base_size = 10)
}

fig1 <- (
  map_panel("Whale") | depth_panel("Whale")
) / (
  map_panel("Wood") | depth_panel("Wood")
) +
  patchwork::plot_annotation(tag_levels = "A")

ggsave(
  "Outputs/Figures/Figure_1_spatial_bathymetric.png",
  fig1,
  width = 13,
  height = 7.25,
  dpi = 450
)

ggsave(
  "Outputs/Figures/Figure_1_spatial_bathymetric.pdf",
  fig1,
  width = 13,
  height = 7.25
)

# Figure 2: research footprint across basin and origin categories.

basin_counts <- site_universe %>%
  count(Fall_Type, Basin) %>%
  mutate(
    Basin = factor(
      Basin,
      levels = c(
        "Unresolved",
        "Mediterranean",
        "Arctic",
        "Southern",
        "Indian",
        "Pacific",
        "Atlantic"
      )
    )
  )

origin_counts <- events %>%
  filter(
    Analysis_Include == 1,
    Fall_Type %in% c("Whale", "Wood")
  ) %>%
  count(Fall_Type, Origin) %>%
  mutate(
    Origin = factor(
      Origin,
      levels = c("Unknown", "Experimental", "Natural")
    )
  )

p1 <- ggplot(
  basin_counts,
  aes(Basin, n, fill = Fall_Type)
) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    title = "Sampling distribution among basins",
    x = NULL,
    y = "Analytical nodes",
    fill = "Fall type"
  )

p2 <- ggplot(
  origin_counts,
  aes(Origin, n, fill = Fall_Type)
) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    title = "Sampling origin",
    x = NULL,
    y = "Included events",
    fill = "Fall type"
  )

fig2 <- (p1 | p2) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(legend.position = "right")

ggsave(
  "Outputs/Figures/Figure_2_research_footprint.png",
  fig2,
  width = 11,
  height = 5.5,
  dpi = 450
)

ggsave(
  "Outputs/Figures/Figure_2_research_footprint.pdf",
  fig2,
  width = 11,
  height = 5.5
)