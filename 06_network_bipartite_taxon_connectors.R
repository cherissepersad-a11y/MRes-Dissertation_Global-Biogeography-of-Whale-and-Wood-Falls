# 06_network_bipartite_taxon_connectors.R
# Taxon-resolved bipartite networks: identify WHICH taxa contribute to links.
source("R/00_setup_and_validate.R")

build_bipartite <- function(long_data, resolution, scope="Combined") {
  x <- long_data %>% distinct(Site_ID,Fall_Type,Taxon) %>%
    mutate(Site_Node=paste(Site_ID,Fall_Type,sep="__"), Taxon_Node=paste0("taxon__",resolution,"__",Taxon))
  if (scope != "Combined") x <- filter(x,Fall_Type==scope)
  g <- igraph::graph_from_data_frame(x %>% transmute(from=Site_Node,to=Taxon_Node),directed=FALSE)
  V(g)$type <- startsWith(V(g)$name,"taxon__")
  list(graph=g,data=x)
}

connector_table <- function(long_data,resolution,scope="Combined") {
  b <- build_bipartite(long_data,resolution,scope); g <- b$graph; x <- b$data
  # Participation across basins is a transparent, data-traceable descriptor of
  # how broadly a taxon spans the sampled geography. It is not gene flow.
  site_meta <- site_universe %>%
    select(Node_ID, Basin, Region)
  x2 <- x %>% 
    left_join(site_meta,by=c("Site_Node"="Node_ID"))
  x2 %>% group_by(Taxon) %>% summarise(
    Sites_Occupied=n_distinct(Site_Node), Basins_Spanned=n_distinct(Basin[!is.na(Basin)]),
    Regions_Spanned=n_distinct(Region[!is.na(Region)]),
    Whale_Sites=n_distinct(Site_Node[Fall_Type=="Whale"]), Wood_Sites=n_distinct(Site_Node[Fall_Type=="Wood"]),
    Fall_Type_Breadth=case_when(Whale_Sites>0 & Wood_Sites>0 ~ "Whale and Wood", Whale_Sites>0 ~ "Whale only", TRUE ~ "Wood only"),
    .groups="drop") %>% arrange(desc(Basins_Spanned),desc(Sites_Occupied)) %>% mutate(Resolution=resolution,Scope=scope)
}

all_connectors <- bind_rows(
  connector_table(species_long,"Species_OTU"), connector_table(genus_long,"Genus"), connector_table(family_long,"Family"))
write_csv(all_connectors,"Outputs/Tables/taxon_connector_summary.csv")

# -------------------------------------------------------------------------
# Figure 6: taxon-resolved bipartite networks
#
# Sites and taxa are represented as the two node classes. Taxon labels are
# restricted using an explicit analytical rule: taxa must occur at >=3
# analytical nodes and span >=2 sampled ocean basins.
#
# These networks represent recorded site-taxon incidence. They are descriptive
# representations of compositional connectivity and are not evidence of direct
# dispersal or gene flow.
# -------------------------------------------------------------------------

plot_bipartite_connectors <- function(
    long_data,
    resolution,
    panel_label,
    scope = "Combined"
) {
  
  b <- build_bipartite(
    long_data,
    resolution,
    scope
  )
  
  g <- b$graph
  
  conn <- connector_table(
    long_data,
    resolution,
    scope
  )
  
  # Connector taxa are defined as spanning >=2 sampled basins and occurring
  # at >=3 analytical nodes. For figure readability, labels are restricted
  # to the 15 highest-ranking connector taxa by basin breadth and site occupancy.
  # All qualifying taxa remain retained in the connector summary table.
  label_taxa <- conn %>%
    filter(
      Basins_Spanned >= 2,
      Sites_Occupied >= 3
    ) %>%
    arrange(
      desc(Basins_Spanned),
      desc(Sites_Occupied),
      Taxon
    ) %>%
    slice_head(n = 15) %>%
    pull(Taxon)
  
  # Recover taxon name directly from the graph node name.
  taxon_prefix <- paste0(
    "taxon__",
    resolution,
    "__"
  )
  
  taxon_name <- ifelse(
    V(g)$type,
    sub(
      paste0("^", taxon_prefix),
      "",
      V(g)$name
    ),
    NA_character_
  )
  
  V(g)$display_label <- ifelse(
    V(g)$type &
      taxon_name %in% label_taxa,
    taxon_name,
    NA_character_
  )
  
  # Fixed seed makes the force-directed layout reproducible.
  set.seed(20260917)
  
  ggraph::ggraph(
    g,
    layout = "fr"
  ) +
    ggraph::geom_edge_link(
      colour = "grey75",
      alpha = 0.10,
      linewidth = 0.20
    ) +
    ggraph::geom_node_point(
      aes(shape = type),
      size = 1.35,
      colour = "black"
    ) +
    ggraph::geom_node_text(
      aes(label = display_label),
      repel = TRUE,
      size = 2.35,
      na.rm = TRUE,
      max.overlaps = Inf
    ) +
    scale_shape_manual(
      values = c(
        `FALSE` = 16,
        `TRUE` = 15
      ),
      labels = c(
        "Site",
        "Taxon"
      ),
      name = "Node type"
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 12,
        face = "bold"
      ),
      legend.position = "right"
    ) +
    labs(
      title = panel_label
    )
}


# Species/OTU -> Genus -> Family ordering is retained across Chapter 2 figures.

p6a <- plot_bipartite_connectors(
  species_long,
  "Species_OTU",
  "Species/OTU"
)

p6b <- plot_bipartite_connectors(
  genus_long,
  "Genus",
  "Genus"
)

p6c <- plot_bipartite_connectors(
  family_long,
  "Family",
  "Family"
)


# Collect the node-type legend once for the entire figure.

fig6 <- (
  p6a | p6b | p6c
) +
  patchwork::plot_layout(
    guides = "collect"
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "right"
  )


ggsave(
  "Outputs/Figures/Figure_6_taxon_resolved_connectivity.png",
  fig6,
  width = 15,
  height = 6,
  dpi = 450
)

ggsave(
  "Outputs/Figures/Figure_6_taxon_resolved_connectivity.pdf",
  fig6,
  width = 15,
  height = 6
)