# ============================================================================
# 10_network_module_basin_summary.R
# Descriptive geographic summary of EXISTING observed network modules.
# Does not rebuild networks, rerun Louvain, or run null permutations.
# ============================================================================

source("R/00_setup_and_validate.R")

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(ggplot2); library(scales)
})

dir.create("Outputs/Tables", recursive=TRUE, showWarnings=FALSE)
dir.create("Outputs/Figures", recursive=TRUE, showWarnings=FALSE)

MODULE_FILE <- "Outputs/Tables/network_observed_module_membership.csv"
if(!file.exists(MODULE_FILE)) stop("Could not find: ", MODULE_FILE)

modules <- read_csv(MODULE_FILE, show_col_types=FALSE)
required <- c("Resolution","Node","Fall_Type","Module","Degree","Strength")
missing <- setdiff(required,names(modules))
if(length(missing)) stop("Module file missing: ",paste(missing,collapse=", "))

# Build stable node key from the final site universe.
if("Node_ID" %in% names(site_universe)){
  site_meta <- site_universe %>% mutate(Node=as.character(Node_ID))
} else if(all(c("Site_ID","Fall_Type") %in% names(site_universe))){
  site_meta <- site_universe %>%
    mutate(Node=paste(as.character(Site_ID),as.character(Fall_Type),sep="__"))
} else stop("Could not construct Node from site_universe.")

# Accept either canonical basin field name.
if("Basin" %in% names(site_meta)){
  site_meta <- site_meta %>% mutate(Basin_Final=as.character(Basin))
} else if("Assigned_Basin" %in% names(site_meta)){
  site_meta <- site_meta %>% mutate(Basin_Final=as.character(Assigned_Basin))
} else stop("Could not identify Basin or Assigned_Basin in site_universe.")

site_meta <- site_meta %>%
  transmute(Node,Basin=trimws(Basin_Final)) %>%
  mutate(Basin=ifelse(is.na(Basin)|Basin=="","Unresolved",Basin)) %>%
  distinct(Node,.keep_all=TRUE)

module_geo <- modules %>%
  left_join(site_meta,by="Node") %>%
  mutate(Basin=ifelse(is.na(Basin)|Basin=="","Unresolved",Basin))

if(nrow(module_geo)!=nrow(modules)) stop("Row count changed after basin join.")

cat("\nMODULE x BASIN DESCRIPTIVE SUMMARY\n")
cat("Module-membership rows:",nrow(module_geo),"\n")
cat("Basin values:\n"); print(sort(unique(module_geo$Basin)))

module_basin_counts <- module_geo %>%
  count(Resolution,Module,Basin,name="Nodes") %>%
  group_by(Resolution,Module) %>%
  mutate(Module_Size=sum(Nodes),Basin_Proportion=Nodes/Module_Size) %>%
  ungroup() %>%
  arrange(Resolution,Module,desc(Nodes),Basin)

write_csv(module_basin_counts,
          "Outputs/Tables/network_module_basin_counts.csv")

module_summary <- module_basin_counts %>%
  group_by(Resolution,Module) %>%
  summarise(
    Module_Size=first(Module_Size),
    Basins_Represented=n_distinct(Basin),
    Dominant_Basin=Basin[which.max(Nodes)],
    Dominant_Basin_Nodes=max(Nodes),
    Dominant_Basin_Proportion=max(Basin_Proportion),
    Basin_Entropy=-sum(Basin_Proportion*log(Basin_Proportion)),
    .groups="drop"
  ) %>%
  arrange(Resolution,desc(Module_Size),Module)

write_csv(module_summary,
          "Outputs/Tables/network_module_basin_summary.csv")
cat("\nModule summary:\n"); print(module_summary,n=Inf)

# Largest modules only in the display; all modules remain in the CSVs.
figure_limits <- c("Species/OTU"=8,"Genus"=8,"Family"=9)

largest_modules <- module_summary %>%
  group_by(Resolution) %>%
  arrange(desc(Module_Size),Module,.by_group=TRUE) %>%
  mutate(Module_Rank=row_number()) %>%
  ungroup() %>%
  filter(Module_Rank <= figure_limits[Resolution]) %>%
  select(Resolution,Module,Module_Size,Module_Rank)

plot_data <- module_basin_counts %>%
  inner_join(largest_modules,
             by=c("Resolution","Module","Module_Size")) %>%
  mutate(Module_Label=paste0("M",Module," (n=",Module_Size,")"))

# Explicit ordering key avoids factor-level collisions across facets.
plot_data <- plot_data %>%
  mutate(Module_Key=paste(Resolution,Module,sep="__")) %>%
  left_join(
    largest_modules %>%
      mutate(Module_Key=paste(Resolution,Module,sep="__")) %>%
      select(Module_Key,Module_Rank),
    by="Module_Key",
    suffix=c("","_order")
  )

# Reorder labels globally by rank while facets remain free_y.
label_order <- plot_data %>%
  distinct(Resolution,Module_Label,Module_Rank) %>%
  arrange(Resolution,desc(Module_Rank)) %>%
  pull(Module_Label) %>% unique()

plot_data$Module_Label <- factor(
  plot_data$Module_Label,
  levels = rev(label_order)
)
plot_data <- plot_data %>%
  mutate(
    Resolution = factor(
      Resolution,
      levels = c(
        "Species/OTU",
        "Genus",
        "Family"
      )
    )
  )

p_module_basin <- ggplot(
  plot_data,
  aes(x=Basin_Proportion,y=Module_Label,fill=Basin)
) +
  geom_col(width=0.78) +
  facet_wrap(~Resolution,scales="free_y",ncol=1) +
  scale_x_continuous(labels=percent_format(accuracy=1),
                     limits=c(0,1),expand=c(0,0)) +
  labs(
    x="Proportion of module nodes",
    y=NULL,
    fill="Ocean basin",
    title="Ocean-basin composition of observed community-similarity modules",
    subtitle=paste(
      "Largest observed modules shown;",
      "module membership was derived from community similarity, not geography"
    )
  ) +
  theme_classic(base_size=11) +
  theme(
    plot.title=element_text(face="bold",size=13),
    plot.subtitle=element_text(size=9.5,margin=margin(b=8)),
    strip.text=element_text(face="bold",size=11),
    legend.position="bottom",
    legend.title=element_text(face="bold"),
    panel.spacing=grid::unit(0.8,"lines")
  )

ggsave("Outputs/Figures/Figure_module_basin_summary.png",
       p_module_basin,width=9,height=10,dpi=450,bg="white")
ggsave("Outputs/Figures/Figure_module_basin_summary.pdf",
       p_module_basin,width=9,height=10,bg="white")

cat(
  "\n============================================================\n",
  "MODULE x BASIN SUMMARY COMPLETED SUCCESSFULLY\n",
  "============================================================\n",
  "Created:\n",
  "  Outputs/Tables/network_module_basin_counts.csv\n",
  "  Outputs/Tables/network_module_basin_summary.csv\n",
  "  Outputs/Figures/Figure_module_basin_summary.png\n",
  "  Outputs/Figures/Figure_module_basin_summary.pdf\n",
  "No networks, Louvain modules or null models were rerun.\n",
  sep=""
)
