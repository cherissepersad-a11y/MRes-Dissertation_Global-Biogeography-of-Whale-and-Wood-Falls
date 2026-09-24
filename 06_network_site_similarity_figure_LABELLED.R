# ============================================================================
# 06_network_site_similarity_figure_LABELLED.R
# Reproducible site-name callouts added to the final observed network figure.
# Labels alter presentation only; analytical networks/results remain unchanged.
# ============================================================================

source("R/00_setup_and_validate.R")

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(igraph); library(ggplot2)
  library(ggraph); library(ggrepel); library(patchwork); library(scales)
})

dir.create("Outputs/Figures", recursive=TRUE, showWarnings=FALSE)
dir.create("Outputs/Tables", recursive=TRUE, showWarnings=FALSE)

edge_files <- c(
  Species_OTU="Outputs/Tables/Species_OTU_Combined_site_similarity_edges.csv",
  Genus="Outputs/Tables/Genus_Combined_site_similarity_edges.csv",
  Family="Outputs/Tables/Family_Combined_site_similarity_edges.csv"
)
obs <- read_csv("Outputs/Tables/network_observed_summary.csv", show_col_types=FALSE)
seeds <- c(Species_OTU=2026092203, Genus=2026092202, Family=2026092201)
pretty_resolution <- c(Species_OTU="Species/OTU", Genus="Genus", Family="Family")
label_limits <- c(Species_OTU=6, Genus=6, Family=5)

build_observed_network <- function(resolution) {
  edges <- read_csv(edge_files[[resolution]], show_col_types=FALSE) %>%
    filter(is.finite(Weight), Weight > 0)
  g <- graph_from_data_frame(edges, directed=FALSE)
  set.seed(seeds[[resolution]])
  cl <- cluster_louvain(g, weights=E(g)$Weight)
  V(g)$Module <- membership(cl)
  V(g)$Fall_Type <- ifelse(grepl("__Whale$",V(g)$name),"Whale",
                           ifelse(grepl("__Wood$",V(g)$name),"Wood","Unknown"))
  expected <- obs %>% filter(Resolution==resolution)
  if(nrow(expected)!=1) stop("Could not identify observed-summary row for ",resolution)
  got_nodes<-vcount(g); got_edges<-ecount(g); got_components<-components(g)$no
  got_modules<-length(unique(membership(cl))); got_q<-modularity(cl)
  cat("\n",pretty_resolution[[resolution]]," validation:\n",
      "  Nodes:      ",got_nodes,"\n","  Edges:      ",got_edges,"\n",
      "  Components: ",got_components,"\n","  Modules:    ",got_modules,"\n",
      "  Modularity: ",format(got_q,digits=10),"\n",sep="")
  if(got_nodes!=expected$Nodes) stop("Node validation failed for ",resolution)
  if(got_edges!=expected$Edges) stop("Edge validation failed for ",resolution)
  if(got_components!=expected$Components) stop("Component validation failed for ",resolution)
  if(got_modules!=expected$Modules) stop("Module validation failed for ",resolution)
  if(!isTRUE(all.equal(got_q,expected$Modularity,tolerance=1e-7)))
    stop("Modularity validation failed for ",resolution)
  list(graph=g,community=cl,q=got_q,modules=got_modules)
}

networks <- lapply(names(edge_files),build_observed_network)
names(networks) <- names(edge_files)
cat("\nAll observed networks reproduced successfully.\n")

module_table <- bind_rows(lapply(names(networks),function(resolution){
  g<-networks[[resolution]]$graph
  tibble(Resolution=pretty_resolution[[resolution]],Node=V(g)$name,
         Fall_Type=V(g)$Fall_Type,Module=as.integer(V(g)$Module),
         Degree=degree(g),Strength=strength(g,weights=E(g)$Weight))
}))
write_csv(module_table,"Outputs/Tables/network_observed_module_membership.csv")

# ---------------------------------------------------------------------------
# Site-name lookup
# ---------------------------------------------------------------------------
if(!"Site_Name" %in% names(site_universe))
  stop("site_universe does not contain Site_Name.")

if("Node_ID" %in% names(site_universe)){
  site_label_lookup <- site_universe %>%
    transmute(Node=as.character(Node_ID),Site_Name=as.character(Site_Name))
} else if(all(c("Site_ID","Fall_Type") %in% names(site_universe))){
  site_label_lookup <- site_universe %>%
    transmute(Node=paste(as.character(Site_ID),as.character(Fall_Type),sep="__"),
              Site_Name=as.character(Site_Name))
} else stop("site_universe must contain Node_ID, or Site_ID + Fall_Type.")

site_label_lookup <- site_label_lookup %>%
  mutate(Site_Name=trimws(Site_Name)) %>%
  filter(!is.na(Node),Node!="",!is.na(Site_Name),Site_Name!="") %>%
  distinct(Node,.keep_all=TRUE)

# ---------------------------------------------------------------------------
# Reproducible label selection
# Strongest node represents each module; Degree breaks ties.
# Representatives from the largest modules receive priority.
# Both fall types are retained where available.
# ---------------------------------------------------------------------------
select_network_labels <- function(resolution,n_labels){
  g<-networks[[resolution]]$graph
  nm<-tibble(
    Node=V(g)$name,Fall_Type=V(g)$Fall_Type,Module=as.integer(V(g)$Module),
    Degree=as.numeric(degree(g)),
    Strength=as.numeric(strength(g,weights=E(g)$Weight))
  ) %>% left_join(site_label_lookup,by="Node")

  sizes<-nm %>% count(Module,name="Module_Size")
  nm<-nm %>% left_join(sizes,by="Module")

  reps<-nm %>%
    filter(!is.na(Site_Name),Site_Name!="") %>%
    group_by(Module) %>%
    arrange(desc(Strength),desc(Degree),Node,.by_group=TRUE) %>%
    slice(1) %>% ungroup() %>%
    arrange(desc(Module_Size),desc(Strength),desc(Degree),Node)

  selected<-reps %>% slice_head(n=n_labels)

  available<-sort(unique(nm$Fall_Type[nm$Fall_Type %in% c("Whale","Wood")]))
  missing<-setdiff(available,unique(selected$Fall_Type))

  for(ft in missing){
    replacement<-nm %>%
      filter(Fall_Type==ft,!is.na(Site_Name),Site_Name!="",
             !Node %in% selected$Node) %>%
      arrange(desc(Strength),desc(Degree),desc(Module_Size),Node) %>%
      slice_head(n=1)
    if(nrow(replacement)==1 && nrow(selected)>0){
      removable<-selected %>%
        group_by(Fall_Type) %>% mutate(Type_Count=n()) %>% ungroup() %>%
        filter(Type_Count>1 | !Fall_Type %in% c("Whale","Wood")) %>%
        arrange(Module_Size,Strength,Degree,Node)
      if(nrow(removable)>0){
        selected<-selected %>% filter(Node!=removable$Node[1]) %>%
          bind_rows(replacement)
      }
    }
  }

  # Avoid repeating the same site name within a panel.
  selected<-selected %>%
    arrange(desc(Module_Size),desc(Strength),desc(Degree),Node) %>%
    distinct(Site_Name,.keep_all=TRUE)

  if(nrow(selected)<n_labels){
    fill<-reps %>%
      filter(!Node %in% selected$Node,!Site_Name %in% selected$Site_Name)
    selected<-bind_rows(selected,fill %>% slice_head(n=n_labels-nrow(selected)))
  }

  selected %>%
    arrange(desc(Module_Size),desc(Strength),desc(Degree),Node) %>%
    mutate(Resolution=resolution,Label=Site_Name) %>%
    select(Resolution,Node,Site_Name,Label,Fall_Type,Module,Module_Size,Degree,Strength)
}

network_label_tables<-lapply(names(networks),function(resolution)
  select_network_labels(resolution,label_limits[[resolution]]))
names(network_label_tables)<-names(networks)
network_labels_output<-bind_rows(network_label_tables)
write_csv(network_labels_output,"Outputs/Tables/network_figure_labelled_sites.csv")

cat("\nReproducibly selected site labels:\n")
print(network_labels_output %>%
        select(Resolution,Site_Name,Fall_Type,Module,Module_Size,Degree,Strength),
      n=Inf)

# ---------------------------------------------------------------------------
# Plotting function
# ---------------------------------------------------------------------------
make_network_panel <- function(resolution,panel_letter){
  g<-networks[[resolution]]$graph
  set.seed(seeds[[resolution]])
  lay<-create_layout(g,layout="fr",weights=E(g)$Weight,niter=1500)
  lay$Module_plot<-factor(lay$Module)

  labels_this_panel<-network_label_tables[[resolution]] %>% select(Node,Label)
  lay<-lay %>% left_join(labels_this_panel,by=c("name"="Node"))

  q_value<-networks[[resolution]]$q
  n_modules<-networks[[resolution]]$modules

  ggraph(lay)+
    geom_edge_link(aes(edge_alpha=Weight,edge_width=Weight),
                   colour="grey55",show.legend=FALSE)+
    scale_edge_alpha(range=c(0.03,0.35))+
    scale_edge_width(range=c(0.05,0.65))+
    geom_node_point(aes(colour=Module_plot,shape=Fall_Type),
                    size=2.15,alpha=0.92,stroke=0.25)+
    ggrepel::geom_text_repel(
      aes(x=x, 
          y=y, 
          label=Label),
      na.rm=TRUE,size=2.4,
      box.padding=0.45,point.padding=0.25,min.segment.length=0,
      max.overlaps=Inf,seed=seeds[[resolution]],
      segment.size=0.25,segment.alpha=0.65,show.legend=FALSE
    )+
    scale_shape_manual(values=c(Whale=16,Wood=17,Unknown=15),name="Fall type")+
    guides(colour="none",
           shape=guide_legend(override.aes=list(size=4,alpha=1)))+
    labs(title=paste0(panel_letter,"  ",pretty_resolution[[resolution]]),
         subtitle=paste0("Q = ",sprintf("%.3f",q_value),
                         "   |   ",n_modules," modules"))+
    theme_void(base_size=11)+
    theme(plot.title=element_text(face="bold",size=12,hjust=0),
          plot.subtitle=element_text(size=9.5,hjust=0,margin=margin(b=4)),
          legend.position="bottom",
          legend.title=element_text(face="bold"),
          plot.margin=margin(12,16,12,16))
}

p_species<-make_network_panel("Species_OTU","A")
p_genus<-make_network_panel("Genus","B")
p_family<-make_network_panel("Family","C")

figure6_labelled<-(
  p_species | p_genus | p_family
)+
  plot_layout(guides="collect",widths=c(1,1,1))+
  plot_annotation(
    title="Observed site-similarity network structure across taxonomic resolutions",
    theme=theme(plot.title=element_text(face="bold",size=14,hjust=0.5))
  )&
  theme(legend.position="bottom")

# New filenames preserve the existing frozen unlabelled figure.
ggsave("Outputs/Figures/Figure_6_site_similarity_network_LABELLED.png",
       figure6_labelled,width=15,height=6.3,dpi=450,bg="white")
ggsave("Outputs/Figures/Figure_6_site_similarity_network_LABELLED.pdf",
       figure6_labelled,width=15,height=6.3,bg="white")

cat(
  "\n============================================================\n",
  "LABELLED NETWORK FIGURE COMPLETED SUCCESSFULLY\n",
  "============================================================\n",
  "Outputs/Figures/Figure_6_site_similarity_network_LABELLED.png\n",
  "Outputs/Figures/Figure_6_site_similarity_network_LABELLED.pdf\n",
  "Outputs/Tables/network_figure_labelled_sites.csv\n",
  "No null-model permutations were rerun.\n",
  sep=""
)
