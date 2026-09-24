# 08_export_gephi_files.R
# Export Gephi-ready site and taxon nodes plus two complementary edge types.
source("R/00_setup_and_validate.R")

export_gephi <- function(long_data,resolution,scope="Combined") {
  x <- long_data %>% distinct(Site_ID,Fall_Type,Taxon) %>% mutate(Node_ID=paste(Site_ID,Fall_Type,sep="__"))
  if (scope!="Combined") x <- filter(x,Fall_Type==scope)
  meta <- site_universe %>% filter(Node_ID %in% x$Node_ID) %>%
    transmute(Id=Node_ID,Label=coalesce(Site_Name,Site_ID),Node_Type="Site",Site_ID,Fall_Type,
              Latitude,Longitude,Basin,Region,Taxonomic_Resolution=resolution)
  taxa <- x %>% count(Taxon,name="Sites_Occupied") %>%
    transmute(Id=paste0("taxon__",resolution,"__",Taxon),Label=Taxon,Node_Type="Taxon",
              Site_ID=NA_character_,Fall_Type=NA_character_,Latitude=NA_real_,Longitude=NA_real_,
              Basin=NA_character_,Region=NA_character_,Taxonomic_Resolution=resolution,Sites_Occupied)
  nodes <- bind_rows(meta,taxa)
  bip <- x %>% transmute(Source=Node_ID,Target=paste0("taxon__",resolution,"__",Taxon),Type="Undirected",Weight=1,Taxon)
  dir.create("Gephi/nodes",recursive=TRUE,showWarnings=FALSE); dir.create("Gephi/edges",recursive=TRUE,showWarnings=FALSE)
  write_csv(nodes,file.path("Gephi/nodes",paste0(resolution,"_",scope,"_bipartite_nodes.csv")))
  write_csv(bip,file.path("Gephi/edges",paste0(resolution,"_",scope,"_bipartite_edges.csv")))
}
for (res in c("Species_OTU","Genus","Family")) {
  long <- if(res=="Species_OTU") species_long else if(res=="Genus") genus_long else family_long
  for(scope in c("Whale","Wood","Combined")) export_gephi(long,res,scope)
}
message("Gephi bipartite files exported. See Documentation/GEphi_step_by_step.md.")
