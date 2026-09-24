# 05_network_site_similarity.R
# Site-similarity networks: compositional association among organic-fall nodes.
source("R/00_setup_and_validate.R")

build_site_similarity <- function(long_data, scope="Combined") {
  x <- long_data %>% distinct(Site_ID,Fall_Type,Taxon) %>% mutate(Node_ID=paste(Site_ID,Fall_Type,sep="__"))
  if (scope != "Combined") x <- filter(x,Fall_Type==scope)
  incidence <- xtabs(~ Node_ID + Taxon, data=x) > 0
  d <- as.matrix(vegan::vegdist(incidence,method="jaccard",binary=TRUE))
  s <- 1-d; diag(s) <- 0
  idx <- which(upper.tri(s) & s>0,arr.ind=TRUE)
  tibble(Source=rownames(s)[idx[,1]],Target=colnames(s)[idx[,2]],Weight=s[idx],Jaccard=s[idx])
}

# IMPORTANT: an edge means two analytical nodes share recorded taxa. It is a
# compositional-similarity edge, not evidence that organisms dispersed directly
# between those locations. Threshold sensitivity and null models are therefore
# required before interpreting modules or connector sites biologically.
for (res in c("Species_OTU","Genus","Family")) {
  long <- get(paste0(tolower(sub("_OTU","",res)),"_long"), inherits=TRUE)
  if (res=="Species_OTU") long <- species_long
  for (scope in c("Whale","Wood","Combined")) {
    e <- build_site_similarity(long,scope)
    write_csv(e,file.path("Outputs","Tables",paste0(res,"_",scope,"_site_similarity_edges.csv")))
  }
}
