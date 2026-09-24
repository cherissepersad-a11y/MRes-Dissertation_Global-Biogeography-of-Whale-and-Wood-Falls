# 10_distance_decay_spearman_permutation_check.R
# Standalone verification of distance-decay Spearman correlations using
# node-label permutations. Source data are read only.

pkgs <- c("readr","dplyr","vegan","geosphere")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if(length(miss)) stop("Install missing package(s): ", paste(miss, collapse=", "))
library(readr); library(dplyr); library(vegan); library(geosphere)

N_PERM <- 999
SEEDS <- c(Species_OTU=2026092301, Genus=2026092302, Family=2026092303)
EXPECTED_RHO <- c(Species_OTU=-0.2098, Genus=-0.1141, Family=-0.0610)
TOL <- 0.001

cat("Select Analysis_SiteUniverse.csv\n"); SITE_FILE <- file.choose()
cat("Select Species_OTU_SitePresence_Matrix.csv\n"); SPECIES_FILE <- file.choose()
cat("Select Genus_SitePresence_Matrix.csv\n"); GENUS_FILE <- file.choose()
cat("Select Family_SitePresence_Matrix.csv\n"); FAMILY_FILE <- file.choose()

OUT <- file.path(dirname(SPECIES_FILE),"distance_decay_permutation_verification")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

sites <- read_csv(SITE_FILE, show_col_types=FALSE)
sp <- read_csv(SPECIES_FILE, show_col_types=FALSE)
ge <- read_csv(GENUS_FILE, show_col_types=FALSE)
fa <- read_csv(FAMILY_FILE, show_col_types=FALSE)

req <- c("Site_ID","Fall_Type","Latitude","Longitude")
if(length(setdiff(req,names(sites)))) stop("Wrong SiteUniverse file selected.")

sites <- sites %>% mutate(
 Site_ID=trimws(as.character(Site_ID)),
 Fall_Type=trimws(as.character(Fall_Type)),
 Node_Key=paste(Site_ID,Fall_Type,sep="__"),
 Latitude_final=suppressWarnings(parse_number(as.character(Latitude))),
 Longitude_final=suppressWarnings(parse_number(as.character(Longitude)))
)
if(anyDuplicated(sites$Node_Key)) stop("Duplicate Site x Fall-Type nodes.")
if(nrow(sites)!=863) stop("Expected 863 analytical nodes; found ",nrow(sites))
cat("Analytical-universe check PASSED: 863 nodes.\n")

prep <- function(df,res){
 if(length(setdiff(c("Site_ID","Fall_Type"),names(df)))) stop(res,": wrong matrix selected.")
 df <- df %>% mutate(
   Site_ID=trimws(as.character(Site_ID)),
   Fall_Type=trimws(as.character(Fall_Type)),
   Node_Key=paste(Site_ID,Fall_Type,sep="__")
 )
 if(anyDuplicated(df$Node_Key)) stop(res,": duplicate node rows.")
 meta <- c("Site_ID","Fall_Type","Node_Key","Site_Name","Latitude","Longitude",
           "Ocean","Region","Verbatim_Basin","Assigned_Basin","Depth_m","Node_Depth","Richness")
 taxa <- setdiff(names(df),intersect(meta,names(df)))
 m <- as.matrix(as.data.frame(lapply(df[,taxa,drop=FALSE], function(x)
   suppressWarnings(as.numeric(as.character(x)))),check.names=FALSE))
 rownames(m) <- df$Node_Key
 if(anyNA(m)) stop(res,": NA/non-numeric values in taxon columns.")
 if(!all(unique(as.vector(m)) %in% c(0,1))) stop(res,": matrix is not binary.")
 if(nrow(m)!=863) stop(res,": expected 863 matrix rows; found ",nrow(m))
 if(length(setdiff(rownames(m),sites$Node_Key))) stop(res,": unknown analytical nodes.")
 m <- m[rowSums(m)>0,,drop=FALSE]
 cat(res,": ",nrow(m)," non-empty nodes; ",ncol(m)," taxa.\n",sep="")
 m
}

spm<-prep(sp,"Species_OTU"); gem<-prep(ge,"Genus"); fam<-prep(fa,"Family")
if(!identical(c(nrow(spm),nrow(gem),nrow(fam)),c(818L,286L,115L)))
 stop("Non-empty node counts do not match verified final values.")
if(!identical(c(ncol(spm),ncol(gem),ncol(fam)),c(1499L,473L,147L)))
 stop("Taxon counts do not match verified final values.")

build <- function(m,res){
 md <- sites[match(rownames(m),sites$Node_Key),,drop=FALSE]
 if(anyNA(md$Node_Key) || !identical(rownames(m),md$Node_Key)) stop(res,": node matching failed.")
 ok <- is.finite(md$Latitude_final) & is.finite(md$Longitude_final)
 m <- m[ok,,drop=FALSE]; md <- md[ok,,drop=FALSE]
 sim <- 1-as.matrix(vegdist(m,method="jaccard",binary=TRUE))
 rownames(sim)<-colnames(sim)<-rownames(m)
 xy <- cbind(md$Longitude_final,md$Latitude_final)
 geo <- distm(xy,fun=distHaversine)/1000
 rownames(geo)<-colnames(geo)<-md$Node_Key
 if(!identical(rownames(geo),rownames(sim))) stop(res,": matrix alignment failed.")
 cat(res,": ",nrow(geo)," coordinate-complete nodes; ",choose(nrow(geo),2)," pairs.\n",sep="")
 list(geo=geo,sim=sim,n=nrow(geo),pairs=choose(nrow(geo),2))
}

spo<-build(spm,"Species_OTU"); geo<-build(gem,"Genus"); fao<-build(fam,"Family")

obs <- function(o){
 u<-upper.tri(o$geo)
 gv<-o$geo[u]; sv<-o$sim[u]; keep<-is.finite(gv)&is.finite(sv)
 list(rho=cor(gv[keep],sv[keep],method="spearman"),upper=u,keep=keep,gv=gv[keep])
}
sobs<-obs(spo); gobs<-obs(geo); fobs<-obs(fao)
got<-c(Species_OTU=sobs$rho,Genus=gobs$rho,Family=fobs$rho)
check<-data.frame(Resolution=names(EXPECTED_RHO),Expected_rho=as.numeric(EXPECTED_RHO),
 Recalculated_rho=as.numeric(got[names(EXPECTED_RHO)]))
check$Absolute_difference<-abs(check$Recalculated_rho-check$Expected_rho)
check$Pass<-check$Absolute_difference<=TOL
print(check,row.names=FALSE)
write_csv(check,file.path(OUT,"distance_decay_rho_reproduction_check.csv"))
if(!all(check$Pass)) stop("Rho reproduction FAILED. Permutations not run.")
cat("Observed-rho reproduction PASSED.\n")

permtest <- function(o,ob,res,seed){
 set.seed(seed); null<-rep(NA_real_,N_PERM)
 for(i in seq_len(N_PERM)){
   p<-sample.int(o$n)
   pv<-o$sim[p,p,drop=FALSE][ob$upper][ob$keep]
   null[i]<-cor(ob$gv,pv,method="spearman")
   if(i%%100==0 || i==N_PERM) cat(res,": ",i,"/",N_PERM,"\n",sep="")
 }
 v<-null[is.finite(null)]; nv<-length(v)
 pl<-(sum(v<=ob$rho)+1)/(nv+1); pu<-(sum(v>=ob$rho)+1)/(nv+1)
 pt<-min(1,2*min(pl,pu))
 list(
  summary=data.frame(Resolution=res,Nodes=o$n,Pairs=o$pairs,Spearman_rho=ob$rho,
    Null_mean_rho=mean(v),Null_SD_rho=sd(v),
    Null_q025=as.numeric(quantile(v,.025)),Null_q975=as.numeric(quantile(v,.975)),
    P_lower=pl,P_upper=pu,P_two_sided=pt,Permutations_requested=N_PERM,
    Permutations_valid=nv,Seed=seed),
  null=data.frame(Resolution=res,Permutation=seq_along(null),Spearman_rho=null)
 )
}

cat("Running Species/OTU permutations...\n")
ps<-permtest(spo,sobs,"Species_OTU",SEEDS[["Species_OTU"]])
cat("Running Genus permutations...\n")
pg<-permtest(geo,gobs,"Genus",SEEDS[["Genus"]])
cat("Running Family permutations...\n")
pf<-permtest(fao,fobs,"Family",SEEDS[["Family"]])

summary_all<-bind_rows(ps$summary,pg$summary,pf$summary)
null_all<-bind_rows(ps$null,pg$null,pf$null)
write_csv(summary_all,file.path(OUT,"distance_decay_spearman_permutation_summary.csv"))
write_csv(null_all,file.path(OUT,"distance_decay_spearman_null_999.csv"))

cat("\nFINAL RESULTS\n"); print(summary_all,row.names=FALSE)
cat("\n*** SCRIPT COMPLETED SUCCESSFULLY ***\n")
cat("Outputs saved to: ",OUT,"\n",sep="")
