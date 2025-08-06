##################################################################################################
# Propagation of evidence compound2compound
##################################################################################################
rm(list = ls())
setwd("~/data/Untitled/Trabajo/R_practice/exp_compound")
library(bnlearn)
library(magrittr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(parallel)
########################################################################################################
# 
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/propagationFunctions.R")
load("data/raw_aggr/t2m_ERA5_monthly_1940_2022_10d.rda")
load("data/raw_aggr/tp_ERA5_monthly_1940_2022_10d.rda")
###########################################
# Standardize
###########################################
limit <- TRUE
############################################
# discretize in b
############################################
b <- 3

############################################################################
# Function to load hciterations of models in a list 
############################################################################
#permused <- 0
#algo <- "tabu"
#score <- "aic"

loadIterationsComp <- function(permused = 0, algo = c("hc","tabu"),score = c("aic","bic"), it = NULL, limit = FALSE,compound = NULL) {
  if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
  if(!is.null(compound)){comp <- paste0("_comp_",compound)} else {comp <-""}
  
  hc_list <- list.files(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_10d",lim,comp,"/perm",permused), full.names = T)
  hc_names <- list.files(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_10d",lim,comp,"/perm",permused))
  hc_names <- gsub(".rda", "", hc_names)
  hc_networks <- lapply(hc_list, function(x){get(load(x))})
  names(hc_networks) <- hc_names
  sizes <- sapply(hc_networks,narcs)
  hc_networks <- hc_networks[order(sizes)]
  return(hc_networks)
}

tabu_bic_0 <- loadIterationsComp(permused = 0,algo = "tabu", score = "aic",limit = limit)

art.compound.dag <- function(g){
  ngrids <- length(nodes(g))/2
  comp.nodes <- paste0("C",1:ngrids)
  g.comp <- g
  for(i in 1:length(comp.nodes)){
    g.comp <- add.node(g.comp,comp.nodes[i])
    g.comp <- set.arc(g.comp, nodes(g)[i],comp.nodes[i])
    g.comp <- set.arc(g.comp, nodes(g)[i+ngrids],comp.nodes[i])
  }
  return(g.comp)
}

tabu_bic_0_art <- lapply(tabu_bic_0,art.compound.dag)
names(tabu_bic_0_art)<- gsub("lim_0","lim_comp_art_0",names(tabu_bic_0_art))

tabu_bic_0_art[[1]]
k<-0 
score <- "aic"
for (x in 1:length(tabu_bic_0_art)){
  assign(names(tabu_bic_0_art)[x],tabu_bic_0_art[[x]])
  save(list = names(tabu_bic_0_art)[x],file = paste0("data/tabuiterations/disc_",b,"_",score,"_t2m_tp_10d_lim_comp_art/perm",k,"/",names(tabu_bic_0_art)[x],".rda"))
}
