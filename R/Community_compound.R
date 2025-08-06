########################################################
rm(list = ls())
setwd("~/data/Untitled/Trabajo/R_practice/exp_compound")
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("data/raw_aggr/t2m_ERA5_monthly_1940_2022_10d.rda")
load("data/raw_aggr/tp_ERA5_monthly_1940_2022_10d.rda")
###########################################
# Standardize
###########################################
limit <- TRUE
if(isTRUE(limit)){
  sel.tp <- subsetGrid(tp_ERA5_monthly_1940_2022_10d,latLim = c(-50,75))
  sel.t2m <- subsetGrid(t2m_ERA5_monthly_1940_2022_10d,latLim = c(-50,75))
} else {
  sel.tp <- tp_ERA5_monthly_1940_2022_10d
  sel.t2m <- t2m_ERA5_monthly_1940_2022_10d
}
tp_ERA5_monthly_1940_2022_10d <- NULL
t2m_ERA5_monthly_1940_2022_10d <- NULL
df.tp <-as.data.frame(TimeCoordsAnom_from_Grid_rms(sel.tp,rms = TRUE))
df.t2m <-as.data.frame(TimeCoordsAnom_from_Grid_rms(sel.t2m,rms = TRUE))
#sel.tp<- NULL
#sel.t2m <- NULL
############################################
# discretize in b
############################################
b <- 3
disc.df.tp <- discretize(df.tp, method = "quantile", breaks = b)
disc.df.t2m <- discretize(df.t2m, method = "quantile", breaks = b)
names(disc.df.t2m) <- paste0(names(disc.df.t2m),".t2m")
names(disc.df.tp) <- paste0(names(disc.df.tp),".tp")
disc.df.t2m.tp <- cbind(disc.df.t2m,disc.df.tp)
df.t2m<- NULL
df.tp <- NULL
disc.df.t2m <- NULL
disc.df.tp <- NULL

compound <- "HD_only"
if(!is.null(compound)){
  ngrids <- ncol(disc.df.t2m.tp)/2
  comp.nodes <- paste0("C",1:ngrids)
  
  compound.data <- function(data,comp.type,node){ 
    data.num <- as.data.frame(sapply(data, function(x) as.numeric(x), simplify = FALSE))
    c<- node
    data[[comp.nodes[c]]]<- numeric(nrow(data))
    
    if (comp.type == "all.distinct"){
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==1 & data.num[[paste0("V",c,".tp")]]==3] <- 13
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==1 & data.num[[paste0("V",c,".tp")]]==1] <- 11
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==3 & data.num[[paste0("V",c,".tp")]]==1] <- 31
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==3 & data.num[[paste0("V",c,".tp")]]==3] <- 33
    } else if (comp.type == "all"){
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==1 & data.num[[paste0("V",c,".tp")]]==3] <- 1
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==1 & data.num[[paste0("V",c,".tp")]]==1] <- 1
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==3 & data.num[[paste0("V",c,".tp")]]==1] <- 1
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==3 & data.num[[paste0("V",c,".tp")]]==3] <- 1
    } else if (comp.type == "HD"){   
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==3 & data.num[[paste0("V",c,".tp")]]==1] <- 1
    } else if (comp.type == "HW"){
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==3 & data.num[[paste0("V",c,".tp")]]==3] <- 1
    } else if (comp.type == "CW"){
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==1 & data.num[[paste0("V",c,".tp")]]==3] <- 1
    } else if (comp.type == "CD"){
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==1 & data.num[[paste0("V",c,".tp")]]==1] <- 1
    }
    data[[comp.nodes[c]]] <- as.factor(data[[comp.nodes[c]]])
    return(data)
  }
  
  
  for(i in 1:(ncol(disc.df.t2m.tp)/2)){
    disc.df.t2m.tp <- compound.data(disc.df.t2m.tp,gsub("_only","",compound),i)
  }
  
  #   ,function(x)compound.data(disc.df.t2m.tp,compound,x))
  # test[sapply(test, is.numeric)] <- lapply(test[,sapply(test, is.numeric)], as.factor)
  # 
  # colnames(test)<- comp.nodes
  # disc.df.t2m.tp <- data.frame(disc.df.t2m.tp,as.data.frame(test))
  # 
}
# 
# test[sapply(test, is.numeric)] <- lapply(test[sapply(test, is.numeric)], as.factor)
# str(mydata)
# 
# 


#####################################################################################
# create directories
#####################################################################################
#k <- 0
#algo<- "tabu"
#score<- "bic"
if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
if(!is.null(compound)){comp <- paste0("_comp_",compound)} else {comp <-""}

#  dir.create(paste0("data/",algo,"iterations/"))
# dir.create(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_10d",lim,comp,"/"))
# dir.create(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_10d",lim,comp,"/perm",k))
# dir.create(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_10d",lim,comp,"/perm",k,"train"))
#  dir.create(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_10d",lim,comp,"/perm",k,"test"))
###########################################################################
# 
###########################################################################
if (!is.null(grep("only",compound))){
  data <- disc.df.t2m.tp[,grep("C",colnames(disc.df.t2m.tp))]
} else {data <- disc.df.t2m.tp}
# data.num <- as.data.frame(sapply(data, function(x) as.numeric(x), simplify = FALSE))
#data <- df[permutations[[k]]]
############################################################################
# Function to load hciterations of models in a list 
############################################################################
#permused <- 0
#algo <- "tabu"
#score <- "aic"

loadIterationsComp <- function(permused = 0, algo = c("hc","tabu"),score = c("aic","bic"), it = NULL, limit = FALSE,compound = NULL,artificial = FALSE) {
  if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
  if(compound =="comp"){comp<- "_comp"} else if(!is.null(compound)){comp <- paste0("_comp_",compound)} else {comp <-""}
  if(isTRUE(artificial)){art <- paste0("_art")} else {art <-""}
  hc_list <- list.files(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_10d",lim,comp,art,"/perm",permused), full.names = T)
  hc_names <- list.files(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_10d",lim,comp,art,"/perm",permused))
  hc_names <- gsub(".rda", "", hc_names)
  
  if(!is.null(it)){
    hc_list <- hc_list[grep(it,hc_list)]
    hc_names <- hc_names[grep(it,hc_names)]
  }
  
  hc_networks <- lapply(hc_list, function(x){get(load(x))})
  names(hc_networks) <- hc_names
  sizes <- sapply(hc_networks,narcs)
  hc_networks <- hc_networks[order(sizes)]
  
  
  hc_networks <- lapply(hc_list, function(x){get(load(x))})
  names(hc_networks) <- hc_names
  sizes <- sapply(hc_networks,narcs)
  hc_networks <- hc_networks[order(sizes)]
  return(hc_networks)
}

Algo <- "tabu"
Score <- "bic"
IT<- "0_1000_1100"
tabu_bic_0_comp_only <- loadIterationsComp(
  permused = 0,algo = Algo, score = Score,
  limit = limit,compound = "HD_only",artificial = FALSE, it = IT)

tabu_aic_0_comp_only_fit<- lapply(tabu_bic_0_comp_only, function(x) bn.fit(x = x, data = data))
#############################################################
#
#############################################################
igraph <- igraph.from.graphNEL(as.graphNEL(tabu_aic_0_comp_only_fit[[1]]))
HD_clusters_bic <- cluster_edge_betweenness(igraph,directed = FALSE)

ncom <- length(levels(factor(HD_clusters_bic$membership)))
mem <- membership(HD_clusters_bic)
memClim <- quantity2clim(mem, "membership", sel.t2m, backperm = NULL)
if(ncom <= 9){
  colRainbow <- brewer.pal(ncom,"Set1")
} else {
  colRainbow <- brewer.pal(9,"Set1")
  colRainbow<- colorRampPalette(colRainbow)(ncom)
}
if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
plotclim <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                                    set.min = NULL, set.max = NULL, 
                                    lonCenter = 0, 
                                    regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                                   # main = paste0("Louvain BN ",whichsizesBN[i]," mod ",round(modularity(commu_louvainsbn[[i]]),3)),
                                    colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                                    lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom)))
plotname <- paste0("figs/HD_impact/HD_communities_optimum")
pdf(paste0(plotname,".pdf"),width = 10, height = 5)
plotclim
dev.off()



##############################################################
#
##############################################################
par(mfrow = c(1,3))
k <-1
plots <- list()

cuts <- c(4,6,8,10,14,19)
cuts <- c(4,20,30,40,50,60)
cuts <- 1:16
# commu_betsbn <- commu_betsbn_invw_int
# commu_betscn <- commu_betscn_unw_int
colRainbow <- brewer.pal(9,"Set1")
colRainbow<- colorRampPalette(colRainbow)(max(cuts))
for(j in 1:length(cuts)){
  i <- cuts[[j]]
  commu_bet12 <- cut_at(HD_clusters_bic,i)
 # g <- bnlist[[k]]
  ncom <- length(levels(factor(commu_bet12)))
  mem <- commu_bet12
  memClim <- quantity2clim(mem, "membership", sel.t2m, backperm = NULL)
  if(ncom <= 9){
    #colRainbow <- brewer.pal(ncom,"Set1")
    colRainbowplot <- colRainbow[1:ncom]
  } else {
    #colRainbow <- brewer.pal(9,"Set1")
    #colRainbow<- colorRampPalette(colRainbow)(ncom)
    colRainbowplot <- colRainbow[1:ncom]
  }
  #if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plots[[j]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                            set.min = NULL, set.max = NULL, 
                            lonCenter = 0, 
                            regions = TRUE,col.regions = colRainbowplot, at = 0:ncom, rev.colors = FALSE, 
                           # main = paste0("Comunities BN ",length(E(g))),
                            colorkey = list(col = colRainbowplot, width = 0.6, at = 0:ncom,
                                            lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
  
}

plots

plotname <- paste0("figs/HD_impact/HD_communities")
pdf(paste0(plotname,".pdf"),width = 20, height = 15)

do.call(grid.arrange,c(plots))
dev.off()
