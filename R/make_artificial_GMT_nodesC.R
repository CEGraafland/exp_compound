##################################################################################################
# Make artificial GMT nodes
##################################################################################################
rm(list = ls())
setwd("~/data/Untitled/Trabajo/R_practice/exp_compound")
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
library(RColorBrewer)
library(gridExtra)
library(visualizeR)
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
# load("data/raw_aggr/t2m_ERA5_monthly_1940_2022_10d.rda")
# load("data/raw_aggr/tp_ERA5_monthly_1940_2022_10d.rda")
########################################################################################################
# Load functions and data
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/propagationFunctions.R")
source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
#load("data/raw_aggr/t2m_ERA5_monthly_1940_2022_10d.rda")
#load("data/raw_aggr/tp_ERA5_monthly_1940_2022_10d.rda")

degrees <- "5d"
mask <- TRUE # Only land or land & sea
limit <- TRUE # include poles?
b <- 3 # How many levels TP & T2M:
k <- 0 # Which permutation of variables later
IT <- paste0(k,"_1700_1800")
IT <- paste0(k,"_2000_2100")
IT<- paste0(k,"_1000_1100")
IT <- paste0(k,"_2100_2200")
IT <- paste0(k,"_4000_4100")
compound <- "HD" # Which type of compound Hot Dry, Cold Dry, Hot Wet, Cold Wet, all together, all.distinct
global_mean_temp <- TRUE
global_mean_temp_cont <- FALSE# Extra node for global_mean_temp?
combi <- "only" # only Compound nodes or also "TPT2M" nodes
algo<- "tabu" # Which algorithm
score<- "aic" # Which score in algorithm
cut.art <- 30
geo <- "auto"
# quantsel <-"betw"

if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
if(!is.null(compound)){comp <- paste0("_comp_",compound)} else {comp <-""}
if(isTRUE(mask)){superf <- "_land"}else {superf <-"_landsea"}
if(isTRUE(global_mean_temp) & !geo == "auto"){n1 <- "_GMT"} else if(isTRUE(global_mean_temp_cont)){n1 <- "_GMTcont"} else{n1 <- ""}
# if (isTRUE(global_mean_temp_cont) & score == "bic"){score <- "bic-cg"} 

####################################
# mask
####################################
degrees <- "5d" #10d
mask <- TRUE
# if(isTRUE(mask)){
#   t2m_ERA5 <- get(load(paste0("data/mask_aggr/t2m_land_ERA5_monthly_1940_2022_",degrees,".rda")))
#   tp_ERA5 <- get(load(paste0("data/mask_aggr/tp_land_ERA5_monthly_1940_2022_",degrees,".rda")))
#   landmask <- get(load(paste0("data/mask_aggr/landmask_",degrees,".rda")))
#   if(isTRUE(limit)){landmask <- subsetGrid(landmask,latLim = c(-50,75))
#   t2m_ERA5 <- subsetGrid(t2m_ERA5,latLim = c(-50,75))
#   tp_ERA5 <-  subsetGrid(tp_ERA5,latLim = c(-50,75)) }
# } else {
#   t2m_ERA5 <- get(load(paste0("data/raw_aggr/t2m_ERA5_monthly_1940_2022_",degrees,".rda")))
#   tp_ERA5 <- get(load(paste0("data/raw_aggr/tp_ERA5_monthly_1940_2022_",degrees,".rda")))
#   if(isTRUE(limit)){
#     t2m_ERA5 <- subsetGrid(t2m_ERA5,latLim = c(-50,75))
#     tp_ERA5 <-  subsetGrid(tp_ERA5,latLim = c(-50,75)) }
# }


# if(!k==0){permutations <-get(load(paste0("data/permutations/permutations_",degrees,lim,superf,".rda")))
# backpermutations <- get(load(paste0("data/permutations/backpermutations_",degrees,lim,superf,".rda")))}
# 
# data <- get(load(file =paste0("data/manipulate_aggr/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".rda")))
# if(!k==0){data <-data[,c(permutations[[k]],(length(permutations[[k]])+1):ncol(data))]}
############################################################################
#
############################################################################
loadIterationsComp <- function(permused = 0, algo = c("hc","tabu"),score = c("aic","bic"), it = NULL) {
  
  hc_list <- list.files(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"/perm",k), full.names = T)
  hc_names <- list.files(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"/perm",k))
  hc_names <- gsub(".rda", "", hc_names)
  
  if(!is.null(it)){
    hc_list <- hc_list[grep(it,hc_list)]
    hc_names <- hc_names[grep(it,hc_names)]
  }
  
  hc_networks <- lapply(hc_list, function(x){get(load(x))})
  names(hc_networks) <- hc_names
  sizes <- sapply(hc_networks,narcs)
  hc_networks <- hc_networks[order(sizes)]
  return(hc_networks)
}


tabu_bic_0 <-loadIterationsComp(permused = k,algo = algo, score = score,it = IT)

if(is.null(geo)){attach.toGMT <- get(load(file = paste0("results/nodes_in_com/df.com.nodes_",degrees,"_disc_",b,lim,comp,combi,superf,n1,"_",IT,"_cut.art_",cut.art,".rda")))
} else if (geo == "totdegree") {igraph.bn <- igraph.from.graphNEL(as.graphNEL(tabu_bic_0[[1]]))
attach.toGMT <- data.frame("nodes" = names(V(igraph.bn)[order(degree(igraph.bn),decreasing = TRUE)][1:cut.art]))
} else if (geo == "totclose") {igraph.bn <- igraph.from.graphNEL(as.graphNEL(tabu_bic_0[[1]]))
attach.toGMT <- data.frame("nodes" = names(V(igraph.bn)[order(closeness(igraph.bn),decreasing = TRUE)][1:cut.art]))
} else if (geo == "totbetw") {igraph.bn <- igraph.from.graphNEL(as.graphNEL(tabu_bic_0[[1]]))
attach.toGMT <- data.frame("nodes" = names(V(igraph.bn)[order(betweenness(igraph.bn),decreasing = TRUE)][1:cut.art]))
} else if (geo == "auto") {
  global_mean_temp <- TRUE
  n1 <- "_GMT"
  load("/data/Untitled/Trabajo/R_practice/exp_compound/data/tabuiterations/ERA5_monthly_1940_2022_5d_disc_3_aic_lim_comp_HDonly_land/perm0/tabu_disc_3_aic_lim_comp_HDonly_land_0_4000_4100i.rda")
  data.k <- get(load(file =paste0("data/manipulate_aggr/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".rda")))
  data <- data.k
  if (algo == "tabu"){
    
    steps <- cut.art 
    start <- tabu_disc_3_aic_lim_comp_HDonly_land_0_4000_4100i
    start.edges <- arcs(start)
    start.igraph <- graph_from_edgelist(start.edges)
    all.edges <-as.matrix(expand.grid(nodes(start),nodes(start),stringsAsFactors = FALSE))
    full.igraph <- graph_from_edgelist(all.edges)
    
    attributes(V(start.igraph))
    attributes(V(full.igraph))
    
    test <- difference(full.igraph,start.igraph)
    
    # 
    # apply(start.edges[1:2,], 1, function(x, want) isTRUE(all.equal(x, want)), all.edges[1,])
    # apply(start.edges,1,function(y) apply(all.edges, 1, function(x, want) isTRUE(all.equal(x, want)), y))
    # blacklist <- which(start.edges %in%all.edges)
    blacklist <- as_edgelist(test)
    start <- add.node(start,"GT")
    berekening <- tabu(data, max.iter = steps, score = "loglik", start = start, blacklist = blacklist)
  }
  tabu_bic_0_art <- list(berekening)
  names(tabu_bic_0_art) <- names(tabu_bic_0)
  
}else {attach.toGMT <- get(load(file = paste0("results/nodes_in_com/df.com.nodes_",degrees,"_disc_",b,lim,comp,combi,superf,n1,"_",IT,"_cut.art_",cut.art,geo,".rda")))
        attach.toGMT <- data.frame("nodes" = unlist(attach.toGMT)) 
        }

art.GMT.dag <- function(g,df.comp.nodes,mode = c("GTC","CGT")){
  g.comp <- g
  g.comp <- add.node(g.comp,"GT")
  for(i in 1:nrow(df.comp.nodes)){
    if (mode == "CGT"){g.comp <- set.arc(g.comp, df.comp.nodes[i,],"GT")
    } else if (mode == "GTC") {g.comp <- set.arc(g.comp,"GT", df.comp.nodes[i,])}
  }
  return(g.comp)
}

art.RMT.dag <- function(g,df.comp.nodes){
  g.comp <- g
  for(i in 1:nrow(df.comp.nodes)){
  g.comp <- add.node(g.comp,paste0("RT",i))}
  
  for(i in 1:nrow(df.comp.nodes)){
    g.comp <- set.arc(g.comp, nodes(g.comp)[df.comp.nodes[i,]],paste0("RT",i))
  }
  return(g.comp)
}

if(is.null(tabu_bic_0_art)){tabu_bic_0_art <- lapply(tabu_bic_0,art.GMT.dag, df.comp.nodes = na.omit(attach.toGMT), mode = "GTC")}


##############################################################################
#
##############################################################################
names(tabu_bic_0_art)<- gsub("lim",paste0("lim_GMT_art",cut.art),names(tabu_bic_0_art))
names(tabu_bic_0_art)

if(geo == "auto"){n1 <- ""}
if(!dir.exists(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_GMTart_",cut.art,geo,"/"))){dir.create(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_GMTart_",cut.art,geo,"/"))}
if(!dir.exists(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_GMTart_",cut.art,geo,"/perm",k))){ dir.create(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_GMTart_",cut.art,geo,"/perm",k))}


for (x in 1:length(tabu_bic_0_art)){
  assign(names(tabu_bic_0_art)[x],tabu_bic_0_art[[x]])
  save(list = names(tabu_bic_0_art)[x],file = paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_GMTart_",cut.art,geo,"/perm",k,"/",names(tabu_bic_0_art)[x],".rda"))
}


##############################################################################
#
##############################################################################
tabu_disc_3_bic_lim_comp_HDonly_land_GMT_0_2100_2200i$nodes$GT$nbr
tabu_disc_3_bic_lim_GMT_art12_comp_HDonly_land_0_2100_2200i$nodes$GT$nbr
save(tabu_disc_3_bic_lim_GMT_art12_comp_HDonly_land_0_2100_2200i,file = paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_GMTart_",cut.art,geo,"/perm",k,"/",names(tabu_bic_0_art)[x],".rda"))


tabu_disc_3_bic_lim_GMT_art12_comp_HDonly_land_0_2100_2200i <- set.arc(tabu_disc_3_bic_lim_GMT_art12_comp_HDonly_land_0_2100_2200i, "C396","GT")
tabu_disc_3_bic_lim_GMT_art12_comp_HDonly_land_0_2100_2200i <- set.arc(tabu_disc_3_bic_lim_GMT_art12_comp_HDonly_land_0_2100_2200i, "C456","GT")
tabu_disc_3_bic_lim_GMT_art12_comp_HDonly_land_0_2100_2200i <- set.arc(tabu_disc_3_bic_lim_GMT_art12_comp_HDonly_land_0_2100_2200i, "C624","GT")

landmask <- get(load(paste0("data/mask_aggr/landmask_",degrees,".rda")))
landmask <- subsetGrid(landmask,latLim = c(-50,75))


textdf <- cbind(attr(sel.t2m.1,"VertexCoords"),1:ngrids)
#textdf$x[textdf$x<0] <-textdf$x[textdf$x<0] + 360
textdf <- as.matrix(textdf)
textlay <- apply(textdf,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

clim.0vector <- quantity2clim(numeric(length=ngrids),"node vector", ref.grid = sel.t2m)
Xvectorname <-attr(clim.0vector$Data,"climatology:fun")
Xvector <-spatialPlot(clim.0vector, main = Xvectorname,  backdrop.theme = "coastline",sp.layout = textlay, lonCenter =0, colorkey = FALSE)

mask_loc <- as.vector(array3Dto2Dmat(redim(landmask,member = FALSE)$Data))
mask_loc_perm <- mask_loc
mask_loc_perm[!is.na(mask_loc)]<- 1:length(mask_loc_perm[!is.na(mask_loc)])

t2m_ERA5 <- get(load(paste0("data/mask_aggr/t2m_land_ERA5_monthly_1940_2022_",degrees,".rda")))
t2m_ERA5 <- subsetGrid(t2m_ERA5,latLim = c(-50,75))
TCA_t2m_ERA5 <- TimeCoordsAnom_from_Grid(t2m_ERA5)

textdf <- cbind(attr(TCA_t2m_ERA5,"VertexCoords"),mask_loc_perm)
#textdf$x[textdf$x<0] <-textdf$x[textdf$x<0] + 360
textdf <- as.matrix(textdf)
textlay <- apply(textdf,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.5,col = "black"))

clim.0vector <- quantity2clim(mask_loc_perm,"node vector", ref.grid = landmask)
Xvectorname <-attr(clim.0vector$Data,"climatology:fun")
Xvector <-spatialPlot(clim.0vector, main = Xvectorname,  backdrop.theme = "coastline",sp.layout = textlay, lonCenter =0, colorkey = FALSE)
Xvector

################################################################
#
################################################################
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

degrees <- "10d"
mask <- FALSE # Only land or land & sea
limit <- TRUE # include poles?
b <- 3 # How many levels TP & T2M:
k <- 0 # Which permutation of variables later
IT <- paste0(k,"_2700_2800")
compound <- "HD" # Which type of compound Hot Dry, Cold Dry, Hot Wet, Cold Wet, all together, all.distinct
global_mean_temp <- FALSE
global_mean_temp_cont <- FALSE# Extra node for global_mean_temp?
combi <- "TPT2M" # only Compound nodes or also "TPT2M" nodes
algo<- "tabu" # Which algorithm
score<- "aic" # Which score in algorithm
global_mean_temp_art <- FALSE
scaleGT <- FALSE
cut.art <- Inf
geo <- "auto"

if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
if(!is.null(compound)){comp <- paste0("_comp_",compound)} else {comp <-""}
if(isTRUE(mask)){superf <- "_land"}else {superf <-"_landsea"}
if(isTRUE(global_mean_temp)){n1 <- "_GMT"} else if(isTRUE(global_mean_temp_cont)){n1 <- "_GMTcont"} else {n1 <- ""}
if(isTRUE(global_mean_temp_art)){n2 <- paste0("_GMTart_",cut.art,geo)} else {n2 <- ""}

if(isTRUE(global_mean_temp_cont) & !isTRUE(global_mean_temp_art) & score == "bic"){score <- "bic-cg"} 
if(isTRUE(global_mean_temp_cont) & !isTRUE(global_mean_temp_art) & score == "aic"){score <- "aic-cg"} 
if (isTRUE(global_mean_temp_art)){n <- n2} else {n <- n1}

if(!k==0){permutations <-get(load(paste0("data/permutations/permutations_",degrees,lim,superf,".rda")))}

data.k <- get(load(file =paste0("data/manipulate_aggr/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".rda")))
data <- data.k
if(scaleGT == TRUE) {data$GT <- scale(data$GT)}
if(!k==0){data <-data[,c(permutations[[k]],(length(permutations[[k]])+1):ncol(data))]}
######################################################################################
#
######################################################################################
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
Score <- "aic"
IT<- "0_2700_2800"
tabu_bic_0_comp_art <- loadIterationsComp(
  permused = 0,algo = algo, score = score,
  limit = limit,compound = "comp",artificial = TRUE, it = IT)

tabu_aic_0_comp_art_fit<- lapply(tabu_bic_0_comp_art, function(x) bn.fit(x = x, data = data))
#############################################################################################
#
#############################################################################################
ngrids <- length(nodes(tabu_bic_0_comp_art$tabu_disc_3_bic_t2m_tp_10d_lim_comp_art_0_1700_1800i))/3
ngrids <- length(nodes(tabu_bic_0_comp_art$tabu_disc_3_aic_t2m_tp_10d_lim_comp_art_0_2700_2800i))/3
if(geo == "geo"){
  attach.toGMT<- data.frame("t2m-nodes" = c(7,47,192,362,466,424,29,163,185,45,266,161,257,406),
                          "tp-nodes" = c(7,47,192,362,466,424,29,163,185,45,266,161,257,406)+ngrids,"C-nodes" = c(7,47,192,362,466,424,29,163,185,45,266,161,257,406)+2*ngrids)

  art.GMT.dag <- function(g,df.comp.nodes){
    
    g.comp <- g
    g.comp <- add.node(g.comp,"GT")
    for(i in 1:nrow(df.comp.nodes)){
      g.comp <- set.arc(g.comp, "GT",nodes(g.comp)[df.comp.nodes[i,1]])
      g.comp <- set.arc(g.comp, "GT",nodes(g.comp)[df.comp.nodes[i,2]])
      g.comp <- set.arc(g.comp, "GT",nodes(g.comp)[df.comp.nodes[i,3]])
    }
    return(g.comp)
  }
  
  tabu_bic_0_art <- lapply(tabu_bic_0_comp_art,art.GMT.dag, df.comp.nodes = attach.toGMT)
  names(tabu_bic_0_art)<- gsub("lim",paste0("lim_GMT_art",cut.art),names(tabu_bic_0_art))
  
  if(!dir.exists(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_GMTart_",cut.art,geo,"/"))){dir.create(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_GMTart_",cut.art,geo,"/"))}
  if(!dir.exists(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_GMTart_",cut.art,geo,"/perm",k))){ dir.create(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_GMTart_",cut.art,geo,"/perm",k))}
  } else if (geo == "auto") {
  global_mean_temp <- TRUE
  n1 <- "_GMT"
  # load("/data/Untitled/Trabajo/R_practice/exp_compound/data/tabuiterations/ERA5_monthly_1940_2022_5d_disc_3_aic_lim_comp_HDonly_land/perm0/tabu_disc_3_aic_lim_comp_HDonly_land_0_4000_4100i.rda")
  data.k <- get(load(file =paste0("data/manipulate_aggr/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".rda")))
  data <- data.k
  if (algo == "tabu"){
    
    steps <- cut.art 
    start <- tabu_bic_0_comp_art$tabu_disc_3_aic_t2m_tp_10d_lim_comp_art_0_2700_2800i
    start.edges <- arcs(start)
    start.igraph <- graph_from_edgelist(start.edges)
    all.edges <-as.matrix(expand.grid(nodes(start),nodes(start),stringsAsFactors = FALSE))
    full.igraph <- graph_from_edgelist(all.edges)
    
    attributes(V(start.igraph))
    attributes(V(full.igraph))
    
    test <- difference(full.igraph,start.igraph)
    
    # 
    # apply(start.edges[1:2,], 1, function(x, want) isTRUE(all.equal(x, want)), all.edges[1,])
    # apply(start.edges,1,function(y) apply(all.edges, 1, function(x, want) isTRUE(all.equal(x, want)), y))
    # blacklist <- which(start.edges %in%all.edges)
    blacklist <- as_edgelist(test)
    start <- add.node(start,"GT")
    if(!cut.art == Inf){
      berekening <- tabu(data, max.iter = steps, score = "loglik", start = start, blacklist = blacklist)
    } else {berekening <- tabu(data, max.iter = steps, score = "aic", start = start, blacklist = blacklist)}
  }   
  tabu_bic_0_art <- list(berekening)
  names(tabu_bic_0_art)<- gsub("lim",paste0("lim_GMT_art",cut.art),names(tabu_bic_0_comp_art))
  
  
  if(!dir.exists(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,"_GMTart_",cut.art,geo,"/"))){dir.create(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,"_GMTart_",cut.art,geo,"/"))}
  if(!dir.exists(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,"_GMTart_",cut.art,geo,"/perm",k))){ dir.create(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,"_GMTart_",cut.art,geo,"/perm",k))}
  
  for (x in 1:length(tabu_bic_0_comp_art)){
    assign(names(tabu_bic_0_art)[x],tabu_bic_0_art[[x]])
    save(list = names(tabu_bic_0_art)[x],file = paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,"_GMTart_",cut.art,geo,"/perm",k,"/",names(tabu_bic_0_art)[x],".rda"))
  }
  }  else if (geo == "autoext") {
    global_mean_temp <- TRUE
    n1 <- "_GMT"
    # load("/data/Untitled/Trabajo/R_practice/exp_compound/data/tabuiterations/ERA5_monthly_1940_2022_5d_disc_3_aic_lim_comp_HDonly_land/perm0/tabu_disc_3_aic_lim_comp_HDonly_land_0_4000_4100i.rda")
    
    data.k <- get(load(file =paste0("data/manipulate_aggr/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".rda")))
    data <- data.k
    if (algo == "tabu"){
      
      steps <- cut.art 
      start <- tabu_bic_0_comp_art$tabu_disc_3_aic_t2m_tp_10d_lim_comp_art_0_2700_2800i
      start.edges <- arcs(start)
      start.igraph <- graph_from_edgelist(start.edges)
      all.edges <-as.matrix(expand.grid(nodes(start),nodes(start),stringsAsFactors = FALSE))
     # HIER VERDER GAAN: OOK COMBINATIES GT naar T2M, GT naar TP, Verbieden
     
       full.igraph <- graph_from_edgelist(all.edges)
      
      attributes(V(start.igraph))
      attributes(V(full.igraph))
      
      test <- difference(full.igraph,start.igraph)
      
      # load("/data/Untitled/Trabajo/R_practice/exp_compound/results/propagation/spatext_tabu_disc_3_aic_lim_comp_art_0_2700_2700/spatext_tabu_aic_0_2700_2800_C_2_C_2.rda")
      load("/data/Untitled/Trabajo/R_practice/exp_compound/results/propagation/spatext_tabu_disc_3_aic_lim_comp_art_0_2700_2800/spatext_tabu_aic_0_2700_2800_C_2_C_2.rda")
      load("/data/Untitled/Trabajo/R_practice/exp_compound/results/propagation/spatext_tabu_disc_3_aic_lim_comp_art_0_2700_2800/spatext_tabu_aic_0_2700_2800_C_2_t2m_3.rda")
      # load(file = paste0("results/propagation/spatext_",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/spatext_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_C_",lvl.ev,".rda"))
     o1 <- order(spatextHtoHD,decreasing = TRUE)
     Vo1<- o1[(cut.art+1):length(o1)]
      
     Vnodes <- nodes(start)[c(Vo1,Vo1+ngrids,Vo1+2*ngrids)]
     GTV.edges.1 <-as.matrix(expand.grid("GT",Vnodes,stringsAsFactors = FALSE))
     GTV.edges.2 <-as.matrix(expand.grid(Vnodes,"GT",stringsAsFactors = FALSE))
     GTV.igraph <- graph_from_edgelist(rbind(GTV.edges.1,GTV.edges.2))
    
      # apply(start.edges[1:2,], 1, function(x, want) isTRUE(all.equal(x, want)), all.edges[1,])
      # apply(start.edges,1,function(y) apply(all.edges, 1, function(x, want) isTRUE(all.equal(x, want)), y))
      # blacklist <- which(start.edges %in%all.edges)
      blacklist.1 <- as_edgelist(test)
      blacklist.2 <- as_edgelist(GTV.igraph)
      blacklist <- rbind(blacklist.1,blacklist.2)
     
      start <- add.node(start,"GT")
      
      berekening <- tabu(data, max.iter = steps, score = "loglik", start = start, blacklist = blacklist)
    }
  tabu_bic_0_art <- list(berekening)
  names(tabu_bic_0_art)<- gsub("lim",paste0("lim_GMT_art",cut.art),names(tabu_bic_0_comp_art))

  
  if(!dir.exists(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,"_GMTart_",cut.art,geo,"/"))){dir.create(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,"_GMTart_",cut.art,geo,"/"))}
  if(!dir.exists(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,"_GMTart_",cut.art,geo,"/perm",k))){ dir.create(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,"_GMTart_",cut.art,geo,"/perm",k))}
  
  for (x in 1:length(tabu_bic_0_comp_art)){
    assign(names(tabu_bic_0_art)[x],tabu_bic_0_art[[x]])
    save(list = names(tabu_bic_0_art)[x],file = paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,"_GMTart_",cut.art,geo,"/perm",k,"/",names(tabu_bic_0_art)[x],".rda"))
  }
  }
# art.GMT.dag <- function(g,df.comp.nodes){
#   g.comp <- g
#   g.comp <- add.node(g.comp,"GT")
#   for(i in 1:nrow(df.comp.nodes)){
#     g.comp <- set.arc(g.comp, nodes(g.comp)[df.comp.nodes[i,]],"GT")
#   }
#   return(g.comp)
# }





tabu_bic_0_art[[1]]$nodes
max(sapply(tabu_bic_0_art[[1]]$nodes, function(x)length(x$parents)))

for (x in 1:length(tabu_bic_0_art)){
  assign(names(tabu_bic_0_art)[x],tabu_bic_0_art[[x]])
  save(list = names(tabu_bic_0_art)[x],file = paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_GMTart_",cut.art,geo,"/perm",k,"/",names(tabu_bic_0_art)[x],".rda"))
}

tabu_disc_3_aic_t2m_tp_10d_lim_GMT_art100_comp_art_0_2700_2800i$nodes$GT
