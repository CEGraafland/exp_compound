##################################################################################################
# Propagation of evidence 
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

degrees <- "5d"
mask <- TRUE # Only land or land & sea
limit <- TRUE # include poles?
b <- 3 # How many levels TP & T2M:
k <- 0 # Which permutation of variables later
# IT <- paste0(k,"_2000_2100")
 IT <- paste0(k,"_2100_2200")
#IT <- paste0(k,"_1000_1100")
# IT <- paste0(k,"_1700_1800")
 IT <- paste0(k,"_4000_4100")
compound <- "HD" # Which type of compound Hot Dry, Cold Dry, Hot Wet, Cold Wet, all together, all.distinct
global_mean_temp <- TRUE 
global_mean_temp_cont <- FALSE # Extra node for global_mean_temp?
combi <- "only" # only Compound nodes or also "TPT2M" nodes
algo<- "tabu" # Which algorithm
score<- "aic" # Which score in algorithm
global_mean_temp_art <- TRUE
cut.art <- 30
scaleGT <- FALSE
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
############################################################################
#
############################################################################
# 
# save(list = paste0(algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_",k,"_",i,"_",j,"i"),
#      file = paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_10d_disc_",b,lim,comp,combi,superf,n1,"/perm",k,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,"_",k,"_",i,"_",j,"i.rda"))

loadIterationsComp <- function(permused = 0, algo = c("hc","tabu"),score = c("aic","bic"), it = NULL) {
  
  hc_list <- list.files(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n,"/perm",k), full.names = T)
  hc_names <- list.files(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n,"/perm",k))
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


tabu_bic_0 <- loadIterationsComp(permused = k,algo = algo, score = score,it = IT)
tabu_bic_0$tabu_disc_3_aic_lim_comp_HDonly_land_0_4000_4100i$nodes$GT
tabu_bic_0$tabu_disc_3_bic_lim_GMT_comp_HDonly_land_0_4000_4100i$nodes$GT


tabu_bic_0_fits<- lapply(tabu_bic_0, function(x) bn.fit(x = x, data = data))

#################################################################################
# Single evidence. V81 postive (+ +)
#################################################################################
fittedbn <- tabu_bic_0_fits[[1]]


PropagationComp2Comp <- function(baysnet, data.c, nodeComp, level.Comp, nodesEvidence, level.ev, perm, add.without = TRUE){
  ngrids<- length(nodes(baysnet))-1
  nodesEvents <- c(nodeComp)
  

  
  valueEvents <- c(paste0("(levels(data.c[,names(baysnet)[nodeComp]])[level.Comp])")
                   #,paste0("(levels(data.c[,names(baysnet)[ngrids+nodeComp]])[level.tp])")
                           )
  
  if(length(nodesEvidence)==1){
    valueEvidence <- c(paste0("(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]])"))
  } else if (length(nodesEvidence)>1){
    valueEvidence <- c(paste0("(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]])"),
                       paste0("(levels(data.c[,names(baysnet)[nodesEvidence[2]]])[level.ev[2]])"))
  }
  

  
  if (is.null(perm)) {nodesEventsRef <- nodesEvents} else {
    nodesEventsRef <- c()
    for (i in 1:length(nodesEvents)){
      nodesEventsRef[i] <- which(perm == nodesEvents[i])
    }
  }
  

  
  if (length(nodesEvidence) == 1) {
    # For LS:
    str2 <- paste0("(", names(baysnet)[nodesEvidence[1]],"==", valueEvidence[1], ")")
    #FOR LW:
    #str2 <- paste0("(", names(baysnet)[nodesEvidence[1]], "=", valueEvidence, ")")
    #str2 <- valueEvidence
    # str2 <- paste0("(", names(baysnet)[nodesEvidence[1]]," == ", valueEvidence[1], ")")
    probname <- paste0("P(",names(baysnet)[nodeComp],"= ",lvl.comp,"|",names(baysnet)[nodesEvidence[1]]," = ", level.ev[1], ")")
  } else if (length(nodesEvidence) > 1){
    proves <- c()
    j <- 1
    for (j in 1:length(nodesEvidence)){
      proves[j]<- paste0(names(baysnet)[nodesEvidence[j]],"== ", valueEvidence[j])
    }
    
    text <- "(("
    for (j in 1:(length(nodesEvidence)-1)){
      text <- paste0(text,proves[j],") & (")
    }
    text <- paste0(text,proves[length(nodesEvidence)],"))")
    str2 <- text
    
    
    probname <- paste0("P(",names(baysnet)[nodeComp],"= ",lvl.comp,"|")
    for (j in 1:(length(nodesEvidence)-1)){
      probname <- paste0(probname,names(baysnet)[nodesEvidence[j]],"=",level.ev[j]," & ")
    }
    probname <- paste0(probname,names(baysnet)[nodesEvidence[length(nodesEvidence)]],"=",level.ev[length(nodesEvidence)],")")
    
  }
  
  
  str <- paste0("(", names(baysnet)[nodesEvents],"==", valueEvents, ")", collapse = " & ")
  #str
  
  #cmd1 = paste0("cpquery(baysnet, ", str, ", ", str2, ", method = ","'ls'",",n = 1000000)")
  cmd1 = paste0("cpquery(baysnet, ", str, ", ", str2, ", method = ","'ls'",")")
  
  cmd3 = paste0("cpquery(baysnet, ", str, ", ", "TRUE", ", method = ","'ls'",")")
  #cmd1
  #cmd3

  envlist <- list(baysnet = baysnet,
                  data.c = data.c,
                  nodeComp = nodeComp,
                  level.Comp = level.Comp,
                  #level.tp = level.tp,
                  nodesEvidence = nodesEvidence,
                  level.ev = level.ev,
                  perm = perm,
                  ngrids = ngrids,
                  nodesEvents = nodesEvents,
                  valueEvents = valueEvents,
                  valueEvidence = valueEvidence,
                  nodesEventsRef = nodesEventsRef,
                  str = str,
                  str2 = str2,
                  cmd1 = cmd1,
                  cmd3 = cmd3)
  list2env(envlist, envir = .GlobalEnv)  
  
  with <- eval(parse(text = cmd1))
  if (isTRUE(add.without)) {without <- eval(parse(text = cmd3))} else {without <- numeric(length = 1) }
  
  attr(with, "probability") <- probname
  attr(without, "probability") <- paste0("P(",names(baysnet)[nodeComp],"= ",lvl.comp,")")
  df <- data.frame(names = paste0(names(baysnet)[nodesEventsRef], collapse = "&"), with = with, without = without)
  return(df)
  
}

# data.c <- data
# perm <- permutations[[k]]
# baysnet <- fittedbn
# #ngrids<- length(nodes(baysnet))/2
# nodeComp <- 1
# level.Comp <- 2
# # nodesEvents <- c(nodeComp)
# 
# # levelEvents <- 1
# #level.tp <- 1
# 
# # FOR LS: 
#  nodesEvidence <- length(nodes(fittedbn))
#  level.ev <- c(4)
# 
# # nodesEvidence <- c(81,648 + 81)
# # level.ev <- c(4,1)
# 
# # perm <- NULL
# # lvl.comp <- 1
# # add.without <- TRUE




########################################################
# compound HD events impact on other compound HD events
########################################################

ngrids <- length(nodes(fittedbn)) - 1
j<- length(nodes(fittedbn)) # THis is GT
dfs <- list()
lvl.comp <- 2
nE <- c(j)
lvl.ev <- c(4)


 for(i in 1:ngrids){
  
   dfs[[i]] <- PropagationComp2Comp(baysnet = fittedbn,
                                   nodeComp = i,
                                   data.c = data,
                                   level.Comp = lvl.comp,
                                   nodesEvidence = nE,
                                   level.ev = lvl.ev,
                                   perm = NULL, add.without = TRUE)
}

if(!dir.exists(paste0("results/propagation/",degrees,"/"))){dir.create(paste0("results/propagation/",degrees,"/"))}
if(!dir.exists(paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/"))){dir.create(paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/"))}

if(length(lvl.ev) ==1){
  assign(paste0("prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
  
}

 








##################################################################################################
# Propagation of evidence 
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

degrees <- "5d"
mask <- TRUE # Only land or land & sea
limit <- TRUE # include poles?
b <- 3 # How many levels TP & T2M:
k <- 0 # Which permutation of variables later
IT <- paste0(k,"_2100_2200")
IT <- paste0(k,"_1000_1100")
IT <- paste0(k,"_1700_1800")
compound <- "HW" # Which type of compound Hot Dry, Cold Dry, Hot Wet, Cold Wet, all together, all.distinct
global_mean_temp <- FALSE 
global_mean_temp_cont <- TRUE# Extra node for global_mean_temp?
combi <- "only" # only Compound nodes or also "TPT2M" nodes
algo<- "tabu" # Which algorithm
score<- "bic" # Which score in algorithm
global_mean_temp_art <- TRUE
cut.art <- 15
scaleGT <- FALSE

if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
if(!is.null(compound)){comp <- paste0("_comp_",compound)} else {comp <-""}
if(isTRUE(mask)){superf <- "_land"}else {superf <-"_landsea"}
if(isTRUE(global_mean_temp)){n1 <- "_GMT"} else if(isTRUE(global_mean_temp_cont)){n1 <- "_GMTcont"} else {n1 <- ""}
if(isTRUE(global_mean_temp_art)){n2 <- "_GMTart"} else {n2 <- ""}

if(isTRUE(global_mean_temp_cont) & !isTRUE(global_mean_temp_art) & score == "bic"){score <- "bic-cg"} 
if(isTRUE(global_mean_temp_cont) & !isTRUE(global_mean_temp_art) & score == "aic"){score <- "aic-cg"} 
if (isTRUE(global_mean_temp_art)){n <- n2} else {n <- n1}

if(!k==0){permutations <-get(load(paste0("data/permutations/permutations_",degrees,lim,superf,".rda")))}

data.k <- get(load(file =paste0("data/manipulate_aggr/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".rda")))
data <- data.k
if(scaleGT == TRUE) {data$GT <- scale(data$GT)}
if(!k==0){data <-data[,c(permutations[[k]],(length(permutations[[k]])+1):ncol(data))]}
############################################################################
#
############################################################################
# 
# save(list = paste0(algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_",k,"_",i,"_",j,"i"),
#      file = paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_10d_disc_",b,lim,comp,combi,superf,n1,"/perm",k,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,"_",k,"_",i,"_",j,"i.rda"))

loadIterationsComp <- function(permused = 0, algo = c("hc","tabu"),score = c("aic","bic"), it = NULL) {
  
  hc_list <- list.files(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n,"/perm",k), full.names = T)
  hc_names <- list.files(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n,"/perm",k))
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



tabu_bic_0 <- loadIterationsComp(permused = k,algo = algo, score = score,it = IT)
tabu_bic_0_fits<- lapply(tabu_bic_0, function(x) bn.fit(x = x, data = data))


fittedbn <- tabu_bic_0_fits[[1]]


#################################################################
# compound HD events impact on other compound HD events CONTUINUE
#################################################################
PropagationGTcont2Comp <- function(baysnet, data.c, nodeComp, level.Comp, nodesEvidence, valueEvidence, perm, add.without = TRUE){
  ngrids<- length(nodes(baysnet))-1
  nodesEvents <- c(nodeComp)
  
  
  with <- numeric(length = length(nodesEvents))
  without <- numeric(length = length(nodesEvents))


  valueEvents <- c(paste0("(levels(data.c[,names(baysnet)[nodeComp]])[level.Comp])")
                   #,paste0("(levels(data.c[,names(baysnet)[ngrids+nodeComp]])[level.tp])")
  )
  
  

  
  
  if (length(nodesEvidence) == 1) {
    str2 <- paste0("(", names(baysnet)[nodesEvidence[1]],">=", valueEvidence[1], ")")
    probname <- paste0("P(",names(baysnet)[nodeComp],"= ",lvl.comp,"|",names(baysnet)[nodesEvidence[1]]," >= ", valueEvidence[1], ")")
  }
  
  ### *** TO DO ***
  if (length(nodesEvidence) > 1){
    proves <- c()
    for (j in 1:length(nodesEvidence)){
      proves[j]<- paste0("V",nodesEvidence[j]," = ", valueEvidence[j])
    }
    
    text <- "list("
    for (j in 1:(length(nodesEvidence)-1)){
      text <- paste0(text,proves[j],",")
    }
    text <- paste0(text,proves[length(nodesEvidence)],")")
    str2 <- text
    
    
    probname <- paste0("P(V ", valueEvent,"|")
    for (j in 1:(length(nodesEvidence)-1)){
      probname <- paste0(probname,proves[j],",")
    }
    probname <- paste0(probname,proves[length(nodesEvidence)],")")
  }
  ####**
  
  
  # if(length(nodesEvidence)==1){
  #   valueEvidence <- c(paste0("(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]])"))
  # } else if (length(nodesEvidence)>1){
  #   valueEvidence <- c(paste0("(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]])"),
  #                      paste0("(levels(data.c[,names(baysnet)[nodesEvidence[2]]])[level.ev[2]])"))
  # }
  # 
  
  
  if (is.null(perm)) {nodesEventsRef <- nodesEvents} else {
    nodesEventsRef <- c()
    for (i in 1:length(nodesEvents)){
      nodesEventsRef[i] <- which(perm == nodesEvents[i])
    }
  }
  
  
  
  # if (length(nodesEvidence) == 1) {
  #   # For LS:
  #   str2 <- paste0("(", names(baysnet)[nodesEvidence[1]],">=", valueEvidence[1], ")")
  #   #FOR LW:
  #   #str2 <- paste0("(", names(baysnet)[nodesEvidence[1]], "=", valueEvidence, ")")
  #   #str2 <- valueEvidence
  #   # str2 <- paste0("(", names(baysnet)[nodesEvidence[1]]," == ", valueEvidence[1], ")")
  #   probname <- paste0("P(",names(baysnet)[nodeComp],"= ",lvl.comp,"|",names(baysnet)[nodesEvidence[1]]," = ", valueEvidence[1], ")")
  # } else if (length(nodesEvidence) > 1){
  #   proves <- c()
  #   j <- 1
  #   for (j in 1:length(nodesEvidence)){
  #     proves[j]<- paste0(names(baysnet)[nodesEvidence[j]],"== ", valueEvidence[j])
  #   }
  #   
  #   text <- "(("
  #   for (j in 1:(length(nodesEvidence)-1)){
  #     text <- paste0(text,proves[j],") & (")
  #   }
  #   text <- paste0(text,proves[length(nodesEvidence)],"))")
  #   str2 <- text
  #   
  #   
  #   probname <- paste0("P(",names(baysnet)[nodeComp],"= ",lvl.comp,"|")
  #   for (j in 1:(length(nodesEvidence)-1)){
  #     probname <- paste0(probname,names(baysnet)[nodesEvidence[j]],"=",valueEvidence[j]," & ")
  #   }
  #   probname <- paste0(probname,names(baysnet)[nodesEvidence[length(nodesEvidence)]],"=",valueEvidence[length(nodesEvidence)],")")
  #   
  # }
  # 
  
  str <- paste0("(", names(baysnet)[nodesEvents],"==", valueEvents, ")", collapse = " & ")
  #str
  
   cmd1 = paste0("cpquery(baysnet, ", str, ", ", str2, ", method = ","'ls'",",n = 100000)")
  #cmd1 = paste0("cpquery(baysnet, ", str, ", ", str2, ", method = ","'lw'",")")
  
  cmd3 = paste0("cpquery(baysnet, ", str, ", ", "TRUE", ", method = ","'ls'",")")
  #cmd1
  #cmd3
  
  envlist <- list(baysnet = baysnet,
                  data.c = data.c,
                  nodeComp = nodeComp,
                  level.Comp = level.Comp,
                  #level.tp = level.tp,
                  nodesEvidence = nodesEvidence,
                  # level.ev = level.ev,
                  perm = perm,
                  ngrids = ngrids,
                  nodesEvents = nodesEvents,
                  valueEvents = valueEvents,
                  valueEvidence = valueEvidence,
                  nodesEventsRef = nodesEventsRef,
                  str = str,
                  str2 = str2,
                  cmd1 = cmd1,
                  cmd3 = cmd3)
  list2env(envlist, envir = .GlobalEnv)  
  
  with <- eval(parse(text = cmd1))
  if (isTRUE(add.without)) {without <- eval(parse(text = cmd3))} else {without <- numeric(length = 1) }
  
  attr(with, "probability") <- probname
  attr(without, "probability") <- paste0("P(",names(baysnet)[nodeComp],"= ",lvl.comp,")")
  df <- data.frame(names = paste0(names(baysnet)[nodesEventsRef], collapse = "&"), with = with, without = without)
  return(df)
  
}


ngrids <- length(nodes(fittedbn)) - 1
j<- length(nodes(fittedbn)) # THis is GT
grados <- 2
lvl.comp <- 2
nE <- c(j)
if(isTRUE(scaleGT)){value.ev <- c(mean(data[["GT"]]) + grados/sd(data.k[["GT"]]))} else {value.ev <- c(mean(data[["GT"]]) + grados)} 
names(value.ev)<- "GT"

dfs <- list()
i <- 1
for(i in 1:ngrids){
  dfs[[i]]<- PropagationGTcont2Comp(baysnet = fittedbn,
                                  nodeComp = i,
                                  data.c = data,
                                  level.Comp = lvl.comp,
                                  nodesEvidence = nE,
                                  valueEvidence = value.ev,
                                  perm = NULL, add.without = TRUE)
}

if(!dir.exists(paste0("results/propagation/",degrees,"/"))){dir.create(paste0("results/propagation/",degrees,"/"))}
if(!dir.exists(paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/"))){dir.create(paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/"))}

if(length(value.ev) ==1){
  assign(paste0("prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",value.ev),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",value.ev),file = paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",value.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",value.ev[1],"_",names(fittedbn)[nE][2],"_",value.ev[2]),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",value.ev[1],"_",names(fittedbn)[nE][2],"_",value.ev[2]),file = paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",value.ev[1],"_",names(fittedbn)[nE][2],"_",value.ev[2],".rda"))
  
}










##################################################################################################
# Propagation of evidence 
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

degrees <- "5d"
mask <- TRUE # Only land or land & sea
limit <- TRUE # include poles?
b <- 3 # How many levels TP & T2M:
k <- 0 # Which permutation of variables later
IT <- paste0(k,"_2100_2200")
IT <- paste0(k,"_1700_1800")
IT <- paste0(k,"_1000_1100")
compound <- "HD" # Which type of compound Hot Dry, Cold Dry, Hot Wet, Cold Wet, all together, all.distinct
global_mean_temp <- FALSE 
global_mean_temp_cont <- TRUE# Extra node for global_mean_temp?
combi <- "only" # only Compound nodes or also "TPT2M" nodes
algo<- "tabu" # Which algorithm
score<- "bic" # Which score in algorithm
global_mean_temp_art <- TRUE
cut.art <- 15
scaleGT <- FALSE

if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
if(!is.null(compound)){comp <- paste0("_comp_",compound)} else {comp <-""}
if(isTRUE(mask)){superf <- "_land"}else {superf <-"_landsea"}
if(isTRUE(global_mean_temp)){n1 <- "_GMT"} else if(isTRUE(global_mean_temp_cont)){n1 <- "_GMTcont"} else {n1 <- ""}
if(isTRUE(global_mean_temp_art)){n2 <- "_GMTart"} else {n2 <- ""}

if(isTRUE(global_mean_temp_cont) & !isTRUE(global_mean_temp_art) & score == "bic"){score <- "bic-cg"} 
if(isTRUE(global_mean_temp_cont) & !isTRUE(global_mean_temp_art) & score == "aic"){score <- "aic-cg"} 
if (isTRUE(global_mean_temp_art)){n <- n2} else {n <- n1}

if(!k==0){permutations <-get(load(paste0("data/permutations/permutations_",degrees,lim,superf,".rda")))}

data.k <- get(load(file =paste0("data/manipulate_aggr/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".rda")))
data <- data.k
if(scaleGT == TRUE) {data$GT <- scale(data$GT)}
if(!k==0){data <-data[,c(permutations[[k]],(length(permutations[[k]])+1):ncol(data))]}
############################################################################
#
############################################################################
# 
# save(list = paste0(algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_",k,"_",i,"_",j,"i"),
#      file = paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_10d_disc_",b,lim,comp,combi,superf,n1,"/perm",k,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,"_",k,"_",i,"_",j,"i.rda"))


loadIterationsComp <- function(permused = 0, algo = c("hc","tabu"),score = c("aic","bic"), it = NULL) {
  
  hc_list <- list.files(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n,"/perm",k), full.names = T)
  hc_names <- list.files(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n,"/perm",k))
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


# IT <- paste0(k,"_700_800") #CD
#IT <- paste0(k,"_1000_1100") #HW
# IT <- paste0(k,"_2100_2200") #CW
# IT <- paste0(k,"_1700_1800") #HD

tabu_bic_0 <- loadIterationsComp(permused = k,algo = algo, score = score,it = IT)
tabu_bic_0_fits<- lapply(tabu_bic_0, function(x) bn.fit(x = x, data = data))

fittedbn <- tabu_bic_0_fits[[1]]

##########################################################################
#
##########################################################################
meth <- "lw"
PropagationGTcont2CompLW <- function(baysnet, data.c, nodeComp, level.Comp, nodesEvidence, level.ev, perm, add.without = TRUE){
  ngrids<- length(nodes(baysnet))-1
  nodesEvents <- c(nodeComp)
  
  
  if (is.null(perm)) {nodesEventsRef <- nodesEvents} else {
    nodesEventsRef <- c()
    for (i in 1:length(nodesEvents)){
      nodesEventsRef[i] <- which(perm == nodesEvents[i])
    }
  }
  
  
  valueEvents <- c(paste0("(levels(data.c[,names(baysnet)[nodeComp]])[level.Comp])")
                   #,paste0("(levels(data.c[,names(baysnet)[ngrids+nodeComp]])[level.tp])")
  )
  
  # if(length(nodesEvidence)==1){
  #   valueEvidence <- c(paste0("(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]])"))
  # } else if (length(nodesEvidence)>1){
  #   valueEvidence <- c(paste0("(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]])"),
  #                      paste0("(levels(data.c[,names(baysnet)[nodesEvidence[2]]])[level.ev[2]])"))
  # }
  
  valueEvidence <- level.ev
  
  if (length(nodesEvidence) == 1) {
    str2 <- paste0("list(", names(baysnet)[nodesEvidence[1]]," = ", valueEvidence[1], ")")
    probname <- paste0("P(C ", valueEvents,"|",names(baysnet)[nodesEvidence[1]]," = ", valueEvidence[1], ")")
  }
  

  str <- paste0("(", names(baysnet)[nodesEvents],"==", valueEvents, ")", collapse = " & ")
  #str
  
  #cmd1 = paste0("cpquery(baysnet, ", str, ", ", str2, ", method = ","'ls'",",n = 1000000)")
  #cmd1 = paste0("cpquery(baysnet, ", str, ", ", str2, ", method = ","'lw'",",debug = TRUE)")
  cmd1 = paste0("cpquery(baysnet, ", str, ", ", str2, ", method = ","'lw'",")")
  
  cmd3 = paste0("cpquery(baysnet, ", str, ", ", "TRUE", ", method = ","'lw'",")")
  #cmd1
  #cmd3
  
  envlist <- list(baysnet = baysnet,
                  data.c = data.c,
                  nodeComp = nodeComp,
                  level.Comp = level.Comp,
                  #level.tp = level.tp,
                  nodesEvidence = nodesEvidence,
                  level.ev = level.ev,
                  perm = perm,
                  ngrids = ngrids,
                  nodesEvents = nodesEvents,
                  valueEvents = valueEvents,
                  valueEvidence = valueEvidence,
                  nodesEventsRef = nodesEventsRef,
                  str = str,
                  str2 = str2,
                  cmd1 = cmd1,
                  cmd3 = cmd3)
  list2env(envlist, envir = .GlobalEnv)  
  
  with <- eval(parse(text = cmd1))
  if (isTRUE(add.without)) {without <- eval(parse(text = cmd3))} else {without <- numeric(length = 1) }
  
  attr(with, "probability") <- probname
  attr(without, "probability") <- paste0("P(",names(baysnet)[nodeComp],"= ",lvl.comp,")")
  df <- data.frame(names = paste0(names(baysnet)[nodesEventsRef], collapse = "&"), with = with, without = without)
  return(df)
  
}

ngrids <- length(nodes(fittedbn)) - 1
j<- length(nodes(fittedbn)) # THis is GT
grados <- 4
lvl.comp <- 2
nE <- c(j)
if(isTRUE(scaleGT)){value.ev <- c(mean(data[["GT"]]) + grados/sd(data.k[["GT"]]))} else {value.ev <- c(mean(data[["GT"]]) + grados)} 
names(value.ev)<- "GT"

dfs <- list()
i <- 1
for(i in 1:ngrids){
  dfs[[i]]<- PropagationGTcont2CompLW(baysnet = fittedbn,
                                    nodeComp = i,
                                    data.c = data,
                                    level.Comp = lvl.comp,
                                    nodesEvidence = nE,
                                    level.ev  = value.ev,
                                    perm = NULL, add.without = TRUE)
}

if(!dir.exists(paste0("results/propagation/",degrees,"/"))){dir.create(paste0("results/propagation/",degrees,"/"))}
if(!dir.exists(paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/"))){dir.create(paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/"))}

if(length(value.ev) ==1){
  assign(paste0("prop",meth,"_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",value.ev),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop",meth,"_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",value.ev),file = paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop",meth,"_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",value.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop",meth,"_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",value.ev[1],"_",names(fittedbn)[nE][2],"_",value.ev[2]),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",value.ev[1],"_",names(fittedbn)[nE][2],"_",value.ev[2]),file = paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop",meth,"_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",value.ev[1],"_",names(fittedbn)[nE][2],"_",value.ev[2],".rda"))
  
}

