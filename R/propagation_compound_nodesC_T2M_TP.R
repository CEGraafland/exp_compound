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
sel.tp.1<- TimeCoordsAnom_from_Grid_rms(sel.tp,rms = TRUE)
sel.t2m.1<-TimeCoordsAnom_from_Grid_rms(sel.t2m,rms = TRUE)
df.tp <-as.data.frame(sel.tp.1)
df.t2m <-as.data.frame(sel.t2m.1)
############################################
# discretize in 3
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

# sel.tp <- tp_ERA5_monthly_1940_2022_10d
# sel.t2m <- t2m_ERA5_monthly_1940_2022_10d
# tp_ERA5_monthly_1940_2022_10d <- NULL
# t2m_ERA5_monthly_1940_2022_10d <- NULL
# sel.tp.1<- TimeCoordsAnom_from_Grid_rms(sel.tp,rms = TRUE)
# sel.t2m.1<-TimeCoordsAnom_from_Grid_rms(sel.t2m,rms = TRUE)
# df.tp <-as.data.frame(sel.tp.1)
# df.t2m <-as.data.frame(sel.t2m.1)
# # sel.tp<- NULL
# # sel.t2m <- NULL
# ############################################
# # discretize in 4
# ############################################
# b <- 3
# disc.df.tp <- discretize(df.tp, method = "quantile", breaks = b)
# disc.df.t2m <- discretize(df.t2m, method = "quantile", breaks = b)
# names(disc.df.t2m) <- paste0(names(disc.df.t2m),".t2m")
# names(disc.df.tp) <- paste0(names(disc.df.tp),".tp")
# disc.df.t2m.tp <- cbind(disc.df.t2m,disc.df.tp)
############################################################################
# Function to load hciterations of models in a list 
############################################################################
#permused <- 0
#algo <- "tabu"
#score <- "aic"


loadIterationsComp <- function(permused = 0, algo = c("hc","tabu"),score = c("aic","bic"), it = NULL, limit = FALSE,compound = NULL,artificial = FALSE) {
  if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
  if(is.null(compound)) {comp <-""} else if(compound =="comp"){comp<- "_comp"} else if(!is.null(compound)){comp <- paste0("_comp_",compound)} 
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
tabu_aic_0 <- loadIterationsComp(
  permused = 0,
  algo = Algo, score = Score,
  limit = limit, compound = NULL, it = IT)


# loadIterationsComp <- function(permused = 0, algo = c("hc","tabu"),score = c("aic","bic"), it = NULL) {
#   
#   hc_list <- list.files(paste0("data/",algo,"iterations/disc_4_",score,"_t2m_tp_10d/perm",permused), full.names = T)
#   hc_names <- list.files(paste0("data/",algo,"iterations/disc_4_",score,"_t2m_tp_10d/perm",permused))
#   hc_names <- gsub(".rda", "", hc_names)
#   
#   if(!is.null(it)){
#     hc_list <- hc_list[grep(it,hc_list)]
#     hc_names <- hc_names[grep(it,hc_names)]
#   }
#   
#   hc_networks <- lapply(hc_list, function(x){get(load(x))})
#   names(hc_networks) <- hc_names
#   sizes <- sapply(hc_networks,narcs)
#   hc_networks <- hc_networks[order(sizes)]
#   return(hc_networks)
# }
# 
# Permused <- 0
# Algo <- "tabu"
# Score <- "aic"
# IT <- "0_2600_2700"
# #tabu_aic_0 <- loadIterationsComp(permused = 0,algo = "tabu", score = "aic",it = "0_1600_1700")
# tabu_aic_0 <- loadIterationsComp(permused = Permused,algo = Algo, score = Score,it = IT)
# 
# #tabu_bic_0 <- loadIterationsComp(permused = 0,algo = "tabu", score = "bic")
# hc_aic_0 <- loadIterationsComp(permused = 0,algo = "hc", score = "aic")
# hc_bic_0 <- loadIterationsComp(permused = 0,algo = "hc", score = "bic")


#tabu_aic_0_fits<- lapply(tabu_aic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))
tabu_aic_0_fits<- lapply(tabu_aic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))

#tabu_bic_0_fits<- lapply(tabu_bic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))
# hc_aic_0_fits<- lapply(hc_aic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))
# hc_bic_0_fits<- lapply(hc_bic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))

#################################################################################
# Single evidence. V81 postive (+ +)
#################################################################################
fitted <- tabu_aic_0_fits[[1]]
# cpquery(fitted, (V82.t2m == levels(disc.df.t2m.tp[,"V82.t2m"])[4]) & (V82.tp == levels(disc.df.t2m.tp[,"V82.tp"])[1]), evidence = 
#           (V82.t2m == levels(disc.df.t2m.tp[,"V82.t2m"])[4]),method = "ls",n = 100000)
# cpquery(fitted, (V82.t2m == levels(disc.df.t2m.tp[,"V82.t2m"])[4]) & (V82.tp == levels(disc.df.t2m.tp[,"V82.tp"])[1]), evidence = TRUE, method = "ls")
# cpquery(fitted, (V82.t2m == levels(disc.df.t2m.tp[,"V82.t2m"])[4]), evidence = TRUE)*cpquery(fitted,  (V82.tp == levels(disc.df.t2m.tp[,"V82.tp"])[1]), evidence = TRUE, method = "lw")


# data.c <- disc.df.t2m.tp
# perm <- NULL
# baysnet <- fitted
#ngrids<- length(nodes(baysnet))/2
# nodeComp <- 20
#nodesEvents <- c(nodeComp,ngrids+nodeComp)

#level.t2m <- 4
#level.tp <- 1

# FOR LS: 
# nodesEvidence <- c(81)
# level.ev <- c(4)

# nodesEvidence <- c(81,648 + 81)
# level.ev <- c(4,1)


PropagationCompSimple <- function(baysnet, data.c, nodeComp, level.t2m,level.tp, nodesEvidence, level.ev, perm){
  ngrids<- length(nodes(baysnet))/2
  nodesEvents <- c(nodeComp,ngrids+nodeComp)
  
  envlist <- list(baysnet = baysnet,
                  data.c = data.c,
                  nodeComp = nodeComp,
                  level.t2m = level.t2m,
                  level.tp = level.tp,
                  nodesEvidence = nodesEvidence,
                  level.ev = level.ev,
                  perm = perm,
                  ngrids = ngrids,
                  nodesEvents = nodesEvents)
  list2env(envlist, envir = parent.frame())
  
  # valueEvents <- c(quote(levels(data.c[,names(baysnet)[nodeComp]])[level.t2m]),
  #quote(levels(data.c[,names(baysnet)[ngrids+nodeComp]])[level.tp]))
  valueEvents <- c(paste0("(levels(data.c[,names(baysnet)[nodeComp]])[level.t2m])"),
                   paste0("(levels(data.c[,names(baysnet)[ngrids+nodeComp]])[level.tp])"))
  
  if(length(nodesEvidence)==1){
    #valueEvidence <- c(quote(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]]))
    valueEvidence <- c(paste0("(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]])"))
  } else if (length(nodesEvidence)>1){
    # valueEvidence <- c(quote(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]]),
    #          quote(levels(data.c[,names(baysnet)[nodesEvidence[2]]])[level.ev[2]]))
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
    probname <- paste0("P(",names(baysnet)[nodeComp],"=", level.t2m," & ",names(baysnet)[ngrids +nodeComp]," = ", level.tp  ,"|",names(baysnet)[nodesEvidence[1]]," = ", level.ev[1], ")")
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
    
    
    probname <- paste0("P(",names(baysnet)[nodeComp],"=", level.t2m," & ",names(baysnet)[ngrids +nodeComp]," = ", level.tp  ,"|")
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
  
  with <- eval(parse(text = cmd1))
  without <- eval(parse(text = cmd3))
  
  attr(with, "probability") <- probname
  attr(without, "probability") <- paste0("P(",names(baysnet)[nodeComp],"=", level.t2m," & ",names(baysnet)[ngrids +nodeComp]," = ", level.tp ,")")
  df <- data.frame(names = paste0(names(baysnet)[nodesEventsRef], collapse = "&"), with = with, without = without)
  return(df)
  
}



dfs <- list()
lvl.t2m <- b
lvl.tp <- b
# nE <- c(81,81+length(nodes(fitted))/2)
# lvl.ev <- c(b,b)
nE <- c(57)
lvl.ev <- c(b)

i <-1
for(i in 1:(length(nodes(fitted))/2)){
dfs[[i]]<- PropagationCompSimple(baysnet = fitted,
                      nodeComp = i,
                      data.c = disc.df.t2m.tp,
                      level.tp = lvl.tp,
                      level.t2m = lvl.t2m,
                      nodesEvidence = nE,
                      level.ev = lvl.ev,
                      perm = NULL)
}

if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
if(length(lvl.ev) ==1){
assign(paste0("prop_",Algo,"_",Score,"_",IT,lim,"_t2m",lvl.t2m,"_tp",lvl.tp,"_t2m",nE,"_",lvl.ev),do.call("rbind",dfs))
save(list = paste0("prop_",Algo,"_",Score,"_",IT,lim,"_t2m",lvl.t2m,"_tp",lvl.tp,"_t2m",nE,"_",lvl.ev),file = paste0("results/propagation/",Algo,"_",Score,"_",IT,lim,"/prop_",Algo,"_",Score,"_",IT,lim,"_t2m",lvl.t2m,"_tp",lvl.tp,"_t2m",nE,"_",lvl.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,lim,"_t2m",lvl.t2m,"_tp",lvl.tp,"_t2m",nE[1],"_",lvl.ev[1],"_tp",nE[1],"_",lvl.ev[2]),do.call("rbind",dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,lim,"_t2m",lvl.t2m,"_tp",lvl.tp,"_t2m",nE[1],"_",lvl.ev[1],"_tp",nE[1],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_",Score,"_",IT,lim,"/prop_",Algo,"_",Score,"_",IT,"_t2m",lvl.t2m,"_tp",lvl.tp,"_t2m",nE[1],"_",lvl.ev[1],"_tp",nE[1],"_",lvl.ev[2],".rda"))
  
}
# WERKT NIET DOOR DE lapply
# lapply(1:648,PropagationCompSimple, baysnet = fitted,data.c = disc.df.t2m.tp,
#                                                                        level.t2m =4,
#                                                                        level.tp = 4,
#                                                                        nodesEvidence = c(81),
#                                                                        level.ev = c(4),
#                                                                        perm = NULL)
# 
#######################################################################################
#
#######################################################################################
PropagationCompoundElemental <- function(baysnet, data.c, nodeComp, level.t2m,level.tp, nodesEvidence, level.ev, perm){
  ngrids<- length(nodes(baysnet))/3
  nodesEvents <- c(2*ngrids+nodeComp)
  
  envlist <- list(baysnet = baysnet,
                  data.c = data.c,
                  nodeComp = nodeComp,
                  level.t2m = level.t2m,
                  level.tp = level.tp,
                  nodesEvidence = nodesEvidence,
                  level.ev = level.ev,
                  perm = perm,
                  ngrids = ngrids,
                  nodesEvents = nodesEvents)
  list2env(envlist, envir = parent.frame())
  
  # valueEvents <- c(quote(levels(data.c[,names(baysnet)[nodeComp]])[level.t2m]),
  #quote(levels(data.c[,names(baysnet)[ngrids+nodeComp]])[level.tp]))
  valueEvents <- c(paste0("(levels(data.c[,names(baysnet)[nodeComp]])[level.t2m])"),
                   paste0("(levels(data.c[,names(baysnet)[ngrids+nodeComp]])[level.tp])"))
  
  if(length(nodesEvidence)==1){
    #valueEvidence <- c(quote(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]]))
    valueEvidence <- c(paste0("(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]])"))
  } else if (length(nodesEvidence)>1){
    # valueEvidence <- c(quote(levels(data.c[,names(baysnet)[nodesEvidence[1]]])[level.ev[1]]),
    #          quote(levels(data.c[,names(baysnet)[nodesEvidence[2]]])[level.ev[2]]))
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
    probname <- paste0("P(",names(baysnet)[nodeComp],"= HD|",names(baysnet)[nodesEvidence[1]]," = ", level.ev[1], ")")
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
    
    
    probname <- paste0("P(",names(baysnet)[nodeComp],"= HD|")
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
  
  with <- eval(parse(text = cmd1))
  without <- eval(parse(text = cmd3))
  
  attr(with, "probability") <- probname
  attr(without, "probability") <- paste0("P(",names(baysnet)[nodeComp],"= HD)")
  df <- data.frame(names = paste0(names(baysnet)[nodesEventsRef], collapse = "&"), with = with, without = without)
  return(df)
  
}

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
if(isTRUE(limit)){
  sel.tp <- subsetGrid(tp_ERA5_monthly_1940_2022_10d,latLim = c(-50,75))
  sel.t2m <- subsetGrid(t2m_ERA5_monthly_1940_2022_10d,latLim = c(-50,75))
} else {
  sel.tp <- tp_ERA5_monthly_1940_2022_10d
  sel.t2m <- t2m_ERA5_monthly_1940_2022_10d
}
tp_ERA5_monthly_1940_2022_10d <- NULL
t2m_ERA5_monthly_1940_2022_10d <- NULL
sel.tp.1<- TimeCoordsAnom_from_Grid_rms(sel.tp,rms = TRUE)
sel.t2m.1<-TimeCoordsAnom_from_Grid_rms(sel.t2m,rms = TRUE)
df.tp <-as.data.frame(sel.tp.1)
df.t2m <-as.data.frame(sel.t2m.1)
############################################
# discretize in 3
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

compound <- "HD"
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
    disc.df.t2m.tp <- compound.data(disc.df.t2m.tp,"HD",i)
  }
}



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
Score <- "aic"
IT<- "0_2700_2800"
tabu_bic_0_comp_art <- loadIterationsComp(
  permused = 0,algo = Algo, score = Score,
  limit = limit,compound = "comp",artificial = TRUE, it = IT)


#################################################################################
# Single evidence. V81 postive (+ +)
#################################################################################
fittedbn <- bn.fit(tabu_bic_0_comp_art[[1]],disc.df.t2m.tp)

#cpquery(fitted, (C82 == levels(disc.df.t2m.tp[,"C82"])[1]) & (V82.tp == levels(disc.df.t2m.tp[,"V82.tp"])[1]), evidence = 
 #          (V82.t2m == levels(disc.df.t2m.tp[,"V82.t2m"])[4]),method = "ls",n = 100000)
# cpquery(fitted, (V82.t2m == levels(disc.df.t2m.tp[,"V82.t2m"])[4]) & (V82.tp == levels(disc.df.t2m.tp[,"V82.tp"])[1]), evidence = TRUE, method = "ls")
# cpquery(fitted, (V82.t2m == levels(disc.df.t2m.tp[,"V82.t2m"])[4]), evidence = TRUE)*cpquery(fitted,  (V82.tp == levels(disc.df.t2m.tp[,"V82.tp"])[1]), evidence = TRUE, method = "lw")
#fitted$C468$prob

#  data.c <- disc.df.t2m.tp
#  rm(data.c)
#  perm <- NULL
#  baysnet <- fitted
# ngrids<- length(bnlearn::nodes(baysnet))/3
# nodeComp <- 346
# level.Comp <- 1
#nodesEvents <- c(nodeComp,ngrids+nodeComp)

#level.t2m <- 4
#level.tp <- 1

# FOR LS:
#nodesEvidence <- c(346+2*ngrids)
#level.ev <- c(2)

# nodesEvidence <- c(81,648 + 81)
# level.ev <- c(4,1)


PropagationComp2Comp <- function(baysnet, data.c, nodeComp, level.Comp, nodesEvidence, level.ev, perm, add.without = TRUE){
  ngrids<- length(nodes(baysnet))/3
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


########################################################
# compound HD events impact on other compound HD events
########################################################
j<- 1
dfs <- list()
lvl.comp <- 2
nE <- c(j+2*ngrids)
lvl.ev <- c(2)


#  for(i in 1:ngrids){
#   dfs[[i]]<- PropagationComp2Comp(baysnet = fittedbn,
#                                    nodeComp = (2*ngrids)+i,
#                                    data.c = disc.df.t2m.tp,
#                                    level.Comp = lvl.comp,
#                                    nodesEvidence = nE,
#                                    level.ev = lvl.ev,
#                                    perm = NULL, add.without = FALSE)
# }
dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                nodeComp = (2*ngrids)+i,
                                                                data.c = disc.df.t2m.tp,
                                                                level.Comp = lvl.comp,
                                                                nodesEvidence = nE,
                                                                level.ev = lvl.ev,
                                                                perm = NULL, add.without = TRUE), mc.cores = 15)


if(length(lvl.ev) ==1){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
  
}


for(j in 359:ngrids){
dfs <- list()
lvl.comp <- 2
nE <- c(j+2*ngrids)
lvl.ev <- c(2)
i <- 1



start_time <- Sys.time()
dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                 nodeComp = (2*ngrids)+i,
                                                                 data.c = disc.df.t2m.tp,
                                                                 level.Comp = lvl.comp,
                                                                 nodesEvidence = nE,
                                                                 level.ev = lvl.ev,
                                                                 perm = NULL, add.without = FALSE), mc.cores = 15)
end_time <- Sys.time()
end_time - start_time

#dfs.df <- do.call(rbind.data.frame,dfs)
#dfs.df$with-dfs.df$without
# 
# dfs2.df <- do.call(rbind.data.frame,dfs2)
# dfs2.df$with-dfs2.df$without


if(length(lvl.ev) ==1){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
  
}

}

########################################################
# compound HD events impact on Temperature H events
########################################################
j<- 1
dfs <- list()
lvl.comp <- 3
nE <- c(j+2*ngrids)
lvl.ev <- c(2)


# for(i in 1:ngrids){
#   dfs[[i]]<- PropagationComp2Comp(baysnet = fittedbn,
#                                   nodeComp = i,
#                                   data.c = disc.df.t2m.tp,
#                                   level.Comp = lvl.comp,
#                                   nodesEvidence = nE,
#                                   level.ev = lvl.ev,
#                                   perm = NULL, add.without = TRUE)
# }
dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                nodeComp = i,
                                                                data.c = disc.df.t2m.tp,
                                                                level.Comp = lvl.comp,
                                                                nodesEvidence = nE,
                                                                level.ev = lvl.ev,
                                                                perm = NULL, add.without = TRUE), mc.cores = 15)


if(length(lvl.ev) ==1){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
  
}

j <- 2
for(j in 383:ngrids){
  dfs <- list()
  lvl.comp <- 3
  nE <- c(j+2*ngrids)
  lvl.ev <- c(2)

  start_time <- Sys.time()
  dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                  nodeComp = i,
                                                                  data.c = disc.df.t2m.tp,
                                                                  level.Comp = lvl.comp,
                                                                  nodesEvidence = nE,
                                                                  level.ev = lvl.ev,
                                                                  perm = NULL, add.without = FALSE), mc.cores = 15)
  end_time <- Sys.time()
  end_time - start_time

  #dfs.df <- do.call(rbind.data.frame,dfs)
  #dfs.df$with-dfs.df$without
  # 
  # dfs2.df <- do.call(rbind.data.frame,dfs2)
  # dfs2.df$with-dfs2.df$without
  
  
  if(length(lvl.ev) ==1){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
  } else if (length(lvl.ev) ==2){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
    
  }
  
}
########################################################
# compound HD events impact on Temperature C events
########################################################
j<- 1
dfs <- list()
lvl.comp <- 1
nE <- c(j+2*ngrids)
lvl.ev <- c(2)


# for(i in 1:ngrids){
#   dfs[[i]]<- PropagationComp2Comp(baysnet = fittedbn,
#                                   nodeComp = i,
#                                   data.c = disc.df.t2m.tp,
#                                   level.Comp = lvl.comp,
#                                   nodesEvidence = nE,
#                                   level.ev = lvl.ev,
#                                   perm = NULL, add.without = TRUE)
# }
dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                nodeComp = i,
                                                                data.c = disc.df.t2m.tp,
                                                                level.Comp = lvl.comp,
                                                                nodesEvidence = nE,
                                                                level.ev = lvl.ev,
                                                                perm = NULL, add.without = TRUE), mc.cores = 15)


if(length(lvl.ev) ==1){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
  
}

for(j in 2:ngrids){
  dfs <- list()
  lvl.comp <- 1
  nE <- c(j+2*ngrids)
  lvl.ev <- c(2)
  
  # start_time <- Sys.time()
  dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                  nodeComp = i,
                                                                  data.c = disc.df.t2m.tp,
                                                                  level.Comp = lvl.comp,
                                                                  nodesEvidence = nE,
                                                                  level.ev = lvl.ev,
                                                                  perm = NULL, add.without = FALSE), mc.cores = 15)
  # end_time <- Sys.time()
  # end_time - start_time
  
  #dfs.df <- do.call(rbind.data.frame,dfs)
  #dfs.df$with-dfs.df$without
  # 
  # dfs2.df <- do.call(rbind.data.frame,dfs2)
  # dfs2.df$with-dfs2.df$without
  
  
  if(length(lvl.ev) ==1){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
  } else if (length(lvl.ev) ==2){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
    
  }
  
}
########################################################
# compound HD events impact on RainFall D events
########################################################
j<- 1
dfs <- list()
lvl.comp <- 1
nE <- c(j+2*ngrids)
lvl.ev <- c(2)


# for(i in 1:ngrids){
#   dfs[[i]]<- PropagationComp2Comp(baysnet = fittedbn,
#                                   nodeComp = (ngrids +i),
#                                   data.c = disc.df.t2m.tp,
#                                   level.Comp = lvl.comp,
#                                   nodesEvidence = nE,
#                                   level.ev = lvl.ev,
#                                   perm = NULL, add.without = TRUE)
# }
dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                nodeComp = (ngrids +i),
                                                                data.c = disc.df.t2m.tp,
                                                                level.Comp = lvl.comp,
                                                                nodesEvidence = nE,
                                                                level.ev = lvl.ev,
                                                                perm = NULL, add.without = TRUE), mc.cores = 15)



if(length(lvl.ev) ==1){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
  
}


for(j in 373:ngrids){
  dfs <- list()
  lvl.comp <- 1
  nE <- c(j+2*ngrids)
  lvl.ev <- c(2)
  
  start_time <- Sys.time()
  dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                  nodeComp = (ngrids +i),
                                                                  data.c = disc.df.t2m.tp,
                                                                  level.Comp = lvl.comp,
                                                                  nodesEvidence = nE,
                                                                  level.ev = lvl.ev,
                                                                  perm = NULL, add.without = FALSE), mc.cores = 15)
  end_time <- Sys.time()
  end_time - start_time
  
  #dfs.df <- do.call(rbind.data.frame,dfs)
  #dfs.df$with-dfs.df$without
  # 
  # dfs2.df <- do.call(rbind.data.frame,dfs2)
  # dfs2.df$with-dfs2.df$without
  
  
  if(length(lvl.ev) ==1){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
  } else if (length(lvl.ev) ==2){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
    
  }
  
}
########################################################
# compound HD events impact on RainFall W events
########################################################
j<- 1
dfs <- list()
lvl.comp <- 3
nE <- c(j+2*ngrids)
lvl.ev <- c(2)


# for(i in 1:ngrids){
#   dfs[[i]]<- PropagationComp2Comp(baysnet = fittedbn,
#                                   nodeComp = (ngrids +i),
#                                   data.c = disc.df.t2m.tp,
#                                   level.Comp = lvl.comp,
#                                   nodesEvidence = nE,
#                                   level.ev = lvl.ev,
#                                   perm = NULL, add.without = TRUE)
# }
dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                nodeComp = (ngrids +i),
                                                                data.c = disc.df.t2m.tp,
                                                                level.Comp = lvl.comp,
                                                                nodesEvidence = nE,
                                                                level.ev = lvl.ev,
                                                                perm = NULL, add.without = TRUE), mc.cores = 15)



if(length(lvl.ev) ==1){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
  
}


for(j in 2:ngrids){
  dfs <- list()
  lvl.comp <- 3
  nE <- c(j+2*ngrids)
  lvl.ev <- c(2)
  
  start_time <- Sys.time()
  dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                  nodeComp = (ngrids +i),
                                                                  data.c = disc.df.t2m.tp,
                                                                  level.Comp = lvl.comp,
                                                                  nodesEvidence = nE,
                                                                  level.ev = lvl.ev,
                                                                  perm = NULL, add.without = FALSE), mc.cores = 15)
  end_time <- Sys.time()
  end_time - start_time
  
  #dfs.df <- do.call(rbind.data.frame,dfs)
  #dfs.df$with-dfs.df$without
  # 
  # dfs2.df <- do.call(rbind.data.frame,dfs2)
  # dfs2.df$with-dfs2.df$without
  
  
  if(length(lvl.ev) ==1){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
  } else if (length(lvl.ev) ==2){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
    
  }
  
}

########################################################
#  D events impact on other compound HD events
########################################################
j<- 1
dfs <- list()
lvl.comp <- 2
nE <- c(j+ngrids)
lvl.ev <- c(1)


#  for(i in 1:ngrids){
#   dfs[[i]]<- PropagationComp2Comp(baysnet = fittedbn,
#                                    nodeComp = (2*ngrids)+i,
#                                    data.c = disc.df.t2m.tp,
#                                    level.Comp = lvl.comp,
#                                    nodesEvidence = nE,
#                                    level.ev = lvl.ev,
#                                    perm = NULL, add.without = FALSE)
# }
dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                nodeComp = (2*ngrids)+i,
                                                                data.c = disc.df.t2m.tp,
                                                                level.Comp = lvl.comp,
                                                                nodesEvidence = nE,
                                                                level.ev = lvl.ev,
                                                                perm = NULL, add.without = TRUE), mc.cores = 15)


if(length(lvl.ev) ==1){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
  
}


for(j in 2:ngrids){
  dfs <- list()
  lvl.comp <- 2
  nE <- c(j+ngrids)
  lvl.ev <- c(1)

  start_time <- Sys.time()
  dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                  nodeComp = (2*ngrids)+i,
                                                                  data.c = disc.df.t2m.tp,
                                                                  level.Comp = lvl.comp,
                                                                  nodesEvidence = nE,
                                                                  level.ev = lvl.ev,
                                                                  perm = NULL, add.without = FALSE), mc.cores = 15)
  end_time <- Sys.time()
  end_time - start_time
  
  #dfs.df <- do.call(rbind.data.frame,dfs)
  #dfs.df$with-dfs.df$without
  # 
  # dfs2.df <- do.call(rbind.data.frame,dfs2)
  # dfs2.df$with-dfs2.df$without
  
  
  if(length(lvl.ev) ==1){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
  } else if (length(lvl.ev) ==2){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
    
  }
  
}

########################################################
#  H events impact on other compound HD events
########################################################
j<- 1
dfs <- list()
lvl.comp <- 2
nE <- c(j)
lvl.ev <- c(3)


#  for(i in 1:ngrids){
#   dfs[[i]]<- PropagationComp2Comp(baysnet = fittedbn,
#                                    nodeComp = (2*ngrids)+i,
#                                    data.c = disc.df.t2m.tp,
#                                    level.Comp = lvl.comp,
#                                    nodesEvidence = nE,
#                                    level.ev = lvl.ev,
#                                    perm = NULL, add.without = FALSE)
# }
dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                nodeComp = (2*ngrids)+i,
                                                                data.c = disc.df.t2m.tp,
                                                                level.Comp = lvl.comp,
                                                                nodesEvidence = nE,
                                                                level.ev = lvl.ev,
                                                                perm = NULL, add.without = TRUE), mc.cores = 15)


if(length(lvl.ev) ==1){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
  
}


for(j in 2:ngrids){
  dfs <- list()
  lvl.comp <- 2
  nE <- c(j)
  lvl.ev <- c(3)
  
  start_time <- Sys.time()
  dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                  nodeComp = (2*ngrids)+i,
                                                                  data.c = disc.df.t2m.tp,
                                                                  level.Comp = lvl.comp,
                                                                  nodesEvidence = nE,
                                                                  level.ev = lvl.ev,
                                                                  perm = NULL, add.without = FALSE), mc.cores = 15)
  end_time <- Sys.time()
  end_time - start_time
  
  #dfs.df <- do.call(rbind.data.frame,dfs)
  #dfs.df$with-dfs.df$without
  # 
  # dfs2.df <- do.call(rbind.data.frame,dfs2)
  # dfs2.df$with-dfs2.df$without
  
  
  if(length(lvl.ev) ==1){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
  } else if (length(lvl.ev) ==2){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
    
  }
  
}

########################################################
# D events impact on RainFall D events
########################################################
j<- 1
dfs <- list()
lvl.comp <- 1
nE <- c(j+ngrids)
lvl.ev <- c(1)


# for(i in 1:ngrids){
#   dfs[[i]]<- PropagationComp2Comp(baysnet = fittedbn,
#                                   nodeComp = (ngrids +i),
#                                   data.c = disc.df.t2m.tp,
#                                   level.Comp = lvl.comp,
#                                   nodesEvidence = nE,
#                                   level.ev = lvl.ev,
#                                   perm = NULL, add.without = TRUE)
# }
dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                nodeComp = (ngrids +i),
                                                                data.c = disc.df.t2m.tp,
                                                                level.Comp = lvl.comp,
                                                                nodesEvidence = nE,
                                                                level.ev = lvl.ev,
                                                                perm = NULL, add.without = TRUE), mc.cores = 15)



if(length(lvl.ev) ==1){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
  
}


for(j in 2:ngrids){
  dfs <- list()
  lvl.comp <- 1
  nE <- c(j+ngrids)
  lvl.ev <- c(1)
  
  start_time <- Sys.time()
  dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                  nodeComp = (ngrids +i),
                                                                  data.c = disc.df.t2m.tp,
                                                                  level.Comp = lvl.comp,
                                                                  nodesEvidence = nE,
                                                                  level.ev = lvl.ev,
                                                                  perm = NULL, add.without = FALSE), mc.cores = 15)
  end_time <- Sys.time()
  end_time - start_time
  
  #dfs.df <- do.call(rbind.data.frame,dfs)
  #dfs.df$with-dfs.df$without
  # 
  # dfs2.df <- do.call(rbind.data.frame,dfs2)
  # dfs2.df$with-dfs2.df$without
  
  
  if(length(lvl.ev) ==1){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
  } else if (length(lvl.ev) ==2){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
    
  }
  
}

########################################################
# H events impact on Temperature H events
########################################################
j<- 1
dfs <- list()
lvl.comp <- 3
nE <- c(j)
lvl.ev <- c(3)


# for(i in 1:ngrids){
#   dfs[[i]]<- PropagationComp2Comp(baysnet = fittedbn,
#                                   nodeComp = i,
#                                   data.c = disc.df.t2m.tp,
#                                   level.Comp = lvl.comp,
#                                   nodesEvidence = nE,
#                                   level.ev = lvl.ev,
#                                   perm = NULL, add.without = TRUE)
# }
dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                nodeComp = i,
                                                                data.c = disc.df.t2m.tp,
                                                                level.Comp = lvl.comp,
                                                                nodesEvidence = nE,
                                                                level.ev = lvl.ev,
                                                                perm = NULL, add.without = TRUE), mc.cores = 15)


if(length(lvl.ev) ==1){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
  
}


for(j in 2:ngrids){
  dfs <- list()
  lvl.comp <- 3
  nE <- c(j)
  lvl.ev <- c(3)
  
  start_time <- Sys.time()
  dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                  nodeComp = i,
                                                                  data.c = disc.df.t2m.tp,
                                                                  level.Comp = lvl.comp,
                                                                  nodesEvidence = nE,
                                                                  level.ev = lvl.ev,
                                                                  perm = NULL, add.without = FALSE), mc.cores = 15)
  end_time <- Sys.time()
  end_time - start_time
  
  #dfs.df <- do.call(rbind.data.frame,dfs)
  #dfs.df$with-dfs.df$without
  # 
  # dfs2.df <- do.call(rbind.data.frame,dfs2)
  # dfs2.df$with-dfs2.df$without
  
  
  if(length(lvl.ev) ==1){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
  } else if (length(lvl.ev) ==2){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
    
  }
  
}
########################################################
# D events impact on Temperature H events
########################################################
j<- 1
dfs <- list()
lvl.comp <- 3
nE <- c(j+ngrids)
lvl.ev <- c(1)


# for(i in 1:ngrids){
#   dfs[[i]]<- PropagationComp2Comp(baysnet = fittedbn,
#                                   nodeComp = i,
#                                   data.c = disc.df.t2m.tp,
#                                   level.Comp = lvl.comp,
#                                   nodesEvidence = nE,
#                                   level.ev = lvl.ev,
#                                   perm = NULL, add.without = TRUE)
# }
dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                nodeComp = i,
                                                                data.c = disc.df.t2m.tp,
                                                                level.Comp = lvl.comp,
                                                                nodesEvidence = nE,
                                                                level.ev = lvl.ev,
                                                                perm = NULL, add.without = TRUE), mc.cores = 15)


if(length(lvl.ev) ==1){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
  
}


for(j in 2:ngrids){
  dfs <- list()
  lvl.comp <- 3
  nE <- c(j+ngrids)
  lvl.ev <- c(1)
  
  start_time <- Sys.time()
  dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                  nodeComp = i,
                                                                  data.c = disc.df.t2m.tp,
                                                                  level.Comp = lvl.comp,
                                                                  nodesEvidence = nE,
                                                                  level.ev = lvl.ev,
                                                                  perm = NULL, add.without = FALSE), mc.cores = 15)
  end_time <- Sys.time()
  end_time - start_time
  
  #dfs.df <- do.call(rbind.data.frame,dfs)
  #dfs.df$with-dfs.df$without
  # 
  # dfs2.df <- do.call(rbind.data.frame,dfs2)
  # dfs2.df$with-dfs2.df$without
  
  
  if(length(lvl.ev) ==1){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
  } else if (length(lvl.ev) ==2){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
    
  }
  
}

########################################################
# H events impact on RainFall D events
########################################################
j<- 1
dfs <- list()
lvl.comp <- 1
nE <- c(j)
lvl.ev <- c(3)


# for(i in 1:ngrids){
#   dfs[[i]]<- PropagationComp2Comp(baysnet = fittedbn,
#                                   nodeComp = (ngrids +i),
#                                   data.c = disc.df.t2m.tp,
#                                   level.Comp = lvl.comp,
#                                   nodesEvidence = nE,
#                                   level.ev = lvl.ev,
#                                   perm = NULL, add.without = TRUE)
# }
dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                nodeComp = (ngrids +i),
                                                                data.c = disc.df.t2m.tp,
                                                                level.Comp = lvl.comp,
                                                                nodesEvidence = nE,
                                                                level.ev = lvl.ev,
                                                                perm = NULL, add.without = TRUE), mc.cores = 15)



if(length(lvl.ev) ==1){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
} else if (length(lvl.ev) ==2){
  assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
  save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
  
}


for(j in 2:ngrids){
  dfs <- list()
  lvl.comp <- 1
  nE <- c(j)
  lvl.ev <- c(3)
  
  start_time <- Sys.time()
  dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                  nodeComp = (ngrids +i),
                                                                  data.c = disc.df.t2m.tp,
                                                                  level.Comp = lvl.comp,
                                                                  nodesEvidence = nE,
                                                                  level.ev = lvl.ev,
                                                                  perm = NULL, add.without = FALSE), mc.cores = 15)
  end_time <- Sys.time()
  end_time - start_time
  
  #dfs.df <- do.call(rbind.data.frame,dfs)
  #dfs.df$with-dfs.df$without
  # 
  # dfs2.df <- do.call(rbind.data.frame,dfs2)
  # dfs2.df$with-dfs2.df$without
  
  
  if(length(lvl.ev) ==1){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
  } else if (length(lvl.ev) ==2){
    assign(paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
    
  }
  
}