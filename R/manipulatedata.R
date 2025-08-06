########################################################
rm(list = ls())
setwd("~/data/Untitled/Trabajo/R_practice/exp_compound")
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("../R/Functions/BasicNetworkFunctions.R")
####################################
# mask
####################################
degrees <- "10d" #5d
mask <- FALSE
if(isTRUE(mask)){
  t2m_ERA5 <- get(load(paste0("data/mask_aggr/t2m_land_ERA5_monthly_1940_2022_",degrees,".rda")))
  tp_ERA5 <- get(load(paste0("data/mask_aggr/tp_land_ERA5_monthly_1940_2022_",degrees,".rda")))
} else {
  t2m_ERA5 <- get(load(paste0("data/raw_aggr/t2m_ERA5_monthly_1940_2022_",degrees,".rda")))
  tp_ERA5 <- get(load(paste0("data/raw_aggr/tp_ERA5_monthly_1940_2022_",degrees,".rda")))
}
####################################
# limit
####################################
limit <- FALSE
if(isTRUE(limit)){
  sel.tp <- subsetGrid(tp_ERA5,latLim = c(-50,75))
  sel.t2m <- subsetGrid( t2m_ERA5,latLim = c(-50,75))
} else {
  sel.tp <- tp_ERA5
  sel.t2m <-  t2m_ERA5
}
tp_ERA5<- NULL
t2m_ERA5  <- NULL
df.tp <-as.data.frame(TimeCoordsAnom_from_Grid_rms(sel.tp,rms = TRUE,mask=mask))
df.t2m <-as.data.frame(TimeCoordsAnom_from_Grid_rms(sel.t2m,rms = TRUE,mask = mask))
sel.tp<- NULL
sel.t2m <- NULL


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
###########################################
# add compound nodes
###########################################
compound <- "HW"
combi <-"TPT2M"

if(!is.null(compound)){
  ngrids <- ncol(disc.df.t2m.tp)/2
  comp.nodes <- paste0("C",1:ngrids)
  
  compound.data <- function(data,comp.type,node){ 
    data.num <- as.data.frame(sapply(data, function(x) as.numeric(x), simplify = FALSE))
    c<- node
    data[[comp.nodes[c]]]<- numeric(nrow(data))
    
    if (comp.type == "all.distinct"){
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==1 & data.num[[paste0("V",c,".tp")]]==b] <- 13
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==1 & data.num[[paste0("V",c,".tp")]]==1] <- 11
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==b & data.num[[paste0("V",c,".tp")]]==1] <- 31
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==b & data.num[[paste0("V",c,".tp")]]==b] <- 33
    } else if (comp.type == "all"){
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==1 & data.num[[paste0("V",c,".tp")]]==b] <- 1
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==1 & data.num[[paste0("V",c,".tp")]]==1] <- 1
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==b & data.num[[paste0("V",c,".tp")]]==1] <- 1
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==b & data.num[[paste0("V",c,".tp")]]==b] <- 1
    } else if (comp.type == "HD"){   
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==b & data.num[[paste0("V",c,".tp")]]==1] <- 1
    } else if (comp.type == "HW"){
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==b & data.num[[paste0("V",c,".tp")]]==b] <- 1
    } else if (comp.type == "CW"){
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==1 & data.num[[paste0("V",c,".tp")]]==b] <- 1
    } else if (comp.type == "CD"){
      data[[comp.nodes[c]]][data.num[[paste0("V",c,".t2m")]]==1 & data.num[[paste0("V",c,".tp")]]==1] <- 1
    }
    data[[comp.nodes[c]]] <- as.factor(data[[comp.nodes[c]]])
    return(data)
  }
  
  
  for(i in 1:(ncol(disc.df.t2m.tp)/2)){
    disc.df.t2m.tp <- compound.data(disc.df.t2m.tp,compound,i)
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
###################################################
# add Global Mean Temperature node
###################################################
global_mean_temp <- FALSE
if(isTRUE(global_mean_temp)){
t2m_world <- get(load("data/raw_aggr/t2m_ERA5_monthly_1940_2022_5d.rda"))
t2m_global_world_mean <-aggregateGrid(t2m_world, aggr.spatial = list(FUN = "mean", na.rm = TRUE))
GT <-data.frame("GT" = as.vector(t2m_global_world_mean$Data))
GT <- discretize(GT, breaks = 4)
}

global_mean_temp_cont <- FALSE
if(isTRUE(global_mean_temp_cont)){
  t2m_world <- get(load("data/raw_aggr/t2m_ERA5_monthly_1940_2022_5d.rda"))
  t2m_global_world_mean <-aggregateGrid(t2m_world, aggr.spatial = list(FUN = "mean", na.rm = TRUE))
  GT <-data.frame("GT" = as.vector(t2m_global_world_mean$Data))
}
###########################################################################
# 
###########################################################################
if(!is.null(compound)& combi == "only"){
  data <- disc.df.t2m.tp[,grep("C",colnames(disc.df.t2m.tp))]
} else if (!is.null(compound)&combi == "TPT2M") {data <- disc.df.t2m.tp
} else {data <- disc.df.t2m.tp}

# data.num <- as.data.frame(sapply(data, function(x) as.numeric(x), simplify = FALSE))
#data <- df[permutations[[k]]]
if(isTRUE(global_mean_temp)){
  data<- cbind(data,GT)
}
if(isTRUE(global_mean_temp_cont)){
  data<- cbind(data,GT)
}

if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
if(!is.null(compound)){comp <- paste0("_comp_",compound)} else {comp <-""}
if(isTRUE(mask)){superf <- "_land"}else {superf <-"_landsea"}
if(isTRUE(global_mean_temp)){n1 <- "_GMT"} else if(isTRUE(global_mean_temp_cont)){n1 <- "_GMTcont"} else{n1 <- ""}


assign(paste0("ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1),data)
save(list =paste0("ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1),
     file =paste0("data/manipulate_aggr/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".rda"))


