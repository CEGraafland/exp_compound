########################################################
rm(list = ls())
setwd("~/data/Untitled/Trabajo/R_practice/exp_compound")
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("../R/Functions/BasicNetworkFunctions.R")
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
k <- 0
algo<- "tabu"
score<- "bic"
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

start <- NULL
steps <- 100
last <- 10000

for (m in 0:(last/steps)) {
  #for (m in 0:0) {
  i <- m*steps
  j <- i+steps

  if (algo == "hc") {
    berekening <- hc(data, max.iter = steps, score = score, start = start)
  } else if (algo == "tabu"){
    berekening <- tabu(data, max.iter = steps, score = score, start = start)
  }
  assign(paste0(algo,"_disc_",b,"_",score,"_t2m_tp_10d",lim,comp,"_",k,"_",i,"_",j,"i"), berekening)

  save(list = paste0(algo,"_disc_",b,"_",score,"_t2m_tp_10d",lim,comp,"_",k,"_",i,"_",j,"i"),
       file = paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_10d",lim,comp,"/perm",k,"/",algo,"_disc_",b,"_",score,"_t2m_tp_10d",lim,comp,"_",k,"_",i,"_",j,"i.rda"))

  if(m==0){
    start <- berekening
  } else if(narcs(berekening) == narcs(start)){
    break
  } else {start <- berekening}
}
#####################################################################################
#
#####################################################################################

tabu_disc_3_bic_t2m_tp_10d_lim_comp_HD_only_0_0_100i

#####################################################################################
# create directories pruebas compound
#####################################################################################
compound_nodes <- c(105,122,
                     371,290,354,372,
                     480,462,
                     81,
                     27,
                     532,
                     154,
                     388,386,
                     637,638)
compound_vnodes <-sapply(compound_nodes,function(x)paste0("V",x))

# for(v_i in compound_vnodes){
# dir.create(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_comp_",v_i,"_10d/"))
# dir.create(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_comp_",v_i,"_10d/perm",k))
# }
# dir.create(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_comp_",compound_vnodes[v_i],"_10d/perm",k,"train"))
#  dir.create(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_comp_",compound_vnodes[v_i],"_10d/perm",k,"test"))
###########################################################################
# Compound add one. 
###########################################################################
data.num <- as.data.frame(sapply(disc.df.t2m.tp, function(x) as.numeric(x), simplify = FALSE))

#v_i <- "V105"
for(h in 3:length(compound_vnodes)){
  v_i <- compound_vnodes[h]
  data <- disc.df.t2m.tp
  data[[paste0(v_i,".comp")]]<- numeric(nrow(data))
  assign(paste0("data.num.",v_i,".2state"),disc.df.t2m.tp)
  
  data[[paste0(v_i,".comp")]][data.num[[paste0(v_i,".t2m")]]==1 & data.num[[paste0(v_i,".tp")]]==3] <- 13
  data[[paste0(v_i,".comp")]][data.num[[paste0(v_i,".t2m")]]==1 & data.num[[paste0(v_i,".tp")]]==1] <- 11
  data[[paste0(v_i,".comp")]][data.num[[paste0(v_i,".t2m")]]==3 & data.num[[paste0(v_i,".tp")]]==1] <- 31
  data[[paste0(v_i,".comp")]][data.num[[paste0(v_i,".t2m")]]==3 & data.num[[paste0(v_i,".tp")]]==3] <- 33
  data[[paste0(v_i,".comp")]] <- as.factor(data[[paste0(v_i,".comp")]])
  
  # assign(paste0("data.num.",v_i,".5state"),data)
  # str(data)
  # str(data$V1.t2m)
  # disc.df.t2m.tp[1,2]
  #data <- df[permutations[[k]]]
  
  start <- NULL
  steps <- 2400
  last <- 2400
  
# for (m in 0:(last/steps)) {
for (m in 0:0) {
  i <- m*steps
  j <- i+steps
  
  if (algo == "hc") {
    berekening <- hc(data, max.iter = steps, score = score, start = start)
  } else if (algo == "tabu"){
    berekening <- tabu(data, max.iter = steps, score = score, start = start)
}
  assign(paste0(algo,"_disc_",b,"_",score,"_t2m_tp_comp_",v_i,"_10d_",k,"_",i,"_",j,"i"), berekening)
  
  save(list = paste0(algo,"_disc_",b,"_",score,"_t2m_tp_comp_",v_i,"_10d_",k,"_",i,"_",j,"i"), 
       file = paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_comp_",v_i,"_10d/perm",k,"/",algo,"_disc_",b,"_",score,"_t2m_tp_comp_",v_i,"_10d_",k,"_",i,"_",j,"i.rda"))
  
  if(m==0){
    start <- berekening
  } else if(narcs(berekening) == narcs(start)){
    break
  } else {start <- berekening}
}
}


########################
#
########################

b<- 3
algo
k <- 0
score <- "bic"
i <- 0
j <- 2400

h <- 4
for(h in 1:length(compound_vnodes)){
  v_i <- compound_vnodes[h]
netw <- get(load(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_comp_",v_i,"_10d/perm",k,"/",algo,"_disc_",b,"_",score,"_t2m_tp_comp_",v_i,"_10d_",k,"_",i,"_",j,"i.rda")))
netr <-remove.node(netw,paste0(v_i,".comp"))
all.equal(netr,tabu_disc_3_bic_t2m_tp_10d_0_2300_2400i)
}

compare(netr,tabu_disc_3_bic_t2m_tp_10d_0_2300_2400i)
netw$nodes$V290.comp
