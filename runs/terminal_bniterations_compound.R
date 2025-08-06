########################################################
rm(list = ls())
setwd("~/data/Untitled/Trabajo/R_practice/exp_compound")
library(bnlearn)
library(transformeR)
library(magrittr)
library(igraph)
source("../R/Functions/BasicNetworkFunctions.R")
######################################################
# 
######################################################
degrees <- "5d"
mask <- TRUE # Only land or land & sea
limit <- TRUE # include poles?
b <- 3 # How many levels TP & T2M:
k <- 0 # Which permutation of variables later
compound <- "HD" # Which type of compound Hot Dry, Cold Dry, Hot Wet, Cold Wet, all together, all.distinct
global_mean_temp <- FALSE
global_mean_temp_cont <- FALSE# Extra node for global_mean_temp?
combi <- "only" # only Compound nodes or also "TPT2M" nodes
algo<- "tabu" # Which algorithm
score<- "aic" # Which score in algorithm

if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
if(!is.null(compound)){comp <- paste0("_comp_",compound)} else {comp <-""}
if(isTRUE(mask)){superf <- "_land"}else {superf <-"_landsea"}
if(isTRUE(global_mean_temp)){n1 <- "_GMT"} else if(isTRUE(global_mean_temp_cont)){n1 <- "_GMTcont"} else{n1 <- ""}
if (isTRUE(global_mean_temp_cont) & score == "bic"){score <- "bic-cg"} 
if (isTRUE(global_mean_temp_cont) & score == "aic"){score <- "aic-cg"} 
if (isTRUE(global_mean_temp_cont) & score == "loglik"){score <- "loglik-cg"} 

if(!k==0){permutations <-get(load(paste0("data/permutations/permutations_",degrees,lim,superf,".rda")))}

data <- get(load(file =paste0("data/manipulate_aggr/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".rda")))
if(!k==0){data <-data[,c(permutations[[k]],(length(permutations[[k]])+1):ncol(data))]}
#################################################################
# Create permutations
#################################################################
# permlength <- length(grep("C",colnames(data)))
# permutations <- list()
# for(i in 1:25){
#   permutations[[i]]<- sample(1:permlength,permlength)
# }
# assign(paste0("permutations_",degrees,lim,superf),permutations)
#save(list = paste0("permutations_",degrees,lim,superf),file = paste0("data/permutations/permutations_",degrees,lim,superf,".rda"))

# ##########################
# 
# permutations <- get(load(file = paste0("data/permutations/permutations_",degrees,lim,superf,".rda")))
# 
# backpermutations <- lapply(permutations, function (x)order(((1:length(x))[x])))
# 
# backpermutations[[1]] <- order(((1:length(permutations[[1]]))[permutations[[1]]]))
# 
# ((1:length(permutations[[1]]))[permutations[[1]]])[backpermutations[[1]]]
# 
# 
# 1:length(permutations[[1]])
# # mapply(FUN =function(x,y)x[y], x = permutations, y = backpermutations)
# assign(paste0("backpermutations_",degrees,lim,superf),backpermutations)
# save(list = paste0("backpermutations_",degrees,lim,superf),file = paste0("data/permutations/backpermutations_",degrees,lim,superf,".rda"))
# #####################################################################################
# create directories
#####################################################################################
 # dir.create(paste0("data/",algo,"iterations/"))
if(!dir.exists(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"/"))){dir.create(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"/"))}
if(!dir.exists(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"/perm",k))){ dir.create(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"/perm",k))}
if(!dir.exists(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"/perm",k,"train"))){ dir.create(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"/perm",k,"train"))}
if(!dir.exists(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"/perm",k,"test"))){ dir.create(paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"/perm",k,"test"))}
##########################################################################
# Full data 
###########################################################################
start <- NULL
steps <- 100
last <- 10000
m<- 0
for (m in 0:(last/steps)) {
  #for (m in 0:0) {
  i <- m*steps
  j <- i+steps
  
  if (algo == "hc") {
    berekening <- hc(data, max.iter = steps, score = score, start = start)
  } else if (algo == "tabu"){
    berekening <- tabu(data, max.iter = steps, score = score, start = start)
  }
  assign(paste0(algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_",k,"_",i,"_",j,"i"), berekening)
 
  save(list = paste0(algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_",k,"_",i,"_",j,"i"),
       file = paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"/perm",k,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,"_",k,"_",i,"_",j,"i.rda"))
  
  if(m==0){
    start <- berekening
  } else if(narcs(berekening) == narcs(start)){
    break
  } else {start <- berekening}
}

##########################################################################
# Train data 
##########################################################################
# set.seed(90)
# samplesize <- nrow(data)
# indTRAIN <- sample(1:samplesize,samplesize/2)
# indTEST <- (1:samplesize)[-indTRAIN]
# assign(paste0("indTRAIN1"),indTRAIN)
# save(list = paste0("indTRAIN1"),file = paste0("data/trainindices/ERA5_monthly_1940_2022_indTRAIN1.rda"))
###########################################################################
start <- NULL
steps <- 100
last <- 10000

load(paste0("data/trainindices/ERA5_monthly_1940_2022_indTRAIN1.rda"))
learndata <- data[indTRAIN1,]

for (m in 0:(last/steps)) {
  #for (m in 0:0) {
  i <- m*steps
  j <- i+steps
  
  if (algo == "hc") {
    berekening <- hc(learndata, max.iter = steps, score = score, start = start)
  } else if (algo == "tabu"){
    berekening <- tabu(learndata, max.iter = steps, score = score, start = start)
  }
  assign(paste0("train1_",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_",k,"_",i,"_",j,"i"), berekening)
  
  save(list = paste0("train1_",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"_",k,"_",i,"_",j,"i"),
       file = paste0("data/",algo,"iterations/ERA5_monthly_1940_2022_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,"/perm",k,"train/train1_",algo,"_disc_",b,"_",score,lim,comp,combi,superf,"_",k,"_",i,"_",j,"i.rda"))
  
  if(m==0){
    start <- berekening
  } else if(narcs(berekening) == narcs(start)){
    break
  } else {start <- berekening}
}