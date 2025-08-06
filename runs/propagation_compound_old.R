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
###########################################################################################
# Create namelist with models + interim data
###########################################################################################
###########################################
# Standardize
###########################################
sel.tp <- tp_ERA5_monthly_1940_2022_10d
sel.t2m <- t2m_ERA5_monthly_1940_2022_10d
tp_ERA5_monthly_1940_2022_10d <- NULL
t2m_ERA5_monthly_1940_2022_10d <- NULL
sel.tp.1<- TimeCoordsAnom_from_Grid_rms(sel.tp,rms = TRUE)
sel.t2m.1<-TimeCoordsAnom_from_Grid_rms(sel.t2m,rms = TRUE)
df.tp <-as.data.frame(sel.tp.1)
df.t2m <-as.data.frame(sel.t2m.1)
# sel.tp<- NULL
# sel.t2m <- NULL
############################################
# discretize in 4
############################################
b <- 4
disc.df.tp <- discretize(df.tp, method = "quantile", breaks = b)
disc.df.t2m <- discretize(df.t2m, method = "quantile", breaks = b)
names(disc.df.t2m) <- paste0(names(disc.df.t2m),".t2m")
names(disc.df.tp) <- paste0(names(disc.df.tp),".tp")
disc.df.t2m.tp <- cbind(disc.df.t2m,disc.df.tp)
############################################################################
# Function to load hciterations of models in a list 
############################################################################

permused <- 0
algo <- "tabu"
score <- "aic"

loadIterationsComp <- function(permused = 0, algo = c("hc","tabu"),score = c("aic","bic"), it = NULL) {
  
  hc_list <- list.files(paste0("data/",algo,"iterations/disc_4_",score,"_t2m_tp_10d/perm",permused), full.names = T)
  hc_names <- list.files(paste0("data/",algo,"iterations/disc_4_",score,"_t2m_tp_10d/perm",permused))
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

tabu_aic_0 <- loadIterationsComp(permused = 0,algo = "tabu", score = "aic",it = "0_1600_1700")
tabu_bic_0 <- loadIterationsComp(permused = 0,algo = "tabu", score = "bic")
hc_aic_0 <- loadIterationsComp(permused = 0,algo = "hc", score = "aic")
hc_bic_0 <- loadIterationsComp(permused = 0,algo = "hc", score = "bic")


tabu_aic_0_fits<- lapply(tabu_aic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))
tabu_bic_0_fits<- lapply(tabu_bic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))
hc_aic_0_fits<- lapply(hc_aic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))
hc_bic_0_fits<- lapply(hc_bic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))

#################################################################################
# Single evidence. V81 postive (+ +)
#################################################################################
nodes(tabu_aic_0_fits$tabu_disc_4_aic_t2m_tp_10d_0_1600_1700i)
cpquery


detectCores()

levels(disc.df.t2m.tp[,"V1.t2m"])

## discrete Bayesian network (it is the same with ordinal nodes).
data(learning.test)
str(learning.test[1,2])
learning.test[1,2]
fitted = tabu_aic_0_fits$tabu_disc_4_aic_t2m_tp_10d_0_1600_1700i
var = names(disc.df.t2m.tp)
# the result should be around 0.025.
cpquery(fitted, (V82.t2m == levels(disc.df.t2m.tp[,"V82.t2m"])[4]) & (V82.tp == levels(disc.df.t2m.tp[,"V82.tp"])[1]), evidence = 
          (V82.t2m == levels(disc.df.t2m.tp[,"V82.t2m"])[4]),method = "ls",n = 100000)
cpquery(fitted, (V82.t2m == levels(disc.df.t2m.tp[,"V82.t2m"])[4]) & (V82.tp == levels(disc.df.t2m.tp[,"V82.tp"])[1]), evidence = TRUE, method = "ls")
cpquery(fitted, (V82.t2m == levels(disc.df.t2m.tp[,"V82.t2m"])[4]), evidence = TRUE)*cpquery(fitted,  (V82.tp == levels(disc.df.t2m.tp[,"V82.tp"])[1]), evidence = TRUE, method = "lw")

path.exists(tabu_aic_0$tabu_disc_4_aic_t2m_tp_10d_0_1600_1700i, "V82.tp", "V82.t2m", direct = TRUE, underlying.graph = TRUE, debug = FALSE)

fitted$V82.t2m
fitted$V82.tp
fitted$V64.t2m
fitted$V64.tp
fitted$V46.t2m
fitted$V46.tp

10000 * log10(nparams(fitted))
# programmatically build a conditional probability query...
var = names(learning.test)

obs = 2
str = paste("(", names(learning.test)[-3], " == '",
            sapply(learning.test[obs, -3]
                   , as.character), "')",
            sep = "", collapse = " & ")
str
str2 = paste("(", var[3], " == '",
             as.character(learning.test[obs, 3]), "')", sep = "")
str2

cmd = paste("cpquery(fitted, ", str2, ", ", str, ")", sep = "")
eval(parse(text = cmd))
# ... but note that predict works better in this particular case.
attr(predict(fitted, "C", learning.test[obs, -3], prob = TRUE), "prob")
# do the same with likelihood weighting.
cpquery(fitted, event = eval(parse(text = str2)),
        evidence = as.list(learning.test[2, -3]), method = "lw")
attr(predict(fitted, "C", learning.test[obs, -3],
             method = "bayes-lw", prob = TRUE), "prob")
# conditional distribution of A given C == "c".
table(cpdist(fitted, "A", (C == "c")))


PropagationExactGeneralPerm(baysnet = tabu_aic_0_fits$tabu_disc_4_aic_t2m_tp_10d_0_0_100i,
                            nodesEvents = 1:3,
                            valueEvent = ">= 1",
                            nodesEvidence = c(81),
                            valueEvidence = c(2),
                            perm = permutations[[whichperm]])

PropagationExactGeneralPerm
function(baysnet, nodesEvents, valueEvent, nodesEvidence, valueEvidence, perm){
  
  baysnet <- fitted
  nodes(baysnet)
  nodesEvents <- c(82,648+82)
  valueEvents <- c(quote(levels(disc.df.t2m.tp[,"V82.t2m"])[4]), quote(levels(disc.df.t2m.tp[,"V82.tp"])[1]))
  nodesEvidence <- c(81)
  valueEvidence <- c(quote(disc.df.t2m.tp[1,"V81.t2m"]))

  valueEvidence <- c(quote(as.list(disc.df.t2m.tp[1,])[81]))
  as.list(learning.test[1,-3])
  # dataperm <- datapermutations[[1]]
  # perm <- permutations[[1]]
  
  # baysnet = fitted
  # nodesEvents = c(298,299)
  # valueEvent = ">=1"
  # nodesEvidence = c(81,280)
  # valueEvidence = c(2,2)
  # perm = permutations[[1]]
  
  
  if (is.null(perm)) {nodesEventsRef <- nodesEvents} 
  else {
    nodesEventsRef <- c()
    for (i in 1:length(nodesEvents)){
      nodesEventsRef[i] <- which(perm == nodesEvents[i])
    }
  }
  
  with <- numeric(length = length(nodesEvents))
  without <- numeric(length = length(nodesEvents))
  
  
  if (length(nodesEvidence) == 1) {
    str2 <- paste0("(", names(baysnet)[nodesEvidence[1]],"==", valueEvidence[1], ")")
    str2 <- paste0("(", names(baysnet)[nodesEvidence[1]], "=", valueEvidence, ")")
    str2 <- valueEvidence
   # str2 <- paste0("(", names(baysnet)[nodesEvidence[1]]," == ", valueEvidence[1], ")")
   # probname <- paste0("P(V ", valueEvent,"|",nodesEvidence[1]," = ", valueEvidence[1], ")")
  }
  
  str(valueEvidence[[1]])
  
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
  
  # str2
  # i <- 2
  i <- 1
  as.list(learning.test[2, -3])
  for(i in 1:length(nodesEvents)) {
    # l <- nodesEvents[i]
    # l
    #l <- nodesEventsRef[i]
    # str <- paste0("(", names(baysnet)[l], ">=", valueEvent, ")")
    str <- paste0("(", names(baysnet)[nodesEvents],"==", valueEvents, ")", collapse = " & ")
    str
    nparams(baysnet)
    cmd1 = paste0("cpquery(baysnet, ", str, ", ", str2, ", method = ","'ls'",",n = 1000000)")
    cmd3 = paste0("cpquery(baysnet, ", str, ", ", "TRUE", ", method = ","'lw'",")")
    cmd1
    cmd3
    

    with[i] <- eval(parse(text = cmd1))
    with[i]
    
    sampling
    
    without[i] <- eval(parse(text = cmd3))
    without[i]
    
    
    # with[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = eval(parse(text = str2)))
    # withcomplement[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = eval(parse(text = str3)))
    # without[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = TRUE)
    
  }
  
  attr(with, "probability") <- probname
  attr(without, "probability") <- paste0("P(V ", valueEvent,")")
  df <- data.frame(names = names(baysnet)[nodesEventsRef], with = with, without = without)
  return(df)
  
}

##############################################################################################


x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(tabu_aic_0_fits)){
  # for (i in 8:length(selection_fits_gcms.future)){
  # for (i in 2:length(selection_fits_gcms.future)){
  
  assign(paste0("prop_",names(hc_gcms.future)[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2_detrended"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms.future[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("prop_",names(hc_gcms.future) [i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2_detrended"),
       file = paste0("results/propagation/perm",whichperm,"/posV81pos_detrended/CMIP5_rcp85_detrended/prop_",names(hc_gcms.future) [i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2_detrended.rda"))
  
}
#################################################################################
# Negative evidence. V81 (+ -)
#################################################################################
i <- 1
x <- seq(0,9000,100)
y <- seq(100,9100,100)
whichperm <- 3
for (i in 1:length(selection_fits_gcms.future)){
  #for (i in 1:length(selection_fits_gcms.future)){
  
  assign(paste0("propneg_",names(hc_gcms.future)[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2_detrended"),
         PropagationExactGeneralPerm(baysnet = selection_fits_gcms.future[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[whichperm]]))
  save(list = paste0("propneg_",names(hc_gcms.future)[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2_detrended"),
       file = paste0("results/propagation/perm",whichperm,"/posV81neg_detrended/CMIP5_rcp85_detrended/propneg_",names(hc_gcms.future)[i],"_hc_",x[modelsize],"_",y[modelsize],"i_V81_equal2_detrended.rda"))
  
}

