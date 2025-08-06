##################################################################################################
# Propagation of evidence TP T2m C GT 
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

degrees <- "10d"
mask <- FALSE # Only land or land & sea
limit <- TRUE # include poles?
b <- 3 # How many levels TP & T2M:
k <- 0 # Which permutation of variables later
IT <- paste0(k,"_2700_2800")
compound <- "HD" # Which type of compound Hot Dry, Cold Dry, Hot Wet, Cold Wet, all together, all.distinct
global_mean_temp <- TRUE
global_mean_temp_cont <- FALSE # Extra node for global_mean_temp?
combi <- "TPT2M" # only Compound nodes or also "TPT2M" nodes
algo<- "tabu" # Which algorithm
score<- "aic" # Which score in algorithm
global_mean_temp_art <- TRUE
cut.art <- Inf
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

if(isTRUE(global_mean_temp_art)){
tabu_bic_0 <- loadIterationsComp(permused = k,algo = algo, score = score,it = IT)
sel.red <- grep("GMT",names(tabu_bic_0))
tabu_bic_0 <- tabu_bic_0[sel.red]}
# } else {
#   tabu_bic_0 <- loadIterationsComp(permused = k,algo = algo, score = score,it = IT)
#   sel.red <- grep("GMT",names(tabu_bic_0))
#   tabu_bic_0 <- tabu_bic_0[sel.red]
# }
tabu_bic_0_fits<- lapply(tabu_bic_0, function(x) bn.fit(x = x, data = data))
fittedbn <- tabu_bic_0_fits[[1]]
#####################################################################################################
#
#####################################################################################################

PropagationComp2Comp <- function(baysnet, data.c, nodeComp, level.Comp, nodesEvidence, level.ev, perm, add.without = TRUE){
  ngrids<- (length(nodes(baysnet))-1)/3
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
# high GT impact on compound HD events
########################################################
# ngrids <- (length(nodes(fittedbn))-1)/3
# j<- 1
# dfs <- list()
# lvl.comp <- 2
# nE <- c(j+3*ngrids)
# lvl.ev <- c(4)
# 
# i <- 1
# for(i in 1:ngrids){
#   dfs[[i]]<-PropagationComp2Comp(baysnet = fittedbn,
#                                   nodeComp = (2*ngrids)+i,
#                                   data.c = data,
#                                   level.Comp = lvl.comp,
#                                   nodesEvidence = nE,
#                                   level.ev = lvl.ev,
#                                   perm = NULL, add.without = TRUE)
# }
# 
# 
# if(length(lvl.ev) ==1){
#   assign(paste0("prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
#   save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
# } else if (length(lvl.ev) ==2){
#   assign(paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
#   save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
#   
# }


  
  ngrids <- (length(nodes(fittedbn))-1)/3
  j<- 1
  dfs <- list()
  lvl.comp <- 2
  nE <- c(j+3*ngrids)
  lvl.ev <- c(1)
  
  
  start_time <- Sys.time()
  dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                  nodeComp = (2*ngrids)+i,
                                                                  data.c = data,
                                                                  level.Comp = lvl.comp,
                                                                  nodesEvidence = nE,
                                                                  level.ev = lvl.ev,
                                                                  perm = NULL, add.without = TRUE), mc.cores = 13)
  end_time <- Sys.time()
  end_time - start_time
  
  #dfs.df <- do.call(rbind.data.frame,dfs)
  #dfs.df$with-dfs.df$without
  # 
  # dfs2.df <- do.call(rbind.data.frame,dfs2)
  # dfs2.df$with-dfs2.df$without
  
  fittedbn$GT$parents
  
  if(!dir.exists(paste0("results/propagation/",degrees,"/"))){dir.create(paste0("results/propagation/",degrees,"/"))}
  if(!dir.exists(paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/"))){dir.create(paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/"))}
  
  if(length(lvl.ev) ==1){
    assign(paste0("prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
  } else if (length(lvl.ev) ==2){
    assign(paste0("prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
    
  }
  
  load(file = paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
  prop_tabu_aic_0_2700_2800_C_2_GT_4
 proplist <-list(prop_tabu_aic_0_2700_2800_C_2_GT_4)
load("data/raw_aggr/t2m_ERA5_monthly_1940_2022_10d.rda")


limit <- TRUE
if(isTRUE(limit)){

  sel.t2m <- subsetGrid(t2m_ERA5_monthly_1940_2022_10d,latLim = c(-50,75))
} else {
  sel.tp <- tp_ERA5_monthly_1940_2022_10d
  sel.t2m <- t2m_ERA5_monthly_1940_2022_10d
}


col.r <- colorRampPalette(brewer.pal(9, "Reds"))(15)
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(15)
col.rb <- c(rev(col.b),'white','white',col.r)


#proplist2 <- proplist
proplist <- lapply(1:length(proplist), function(x) {proplist[[x]]$without <- proplist[[1]]$without ; return(proplist[[x]])})


#propdiflist <- lapply(1:length(proplist), function(x) {proplist[[x]]$without <- proplist[[1]]$without ; return(proplist[[x]])})

climlist <- lapply(proplist, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlist <- lapply(proplist, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
#climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = sel.t2m))
withoutvector <- proplist[[1]]$without
withoutclim <- quantity2clim(withoutvector, paste0(attr(withoutvector, "probability")),ref.grid = sel.t2m)

# climlistt2m3 <- lapply(proplistt2m3, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
# climwithlistt2m3 <- lapply(proplistt2m3, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
# #climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = sel.t2m))
# withoutvectort2m3 <- proplistt2m3[[1]]$without
# withoutclimt2m3 <- quantity2clim(withoutvectort2m3, paste0(attr(withoutvectort2m3, "probability")),ref.grid = sel.t2m)
#
# climlisttp1 <- lapply(proplisttp1, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
# climwithlisttp1 <- lapply(proplisttp1, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
# #climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = sel.t2m))
# withoutvectortp1 <- proplisttp1[[1]]$without
# withoutclimtp1 <- quantity2clim(withoutvectortp1, paste0(attr(withoutvectortp1, "probability")),ref.grid = sel.t2m)
#
# climlistt2m1 <- lapply(proplistt2m1, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
# climwithlistt2m1 <- lapply(proplistt2m1, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
# #climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = sel.t2m))
# withoutvectort2m1 <- proplistt2m1[[1]]$without
# withoutclimt2m1 <- quantity2clim(withoutvectort2m1, paste0(attr(withoutvectort2m1, "probability")),ref.grid = sel.t2m)
#
# climlisttp3 <- lapply(proplisttp3, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
# climwithlisttp3 <- lapply(proplisttp3, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
# #climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = sel.t2m))
# withoutvectortp3 <- proplisttp3[[1]]$without
# withoutclimtp3 <- quantity2clim(withoutvectortp3, paste0(attr(withoutvectortp3, "probability")),ref.grid = sel.t2m)
#



enspropdifsclims <- bindGrid(climlist, dimension = c("member"))
enspropdifsclims$Members <- cnodes
enspropwithsclims <- bindGrid(climwithlist, dimension = c("member"))
enspropwithsclims$Members <- cnodes
#enspropwithoutsclims <- bindGrid(climwithoutlist, dimension = c("member"))
#enspropwithoutsclims$Members <- cnodes

#namesmembers <- apply(lvls,MARGIN = 1, function(x)paste0("lvl.t2m = ",x[1]," & lvl.tp = ",x[2]))
namesmembers <- cnodes

# plotname <- paste0("figs/evidence_propagation/propdif_",Algo,"_",Score,"_",IT,"_t2m_tp_4lvls_t2m_",nE,"_",lvl.ev,".pdf")
# pdf(plotname, width = 10, height = 7)
# plotname <- paste0("figs/evidence_propagation/propdif_",Algo,"_",Score,"_",IT,"_t2m_tp_4lvls_t2m_",nE,"_",lvl.ev,".png")
# png(plotname, width = 20, height = 15,units = "cm", res = 150)



spatialPlot(enspropdifsclims,as.table = TRUE,
            #names.attr = namesmembers,
            backdrop.theme = "coastline", lonCenter =0, main = list(paste0("P(C = HD & V.tp = lvl.tp|",names(fitted)[nE]," = ",lvl.ev,") - P(V.t2m = lvl.t2m & V.tp = lvl.tp)")),
              at = seq(-0.2,0.2,0.025),
            region = TRUE,
            col.regions= col.rb,
            set.max = 0.4,
            rev.colors = TRUE,
            colorkey = list(width = 0.6, lables = list(cex = 0.5)))



########################################################
# compound HD events impact on Temperature H events
########################################################


# # for(i in 1:ngrids){
# #   dfs[[i]]<- PropagationComp2Comp(baysnet = fittedbn,
# #                                   nodeComp = i,
# #                                   data.c = disc.df.t2m.tp,
# #                                   level.Comp = lvl.comp,
# #                                   nodesEvidence = nE,
# #                                   level.ev = lvl.ev,
# #                                   perm = NULL, add.without = TRUE)
# # }
# dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
#                                                                 nodeComp = i,
#                                                                 data.c = data,
#                                                                 level.Comp = lvl.comp,
#                                                                 nodesEvidence = nE,
#                                                                 level.ev = lvl.ev,
#                                                                 perm = NULL, add.without = TRUE), mc.cores = 15)
# 
# 
# if(length(lvl.ev) ==1){
#   assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
#   save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
# } else if (length(lvl.ev) ==2){
#   assign(paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
#   save(list = paste0("prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
#   
# }

  ngrids <- (length(nodes(fittedbn))-1)/3
  j<- 1
  dfs <- list()
  lvl.comp <- 3
  nE <- c(j+3*ngrids)
  lvl.ev <- c(4)
  
  start_time <- Sys.time()
  dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                  nodeComp = i,
                                                                  data.c = data,
                                                                  level.Comp = lvl.comp,
                                                                  nodesEvidence = nE,
                                                                  level.ev = lvl.ev,
                                                                  perm = NULL, add.without = FALSE), mc.cores = 15)
  end_time <- Sys.time()
  end_time - start_time

  
  if(!dir.exists(paste0("results/propagation/",degrees,"/"))){dir.create(paste0("results/propagation/",degrees,"/"))}
  if(!dir.exists(paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/"))){dir.create(paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/"))}
  
  if(length(lvl.ev) ==1){
    assign(paste0("prop_",algo,"_",score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",algo,"_",score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop_",algo,"_",score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
  } else if (length(lvl.ev) ==2){
    assign(paste0("prop_",algo,"_",score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",algo,"_",score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop_",algo,"_",score,"_",IT,"_t2m_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
    
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
ngrids <- (length(nodes(fittedbn))-1)/3
j<- 1
dfs <- list()
lvl.comp <- 1
nE <- c(j+3*ngrids)
lvl.ev <- c(4)

  
  start_time <- Sys.time()
  dfs <- mclapply(1:ngrids,FUN = function(i) PropagationComp2Comp(baysnet = fittedbn,
                                                                  nodeComp = (ngrids +i),
                                                                  data.c = data,
                                                                  level.Comp = lvl.comp,
                                                                  nodesEvidence = nE,
                                                                  level.ev = lvl.ev,
                                                                  perm = NULL, add.without = FALSE), mc.cores = 15)
  end_time <- Sys.time()
  end_time - start_time
  
  
  if(!dir.exists(paste0("results/propagation/",degrees,"/"))){dir.create(paste0("results/propagation/",degrees,"/"))}
  if(!dir.exists(paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/"))){dir.create(paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/"))}
  
  if(length(lvl.ev) ==1){
    assign(paste0("prop_",algo,"_",score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",algo,"_",score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev),file = paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop_",algo,"_",score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda"))
  } else if (length(lvl.ev) ==2){
    assign(paste0("prop_",algo,"_",score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),do.call(rbind.data.frame,dfs))
    save(list = paste0("prop_",algo,"_",score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2]),file = paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop_",algo,"_",score,"_",IT,"_tp_",lvl.comp,"_",names(fittedbn)[nE][1],"_",lvl.ev[1],"_",names(fittedbn)[nE][2],"_",lvl.ev[2],".rda"))
    
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




