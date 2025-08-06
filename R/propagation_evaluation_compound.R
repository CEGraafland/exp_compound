#############################################################################
# Propgation of evidence
#############################################################################
setwd("~/data/Untitled/Trabajo/R_practice/exp_compound")
rm(list=ls())
library(bnlearn)
library(RColorBrewer)
library(magrittr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(visualizeR)
########################################################################################################
# Load functions and data
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/propagationFunctions.R")
source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
#load("data/raw_aggr/t2m_ERA5_monthly_1940_2022_10d.rda")
#load("data/raw_aggr/tp_ERA5_monthly_1940_2022_10d.rda")

degrees <- "10d"
mask <- FALSE # Only land or land & sea
limit <- TRUE # include poles?
b <- 3 # How many levels TP & T2M:
k <- 0 # Which permutation of variables later
IT <- paste0(k,"_2100_2200")
# IT <- paste0(k,"_1700_1800")
 IT <- paste0(k,"_4000_4100")
 IT <- paste0(k,"_2700_2800")
compound <- "HD" # Which type of compound Hot Dry, Cold Dry, Hot Wet, Cold Wet, all together, all.distinct
global_mean_temp <- TRUE
global_mean_temp_cont <- FALSE# Extra node for global_mean_temp?
combi <- "TPT2M" # only Compound nodes or also "TPT2M" nodes
algo<- "tabu" # Which algorithm
score<- "aic" # Which score in algorithm
global_mean_temp_art <- TRUE
cut.art <- Inf
scaleGT <- FALSE
geo <-"auto"

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

####################################
# mask
####################################
if(isTRUE(mask)){
  # t2m_ERA5 <- get(load(paste0("data/mask_aggr/t2m_land_ERA5_monthly_1940_2022_",degrees,".rda")))
  # tp_ERA5 <- get(load(paste0("data/mask_aggr/tp_land_ERA5_monthly_1940_2022_",degrees,".rda")))
  landmask <- get(load(paste0("data/mask_aggr/landmask_",degrees,".rda")))
  if(isTRUE(limit)){landmask <- subsetGrid(landmask,latLim = c(-50,75))}
  } else {
  t2m_ERA5 <- get(load(paste0("data/raw_aggr/t2m_ERA5_monthly_1940_2022_",degrees,".rda")))
  tp_ERA5 <- get(load(paste0("data/raw_aggr/tp_ERA5_monthly_1940_2022_",degrees,".rda")))
  if(isTRUE(limit)){t2m_ERA5 <- subsetGrid(t2m_ERA5,latLim = c(-50,75))}
  
  }

if(!k==0){permutations <-get(load(paste0("data/permutations/permutations_",degrees,lim,superf,".rda")))
backpermutations <- get(load(paste0("data/permutations/backpermutations_",degrees,lim,superf,".rda")))}

data.k <- get(load(file =paste0("data/manipulate_aggr/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".rda")))
 data <- data.k
if(isTRUE(scaleGT)){data$GT <- scale(data$GT)}
if(!k==0){data <-data[,c(permutations[[k]],(length(permutations[[k]])+1):ncol(data))]}
############################################################################
#
############################################################################
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
 
# 
# IT <- paste0(k,"_700_800") #CD
# IT <- paste0(k,"_1000_1100") #HW
# IT <- paste0(k,"_2100_2200") #CW
# IT <- paste0(k,"_1700_1800") #HD 
# IT <- paste0(k,"_4000_4100") #HD 

tabu_bic_0 <- loadIterationsComp(permused = k,algo = algo, score = score,it = IT)
if(isTRUE(global_mean_temp_art)){
  sel.red <- grep("GMT",names(tabu_bic_0))
  tabu_bic_0 <- tabu_bic_0[sel.red]
}
# 
# tabu_bic_0$tabu_disc_3_aic_lim_comp_HDonly_land_0_4000_4100i$nodes$GT
# tabu_bic_0$tabu_disc_3_aic_lim_GMT_art40_comp_HDonly_land_0_4000_4100i$nodes$GT
# tabu_bic_0$tabu_disc_3_aic_lim_comp_HDonly_land_0_4000_4100i$nodes$GT
hist(sapply(tabu_bic_0$tabu_disc_3_aic_lim_comp_HDonly_land_0_4000_4100i$nodes,function(x)length(x$children)+length(x$parents)))
tabu_bic_0_fits<- lapply(tabu_bic_0, function(x) bn.fit(x = x, data = data))

tabu_bic_0_fits$tabu_disc_3_aic_lim_GMT_art80_comp_HDonly_land_0_4000_4100i$GT
tabu_bic_0$tabu_disc_3_aic_lim_GMT_art80_comp_HDonly_land_0_4000_4100i$nodes$GT
tabu_bic_0_fits$tabu_disc_3_aic_lim_GMT_art80_comp_HDonly_land_0_4000_4100i$C914$children
#################################################################################
# Single evidence. V81 postive (+ +)
#################################################################################
fittedbn <- tabu_bic_0_fits[[1]]
if(combi == "only"){ngrids <- length(nodes(fittedbn)) - 1} else if (combi == "TPT2M") {ngrids <- (length(nodes(fittedbn)) - 1)/3}
j<- length(nodes(fittedbn))
lvl.comp <- 2 # C 2 TP 1 T2M3
nE <- c(j)



meth <- "lw"
meth <- ""

col.r <- colorRampPalette(brewer.pal(9, "Reds"))(15)
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(15)
col.rb <- c(rev(col.b),'white','white',col.r)

proplist_0 <- list()
i <- 1

mask_loc <- as.vector(array3Dto2Dmat(redim(landmask,member = FALSE)$Data))



####################################################################
# First to quarter quantile  
####################################################################
i <- 1
meth <- ""
quantiles <- 1:4
str(proplist_0)

for(i in 1:length(quantiles)){
lvl.ev <- quantiles[i]
proplist_0[[i]]<- get(load(paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop",meth,"_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda")))
}
if(!k==0){proplist_1 <- lapply(proplist_0,function(x)x[backpermutations[[k]],])} else {proplist_1 <- proplist_0}


max(proplist_0[[1]]$with)
min(proplist_0[[1]]$with)

if(isTRUE(mask)){
  for(i in 1:length(proplist_1)){
    proplist_1[[i]][nrow(proplist_1[[i]])+1,] <- NA}
  mask_loc_perm <- mask_loc
  mask_loc_perm[!is.na(mask_loc)]<- 1:nrow(proplist_0[[1]])
  mask_loc_perm[is.na(mask_loc_perm)] <- j
  proplist <- proplist_1
  for(i in 1:length(proplist_1)){
    proplist[[i]] <-proplist_1[[i]][mask_loc_perm,]}
} else {proplist <- proplist_1}

if(isTRUE(mask)){
  # climlist1 <- lapply(proplist, function(x) quantity2clim(sign(x$with - x$without)*abs(x$with - x$without), paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = landmask))
  climlist <- lapply(proplist, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = landmask))
  # all.equal(climlist1,climlist)
  climwithlist <- lapply(proplist, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = landmask))
  climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = landmask))
} else {
  climlist <- lapply(proplist, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = t2m_ERA5))
  climwithlist <- lapply(proplist, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = t2m_ERA5))
  climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = t2m_ERA5))
}


climlist <- lapply(climlist,redim,member = TRUE)
enspropdifsclims <- bindGrid(climlist, dimension = c("member"),skip.temporal.check = TRUE)
# enspropdifsclims$Members <- apply(lvls,MARGIN = 2,as.character)
enspropwithsclims <- bindGrid(climwithlist, dimension = c("member"),skip.temporal.check = TRUE)
#enspropwithsclims$Members <- apply(lvls,MARGIN = 2,as.character)
enspropwithoutsclims <- bindGrid(climwithoutlist, dimension = c("member"),skip.temporal.check = TRUE)
#enspropwithoutsclims$Members <- apply(lvls,MARGIN = 2,as.character)

namesmembers <- paste0("q = ",quantiles)

plotname <- paste0("figs/evidence_propagation/GT/propdif_",meth,algo,"_",score,"_disc_",b,lim,comp,combi,superf,n,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_quantiles.pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/evidence_propagation/GT/propdif_",meth,algo,"_",score,"_disc_",b,lim,comp,combi,superf,n,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_quantiles.png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)
spatialPlot(enspropdifsclims,as.table = TRUE,
            names.attr = namesmembers,
            backdrop.theme = "coastline", 
            lonCenter = 0, 
            main = list(paste0(IT," P(",compound," = 2 |",names(fittedbn)[nE]," = q) - P(",compound," = 2) ",cut.art,geo)), 
                  at = seq(-0.2,0.2,0.4/32),
            region = TRUE, 
            col.regions= col.rb,
             set.max = 0.2, 
            set.min = -0.2,
            rev.colors = TRUE,
            colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()

nparams(fittedbn, effective = TRUE)

plotname <- paste0("figs/evidence_propagation/GT/propwith_",meth,"_",algo,"_",score,"_disc_",b,"_",lim,comp,combi,superf,n1,n2,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,"_wlvl_",paste0(gradosvec,collapse = "_"),".pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/evidence_propagation/GT/propwith_",meth,"_",algo,"_",score,"_disc_",b,"_",lim,comp,combi,superf,n1,n2,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,"_wlvl_",paste0(gradosvec,collapse = "_"),".png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)

spatialPlot(enspropwithoutsclims,as.table = TRUE,
            names.attr = namesmembers,
            backdrop.theme = "coastline", lonCenter = 0, main = list(paste0("P(V.t2m = lvl.t2m & V.tp = lvl.tp|",names(fitted)[nE]," = ",lvl.ev,")")), 
            #at = seq(0,max(enspropwithsclims$Data),0.05),
            #region = TRUE, col.regions= col.r,
            at = seq(0,0.4,0.025),
            region = TRUE, 
            col.regions= c('white',col.r),
            set.max = 0.4, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()

###############################################################
# quartile quantile. 
###############################################################
i <- 1
lvl.ev <- c(4)
proplist_0[[i]]<- get(load(paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop",meth,"_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda")))
if(!k==0){proplist_1 <- lapply(proplist_0,function(x)x[backpermutations[[k]],])} else {proplist_1 <- proplist_0}

max(proplist_0[[1]]$with)
proplist_0[[1]]$without

if(isTRUE(mask)){
for(i in 1:length(proplist_1)){
  proplist_1[[i]][nrow(proplist_1[[i]])+1,] <- NA}
mask_loc_perm <- mask_loc
mask_loc_perm[!is.na(mask_loc)]<- 1:nrow(proplist_0[[1]])
mask_loc_perm[is.na(mask_loc_perm)] <- j
proplist <- proplist_1
for(i in 1:length(proplist_1)){
  proplist[[i]] <-proplist_1[[i]][mask_loc_perm,]}
} else {proplist <- proplist_1}

if(isTRUE(mask)){
# climlist1 <- lapply(proplist, function(x) quantity2clim(sign(x$with - x$without)*abs(x$with - x$without), paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = landmask))
climlist <- lapply(proplist, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = landmask))
# all.equal(climlist1,climlist)
climwithlist <- lapply(proplist, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = landmask))
climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = landmask))
} else {
climlist <- lapply(proplist, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = t2m_ERA5))
climwithlist <- lapply(proplist, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = t2m_ERA5))
climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = t2m_ERA5))
}


climlist <- lapply(climlist,redim,member = TRUE)
enspropdifsclims <- bindGrid(climlist, dimension = c("member"),skip.temporal.check = TRUE)
# enspropdifsclims$Members <- apply(lvls,MARGIN = 2,as.character)
enspropwithsclims <- bindGrid(climwithlist, dimension = c("member"),skip.temporal.check = TRUE)
#enspropwithsclims$Members <- apply(lvls,MARGIN = 2,as.character)
enspropwithoutsclims <- bindGrid(climwithoutlist, dimension = c("member"),skip.temporal.check = TRUE)
#enspropwithoutsclims$Members <- apply(lvls,MARGIN = 2,as.character)

namesmembers <- paste0("q = ",lvl.ev)

plotname <- paste0("figs/evidence_propagation/GT/propdif_",meth,algo,"_",score,"_disc_",b,lim,comp,combi,superf,n1,n2,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/evidence_propagation/GT/propdif_",meth,algo,"_",score,"_disc_",b,lim,comp,combi,superf,n1,n2,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)
spatialPlot(enspropdifsclims,as.table = TRUE,
            names.attr = namesmembers,
            backdrop.theme = "coastline", 
            lonCenter = 0, 
            main = list(paste0(IT," P(",compound," = 2 |",names(fittedbn)[nE]," = ",lvl.ev,") - P(",compound," = 2)")), 
             # at = seq(-,0.06,0.01),
            region = TRUE, 
            col.regions= col.rb,
            set.max = 1, 
            rev.colors = TRUE,
            colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()

plotname <- paste0("figs/evidence_propagation/GT/propwith_",meth,"_",algo,"_",score,"_disc_",b,"_",lim,comp,combi,superf,n1,n2,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,"_wlvl_",paste0(gradosvec,collapse = "_"),".pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/evidence_propagation/GT/propwith_",meth,"_",algo,"_",score,"_disc_",b,"_",lim,comp,combi,superf,n1,n2,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,"_wlvl_",paste0(gradosvec,collapse = "_"),".png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)

spatialPlot(enspropwithsclims,as.table = TRUE,
            names.attr = namesmembers,
            backdrop.theme = "coastline", lonCenter = 180, main = list(paste0("P(V.t2m = lvl.t2m & V.tp = lvl.tp|",names(fitted)[nE]," = ",lvl.ev,")")), 
            #at = seq(0,max(enspropwithsclims$Data),0.05),
            #region = TRUE, col.regions= col.r,
            at = seq(0,0.4,0.025),
            region = TRUE, 
            col.regions= c('white',col.r),
            set.max = 0.4, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()

####################################################################
# Warming levels 
####################################################################
gradosvec <- c(2,3,4)
for(i in 1:length(gradosvec)){
  grados <- gradosvec[i]
  if(!isTRUE(scaleGT)){
    lvl.ev <- c(mean(data[["GT"]]) + grados)
  } else {lvl.ev <- c(mean(data[["GT"]]) + grados/sd(data.k[["GT"]]))}
  proplist_0[[i]]<- get(load(paste0("results/propagation/",degrees,"/",algo,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2,"_",IT,"/prop",meth,"_",algo,"_",score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda")))
}

if(!k==0){proplist_1 <- lapply(proplist_0,function(x)x[backpermutations[[k]],])} else {proplist_1 <- proplist_0}
for(i in 1:length(proplist_1)){
  proplist_1[[i]][nrow(proplist_1[[i]])+1,] <- NA}
mask_loc_perm <- mask_loc
mask_loc_perm[!is.na(mask_loc)]<- 1:nrow(proplist_0[[1]])
mask_loc_perm[is.na(mask_loc_perm)] <- j
proplist <- proplist_1
for(i in 1:length(proplist_1)){
  proplist[[i]] <-proplist_1[[i]][mask_loc_perm,]}

# climlist1 <- lapply(proplist, function(x) quantity2clim(sign(x$with - x$without)*abs(x$with - x$without), paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = landmask))
climlist <- lapply(proplist, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = landmask))
# all.equal(climlist1,climlist)
climwithlist <- lapply(proplist, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = landmask))
climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = landmask))

climlist <- lapply(climlist,redim,member = TRUE)
enspropdifsclims <- bindGrid(climlist, dimension = c("member"),skip.temporal.check = TRUE)
# enspropdifsclims$Members <- apply(lvls,MARGIN = 2,as.character)
enspropwithsclims <- bindGrid(climwithlist, dimension = c("member"),skip.temporal.check = TRUE)
#enspropwithsclims$Members <- apply(lvls,MARGIN = 2,as.character)
enspropwithoutsclims <- bindGrid(climwithoutlist, dimension = c("member"),skip.temporal.check = TRUE)
#enspropwithoutsclims$Members <- apply(lvls,MARGIN = 2,as.character)

namesmembers <- paste0("K = ",gradosvec)

plotname <- paste0("figs/evidence_propagation/GT/propdif_",meth,"_",algo,"_",score,"_disc_",b,"_",lim,comp,combi,superf,n1,n2,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,"_wlvl_",paste0(gradosvec,collapse = "_"),".pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/evidence_propagation/GT/propdif_",meth,"_",algo,"_",score,"_disc_",b,"_",lim,comp,combi,superf,n1,n2,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,"_wlvl_",paste0(gradosvec,collapse = "_"),".png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)
spatialPlot(enspropdifsclims,as.table = TRUE,
            names.attr = namesmembers,
            backdrop.theme = "coastline", 
            lonCenter = 0, 
            main = list(paste0(IT," P(",compound," = 2 |",names(fittedbn)[nE]," = ",round(mean(data[["GT"]]))," + K) - P(",compound," = 2)")), 
            at = seq(-0.25,0.25,0.025),
            region = TRUE, 
            col.regions= col.rb,
            set.max = 0.25,
            set.min = -0.25,
            rev.colors = TRUE,
            colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()

plotname <- paste0("figs/evidence_propagation/GT/propwith_",meth,"_",algo,"_",score,"_disc_",b,"_",lim,comp,combi,superf,n1,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,"_wlvl_",paste0(gradosvec,collapse = "_"),".pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/evidence_propagation/GT/propwith_",meth,"_",algo,"_",score,"_disc_",b,"_",lim,comp,combi,superf,n1,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,"_wlvl_",paste0(gradosvec,collapse = "_"),".png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)

spatialPlot(enspropwithsclims,as.table = TRUE,
            names.attr = namesmembers,
            backdrop.theme = "coastline", lonCenter = 180, main = list(paste0("P(V.t2m = lvl.t2m & V.tp = lvl.tp|",names(fitted)[nE]," = ",lvl.ev,")")), 
            #at = seq(0,max(enspropwithsclims$Data),0.05),
            #region = TRUE, col.regions= col.r,
            at = seq(0,0.4,0.025),
            region = TRUE, 
            col.regions= c('white',col.r),
            set.max = 0.4, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()

# plotname <- paste0("figs/evidence_propagation/propwithout_",Algo,"_",Score,"_",IT,"_t2m_tp_4lvls_t2m_",nE,"_",lvl.ev,".pdf")
# pdf(plotname, width = 10, height = 7)
# plotname <- paste0("figs/evidence_propagation/propwithout_",Algo,"_",Score,"_",IT,"_t2m_tp_4lvls_t2m_",nE,"_",lvl.ev,".png")
# png(plotname, width = 20, height = 15,units = "cm", res = 150)
# 
# 
# spatialPlot(enspropwithoutsclims,as.table = TRUE,names.attr = namesmembers,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0("P(V.t2m = lvl.t2m & V.tp = lvl.tp)")), 
#             #at = seq(0,max(enspropwithsclims$Data),0.05),
#             #region = TRUE, col.regions= col.r,
#             at = seq(0,0.4,0.025),
#             region = TRUE, 
#             col.regions= c('white',col.r),
#             set.max = 0.4, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
# 
# dev.off()
# 

