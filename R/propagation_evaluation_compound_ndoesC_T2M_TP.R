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


# sel.tp<- NULL
# sel.t2m <- NULL
############################################
# discretize in b
############################################
b <- 3
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


loadIterationsComp <- function(permused = 0, algo = c("hc","tabu"),score = c("aic","bic"), it = NULL, limit = FALSE,compound = NULL,artificial = FALSE) {
  if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
  if (is.null(compound)){comp <-""} else if(compound =="comp"){comp<- "_comp"} else if(!is.null(compound)){comp <- paste0("_comp_",compound)} 
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

Permused <- 0
Algo <- "tabu"
Score <- "aic"
IT <- "0_1600_1700"
IT <- "0_2600_2700"
IT <- "0_2700_2800"
tabu_aic_0 <- loadIterationsComp(permused = Permused,algo = Algo , score = Score,it = IT, limit = TRUE)
# tabu_bic_0 <- loadIterationsComp(permused = 0,algo = "tabu", score = "bic")
# hc_aic_0 <- loadIterationsComp(permused = 0,algo = "hc", score = "aic")
# hc_bic_0 <- loadIterationsComp(permused = 0,algo = "hc", score = "bic")


tabu_aic_0_fits<- lapply(tabu_aic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))
# tabu_bic_0_fits<- lapply(tabu_bic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))
# hc_aic_0_fits<- lapply(hc_aic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))
# hc_bic_0_fits<- lapply(hc_bic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))
nparams(tabu_aic_0_fits$tabu_disc_3_aic_t2m_tp_10d_lim_0_2700_2800i,effective = TRUE)
#################################################################################
# Single evidence. V81 postive (+ +)
#################################################################################
# load("/data/Untitled/Trabajo/R_practice/exp_compound/results/propagation/prop_t2m4_tp4_t2m81_4.rda")
# load("/data/Untitled/Trabajo/R_practice/exp_compound/results/propagation/prop_t2m1_tp1_t2m81_4.rda")
# load("/data/Untitled/Trabajo/R_practice/exp_compound/results/propagation/prop_t2m4_tp1_t2m81_4.rda")
# load("/data/Untitled/Trabajo/R_practice/exp_compound/results/propagation/prop_t2m1_tp4_t2m81_4.rda")

# lvl.t2m <- 1
# lvl.tp <- 1
nE <- 57
lvl.ev <- b
lvls <- rbind(c(b,b),c(1,1),c(b,1),c(1,b))

proplist <- list()
for(i in (1:nrow(lvls)) ){
  if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
lvl.t2m <- lvls[i,1]
lvl.tp <- lvls[i,2]
proplist[[i]]<- get(load(paste0("results/propagation/",Algo,"_",Score,"_",IT,lim,"/prop_",Algo,"_",Score,"_",IT,lim,"_t2m",lvl.t2m,"_tp",lvl.tp,"_t2m",nE,"_",lvl.ev,".rda")))
}


fitted <- tabu_aic_0_fits[[1]]
col.r <- colorRampPalette(brewer.pal(9, "Reds"))(15)
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(15)
col.rb <- c(rev(col.b),'white','white',col.r)

proplist[[1]]$with
climlist <- lapply(proplist, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithindlist <- lapply(proplist, function(x) quantity2clim(x$with-(1/b*1/b), paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
climwithlist <- lapply(proplist, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = sel.t2m))
climwithoutindlist <- lapply(proplist, function(x) quantity2clim(x$without-(1/b*1/b), paste0(attr(x$without, "probability")),ref.grid = sel.t2m))


enspropdifsclims <- bindGrid(climlist, dimension = c("member"))
enspropdifsclims$Members <- apply(lvls,MARGIN = 2,as.character)
enspropwithsclims <- bindGrid(climwithlist, dimension = c("member"))
enspropwithsclims$Members <- apply(lvls,MARGIN = 2,as.character)
enspropwithindsclims <- bindGrid(climwithindlist, dimension = c("member"))
enspropwithindsclims$Members <- apply(lvls,MARGIN = 2,as.character)
enspropwithoutsclims <- bindGrid(climwithoutlist, dimension = c("member"))
enspropwithoutsclims$Members <- apply(lvls,MARGIN = 2,as.character)
enspropwithoutindsclims <- bindGrid(climwithoutindlist, dimension = c("member"))
enspropwithoutindsclims$Members <- apply(lvls,MARGIN = 2,as.character)


namesmembers <- apply(lvls,MARGIN = 1, function(x)paste0("lvl.t2m = ",x[1]," & lvl.tp = ",x[2]))

plotname <- paste0("figs/evidence_propagation/propdif_",Algo,"_",Score,"_",IT,lim,"_t2m_tp_",b,"lvls_t2m_",nE,"_",lvl.ev,".pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/evidence_propagation/propdif_",Algo,"_",Score,"_",IT,lim,"_t2m_tp_",b,"lvls_t2m_",nE,"_",lvl.ev,".png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)

spatialPlot(enspropdifsclims,as.table = TRUE,names.attr = namesmembers,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0("P(V.t2m = lvl.t2m & V.tp = lvl.tp|",names(fitted)[nE]," = ",lvl.ev,") - P(V.t2m = lvl.t2m & V.tp = lvl.tp)")), 
            at = seq(-0.4,0.4,0.025),
            region = TRUE, 
            col.regions= col.rb,
            set.max = 0.4, 
            rev.colors = TRUE,
            colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()

plotname <- paste0("figs/evidence_propagation/propwith_",Algo,"_",Score,"_",IT,lim,"_t2m_tp_",b,"lvls_t2m_",nE,"_",lvl.ev,"_respind.pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/evidence_propagation/propwith_",Algo,"_",Score,"_",IT,lim,"_t2m_tp_",b,"lvls_t2m_",nE,"_",lvl.ev,"_respind.png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)

# this one for with - (1/b) * (1/b)
spatialPlot(enspropwithindsclims,as.table = TRUE,names.attr = namesmembers,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0("P(V.t2m = lvl.t2m & V.tp = lvl.tp|",names(fitted)[nE]," = ",lvl.ev,")")), 
            #at = seq(0,max(enspropwithsclims$Data),0.05),
            #region = TRUE, col.regions= col.r,
            at = seq(-0.2,0.2,0.025),
            region = TRUE, 
            col.regions= c(rev(col.b),"white",'white',col.r),
            set.max = 0.2, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()

plotname <- paste0("figs/evidence_propagation/propwith_",Algo,"_",Score,"_",IT,lim,"_t2m_tp_",b,"lvls_t2m_",nE,"_",lvl.ev,".pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/evidence_propagation/propwith_",Algo,"_",Score,"_",IT,lim,"_t2m_tp_",b,"lvls_t2m_",nE,"_",lvl.ev,".png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)

# this one for without -0.0625
spatialPlot(enspropwithsclims,as.table = TRUE,names.attr = namesmembers,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0("P(V.t2m = lvl.t2m & V.tp = lvl.tp|",names(fitted)[nE]," = ",lvl.ev,")")), 
            #at = seq(0,max(enspropwithsclims$Data),0.05),
            #region = TRUE, col.regions= col.r,
            at = seq(0,0.4,0.025),
            region = TRUE, 
            col.regions= c('white',col.r),
            set.max = 0.4, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()

plotname <- paste0("figs/evidence_propagation/propwithout_",Algo,"_",Score,"_",IT,lim,"_t2m_tp_",b,"lvls_t2m_",nE,"_",lvl.ev,".pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/evidence_propagation/propwithout_",Algo,"_",Score,"_",IT,lim,"_t2m_tp_",b,"lvls_t2m_",nE,"_",lvl.ev,".png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)


spatialPlot(enspropwithoutsclims,as.table = TRUE,names.attr = namesmembers,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0("P(V.t2m = lvl.t2m & V.tp = lvl.tp)")), 
            #at = seq(0,max(enspropwithsclims$Data),0.05),
            #region = TRUE, col.regions= col.r,
            at = seq(0,0.4,0.025),
            region = TRUE, 
            col.regions= c('white',col.r),
            set.max = 0.4, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()

plotname <- paste0("figs/evidence_propagation/propwithout_",Algo,"_",Score,"_",IT,lim,"_t2m_tp_",b,"lvls_t2m_",nE,"_",lvl.ev,"_respind.pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/evidence_propagation/propwithout_",Algo,"_",Score,"_",IT,lim,"_t2m_tp_",b,"lvls_t2m_",nE,"_",lvl.ev,"_respind.png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)


spatialPlot(enspropwithoutindsclims,as.table = TRUE,names.attr = namesmembers,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0("P(V.t2m = lvl.t2m & V.tp = lvl.tp)")), 
            #at = seq(0,max(enspropwithsclims$Data),0.05),
            #region = TRUE, col.regions= col.r,
            at = seq(-0.2,0.2,0.025),
            region = TRUE, 
            col.regions= c(rev(col.b),"white",'white',col.r), set.min = -0.2,
            set.max = 0.2, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()



prop_dif <- mapply(function(prop,griddata) quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), griddata),prop = propfiles, griddata = list(cubicsets[[1]]), SIMPLIFY = FALSE)
prop444 <- prop_t2m4_tp4_t2m81_4
propclim444<- quantity2clim(prop444$with - prop444$without, paste0("P(V.t2m = lvl.t2m & V.tp = lvl.tp|",names(fitted)[nE]," = ",lvl.ev,") - P(V.t2m = lvl.t2m & V.tp = lvl.tp)"),ref.grid = sel.t2m)
spatialPlot(propclim444, as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, at = seq(0.00,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

prop_t2m1_tp1_t2m81_4
prop114 <- prop_t2m1_tp1_t2m81_4
propclim114<- quantity2clim(prop114$with - prop114$without, paste0(attr(prop114$with, "probability"),"-", attr(prop114$without, "probability")),ref.grid = sel.t2m)
spatialPlot(propclim114, as.table = TRUE,backdrop.theme = "coastline", lonCenter = 180, at = seq(0.00,0.85,0.01),region = TRUE, col.regions= col.r,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))


#################################################################################
# Double evidence. V81 T2m and TP postive (+ +)
#################################################################################
fitted <- tabu_aic_0_fits[[1]]
lvl.t2m <- 1
lvl.tp <- 1
nE <- c(81,81+length(nodes(fitted))/2)
lvl.ev <- c(4,4)

lvls <- rbind(c(4,4),c(1,1),c(4,1),c(1,4))

proplist <- list()
for(i in (1:nrow(lvls)) ){
  lvl.t2m <- lvls[i,1]
  lvl.tp <- lvls[i,2]
  proplist[[i]]<- get(load(paste0("results/propagation/",Algo,"_",Score,"_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m",lvl.t2m,"_tp",lvl.tp,"_t2m",nE[1],"_",lvl.ev[1],"_tp",nE[1],"_",lvl.ev[2],".rda")))
}



col.r <- colorRampPalette(brewer.pal(9, "Reds"))(15)
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(15)
col.rb <- c(rev(col.b),'white','white',col.r)

proplist[[1]]$with
climlist <- lapply(proplist, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlist <- lapply(proplist, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = sel.t2m))



enspropdifsclims <- bindGrid(climlist, dimension = c("member"))
enspropdifsclims$Members <- apply(lvls,MARGIN = 2,as.character)
enspropwithsclims <- bindGrid(climwithlist, dimension = c("member"))
enspropwithsclims$Members <- apply(lvls,MARGIN = 2,as.character)
enspropwithoutsclims <- bindGrid(climwithoutlist, dimension = c("member"))
enspropwithoutsclims$Members <- apply(lvls,MARGIN = 2,as.character)



namesmembers <- apply(lvls,MARGIN = 1, function(x)paste0("lvl.t2m = ",x[1]," & lvl.tp = ",x[2]))

plotname <- paste0("figs/evidence_propagation/propdif_",Algo,"_",Score,"_",IT,"_t2m_tp_4lvls_t2m_",nE[1],"_",lvl.ev[1],"_tp_",nE[1],"_",lvl.ev[2],".pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/evidence_propagation/propdif_",Algo,"_",Score,"_",IT,"_t2m_tp_4lvls_t2m_",nE[1],"_",lvl.ev[1],"_tp_",nE[1],"_",lvl.ev[2],".png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)

spatialPlot(enspropdifsclims,as.table = TRUE,names.attr = namesmembers,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0("P(V.t2m = lvl.t2m & V.tp = lvl.tp|",names(fitted)[nE[1]]," = ",lvl.ev[1]," & ",names(fitted)[nE[2]]," = ",lvl.ev[2],") - P(V.t2m = lvl.t2m & V.tp = lvl.tp)")), 
            at = seq(-0.4,0.4,0.025),
            region = TRUE, 
            col.regions= col.rb,
            set.max = 0.4, 
            rev.colors = TRUE,
            colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()

plotname <- paste0("figs/evidence_propagation/propwith_",Algo,"_",Score,"_",IT,"_t2m_tp_4lvls_t2m_",nE[1],"_",lvl.ev[1],"_tp_",nE[1],"_",lvl.ev[2],".pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/evidence_propagation/propwith_",Algo,"_",Score,"_",IT,"_t2m_tp_4lvls_t2m_",nE[1],"_",lvl.ev[1],"_tp_",nE[1],"_",lvl.ev[2],".png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)

spatialPlot(enspropwithsclims,as.table = TRUE,names.attr = namesmembers,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0("P(V.t2m = lvl.t2m & V.tp = lvl.tp|",names(fitted)[nE[1]]," = ",lvl.ev[1]," & ",names(fitted)[nE[2]]," = ",lvl.ev[2],")")), 
            #at = seq(0,max(enspropwithsclims$Data),0.05),
            #region = TRUE, col.regions= col.r,
            at = seq(0,0.4,0.025),
            region = TRUE, 
            col.regions= c('white',col.r),
            set.max = 0.4, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()

plotname <- paste0("figs/evidence_propagation/propwithout_",Algo,"_",Score,"_",IT,"_t2m_tp_4lvls_t2m_",nE[1],"_",lvl.ev[1],"_tp_",nE[1],"_",lvl.ev[2],".pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/evidence_propagation/propwithout_",Algo,"_",Score,"_",IT,"_t2m_tp_4lvls_t2m_",nE[1],"_",lvl.ev[1],"_tp_",nE[1],"_",lvl.ev[2],".png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)


spatialPlot(enspropwithoutsclims,as.table = TRUE,names.attr = namesmembers,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0("P(V.t2m = lvl.t2m & V.tp = lvl.tp)")), 
            #at = seq(0,max(enspropwithsclims$Data),0.05),
            #region = TRUE, col.regions= col.r,
            at = seq(0,0.4,0.025),
            region = TRUE, 
            col.regions= c('white',col.r),
            set.max = 0.4, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()

#############################################################################
# Propgation of evidence compound nodes
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
# discretize in 4
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

tabu_aic_0_comp_art_fit<- lapply(tabu_bic_0_comp_art, function(x) bn.fit(x = x, data = disc.df.t2m.tp))
# nparams(tabu_aic_0_comp_art_fit[[1]])
# tabu_aic_0_comp_art_fit$tabu_disc_3_bic_t2m_tp_10d_lim_comp_art_0_1700_1800i$C468
# 1404*1404
#################################################################################
# Single evidence. Ci =  2
#################################################################################
fittedbn <- tabu_aic_0_comp_art_fit$tabu_disc_3_aic_t2m_tp_10d_lim_comp_art_0_2700_2800i
proplist <- proplistt2m3 <- proplistt2m1 <- proplisttp1 <- proplisttp3 <- list()
lvl.comp <- 2
lvl.t2m3 <- lvl.tp3 <- 3
lvl.t2m1 <- lvl.tp1 <- 1
lvl.ev <- c(2)
cnodes <- paste0("C",1:ngrids)
for(j in (1:ngrids) ){
  nE <- c(j+2*ngrids)
  proplist[[j]]<- get(load(
    file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.ev,".rda")))
  proplistt2m3[[j]]<- get(load(
    file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.t2m3,"_",names(fittedbn)[nE],"_",lvl.ev,".rda")))
  # proplistt2m1[[j]]<- get(load(
    # file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.t2m1,"_",names(fittedbn)[nE],"_",lvl.ev,".rda")))
  proplisttp1[[j]]<- get(load(
    file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.tp1,"_",names(fittedbn)[nE],"_",lvl.ev,".rda")))
  # proplisttp3[[j]]<- get(load(
    # file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.tp3,"_",names(fittedbn)[nE],"_",lvl.ev,".rda")))
  
  }


col.r <- colorRampPalette(brewer.pal(9, "Reds"))(15)
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(15)
col.rb <- c(rev(col.b),'white','white',col.r)


#proplist2 <- proplist
proplist <- lapply(1:length(proplist), function(x) {proplist[[x]]$without <- proplist[[1]]$without ; return(proplist[[x]])})
proplistt2m3 <- lapply(1:length(proplistt2m3), function(x) {proplistt2m3[[x]]$without <- proplistt2m3[[1]]$without ; return(proplistt2m3[[x]])})
proplisttp1 <- lapply(1:length(proplisttp1), function(x) {proplisttp1[[x]]$without <- proplisttp1[[1]]$without ; return(proplisttp1[[x]])})
proplistt2m1 <- lapply(1:length(proplistt2m1), function(x) {proplistt2m1[[x]]$without <- proplistt2m1[[1]]$without ; return(proplistt2m1[[x]])})
proplisttp3 <- lapply(1:length(proplisttp3), function(x) {proplisttp3[[x]]$without <- proplisttp3[[1]]$without ; return(proplisttp3[[x]])})


#propdiflist <- lapply(1:length(proplist), function(x) {proplist[[x]]$without <- proplist[[1]]$without ; return(proplist[[x]])})

climlist <- lapply(proplist, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlist <- lapply(proplist, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
#climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = sel.t2m))
withoutvector <- proplist[[1]]$without 
withoutclim <- quantity2clim(withoutvector, paste0(attr(withoutvector, "probability")),ref.grid = sel.t2m)

climlistt2m3 <- lapply(proplistt2m3, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlistt2m3 <- lapply(proplistt2m3, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
#climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = sel.t2m))
withoutvectort2m3 <- proplistt2m3[[1]]$without 
withoutclimt2m3 <- quantity2clim(withoutvectort2m3, paste0(attr(withoutvectort2m3, "probability")),ref.grid = sel.t2m)

climlisttp1 <- lapply(proplisttp1, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlisttp1 <- lapply(proplisttp1, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
#climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = sel.t2m))
withoutvectortp1 <- proplisttp1[[1]]$without 
withoutclimtp1 <- quantity2clim(withoutvectortp1, paste0(attr(withoutvectortp1, "probability")),ref.grid = sel.t2m)

climlistt2m1 <- lapply(proplistt2m1, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlistt2m1 <- lapply(proplistt2m1, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
#climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = sel.t2m))
withoutvectort2m1 <- proplistt2m1[[1]]$without 
withoutclimt2m1 <- quantity2clim(withoutvectort2m1, paste0(attr(withoutvectort2m1, "probability")),ref.grid = sel.t2m)

climlisttp3 <- lapply(proplisttp3, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlisttp3 <- lapply(proplisttp3, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
#climwithoutlist <- lapply(proplist, function(x) quantity2clim(x$without, paste0(attr(x$without, "probability")),ref.grid = sel.t2m))
withoutvectortp3 <- proplisttp3[[1]]$without 
withoutclimtp3 <- quantity2clim(withoutvectortp3, paste0(attr(withoutvectortp3, "probability")),ref.grid = sel.t2m)




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


textdf <- cbind(attr(sel.t2m.1,"VertexCoords"),1:ngrids)
#textdf$x[textdf$x<0] <-textdf$x[textdf$x<0] + 360
textdf <- as.matrix(textdf)
textlay <- apply(textdf,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

clim.0vector <- quantity2clim(numeric(length=ngrids),"node vector", ref.grid = sel.t2m)
Xvectorname <-attr(clim.0vector$Data,"climatology:fun")
Xvector <-spatialPlot(clim.0vector, main = Xvectorname,  backdrop.theme = "coastline",sp.layout = textlay, lonCenter =0, colorkey = FALSE)

# KIJK IUT
 # spatialPlot(enspropdifsclims,as.table = TRUE,names.attr = namesmembers,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0("P(C = HD & V.tp = lvl.tp|",names(fitted)[nE]," = ",lvl.ev,") - P(V.t2m = lvl.t2m & V.tp = lvl.tp)")), 
 #            at = seq(-0.4,0.4,0.025),
 #            region = TRUE, 
 #            col.regions= col.rb,
 #            set.max = 0.4, 
 #            rev.colors = TRUE,
 #            colorkey = list(width = 0.6, lables = list(cex = 0.5)))



#dev.off()
# spatialPlot(withoutclim,as.table = TRUE,names.attr = namesmembers,backdrop.theme = "coastline", lonCenter = 0, 
#             main = list(paste0("P(C = HD)")), 
#             #at = seq(0,0.4,0.025),
#             region = TRUE, 
#             col.regions= col.r,
#             set.max = 0.4, 
#             rev.colors = TRUE,
#             colorkey = list(width = 0.6, lables = list(cex = 0.5)))

gridx <- 254
pC2 <- spatialPlot(climlist[[gridx]],
                   #as.table = TRUE,names.attr = namesmembers,
                   backdrop.theme = "coastline", lonCenter = 0, 
                   main = list(paste0("P(V.Comp.HD = 2 (i.e. HotDry = TRUE) |",names(fittedbn)[2*ngrids+gridx]," = ",lvl.ev," (HotDry)) - P(V.Comp.HD = 2 (TRUE)")), 
                   at = seq(-0.4,0.4,0.025),
            region = TRUE, 
            col.regions= col.rb,
            set.max = 0.4, 
            rev.colors = TRUE,
            colorkey = list(width = 0.6, lables = list(cex = 0.5)))

pTP1 <-spatialPlot(climlisttp1[[gridx]],
                   #as.table = TRUE,names.attr = namesmembers,
                   backdrop.theme = "coastline", lonCenter = 0, 
                   main = list(paste0("P(V.tp = 1 (Dry)|",names(fittedbn)[2*ngrids+gridx]," = ",lvl.ev," (HotDry)) - P(V.tp = 1 (Dry)")), 
                   at = seq(-0.4,0.4,0.025),
            region = TRUE, 
            col.regions= col.rb,
            set.max = 0.4, 
            rev.colors = TRUE,
            colorkey = list(width = 0.6, lables = list(cex = 0.5)))
pT2M3 <- spatialPlot(climlistt2m3[[gridx]],
                     as.table = TRUE,names.attr = namesmembers,
                     backdrop.theme = "coastline", 
                     lonCenter = 0, 
                     main = list(paste0("P(V.t2m = 3 (Hot)|",names(fittedbn)[2*ngrids+gridx]," = ",lvl.ev," (HotDry)) - P(V.t2m = 3 (Hot)")), 
                     at = seq(-0.4,0.4,0.025),
            region = TRUE, 
            col.regions= col.rb,
            set.max = 0.4, 
            rev.colors = TRUE,
            colorkey = list(width = 0.6, lables = list(cex = 0.5)))
pT2M1 <- spatialPlot(climlistt2m1[[gridx]],
                     #as.table = TRUE,names.attr = namesmembers,
                     backdrop.theme = "coastline", lonCenter = 0,
                     main = list(paste0("P(V.t2m = 1 (Cold)|",names(fittedbn)[2*ngrids+gridx]," = ",lvl.ev," (HotDry)) - P(V.t2m = 1 (Cold)")), 
                     at = seq(-0.4,0.4,0.025),
            region = TRUE, 
            col.regions= col.rb,
            set.max = 0.4, 
            rev.colors = TRUE,
            colorkey = list(width = 0.6, lables = list(cex = 0.5)))
pTP3 <- spatialPlot(climlisttp3[[gridx]],
                    #as.table = TRUE,names.attr = namesmembers,
                    backdrop.theme = "coastline", lonCenter = 0, 
                    main = list(paste0("P(V.tp = 3 (Wet)|",names(fittedbn)[2*ngrids+gridx]," = ",lvl.ev," (HotDry)) - P(V.tp3 = 3 (Wet)")), 
            at = seq(-0.40,0.4,0.025),
            region = TRUE, 
            col.regions= col.rb,
            set.max = 0.4, 
            rev.colors = TRUE,
            colorkey = list(width = 0.6, lables = list(cex = 0.5)))

plotname <- paste0("figs/HD_impact/HD_impact_simple_",IT,"_C",gridx)
pdf(paste0(plotname,".pdf"),width = 20, height = 10)
plotlist<- list(pC2,Xvector,pT2M3,pTP1)
do.call(grid.arrange,c(plotlist,nrow = 2))
dev.off()

plotname <- paste0("figs/HD_impact/HD_impact_simple_",IT,"_C",gridx,".png")
png(plotname, width = 50, height = 20,units = "cm", res = 300)
plotlist<- list(pC2,Xvector,pT2M3,pTP1)
do.call(grid.arrange,c(plotlist,nrow = 2))
dev.off()

plotname <- paste0("figs/HD_impact/HD_impact_all_C",gridx)
pdf(paste0(plotname,".pdf"),width = 20, height = 15)
plotlist<- list(pC2,Xvector,pT2M3,pTP1,pT2M1,pTP3)
do.call(grid.arrange,c(plotlist,nrow = 3))
dev.off()


# Hier bezig

plotname <- paste0("figs/HD_impact/HD_impact_to_HD_C",gridx)
png(paste0(plotname,".png"),width = 20, height = 5, units = "cm", res = 300)

spatialPlot(climwithlist[[gridx]],
            #as.table = TRUE,names.attr = namesmembers,
            backdrop.theme = "coastline", lonCenter = 0, 
            main = list(paste0("P(V.Comp.HD = 2 (i.e. HotDry = TRUE) |",names(fittedbn)[2*ngrids+gridx]," = ",lvl.ev," (HotDry))")), 
            at = seq(0,0.4,0.025),
            region = TRUE, 
            col.regions= c("white",col.r),
            set.max = 0.4, 
            rev.colors = TRUE,
            colorkey = list(width = 0.6, lables = list(cex = 0.5)))

dev.off()


#######################################################################
# Sensitive nodes analysis
#######################################################################
posdiffun <- function(x) {if (x >= 0) {x} else {0}}

sum(abs(proplist[[3]]$with-proplist[[3]]$without) >= 0.1)


proplist

HDtoHDextension <- sapply(proplist, function(x){sum(abs(x$with-x$without)>= 0.05)})
HDtoHDimpact <- sapply(proplist, function(x){mean(abs(x$with-x$without))})
HDtoHDimpactclim <- quantity2clim(HDtoHDimpact, "HD-2-HD impact",ref.grid = sel.t2m)
HDtoHDname <-attr(HDtoHDimpactclim$Data,"climatology:fun")

assign("spatextHDtoHD",HDtoHDextension)
save(spatextHDtoHD, file = paste0("results/propagation/spatext_",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/spatext_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_C_",lvl.ev,".rda")
)


textHD <- HDtoHDextension
textHD <- round(HDtoHDimpact*1000,0)
textdfHD <- cbind(attr(sel.t2m.1,"VertexCoords"),textHD)
#textdf$x[textdf$x<0] <- textdf$x[textdf$x<0] + 360
textdfHD <- as.matrix(textdfHD)
textlayHD <- apply(textdfHD,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

HDtoHDimpactclim$Data
HD <- spatialPlot(HDtoHDimpactclim,backdrop.theme = "coastline", lonCenter = 0, 
                  sp.layout = textlayHD,
            rev.colors = FALSE,
            main = HDtoHDname,
            color.theme = "Purples",
            #at = c(0.15,0.25,0.0125),
            #set.min = 0.15,
            colorkey = list(width = 0.6, lables = list(cex = 1)))
HD
dev.off()
HDtoHextension <- sapply(proplistt2m3, function(x){sum(abs(x$with-x$without)>= 0.05)})
HDtoHimpact <- sapply(proplistt2m3, function(x){mean(abs(x$with-x$without))})
HDtoHimpactclim <- quantity2clim(HDtoHimpact, "HD-2-H impact",ref.grid = sel.t2m)
HDtoHname <-attr(HDtoHimpactclim$Data,"climatology:fun")



textH <- HDtoHextension
textH <- round(HDtoHimpact*1000,0)
textdfH <- cbind(attr(sel.t2m.1,"VertexCoords"),textH)
#textdf$x[textdf$x<0] <- textdf$x[textdf$x<0] + 360
textdfH <- as.matrix(textdfH)
textlayH <- apply(textdfH,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

H <- spatialPlot(HDtoHimpactclim,backdrop.theme = "coastline", lonCenter = 0, sp.layout = textlayH,
            rev.colors = FALSE, main = HDtoHname,
            color.theme = "Reds",
            colorkey = list(width = 0.6, lables = list(cex = 1)))


HDtoDextension <- sapply(proplisttp1, function(x){sum(abs(x$with-x$without)>= 0.05)})

HDtoDimpact <- sapply(proplisttp1, function(x){mean(abs(x$with-x$without))})
HDtoDimpactclim <- quantity2clim(HDtoDimpact, "HD-2-D impact",ref.grid = sel.t2m)
HDtoDname <-attr(HDtoDimpactclim$Data,"climatology:fun")

textD <- HDtoDextension
textD <- round(HDtoDimpact*1000,0)
textdfD <- cbind(attr(sel.t2m.1,"VertexCoords"),textD)
#textdf$x[textdf$x<0] <- textdf$x[textdf$x<0] + 360
textdfD <- as.matrix(textdfD)
textlayD <- apply(textdfD,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

D <-spatialPlot(HDtoDimpactclim,backdrop.theme = "coastline", lonCenter = 0, sp.layout = textlayD,
            rev.colors = FALSE, main = HDtoDname,
            color.theme = "Greys",
            colorkey = list(width = 0.6, lables = list(cex = 1)))


HDtoCimpact <- sapply(proplistt2m1, function(x){mean(abs(x$with-x$without))})
HDtoCimpactclim <- quantity2clim(HDtoCimpact, "HD-2-C impact",ref.grid = sel.t2m)
HDtoCname <-attr(HDtoCimpactclim$Data,"climatology:fun")

textC <- round(HDtoCimpact*1000,0)
textdfC <- cbind(attr(sel.t2m.1,"VertexCoords"),textC)
#textdf$x[textdf$x<0] <- textdf$x[textdf$x<0] + 360
textdfC <- as.matrix(textdfC)
textlayC <- apply(textdfC,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

C <- spatialPlot(HDtoCimpactclim,backdrop.theme = "coastline", lonCenter = 0,sp.layout = textlayC,
                 rev.colors = FALSE, main = HDtoCname,
                 color.theme = "Blues",
                 colorkey = list(width = 0.6, lables = list(cex = 0.5)))


HDtoWimpact <- sapply(proplisttp3, function(x){mean(abs(x$with-x$without))})
HDtoWimpactclim <- quantity2clim(HDtoWimpact, "HD-2-W impact",ref.grid = sel.t2m)
HDtoWname <-attr(HDtoWimpactclim$Data,"climatology:fun")

textW <- round(HDtoWimpact*1000,0)
textdfW <- cbind(attr(sel.t2m.1,"VertexCoords"),textW)
#textdf$x[textdf$x<0] <- textdf$x[textdf$x<0] + 360
textdfW <- as.matrix(textdfW)
textlayW <- apply(textdfW,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

W<-spatialPlot(HDtoWimpactclim,backdrop.theme = "coastline", lonCenter = 0,sp.layout = textlayW,
                rev.colors = FALSE, main = HDtoWname,
                color.theme = "Greens",
                colorkey = list(width = 0.6, lables = list(cex = 0.5)))



plotname <- paste0("figs/HD_impact/HD_impact_simple_",IT)
pdf(paste0(plotname,".pdf"),width = 20, height = 10)
png(paste0(plotname,".png"),units = "cm",width = 40, height = 40,res = 600)
grid.arrange(HD,Xvector,H,D)
dev.off()

plotname <- "figs/HD_impact/HD_impact_all"
pdf(paste0(plotname,".pdf"),width = 20, height = 10)
png(paste0(plotname,".png"))
grid.arrange(HD,Xvector,H,D,C,W)
dev.off()


#############################################################################
# Propgation of evidence compound, t2m, and tp nodes
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
# discretize in 4
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

tabu_aic_0_comp_art_fit<- lapply(tabu_bic_0_comp_art, function(x) bn.fit(x = x, data = disc.df.t2m.tp))
# nparams(tabu_aic_0_comp_art_fit[[1]])
# tabu_aic_0_comp_art_fit$tabu_disc_3_bic_t2m_tp_10d_lim_comp_art_0_1700_1800i$C468
# 1404*1404
#################################################################################
# Single evidence. Ci =  2
#################################################################################
fittedbn <- tabu_aic_0_comp_art_fit$tabu_disc_3_aic_t2m_tp_10d_lim_comp_art_0_2700_2800i
proplistC2C2 <- proplistC2t2m3 <- proplistC2tp1 <- proplisttp1C2 <- proplisttp1t2m3 <- proplisttp1tp1<- proplistt2m3C2 <- proplistt2m3t2m3 <- proplistt2m3tp1<- list()
lvl.comp <- 2
lvl.t2m3 <-  3
lvl.t2m1 <-  1
lvl.tp1 <- 1
lv.tp3 <- 3
# load for evidence variable C
cnodes <- paste0("C",1:ngrids)
for(j in (1:ngrids) ){
  nE <- c(j+2*ngrids)
  proplistC2C2[[j]]<- get(load(
    file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.comp,".rda")))
  proplistC2t2m3[[j]]<- get(load(
    file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.t2m3,"_",names(fittedbn)[nE],"_",lvl.comp,".rda")))
  proplistC2tp1[[j]]<- get(load(
    file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.tp1,"_",names(fittedbn)[nE],"_",lvl.comp,".rda")))
}
# load for evidence variable tp1
for(j in (1:ngrids) ){
  nE <- c(j+ngrids)
  proplisttp1C2[[j]]<- get(load(
    file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.tp1,".rda")))
  proplisttp1t2m3[[j]]<- get(load(
    file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.t2m3,"_",names(fittedbn)[nE],"_",lvl.tp1,".rda")))
  proplisttp1tp1[[j]]<- get(load(
    file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.tp1,"_",names(fittedbn)[nE],"_",lvl.tp1,".rda")))
}
# load for evidence variable t2m3
for(j in (1:ngrids) ){
  nE <- c(j)
  proplistt2m3C2[[j]]<- get(load(
    file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_",names(fittedbn)[nE],"_",lvl.t2m3,".rda")))
  proplistt2m3t2m3[[j]]<- get(load(
    file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_t2m_",lvl.t2m3,"_",names(fittedbn)[nE],"_",lvl.t2m3,".rda")))
  proplistt2m3tp1[[j]]<- get(load(
    file = paste0("results/propagation/",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/prop_",Algo,"_",Score,"_",IT,"_tp_",lvl.tp1,"_",names(fittedbn)[nE],"_",lvl.t2m3,".rda")))
}

#Add without vector to all proplists
# proplistC2C2 <- proplistC2t2m3 <- proplistC2tp1 <- proplisttp1C2 <- proplisttp1t2m3 <- proplisttp1tp1<- proplistt2m3C2 <- proplistt2m3t2m3 <- proplistt2m3tp1<- list()

proplistC2C2 <- lapply(1:length(proplistC2C2), function(x) {proplistC2C2[[x]]$without <- proplistC2C2[[1]]$without ; return(proplistC2C2[[x]])})
proplistC2t2m3 <- lapply(1:length(proplistC2t2m3), function(x) {proplistC2t2m3[[x]]$without <- proplistC2t2m3[[1]]$without ; return(proplistC2t2m3[[x]])})
proplistC2tp1 <- lapply(1:length(proplistC2tp1), function(x) {proplistC2tp1[[x]]$without <- proplistC2tp1[[1]]$without ; return(proplistC2tp1[[x]])})
proplisttp1C2 <- lapply(1:length(proplisttp1C2), function(x) {proplisttp1C2[[x]]$without <- proplisttp1C2[[1]]$without ; return(proplisttp1C2[[x]])})
proplisttp1t2m3 <- lapply(1:length(proplisttp1t2m3), function(x) {proplisttp1t2m3[[x]]$without <- proplisttp1t2m3[[1]]$without ; return(proplisttp1t2m3[[x]])})
proplisttp1tp1 <- lapply(1:length(proplisttp1tp1), function(x) {proplisttp1tp1[[x]]$without <- proplisttp1tp1[[1]]$without ; return(proplisttp1tp1[[x]])})
proplistt2m3C2 <- lapply(1:length(proplistt2m3C2), function(x) {proplistt2m3C2[[x]]$without <- proplistt2m3C2[[1]]$without ; return(proplistt2m3C2[[x]])})
proplistt2m3t2m3 <- lapply(1:length(proplistt2m3t2m3), function(x) {proplistt2m3t2m3[[x]]$without <- proplistt2m3t2m3[[1]]$without ; return(proplistt2m3t2m3[[x]])})
proplistt2m3tp1 <- lapply(1:length(proplistt2m3tp1), function(x) {proplistt2m3tp1[[x]]$without <- proplistt2m3tp1[[1]]$without ; return(proplistt2m3tp1[[x]])})


#propdiflist <- lapply(1:length(proplist), function(x) {proplist[[x]]$without <- proplist[[1]]$without ; return(proplist[[x]])})

climlistC2C2 <- lapply(proplistC2C2, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlistC2C2 <- lapply(proplistC2C2, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
withoutvectorC2C2 <- proplistC2C2[[1]]$without 
withoutclimC2C2 <- quantity2clim(withoutvectorC2C2, paste0(attr(withoutvectorC2C2, "probability")),ref.grid = sel.t2m)

climlistC2t2m3 <- lapply(proplistC2t2m3, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlistC2t2m3 <- lapply(proplistC2t2m3, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
withoutvectorC2t2m3 <- proplistC2t2m3[[1]]$without 
withoutclimC2t2m3 <- quantity2clim(withoutvectorC2t2m3, paste0(attr(withoutvectorC2t2m3, "probability")),ref.grid = sel.t2m)

climlistC2tp1 <- lapply(proplistC2tp1, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlistC2tp1 <- lapply(proplistC2tp1, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
withoutvectorC2tp1 <- proplistC2tp1[[1]]$without 
withoutclimC2tp1 <- quantity2clim(withoutvectorC2tp1, paste0(attr(withoutvectorC2tp1, "probability")),ref.grid = sel.t2m)

climlisttp1C2 <- lapply(proplisttp1C2, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlisttp1C2 <- lapply(proplisttp1C2, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
withoutvectortp1C2 <- proplisttp1C2[[1]]$without 
withoutclimtp1C2 <- quantity2clim(withoutvectortp1C2, paste0(attr(withoutvectortp1C2, "probability")),ref.grid = sel.t2m)

climlisttp1t2m3 <- lapply(proplisttp1t2m3, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlisttp1t2m3 <- lapply(proplisttp1t2m3, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
withoutvectortp1t2m3 <- proplisttp1t2m3[[1]]$without 
withoutclimtp1t2m3 <- quantity2clim(withoutvectortp1t2m3, paste0(attr(withoutvectortp1t2m3, "probability")),ref.grid = sel.t2m)

climlisttp1tp1 <- lapply(proplisttp1tp1, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlisttp1tp1 <- lapply(proplisttp1tp1, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
withoutvectortp1tp1 <- proplisttp1tp1[[1]]$without 
withoutclimtp1tp1 <- quantity2clim(withoutvectortp1tp1, paste0(attr(withoutvectortp1tp1, "probability")),ref.grid = sel.t2m)

climlistt2m3C2 <- lapply(proplistt2m3C2, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlistt2m3C2 <- lapply(proplistt2m3C2, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
withoutvectort2m3C2 <- proplistt2m3C2[[1]]$without 
withoutclimt2m3C2 <- quantity2clim(withoutvectort2m3C2, paste0(attr(withoutvectort2m3C2, "probability")),ref.grid = sel.t2m)

climlistt2m3t2m3 <- lapply(proplistt2m3t2m3, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlistt2m3t2m3 <- lapply(proplistt2m3t2m3, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
withoutvectort2m3t2m3 <- proplistt2m3t2m3[[1]]$without 
withoutclimt2m3t2m3 <- quantity2clim(withoutvectort2m3t2m3, paste0(attr(withoutvectort2m3t2m3, "probability")),ref.grid = sel.t2m)

climlistt2m3tp1 <- lapply(proplistt2m3tp1, function(x) quantity2clim(x$with - x$without, paste0(attr(x$with, "probability"),"-", attr(x$without, "probability")),ref.grid = sel.t2m))
climwithlistt2m3tp1 <- lapply(proplistt2m3tp1, function(x) quantity2clim(x$with, paste0(attr(x$with, "probability")),ref.grid = sel.t2m))
withoutvectort2m3tp1 <- proplistt2m3tp1[[1]]$without 
withoutclimt2m3tp1 <- quantity2clim(withoutvectort2m3tp1, paste0(attr(withoutvectort2m3tp1, "probability")),ref.grid = sel.t2m)

# #########################################################
# 
# enspropdifsclims <- bindGrid(climlist, dimension = c("member"))
# enspropdifsclims$Members <- cnodes
# enspropwithsclims <- bindGrid(climwithlist, dimension = c("member"))
# enspropwithsclims$Members <- cnodes
# #enspropwithoutsclims <- bindGrid(climwithoutlist, dimension = c("member"))
#enspropwithoutsclims$Members <- cnodes

#namesmembers <- apply(lvls,MARGIN = 1, function(x)paste0("lvl.t2m = ",x[1]," & lvl.tp = ",x[2]))
namesmembers <- cnodes

col.r <- colorRampPalette(brewer.pal(9, "Reds"))(15)
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(15)
col.p <- colorRampPalette(brewer.pal(9, "Purples"))(15)
col.g <- colorRampPalette(brewer.pal(9, "Greys"))(15)
col.gr <- colorRampPalette(brewer.pal(9, "Greens"))(15)
col.rb <- c(rev(col.b),'white','white',col.r)

# plotname <- paste0("figs/evidence_propagation/propdif_",Algo,"_",Score,"_",IT,"_t2m_tp_4lvls_t2m_",nE,"_",lvl.ev,".pdf")
# pdf(plotname, width = 10, height = 7)
# plotname <- paste0("figs/evidence_propagation/propdif_",Algo,"_",Score,"_",IT,"_t2m_tp_4lvls_t2m_",nE,"_",lvl.ev,".png")
# png(plotname, width = 20, height = 15,units = "cm", res = 150)


textdf <- cbind(attr(sel.t2m.1,"VertexCoords"),1:ngrids)
#textdf$x[textdf$x<0] <-textdf$x[textdf$x<0] + 360
textdf <- as.matrix(textdf)
textlay <- apply(textdf,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

clim.0vector <- quantity2clim(numeric(length=ngrids),"node vector", ref.grid = sel.t2m)
Xvectorname <-attr(clim.0vector$Data,"climatology:fun")
Xvector <-spatialPlot(clim.0vector, main = Xvectorname,  backdrop.theme = "coastline",sp.layout = textlay, lonCenter =0, colorkey = FALSE)

# KIJK IUT
# spatialPlot(enspropdifsclims,as.table = TRUE,names.attr = namesmembers,backdrop.theme = "coastline", lonCenter = 180, main = list(paste0("P(C = HD & V.tp = lvl.tp|",names(fitted)[nE]," = ",lvl.ev,") - P(V.t2m = lvl.t2m & V.tp = lvl.tp)")), 
#            at = seq(-0.4,0.4,0.025),
#            region = TRUE, 
#            col.regions= col.rb,
#            set.max = 0.4, 
#            rev.colors = TRUE,
#            colorkey = list(width = 0.6, lables = list(cex = 0.5)))



#dev.off()
# spatialPlot(withoutclim,as.table = TRUE,names.attr = namesmembers,backdrop.theme = "coastline", lonCenter = 0, 
#             main = list(paste0("P(C = HD)")), 
#             #at = seq(0,0.4,0.025),
#             region = TRUE, 
#             col.regions= col.r,
#             set.max = 0.4, 
#             rev.colors = TRUE,
#             colorkey = list(width = 0.6, lables = list(cex = 0.5)))

gridx <- 408
res <- 0.2
pC2C2 <- spatialPlot(climlistC2C2[[gridx]],
                   #as.table = TRUE,names.attr = namesmembers,
                   backdrop.theme = "coastline", lonCenter = 0, 
                   main = list(paste0("P(C = 2 (HotDry) |",names(fittedbn)[2*ngrids+gridx]," = ",lvl.comp," (HotDry)) - P(V.C = 2 (HotDry)")), 
                   at = seq(-res,res,res/16),
                   region = TRUE, 
                   col.regions= col.rb,
                   set.max = res, 
                   rev.colors = TRUE,
                   colorkey = list(width = 0.6, lables = list(cex = 0.5)))

pC2TP1 <-spatialPlot(climlistC2tp1[[gridx]],
                   #as.table = TRUE,names.attr = namesmembers,
                   backdrop.theme = "coastline", lonCenter = 0, 
                   main = list(paste0("P(V.tp = 1 (Dry)|",names(fittedbn)[2*ngrids+gridx]," = ",lvl.comp," (HotDry)) - P(V.tp = 1 (Dry)")), 
                   at = seq(-res,res,res/16),
                   region = TRUE, 
                   col.regions= col.rb,
                   set.max = res, 
                   rev.colors = TRUE,
                   colorkey = list(width = 0.6, lables = list(cex = 0.5)))
pC2T2M3 <- spatialPlot(climlistC2t2m3[[gridx]],
                     # as.table = TRUE,names.attr = namesmembers,
                     backdrop.theme = "coastline", 
                     lonCenter = 0, 
                     main = list(paste0("P(V.t2m = 3 (Hot)|",names(fittedbn)[2*ngrids+gridx]," = ",lvl.comp," (HotDry)) - P(V.t2m = 3 (Hot))")), 
                     at = seq(-res,res,res/16),
                     region = TRUE, 
                     col.regions= col.rb,
                     set.max = res, 
                     rev.colors = TRUE,
                     colorkey = list(width = 0.6, lables = list(cex = 0.5)))
pT2M3C2 <- spatialPlot(climlistt2m3C2[[gridx]],
                     #as.table = TRUE,names.attr = namesmembers,
                     backdrop.theme = "coastline", lonCenter = 0,
                     main = list(paste0("P(C = 2 (HotDry)|",names(fittedbn)[gridx]," = ",lvl.t2m3," (Hot)) - P(C2 = 2 (HotDry))")), 
                     at = seq(-res,res,res/16),
                     region = TRUE, 
                     col.regions= col.rb,
                     set.max = res, 
                     rev.colors = TRUE,
                     colorkey = list(width = 0.6, lables = list(cex = 0.5)))
pT2M3T2M3 <- spatialPlot(climlistt2m3t2m3[[gridx]],
                    #as.table = TRUE,names.attr = namesmembers,
                    backdrop.theme = "coastline", lonCenter = 0, 
                    main = list(paste0("P(V.t2m = 3 (Hot)|",names(fittedbn)[gridx]," = ",lvl.t2m3," (Hot)) - P(V.t2m = 3 (hot)")), 
                    at = seq(-res,res,res/16),
                    region = TRUE, 
                    col.regions= col.rb,
                    set.max = res, 
                    rev.colors = TRUE,
                    colorkey = list(width = 0.6, lables = list(cex = 0.5)))

pT2M3TP1 <- spatialPlot(climlistt2m3tp1[[gridx]],
                         #as.table = TRUE,names.attr = namesmembers,
                         backdrop.theme = "coastline", lonCenter = 0, 
                         main = list(paste0("P(V.tp = 1 (Dry)|",names(fittedbn)[gridx]," = ",lvl.t2m3," (Hot)) - P(V.tp = 1 (Dry))")), 
                         at = seq(-res,res,res/16),
                         region = TRUE, 
                         col.regions= col.rb,
                         set.max = res, 
                         rev.colors = TRUE,
                         colorkey = list(width = 0.6, lables = list(cex = 0.5)))

pTP1C2 <- spatialPlot(climlisttp1C2[[gridx]],
                       #as.table = TRUE,names.attr = namesmembers,
                       backdrop.theme = "coastline", lonCenter = 0,
                       main = list(paste0("P(C = 2 (HotDry)|",names(fittedbn)[ngrids+ gridx]," = ",lvl.tp1," (Dry)) - P(C2 = 2 (HotDry)) ")), 
                       at = seq(-res,res,res/16),
                       region = TRUE, 
                       col.regions= col.rb,
                       set.max = res, 
                       rev.colors = TRUE,
                       colorkey = list(width = 0.6, lables = list(cex = 0.5)))
pTP1T2M3 <- spatialPlot(climlisttp1t2m3[[gridx]],
                         #as.table = TRUE,names.attr = namesmembers,
                         backdrop.theme = "coastline", lonCenter = 0, 
                         main = list(paste0("P(V.t2m = 3 (Hot)|",names(fittedbn)[ngrids +gridx]," = ",lvl.tp1," (Dry)) - P(V.t2m = 3 (hot)")), 
                         at = seq(-res,res,res/16),
                         region = TRUE, 
                         col.regions= col.rb,
                         set.max = res, 
                         rev.colors = TRUE,
                         colorkey = list(width = 0.6, lables = list(cex = 0.5)))

pTP1TP1 <- spatialPlot(climlisttp1tp1[[gridx]],
                        #as.table = TRUE,names.attr = namesmembers,
                        backdrop.theme = "coastline", lonCenter = 0, 
                        main = list(paste0("P(V.tp = 1 (Dry)|",names(fittedbn)[ngrids +gridx]," = ",lvl.tp1," (Dry)) - P(V.tp = 1 (Dry))")), 
                        at = seq(-res,res,res/16),
                        region = TRUE, 
                        col.regions= col.rb,
                        set.max = res, 
                        rev.colors = TRUE,
                        colorkey = list(width = 0.6, lables = list(cex = 0.5)))


plotname <- paste0("figs/grid_impact/H_D_HD_impact_simple_",IT,"_res",res,"_C",gridx)
pdf(paste0(plotname,".pdf"),width = 20, height = 10)
plotlist<- list(pC2C2,pT2M3C2,pTP1C2,pC2TP1,pT2M3TP1,pTP1TP1,pC2T2M3,pT2M3T2M3,pTP1T2M3 )
do.call(grid.arrange,c(plotlist,nrow = 3))
dev.off()

plotname <- paste0("figs/grid_impact/H_D_HD_impact_simple_",IT,"_res",res,"_C",gridx,".png")
png(plotname, width = 50, height = 20,units = "cm", res = 300)
plotlist<- list(pC2C2,pT2M3C2,pTP1C2,pC2TP1,pT2M3TP1,pTP1TP1,pC2T2M3,pT2M3T2M3,pTP1T2M3 )
do.call(grid.arrange,c(plotlist,nrow = 3))
dev.off()

# plotname <- paste0("figs/HD_impact/HD_impact_all_C",gridx)
# pdf(paste0(plotname,".pdf"),width = 20, height = 15)
# plotlist<- list(pC2,Xvector,pT2M3,pTP1,pT2M1,pTP3)
# do.call(grid.arrange,c(plotlist,nrow = 3))
# dev.off()


# Hier bezig
# 
# plotname <- paste0("figs/HD_impact/HD_impact_to_HD_C",gridx)
# png(paste0(plotname,".png"),width = 20, height = 5, units = "cm", res = 300)
# 
# spatialPlot(climwithlist[[gridx]],
#             #as.table = TRUE,names.attr = namesmembers,
#             backdrop.theme = "coastline", lonCenter = 0, 
#             main = list(paste0("P(V.Comp.HD = 2 (i.e. HotDry = TRUE) |",names(fittedbn)[2*ngrids+gridx]," = ",lvl.ev," (HotDry))")), 
#             at = seq(0,0.4,0.025),
#             region = TRUE, 
#             col.regions= c("white",col.r),
#             set.max = 0.4, 
#             rev.colors = TRUE,
#             colorkey = list(width = 0.6, lables = list(cex = 0.5)))
# 
# dev.off()


#######################################################################
# Sensitive nodes analysis
#######################################################################
posdiffun <- function(x) {if (x >= 0) {x} else {0}}

attr(sel.tp.1,"VertexCoords")
haversine
griddists <- apply(attr(sel.tp.1,"VertexCoords"),MARGIN = 1, FUN = function(x) apply(attr(sel.tp.1,"VertexCoords"),MARGIN = 1, FUN = function(y) haversine(x[1],x[2],y[1],y[2])))
griddists[1,]%*%HDtoHDdif[,1]/4683530
HDtoHDdif[,2]%*%griddists[2,]
sum(HDtoHDdif[,2]*griddists[2,])
HDtoHDdif[3,2]*griddists[3,2]
proplistC2C2[[2]]

p.ch <- 0.1

HDtoHDextension <- sapply(proplistC2C2, function(x){sum(abs(x$with-x$without)>= 0.05)})
HDtoHDextension_a <- sapply(proplistC2C2, function(x){a <- x$with-x$without; a[a<= 0.05]<- 0; a[a>= 0.05]<- 1;return(a)})
HDtoHDimpactweighted <-diag(griddists%*%HDtoHDextension_a)/HDtoHDextension

HDtoHDdif <- sapply(proplistC2C2,function(x,y){abs(x$with-x$without)})
HDtoHDimpactweighted <- diag(griddists%*%HDtoHDdif)/nrow(griddists)
HDtoHDimpactwclim<-quantity2clim(HDtoHDimpactweighted, "HD-2-HD weighted impact",ref.grid = sel.t2m)

# HDtoHDdif%*%griddists
# 
# griddists[5,5]
# max(griddists)

HDtoHDimpact <- sapply(proplistC2C2, function(x){mean(abs(x$with-x$without))})
HDtoHDimpactclim <- quantity2clim(HDtoHDimpact, "HD-2-HD impact",ref.grid = sel.t2m)
HDtoHDname <-attr(HDtoHDimpactclim$Data,"climatology:fun")

assign("spatextHDtoHD",HDtoHDextension)
save(spatextHDtoHD, file = paste0("results/propagation/spatext_",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/spatext_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_C_",lvl.ev,".rda")
)

pi*6300*2/2

textHD <- HDtoHDextension
textHD <- round(HDtoHDimpact*1000,0)
textdfHD <- cbind(attr(sel.t2m.1,"VertexCoords"),textHD)
#textdf$x[textdf$x<0] <- textdf$x[textdf$x<0] + 360
textdfHD <- as.matrix(textdfHD)
textlayHD <- apply(textdfHD,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

HDtoHDimpactclim$Data
HD <- spatialPlot(HDtoHDimpactclim,backdrop.theme = "coastline", lonCenter = 0, 
                  # sp.layout = textlayHD,
                  rev.colors = FALSE,
                  main = HDtoHDname,
                  regions = TRUE, col.regions = col.p,
                  at = seq(0,max(HDtoHDimpact),max(HDtoHDimpact)/15),
                  set.max = max(HDtoHDimpact),
                  colorkey = list(width = 0.6, lables = list(cex = 1)))
HDw <- spatialPlot(HDtoHDimpactwclim,backdrop.theme = "coastline", lonCenter = 0, 
                  # sp.layout = textlayHD,
                  rev.colors = FALSE,
                  main = HDtoHDname,
                  regions = TRUE, col.regions = col.p,
                  at = seq(0,max(HDtoHDimpactweighted),max(HDtoHDimpactweighted)/15),
                  set.max = max(HDtoHDimpactweighted),
                  colorkey = list(width = 0.6, lables = list(cex = 1)))
grid.arrange(HDw,HD)
HDw

dev.off()

spatialPlot(HDtoHDextension,backdrop.theme = "coastline", lonCenter = 0, 
            sp.layout = textlayHD,
            rev.colors = FALSE,
            main = HDtoHDname,
            regions = TRUE, col.regions = col.p,
            # at = seq(0,0.060,0.060/15),
            # set.max = 0.060,
            colorkey = list(width = 0.6, lables = list(cex = 1)))


HDtoHextension <- sapply(proplistC2t2m3, function(x){sum(abs(x$with-x$without)>= 0.05)})
HDtoHimpact <- sapply(proplistC2t2m3, function(x){mean(abs(x$with-x$without))})
HDtoHimpactclim <- quantity2clim(HDtoHimpact, "HD-2-H impact",ref.grid = sel.t2m)
HDtoHname <-attr(HDtoHimpactclim$Data,"climatology:fun")



textH <- HDtoHextension
textH <- round(HDtoHimpact*1000,0)
textdfH <- cbind(attr(sel.t2m.1,"VertexCoords"),textH)
#textdf$x[textdf$x<0] <- textdf$x[textdf$x<0] + 360
textdfH <- as.matrix(textdfH)
textlayH <- apply(textdfH,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

H <- spatialPlot(HDtoHimpactclim,backdrop.theme = "coastline", lonCenter = 0, sp.layout = textlayH,
                 rev.colors = FALSE, main = HDtoHname,
                 regions = TRUE, col.regions = col.r,
                 at = seq(0,0.16,0.16/15),
                  set.max = 0.16,
                 colorkey = list(width = 0.6, lables = list(cex = 1)))
H

HDtoDextension <- sapply(proplistC2tp1, function(x){sum(abs(x$with-x$without)>= 0.05)})

HDtoDimpact <- sapply(proplistC2tp1, function(x){mean(abs(x$with-x$without))})
HDtoDimpactclim <- quantity2clim(HDtoDimpact, "HD-2-D impact",ref.grid = sel.t2m)
HDtoDname <-attr(HDtoDimpactclim$Data,"climatology:fun")

textD <- HDtoDextension
textD <- round(HDtoDimpact*1000,0)
textdfD <- cbind(attr(sel.t2m.1,"VertexCoords"),textD)
#textdf$x[textdf$x<0] <- textdf$x[textdf$x<0] + 360
textdfD <- as.matrix(textdfD)
textlayD <- apply(textdfD,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

D <-spatialPlot(HDtoDimpactclim,backdrop.theme = "coastline", lonCenter = 0, sp.layout = textlayD,
                rev.colors = FALSE, main = HDtoDname,
                regions = TRUE, col.regions = col.g,
                at = seq(0,0.06,0.06/15),
               # set.max = 0.04,
                colorkey = list(width = 0.6, lables = list(cex = 1)))
D

HtoHDdif <- sapply(proplistt2m3C2,function(x,y){abs(x$with-x$without)})
HtoHDimpactweighted <- diag(griddists%*%HtoHDdif)/colSums(griddists)
HtoHDimpactwclim<-quantity2clim(HtoHDimpactweighted, "H-2-HD weighted impact",ref.grid = sel.t2m)


HtoHDextension <- sapply(proplistt2m3C2, function(x){sum(abs(x$with-x$without)>= 0.05)})
HtoHDextension_a <- sapply(proplistt2m3C2, function(x){a <- x$with-x$without; a[a<= 0.05]<- 0; a[a>= 0.05]<- 1;return(a)})
HtoHDimpactweighted <-diag(griddists%*%HtoHDextension_a)/HtoHDextension
HtoHDimpactwclim<-quantity2clim(HtoHDimpactweighted, "H-2-HD spatial extension",ref.grid = sel.t2m)
HtoHDwname <- attr(HtoHDimpactwclim$Data,"climatology:fun")



HtoHDimpact <- sapply(proplistt2m3C2, function(x){mean(abs(x$with-x$without))})
HtoHDimpactclim <- quantity2clim(HtoHDimpact, "H-2-HD impact",ref.grid = sel.t2m)
HtoHDname <-attr(HtoHDimpactclim$Data,"climatology:fun")

assign("spatextHtoHD",HtoHDextension)
save(spatextHtoHD, file = paste0("results/propagation/spatext_",Algo,"_disc_",b,"_",Score,"_lim_comp_art_",IT,"/spatext_",Algo,"_",Score,"_",IT,"_C_",lvl.comp,"_t2m_",lvl.t2m3,".rda")
)


textHtoHD <- HtoHDextension
textHtoHD <- round(HtoHDimpact*1000,0)
textdfHtoHD <- cbind(attr(sel.t2m.1,"VertexCoords"),textHtoHD)
#textdf$x[textdf$x<0] <- textdf$x[textdf$x<0] + 360
textdfHtoHD <- as.matrix(textdfHtoHD)
textlayHtoHD <- apply(textdfHtoHD,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

HtoHDimpactclim$Data
HtoHD <- spatialPlot(HtoHDimpactclim,backdrop.theme = "coastline", lonCenter = 0, 
                  # sp.layout = textlayHtoHD,
                  rev.colors = FALSE,
                  main = HtoHDname,
                  regions = TRUE, col.regions = col.p,
                  at = seq(0,0.060,0.060/15),
                  set.max = 0.060,
                  colorkey = list(width = 0.6, lables = list(cex = 1)))
HtoHDw <- spatialPlot(HtoHDimpactwclim,backdrop.theme = "coastline", lonCenter = 0, 
                               # sp.layout = textlayHtoHD,
                               rev.colors = FALSE,
                               main = HtoHDwname,
                               regions = TRUE, col.regions = col.p,
                                at = seq(0,max(HtoHDimpactweighted),max(HtoHDimpactweighted)/15),
                               # at = seq(0,0.060,0.060/15),
                               #  set.max = 0.060,         
                               #  #set.max = 0.060,
                               colorkey = list(width = 0.6, lables = list(cex = 1)))
HtoHD
HtoHDw
grid.arrange(HtoHD,HtoHDw)


HtoHextension <- sapply(proplistt2m3t2m3, function(x){sum(abs(x$with-x$without)>= 0.05)})
HtoHimpact <- sapply(proplistt2m3t2m3, function(x){mean(abs(x$with-x$without))})
HtoHimpactclim <- quantity2clim(HtoHimpact, "H-2-H impact",ref.grid = sel.t2m)
HtoHname <-attr(HtoHimpactclim$Data,"climatology:fun")



textHtoH <- HtoHextension
textHtoH <- round(HtoHimpact*1000,0)
textdfHtoH <- cbind(attr(sel.t2m.1,"VertexCoords"),textHtoH)
#textdf$x[textdf$x<0] <- textdf$x[textdf$x<0] + 360
textdfHtoH <- as.matrix(textdfHtoH)
textlayHtoH <- apply(textdfHtoH,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

HtoH <- spatialPlot(HtoHimpactclim,backdrop.theme = "coastline", lonCenter = 0, sp.layout = textlayHtoH,
                 rev.colors = FALSE, main = HtoHname,
                 regions = TRUE, col.regions = col.r,
                 at = seq(0,0.16,0.16/15),
                 set.max = 0.16,
                 colorkey = list(width = 0.6, lables = list(cex = 1)))
HtoH




HtoDextension <- sapply(proplistt2m3tp1, function(x){sum(abs(x$with-x$without)>= 0.05)})

HtoDimpact <- sapply(proplistt2m3tp1, function(x){mean(abs(x$with-x$without))})
HtoDimpactclim <- quantity2clim(HtoDimpact, "H-2-D impact",ref.grid = sel.t2m)
HtoDname <-attr(HtoDimpactclim$Data,"climatology:fun")

textHtoD <- HtoDextension
textHtoD <- round(HtoDimpact*1000,0)
textdfHtoD <- cbind(attr(sel.t2m.1,"VertexCoords"),textHtoD)
#textdf$x[textdf$x<0] <- textdf$x[textdf$x<0] + 360
textdfHtoD <- as.matrix(textdfHtoD)
textlayHtoD <- apply(textdfHtoD,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

HtoD <-spatialPlot(HtoDimpactclim,backdrop.theme = "coastline", lonCenter = 0, sp.layout = textlayHtoD,
                rev.colors = FALSE, main = HtoDname,
                regions = TRUE, col.regions = col.g,
                at = seq(0,0.06,0.06/15),
                # set.max = 0.04,
                colorkey = list(width = 0.6, lables = list(cex = 1)))
HtoD
D

DtoHDdif <- sapply(proplisttp1C2,function(x,y){abs(x$with-x$without)})
DtoHDimpactweighted <- colMeans(griddists%*%DtoHDdif)
DtoHDimpactwclim<-quantity2clim(DtoHDimpactweighted, "D-2-HD weighted impact",ref.grid = sel.t2m)



DtoHDextension <- sapply(proplisttp1C2, function(x){sum(abs(x$with-x$without)>= 0.1)})
DtoHDimpact <- sapply(proplisttp1C2, function(x){mean(abs(x$with-x$without))})
DtoHDimpactclim <- quantity2clim(DtoHDimpact, "D-2-HD impact",ref.grid = sel.t2m)
DtoHDname <-attr(DtoHDimpactclim$Data,"climatology:fun")

textDtoHD <- DtoHDextension
textDtoHD <- round(DtoHDimpact*1000,0)
textdfDtoHD <- cbind(attr(sel.t2m.1,"VertexCoords"),textDtoHD)
#textdf$x[textdf$x<0] <- textdf$x[textdf$x<0] + 360
textdfDtoHD <- as.matrix(textdfDtoHD)
textlayDtoHD <- apply(textdfDtoHD,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

DtoHD <-spatialPlot(DtoHDimpactclim,backdrop.theme = "coastline", lonCenter = 0, sp.layout = textlayDtoHD,
                   rev.colors = FALSE, main = DtoHDname,
                   regions = TRUE, col.regions = col.p,
                   at = seq(0,0.060,0.060/15),
                   set.max = 0.060,
                   colorkey = list(width = 0.6, lables = list(cex = 1)))
DtoHDw <-spatialPlot(DtoHDimpactwclim,backdrop.theme = "coastline", lonCenter = 0, sp.layout = textlayDtoHD,
                        rev.colors = FALSE, main = DtoHDname,
                        regions = TRUE, col.regions = col.p,
                        # at = seq(0,40,40/15),
                        #set.max = 0.060,
                        colorkey = list(width = 0.6, lables = list(cex = 1)))

DtoHD
DtoHDw

DtoHextension <- sapply(proplisttp1t2m3, function(x){sum(abs(x$with-x$without)>= 0.05)})

DtoHimpact <- sapply(proplisttp1t2m3, function(x){mean(abs(x$with-x$without))})
DtoHimpactclim <- quantity2clim(DtoHimpact, "D-2-H impact",ref.grid = sel.t2m)
DtoHname <-attr(DtoHimpactclim$Data,"climatology:fun")

textDtoH <- DtoHextension
textDtoH <- round(DtoHimpact*1000,0)
textdfDtoH <- cbind(attr(sel.t2m.1,"VertexCoords"),textDtoH)
#textdf$x[textdf$x<0] <- textdf$x[textdf$x<0] + 360
textdfDtoH <- as.matrix(textdfDtoH)
textlayDtoH <- apply(textdfDtoH,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

DtoH <-spatialPlot(DtoHimpactclim,backdrop.theme = "coastline", lonCenter = 0, sp.layout = textlayDtoH,
                   rev.colors = FALSE, main = DtoHname,
                   regions = TRUE, col.regions = col.r,
                   at = seq(0,0.16,0.16/15),
                   set.max = 0.16,
                   colorkey = list(width = 0.6, lables = list(cex = 1)))
DtoH
H

DtoDextension <- sapply(proplisttp1tp1, function(x){sum(abs(x$with-x$without)>= 0.05)})

DtoDimpact <- sapply(proplisttp1tp1, function(x){mean(abs(x$with-x$without))})
DtoDimpactclim <- quantity2clim(DtoDimpact, "D-2-D impact",ref.grid = sel.t2m)
DtoDname <-attr(DtoDimpactclim$Data,"climatology:fun")

textDtoD <- DtoDextension
textDtoD <- round(DtoDimpact*1000,0)
textdfDtoD <- cbind(attr(sel.t2m.1,"VertexCoords"),textDtoD)
#textdf$x[textdf$x<0] <- textdf$x[textdf$x<0] + 360
textdfDtoD <- as.matrix(textdfDtoD)
textlayDtoD <- apply(textdfDtoD,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

DtoD <-spatialPlot(DtoDimpactclim,backdrop.theme = "coastline", lonCenter = 0, sp.layout = textlayDtoD,
                rev.colors = FALSE, main = DtoDname,
                regions = TRUE, col.regions = col.g,
                at = seq(0,0.06,0.06/15),
                 set.max = 0.05,
                colorkey = list(width = 0.6, lables = list(cex = 1)))
DtoD
D

plotname <- paste0("figs/H_impact/H_impact_simple_",IT)
pdf(paste0(plotname,".pdf"),width = 20, height = 10)
png(paste0(plotname,".png"),units = "cm",width = 40, height = 40,res = 600)
grid.arrange(HtoHD,Xvector,HtoH,HtoD)
dev.off()

plotname <- paste0("figs/D_impact/D_impact_simple_",IT)
pdf(paste0(plotname,".pdf"),width = 20, height = 10)
png(paste0(plotname,".png"),units = "cm",width = 40, height = 40,res = 600)
grid.arrange(DtoHD,Xvector,DtoH,DtoD)
dev.off()

grid.arrange(HD,HtoHD,DtoHD)
grid.arrange(D,HtoD,DtoD)
grid.arrange(H,HtoH,DtoH)

plotname<- paste0("figs/all_var_impact/all_var_impact_",IT)
pdf(paste0(plotname,".pdf"),width = 20, height = 10)
grid.arrange(HD,HtoHD,DtoHD,D,HtoD,DtoD,H,HtoH,DtoH)
dev.off()

plotname<- paste0("figs/all_var_impact/all_var_impact_",IT,"_weighted")
pdf(paste0(plotname,".pdf"),width = 20, height = 10)
grid.arrange(HD,HDw,HtoHD,HtoHDw,DtoHD,DtoHDw)
dev.off()

##################################################
# Plot specific grid and all together
##################################################

plotname <- paste0("figs/all_var_impact/all_var_impact_",IT,"_V",gridx)
pdf(paste0(plotname,".pdf"),width = 13, height = 30)
plotlist<- list(pC2C2,HD,pT2M3C2,HtoHD,pTP1C2,DtoHD,pC2TP1,D,pT2M3TP1,HtoD,pTP1TP1,DtoD,pC2T2M3,H,pT2M3T2M3,HtoH,pTP1T2M3,DtoH)
do.call(grid.arrange,c(plotlist,ncol = 2 ))
dev.off()

#################################################
# Plot spatial extension and impact (without text numbers)
#################################################
attr(sel.tp.1,"VertexCoords")
haversine
griddists <- apply(attr(sel.tp.1,"VertexCoords"),MARGIN = 1, FUN = function(x) apply(attr(sel.tp.1,"VertexCoords"),MARGIN = 1, FUN = function(y) haversine(x[1],x[2],y[1],y[2])))
griddists[1,]%*%HDtoHDdif[,1]/4683530
HDtoHDdif[,2]%*%griddists[2,]
sum(HDtoHDdif[,2]*griddists[2,])
HDtoHDdif[3,2]*griddists[3,2]
proplistC2C2[[2]]

p.ch <- 0.1

# number of variables affected by prob change more than 0.05 counted for each gridbox
HDtoHDextension <- sapply(proplistC2C2, function(x){sum(abs(x$with-x$without)>= p.ch)})
# afected medium distance of minumal changed variables.
HDtoHDextension_a <- sapply(proplistC2C2, function(x){a <- x$with-x$without; a[a<= p.ch]<- 0; a[a>= p.ch]<- 1;return(a)})
HDtoHDimpactweighted <-diag(griddists%*%HDtoHDextension_a)/HDtoHDextension
# HDtoHDdif <- sapply(proplistC2C2,function(x,y){abs(x$with-x$without)})
# HDtoHDimpactweighted <- diag(griddists%*%HDtoHDdif)/nrow(griddists)
HDtoHDimpactwclim <- quantity2clim(HDtoHDimpactweighted, "HD-2-HD spatial extension",ref.grid = sel.t2m)
HDtoHDwname <- attr(HDtoHDimpactwclim$Data,"climatology:fun")

HDtoHDimpact <- sapply(proplistC2C2, function(x){mean(abs(x$with-x$without))})
HDtoHDimpactclim <- quantity2clim(HDtoHDimpact, "HD-2-HD impact",ref.grid = sel.t2m)
HDtoHDname <-attr(HDtoHDimpactclim$Data,"climatology:fun")

HD <- spatialPlot(HDtoHDimpactclim,backdrop.theme = "coastline", lonCenter = 0, 
                  # sp.layout = textlayHD,
                  rev.colors = FALSE,
                  main = HDtoHDname,
                  regions = TRUE, col.regions = col.p,
                  at = seq(0,max(HDtoHDimpact),max(HDtoHDimpact)/15),
                  set.max = max(HDtoHDimpact),
                  colorkey = list(width = 0.6, lables = list(cex = 1)))
HDw <- spatialPlot(HDtoHDimpactwclim,backdrop.theme = "coastline", lonCenter = 0, 
                   # sp.layout = textlayHD,
                   rev.colors = FALSE,
                   main = HDtoHDwname,
                   regions = TRUE, col.regions = col.g,
                   at = seq(0,max(HDtoHDimpactweighted),max(HDtoHDimpactweighted)/15),
                   set.max = max(HDtoHDimpactweighted),
                   colorkey = list(width = 0.6, lables = list(cex = 1)))


HtoHDdif <- sapply(proplistt2m3C2,function(x,y){abs(x$with-x$without)})
# HtoHDimpactweighted <- diag(griddists%*%HtoHDdif)/colSums(griddists)
# HtoHDimpactwclim<-quantity2clim(HtoHDimpactweighted, "H-2-HD weighted impact",ref.grid = sel.t2m)


HtoHDextension <- sapply(proplistt2m3C2, function(x){sum(abs(x$with-x$without)>= p.ch)})
HtoHDextension_a <- sapply(proplistt2m3C2, function(x){a <- x$with-x$without; a[a<= p.ch]<- 0; a[a>= p.ch]<- 1;return(a)})
HtoHDimpactweighted <-diag(griddists%*%HtoHDextension_a)/HtoHDextension
HtoHDimpactwclim<-quantity2clim(HtoHDimpactweighted, "H-2-HD spatial extension",ref.grid = sel.t2m)
HtoHDwname <- attr(HtoHDimpactwclim$Data,"climatology:fun")



HtoHDimpact <- sapply(proplistt2m3C2, function(x){mean(abs(x$with-x$without))})
HtoHDimpactclim <- quantity2clim(HtoHDimpact, "H-2-HD impact",ref.grid = sel.t2m)
HtoHDname <-attr(HtoHDimpactclim$Data,"climatology:fun")

HtoHD <- spatialPlot(HtoHDimpactclim,backdrop.theme = "coastline", lonCenter = 0, 
                     # sp.layout = textlayHtoHD,
                     rev.colors = FALSE,
                     main = HtoHDname,
                     regions = TRUE, col.regions = col.p,
                     at = seq(0,max(HDtoHDimpact),max(HDtoHDimpact)/15),
                     set.max = max(HDtoHDimpact),
                     # at = seq(0,0.060,0.060/15),
                     # set.max = 0.060,
                     colorkey = list(width = 0.6, lables = list(cex = 1)))
HtoHDw <- spatialPlot(HtoHDimpactwclim,backdrop.theme = "coastline", lonCenter = 0, 
                      # sp.layout = textlayHtoHD,
                      rev.colors = FALSE,
                      main = HtoHDwname,
                      regions = TRUE, col.regions = col.g,
                      at = seq(0,max(HDtoHDimpactweighted),max(HDtoHDimpactweighted)/15),
                      set.max = max(HDtoHDimpactweighted),
                     
                      #at = seq(0,max(HtoHDimpactweighted),max(HtoHDimpactweighted)/15),
                      # at = seq(0,0.060,0.060/15),
                      #  set.max = 0.060,         
                      #  #set.max = 0.060,
                      colorkey = list(width = 0.6, lables = list(cex = 1)))


DtoHDdif <- sapply(proplisttp1C2,function(x,y){abs(x$with-x$without)})
DtoHDextension <- sapply(proplisttp1C2, function(x){sum(abs(x$with-x$without)>= p.ch)})
# afected medium distance of minumal changed variables.
DtoHDextension_a <- sapply(proplisttp1C2, function(x){a <- x$with-x$without; a[a<= p.ch]<- 0; a[a>= p.ch]<- 1;return(a)})
DtoHDimpactweighted <-diag(griddists%*%DtoHDextension_a)/DtoHDextension
# HDtoHDdif <- sapply(proplistC2C2,function(x,y){abs(x$with-x$without)})
# HDtoHDimpactweighted <- diag(griddists%*%HDtoHDdif)/nrow(griddists)
DtoHDimpactwclim <- quantity2clim(DtoHDimpactweighted, "D-2-HD spatial extension",ref.grid = sel.t2m)
DtoHDwname <- attr(DtoHDimpactwclim$Data,"climatology:fun")



DtoHDextension <- sapply(proplisttp1C2, function(x){sum(abs(x$with-x$without)>= 0.1)})
DtoHDimpact <- sapply(proplisttp1C2, function(x){mean(abs(x$with-x$without))})
DtoHDimpactclim <- quantity2clim(DtoHDimpact, "D-2-HD impact",ref.grid = sel.t2m)
DtoHDname <-attr(DtoHDimpactclim$Data,"climatology:fun")

DtoHD <-spatialPlot(DtoHDimpactclim,backdrop.theme = "coastline", lonCenter = 0, 
                    # sp.layout = textlayDtoHD,
                    rev.colors = FALSE, main = DtoHDname,
                    regions = TRUE, col.regions = col.p,
                    # at = seq(0,0.060,0.060/15),
                    # set.max = 0.060,
                    at = seq(0,max(HDtoHDimpact),max(HDtoHDimpact)/15),
                    set.max = max(HDtoHDimpact),
                    colorkey = list(width = 0.6, lables = list(cex = 1)))
DtoHDw <-spatialPlot(DtoHDimpactwclim,backdrop.theme = "coastline", lonCenter = 0, 
                     # sp.layout = textlayDtoHD,
                     rev.colors = FALSE, main = DtoHDwname,
                     regions = TRUE, col.regions = col.g,
                     at = seq(0,max(HDtoHDimpactweighted),max(HDtoHDimpactweighted)/15),
                     set.max = max(HDtoHDimpactweighted),
                     #at = seq(0,40,40/15),
                     #set.max = 0.060,
                     colorkey = list(width = 0.6, lables = list(cex = 1)))

plotname <- paste0("figs/all_var_impact/all_var_impact_spatext_",p.ch,"_toHD_",IT)
pdf(paste0(plotname,".pdf"),width = 13, height = 10)
grid.arrange(HD,HDw,HtoHD,HtoHDw,DtoHD,DtoHDw)
dev.off()

which(griddists == max(griddists),arr.ind =TRUE)
pi*6371*2

