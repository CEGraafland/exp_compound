##################################################################################################
# Event apreciation
##################################################################################################
rm(list = ls())
setwd("~/data/Untitled/Trabajo/R_practice/exp_compound")
library(bnlearn)
library(magrittr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(parallel)
library(magrittr)
library(visualizeR)
library('RColorBrewer')
library(gridExtra)

source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")

########################################################################################################
# 
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/propagationFunctions.R")

degrees <- "10d"
mask <- FALSE # Only land or land & sea
limit <- TRUE # include poles?
b <- 3  # How many levels TP & T2M:
k <- 0 # Which permutation of variables later
# IT <- paste0(k,"_2000_2100")
# IT <- paste0(k,"_2100_2200")
# # IT <- paste0(k,"_1000_1100")
# IT <- paste0(k,"_1700_1800")
compound <- "HD" # Which type of compound Hot Dry, Cold Dry, Hot Wet, Cold Wet, all together, all.distinct
global_mean_temp <- TRUE
global_mean_temp_cont <- FALSE# Extra node for global_mean_temp?
combi <- "TPT2M" # "only" Compound nodes or also "TPT2M" nodes
# algo<- "tabu" # Which algorithm
# score<- "bic" # Which score in algorithm
# global_mean_temp_art <- FALSE
# cut.art <- ""
 scaleGT <- FALSE
# geo <- ""

if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
if(!is.null(compound)){comp <- paste0("_comp_",compound)} else {comp <-""}
if(isTRUE(mask)){superf <- "_land"}else {superf <-"_landsea"}
if(isTRUE(global_mean_temp)){n1 <- "_GMT"} else if(isTRUE(global_mean_temp_cont)){n1 <- "_GMTcont"} else {n1 <- ""}
# if(isTRUE(global_mean_temp_art)){n2 <- paste0("_GMTart_",cut.art,geo)} else {n2 <- ""}

# if(isTRUE(global_mean_temp_cont) & !isTRUE(global_mean_temp_art) & score == "bic"){score <- "bic-cg"} 
# if(isTRUE(global_mean_temp_cont) & !isTRUE(global_mean_temp_art) & score == "aic"){score <- "aic-cg"} 
# if (isTRUE(global_mean_temp_art)){n <- n2} else {n <- n1}

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
 
 
if(!k==0){permutations <-get(load(paste0("data/permutations/permutations_",degrees,lim,superf,".rda")))}

data.k <- get(load(file =paste0("data/manipulate_aggr/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".rda")))
data <- data.k
if(scaleGT == TRUE) {data$GT <- scale(data$GT)}
if(!k==0){data <-data[,c(permutations[[k]],(length(permutations[[k]])+1):ncol(data))]}


t2m_ERA5 <- get(load(paste0("data/raw_aggr/t2m_ERA5_monthly_1940_2022_",degrees,".rda")))
tp_ERA5 <- get(load(paste0("data/raw_aggr/tp_ERA5_monthly_1940_2022_",degrees,".rda")))
sel.tp <- subsetGrid(tp_ERA5,latLim = c(-50,75))
sel.t2m <- subsetGrid( t2m_ERA5,latLim = c(-50,75))

#####################################################################
# Compound & GT discrete & lvls 1:4
#####################################################################
dataRMS2 <- data
lvls <- 1:4
dif.extremes.list <- list()
for(i in lvls){
Big2ind <-  which(dataRMS2$GT == levels(dataRMS2$GT)[i])
c <- sapply(dataRMS2[Big2ind,], FUN = as.numeric)
d <- sapply(dataRMS2, FUN = as.numeric)
extremes.C <- apply(c[,grep("C",colnames(c))],MARGIN = 2,FUN = function(x)sum(x==2)/length(x))
we.extremes.C <- apply(d[,grep("C",colnames(c))],MARGIN = 2,FUN = function(x)sum(x==2)/length(x))
dif.extremes.C <- extremes.C - we.extremes.C
dif.extremes.list[[i]]<- dif.extremes.C
}



col.blue <- rev(brewer.pal(9,"Blues"))
col.blue
col.red <- brewer.pal(9,"Reds")

col.div <- c(col.blue, col.red)
col.div[4:5] <- "white"
col.b <- c(col.blue,col.red)
length(col.b)

col.r <- colorRampPalette(brewer.pal(9, "Reds"))(15)
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(15)
col.rb <- c(rev(col.b),'white','white',col.r)

if(mask == TRUE){
mask_loc <- as.vector(array3Dto2Dmat(redim(landmask,member = FALSE)$Data))
list1 <- list2 <-dif.extremes.list
for(i in 1:length(list1)){
  list1[[i]][length(list1[[i]])+1] <- NA}
mask_loc_perm <- mask_loc
mask_loc_perm[!is.na(mask_loc)]<- 1:length(list2[[1]])
#mask_loc_perm[is.na(mask_loc_perm)] <- j
list3 <-list1
for(i in 1:length(list1)){
  list3[[i]] <-list1[[i]][mask_loc_perm]}
} else {list3 <- dif.extremes.list}

attributes(sel.t2m$Data)
climBig2 <- lapply(list3,quantity2clim,what = "Big2",ref.grid = sel.t2m)
namesmembers <- paste0("Rel.Freq q = ",lvls)

climBig2.ens <- makeMultiGrid(climBig2)
climBig2.ens$Members <- namesmembers
figdata <- spatialPlot(climBig2.ens,as.table = TRUE,names.attr = namesmembers,lonCenter = 0, backdrop.theme = "coastline", rev.colors = TRUE, 
                       main = paste0("relative frequency ",compound),
                       col.regions = col.rb,set.min = -0.2, set.max = 0.2,
                       at = seq(-0.2,0.2,0.4/32), cex = 0.5,
                       colorkey = list(width = 0.6, lables = list(cex = 0.5)))
figdata

plotname <- paste0("figs/dataevents/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,"_quantiles.pdf")
pdf(plotname, width = 10, height = 5)
plotname <- paste0("figs/dataevents/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,"_quantiles.png")
png(plotname, width = 20, height = 10,units = "cm", res = 300)
figdata
dev.off()


#####################################################################
#
#####################################################################
##################################################################################################
# Event apreciation
##################################################################################################
rm(list = ls())
setwd("~/data/Untitled/Trabajo/R_practice/exp_compound")
library(bnlearn)
library(magrittr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(parallel)
library(magrittr)
library(visualizeR)
library('RColorBrewer')
library(gridExtra)

source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")

########################################################################################################
# 
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/propagationFunctions.R")

degrees <- "10d"
mask <- FALSE # Only land or land & sea
limit <- FALSE # include poles?
b <- 3 # How many levels TP & T2M:
k <- 0 # Which permutation of variables later
compound <- "HD" # Which type of compound Hot Dry, Cold Dry, Hot Wet, Cold Wet, all together, all.distinct
global_mean_temp <- FALSE
global_mean_temp_cont <- FALSE# Extra node for global_mean_temp?
combi <- "TPT2M" # "only" Compound nodes or also "TPT2M" nodes
# algo<- "tabu" # Which algorithm
# score<- "bic" # Which score in algorithm
# global_mean_temp_art <- FALSE
# cut.art <- ""
scaleGT <- FALSE
# geo <- ""

if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
if(!is.null(compound)){comp <- paste0("_comp_",compound)} else {comp <-""}
if(isTRUE(mask)){superf <- "_land"}else {superf <-"_landsea"}
if(isTRUE(global_mean_temp)){n1 <- "_GMT"} else if(isTRUE(global_mean_temp_cont)){n1 <- "_GMTcont"} else {n1 <- ""}
# if(isTRUE(global_mean_temp_art)){n2 <- paste0("_GMTart_",cut.art,geo)} else {n2 <- ""}

# if(isTRUE(global_mean_temp_cont) & !isTRUE(global_mean_temp_art) & score == "bic"){score <- "bic-cg"} 
# if(isTRUE(global_mean_temp_cont) & !isTRUE(global_mean_temp_art) & score == "aic"){score <- "aic-cg"} 
# if (isTRUE(global_mean_temp_art)){n <- n2} else {n <- n1}

if(!k==0){permutations <-get(load(paste0("data/permutations/permutations_",degrees,lim,superf,".rda")))}
data.k <- get(load(file =paste0("data/manipulate_aggr/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".rda")))
data <- data.k
if(scaleGT == TRUE) {data$GT <- scale(data$GT)}
if(!k==0){data <-data[,c(permutations[[k]],(length(permutations[[k]])+1):ncol(data))]}


t2m_ERA5 <- get(load(paste0("data/raw_aggr/t2m_ERA5_monthly_1940_2022_",degrees,".rda")))
tp_ERA5 <- get(load(paste0("data/raw_aggr/tp_ERA5_monthly_1940_2022_",degrees,".rda")))
sel.tp <- subsetGrid(tp_ERA5,latLim = c(-50,75))
sel.t2m <- subsetGrid( t2m_ERA5,latLim = c(-50,75))

dataRMS2 <- data
x <- 81
var <- "t2m"
condvar <- paste0("V",x,".",var)
condlvl <- b

# condlvl <- 2
# condvar <- paste0("C",x)
# Big2ind <-  which(dataRMS2$GT == levels(dataRMS2$GT)[4])
# Big2ind <-  which(GT$GT == levels(GT$GT)[4])
# Big2ind <- 1:length(GT$GT)

Big2ind <-  which(dataRMS2[,condvar] == levels(dataRMS2[,condvar])[condlvl])

# event1 <- 1:5
# event2 <- 128:130
# 
# dataRMS2[,]
c <- sapply(dataRMS2[Big2ind,], FUN = as.numeric)
d <- sapply(dataRMS2, FUN = as.numeric)
extremes.T2M <- apply(c[,grep("C",colnames(c))],MARGIN = 2,FUN = function(x)sum(x==2)/length(x))
we.extremes.T2M <- apply(d[,grep("C",colnames(c))],MARGIN = 2,FUN = function(x)sum(x==2)/length(x))
dif.extremes.T2M <- extremes.T2M - we.extremes.T2M

# extremes.T2M <- apply(c[,grep("t2m",colnames(c))],MARGIN = 2,FUN = function(x)sum(x==1)/length(x))
# extremes.T2M <- apply(c[,grep("tp",colnames(c))],MARGIN = 2,FUN = function(x)sum(x==1)/length(x))



# 
# extremes <- colMeans(c)
# extremes.T2M <-extremes[grep("C",names(extremes))]


col.blue <- rev(brewer.pal(9,"Blues"))
col.blue
col.red <- brewer.pal(9,"Reds")

col.div <- c(col.blue, col.red)
col.div[4:5] <- "white"
col.b <- c(col.blue,col.red)
length(col.b)

col.r <- colorRampPalette(brewer.pal(9, "Reds"))(15)
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(15)
col.rb <- c(rev(col.b),'white','white',col.r)



length(extremes.T2M)
dim(sel.t2m$Data)
13*36

min(dif.extremes.T2M)
sel.t2m$Data


plotname <- paste0("figs/dataevents/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/dataevents/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)
climBig2 <- quantity2clim(we.extremes.T2M-((1/b)*(1/b)), what = "Big 2",t2m_ERA5)
figdata <- spatialPlot(climBig2, 
                       main = paste0("frequency ",compound," TRUE - independent frequency (",round(((1/b)*(1/b)),2),")"),
                       lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
                       col.regions = col.rb,set.min = -0.4, set.max = 0.4,
                       at = seq(-0.2,0.2,0.4/15))
figdata
dev.off()

assign(paste0("event_ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1),figdata)
save(list = paste0("event_ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1),file = paste0("results/dataevents/event_ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".rda"))

plotname <- paste0("figs/dataevents/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,"_",condvar,"_",condlvl,".pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/dataevents/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,"_",condvar,"_",condlvl,".png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)
climBig2 <- quantity2clim(extremes.T2M-((1/b)*(1/b)), what = "Big 2",t2m_ERA5)
figdata <- spatialPlot(climBig2, 
            # main = paste0("relative frequency ",compound," TRUE ", condvar," = ",condlvl, "- independent frequency (0.0625)"),
             main = paste0("frequency ",compound," TRUE ", condvar," = ",condlvl, " - independent frequency (",round(((1/b)*(1/b)),2),")"),
            lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.rb,set.min = -0.2, set.max = 0.2,
            at = seq(-0.2,0.2,0.4/15))
figdata
dev.off()

assign(paste0("event_ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,"_",condvar,"_",condlvl),figdata)
save(list = paste0("event_ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,"_",condvar,"_",condlvl),file = paste0("results/dataevents/event_ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,"_",condvar,"_",condlvl,".rda"))


figdata <- spatialPlot(climBig2, 
                       # main = paste0("relative frequency ",compound," TRUE ", condvar," = ",condlvl),
                       lonCenter = 0, backdrop.theme = "coastline", rev.colors = TRUE, 
                       col.regions = rev(col.b), set.min = -0.4,
                       at = seq(-0.4,0,0.4/15))
figdata

##################################################################
#
##################################################################
fig.HD.V81 <-get(load("/data/Untitled/Trabajo/R_practice/exp_compound/results/dataevents/event_ERA5_monthly_1940_2022_10d_disc_3_comp_HDTPT2M_landsea_V81.t2m_3.rda"))
fig.HD <- get(load("/data/Untitled/Trabajo/R_practice/exp_compound/results/dataevents/event_ERA5_monthly_1940_2022_10d_disc_3_comp_HDTPT2M_landsea.rda"))
fig.HW.V81 <- get(load("/data/Untitled/Trabajo/R_practice/exp_compound/results/dataevents/event_ERA5_monthly_1940_2022_10d_disc_3_comp_HWTPT2M_landsea_V81.t2m_3.rda"))
fig.HW <- get(load("/data/Untitled/Trabajo/R_practice/exp_compound/results/dataevents/event_ERA5_monthly_1940_2022_10d_disc_3_comp_HWTPT2M_landsea.rda"))

plotname <- paste0("figs/dataevents/event_ERA5_monthly_1940_2022_10d_disc_3_comp_HDandHW_landsea_V81.t2m_3",".pdf")
pdf(plotname,width = 12, height = 7)
grid.arrange(fig.HW,fig.HW.V81,fig.HD,fig.HD.V81)
dev.off()


# 
# Big2ind <- which(eval(parse(text = paste0("dataRMS2$",node," >=",value))))
# Big2ind <- which(eval(parse(text = paste0("dataRMS2$",node," <=",-value))))
# Big2ind <-  which(dataRMS2$V81 > 1.5 & dataRMS2$V227 > 1)
# Big2ind <-  which(dataRMS2$V81 > 0.5)
# events <- list()
# events[[1]] <- c(Big2ind[1])
# m <- 1
# if (!length(Big2ind)==1){
#   for (k in 2:length(Big2ind)){
#     if ((Big2ind[k]-1) == Big2ind[k-1]){events[[m]] <- append(events[[m]], Big2ind[k])
#     } else {events[[m+1]] <- c(Big2ind[k])
#     m <- m+1}
#   }
# } 
# plots <- list()
# for(i in 1:length(events)){
#   extremes <- colMeans(dataRMS2[events[[i]],])
#   
#   climBig2 <- quantity2clim(extremes, what = "Big 2",tas_interim_10d_akima_cubic)
#   plots[[i]] <- spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
#                             col.regions = col.div,set.min = -2.5, set.max = 2.5,
#                             at = seq(-3,3,1),
#                             main = list(paste0(length(events[[i]]))))
# }
# 
# do.call(grid.arrange,c(plots))

global_mean_temp <- TRUE
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


########################################################################################
# two evidences
########################################################################################
dataRMS2 <- data
x <- c(81,81)
var <- c("t2m","tp")
condvar <- paste0("V",x,".",var)
condlvl <- c(4,4)

# condlvl <- 2
# condvar <- paste0("C",x)
# Big2ind <-  which(dataRMS2$GT == levels(dataRMS2$GT)[4])
# Big2ind <-  which(GT$GT == levels(GT$GT)[4])
# Big2ind <- 1:length(GT$GT)

Big2ind <-  which(dataRMS2[,condvar[1]] == levels(dataRMS2[,condvar[1]])[condlvl[1]]&dataRMS2[,condvar[2]] == levels(dataRMS2[,condvar[2]])[condlvl[2]])

# event1 <- 1:5
# event2 <- 128:130
# 
# dataRMS2[,]
c <- sapply(dataRMS2[Big2ind,], FUN = as.numeric)
d <- sapply(dataRMS2, FUN = as.numeric)
extremes.T2M <- apply(c[,grep("C",colnames(c))],MARGIN = 2,FUN = function(x)sum(x==2)/length(x))
we.extremes.T2M <- apply(d[,grep("C",colnames(c))],MARGIN = 2,FUN = function(x)sum(x==2)/length(x))
dif.extremes.T2M <- extremes.T2M - we.extremes.T2M - 0.0625

# extremes.T2M <- apply(c[,grep("t2m",colnames(c))],MARGIN = 2,FUN = function(x)sum(x==1)/length(x))
# extremes.T2M <- apply(c[,grep("tp",colnames(c))],MARGIN = 2,FUN = function(x)sum(x==1)/length(x))



# 
# extremes <- colMeans(c)
# extremes.T2M <-extremes[grep("C",names(extremes))]


col.blue <- rev(brewer.pal(9,"Blues"))
col.blue
col.red <- brewer.pal(9,"Reds")

col.div <- c(col.blue, col.red)
col.div[4:5] <- "white"
col.b <- c(col.blue,col.red)
length(col.b)

col.r <- colorRampPalette(brewer.pal(9, "Reds"))(15)
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(15)
col.rb <- c(rev(col.b),'white','white',col.r)

plotname <- paste0("figs/dataevents/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,"_",condvar[1],"_",condlvl[1],"_",condvar[2],"_",condlvl[2],".pdf")
pdf(plotname, width = 10, height = 7)
plotname <- paste0("figs/dataevents/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,"_",condvar[1],"_",condlvl[1],"_",condvar[2],"_",condlvl[2],".png")
png(plotname, width = 20, height = 15,units = "cm", res = 150)
climBig2 <- quantity2clim(dif.extremes.T2M, what = "Big 2",t2m_ERA5)
figdata <- spatialPlot(climBig2, 
                       main = paste0("relative frequency ",compound," TRUE ", condvar[1]," = ",condlvl[1]," & ", condvar[2]," = ",condlvl[2]),
                       lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
                       col.regions = col.rb,set.min = -0.4, set.max = 0.4
                       ,at = seq(-0.4,0.4,0.8/15)
                      )
figdata
dev.off()

###############################################################################
#
###############################################################################
spatialPlot(climatology(sel.t2m), backdrop.theme = "coastline",color.theme = "RdBu",rev.colors = TRUE )
spatialPlot(climatology(sel.tp), backdrop.theme = "coastline",color.theme = "BrBG",rev.colors = FALSE)
subsetGrid(sel.t2m, lonLim = c(-18),latLim =c(46) )
madagaskar.t2m <- subsetGrid(sel.t2m, lonLim = c(-18),latLim =c(46) )
quantile(madagaskar.t2m$Data)[4]
temporalPlot(subsetGrid(sel.t2m, lonLim = c(-18),latLim =c(46) ))
temporalPlot(subsetGrid(sel.tp, lonLim = c(-18),latLim =c(46) ))
add_line(x = c(0,82*12),y = c(289,289))
lines
brewer.pal.info
