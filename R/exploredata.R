#############################################################
# explore data
#############################################################
#############################################################
rm(list = ls())
setwd("~/data/Untitled/Trabajo/R_practice/exp_compound")
#############################################################
# source("../../creds")
library(devtools)
library(loadeR)
library(visualizeR)
library(bnlearn)
#options(java.parameters = "-Xmx8000m")
#rJava::.jinit(parameters = "-Xmx1000m")
library(transformeR)
# library(visualizeR)
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
library(RColorBrewer)

####################################################################
# tp
####################################################################
var <- "tp"
if (var =="t2m"){aggr.fun <- "mean"} else if (var == "tp"){aggr.fun <- "sum"}
files <- list.files(paste0("/data/Untitled/Trabajo/R_practice/exp_compound/data/",var,"_pyth"), full.names = T)
names.ext <- list.files(paste0("/data/Untitled/Trabajo/R_practice/exp_compound/data/",var,"_pyth"))
names <- gsub(".rda", "",names.ext)
gridList <- lapply(files, function(x){get(load(x))})
names(gridList) <- names
tp_ERA5_monthly_1940_2022_5d <- bindGrid(gridList,dimension = "time")
# if ... > 90 <-90 then lonLim <- ... 
tp_ERA5_monthly_1940_2022_5d <- subsetGrid(tp_ERA5_monthly_1940_2022_5d,latLim = c(-87.625,87.375))
tp_ERA5_monthly_1940_2022_5d_detr <- detrendGrid(tp_ERA5_monthly_1940_2022_5d)
tp_ERA5_monthly_1940_2022_10d <- upscaleGrid(tp_ERA5_monthly_1940_2022_5d, times = 2,aggr.fun = list(FUN = aggr.fun,na.rm = TRUE))

save(tp_ERA5_monthly_1940_2022_10d,file = "data/raw_aggr/tp_ERA5_monthly_1940_2022_10d.rda")
save(tp_ERA5_monthly_1940_2022_5d,file = "data/raw_aggr/tp_ERA5_monthly_1940_2022_5d.rda")

###################################################################
# t2m
###################################################################
var <- "t2m"
files <- list.files(paste0("/data/Untitled/Trabajo/R_practice/exp_compound/data/",var,"_pyth"), full.names = T)
names.ext <- list.files(paste0("/data/Untitled/Trabajo/R_practice/exp_compound/data/",var,"_pyth"))
names <- gsub(".rda", "",names.ext)
gridList <- lapply(files, function(x){get(load(x))})
names(gridList) <- names
t2m_ERA5_monthly_1940_2022_5d <- bindGrid(gridList,dimension = "time")
# if ... > 90 <-90 then lonLim <- ... 
t2m_ERA5_monthly_1940_2022_5d <- subsetGrid(t2m_ERA5_monthly_1940_2022_5d,latLim = c(-87.625,87.375))

if (var =="t2m"){aggr.fun <- "mean"} else if (var == "tp"){aggr.fun <- "sum"}
t2m_ERA5_monthly_1940_2022_10d <- upscaleGrid(t2m_ERA5_monthly_1940_2022_5d, times = 2,aggr.fun = list(FUN = aggr.fun,na.rm = TRUE))
t2m_ERA5_monthly_1940_2022_5d_detr <- detrendGrid(t2m_ERA5_monthly_1940_2022_5d)
save(t2m_ERA5_monthly_1940_2022_10d,file = "data/raw_aggr/t2m_ERA5_monthly_1940_2022_10d.rda")
save(t2m_ERA5_monthly_1940_2022_5d,file = "data/raw_aggr/t2m_ERA5_monthly_1940_2022_5d.rda")

##############################################################
# combine t2m tp
##############################################################
t2m_tp_monthly_1940_2022_5d<- bindGrid(t2m_ERA5_monthly_1940_2022_5d,tp_ERA5_monthly_1940_2022_5d,dimension = "member")
t2m_tp_monthly_1940_2022_5d_detr <- bindGrid(t2m_ERA5_monthly_1940_2022_5d_detr,tp_ERA5_monthly_1940_2022_5d_detr,dimension = "member")

###############################################################################
#
###############################################################################
display.brewer.all()

tp_ERA5_monthly_1940_2022_10d_scaled  <- scaleGrid(tp_ERA5_monthly_1940_2022_10d,type = "center", time.frame = "monthly")
t2m_ERA5_monthly_1940_2022_10d_scaled <- scaleGrid(t2m_ERA5_monthly_1940_2022_10d, type = "center",time.frame = "monthly")



# plotname <- "figs/"
# pdf()
spatialPlot(climatology(tp_ERA5_monthly_1940_2022_10d_scaled),color.theme = "BrBG" ,rev.colors = FALSE, lonCenter = 180,backdrop.theme = "coastline")
spatialPlot(climatology(t2m_ERA5_monthly_1940_2022_10d_scaled),rev.colors = TRUE, lonCenter = 180,backdrop.theme = "coastline")

temporalPlot(t2m_ERA5_monthly_1940_2022_5d)
temporalPlot(t2m_ERA5_monthly_1940_2022_5d_detr)
temporalPlot(tp_ERA5_monthly_1940_2022_5d)
temporalPlot(tp_ERA5_monthly_1940_2022_5d_detr)
spatialPlot(climatology(t2m_ERA5_monthly_1940_2022_5d_detr),rev.colors = TRUE, lonCenter = 180,backdrop.theme = "coastline")
spatialPlot(climatology(tp_ERA5_monthly_1940_2022_5d_detr),rev.colors = TRUE, lonCenter = 180,backdrop.theme = "coastline")

##########################################
# explore cutting worldmap
##########################################
load("data/raw_aggr/tp_ERA5_monthly_1940_2022_10d.rda")
load("data/raw_aggr/t2m_ERA5_monthly_1940_2022_10d.rda")
test <- subsetGrid(tp_ERA5_monthly_1940_2022_10d,latLim = c(-50,75))
spatialPlot(climatology(test),rev.colors = TRUE, lonCenter = 180,backdrop.theme = "coastline")
str(test)
test$xyCoords$y


###########################################
# distribution tests
###########################################
sel.tp <- tp_ERA5_monthly_1940_2022_10d
sel.t2m <- t2m_ERA5_monthly_1940_2022_10d
df.tp <-as.data.frame(TimeCoordsAnom_from_Grid_rms(sel.tp,rms = TRUE))
df.t2m <-as.data.frame(TimeCoordsAnom_from_Grid_rms(sel.t2m,rms = TRUE))

sel.node <- "V600"
ks.test(x = df.tp[,sel.node], y = "pbeta",shape1 = 2,shape2 = 3)
hist(df.tp[,sel.node])
d <-density(df.t2m[,sel.node])
plot(d)
shapiro.test(df.tp[,sel.node])
chisq.test(df.tp[,sel.node],rescale.p = TRUE)
qqnorm(y = df.tp[,sel.node])
############################################
# discretize
############################################
disc.df.tp <- discretize(df.tp, method = "quantile", breaks = 5)
disc.df.t2m <- discretize(df.t2m, method = "quantile", breaks = 5)
names(disc.df.t2m) <- paste0(names(disc.df.t2m),".t2m")
names(disc.df.tp) <- paste0(names(disc.df.tp),".tp")
disc.df.t2m.tp <- cbind(disc.df.t2m,disc.df.tp)

###########################################
#
###########################################
bn <- hc(disc.df.t2m.tp,max.iter = 10)
bn

gc()

plot(disc.df.t2m,disc.df.tp)
plot(disc.df.t2m)
hist(disc.df.t2m)
disc.df.t2m[1:3,1]

cor(as.numeric(disc.df.t2m[,1]),as.numeric(disc.df.tp[,1]))
cor(as.numeric(disc.df.t2m[,1]),as.numeric(disc.df.tp[,1]))
plot(as.numeric(disc.df.t2m[,1]),as.numeric(disc.df.tp[,1]))
cor(as.numeric(disc.df.t2m[,1]),as.numeric(disc.df.t2m[,2]))
n.t2m <-sapply(names(disc.df.t2m),FUN = function(x)as.numeric(disc.df.t2m[[x]]))
n.tp <-sapply(names(disc.df.tp),FUN = function(x)as.numeric(disc.df.tp[[x]]))

clim.grid.cor.t2m.tp <- quantity2clim(diag(cor.t2m.tp), "cor tp t2m",sel.tp)
spatialPlot((clim.grid.cor.t2m.tp),rev.colors = TRUE, lonCenter = 180,backdrop.theme = "coastline")

cor(as.vector(n.t2m),as.vector(n.tp), method = "spearman")
cor.t2m.tp <- cor(n.t2m,n.tp)
plot(diag(cor.t2m.tp))
rm(list = ls())
gc()
dev.off()

################################################################
#
################################################################

