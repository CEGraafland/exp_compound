########################################################
rm(list = ls())
setwd("~/data/Untitled/Trabajo/R_practice/exp_compound")
library("abind")
library(loadeR)
library(visualizeR)
library(transformeR)
########################################################
# load data and mask
########################################################
landmask <- loadGridData("data/mask_W5E5v2.0.nc", var = "mask")
load("data/raw_aggr/t2m_ERA5_monthly_1940_2022_10d.rda")
load("data/raw_aggr/tp_ERA5_monthly_1940_2022_10d.rda")
########################################################
# interpolate mask to 10 data
########################################################
spatialPlot(climatology(landmask))
landmask_up10d <- upscaleGrid(landmask, times = 20)
landmask_comp <- interpGrid(redim(landmask_up10d,drop = TRUE),new.coordinates = getCoordinates(t2m_ERA5_monthly_1940_2022_10d))
landmask_comp$Data[landmask_comp$Data==0]<-NA
spatialPlot(landmask_comp,backdrop.theme = "coastline")

t2m_ERA5_monthly_1940_2022_10d<- redim(t2m_ERA5_monthly_1940_2022_10d,drop = TRUE)

timeslices <- dim(t2m_ERA5_monthly_1940_2022_10d$Data)[1]
masks<- rep(list(landmask_comp$Data),timeslices)
landmask_time <- abind(masks,along = 0)
t2m_land_ERA5_monthly_1940_2022_10d <- gridArithmetics(t2m_ERA5_monthly_1940_2022_10d,landmask_time)
spatialPlot(climatology(t2m_land_ERA5_monthly_1940_2022_10d))

tp_ERA5_monthly_1940_2022_10d<- redim(tp_ERA5_monthly_1940_2022_10d,drop = TRUE)
tp_land_ERA5_monthly_1940_2022_10d <- gridArithmetics(tp_ERA5_monthly_1940_2022_10d,landmask_time)
spatialPlot(climatology(tp_land_ERA5_monthly_1940_2022_10d))

landmask_10d <- landmask_comp
save(landmask_10d,file ="data/mask_aggr/landmask_10d.rda")
save(t2m_land_ERA5_monthly_1940_2022_10d,file ="data/mask_aggr/t2m_land_ERA5_monthly_1940_2022_10d.rda")
save(tp_land_ERA5_monthly_1940_2022_10d,file ="data/mask_aggr/tp_land_ERA5_monthly_1940_2022_10d.rda")


########################################################
# load data and mask
########################################################
landmask <- loadGridData("data/mask_W5E5v2.0.nc", var = "mask")
load("data/raw_aggr/t2m_ERA5_monthly_1940_2022_5d.rda")
load("data/raw_aggr/tp_ERA5_monthly_1940_2022_5d.rda")
########################################################
# interpolate mask to 5 data
########################################################
spatialPlot(climatology(landmask))
landmask_up5d <- upscaleGrid(landmask, times = 10)
landmask_comp <- interpGrid(redim(landmask_up5d,drop = TRUE),new.coordinates = getCoordinates(t2m_ERA5_monthly_1940_2022_5d))
landmask_comp$Data[landmask_comp$Data==0]<-NA
spatialPlot(landmask_comp,backdrop.theme = "coastline")

t2m_ERA5_monthly_1940_2022_5d<- redim(t2m_ERA5_monthly_1940_2022_5d,drop = TRUE)
timeslices <- dim(t2m_ERA5_monthly_1940_2022_5d$Data)[1]
masks<- rep(list(landmask_comp$Data),timeslices)
landmask_time <- abind(masks,along = 0)
t2m_land_ERA5_monthly_1940_2022_5d <- gridArithmetics(t2m_ERA5_monthly_1940_2022_5d,landmask_time)
spatialPlot(climatology(t2m_land_ERA5_monthly_1940_2022_5d))

tp_ERA5_monthly_1940_2022_5d<- redim(tp_ERA5_monthly_1940_2022_5d,drop = TRUE)
tp_land_ERA5_monthly_1940_2022_5d <- gridArithmetics(tp_ERA5_monthly_1940_2022_5d,landmask_time)
spatialPlot(climatology(tp_land_ERA5_monthly_1940_2022_5d))

landmask_5d <- landmask_comp
save(landmask_5d,file ="data/mask_aggr/landmask_5d.rda")
save(t2m_land_ERA5_monthly_1940_2022_5d,file ="data/mask_aggr/t2m_land_ERA5_monthly_1940_2022_5d.rda")
save(tp_land_ERA5_monthly_1940_2022_5d,file ="data/mask_aggr/tp_land_ERA5_monthly_1940_2022_5d.rda")
