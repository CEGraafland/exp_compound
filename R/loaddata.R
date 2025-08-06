###################################################################################
# load ERA5 tp y t2m
###################################################################################
rm(list = ls())
setwd("~/data/Untitled/Trabajo/R_practice/exp_compound")
# source("../../creds")
library(devtools)
library(loadeR)
library(visualizeR)
library(bnlearn)
options(java.parameters = "-Xmx8000m")
#rJava::.jinit(parameters = "-Xmx1000m")
library(transformeR)
# library(visualizeR)
####################################################################################
# ERA5 t2m y tp 
# 1981 t/m 2010
####################################################################################
 year <- 1981 
# dir.create(paste0("/data/Untitled/Trabajo/R_practice/exp_compound/data/",var))

for(year in 1995:2010){
  var <- "t2m"
  if (var =="t2m"){aggr.fun <- "mean"} else if (var == "tp"){aggr.fun <- "sum"}
  files <- list.files(paste0("/oceano/gmeteo/WORK/DATA/C3S-CDS/ERA5/day/",var,"/",year), full.names = T)
  names.ext <- list.files(paste0("/oceano/gmeteo/WORK/DATA/C3S-CDS/ERA5/day/",var,"/",year))
  names <- gsub(".nc", "",names.ext)
  gridList.daily.year <- lapply(files, function(x){loadGridData(x,var=var)})# Aqui hay problemas
  names(gridList.daily.year) <- names
  grid.daily.year <- bindGrid(gridList.daily.year,dimension = "time")
  gridList.daily.year<- NULL
  gc()
  grid.monthly.year <- aggregateGrid(grid.daily.year,aggr.m = list(FUN = aggr.fun))
  grid.daily.year <- NULL
  gc()
  grid.monthly.year.5d <- upscaleGrid(grid.monthly.year, times = 20,aggr.fun = list(FUN = aggr.fun,na.rm = TRUE))
  grid.monthly.year <- NULL
  gc()
  assign(paste0(var,"_ERA5_monthly_",year),grid.monthly.year.5d)
  save(list = paste0(var,"_ERA5_monthly_",year),file = paste0("/data/Untitled/Trabajo/R_practice/exp_compound/data/",var,"/",var,"_ERA5_monthly_",year,".rda"))
  rm(list = ls())
  gc()
}

loadGridData("/oceano/gmeteo/WORK/DATA/C3S-CDS/ERA5/day/t2m/1981/t2m_ERA5_day_198101.nc", var = "t2m")
############################################################################
# directamente aggregado con python
############################################################################

###################################################################
# t2m
###################################################################

evenyears <- (1940:2021)[c(TRUE,FALSE)]
odyears<-(1940:2021)[c(FALSE,TRUE)]
pairyears <- cbind(evenyears,odyears)
var <- "t2m"
if (var =="t2m"){aggr.fun <- "mean"} else if (var == "tp"){aggr.fun <- "sum"}
i <- 3
for(i in 3:nrow(pairyears)){
  years <- pairyears[i,]
  # years <- 2022
  grid.monthly <- loadGridData(paste0("/oceano/gmeteo/WORK/DATA/C3S-CDS/ERA5/mon/era5_mon_",var,"_1940_2022.nc"), var ="t2m",years = years)
  grid.monthly.5d <- upscaleGrid(grid.monthly, times = 20,aggr.fun = list(FUN = aggr.fun,na.rm = TRUE))
  grid.monthly<- NULL
  gc()
  for (year in years) {
    assign(paste0(var,"_ERA5_monthly_",year), subsetGrid(grid.monthly.5d,years = year))
  }
  lapply(years,function(year) save(list = paste0(var,"_ERA5_monthly_",year),file = paste0("/data/Untitled/Trabajo/R_practice/exp_compound/data/",var,"_pyth/",var,"_ERA5_monthly_",year,".rda")))
  for (year in years) {
    assign(paste0(var,"_ERA5_monthly_",year), NULL)
  }
  gc()
}

###################################################################
# tp
###################################################################

evenyears <- (1940:2021)[c(TRUE,FALSE)]
odyears<-(1940:2021)[c(FALSE,TRUE)]
pairyears <- cbind(evenyears,odyears)
var <- "tp"
if (var =="t2m"){aggr.fun <- "mean"} else if (var == "tp"){aggr.fun <- "sum"}

for(i in 2:nrow(pairyears)){
  years <- pairyears[i,]
  years <- 2022
  grid.monthly <- loadGridData(paste0("/oceano/gmeteo/WORK/DATA/C3S-CDS/ERA5/mon/era5_mon_",var,"_1940_2022.nc"), var =var,years = years)
  grid.monthly.5d <- upscaleGrid(grid.monthly, times = 20,aggr.fun = list(FUN = aggr.fun,na.rm = TRUE))
  grid.monthly<- NULL
  gc()
  for (year in years) {
    assign(paste0(var,"_ERA5_monthly_",year), subsetGrid(grid.monthly.5d,years = year))
  }
  lapply(years,function(year) save(list = paste0(var,"_ERA5_monthly_",year),file = paste0("/data/Untitled/Trabajo/R_practice/exp_compound/data/",var,"_pyth/",var,"_ERA5_monthly_",year,".rda")))
  for (year in years) {
    assign(paste0(var,"_ERA5_monthly_",year), NULL)
  }

}

