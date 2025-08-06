#########################################################################################
# Template calcular distancias y plotear diferentes timesteps redes bayesianas dinamicas.
#########################################################################################
rm(list = ls())
setwd("~/data/Untitled/Trabajo/R_practice/exp_compound")
source("../R/Functions/BasicNetworkFunctions.R")
###########################################################################################
# Funcion para luego calcular distancias: 
# puedes guardarlo en un script apartado y luego cargarlo con 'source'
# Tambien puedes utilizar la funcion que has utilizado antes que no era haversine
###########################################################################################
haversine <- function(x1Lat,x1Lon,x2Lat,x2Lon){
  earthR <- 6371 #using mean radius
  mLat <- as.double(x1Lat)
  bLat <- as.double(x2Lat)
  mLong <- as.double(x1Lon)
  bLong <- as.double(x2Lon)
  changeLat <- (mLat - bLat)/180*pi
  changeLong <- (mLong - bLong)/180*pi
  a <- sin(changeLat/2) * sin(changeLat/2) + cos((mLat)/180*pi) * 
    cos(bLat/180*pi) * sin(changeLong/2) * sin(changeLong/2)
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  distKm <- earthR * c
  distMi <- as.double(distKm * 0.621371192)
  output <- c(x1Lat,x1Lon,x2Lat,x2Lon,distKm,distMi)
  return(output[5])}     

######################################################################
# descargar packages
######################################################################
library(bnlearn)
library(transformeR)
library(visualizeR)
library(magrittr)
library(igraph)
library(broom)
library(maps)
library(geosphere)
library(rgeos)
library(gridExtra)
library(ggplot2)
require(data.table)  

##################################################################
# Load Redes Bayesianas
##################################################################
limit <- TRUE 
if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
b <- 3
# load("/data/Untitled/Trabajo/R_practice/TFG_Celia/results/wg2step.rda")
# load("/data/Untitled/Trabajo/R_practice/TFG_Celia/wg.rda")
# load("~/data/Untitled/Trabajo/R_practice/TFG_Celia/results/RBD_hc_interim_discrete_wg3.Rda")
# load("~/data/Untitled/Trabajo/R_practice/TFG_Celia/results/RBD_hc_interim_discrete_wgf.Rda")
#load("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_10d/perm0/",algo,"_disc_",b,"_",score,"_t2m_tp_10d_0_2600_2700i.rda")
#load("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_10d/perm0/",algo,"_disc_",b,"_",score,"_t2m_tp_10d_0_2600_2700i.rda")
algo<-"tabu"
score <- "aic"
load("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_10d/perm0/",algo,"_disc_",b,"_",score,"_t2m_tp_10d_0_1600_1700i.rda")
load("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_comp_V105_10d/perm0/",algo,"_disc_",b,"_",score,"_t2m_tp_comp_V105_10d_0_3600_3700i.rda")
load("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_10d",lim,"/perm0/",algo,"_disc_",b,"_",score,"_t2m_tp_10d",lim,"_0_1200_1300i.rda")


load("data/raw_aggr/t2m_ERA5_monthly_1940_2022_10d.rda")
load("data/raw_aggr/tp_ERA5_monthly_1940_2022_10d.rda")

###########################################
#
###########################################




###########################################
# Standardize
###########################################

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

disc.df.tp <- discretize(df.tp, method = "quantile", breaks = b)
disc.df.t2m <- discretize(df.t2m, method = "quantile", breaks = b)
names(disc.df.t2m) <- paste0(names(disc.df.t2m),".t2m")
names(disc.df.tp) <- paste0(names(disc.df.tp),".tp")
disc.df.t2m.tp <- cbind(disc.df.t2m,disc.df.tp)
# df.t2m<- NULL
# df.tp <- NULL
# disc.df.t2m <- NULL
# disc.df.tp <- NULL
gc()
##########################################################################
# Preparacion histograma (y plots)
##########################################################################
# Extraer algunas caracterisitas
dag <- hc_disc_4_aic_t2m_tp_10d_0_2600_2700i
dag <- tabu_disc_4_bic_t2m_tp_10d_0_1600_1700i
dag <- tabu_disc_3_aic_t2m_tp_comp_V105_10d_0_3600_3700i
dag <- tabu_disc_3_bic_t2m_tp_10d_lim_0_1200_1300i

dag <- get(load("/data/Untitled/Trabajo/R_practice/exp_compound/data/tabuiterations/disc_3_aic_t2m_tp_10d_lim/perm0/tabu_disc_3_aic_t2m_tp_10d_lim_0_2700_2800i.rda"))

#data.network <- sel.t2m.1
#rm(data.network)
x.t2m<- attr(sel.t2m.1, "Xcoords", exact = FALSE)
x.tp <- attr(sel.tp.1,"Xcoords", exact = FALSE)

y.t2m <- attr(sel.t2m.1, "Ycoords", exact = FALSE)
y.tp <- attr(sel.tp.1, "Ycoords", exact = FALSE)

# x <- (x+360)%%360
p.t2m <- expand.grid(y.t2m, x.t2m)[2:1]
p.tp <- expand.grid(y.tp,x.tp)[2:1]
p <- rbind(p.t2m,p.tp)
colnames(p)<-c("x","y")
arcsbn <- arcs(dag)
epochs <- 2

# Convertir hacia objeto igraph
igraph <- igraph.from.graphNEL(as.graphNEL(dag))
nodenames <- nodes(dag) 
row.names(p)<- nodenames
# Make edgelist in igraph class by edgeindices adjacancy matrix
# for selecting edgeindices when applying Haversine
adjmat <- as.matrix(as_adjacency_matrix(igraph))
edgeindices <- which(adjmat == 1, arr.ind = TRUE)
fromV <- nodenames[edgeindices[,1]]
toV <- nodenames[edgeindices[,2]]
fromtoV <- cbind(fromV,toV)

# Make edgelist in igraph class method 2
indsame <- c()
for (i in 1:nrow(arcsbn)){
  int <- intersect(which(as_edgelist(igraph)[,1] == arcsbn[i,1]),
                   which(as_edgelist(igraph)[,2] == arcsbn[i,2]))
  indsame[i] <- int
}

# Make coordinates for every variable
p
longitude <- p[,"x"]
lattitude <- p[,"y"]

# estimate distance of all edges in igraph-edgelist.
distances <- c()

for (i in 1:nrow(edgeindices)){
  x1Lat <- lattitude[edgeindices[i,1]]
  x1Lon <- longitude[edgeindices[i,1]]
  x2Lat <- lattitude[edgeindices[i,2]]
  x2Lon <- longitude[edgeindices[i,2]]
  
  disti <- haversine(x1Lat,x1Lon,x2Lat,x2Lon)
  distances[i] <- disti}

# Make dataframe with departing variable, end variable and distance
arcdistances <- cbind(unname(fromtoV), round(distances, digits = 3)) 
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

abreviations <- paste0(sub('.*\\.', '', fromV),sub('.*\\.', '', toV))

# Utilizando este data.frame puedes hacer el histograma:
arcdistances.df <- data.frame(fromV,toV,distances,abreviations)
arcdistances.df[arcdistances.df$distances==0,]
######################################################################
# Continuacion para los plots
######################################################################
# identify which indices in igraph edgelist correspond to indices in bn
indsame2 <- c()
nrow(as_edgelist(igraph))
for (i in 1:nrow(as_edgelist(igraph))){
  int <- intersect(which(arcdistances[,1] == as_edgelist(igraph)[i,1]),
                   which(arcdistances[,2] == as_edgelist(igraph)[i,2]))
  indsame2[i] <- int
}

# permute distance vector to find belong distance vector for igraph object
newdistances <- distances[indsame2]
E(igraph)$distances <- newdistances

edgelist <- as_edgelist(igraph) 
gridpoints <- row.names(p)# AQUI 

dists <- edge.attributes(igraph)$distances
distsv <- round(c(0,max(dists)/4,max(dists)/2,max(dists)/4*3, max(dists)))
distsc <- as.character(distsv)

col.1 <- adjustcolor("light blue", alpha=0.8)
col.2 <- adjustcolor("navy blue", alpha=0.8)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.col <- edge.pal(100)

arcslist <- list()

for(i in 1:nrow(edgelist))  {
  node1 <- gridpoints[gridpoints == edgelist[i,1]]
  
  node2 <- gridpoints[gridpoints == edgelist[i,2]]
  # Hier er tussen proppen of het groter dan 90 is of kleiner lattitude
  if(!(p[node1,"x"] == p[node2,"x"] && p[node1,"y"] ==p[node2,"y"])){
  arc <- gcIntermediate(c(p[node1,"x"], p[node1,"y"]), 
                        c(p[node2,"x"], p[node2,"y"]),
                        n=100, addStartEnd=TRUE, breakAtDateLine = FALSE)
  } else {arc <- matrix(c(p[node1,"x"], p[node1,"y"]),nrow = 1, ncol =2,dimnames = list(NULL,c('lon','lat')))}
  
  edge.ind <- ceiling(edge_attr(igraph,"distances",E(igraph)[i])*100 / max(dists))
  if(grepl("tp",node1) && grepl("tp",node2)){
    edge.time <- 1
  } else if (grepl("t2m",node1) && grepl("t2m",node2)) {
    edge.time <- 2
  } else if (grepl("T2",node1) && grepl("T2",node2)) {
    edge.time <- 5
  } else if (grepl("t2m",node1) && grepl("tp",node2)){
    edge.time <- 3
  } else if (grepl("tp",node1) && grepl("t2m",node2)){
    edge.time <- 4
  } else if (grepl("T0",node1) && grepl("T2",node2)){
  edge.time <- 6
  } else if (grepl("T2",node1) && grepl("T0",node2)){
  edge.time <- 7
  } else if (grepl("T1",node1) && grepl("T2",node2)){
    edge.time <- 8
  } else if (grepl("T2",node1) && grepl("T1",node2)){
    edge.time <- 9
  }
  
  col <- as.character(edge.col[edge.ind])
  negs <- which(arc[,'lon']<0)
  
  if(!is.null(negs)){arc[negs] <- (arc[negs]+360)%%360}
  
  arc <- cbind(arc,edge.ind,edge.time)
  arcslist[[i]]<- arc
  
}

namesarcs <- character()
for(i in 1:length(arcslist)){namesarcs[i] <- paste0("arc",i)}
names(arcslist) <- namesarcs

mapworld2 <- map("world2", fill = TRUE,plot = FALSE)
world_pol_df <- tidy(mapworld2, IDs = "region")

mapworld2$orientation

do.call(rbind.data.frame, lapply(names(arcslist), function(x) {
  cbind.data.frame(route=x, arcslist[[x]], stringsAsFactors=FALSE)
})) -> all_paths



# table.all_paths <- table(all_paths$route)
# names.single <- names(which(table.all_paths ==1))
# str(names.single)
# str(all_paths$route)
# all_paths[all_paths["route"] == names.single,]
# all_paths[all_paths["route"][names.single]]
# all_paths[all_paths["edge.ind"] == 0,]$route
# names.single
times <- c("tptp","t2mt2m","t2mtp","tpt2m","T2T2","T0T2","T2T0","T1T2","T2T1")
times.indices<- 1:9

names(times.indices) <- times

splitplot <- function(x) ggplot() +
  geom_polygon(data = world_pol_df, mapping = aes(x = long, y = lat, group = group)) +
  geom_line(data=all_paths[all_paths$edge.time==times.indices[x],], 
            aes(x=lon, y=lat, group=route, col = edge.ind),
            arrow = arrow(angle = 30, length = unit(0.02, "inches"),
                          ends = "last", type = "closed")) +
  scale_color_gradient(low = col.1, high = col.2,limits = c(0,100), breaks = c(0,25,50,75,100), guide = guide_colorbar(title = "Distance"), labels = distsc) +
  geom_point(data = all_paths[all_paths$edge.time==times.indices[x]&all_paths$edge.ind==0,],aes(x=lon, y=lat, group=route,col = edge.ind))+
  ggtitle(paste0("BN: ",times[x]," ",length(unique(all_paths[all_paths$edge.time == x,]$route)),"/",nrow(edgelist))) +
  theme(plot.title = element_text(hjust = 0.5))

plotname <- paste0("figs/",algo,"_",score,"_",nrow(edgelist),"_t2m_tp_windows",lim)

pdf(paste0(plotname,".pdf"), height = 7, width = 10)
png(paste0(plotname,".png"), width = 25, height = 20,units = "cm", res = 300)

splitplots <- lapply(times.indices[1:(epochs*epochs)],splitplot)
do.call("grid.arrange",c(list(allplot),splitplots))

dev.off()
####################################################################
#
###################################################################
# all.plot <- function(x) 
  
allplot<-   ggplot() +
  geom_polygon(data = world_pol_df, mapping = aes(x = long, y = lat, group = group)) +
  geom_line(data=all_paths, 
            aes(x=lon, y=lat, group=route, col = edge.ind),
            arrow = arrow(angle = 30, length = unit(0.02, "inches"),
                          ends = "last", type = "closed")) +
  scale_color_gradient(low = col.1, high = col.2,limits = c(0,100), breaks = c(0,25,50,75,100), guide = guide_colorbar(title = "Distance"), labels = distsc) +
  geom_point(data = all_paths[all_paths$edge.ind==0,],aes(x=lon, y=lat, group=route,col = edge.ind))+
  ggtitle(paste0("BN:  ",nrow(edgelist))) +
  theme(plot.title = element_text(hjust = 0.5))

combiplot <-  ggplot() +
  geom_polygon(data = world_pol_df, mapping = aes(x = long, y = lat, group = group)) +
  geom_line(data=all_paths[all_paths$edge.time==times.indices[3]|all_paths$edge.time==times.indices[4],], 
            aes(x=lon, y=lat, group=route, col = edge.ind),
            arrow = arrow(angle = 30, length = unit(0.02, "inches"),
                          ends = "last", type = "closed")) +
  scale_color_gradient(low = col.1, high = col.2,limits = c(0,100), breaks = c(0,25,50,75,100), guide = guide_colorbar(title = "Distance"), labels = distsc) +
  geom_point(data = all_paths[(all_paths$edge.time==times.indices[3]|all_paths$edge.time==times.indices[4])&all_paths$edge.ind==0,],aes(x=lon, y=lat, group=route,col = edge.ind))+
  ggtitle(paste0("BN: ",times[3],"&",times[4],"",length(unique(all_paths[all_paths$edge.time == 3|all_paths$edge.time == 4,]$route)),"/",nrow(edgelist))) +
  theme(plot.title = element_text(hjust = 0.5))

plotname <- paste0("figs/",algo,"_",score,"_",nrow(edgelist),"_t2m_tp_all_split_combi_windows",lim)

pdf(paste0(plotname,".pdf"), height = 7, width = 10)
png(paste0(plotname,".png"), width = 25, height = 20,units = "cm", res = 300)

splitplots <- lapply(times.indices[1:2],splitplot)
do.call("grid.arrange",c(list(allplot),splitplots,list(combiplot)))

dev.off()
