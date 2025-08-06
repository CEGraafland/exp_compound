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


degrees <- "10d"
mask <- FALSE # Only land or land & sea
limit <- TRUE # include poles?
b <- 3 # How many levels TP & T2M:
k <- 0 # Which permutation of variables later
IT <- paste0(k,"_4000_4100")
IT <- paste0(k,"_2700_2800")
compound <- "HD" # Which type of compound Hot Dry, Cold Dry, Hot Wet, Cold Wet, all together, all.distinct
global_mean_temp <-TRUE 
global_mean_temp_cont <- FALSE# Extra node for global_mean_temp?
combi <- "TPT2M" # only Compound nodes or also "TPT2M" nodes
algo<- "tabu" # Which algorithm
score<- "aic" # Which score in algorithm
global_mean_temp_art <- TRUE
cut.art <- 30
geo <- "auto"

if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
if(!is.null(compound)){comp <- paste0("_comp_",compound)} else {comp <-""}
if(isTRUE(mask)){superf <- "_land"}else {superf <-"_landsea"}
if(isTRUE(global_mean_temp)){n1 <- "_GMT"} else if(isTRUE(global_mean_temp_cont)){n1 <- "_GMTcont"} else {n1 <- ""}
if(isTRUE(global_mean_temp_art)){n2 <- paste0("_GMTart_",cut.art,geo)} else {n2 <- ""}

if(isTRUE(global_mean_temp_cont) & !isTRUE(global_mean_temp_art) & score == "bic"){score <- "bic-cg"} 
if(isTRUE(global_mean_temp_cont) & !isTRUE(global_mean_temp_art) & score == "aic"){score <- "aic-cg"} 
if (isTRUE(global_mean_temp_art)){n <- n2} else {n <- n1}
####################################
# mask
####################################
degrees <- "10d" #10d
mask <- FALSE
if(isTRUE(mask)){
  t2m_ERA5 <- get(load(paste0("data/mask_aggr/t2m_land_ERA5_monthly_1940_2022_",degrees,".rda")))
  tp_ERA5 <- get(load(paste0("data/mask_aggr/tp_land_ERA5_monthly_1940_2022_",degrees,".rda")))
  landmask <- get(load(paste0("data/mask_aggr/landmask_",degrees,".rda")))
  if(isTRUE(limit)){landmask <- subsetGrid(landmask,latLim = c(-50,75))
  t2m_ERA5 <- subsetGrid(t2m_ERA5,latLim = c(-50,75))
  tp_ERA5 <-  subsetGrid(tp_ERA5,latLim = c(-50,75)) }
} else {
  t2m_ERA5 <- get(load(paste0("data/raw_aggr/t2m_ERA5_monthly_1940_2022_",degrees,".rda")))
  tp_ERA5 <- get(load(paste0("data/raw_aggr/tp_ERA5_monthly_1940_2022_",degrees,".rda")))
  if(isTRUE(limit)){
  t2m_ERA5 <- subsetGrid(t2m_ERA5,latLim = c(-50,75))
  tp_ERA5 <-  subsetGrid(tp_ERA5,latLim = c(-50,75)) }
  }

data <- get(load(file =paste0("data/manipulate_aggr/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".rda")))
##################################################################
# Load Redes Bayesianas
##################################################################
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

dag <- tabu_disc_3_bic_lim_GMT_art12_comp_HDonly_land_0_2100_2200i
dag <- loadIterationsComp(permused = k,algo = algo, score = score,it = IT)[[1]]
# tabu_bic_0_fits<- lapply(tabu_bic_0, function(x) bn.fit(x = x, data = data))
nodes(dag)
which(arcs(dag)[,1]=="GT")
arcs(dag)[which(arcs(dag)[,1]=="GT"),]
##########################################################################
# Preparacion histograma (y plots)
##########################################################################
#rm(data.network)

x.template<- t2m_ERA5$xyCoords$x
y.template<- t2m_ERA5$xyCoords$y

# x <- (x+360)%%360
keep <- as.vector(landmask$Data)
keepsel <- as.logical(keep)
keepsel[is.na(keep)] <- FALSE
p.template <- expand.grid(y.template, x.template)[2:1]
p.template.sel <-p.template[keepsel,]
p.extra <- expand.grid(0,-180)[2:1]
p <- rbind(p.template.sel,p.extra)
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
  if(grepl("C",node1) && grepl("C",node2)){
    edge.time <- 1
  } else if (grepl("GT",node1) && grepl("C",node2)) {
    edge.time <- 2
  } else if (grepl("T2",node1) && grepl("T2",node2)) {
    edge.time <- 5
  } else if (grepl("C",node1) && grepl("GT",node2)){
    edge.time <- 2
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

mapworld2 <- map("world2", fill = TRUE,plot = FALSE, col = "grey")
world_pol_df <- tidy(mapworld2, IDs = "region")

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
times <- c("CC","GT-C","C-GT","tpt2m","T2T2","T0T2","T2T0","T1T2","T2T1")
times.indices<- 1:9
names(times.indices) <- times
x<-1
splitplot <- function(x) ggplot() +
  geom_polygon(data = world_pol_df, mapping = aes(x = long, y = lat, group = group)) +
  geom_line(data=all_paths[all_paths$edge.time==times.indices[x],], 
            aes(x=lon, y=lat, group=route, col = edge.ind),
            arrow = arrow(angle = 30, length = unit(0.02, "inches"),
                          ends = "last", type = "closed")) +
  scale_color_gradient(low = col.1, high = col.2,limits = c(0,100), breaks = c(0,25,50,75,100), guide = guide_colorbar(title = "Distance"), labels = distsc) +
  geom_point(data = all_paths[all_paths$edge.time==times.indices[x]&all_paths$edge.ind==0,],aes(x=lon, y=lat, group=route,col = edge.ind))+
  ggtitle(paste0("BN: ",times[x]," ",nrow(edgelist))) +
  theme(plot.title = element_text(hjust = 0.5))

plotname <- paste0("figs/",algo,"_",score,"_",nrow(edgelist),"_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2)

pdf(paste0(plotname,".pdf"), height = 20, width = 10)
png(paste0(plotname,".png"), width = 25, height = 20,units = "cm", res = 300)

splitplots <- lapply(times.indices[1:3],splitplot)
do.call("grid.arrange",splitplots)
splitplots[[2]]

dev.off()



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


degrees <- "10d"
mask <- FALSE # Only land or land & sea
limit <- TRUE # include poles?
b <- 3 # How many levels TP & T2M:
k <- 0 # Which permutation of variables later
IT <- paste0(k,"_4000_4100")
IT <- paste0(k,"_2700_2800")
compound <- "HD" # Which type of compound Hot Dry, Cold Dry, Hot Wet, Cold Wet, all together, all.distinct
global_mean_temp <-TRUE 
global_mean_temp_cont <- FALSE# Extra node for global_mean_temp?
combi <- "TPT2M" # only Compound nodes or also "TPT2M" nodes
algo<- "tabu" # Which algorithm
score<- "aic" # Which score in algorithm
global_mean_temp_art <- TRUE
cut.art <- 30
geo <- "auto"

if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
if(!is.null(compound)){comp <- paste0("_comp_",compound)} else {comp <-""}
if(isTRUE(mask)){superf <- "_land"}else {superf <-"_landsea"}
if(isTRUE(global_mean_temp)){n1 <- "_GMT"} else if(isTRUE(global_mean_temp_cont)){n1 <- "_GMTcont"} else {n1 <- ""}
if(isTRUE(global_mean_temp_art)){n2 <- paste0("_GMTart_",cut.art,geo)} else {n2 <- ""}

if(isTRUE(global_mean_temp_cont) & !isTRUE(global_mean_temp_art) & score == "bic"){score <- "bic-cg"} 
if(isTRUE(global_mean_temp_cont) & !isTRUE(global_mean_temp_art) & score == "aic"){score <- "aic-cg"} 
if (isTRUE(global_mean_temp_art)){n <- n2} else {n <- n1}
####################################
# mask
####################################
degrees <- "10d" #10d
mask <- FALSE
if(isTRUE(mask)){
  t2m_ERA5 <- get(load(paste0("data/mask_aggr/t2m_land_ERA5_monthly_1940_2022_",degrees,".rda")))
  tp_ERA5 <- get(load(paste0("data/mask_aggr/tp_land_ERA5_monthly_1940_2022_",degrees,".rda")))
  landmask <- get(load(paste0("data/mask_aggr/landmask_",degrees,".rda")))
  if(isTRUE(limit)){landmask <- subsetGrid(landmask,latLim = c(-50,75))
  t2m_ERA5 <- subsetGrid(t2m_ERA5,latLim = c(-50,75))
  tp_ERA5 <-  subsetGrid(tp_ERA5,latLim = c(-50,75)) }
} else {
  t2m_ERA5 <- get(load(paste0("data/raw_aggr/t2m_ERA5_monthly_1940_2022_",degrees,".rda")))
  tp_ERA5 <- get(load(paste0("data/raw_aggr/tp_ERA5_monthly_1940_2022_",degrees,".rda")))
  if(isTRUE(limit)){
    t2m_ERA5 <- subsetGrid(t2m_ERA5,latLim = c(-50,75))
    tp_ERA5 <-  subsetGrid(tp_ERA5,latLim = c(-50,75)) }
}

data <- get(load(file =paste0("data/manipulate_aggr/ERA5_monthly_1940_2022_",degrees,"_disc_",b,lim,comp,combi,superf,n1,".rda")))
##################################################################
# Load Redes Bayesianas
##################################################################
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


dag <- loadIterationsComp(permused = k,algo = algo, score = score,it = IT)[[1]]
# tabu_bic_0_fits<- lapply(tabu_bic_0, function(x) bn.fit(x = x, data = data))
nodes(dag)
which(arcs(dag)[,1]=="GT")
arcs(dag)[which(arcs(dag)[,1]=="GT"),]
##########################################################################
# Preparacion histograma (y plots)
##########################################################################
#rm(data.network)

x.template<- t2m_ERA5$xyCoords$x
y.template<- t2m_ERA5$xyCoords$y

# x <- (x+360)%%360
# keep <- as.vector(landmask$Data)
# keepsel <- as.logical(keep)
# keepsel[is.na(keep)] <- FALSE
p.template <- expand.grid(y.template, x.template)[2:1]
p.template.sel <- rbind(p.template,p.template,p.template)
# p.template.sel <-p.template[keepsel,]
p.extra <- expand.grid(0,-180)[2:1]
p <- rbind(p.template.sel,p.extra)
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
i <- 1
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
  if(grepl("C",node1) && grepl("C",node2)){
    edge.time <- 1
  } else if (grepl("GT",node1) && grepl("C",node2)) {
   edge.time <- 2  
  } else if (grepl("C",node1) && grepl("GT",node2)){
     edge.time <- 2
  } else if (grepl("GT",node1) && grepl("t2m",node2)) {
    edge.time <- 3
  } else if (grepl("GT",node1) && grepl("tp",node2)){
    edge.time <- 3
  } else if (grepl("tp",node1) && grepl("GT",node2)){
    edge.time <- 3
  } else if (grepl("t2m",node1) && grepl("GT",node2)){
    edge.time <- 3
  } else if (grepl("t2m",node1) && grepl("tp",node2)){
    edge.time <- 4
  } else if (grepl("t2m",node1) && grepl("t2m",node2)){
    edge.time <- 4
  } else if (grepl("tp",node1) && grepl("tp",node2)){
    edge.time <- 4  
  } else if (grepl("tp",node1) && grepl("t2m",node2)){
    edge.time <- 4  
  } else if (grepl("tp",node1) && grepl("C",node2)) {
    edge.time <- 5  
  } else if (grepl("t2m",node1) && grepl("C",node2)){
    edge.time <- 5  
  } else if (grepl("C",node1) && grepl("tp",node2)) {
      edge.time <- 5  
  } else if (grepl("C",node1) && grepl("t2m",node2)){
      edge.time <- 5
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

mapworld2 <- map("world2", fill = TRUE,plot = FALSE, col = "grey")
world_pol_df <- tidy(mapworld2, IDs = "region")

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
times <- c("CC","GT-C","GT-tp-t2m","tpt2m","C-tp-t2m")
times.indices<- 1:5
names(times.indices) <- times
x<-4
splitplot<- function(x) ggplot() +
  geom_polygon(data = world_pol_df, mapping = aes(x = long, y = lat, group = group)) +
  geom_line(data=all_paths[all_paths$edge.time==times.indices[x],], 
            aes(x=lon, y=lat, group=route, col = edge.ind),
            arrow = arrow(angle = 30, length = unit(0.02, "inches"),
                          ends = "last", type = "closed")) +
  scale_color_gradient(low = col.1, high = col.2,limits = c(0,100), breaks = c(0,25,50,75,100), guide = guide_colorbar(title = "Distance"), labels = distsc) +
   geom_point(data = all_paths[all_paths$edge.time==times.indices[x]&all_paths$edge.ind==0,],aes(x=lon, y=lat, group=route,col = edge.ind))+
  ggtitle(paste0("BN: ",times[x]," ",length(unique(all_paths[all_paths$edge.time == times.indices[x],]$route)),"/",length(unique(all_paths[all_paths$edge.time == times.indices[4],]$route)))) +
  theme(plot.title = element_text(hjust = 0.5))

plotname <- paste0("figs/",algo,"_",score,"_",length(unique(all_paths[all_paths$edge.time == times.indices[4],]$route)),"_",degrees,"_disc_",b,"_",score,lim,comp,combi,superf,n1,n2)

pdf(paste0(plotname,".pdf"), height = 20, width = 10)
png(paste0(plotname,".png"), width = 25, height = 20,units = "cm", res = 300)

splitplots <- lapply(times.indices[1:5],splitplot)
do.call("grid.arrange",splitplots[3:4])
splitplots[[2]]

dev.off()
