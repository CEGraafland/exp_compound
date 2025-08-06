rm(list = ls())
setwd("~/data/Untitled/Trabajo/R_practice/exp_compound")
library(bnlearn)
library(magrittr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(sparsebn)
library(sparsebnUtils)
library(RColorBrewer)
library(gaussDiff)
library("copula")
library("bnmonitor")
library(visualizeR)
########################################################################################################
# Model Evaluation Compound
########################################################################################################
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
#load("../Data/tas_ncep_10d.rda")
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

permused <- 0
algo <- "tabu"
score <- "aic"

loadIterationsComp <- function(permused = 0, algo = c("hc","tabu"),score = c("aic","bic"), it = NULL, limit = FALSE,compound = NULL) {
  if(isTRUE(limit)){lim <- "_lim"} else {lim <-""}
  if(!is.null(compound)){comp <- paste0("_comp_",compound)} else {comp <-""}
  
  hc_list <- list.files(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_10d",lim,comp,"/perm",permused), full.names = T)
  hc_names <- list.files(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_10d",lim,comp,"/perm",permused))
  hc_names <- gsub(".rda", "", hc_names)
  hc_networks <- lapply(hc_list, function(x){get(load(x))})
  names(hc_networks) <- hc_names
sizes <- sapply(hc_networks,narcs)
  hc_networks <- hc_networks[order(sizes)]
   return(hc_networks)
 }

###############################################################
# evalutation tabu versus hc aic versus bic
###############################################################
tabu_aic_0 <- loadIterationsComp(permused = 0,algo = "tabu", score = "aic",limit = limit)
tabu_bic_0 <- loadIterationsComp(permused = 0,algo = "tabu", score = "bic",limit = limit)
hc_aic_0 <- loadIterationsComp(permused = 0,algo = "hc", score = "aic")
hc_bic_0 <- loadIterationsComp(permused = 0,algo = "hc", score = "bic")


tabu_aic_0_fits<- lapply(tabu_aic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp,debug = TRUE))
tabu_bic_0_fits<- lapply(tabu_bic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))
hc_aic_0_fits<- lapply(hc_aic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))
hc_bic_0_fits<- lapply(hc_bic_0, function(x) bn.fit(x = x, data = disc.df.t2m.tp))

lo_tabu_aic_0 <- sapply(X = tabu_aic_0_fits,FUN =logLik,data = disc.df.t2m.tp)
lo_tabu_bic_0 <- sapply(X = tabu_bic_0_fits,FUN =logLik,data = disc.df.t2m.tp)
lo_hc_aic_0 <-sapply(X = hc_aic_0_fits,FUN =logLik,data = disc.df.t2m.tp)
lo_hc_bic_0 <- sapply(X = hc_bic_0_fits,FUN =logLik,data = disc.df.t2m.tp)

sta0 <-sapply(FUN = narcs, tabu_aic_0)
stb0 <-sapply(FUN = narcs, tabu_bic_0)
sha0 <-sapply(FUN = narcs, hc_aic_0)
shb0 <-sapply(FUN = narcs, hc_bic_0)

plotname <- "figs/loglik_tabu_hc_bic_aic_compound"
pdf(paste0(plotname,".pdf"))
png(paste0(plotname,".png"))

plot(sta0,lo_tabu_aic_0,col = rainbow(4)[1],pch = 16, xlab = "|E|",ylab = c("Loglik(x|model)"))
points(stb0,lo_tabu_bic_0,col = rainbow(4)[2],pch = 16)
points(sha0,lo_hc_aic_0,col = rainbow(4)[3],pch = 16)
points(shb0,lo_hc_bic_0,col = rainbow(4)[4],pch = 16)

legend("bottomright",legend = c("tabu aic","tabu bic","hc aic","hc bic"), fill = rainbow(4),cex = 0.7)
dev.off()

#####################################################
# Evaluation artificial compound versus learned compound
#####################################################
tabu_bic_0_HD_learn <- loadIterationsComp(permused = 0,algo = "tabu", score = "bic",limit = limit,compound = "HD")
tabu_bic_0 <- loadIterationsComp(permused = 0,algo = "tabu", score = "bic",limit = limit)

art.compound.dag <- function(g){
  ngrids <- length(nodes(g))/2
  comp.nodes <- paste0("C",1:ngrids)
  g.comp <- g
  for(i in 1:length(comp.nodes)){
    g.comp <- add.node(g.comp,comp.nodes[i])
    g.comp <- set.arc(g.comp, nodes(g)[i],comp.nodes[i])
    g.comp <- set.arc(g.comp, nodes(g)[i+ngrids],comp.nodes[i])
  }
  return(g.comp)
}

tabu_bic_0_art <- lapply(tabu_bic_0,art.compound.dag)

tabu_bic_0_HD_learn_fits<- lapply(tabu_bic_0_HD_learn, function(x) bn.fit(x = x, data = disc.df.t2m.tp))
tabu_bic_0_HD_art_fits <- lapply(tabu_bic_0_art, function(x) bn.fit(x = x, data = disc.df.t2m.tp))

lo_tabu_bic_0_HD_learn <- sapply(X = tabu_bic_0_HD_learn_fits,FUN =logLik,data = disc.df.t2m.tp)
lo_tabu_bic_0_art_HD <- sapply(X = tabu_bic_0_HD_art_fits,FUN =logLik,data = disc.df.t2m.tp)

stb0_HD_learn <-sapply(FUN = narcs, tabu_bic_0_HD_learn)
stb0_HD_art <-sapply(FUN = narcs, tabu_bic_0_art)

plotname <- "figs/loglik_tabu_bic_art_learn_compound"
pdf(paste0(plotname,".pdf"))
png(paste0(plotname,".png"))

plot(stb0_HD_learn,lo_tabu_bic_0_HD_learn,col = rainbow(2)[1],pch = 16, xlab = "|E|",ylab = c("Loglik(x|model)"))
points(stb0_HD_art,lo_tabu_bic_0_art_HD,col = rainbow(2)[2],pch = 16)

legend("bottomright",legend = c("tabu bic comp learned","tabu bic comp art"), fill = rainbow(2),cex = 0.7)
dev.off()


##################################################
#
##################################################
textdf <- cbind(attr(sel.t2m.1,"VertexCoords"),1:ngrids)
#textdf$x[textdf$x<0] <-textdf$x[textdf$x<0] + 360
textdf <- as.matrix(textdf)
textlay <- apply(textdf,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "black"))

clim.0vector <- quantity2clim(numeric(length=ngrids),"zero vector", ref.grid = sel.t2m)
spatialPlot(clim.0vector, backdrop.theme = "coastline",sp.layout = textlay, lonCenter =0, colorkey = FALSE)
#####################################################
#
#####################################################
art <- TRUE
narcs(tabu_bic_0_art$tabu_disc_3_bic_t2m_tp_10d_lim_0_1100_1200i)
tabu
if(isTRUE(art)){g <-tabu_bic_0_art$tabu_disc_3_bic_t2m_tp_10d_lim_0_1100_1200i}

g <-tabu_bic_0_art$tabu_disc_3_bic_t2m_tp_10d_lim_0_1100_1200i
g <-tabu_bic_0_HD_learn$tabu_disc_3_bic_t2m_tp_10d_lim_comp_HD_0_1100_1200i
gc()
# virt.nodes <- gsub("V","X",nodes(g))
# g.virt <- g
# for(i in 1:length(virt.nodes)){
#   g.virt <- add.node(g.virt,virt.nodes[i])
#   g.virt <- set.arc(g.virt, virt.nodes[i],nodes(g)[i])
# }
# inset <- sapply(virt.nodes, function(x) dsep(g.virt, x ,nodes(g)[300],nodes(g)[303]))
# inset.det <- sapply(nodes(g), function(x) dsep(g.virt, x ,nodes(g)[300],nodes(g)[303]))
# 
# depset <- which(inset == FALSE)
# depset.det <- which(inset.det == FALSE)
# 
# length(depset)
# length(depset.det)
# depset %in%depset.det
# 
# 
# depset[39]
# 
# all.equal(depset, depset.det)
# 
# g$nodes[300]
# g$nodes[303]

# ngrids <- length(nodes(g))/2
# comp.nodes <- paste0("C",1:ngrids)
# g.comp <- g
# for(i in 1:length(comp.nodes)){
#   g.comp <- add.node(g.comp,comp.nodes[i])
#   g.comp <- set.arc(g.comp, nodes(g)[i],comp.nodes[i])
#   g.comp <- set.arc(g.comp, nodes(g)[i+ngrids],comp.nodes[i])
# }

#virt.nodes <- gsub("V","X",nodes(g))
virt.nodes <- gsub("V","X",nodes(g)[grep("V",nodes(g))])
gc()
g.virt <- g
for(i in 1:length(virt.nodes)){
  g.virt <- add.node(g.virt,virt.nodes[i])
  g.virt <- set.arc(g.virt, virt.nodes[i],nodes(g)[i])
}

z<- 346
inset <- sapply(virt.nodes, function(x) dsep(g.virt, x ,comp.nodes[z]))
inset.det <- sapply(nodes(g), function(x) dsep(g.virt, x ,comp.nodes[z]))
# inset.det2 <-sapply(nodes(g), function(x) dsep(g.comp, x ,comp.nodes[z]))

gc()

depset <- which(inset == FALSE)
depset.det <- which(inset.det == FALSE)

t2ms <- which(depset<=ngrids)
t2ms <- as.numeric(!inset[1:ngrids])
t2ms[z]<- 2

tps <- which(depset>ngrids)
tps <- as.numeric(!inset[(ngrids+1):(2*ngrids)])
tps[z]<- 2

tps.det <- which(depset.det>ngrids)
tps.det <- as.numeric(!inset.det[(ngrids+1):(2*ngrids)])
tps.det[z]<- 2

t2ms.det <- which(depset.det<=ngrids)
t2ms.det <- as.numeric(!inset.det[1:ngrids])
t2ms.det[z]<- 2

clim.prop.depset.t2ms <- quantity2clim(t2ms,"prob dependence set", ref.grid = sel.t2m)
spatialPlot(clim.prop.depset.t2ms, backdrop.theme = "coastline")

clim.prop.depset.det.t2ms <- quantity2clim(t2ms.det,"det dependence set", ref.grid = sel.t2m)
spatialPlot(clim.prop.depset.det.t2ms, backdrop.theme = "coastline")

clim.prop.depset.tps <- quantity2clim(tps,"prob dependence set", ref.grid = sel.tp)
spatialPlot(clim.prop.depset.tps, backdrop.theme = "coastline")

clim.prop.depset.det.tps <- quantity2clim(tps.det,"det dependence set", ref.grid = sel.tp)
spatialPlot(clim.prop.depset.det.tps, backdrop.theme = "coastline")

clim.props <- bindGrid(clim.prop.depset.t2ms,clim.prop.depset.tps,clim.prop.depset.det.t2ms,clim.prop.depset.det.tps, dimension = "member")
namesmembers <- c(paste0("C",z,"~ P(Vi.t2m|parents(Vi.t2m))"), paste0("C",z,"~ P(Vi.tp|parents(Vi.tp))"),paste0("C",z,"~ Vi.t2m"), paste0("C",z,"~ Vi.tp"))

spatialPlot(clim.props, backdrop.theme = "coastline",names.attr = namesmembers, as.table = TRUE, color.theme = "YlGnBu")
brewer.pal.info
#####################################################
# pruebas analisis de sensitividad.
#####################################################
fitted <- tabu_aic_0_fits$tabu_disc_3_aic_t2m_tp_10d_lim_0_900_1000i
param_node <- "V81.t2m"
lvl.param_node <- 3
lvl.param_padres <- 3
padres <- bnlearn::parents(fitted,param_node)
valor_padres<- sapply(disc.df.t2m.tp[padres],FUN = function(x)levels(x))
if(length(padres)==1){vp <- valor_padres[lvl.param_padres]} else {vp <- vp <- valor_padres[lvl.param_padres,]}


query_node <-"V99.tp"
lvl.query_node <- 1
query_ev <- "V99.t2m"
lvl.query_ev <- 3

send_dic <- sensitivity(fitted, 
                        interest_node = query_node, interest_node_value = levels(disc.df.t2m.tp[[query_node]])[[lvl.query_node]],
                        evidence_nodes =query_ev, evidence_states =levels(disc.df.t2m.tp[[query_ev]])[[lvl.query_ev]], 
                        node = param_node, 
                        value_node = levels(disc.df.t2m.tp[[param_node]])[[lvl.param_node]], 
                        value_parents = vp,
                        new_value = "all",
                        covariation = "proportional")
send_dic$plot <- send_dic$plot + ggtitle(paste0("Sensitivity of P(",query_node,"=",lvl.query_node,"|",query_ev,"=",lvl.query_ev,") under new_values of P(",param_node,"=",lvl.param_node,"|",padres," = ",lvl.param_padres,")"))

plot(send_dic)

write.net("data/hciterations/tabu_disc_4_aic_t2m_tp_10d_0_1200_1300i.net",tabu_aic_0_fits$tabu_disc_4_aic_t2m_tp_10d_0_1200_1300i)

cpquery
#################################################################
#
#################################################################
###########################################
# Standardize
###########################################
sel.tp <- tp_ERA5_monthly_1940_2022_10d
sel.t2m <- t2m_ERA5_monthly_1940_2022_10d
tp_ERA5_monthly_1940_2022_10d <- NULL
t2m_ERA5_monthly_1940_2022_10d <- NULL
df.tp <-as.data.frame(TimeCoordsAnom_from_Grid_rms(sel.tp,rms = TRUE))
df.t2m <-as.data.frame(TimeCoordsAnom_from_Grid_rms(sel.t2m,rms = TRUE))
sel.tp<- NULL
sel.t2m <- NULL
############################################
# discretize in b
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


############################################################################
# Function to load hciterations of models in a list 
############################################################################

permused <- 0
algo <- "tabu"
score <- "aic"
node.comp <-"V105"

loadIterationsComp <- function(permused = 0, algo = c("hc","tabu"),score = c("aic","bic"), it = NULL,node.comp = NULL) {
  if(is.null(node.comp)){comp.ind <- ""} else {comp.ind <- paste0("comp_",node.comp,"_")}
  hc_list <- list.files(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_",comp.ind,"10d/perm",permused), full.names = T)
  hc_names <- list.files(paste0("data/",algo,"iterations/disc_",b,"_",score,"_t2m_tp_",comp.ind,"10d/perm",permused))
  hc_names <- gsub(".rda", "", hc_names)
  hc_networks <- lapply(hc_list, function(x){get(load(x))})
  names(hc_networks) <- hc_names
  sizes <- sapply(hc_networks,narcs)
  hc_networks <- hc_networks[order(sizes)]
  return(hc_networks)
}


tabu_aic_0_V105 <- loadIterationsComp(permused = 0,algo = "tabu", score = "aic",node.comp = "V105")


v_i <- "V105"
data <- disc.df.t2m.tp
data.num <- as.data.frame(sapply(disc.df.t2m.tp, function(x) as.numeric(x), simplify = FALSE))

data[[paste0(v_i,".comp")]]<- numeric(nrow(data))
assign(paste0("data.num.",v_i,".2state"),disc.df.t2m.tp)

data[[paste0(v_i,".comp")]][data.num[[paste0(v_i,".t2m")]]==1 & data.num[[paste0(v_i,".tp")]]==3] <- 13
data[[paste0(v_i,".comp")]][data.num[[paste0(v_i,".t2m")]]==1 & data.num[[paste0(v_i,".tp")]]==1] <- 11
data[[paste0(v_i,".comp")]][data.num[[paste0(v_i,".t2m")]]==3 & data.num[[paste0(v_i,".tp")]]==1] <- 31
data[[paste0(v_i,".comp")]][data.num[[paste0(v_i,".t2m")]]==3 & data.num[[paste0(v_i,".tp")]]==3] <- 33
data[[paste0(v_i,".comp")]] <- as.factor(data[[paste0(v_i,".comp")]])


tabu_aic_0_V105_fits<- lapply(tabu_aic_0_V105, function(x) bn.fit(x = x, data = data))


lo_tabu_aic_0_V105 <- sapply(X = tabu_aic_0_V105_fits,FUN =logLik,data = data)

sta0 <-sapply(FUN = narcs, tabu_aic_0_V105)

# plotname <- "figs/loglik_tabu_hc_bic_aic_compound"
# pdf(paste0(plotname,".pdf"))
# png(paste0(plotname,".png"))

plot(sta0,lo_tabu_aic_0_V105,col = rainbow(4)[1],pch = 16, xlab = "|E|",ylab = c("Loglik(x|model)"))

# legend("bottomright",legend = c("tabu aic","tabu bic","hc aic","hc bic"), fill = rainbow(4),cex = 0.7)
# dev.off()

fitted <- tabu_aic_0_V105_fits$tabu_disc_3_aic_t2m_tp_comp_V105_10d_0_2000_2100i
param_node <- "V82.t2m"
padres <- bnlearn::parents(fitted,param_node)
valor_padres<- sapply(data[padres],FUN = function(x)levels(x))
lvl.param_node <- 3
lvl.param_padres <- 3
if(length(padres)==1){vp <- valor_padres[lvl.param_padres]} else {vp <- vp <- valor_padres[lvl.param_padres,]}



query_node <-"V105.comp"
lvl.query_node <- 2
query_ev <- NULL
# query_ev <- "V99.t2m"
# lvl.query_ev <- 3
if(is.null(query_ev)){es <- NULL} else {es <- levels(data[[query_ev]])[[lvl.query_ev]]}


send_dic <- sensitivity(fitted, 
                        interest_node = query_node, interest_node_value = levels(data[[query_node]])[[lvl.query_node]],
                        evidence_nodes =query_ev, evidence_states =es, 
                        node = param_node, 
                        value_node = levels(data[[param_node]])[[lvl.param_node]], 
                        value_parents = vp,
                        new_value = "all",
                        covariation = "proportional")
send_dic$plot <- send_dic$plot + ggtitle(paste0("Sensitivity of P(",query_node,"=",lvl.query_node,"|",query_ev,"=",lvl.query_ev,") under new_values of P(",param_node,"=",lvl.param_node,"|",padres," = ",lvl.param_padres,")"))

plot(send_dic)
