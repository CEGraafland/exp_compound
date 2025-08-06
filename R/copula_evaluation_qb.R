#############################################################
# Copula evaluation q = 3
#############################################################
rm(list = ls())
setwd("~/data/Untitled/Trabajo/R_practice/exp_compound")
#############################################################
library("copula")
library("VineCopula")
#library("bnmonitor")
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
###############################################################
#
###############################################################
load("data/raw_aggr/tp_ERA5_monthly_1940_2022_10d.rda")
#save(tp_ERA5_monthly_1940_2022_5d,file = "data/raw_aggr/tp_ERA5_monthly_1940_2022_5d.rda")
load("data/raw_aggr/t2m_ERA5_monthly_1940_2022_10d.rda")
#save(t2m_ERA5_monthly_1940_2022_5d,file = "data/raw_aggr/t2m_ERA5_monthly_1940_2022_5d.rda")
###########################################
# Standardize
###########################################
sel.tp <- tp_ERA5_monthly_1940_2022_10d
sel.t2m <- t2m_ERA5_monthly_1940_2022_10d
tp_ERA5_monthly_1940_2022_10d <- NULL
t2m_ERA5_monthly_1940_2022_10d <- NULL
sel.tp.1<- TimeCoordsAnom_from_Grid_rms(sel.tp,rms = TRUE)
sel.t2m.1<-TimeCoordsAnom_from_Grid_rms(sel.t2m,rms = TRUE)
df.tp <-as.data.frame(sel.tp.1)
df.t2m <-as.data.frame(sel.t2m.1)
# sel.tp<- NULL
# sel.t2m <- NULL

########################################
# df
########################################
b <- 3
quantile(df.t2m$V81,probs = seq(0,1,1/3))[b]
ind.V81.t2m.qb <- which(df.t2m$V81>=quantile(df.t2m$V81)[b])
df.t2m.V81.t2m.qb <- df.t2m[ind.V81.t2m.qb,]
df.tp.V81.t2m.qb <- df.tp[ind.V81.t2m.qb,]

# random quantile:
set.seed(1)
sample.ind <- sample(1:nrow(df.t2m),nrow(df.t2m)/b)
random.qb.sample.ind <- sample.ind
df.t2m.random.qb <- df.t2m[random.qb.sample.ind,]
df.tp.random.qb <- df.tp[random.qb.sample.ind,]
########################################
# psuedo observaciones 
########################################
p.df.tp<- pobs(df.tp)
p.df.t2m <- pobs(df.t2m)

p.df.tp.V81.t2m.qb<- pobs(df.tp.V81.t2m.qb)
p.df.t2m.V81.t2m.qb <- pobs(df.t2m.V81.t2m.qb)

p.df.t2m.random.qb <- pobs(df.t2m.random.qb)
p.df.tp.random.qb <- pobs(df.tp.random.qb)

#bi2<- BiCopSelect(p.df.tp[,2], p.df.t2m[,2], familyset = NA)
# plot(BiCops[[5]])
# BiCops[[646]]$p.value.indeptest
#dev.off()

x_test <- 82
plot(p.df.tp[,x_test],p.df.t2m[,x_test])
plot(p.df.tp.V81.t2m.qb[,x_test],p.df.t2m.V81.t2m.qb[,x_test])
plot(df.tp[,x_test],df.t2m[,x_test])
plot(df.tp.V81.t2m.qb[,x_test],df.t2m.V81.t2m.qb[,x_test],ylim = c(min(df.t2m[,x_test]),max(df.t2m[,x_test])))


#########################################
# Estimate bivariate copulas
#########################################
# Without inicial ind. test
# BiCops <- lapply(1:ncol(p.df.tp),FUN = function(i)BiCopSelect(p.df.tp[,i], p.df.t2m[,i], familyset = NA))
# With inicial ind.test
# BiCops.ind <- lapply(1:ncol(p.df.tp),FUN = function(i)BiCopSelect(p.df.tp[,i], p.df.t2m[,i], familyset = NA,indeptest = TRUE))
# assign(paste0("BiCops.ind.q",b),BiCops.ind) 
# save(list = paste0("BiCops.ind.q",b),file = paste0("results/copulas/BiCops.ind.q",b,".rda"))
BiCops.ind <- get(load(file = paste0("results/copulas/BiCops.ind.q",b,".rda")))

# BiCops.V81.t2m.ind <- lapply(1:ncol(p.df.tp.V81.t2m.qb),FUN = function(i)BiCopSelect(p.df.tp.V81.t2m.qb[,i], p.df.t2m.V81.t2m.qb[,i], familyset = NA,indeptest = TRUE))
# assign(paste("BiCops.V81.t2m.q",b,"ind"),BiCops.V81.t2m.ind)
# save(list = paste("BiCops.V81.t2m.q",b,"ind"),file = paste0("results/copulas/BiCops.V81.t2m.q",b,".ind.rda"))
BiCops.V81.t2m.ind <- get(load(file = paste0("results/copulas/BiCops.V81.t2m.q",b,".ind.rda")))

# BiCops.random.ind <- lapply(1:ncol(p.df.tp.random.qb),FUN = function(i)BiCopSelect(p.df.tp.random.qb[,i], p.df.t2m.random.qb[,i], familyset = NA,indeptest = TRUE))
# assign(paste("BiCops.random.q",b,".ind"),BiCops.random.ind)
# save(list = paste("BiCops.V81.t2m.random.q",b,".ind"),file = paste0("results/copulas/BiCops.random.q",b,".ind.rda"))
BiCops.random.ind <- get(load(file = paste0("results/copulas/BiCops.random.q",b,".ind.rda")))
#BiCops.randomq4.ind <- BiCops.V81.t2m.randomq4.ind
#save(BiCops.randomq4.ind,file = "results/copulas/BiCops.randomq4.ind.rda")


# Only ind.tests
IndTests<- lapply(1:ncol(p.df.tp),FUN = function(i)BiCopIndTest(p.df.tp[,i], p.df.t2m[,i]))
IndTests.V81.t2m<- lapply(1:ncol(p.df.tp.V81.t2m.qb),FUN = function(i)BiCopIndTest(p.df.tp.V81.t2m.qb[,i], p.df.t2m.V81.t2m.qb[,i]))
IndTests.random<- lapply(1:ncol(p.df.tp.random.qb),FUN = function(i)BiCopIndTest(p.df.tp.random.qb[,i], p.df.t2m.random.qb[,i]))
##################################
# extraer caracteristicas copulas
##################################
BiCops.V81.t2m.ind[[621]]$par
BiCopsFam <- sapply(BiCops.ind,function(x) x$family)
BiCopsFam.V81.t2m<- sapply(BiCops.V81.t2m.ind,function(x) x$family)
BiCopsFam.random<- sapply(BiCops.random.ind,function(x) x$family)


gaus.min <- sapply(BiCops.ind,function(x) x$family ==1 && x$par<0)
BiCopsFam[gaus.min]<- 301
t.min <- sapply(BiCops.ind,function(x) x$family ==2 && x$par<0)
BiCopsFam[t.min]<- 302
frank.min <- sapply(BiCops.ind,function(x) x$family ==5 && x$par<0)
BiCopsFam[frank.min]<- 305

gaus.min.V81.t2m <- sapply(BiCops.V81.t2m.ind,function(x) x$family ==1 && x$par<0)
BiCopsFam.V81.t2m[gaus.min.V81.t2m]<- 301
t.min.V81.t2m <- sapply(BiCops.V81.t2m.ind,function(x) x$family ==2 && x$par<0)
BiCopsFam.V81.t2m[t.min.V81.t2m]<- 302
frank.min.V81.t2m <- sapply(BiCops.V81.t2m.ind,function(x) x$family ==5 && x$par<0)
BiCopsFam.V81.t2m[frank.min.V81.t2m]<- 305

gaus.min.random <- sapply(BiCops.random.ind,function(x) x$family ==1 && x$par<0)
BiCopsFam.random[gaus.min.V81.t2m]<- 301
t.min.random <- sapply(BiCops.random.ind,function(x) x$family ==2 && x$par<0)
BiCopsFam.random[t.min.V81.t2m]<- 302
frank.min.random <- sapply(BiCops.random.ind,function(x) x$family ==5 && x$par<0)
BiCopsFam.random[frank.min.V81.t2m]<- 305

BiCops.tail.low <- sapply(BiCops.ind,function(x) BiCopPar2TailDep(obj = x)$lower)
BiCops.tail.up <- sapply(BiCops.ind,function(x) BiCopPar2TailDep(obj = x)$upper)
BiCops.emptau <- sapply(BiCops.ind,function(x) x$emptau)
BiCops.V81.t2m.emptau <- sapply(BiCops.V81.t2m.ind,function(x) x$emptau)
BiCops.theortau <- sapply(BiCops.ind,function(x) BiCopPar2Tau(family = x$family, par = x$par, par2 = x$par2))

IndTestsP <- sapply(IndTests,function(x) x$p.value)
IndTestsP005 <- sapply(IndTests,function(x) as.numeric(x$p.value <= 0.05))
IndTestsP.V81.t2m <- sapply(IndTests.V81.t2m,function(x) x$p.value)
IndTestsP005.V81.t2m<- sapply(IndTests.V81.t2m,function(x) as.numeric(x$p.value <= 0.05))

# Small-> reject 0 of independence -> dependent  -- TRUE = 1 = dependent
# high -> acept 0 of independence -> indepdent -- FALSE = 0 = independent

clim.IndTestsP <- quantity2clim(IndTestsP005, what = "p.value Ind test",sel.tp)
clim.BiCopsFam <- quantity2clim(BiCopsFam, what = "Copula Family",sel.tp)
clim.tail.low <- quantity2clim(BiCops.tail.low, what = "lower tail dependency",sel.tp)
clim.tail.up <- quantity2clim(BiCops.tail.up, what = "upper tail dependency",sel.tp)
clim.emptau <- quantity2clim(BiCops.emptau, what = "emperical kendall tau",sel.tp)
clim.theortau <- quantity2clim(BiCops.theortau, what = "theortical kendall tau", sel.tp)

clim.IndTestsP.V81.t2m <- quantity2clim(IndTestsP005.V81.t2m, what = "p.value Ind test",sel.tp)
clim.BiCopsFam.V81.t2m <- quantity2clim(BiCopsFam.V81.t2m, what = "Copula Family",sel.tp)
clim.emptau.V81.t2m.ind <- quantity2clim(BiCops.V81.t2m.emptau, what = "emperical kendall tau",sel.tp)

clim.BiCopsFam.random <-quantity2clim(BiCopsFam.random, what = "Copula Family",sel.tp)

textdf <- cbind(attr(sel.t2m.1,"VertexCoords"),BiCopsFam)
textdf$x[textdf$x<0] <-textdf$x[textdf$x<0] + 360
textdf <- as.matrix(textdf)
textlay <- apply(textdf,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.5,col = "green"))


RColorBrewer::brewer.pal(11,"Spectral")

col.r <- colorRampPalette(brewer.pal(9, "Reds"))(15)
col.b <- colorRampPalette(brewer.pal(9, "Blues"))(15)
col.rb <- c(rev(col.b),'white','white',col.r)

asym.pos.low.tail <- c(3,14,16,114,204)
asym.pos.up.tail <- c(4,6,13,104,214)
symm.pos.two.tails <- c(10,20,7,9,17,19)

col.asym.pos.low.tail <-colorRampPalette(c("orange","darkred"))(7)[3:7]
col.asym.pos.up.tail<- colorRampPalette(c("orange","red"))(6)[2:6]
col.symm.pos.two.tails<- colorRampPalette(c("red","darkred"))(8)[2:7]

#symm.pos.no.tail <- c(1,2)
symm.ind<- c(0)
#symm.neg.no.tail <- c(5)
symm.pos.no.tail <- c(1,2,5)
symm.neg.no.tail <- c(301,302,305)

#col.symm.pos.no.tail <- colorRampPalette(c("white","lightcoral"))(5)[2:3]
col.symm.pos.no.tail <- colorRampPalette(c("white","lightcoral"))(10)[c(2,3,4)]
col.symm.ind <- "white"
#col.symm.neg.no.tail <- colorRampPalette(c("white","blue"))(5)[2]
col.symm.neg.no.tail <- colorRampPalette(c("white","blue"))(10)[c(2,3,4)]

asym.neg.left.tail <- c(24,33,26,124,234)
symm.neg.two.tails <- c(30,40,27,29,37,39)
asym.neg.right.tail <- c(23,34,36,134,224)

col.asym.neg.left.tail <-colorRampPalette(c("lightblue","blue"))(6)[2:6]
col.asym.neg.right.tail<-colorRampPalette(c("lightblue","purple"))(6)[2:6]
col.symm.neg.two.tails<- colorRampPalette(c("blue","purple"))(8)[2:7]


order.tails <- c(rev(asym.pos.low.tail),
                 rev(symm.pos.two.tails),
                 rev(asym.pos.up.tail),
                 rev(symm.pos.no.tail),
                 symm.ind,
                 symm.neg.no.tail,
                 asym.neg.left.tail,
                 symm.neg.two.tails,
                 asym.neg.right.tail)
col.order.tails <- c(rev(col.asym.pos.low.tail),
                     rev(col.symm.pos.two.tails),
                     rev(col.asym.pos.up.tail),
                     rev(col.symm.pos.no.tail),
                     col.symm.ind,
                     col.symm.neg.no.tail,
                     col.asym.neg.left.tail,
                     col.symm.neg.two.tails,
                     col.asym.neg.right.tail)

df.col.order.tails <- data.frame(order.tails,col.order.tails)
length(order.tails)
length(col.order.tails)

unique.fams <- unique(BiCopsFam)[order(unique(BiCopsFam))]
col.fams <- colorRampPalette(c("red","orange","yellow","green","blue"))(length(unique.fams))
col.fams <- c("white",col.fams)

fams1 <-order.tails[order(order.tails)]
cols1<- col.order.tails[order(order.tails)]


plotname <- paste0("figs/copula/cop.families.global.q",b,".png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
fig.BiCopsFam <- spatialPlot(clim.BiCopsFam, 
            backdrop.theme = "coastline", 
            lonCenter = 180, 
            #at = c(-1,unique(BiCopsFam)[order(unique(BiCopsFam))]),
            #col.regions = col.fams,
            #colorkey = list(at =1:length(unique.fams), labels = as.character(unique.fams)),
            at = c(-1,fams1,306),
            col.regions = cols1,
            colorkey = list(at =1:(length(order.tails)+1), 
                            labels = as.character(c(fams1,306))),
            sp.layout = textlay)
fig.BiCopsFam
dev.off()

plotname <- paste0("figs/copula/ind.test.global.q",b,".png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.IndTestsP,
            backdrop.theme = "coastline", 
            lonCenter = 180)
dev.off()

plotname <- paste0("figs/copula/upper.dep.q",b,".png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.tail.up,
            backdrop.theme = "coastline", 
            color.theme = "Reds",
            lonCenter = 180)
dev.off()

plotname <- paste0("figs/copula/lower.dep.q",b,".png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.tail.low,
            backdrop.theme = "coastline", 
            color.theme = "Blues",
            lonCenter = 180)
dev.off()

plotname <- paste0("figs/copula/emptau.global.q",b,".png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
fig.emptau <- spatialPlot(clim.emptau,
            backdrop.theme = "coastline", rev.colors = TRUE,
            col.regions = col.rb,
            lonCenter = 180)
fig.emptau
dev.off()

plotname <- paste0("figs/copula/theortau.global.q",b,".png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.theortau, 
            backdrop.theme = "coastline", 
            col.regions = col.rb,
            lonCenter = 180)
dev.off()
#########################################
#
#########################################
textdf.V81.t2m <- cbind(attr(sel.t2m.1,"VertexCoords"),BiCopsFam.V81.t2m)
textdf.V81.t2m$x[textdf.V81.t2m$x<0] <-textdf.V81.t2m$x[textdf.V81.t2m$x<0] + 360
textdf.V81.t2m <- as.matrix(textdf.V81.t2m)
textlay.V81.t2m <- apply(textdf.V81.t2m,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.5,col = "green"))

plotname <- paste0("figs/copula/cop.families.global.V81.t2m.q",b,".png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
fig.BiCopsFam.V81.t2m <- spatialPlot(clim.BiCopsFam.V81.t2m, 
            backdrop.theme = "coastline", 
            lonCenter = 180, 
            # at = c(-1,unique(BiCopsFam)[order(unique(BiCopsFam))]),
            # col.regions = col.fams,
            # colorkey = list(at =1:length(unique.fams), labels = as.character(unique.fams)),
            at = c(-1,fams1,306),
            col.regions = cols1,
            colorkey = list(at =1:(length(order.tails)+1), 
                            labels = as.character(c(fams1,306))),
            sp.layout = textlay.V81.t2m)
fig.BiCopsFam.V81.t2m
dev.off()

plotname <- paste0("figs/copula/ind.test.global.V81.t2m.q",b,".png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.IndTestsP.V81.t2m,
            backdrop.theme = "coastline", 
            lonCenter = 180)
dev.off()

plotname <- paste0("figs/copula/emptau.global.V81.t2m.q",b,".png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
fig.emptau.V81.t2m <- spatialPlot(clim.emptau.V81.t2m.ind,
            rev.colors = TRUE,
            col.regions = col.rb,
            backdrop.theme = "coastline", 
            lonCenter = 180)
fig.emptau.V81.t2m
dev.off()



textdf.random <- cbind(attr(sel.t2m.1,"VertexCoords"),BiCopsFam.random)
textdf.random$x[textdf.random$x<0] <-textdf.random$x[textdf.random$x<0] + 360
textdf.random <- as.matrix(textdf.random)
textlay.random <- apply(textdf.random,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "green"))

plotname <- paste0("figs/copula/cop.families.global.random.q",b,".png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.BiCopsFam.random, 
            backdrop.theme = "coastline", 
            lonCenter = 180, 
            # at = c(-1,unique(BiCopsFam)[order(unique(BiCopsFam))]),
            # col.regions = col.fams,
            # colorkey = list(at =1:length(unique.fams), labels = as.character(unique.fams)),
            at = c(-1,fams1,306),
            col.regions = cols1,
            colorkey = list(at =1:(length(order.tails)+1), 
                            labels = as.character(c(fams1,306))),
            sp.layout = textlay.random)
dev.off()


########################################
#
########################################
plotname <- paste0("figs/copula/cop.families.emptau.normal.and.condV81.q",b,".pdf")
pdf(plotname,width = 12, height = 7)
grid.arrange(fig.emptau,fig.emptau.V81.t2m,fig.BiCopsFam,fig.BiCopsFam.V81.t2m)
dev.off()

#########################################
# Illustracion Copulas. Yet to adatpt to b general inestad of b = 4
#########################################
which(BiCopsFam == 104)
box.fam <- 89
par(mfrow = c(1,2))
#p <- contour(BiCops.ind[[box.fam]])
q <- plot(BiCops.ind[[box.fam]])

plotname <- paste0("figs/copula/ex.V",box.fam,"_",BiCopsFam[[box.fam]],"_",BiCops.ind[[box.fam]]$familyname,".png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
par(mfrow=c(1,2))
BiCopMetaContour(family = 10)
contour(BiCops.ind[[box.fam]])

print(q, position = c(0.5, 0, 1, 1), more = TRUE)
title(paste0("Box :V",box.fam,"\nFamily ",BiCopsFam[[box.fam]],": ",BiCops.ind[[box.fam]]$familyname))
dev.off()


plotfams <- function(x) {
  contour(BiCops.ind[[which(BiCopsFam == unique.fams[x])[1]]])
  text(x=2.7,y = 3,labels = paste0(round(BiCops.ind[[which(BiCopsFam == unique.fams[x])[1]]]$tau,2)))
  title(paste0(unique.fams[x],": ",gsub("degrees",BiCopName(unique.fams[x],short = FALSE),replacement = "")))
}

nrow <- floor(sqrt(length(unique.fams)))
plotname <- "figs/copula/families.png"
png(plotname, width = 30, height = 30, units = "cm", res = 150)
par(mfrow = c(nrow+1,nrow+1),mar = c(1.5, 1.5, 2, 0.5) + 0.1)
lapply(1:length(unique.fams),plotfams)
dev.off()
