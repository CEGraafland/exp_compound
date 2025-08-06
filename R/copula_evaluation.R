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
quantile(df.t2m$V81)[4]
ind.V81.t2m.q4 <- which(df.t2m$V81>=quantile(df.t2m$V81)[4])
df.t2m.V81.t2m.q4 <- df.t2m[ind.V81.t2m.q4,]
df.tp.V81.t2m.q4 <- df.tp[ind.V81.t2m.q4,]

# random quantile:
set.seed(1)
sample.ind <- sample(1:nrow(df.t2m),nrow(df.t2m)/4)
random.q.sample.ind <- sample.ind
save(random.q.sample.ind, file = "results/copulas/random.q.sample.ind")
df.t2m.random.q <- df.t2m[random.q.sample.ind,]
df.tp.random.q <- df.tp[random.q.sample.ind,]
########################################
# psuedo observaciones 
########################################
p.df.tp<- pobs(df.tp)
p.df.t2m <- pobs(df.t2m)

p.df.tp.V81.t2m.q4<- pobs(df.tp.V81.t2m.q4)
p.df.t2m.V81.t2m.q4 <- pobs(df.t2m.V81.t2m.q4)

p.df.t2m.random.q <- pobs(df.t2m.random.q)
p.df.tp.random.q <- pobs(df.tp.random.q)

#bi2<- BiCopSelect(p.df.tp[,2], p.df.t2m[,2], familyset = NA)
# plot(BiCops[[5]])
# BiCops[[646]]$p.value.indeptest
dev.off()

x_test <- 82
 plot(p.df.tp[,x_test],p.df.t2m[,x_test])
 plot(p.df.tp.V81.t2m.q4[,x_test],p.df.t2m.V81.t2m.q4[,x_test])
 plot(df.tp[,x_test],df.t2m[,x_test])
 plot(df.tp.V81.t2m.q4[,x_test],df.t2m.V81.t2m.q4[,x_test],ylim = c(min(df.t2m[,x_test]),max(df.t2m[,x_test])))
 
#########################################
# Estimate bivariate copulas
#########################################
# Without inicial ind. test
# BiCops <- lapply(1:ncol(p.df.tp),FUN = function(i)BiCopSelect(p.df.tp[,i], p.df.t2m[,i], familyset = NA))
# With inicial ind.test
# BiCops.ind <- lapply(1:ncol(p.df.tp),FUN = function(i)BiCopSelect(p.df.tp[,i], p.df.t2m[,i], familyset = NA,indeptest = TRUE))
# save(BiCops.ind,file = "results/copulas/BiCops.ind.rda")
load(file = "results/copulas/BiCops.ind.rda")

#BiCops.V81.t2m.q4.ind <- lapply(1:ncol(p.df.tp.V81.t2m.q4),FUN = function(i)BiCopSelect(p.df.tp.V81.t2m.q4[,i], p.df.t2m.V81.t2m.q4[,i], familyset = NA,indeptest = TRUE))
#save(BiCops.V81.t2m.q4.ind,file = "results/copulas/BiCops.V81.t2m.q4.ind.rda")
load(file = "results/copulas/BiCops.V81.t2m.q4.ind.rda")



#BiCops.V81.t2m.randomq4.ind <- lapply(1:ncol(p.df.tp.random.q),FUN = function(i)BiCopSelect(p.df.tp.random.q[,i], p.df.t2m.random.q[,i], familyset = NA,indeptest = TRUE))
#save(BiCops.V81.t2m.randomq4.ind,file = "results/copulas/BiCops.V81.t2m.randomq4.ind.rda")
load(file = "results/copulas/BiCops.randomq4.ind.rda")
#BiCops.randomq4.ind <- BiCops.V81.t2m.randomq4.ind
#save(BiCops.randomq4.ind,file = "results/copulas/BiCops.randomq4.ind.rda")


# Only ind.tests
IndTests<- lapply(1:ncol(p.df.tp),FUN = function(i)BiCopIndTest(p.df.tp[,i], p.df.t2m[,i]))
IndTests.V81.t2m.q4<- lapply(1:ncol(p.df.tp.V81.t2m.q4),FUN = function(i)BiCopIndTest(p.df.tp.V81.t2m.q4[,i], p.df.t2m.V81.t2m.q4[,i]))
IndTests.random.q <- lapply(1:ncol(p.df.tp.random.q),FUN = function(i)BiCopIndTest(p.df.tp.random.q[,i], p.df.t2m.random.q[,i]))
##################################
# extraer caracteristicas copulas
##################################
BiCops.V81.t2m.q4.ind[[621]]$par
BiCopsFam <- sapply(BiCops.ind,function(x) x$family)
BiCopsFam.V81.t2m.q4 <- sapply(BiCops.V81.t2m.q4.ind,function(x) x$family)
BiCopsFam.random.q <- sapply(BiCops.randomq4.ind,function(x) x$family)


gaus.min <- sapply(BiCops.ind,function(x) x$family ==1 && x$par<0)
BiCopsFam[gaus.min]<- 301
t.min <- sapply(BiCops.ind,function(x) x$family ==2 && x$par<0)
BiCopsFam[t.min]<- 302
frank.min <- sapply(BiCops.ind,function(x) x$family ==5 && x$par<0)
BiCopsFam[frank.min]<- 305

gaus.min.V81.t2m.q4 <- sapply(BiCops.V81.t2m.q4.ind,function(x) x$family ==1 && x$par<0)
BiCopsFam.V81.t2m.q4[gaus.min.V81.t2m.q4]<- 301
t.min.V81.t2m.q4 <- sapply(BiCops.V81.t2m.q4.ind,function(x) x$family ==2 && x$par<0)
BiCopsFam.V81.t2m.q4[t.min.V81.t2m.q4]<- 302
frank.min.V81.t2m.q4 <- sapply(BiCops.V81.t2m.q4.ind,function(x) x$family ==5 && x$par<0)
BiCopsFam.V81.t2m.q4[frank.min.V81.t2m.q4]<- 305

gaus.min.randomq4 <- sapply(BiCops.randomq4.ind,function(x) x$family ==1 && x$par<0)
BiCopsFam.random.q[gaus.min.V81.t2m.q4]<- 301
t.min.randomq4 <- sapply(BiCops.randomq4.ind,function(x) x$family ==2 && x$par<0)
BiCopsFam.random.q[t.min.V81.t2m.q4]<- 302
frank.min.randomq4 <- sapply(BiCops.randomq4.ind,function(x) x$family ==5 && x$par<0)
BiCopsFam.random.q[frank.min.V81.t2m.q4]<- 305

BiCops.tail.low <- sapply(BiCops.ind,function(x) BiCopPar2TailDep(obj = x)$lower)
BiCops.tail.up <- sapply(BiCops.ind,function(x) BiCopPar2TailDep(obj = x)$upper)
BiCops.emptau <- sapply(BiCops.ind,function(x) x$emptau)
BiCops.V81.t2m.q4.emptau <- sapply(BiCops.V81.t2m.q4.ind,function(x) x$emptau)
BiCops.theortau <- sapply(BiCops.ind,function(x) BiCopPar2Tau(family = x$family, par = x$par, par2 = x$par2))
cor(BiCops.tail.low,we.extremes.T2M)
IndTestsP <- sapply(IndTests,function(x) x$p.value)
IndTestsP005 <- sapply(IndTests,function(x) as.numeric(x$p.value <= 0.05))

IndTestsP.V81.t2m.q4 <- sapply(IndTests.V81.t2m.q4,function(x) x$p.value)
IndTestsP005.V81.t2m.q4 <- sapply(IndTests.V81.t2m.q4,function(x) as.numeric(x$p.value <= 0.05))

# Small-> reject 0 of independence -> dependent  -- TRUE = 1 = dependent
# high -> acept 0 of independence -> indepdent -- FALSE = 0 = independent

clim.IndTestsP <- quantity2clim(IndTestsP005, what = "p.value Ind test",sel.tp)
clim.BiCopsFam <- quantity2clim(BiCopsFam, what = "Copula Family",sel.tp)
clim.tail.low <- quantity2clim(BiCops.tail.low, what = "lower tail dependency",sel.tp)
clim.tail.up <- quantity2clim(BiCops.tail.up, what = "upper tail dependency",sel.tp)
clim.emptau <- quantity2clim(BiCops.emptau, what = "emperical kendall tau",sel.tp)
clim.theortau <- quantity2clim(BiCops.theortau, what = "theortical kendall tau", sel.tp)

clim.IndTestsP.V81.t2m.q4 <- quantity2clim(IndTestsP005.V81.t2m.q4, what = "p.value Ind test",sel.tp)
clim.BiCopsFam.V81.t2m.q4 <- quantity2clim(BiCopsFam.V81.t2m.q4, what = "Copula Family",sel.tp)
clim.emptau.V81.t2m.q4.ind <- quantity2clim(BiCops.V81.t2m.q4.emptau, what = "emperical kendall tau",sel.tp)

clim.BiCopsFam.random.q <-quantity2clim(BiCopsFam.random.q, what = "Copula Family",sel.tp)

textdf <- cbind(attr(sel.t2m.1,"VertexCoords"),BiCopsFam)
textdf$x[textdf$x<0] <-textdf$x[textdf$x<0] + 360
textdf <- as.matrix(textdf)
textlay <- apply(textdf,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "green"))


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


plotname <- paste0("figs/copula/cop.families.global.png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.BiCopsFam, 
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
dev.off()

plotname <- paste0("figs/copula/ind.test.global.png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.IndTestsP,
            backdrop.theme = "coastline", 
            lonCenter = 180)
dev.off()

plotname <- paste0("figs/copula/upper.dep.png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.tail.up,
            backdrop.theme = "coastline", 
            color.theme = "Reds",
            lonCenter = 180)
dev.off()

plotname <- paste0("figs/copula/lower.dep.png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.tail.low,
            backdrop.theme = "coastline", 
            color.theme = "Blues",
            lonCenter = 180)
dev.off()

plotname <- paste0("figs/copula/emptau.global.png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.emptau,
            backdrop.theme = "coastline", rev.colors = TRUE,
            col.regions = col.rb,
            lonCenter = 180)
dev.off()

plotname <- paste0("figs/copula/theortau.global.png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.theortau, 
            backdrop.theme = "coastline", 
            col.regions = col.rb,
            lonCenter = 180)
dev.off()
#########################################
#
#########################################
textdf.V81.t2m.q4 <- cbind(attr(sel.t2m.1,"VertexCoords"),BiCopsFam.V81.t2m.q4)
textdf.V81.t2m.q4$x[textdf.V81.t2m.q4$x<0] <-textdf.V81.t2m.q4$x[textdf.V81.t2m.q4$x<0] + 360
textdf.V81.t2m.q4 <- as.matrix(textdf.V81.t2m.q4)
textlay.V81.t2m.q4 <- apply(textdf.V81.t2m.q4,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "green"))

plotname <- paste0("figs/copula/cop.families.global.V81.t2m.q4.png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.BiCopsFam.V81.t2m.q4, 
            backdrop.theme = "coastline", 
            lonCenter = 180, 
            # at = c(-1,unique(BiCopsFam)[order(unique(BiCopsFam))]),
            # col.regions = col.fams,
            # colorkey = list(at =1:length(unique.fams), labels = as.character(unique.fams)),
            at = c(-1,fams1,306),
            col.regions = cols1,
            colorkey = list(at =1:(length(order.tails)+1), 
                            labels = as.character(c(fams1,306))),
            sp.layout = textlay.V81.t2m.q4)
dev.off()

plotname <- paste0("figs/copula/ind.test.global.V81.t2m.q4.png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.IndTestsP.V81.t2m.q4,
            backdrop.theme = "coastline", 
            lonCenter = 180)
dev.off()

plotname <- paste0("figs/copula/emptau.global.V81.t2m.q4.png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.emptau.V81.t2m.q4.ind,
            rev.colors = TRUE,
            col.regions = col.rb,
            backdrop.theme = "coastline", 
            lonCenter = 180)
dev.off()



textdf.random.q <- cbind(attr(sel.t2m.1,"VertexCoords"),BiCopsFam.random.q)
textdf.random.q$x[textdf.random.q$x<0] <-textdf.random.q$x[textdf.random.q$x<0] + 360
textdf.random.q <- as.matrix(textdf.random.q)
textlay.random.q <- apply(textdf.random.q,MARGIN = 1,FUN = function(x) list("sp.text",x[-3],as.character(x[3]),cex = 0.75,col = "green"))

plotname <- paste0("figs/copula/cop.families.global.random.q.png")
png(plotname,width = 18, height = 10,units = "cm", res = 150)
spatialPlot(clim.BiCopsFam.random.q, 
            backdrop.theme = "coastline", 
            lonCenter = 180, 
            # at = c(-1,unique(BiCopsFam)[order(unique(BiCopsFam))]),
            # col.regions = col.fams,
            # colorkey = list(at =1:length(unique.fams), labels = as.character(unique.fams)),
            at = c(-1,fams1,306),
            col.regions = cols1,
            colorkey = list(at =1:(length(order.tails)+1), 
                            labels = as.character(c(fams1,306))),
            sp.layout = textlay.random.q)
dev.off()




#########################################
# Illustracion Copulas.
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
##########################
# template text in spplot
##########################

sp.label <- function(x, label) {list("sp.text", coordinates(x), label,cex=0.5, col="green")}
NUMB.sp.label <- function(x) {sp.label(x, as.vector(x@data$NUMB))}
make.NUMB.sp.label <- function(x) {do.call("list", NUMB.sp.label(x))}
        

sp.label(clim.BiCopsFam)
sp.label <- function(x) {list("sp.text", expand.grid(x$xyCoords$y, x$xyCoords$x)[2:1], label = as.character(as.vector(redim(x,drop = TRUE)$Data)),cex=0.5, col="green")}
make.NUMB.sp.label <- function(x) {do.call("list", sp.label(x))}    
make.NUMB.sp.label(clim.BiCopsFam)

coordinates(clim.BiCopsFam$xyCoords)
sp.label(coordinates(clim.BiCopsFam$xyCoords),)
            
mapply(FUN = function(x,y)BiCopSelect(x,y,familyset = NA),x = p.df.tp,y = p.df.t2m)
apply


mapply(function(x, y) seq_len(x) + y,
       c(a =  1, b = 2, c = 3),  # names from first
       c(A = 10, B = 0, C = -10))
mapply(rep, 1:4, 4:1)
