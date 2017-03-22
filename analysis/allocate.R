R --vanilla

## Initiate and go to R-project

library(raster)
library(rgdal)
library(rgeos)



p_wd <- normalizePath("P:\\NEC04717_WessexBESS\\Data\\WESSEX_BESS_GIS\\WESSEX_BESS_50KmCLIP")
setwd(normalizePath("P:\\NEC04717_WessexBESS\\Data\\WESSEX_BESS_GIS\\WESSEX_BESS_50KmCLIP"))

lcmap <- readOGR(file.path(p_wd),"LCM2007_50kCLIP")

box_crp <- extent(lcmap)-100000
lcmap_crp <- crop(lcmap,box_crp)

r <- raster(extent(lcmap_crp))
projection(r) <- proj4string(lcmap_crp)
res(r) <- 100

lcmap_r <- rasterize(lcmap_crp,r,field='objectid')
lcmap.dt <- data.table::data.table(rasterToPoints(lcmap_r))


layer_array <- array(NA,dim=c(dim(as.matrix(lcmap_r)),length(unique(values(lcmap_r)))))

for (i in 1:dim(layer_array)[3]){
prob.layer <- data.table::data.table(layer=unique(lcmap.dt$layer), prob=sample(seq(from=0,to=1,by=0.01),length(unique(lcmap.dt$layer)),replace=TRUE))
val.layer <- data.table::data.table(layer=values(lcmap_r))
ord_prob <- merge(val.layer,prob.layer,by="layer",sort=FALSE)

lcmap_r2 <- lcmap_r
lcmap_r2[] <- ord_prob$prob

layer_array[,,i] <- as.matrix(lcmap_r2)
}

r1 <- raster(layer_array[,,1])
extent(r1) <- extent(lcmap_r2)
projection(r1) <- proj4string(lcmap_crp)
plot(r1)
plot(lcmap_crp,add=TRUE)

dev.new()
r1 <- raster(layer_array[,,2])
extent(r1) <- extent(lcmap_r2)
projection(r1) <- proj4string(lcmap_crp)
plot(r1)
plot(lcmap_crp,add=TRUE)

dev.new()
lcmap_r <- rasterize(lcmap_crp,r,field='intcode')
plot(lcmap_r,col=RColorBrewer::brewer.pal(length(unique(values(lcmap_r))),'Set3'))
plot(lcmap_crp,add=TRUE)

lcmap_lu <- as.matrix(lcmap_r)

unique(values(lcmap_r))[1]

lcmap.l <- lcmap_r; lcmap.l[]<- 1
lcmap.l[lcmap_r==unique(values(lcmap_r))[1]] <- 0


prob.chg <- as.matrix(lcmap.l)*layer_array[,,1]

prob.chg_r <- raster(prob.chg)
extent(prob.chg_r ) <- extent(lcmap_r2)
projection(prob.chg_r ) <- proj4string(lcmap_crp)

plot(prob.chg_r)
extr.prob.map <- extract(prob.chg_r,lcmap_crp,fun=mean)

lcmap_crp@data$prob <-extr.prob.map
plot(lcmap_crp,col=lcmap_crp@data$prob)

colfunc<-RColorBrewer::brewer.pal(11,'Spectral')

plot(prob.chg_r,col=colfunc[findInterval(lcmap_crp@data$prob,seq(from=0,to=1,by=0.1))])
dev.new()
plot(lcmap_crp,col=colfunc[findInterval(lcmap_crp@data$prob,seq(from=0,to=1,by=0.1))])

which(lcmap_crp@data$prob==max(lcmap_crp@data$prob,na.rm=TRUE))


dev.new()
plot(lcmap_crp,col=colfunc[findInterval(lcmap_crp@data$prob,seq(from=0,to=1,by=0.1))])
plot(lcmap_crp[which(lcmap_crp@data$prob==max(lcmap_crp@data$prob,na.rm=TRUE)),],col=colfunc[1],add=TRUE)







lcmap.dt <- merge(lcmap.dt,prob.layer,by="layer")


plot(lcmap_r)
plot(lcmap_crp,add=TRUE)

length(unique(data.frame(lcmap_crp)$fieldcode))
length((data.frame(lcmap_crp)$fieldcode))