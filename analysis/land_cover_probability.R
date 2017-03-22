
##
## Analysis and prediction of Land Use/Land Cover in WessexBESS
## Created 6 March 2017
## Reto Schmucki
##

## C:\Users\RETOSCHM\OneDrive - Natural Environment Research Council\WessexLULCchange\analysis\land_cover_probability.R

##=======================================================================
## NOTE: 
## Many layers are extracted from from the resulting dataset.
##
## C:\Users\RETOSCHM\OneDrive - Natural Environment Research Council\Wessex_BESSgithub\WessexBess\Build WessexBESS database ## in PostgreSQL.R
##
##=======================================================================

R --vanilla

## Initiate and go to R-project

library(raster)
library(rgdal)
library(INLA)
library(geostatsp)
library(rgeos)
library(RgoogleMaps)
library(PBSmapping)
library(ggplot2)


project_name <- "WessexLULCchange"

if (dir.exists(normalizePath(file.path(Sys.getenv("HOME"),"OneDrive - Natural Environment Research Council",project_name)))) {
	cat(paste('You are all set to work on',project_name,'project!','\n'))
	} else { source(normalizePath(file.path(Sys.getenv("HOME"),"OneDrive - Natural Environment Research Council",project_name,"set_project.R")))
}

p_wd <- normalizePath(file.path(Sys.getenv("HOME"),"OneDrive - Natural Environment Research Council",project_name))
setwd(p_wd)

## Get the original data

file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/wessex_grid.shp",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/wessex_grid.dbf",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/wessex_grid.prj",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/wessex_grid.shx",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)

file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/spta_boundary/spta_bndry.shp",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/spta_boundary/spta_bndry.dbf",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/spta_boundary/spta_bndry.prj",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/spta_boundary/spta_bndry.shx",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)

file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/Proposed_BESS_Boundaryy.shp",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/Proposed_BESS_Boundary.dbf",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/Proposed_BESS_Boundary.prj",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/Proposed_BESS_Boundary.shx",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)

file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/WoodyLinearFeatures50kCLIP.shp",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/WoodyLinearFeatures50kCLIP.dbf",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/WoodyLinearFeatures50kCLIP.prj",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/WoodyLinearFeatures50kCLIP.shx",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)

file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/WorldHeritageSite50kCLIP.shp",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/WorldHeritageSite50kCLIP.dbf",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/WorldHeritageSite50kCLIP.prj",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/WorldHeritageSite50kCLIP.shx",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)

file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/OSroads_50kCLIP.shp",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/OSroads_50kCLIP.dbf",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/OSroads_50kCLIP.prj",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/OSroads_50kCLIP.shx",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)

file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/OSrivers_50kCLIP.shp",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/OSrivers_50kCLIP.dbf",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/OSrivers_50kCLIP.prj",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/OSrivers_50kCLIP.shx",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)

file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/RamsarSite_50kCLIP.shp",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/RamsarSite_50kCLIP.dbf",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/RamsarSite_50kCLIP.prj",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/RamsarSite_50kCLIP.shx",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)

file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/SSSI_50kCLIP.shp",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/SSSI_50kCLIP.dbf",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/SSSI_50kCLIP.prj",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/SSSI_50kCLIP.shx",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)

file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/StampDudleyEdited_50kCLIP.shp",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/StampDudleyEdited_50kCLIP.dbf",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/StampDudleyEdited_50kCLIP.prj",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/StampDudleyEdited_50kCLIP.shx",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
StampDudleyEdited_50kCLIP

file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/land_cover_map_30m.tif",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/elev_dem30mbilinear.tif",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/slope30m.tif",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)
file.copy("P:/NEC04717_WessexBESS/Data/WESSEX_BESS_GIS/WESSEX_BESS_50KmCLIP/wrb_soil_30m.tif",file.path(p_wd,"data"), overwrite = TRUE, recursive = FALSE,
          copy.mode = TRUE, copy.date = FALSE)


grid_crop <- readOGR(file.path(p_wd,"data"),"wessex_grid")
wessex_boundary <- readOGR(file.path(p_wd,"data"),"Proposed_BESS_Boundary")
StampDudley <- readOGR(file.path(p_wd,"data"),"StampDudleyEdited_50kCLIP")

land_use30m <- raster(file.path(p_wd,"data","land_cover_map_30m.tif"))

elev_dem30m <- raster(file.path(p_wd,"data","elev_dem30mbilinear.tif"))
slope_dem30m <- raster(file.path(p_wd,"data","slope30m.tif"))
soil_wrb <- raster(file.path(p_wd,"data","wrb_soil_30m.tif"))

wood_linear <- readOGR(file.path(p_wd,"data"),"WoodyLinearFeatures50kCLIP")
road <- readOGR(file.path(p_wd,"data"),"OSroads_50kCLIP")
river <- readOGR(file.path(p_wd,"data"),"OSrivers_50kCLIP")

spta_bndry <- readOGR(file.path(p_wd,"data"),"spta_bndry")
world_heritage <- readOGR(file.path(p_wd,"data"),"WorldHeritageSite50kCLIP")
ramsar <- readOGR(file.path(p_wd,"data"),"RamsarSite_50kCLIP")
sssi <- readOGR(file.path(p_wd,"data"),"SSSI_50kCLIP")

road <- road[as.character(road@data$class) %in% c('A Road','B Road','Motorway'),]

river_mask <- mask(land_use30m, gBuffer(river,width=50))
river_mask[!is.na(river_mask)] <- 1

road_mask <- mask(land_use30m,gBuffer(road,width=50))
road_mask[!is.na(road_mask)] <- 1

projection(spta_bndry) <- crs(land_use30m)
spta_bndry_mask <- mask(land_use30m, spta_bndry)
spta_bndry_mask[!is.na(spta_bndry_mask)] <- 1

projection(world_heritage) <- crs(land_use30m)
world_heritage_mask <- mask(land_use30m, world_heritage)
world_heritage_mask[!is.na(world_heritage_mask)] <- 1

projection(sssi) <- crs(land_use30m)
sssi_mask <- mask(land_use30m, sssi)
sssi_mask[!is.na(sssi_mask)] <- 1

projection(ramsar) <- crs(land_use30m)
ramsar_mask <- mask(land_use30m, ramsar)
ramsar_mask[!is.na(ramsar_mask)] <- 1

if (!file.exists(file.path(p_wd,"data","urban_poly.shp"))){
urban_poly <- rasterToPolygons(land_use30m,fun=function(x){x==35},dissolve=TRUE)
urban_poly <- disaggregate(urban_poly)
writeOGR(urban_poly,dsn=file.path(p_wd,"data"),layer="urban_poly",driver="ESRI Shapefile")
}else{
urban_poly <- readOGR(file.path(p_wd,"data"),"urban_poly")
}

urban_poly_subset <-urban_poly[area(urban_poly) >= 250000,] # larger than 500 X 500 meters
projection(urban_poly_subset) <- crs(land_use30m) 
urban_mask <- mask(land_use30m, urban_poly_subset)

total_mask <- stackApply(stack(river_mask, road_mask, spta_bndry_mask, world_heritage_mask, ramsar_mask, sssi_mask, urban_mask),1,fun=sum,na.rm=TRUE) 
total_mask[total_mask>0] <- 100
total_mask[total_mask==0] <- NA

save.image(file="data/lulc.RData")

load("data/lulc.RData")

png(filename = "figs/wessex_map.png",,width = 7, height = 7, units = 'in', res = 300,bg = "transparent")
b <- as.data.frame(t(bbox(land_use30m)))
coordinates(b) <- ~ s1 + s2
crs(b)<- sp::CRS("+init=epsg:27700")
b <- spTransform(b,sp::CRS("+init=epsg:4326"))
lat <- coordinates(b)[,2]
long <- coordinates(b)[,1]
wessexMap <- MapBackground(lat=lat,lon=long,zoom=1)
wessex_limit <- maptools::SpatialPolygons2PolySet(spTransform(wessex_boundary,sp::CRS("+init=epsg:4326")))
bb <- as(extent(land_use30m),'SpatialPolygons')
crs(bb) <- sp::CRS("+init=epsg:27700")
wessex_wide <- maptools::SpatialPolygons2PolySet(spTransform(bb,sp::CRS("+init=epsg:4326")))

PlotPolysOnStaticMap(wessexMap,wessex_limit,col=rgb(1,0,1,0.1),add=FALSE)
PlotPolysOnStaticMap(wessexMap,wessex_wide,col=rgb(0,0,0,0),border=rgb(1,0,0,1),lwd=1.5,lty=2,add=TRUE)
dev.off()

png(filename = "figs/wessex_lcover.png",width = 7, height = 7, units = 'in', res = 300,bg = "transparent") ##,res = 300)
plot(land_use30m,col=c(rep("wheat3",9),grDevices::terrain.colors(34)[10:34],"grey50"),)
dev.off()


png(filename = "figs/wessex_mask.png",width = 7, height = 7, units = 'in', res = 300,bg = "transparent")
plot(land_use30m,col=c(rep("wheat3",9),grDevices::terrain.colors(34)[10:34],"grey50"),xlab="longitude",ylab="latitude",legend=FALSE)
plot(river,col='darkblue',add=TRUE)
plot(road,col='black',add=TRUE)
plot(spta_bndry,lwd=1,lty=2,border='blue',col=rgb(0,0,1,0.2),add=TRUE)
plot(urban_poly_subset,border=rgb(1,0,1,0.7),col=rgb(1,0,1,0.7),lwd=0.5,add=TRUE)
plot(world_heritage,border=FALSE,col=rgb(0.4,0.3,1,0.8),add=TRUE)
plot(ramsar,border=FALSE,col=rgb(0,1,1,0.7),add=TRUE)
dev.off()

png(filename = "figs/wessex_maskgrey.png",width = 7, height = 7, units = 'in', res = 300,bg = "transparent")
plot(land_use30m,col=c(rep("wheat3",9),grDevices::terrain.colors(34)[10:34],"grey50"),xlab="longitude",ylab="latitude",legend=FALSE)
plot(total_mask,col='grey10',add=TRUE,border=FALSE,legend=FALSE)
dev.off()

png(filename = "figs/wessex_maskDudley.png",width = 14, height = 7, units = 'in', res = 300,bg = "transparent")
par(mfrow=c(1,2))
plot(land_use30m,col=c(rep("wheat3",9),grDevices::terrain.colors(34)[10:34],"grey50"),xlab="longitude",ylab="latitude",legend=FALSE,axes=FALSE,bty="n")
plot(total_mask,col='grey10',add=TRUE,border=FALSE,legend=FALSE)
plot(StampDudley,col=c(rep("wheat3",9),grDevices::terrain.colors(34)[10:34],"grey50")[c(15,35,33,5,35,12,21,11)][round(as.numeric(StampDudley@data$GRIDCODE))],border=FALSE,xlab="longitude",ylab="latitude")
plot(total_mask,col='grey10',add=TRUE,border=FALSE,legend=FALSE)
dev.off()



##crop one area
region <- c(22:27,32:37,42:47,52:57,62:67,72:77)

lu_crop <- crop(land_use30m,grid_crop[grid_crop@data$FID %in% region,])
elev_crop <- crop(elev_dem30m,grid_crop[grid_crop@data$FID %in% region,])
slope_crop <- crop(slope_dem30m,grid_crop[grid_crop@data$FID %in% region,])
soil_crop <- crop(soil_wrb,grid_crop[grid_crop@data$FID %in% region,])
mask_crop <- crop(total_mask,grid_crop[grid_crop@data$FID %in% region,])

agri <- lu_crop
agri[agri %in% c(1:9)] <- 1
agri[agri > 9] <- 2

agri_nomask <- stackApply(stack(agri,mask_crop),1,fun=sum,na.rm=TRUE)
agri_nomask.pt <- rasterToPoints(agri_nomask,spatial=TRUE)
agri_nomask1.pt <- agri_nomask.pt[agri_nomask.pt@data$index_1==1,]
agri_nomask0.pt <- agri_nomask.pt[agri_nomask.pt@data$index_1==2,]

agri.pt <- agri_nomask1.pt[sample(c(1:dim(agri_nomask1.pt)[1]),500,replace=FALSE),]
projection(agri.pt) <- crs(lu_crop)

nonagri.pt <- agri_nomask0.pt[sample(c(1:dim(agri_nomask0.pt)[1]),500,replace=FALSE),]
projection(nonagri.pt) <- crs(lu_crop)

agri.lc.pt <- bind(agri.pt,nonagri.pt)

png(filename = "figs/inla_sampling.png",width = 7, height = 7, units = 'in', res = 300,bg = "transparent")
plot(lu_crop,col=c(rep("wheat3",9),grDevices::terrain.colors(34)[10:34],"grey50"),xlab="longitude",ylab="latitude",legend=FALSE,)
image(mask_crop,col='grey50',add=TRUE)
plot(agri.pt,col='red',pch=19,cex=0.7,add=TRUE)
plot(nonagri.pt,col='blue',pch=19,cex=0.7,add=TRUE)
dev.off()

png(filename = "figs/inla_variable.png",width = 14, height = 14, units = 'in', res = 300,bg = "transparent")
par(mfrow=c(2,2))
plot(lu_crop,col=c(rep("wheat3",9),grDevices::terrain.colors(34)[10:34],"grey50"),xlab="longitude",ylab="latitude",legend=FALSE,)
image(mask_crop,col='grey50',add=TRUE)
plot(agri.pt,col='red',pch=19,cex=0.7,add=TRUE)
plot(nonagri.pt,col='blue',pch=19,cex=0.7,add=TRUE)
plot(soil_crop,col=colorRamps::primary.colors(15),xlab="longitude",ylab="latitude",legend=FALSE,main='soil')
plot(elev_crop,xlab="longitude",ylab="latitude",main='elevation')
plot(slope_crop,xlab="longitude",ylab="latitude",main='slope')
dev.off()

set.seed(6543)
coords <- agri.lc.pt[sample(1:1000,250,replace=FALSE),]
sp.mesh <- inla.mesh.2d(loc=coords,cutoff=500,max.edge=c(1500,5000))
plot(sp.mesh)

agri.spde <- inla.spde2.matern(mesh=sp.mesh,alpha=2)

est.index <- sample(c(1:length(agri.lc.pt)),length(agri.lc.pt),replace=FALSE)[1:round(0.7*length(agri.lc.pt))]
val.index <- c(1:length(agri.lc.pt))[!(1:length(agri.lc.pt) %in% est.index)]

est.coords <- agri.lc.pt@coords[est.index,]
val.coords <- agri.lc.pt@coords[val.index,]

agri.est <- rep(c(1,0),c(500,500))[est.index]
agri.val <- rep(c(1,0),c(500,500))[val.index]

soil.est <- extract(soil_crop,agri.lc.pt[est.index,])
soil.val <- extract(soil_crop,agri.lc.pt[val.index,])

elev.est <- extract(elev_crop,agri.lc.pt[est.index,])
elev.val <- extract(elev_crop,agri.lc.pt[val.index,])

slope.est <- extract(slope_crop,agri.lc.pt[est.index,])
slope.val <- extract(slope_crop,agri.lc.pt[val.index,])

agri.spde <- inla.spde2.matern(mesh=sp.mesh,alpha=2)

A.est <- inla.spde.make.A(mesh=sp.mesh, loc=est.coords)
A.val <- inla.spde.make.A(mesh=sp.mesh, loc=val.coords)

A.pred <- inla.spde.make.A(mesh=sp.mesh)

s.index <- inla.spde.make.index(name="spatial.field",
								n.spde=agri.spde$n.spde)

stack.est <- inla.stack(data=list(agri=agri.est),
					A= list(A.est,1,1,1),
					effects= list(c(s.index, list(Intercept=1)),list(soil=factor(soil.est)),
								list(elevation=elev.est),list(slope=slope.est)),
					tag="est")

stack.val <- inla.stack(data=list(agri=NA),
					A= list(A.val,1,1,1),
					effects= list(c(s.index, list(Intercept=1)),list(soil=factor(soil.val)),
								list(elevation=elev.val),list(slope=slope.val)),
					tag='val')

stack.pred <- inla.stack(data=list(agri=NA),
					A= list(A.pred),
					effects= list(c(s.index, list(Intercept=1))),
					tag='pred')

stack.join <- inla.stack(stack.est,stack.val,stack.pred)

formula <- agri ~ -1 + Intercept + soil + elevation + slope + f(spatial.field, model=spde)
formula.2 <- agri ~ -1 + Intercept + soil + slope + f(spatial.field, model=spde)
formula.3 <- agri ~ -1 + Intercept + f(spatial.field, model=spde)

agri.out <- inla(formula.1,
			data=inla.stack.data(stack.join,spde=agri.spde),
			family="binomial",Ntrials=1,
			control.predictor=list(A=inla.stack.A(stack.join), compute=TRUE),
			control.fixed = list(expand.factor.strategy='inla'),
			control.compute=list(dic=TRUE))

agri.out2 <- inla(formula.2,
			data=inla.stack.data(stack.join,spde=agri.spde),
			family="binomial",Ntrials=1,
			control.predictor=list(A=inla.stack.A(stack.join), compute=TRUE),
			control.fixed = list(expand.factor.strategy='inla'),
			control.compute=list(dic=TRUE))

agri.out3 <- inla(formula.3,
			data=inla.stack.data(stack.join,spde=agri.spde),
			family="binomial",Ntrials=1,
			control.predictor=list(A=inla.stack.A(stack.join), compute=TRUE),
			control.fixed = list(expand.factor.strategy='inla'),
			control.compute=list(dic=TRUE))


save.image(file="data/lulc.RData")

agri.out <- agri.out3

index.val <- inla.stack.index(stack.join,'val')$data
post.agri.val_mean <- agri.out$summary.linear.predictor[index.val,'mean']
post.agri.val_sd <- agri.out$summary.linear.predictor[index.val,'sd']

invlogit <- function(x) {
						invl <- exp(x)/(1+exp(x))
						 return(invl)
}

t.post.agri.val_mean <- invlogit(post.agri.val_mean)

validation <- data.frame(obs=as.factor(agri.val),est.prob=t.post.agri.val_mean)

X <- validation$est.prob[validation$obs==1]
O <- validation$est.prob[validation$obs==0]

png(filename = "figs/inla_predict_spatialfield.png",width = 14, height = 7, units = 'in', res = 300,bg = "transparent")
par(mfrow=c(1,2))
hist(X, prob=TRUE, col=rgb(0,0,1,0.2),main='agricultural')
lines(density(X, adjust=2), lty=1, lwd=2, col='blue')
hist(O, prob=TRUE, col=rgb(1,0.3,0,0.5),main='non agricultural')
lines(density(O, adjust=2), lty=1, lwd=2, col=rgb(1,0,0,1))
dev.off()

index.pred <- inla.stack.index(stack.join,'pred')$data
post.mean.pred <- agri.out3$summary.linear.predictor[index.pred,'mean']
post.mean.pred <- invlogit(post.mean.pred)

seq.x.grid <- seq(from=extent(slope_crop)[1],to=extent(slope_crop)[2],by=1000)
seq.y.grid <- seq(from=extent(slope_crop)[3],to=extent(slope_crop)[4],by=1000)

pred.grid <- as.matrix(expand.grid(x=seq.x.grid,y=seq.y.grid))

proj.grid <- inla.mesh.projector(sp.mesh,
				xlim=range(pred.grid[,1]),
				ylim=range(pred.grid[,2]),
				dims=c(length(seq.x.grid),length(seq.y.grid)))

post.mean.pred.grid <- inla.mesh.project(proj.grid,post.mean.pred)

png(filename = "figs/inla_spatialfield.png",width = 7, height = 7, units = 'in', res = 300,bg = "transparent")
latent.field <- raster(as.matrix(post.mean.pred.grid))
plot(latent.field,axes=FALSE,xlab="longitude",ylab="latitude")
dev.off()



round(agri.out$summary.fixed,3)