# BASE --------------------------------------------------------------------


# multispectral imagery processing
# vegetation indices calculation

# sr
# ndvi
# gndvi
# savi
# evi
# mcari
# ndre
# exg
# exr
# exgexr

# CRS 4326 : WGS 84

# FILE SYSTEM:

# 'red' == source band rasters (dates)
# 'shp' == source vectors for zonal statistics
# 'VIs' == vegetation indices rasters stored (dates)
# 'final tabs' == final tables as output of zonal statistics (dates)

# PACKAGES:

# list.of.packages <- c("sf","raster","mapview","whitebox","randomcoloR","leaflet","Rcpp")
# install.packages(list.of.packages)
# install.packages("devtools")
# install.packages('raster', repos='https://rspatial.r-universe.dev')
# library(devtools)
# install_github("r-spatial/sf")

library(sf)
library(raster)
library(mapview)
library(randomcoloR)
library(leaflet)
require(tidyverse)


# VARIANTY ----------------------------------------------------------------

# datum <- "m20210521"
# datum <- "v20210521"
datum <- "v20210609"

if(startsWith(datum, "m")){
  pole <- "male"
  } else if(startsWith(datum, "v")){
    pole <- "velke"
  } else{
    print("ZEROOO")
  }
    

# LOAD --------------------------------------------------------------------


cesta <- paste0("red/", datum,"/") # změnit složku podle data

red <- raster(paste0(cesta, "result_Red.tif"))
green <- raster(paste0(cesta, "result_Green.tif"))
blue <- raster(paste0(cesta, "result_Blue.tif"))
nir <- raster(paste0(cesta, "result_NIR.tif"))
edge <- raster(paste0(cesta, "result_RedEdge.tif"))

# check load and crs
plot(blue)
blue@crs
graphics.off()


# INDICES CALCULATION -----------------------------------------------------


# SR 
# Simple Ratio
# NIR/R

sr <- nir/red
writeRaster(sr, filename = paste0("VIs/", datum, "_sr.tif"), format="GTiff", overwrite=TRUE)

# NDVI 
# Normalized Difference Vegetation Index
# (NIR - R)/(NIR + R)

ndvi <- (nir-red)/(nir+red)
writeRaster(ndvi, filename = paste0("VIs/", datum, "_ndvi.tif"), format="GTiff", overwrite=TRUE)

# GNDVI 
# Green Normalized Difference Vegetation Index
# (NIR - Green)/(NIR + Green)

gndvi <- (nir - green)/(nir + green)
writeRaster(gndvi, filename = paste0("VIs/", datum, "_gndvi.tif"), format="GTiff", overwrite=TRUE)

# SAVI 
# Soil Adjusted Vegetation Index
# (1.5*(nir-red))/(nir+red+0.5)

savi <- (1.5*(nir-red))/(nir+red+0.5)
writeRaster(savi, filename = paste0("VIs/", datum, "_savi.tif"), format="GTiff", overwrite=TRUE)

# EVI 
# Enhanced Vegetation Index
# 2.5*(NIR-R)/(NIR + 6*R - 7.5*B + 1)

evi <- 2.5*(nir-red)/(nir + 6*red - 7.5*blue + 1)
writeRaster(evi, filename = paste0("VIs/", datum, "_evi.tif"), format="GTiff", overwrite=TRUE)

# MCARI 
# Modified Chlorophyll Absorption Ratio Index
# ((NIR-R)-0.2*(NIR-G))*(NIR/R)

mcari <- ((nir-red)-0.2*(nir-green))*(nir/red)
writeRaster(mcari, filename = paste0("VIs/", datum, "_mcari.tif"), format="GTiff", overwrite=TRUE)

# NDRE 
# Normalized Difference Red Edge Index
# (NIR - RE)/(NIR + RE)

ndre <- (nir - edge)/(nir + edge)
writeRaster(ndre, filename = paste0("VIs/", datum, "_ndre.tif"), format="GTiff", overwrite=TRUE)

# ExG 
# Excess Green Index
# 2*g-r-b

exg <- 2*(green/(green+blue+red))-(red/(green+blue+red))-(blue/(green+blue+red))
writeRaster(exg, filename = paste0("VIs/", datum, "_exg.tif"), format="GTiff", overwrite=TRUE)

# ExR 
# Excess Red Index
# 1.4(g-b)

exr <- 1.4*((green/(green+blue+red))-(blue/(green+blue+red)))
writeRaster(exr, filename = paste0("VIs/", datum, "_exr.tif"), format="GTiff", overwrite=TRUE)

# ExG-ExR 
# difference of ExG and ExR

exgexr <- exg-exr
writeRaster(exgexr, filename = paste0("VIs/", datum, "_exgexr.tif"), format="GTiff", overwrite=TRUE)

# REMOVE MS BANDS ---------------------------------------------------------

rm("blue", "edge", "green", "nir", "red")


# ZONE SHP LOAD -----------------------------------------------------------

zone <- if(pole == "male"){
  st_read("shp/male.shp")       # n = 1920
} else if (pole == "velke"){
  st_read("shp/velke.shp")      # n = 6137
} else{
  print("ZEROOO")
}

# plot(zone)

zone <- zone %>% 
  select(2,3)


# ZONAL STATISTICS --------------------------------------------------------

# stacked <- stack(sr, ndvi, gndvi, savi, evi, mcari, ndre, exg, exr, exgexr)
# ex <- raster::extract(stacked, zone, fun ='mean', na.rm=TRUE, df=TRUE, weights = TRUE) # https://stackoverflow.com/questions/49936943/calculate-zonal-statistics-in-r-as-in-gis
# TAKES HOURS !!!


# ZONAL STATISTICS WO STACK -----------------------------------------------


zs_sr <- raster::extract(sr, zone, fun ='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
names(zs_sr)[2] <- 'sr'
zs_sr$ID <- as.numeric(zs_sr$ID)

zs_ndvi <- raster::extract(ndvi, zone, fun ='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
names(zs_ndvi)[2] <- 'ndvi'
zs_ndvi$ID <- as.numeric(zs_ndvi$ID)

zs_gndvi <- raster::extract(gndvi, zone, fun ='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
names(zs_gndvi)[2] <- 'gndvi'
zs_gndvi$ID <- as.numeric(zs_gndvi$ID)

zs_savi <- raster::extract(savi, zone, fun ='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
names(zs_savi)[2] <- 'savi'
zs_savi$ID <- as.numeric(zs_savi$ID)

zs_evi <- raster::extract(evi, zone, fun ='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
names(zs_evi)[2] <- 'evi'
zs_evi$ID <- as.numeric(zs_evi$ID)

zs_mcari <- raster::extract(mcari, zone, fun ='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
names(zs_mcari)[2] <- 'mcari'
zs_mcari$ID <- as.numeric(zs_mcari$ID)

zs_ndre <- raster::extract(ndre, zone, fun ='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
names(zs_ndre)[2] <- 'ndre'
zs_ndre$ID <- as.numeric(zs_ndre$ID)

zs_exg <- raster::extract(exg, zone, fun ='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
names(zs_exg)[2] <- 'exg'
zs_exg$ID <- as.numeric(zs_exg$ID)

zs_exr <- raster::extract(exr, zone, fun ='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
names(zs_exr)[2] <- 'exr'
zs_exr$ID <- as.numeric(zs_exr$ID)

zs_exgexr <- raster::extract(exgexr, zone, fun ='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
names(zs_exgexr)[2] <- 'exgexr'
zs_exgexr$ID <- as.numeric(zs_exgexr$ID)

# REMOVE VIs AND CHECK NAs--------------------------------------------------


apply(zs_sr, 2, function(x) any(is.na(x)))
apply(zs_ndvi, 2, function(x) any(is.na(x)))
apply(zs_gndvi, 2, function(x) any(is.na(x)))
apply(zs_savi, 2, function(x) any(is.na(x)))
apply(zs_evi, 2, function(x) any(is.na(x)))
apply(zs_mcari, 2, function(x) any(is.na(x)))
apply(zs_ndre, 2, function(x) any(is.na(x)))
apply(zs_exg, 2, function(x) any(is.na(x)))
apply(zs_exr, 2, function(x) any(is.na(x)))
apply(zs_exgexr, 2, function(x) any(is.na(x)))


rm("sr", "ndvi", "gndvi", "savi", "evi", "mcari", "ndre", "exg", "exr", "exgexr")


# FINAL TABLE -------------------------------------------------------------


require(plyr)
tab_fin  <- join_all(list(zs_sr, zs_ndvi, zs_gndvi, 
                    zs_savi, zs_evi, zs_mcari, zs_ndre,
                    zs_exg, zs_exr, zs_exgexr), by = 'ID', type = 'full')


# *** WR *** WR *** WR ** WR ***
write.csv(tab_fin, file = paste0("final_tabs/", datum, ".csv"), row.names=FALSE)
# *** WR *** WR *** WR ** WR ***                
