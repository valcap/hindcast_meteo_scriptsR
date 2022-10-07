# Purpose

# Load libraries
library(raster)
library(rgdal)
library(elevatr)
library(tidyverse)

# Global variables
my_working_dir  <- paste0("/home/capecchi/art_hindcast_meteo/DATA")
my_figs_dir     <- paste0("/home/capecchi/art_hindcast_meteo/figs")
my_output_dir   <- paste0("/home/capecchi/art_hindcast_meteo/res")
my_shp_dir      <- paste0("/home/capecchi/Sync/shp")

# Working directory
setwd(my_working_dir)

# Load the point locations
my_anag_file     <- paste0(my_working_dir,"/points050_arcis_extent.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))
world_adm0       <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/ne_10m_coastline")),layer="ne_10m_coastline")
ita_adm0         <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/ITA_adm")),layer="ITA_adm0")
extent_ita       <- ita_adm0[ita_adm0$NAME_ENGLI=="Italy",]
extent_eobs      <- extent(pointCoordinates)
elevation        <- get_elev_raster(extent_ita, z = 8)
elevation_crop   <- crop(elevation,extent_eobs)

#####
setEPS()
postscript(paste0(my_figs_dir,"/virtual_rain-gauges_ARCIS_loc.eps"))
plot(pointCoordinates, main=paste0("Virtual ARCIS weather stations\n","Total number=",length(pointCoordinates)), add=FALSE)
my_col <- terrain.colors(15)
plot(elevation_crop, add=TRUE, breaks = c(-2,5,50,100,200,300,400,500,750,1000,1500,2000,3000), col=my_col, legend=FALSE)
plot(pointCoordinates, main=paste0("Virtual ARCIS weather stations\n","Total number=",length(pointCoordinates)), add=TRUE)
plot(world_adm0, border="black", lwd=1, main="", add=TRUE)
dev.off()

#####
eobs_topo  <- as.tibble(raster::extract(elevation_crop, pointCoordinates))[,1]
eobs_topo2 <- eobs_topo %>% filter(between(value,0,4000))

setEPS()
postscript(paste0(my_figs_dir,"/virtual_rain-gauges_ARCIS_dist.eps"))
hist(eobs_topo2$value,breaks=15, xlab="Altitude[m]", ylab="Count", main="")
dev.off()

