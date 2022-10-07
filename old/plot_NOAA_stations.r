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
my_anag_file     <- paste0(my_working_dir,"/NOAA/stations_synop_BIGGER-DOMAIN_LOCALE.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))
world_adm0       <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/ne_10m_coastline")),layer="ne_10m_coastline")
ita_adm0         <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/ITA_adm")),layer="ITA_adm0")
extent_ita       <- ita_adm0[ita_adm0$NAME_ENGLI=="Italy",]
extent_noaa      <- extent(pointCoordinates)
elevation        <- get_elev_raster(extent_ita, z = 8)
elevation_crop   <- crop(elevation,extent_noaa)

#####
setEPS()
postscript(paste0(my_figs_dir,"/noaa_weather_stations.eps"))
my_title=paste0("Synop weather stations\n","Total number=",length(pointCoordinates))
plot(pointCoordinates, main=my_title, add=FALSE)
my_col <- terrain.colors(15)
plot(elevation_crop, add=TRUE, breaks = c(5,50,100,200,300,400,500,750,1000,1500,2000,3000), col=my_col)
plot(pointCoordinates, main=my_title, add=TRUE)
plot(world_adm0, border="black", lwd=1, main="", add=TRUE)
dev.off()

#####
noaa_topo  <- as.tibble(raster::extract(elevation_crop, pointCoordinates))[,1]
noaa_topo2 <- noaa_topo %>% filter(between(value,0,4000))

setEPS()
postscript(paste0(my_figs_dir,"/noaa_weather_stations_topo.eps"))
hist(noaa_topo2$value,breaks=15, xlab="Altitude[m]", ylab="Count", main="")
dev.off()

