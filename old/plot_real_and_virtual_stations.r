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
# Extent of the MOLOCH model
#minlon          <- 2.4
minlon          <- 4.4
maxlon          <- 19.875
minlat          <- 34.2123
#maxlat          <- 49.6498
maxlat          <- 48.6498

# Working directory
setwd(my_working_dir)

my_extent_moloch <- as(extent(minlon,maxlon,minlat,maxlat), "SpatialPolygons")
# Load the NOAA locations
my_anag_file     <- paste0(my_working_dir,"/NOAA/stations_synop_BIGGER-DOMAIN_LOCALE.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates_noaa <- anag_full
coordinates(pointCoordinates_noaa)= ~ LONGITUDE+LATITUDE
world_adm0       <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/ne_10m_coastline")),layer="ne_10m_coastline")
ita_adm0         <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/ITA_adm")),layer="ITA_adm0")
coasts_lowres    <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/gshhg-shp-2.3.7/GSHHS_shp/c")),layer="GSHHS_c_L1")
coasts_medres    <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/gshhg-shp-2.3.7/GSHHS_shp/l")),layer="GSHHS_l_L1")
coasts_higres    <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/gshhg-shp-2.3.7/GSHHS_shp/i")),layer="GSHHS_i_L1")
extent_ita       <- ita_adm0[ita_adm0$NAME_ENGLI=="Italy",]
extent_noaa      <- extent(pointCoordinates_noaa)
elevation        <- get_elev_raster(extent_ita, z = 4)
elevation_moloch <- crop(elevation,my_extent_moloch)
# Load the E-OBS locations
my_anag_file     <- paste0(my_working_dir,"/E-OBS/points050_moloch_extent.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates_eobs <- anag_full
coordinates(pointCoordinates_eobs)= ~ LONGITUDE+LATITUDE
#world_adm0       <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/ne_10m_coastline")),layer="ne_10m_coastline")
#ita_adm0         <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/ITA_adm")),layer="ITA_adm0")
#extent_ita       <- ita_adm0[ita_adm0$NAME_ENGLI=="Italy",]
#extent_eobs      <- extent(pointCoordinates_eobs)
#elevation        <- get_elev_raster(extent_ita, z = 4)

#####
# IL SALVATAGGIO E' FATTO CON Rstudio PER AVERE UN'IMMAGINE BEN FORMATTATA
setEPS()
postscript(paste0(my_figs_dir,"/weather_stations.eps"))
par(mfrow=c(2,2))

my_title=paste0("(a)\n Virtual weather stations\n","Total number=",length(pointCoordinates_eobs))
#my_col <- colorRampPalette(c("white","black"))
#my_breaks <- c(5,25,50,75,100,125,200,250,300,350,400,450,500,550,600,650,700,750,800,900,1000,1200,1500,1700,2000,2500,3000)
my_breaks <- c(5,50,500,750,1000,1500,2000,3000)
plot(pointCoordinates_eobs, main=my_title, add=FALSE, pch=20, cex=0.5, xlim=c(minlon,maxlon), ylim=c(minlat,maxlat))
plot(elevation_moloch, add=TRUE, breaks = my_breaks, col=terrain.colors(10), xlim=c(minlon,maxlon), ylim=c(minlat,maxlat), legend=TRUE, legend.args = list(text = 'Elevation (m)', side = 1, font = 2, line = 2.5, cex = 0.8))
plot(pointCoordinates_eobs, main=my_title, add=TRUE, pch=20, cex=0.5, xlim=c(minlon,maxlon), ylim=c(minlat,maxlat))
plot(coasts_lowres, border="black", lwd=1, main="", add=TRUE, xlim=c(minlon,maxlon), ylim=c(minlat,maxlat))

my_title=paste0("(b)\n Wind weather stations\n","Total number=",length(pointCoordinates_noaa))
plot(pointCoordinates_noaa, main=my_title, add=FALSE, pch=20, cex=0.5, xlim=c(minlon,maxlon), ylim=c(minlat,maxlat))
plot(elevation_moloch, add=TRUE, breaks = my_breaks, col=terrain.colors(10), xlim=c(minlon,maxlon), ylim=c(minlat,maxlat), legend=FALSE)
plot(pointCoordinates_noaa, main=my_title, add=TRUE, pch=20, cex=0.5, xlim=c(minlon,maxlon), ylim=c(minlat,maxlat))
plot(coasts_lowres, border="black", lwd=1, main="", add=TRUE, xlim=c(minlon,maxlon), ylim=c(minlat,maxlat))
#dev.off()

######################
#elevation_moloch[elevation_moloch<0] <- NA
#summary(elevation_moloch)
#hist(elevation_moloch)
######################


#####

noaa_topo  <- as_tibble(raster::extract(elevation_moloch, pointCoordinates_noaa))[,1]
noaa_topo2 <- noaa_topo %>% filter(between(value,0,2500))
eobs_topo  <- as_tibble(raster::extract(elevation_moloch, pointCoordinates_eobs))[,1]
eobs_topo2 <- eobs_topo %>% filter(between(value,0,2500))

# IL SALVATAGGIO E' FATTO CON Rstudio PER AVERE UN'IMMAGINE BEN FORMATTATA
#setEPS()
#postscript(paste0(my_figs_dir,"/stations_hist.eps"))
#par(mfrow=c(1,2))
hist(eobs_topo2$value,breaks=15, xlab="Elevation[m]", ylab="Count", main="(c)\n Virtual Weather Stations\nElevation distribution")
hist(noaa_topo2$value,breaks=15, xlab="Elevation[m]", ylab="Count", main="(d)\n Wind Weather Stations\nElevation distribution")
dev.off()

