# Purpose
print(paste0(""))

# Load libraries
library(raster)
library(ncdf4)
library(rgdal)
library(graphics)
library(elevatr)
library(tidyverse)

# Global variables
u               <- 'capecchi'
my_working_dir  <- paste0("/home/",u,"/art_hindcast_meteo/")
my_figs_dir     <- paste0("/home/",u,"/art_hindcast_meteo/figs")
my_output_dir   <- paste0("/home/",u,"/art_hindcast_meteo/res")
my_shp_dir      <- paste0("/home/",u,"/Sync/shp")

# Shapefiles
my_italy <- readOGR(
  dsn= paste0(my_shp_dir,"/") ,
  layer="ITA_adm0",
  verbose=FALSE
)
my_coastline <- readOGR( 
  dsn= paste0(my_shp_dir,"/ne_10m_coastline/") , 
  layer="ne_10m_coastline",
  verbose=FALSE
)
my_coastline_lr <- readOGR(
  dsn= paste0(my_shp_dir,"/gshhg-shp-2.3.7/GSHHS_shp/l/") ,
  layer="GSHHS_l_L1",
  verbose=FALSE
)

# ALPS-OVEST (NW in Brunetti et al 2001)
lon1ALPO=8.00
lon2ALPO=11.0
lat1ALPO=45.7
lat2ALPO=46.5
# ALPS-EST (NEN in Brunetti et al 2001)
lon1ALPE=12.5
lon2ALPE=13.5
lat1ALPE=46.0
lat2ALPE=46.7
# PIANURA PADANA (NES in Brunetti et al 2001)
lon1PAD=8.00
lon1PAD=10.00
lon2PAD=12.5
lat1PAD=44.7
lat2PAD=45.5
# TIRRENO (CE in Brunetti et al 2001)
lon1TIR=10.5
lon2TIR=11.8
lat1TIR=42.5
lat2TIR=43.8
# SARDEGNA (IS in Brunetti et al 2001)
lon1SAR=8.5
lon2SAR=9.5
lat1SAR=39.0
lat2SAR=41.0
# SUD (SO in Brunetti et al 2001)
lon1SUD=14.0
lon2SUD=16.0
lat1SUD=40.8
lat2SUD=41.8
# MIL (as in Reder et al 2022)
lon1MIL=8.75
lon2MIL=9.75
lat1MIL=45.0
lat2MIL=45.75

# Working directory
setwd(my_working_dir)

# Get the extent of the GRIPHO dataset
file_obse        <- paste0("DATA/GRIPHO/tp/gripho-v2-3km_year_tp_AVG2001-2016.nc")
obse_orig        <- brick(file_obse)
my_extent_gripho <- extent(obse_orig)

# Loading observational dataset
my_anag_file     <- paste0(my_working_dir,"/DATA/points025_arcis_extent.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates_arcis <- anag_full
coordinates(pointCoordinates_arcis) = ~ LONGITUDE+LATITUDE
my_anag_file     <- paste0(my_working_dir,"/DATA/points025_moloch_extent.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates_gripho <- anag_full
coordinates(pointCoordinates_gripho) = ~ LONGITUDE+LATITUDE
#print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates_arcis)))

# ...and plot them on a map
setEPS()
postscript(paste0(my_figs_dir,"/virtual_rain-gauges.eps"))
#plot(extent(my_extent_gripho), main=paste0("Crosses and closed circles: GRIPHO/E-OBS virtual rain-gauges (", length(pointCoordinates_gripho),")\n","Crosses: ARCIS virtual rain-gauges (", length(pointCoordinates_arcis),")"), add=F, col='white', lwd=1, xlab="Longitude [degree]", ylab="Latitude [degree]")
plot(extent(my_extent_gripho), main="", add=F, col='white', lwd=1, xlab="Longitude [degree]", ylab="Latitude [degree]")
plot(pointCoordinates_gripho, pch=19, cex=0.3, add=T, col="black")
plot(pointCoordinates_arcis, pch=4, cex=0.7, add=T ,col="black")
plot(my_coastline_lr, lwd=1.5, border="black", add=TRUE)
rect(lon1SUD,lat1SUD,lon2SUD,lat2SUD,col=NA,border="red")
text(lon2SUD+1.,lat2SUD+0.2, "South Italy",col="red",cex=1.1,)
rect(lon1TIR,lat1TIR,lon2TIR,lat2TIR,col=NA,border="red")
text(lon1TIR+0.5,lat1TIR-0.6, paste0("Tyrrhenian\nCoast"),col="red",cex=1.1,)
rect(lon1ALPO,lat1ALPO,lon2ALPO,lat2ALPO,col=NA,border="red")
text(lon1ALPO+1.,lat2ALPO+0.2, "Western Alps",col="red",cex=1.1,)
rect(lon1ALPE,lat1ALPE,lon2ALPE,lat2ALPE,col=NA,border="red")
text(lon2ALPE+1.4,(lat2ALPE+lat1ALPE)/2, "Eastern Alps",col="red",cex=1.1,)
rect(lon1PAD,lat1PAD,lon2PAD,lat2PAD,col=NA,border="red")
text(lon2PAD+1.,(lat2PAD+lat1PAD)/2, "Po valley",col="red",cex=1.1,)
rect(lon1SAR,lat1SAR,lon2SAR,lat2SAR,col=NA,border="red")
text(lon2SAR+1.1,(lat2SAR+lat1SAR)/2, "Sardinia",col="red",cex=1.1,)
dev.off()

# Histograms of elevations
ita_adm0         <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/ITA_adm")),layer="ITA_adm0")
extent_ita       <- ita_adm0[ita_adm0$NAME_ENGLI=="Italy",]
elevation        <- get_elev_raster(extent_ita, z = 8)
elevation_crop   <- crop(elevation,extent_ita)
#####
gripho_topo  <- as.tibble(raster::extract(elevation_crop, pointCoordinates_gripho))[,1]
gripho_topo2 <- gripho_topo %>% filter(between(value,0,4000))
arcis_topo  <- as.tibble(raster::extract(elevation_crop, pointCoordinates_arcis))[,1]
arcis_topo2 <- arcis_topo %>% filter(between(value,0,4000))

actual_gripho_topo <- as.tibble(raster::extract(elevation_crop,extent(pointCoordinates_gripho)))[,1]
actual_gripho_topo2 <- actual_gripho_topo %>% filter(between(value,0,4000))
actual_arcis_topo <- as.tibble(raster::extract(elevation_crop,extent(pointCoordinates_arcis)))[,1]
actual_arcis_topo2 <- actual_arcis_topo %>% filter(between(value,0,4000))

# Plot histograms of GRIPHO and ARCIS virtual rain-gauges
# plus GRIPHO ad ARCIS domain actual data
setEPS()
postscript(paste0(my_figs_dir,"/virtual_rain-gauges_topo.eps"))
par(mfrow=c(2,2))
hist(gripho_topo2$value,breaks=15, xlab="Altitude[m]", ylab="Count", main="Topography of\nGRIPHO/E-OBS virtual rain-gauges")
hist(arcis_topo2$value,breaks=15, xlab="Altitude[m]", ylab="Count", main="Topography of\nARCIS virtual rain-gauges")
hist(actual_gripho_topo2$value,breaks=15, xlim=c(0,2500), xlab="Altitude[m]", ylab="Count", main="Topography of\nGRIPHO/E-OBS domaim")
hist(actual_arcis_topo2$value,breaks=15, xlim=c(0,2500), xlab="Altitude[m]", ylab="Count", main="Topography of\nARCIS domaim")
dev.off()

# Byebye
print(paste0("Byebye"))

# Get average orography out of the six areas
summary(raster::extract(elevation,extent(lon1ALPO,lon2ALPO,lat1ALPO,lat2ALPO)))[4]
summary(raster::extract(elevation,extent(lon1ALPE,lon2ALPE,lat1ALPE,lat2ALPE)))[4]
summary(raster::extract(elevation,extent(lon1PAD,lon2PAD,lat1PAD,lat2PAD)))[4]
summary(raster::extract(elevation,extent(lon1TIR,lon2TIR,lat1TIR,lat2TIR)))[4]
summary(raster::extract(elevation,extent(lon1SUD,lon2SUD,lat1SUD,lat2SUD)))[4]
summary(raster::extract(elevation,extent(lon1SAR,lon2SAR,lat1SAR,lat2SAR)))[4]

# quit()
