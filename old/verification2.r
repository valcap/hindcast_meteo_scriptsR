# Purpose
print(paste0("RUNNING CONTINGENCY TABLE AND RELATED SCORES"))

# Load libraries
library(raster)
library(rgdal)
library(ncdf4)
library(verification)
library(Metrics)
library(plotrix)

# Global variables
my_working_dir  <- paste0("~/art_hindcast_meteo/DATA")
my_figs_dir     <- paste0("~/art_hindcast_meteo/figs")
my_output_dir   <- paste0("~/art_hindcast_meteo/res")
my_shp_dir      <- paste0("~/Sync/shp")
#my_anag_file <- paste0("ARCIS/Arcis_Stazioni_Pluviometriche.csv")
my_anag_file    <- paste0("NOAA/stations_synop_BIGGER-DOMAIN.csv")
my_figs_dir     <- paste0("../figs")
my_eobs_long    <- NULL
my_era5_long    <- NULL
my_molo_long    <- NULL
my_bola_long    <- NULL
vec_thres       <- c(2.5,5,10,20,30,40,50,60,70,80,90,100)
vec_thres       <- 50
TD_normal       <- FALSE
bool_plot_daily <- FALSE
max_plot_daily  <- 200
vec_years       <- c(1982,2015)
vec_eobs_pod    <- NULL
vec_era5_pod    <- NULL
vec_molo_pod    <- NULL
vec_bola_pod    <- NULL
vec_eobs_far    <- NULL
vec_era5_far    <- NULL
vec_molo_far    <- NULL
vec_bola_far    <- NULL

# Working directory
setwd(my_working_dir)

# Load the point locations
anag_full <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))
print(paste0(""))
# Load shapefiles && plot point locations
ITA_adm0 <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/ITA_adm")),layer="ITA_adm0")
ITA_adm1 <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/ITA_adm")),layer="ITA_adm1")
setEPS()
postscript(paste0(my_figs_dir,"/rain-gauges.eps"))
plot(pointCoordinates, main=paste0("# of rain-gauges ",length(pointCoordinates)))
#plot(ITA_adm1, border="red", lwd=1, main="", add=TRUE)
plot(ITA_adm0, border="black", lwd=1, main="", add=TRUE)
plot(extent(pointCoordinates), add=TRUE, col='red', lwd=1)
dev.off()

###########################
# MAIN LOOP OVER THRESHOLDS
for (thres in vec_thres) {

#################
# LOOP OVER YEARS
for (j in seq(1,length(vec_years))) {
  my_year <- vec_years[j]
  print(paste0("___Running year ",my_year,"___"))

###
  print(paste0("+++Reading and resampling files..."))
  # Raster files definition
  file_arcis  <- paste0("ARCIS/tp/ARCIS_",my_year,".nc")
  file_eobs   <- paste0("E-OBS/tp/E-OBS_day_tp_",my_year,".nc")
  file_era5   <- paste0("ERA5-Land/tp/ERA5-Land_",my_year,"-daysum_tp.nc")
  file_moloch <- paste0("MOLOCH/moloch_",my_year,"_daysum_apcp.nc")
  file_bolam  <- paste0("BOLAM/bolam_",my_year,"_daysum_apcp.nc") 
  # Allocate rasters
  arcis_orig  <- brick(file_arcis)
  eobs_orig   <- brick(file_eobs)
  era5_orig   <- brick(file_era5)
  moloch_orig <- brick(file_moloch)
  bolam_orig  <- brick(file_bolam)
###

###
  # Get extent on the basis of the point coordinates
  my_extent       <- extent(pointCoordinates)
  # Crop E-OBS grid on my_extent and regrid (if needed)
  eobs_cropped    <- crop(eobs_orig, my_extent)
#  eobs_resam      <- resample(eobs_reduced, eobs_orig, method="ngb")
  eobs_resam      <- eobs_cropped
  # Crop ERA5-Land grid on my_extent and regrid (if needed)
  era5_cropped    <- crop(era5_orig, my_extent)
  era5_resam      <- resample(era5_cropped, eobs_orig, method="bilinear")
  # Crop MOLOCH grid on my_extent and regrid (if needed)
  moloch_cropped  <- crop(moloch_orig, my_extent)
  moloch_resam    <- resample(moloch_cropped, eobs_orig, method="bilinear")
  # Crop BOLAM grid on my_extent and regrid (if needed)
  bolam_cropped   <- crop(bolam_orig, my_extent)
  bolam_resam     <- resample(bolam_cropped, eobs_orig, method="bilinear")
  print(paste0("+++OK"))
###

###
  # seleziono solo alcuni layers sulla finestra temporale comune
  arcis_new        <- subset(arcis_orig, 2:365)
  eobs_resam_new   <- subset(eobs_resam, 2:365)
  era5_resam_new   <- subset(era5_resam, 2:365)
  moloch_resam_new <- subset(moloch_resam, 1:364)
  bolam_resam_new  <- subset(bolam_resam, 1:364)
  my_raster_names         <- seq(1,364)
  names(arcis_new)        <- my_raster_names
  names(eobs_resam_new)   <- my_raster_names
  names(era5_resam_new)   <- my_raster_names
  names(moloch_resam_new) <- my_raster_names
  names(bolam_resam_new)  <- my_raster_names
###

###
  if ( max_plot_daily == TRUE ) {
    print(paste0("+++Plotting data..."))
    # Plot raw data
    for (i in seq(1,max_num_plot)) {
      jday=i+1
      print(paste0("___Plot tp ",my_year," and jday",jday,"___"))
      png(paste0("../figs/tp_",my_year,"_jday",jday,".png"))
      par(mfrow=c(2,2))
#      plot(arcis_new[[i]], main="ARCIS", zlim=c(0,50))
#      plot(ITA_adm0, border="black", lwd=1, add=TRUE)
      plot(eobs_resam_new[[i]], main="E-OBS", zlim=c(0,50))
      plot(ITA_adm0, border="black", lwd=1, add=TRUE)
      plot(era5_resam_new[[i]]*100, main="ERA5-Land", zlim=c(0,50))
      plot(ITA_adm0, border="black", lwd=1, add=TRUE)
      plot(moloch_resam_new[[i]], main="MOLOCH", zlim=c(0,50))
      plot(ITA_adm0, border="black", lwd=1, add=TRUE)
      plot(bolam_resam_new[[i]], main="BOLAM", zlim=c(0,50))
      plot(ITA_adm0, border="black", lwd=1, add=TRUE)
      dev.off()
    }
    print(paste0("+++OK"))
  }
###
  
###
  # Extracting data
  print(paste0("+++Extracting data..."))
  # soluzione più veloce (forse)
  # estrarre tutti i punti e per tutti i layers (cioè per tutti i time step)
  arci_values <- as.data.frame(extract(arcis_new, pointCoordinates))
  eobs_values <- as.data.frame(extract(eobs_resam_new, pointCoordinates))
  era5_values <- as.data.frame(extract(era5_resam_new, pointCoordinates))
  molo_values <- as.data.frame(extract(moloch_resam_new, pointCoordinates))
  bola_values <- as.data.frame(extract(bolam_resam_new, pointCoordinates))
  my_arci <- NULL
  my_eobs <- NULL
  my_era5 <- NULL
  my_molo <- NULL
  my_bola <- NULL
  # Da data frame a vector
  for (i in 1:364) {
    my_arci <- c(my_arci,arci_values[,i]) 
    my_eobs <- c(my_eobs,eobs_values[,i]) 
    my_era5 <- c(my_era5,era5_values[,i])
    my_molo <- c(my_molo,molo_values[,i])
    my_bola <- c(my_bola,bola_values[,i])
  }
  print(paste0("+++OK"))

###
  print(paste0("+++Creating contingency table..."))
  # Create binary vector for E-OBS
  my_eobs_thres <- replace(my_eobs, my_eobs <= thres,0)
  my_eobs_thres <- replace(my_eobs_thres, my_eobs > thres,1)
  # Create binary vector for ERA5-Land
  my_era5_thres <- replace(my_era5*100, my_era5*100 <= thres,0)
  my_era5_thres <- replace(my_era5_thres,    my_era5*100 > thres,1)
  # Create binary vector for MOLOCH
  my_molo_thres <- replace(my_molo, my_molo <= thres,0)
  my_molo_thres <- replace(my_molo_thres, my_molo > thres,1)
  # Create binary vector for BOLAM
  my_bola_thres <- replace(my_bola, my_bola <= thres,0)
  my_bola_thres <- replace(my_bola_thres, my_bola > thres,1)
  # Plot performance diagram
  my_tab_era5 <- table.stats(my_eobs_thres, my_era5_thres)
  my_tab_molo <- table.stats(my_eobs_thres, my_molo_thres)
  my_tab_bola <- table.stats(my_eobs_thres, my_bola_thres)
  vec_era5_pod <- c(vec_era5_pod,my_tab_era5$POD)
  vec_molo_pod <- c(vec_molo_pod,my_tab_molo$POD)
  vec_bola_pod <- c(vec_bola_pod,my_tab_bola$POD)
  vec_era5_far <- c(vec_era5_far,my_tab_era5$FAR)
  vec_molo_far <- c(vec_molo_far,my_tab_molo$FAR)
  vec_bola_far <- c(vec_bola_far,my_tab_bola$FAR)
  print(paste0("+++OK"))
###  

###  
  # Empty cache
  rm(my_eobs_thres, my_eobs, eobs_values, eobs_resam_new, eobs_resam, eobs_orig, eobs_cropped)
  rm(my_era5_thres, my_era5, era5_values, era5_resam_new, era5_resam, era5_orig, era5_cropped)
  rm(my_molo_thres, my_molo, molo_values, moloch_resam_new, moloch_resam, moloch_orig, moloch_cropped)
  rm(my_bola_thres, my_bola, bola_values, bolam_resam_new, bolam_resam, bolam_orig, bolam_cropped)
  my_eobs_thres<-NULL; my_eobs<-NULL; eobs_values<- NULL
  eobs_cropped<-NULL; eobs_orig<-NULL; eobs_resam<-NULL; eobs_resam_new<-NULL
  my_era5_thres<-NULL; my_era5<-NULL; era5_values<-NULL; my_tab_era5<-NULL
  era5_cropped<-NULL; era5_orig<-NULL; era5_resam<-NULL; era5_resam_new<-NULL
  my_molo_thres<-NULL; my_molo<-NULL; molo_values<-NULL; my_tab_molo<-NULL
  molo_cropped<-NULL; molo_orig<-NULL; molo_resam<-NULL; molo_resam_new<-NULL
  my_bola_thres<-NULL; my_bola<-NULL; bola_values<-NULL; my_tab_bola<-NULL
  bola_cropped<-NULL; bola_orig<-NULL; bola_resam<-NULL; bola_resam_new<-NULL
###  
  
  # End of my_year loop
  print(paste0("___End of ",my_year,"___"))
}

# Plot final performance diagram
setEPS()
postscript(paste0(my_figs_dir,"/Perfor_Diag_",thres,".eps"))
performance.diagram(main = paste0("Precipitation ",thres," mm/24-hour"))
points(1 - mean(vec_era5_far), mean(vec_era5_pod), col = "red", cex = 2.0, pch=19)
points(1 - mean(vec_molo_far), mean(vec_molo_pod), col = "orange", cex = 2.0, pch=19)
points(1 - mean(vec_bola_far), mean(vec_bola_pod), col = "blue", cex = 2.0, pch=19)
# Draw the lengend on the plot and save it
legend("topleft",legend=c("ERA5-Land","MOLOCH","BOLAM"),pch=19,col=c("red","orange","blue"))
dev.off()

# End of thres loop
}

# Byebye
print(paste0("Byebye"))
#quit()

