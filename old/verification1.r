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
bool_plot_daily <- FALSE
my_working_dir  <- paste0("/home/capecchi/art_hindcast_meteo/DATA")
my_figs_dir     <- paste0("/home/capecchi/art_hindcast_meteo/figs")
my_output_dir   <- paste0("/home/capecchi/art_hindcast_meteo/res")
my_shp_dir      <- paste0("/home/capecchi/Sync/shp")
#my_anag_file <- paste0("ARCIS/Arcis_Stazioni_Pluviometriche.csv")
my_anag_file    <- paste0("NOAA/stations_synop_BIGGER-DOMAIN.csv")
my_figs_dir     <- paste0("../figs")
my_eobs_long    <- NULL
my_era5_long    <- NULL
my_molo_long    <- NULL
my_bola_long    <- NULL
vec_thres       <- c(1,2,5,10,20,30,40,50,60,70,80,90,100)
vec_thres       <- c(1,2,5,10)
TD_normal       <- FALSE

#
setEPS()
postscript(paste0(my_figs_dir,"/Taylor_Diag.eps"))

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

# MAIN LOOP OVER YEARS
for (my_year in c(1982,2015)) {
#for (my_year in c(1982)) {
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
  if ( bool_plot_daily == TRUE ) {
    # Plot raw data
    for (i in seq(1,100)) {
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
###

###  
  # Store data for final Taylor diagram
  my_eobs_long <- c(my_eobs_long,my_eobs)
  my_era5_long <- c(my_era5_long,my_era5)
  my_molo_long <- c(my_molo_long,my_molo)
  my_bola_long <- c(my_bola_long,my_bola)  
###

###
  # Plot points in the Taylor diagram
  print(paste0("+++Plotting data..."))
  if ( (my_year == 1981) || (my_year == 1982) ) {
    taylor.diagram(my_eobs,my_era5*100,add=F,col="red",ref.sd=TRUE,normalize=TD_normal,main="",grad.corr.lines=c(0.4,0.6,0.8))
  } else {
    taylor.diagram(my_eobs,my_era5*100,add=T,col="red",ref.sd=TRUE,normalize=TD_normal)
  }
  taylor.diagram(my_eobs,my_molo,add=T,col="orange",ref.sd=TRUE,normalize=TD_normal)
  taylor.diagram(my_eobs,my_bola,add=T,col="blue",ref.sd=TRUE,normalize=TD_normal)
###
  
  # End of my_year loop
  print(paste0("___End of ",my_year,"___"))
}
# Draw the lengend on the plot and save it
legend("topright",legend=c("ERA5-Land","MOLOCH","BOLAM"),pch=19,col=c("red","orange","blue"))
dev.off()
 
# Standard statistical verification scores for continuos variables
print(paste0("+++Calculating standard statistical verification scores..."))
# Correlation
cor_era5 <- cor(my_era5_long*100,my_eobs_long, use = "complete.obs")
cor_molo <- cor(my_molo_long,my_eobs_long, use = "complete.obs")
cor_bola <- cor(my_bola_long,my_eobs_long, use = "complete.obs")
# RMSE
rmse_era5 <- sqrt( mean( (my_era5_long*100-my_eobs_long)^2 , na.rm = TRUE ) ) 
rmse_molo <- sqrt( mean( (my_molo_long-my_eobs_long)^2 , na.rm = TRUE ) ) 
rmse_bola <- sqrt( mean( (my_bola_long-my_eobs_long)^2 , na.rm = TRUE ) ) 
# (multiplicative) bias
bias_era5 <- (mean(my_era5_long*100, na.rm = TRUE))/(mean(my_eobs_long, na.rm = TRUE))
bias_molo <- (mean(my_molo_long, na.rm = TRUE))/(mean(my_eobs_long, na.rm = TRUE))
bias_bola <- (mean(my_bola_long, na.rm = TRUE))/(mean(my_eobs_long, na.rm = TRUE))
# ME, MSE, MAE
verify_era5 <- verify(my_eobs_long,my_era5_long*100, frcst.type = "cont", obs.type = "cont",)
verify_molo <- verify(my_eobs_long,my_molo_long, frcst.type = "cont", obs.type = "cont",)
verify_bola <- verify(my_eobs_long,my_bola_long, frcst.type = "cont", obs.type = "cont",)
# Write to output file
df <- data.frame(c(cor_era5,rmse_era5,bias_era5,verify_era5$ME,verify_era5$MAE), c(cor_molo,rmse_molo,bias_molo,verify_molo$ME,verify_molo$MAE), c(cor_bola,rmse_bola,bias_bola,verify_bola$ME,verify_bola$MAE))
write.table(df, paste0(my_output_dir,"/verification_scores.csv"), col.names = c("ERA5-Land","MOLOCH","BOLAM"), row.names = c("Correlation","RMSE","Bias","Mean Error","Mean Absolute Error"), sep=";")
print(paste0("+++OK"))

# Set output EPS file
print(paste0("+++Plotting the summary plot..."))
setEPS()
postscript(paste0(my_figs_dir,"/Taylor_Diag_final.eps"))
taylor.diagram(my_eobs_long,my_era5_long*100,add=F, col="red",   ref.sd=TRUE, normalize=TD_normal, main="", pcex=1.5, grad.corr.lines=c(0.4,0.5,0.6,0.7,0.8))
taylor.diagram(my_eobs_long,my_molo_long,    add=T, col="orange",ref.sd=TRUE, normalize=TD_normal, pcex=1.5)
taylor.diagram(my_eobs_long,my_bola_long,    add=T, col="blue",  ref.sd=TRUE, normalize=TD_normal, pcex=1.5)
lpos<-sd(my_eobs_long)
legend("topright",legend=c("ERA5-Land","MOLOCH","BOLAM"),pch=19,col=c("red","orange","blue"), cex=1)
dev.off()
print(paste0("Ok"))

# Byebye
print(paste0("Byebye"))
quit()

