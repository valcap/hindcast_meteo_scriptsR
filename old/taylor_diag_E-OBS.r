# Purpose
print(paste0("RUNNING TAYLOR DIAGRAMS"))

# Load libraries
library(raster)
library(ncdf4)
library(plotrix)
library(verification)
library(Metrics)

# Working directory
setwd("/OCEANASTORE/progetti/spitbran/")

# Load the point locations
anag_file <- paste0("LaMMA/stations_synop.csv")
anag_full <- read.csv(anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))
print(paste0(""))

# Set output EPS file
setEPS()
postscript(paste0("work/Taylor_Diag_normalized.eps"))

# Initialize variables for final Taylor diagram
my_eobs_long <- NULL 
my_era5_long <- NULL 
my_molo_long <- NULL
my_bola_long <- NULL

# MAIN LOOP OVER YEARS in `seq 1979 2019`
for (my_year in seq(1981,2019,1)) {
#for (my_year in seq(1981,1982,1)) {
  print(paste0("___Running year ",my_year,"___"))

  print(paste0("+++Reading and resampling files..."))

###
  # Raster files definition
  file_eobs   <- paste0("E-OBS/tp/E-OBS_day_tp_",my_year,".nc")
  file_era5   <- paste0("ERA5-Land/tp/ERA5-Land_",my_year,"-daysum_tp.nc")
  file_moloch <- paste0("work/moloch/moloch_",my_year,"_daysum_apcp.nc")
  file_bolam  <- paste0("work/bolam/bolam_",my_year,"_daysum_apcp.nc") 
  
  # Allocate rasters
  eobs_orig   <- stack(file_eobs)
  era5_orig   <- stack(file_era5)
  moloch_orig <- stack(file_moloch)
  bolam_orig  <- stack(file_bolam)
###

###
  # Get MOLOCH extent
  extent_moloch <- extent(moloch_orig)
  # Crop E-OBS grid on the MOLOCH extent
  eobs_reduced  <- crop(eobs_orig,extent_moloch)
  # Crop ERA5 grid on the MOLOCH extent and regrid on the E-OBS grid spacing
  era5_reduced  <- crop(era5_orig,extent_moloch)
  era5_resam    <- resample(era5_reduced, eobs_reduced, method="bilinear")
  # Regrid the MOLOCH grid on the E-OBS grid spacing
  moloch_resam  <- resample(moloch_orig, eobs_reduced, method="bilinear")
  # Crop BOLAM grid on the MOLOCH extent and regrid on the E-OBS grid spacing
  bolam_reduced <- crop(bolam_orig,extent_moloch)
  bolam_resam   <- resample(bolam_reduced, eobs_reduced, method="bilinear")
  print(paste0("+++OK"))
###

###
  # seleziono solo alcuni layers su una finestra temporale comune
  eobs_reduced_new <- subset(eobs_reduced, 2:365)
  era5_reduced_new <- subset(era5_reduced, 2:365)
  moloch_resam_new <- subset(moloch_resam, 1:364)
  bolam_resam_new <- subset(bolam_resam, 1:364)
  my_raster_names     <- seq(1,364)
  names(eobs_reduced_new) <- my_raster_names
  names(era5_reduced_new) <- my_raster_names
  names(moloch_resam_new) <- my_raster_names
  names(bolam_resam_new) <- my_raster_names
###

###
  # Extracting data
  print(paste0("+++Extracting data..."))
  # soluzione più veloce (forse)
  # estrarre tutti i punti e per tutti i layers (cioè per tutti i time step)
  eobs_values <- as.data.frame(extract(eobs_reduced_new, pointCoordinates))
  era5_values <- as.data.frame(extract(era5_reduced_new, pointCoordinates))
  molo_values <- as.data.frame(extract(moloch_resam_new, pointCoordinates))
  bola_values <- as.data.frame(extract(bolam_resam_new, pointCoordinates))
  my_eobs <- NULL
  my_era5 <- NULL
  my_molo <- NULL
  my_bola <- NULL
  for (i in 1:364) {
    my_eobs <- c(my_eobs,eobs_values[,i]) 
    my_era5 <- c(my_era5,era5_values[,i])
    my_molo <- c(my_molo,molo_values[,i])
    my_bola <- c(my_bola,bola_values[,i])
  }
  print(paste0("+++OK"))
###
  
###
  # Plot points in the Taylor diagram
  print(paste0("+++Plotting data..."))
  if ( (my_year == 1981) || (my_year == 1989) ) {
    taylor.diagram(my_eobs,my_era5*100,add=F,col="red",ref.sd=TRUE,normalize=T,main="",grad.corr.lines=c(0.4,0.6,0.8))
  } else {
    taylor.diagram(my_eobs,my_era5*100,add=T,col="red",ref.sd=TRUE,normalize=T)
  }
  taylor.diagram(my_eobs,my_molo,add=T,col="orange",ref.sd=TRUE,normalize=T)
  taylor.diagram(my_eobs,my_bola,add=T,col="blue",ref.sd=TRUE,normalize=T)
  print(paste0("+++OK"))
###

###
  # Store data for final Taylor diagram
  my_eobs_long <- c(my_eobs_long,my_eobs)
  my_era5_long <- c(my_era5_long,my_era5)
  my_molo_long <- c(my_molo_long,my_molo)
  my_bola_long <- c(my_bola_long,my_bola)  
###

###
  # End of my_year loop
  print(paste0("___End of ",my_year,"___"))
###
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
df <- data.frame(c(cor_era5,rmse_era5,bias_era5,verify_era5$ME,verify_era5$MAE,verify_era5$MSE), c(cor_molo,rmse_molo,bias_molo,verify_molo$ME,verify_molo$MAE,verify_molo$MSE), c(cor_bola,rmse_bola,bias_bola,verify_bola$ME,verify_bola$MAE,verify_bola$MSE))
write.table(df, "work/verification_scores.csv", col.names = c("ERA5-Land","MOLOCH","BOLAM"), row.names = c("Correlation","RMSE","Bias","Mean Error","Mean Absolute Error","Mean Square Error"), sep=";")
print(paste0("+++OK"))

# Set output EPS file
print(paste0("+++Plotting the summary plot..."))
setEPS()
postscript(paste0("work/Taylor_Diag_normalized_final.eps"))
taylor.diagram(my_eobs_long,my_era5_long*100,add=F,col="red",   ref.sd=TRUE,normalize=T, main="", pcex=1.5,grad.corr.lines=c(0.4,0.5,0.6,0.7,0.8))
taylor.diagram(my_eobs_long,my_molo_long,    add=T,col="orange",ref.sd=TRUE,normalize=T, pcex=1.5)
taylor.diagram(my_eobs_long,my_bola_long,    add=T,col="blue",  ref.sd=TRUE,normalize=T, pcex=1.5)
lpos<-sd(my_eobs_long)
legend("topright",legend=c("ERA5-Land","MOLOCH","BOLAM"),pch=19,col=c("red","orange","blue"), cex=1)
dev.off()
print(paste0("Ok"))

# Byebye
print(paste0("Byebye"))
quit()
