# Purpose
print(paste0("RUNNING PERFORMANCE DIAGRAMS"))

# Load libraries
library(raster)
library(ncdf4)
library(verification)

# Working directory
setwd("/OCEANASTORE/progetti/spitbran/")

# Load the point locations
anag_file <- paste0("LaMMA/stations_synop.csv")
anag_full <- read.csv(anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))
print(paste0(""))

# Initialize the performance diagram plot
#vector_thres <- c(1,2,5,10)    | provare questi
#vector_thres <- c(20,30,40,50) | due set di thresholds
vector_thres <- c(5,10,20,40)
#for (i in c(1,2,3,4)) {
#  my_pch <- 14+i
#  thres <- vector_thres[i]
#  print(paste0("---------Running threshold ",thres,"m/24-hour---------"))
#  setEPS()
#  postscript(paste0("work/Perfor_Diag_Thres-",thres,"mm.eps"))
#  performance.diagram(main = paste0("Precipitation threshold ",thres," mm/24-hour"))

  # Initialize variables for final performance diagram
  my_eobs_long <- NULL 
  my_era5_long <- NULL 
  my_molo_long <- NULL
  my_bola_long <- NULL

  # MAIN LOOP OVER YEARS in `seq 1981 2019`
#  for (my_year in seq(1981,2019,1)) {
  for (my_year in seq(1981,2000,1)) {
    print(paste0("___Running year ",my_year,"___"))
  
    # Raster files definition
    file_eobs   <- paste0("E-OBS/tp/E-OBS_day_tp_",my_year,".nc")
    file_era5   <- paste0("ERA5-Land/tp/ERA5-Land_",my_year,"-daysum_tp.nc")
    file_moloch <- paste0("work/moloch/moloch_",my_year,"_daysum_apcp.nc")
    file_bolam  <- paste0("work/bolam/bolam_",my_year,"_daysum_apcp.nc") 
    
    print(paste0("+++Reading and resampling files..."))
    # Allocate rasters
    eobs_orig   <- stack(file_eobs)
    era5_orig   <- stack(file_era5)
    moloch_orig <- stack(file_moloch)
    bolam_orig  <- stack(file_bolam)
  
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
  
    # Re-define the layer names of each raster
         # BOLAM & MOLOCH have a 2001 missing day
         # it would be annoying to handle this exception.
         # Skipping it is faster...
         if ((my_year == "2001")) { next }
    if ( (my_year == "2020") || (my_year == "2016") || (my_year == "2012")
    || (my_year == "2008") || (my_year == "2004") || (my_year == "2000")
    || (my_year == "1996") || (my_year == "1992") || (my_year == "1988")
    || (my_year == "1984") || (my_year == "1980") ) {
       my_raster_names     <- seq(1,366)
       names(eobs_reduced) <- my_raster_names
       my_raster_names     <- seq(1,367)
       names(era5_resam)   <- my_raster_names
       names(era5_reduced) <- my_raster_names
       my_raster_names     <- seq(1,366)
       names(moloch_resam) <- my_raster_names
       names(moloch_orig)  <- my_raster_names
       names(bolam_resam)  <- my_raster_names
       names(bolam_orig)   <- my_raster_names
    } else {
       my_raster_names     <- seq(1,365)
       names(eobs_reduced) <- my_raster_names
       if ( (my_year == "1989") || (my_year == "1999") || (my_year == "2009") || (my_year == "2019")) {
         my_raster_names     <- seq(1,365)
       } else {
         my_raster_names     <- seq(1,366)
       }
       names(era5_resam)   <- my_raster_names
       names(era5_reduced) <- my_raster_names
       my_raster_names     <- seq(1,365)
       names(moloch_resam) <- my_raster_names
       names(moloch_orig)  <- my_raster_names
       names(bolam_resam)  <- my_raster_names
       names(bolam_orig)   <- my_raster_names
    }
  
    # Extracting data
    print(paste0("+++Extracting data..."))
    my_eobs <- NULL
    my_era5 <- NULL
    my_molo <- NULL
    my_bola <- NULL
    for (i in 2:365) {
      # Get E-OBS Values for a specific time step
      eobs_values <- extract(eobs_reduced[[i]], pointCoordinates)
      # Get ERA5 values for a specific time step
      era5_values <- extract(era5_resam[[i]], pointCoordinates)
      # Get MOLOCH values for a specific time step
      molo_values <- extract(moloch_resam[[i-1]], pointCoordinates)
      # Get BOLAM values for a specific time step
      bola_values <- extract(bolam_resam[[i-1]], pointCoordinates)
      my_eobs <- c(my_eobs,eobs_values)
      my_era5 <- c(my_era5,era5_values)
      my_molo <- c(my_molo,molo_values)
      my_bola <- c(my_bola,bola_values)
    }
    print(paste0("+++OK"))

#    # Plot points in the performance diagram
#    print(paste0("+++Plotting data..."))
#    # Create binary vector for E-OBS
#    my_eobs_thres <- replace(my_eobs, my_eobs <= thres,0)
#    my_eobs_thres <- replace(my_eobs_thres, my_eobs > thres,1)
#    # Create binary vector for ERA5-Land
#    my_era5_thres <- replace(my_era5*100, my_era5*100 <= thres,0)
#    my_era5_thres <- replace(my_era5_thres, my_era5*100 > thres,1)
#    # Create binary vector for MOLOCH
#    my_molo_thres <- replace(my_molo, my_molo <= thres,0)
#    my_molo_thres <- replace(my_molo_thres, my_molo > thres,1)
#    # Create binary vector for BOLAM
#    my_bola_thres <- replace(my_bola, my_bola <= thres,0)
#    my_bola_thres <- replace(my_bola_thres, my_bola > thres,1)
#    # Plot performance diagram
#    my_tab_era5 <- table.stats(my_eobs_thres, my_era5_thres)
#    my_tab_molo <- table.stats(my_eobs_thres, my_molo_thres)
#    my_tab_bola <- table.stats(my_eobs_thres, my_bola_thres)
#    points(1 - my_tab_era5$FAR, my_tab_era5$POD, col = "red", cex = 1.5, pch=my_pch, add=T)
#    points(1 - my_tab_molo$FAR, my_tab_molo$POD, col = "orange", cex = 1.5, pch=my_pch, add=T)
#    points(1 - my_tab_bola$FAR, my_tab_bola$POD, col = "blue", cex = 1.5, pch=my_pch, add=T)
#    print(paste0("+++OK"))
  
    # Store data for final Taylor diagram
    my_eobs_long <- c(my_eobs_long,my_eobs)
    my_era5_long <- c(my_era5_long,my_era5)
    my_molo_long <- c(my_molo_long,my_molo)
    my_bola_long <- c(my_bola_long,my_bola)  
  
    # End of my_year loop
    print(paste0("___End of ",my_year,"___"))
  } 
#  # Draw the legend on the plot and save it
#  legend("topleft",legend=c("ERA5-Land","MOLOCH","BOLAM"),pch=my_pch,col=c("red","orange","blue"), bg='azure',cex = 1.5)
#  dev.off()

#  print(paste0("---------End of threshold ",thres,"m/24-hour---------"))
#}

# Set output EPS file
setEPS()
postscript(paste0("work/Perfor_Diag_final.eps"))
performance.diagram(main = paste0("Precipitation mm/24-hour"))
for (i in c(1,2,3,4)) {
  my_pch <- 14+i
  thres <- vector_thres[i]
  # Create binary vector for E-OBS
  my_eobs_thres <- replace(my_eobs_long, my_eobs_long <= thres,0)
  my_eobs_thres <- replace(my_eobs_thres, my_eobs_long > thres,1)
  # Create binary vector for ERA5-Land
  my_era5_thres <- replace(my_era5_long*100, my_era5_long*100 <= thres,0)
  my_era5_thres <- replace(my_era5_thres,    my_era5_long*100 > thres,1)
  # Create binary vector for MOLOCH
  my_molo_thres <- replace(my_molo_long, my_molo_long <= thres,0)
  my_molo_thres <- replace(my_molo_thres, my_molo_long > thres,1)
  # Create binary vector for BOLAM
  my_bola_thres <- replace(my_bola_long, my_bola_long <= thres,0)
  my_bola_thres <- replace(my_bola_thres, my_bola_long > thres,1)
  # Plot performance diagram
  my_tab_era5 <- table.stats(my_eobs_thres, my_era5_thres)
  my_tab_molo <- table.stats(my_eobs_thres, my_molo_thres)
  my_tab_bola <- table.stats(my_eobs_thres, my_bola_thres)
  points(1 - my_tab_era5$FAR, my_tab_era5$POD, col = "red", cex = 2.0, pch=my_pch, add=T)
  points(1 - my_tab_molo$FAR, my_tab_molo$POD, col = "orange", cex = 2.0, pch=my_pch, add=T)
  points(1 - my_tab_bola$FAR, my_tab_bola$POD, col = "blue", cex = 2.0, pch=my_pch, add=T)
}
# Draw the legend on the plot and save it
legend("topleft", legend=c("ERA5-Land","MOLOCH","BOLAM"), pch=15, col=c("red","orange","blue"), bg='azure')
lab<-NULL
for (i in c(1,2,3,4)) {
  thres <- vector_thres[i]
  lab[i] <- paste0(thres," mm/24-hour")
}
legend("topright", legend=c(lab[1],lab[2],lab[3],lab[4]), pch=c(15,16,17,18), col=c("black","black","black","black"), bg='azure')
dev.off()
  
# Byebye
quit()

