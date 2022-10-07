# Load libraries
library(raster)
library(ncdf4)
library("verification")
library(Metrics)
library(plotrix)
library("lubridate")
library(dplyr)
library(rgdal)
# add here any other packages you need

# Global variables
u               <- 'capecchi'
my_working_dir  <- paste0("/home/",u,"/art_hindcast_meteo")
my_figs_dir     <- paste0(my_working_dir,"/figs")
# Intervallo comune a GRIPHO, ARCIS and E-OBS
ini_year        <- 2001
end_year        <- 2015
vec_years       <- seq(ini_year,end_year)
sf_era5         <- 1000
vec_thres       <- c(1,2,5,10,20,50)
vec_labels      <- c("(a)","(b)","(c)","(d)","(e)","(f)")

# Working directory
setwd(my_working_dir)

obse_dataset    <- "GRIPHO"
# Load the point locations...
my_anag_file     <- paste0(my_working_dir,"/DATA/points025_italy_random-subset2.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
######################
# MAIN LOOP OVER YEARS
for (my_year in vec_years) {
  print(paste0("___Running year ",my_year,"___"))
  if (my_year==2009) {
    next
  }
  my_vector_obse  <- NULL
  my_vector_era5  <- NULL
  my_vector_bola  <- NULL
  my_vector_molo  <- NULL
  
###
  # Loading raster file names
  print(paste0("+++Reading raster files..."))
  # Raster files definition
  file_obse   <- paste0("DATA/GRIPHO/tp/gripho-v2-3km_day_tp_",my_year,".nc")
  file_era5   <- paste0("DATA/ERA5-Land/tp/ERA5-Land_day_tp_",my_year,".nc")
  file_molo   <- paste0("DATA/MOLOCH/tp/moloch_day_tp_",my_year,".nc")
  file_bola   <- paste0("DATA/BOLAM/tp/bolam_day_tp_",my_year,".nc") 
  # Allocate rasters
  obse_orig   <- NULL
  era5_orig   <- NULL
  molo_orig   <- NULL
  bola_orig   <- NULL
  obse_orig   <- brick(file_obse)
  era5_orig   <- brick(file_era5)
  molo_orig   <- brick(file_molo)
  bola_orig   <- brick(file_bola)
  print(paste0("+++OK"))
###

###
  print(paste0("+++Subsetting and adjusting name conventions..."))
  obse_subset           <- subset(obse_orig,2:365)
  era5_subset           <- subset(era5_orig,2:365)
  molo_subset           <- subset(molo_orig,1:364)
  bola_subset           <- subset(bola_orig,1:364)
  ##################################################################
  # ATTENZIONE, POICHE' ASSUMIAMO LA CONVENZIONE DI ERA5 PER LA DATA
  # %y-%m-%d SI RIFERISCE AL GIORNO PRECEDENTE (CIOE' %(d-1)
  names(obse_subset)     <- as.Date(names(era5_subset), format="X%Y.%m.%d")
  names(era5_subset)     <- as.Date(names(era5_subset), format="X%Y.%m.%d")
  names(molo_subset)     <- as.Date(names(era5_subset), format="X%Y.%m.%d")
  names(bola_subset)     <- as.Date(names(era5_subset), format="X%Y.%m.%d")
  ##################################################################
###
  
###
  # Get extent on the basis of the point coordinates
  my_extent           <- extent(pointCoordinates)
  print(paste0("+++Cropping extents and resampling rasters..."))
  # Crop E-OBS grid on my_extent and regrid (if needed)
  obse_crop           <- crop(obse_subset, my_extent)
  obse_resam          <- as.integer(obse_crop)
  # Crop ERA5-Land grid on my_extent and regrid (if needed)
  era5_crop           <- crop(era5_subset, my_extent)
  era5_resam          <- as.integer(era5_crop*sf_era5)
#  era5_resam          <- resample(era5_crop*sf_era5, obse_crop, method="ngb")
  # Crop BOLAM grid on my_extent and regrid (if needed)
  bola_crop           <- crop(bola_subset, my_extent)
  bola_resam          <- as.integer(bola_crop)
#  bola_resam          <- resample(as.integer(bola_crop), era5_crop, method="bilinear")
  # Crop MOLOCH grid on my_extent and regrid (if needed)
  molo_crop           <- crop(molo_subset, my_extent)
  molo_resam          <- as.integer(molo_crop)
#  molo_resam          <- resample(as.integer(molo_crop), obse_crop, method="bilinear")
  print(paste0("+++OK"))
###

###
  # Extracting data
  print(paste0("+++Extracting data..."))
  obse_values <- extract(as.integer(obse_resam), pointCoordinates, df=TRUE)
  era5_values <- extract(as.integer(era5_resam), pointCoordinates, df=TRUE)
  molo_values <- extract(as.integer(molo_resam), pointCoordinates, df=TRUE)
  bola_values <- extract(as.integer(bola_resam), pointCoordinates, df=TRUE)
  # Re-organizing data in a single column vector
  for (i in seq(2,ncol(obse_values))) {
    my_vector_obse <- c(my_vector_obse, as.vector(obse_values[,i]))
  }
  for (i in seq(2,ncol(era5_values))) {
      my_vector_era5 <- c(my_vector_era5, as.vector(era5_values[,i]))
  }
  for (i in seq(2,ncol(bola_values))) {
      my_vector_bola <- c(my_vector_bola, as.vector(bola_values[,i]))
  }
  for (i in seq(2,ncol(molo_values))) {
      my_vector_molo <- c(my_vector_molo, as.vector(molo_values[,i]))
  }
  print(paste0("+++OK"))
###

###
# PERFORMANCE DIAGRAMS
  print(paste0("+++Producing performance diagram..."))
  setEPS()
  file_taylor_eps <- paste0(my_figs_dir,"/tp-daily_PerfDiag",my_year,"_vs_",obse_dataset,".eps")
  postscript(file_taylor_eps, horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 9)
  par(mfrow=c(2,3))
  for (i in seq(1,6)) {  
    thres  <- vec_thres[i]  
    my_lab <- vec_labels[i]  
    obse_binary          <- replace(my_vector_obse, my_vector_obse <= thres, 0)
    obse_binary          <- replace(obse_binary, obse_binary  > thres, 1)
    era5_binary          <- replace(my_vector_era5, my_vector_era5 <= thres, 0)
    era5_binary          <- replace(era5_binary, era5_binary  > thres, 1)
    molo_binary          <- replace(my_vector_molo, my_vector_molo <= thres, 0)
    molo_binary          <- replace(molo_binary, molo_binary  > thres, 1)
    bola_binary          <- replace(my_vector_bola, my_vector_bola <= thres, 0)
    bola_binary          <- replace(bola_binary, bola_binary  > thres, 1)
    performance.diagram(main = paste0(my_lab,' Threshold = ',thres,' mm/day'), cex.main=0.99, cex.lab=1.1, cex.axis=0.99)
    my_table_stats     <- table.stats(obse_binary,era5_binary)
    points(1 - my_table_stats$FAR, my_table_stats$POD, col = "red", cex = 2.5, pch=16)
    print(paste0(my_table_stats$FAR,my_table_stats$POD,sep=" "))
    my_table_stats     <- table.stats(obse_binary,bola_binary)
    points(1 - my_table_stats$FAR, my_table_stats$POD, col = "blue", cex = 2.5, pch=16)
    my_table_stats     <- table.stats(obse_binary,molo_binary)
    points(1 - my_table_stats$FAR, my_table_stats$POD, col = "orange", cex = 2.5, pch=16)
    obse_binary<-NULL; era5_binary<-NULL; molo_binary<-NULL; bola_binary<-NULL
    # Draw the legend on the plot
    if (thres==1) {
      legend("topleft",legend=c("ERA5-Land","MOLOCH","BOLAM"), pch=19, col=c("red","orange","blue"), bg="white")
    }
  }
  dev.off()
  ###
  # End of my_year loop
  print(paste0("___End of ",my_year,"___"))
  # END OF MAIN LOOP OVER YEARS
#############################
}
print(paste0("+++OK ",obse_dataset))
###


quit()






###############################################################
###############################################################
###############################################################
obse_dataset    <- "ARCIS"
# Load the point locations...
my_anag_file     <- paste0(my_working_dir,"/DATA/points025_arcis_extent.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
######################
# MAIN LOOP OVER YEARS
for (my_year in vec_years) {
  print(paste0("___Running year ",my_year,"___"))
  if (my_year==2009) {
    next
  }
  my_vector_obse  <- NULL
  my_vector_era5  <- NULL
  my_vector_bola  <- NULL
  my_vector_molo  <- NULL
  
  ###
  # Loading raster file names
  print(paste0("+++Reading raster files..."))
  # Raster files definition
  file_obse   <- paste0("DATA/ARCIS/tp/ARCIS_day_tp_",my_year,".nc")
  file_era5   <- paste0("DATA/ERA5-Land/tp/ERA5-Land_day_tp_",my_year,".nc")
  file_molo   <- paste0("DATA/MOLOCH/tp/moloch_day_tp_",my_year,".nc")
  file_bola   <- paste0("DATA/BOLAM/tp/bolam_day_tp_",my_year,".nc") 
  # Allocate rasters
  obse_orig   <- NULL
  era5_orig   <- NULL
  molo_orig   <- NULL
  bola_orig   <- NULL
  obse_orig   <- brick(file_obse)
  era5_orig   <- brick(file_era5)
  molo_orig   <- brick(file_molo)
  bola_orig   <- brick(file_bola)
  print(paste0("+++OK"))
  ###
  
  ###
  print(paste0("+++Subsetting and adjusting name conventions..."))
  obse_subset           <- subset(obse_orig,2:365)
  era5_subset           <- subset(era5_orig,2:365)
  molo_subset           <- subset(molo_orig,1:364)
  bola_subset           <- subset(bola_orig,1:364)
  ##################################################################
  # ATTENZIONE, POICHE' ASSUMIAMO LA CONVENZIONE DI ERA5 PER LA DATA
  # %y-%m-%d SI RIFERISCE AL GIORNO PRECEDENTE (CIOE' %(d-1)
  names(obse_subset)     <- as.Date(names(era5_subset), format="X%Y.%m.%d")
  names(era5_subset)     <- as.Date(names(era5_subset), format="X%Y.%m.%d")
  names(molo_subset)     <- as.Date(names(era5_subset), format="X%Y.%m.%d")
  names(bola_subset)     <- as.Date(names(era5_subset), format="X%Y.%m.%d")
  ##################################################################
  ###
  
  ###
  # Get extent on the basis of the point coordinates
  my_extent           <- extent(pointCoordinates)
  print(paste0("+++Cropping extents and resampling rasters..."))
  # Crop E-OBS grid on my_extent and regrid (if needed)
  obse_crop           <- crop(obse_subset, my_extent)
  obse_resam          <- as.integer(obse_crop)
  # Crop ERA5-Land grid on my_extent and regrid (if needed)
  era5_crop           <- crop(era5_subset, my_extent)
  era5_resam          <- as.integer(era5_crop*sf_era5)
  #  era5_resam          <- resample(era5_crop*sf_era5, obse_crop, method="ngb")
  # Crop BOLAM grid on my_extent and regrid (if needed)
  bola_crop           <- crop(bola_subset, my_extent)
  bola_resam          <- as.integer(bola_crop)
  #bola_resam          <- resample(as.integer(bola_crop), era5_crop, method="bilinear")
  # Crop MOLOCH grid on my_extent and regrid (if needed)
  molo_crop           <- crop(molo_subset, my_extent)
  molo_resam          <- as.integer(molo_crop)
  #molo_resam          <- resample(as.integer(molo_crop), obse_crop, method="bilinear")
  print(paste0("+++OK"))
  ###
  
  ###
  # Extracting data
  print(paste0("+++Extracting data..."))
  obse_values <- extract(as.integer(obse_resam), pointCoordinates, df=TRUE)
  era5_values <- extract(as.integer(era5_resam), pointCoordinates, df=TRUE)
  molo_values <- extract(as.integer(molo_resam), pointCoordinates, df=TRUE)
  bola_values <- extract(as.integer(bola_resam), pointCoordinates, df=TRUE)
  # Re-organizing data in a single column vector
  for (i in seq(2,ncol(obse_values))) {  
    my_vector_obse <- c(my_vector_obse, as.vector(obse_values[,i]))
  }
  for (i in seq(2,ncol(era5_values))) {  
    my_vector_era5 <- c(my_vector_era5, as.vector(era5_values[,i]))
  }
  for (i in seq(2,ncol(bola_values))) {
    my_vector_bola <- c(my_vector_bola, as.vector(bola_values[,i]))
  }
  for (i in seq(2,ncol(molo_values))) {
    my_vector_molo <- c(my_vector_molo, as.vector(molo_values[,i]))
  }
  print(paste0("+++OK"))
  ###
  
  ###
  # End of my_year loop
  print(paste0("___End of ",my_year,"___"))
# END OF MAIN LOOP OVER YEARS
#############################

###
# PERFORMANCE DIAGRAMS
  print(paste0("+++Producing performance diagram..."))
  setEPS()
  file_taylor_eps <- paste0(my_figs_dir,"/tp-daily_PerfDiag",my_year,"_vs_",obse_dataset,".eps")
  postscript(file_taylor_eps, horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 9)
  par(mfrow=c(2,3))
  for (i in seq(1,6)) {  
    thres  <- vec_thres[i]  
    my_lab <- vec_labels[i]  
    obse_binary          <- replace(my_vector_obse, my_vector_obse <= thres, 0)
    obse_binary          <- replace(obse_binary, obse_binary  > thres, 1)
    era5_binary          <- replace(my_vector_era5, my_vector_era5 <= thres, 0)
    era5_binary          <- replace(era5_binary, era5_binary  > thres, 1)
    molo_binary          <- replace(my_vector_molo, my_vector_molo <= thres, 0)
    molo_binary          <- replace(molo_binary, molo_binary  > thres, 1)
    bola_binary          <- replace(my_vector_bola, my_vector_bola <= thres, 0)
    bola_binary          <- replace(bola_binary, bola_binary  > thres, 1)
    performance.diagram(main = paste0(my_lab,' Threshold = ',thres,' mm/day'), cex.main=0.99, cex.lab=1.1, cex.axis=0.99)
    my_table_stats     <- table.stats(obse_binary,era5_binary)
    points(1 - my_table_stats$FAR, my_table_stats$POD, col = "red", cex = 2.5, pch=16)
    my_table_stats     <- table.stats(obse_binary,bola_binary)
    points(1 - my_table_stats$FAR, my_table_stats$POD, col = "blue", cex = 2.5, pch=16)
    my_table_stats     <- table.stats(obse_binary,molo_binary)
    points(1 - my_table_stats$FAR, my_table_stats$POD, col = "orange", cex = 2.5, pch=16)
    obse_binary<-NULL; era5_binary<-NULL; molo_binary<-NULL; bola_binary<-NULL
  # Draw the legend on the plot
    if (thres==1) {
      legend("topleft",legend=c("ERA5-Land","MOLOCH","BOLAM"), pch=19, col=c("red","orange","blue"), bg="white")
    }
  }
dev.off()
}
print(paste0("+++OK"))
###

###############################################################
###############################################################
###############################################################
obse_dataset    <- "E-OBS"
# Load the point locations...
my_anag_file     <- paste0(my_working_dir,"/DATA/points025_italy_random-subset2.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
my_vector_obse  <- NULL
my_vector_era5  <- NULL
my_vector_bola  <- NULL
my_vector_molo  <- NULL
######################
# MAIN LOOP OVER YEARS
for (my_year in vec_years) {
  print(paste0("___Running year ",my_year,"___"))
  if (my_year==2009) {
    next
  }
  my_vector_obse  <- NULL
  my_vector_era5  <- NULL
  my_vector_bola  <- NULL
  my_vector_molo  <- NULL
  
  ###
  # Loading raster file names
  print(paste0("+++Reading raster files..."))
  # Raster files definition
  file_obse   <- paste0("DATA/E-OBS/tp/E-OBS_day_tp_",my_year,".nc")
  file_era5   <- paste0("DATA/ERA5-Land/tp/ERA5-Land_day_tp_",my_year,".nc")
  file_molo   <- paste0("DATA/MOLOCH/tp/moloch_day_tp_",my_year,".nc")
  file_bola   <- paste0("DATA/BOLAM/tp/bolam_day_tp_",my_year,".nc") 
  # Allocate rasters
  obse_orig   <- NULL
  era5_orig   <- NULL
  molo_orig   <- NULL
  bola_orig   <- NULL
  obse_orig   <- brick(file_obse)
  era5_orig   <- brick(file_era5)
  molo_orig   <- brick(file_molo)
  bola_orig   <- brick(file_bola)
  print(paste0("+++OK"))
  ###
  
  ###
  print(paste0("+++Subsetting and adjusting name conventions..."))
  obse_subset           <- subset(obse_orig,2:365)
  era5_subset           <- subset(era5_orig,2:365)
  molo_subset           <- subset(molo_orig,1:364)
  bola_subset           <- subset(bola_orig,1:364)
  ##################################################################
  # ATTENZIONE, POICHE' ASSUMIAMO LA CONVENZIONE DI ERA5 PER LA DATA
  # %y-%m-%d SI RIFERISCE AL GIORNO PRECEDENTE (CIOE' %(d-1)
  names(obse_subset)     <- as.Date(names(era5_subset), format="X%Y.%m.%d")
  names(era5_subset)     <- as.Date(names(era5_subset), format="X%Y.%m.%d")
  names(molo_subset)     <- as.Date(names(era5_subset), format="X%Y.%m.%d")
  names(bola_subset)     <- as.Date(names(era5_subset), format="X%Y.%m.%d")
  ##################################################################
  ###
  
  ###
  # Get extent on the basis of the point coordinates
  my_extent           <- extent(pointCoordinates)
  print(paste0("+++Cropping extents and resampling rasters..."))
  # Crop E-OBS grid on my_extent and regrid (if needed)
  obse_crop           <- crop(obse_subset, my_extent)
  obse_resam          <- as.integer(obse_crop)
  # Crop ERA5-Land grid on my_extent and regrid (if needed)
  era5_crop           <- crop(era5_subset, my_extent)
  era5_resam          <- as.integer(era5_crop*sf_era5)
  #  era5_resam          <- resample(era5_crop*sf_era5, obse_crop, method="ngb")
  # Crop BOLAM grid on my_extent and regrid (if needed)
  bola_crop           <- crop(bola_subset, my_extent)
  #  bola_resam          <- as.integer(bola_crop)
  bola_resam          <- resample(as.integer(bola_crop), obse_resam, method="bilinear")
  # Crop MOLOCH grid on my_extent and regrid (if needed)
  molo_crop           <- crop(molo_subset, my_extent)
  #  molo_resam          <- as.integer(molo_crop)
  molo_resam          <- resample(as.integer(molo_crop), obse_resam, method="bilinear")
  print(paste0("+++OK"))
  ###
  
  ###
  # Extracting data
  print(paste0("+++Extracting data..."))
  obse_values <- extract(as.integer(obse_resam), pointCoordinates, df=TRUE)
  era5_values <- extract(as.integer(era5_resam), pointCoordinates, df=TRUE)
  molo_values <- extract(as.integer(molo_resam), pointCoordinates, df=TRUE)
  bola_values <- extract(as.integer(bola_resam), pointCoordinates, df=TRUE)
  # Re-organizing data in a single column vector
  #  for (i in seq(2,ncol(obse_values))) {
  for (i in seq(2,300)) {  
    my_vector_obse <- c(my_vector_obse, as.vector(obse_values[,i]))
  }
  #  for (i in seq(2,ncol(era5_values))) {
  for (i in seq(2,300)) {  
    my_vector_era5 <- c(my_vector_era5, as.vector(era5_values[,i]))
  }
  #  for (i in seq(2,ncol(bola_values))) {
  for (i in seq(2,300)) {  
    my_vector_bola <- c(my_vector_bola, as.vector(bola_values[,i]))
  }
  #  for (i in seq(2,ncol(molo_values))) {
  for (i in seq(2,300)) {  
    my_vector_molo <- c(my_vector_molo, as.vector(molo_values[,i]))
  }
  print(paste0("+++OK"))
  ###
  
  ###
  # End of my_year loop
  print(paste0("___End of ",my_year,"___"))
# END OF MAIN LOOP OVER YEARS
#############################

###
# PERFORMANCE DIAGRAMS
print(paste0("+++Producing performance diagram..."))
setEPS()
file_taylor_eps <- paste0(my_figs_dir,"/tp-daily_PerfDiag",my_year,"_vs_",obse_dataset,".eps")
postscript(file_taylor_eps, horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 9)
par(mfrow=c(2,3))
for (i in seq(1,6)) {  
  thres  <- vec_thres[i]  
  my_lab <- vec_labels[i]  
  obse_binary          <- replace(my_vector_obse, my_vector_obse <= thres, 0)
  obse_binary          <- replace(obse_binary, obse_binary  > thres, 1)
  era5_binary          <- replace(my_vector_era5, my_vector_era5 <= thres, 0)
  era5_binary          <- replace(era5_binary, era5_binary  > thres, 1)
  molo_binary          <- replace(my_vector_molo, my_vector_molo <= thres, 0)
  molo_binary          <- replace(molo_binary, molo_binary  > thres, 1)
  bola_binary          <- replace(my_vector_bola, my_vector_bola <= thres, 0)
  bola_binary          <- replace(bola_binary, bola_binary  > thres, 1)
  performance.diagram(main = paste0(my_lab,' Threshold = ',thres,' mm/day'), cex.main=0.99, cex.lab=1.1, cex.axis=0.99)
  my_table_stats     <- table.stats(obse_binary,era5_binary)
  points(1 - my_table_stats$FAR, my_table_stats$POD, col = "red", cex = 2.5, pch=16)
  my_table_stats     <- table.stats(obse_binary,bola_binary)
  points(1 - my_table_stats$FAR, my_table_stats$POD, col = "blue", cex = 2.5, pch=16)
  my_table_stats     <- table.stats(obse_binary,molo_binary)
  points(1 - my_table_stats$FAR, my_table_stats$POD, col = "orange", cex = 2.5, pch=16)
  obse_binary<-NULL; era5_binary<-NULL; molo_binary<-NULL; bola_binary<-NULL
  # Draw the legend on the plot
  if (thres==1) {
    legend("topleft",legend=c("ERA5-Land","MOLOCH","BOLAM"), pch=19, col=c("red","orange","blue"), bg="white")
  }
}
dev.off()
}
print(paste0("+++OK"))
###

# Byebye
print(paste0("Byebye"))
#quit()
