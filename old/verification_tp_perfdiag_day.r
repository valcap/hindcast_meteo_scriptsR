# Load libraries
library(raster)
library(ncdf4)
library(verification)
library(Metrics)
library(plotrix)
library(lubridate)
library(dplyr)
# add here any other packages you need

# Global variables
my_working_dir  <- paste0("/OCEANASTORE/progetti/spitbran")
my_figs_dir     <- paste0("/OCEANASTORE/progetti/spitbran/work/figs")
my_output_dir   <- paste0("/OCEANASTORE/progetti/spitbran/work/res")
my_shp_dir      <- paste0("/OCEANASTORE/progetti/spitbran/work/shp")
ini_year        <- 1981
end_year        <- 2019
vec_years       <- seq(ini_year,end_year)
sf_era5         <- 1000

# Working directory
setwd(my_working_dir)

# Load the point locations...
my_anag_file     <- paste0(my_working_dir,"/E-OBS/points050_moloch_extent.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))
print(paste0(""))

######################
# MAIN LOOP OVER YEARS
for (my_year in vec_years) {
  print(paste0("___Running year ",my_year,"___"))
  if (my_year==2001) {
    next
  }

###
  # Loading raster file names
  print(paste0("+++Reading raster files..."))
  # Raster files definition
  file_eobs   <- paste0("E-OBS/tp/E-OBS_day_tp_",my_year,".nc")
  file_era5   <- paste0("ERA5-Land/tp/ERA5-Land_day_tp_",my_year,".nc")
  file_molo   <- paste0("work/moloch/tp/moloch_day_tp_",my_year,".nc")
  file_bola   <- paste0("work/bolam/tp/bolam_day_tp_",my_year,".nc") 
  # Allocate rasters
  eobs_orig   <- NULL
  era5_orig   <- NULL
  molo_orig   <- NULL
  bola_orig   <- NULL
  eobs_orig   <- brick(file_eobs)
  era5_orig   <- brick(file_era5)
  molo_orig   <- brick(file_molo)
  bola_orig   <- brick(file_bola)
  print(paste0("+++OK"))
###

###
  # Get extent on the basis of the point coordinates
  print(paste0("+++Cropping extents and resampling rasters..."))
  my_extent           <- extent(pointCoordinates)
  # Crop E-OBS grid on my_extent and regrid (if needed)
  eobs_crop           <- crop(eobs_orig, my_extent)
  eobs_resam          <- eobs_crop
  # Crop ERA5-Land grid on my_extent and regrid (if needed)
  era5_crop           <- crop(era5_orig, my_extent)
  era5_resam          <- resample(era5_crop, eobs_crop, method="ngb")
  # Crop BOLAM grid on my_extent and regrid (if needed)
  bola_crop           <- crop(bola_orig, my_extent)
  bola_resam          <- resample(bola_crop, eobs_crop, method="bilinear")
  # Crop MOLOCH grid on my_extent and regrid (if needed)
  molo_crop           <- crop(molo_orig, my_extent)
  molo_resam          <- resample(molo_crop, eobs_crop, method="bilinear")
  print(paste0("+++OK"))
###

###
  # Select only data from the current year (except the new year's day)
  print(paste0("+++Adjusting name conventions..."))
  names(eobs_resam)     <- as.Date(names(eobs_resam),format="X%Y.%m.%d")
  names(era5_resam)     <- as.Date(names(era5_orig),format="X%Y.%m.%d")
  names(bola_resam)     <- as.Date(names(bola_resam),format="X%Y.%m.%d.%H.%M.%S")
  names(molo_resam)     <- names(bola_resam)
  print(paste0("+++...and selecting data from the current year..."))
  # E-OBS
  r <- setZ(eobs_resam, as.Date(names(eobs_resam),format="X%Y.%m.%d"))
  eobs_ok <- subset(r, which(getZ(r) >= paste0(my_year,'-01-02') & (getZ(r) <= paste0(my_year,'-12-31'))))
  # ERA5
  r <- setZ(era5_resam, as.Date(names(era5_resam),format="X%Y.%m.%d"))
  era5_ok <- subset(r, which(getZ(r) >= paste0(my_year,'-01-02') & (getZ(r) <= paste0(my_year,'-12-31'))))
  # MOLOCH
  r <- setZ(molo_resam, as.Date(names(molo_resam),format="X%Y.%m.%d"))
  molo_ok <- subset(r, which(getZ(r) >= paste0(my_year,'-01-02') & (getZ(r) <= paste0(my_year,'-12-31'))))
  # BOLAM
  r <- setZ(bola_resam, as.Date(names(bola_resam),format="X%Y.%m.%d"))
  bola_ok <- subset(r, which(getZ(r) >= paste0(my_year,'-01-02') & (getZ(r) <= paste0(my_year,'-12-31'))))
  print(paste0("+++OK"))
###

###
  # Extracting data
  print(paste0("+++Extracting data..."))
  eobs_values <- extract(eobs_ok, pointCoordinates, df=TRUE)
  era5_values <- extract(era5_ok*sf_era5, pointCoordinates, df=TRUE)
  molo_values <- extract(molo_ok, pointCoordinates, df=TRUE)
  bola_values <- extract(bola_ok, pointCoordinates, df=TRUE)
  # Re-organizing data in a single column vector
  my_vector_eobs <- NULL; my_vector_era5 <- NULL; my_vector_bola <- NULL; my_vector_molo <- NULL
  for (i in seq(2,ncol(eobs_values))) {
    my_vector_eobs <- c(my_vector_eobs, as.vector(eobs_values[,i]))
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
  print(paste0("+++Producing performance diagrams..."))
  setEPS()
  file_taylor_eps <- paste0(my_figs_dir,"/tp_PerfDiag_",my_year,"-day.eps")
  postscript(file_taylor_eps, horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 9)
#  postscript(file_taylor_eps)
  par(mfrow=c(2,3))
  for (thres in c(1,5,10,20,30,50)) {  # thresholds from Bonanno et al (2019)
    eobs_binary          <- replace(my_vector_eobs, my_vector_eobs <= thres, 0)
    eobs_binary          <- replace(eobs_binary, eobs_binary  > thres, 1)
    era5_binary          <- replace(my_vector_era5, my_vector_era5 <= thres, 0)
    era5_binary          <- replace(era5_binary, era5_binary  > thres, 1)
    molo_binary          <- replace(my_vector_molo, my_vector_molo <= thres, 0)
    molo_binary          <- replace(molo_binary, molo_binary  > thres, 1)
    bola_binary          <- replace(my_vector_bola, my_vector_bola <= thres, 0)
    bola_binary          <- replace(bola_binary, bola_binary  > thres, 1)
    performance.diagram(main = paste0('Threshold = ',thres,' mm/24-hour'), cex.main=0.99, cex.lab=1.1, cex.axis=0.99)
    my_table_stats     <- table.stats(eobs_binary,era5_binary)
    points(1 - my_table_stats$FAR, my_table_stats$POD, col = "red", cex = 2.5, pch=16)
    my_table_stats     <- table.stats(eobs_binary,bola_binary)
    points(1 - my_table_stats$FAR, my_table_stats$POD, col = "blue", cex = 2.5, pch=16)
    my_table_stats     <- table.stats(eobs_binary,molo_binary)
    points(1 - my_table_stats$FAR, my_table_stats$POD, col = "orange", cex = 2.5, pch=16)
    eobs_binary<-NULL; era5_binary<-NULL; molo_binary<-NULL; bola_binary<-NULL
    # Draw the legend on the plot
    if (thres==1) {
      legend("topleft",legend=c("ERA5-Land","MOLOCH","BOLAM"), pch=19, col=c("red","orange","blue"), bg="white")
    }
  }
  dev.off()
  print(paste0("+++OK"))
###

###
  # End of my_year loop
  print(paste0("___End of ",my_year,"___"))
}
# END OF MAIN LOOP OVER YEARS
#############################

# Byebye
print(paste0("Byebye"))
quit()

