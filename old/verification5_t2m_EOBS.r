# Purpose

# Load libraries
library(raster)
library(ncdf4)
library(verification)
library(Metrics)
library(plotrix)
library(lubridate)
library(dplyr)

# Global variables
my_working_dir  <- paste0("/OCEANASTORE/progetti/spitbran")
my_figs_dir     <- paste0("/OCEANASTORE/progetti/spitbran/work/figs")
my_output_dir   <- paste0("/OCEANASTORE/progetti/spitbran/work/res")
my_shp_dir      <- paste0("/OCEANASTORE/progetti/spitbran/work/shp")
TD_normal       <- TRUE
ini_year        <- 1981
end_year        <- 2019
vec_years       <- seq(ini_year,end_year)
add_offset      <- 273.15
var             <- 'tn'
my_rmse <- function(x) {  
  sqrt( mean( (x)^2 , na.rm = TRUE ) )
  }


# Working directory
setwd(my_working_dir)

# Load the point locations
my_anag_file     <- paste0(my_working_dir,"/E-OBS/points050_moloch_extent.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))

######################
# MAIN LOOP OVER YEARS
for (my_year in vec_years) {
  if (my_year==2001) {
    print(paste0("___Skipping ",my_year,"___"))
    next
  }
  print(paste0("___Running year ",my_year,"___"))

###
  # Loading raster file names
  print(paste0("+++Reading raster files..."))
  # Raster files definition
  file_eobs   <- paste0("E-OBS/t2m/E-OBS_day_",var,"_",my_year,".nc")
  file_era5   <- paste0("ERA5-Land/t2m/ERA5-Land_",my_year,"-daystat_",var,".nc")
  file_molo   <- paste0("work/moloch/moloch_",my_year,"_daystat_",var,".nc")
  file_bola   <- paste0("work/bolam/bolam_",my_year,"_daystat_",var,".nc") 

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
#  era5_resam          <- resample(era5_crop, eobs_crop, method="ngb")
  era5_resam          <- era5_crop
  # Crop BOLAM grid on my_extent and regrid (if needed)
  bola_crop           <- crop(bola_orig, my_extent)
#  bola_resam          <- resample(bola_crop, eobs_crop, method="bilinear")
  bola_resam          <- bola_crop
  # Crop MOLOCH grid on my_extent and regrid (if needed)
  molo_crop           <- crop(molo_orig, my_extent)
#  molo_resam          <- resample(molo_crop, eobs_crop, method="bilinear")
  molo_resam          <- molo_crop
  print(paste0("+++OK"))
###

###
  # select only data from the current year
  print(paste0("+++Adjusting name conventions..."))
  names(eobs_resam)     <- as.Date(names(eobs_resam),format="X%Y.%m.%d")
  names(era5_resam)     <- as.Date(names(era5_resam),format="X%Y.%m.%d.%H.%M.%S")
  names(bola_resam)     <- as.Date(names(bola_resam),format="X%Y.%m.%d.%H.%M.%S")
  names(molo_resam)     <- names(bola_resam)
  # select only data from the current year
  print(paste0("+++selecting data from the current year..."))
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
  eobs_values    <- as.data.frame(extract(eobs_ok, pointCoordinates))
  era5_values    <- as.data.frame(extract(era5_ok-add_offset, pointCoordinates))
  molo_values    <- as.data.frame(extract(molo_ok-add_offset, pointCoordinates))
  bola_values    <- as.data.frame(extract(bola_ok-add_offset, pointCoordinates))
  print(paste0("+++OK"))
###

diff_era5_eobs <- (era5_values)-(eobs_values)
pippo          <- apply(diff_era5_eobs,1,my_rmse)
rmse_era5      <- mean(pippo,na.rm = TRUE)
bias_era5      <- mean(apply(era5_values,1,mean),na.rm=T)/mean(apply(eobs_values,1,mean),na.rm=T)
my_corr_mean <- NULL
for (my_ind in seq(1,length(pointCoordinates))) {
  if (all(is.na(t(eobs_values[my_ind,]))) == TRUE ) {next}
  if (all(is.na(t(era5_values[my_ind,]))) == TRUE ) {next}
  tmp          <- cor(t(eobs_values[my_ind,]), t(era5_values[my_ind,]), use = "complete.obs")
  my_corr_mean <- c(my_corr_mean,tmp)
}
corr_era5      <- mean(my_corr_mean, na.rm = TRUE)

diff_molo_eobs <- (molo_values)-(eobs_values)
pippo          <- apply(diff_molo_eobs,1,my_rmse)
rmse_molo      <- mean(pippo,na.rm = TRUE)
bias_molo      <- mean(apply(molo_values,1,mean),na.rm=T)/mean(apply(eobs_values,1,mean),na.rm=T)
my_corr_mean <- NULL
for (my_ind in seq(1,length(pointCoordinates))) {
  if (all(is.na(t(eobs_values[my_ind,]))) == TRUE ) {next}
  if (all(is.na(t(molo_values[my_ind,]))) == TRUE ) {next}
  tmp          <- cor(t(eobs_values[my_ind,]), t(molo_values[my_ind,]), use = "complete.obs")
  my_corr_mean <- c(my_corr_mean,tmp)
}
corr_molo      <- mean(my_corr_mean, na.rm = TRUE)


diff_bola_eobs <- (bola_values)-(eobs_values)
pippo          <- apply(diff_bola_eobs,1,my_rmse)
rmse_bola      <- mean(pippo,na.rm = TRUE)
bias_bola      <- mean(apply(bola_values,1,mean),na.rm=T)/mean(apply(eobs_values,1,mean),na.rm=T)
my_corr_mean <- NULL
for (my_ind in seq(1,length(pointCoordinates))) {
  if (all(is.na(t(eobs_values[my_ind,]))) == TRUE ) {next}
  if (all(is.na(t(bola_values[my_ind,]))) == TRUE ) {next}
  tmp          <- cor(t(eobs_values[my_ind,]), t(bola_values[my_ind,]), use = "complete.obs")
  my_corr_mean <- c(my_corr_mean,tmp)
}
corr_bola      <- mean(my_corr_mean, na.rm = TRUE)

print(paste("VER_CORR",my_year,sprintf("%.2f",corr_era5),sprintf("%.2f",corr_molo),sprintf("%.2f",corr_bola),sep=";"))
print(paste("VER_RMSE",my_year,sprintf("%.2f",rmse_era5),sprintf("%.2f",rmse_molo),sprintf("%.2f",rmse_bola),sep=";"))
print(paste("VER_BIAS",my_year,sprintf("%.2f",bias_era5),sprintf("%.2f",bias_molo),sprintf("%.2f",bias_bola),sep=";"))

###
  # End of my_year loop
  print(paste0("___End of ",my_year,"___"))
}
# END OF MAIN LOOP OVER YEARS
#############################

# Quit
quit()
