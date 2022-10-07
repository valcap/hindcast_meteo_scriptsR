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
eobs_values_ts  <- NULL 
era5_values_ts  <- NULL
molo_values_ts  <- NULL
bola_values_ts  <- NULL
cutoff_NA       <- 50

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
  # select only data from the current year (except the new year day)
  print(paste0("+++Adjusting name conventions..."))
  names(eobs_resam)     <- as.Date(names(eobs_resam),format="X%Y.%m.%d")
  names(era5_resam)     <- as.Date(names(era5_orig),format="X%Y.%m.%d")
  names(bola_resam)     <- as.Date(names(bola_resam),format="X%Y.%m.%d.%H.%M.%S")
  names(molo_resam)     <- names(bola_resam)
  # select only data from the current year
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
  # Creating yearly accumulated rainfall
#  print(paste0("+++Creating yearly accumulated rainfall"))
#  index_years          <- format(as.Date(names(molo_ok), format = "X%Y.%m.%d"), format = "%Y")
#  eobs_ok_year         <- stackApply(eobs_ok, index_years, fun = sum)
#  names(eobs_ok_year)  <- my_year
#  era5_ok_year         <- stackApply(era5_ok, index_years, fun = sum)
#  names(era5_ok_year)  <- my_year
#  molo_ok_year         <- stackApply(molo_ok, index_years, fun = sum)
#  names(molo_ok_year)  <- my_year
#  bola_ok_year         <- stackApply(bola_ok, index_years, fun = sum)
#  names(bola_ok_year)  <- my_year
#  print(paste0("+++OK"))
###

###
  # Extracting data
  print(paste0("+++Extracting data..."))
  eobs_values <- as.data.frame(extract(eobs_ok_year, pointCoordinates))[,1]
  era5_values <- as.data.frame(extract(era5_ok_year*sf_era5, pointCoordinates))[,1]
  molo_values <- as.data.frame(extract(molo_ok_year, pointCoordinates))[,1]
  bola_values <- as.data.frame(extract(bola_ok_year, pointCoordinates))[,1]
  #
  eobs_values <- replace(eobs_values,eobs_values<cutoff_NA, NA)
  era5_values <- replace(era5_values,era5_values<cutoff_NA, NA)
  molo_values <- replace(molo_values,molo_values<cutoff_NA, NA)
  bola_values <- replace(bola_values,bola_values<cutoff_NA, NA)
  print(paste0("+++OK"))
###

###
  # Standard statistical verification scores for continuos variables
  print(paste0("+++Calculating standard statistical verification scores..."))
  # Correlation
  cor_era5 <- cor(era5_values, eobs_values, use = "complete.obs")
  cor_molo <- cor(molo_values, eobs_values, use = "complete.obs")
  cor_bola <- cor(bola_values, eobs_values, use = "complete.obs")
  # RMSE
  rmse_era5 <- sqrt( mean( (era5_values-eobs_values)^2 , na.rm = TRUE ) ) 
  rmse_molo <- sqrt( mean( (molo_values-eobs_values)^2 , na.rm = TRUE ) ) 
  rmse_bola <- sqrt( mean( (bola_values-eobs_values)^2 , na.rm = TRUE ) ) 
  # (multiplicative) bias
  bias_era5 <- (mean(era5_values, na.rm = TRUE))/(mean(eobs_values, na.rm = TRUE))
  bias_molo <- (mean(molo_values, na.rm = TRUE))/(mean(eobs_values, na.rm = TRUE))
  bias_bola <- (mean(bola_values, na.rm = TRUE))/(mean(eobs_values, na.rm = TRUE))
  # ME, MSE, MAE
  verify_era5 <- verify(eobs_values,era5_values, frcst.type = "cont", obs.type = "cont",)
  verify_molo <- verify(eobs_values,molo_values, frcst.type = "cont", obs.type = "cont",)
  verify_bola <- verify(eobs_values,bola_values, frcst.type = "cont", obs.type = "cont",)
  # Write to output file
  df <- data.frame(c(my_year,cor_era5,rmse_era5,bias_era5,verify_era5$ME,verify_era5$MAE), c(my_year,cor_molo,rmse_molo,bias_molo,verify_molo$ME,verify_molo$MAE), c(my_year,cor_bola,rmse_bola,bias_bola,verify_bola$ME,verify_bola$MAE))
  if (my_year == ini_year) {
    write.table(df, file_scores_csv, col.names = c("ERA5-Land","MOLOCH","BOLAM"), row.names = c("Year","Correlation","RMSE","Bias","Mean Error","Mean Absolute Error"), sep=";", append = FALSE)
  } else {
    write.table(df, file_scores_csv, col.names = c("ERA5-Land","MOLOCH","BOLAM"), row.names = c("Year","Correlation","RMSE","Bias","Mean Error","Mean Absolute Error"), sep=";", append = TRUE)
  }
  print(paste("VER_CORR",my_year,sprintf("%.2f",cor_era5),sprintf("%.2f",cor_molo),sprintf("%.2f",cor_bola),sep=";"))
  print(paste("VER_RMSE",my_year,sprintf("%.2f",rmse_era5),sprintf("%.2f",rmse_molo),sprintf("%.2f",rmse_bola),sep=";"))
  print(paste("VER_BIAS",my_year,sprintf("%.2f",bias_era5),sprintf("%.2f",bias_molo),sprintf("%.2f",bias_bola),sep=";"))
  print(paste("VER_ME",my_year,sprintf("%.2f",verify_era5$ME),sprintf("%.2f",verify_molo$ME),sprintf("%.2f",verify_bola$ME),sep=";"))

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

