# Purpose
print(paste0("Purpose:...TODO..."))

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
TD_normal       <- FALSE
ini_year        <- 1981
end_year        <- 2019
vec_years       <- seq(ini_year,end_year)
sf_era5         <- 100
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
# ...and plot them on a map
setEPS()
postscript(paste0(my_figs_dir,"/virtual_rain-gauges_EOBS.eps"))
plot(pointCoordinates, main=paste0("# of (virtual) rain-gauges ",length(pointCoordinates)))
plot(extent(pointCoordinates), add=TRUE, col='red', lwd=1)
dev.off()

######################
# MAIN LOOP OVER YEARS
for (my_year in vec_years) {
  print(paste0("___Running year ",my_year,"___"))

###
  # Loading raster file names
  print(paste0("+++Reading raster files..."))
  # Raster files definition
  file_eobs   <- paste0("E-OBS/tp/E-OBS_day_tp_",my_year,".nc")
  file_era5   <- paste0("ERA5-Land/tp/ERA5-Land_",my_year,"-daysum_tp.nc")
  file_molo   <- paste0("work/moloch/moloch_",my_year,"_daysum_apcp.nc")
  file_bola   <- paste0("work/bolam/bolam_",my_year,"_daysum_apcp.nc") 
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
  # select only data from the current year
  print(paste0("+++Adjusting name conventions..."))
  names(eobs_resam)     <- as.Date(names(eobs_resam),format="X%Y.%m.%d")
  if ((my_year == '1990') | (my_year == '2000') | (my_year == '2010')) { 
    names(era5_resam)   <- as.Date(names(era5_resam),format="X%Y%m%d.4791667")
  } else {
    names(era5_resam)   <- as.Date(names(era5_resam),format="X%Y.%m.%d.%H.%M.%S")
  }
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
  # Creating yearly data
  print(paste0("+++Creating yearly data"))
  index_years          <- format(as.Date(names(molo_ok), format = "X%Y.%m.%d"), format = "%Y")
  eobs_ok_year         <- stackApply(eobs_ok, index_years, fun = sum)
  names(eobs_ok_year)  <- my_year
  era5_ok_year         <- stackApply(era5_ok, index_years, fun = sum)
  names(era5_ok_year)  <- my_year
  molo_ok_year         <- stackApply(molo_ok, index_years, fun = sum)
  names(molo_ok_year)  <- my_year
  bola_ok_year         <- stackApply(bola_ok, index_years, fun = sum)
  names(bola_ok_year)  <- my_year
  print(paste0("+++OK"))
###

###
  # Extracting data
  print(paste0("+++Extracting data..."))
  eobs_values <- as.data.frame(extract(eobs_ok_year, pointCoordinates))
  era5_values <- as.data.frame(extract(era5_ok_year*sf_era5, pointCoordinates))
  molo_values <- as.data.frame(extract(molo_ok_year, pointCoordinates))
  bola_values <- as.data.frame(extract(bola_ok_year, pointCoordinates))
  eobs_values_ts <- c(eobs_values_ts,eobs_values[,1])
  era5_values_ts <- c(era5_values_ts,era5_values[,1])
  molo_values_ts <- c(molo_values_ts,molo_values[,1])
  bola_values_ts <- c(bola_values_ts,bola_values[,1])
  print(paste0("+++OK"))
###

###
  # End of my_year loop
  print(paste0("___End of ",my_year,"___"))
}
# END OF MAIN LOOP OVER YEARS
#############################


eobs_values_ts <- replace(eobs_values_ts,eobs_values_ts<cutoff_NA, NA)
era5_values_ts <- replace(era5_values_ts,era5_values_ts<cutoff_NA, NA)
molo_values_ts <- replace(molo_values_ts,molo_values_ts<cutoff_NA, NA)
bola_values_ts <- replace(bola_values_ts,bola_values_ts<cutoff_NA, NA)

# Standard statistical verification scores for continuos variables
print(paste0("+++Calculating standard statistical verification scores..."))
# Correlation
cor_era5 <- cor(era5_values_ts, eobs_values_ts, use = "complete.obs")
cor_molo <- cor(molo_values_ts, eobs_values_ts, use = "complete.obs")
cor_bola <- cor(bola_values_ts, eobs_values_ts, use = "complete.obs")
# RMSE
rmse_era5 <- sqrt( mean( (era5_values_ts-eobs_values_ts)^2 , na.rm = TRUE ) ) 
rmse_molo <- sqrt( mean( (molo_values_ts-eobs_values_ts)^2 , na.rm = TRUE ) ) 
rmse_bola <- sqrt( mean( (bola_values_ts-eobs_values_ts)^2 , na.rm = TRUE ) ) 
# (multiplicative) bias
bias_era5 <- (mean(era5_values_ts, na.rm = TRUE))/(mean(eobs_values_ts, na.rm = TRUE))
bias_molo <- (mean(molo_values_ts, na.rm = TRUE))/(mean(eobs_values_ts, na.rm = TRUE))
bias_bola <- (mean(bola_values_ts, na.rm = TRUE))/(mean(eobs_values_ts, na.rm = TRUE))
# ME, MSE, MAE
verify_era5 <- verify(eobs_values_ts,era5_values_ts, frcst.type = "cont", obs.type = "cont",)
verify_molo <- verify(eobs_values_ts,molo_values_ts, frcst.type = "cont", obs.type = "cont",)
verify_bola <- verify(eobs_values_ts,bola_values_ts, frcst.type = "cont", obs.type = "cont",)
# Write to output file
file_scores_csv <- paste0(my_output_dir,"/tp_scores_vs_EOBS.csv")
df <- data.frame(c(cor_era5,rmse_era5,bias_era5,verify_era5$ME,verify_era5$MAE), c(cor_molo,rmse_molo,bias_molo,verify_molo$ME,verify_molo$MAE), c(cor_bola,rmse_bola,bias_bola,verify_bola$ME,verify_bola$MAE))
write.table(df, file_scores_csv, col.names = c("ERA5-Land","MOLOCH","BOLAM"), row.names = c("Correlation","RMSE","Bias","Mean Error","Mean Absolute Error"), sep=";")
print(paste0("+++OK"))

# Byebye
print(paste0("Byebye"))
quit()
