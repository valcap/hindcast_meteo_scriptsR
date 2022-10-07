# Purpose
print(paste0("Calculating standard statistical verification scores"))

# Load libraries
library(raster)
library(ncdf4)
library(verification)
library(Metrics)
library(plotrix)
library(lubridate)
library(dplyr)
library(rgdal)

# Global variables
u               <- 'valcap'
my_working_dir  <- paste0("/home/",u,"/art_hindcast_meteo/")
my_figs_dir     <- paste0("/home/",u,"/art_hindcast_meteo/figs")
my_output_dir   <- paste0("/home/",u,"/art_hindcast_meteo/res")
my_shp_dir      <- paste0("/home/",u,"/Sync/shp")
ini_year        <- 1981
end_year        <- 2019
vec_years       <- seq(ini_year,end_year)
sf_era5         <- 1000
cutoff_NA       <- 50
my_rmse <- function(x) {  
  sqrt( mean( (x)^2 , na.rm = TRUE ) )
}
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

# Working directory
setwd(my_working_dir)

######################
print(paste("VER_SCORE","YEAR","ERA5-Land","MOLOCH","BOLAM",sep=";"),quote=FALSE)
# MAIN LOOP OVER YEARS
for (my_year in vec_years) {
  print(paste0("___Running year ",my_year,"___"))

###
  # Loading observational dataset
  my_anag_file     <- paste0(my_working_dir,"/DATA/SCIA/SCIA_TotAnnPre_",my_year,".csv")
  anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
  pointCoordinates <- anag_full
  coordinates(pointCoordinates)= ~ lon+lat
  print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))
  print(paste0(""))
  # ...and plot them on a map
  setEPS()
  postscript(paste0(my_figs_dir,"/SCIA_rain-gauges_",my_year,".eps"))
  plot(extent(my_italy), main=paste0(my_year," - Number of rain-gauges ",length(pointCoordinates)), add=F, col='white', lwd=1, xlab="Longitude [degree]", ylab="Latitude [degree]")
  plot(pointCoordinates, pch=1, cex=0.5, add=T)
  plot(my_coastline, lwd=1.5, border="black",add=TRUE)
  dev.off()

###
  # Loading raster file names
  print(paste0("+++Reading raster files..."))
  # Raster files definition
  file_era5   <- paste0("DATA/ERA5-Land/tp/ERA5-Land_year_tp_",my_year,".nc")
  file_molo   <- paste0("DATA/MOLOCH/tp/moloch_year_tp_",my_year,".nc")
  file_bola   <- paste0("DATA/BOLAM/tp/bolam_year_tp_",my_year,".nc") 
  # Allocate rasters
  era5_orig   <- NULL
  molo_orig   <- NULL
  bola_orig   <- NULL
  era5_orig   <- brick(file_era5)
  molo_orig   <- brick(file_molo)
  bola_orig   <- brick(file_bola)
  print(paste0("+++OK"))
###

###
  # Get extent on the basis of the point coordinates
  print(paste0("+++Cropping extents and resampling rasters..."))
  my_extent           <- extent(molo_orig)
  # Crop ERA5-Land grid on my_extent
  era5_crop           <- crop(era5_orig, my_extent)
  era5_resam          <- era5_crop*sf_era5
  # Crop BOLAM grid on my_extent and regrid
  bola_crop           <- crop(bola_orig, my_extent)
  bola_resam          <- resample(bola_crop, era5_crop, method="bilinear")
  # Crop MOLOCH grid on my_extent and regrid
  molo_crop           <- crop(molo_orig, my_extent)
  molo_resam          <- resample(molo_crop, era5_crop, method="bilinear")
  print(paste0("+++OK"))
###

###
  # select only data from the current year
  print(paste0("+++Adjusting name conventions..."))
  names(era5_resam)     <- as.Date(names(era5_resam),format="X%Y.%m.%d")
  names(bola_resam)     <- as.Date(names(bola_resam),format="X%Y.%m.%d.%H.%M.%S")
  names(molo_resam)     <- names(bola_resam)
  print(paste0("+++OK"))
###

###
  # Extracting data
  print(paste0("+++Extracting data..."))
  era5_values <- as.data.frame(extract(era5_resam, pointCoordinates))[,1]
  molo_values <- as.data.frame(extract(molo_resam, pointCoordinates))[,1]
  bola_values <- as.data.frame(extract(bola_resam, pointCoordinates))[,1]
  #
  era5_values <- replace(era5_values,era5_values<cutoff_NA, NA)
  molo_values <- replace(molo_values,molo_values<cutoff_NA, NA)
  bola_values <- replace(bola_values,bola_values<cutoff_NA, NA)
  print(paste0("+++OK"))
  # Observed values as data frame
  obse_values <- pointCoordinates$val
###

###
  # Standard statistical verification scores for continuos variables
  print(paste0("+++Calculating standard statistical verification scores..."))
  #ERA5
  corr_era5      <- cor(era5_values,obse_values, use = "complete.obs")
  rmse_era5      <- my_rmse(era5_values-obse_values)
  mbias_era5     <- mean(era5_values,na.rm=T)/mean(obse_values,na.rm=T)
  obias_era5     <- mean(era5_values,na.rm=T)-mean(obse_values,na.rm=T)
  me_era5        <- mean(era5_values-obse_values,na.rm=T)
  mae_era5       <- mean(abs(era5_values-obse_values),na.rm=T)
  #MOLOCH
  corr_molo      <- cor(molo_values,obse_values, use = "complete.obs")
  rmse_molo      <- my_rmse(molo_values-obse_values)
  mbias_molo     <- mean(molo_values,na.rm=T)/mean(obse_values,na.rm=T)
  obias_molo     <- mean(molo_values,na.rm=T)-mean(obse_values,na.rm=T)
  me_molo        <- mean(molo_values-obse_values,na.rm=T)
  mae_molo       <- mean(abs(molo_values-obse_values),na.rm=T)
  #BOLAM
  corr_bola      <- cor(bola_values,obse_values, use = "complete.obs")
  rmse_bola      <- my_rmse(bola_values-obse_values)
  mbias_bola     <- mean(bola_values,na.rm=T)/mean(obse_values,na.rm=T)
  obias_bola     <- mean(bola_values,na.rm=T)-mean(obse_values,na.rm=T)
  me_bola        <- mean(bola_values-obse_values,na.rm=T)
  mae_bola       <- mean(abs(bola_values-obse_values),na.rm=T)
  
  print(paste("VER_CORR",my_year,sprintf("%.2f",corr_era5),sprintf("%.2f",corr_molo),sprintf("%.2f",corr_bola),sep=";"),quote=FALSE)
  print(paste("VER_RMSE",my_year,sprintf("%.2f",rmse_era5),sprintf("%.2f",rmse_molo),sprintf("%.2f",rmse_bola),sep=";"),quote=FALSE)
  print(paste("VER_MBIAS",my_year,sprintf("%.2f",mbias_era5),sprintf("%.2f",mbias_molo),sprintf("%.2f",mbias_bola),sep=";"),quote=FALSE)
  print(paste("VER_OBIAS",my_year,sprintf("%.2f",obias_era5),sprintf("%.2f",obias_molo),sprintf("%.2f",obias_bola),sep=";"),quote=FALSE)
  print(paste("VER_ME",my_year,sprintf("%.2f",me_era5),sprintf("%.2f",me_molo),sprintf("%.2f",me_bola),sep=";"),quote=FALSE)
  print(paste("VER_MAE",my_year,sprintf("%.2f",mae_era5),sprintf("%.2f",mae_molo),sprintf("%.2f",mae_bola),sep=";"),quote=FALSE)
  
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
# quit()
