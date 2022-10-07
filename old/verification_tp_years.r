# Load libraries
library(raster)
library(ncdf4)
library(verification)
library(Metrics)
library(plotrix)
library(lubridate)
library(dplyr)

# Global variables
u               <- 'capecchi'
my_working_dir  <- paste0("/home/",u,"/art_hindcast_meteo/")
my_figs_dir     <- paste0(my_working_dir,"/figs")
my_output_dir   <- paste0(my_working_dir,"/res")
my_shp_dir      <- paste0("/home/",u,"Sync/shp")
TD_normal       <- FALSE
ini_year        <- 1981
end_year        <- 2019
vec_years       <- seq(ini_year,end_year)
sf_era5         <- 1000
my_rmse <- function(x) {
  sqrt( mean( (x)^2 , na.rm = TRUE ) )
}
my_me <- function(x) {
  mean( x , na.rm = TRUE )
}
my_mbias <- function(x,y) {
  mean( x , na.rm = TRUE )/mean( y , na.rm = TRUE )
}
my_obias <- function(x,y) {
  mean( x , na.rm = TRUE ) - mean( y , na.rm = TRUE )
}
my_corr <- function(x,y) {
  cor(x,y,use="complete.obs")
}



# Working directory
setwd(my_working_dir)

###
# Load the point locations
my_anag_file     <- paste0(my_working_dir,"/DATA/E-OBS/points050_moloch_extent.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))
print(paste0(""))
###  

###
# Open raster files
file_eobs   <- paste0("DATA/E-OBS/tp/E-OBS_year_tp_years.nc")
file_era5   <- paste0("DATA/ERA5-Land/tp/ERA5-Land_year_tp_years.nc")
file_molo   <- paste0("DATA/MOLOCH/tp/moloch_year_tp_years.nc")
file_bola   <- paste0("DATA/BOLAM/tp/bolam_year_tp_years.nc")
eobs_orig   <- brick(file_eobs)
era5_orig   <- brick(file_era5)
molo_orig   <- brick(file_molo)
bola_orig   <- brick(file_bola)
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
era5_resam          <- resample(era5_crop, eobs_crop, method="bilinear")
# Crop BOLAM grid on my_extent and regrid (if needed)
bola_crop           <- crop(bola_orig, my_extent)
bola_resam          <- resample(bola_crop, eobs_crop, method="bilinear")
# Crop MOLOCH grid on my_extent and regrid (if needed)
molo_crop           <- crop(molo_orig, my_extent)
molo_resam          <- resample(molo_crop, eobs_crop, method="bilinear")
###

###
print(paste0("+++Adjusting name conventions..."))
names(eobs_resam)     <- as.Date(names(eobs_resam),format="X%Y.%m.%d")
names(era5_resam)     <- as.Date(names(era5_resam),format="X%Y.%m.%d.%H.%M.%S")
names(bola_resam)     <- as.Date(names(bola_resam),format="X%Y.%m.%d.%H.%M.%S")
names(molo_resam)     <- names(bola_resam)
# select only data from 1981 to 2019
print(paste0("+++selecting data from 1981 to 2019..."))
# E-OBS
r <- setZ(eobs_resam, as.Date(names(eobs_resam),format="X%Y.%m.%d"))
eobs_ok <- subset(r, which(getZ(r) >= paste0(ini_year,'-01-02') & (getZ(r) <= paste0(end_year,'-12-31'))))
index_years <- format(as.Date(names(eobs_ok), format = "X%Y.%m.%d"), format = "%Y")
names(eobs_ok)<-index_years
# ERA5
r <- setZ(era5_resam, as.Date(names(era5_resam),format="X%Y.%m.%d"))
era5_ok <- subset(r, which(getZ(r) >= paste0(ini_year,'-01-02') & (getZ(r) <= paste0(end_year,'-12-31'))))
names(era5_ok)<-index_years
# MOLOCH
r <- setZ(molo_resam, as.Date(names(molo_resam),format="X%Y.%m.%d"))
molo_ok <- subset(r, which(getZ(r) >= paste0(ini_year,'-01-02') & (getZ(r) <= paste0(end_year,'-12-31'))))
names(molo_ok)<-index_years
# BOLAM
r <- setZ(bola_resam, as.Date(names(bola_resam),format="X%Y.%m.%d"))
bola_ok <- subset(r, which(getZ(r) >= paste0(ini_year,'-01-02') & (getZ(r) <= paste0(end_year,'-12-31'))))
names(bola_ok)<-index_years
print(paste0("+++OK"))
###

###
# Extracting data
print(paste0("+++Extracting data..."))
eobs_values <- as.data.frame(extract(eobs_ok, pointCoordinates))
eobs_values_mean <- mean(unlist(eobs_values),na.rm=T)
era5_values <- as.data.frame(extract(era5_ok*sf_era5, pointCoordinates))
molo_values <- as.data.frame(extract(molo_ok, pointCoordinates))
bola_values <- as.data.frame(extract(bola_ok, pointCoordinates))
###

###
era5_rmse_years  <- apply(eobs_values-era5_values,2,my_rmse)
molo_rmse_years  <- apply(eobs_values-molo_values,2,my_rmse)
bola_rmse_years  <- apply(eobs_values-bola_values,2,my_rmse)
era5_me_years    <- apply(era5_values-eobs_values,2,my_me)
molo_me_years    <- apply(molo_values-eobs_values,2,my_me)
bola_me_years    <- apply(bola_values-eobs_values,2,my_me)
era5_mbias_years <- apply(era5_values/eobs_values,2,my_me)
molo_mbias_years <- apply(molo_values/eobs_values,2,my_me)
bola_mbias_years <- apply(bola_values/eobs_values,2,my_me)
###
  # GIVE IT A TRY WITH THE DISO INDEX...
eobs_mean_years  <- apply(eobs_values,2,my_me)
era5_corr_year   <- apply(as.data.frame(cor(era5_values,eobs_values,use="complete.obs")),2,my_me)
molo_corr_year   <- apply(as.data.frame(cor(molo_values,eobs_values,use="complete.obs")),2,my_me)
bola_corr_year   <- apply(as.data.frame(cor(bola_values,eobs_values,use="complete.obs")),2,my_me)

era5_diso_years  <- sqrt((era5_corr_year)^2+(era5_me_years/eobs_mean_years)^2+(era5_rmse_years/eobs_mean_years)^2)
molo_diso_years  <- sqrt((molo_corr_year)^2+(molo_me_years/eobs_mean_years)^2+(molo_rmse_years/eobs_mean_years)^2)
bola_diso_years  <- sqrt((bola_corr_year)^2+(bola_me_years/eobs_mean_years)^2+(bola_rmse_years/eobs_mean_years)^2)
###

print("X_RMSE")
print(paste("X","ERA5-Land","MOLOCH","BOLAM",sep="        "))
print(cbind(as.data.frame(era5_rmse_years),as.data.frame(molo_rmse_years),as.data.frame(bola_rmse_years)))
print("X_ME")
print(paste("X","ERA5-Land","MOLOCH","BOLAM",sep="        "))
print(cbind(as.data.frame(era5_me_years),as.data.frame(molo_me_years),as.data.frame(bola_me_years)))
print("X_MBIAS")
print(paste("X","ERA5-Land","MOLOCH","BOLAM",sep="        "))
print(cbind(as.data.frame(era5_mbias_years),as.data.frame(molo_mbias_years),as.data.frame(bola_mbias_years)))
print("X_DISO")
print(paste("X","ERA5-Land","MOLOCH","BOLAM",sep="        "))
print(cbind(as.data.frame(era5_diso_years),as.data.frame(molo_diso_years),as.data.frame(bola_diso_years)))

print(paste0("OK"))


# lanciare il seguente comando
# Rscript QUESTO_SCRIPT.r | grep X | sed 's/ \{1,\}/;/g' | sed -e "s/X//g"

