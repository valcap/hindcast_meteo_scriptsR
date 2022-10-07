# Load libraries
library(raster)
library(ncdf4)
library(verification)
library(Metrics)
library(plotrix)
library(lubridate)
library(dplyr)
library(rgdal)
library(RColorBrewer)
# add here any other packages you need

# Global variables
my_working_dir  <- paste0("~/art_hindcast_meteo/")
my_figs_dir     <- paste0("~/art_hindcast_meteo/figs")
my_output_dir   <- paste0("~/art_hindcast_meteo/res")
my_shp_dir      <- paste0("~/Sync/shp")
ini_year        <- 1981
end_year        <- 2019
vec_years       <- seq(ini_year,end_year)
sf_era5         <- 1000
my_pal=brewer.pal(n = 6, name = "RdBu")

# Working directory
setwd(my_working_dir)

# Load the point locations...
my_anag_file     <- paste0(my_working_dir,"/DATA/points050_moloch_extent.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))
print(paste0(""))

###
  # Loading raster file names
  print(paste0("+++Reading raster files..."))
  # Raster files definition
  file_eobs   <- paste0("DATA/E-OBS/tp/E-OBS_year_tp_years.nc")
  file_era5   <- paste0("DATA/ERA5-Land/tp/ERA5-Land_year_tp_years.nc")
  file_molo   <- paste0("DATA/MOLOCH/tp/moloch_year_tp_years.nc")
  file_bola   <- paste0("DATA/BOLAM/tp/bolam_year_tp_years.nc") 
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
  print(paste0("+++Cropping rasters..."))
#  my_extent           <- extent(pointCoordinates)
  my_extent           <- extent(molo_orig) 
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

  index_years          <- format(as.Date(names(eobs_crop), format = "X%Y.%m.%d"), format = "%Y")
  index_years_era5     <- format(as.Date(names(era5_crop), format = "X%Y.%m.%d.%H.%M.%S"), format = "%Y")

###
  # Adjusting names
  print(paste0("+++Adjusting name conventions..."))
  names(eobs_resam)     <- as.Date(names(eobs_resam),format="X%Y.%m.%d")
  names(era5_resam)     <- as.Date(names(era5_resam),format="X%Y.%m.%d")
  names(bola_resam)     <- as.Date(names(bola_resam),format="X%Y.%m.%d.%H.%M.%S")
  names(molo_resam)     <- names(bola_resam)
  # select only data from 1981
  print(paste0("+++...and selecting data from 1981..."))
  # E-OBS
  r <- setZ(eobs_resam, as.Date(names(eobs_resam),format="X%Y.%m.%d"))
  eobs_ok <- subset(r, which(getZ(r) > '1981-01-02'))
  # ERA5
  r <- setZ(era5_resam, as.Date(names(era5_resam),format="X%Y.%m.%d"))
  era5_ok <- subset(r, which(getZ(r) > '1981-01-02'))
  # MOLOCH
  r <- setZ(molo_resam, as.Date(names(molo_resam),format="X%Y.%m.%d"))
  molo_ok <- subset(r, which(getZ(r) > '1981-01-02'))
  # BOLAM
  r <- setZ(bola_resam, as.Date(names(bola_resam),format="X%Y.%m.%d"))
  bola_ok <- subset(r, which(getZ(r) > '1981-01-02'))
  print(paste0("+++OK"))
###

###
  # Plotting rasters
  shp_file <- readOGR(paste0(my_shp_dir,"/gshhg-shp-2.3.7/GSHHS_shp/c/GSHHS_c_L1.shp"))
  add_shp=function(){plot(shp_file, bg="transparent", add=TRUE)}
# E-OBS
  eobs_mask <- mask(x = eobs_ok, mask = shp_file)
  # ERA5-Land
  era5_eps <- paste0(my_figs_dir,"/era5_tp-ano_year.eps")
  setEPS()
  postscript(era5_eps, horizontal = FALSE, onefile = FALSE, paper = "special", height = 15, width = 15)
  era5_mask <- mask(x = era5_ok, mask = shp_file)
  plot((((era5_mask*sf_era5)-eobs_mask)/eobs_mask)*100, main=index_years_era5, maxnl=length(index_years_era5), nc=6, nr=7, col = my_pal, zlim=c(-100,100), axes=FALSE, addfun=add_shp, legend.width=3, legend=T, box=F)
  dev.off()
# MOLOCH
  molo_eps <- paste0(my_figs_dir,"/molo_tp-ano_year.eps")
  setEPS()
  postscript(molo_eps, horizontal = FALSE, onefile = FALSE, paper = "special", height = 15, width = 15)
  molo_mask <- mask(x = molo_ok, mask = shp_file)
  plot(((molo_mask-eobs_mask)/eobs_mask)*100, main=index_years_era5, maxnl=length(index_years_era5), nc=6, nr=7, col = my_pal, zlim=c(-100,100), axes=FALSE, addfun=add_shp, legend.width=3, legend=T, box=F)
  dev.off()
# BOLAM
  bola_eps <- paste0(my_figs_dir,"/bola_tp-ano_year.eps")
  setEPS()
  postscript(bola_eps, horizontal = FALSE, onefile = FALSE, paper = "special", height = 15, width = 15)
  bola_mask <- mask(x = bola_ok, mask = shp_file)
  plot(((bola_mask-eobs_mask)/eobs_mask)*100, main=index_years_era5, maxnl=length(index_years_era5), nc=6, nr=7, col = my_pal, zlim=c(-100,100), axes=FALSE, addfun=add_shp, legend.width=3, legend=T, box=F)
  dev.off()
  
# Byebye
print(paste0("Byebye"))
quit()

