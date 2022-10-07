# Load libraries
library(raster)
library(ncdf4)
# library(verification)
# library(Metrics)
# library(plotrix)
# library(lubridate)
# library(dplyr)
library(rgdal)
library(RColorBrewer)
# add here any other packages you need

# Global variables
u               <- 'capecchi'
my_working_dir  <- paste0("/home/",u,"/art_hindcast_meteo/DATA")
my_figs_dir     <- paste0("/home/",u,"/art_hindcast_meteo/figs")
my_output_dir   <- paste0("/home/",u,"/art_hindcast_meteo/res")
my_shp_dir      <- paste0("/home/",u,"/Sync/shp")
ini_year        <- 1981
end_year        <- 2015
vec_years       <- seq(ini_year,end_year)
sf_era5         <- 1000
my_pal          <- brewer.pal(n = 6, name = "RdBu")
my_limits       <- c(300,3000)
#############################
ref_dataset     <- 'ARCIS'
bool_arci       <- TRUE
#############################
bool_eobs       <- !bool_arci

# Working directory
setwd(my_working_dir)

# Adding shapefile
if (bool_arci == TRUE) {
  shp_file <- readOGR(paste0(my_shp_dir,"/gshhg-shp-2.3.7/GSHHS_shp/c/GSHHS_c_L1.shp"))
} else {
  shp_file <- readOGR(paste0(my_shp_dir,"/gshhg-shp-2.3.7/GSHHS_shp/c/GSHHS_c_L1.shp"))
}
add_shp=function(){plot(shp_file, bg="transparent", add=TRUE)}

for (my_year in vec_years) {

###
  # Loading raster file names
  print(paste0("+++Reading raster files..."))
  # Raster files definition
  file_arci   <- paste0("ARCIS/tp/ARCIS_year_tp_",my_year,".nc")
  file_eobs   <- paste0("E-OBS/tp/E-OBS_year_tp_",my_year,".nc")
  file_era5   <- paste0("ERA5-Land/tp/ERA5-Land_year_tp_",my_year,".nc")
  file_molo   <- paste0("MOLOCH/tp/moloch_year_tp_",my_year,".nc")
  file_bola   <- paste0("BOLAM/tp/bolam_year_tp_",my_year,".nc") 
  # Allocate rasters
  arci_orig   <- NULL
  eobs_orig   <- NULL
  era5_orig   <- NULL
  molo_orig   <- NULL
  bola_orig   <- NULL
  arci_orig   <- brick(file_arci)
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
  if (bool_arci == TRUE) {
    my_extent <- extent(arci_orig)
  } else {
    my_extent <- extent(molo_orig)
  }
  # Crop E-OBS grid on my_extent and regrid (if needed)
  eobs_crop           <- crop(eobs_orig, my_extent)
  eobs_resam          <- eobs_crop
  # Crop ERA5-Land grid on my_extent and regrid (if needed)
  era5_crop           <- crop(era5_orig*sf_era5, my_extent)
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
  # Plotting rasters breaks = c(5,50,100,200,300,400,500,750,1000,1500,2000,3000)
  my_eps <- paste0(my_figs_dir,"/cfr_year_tp_",my_year,"_vs_",ref_dataset,".eps")
  setEPS()  
  postscript(my_eps)
  par(mfrow=c(2,2))
  # ARCIS
  if (bool_arci==TRUE) {
    arci_mask <- mask(x = arci_orig, mask = shp_file)
    plot(arci_mask, main="ARCIS", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, addfun=add_shp, legend.width=3, legend=T, box=F)
  }
  # E-OBS
  if (bool_eobs==TRUE) {
    eobs_mask <- mask(x = eobs_crop, mask = shp_file)
    plot(eobs_mask, main="E-OBS", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, addfun=add_shp, legend.width=3, legend=T, box=F)
  }
  # ERA5-Land
  era5_mask <- mask(x = era5_crop, mask = shp_file)
  plot(era5_mask, main="ERA5-Land", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
  # MOLOCH
  molo_mask <- mask(x = molo_crop, mask = shp_file)
  plot(molo_mask, main="MOLOCH", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
  # BOLAM
  bola_mask <- mask(x = bola_crop, mask = shp_file)
  plot(bola_mask, main="BOLAM", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
  mtext(paste0(my_year), # Add main title
        side = 3,
        line = - 2,
        outer = TRUE)
  dev.off()
}

# Now plot average PCP
file_arci   <- paste0("ARCIS/tp/ARCIS_year_tp_AVG.nc")
file_eobs   <- paste0("E-OBS/tp/E-OBS_year_tp_AVG.nc")
file_era5   <- paste0("ERA5-Land/tp/ERA5-Land_year_tp_AVG.nc")
file_molo   <- paste0("MOLOCH/tp/moloch_year_tp_AVG.nc")
file_bola   <- paste0("BOLAM/tp/bolam_year_tp_AVG.nc") 
# Allocate rasters
arci_orig   <- NULL
eobs_orig   <- NULL
era5_orig   <- NULL
molo_orig   <- NULL
bola_orig   <- NULL
arci_orig   <- brick(file_arci)
eobs_orig   <- brick(file_eobs)
era5_orig   <- brick(file_era5)
molo_orig   <- brick(file_molo)
bola_orig   <- brick(file_bola)

###
# Get extent on the basis of the point coordinates
print(paste0("+++Cropping rasters..."))
#  my_extent           <- extent(pointCoordinates)
if (bool_arci == TRUE) {
  my_extent <- extent(arci_orig)
} else {
  my_extent <- extent(molo_orig)
}
# Crop E-OBS grid on my_extent and regrid (if needed)
eobs_crop           <- crop(eobs_orig, my_extent)
eobs_resam          <- eobs_crop
# Crop ERA5-Land grid on my_extent and regrid (if needed)
era5_crop           <- crop(era5_orig*sf_era5, my_extent)
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
# Adding shapefile
if (bool_arci == TRUE) {
  shp_file <- readOGR(paste0(my_shp_dir,"/gshhg-shp-2.3.7/GSHHS_shp/i/GSHHS_i_L1.shp"))
} else {
  shp_file <- readOGR(paste0(my_shp_dir,"/gshhg-shp-2.3.7/GSHHS_shp/i/GSHHS_i_L1.shp"))
}
add_shp=function(){plot(shp_file, bg="transparent", add=TRUE)}

# Plotting rasters breaks = c(5,50,100,200,300,400,500,750,1000,1500,2000,3000)
my_eps <- paste0(my_figs_dir,"/cfr_year_tp_AVG_vs_",ref_dataset,".eps")
setEPS()  
postscript(my_eps)
par(mfrow=c(2,2))
# ARCIS
if (bool_arci==TRUE) {
  arci_mask <- mask(x = arci_orig, mask = shp_file)
  plot(arci_mask, main="ARCIS", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, addfun=add_shp, legend.width=3, legend=T, box=F)
}
# E-OBS
if (bool_eobs==TRUE) {
  eobs_mask <- mask(x = eobs_crop, mask = shp_file)
  plot(eobs_mask, main="E-OBS", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, addfun=add_shp, legend.width=3, legend=T, box=F)
}
# ERA5-Land
era5_mask <- mask(x = era5_crop, mask = shp_file)
plot(era5_mask, main="ERA5-Land", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
# MOLOCH
molo_mask <- mask(x = molo_crop, mask = shp_file)
plot(molo_mask, main="MOLOCH", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
# BOLAM
bola_mask <- mask(x = bola_crop, mask = shp_file)
plot(bola_mask, main="BOLAM", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
mtext(paste0("AVERAGE ANNUAL PRECIPITATION"), # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

# Byebye
print(paste0("Byebye"))
# quit()
