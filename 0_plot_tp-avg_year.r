# Load libraries
library(raster)
library(ncdf4)
library(rgdal)
library(RColorBrewer)
library(viridis)
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

# Working directory
setwd(my_working_dir)

### 
# Adding shapefile
my_italy <- readOGR(
  dsn= paste0(my_shp_dir,"/") ,
  layer="ITA_adm0",
  verbose=FALSE)
my_coastline <- readOGR( 
  dsn= paste0(my_shp_dir,"/ne_10m_coastline/") , 
  layer="ne_10m_coastline",
  verbose=FALSE)
my_europe_adm0 <- readOGR( 
  dsn= paste0(my_shp_dir,"/gadm_v1_lev0_shp/") , 
  layer="europe_adm1_lev0",
  verbose=FALSE)
my_coastline_lr <- readOGR(
  dsn= paste0(my_shp_dir,"/gshhg-shp-2.3.7/GSHHS_shp/l/") ,
  layer="GSHHS_l_L1",
  verbose=FALSE)
add_shp=function(){plot(my_coastline_lr, bg="transparent", add=TRUE, lwd=0.5)}

#############################
ref_dataset     <- 'ARCIS'
#############################
# Average PCP
file_obse   <- paste0("ARCIS/tp/ARCIS_year_tp_AVG1981-2015.nc")
file_era5   <- paste0("ERA5-Land/tp/ERA5-Land_year_tp_AVG1981-2015.nc")
file_molo   <- paste0("MOLOCH/tp/moloch_year_tp_AVG1981-2015.nc")
file_bola   <- paste0("BOLAM/tp/bolam_year_tp_AVG1981-2015.nc") 
# Allocate rasters
obse_orig   <- NULL
era5_orig   <- NULL
molo_orig   <- NULL
bola_orig   <- NULL
obse_orig   <- brick(file_obse)
era5_orig   <- brick(file_era5)
molo_orig   <- brick(file_molo)
bola_orig   <- brick(file_bola)

###
# Get extent on the basis of the point coordinates
my_extent <- extent(obse_orig)
# Crop E-OBS grid on my_extent and regrid (if needed)
print(paste0("+++Cropping rasters..."))
obse_crop           <- obse_orig
obse_resam          <- obse_crop
# Crop ERA5-Land grid on my_extent and regrid (if needed)
era5_crop           <- crop(era5_orig*sf_era5, my_extent)
era5_resam          <- era5_crop
# Crop BOLAM grid on my_extent and regrid (if needed)
bola_crop           <- crop(bola_orig, my_extent)
bola_resam          <- bola_crop
# Crop MOLOCH grid on my_extent and regrid (if needed)
molo_crop           <- crop(molo_orig, my_extent)
molo_resam          <- molo_crop
print(paste0("+++OK"))
###

# Plotting rasters breaks = c(5,50,100,200,300,400,500,750,1000,1500,2000,3000)
my_eps <- paste0(my_figs_dir,"/tp_YearMean_vs_",ref_dataset,".eps")
setEPS()  
postscript(my_eps)
par(mfrow = c(2, 2))
# ARCIS
plot(obse_crop, main="ARCIS", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, legend.width=3, legend=T, box=F, addfun=add_shp)
# ERA5-Land
plot(era5_crop, main="ERA5-Land", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, legend.width=3, legend=F, box=F, addfun=add_shp)
# MOLOCH
plot(molo_crop, main="MOLOCH", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, legend.width=3, legend=F, box=F, addfun=add_shp)
# BOLAM
plot(bola_crop, main="BOLAM", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, legend.width=3, legend=F, box=F, add=F, addfun=add_shp)
mtext(paste0("1981-2015"), # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

#############################
ref_dataset     <- 'GRIPHO'
#############################
# Average PCP
file_obse   <- paste0("GRIPHO/tp/gripho-v2-3km_year_tp_AVG2001-2016.nc")
file_era5   <- paste0("ERA5-Land/tp/ERA5-Land_year_tp_AVG2001-2016.nc")
file_molo   <- paste0("MOLOCH/tp/moloch_year_tp_AVG2001-2016.nc")
file_bola   <- paste0("BOLAM/tp/bolam_year_tp_AVG2001-2016.nc") 
# Allocate rasters
obse_orig   <- NULL
era5_orig   <- NULL
molo_orig   <- NULL
bola_orig   <- NULL
obse_orig   <- brick(file_obse)
era5_orig   <- brick(file_era5)
molo_orig   <- brick(file_molo)
bola_orig   <- brick(file_bola)

###
# Get extent on the basis of the point coordinates
my_extent <- extent(obse_orig)
my_extent_gripho <- my_extent
# Crop E-OBS grid on my_extent and regrid (if needed)
print(paste0("+++Cropping rasters..."))
obse_crop           <- obse_orig
obse_resam          <- obse_crop
# Crop ERA5-Land grid on my_extent and regrid (if needed)
era5_crop           <- crop(era5_orig*sf_era5, my_extent)
era5_resam          <- era5_crop
# Crop BOLAM grid on my_extent and regrid (if needed)
bola_crop           <- crop(bola_orig, my_extent)
bola_resam          <- bola_crop
# Crop MOLOCH grid on my_extent and regrid (if needed)
molo_crop           <- crop(molo_orig, my_extent)
molo_resam          <- molo_crop
print(paste0("+++OK"))
###

# Plotting rasters breaks = c(5,50,100,200,300,400,500,750,1000,1500,2000,3000)
my_eps <- paste0(my_figs_dir,"/tp_YearMean_vs_",ref_dataset,".eps")
setEPS()  
postscript(my_eps)
par(mfrow = c(2, 2))
# E-OBS
plot(obse_crop, main="GRIPHO", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, legend.width=3, legend=T, box=F, addfun=add_shp)
# ERA5-Land
plot(era5_crop, main="ERA5-Land", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, legend.width=3, legend=F, box=F, addfun=add_shp)
# MOLOCH
plot(molo_crop, main="MOLOCH", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, legend.width=3, legend=F, box=F, addfun=add_shp)
# BOLAM
plot(bola_crop, main="BOLAM", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, legend.width=3, legend=F, box=F, add=F, addfun=add_shp)
mtext(paste0("2001-2016"), # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

#############################
ref_dataset     <- 'E-OBS'
#############################
# Average PCP
file_obse   <- paste0("E-OBS/tp/E-OBS_year_tp_AVG1981-2019.nc")
file_era5   <- paste0("ERA5-Land/tp/ERA5-Land_year_tp_AVG1981-2019.nc")
file_molo   <- paste0("MOLOCH/tp/moloch_year_tp_AVG1981-2019.nc")
file_bola   <- paste0("BOLAM/tp/bolam_year_tp_AVG1981-2019.nc") 
# Allocate rasters
obse_orig   <- NULL
era5_orig   <- NULL
molo_orig   <- NULL
bola_orig   <- NULL
obse_orig   <- brick(file_obse)
era5_orig   <- brick(file_era5)
molo_orig   <- brick(file_molo)
bola_orig   <- brick(file_bola)

###
# Get extent on the basis of the point coordinates
my_extent <- my_extent_gripho
# Crop E-OBS grid on my_extent and regrid (if needed)
print(paste0("+++Cropping rasters..."))
obse_crop           <- crop(obse_orig, my_extent)
obse_resam          <- obse_crop
# Crop ERA5-Land grid on my_extent and regrid (if needed)
era5_crop           <- crop(era5_orig*sf_era5, my_extent)
era5_resam          <- era5_crop
# Crop BOLAM grid on my_extent and regrid (if needed)
bola_crop           <- crop(bola_orig, my_extent)
bola_resam          <- bola_crop
# Crop MOLOCH grid on my_extent and regrid (if needed)
molo_crop           <- crop(molo_orig, my_extent)
molo_resam          <- molo_crop
print(paste0("+++OK"))
###

# Plotting rasters breaks = c(5,50,100,200,300,400,500,750,1000,1500,2000,3000)
my_eps <- paste0(my_figs_dir,"/tp_YearMean_vs_",ref_dataset,".eps")
setEPS()  
postscript(my_eps)
par(mfrow = c(2, 2))
# E-OBS
plot(obse_crop, main="E-OBS", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, legend.width=3, legend=T, box=F, addfun=add_shp)
# ERA5-Land
plot(era5_crop, main="ERA5-Land", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, legend.width=3, legend=F, box=F, addfun=add_shp)
# MOLOCH
plot(molo_crop, main="MOLOCH", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, legend.width=3, legend=F, box=F, addfun=add_shp)
# BOLAM
plot(bola_crop, main="BOLAM", col = rev(viridis(100)),zlim=my_limits, axes=FALSE, legend.width=3, legend=F, box=F, add=F, addfun=add_shp)
mtext(paste0("1981-2019"), # Add main title
      side = 3,
      line = - 2,
      outer = TRUE)
dev.off()

# Byebye
print(paste0("Byebye"))
# quit()
