# Load libraries
library(raster)
library(ncdf4)
library(rgdal)
library(RColorBrewer)
library(viridis)
library(elevatr)
# add here any other packages you need

# Global variables
u               <- 'capecchi'
my_working_dir  <- paste0("/home/",u,"/art_hindcast_meteo/DATA")
my_figs_dir     <- paste0("/home/",u,"/art_hindcast_meteo/figs")
my_output_dir   <- paste0("/home/",u,"/art_hindcast_meteo/res")
my_shp_dir      <- paste0("/home/",u,"/Sync/shp")
topo_cutoff     <- 1
sf_era5         <- 1000
my_pal          <- brewer.pal(n = 20, name = "RdBu")
my_limits       <- c(300,3000)
my_limits_perc  <- c(-100,100)

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
my_coastline_lr <- readOGR(
  dsn= paste0(my_shp_dir,"/gshhg-shp-2.3.7/GSHHS_shp/l/") ,
  layer="GSHHS_l_L1",
  verbose=FALSE)
add_shp=function(){plot(my_coastline_lr, bg="transparent", add=TRUE,lwd=0.5)}

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
# Get extent on the basis of the smallest dataset (ARCIS)
my_extent <- extent(obse_orig)
# Crop ERA5-Land grid on my_extent and regrid (if needed)
era5_crop           <- crop(era5_orig*sf_era5, my_extent)
era5_resam          <- resample(era5_crop, era5_crop, method="bilinear")
era5_mask           <- mask(x = era5_resam, mask = my_italy)
# Crop E-OBS grid on my_extent and regrid (if needed)
print(paste0("+++Cropping rasters..."))
obse_crop           <- obse_orig
obse_resam          <- resample(obse_crop, era5_crop, method="bilinear")
obse_mask           <- mask(x = obse_resam, mask = my_italy)
obse_mask2           <- mask(x = obse_crop, mask = my_italy)
# Crop BOLAM grid on my_extent and regrid (if needed)
bola_crop           <- crop(bola_orig, my_extent)
bola_resam          <- resample(bola_crop, era5_crop, method="bilinear")
bola_mask           <- mask(x = bola_resam, mask = my_italy)
# Crop MOLOCH grid on my_extent and regrid (if needed)
molo_crop           <- crop(molo_orig, my_extent)
molo_resam          <- resample(molo_crop, era5_crop, method="bilinear")
molo_mask           <- mask(x = molo_resam, mask = my_italy)
print(paste0("+++OK"))
###

###
# Calculate bias only for (elevation gt topo_cutoff m)
ita_adm0         <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/ITA_adm")),layer="ITA_adm0")
extent_ita       <- ita_adm0[ita_adm0$NAME_ENGLI=="Italy",]
elevation        <- get_elev_raster(extent_ita, z = 4)
elevation_crop   <- crop(elevation,my_extent)
elev_cutoff      <- elevation_crop > topo_cutoff
elev_cutoff_resam <- resample(elev_cutoff, era5_crop, method="bilinear")
# ERA5
r <- (((era5_resam-obse_resam)/obse_resam)*elev_cutoff_resam)*100
m_era5 <- sprintf("%.1f",cellStats(r, 'mean'))
r <- NULL
# MOLOCH
r <- (((molo_resam-obse_resam)/obse_resam)*elev_cutoff_resam)*100
m_molo <- sprintf("%.1f",cellStats(r, 'mean'))
r <- NULL
# BOLAM
r <- (((bola_resam-obse_resam)/obse_resam)*elev_cutoff_resam)*100
m_bola <- sprintf("%.1f",cellStats(r, 'mean'))
r <- NULL

# Now plot
my_eps <- paste0(my_figs_dir,"/tp_YearAno_vs_",ref_dataset,".eps")
setEPS()  
postscript(my_eps)
par(mfrow = c(2, 2))
# ARCIS
plot(obse_mask2, main="(a) ARCIS [mm/year]", col = rev(viridis(20)), zlim=my_limits, axes=FALSE, legend.width=3, legend=T, box=F, addfun=add_shp)
# ERA5-Land
plot(((era5_mask-obse_mask)/obse_mask)*100, main="(b) ERA5-Land bias [%]", col = my_pal, zlim=my_limits_perc, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
legend("topleft", legend=paste0("bias=",m_era5,"%"),cex=0.8)
# MOLOCH
plot(((molo_mask-obse_mask)/obse_mask)*100, main="(c) MOLOCH bias [%]", col = my_pal, zlim=my_limits_perc, axes=FALSE, addfun=add_shp, legend.width=3, legend=T, box=F)
legend("topleft", legend=paste0("bias=",m_molo,"%"),cex=0.8)
# BOLAM
plot(((bola_mask-obse_mask)/obse_mask)*100, main="(d) BOLAM bias [%]", col = my_pal, zlim=my_limits_perc, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
legend("topleft", legend=paste0("bias=",m_bola,"%"),cex=0.8)
dev.off()

#############################
ref_dataset     <- 'GRIPHO'
#############################
# Average PCP
file_obse   <- paste0("GRIPHO/tp/gripho-v2-3km_year_tp_AVG2001-2016.nc")
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
# Get extent on the basis of the smallest dataset (GRIPHO)
my_extent           <- extent(obse_orig)
my_extent_gripho    <- extent(obse_orig)
print(paste0("+++Cropping rasters..."))
# Crop ERA5-Land grid on my_extent and regrid (if needed)
era5_crop           <- crop(era5_orig*sf_era5, my_extent)
era5_resam          <- era5_crop
era5_mask           <- mask(x = era5_resam, mask = my_italy)
# Crop E-OBS grid on my_extent and regrid (if needed)
obse_crop           <- crop(obse_orig, my_extent)
obse_resam          <- resample(obse_crop, era5_resam, method="bilinear")
obse_mask           <- mask(x = obse_resam, mask = my_italy)
obse_mask2           <- mask(x = obse_crop, mask = my_italy)
# Crop BOLAM grid on my_extent and regrid (if needed)
bola_crop           <- crop(bola_orig, my_extent)
bola_resam          <- resample(bola_crop, era5_resam, method="bilinear")
bola_mask           <- mask(x = bola_resam, mask = my_italy)
# Crop MOLOCH grid on my_extent and regrid (if needed)
molo_crop           <- crop(molo_orig, my_extent)
molo_resam          <- resample(molo_crop, era5_resam, method="bilinear")
molo_mask           <- mask(x = molo_resam, mask = my_italy)
print(paste0("+++OK"))
###

###
# Calculate bias only for (elevation gt topo_cutoff m)
ita_adm0         <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/ITA_adm")),layer="ITA_adm0")
extent_ita       <- ita_adm0[ita_adm0$NAME_ENGLI=="Italy",]
elevation        <- get_elev_raster(extent_ita, z = 4)
elevation_crop   <- crop(elevation,my_extent)
elev_cutoff      <- elevation_crop > topo_cutoff
elev_cutoff_resam <- resample(elev_cutoff, era5_crop, method="bilinear")
# ERA5
r <- (((era5_resam-obse_resam)/obse_resam)*elev_cutoff_resam)*100
m_era5 <- sprintf("%.1f",cellStats(r, 'mean'))
r <- NULL
# MOLOCH
r <- (((molo_resam-obse_resam)/obse_resam)*elev_cutoff_resam)*100
m_molo <- sprintf("%.1f",cellStats(r, 'mean'))
r <- NULL
# BOLAM
r <- (((bola_resam-obse_resam)/obse_resam)*elev_cutoff_resam)*100
m_bola <- sprintf("%.1f",cellStats(r, 'mean'))
r <- NULL

# Now plot
my_eps <- paste0(my_figs_dir,"/tp_YearAno_vs_",ref_dataset,".eps")
setEPS()  
postscript(my_eps)
par(mfrow = c(2, 2))
# E-OBS
plot(obse_mask2, main="(a) GRIPHO [mm/year]", col = rev(viridis(20)), zlim=my_limits, axes=FALSE, legend.width=3, legend=T, box=F, addfun=add_shp)
# ERA5-Land
plot(((era5_mask-obse_mask)/obse_mask)*100, main="(b) ERA5-Land bias [%]", col = my_pal, zlim=my_limits_perc, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
legend("topright", legend=paste0("bias=",m_era5,"%"),cex=0.8)
# MOLOCH
plot(((molo_mask-obse_mask)/obse_mask)*100, main="(c) MOLOCH bias [%]", col = my_pal, zlim=my_limits_perc, axes=FALSE, addfun=add_shp, legend.width=3, legend=T, box=F)
legend("topright", legend=paste0("bias=",m_molo,"%"),cex=0.8)
# BOLAM
plot(((bola_mask-obse_mask)/obse_mask)*100, main="(d) BOLAM bias [%]", col = my_pal, zlim=my_limits_perc, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
legend("topright", legend=paste0("bias=",m_bola,"%"),cex=0.8)
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
# Get extent on the basis of the smallest dataset (MOLOCH)
my_extent           <- my_extent_gripho
# Crop E-OBS grid on my_extent and regrid (if needed)
print(paste0("+++Cropping rasters..."))
obse_crop           <- crop(obse_orig, my_extent)
obse_resam          <- obse_crop
obse_mask           <- mask(x = obse_resam, mask = my_italy)
# Crop ERA5-Land grid on my_extent and regrid (if needed)
era5_crop           <- crop(era5_orig*sf_era5, my_extent)
era5_resam          <- resample(era5_crop, obse_crop, method="bilinear")
era5_mask           <- mask(x = era5_resam, mask = my_italy)
# Crop BOLAM grid on my_extent and regrid (if needed)
bola_crop           <- crop(bola_orig, my_extent)
bola_resam          <- resample(bola_crop, obse_crop, method="bilinear")
bola_mask           <- mask(x = bola_resam, mask = my_italy)
# Crop MOLOCH grid on my_extent and regrid (if needed)
molo_crop           <- crop(molo_orig, my_extent)
molo_resam          <- resample(molo_crop, obse_crop, method="bilinear")
molo_mask           <- mask(x = molo_resam, mask = my_italy)
print(paste0("+++OK"))
###

###
# Calculate bias only for (elevation gt topo_cutoff m)
ita_adm0         <- readOGR(dsn=path.expand(paste0(my_shp_dir,"/ITA_adm")),layer="ITA_adm0")
extent_ita       <- ita_adm0[ita_adm0$NAME_ENGLI=="Italy",]
elevation        <- get_elev_raster(extent_ita, z = 4)
elevation_crop   <- crop(elevation,my_extent)
elev_cutoff      <- elevation_crop > topo_cutoff
elev_cutoff_resam <- resample(elev_cutoff, obse_crop, method="bilinear")
# ERA5
r <- (((era5_resam-obse_resam)/obse_resam)*elev_cutoff_resam)*100
m_era5 <- sprintf("%.1f",cellStats(r, 'mean'))
r <- NULL
# MOLOCH
r <- (((molo_resam-obse_resam)/obse_resam)*elev_cutoff_resam)*100
m_molo <- sprintf("%.1f",cellStats(r, 'mean'))
r <- NULL
# BOLAM
r <- (((bola_resam-obse_resam)/obse_resam)*elev_cutoff_resam)*100
m_bola <- sprintf("%.1f",cellStats(r, 'mean'))
r <- NULL


# Now plot
my_eps <- paste0(my_figs_dir,"/tp_YearAno_vs_",ref_dataset,".eps")
setEPS()  
postscript(my_eps)
par(mfrow = c(2, 2))
# E-OBS
plot(obse_mask, main="(a) E-OBS [mm/year]", col = rev(viridis(20)), zlim=my_limits, axes=FALSE, legend.width=3, legend=T, box=F, addfun=add_shp)
# ERA5-Land
plot(((era5_mask-obse_mask)/obse_mask)*100, main="(b) ERA5-Land bias [%]", col = my_pal, zlim=my_limits_perc, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
legend("topright", legend=paste0("bias=",m_era5,"%"),cex=0.8)
# MOLOCH
plot(((molo_mask-obse_mask)/obse_mask)*100, main="(c) MOLOCH bias [%]", col = my_pal, zlim=my_limits_perc, axes=FALSE, addfun=add_shp, legend.width=3, legend=T, box=F)
legend("topright", legend=paste0("bias=",m_molo,"%"),cex=0.8)
# BOLAM
plot(((bola_mask-obse_mask)/obse_mask)*100, main="(d) BOLAM bias [%]", col = my_pal, zlim=my_limits_perc, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
legend("topright", legend=paste0("bias=",m_bola,"%"),cex=0.8)
dev.off()

# Byebye
print(paste0("Byebye"))
# quit()
