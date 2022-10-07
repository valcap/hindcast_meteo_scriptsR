# Load libraries
library(raster)
library(ncdf4)
library(rgdal)
library(RColorBrewer)
library(viridis)
library(elevatr)
library(rayshader)
# add here any other packages you need

# Global variables
u               <- 'capecchi'
my_working_dir  <- paste0("/home/",u,"/art_hindcast_meteo/DATA")
my_figs_dir     <- paste0("/home/",u,"/art_hindcast_meteo/figs")
my_output_dir   <- paste0("/home/",u,"/art_hindcast_meteo/res")
my_shp_dir      <- paste0("/home/",u,"/Sync/shp")
my_pal          <- brewer.pal(n = 11, name = "RdBu")
yini            <- 2002
yend            <- 2016
my_limits_perc  <- c(-100,100)
my_seas         <- 'DJF'

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
ref_dataset     <- 'GRIPHO'
#############################
for (pctl in c(90,95,98,99)) {
if (pctl==90) {my_max = 10}
if (pctl==95) {my_max = 25}
if (pctl==98) {my_max = 40}
if (pctl==99) {my_max = 50}
my_limits   <- c(1,my_max)
# pctl th percentile of hourly mean precipitation
file_obse   <- paste0("GRIPHO/tp/gripho-v2-3km_wethourpctl",pctl,"_tp_",yini,"-",yend,"_",my_seas,".nc")
file_era5   <- paste0("ERA5-Land/tp/ERA5-Land_wethourpctl",pctl,"_tp_",yini,"-",yend,"_",my_seas,".nc")
file_bola   <- paste0("BOLAM/tp/bolam_wethourpctl",pctl,"_tp_",yini,"-",yend,"_",my_seas,".nc") 
file_molo   <- paste0("MOLOCH/tp/moloch_wethourpctl",pctl,"_tp_",yini,"-",yend,"_",my_seas,".nc")
# Allocate rasters
obse_orig   <- NULL
era5_orig   <- NULL
molo_orig   <- NULL
bola_orig   <- NULL
obse_orig   <- brick(file_obse)
era5_orig   <- brick(file_era5)
bola_orig   <- brick(file_bola)
molo_orig   <- brick(file_molo)

###
# Get extent on the basis of the smallest dataset (ARCIS)
my_extent <- extent(obse_orig)
# Crop ERA5-Land grid on my_extent and regrid (if needed)
era5_crop           <- crop(era5_orig, my_extent)
era5_resam          <- resample(era5_crop, era5_crop, method="bilinear")
era5_mask           <- mask(x = era5_resam, mask = my_italy)
#  
print(paste0("+++Cropping rasters..."))
obse_crop           <- obse_orig
obse_resam          <- resample(obse_crop, era5_crop, method="bilinear")
obse_mask           <- mask(x = obse_crop, mask = my_italy)
# Crop BOLAM grid on my_extent and regrid (if needed)
bola_crop           <- crop(bola_orig, my_extent)
bola_resam          <- resample(bola_crop, era5_crop, method="bilinear")
bola_mask           <- mask(x = bola_crop, mask = my_italy)
# Crop MOLOCH grid on my_extent and regrid (if needed)
molo_crop           <- crop(molo_orig, my_extent)
molo_resam          <- resample(molo_crop, era5_crop, method="bilinear")
molo_mask           <- mask(x = molo_crop, mask = my_italy)
print(paste0("+++OK"))
###

# ERA5
r <- (((era5_resam-obse_resam)/obse_resam))*100
m_era5 <- sprintf("%.1f",cellStats(r, 'mean'))
max_era5 <- sprintf("%.1f",cellStats(r, 'max'))
r <- NULL
# MOLOCH
r <- (((molo_resam-obse_resam)/obse_resam))*100
m_molo <- sprintf("%.1f",cellStats(r, 'mean'))
max_molo <- sprintf("%.1f",cellStats(r, 'max'))
r <- NULL
# BOLAM
r <- (((bola_resam-obse_resam)/obse_resam))*100
m_bola <- sprintf("%.1f",cellStats(r, 'mean'))
max_bola <- sprintf("%.1f",cellStats(r, 'max'))
r <- NULL
#print(paste0(max_era5," ",max_molo," ",max_bola))

###
# Now plot
my_eps <- paste0(my_figs_dir,"/tp_wethourpctl",pctl,"_vs_",ref_dataset,"_",my_seas,".eps")
setEPS()  
postscript(my_eps)
par(mfrow = c(2, 2))
# GRIPHO
plot(obse_mask, main=paste0("(a) GRIPHO [mm/1-hour]\n",my_seas), col = rev(viridis(100)), zlim=my_limits, axes=FALSE, legend.width=3, legend=T, box=F, addfun=add_shp)
# ERA5-Land
plot(era5_mask, main="(b) ERA5-Land", col = rev(viridis(100)), zlim=my_limits, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
legend("topright", legend=paste0("bias=",m_era5,"%"),cex=0.8)
# MOLOCH
r <- ((molo_resam-obse_resam)/obse_resam)*100
m <- sprintf("%.1f",cellStats(r, 'mean'))
plot(molo_mask, main="(c) MOLOCH", col = rev(viridis(100)), zlim=my_limits, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
legend("topright", legend=paste0("bias=",m_molo,"%"),cex=0.8)
# BOLAM
r <- ((bola_resam-obse_resam)/obse_resam)*100
m <- sprintf("%.1f",cellStats(r, 'mean'))
plot(bola_mask, main="(d) BOLAM", col = rev(viridis(100)), zlim=my_limits, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
legend("topright", legend=paste0("bias=",m_bola,"%"),cex=0.8)
dev.off()

###
# Now plot
my_eps <- paste0(my_figs_dir,"/tp_wethourpctl",pctl,"_vs_",ref_dataset,"_",my_seas,"_ano.eps")
setEPS()
postscript(my_eps)
par(mfrow = c(2, 2))
# GRIPHO
plot(obse_mask, main=paste0("(a) GRIPHO [mm/1-hour]\n",my_seas), col = rev(viridis(20)), zlim=my_limits, axes=FALSE, legend.width=3, legend=T, box=F, addfun=add_shp)
# ERA5-Land
plot((((era5_resam-obse_resam)/obse_resam))*100, main="(b) ERA5-Land [%]", col = my_pal, zlim=my_limits_perc, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
# plot (rasterToContour((((era5_resam-obse_resam)/obse_resam))*100), levels=c(-10), add=TRUE)
legend("topright", legend=paste0("bias=",m_era5,"%"),cex=0.8)
# MOLOCH
r <- ((molo_resam-obse_resam)/obse_resam)*100
m <- sprintf("%.1f",cellStats(r, 'mean'))
plot((((molo_resam-obse_resam)/obse_resam))*100, main="(c) MOLOCH [%]", col = my_pal, zlim=my_limits_perc, axes=FALSE, addfun=add_shp, legend.width=3, legend=T, box=F)
legend("topright", legend=paste0("bias=",m_molo,"%"),cex=0.8)
# BOLAM
r <- ((bola_resam-obse_resam)/obse_resam)*100
m <- sprintf("%.1f",cellStats(r, 'mean'))
plot((((bola_resam-obse_resam)/obse_resam))*100, main="(d) BOLAM [%]", col = my_pal, zlim=my_limits_perc, axes=FALSE, addfun=add_shp, legend.width=3, legend=F, box=F)
legend("topright", legend=paste0("bias=",m_bola,"%"),cex=0.8)
dev.off()
}

# Byebye
print(paste0("Byebye"))
# quit()
