# Load libraries
library(raster)
library(tabularaster)
library(tidyr)
# add here any other packages you need

# Global variables
u               <- 'capecchi'
my_working_dir  <- paste0("/home/",u,"/art_hindcast_meteo/DATA")
my_figs_dir     <- paste0("/home/",u,"/art_hindcast_meteo/figs")
my_output_dir   <- paste0("/home/",u,"/art_hindcast_meteo/res")
my_shp_dir      <- paste0("/home/",u,"/Sync/shp")
my_pal          <- brewer.pal(n = 10, name = "RdBu")
yini            <- 2002
yend            <- 2016
pctl            <- 90

# Working directory
setwd(my_working_dir)

#############################
ref_dataset     <- 'GRIPHO'
#############################
# pctl th percentile of hourly mean precipitation
#file_obse   <- paste0("GRIPHO/tp/gripho-v2-3km_wethourpctl",pctl,"_tp_",yini,"-",yend,"_JJA.nc")
#file_era5   <- paste0("ERA5-Land/tp/ERA5-Land_wethourpctl",pctl,"_tp_",yini,"-",yend,"_JJA.nc")
#file_bola   <- paste0("BOLAM/tp/bolam_wethourpctl",pctl,"_tp_",yini,"-",yend,"_JJA.nc") 
#file_molo   <- paste0("MOLOCH/tp/moloch_wethourpctl",pctl,"_tp_",yini,"-",yend,"_JJA.nc")
for (pctl in c(90,95,99)) {
  if (pctl==90) {my_max = 50; my_ylim=0.1}
  if (pctl==95) {my_max = 75; my_ylim=0.1}
  if (pctl==99) {my_max = 175; my_ylim=0.05}
  my_limits       <- c(1,my_max)
file_obse   <- paste0("GRIPHO/tp/gripho-v2-3km_wetdaypctl",pctl,"_tp_",yini,"-",yend,"_JJA.nc")
file_era5   <- paste0("ERA5-Land/tp/ERA5-Land_wetdaypctl",pctl,"_tp_",yini,"-",yend,"_JJA.nc")
file_bola   <- paste0("BOLAM/tp/bolam_wetdaypctl",pctl,"_tp_",yini,"-",yend,"_JJA.nc")
file_molo   <- paste0("MOLOCH/tp/moloch_wetdaypctl",pctl,"_tp_",yini,"-",yend,"_JJA.nc")
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
print(paste0("+++Cropping rasters..."))
# Crop ERA5-Land grid on my_extent and regrid (if needed)
era5_crop           <- crop(era5_orig, my_extent)
era5_tbl            <- as_tibble(era5_crop*1000)
#  Observations
obse_crop           <- obse_orig
obse_tbl            <- as_tibble(obse_crop)
# Crop BOLAM grid on my_extent and regrid (if needed)
bola_crop           <- crop(bola_orig, my_extent)
bola_tbl            <- as_tibble(bola_crop)
# Crop MOLOCH grid on my_extent and regrid (if needed)
molo_crop           <- crop(molo_orig, my_extent)
molo_tbl            <- as_tibble(molo_crop)
print(paste0("+++OK"))
###

###
# Now plot
my_eps <- paste0(my_figs_dir,"/tp_hist_wetday",pctl,"_vs_",ref_dataset,".eps")
setEPS()  
postscript(my_eps)
par(mfrow = c(2, 2))
# GRIPHO
hist(obse_tbl$cellvalue, main="GRIPHO", col = "gray", breaks="Scott", xlim=my_limits, freq=F, xlab="GRIPHO", ylim=c(0,my_ylim))
# ERA5-Land
hist(era5_tbl$cellvalue, main="ERA5-Land", col = "gray", breaks="Scott", xlim=my_limits, freq=F, xlab="ERA5-Land", ylim=c(0,my_ylim))
# MOLOCH
hist(molo_tbl$cellvalue, main="MOLOCH", col = "gray", breaks="Scott", xlim=my_limits, freq=F, xlab="MOLOCH", ylim=c(0,my_ylim))
# BOLAM
hist(bola_tbl$cellvalue, main="BOLAM", col = "gray", breaks="Scott", xlim=my_limits, freq=F, xlab="BOLAM", ylim=c(0,my_ylim))
dev.off()

###
# Now plot
my_eps <- paste0(my_figs_dir,"/tp_cdf_wetday",pctl,"_vs_",ref_dataset,".eps")
setEPS()  
postscript(my_eps)
# GRIPHO
p1<-ecdf(obse_tbl$cellvalue)
plot(p1,col="black",xlim=c(0,my_max), main="TODO", xlab="TODO", ylab="TODO")
# ERA5-Land
p2<-ecdf(era5_tbl$cellvalue)
plot(p2,col="red",add=TRUE)
# MOLOCH
p3<-ecdf(molo_tbl$cellvalue)
plot(p3,col="orange",add=TRUE)
# BOLAM
p4<-ecdf(bola_tbl$cellvalue)
plot(p4,col="blue",add=TRUE)
dev.off()
}

# Byebye
print(paste0("Byebye"))
# quit()
