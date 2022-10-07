# Load libraries
library(raster)
library(ncdf4)
library(verification)
library(viridis)
library(rgdal)
library(tabularaster)
library(SpatialVx)
library(zoo)

# Global variables
u               <- 'capecchi'
my_working_dir  <- paste0("/home/",u,"/art_hindcast_meteo")
my_figs_dir     <- paste0(my_working_dir,"/figs")
my_output_dir   <- paste0(my_working_dir,"/res")
my_shp_dir      <- paste0("/home/",u,"/Sync/shp")
my_limits       <- c(0,450)

# Working directory
setwd(my_working_dir)

# Shapefiles
my_coastline_lr <- readOGR(
  dsn= paste0(my_shp_dir,"/gshhg-shp-2.3.7/GSHHS_shp/l/") ,
  layer="GSHHS_l_L1",
  verbose=FALSE)
add_shp=function(){plot(my_coastline_lr, bg="transparent", add=TRUE,lwd=0.5)}

setEPS()
file_taylor_eps <- paste0(my_figs_dir,"/tp_FSS.eps")
postscript(file_taylor_eps, horizontal = FALSE, onefile = FALSE,
           paper = "special", height = 5, width =12)
par(mfrow = c(1, 2))
# ##########################################################################################
# 5/Nov/1994
obse_dataset="ARCIS"
vec_thres       <- seq(1,250,1)
study_case      <- "PIE1994"
# Load observations
if (obse_dataset=="GRIPHO"){
   file_obse     <- paste0(my_working_dir,"/DATA/",obse_dataset,"/tp/gripho-v2-3km_20111025_tp.nc")
} else {
   file_obse     <- paste0(my_working_dir,"/DATA/",obse_dataset,"/tp/ARCIS_19941105_tp.nc")
}
file_era5     <- paste0(my_working_dir,"/DATA/ERA5-Land/tp/ERA5-Land_19941105_tp.grib2")
file_molo     <- paste0(my_working_dir,"/DATA/MOLOCH/tp/moloch_19941105-daysum_tp.nc")
file_bola     <- paste0(my_working_dir,"/DATA/BOLAM/tp/bolam_19941105-daysum_tp.nc")
obse_orig     <- brick(file_obse)
era5_orig     <- brick(file_era5)
molo_orig     <- brick(file_molo)
bola_orig     <- brick(file_bola)

# Load the point locations
my_anag_file     <- paste0("/home/capecchi/casi_studio/PIE94/data/obs_giornalieri/COD-ALL_cum_P1105.txt")
anag_full        <- read.csv(my_anag_file, sep=";", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))
print(paste0(""))
###

###
# Get extent on the basis of the point coordinates
print(paste0("+++Cropping extents and resampling rasters..."))
my_extent           <- extent(pointCoordinates)
# ERA5
# Crop ERA5-Land grid on my_extent and regrid (if needed)
obse_crop                <- crop(obse_orig, my_extent)
era5_crop                <- crop(era5_orig, my_extent)
molo_crop                <- crop(molo_orig, my_extent)
bola_crop                <- crop(bola_orig, my_extent)
# Resampling
obse_resam2era5 <- resample(obse_crop, era5_crop, method="ngb")
obse_resam2molo <- resample(obse_crop, molo_crop, method="ngb")
obse_resam2bola <- resample(obse_crop, bola_crop, method="ngb")
#
era5_resam2obse <- resample(era5_crop, obse_crop, method="ngb")
era5_resam2molo <- resample(era5_crop, molo_crop, method="ngb")
era5_resam2bola <- resample(era5_crop, bola_crop, method="ngb")
#
molo_resam2obse <- resample(molo_crop, obse_crop, method="ngb")
molo_resam2era5 <- resample(molo_crop, era5_crop, method="ngb")
molo_resam2bola <- resample(molo_crop, bola_crop, method="ngb")
#
bola_resam2obse <- resample(bola_crop, obse_crop, method="ngb")
bola_resam2era5 <- resample(bola_crop, era5_crop, method="ngb")
bola_resam2molo <- resample(bola_crop, molo_crop, method="ngb")

# # # plot
## setEPS()
## file_eps <- paste0(my_figs_dir,"/casestudy_",obse_dataset,"_",study_case,".eps")
## postscript(file_eps)
## par(mfrow=c(2,2))
## plot(obse_resam2molo, main=paste0("(a) ",obse_dataset," [mm/day]"), col = rev(viridis(100)), zlim=my_limits, axes=FALSE, legend.width=3, legend=T, box=T, addfun=add_shp)
## plot(era5_resam2molo*1000, main="(b) ERA5-Land [mm/day]", col = rev(viridis(100)), zlim=my_limits, axes=FALSE, legend.width=3, legend=F, box=T, addfun=add_shp)
## plot(molo_crop, main="(c) MOLOCH [mm/day]", col = rev(viridis(100)), zlim=my_limits, axes=FALSE, legend.width=3, legend=F, box=T, addfun=add_shp)
## plot(bola_resam2molo, main="(d) BOLAM [mm/day]", col = rev(viridis(100)), zlim=my_limits, axes=FALSE, legend.width=3, legend=F, box=T, addfun=add_shp)
## dev.off()

# Create matrix from raster
obse_tibble        <-  tabularaster::as_tibble(obse_resam2bola)
era5_tibble        <-  tabularaster::as_tibble(era5_resam2bola*1000)
molo_tibble        <-  tabularaster::as_tibble(molo_resam2bola)
bola_tibble        <-  tabularaster::as_tibble(bola_crop)
obse_matrix        <-  matrix(obse_tibble$cellvalue ,nrow=nrow(bola_crop),ncol=ncol(bola_crop))
era5_matrix        <-  matrix(era5_tibble$cellvalue ,nrow=nrow(bola_crop),ncol=ncol(bola_crop))
molo_matrix        <-  matrix(molo_tibble$cellvalue ,nrow=nrow(bola_crop),ncol=ncol(bola_crop))
bola_matrix        <-  matrix(bola_tibble$cellvalue ,nrow=nrow(bola_crop),ncol=ncol(bola_crop))

# Create binary matrix and calculate FSS
print(paste0("THRES;FSS_ERA5;FSS_MOLO;FSS_BOLA"))
fss_era5<-NULL
fss_molo<-NULL
fss_bola<-NULL
ets_era5<-NULL
ets_molo<-NULL
ets_bola<-NULL
for (i in vec_thres) {
   thres  <- i
   obse_binary          <- replace(obse_matrix, obse_matrix <= thres, 0)
   obse_binary          <- as.matrix(replace(obse_binary, obse_binary  > thres, 1))
   era5_binary          <- replace(era5_matrix, era5_matrix <= thres, 0)
   era5_binary          <- as.matrix(replace(era5_binary, era5_binary  > thres, 1))
   molo_binary          <- replace(molo_matrix, molo_matrix <= thres, 0)
   molo_binary          <- as.matrix(replace(molo_binary, molo_binary  > thres, 1))
   bola_binary          <- replace(bola_matrix, bola_matrix <= thres, 0)
   bola_binary          <- as.matrix(replace(bola_binary, bola_binary  > thres, 1))
   fss_era5[i] <- fss(obs = obse_binary, pred = era5_binary, w=2, fun=max)
   fss_molo[i] <- fss(obs = obse_binary, pred = molo_binary, w=2, fun=max)
   fss_bola[i] <- fss(obs = obse_binary, pred = bola_binary, w=2, fun=max)
   # ets_era5[i] <- MinCvg2dfun(era5_binary,obse_binary)$ets
   # ets_molo[i] <- MinCvg2dfun(molo_binary,obse_binary)$ets
   # ets_bola[i] <- MinCvg2dfun(bola_binary,obse_binary)$ets
   # fss_era5[i] <- fss2dfun(era5_binary, obse_binary)
   # fss_molo[i] <- fss2dfun(molo_binary, obse_binary)
   # fss_bola[i] <- fss2dfun(bola_binary, obse_binary)
   # print(paste0("THRES",thres,";",sprintf("%.3f",fss_era5[i]),";",sprintf("%.3f",fss_molo),";",sprintf("%.3f",fss_bola)))
}
# setEPS()
# file_taylor_eps <- paste0(my_figs_dir,"/tp_FSS_vs_",obse_dataset,"_",study_case,".eps")
# postscript(file_taylor_eps, horizontal = FALSE, onefile = FALSE,
#            paper = "special", height = 6, width = 9)
b <- mean(obse_binary,na.rm=T) ## base rate
fss.uniform  <- 0.5 + b/2
fss.random   <- b
plot(zoo::rollmean(fss_era5, k = 5, fill = NA), ylim = c(0,1), main = study_case, ylab = "FSS", xlab = "Precipitation threshold [mm/day]",
     type = "l", lty = 1, axes = T , lwd = 2, col="red")
abline(h = c(fss.uniform), lty = 2)
abline(v = c(50,100,125,150,175,200,225,250), lty = 2, lwd=0.5, col="gray")
lines(zoo::rollmean(fss_bola, k = 5, fill = NA), lty = 1, lwd = 2, col="blue")
lines(zoo::rollmean(fss_molo, k = 5, fill = NA), lty = 1, lwd = 2, col="orange")
legend("bottomleft", legend=c("ERA5-Land","MOLOCH","BOLAM"),
       lwd=c(2,2,2), col=c("red","orange","blue"), bg="white")

# setEPS()
# file_taylor_eps <- paste0(my_figs_dir,"/tp_ETS_vs_",obse_dataset,"_",study_case,".eps")
# postscript(file_taylor_eps, horizontal = FALSE, onefile = FALSE,
#            paper = "special", height = 6, width = 9)
# plot(ets_era5, ylim = c(0,1), main = study_case, ylab = "ETS", xlab = "Precipitation threshold [mm/day]",
#      type = "l", lty = 1, axes = T , lwd = 2, col="red")
# lines(ets_bola, lty = 1, lwd = 2, col="blue")
# lines(ets_molo, lty = 1, lwd = 2, col="orange")
# legend("topright", legend=c("ERA5-Land","MOLOCH","BOLAM"),
#        lwd=c(2,2,2), col=c("red","orange","blue"), bg="white")
# dev.off()

# # ###1# ###1# ###1# ###1# ###1# ###1# ###1# ###1# ###1# ###1# ###1# ###1# ###1# ###1# ###1
 # 25/Oct/2011
 vec_thres       <- seq(0,200,1)
 study_case      <- "CT2011"
 obse_dataset    <- "ARCIS"
 # Load observations
 if (obse_dataset=="GRIPHO"){
   file_obse     <- paste0(my_working_dir,"/DATA/",obse_dataset,"/tp/gripho-v2-3km_20111025_tp.nc")
 } else {
   file_obse     <- paste0(my_working_dir,"/DATA/",obse_dataset,"/tp/ARCIS_20111025_tp.nc")
 }
 file_era5     <- paste0(my_working_dir,"/DATA/ERA5-Land/tp/ERA5-Land_20111025_tp.grib2")
 file_molo     <- paste0(my_working_dir,"/DATA/MOLOCH/tp/moloch_20111024-daysum_tp.nc")
 file_bola     <- paste0(my_working_dir,"/DATA/BOLAM/tp/bolam_20111024-daysum_tp.nc")
 obse_orig     <- brick(file_obse)
 era5_orig     <- brick(file_era5)
 molo_orig     <- brick(file_molo)
 bola_orig     <- brick(file_bola)

 # Load the point locations
 my_anag_file     <- paste0("/home/capecchi/ECMWF_SP/UEF2019/data/obs/PLUV_25oct2011.txt")
 anag_full        <- read.csv(my_anag_file, sep=";", header=TRUE)
 pointCoordinates <- anag_full
 coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
 print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))
 print(paste0(""))
 ###

 ###
 # Get extent on the basis of the point coordinates
 print(paste0("+++Cropping extents and resampling rasters..."))
 my_extent           <- extent(pointCoordinates)
 # ERA5
 # Crop ERA5-Land grid on my_extent and regrid (if needed)
 era5_crop                <- crop(era5_orig, my_extent)
 obse_crop                <- crop(obse_orig, my_extent)
 molo_crop                <- crop(molo_orig, my_extent)
 bola_crop                <- crop(bola_orig, my_extent)
 # Resampling
 obse_resam2era5 <- resample(obse_crop, era5_crop, method="ngb")
 obse_resam2molo <- resample(obse_crop, molo_crop, method="ngb")
 obse_resam2bola <- resample(obse_crop, bola_crop, method="ngb")
 #
 era5_resam2obse <- resample(era5_crop, obse_crop, method="ngb")
 era5_resam2molo <- resample(era5_crop, molo_crop, method="ngb")
 era5_resam2bola <- resample(era5_crop, bola_crop, method="ngb")
 #
 molo_resam2obse <- resample(molo_crop, obse_crop, method="ngb")
 molo_resam2era5 <- resample(molo_crop, era5_crop, method="ngb")
 molo_resam2bola <- resample(molo_crop, bola_crop, method="ngb")
 #
 bola_resam2obse <- resample(bola_crop, obse_crop, method="ngb")
 bola_resam2era5 <- resample(bola_crop, era5_crop, method="ngb")
 bola_resam2molo <- resample(bola_crop, molo_crop, method="ngb")

# plot
 #setEPS()
 #file_eps <- paste0(my_figs_dir,"/casestudy_",obse_dataset,"_",study_case,".eps")
 #postscript(file_eps)
 #par(mfrow=c(2,2))
 #plot(obse_resam2molo, main=paste0("(a) ",obse_dataset," [mm/day]"), col = rev(viridis(100)), zlim=my_limits, axes=FALSE, legend.width=3, legend=T, box=T, addfun=add_shp)
 #plot(era5_resam2molo*1000, main="(b) ERA5-Land [mm/day]", col = rev(viridis(100)), zlim=my_limits, axes=FALSE, legend.width=3, legend=F, box=T, addfun=add_shp)
 #plot(molo_crop, main="(c) MOLOCH [mm/day]", col = rev(viridis(100)), zlim=my_limits, axes=FALSE, legend.width=3, legend=F, box=T, addfun=add_shp)
 #plot(bola_resam2molo, main="(d) BOLAM [mm/day]", col = rev(viridis(100)), zlim=my_limits, axes=FALSE, legend.width=3, legend=F, box=T, addfun=add_shp)
 #dev.off()

 # Create matrix from raster
 obse_tibble        <-  tabularaster::as_tibble(obse_resam2bola)
 era5_tibble        <-  tabularaster::as_tibble(era5_resam2bola*1000)
 molo_tibble        <-  tabularaster::as_tibble(molo_resam2bola)
 bola_tibble        <-  tabularaster::as_tibble(bola_crop)
 obse_matrix        <-  matrix(obse_tibble$cellvalue ,nrow=nrow(bola_crop),ncol=ncol(bola_crop))
 era5_matrix        <-  matrix(era5_tibble$cellvalue ,nrow=nrow(bola_crop),ncol=ncol(bola_crop))
 molo_matrix        <-  matrix(molo_tibble$cellvalue ,nrow=nrow(bola_crop),ncol=ncol(bola_crop))
 bola_matrix        <-  matrix(bola_tibble$cellvalue ,nrow=nrow(bola_crop),ncol=ncol(bola_crop))

 # Create binary matrix and calculate FSS
 print(paste0("THRES;FSS_ERA5;FSS_MOLO;FSS_BOLA"))
 fss_era5<-NULL
 fss_molo<-NULL
 fss_bola<-NULL
 ets_era5<-NULL
 ets_molo<-NULL
 ets_bola<-NULL
 for (i in vec_thres) {
   thres  <- i
   obse_binary          <- replace(obse_matrix, obse_matrix <= thres, 0)
   obse_binary          <- as.matrix(replace(obse_binary, obse_binary  > thres, 1))
   era5_binary          <- replace(era5_matrix, era5_matrix <= thres, 0)
   era5_binary          <- as.matrix(replace(era5_binary, era5_binary  > thres, 1))
   molo_binary          <- replace(molo_matrix, molo_matrix <= thres, 0)
   molo_binary          <- as.matrix(replace(molo_binary, molo_binary  > thres, 1))
   bola_binary          <- replace(bola_matrix, bola_matrix <= thres, 0)
   bola_binary          <- as.matrix(replace(bola_binary, bola_binary  > thres, 1))
   fss_era5[i] <- fss(obs = obse_binary, pred = era5_binary, w=2, fun=max)
   fss_molo[i] <- fss(obs = obse_binary, pred = molo_binary, w=2, fun=max)
   fss_bola[i] <- fss(obs = obse_binary, pred = bola_binary, w=2, fun=max)
   # ets_era5[i] <- MinCvg2dfun(era5_binary,obse_binary)$ets
   # ets_molo[i] <- MinCvg2dfun(molo_binary,obse_binary)$ets
   # ets_bola[i] <- MinCvg2dfun(bola_binary,obse_binary)$ets

 }
 # setEPS()
 # file_taylor_eps <- paste0(my_figs_dir,"/tp_FSS_vs_",obse_dataset,"_",study_case,".eps")
 # postscript(file_taylor_eps, horizontal = FALSE, onefile = FALSE,
 #            paper = "special", height = 6, width = 9)
 b <- mean(obse_binary,na.rm=T) ## base rate for maximum threshold
 fss.uniform  <- 0.5 + b/2
# fss.random   <- b
# plot(fss_era5, ylim = c(0,1), main = study_case, ylab = "FSS", xlab = "Precipitation threshold [mm/day]",
#      type = "l", lty = 1, axes = T , lwd = 2, col="red")
 plot(zoo::rollmean(fss_era5, k = 5, fill = NA), ylim = c(0,1), main = study_case, ylab = "FSS", xlab = "Precipitation threshold [mm/day]",
      type = "l", lty = 1, axes = T , lwd = 2, col="red")
 abline(h = c(fss.uniform), lty = 2)
 abline(v = c(50,100,125,150,175,200,225,250), lty = 2, lwd=0.5, col="gray")
 lines(zoo::rollmean(fss_bola, k = 5, fill = NA), lty = 1, lwd = 2, col="blue")
 lines(zoo::rollmean(fss_molo, k = 5, fill = NA), lty = 1, lwd = 2, col="orange")
 legend("bottomleft", legend=c("ERA5-Land","MOLOCH","BOLAM"),
        lwd=c(2,2,2), col=c("red","orange","blue"), bg="white")
 dev.off()

# setEPS()
# file_taylor_eps <- paste0(my_figs_dir,"/tp_ETS_vs_",obse_dataset,"_",study_case,".eps")
# postscript(file_taylor_eps, horizontal = FALSE, onefile = FALSE,
#            paper = "special", height = 6, width = 9)
# plot(ets_era5, ylim = c(0,1), main = study_case, ylab = "ETS", xlab = "Precipitation threshold [mm/day]",
#      type = "l", lty = 1, axes = T , lwd = 2, col="red")
# lines(ets_bola, lty = 1, lwd = 2, col="blue")
# lines(ets_molo, lty = 1, lwd = 2, col="orange")
# legend("topright", legend=c("ERA5-Land","MOLOCH","BOLAM"),
#        lwd=c(2,2,2), col=c("red","orange","blue"), bg="white")
# dev.off()

# byebye
#quit()
