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
my_working_dir  <- paste0("/home/",u,"/art_hindcast_meteo")
my_figs_dir     <- paste0(my_working_dir,"/figs")
my_output_dir   <- paste0(my_working_dir,"/res")
my_shp_dir      <- paste0("/home/",u,"Sync/shp")
TD_normal       <- FALSE
sf_era5         <- 1000

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
file_eobs   <- paste0("DATA/E-OBS/tp/E-OBS_yearsum_tp_AVG.nc")
file_era5   <- paste0("DATA/ERA5-Land/tp/ERA5-Land_year_tp_AVG.nc")
file_molo   <- paste0("DATA/MOLOCH/tp/moloch_year_apcp_AVG.nc")
file_bola   <- paste0("DATA/BOLAM/tp/bolam_year_apcp_AVG.nc")
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
names(eobs_resam)     <- "pcp"
names(era5_resam)     <- "pcp"
names(molo_resam)     <- "pcp"
names(bola_resam)     <- "pcp"
print(paste0("+++OK"))
###

###
# Extract values
eobs_values <- as.data.frame(extract(eobs_resam, pointCoordinates))
era5_values <- as.data.frame(extract(era5_resam*1000, pointCoordinates))
molo_values <- as.data.frame(extract(molo_resam, pointCoordinates))
bola_values <- as.data.frame(extract(bola_resam, pointCoordinates))
eobs_values_ts <- c(eobs_values[,1])
era5_values_ts <- c(era5_values[,1])
molo_values_ts <- c(molo_values[,1])
bola_values_ts <- c(bola_values[,1])
###

# E-OBS
my_sd=sprintf("%.2f",sd(eobs_values_ts,na.rm=T))
print(paste("E-OBS",my_sd,sep=";"))
print(paste("VER_DATASET","ST.DEV","CORR","C-RMSE","RMSE","ME","MBIAS","OBIAS",sep=" & "))

# ERA5-Land
my_sd=sprintf("%.0f",sd(era5_values_ts,na.rm=T))
my_cor=sprintf("%.2f",cor(eobs_values_ts, era5_values_ts,use="complete.obs"))
my_ver=verify(obs=eobs_values_ts, pred=era5_values_ts,frcst.type="cont",obs.type="cont")
my_rmse=as.numeric(sprintf("%.0f",sqrt(my_ver$MSE)))
my_me=sprintf("%.0f",(my_ver$ME))
my_mbias=sprintf("%.2f",mean(era5_values_ts,na.rm=T)/mean(eobs_values_ts,na.rm=T))
my_overall_bias=as.numeric(sprintf("%.0f",mean(era5_values_ts,na.rm=T)-mean(eobs_values_ts,na.rm=T)))
my_crmse=as.numeric(sprintf("%.0f",sqrt((my_rmse*my_rmse)-(my_overall_bias*my_overall_bias))))
print(paste("VER_ERA5-Land",my_sd,my_cor,my_crmse,my_rmse,my_me,my_mbias,my_overall_bias,sep=" & "))

# MOLOCH
my_sd=sprintf("%.0f",sd(molo_values_ts,na.rm=T))
my_cor=sprintf("%.2f",cor(eobs_values_ts, molo_values_ts,use="complete.obs"))
my_ver=verify(obs=eobs_values_ts, pred=molo_values_ts,frcst.type="cont",obs.type="cont")
my_rmse=as.numeric(sprintf("%.0f",sqrt(my_ver$MSE)))
my_me=sprintf("%.0f",(my_ver$ME))
my_mbias=sprintf("%.2f",mean(molo_values_ts,na.rm=T)/mean(eobs_values_ts,na.rm=T))
my_overall_bias=as.numeric(sprintf("%.0f",mean(molo_values_ts,na.rm=T)-mean(eobs_values_ts,na.rm=T)))
my_crmse=as.numeric(sprintf("%.0f",sqrt((my_rmse*my_rmse)-(my_overall_bias*my_overall_bias))))
print(paste("VER_MOLOCH",my_sd,my_cor,my_crmse,my_rmse,my_me,my_mbias,my_overall_bias,sep=" & "))

# BOLAM
my_sd=sprintf("%.0f",sd(bola_values_ts,na.rm=T))
my_cor=sprintf("%.2f",cor(eobs_values_ts, bola_values_ts,use="complete.obs"))
my_ver=verify(obs=eobs_values_ts, pred=bola_values_ts,frcst.type="cont",obs.type="cont")
my_rmse=as.numeric(sprintf("%.0f",sqrt(my_ver$MSE)))
my_me=sprintf("%.0f",(my_ver$ME))
my_mbias=sprintf("%.2f",mean(bola_values_ts,na.rm=T)/mean(eobs_values_ts,na.rm=T))
my_overall_bias=as.numeric(sprintf("%.0f",mean(bola_values_ts,na.rm=T)-mean(eobs_values_ts,na.rm=T)))
my_crmse=as.numeric(sprintf("%.0f",sqrt((my_rmse*my_rmse)-(my_overall_bias*my_overall_bias))))
print(paste("VER_BOLAM",my_sd,my_cor,my_crmse,my_rmse,my_me,my_mbias,my_overall_bias,sep=" & "))

###
# Taylor Diagram
print(paste0("Plotting Taylor diagrams"))
setEPS()
file_taylor_eps <- paste0(my_figs_dir,"/tp_Taylor_Diag_vs_EOBS_norm",TD_normal,"_Year.eps")
postscript(file_taylor_eps)
par(mfrow=c(2,2))
obs<-eobs_values_ts

# ERA5
pre <- era5_values_ts
taylor.diagram(obs, pre, add=F, col="red", ref.sd=TRUE, normalize=TD_normal, main="(a) Total annual precipitation", grad.corr.lines=c(0.7,0.8), ngamma=5, pch=19, pcex=1.5)

# MOLOCH
pre <- molo_values_ts
taylor.diagram(obs, pre, add=T, col="orange", ref.sd=TRUE, normalize=TD_normal, pch=19, pcex=1.5)

# BOLAM
pre <- bola_values_ts
taylor.diagram(obs, pre, add=T, col="blue", ref.sd=TRUE, normalize=TD_normal, pch=19, pcex=1.5)

# Draw the legend on the plot and save it
legend("topright",legend=c("ERA5-Land","MOLOCH","BOLAM"), pch=19, col=c("red","orange","blue"),cex=0.8,bty = "n")
###

# Scatter plots
plot(eobs_values_ts, era5_values_ts,
     xlim=c(250,3000) , ylim=c(250,3000),
     pch=19,
     cex=1,
     col="red",
     xlab="E-OBS",
     ylab="ERA5-Land",
     main="(b) ERA5-Land"
)
lines(x = c(250,3000), y = c(250,3000))
#abline(coef = c(0,1))

plot(eobs_values_ts, molo_values_ts,
     xlim=c(250,3000) , ylim=c(250,3000),
     pch=19,
     cex=1,
     col="orange",
     xlab="E-OBS", 
     ylab="MOLOCH",
     main="(c) MOLOCH"
)
lines(x = c(250,3000), y = c(250,3000))

plot(eobs_values_ts, bola_values_ts,
     xlim=c(250,3000) , ylim=c(250,3000),
     pch=19,
     cex=1,
     col="blue",
     xlab="E-OBS",
     ylab="BOLAM",
     main="(d) BOLAM"
)
lines(x = c(250,3000), y = c(250,3000))

dev.off()
print(paste0("OK"))

