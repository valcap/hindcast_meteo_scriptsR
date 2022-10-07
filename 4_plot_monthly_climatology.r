# Load libraries
library(raster)
library(ncdf4)
library(hydroGOF)

##########################################
# General settings
u                <- 'capecchi'
my_working_dir   <- paste0("/home/",u,"/art_hindcast_meteo")
my_figs_dir      <- paste0(my_working_dir,"/figs")
my_output_dir    <- paste0(my_working_dir,"/res")
my_shp_dir       <- paste0("/home/",u,"Sync/shp")
sf_era5          <- 1000

# Set working directory
setwd(my_working_dir)

#########################################
# Load data 2001-2016 (GRIPHO)
file_grip        <- paste0("DATA/GRIPHO/tp/gripho-v2-3km_mon_tp_AVG2001-2016.nc")
file_era5        <- paste0("DATA/ERA5-Land/tp/ERA5-Land_mon_tp_AVG2001-2016.nc")
file_molo        <- paste0("DATA/MOLOCH/tp/moloch_mon_tp_AVG2001-2016.nc")
file_bola        <- paste0("DATA/BOLAM/tp/bolam_mon_tp_AVG2001-2016.nc")
grip_data        <- as.integer(brick(file_grip))
era5_data        <- brick(file_era5)
molo_data        <- as.integer(brick(file_molo))
bola_data        <- as.integer(brick(file_bola))

# Load point locations (over which climatology is computed)
my_anag_file     <- paste0("DATA/points025_italy_extent.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE

# Extract data and calculates mean over the columns
grip_values      <- as.data.frame(extract(grip_data, pointCoordinates))
grip_mon_mean    <- colMeans(grip_values,na.rm=TRUE)
era5_values      <- as.data.frame(extract(as.integer(era5_data*sf_era5), pointCoordinates))
era5_mon_mean    <- colMeans(era5_values,na.rm=TRUE)
molo_values      <- as.data.frame(extract(molo_data, pointCoordinates))
molo_mon_mean    <- colMeans(molo_values,na.rm=TRUE)
bola_values      <- as.data.frame(extract(bola_data, pointCoordinates))
bola_mon_mean    <- colMeans(bola_values,na.rm=TRUE)

setEPS()
postscript(paste0(my_figs_dir,"/tp_MonthlyMean_vs_ALL.eps"))#, width=7, height=4.5)
par(mfrow=c(3,1))
par(mar = c(3, 3, 3, 8), xpd = TRUE)
# Now plot
#setEPS()
#postscript(paste0(my_figs_dir,"/tp_MonthlyMean_vs_GRIPHO.eps"), width=7, height=4.5)
my_min           <- as.integer(min(era5_mon_mean,molo_mon_mean,grip_mon_mean,bola_mon_mean)-10)
my_max           <- as.integer(max(era5_mon_mean,molo_mon_mean,grip_mon_mean,bola_mon_mean)+10)
# my_min           <- 40
# my_max           <- 160
my_months        <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
plot(grip_mon_mean,type="l", col="black", lwd=4, ylab="Average monthly precipitation [mm]", xlab="Month", xaxt="n", ylim=c(my_min,my_max), main="GRIPHO - 2001-2016")
axis(1, at=seq(1,12),labels=my_months, col.axis="black", las=2)
lines(era5_mon_mean, col="red", lwd=2.5)
points(era5_mon_mean, col="red", lwd=2.5, pch=19, cex=1.5)
lines(molo_mon_mean, col="orange", lwd=2.5)
points(molo_mon_mean, col="orange", lwd=2.5, pch=19, cex=1.5)
lines(bola_mon_mean, col="blue", lwd=2.5)
points(bola_mon_mean, col="blue", lwd=2.5, pch=19, cex=1.5)
legend("topleft",legend=c("GRIPHO","ERA5-Land","MOLOCH","BOLAM"), lty=1, pch=c(NA,19,19,19), col=c("black","red","orange","blue"), bg="white", lwd=c(2.5,1.5,1.5,1.5))
#dev.off()

# KGE with the hydroGOF package
kge_era5 <- hydroGOF::KGE(era5_mon_mean, grip_mon_mean)
kge_molo <- hydroGOF::KGE(molo_mon_mean, grip_mon_mean)
kge_bola <- hydroGOF::KGE(bola_mon_mean, grip_mon_mean)
print(paste("KGE","GRIPHO","ERA5-Land",sprintf("%.2f",kge_era5),"MOLOCH",sprintf("%.2f",kge_molo),"BOLAM",sprintf("%.2f",kge_bola),sep=";"),quote=FALSE)

#########################################
# Load data 1981-2015 (ARCIS)
file_arci        <- paste0("DATA/ARCIS/tp/ARCIS_mon_tp_AVG1981-2015.nc")
file_era5        <- paste0("DATA/ERA5-Land/tp/ERA5-Land_mon_tp_AVG1981-2015.nc")
file_molo        <- paste0("DATA/MOLOCH/tp/moloch_mon_tp_AVG1981-2015.nc")
file_bola        <- paste0("DATA/BOLAM/tp/bolam_mon_tp_AVG1981-2015.nc")
arci_data        <- as.integer(brick(file_arci))
era5_data        <- brick(file_era5)
molo_data        <- as.integer(brick(file_molo))
bola_data        <- as.integer(brick(file_bola))

# Load point locations (over which climatology is computed)
my_anag_file     <- paste0("DATA/points025_arcis_extent.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE

# Extract data and calculates mean over the columns
arci_values      <- as.data.frame(extract(arci_data, pointCoordinates))
arci_mon_mean    <- colMeans(arci_values,na.rm=TRUE)
era5_values      <- as.data.frame(extract(as.integer(era5_data*sf_era5), pointCoordinates))
era5_mon_mean    <- colMeans(era5_values,na.rm=TRUE)
molo_values      <- as.data.frame(extract(molo_data, pointCoordinates))
molo_mon_mean    <- colMeans(molo_values,na.rm=TRUE)
bola_values      <- as.data.frame(extract(bola_data, pointCoordinates))
bola_mon_mean    <- colMeans(bola_values,na.rm=TRUE)

# Now plot
#setEPS()
#postscript(paste0(my_figs_dir,"/tp_MonthlyMean_vs_ARCIS.eps"), width=7, height=4.5)
my_min           <- as.integer(min(era5_mon_mean,molo_mon_mean,arci_mon_mean,bola_mon_mean)-10)
my_max           <- as.integer(max(era5_mon_mean,molo_mon_mean,arci_mon_mean,bola_mon_mean)+10)
# my_min           <- 40
# my_max           <- 140
my_months        <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
plot(arci_mon_mean,type="l", col="black", lwd=4, ylab="Average monthly precipitation [mm]", xlab="Month", xaxt="n", ylim=c(my_min,my_max), main="ARCIS - 1981-2015")
axis(1, at=seq(1,12),labels=my_months, col.axis="black", las=2)
lines(era5_mon_mean, col="red", lwd=2.5)
points(era5_mon_mean, col="red", lwd=2.5, pch=19, cex=1.5)
lines(molo_mon_mean, col="orange", lwd=2.5)
points(molo_mon_mean, col="orange", lwd=2.5, pch=19, cex=1.5)
lines(bola_mon_mean, col="blue", lwd=2.5)
points(bola_mon_mean, col="blue", lwd=2.5, pch=19, cex=1.5)
legend("topleft",legend=c("ARCIS","ERA5-Land","MOLOCH","BOLAM"), lty=1, pch=c(NA,19,19,19), col=c("black","red","orange","blue"), bg="white", lwd=c(2.5,1.5,1.5,1.5))
#dev.off()

# KGE with the hydroGOF package
kge_era5 <- hydroGOF::KGE(era5_mon_mean, arci_mon_mean)
kge_molo <- hydroGOF::KGE(molo_mon_mean, arci_mon_mean)
kge_bola <- hydroGOF::KGE(bola_mon_mean, arci_mon_mean)
print(paste("KGE","ARCIS","ERA5-Land",sprintf("%.2f",kge_era5),"MOLOCH",sprintf("%.2f",kge_molo),"BOLAM",sprintf("%.2f",kge_bola),sep=";"),quote=FALSE)

#########################################
# Load data 1981-2019 (E-OBS)
file_eobs        <- paste0("DATA/E-OBS/tp/E-OBS_mon_tp_AVG1981-2019.nc")
file_era5        <- paste0("DATA/ERA5-Land/tp/ERA5-Land_mon_tp_AVG1981-2019.nc")
file_molo        <- paste0("DATA/MOLOCH/tp/moloch_mon_tp_AVG1981-2019.nc")
file_bola        <- paste0("DATA/BOLAM/tp/bolam_mon_tp_AVG1981-2019.nc")
eobs_data        <- as.integer(brick(file_eobs))
era5_data        <- brick(file_era5)
molo_data        <- as.integer(brick(file_molo))
bola_data        <- as.integer(brick(file_bola))

# Load point locations (over which climatology is computed)
my_anag_file     <- paste0("DATA/points025_italy_extent.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE

# Extract data and calculates mean over the columns
eobs_values      <- as.data.frame(extract(eobs_data, pointCoordinates))
eobs_mon_mean    <- colMeans(eobs_values,na.rm=TRUE)
era5_values      <- as.data.frame(extract(as.integer(era5_data*sf_era5), pointCoordinates))
era5_mon_mean    <- colMeans(era5_values,na.rm=TRUE)
molo_values      <- as.data.frame(extract(molo_data, pointCoordinates))
molo_mon_mean    <- colMeans(molo_values,na.rm=TRUE)
bola_values      <- as.data.frame(extract(bola_data, pointCoordinates))
bola_mon_mean    <- colMeans(bola_values,na.rm=TRUE)

# Now plot
#setEPS()
#postscript(paste0(my_figs_dir,"/tp_MonthlyMean_vs_E-OBS.eps"), width=7, height=4.5)
my_min           <- as.integer(min(era5_mon_mean,molo_mon_mean,eobs_mon_mean,bola_mon_mean)-10)
my_max           <- as.integer(max(era5_mon_mean,molo_mon_mean,eobs_mon_mean,bola_mon_mean)+10)
my_months        <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
# my_min           <- 40
# my_max           <- 140
plot(eobs_mon_mean,type="l", col="black", lwd=4, ylab="Average monthly precipitation [mm]", xlab="Month", xaxt="n", ylim=c(my_min,my_max), main="E-OBS - 1981-2019")
axis(1, at=seq(1,12),labels=my_months, col.axis="black", las=2)
lines(era5_mon_mean, col="red", lwd=2.5)
points(era5_mon_mean, col="red", lwd=2.5, pch=19, cex=1.5)
lines(molo_mon_mean, col="orange", lwd=2.5)
points(molo_mon_mean, col="orange", lwd=2.5, pch=19, cex=1.5)
lines(bola_mon_mean, col="blue", lwd=2.5)
points(bola_mon_mean, col="blue", lwd=2.5, pch=19, cex=1.5)
legend("topleft",legend=c("E-OBS","ERA5-Land","MOLOCH","BOLAM"), lty=1, pch=c(NA,19,19,19), col=c("black","red","orange","blue"), bg="white", lwd=c(2.5,1.5,1.5,1.5))
dev.off()

kge_era5 <- hydroGOF::KGE(era5_mon_mean, eobs_mon_mean)
kge_molo <- hydroGOF::KGE(molo_mon_mean, eobs_mon_mean)
kge_bola <- hydroGOF::KGE(bola_mon_mean, eobs_mon_mean)
print(paste("KGE","E-OBS","ERA5-Land",sprintf("%.2f",kge_era5),"MOLOCH",sprintf("%.2f",kge_molo),"BOLAM",sprintf("%.2f",kge_bola),sep=";"),quote=FALSE)

#quit()
