# Purpose

# Load libraries
library(raster)
library(ncdf4)
library(verification)
library(Metrics)
library(plotrix)
library(lubridate)
library(dplyr)

# Global variables
my_working_dir  <- paste0("/OCEANASTORE/progetti/spitbran")
my_figs_dir     <- paste0("/OCEANASTORE/progetti/spitbran/work/figs")
my_output_dir   <- paste0("/OCEANASTORE/progetti/spitbran/work/res")
my_shp_dir      <- paste0("/OCEANASTORE/progetti/spitbran/work/shp")
TD_normal       <- FALSE
ini_year        <- 1981
end_year        <- 2019
vec_years       <- seq(ini_year,end_year)
sf_era5         <- 100
cutoff_NA       <- 50
eobs_values_MAM <- NULL 
era5_values_MAM <- NULL
molo_values_MAM <- NULL
bola_values_MAM <- NULL
eobs_values_JJA <- NULL 
era5_values_JJA <- NULL
molo_values_JJA <- NULL
bola_values_JJA <- NULL
eobs_values_SON <- NULL 
era5_values_SON <- NULL
molo_values_SON <- NULL
bola_values_SON <- NULL
eobs_values_DJF <- NULL 
era5_values_DJF <- NULL
molo_values_DJF <- NULL
bola_values_DJF <- NULL

# Working directory
setwd(my_working_dir)

# Load the point locations
my_anag_file     <- paste0(my_working_dir,"/E-OBS/points050_moloch_extent.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))
print(paste0(""))
setEPS()
postscript(paste0(my_figs_dir,"/virtual_rain-gauges_EOBS.eps"))
plot(pointCoordinates, main=paste0("Number of virtual rain-gauges ",length(pointCoordinates)))
plot(extent(pointCoordinates), add=TRUE, col='red', lwd=1)
dev.off()

######################
# MAIN LOOP OVER YEARS
for (my_year in vec_years) {
  print(paste0("___Running year ",my_year,"___"))

###
  # Loading raster file names
  print(paste0("+++Reading raster files..."))
  # Raster files definition
  file_eobs   <- paste0("E-OBS/tp/E-OBS_day_tp_",my_year,".nc")
  file_era5   <- paste0("ERA5-Land/tp/ERA5-Land_",my_year,"-daysum_tp.nc")
  file_molo   <- paste0("work/moloch/moloch_",my_year,"_daysum_apcp.nc")
  file_bola   <- paste0("work/bolam/bolam_",my_year,"_daysum_apcp.nc") 
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
  print(paste0("+++Cropping extents and resampling rasters..."))
  my_extent           <- extent(pointCoordinates)
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

###
  # select only data from the current year
  print(paste0("+++Adjusting name conventions..."))
  names(eobs_resam)     <- as.Date(names(eobs_resam),format="X%Y.%m.%d")
  if ((my_year == '1990') | (my_year == '2000') | (my_year == '2010')) { 
    names(era5_resam)   <- as.Date(names(era5_resam),format="X%Y%m%d.4791667")
  } else {
    names(era5_resam)   <- as.Date(names(era5_resam),format="X%Y.%m.%d.%H.%M.%S")
  }
  names(bola_resam)     <- as.Date(names(bola_resam),format="X%Y.%m.%d.%H.%M.%S")
  names(molo_resam)     <- names(bola_resam)
  # select only data from the current year
  print(paste0("+++selecting data from the current year..."))
  # E-OBS
  r <- setZ(eobs_resam, as.Date(names(eobs_resam),format="X%Y.%m.%d"))
  eobs_ok <- subset(r, which(getZ(r) >= paste0(my_year,'-01-02') & (getZ(r) <= paste0(my_year,'-12-31'))))
  # ERA5
  r <- setZ(era5_resam, as.Date(names(era5_resam),format="X%Y.%m.%d"))
  era5_ok <- subset(r, which(getZ(r) >= paste0(my_year,'-01-02') & (getZ(r) <= paste0(my_year,'-12-31'))))
  # MOLOCH
  r <- setZ(molo_resam, as.Date(names(molo_resam),format="X%Y.%m.%d"))
  molo_ok <- subset(r, which(getZ(r) >= paste0(my_year,'-01-02') & (getZ(r) <= paste0(my_year,'-12-31'))))
  # BOLAM
  r <- setZ(bola_resam, as.Date(names(bola_resam),format="X%Y.%m.%d"))
  bola_ok <- subset(r, which(getZ(r) >= paste0(my_year,'-01-02') & (getZ(r) <= paste0(my_year,'-12-31'))))
  print(paste0("+++OK"))
###

###
  # Creating monthly data 
  print(paste0("+++Creating monthly data"))
  index_months         <- format(as.Date(names(molo_ok), format = "X%Y.%m.%d"), format = "%m")
  eobs_ok_month        <- stackApply(eobs_ok, index_months, fun = sum)
  names(eobs_ok_month) <- month.abb
  era5_ok_month        <- stackApply(era5_ok, index_months, fun = sum)
  names(era5_ok_month) <- month.abb
  molo_ok_month        <- stackApply(molo_ok, index_months, fun = sum)
  names(molo_ok_month) <- month.abb
  bola_ok_month        <- stackApply(bola_ok, index_months, fun = sum)
  names(bola_ok_month) <- month.abb
  print(paste0("+++OK"))
###

###
  # Creating Season data 
  print(paste0("+++Creating Season data"))
  # EOBS
  eobs_ok_MAM <- subset(eobs_ok_month,3:5)
  eobs_ok_JJA <- subset(eobs_ok_month,6:8)
  eobs_ok_SON <- subset(eobs_ok_month,9:11)
  eobs_ok_JAN <- subset(eobs_ok_month,1)
  eobs_ok_FEB <- subset(eobs_ok_month,2)
  eobs_ok_DEC <- subset(eobs_ok_month,12)
  eobs_ok_DJF <- brick(eobs_ok_JAN,eobs_ok_FEB,eobs_ok_DEC)
  # ERA5
  era5_ok_MAM <- subset(era5_ok_month,3:5)
  era5_ok_JJA <- subset(era5_ok_month,6:8)
  era5_ok_SON <- subset(era5_ok_month,9:11)
  era5_ok_JAN <- subset(era5_ok_month,1)
  era5_ok_FEB <- subset(era5_ok_month,2)
  era5_ok_DEC <- subset(era5_ok_month,12)
  era5_ok_DJF <- brick(era5_ok_JAN,era5_ok_FEB,era5_ok_DEC)
  # MOLOCH
  molo_ok_MAM <- subset(molo_ok_month,3:5)
  molo_ok_JJA <- subset(molo_ok_month,6:8)
  molo_ok_SON <- subset(molo_ok_month,9:11)
  molo_ok_JAN <- subset(molo_ok_month,1)
  molo_ok_FEB <- subset(molo_ok_month,2)
  molo_ok_DEC <- subset(molo_ok_month,12)
  molo_ok_DJF <- brick(molo_ok_JAN,molo_ok_FEB,molo_ok_DEC)
  # BOLAM
  bola_ok_MAM <- subset(bola_ok_month,3:5)
  bola_ok_JJA <- subset(bola_ok_month,6:8)
  bola_ok_SON <- subset(bola_ok_month,9:11)
  bola_ok_JAN <- subset(bola_ok_month,1)
  bola_ok_FEB <- subset(bola_ok_month,2)
  bola_ok_DEC <- subset(bola_ok_month,12)
  bola_ok_DJF <- brick(bola_ok_JAN,bola_ok_FEB,bola_ok_DEC)
  print(paste0("+++OK"))
###

###
  # Extracting data
  print(paste0("+++Extracting data MAM..."))
  eobs_values <- as.data.frame(extract(eobs_ok_MAM, pointCoordinates))
  era5_values <- as.data.frame(extract(era5_ok_MAM*sf_era5, pointCoordinates))
  molo_values <- as.data.frame(extract(molo_ok_MAM, pointCoordinates))
  bola_values <- as.data.frame(extract(bola_ok_MAM, pointCoordinates))
  eobs_values_MAM <- c(eobs_values_MAM,eobs_values[,1])
  era5_values_MAM <- c(era5_values_MAM,era5_values[,1])
  molo_values_MAM <- c(molo_values_MAM,molo_values[,1])
  bola_values_MAM <- c(bola_values_MAM,bola_values[,1])
  print(paste0("+++OK"))
###
###
  # Extracting data
  print(paste0("+++Extracting data JJA..."))
  eobs_values <- as.data.frame(extract(eobs_ok_JJA, pointCoordinates))
  era5_values <- as.data.frame(extract(era5_ok_JJA*sf_era5, pointCoordinates))
  molo_values <- as.data.frame(extract(molo_ok_JJA, pointCoordinates))
  bola_values <- as.data.frame(extract(bola_ok_JJA, pointCoordinates))
  eobs_values_JJA <- c(eobs_values_JJA,eobs_values[,1])
  era5_values_JJA <- c(era5_values_JJA,era5_values[,1])
  molo_values_JJA <- c(molo_values_JJA,molo_values[,1])
  bola_values_JJA <- c(bola_values_JJA,bola_values[,1])
  print(paste0("+++OK"))
###
###
  # Extracting data
  print(paste0("+++Extracting data SON..."))
  eobs_values <- as.data.frame(extract(eobs_ok_SON, pointCoordinates))
  era5_values <- as.data.frame(extract(era5_ok_SON*sf_era5, pointCoordinates))
  molo_values <- as.data.frame(extract(molo_ok_SON, pointCoordinates))
  bola_values <- as.data.frame(extract(bola_ok_SON, pointCoordinates))
  eobs_values_SON <- c(eobs_values_SON,eobs_values[,1])
  era5_values_SON <- c(era5_values_SON,era5_values[,1])
  molo_values_SON <- c(molo_values_SON,molo_values[,1])
  bola_values_SON <- c(bola_values_SON,bola_values[,1])
  print(paste0("+++OK"))
###
###
  # Extracting data
  print(paste0("+++Extracting data DJF..."))
  eobs_values <- as.data.frame(extract(eobs_ok_DJF, pointCoordinates))
  era5_values <- as.data.frame(extract(era5_ok_DJF*sf_era5, pointCoordinates))
  molo_values <- as.data.frame(extract(molo_ok_DJF, pointCoordinates))
  bola_values <- as.data.frame(extract(bola_ok_DJF, pointCoordinates))
  eobs_values_DJF <- c(eobs_values_DJF,eobs_values[,1])
  era5_values_DJF <- c(era5_values_DJF,era5_values[,1])
  molo_values_DJF <- c(molo_values_DJF,molo_values[,1])
  bola_values_DJF <- c(bola_values_DJF,bola_values[,1])
  print(paste0("+++OK"))
###

###
  # End of my_year loop
  print(paste0("___End of ",my_year,"___"))
}
# END OF MAIN LOOP OVER YEARS
#############################

eobs_values_MAM <- replace(eobs_values_MAM,eobs_values_MAM<cutoff_NA,NA)
era5_values_MAM <- replace(era5_values_MAM,era5_values_MAM<cutoff_NA,NA) 
molo_values_MAM <- replace(molo_values_MAM,molo_values_MAM<cutoff_NA,NA)
bola_values_MAM <- replace(bola_values_MAM,bola_values_MAM<cutoff_NA,NA)
eobs_values_JJA <- replace(eobs_values_JJA,eobs_values_JJA<cutoff_NA,NA)
era5_values_JJA <- replace(era5_values_JJA,era5_values_JJA<cutoff_NA,NA) 
molo_values_JJA <- replace(molo_values_JJA,molo_values_JJA<cutoff_NA,NA)
bola_values_JJA <- replace(bola_values_JJA,bola_values_JJA<cutoff_NA,NA)
eobs_values_SON <- replace(eobs_values_SON,eobs_values_SON<cutoff_NA,NA)
era5_values_SON <- replace(era5_values_SON,era5_values_SON<cutoff_NA,NA) 
molo_values_SON <- replace(molo_values_SON,molo_values_SON<cutoff_NA,NA)
bola_values_SON <- replace(bola_values_SON,bola_values_SON<cutoff_NA,NA)
eobs_values_DJF <- replace(eobs_values_DJF,eobs_values_DJF<cutoff_NA,NA)
era5_values_DJF <- replace(era5_values_DJF,era5_values_DJF<cutoff_NA,NA) 
molo_values_DJF <- replace(molo_values_DJF,molo_values_DJF<cutoff_NA,NA)
bola_values_DJF <- replace(bola_values_DJF,bola_values_DJF<cutoff_NA,NA)
###
print(paste0("Plotting Taylor diagrams MAM"))
# Set output eps file for taylor diagram
setEPS()
file_taylor_eps <- paste0(my_figs_dir,"/tp_Taylor_Diag_vs_EOBS_norm",TD_normal,"_MAM.eps")
postscript(file_taylor_eps)
# Taylor Diagrams
obs<-eobs_values_MAM
# ERA5
pre <- era5_values_MAM
taylor.diagram(obs, pre, add=F, col="red", ref.sd=TRUE, normalize=TD_normal, main="Total precipitation - MAM", grad.corr.lines=c(0.4,0.6,0.8), ngamma=5, pch=19, pcex=1.5)
# MOLOCH
pre <- molo_values_MAM
taylor.diagram(obs, pre, add=T, col="orange", ref.sd=TRUE, normalize=TD_normal, pch=19, pcex=1.5)
# BOLAM
pre <- bola_values_MAM
taylor.diagram(obs, pre, add=T, col="blue", ref.sd=TRUE, normalize=TD_normal, pch=19, pcex=1.5)
# Draw the legend on the plot and save it
legend("topright",legend=c("ERA5-Land","MOLOCH","BOLAM"), pch=19, col=c("red","orange","blue"))
dev.off()
print(paste0("OK"))
###
print(paste0("Plotting Taylor diagrams JJA"))
# Set output eps file for taylor diagram
setEPS()
file_taylor_eps <- paste0(my_figs_dir,"/tp_Taylor_Diag_vs_EOBS_norm",TD_normal,"_JJA.eps")
postscript(file_taylor_eps)
# Taylor Diagrams
obs<-eobs_values_JJA
# ERA5
pre <- era5_values_JJA
taylor.diagram(obs, pre, add=F, col="red", ref.sd=TRUE, normalize=TD_normal, main="Total precipitation - JJA", grad.corr.lines=c(0.4,0.6,0.8), ngamma=5, pch=19, pcex=1.5)
# MOLOCH
pre <- molo_values_JJA
taylor.diagram(obs, pre, add=T, col="orange", ref.sd=TRUE, normalize=TD_normal, pch=19, pcex=1.5)
# BOLAM
pre <- bola_values_JJA
taylor.diagram(obs, pre, add=T, col="blue", ref.sd=TRUE, normalize=TD_normal, pch=19, pcex=1.5)
# Draw the legend on the plot and save it
legend("topright",legend=c("ERA5-Land","MOLOCH","BOLAM"), pch=19, col=c("red","orange","blue"))
dev.off()
print(paste0("OK"))
###
print(paste0("Plotting Taylor diagrams SON"))
# Set output eps file for taylor diagram
setEPS()
file_taylor_eps <- paste0(my_figs_dir,"/tp_Taylor_Diag_vs_EOBS_norm",TD_normal,"_SON.eps")
postscript(file_taylor_eps)
# Taylor Diagrams
obs<-eobs_values_SON
# ERA5
pre <- era5_values_SON
taylor.diagram(obs, pre, add=F, col="red", ref.sd=TRUE, normalize=TD_normal, main="Total precipitation - SON", grad.corr.lines=c(0.4,0.6,0.8), ngamma=5, pch=19, pcex=1.5)
# MOLOCH
pre <- molo_values_SON
taylor.diagram(obs, pre, add=T, col="orange", ref.sd=TRUE, normalize=TD_normal, pch=19, pcex=1.5)
# BOLAM
pre <- bola_values_SON
taylor.diagram(obs, pre, add=T, col="blue", ref.sd=TRUE, normalize=TD_normal, pch=19, pcex=1.5)
# Draw the legend on the plot and save it
legend("topright",legend=c("ERA5-Land","MOLOCH","BOLAM"), pch=19, col=c("red","orange","blue"))
dev.off()
print(paste0("OK"))
###
print(paste0("Plotting Taylor diagrams DJF"))
# Set output eps file for taylor diagram
setEPS()
file_taylor_eps <- paste0(my_figs_dir,"/tp_Taylor_Diag_vs_EOBS_norm",TD_normal,"_DJF.eps")
postscript(file_taylor_eps)
# Taylor Diagrams
obs<-eobs_values_DJF
# ERA5
pre <- era5_values_DJF
taylor.diagram(obs, pre, add=F, col="red", ref.sd=TRUE, normalize=TD_normal, main="Total precipitation - DJF", grad.corr.lines=c(0.4,0.6,0.8), ngamma=5, pch=19, pcex=1.5)
# MOLOCH
pre <- molo_values_DJF
taylor.diagram(obs, pre, add=T, col="orange", ref.sd=TRUE, normalize=TD_normal, pch=19, pcex=1.5)
# BOLAM
pre <- bola_values_DJF
taylor.diagram(obs, pre, add=T, col="blue", ref.sd=TRUE, normalize=TD_normal, pch=19, pcex=1.5)
# Draw the legend on the plot and save it
legend("topright",legend=c("ERA5-Land","MOLOCH","BOLAM"), pch=19, col=c("red","orange","blue"))
dev.off()
print(paste0("OK"))

# Quit
quit()

