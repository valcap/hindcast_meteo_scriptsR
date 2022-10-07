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
eobs_values_ts  <- NULL 
era5_values_ts  <- NULL
molo_values_ts  <- NULL
bola_values_ts  <- NULL
add_offset      <- 273.15
var             <- 'tx'

# Working directory
setwd(my_working_dir)

# Load the point locations
my_anag_file     <- paste0(my_working_dir,"/E-OBS/points050_moloch_extent.csv")
anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
pointCoordinates <- anag_full
coordinates(pointCoordinates)= ~ LONGITUDE+LATITUDE
print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))
#setEPS()
#postscript(paste0(my_figs_dir,"/virtual_rain-gauges_EOBS.eps"))
#plot(pointCoordinates, main=paste0("# of rain-gauges ",length(pointCoordinates)))
#plot(extent(pointCoordinates), add=TRUE, col='red', lwd=1)
#dev.off()

######################
# MAIN LOOP OVER YEARS
for (my_year in vec_years) {
  if (my_year==2001) {
    print(paste0("___Skipping ",my_year,"___"))
    next
  }
  print(paste0("___Running year ",my_year,"___"))

###
  # Loading raster file names
  print(paste0("+++Reading raster files..."))
  # Raster files definition
  file_eobs   <- paste0("E-OBS/t2m/E-OBS_day_",var,"_",my_year,".nc")
  file_era5   <- paste0("ERA5-Land/t2m/ERA5-Land_",my_year,"-daystat_",var,".nc")
  file_molo   <- paste0("work/moloch/moloch_",my_year,"_daystat_",var,".nc")
  file_bola   <- paste0("work/bolam/bolam_",my_year,"_daystat_",var,".nc") 

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
#  era5_resam          <- resample(era5_crop, eobs_crop, method="ngb")
  era5_resam          <- era5_crop
  # Crop BOLAM grid on my_extent and regrid (if needed)
  bola_crop           <- crop(bola_orig, my_extent)
#  bola_resam          <- resample(bola_crop, eobs_crop, method="bilinear")
  bola_resam          <- bola_crop
  # Crop MOLOCH grid on my_extent and regrid (if needed)
  molo_crop           <- crop(molo_orig, my_extent)
#  molo_resam          <- resample(molo_crop, eobs_crop, method="bilinear")
  molo_resam          <- molo_crop
  print(paste0("+++OK"))
###

###
  # select only data from the current year
  print(paste0("+++Adjusting name conventions..."))
  names(eobs_resam)     <- as.Date(names(eobs_resam),format="X%Y.%m.%d")
  names(era5_resam)     <- as.Date(names(era5_resam),format="X%Y.%m.%d.%H.%M.%S")
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
  # Extracting data
  print(paste0("+++Extracting data..."))
  eobs_values    <- as.vector(extract(eobs_ok, pointCoordinates))
  era5_values    <- as.vector(extract(era5_ok-add_offset, pointCoordinates))
  molo_values    <- as.vector(extract(molo_ok-add_offset, pointCoordinates))
  bola_values    <- as.vector(extract(bola_ok-add_offset, pointCoordinates))
  eobs_values_ts <- c(eobs_values_ts,eobs_values)
  era5_values_ts <- c(era5_values_ts,era5_values)
  molo_values_ts <- c(molo_values_ts,molo_values)
  bola_values_ts <- c(bola_values_ts,bola_values)
  print(paste0("+++OK"))
###

###
  # End of my_year loop
  print(paste0("___End of ",my_year,"___"))
}
# END OF MAIN LOOP OVER YEARS
#############################

###
print(paste0("Plotting Taylor diagrams"))
# Set output eps file for taylor diagram
setEPS()
file_taylor_eps <- paste0(my_figs_dir,"/",var,"_Taylor_Diag_vs_EOBS_norm",TD_normal,".eps")
postscript(file_taylor_eps)
# Taylor Diagrams
my_title<-NULL
if (var == 'tx') {
  my_title<-'2-metre Maximum Temperature'
}
if (var == 'tn') {
  my_title<-'2-metre Minimum Temperature'
}
obs <- eobs_values_ts
# ERA5
pre <- era5_values_ts
taylor.diagram(obs, pre, add=F, col="red", ref.sd=TRUE, normalize=TD_normal, main=my_title, grad.corr.lines=c(0.4,0.6,0.8), ngamma=5, pch=19, pcex=1.5)
# MOLOCH
pre <- molo_values_ts
taylor.diagram(obs, pre, add=T, col="orange", ref.sd=TRUE, normalize=TD_normal, pch=19, pcex=1.5)
# BOLAM
pre <- bola_values_ts
taylor.diagram(obs, pre, add=T, col="blue", ref.sd=TRUE, normalize=TD_normal, pch=19, pcex=1.5)
# Draw the legend on the plot and save it
legend("topright",legend=c("ERA5-Land","MOLOCH","BOLAM"), pch=19, col=c("red","orange","blue"))
dev.off()
print(paste0("OK"))

# Quit
quit()
