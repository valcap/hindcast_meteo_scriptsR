# Load libraries
library(raster)
library(ncdf4)
library(verification)
library(Metrics)
library(plotrix)
library(lubridate)
library(dplyr)
library(hydroGOF)
library(greenbrown)

# Global variables
u               <- 'capecchi'
my_working_dir  <- paste0("/home/",u,"/art_hindcast_meteo")
my_figs_dir     <- paste0(my_working_dir,"/figs")
my_output_dir   <- paste0(my_working_dir,"/res")
my_shp_dir      <- paste0("/home/",u,"Sync/shp")
TD_normal       <- FALSE
sf_era5         <- 1000
max_sample_size <- 500
ini_year        <- 1981
end_year        <- 2019
vec_years       <- seq(ini_year,end_year)
TD_normal       <- 'FALSE'
my_vector_obse  <- NULL 
my_vector_era5  <- NULL
my_vector_molo  <- NULL
my_vector_bola  <- NULL

# Working directory
setwd(my_working_dir)

# Taylor Diagram
print(paste0("Plotting Taylor diagrams"))
setEPS()
file_taylor_eps <- paste0(my_figs_dir,"/tp_Taylor_Diag_vs_SCIA_years.eps")
postscript(file_taylor_eps)

######################
# MAIN LOOP OVER YEARS
for (my_year in vec_years) {
  print(paste0("___Running year ",my_year,"___"))
  if (my_year==2001) {next}

###
  # Load the point locations...
  my_anag_file     <- paste0(my_working_dir,"/DATA/SCIA/SCIA_TotAnnPre_",my_year,".csv")
  anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
  pointCoordinates <- anag_full
  coordinates(pointCoordinates)= ~ lon+lat
###

###
  # Loading raster file names
  print(paste0("+++Reading raster files..."))
  # Raster files definition
  file_era5   <- paste0("DATA/ERA5-Land/tp/ERA5-Land_year_tp_",my_year,".nc")
  file_molo   <- paste0("DATA/MOLOCH/tp/moloch_year_tp_",my_year,".nc")
  file_bola   <- paste0("DATA/BOLAM/tp/bolam_year_tp_",my_year,".nc") 
  # Allocate rasters
  era5_orig   <- NULL
  molo_orig   <- NULL
  bola_orig   <- NULL
  era5_orig   <- brick(file_era5)
  molo_orig   <- brick(file_molo)
  bola_orig   <- brick(file_bola)
  print(paste0("+++OK"))
###

###
  # Get extent on the basis of the point coordinates
  print(paste0("+++Cropping extents and resampling rasters..."))
  my_extent           <- extent(molo_orig)
  # Crop ERA5-Land grid on my_extent and regrid (if needed)
  era5_crop           <- crop(era5_orig, my_extent)
  era5_resam          <- era5_crop*sf_era5
  # Crop BOLAM grid on my_extent and regrid (if needed)
  bola_crop           <- crop(bola_orig, my_extent)
#  bola_resam          <- resample(bola_crop, era5_crop, method="bilinear")
  bola_resam          <- bola_crop
  # Crop MOLOCH grid on my_extent and regrid (if needed)
  molo_crop           <- crop(molo_orig, my_extent)
# molo_resam          <- resample(molo_crop, era5_crop, method="bilinear")
  molo_resam          <- molo_crop
  print(paste0("+++OK"))
###

###
  # select only data from the current year
  print(paste0("+++Adjusting name conventions..."))
  names(era5_resam)     <- as.Date(names(era5_orig), format="X%Y.%m.%d.%H.%M.%S")
  names(bola_resam)     <- as.Date(names(bola_resam),format="X%Y.%m.%d.%H.%M.%S")
  names(molo_resam)     <- names(bola_resam)
  print(paste0("+++OK"))
###

###
  # Extracting data
  print(paste0("+++Extracting data..."))
  # Selecting a subset of pointCoordinates to avoid overflow over the time series vec_years
#  sample_size <- min(max_sample_size,length(pointCoordinates$val))
#  my_sample_coordinates <- pointCoordinates[sample(1:length(pointCoordinates),sample_size),]
  era5_values <- as.integer(as.data.frame(extract(era5_resam, pointCoordinates))[,1])
  molo_values <- as.integer(as.data.frame(extract(molo_resam, pointCoordinates))[,1])
  bola_values <- as.integer(as.data.frame(extract(bola_resam, pointCoordinates))[,1])
  print(paste0("+++OK"))
  # Observed values as data frame
  obse_values <- as.integer(pointCoordinates$val)
###

###
  # Calculating standard skill scores
  # SCIA
  my_sd=sprintf("%.2f",sd(obse_values,na.rm=T))
  print(paste("SCIA",my_sd,sep=";"))
  print(paste("VER_DATASET","ST.DEV","CORR","C-RMSE","RMSE","ME","MBIAS","OBIAS",sep=" & "))
  # ERA5-Land
  my_sd=sprintf("%.0f",sd(era5_values,na.rm=T))
  my_cor=sprintf("%.2f",cor(obse_values, era5_values,use="complete.obs"))
  my_ver=verify(obs=obse_values, pred=era5_values,frcst.type="cont",obs.type="cont")
  my_rmse=as.numeric(sprintf("%.0f",sqrt(my_ver$MSE)))
  my_me=sprintf("%.0f",(my_ver$ME))
  my_mbias=sprintf("%.2f",mean(era5_values,na.rm=T)/mean(obse_values,na.rm=T))
  my_overall_bias=as.numeric(sprintf("%.0f",mean(era5_values,na.rm=T)-mean(obse_values,na.rm=T)))
  my_crmse=as.numeric(sprintf("%.0f",sqrt((my_rmse*my_rmse)-(my_overall_bias*my_overall_bias))))
  print(paste("VER_ERA5-Land",my_sd,my_cor,my_crmse,my_rmse,my_me,my_mbias,my_overall_bias,sep=" & "))
  # MOLOCH
  my_sd=sprintf("%.0f",sd(molo_values,na.rm=T))
  my_cor=sprintf("%.2f",cor(obse_values, molo_values,use="complete.obs"))
  my_ver=verify(obs=obse_values, pred=molo_values,frcst.type="cont",obs.type="cont")
  my_rmse=as.numeric(sprintf("%.0f",sqrt(my_ver$MSE)))
  my_me=sprintf("%.0f",(my_ver$ME))
  my_mbias=sprintf("%.2f",mean(molo_values,na.rm=T)/mean(obse_values,na.rm=T))
  my_overall_bias=as.numeric(sprintf("%.0f",mean(molo_values,na.rm=T)-mean(obse_values,na.rm=T)))
  my_crmse=as.numeric(sprintf("%.0f",sqrt((my_rmse*my_rmse)-(my_overall_bias*my_overall_bias))))
  print(paste("VER_MOLOCH",my_sd,my_cor,my_crmse,my_rmse,my_me,my_mbias,my_overall_bias,sep=" & "))
  # BOLAM
  my_sd=sprintf("%.0f",sd(bola_values,na.rm=T))
  my_cor=sprintf("%.2f",cor(obse_values, bola_values,use="complete.obs"))
  my_ver=verify(obs=obse_values, pred=bola_values,frcst.type="cont",obs.type="cont")
  my_rmse=as.numeric(sprintf("%.0f",sqrt(my_ver$MSE)))
  my_me=sprintf("%.0f",(my_ver$ME))
  my_mbias=sprintf("%.2f",mean(bola_values,na.rm=T)/mean(obse_values,na.rm=T))
  my_overall_bias=as.numeric(sprintf("%.0f",mean(bola_values,na.rm=T)-mean(obse_values,na.rm=T)))
  my_crmse=as.numeric(sprintf("%.0f",sqrt((my_rmse*my_rmse)-(my_overall_bias*my_overall_bias))))
  print(paste("VER_BOLAM",my_sd,my_cor,my_crmse,my_rmse,my_me,my_mbias,my_overall_bias,sep=" & "))
###

###
  # Storing data
  my_vector_obse <- as.vector(c(my_vector_obse,obse_values))
  my_vector_era5 <- as.vector(c(my_vector_era5,era5_values))
  my_vector_molo <- as.vector(c(my_vector_molo,molo_values))
  my_vector_bola <- as.vector(c(my_vector_bola,bola_values))
  print(paste("LENGTH",length(my_vector_obse),length(my_vector_era5),length(my_vector_molo),length(my_vector_bola),sep=" "))
###

  # End of my_year loop
  print(paste0("___End of ",my_year,"___"))
###
}

# Calculates KGE goodness-of-fit measure (Knoben et al 2019)
# Per qualche motivo a me ignoto, non fornisce i risultati delle librerie (che assumiamo invece corretti)
# kge_era5 <- 1-sqrt((cor(my_vector_obse,my_vector_era5,use="complete.obs")-1)^2+((sd(my_vector_era5,na.rm=T))/(sd(my_vector_obse,na.rm=T))-1)^2+((mean(my_vector_era5,na.rm=T))/(sd(my_vector_obse,na.rm=T))-1)^2)
# kge_molo <- 1-sqrt((cor(my_vector_obse,my_vector_molo,use="complete.obs")-1)^2+((sd(my_vector_molo,na.rm=T))/(sd(my_vector_obse,na.rm=T))-1)^2+((mean(my_vector_molo,na.rm=T))/(sd(my_vector_obse,na.rm=T))-1)^2)
# kge_bola <- 1-sqrt((cor(my_vector_obse,my_vector_bola,use="complete.obs")-1)^2+((sd(my_vector_bola,na.rm=T))/(sd(my_vector_obse,na.rm=T))-1)^2+((mean(my_vector_bola,na.rm=T))/(sd(my_vector_obse,na.rm=T))-1)^2)
# kge_era5
# kge_molo
# kge_bola

# KGE with the greenbrown package
# greenbrown::KGE(my_vector_era5, my_vector_obse)
# greenbrown::KGE(my_vector_molo, my_vector_obse)
# greenbrown::KGE(my_vector_bola, my_vector_obse)
# KGE with the hydroGOF package
hydroGOF::KGE(my_vector_era5, my_vector_obse)
hydroGOF::KGE(my_vector_molo, my_vector_obse)
hydroGOF::KGE(my_vector_bola, my_vector_obse)


# quit()
# 
# ##########################################################################
# ##########################################################################
# # Scatter plots
# plot(eobs_values_ts, era5_values_ts,
#      xlim=c(250,3000) , ylim=c(250,3000),
#      pch=19,
#      cex=1,
#      col="red",
#      xlab="E-OBS",
#      ylab="ERA5-Land",
#      main="(b) ERA5-Land"
# )
# lines(x = c(250,3000), y = c(250,3000))
# #abline(coef = c(0,1))
# 
# plot(eobs_values_ts, molo_values_ts,
#      xlim=c(250,3000) , ylim=c(250,3000),
#      pch=19,
#      cex=1,
#      col="orange",
#      xlab="E-OBS", 
#      ylab="MOLOCH",
#      main="(c) MOLOCH"
# )
# lines(x = c(250,3000), y = c(250,3000))
# 
# plot(eobs_values_ts, bola_values_ts,
#      xlim=c(250,3000) , ylim=c(250,3000),
#      pch=19,
#      cex=1,
#      col="blue",
#      xlab="E-OBS",
#      ylab="BOLAM",
#      main="(d) BOLAM"
# )
# lines(x = c(250,3000), y = c(250,3000))
# 
# dev.off()
# print(paste0("OK"))
# ##########################################################################
# ##########################################################################
