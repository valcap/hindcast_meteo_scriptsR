# Purpose
print(paste0("Calculating performance diagrams for total annual precipitation"))

# Load libraries
library(raster)
library(ncdf4)
library(verification)
library(Metrics)
library(plotrix)
library(lubridate)
library(dplyr)
# add here any other packages you need

# Global variables
u               <- 'capecchi'
my_working_dir  <- paste0("/home/",u,"/art_hindcast_meteo")
my_figs_dir     <- paste0(my_working_dir,"/figs")
ini_year        <- 1983
end_year        <- 2019
# ATT: per motivi di integer overflow (vale a dire troppe stazioni e serie temporale troppo lunga)
#      non \'e possibile fare la verifica su tutto il periodo, ma bisogna scegliere tra queste due opzioni:
#      1. si riduce il numero di stazioni
#      2. si riduce il numero di anni (opzione scelta, serie 1983-2019, 37 anni)
vec_years       <- seq(ini_year,end_year)
sf_era5         <- 1000
my_vector_scia  <- NULL
my_vector_era5  <- NULL
my_vector_bola  <- NULL
my_vector_molo  <- NULL
# PERCENTILI c(0.25,0.33,0.50,0.66,0.75,0.90)
my_vector_thres <- c(675,740,875,1040,1170,1550)
my_vector_perce <- c("25th Percentile","33th Percentile","50th Percentile","66th Percentile","75th Percentile","90th Percentile")
prog            <- 1
# PERCENTILI c(0.10,0.25,0.50,0.75,0.90,0.95)
# my_vector_thres <- c(530,675,875,1170,1550,1800)
# my_vector_perce <- c("10th Percentile","25th Percentile","50th Percentile","75th Percentile","90th Percentile","95th Percentile")
# prog            <- 2
# PERCENTILI c(0.66,0.75,0.90,0.95,0.98,0.99)
# my_vector_thres <- c(1040,1170,1550,1800,2115,2325)
# my_vector_perce <- c("66th Percentile","75th Percentile","90th Percentile","95th Percentile","98th Percentile","99th Percentile")
# prog            <- 3

# Working directory
setwd(my_working_dir)

######################
# MAIN LOOP OVER YEARS
for (my_year in vec_years) {
  print(paste0("___Running year ",my_year,"___"))
  if (my_year==2001) {
    next
  }
  # Load the point locations...
  my_anag_file     <- paste0(my_working_dir,"/DATA/SCIA/SCIA_TotAnnPre_",my_year,".csv")
  anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
  pointCoordinates <- anag_full
  coordinates(pointCoordinates)= ~ lon+lat
  print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))
  print(paste0(""))

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
#  molo_resam          <- resample(molo_crop, era5_crop, method="bilinear")
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
  # Oh yeah
  print(paste0("+++Extracting data..."))
  era5_values <- as.data.frame(extract(era5_resam, pointCoordinates))[,1]
  molo_values <- as.data.frame(extract(molo_resam, pointCoordinates))[,1]
  bola_values <- as.data.frame(extract(bola_resam, pointCoordinates))[,1]
  #
#  era5_values <- replace(era5_values,era5_values<cutoff_NA, NA)
#  molo_values <- replace(molo_values,molo_values<cutoff_NA, NA)
#  bola_values <- replace(bola_values,bola_values<cutoff_NA, NA)
  print(paste0("+++OK"))
  # Observed values as data frame
  obse_values <- pointCoordinates$val
 
###
  # Storing data
  my_vector_scia <- as.vector(c(my_vector_scia,obse_values))
  my_vector_era5 <- as.vector(c(my_vector_era5,era5_values))
  my_vector_molo <- as.vector(c(my_vector_molo,molo_values))
  my_vector_bola <- as.vector(c(my_vector_bola,bola_values))
  print(paste("LENGTH",length(my_vector_scia),length(my_vector_era5),length(my_vector_molo),length(my_vector_bola),sep=" "))
###

###
  # End of my_year loop
  print(paste0("___End of ",my_year,"___"))
}
# END OF MAIN LOOP OVER YEARS
#############################

###
# PERFORMANCE DIAGRAMS
print(paste0("+++Producing performance diagram..."))
setEPS()
file_taylor_eps <- paste0(my_figs_dir,"/tp_PerfDiag_year_SCIA",prog,".eps")
postscript(file_taylor_eps, horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 9)
par(mfrow=c(2,3))
i<-1
for (thres in my_vector_thres) {
  
  scia_binary          <- replace(my_vector_scia, my_vector_scia <= thres, 0)
  scia_binary          <- replace(scia_binary, scia_binary  > thres, 1)
  era5_binary          <- replace(my_vector_era5, my_vector_era5 <= thres, 0)
  era5_binary          <- replace(era5_binary, era5_binary  > thres, 1)
  molo_binary          <- replace(my_vector_molo, my_vector_molo <= thres, 0)
  molo_binary          <- replace(molo_binary, molo_binary  > thres, 1)
  bola_binary          <- replace(my_vector_bola, my_vector_bola <= thres, 0)
  bola_binary          <- replace(bola_binary, bola_binary  > thres, 1)
  performance.diagram(main = paste0('Threshold = ',thres,' mm/year',"\n",my_vector_perce[i]), cex.main=0.99, cex.lab=1.1, cex.axis=0.99)
  my_table_stats     <- table.stats(scia_binary,era5_binary)
  points(1 - my_table_stats$FAR, my_table_stats$POD, col = "red", cex = 2.5, pch=16)
  my_table_stats     <- NULL
  my_table_stats     <- table.stats(scia_binary,bola_binary)
  points(1 - my_table_stats$FAR, my_table_stats$POD, col = "blue", cex = 2.5, pch=16)
  my_table_stats     <- NULL
  my_table_stats     <- table.stats(scia_binary,molo_binary)
  points(1 - my_table_stats$FAR, my_table_stats$POD, col = "orange", cex = 2.5, pch=16)
  my_table_stats     <- NULL
  scia_binary<-NULL; era5_binary<-NULL; molo_binary<-NULL; bola_binary<-NULL
  # Draw the legend on the plot
  if ((my_vector_perce[i]=="25th Percentile") && (prog==1)) {
    legend("topleft",legend=c("ERA5-Land","MOLOCH","BOLAM"), pch=19, col=c("red","orange","blue"), bg="white")
  }
  if ((my_vector_perce[i]=="10th Percentile") && (prog==2)) {
    legend("topleft",legend=c("ERA5-Land","MOLOCH","BOLAM"), pch=19, col=c("red","orange","blue"), bg="white")
  }
  if ((my_vector_perce[i]=="66th Percentile") && (prog==3)) {
    legend("topleft",legend=c("ERA5-Land","MOLOCH","BOLAM"), pch=19, col=c("red","orange","blue"), bg="white")
  }
  i<-i+1
}
dev.off()
print(paste0("+++OK"))
###

# Byebye
print(paste0("Byebye"))

#quit()

