# Load libraries
library("hydroGOF")
library("lubridate")

#########################################################
# General settings
u                <- 'capecchi'
my_working_dir   <- paste0("/home/",u,"/art_hindcast_meteo")
my_figs_dir      <- paste0(my_working_dir,"/figs")
my_output_dir    <- paste0(my_working_dir,"/res")
my_shp_dir       <- paste0("/home/",u,"Sync/shp")
sf_era5          <- 1000
my_vec_zone      <- c("ALPO","ALPE","PAD","TIR","SAR","SUD")
my_vec_title     <- c("(a) Western Alps","(b) Eastern Alps","(c) Po Valley","(d) Tyrrhenian Coast","(e) Sardinia","(f) Central-South Italy")

# Set working directory
setwd(my_output_dir)

for (my_season in c("JJA")) {
# Now plot
setEPS()
postscript(paste0(my_figs_dir,"/tp_DailyCycle-",my_season,"_vs_GRIPHO.eps"), width=7, height=4.5)
par(mfrow=c(2,3))
for (i in seq(1,length(my_vec_zone))) {
  zone <- my_vec_zone[i]
  # Set file naming
  grip_file <- paste0("gripho-v2-3km_",zone,"_",my_season,".csv")
  era5_file <- paste0("ERA5-Land_",zone,"_",my_season,".csv")
  molo_file <- paste0("moloch_",zone,"_",my_season,".csv")
  bola_file <- paste0("bolam_",zone,"_",my_season,".csv")

  # Read data file
  grip_data <- read.csv(grip_file,header=T,col.name="tp_cycle")
  era5_data <- read.csv(era5_file,header=T,col.name="tp_cycle")
  molo_data <- read.csv(molo_file,header=T,col.name="tp_cycle")
  bola_data <- read.csv(bola_file,header=T,col.name="tp_cycle")
#  # Now plot
#  setEPS()
#  postscript(paste0(my_figs_dir,"/tp_DailyCycle_vs_GRIPHO_",zone,".eps"), width=7, height=4.5)
  my_min           <- (min(grip_data,molo_data,bola_data))
  my_min           <- (min(grip_data,molo_data,bola_data,era5_data*1000))
  my_max           <- (max(grip_data,molo_data,bola_data,era5_data*1000))
  my_min           <- 0
  if ((zone=='ALPO') || (zone=='ALPE')) {
    my_max=0.5
  } else {
    my_max=0.2
  }
  my_hours         <- c(0,2,4,6,8,10,12,14,16,18,20,22)
  my_hours         <- c(1,3,5,7,9,11,13,15,17,19,21,23)
  plot(grip_data$tp_cycle,type="l", col="black", lwd=1, ylab="Precipitation [mm/h]", xlab="Hour [UTC]", ylim=c(my_min,my_max), main=my_vec_title[i], xaxt="n")
  axis(1, at=seq(1, 24, by = 2),labels=my_hours, col.axis="black", las=1, cex.axis=0.95)
  lines(era5_data$tp_cycle*1000, col="red", lwd=1.0)
#  points(era5_data$tp_cycle*1000, col="red", lwd=0.5, pch=19, cex=0.5)
  lines(molo_data$tp_cycle, col="orange", lwd=1.0)
#  points(molo_data$tp_cycle, col="orange", lwd=0.5, pch=19, cex=0.5)
  lines(bola_data$tp_cycle, col="blue", lwd=1.0)
#  points(bola_data$tp_cycle, col="blue", lwd=0.5, pch=19, cex=0.5)
  abline(v = c(9, 13, 17), col = c("gray","gray","gray"),
         lty = c(2, 2, 2), lwd = c(0.5, 0.5, 0.5))
  if (zone=="ALPO") {
#    legend("topleft",legend=c("GRIPHO","ERA5-Land","MOLOCH","BOLAM"), lty=1, pch=c(NA,19,19,19), col=c("black","red","orange","blue"), bg="white", lwd=c(1.0,1.0,1.0,1.0),cex=0.7)
    legend("topleft",legend=c("GRIPHO","ERA5-Land","MOLOCH","BOLAM"), lty=1, pch=c(NA,NA,NA,NA), col=c("black","red","orange","blue"), bg="white", lwd=c(1.0,1.0,1.0,1.0),cex=0.7)
  }
#  plot(grip_data$tp_cycle,type="l", col="black", lwd=4, ylab="Precipitation [mm/h]", xlab="Hour [UTC]", ylim=c(my_min,my_max), main=my_vec_title[i], xaxt="n")
#  axis(1, at=seq(1, 24, by = 2),labels=my_hours, col.axis="black", las=1, cex.axis=0.95)
#  lines(era5_data$tp_cycle*1000, col="red", lwd=2.5)
#  points(era5_data$tp_cycle*1000, col="red", lwd=2.5, pch=19, cex=1.5)
#  lines(molo_data$tp_cycle, col="orange", lwd=2.5)
#  points(molo_data$tp_cycle, col="orange", lwd=2.5, pch=19, cex=1.5)
#  lines(bola_data$tp_cycle, col="blue", lwd=2.5)
#  points(bola_data$tp_cycle, col="blue", lwd=2.5, pch=19, cex=1.5)
#  legend("topleft",legend=c("GRIPHO","ERA5-Land","MOLOCH","BOLAM"), lty=1, pch=c(NA,19,19,19), col=c("black","red","orange","blue"), bg="white", lwd=c(2.5,1.5,1.5,1.5))
#  dev.off()

  #
  # KGE with the hydroGOF package
  # method="2009" Gupta et al. 2009 - Gupta, H. V., Kling, H., Yilmaz, K. K., & Martinez, G. F. (2009).
  #               Decomposition of the mean squared error and NSE performance criteria: Implications for improving hydrological modelling. Journal of hydrology, 377(1-2), 80-91
  # method="2012" Kling et al. 2012 - Kling, H., Fuchs, M., & Paulin, M. (2012). Runoff conditions in the upper Danube basin under an ensemble of climate change scenarios. Journal of Hydrology, 424, 264-277
  kge_era5 <- hydroGOF::KGE(era5_data$tp_cycle*1000, grip_data$tp_cycle,method="2009")
  kge_molo <- hydroGOF::KGE(molo_data$tp_cycle, grip_data$tp_cycle,method="2009")
  kge_bola <- hydroGOF::KGE(bola_data$tp_cycle, grip_data$tp_cycle,method="2009")
  print(paste("KGE",zone,"ERA5-Land",sprintf("%.2f",kge_era5),"MOLOCH",sprintf("%.2f",kge_molo),"BOLAM",sprintf("%.2f",kge_bola),sep=";"),quote=FALSE)

}
dev.off()
print(paste0(""))
}
# ciao
# quit()
