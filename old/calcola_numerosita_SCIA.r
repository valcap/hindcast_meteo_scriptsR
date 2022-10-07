# Purpose
print(paste0("Calculating performance diagrams for total annual precipitation"))

# Load libraries
library(raster)
# add here any other packages you need

# Global variables
u               <- 'valcap'
my_working_dir  <- paste0("/home/",u,"/art_hindcast_meteo")
my_figs_dir     <- paste0(my_working_dir,"/figs")
ini_year        <- 1981
end_year        <- 2019
vec_years       <- seq(ini_year,end_year)
vec_stations    <- NULL

# Working directory
setwd(my_working_dir)

######################
# MAIN LOOP OVER YEARS
for (my_year in vec_years) {

# Load the point locations...
  my_anag_file     <- paste0(my_working_dir,"/DATA/SCIA/SCIA_TotAnnPre_",my_year,".csv")
  anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
  pointCoordinates <- anag_full
  coordinates(pointCoordinates)= ~ lon+lat
  n_stat           <- as.integer(length(pointCoordinates))
  vec_stations     <- as.vector(c(vec_stations,n_stat))
}
# END OF MAIN LOOP OVER YEARS
#############################

# Plot
setEPS()
postscript(paste0(my_figs_dir,"/SCIA_rain-gauges.eps"))
plot(vec_stations~vec_years,type="p",pch=19,ylim=c(1500,5000), ylab="NUMBER OF RAIN-GAUGES", xlab="YEAR",main="SCIA", xaxt="n")
axis(side=1, at=seq(ini_year,end_year,5))
lines(vec_stations~vec_years)
dev.off()

# Byebye
print(paste0("Byebye"))
#quit()

