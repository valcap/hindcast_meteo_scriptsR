# Purpose
print(paste0("Calculating performance diagrams for total annual precipitation"))

# Load libraries
library(raster)
# add here any other packages you need

# Global variables
u               <- 'capecchi'
my_working_dir  <- paste0("/home/",u,"/art_hindcast_meteo")
my_figs_dir     <- paste0(my_working_dir,"/figs")
ini_year        <- 1979
end_year        <- 2020
vec_years       <- seq(ini_year,end_year)
my_vector_perc  <- c(0.66,0.75,0.90,0.95,0.98,0.99)
#my_vector_perc  <- c(0.10,0.25,0.50,0.75,0.90,0.95)
#my_vector_perc  <- c(0.25,0.33,0.50,0.66,0.75,0.90)
my_matrix       <- NULL

# Working directory
setwd(my_working_dir)

######################
# MAIN LOOP OVER YEARS
for (my_year in vec_years) {
  print(paste0("___Running year ",my_year,"___"))

# Load the point locations...
  my_anag_file     <- paste0(my_working_dir,"/DATA/SCIA/SCIA_TotAnnPre_",my_year,".csv")
  anag_full        <- read.csv(my_anag_file, sep=",", header=TRUE)
  pointCoordinates <- anag_full
  coordinates(pointCoordinates)= ~ lon+lat
  print(paste0("Total number of points (rain-gauge locations) is ",length(pointCoordinates)))
  my_perc <- as.integer(as.vector(quantile(pointCoordinates$val, probs = my_vector_perc)))
  my_matrix <- cbind(my_matrix,my_perc)
#  print(paste0("P ",my_perc))
  my_perc <- NULL
  print(paste0("___End of ",my_year,"___"))
}
# END OF MAIN LOOP OVER YEARS
#############################
#print(paste0("P ",as.integer(rowMeans(my_matrix))))
print(as.integer(rowMeans(my_matrix)))

# Byebye
print(paste0("Byebye"))
quit()

