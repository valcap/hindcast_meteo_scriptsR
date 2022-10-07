# Load libraries
library(raster)
library(ncdf4)
library(RColorBrewer)
library(viridis)
library(rgdal)
# library(verification)
# library(Metrics)
# library(plotrix)
# library(lubridate)
# library(dplyr)
# add here any other packages you need

# Global variables
u               <- 'capecchi'
my_working_dir  <- paste0("/home/",u,"/art_hindcast_meteo")
my_figs_dir     <- paste0(my_working_dir,"/figs")
my_shp_dir      <- paste0("/home/",u,"/Sync/shp")
my_pal          <- brewer.pal(n = 2, name = "RdBu")
my_limits       <- c(0,3000)

# Working directory
setwd(my_working_dir)

my_coastline_lr <- readOGR(
  dsn= paste0(my_shp_dir,"/gshhg-shp-2.3.7/GSHHS_shp/l/") ,
  layer="GSHHS_l_L1",
  verbose=FALSE)
add_shp=function(){plot(my_coastline_lr, bg="transparent", add=TRUE,lwd=0.5)}

# Allocate raster
obse_dataset    <- "ARCIS"
#file_obse       <- paste0("DATA/",obse_dataset,"/tp/gripho-v2-3km_year_tp_AVG2001-2016.nc")
file_obse       <- paste0("DATA/",obse_dataset,"/tp/ARCIS_year_tp_AVG1981-2015.nc")
obse_orig       <- NULL
obse_orig       <- brick(file_obse)

###
# Plot percentiles
print(paste0("+++Producing tp percentile plot..."))
setEPS()
file_taylor_eps <- paste0(my_figs_dir,"/tp_",obse_dataset,"_percentile-mm.eps")
postscript(file_taylor_eps, horizontal = FALSE, onefile = FALSE, paper = "special", height = 6, width = 9)
par(mfrow=c(2,3))
for (thres in c(700,900,1100,1400,1600,1800)) {
  obse_cutoff      <- obse_orig > thres
  if (thres==700) {plot_legend=TRUE} else {plot_legend=FALSE}
  plot(obse_cutoff, main=paste0("Mean annual precipitation > ",thres," mm"), col = c("white","blue"), axes=FALSE, legend.width=3, legend=plot_legend, box=F, addfun=add_shp)
}
dev.off()
print(paste0("+++OK"))
###

# Byebye
print(paste0("Byebye"))
#quit()

