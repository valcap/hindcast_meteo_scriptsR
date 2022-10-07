#!/usr/bin/env Rscript

# Load libraries
library(raster)
library(ncdf4)
library(udunits2)
library(HeatStress)
library(ClimInd)
library(biometeoR)
library(future.apply)
plan(multisession) 
library(foreach)
library(doParallel)

# Global variables
my_working_dir  <- paste0("/OCEANASTORE/progetti/spitbran/test")
#my_figs_dir     <- paste0("/OCEANASTORE/progetti/spitbran/work/figs")
#my_output_dir   <- paste0("/OCEANASTORE/progetti/spitbran/work/res")
#my_shp_dir      <- paste0("/OCEANASTORE/progetti/spitbran/work/shp")
#ini_year        <- 1981
#end_year        <- 2019
#vec_years       <- seq(ini_year,end_year)
my_year         <- 2019
prefix_fname    <- 'moloch'
my_original_raster_dir <- paste0("/OCEANASTORE/progetti/spitbran/ERA5/winds/MOLOCH/",my_year)


# File names (moloch_2018_dp2m.grib2 in /OCEANASTORE/progetti/spitbran/ERA5/winds/MOLOCH/2018/)
original_g2_files <- c(paste0(my_original_raster_dir,"/",prefix_fname,"_",my_year,"_t2m.grib2"),
                       paste0(my_original_raster_dir,"/",prefix_fname,"_",my_year,"_dp2m.grib2"),
                       paste0(my_original_raster_dir,"/",prefix_fname,"_",my_year,"_dswrf.grib2"),
                       paste0(my_original_raster_dir,"/",prefix_fname,"_",my_year,"_rh2m.grib2"),
                       paste0(my_original_raster_dir,"/",prefix_fname,"_",my_year,"_wind10m.grib2") 
                       )
g2_files <- c(paste0("test_",my_year,"_t2m.grib2"),
              paste0("test_",my_year,"_dp2m.grib2"),
              paste0("test_",my_year,"_dswrf.grib2"),
              paste0("test_",my_year,"_rh2m.grib2"),
              paste0("test_",my_year,"_wind10m.grib2") 
              )
nc_files <- gsub(".grib2",".nc",g2_files)

# Working directory
setwd(my_working_dir)

# MAIN LOOP OVER DAY of YEAR
#for (my_day in vec_days ) {
#...
#}
my_day <- paste0(my_year,"0131")
print(paste0("+++Working on ",my_day,"+++"))
for ( i in 1:length(g2_files) ) {
  system(paste0("cdo seldate,",my_day," ",original_g2_files[i]," ",g2_files[i]))
}

        
# Convert grib2 files to NetCDF format
print(paste0("+++Convert grib2 files to NetCDF format+++"))
for ( i in 1:length(g2_files) ) {
  system(paste0("cdo -f nc copy ", g2_files[i]," ", nc_files[i]))
  system(paste0("rm -f ", my_working_dir,"/",g2_files[i]))
}
print(paste0("+++OK+++"))

# Create wind speed from u&v scalars
print(paste0("+++Create wind speed from u&v scalars+++"))
system(paste0("cdo expr,'ws=sqrt(10u*10u+10v*10v)' ", nc_files[5], " ", paste0("test_",my_year,"_ws10m.nc")))
system(paste0("rm -f ", my_working_dir,"/",nc_files[5]))
nc_files[5]=paste0("test_",my_year,"_ws10m.nc")
print(paste0("+++OK+++"))

# Read files data and features
print(paste0("+++Reading input raster files+++"))
brickdataras=sapply(nc_files,brick)
brickdata=lapply(brickdataras,as.data.frame)
files_data=nc_open(nc_files[1])
basetime=as.POSIXct(gsub("minutes since ","",ncatt_get(files_data,"time")$units))
time=ncvar_get(files_data,"time")
newtime=basetime+time*60
nc_close(files_data);
print(paste0("+++OK+++"))

# Function: calculate mean radiant temperature from XXX, YYY, etc...
mrt_globe=function(t, tg, wind, diam_globe= 0.05, emis_globe= 0.97) 
{    stefanb = 0.0000000567;
     (((tg + 273.15)**4) + ((1.1 * (10**8) * (wind**0.6)) /(emis_globe*(diam_globe**0.4))) * (tg - t)**0.25) - 273.15;
}
 
# Function: reduce wind speed from 10-metre to 2-metre above surface 
# Roughly....reducewind=0.75*input
reducewind=function(x,ref=10,fin=2) (x*1/(log(ref/0.01)/log(fin/0.01)))

# utci vento a 10 m 
# wbgt vento a 2  m (see Brode et al, 2011)

########################################################################################################################################################
# hic sunt leones... 

outfilename=paste0("test_",format(newtime[1],"%Y%m%d"),".rds")
outfilenametime=paste0("times_test_",format(newtime[1],"%Y%m%d"),".rds")
res=list()

print(paste0("+++Working on sub-daily data+++"))
numCores<-detectCores()
registerDoParallel(numCores)

foreach (j=1:24, .combine=rbind) %dopar% {
#for ( j in 1:2) {
  print(paste0("++++++Working on time step",j,"+++"))

datafile=data.frame(lonvect=coordinates(brickdataras[[1]])[,1],
                    latvect=coordinates(brickdataras[[1]])[,2],
                    T_fhg=as.numeric(brickdata[[1]][,j])-273.15,
                    TD_fhg=as.numeric(brickdata[[2]][,j])-273.15)


datafile$RH_fhg=as.numeric(brickdata[[4]][,j])
datafile$dswrf_sfc=as.numeric(brickdata[[3]][,j])
datafile$wind_speed=as.numeric(brickdata[[5]][,j])
datafile$wind_speed2m=as.numeric(reducewind(brickdata[[5]][,j]))
datafile$albedo=rep(0.3,length(datafile$lonvect))
datafile$tglob=rep(NA,length(datafile$lonvect))
datafile$wbgtsha_var=rep(NA,length(datafile$lonvect))
datafile$wbgtsun_var=rep(NA,length(datafile$lonvect))

idbolam=1:nrow(datafile)
posixtimefile=as.POSIXct(newtime[j],format="%Y%m%d",tz="GMT")
datafile$timeh=rep(posixtimefile,length(datafile$lonvect))
wbgt.outdoors.ls<- lapply(idbolam,function(x) wbgt.Liljegren(datafile$T_fhg[x], 
                                                             datafile$TD_fhg[x],  
                                                             datafile$wind_speed2m[x], 
                                                             datafile$dswrf_sfc[x],  
                                                             dates=datafile$timeh[x],
                                                             lon=datafile$lonvect[x],
                                                             lat=datafile$latvect[x]))

tglob=unlist(lapply(wbgt.outdoors.ls,function(x) x$Tg))
datafile$wbgtsha_var=unlist(lapply(wbgt.outdoors.ls,function(x) x$Tnwb))*0.7+0.3*datafile$T_fhg
datafile$wbgtsun_var=unlist(lapply(wbgt.outdoors.ls,function(x) x$data))
datafile$MRT=mrt_globe(datafile$T_fhg,tglob,datafile$wind_speed2m,150)    # 150 mm
datafile$utci_var=future_sapply(seq_along(datafile$T_fhg),function(x) {biometeoR::utci(datafile$T_fhg[x],datafile$RH_fhg[x],datafile$wind_speed[x],datafile$MRT[x])})
datafile$utci_shade_var=future_sapply(seq_along(datafile$T_fhg),function(x) {biometeoR::utci(datafile$T_fhg[x],datafile$RH_fhg[x],datafile$wind_speed[x],datafile$T_fhg[x])})
res[[i]]=datafile[,c("wbgtsun_var","wbgtsha_var","utci_var","utci_shade_var","wind_speed","dswrf_sfc","tglob","RH_fhg","T_fhg")]
}
print(paste0("+++OK+++"))

# OK fino a qui
print(paste0("+++Saving biometeo data+++"))
saveRDS(res,file=outfilename)
saveRDS(newtime,outfilenametime)
#file.remove(dir(pattern(".nc")))
print(paste0("+++OK+++"))

stopImplicitCluster()
quit()

