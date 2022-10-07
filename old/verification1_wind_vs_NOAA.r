# Import libraries
library(openair)
library(worldmet)
library(raster)
library(ncdf4)
library(rWind)
library(verification)
library(plotrix)
library(lubridate)
library(timetk)
library(tidyverse)
library(padr)
library(numform)
library(RCurl)

# Global variables
my_working_dir  <- paste0("/OCEANASTORE/progetti/spitbran")
setwd(my_working_dir)
my_figs_dir     <- paste0("/OCEANASTORE/progetti/spitbran/work/figs")
my_output_dir   <- paste0("/OCEANASTORE/progetti/spitbran/work/res")
my_moloch_dir   <- paste0("/OCEANASTORE/progetti/spitbran/work/moloch")
my_NOAA_dir     <- paste0("/OCEANASTORE/progetti/spitbran/NOAA/global-hourly")
my_shp_dir      <- paste0("/OCEANASTORE/progetti/spitbran/work/shp")
my_anag_file    <- paste0("NOAA/stations_synop_BIGGER-DOMAIN_LOCALE.csv")
my_anag_data    <- read.csv(my_anag_file, sep=",", header=TRUE)
vector_years    <- seq(2018,2019)
append_results  <- TRUE

# Set-up output csv file and write the header
file_scores_csv <- paste0(my_output_dir,"/verification_scores_wind.csv")
df              <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(df)    <- c("YEAR","STATION_CODE","STATION_NAME","STATION_LAT","STATION_LON","STATION_ELE","INI_DATE","END_DATE","CORR","RMSE","BIAS","ME","MAE")
write.table(df, file=file_scores_csv, sep=";", append=append_results, quote = FALSE)

# Loop over year
for (my_year in vector_years) {
###
  # 1. Preliminary stuff about leap years
  print(paste0("___START year ",my_year,"___"))
  if (leap_year(my_year)==TRUE) {
    print(paste0(my_year," is a leap year"))
    nts<-8760
  } else {
    nts<-8736
  }
  my_ini_date       <- as.POSIXct(paste0(my_year,'-01-02 00:00:00'), format="%Y-%m-%d %H:%M:%S", tz="UTC")
  my_end_date       <- as.POSIXct(paste0(my_year,'-12-31 23:00:00'), format="%Y-%m-%d %H:%M:%S", tz="UTC")
###
  # 2. Check if files exists and load model data
  file_moloch_10u   <- paste0(my_moloch_dir,"/moloch_",my_year,"_10u.nc")
  file_moloch_10v   <- paste0(my_moloch_dir,"/moloch_",my_year,"_10v.nc")
  if ((file.exists(file_moloch_10u)==FALSE) || (file.exists(file_moloch_10v)==FALSE)) {
    print("ops missing MOLOCH 10u or 10v file (or both)")
    quit()
  }
  # Start the clock!
  stm               <- proc.time()
  print(paste0("___Loading data from MOLOCH 10u file___"))
  raster_moloch_10u <- brick(file_moloch_10u)
  print(paste0("___Loading data from MOLOCH 10v file___"))
  raster_moloch_10v <- brick(file_moloch_10v)
  # Stop the clock!
  etm               <- proc.time() - stm
  print(paste0("___OK in ",as.character(etm[3])," seconds"))
###
  # 3. Loop over codes
  num_valid_data     <- 0
  for (my_obse_code in c(as.vector(my_anag_data$STATION_ID))) {
    my_obse_code     <- f_pad_zero(my_obse_code,width=11)
    # some codes, in some years, have problems...
    if ((my_obse_code == '16266099999') || (my_obse_code == '11185099999')) {next}
    if (my_year == '1983') {
      if ((my_obse_code == '11185099999')) {next}
    }
    if (my_year == '1984') {
      if ((my_obse_code == '06632099999')) {next}
    }
    if (my_year == '1985') {
      if ((my_obse_code == '16234099999') || (my_obse_code == '11204099999') || (my_obse_code == '16542099999')) {
        next
      }
    }
    if (my_year == '1990') {
      if ((my_obse_code == '11220099999')) {next}
    }
    if (my_year == '1992') {
      if ((my_obse_code == '07768099999')) {next}
    }
    if (my_year == '1994') {
      if ((my_obse_code == '07752099999') || (my_obse_code == '13348099999') || (my_obse_code == '07677099999') || (my_obse_code == '11229099999')) {
        next
      }
    }
    if (my_year == '1995') {
      if ((my_obse_code == '16434099999') || (my_obse_code == '16232099999')  || (my_obse_code == '07688099999')) {
        next
      }
    }
    if (my_year == '1996') {
      if ((my_obse_code == '16294099999') || (my_obse_code == '16214099999')  || (my_obse_code == '16310099999')) {
        next
      }
    }
    if (my_year == '1997') {
      if ((my_obse_code == '07593099999')  || (my_obse_code == '16450099999')) {
        next
      }
    }
    if (my_year == '1998') {
      if ((my_obse_code == '07768099999')  || (my_obse_code == '60714099999')) {
        next
      }
    }
    if (my_year == '1999') {
      if ((my_obse_code == '07593099999')   || (my_obse_code == '16294099999')  || (my_obse_code == '16541099999') || (my_obse_code == '16337099999') || (my_obse_code == '16542099999') || (my_obse_code == '16224099999') || (my_obse_code == '16172099999')) {
        next
      }
    }
    if (my_year == '2000') {
      if ((my_obse_code == '07768099999')) {next}
    }
    if (my_year == '2001') {
      if ((my_obse_code == '16052099999')) {next}
    }
    if (my_year == '2002') {
      if ((my_obse_code == '07688099999')) {next}
    }
    if (my_year == '2003') {
      if ((my_obse_code == '11210099999') || (my_obse_code == '11140099999')  || (my_obse_code == '07688099999')) {
        next
      }
    }
    if (my_year == '2004') {
      if ((my_obse_code == '16260099999') || (my_obse_code == '11140099999')  || (my_obse_code == '07688099999') || (my_obse_code == '11220099999')) {# || (my_obse_code == '16542099999') || (my_obse_code == '16224099999') || (my_obse_code == '16172099999')) {
        next
      }
    }
    if (my_year == '2005') {
      if ((my_obse_code == '16364099999') || (my_obse_code == '11210099999')  || (my_obse_code == '07688099999') || (my_obse_code == '16095099999')) {# || (my_obse_code == '16542099999') || (my_obse_code == '16224099999') || (my_obse_code == '16172099999')) {
        next
      }
    }
    if (my_year == '2006') {
      if ((my_obse_code == '16364099999') || (my_obse_code == '11214099999')  || (my_obse_code == '16450099999') || (my_obse_code == '11220099999')) {# || (my_obse_code == '16542099999') || (my_obse_code == '16224099999') || (my_obse_code == '16172099999')) {
        next
      }
    }
    if (my_year == '2007') {
      if ((my_obse_code == '16364099999') || (my_obse_code == '16260099999')  || (my_obse_code == '07688099999') || (my_obse_code == '11220099999')) {# || (my_obse_code == '16542099999') || (my_obse_code == '16224099999') || (my_obse_code == '16172099999')) {
        next
      }
    }
    if (my_year == '2008') {
      if ((my_obse_code == '11119099999') || (my_obse_code == '16129099999')  || (my_obse_code == '16364099999') || (my_obse_code == '11214099999') || (my_obse_code == '11140099999') || (my_obse_code == '11220099999')) {# || (my_obse_code == '16172099999')) {
        next
      }
    }
    if (my_year == '2009') {
      if ((my_obse_code == '11119099999') || (my_obse_code == '16416099999')  || (my_obse_code == '16334099999') || (my_obse_code == '16129099999') || (my_obse_code == '16364099999') || (my_obse_code == '11214099999') || (my_obse_code == '11140099999') || (my_obse_code == '07688099999') || (my_obse_code == '11220099999') || (my_obse_code == '11241099999')) {
        next
      }
    }
    if (my_year == '2010') {
      if ((my_obse_code == '11119099999') || (my_obse_code == '11175099999')  || (my_obse_code == '07688099999')) {# || (my_obse_code == '16129099999') || (my_obse_code == '16364099999') || (my_obse_code == '11214099999')) {# || (my_obse_code == '16172099999')) {
        next
      }
    }
#    if (my_year == '2011') {    }
    if (my_year == '2012') {
      if ((my_obse_code == '16416099999') || (my_obse_code == '16364099999')  || (my_obse_code == '11161099999') || (my_obse_code == '16469099999')) {# || (my_obse_code == '16364099999') || (my_obse_code == '11214099999')) {# || (my_obse_code == '16172099999')) {
        next
      }
    }
    if (my_year == '2018') {
      if ((my_obse_code == '16121099999')) {# || (my_obse_code == '16364099999')  || (my_obse_code == '11161099999') || (my_obse_code == '16469099999')) {# || (my_obse_code == '16364099999') || (my_obse_code == '11214099999')) {# || (my_obse_code == '16172099999')) {
        next
      }
    }

#####
    # 3.1 (Down- and )Load observed data from
    # https://www.ncei.noaa.gov/data/global-hourly/access/2019/01007099999.csv
    print(paste0("+++START code ",my_obse_code,""))
    local_rds_file   <- paste0(my_NOAA_dir,"/",my_year,"/NOAA_global-hourly_",my_obse_code,".rds")
    if (file.exists(local_rds_file)) {
        my_obse_data       <- readRDS(local_rds_file)
        print(paste0("++++++OK (local file)"))
    } else {
      remote_csv_file      <- paste0("https://www.ncei.noaa.gov/data/global-hourly/access/",my_year,"/",my_obse_code,".csv")
      file_exist           <- url.exists(remote_csv_file)
      if (file_exist == FALSE) {
        print(paste0("+++Remote file is missing"))
        next
      } else {
        print(paste0("+++Remote file exists"))
        num_valid_data     <- num_valid_data + 1
        print(paste0("++++++Downloading data"))
        stm                <- proc.time()
        my_obse_data       <- importNOAA(code = as.character(my_obse_code), year = my_year, quiet = TRUE, n.cores = 1)
        saveRDS(my_obse_data, file=local_rds_file)
        etm                <- proc.time() - stm
        print(paste0("++++++OK in ",as.character(etm[3])," seconds"))
      }
    }
#####
    # 3.2 Get time resolution of the imported data
    print(paste0("++++++Filling missing data with NA"))
    my_obse_data_timeres          <- my_obse_data$date %>% get_interval()
    if (my_obse_data_timeres == "hour") {
      # fill the gap (hour)
      my_obse_data_hour_filled    <- my_obse_data %>% pad('hour', by='date', start_val=my_ini_date, end_val=my_end_date)
      # subset data and remove duplicates (if any)
      my_obse_data_hour_filled_OK <- filter(my_obse_data_hour_filled, date >= my_ini_date & date <= my_end_date) %>% distinct(date, .keep_all = TRUE)
    } else {
      # if time resolution is not "hour" => create a new column with only hour time step named "date_hour"
      my_obse_data_hour           <- my_obse_data %>% thicken('hour')
      # fill the gap (hour)
      my_obse_data_hour_filled    <- my_obse_data_hour %>% pad('hour', by='date_hour', start_val=my_ini_date, end_val=my_end_date)
      # subset data and remove duplicates (if any)
      my_obse_data_hour_filled_OK <- filter(my_obse_data_hour_filled, date_hour >= my_ini_date & date_hour <= my_end_date) %>% distinct(date_hour, .keep_all = TRUE)
    }
    print(paste0("++++++OK"))
#####
    # 3.3 Re-arrange observed data
    print(paste0("++++++Re-arranging data"))
    obse_ws         <- my_obse_data_hour_filled_OK$ws
    obse_wd         <- my_obse_data_hour_filled_OK$wd
    prog <- 1; obse_code <- NA
    while (is.na(obse_code)) {
      obse_code     <- my_obse_data_hour_filled_OK$code[prog]
      prog <- prog  + 1 
    }
    prog <- 1; obse_stat <- NA
    while (is.na(obse_stat)) {
      obse_stat     <- my_obse_data_hour_filled_OK$station[prog]
      prog <- prog  + 1 
    }
    prog <- 1; obse_lon <- NA
    while (is.na(obse_lon)) {
      obse_lon      <- my_obse_data_hour_filled_OK$longitude[prog]
      prog <- prog  + 1 
    }
    prog <- 1; obse_lat <- NA
    while (is.na(obse_lat)) {
      obse_lat      <- my_obse_data_hour_filled_OK$latitude[prog]
      prog <- prog  + 1 
    }
    prog <- 1; obse_ele <- NA
    while (is.na(obse_ele)) {
      obse_ele      <- my_obse_data_hour_filled_OK$elev[prog]
      prog <- prog  + 1 
    }
    print(paste0("++++++OK"))
#####
    # 4. Extract wind model data @ (obse_lon,obse_lat)
    print(paste0("++++++Extracting data from MOLOCH files"))
    stm             <- proc.time()
    molo_10u        <- as.data.frame(raster::extract(raster_moloch_10u, SpatialPoints(cbind(obse_lon,obse_lat))))
    molo_10v        <- as.data.frame(raster::extract(raster_moloch_10v, SpatialPoints(cbind(obse_lon,obse_lat))))
    etm             <- proc.time() - stm
    print(paste0("++++++OK in ",as.character(etm[3])," seconds"))
#####
    # 5. Model data: from zonal and meridional wind to wind speed and direction using the rWind package
    data_from_rWind <- uv2ds(t(molo_10u),t(molo_10v))
    molo_wd         <- data_from_rWind[,1]
    molo_wd         <- molo_wd[1:nts]
    molo_ws         <- data_from_rWind[,2]
    molo_ws         <- molo_ws[1:nts]
    my_molo_data    <- as.data.frame(cbind(molo_ws,molo_wd))
#####
    # 6. Calculate standard verification scores
    #    such as r, RMSE, (multiplicative) bias, ME, MAE
    print(paste0("++++++Calculating standard verification scores for wind speed"))
    cor_molo        <- sprintf("%.3f",cor(molo_ws,obse_ws, use = "complete.obs"))
    rmse_molo       <- sprintf("%.3f",sqrt( mean( (molo_ws-obse_ws)^2 , na.rm = TRUE ) ))
    mbias_molo      <- sprintf("%.3f",(mean(molo_ws, na.rm = TRUE))/(mean(obse_ws, na.rm = TRUE)))
    verify_molo     <- verify(obse_ws,molo_ws, frcst.type = "cont", obs.type = "cont",)
    me_molo         <- sprintf("%.3f",verify_molo$ME)
    mae_molo        <- sprintf("%.3f",verify_molo$MAE)
    print(paste0("++++++OK"))
#####
    # 7. Write the standard verification scores to the csv output file
    print(paste0("++++++Saving data to ",file_scores_csv,""))
    df <- data.frame(my_year,obse_code,obse_stat,obse_lat,obse_lon,obse_ele,my_ini_date,my_end_date,cor_molo,rmse_molo,mbias_molo,me_molo,mae_molo)
    write.table(df, file=file_scores_csv, col.names=F, row.names=FALSE, sep=";", append=TRUE, quote=FALSE)
    print(paste0("++++++OK"))
    print(paste0("+++END code ",my_obse_code,""))
  # End of loop over codes
  }
  print(paste0("___END year ",my_year,"___"))

# End of loop over years
}

#
quit()

