library(openair)
library(worldmet)
# Original data
# https://www.ncei.noaa.gov/data/global-summary-of-the-day/access/2019/

# Load observed data
firenze_met <- importNOAA(code = "161700-99999", year = 2019)
# Plot windorose
windRose(dublin_met)

# Dati MOLOCH (fare loop su ogni anno)
# 1. da wind10m a 10u && 10v
grib_copy -w shortName=10u /OCEANASTORE/progetti/spitbran/ERA5/winds/MOLOCH/2019/moloch_2019_wind10m.grib2 /OCEANASTORE/progetti/spitbran/work/moloch/moloch_2019_10u.grib2
grib_copy -w shortName=10v /OCEANASTORE/progetti/spitbran/ERA5/winds/MOLOCH/2019/moloch_2019_wind10m.grib2 /OCEANASTORE/progetti/spitbran/work/moloch/moloch_2019_10v.grib2
# 2. convertire da grib2 a nc (facoltativo???)
cdo -f nc copy /OCEANASTORE/progetti/spitbran/work/moloch/moloch_2019_10u.grib2 /OCEANASTORE/progetti/spitbran/work/moloch/moloch_2019_10u.nc
cdo -f nc copy /OCEANASTORE/progetti/spitbran/work/moloch/moloch_2019_10v.grib2 /OCEANASTORE/progetti/spitbran/work/moloch/moloch_2019_10v.nc
rm -f /OCEANASTORE/progetti/spitbran/work/moloch/moloch_2019_10u.grib2 /OCEANASTORE/progetti/spitbran/work/moloch/moloch_2019_10v.grib2
# 3. importare i file nc così costruiti
file_moloch_10u <- paste0("/OCEANASTORE/progetti/spitbran/work/moloch/moloch_2019_10u.nc")
file_moloch_10v <- paste0("/OCEANASTORE/progetti/spitbran/work/moloch/moloch_2019_10v.nc")
