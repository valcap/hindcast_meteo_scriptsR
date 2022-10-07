library(dplyr)
library(ggplot2)
library(sf)
library(gstat)
library(stars)

## a) Defining the projection of the maps
proj_lambert <- "+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

## b) Map Canada selected provinces.
canada.map <- st_read("C:/Users/rotundjo/Desktop/lpr_000b16a_e/lpr_000b16a_e.shp") %>% 
  filter(PRNAME == "Manitoba" |
           PRNAME == "Alberta" |
           PRNAME == "Saskatchewan") %>%
  st_transform(proj_lambert )

## c) Data set
data.yield <- read_csv("C:/Users/rotundjo/Desktop/yield.data.csv")
data.yield.coords <- data.yield %>% select(lat, long)
data.yield.data <- data.yield %>% select(av.yield)

#### Creating LAT LONG SpatialPoints
coordinates(data.yield.coords) = c("long", "lat")
proj4string(data.yield.coords) <- CRS("+init=epsg:4326") # 4326 is the EPSG identifier of WGS84.

### Converting the LAT LONG to the lambert CRS shown above
lambert_data.yield.coords = spTransform(data.yield.coords, proj_lambert)
lambert_data.yield.coords = as.data.frame(lambert_data.yield.coords)

#### Adding the Type column back to the data frame, with the new polar coordinates
lambert_data = cbind(lambert_data.yield.coords, data.yield.data)


coordinates(data.yield) = c("long", "lat")
proj4string(data.yield) <- CRS("+init=epsg:4326")
data.yield = spTransform(data.yield, proj_lambert)
data.yield.sf <- st_as_sf(data.yield, coords = c("long", "lat"), crs = proj_lambert)


## d) Variogram fit
v = variogram(av.yield~1, data.yield.sf)
plot(v, plot.numbers = TRUE)

v.m = fit.variogram(v, vgm("Sph", psill = 3, range = 3000000, nugget = 1))
plot(v, v.m, plot.numbers = FALSE)

## e) Building a empty grid
st_bbox(canada.map) %>%
  st_as_stars(dx = 1000) %>%
  st_set_crs(proj_lambert) %>%
  st_crop(canada.map) -> grd

## f) Kriging interpolation
k = krige(av.yield~1, data.yield.sf, grd, v.m)
