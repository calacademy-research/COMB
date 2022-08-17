#read points summarize Veg plot and SALO data for 11.3 m^2 comparison

# 0. Getting started
# 
# [] build new script to fix shape file for ALL points
# install.packages("rgdal")
# library(rgdal)
library(sf)
library(here)
library(tidyverse)

# read in wildlife points
#[ ] update wild_points to include plant plots, see e-mail 2022-06-28 Sarah Jacobes
#Import
wild_points <- sf::read_sf(here("spatial", "input", "shapefiles", "WildlifePoints.shp")) # crs not included


names(wild_points) <- c(
  "point_d", "Cpls_Wt", "VEG_CSE", "AVIAN_S", "WHR_TSD", "SZ_DNS2",
  "Treatment", "geometry"
)

# get the right coordinate reference system & transform
st_crs(wild_points) <- 26910 # were collected using NAD83 coordinate
wild_points <- st_transform(wild_points, crs(study_area)) # transformed crs

# are points in the fire boundary area?
wild_points$inside_fire_boundary <- as.vector(st_intersects(fire_boundary, wild_points, sparse = FALSE))

#fix with "plotID_UTM.csv" that is now in the input directory (so running "read_points_output_data.R" creates it)
plotID_UTM <- read_csv(here("spatial", "input", "shapefiles", "plotID_UTM.csv"))
  
# plotID_UTM$plotID_av are the bird_points
# plotID_UTM$plotID_veg are the veg_points

#building a bigger table
new_wild_points <- left_join(plotID_UTM, wild_points, by = c("plotID_av" = "point_d"), keep = TRUE)

#replace the 87th point (avian # 1072) that has an empty geometry see: 
# https://gis.stackexchange.com/questions/244756/edit-sf-point-conditionally 
#scroll all the way down

#fix the 87th point and set the CRS
new_point <- st_point(c(743246.3, 4287041.4)) %>% 
  st_sfc(crs = 32610) 

#conditionally replace it
new_wild_points <- new_wild_points %>% mutate(geometry = st_sfc(ifelse(new_wild_points$plotID_av==1072, st_geometry(new_point), geometry)))

#could do for all vegetation points (hold off until edited [ ])

#write out fixed shapefile
sf::st_write(new_wild_points, dsn = here("spatial", "input", "shapefiles", "WildlifePointsFixed.shp"), append = FALSE)

# 
# 1. Vegetation plots
# 
# Summarize key variables (Sarah)
# 
# 2. Satellite data
# 
# [] still need to add in 11.3 m^2 (extract and summarize) for some variables ...
# 
# 3. Comparisons