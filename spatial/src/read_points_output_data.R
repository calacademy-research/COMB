#### Sample raster data around the points ####
# read_points_output_data.R
# Created by: Durrell Kapan
# Created on: 5 April 2022
# modified from Caples_Spatial_EDA.md

#### Load Libraries ####

#load libraries
library(tidyverse)
library(broom)
library(circular)
library(clipr)
library(elevatr)
library(exactextractr)
library(forcats)
library(GGally)
library(ggpubr)
library(googlesheets4)
library(leaflet) # for interactive mapping (TBD)
library(plotKML)
library(purrr)
library(stringr)
library(cowplot)
library(RStoolbox)
library(sf)
library(raster)
library(rgdal)

library(here)
library(googledrive) #connect to google drive
library(data.table)

#load functions
source(here("comb_functions.R"))

#------------------------------------------------------------
# Download/organize input files from google drive if not already downloaded
#------------------------------------------------------------
# EITHER
# run read_adjust_combine_rasters.R
# source(here("/spatial/src/read_adjust_combine_rasters.R"))
# -OR-
# **Import raster brick**

# Creating the folder within inputs that contains the symlinks 
# that point to the output directory 
if (dir.exists(here("spatial/output/rasters/")) == F) {
  dir.create(here("spatial/output/rasters/"))
}
#run to ensure raster brick is loaded
drive_sync(here("spatial","output","rasters"), drive_folder = drive_ls("https://drive.google.com/drive/folders/1sbgR_OMtK-Hq6P6lVBFIK8xbcjmJDQVV")$id[1])
#read raster

if (exists(x = "canopy_fuel_nbr_dem_RAVG_LIDAR") == F) {
canopy_fuel_nbr_dem_RAVG_LIDAR <- stack(here("spatial","output","rasters","canopy_fuel_nbr_dem_RAVG_LIDAR.grd"))
}

#------------------------------------------------------------
# **Import shape files**
#   
# Study area file
# -   study_area
# -   Fire boundary file
# -   fire_boundary
# -   Sampling points
# -   wild_points
# everything to be placed into WGS 84 / UTM zone 10N = crs = 32610
#------------------------------------------------------------

#get the data
shapedir <- here("spatial","input","shapefiles")
google_file_names <- "https://drive.google.com/drive/folders/15uXx-I1U9cDTzoI0ptNv_wpuDYju_IRv"
drive_sync(local_dir = shapedir, drive_folder = google_file_names)

#read in study area
study_area <- sf::read_sf(here("spatial","input","shapefiles","study_boundary.shp"))

#read in fire boundary
fire_boundary <- sf::read_sf(here("spatial","input","shapefiles","ca3872412014620191010_20181118_20191118_burn_bndy.shp"))
fire_boundary<-st_transform(fire_boundary, crs(study_area)) #transformed crs

#read in wildlife points 
wild_points <- sf::read_sf(here("spatial","input","shapefiles","WildlifePoints.shp")) #crs not included
names(wild_points)<-c("point_d", "Cpls_Wt", "VEG_CSE", "AVIAN_S", "WHR_TSD", "SZ_DNS2", 
                      "Treatment", "geometry")

#get the right coordinate reference system & transform
st_crs(wild_points) = 26910 #were collected using NAD83 coordinate 
wild_points<-st_transform(wild_points, crs(study_area)) #transformed crs


buffer<-400
inner_boundary<-as.vector(extent(canopy_imgStack)+c(buffer,-buffer,buffer,-buffer))
inner_boundary<-rbind(inner_boundary[c(1,3)], 
                      inner_boundary[c(2,4)],
                      inner_boundary[c(1,4)],
                      inner_boundary[c(2,3)])

inner_boundary<-rbind(inner_boundary,c(744000,4289000))

## 50 & 100 meter buffer around sampling points = 1 & 4ha plots: 
##(100m in all four directions from point center, 200m ditto)
wldf_50 <- wild_points  %>% sf::st_buffer(50,endCapStyle = "SQUARE") 
wldf_100 <- wild_points  %>% sf::st_buffer(100,endCapStyle = "SQUARE") 

#------------------------------------------------------------
# -   select these for each radius
# 
# -   (50, 100m radii squares = 1ha, 4ha)
# 
# -   extract all data (raw) for and summarize for these radii
# 
# -   delint variable names ...
#------------------------------------------------------------

#extract variables for 1 & 4ha plots 
#for now missing standard error (standard deviation divided by the square root of the sample size) "1" not subtracted...
#[ ]eventually build custom extractor functions e.g. stderr <- function(x) sd(x)/sqrt(length(x)) #or depending on your flavor -1.

#[ ] fix so we can make it depend upon
# canopy_fuel_nbr_dem_RAVG_LIDAR$[1:21]
# e.g. canopy_fuel_nbr_dem_RAVG_LIDAR$[1:21] #==canopy_imgStack

extract_canopy_1ha <- exactextractr::exact_extract(canopy_imgStack, wldf_50, c("mean","median","min","max","count"))
extract_canopy_4ha <- exactextractr::exact_extract(canopy_imgStack, wldf_100, c("mean","median","min","max","count"))

dput(colnames(extract_canopy_1ha)) #adjust names [ ] adopt standard naming NEXT

#variable pattern to hand-code for now

colnames(extract_canopy_1ha)<-
  c("mean.CanopyBaseHeight_2018_1ha", "mean.CanopyBaseHeight_2019_1ha", "mean.CanopyBaseHeight_2020_1ha", 
    "mean.CanopyBulkDensity_2018_1ha", "mean.CanopyBulkDensity_2019_1ha", "mean.CanopyBulkDensity_2020_1ha", 
    "mean.CanopyLayerCount_2018_1ha", "mean.CanopyLayerCount_2019_1ha", "mean.CanopyLayerCount_2020_1ha", 
    "mean.CanopyCover_2018_1ha", "mean.CanopyCover_2019_1ha", "mean.CanopyCover_2020_1ha", 
    "mean.CanopyHeight_2018_1ha", "mean.CanopyHeight_2019_1ha", 
    "mean.CanopyHeight_2020_1ha", "median.CanopyBaseHeight_2018_1ha", 
    "median.CanopyBaseHeight_2019_1ha", "median.CanopyBaseHeight_2020_1ha", 
    "median.CanopyBulkDensity_2018_1ha", "median.CanopyBulkDensity_2019_1ha", 
    "median.CanopyBulkDensity_2020_1ha", "median.CanopyLayerCount_2018_1ha", 
    "median.CanopyLayerCount_2019_1ha", "median.CanopyLayerCount_2020_1ha", 
    "median.CanopyCover_2018_1ha", "median.CanopyCover_2019_1ha", 
    "median.CanopyCover_2020_1ha", "median.CanopyHeight_2018_1ha", 
    "median.CanopyHeight_2019_1ha", "median.CanopyHeight_2020_1ha", 
    "min.CanopyBaseHeight_2018_1ha", "min.CanopyBaseHeight_2019_1ha", "min.CanopyBaseHeight_2020_1ha", 
    "min.CanopyBulkDensity_2018_1ha", "min.CanopyBulkDensity_2019_1ha", "min.CanopyBulkDensity_2020_1ha", 
    "min.CanopyLayerCount_2018_1ha", "min.CanopyLayerCount_2019_1ha", "min.CanopyLayerCount_2020_1ha", 
    "min.CanopyCover_2018_1ha", "min.CanopyCover_2019_1ha", "min.CanopyCover_2020_1ha", 
    "min.CanopyHeight_2018_1ha", "min.CanopyHeight_2019_1ha", "min.CanopyHeight_2020_1ha", 
    "max.CanopyBaseHeight_2018_1ha", "max.CanopyBaseHeight_2019_1ha", "max.CanopyBaseHeight_2020_1ha", 
    "max.CanopyBulkDensity_2018_1ha", "max.CanopyBulkDensity_2019_1ha", "max.CanopyBulkDensity_2020_1ha", 
    "max.CanopyLayerCount_2018_1ha", "max.CanopyLayerCount_2019_1ha", "max.CanopyLayerCount_2020_1ha", 
    "max.CanopyCover_2018_1ha", "max.CanopyCover_2019_1ha", "max.CanopyCover_2020_1ha", 
    "max.CanopyHeight_2018_1ha", "max.CanopyHeight_2019_1ha", "max.CanopyHeight_2020_1ha", 
    "count.CanopyBaseHeight_2018_1ha", "count.CanopyBaseHeight_2019_1ha", "count.CanopyBaseHeight_2020_1ha", 
    "count.CanopyBulkDensity_2018_1ha", "count.CanopyBulkDensity_2019_1ha", 
    "count.CanopyBulkDensity_2020_1ha", "count.CanopyLayerCount_2018_1ha", 
    "count.CanopyLayerCount_2019_1ha", "count.CanopyLayerCount_2020_1ha", "count.CanopyCover_2018_1ha", 
    "count.CanopyCover_2019_1ha", "count.CanopyCover_2020_1ha", 
    "count.CanopyHeight_2018_1ha", "count.CanopyHeight_2019_1ha", 
    "count.CanopyHeight_2020_1ha")

colnames(extract_canopy_4ha)<-
  c("mean.CanopyBaseHeight_2018_4ha", "mean.CanopyBaseHeight_2019_4ha", "mean.CanopyBaseHeight_2020_4ha", 
    "mean.CanopyBulkDensity_2018_4ha", "mean.CanopyBulkDensity_2019_4ha", "mean.CanopyBulkDensity_2020_4ha", 
    "mean.CanopyLayerCount_2018_4ha", "mean.CanopyLayerCount_2019_4ha", "mean.CanopyLayerCount_2020_4ha", 
    "mean.CanopyCover_2018_4ha", "mean.CanopyCover_2019_4ha", "mean.CanopyCover_2020_4ha", 
    "mean.CanopyHeight_2018_4ha", "mean.CanopyHeight_2019_4ha", 
    "mean.CanopyHeight_2020_4ha", "median.CanopyBaseHeight_2018_4ha", 
    "median.CanopyBaseHeight_2019_4ha", "median.CanopyBaseHeight_2020_4ha", 
    "median.CanopyBulkDensity_2018_4ha", "median.CanopyBulkDensity_2019_4ha", 
    "median.CanopyBulkDensity_2020_4ha", "median.CanopyLayerCount_2018_4ha", 
    "median.CanopyLayerCount_2019_4ha", "median.CanopyLayerCount_2020_4ha", 
    "median.CanopyCover_2018_4ha", "median.CanopyCover_2019_4ha", 
    "median.CanopyCover_2020_4ha", "median.CanopyHeight_2018_4ha", 
    "median.CanopyHeight_2019_4ha", "median.CanopyHeight_2020_4ha", 
    "min.CanopyBaseHeight_2018_4ha", "min.CanopyBaseHeight_2019_4ha", "min.CanopyBaseHeight_2020_4ha", 
    "min.CanopyBulkDensity_2018_4ha", "min.CanopyBulkDensity_2019_4ha", "min.CanopyBulkDensity_2020_4ha", 
    "min.CanopyLayerCount_2018_4ha", "min.CanopyLayerCount_2019_4ha", "min.CanopyLayerCount_2020_4ha", 
    "min.CanopyCover_2018_4ha", "min.CanopyCover_2019_4ha", "min.CanopyCover_2020_4ha", 
    "min.CanopyHeight_2018_4ha", "min.CanopyHeight_2019_4ha", "min.CanopyHeight_2020_4ha", 
    "max.CanopyBaseHeight_2018_4ha", "max.CanopyBaseHeight_2019_4ha", "max.CanopyBaseHeight_2020_4ha", 
    "max.CanopyBulkDensity_2018_4ha", "max.CanopyBulkDensity_2019_4ha", "max.CanopyBulkDensity_2020_4ha", 
    "max.CanopyLayerCount_2018_4ha", "max.CanopyLayerCount_2019_4ha", "max.CanopyLayerCount_2020_4ha", 
    "max.CanopyCover_2018_4ha", "max.CanopyCover_2019_4ha", "max.CanopyCover_2020_4ha", 
    "max.CanopyHeight_2018_4ha", "max.CanopyHeight_2019_4ha", "max.CanopyHeight_2020_4ha", 
    "count.CanopyBaseHeight_2018_4ha", "count.CanopyBaseHeight_2019_4ha", "count.CanopyBaseHeight_2020_4ha", 
    "count.CanopyBulkDensity_2018_4ha", "count.CanopyBulkDensity_2019_4ha", 
    "count.CanopyBulkDensity_2020_4ha", "count.CanopyLayerCount_2018_4ha", 
    "count.CanopyLayerCount_2019_4ha", "count.CanopyLayerCount_2020_4ha", "count.CanopyCover_2018_4ha", 
    "count.CanopyCover_2019_4ha", "count.CanopyCover_2020_4ha", 
    "count.CanopyHeight_2018_4ha", "count.CanopyHeight_2019_4ha", 
    "count.CanopyHeight_2020_4ha")

canopy<-cbind(extract_canopy_1ha,extract_canopy_4ha)

extract_fuel_1ha <- exactextractr::exact_extract(fuel_imgStack, wldf_50, c("mean","median","min","max","count"))

extract_fuel_4ha <- exactextractr::exact_extract(fuel_imgStack, wldf_100, c("mean","median","min","max","count"))

dput(colnames(extract_fuel_1ha)) #adjust names [ ] adopt standard naming NEXT

colnames(extract_fuel_1ha)<-c("mean.LadderFuelDensity_2018_1ha", "mean.LadderFuelDensity_2019_1ha", 
                              "mean.LadderFuelDensity_2020_1ha", "mean.SurfaceFuels_2018_1ha", "mean.SurfaceFuels_2019_1ha", 
                              "mean.SurfaceFuels_2020_1ha", "median.LadderFuelDensity_2018_1ha", "median.LadderFuelDensity_2019_1ha", 
                              "median.LadderFuelDensity_2020_1ha", "median.SurfaceFuels_2018_1ha", "median.SurfaceFuels_2019_1ha", 
                              "median.SurfaceFuels_2020_1ha", "min.LadderFuelDensity_2018_1ha", "min.LadderFuelDensity_2019_1ha", 
                              "min.LadderFuelDensity_2020_1ha", "min.SurfaceFuels_2018_1ha", "min.SurfaceFuels_2019_1ha", 
                              "min.SurfaceFuels_2020_1ha", "max.LadderFuelDensity_2018_1ha", "max.LadderFuelDensity_2019_1ha", 
                              "max.LadderFuelDensity_2020_1ha", "max.SurfaceFuels_2018_1ha", "max.SurfaceFuels_2019_1ha", 
                              "max.SurfaceFuels_2020_1ha", "count.LadderFuelDensity_2018_1ha", "count.LadderFuelDensity_2019_1ha", 
                              "count.LadderFuelDensity_2020_1ha", "count.SurfaceFuels_2018_1ha", "count.SurfaceFuels_2019_1ha", 
                              "count.SurfaceFuels_2020_1ha")

colnames(extract_fuel_4ha)<-c("mean.LadderFuelDensity_2018_4ha", "mean.LadderFuelDensity_2019_4ha", 
                              "mean.LadderFuelDensity_2020_4ha", "mean.SurfaceFuels_2018_4ha", "mean.SurfaceFuels_2019_4ha", 
                              "mean.SurfaceFuels_2020_4ha", "median.LadderFuelDensity_2018_4ha", "median.LadderFuelDensity_2019_4ha", 
                              "median.LadderFuelDensity_2020_4ha", "median.SurfaceFuels_2018_4ha", "median.SurfaceFuels_2019_4ha", 
                              "median.SurfaceFuels_2020_4ha", "min.LadderFuelDensity_2018_4ha", "min.LadderFuelDensity_2019_4ha", 
                              "min.LadderFuelDensity_2020_4ha", "min.SurfaceFuels_2018_4ha", "min.SurfaceFuels_2019_4ha", 
                              "min.SurfaceFuels_2020_4ha", "max.LadderFuelDensity_2018_4ha", "max.LadderFuelDensity_2019_4ha", 
                              "max.LadderFuelDensity_2020_4ha", "max.SurfaceFuels_2018_4ha", "max.SurfaceFuels_2019_4ha", 
                              "max.SurfaceFuels_2020_4ha", "count.LadderFuelDensity_2018_4ha", "count.LadderFuelDensity_2019_4ha", 
                              "count.LadderFuelDensity_2020_4ha", "count.SurfaceFuels_2018_4ha", "count.SurfaceFuels_2019_4ha", 
                              "count.SurfaceFuels_2020_4ha")

fuels<-cbind(extract_fuel_1ha,extract_fuel_4ha)

#nbr
extract_nbr_1ha <- exactextractr::exact_extract(nbr_imgStack, wldf_50, c("mean","median","min","max","count"))
extract_nbr_4ha <- exactextractr::exact_extract(nbr_imgStack, wldf_100, c("mean","median","min","max","count"))

dput(colnames(extract_nbr_1ha)) #adjust names

colnames(extract_nbr_1ha)<-
  c("mean.Caldor_dNBR_Nov20_Oct21_1ha", "mean.Caldor_dNBR2_Nov20_Oct21_1ha", 
    "mean.Caldor_NBR_Oct21_1ha", "mean.Caples_dNBR_Nov18_Nov19_1ha", "mean.Caples_NBR_Nov18_1ha", 
    "mean.Caples_NBR_Nov19_1ha", "mean.Caples_NBR_Nov20_1ha", "mean.Caples_NBR_Oct21_NBR_1ha", 
    "median.Caldor_dNBR_Nov20_Oct21_1ha", "median.Caldor_dNBR2_Nov20_Oct21_1ha", 
    "median.Caldor_NBR_Oct21_1ha", "median.Caples_dNBR_Nov18_Nov19_1ha", 
    "median.Caples_NBR_Nov18_1ha", "median.Caples_NBR_Nov19_1ha", "median.Caples_NBR_Nov20_1ha", 
    "median.Caples_NBR_Oct21_NBR_1ha", "min.Caldor_dNBR_Nov20_Oct21_1ha", 
    "min.Caldor_dNBR2_Nov20_Oct21_1ha", "min.Caldor_NBR_Oct21_1ha", "min.Caples_dNBR_Nov18_Nov19_1ha", 
    "min.Caples_NBR_Nov18_1ha", "min.Caples_NBR_Nov19_1ha", "min.Caples_NBR_Nov20_1ha", 
    "min.Caples_NBR_Oct21_NBR_1ha", "max.Caldor_dNBR_Nov20_Oct21_1ha", "max.Caldor_dNBR2_Nov20_Oct21_1ha", 
    "max.Caldor_NBR_Oct21_1ha", "max.Caples_dNBR_Nov18_Nov19_1ha", "max.Caples_NBR_Nov18_1ha", 
    "max.Caples_NBR_Nov19_1ha", "max.Caples_NBR_Nov20_1ha", "max.Caples_NBR_Oct21_NBR_1ha", 
    "count.Caldor_dNBR_Nov20_Oct21_1ha", "count.Caldor_dNBR2_Nov20_Oct21_1ha", 
    "count.Caldor_NBR_Oct21_1ha", "count.Caples_dNBR_Nov18_Nov19_1ha", "count.Caples_NBR_Nov18_1ha", 
    "count.Caples_NBR_Nov19_1ha", "count.Caples_NBR_Nov20_1ha", "count.Caples_NBR_Oct21_NBR_1ha"
  )

colnames(extract_nbr_4ha)<-
  c("mean.Caldor_dNBR_Nov20_Oct21_4ha", "mean.Caldor_dNBR2_Nov20_Oct21_4ha", 
    "mean.Caldor_NBR_Oct21_4ha", "mean.Caples_dNBR_Nov18_Nov19_4ha", "mean.Caples_NBR_Nov18_4ha", 
    "mean.Caples_NBR_Nov19_4ha", "mean.Caples_NBR_Nov20_4ha", "mean.Caples_NBR_Oct21_NBR_4ha", 
    "median.Caldor_dNBR_Nov20_Oct21_4ha", "median.Caldor_dNBR2_Nov20_Oct21_4ha", 
    "median.Caldor_NBR_Oct21_4ha", "median.Caples_dNBR_Nov18_Nov19_4ha", 
    "median.Caples_NBR_Nov18_4ha", "median.Caples_NBR_Nov19_4ha", "median.Caples_NBR_Nov20_4ha", 
    "median.Caples_NBR_Oct21_NBR_4ha", "min.Caldor_dNBR_Nov20_Oct21_4ha", 
    "min.Caldor_dNBR2_Nov20_Oct21_4ha", "min.Caldor_NBR_Oct21_4ha", "min.Caples_dNBR_Nov18_Nov19_4ha", 
    "min.Caples_NBR_Nov18_4ha", "min.Caples_NBR_Nov19_4ha", "min.Caples_NBR_Nov20_4ha", 
    "min.Caples_NBR_Oct21_NBR_4ha", "max.Caldor_dNBR_Nov20_Oct21_4ha", "max.Caldor_dNBR2_Nov20_Oct21_4ha", 
    "max.Caldor_NBR_Oct21_4ha", "max.Caples_dNBR_Nov18_Nov19_4ha", "max.Caples_NBR_Nov18_4ha", 
    "max.Caples_NBR_Nov19_4ha", "max.Caples_NBR_Nov20_4ha", "max.Caples_NBR_Oct21_NBR_4ha", 
    "count.Caldor_dNBR_Nov20_Oct21_4ha", "count.Caldor_dNBR2_Nov20_Oct21_4ha", 
    "count.Caldor_NBR_Oct21_4ha", "count.Caples_dNBR_Nov18_Nov19_4ha", "count.Caples_NBR_Nov18_4ha", 
    "count.Caples_NBR_Nov19_4ha", "count.Caples_NBR_Nov20_4ha", "count.Caples_NBR_Oct21_NBR_4ha"
  )

nbr<-cbind(extract_nbr_1ha,extract_nbr_4ha)

#elevation too
extract_elevation_1ha<-exactextractr::exact_extract(dem_imgStack[[1]], wldf_50, c("mean","median","min","max","count"))
extract_elevation_4ha<-exactextractr::exact_extract(dem_imgStack[[1]], wldf_100, c("mean","median","min","max","count"))

colnames(extract_elevation_1ha)<-c("mean.Elevation_2020_1ha", "median.Elevation_2020_1ha", "min.Elevation_2020_1ha", "max.Elevation_NA_1ha", "count.Elevation_NA_1ha")
colnames(extract_elevation_4ha)<-c("mean.Elevation_2020_4ha", "median.Elevation_2020_4ha", "min.Elevation_2020_4ha", "max.Elevation_NA_4ha", "count.Elevation_NA_4ha")

elevation<-cbind(extract_elevation_1ha, extract_elevation_4ha)

#binary extraction easiest
#define binary aspect Raster
aspect.bi<-dem_imgStack$aspect
#extract values for Pat's binary categorization: aspect>315 | aspect<135
#aspect.bi <- getValues(asp_bi_imgStack) 
aspect.bi <- (aspect.bi > 315 | aspect.bi < 135)
#values(asp_bi_imgStack) <- aspect.bi

#binary at 50
extract_aspect_bi_1ha <- exactextractr::exact_extract(aspect.bi, wldf_50, c("sum", "count"))
Perc.NtoE_1ha <- extract_aspect_bi_1ha$sum/extract_aspect_bi_1ha$count
Perc.NtoE_1ha <- as.data.frame(Perc.NtoE_1ha)
Perc.NtoE_1ha$Binary.NtoE_1ha <- ifelse (Perc.NtoE_1ha$Perc.NtoE_1ha >= 0.5, 1, 0)

#binary at 100
extract_aspect_bi_4ha <- exactextractr::exact_extract(aspect.bi, wldf_100, c("sum", "count"))
Perc.NtoE_4ha <- extract_aspect_bi_4ha$sum/extract_aspect_bi_4ha$count
Perc.NtoE_4ha <- as.data.frame(Perc.NtoE_4ha)
Perc.NtoE_4ha$Binary.NtoE_2020_4ha <- ifelse (Perc.NtoE_4ha$Perc.NtoE_4ha >= 0.5, 1, 0)

# dput(colnames(Perc.NtoE_1ha)) #add dummy year [change to 2018 for immutable variables]

colnames(Perc.NtoE_1ha)<-c("Perc.NtoE_2020_1ha", "Binary.NtoE_2020_1ha")
colnames(Perc.NtoE_4ha)<-c("Perc.NtoE_2020_4ha", "Binary.NtoE_2020_4ha")

aspect.NE<-cbind(Perc.NtoE_1ha,Perc.NtoE_4ha)

#get metadata columns
metadatavars<-wldf_50[,1:7]
st_geometry(metadatavars) <- NULL
#collate into wide dataset
#Put all together for 1ha 4ha 

wide_forest_variables<-cbind(metadatavars, canopy, fuels, nbr, elevation, aspect.NE)

#------------------------------------------------------------
# -   Make clean tall dataset
#------------------------------------------------------------

wide_forest_variables %>% #make into a 'long or tall' dataset
  pivot_longer(
    cols = contains("ha"),
    names_to = c("sum_fn", "var", "Year", "scale"),
    # min.Elevation_NA_4ha
    names_pattern = "(.*)\\.(.*)_(.*)_(.*)",
    values_to = "value",
    values_drop_na = TRUE
  ) %>% 
  mutate(point_d=as.factor(point_d)) %>% 
  # mutate(Treatment=Tretmnt) %>% 
  separate(col = SZ_DNS2, into=c("Size","Density"), sep=1) %>% 
  mutate(Density = recode_factor(Density, SP = "Sparse", M = "Moderate", D = "Dense")) %>%
  dplyr::select(point_d, sum_fn, var, Year, scale, value, 
                Cpls_Wt, Treatment, Size, Density)-> tall_forest_variables

NotIn <- function(x,y) !(x %in% y)

tall_forest_variables %>%
  dplyr::select(point_d, sum_fn, var, Year, scale, value) %>%
  filter(sum_fn %in% c("max", "mean","Perc","Binary"), scale=='4ha', NotIn(var, c("dNBR","NA"))) %>% 
  pivot_wider(names_from = sum_fn:scale,
              #names_glue ="{var}_{.value}", #printf('{%s}_{%s}_{%s', sum_fn, scale, Year
              values_from = value
  ) -> wide_forest_variables_mean_perc_bin_4ha

# write_clip(wide_forest_variables_mean_perc_bin_4ha) #[ ] should point directly to output (.csv and equivalent google sheet) [ ] FIX THIS

tall_forest_variables %>%
  dplyr::select(point_d, sum_fn, var, Year, scale, value) %>%
  filter(sum_fn %in% c("max", "mean","Perc","Binary"), scale=='1ha', NotIn(var, c("dNBR","NA"))) %>% 
  pivot_wider(names_from = sum_fn:scale,
              #names_glue ="{var}_{.value}", #printf('{%s}_{%s}_{%s', sum_fn, scale, Year
              values_from = value
  ) -> wide_forest_variables_mean_perc_bin_1ha

#------------------------------------------------------------
# -   conversion of continuous circular data is more complex
# -   [X] FIX not working 2021-10-28, 
#     [ ] would be good to do in previous script
# -   Export final converted data
#------------------------------------------------------------

#aspect image stack (in degrees) (non-converted raster)
asp_imgStack<-dem_imgStack$aspect
extract_asp <- exactextractr::exact_extract(asp_imgStack, wldf_100)
extract_asp[[1]]$value

#convert degrees to a circular object and then to radians 
#build a function to simplify converting degrees to radians with zero degrees/radians = North 
con_circle<-function(degs) {
  conversion.circular(       #do conversion to radians
    circular(degs, units = "degrees", zero = pi/2, rotation="clock"), #of a circular object with degrees with N as zero
    units="radians") #circular frame is inherited
}

#now get mean or other circular stats from converted values (in radians, 0 up top, clockwise compass)
#define aspect Raster
asp_rad_imgStack<-asp_imgStack
#calculate the updated values in radians
circ_vals <- as.vector(con_circle(getValues(asp_rad_imgStack)))
#apply to raster, replacing degrees with the radians
values(asp_rad_imgStack) <- circ_vals

#checking values line up, at record 88888 (which is conveniently 1/2 a circle = ~180 degrees)
asp_imgStack$aspect[88888] #180.3275
asp_rad_imgStack$aspect[88888] #3.147309 !

#sanity plots
# plot(asp_imgStack)
# plot(asp_rad_imgStack)

#recall these are NOT in circular form (though they were calculated with it)!

#extract radian data for points
extract_rad_asp_100 <- exactextractr::exact_extract(asp_rad_imgStack, wldf_100)

#need to calculate the weighted mean for each value, perhaps do the WM then convert to circular, or just do by hand:
extract_rad_asp_100[[1]]$value%*%extract_rad_asp_100[[1]]$coverage_fraction/length(extract_rad_asp_100[[1]]$value)
#          [,1]
# [1,] 4.417106
# looks great!

#this can be calculated as follows (test by hand below) already in radians just need circular specs
circular(weighted.mean(extract_rad_asp_100[[1]]$value, extract_rad_asp_100[[1]]$coverage_fraction),units = 'radians', zero=pi/2, rotation="clock")
#[1] 4.875364 -- works so ...

#wasn't able to get reasonable numbers out of the
## circular::weighted.mean() == weighted.mean.circular() function [ ]#worth double checking why
#make the weighted mean circular function 
wm<-function(.x){circular(stats::weighted.mean(x=.x$value, w=.x$coverage_fraction), units = 'radians', zero=pi/2, rotation="clock")}

purrr::map(extract_rad_asp_100, wm)[[1]]

# Circular Data: 
# Type = angles 
# Units = radians 
# Template = none 
# Modulo = asis 
# Zero = 1.570796 
# Rotation = clock 
# [1] 4.875364

#it works!

#use map to do all the weighted means at 1x 
wmDirVals_100<-purrr::map(extract_rad_asp_100,wm)
as.data.frame(wmDirVals_100) %>% t() -> wmDirVals_100 #sloppy, had to transpose (t())

#add to the dataset
wide_forest_variables_mean_perc_bin_4ha$mean_aspCirc_2020_4ha<-wmDirVals_100[,1]

# final test
wide_forest_variables_mean_perc_bin_4ha$mean_aspCirc_2020_4ha[1]
# structure.4.87536440614317..circularp...list.type....angles... 
#                                                       4.875364 

extract_rad_asp_50 <- exactextractr::exact_extract(asp_rad_imgStack, wldf_50)

#use map to do all the weighted means at 1x 
wmDirVals_50<-purrr::map(extract_rad_asp_50,wm)
as.data.frame(wmDirVals_50) %>% t() -> wmDirVals_50

wide_forest_variables_mean_perc_bin_1ha$mean_aspCirc_2020_1ha<-wmDirVals_50[,1]

#convert back for sanity (and readability)
con_circular_deg <- function(rads) {conversion.circular(circular(rads, units = "radians", zero=pi/2, rotation="clock"),units="degrees")} #not used?

con_circular_deg(wide_forest_variables_mean_perc_bin_1ha$mean_aspCirc_2020_1ha)[[1]]
#[1] 302.508, yay looks good and very close to the median

#------------------------------------------------------------
# Final all-at-once pulling of these variables see 
# https://github.com/calacademy-research/COMB/issues/13
#------------------------------------------------------------
# #The 16 standard variables, NOT USED
# std_var <- names(canopy_fuel_nbr_dem_RAVG_LIDAR)[c(1:3, 7:15, 25:27, 30,32)]
# # [ ] need to calculate the circular weighted mean aspect (var 31) separately see above
# 
# std_layers <- subset(canopy_fuel_nbr_dem_RAVG_LIDAR, std_var)
# 
# dput(names(std_layers))
# fixnames <- c("CanopyBaseHeight_2018", "CanopyBaseHeight_2019", "CanopyBaseHeight_2020", 
#   "CanopyLayerCount_2018", "CanopyLayerCount_2019", "CanopyLayerCount_2020", 
#   "CanopyCover_2018", "CanopyCover_2019", "CanopyCover_2020", 
#   "CanopyHeight_2018", "CanopyHeight_2019", "CanopyHeight_2020", 
#   "CaplesdNBR_Nov18Nov19", "CaplesNBR_Nov18", "CaplesNBR_Nov19", 
#   "elevation_181920", "RAVGdnbr_2018111820191118"
# )
# 
# #fix the names
# names(std_layers) <- fixnames
# 
# #then do the standard variable extraction
# extract_std_var_1ha <- exactextractr::exact_extract(std_layers, wldf_50, c("mean","median","min","max","count"))
# extract_std_var_4ha <- exactextractr::exact_extract(std_layers, wldf_100, c("mean","median","min","max","count"))

#And add in the 
#special_vars <- c("RAVGdnbr_2018111820191118", LargeTreeHeightFraction","LargeTreeCoverFraction")

#RAVG
RAVG_var <- names(canopy_fuel_nbr_dem_RAVG_LIDAR)[c(37)]
RAVG_layers <- subset(canopy_fuel_nbr_dem_RAVG_LIDAR, RAVG_var)
names(RAVG_layers) <- c("RAVGrdnbrcbi4_20182019")

extract_RAVG_var_1ha <- exactextractr::exact_extract(RAVG_layers, wldf_50, c("mean","median","min","max","count"))
extract_RAVG_var_4ha <- exactextractr::exact_extract(RAVG_layers, wldf_100, c("mean","median","min","max","count"))

mean_RAVGcbi4_20182019_1ha <- extract_RAVG_var_1ha$mean
min_RAVGcbi4_20182019_1ha <- extract_RAVG_var_1ha$min
max_RAVGcbi4_20182019_1ha <- extract_RAVG_var_1ha$max

mean_RAVGcbi4_20182019_4ha <- extract_RAVG_var_4ha$mean
min_RAVGcbi4_20182019_4ha <- extract_RAVG_var_4ha$min
max_RAVGcbi4_20182019_4ha <- extract_RAVG_var_4ha$max


#split into named variables

#Large trees
#LargeTreeHeightFraction == the fraction of x_ha of CanopyHeight > 28
#function what % of trees are > cutoff
Perc_LargeTreeHeight_fn <- function(.canopy_ht, lg_tree_cut = 28){
  as.numeric(.canopy_ht$value >= lg_tree_cut)%*%.canopy_ht$coverage_fraction/length(.canopy_ht$coverage_fraction)
}

#hard code each year as function doesn't play nice with different layers at 1x
Perc_LargeTreeHeight_2018_1ha <- std_layers[[10]] %>%
exactextractr::exact_extract(., wldf_50, Perc_LargeTreeHeight_fn, summarize_df = TRUE)
Perc_LargeTreeHeight_2019_1ha <- std_layers[[11]] %>%
  exactextractr::exact_extract(., wldf_50, Perc_LargeTreeHeight_fn, summarize_df = TRUE)
Perc_LargeTreeHeight_2020_1ha <- std_layers[[12]] %>%
  exactextractr::exact_extract(., wldf_50, Perc_LargeTreeHeight_fn, summarize_df = TRUE)

Perc_LargeTreeHeight_2018_4ha <- std_layers[[10]] %>%
  exactextractr::exact_extract(., wldf_100, Perc_LargeTreeHeight_fn, summarize_df = TRUE)
Perc_LargeTreeHeight_2019_4ha <- std_layers[[11]] %>%
  exactextractr::exact_extract(., wldf_100, Perc_LargeTreeHeight_fn, summarize_df = TRUE)
Perc_LargeTreeHeight_2020_4ha <- std_layers[[12]] %>%
  exactextractr::exact_extract(., wldf_100, Perc_LargeTreeHeight_fn, summarize_df = TRUE)

#Perc_LargeTreeCover == Actual canopy cover of pixels where CanopyHeight >28
Perc_LargeTreeCover <- function(.canopy_ht, .canopy_cov, lg_tree_cut = 28){
  .canopy_cov$value[.canopy_ht$value >= lg_tree_cut]%*%.canopy_cov$coverage_fraction[.canopy_ht$value >= lg_tree_cut]/length(.canopy_cov$coverage_fraction[.canopy_ht$value >= lg_tree_cut])
}

#test
# canopy_ht <- exactextractr::exact_extract(std_layers$CanopyHeight_2018, wldf_50)
# canopy_cov <- exactextractr::exact_extract(std_layers$CanopyCover_2018, wldf_50)

#works OK for the 'by hand' function on parallel elements of two lists
#canopy_cov[[3]]$value[canopy_ht[[3]]$value >= 28]%*%canopy_cov[[3]]$coverage_fraction[canopy_ht[[3]]$value >= 28]/length(canopy_cov[[3]]$coverage_fraction[canopy_ht[[3]]$value >= 28])
#equivalent to function
#LargeTreeCoverFraction(canopy_ht[[3]], canopy_cov[[3]])

# [ ] not sure why map function doesn't work, fix later
# purrr::map2_chr(.x = canopy_ht, .y = canopy_cov, .f = LargeTreeCoverFraction(.canopy_ht = .canopy_ht, .canopy_cov = .canopy_cov))

#hardcode for each year 2018:2020
# [ ] bad practice, having problems with mapX family of fcns for these raster extraction lists
#for year, 2018:
canopy_ht_2018_1ha <- exactextractr::exact_extract(std_layers$CanopyHeight_2018, wldf_50)
canopy_cov_2018_1ha <- exactextractr::exact_extract(std_layers$CanopyCover_2018, wldf_50)

Perc_LargeTreeCover_2018_1ha <- 1:length(canopy_ht_2018_1ha)
for(i in 1:length(canopy_ht_2018_1ha)) Perc_LargeTreeCover_2018_1ha[i] <- Perc_LargeTreeCover(canopy_ht_2018_1ha[[i]], canopy_cov_2018_1ha[[i]])

#for year, 2019:
canopy_ht_2019_1ha <- exactextractr::exact_extract(std_layers$CanopyHeight_2019, wldf_50)
canopy_cov_2019_1ha <- exactextractr::exact_extract(std_layers$CanopyCover_2019, wldf_50)

Perc_LargeTreeCover_2019_1ha <- 1:length(canopy_ht_2019_1ha)
for(i in 1:length(canopy_ht_2019_1ha)) Perc_LargeTreeCover_2019_1ha[i] <- Perc_LargeTreeCover(canopy_ht_2019_1ha[[i]], canopy_cov_2019_1ha[[i]])

#for year, 2020:
canopy_ht_2020_1ha <- exactextractr::exact_extract(std_layers$CanopyHeight_2020, wldf_50)
canopy_cov_2020_1ha <- exactextractr::exact_extract(std_layers$CanopyCover_2020, wldf_50)

Perc_LargeTreeCover_2020_1ha <- 1:length(canopy_ht_2020_1ha)
for(i in 1:length(canopy_ht_2020_1ha)) Perc_LargeTreeCover_2020_1ha[i] <- Perc_LargeTreeCover(canopy_ht_2020_1ha[[i]], canopy_cov_2020_1ha[[i]])

#4ha
#foryear, 2018:
canopy_ht_2018_4ha <- exactextractr::exact_extract(std_layers$CanopyHeight_2018, wldf_50)
canopy_cov_2018_4ha <- exactextractr::exact_extract(std_layers$CanopyCover_2018, wldf_50)

Perc_LargeTreeCover_2018_4ha <- 1:length(canopy_ht_2018_4ha)
for(i in 1:length(canopy_ht_2018_4ha)) Perc_LargeTreeCover_2018_4ha[i] <- Perc_LargeTreeCover(canopy_ht_2018_4ha[[i]], canopy_cov_2018_4ha[[i]])

#for year, 2019:
canopy_ht_2019_4ha <- exactextractr::exact_extract(std_layers$CanopyHeight_2019, wldf_50)
canopy_cov_2019_4ha <- exactextractr::exact_extract(std_layers$CanopyCover_2019, wldf_50)

Perc_LargeTreeCover_2019_4ha <- 1:length(canopy_ht_2019_4ha)
for(i in 1:length(canopy_ht_2019_4ha)) Perc_LargeTreeCover_2019_4ha[i] <- Perc_LargeTreeCover(canopy_ht_2019_4ha[[i]], canopy_cov_2019_4ha[[i]])

#for year, 2020:
canopy_ht_2020_4ha <- exactextractr::exact_extract(std_layers$CanopyHeight_2020, wldf_50)
canopy_cov_2020_4ha <- exactextractr::exact_extract(std_layers$CanopyCover_2020, wldf_50)

Perc_LargeTreeCover_2020_4ha <- 1:length(canopy_ht_2020_4ha)
for(i in 1:length(canopy_ht_2020_4ha)) Perc_LargeTreeCover_2020_4ha[i] <- Perc_LargeTreeCover(canopy_ht_2020_4ha[[i]], canopy_cov_2020_4ha[[i]])

#------------------------------------------------------------
# to make final objects "wide1havars" add into 
# wide_forest_variables_mean_perc_bin_1ha 
# and "wide4havars"
# wide_forest_variables_mean_perc_bin_4ha
#------------------------------------------------------------
wide1havars <- cbind(wide_forest_variables_mean_perc_bin_1ha,
      mean_RAVGcbi4_20182019_1ha,
      min_RAVGcbi4_20182019_1ha,
      max_RAVGcbi4_20182019_1ha,
      Perc_LargeTreeHeight_2018_1ha,
      Perc_LargeTreeHeight_2019_1ha,
      Perc_LargeTreeHeight_2020_1ha,
      Perc_LargeTreeCover_2018_1ha,
      Perc_LargeTreeCover_2019_1ha,
      Perc_LargeTreeCover_2020_1ha)

wide4havars <- cbind(wide_forest_variables_mean_perc_bin_4ha,
      mean_RAVGcbi4_20182019_4ha,
      min_RAVGcbi4_20182019_4ha,
      max_RAVGcbi4_20182019_4ha,
      Perc_LargeTreeHeight_2018_4ha,
      Perc_LargeTreeHeight_2019_4ha,
      Perc_LargeTreeHeight_2020_4ha,
      Perc_LargeTreeCover_2018_4ha,
      Perc_LargeTreeCover_2019_4ha,
      Perc_LargeTreeCover_2020_4ha)

#------------------------------------------------------------
# output wide forest variables for Mary & Angela
# 1. wide_forest_variables_mean_perc_bin_1ha, wide_forest_variables_mean_perc_bin_4ha
# -> /spatial/output/tables
# drive_sync(here("spatial","output","tables"), drive_folder = drive_ls("FIX")$id[2][1])
#------------------------------------------------------------

# Write diff to output folder
if (dir.exists(here("spatial/output/tables")) == F) {
  dir.create(here("spatial/output/tables"))
}

wide1havars_filename <- paste0("wide1havars", ".csv")
fwrite(wide1havars, here("spatial", "output", "tables", wide1havars_filename))

wide4havars_filename <- paste0("wide4havars", ".csv")
fwrite(wide4havars, here("spatial", "output", "tables", wide4havars_filename))

#sync up to google drive, if not already, 'sunk' :)
drive_sync(here("spatial","output","tables"), drive_folder = as_id(drive_ls("https://drive.google.com/drive/folders/1vNSUlaz34OA4uYNDjOkT_lToY60y0BPo")))
#currently ^ not working for upload, [ ]?'

#------------------------------------------------------------
# Make final object tall for analysis
#------------------------------------------------------------
wide1havars %>% #make into a 'long or tall' dataset
  pivot_longer(
    cols = contains("ha"),
    names_to = c("sum_fn", "var", "Year", "scale"),
    # min.Elevation_NA_4ha
    names_pattern = "(.*)\\.(.*)_(.*)_(.*)",
    values_to = "value",
    values_drop_na = TRUE
  ) %>% tail()
  mutate(point_d=as.factor(point_d)) %>% 
  # mutate(Treatment=Tretmnt) %>% 
  separate(col = SZ_DNS2, into=c("Size","Density"), sep=1) %>% 
  mutate(Density = recode_factor(Density, SP = "Sparse", M = "Moderate", D = "Dense")) %>%
  dplyr::select(point_d, sum_fn, var, Year, scale, value, 
                Cpls_Wt, Treatment, Size, Density)-> tall_forest_variables

