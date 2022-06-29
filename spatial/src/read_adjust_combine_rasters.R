#### Fix the rasters: read_adjust_combine_rasters.R ####
# read_adjust_combine_rasters.R
# Created by: Durrell Kapan
# Created on: 5 April 2022
# modified from Caples_Spatial_EDA.md
# see that file for details on creating inputs

#### Load Libraries ####

# load libraries
library(here)
library(tidyverse)
library(googledrive) # connect to google drive
library(purrr)
library(raster)
library(rgdal)

# load functions
source(here("comb_functions.R"))

#------------------------------------------------------------
# Download/organize input files from google drive if not already downloaded
#------------------------------------------------------------
# this might be a bit buggy, try it out @MARY
# define rastertype subdirectories
rastertypes <- c("SALO", "NBR", "RAVG", "LIDAR", "DEM")
# define base directory for input rasters
rasterdir <- here("spatial", "input", "rasters")
# list of directories to use
rasterdirs <- here(rasterdir, rastertypes[])
# add in 'super' directories
rasterdirs <- c("/Users/dkapan/GitHub/COMB/spatial/input/", "/Users/dkapan/GitHub/COMB/spatial/input/rasters/", rasterdirs)
# Creating the folder within inputs that contains the raw_files
map(.x = rasterdirs, .f = function(.x) {
  if (dir.exists(.x) == F) {
    dir.create(.x)
  }
})
# wasn't able to write ^ as anonymous function, too tired.
# Getting list of files that need to be downloaded
google_raw_file_names <- drive_ls("https://drive.google.com/drive/folders/1fN22kGBoWF03p1cNZ7h_MbJf0X-CSyda", recursive = TRUE)
local_file_names <- list.files(rasterdir, pattern = "*", recursive = TRUE)
# Filtering google files that are not already downloaded
needed_google_raw_file_names <- google_raw_file_names %>% filter(!(google_raw_file_names$name %in% local_file_names))
# Downloading raw_data csv files from google drive if not already downloaded
map2(
  needed_google_raw_file_names$name, needed_google_raw_file_names$id,
  ~ drive_sync(local_dir = here(rasterdir, .x), drive_folder = as_id(.y))
)
# if you hard-coded solution it would look like this
# first go into each rasterdir
# rasterdirs[1]
# [1] "/Users/dkapan/GitHub/Caples_Spatial/spatial/input/rasters/SALO"
# don't forget to match by hand to the appropriate google_drive input directory, e.g.
# drive_sync(local_dir = here(rasterdirs[1]), drive_folder = as_dribble("https://drive.google.com/drive/folders/1SKzKZNMYpDmzcUCACOjtEJCLybHN6wgQ"))
# but the overall code above seems to work, so no need! :)

rasterdir <- rasterdirs[2 + 1] # SALO the 2+ skips the first two 'super directories'
canopy_files <- list.files(rasterdir, pattern = "Canopy", full.names = TRUE)
canopy_imgStack <- stack(canopy_files)

fuel_files <- list.files(rasterdir, pattern = "Fuel", full.names = TRUE)
fuel_imgStack <- stack(fuel_files)

# NBR from GEE code
rasterdir <- rasterdirs[2 + 2] # NBR
nbr_files <- list.files(rasterdir, pattern = "NBR", full.names = TRUE)
nbr_imgStack <- stack(nbr_files[1:8])
# fix CRS by re-projection
nbr_imgStack <- projectRaster(nbr_imgStack, crs = crs(canopy_imgStack))

rasterdir <- rasterdirs[2 + 3] # RAVG
# RAVG from BE Caples
RAVG_Caples_files <- list.files(rasterdir, pattern = "ca3872412014620191010_20181118_20191118", full.names = TRUE)
RAVG_Caples_imgStack <- stack(RAVG_Caples_files)

# RAVG from BE Caldor
RAVG_Caldor_files <- list.files(rasterdir, pattern = "ca3858612053820210815_20201011_20211016", full.names = TRUE)
RAVG_Caldor_imgStack <- stack(RAVG_Caldor_files)

# same CRS (reproject below)
identical(crs(RAVG_Caples_imgStack), crs(RAVG_Caldor_imgStack))

# fix CRS by reprojection
RAVG_Caples_imgStack <- projectRaster(RAVG_Caples_imgStack, crs = crs(canopy_imgStack))
RAVG_Caldor_imgStack <- projectRaster(RAVG_Caldor_imgStack, crs = crs(canopy_imgStack))

extent(RAVG_Caldor_imgStack)
extent(RAVG_Caples_imgStack)

# get extent for one, set for both
RAVG_extent <- extent(RAVG_Caples_imgStack)

# crop and set extents
RAVG_Caldor_imgStack <- crop(RAVG_Caldor_imgStack, RAVG_extent)
RAVG_Caldor_imgStack <- setExtent(RAVG_Caldor_imgStack, RAVG_extent, keepres = TRUE)

# stack the stacks!
RAVG_imgStack <- stack(RAVG_Caples_imgStack, RAVG_Caldor_imgStack)

# LIDAR data !!!
# [ ] try to update with higher resolution data
rasterdir <- rasterdirs[2 + 4] # LIDAR
# (same as rasterdir  <- here("spatial","input","Rasters","LIDAR"))

# LIDAR from Becky Estes
LIDAR_files <- list.files(rasterdir, pattern = "*.tif", full.names = TRUE)

# DEM from GEE code
rasterdir <- rasterdirs[2 + 5] # DEM
# alternatively rasterdir  <- here("spatial","input","Rasters")
dem_files <- list.files(rasterdir, pattern = "dem", full.names = TRUE)
dem_imgStack <- stack(dem_files)

# fix CRS by reprojection
dem_imgStack <- projectRaster(dem_imgStack, crs = crs(canopy_imgStack))

# note odd extents (not the same between the three groups below and error on PROJ4 when looking at these RASTERS
# [ ] check to make sure everything lines up at the end)
# use raster::raster(LIDAR_files[1]) to see error
# Warning messages:
# 1: In showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj = prefer_proj) :
# Discarded datum NAD83 (National Spatial Reference System 2011) in Proj4 definition

# hack method to fix LIDAR extents
# stack the files with same extent

LIDAR_int_134 <- raster::stack(LIDAR_files[c(1, 3:4)])
LIDAR_int_2 <- raster::stack(LIDAR_files[2])
LIDAR_int_58 <- raster::stack(LIDAR_files[5:8])

# fix CRS by reprojection
LIDAR_int_134 <- projectRaster(LIDAR_int_134, crs = crs(canopy_imgStack))
LIDAR_int_2 <- projectRaster(LIDAR_int_2, crs = crs(canopy_imgStack))
LIDAR_int_58 <- projectRaster(LIDAR_int_58, crs = crs(canopy_imgStack))

#check extents
extent(LIDAR_int_134)
extent(LIDAR_int_2)
extent(LIDAR_int_58)

# fix extent by 'clipping' using extent.
# set extents to be the max of mins and the mins of maxs so we can stack()
# get extent from 'smallest'
LIDAR_extent <- extent(LIDAR_int_2)

# crop and set extents
LIDAR_int_134 <- crop(LIDAR_int_134, LIDAR_extent)
LIDAR_int_134 <- setExtent(LIDAR_int_134, LIDAR_extent, keepres = TRUE)

LIDAR_int_58 <- crop(LIDAR_int_58, LIDAR_extent)
LIDAR_int_58 <- setExtent(LIDAR_int_58, LIDAR_extent, keepres = TRUE)

#stack
LIDAR_imgStack <- raster::stack(LIDAR_int_134, LIDAR_int_2, LIDAR_int_58)

#fix name for layer 2
names(LIDAR_imgStack)[[4]]<-"Canopy_Rugosity"

names(LIDAR_imgStack)

# **Do calculations**
#
#   -   Digital Elevation Model
#
# -   calculate Aspect

dem_imgStack$aspect <- terrain(dem_imgStack$elevation, opt = "aspect", neighbors = 8, unit = "degrees", filename = here("spatial", "output", "rasters", "caples_aspect.tif"), overwrite = TRUE) # note the aspect must be translated to circular coordinates if using summary functions see points script

# **Do analysis**
#
#   -   Merge rasters\

# merge both stacks with different extents, inspect them first
extent(canopy_imgStack)
extent(fuel_imgStack)
extent(nbr_imgStack)
extent(dem_imgStack)
extent(RAVG_imgStack)
extent(LIDAR_imgStack)

# merge rasters to match
# first downsample NBR and DEM to get same 10m resolution then merge
nbr_imgStack_10m <- resample(nbr_imgStack, canopy_imgStack, method = "bilinear")
dem_imgStack_10m <- resample(dem_imgStack, canopy_imgStack, method = "bilinear")
RAVG_imgStack_10m <- resample(RAVG_imgStack, canopy_imgStack, method = "bilinear")
LIDAR_imgStack_10m <- resample(LIDAR_imgStack, canopy_imgStack, method = "bilinear")

# can follow this code to crop before merging
# get extents
study_area_extent <- extent(canopy_imgStack)

# crop and set extents
nbr_imgStack_10m <- crop(nbr_imgStack_10m, study_area_extent)
nbr_imgStack_10m <- setExtent(nbr_imgStack_10m, study_area_extent, keepres = TRUE)

dem_imgStack_10m <- crop(dem_imgStack_10m, study_area_extent)
dem_imgStack_10m <- setExtent(dem_imgStack_10m, study_area_extent, keepres = TRUE)

RAVG_imgStack_10m <- crop(RAVG_imgStack_10m, study_area_extent)
RAVG_imgStack_10m <- setExtent(RAVG_imgStack_10m, study_area_extent, keepres = TRUE)

LIDAR_imgStack_10m <- crop(LIDAR_imgStack_10m, study_area_extent)
LIDAR_imgStack_10m <- setExtent(LIDAR_imgStack_10m, study_area_extent, keepres = TRUE)

# now working, but not sure I like upscaling the ~9m pixels to 10m ... only for exploratory data analysis
canopy_fuel_nbr_dem_RAVG_LIDAR <- raster::stack(canopy_imgStack, fuel_imgStack, nbr_imgStack_10m, dem_imgStack_10m, RAVG_imgStack_10m, LIDAR_imgStack_10m, orig = FALSE, tolerance = 0.1)
# [ ] check this doesn't corrupt everything ??? never did check, but is only adding a few decimeters to any pixel

# make a place to output the raster brick
outrasterdir <- c("/Users/dkapan/GitHub/COMB/spatial/output/")
# Creating the folder within spatial to contains the processed raster brick
map(.x = outrasterdir, .f = function(.x) {
  if (dir.exists(.x) == F) {
    dir.create(.x)
  }
})

# write to hard drive
writeRaster(canopy_fuel_nbr_dem_RAVG_LIDAR, here("spatial", "output", "rasters", "canopy_fuel_nbr_dem_RAVG_LIDAR.grd"), format = "raster", overwrite = TRUE)

# sync up to google drive, if not already, 'sunk' :)
drive_sync(here("spatial", "output", "rasters"), drive_folder = drive_ls("https://drive.google.com/drive/folders/1sbgR_OMtK-Hq6P6lVBFIK8xbcjmJDQVV")$id[1])
# [ ]final drive sync doesn't work:
# Error in curl::curl_fetch_memory(url, handle = handle) :
# Error in the HTTP2 framing layer
