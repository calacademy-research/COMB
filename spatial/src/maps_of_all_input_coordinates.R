#QUALITY CHECKING OF ALL LOCATION DATA SHARED WITH KAPAN 2022-09-14
# maps_of_all_input_coordinates.R
# 0. Getting started
# 
# install.packages("rgdal")
# install.packages("mapview")
library(rgdal)
library(sf)
library(here)
library(tidyverse)
library(raster)
library(ggpubr)
library(leaflet)
library(mapview)

#1. Read in wildlife points (the original!)
#Import original wildlife points as a shape file (all on the 74xx46.3 UTM_E & 428xx41 UTM_N grids -- so == bird points)
#These were likely shared by e-mail but are consistent with all 'bird-points'.
WildlifePoints.shp <- sf::read_sf(here("spatial", "input", "shapefiles", "WildlifePoints.shp")) %>% # crs not included
  st_set_crs(26910) %>% #were collected using NAD83 coordinate (26910) coordinate reference system
  st_transform(crs(study_area)) #, aoi = study_area$geometry) #defined in other script
#[ ] NOTE PROBLEM NO TRANSFORMATION IS HAPPENING BETWEEN NAD83 AND WGS84 NEED TO FIGURE OUT WHAT'S UP

#FIX NAMES
names(WildlifePoints.shp) <- c(
  "point_d", "Cpls_Wt", "VEG_CSE", "AVIAN_S", "WHR_TSD", "SZ_DNS2",
  "Treatment", "geometry"
)

#save NAD83crs for later 
st_crs(WildlifePoints.shp) -> NAD83crs

# are points in the fire boundary area?
WildlifePoints.shp$inside_fire_boundary <- as.vector(st_intersects(fire_boundary, WildlifePoints.shp, sparse = FALSE))

##[ ] update wild_points to include plant plots, see e-mail 2022-06-28 Sarah Jacobs
#Now adding in a 'fix' from Sarha with "plotID_UTM.csv" that is now in the input directory (so running "read_points_output_data.R" creates it)
plotID_UTM <- read_csv(here("spatial", "input", "shapefiles", "plotID_UTM.csv"))  #from Sarah Jacobs 2022-06-28 see email
  
# plotID_UTM$plotID_av are the bird_points
# plotID_UTM$plotID_veg are the veg_points
# Northing & Easting from plotID_UTM are wierd 
# but is a good join table for plot numbers to names...

#building a bigger table
new_wild_points <- left_join(plotID_UTM, WildlifePoints.shp, by = c("plotID_av" = "point_d"), keep = TRUE)

#remake into sf add back crs
new_wild_points %>% 
  st_as_sf() %>%
  st_set_crs(crs(study_area)) -> new_wild_points # transformed crs

#replace the 87th point (avian # 1072) that has an empty geometry see: 
# https://gis.stackexchange.com/questions/244756/edit-sf-point-conditionally 
#scroll all the way down

#fix the 87th point and set the CRS
new_point <- st_point(c(743246.3, 4287041.4)) %>% 
  st_sfc(crs = 26910) 

new_point <- st_transform(new_point, crs(study_area)) # transformed crs

#conditionally replace it
new_wild_points %>% 
  mutate(geometry = st_sfc(ifelse(new_wild_points$plotID_av==1072, st_geometry(new_point), geometry))) %>%
  st_as_sf() %>%
  st_set_crs(crs(study_area)) -> new_wild_points_WildlifePoints.shp # transformed crs, locations of WildlifePoints.shp

#add in coordinates to regular data
new_wild_points_WildlifePoints.shp %>%
  cbind(., st_coordinates(.)) %>%
  mutate(Easting = X, Northing = Y) %>%
  dplyr::select(-X, -Y) -> new_wild_points_WildlifePoints.shp

#could do for all vegetation points (hold off until edited [ ])

crs(new_wild_points_WildlifePoints.shp)

#write out fixed shapefile
sf::st_write(new_wild_points_WildlifePoints.shp, dsn = here("spatial", "input", "shapefiles", "WildlifePointsFixed.shp"), append = FALSE)

new_wild_points_WildlifePoints.shp_df <- cbind(st_drop_geometry(new_wild_points_WildlifePoints.shp), st_coordinates(new_wild_points_WildlifePoints.shp))

# alternatively, add in sarah's geometry
new_wild_points_WildlifePoints.shp_df %>%
  mutate(Northing2 = Northing, Easting2 = Easting) %>%
  st_as_sf(., coords = c("Easting", "Northing")) -> new_wild_points_plotID_UTM #now has locations of plot_ID_UTM

#"FIX" it [ ] fix this step since it's not changing anythign
  st_crs(new_wild_points_plotID_UTM) <- 26910 
  st_transform(new_wild_points_plotID_UTM, crs(study_area))  #make new geometry

#1.5 also do the same for FinalCaplesMonitoringPlots2022 which I am not sure where they came from?
FinalCaplesMonitoringPlots2022 <- sf::read_sf(here("spatial", "input", "shapefiles", "Monitoring2022", "FinalCaplesMonitoringPlots2022.shp")) %>% # crs not included
  st_set_crs(26910) %>% #were collected using NAD83 coordinate (26910) coordinate reference system
  st_transform(crs(study_area)) #, aoi = study_area$geometry) #defined in other script

# 2. Vegetation Data
# 
# Summarize key variables (Sarah)
# using output 'Caples_PlotData_20220225'
# so must source("Caples_Tidy_20220221.R" to make it work)

Caples_PlotData_20220225 %>% 
  st_as_sf(., coords = c("Easting", "Northing"), crs = NAD83crs) %>% #define crs
  st_transform(., crs(study_area)) -> Caples_PlotData_20220225_c #change crs

#add in columns for the coordinates for display
Caples_PlotData_20220225_c <- cbind(Caples_PlotData_20220225_c, st_coordinates(Caples_PlotData_20220225_c))

#select out points for comparisons
Caples_PlotData_20220225_c %>%
  filter(SamplingTime == 0) %>%
  distinct() -> veg_plot_points_for_comparison_Caples_PlotData_20220225

# 3. Bring in field GPS data (from units and iPhone)
install.packages("gpx")
library(gpx)

# get the data
shapedir <- here("spatial", "input", "shapefiles","field_validation_points", "survey_point_by_cluster")
google_file_names <- "https://drive.google.com/drive/folders/1vU0LSMa7dDhfy9dakarnDogw-spo1iB2"
drive_sync(local_dir = shapedir, drive_folder = google_file_names)

gpx_files[1] <- basename(system(paste0("find ", shapedir, " -mindepth 1 -maxdepth 1 ! -type l"), intern = TRUE))

read_gpx(paste0(shapedir,"/",gpx_files[1]))


#view each one
new_wild_points_WildlifePoints.shp %>% 
  mutate(file_name = "new_wild_points_WildlifePoints.shp") %>%
  dplyr::select(file_name,
                plotID_av,
                plotID_veg,
                WHR_TSD,
                SZ_DNS2,
                Easting,
                Northing,
                geometry)  %>%
  mapView(zcol = "plotID_av") +
new_wild_points_plotID_UTM %>% 
  mutate(file_name = "new_wild_points_plotID_UTM") %>%
  dplyr::select(file_name,
                plotID_av,
                plotID_veg,
                WHR_TSD,
                SZ_DNS2,
                Easting = Easting2,
                Northing = Northing2,
                geometry)  %>%
  mapView(zcol = "plotID_veg") +
FinalCaplesMonitoringPlots2022 %>% 
  mutate(file_name = "FinalCaplesMonitoringPlots2022") %>%
  dplyr::select(file_name,
                Avian_Point = Avian_Poin, 
                CSE_ID, 
                RedFir_ID, 
                Notes, 
                Easting = UTM_E,
                Northing = UTM_N,
                geometry) %>%
  mapView(zcol = "Avian_Point") +
veg_plot_points_for_comparison_Caples_PlotData_20220225 %>% 
  mutate(file_name = "Caples_PlotData_20220225") %>%
  dplyr::select(file_name,
                PlotID, 
                Veg_Type, 
                Elevation, 
                Burned_Caples, 
                Caples_Severity, 
                Caples_Severity_Class, 
                Total_Tree_Cover, 
                Total_Tree_Ht, 
                Easting = X,
                Northing = Y,
                geometry) %>%
  mapView(zcol = "PlotID") +
fire_boundary %>%
  mapView(alpha.regions = 0) #+
# vwf_11.37 %>%
#   mapView(alpha.regions = .3)

# not working, need to TILE it [ ]
# canopy_fuel_nbr_dem_RAVG_LIDAR$CaplesCanopyCover2018 %>%
#   viewRGB()


# 3. Merge new_wild_points and vegetation data based on discussion
# calculate things that are messed up ... and FIX [ ]

#drop geometry from both:
FinalCaplesMonitoringPlots2022_df <- cbind(st_drop_geometry(FinalCaplesMonitoringPlots2022), st_coordinates(FinalCaplesMonitoringPlots2022))
new_wild_points_df <- cbind(st_drop_geometry(new_wild_points), st_coordinates(new_wild_points))
veg_plot_points_for_comparison_df <- cbind(st_drop_geometry(veg_plot_points_for_comparison_Caples_PlotData_20220225), st_coordinates(veg_plot_points_for_comparison_Caples_PlotData_20220225))

veg_plot_points_for_comparison_df %>% setNames(make.names(names(.), unique = TRUE)) -> veg_plot_points_for_comparison_df

#merge data based on shared key(s)
FinalCaplesMonitoringPlots2022_df %>%
  left_join(new_wild_points_df, c("CSE_ID" = "plotID_veg"), keep = TRUE) %>% 
  mutate(CSE_RF_ID = ifelse(is.na(RedFir_ID), CSE_ID, RedFir_ID)) %>% 
  left_join(veg_plot_points_for_comparison_df, by = c("CSE_RF_ID" = "PlotID"),  keep = TRUE) %>% 
    dplyr::select(UTM_N_FCMP22 = UTM_N, UTM_E_FCMP22 = UTM_E, 
                FCMP22_Y.y = Y.x, FCMP22_X.x = X.x, 
                UTM_N_nwp_Northing = Northing, UTM_N_nwp_Easting = Easting,
                nwp_Y.y = Y.y, nwp_X.y = X.y,
                vppfc_Y = Y, vppfc_X = X, CSE_ID:RedFir_ID, Cpls_Wt:inside_fire_boundary, PV_ID:WYMO) %>% #yikes, they are all over the place!
  filter(!is.na(UTM_N_nwp_Easting)) %>%
  st_as_sf(., coords = c("UTM_N_nwp_Easting", "UTM_N_nwp_Northing")) -> fcmp22data #make first two into the geometry

st_crs(fcmp22data) <- 26910
fcmp22data <- st_transform(fcmp22data, crs(study_area))

# 4. Comparisons
# a. Satellite data
# using processed data after running read_adjust_combine_rasters.R
# add in 11.3 m^2 (extract and summarize) for some variables ...
#USE ^^^ to see how far off everything is

#extract data from each source and do the correlations
#start with FinalCaplesMonitoringPlots2022

vwf_11.37 <- fcmp22data  %>% sf::st_buffer(11.37,endCapStyle = "ROUND") 

canopy_fuel_nbr_dem_RAVG_LIDAR %>%
  exactextractr::exact_extract(., vwf_11.37, c("mean", "median", "min", "max", "count")) -> canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m

canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m %>% dim()

dim(fcmp22data)[[1]] == dim(canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m)[[1]]

#[ ] need a validation step here to make sure it is lined up (too bad the extractr doesn't give rownames for objects!)

cbind(fcmp22data,canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m) -> veg_sat_comp

veg_sat_comp$geometry

veg_sat_comp$UTM_N <- sf::st_coordinates(veg_sat_comp$geometry)[,2]
veg_sat_comp$UTM_E <- sf::st_coordinates(veg_sat_comp$geometry)[,1]

veg_sat_comp$Caples_Severity_Class_dRAVG = cut(veg_sat_comp$mean.ca3872412014620191010_20181118_20191118_rdnbr_cbi4, breaks = c(-1, 0.5, 2, 3, 4), labels=c(0,1,2,3))

#adjust levels for severity class
# dput(unique(veg_sat_comp$Caples_Severity_Class))
veg_sat_comp$Caples_Severity_Class <- ordered(veg_sat_comp$Caples_Severity_Class, levels = c("High", "Mod", "Low", "Unburned", NA))

# #analysis of tree height
# #histograms to see if they make sense
# veg_sat_comp$MaxTreeHeight %>% hist()
# veg_sat_comp$mean.CaplesCanopyHeight2018 %>% hist()

#linear model of one versus the other 
reg<-lm(formula = max.CaplesCanopyHeight2018 ~ MaxTreeHeight + 0,
        data=veg_sat_comp)                      

#get intercept and slope value
coeff<-coefficients(reg)          
slope<-coeff[1]
# slope<- coeff[2]
# intercept<-coeff[1]

reg2<-lm(formula = max.CaplesCanopyHeight2018 ~ MaxTreeHeight,
         data=veg_sat_comp)                      

#get intercept and slope value
coeff<-coefficients(reg2)          
intercept2<-coeff[1]
slope2<-coeff[2]

veg_sat_comp %>%
  filter(!is.na(MaxTreeHeight)) %>%
  filter(Avian_Poin > 0) %>%
  # filter()
  ggplot(.) +
    geom_point(mapping = aes(x = MaxTreeHeight, y = max.CaplesCanopyHeight2018)) +
    geom_abline(mapping = aes(intercept = 0, slope = 1)) + 
    geom_abline(slope = slope, color="red", 
              linetype="dashed", size=1.5) + 
  geom_text(aes(x = MaxTreeHeight, y = jitter(max.CaplesCanopyHeight2018), label = PlotID, hjust = -.25, vjust = -.35)) +
  stat_cor(aes(x = MaxTreeHeight, y = max.CaplesCanopyHeight2018, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 0, label.y = 35) +
  #facet_wrap(~Caples_Severity_Class) +
    #geom_abline(intercept = intercept2, slope = slope2, color="blue", 
    #          linetype="dotted", size=1.5) + 
  #  facet_wrap(~Caples_Severity_Class) +
    labs(y="max Canopy Height before fire (2018 satellite)", x = "Max Tree Height before fire (ground)", title="Correlation between on-the-ground and satellite measurements\n height for vegetation plot (~ 400 meters squared) veg points") 
  

summary(reg)  
summary(reg2)  

anova(reg,reg2)

###REDO ABOVE FOR BIRD POINT LOCS (also not correct)
#merge data based on shared key(s)
FinalCaplesMonitoringPlots2022_df %>%
  left_join(new_wild_points_df, c("CSE_ID" = "plotID_veg"), keep = TRUE) %>% 
  mutate(CSE_RF_ID = ifelse(is.na(RedFir_ID), CSE_ID, RedFir_ID)) %>% 
  left_join(veg_plot_points_for_comparison_df, by = c("CSE_RF_ID" = "PlotID"),  keep = TRUE) %>% 
  dplyr::select(UTM_N_FCMP22 = UTM_N, UTM_E_FCMP22 = UTM_E, 
                FCMP22_Y.y = Y.x, FCMP22_X.x = X.x, 
                UTM_N_nwp_Northing = Northing, UTM_N_nwp_Easting = Easting,
                nwp_Y.y = Y.y, nwp_X.y = X.y,
                vppfc_Y = Y, vppfc_X = X, CSE_ID:RedFir_ID, Cpls_Wt:inside_fire_boundary, PV_ID:WYMO) %>% #yikes, they are all over the place!
  filter(!is.na(UTM_N_nwp_Easting)) %>%
  st_as_sf(., coords = c("nwp_X.y", "nwp_Y.y")) -> fcmp22data #make first two into the geometry

st_crs(fcmp22data) <- 26910
fcmp22data <- st_transform(fcmp22data, crs(study_area))

# 4. Comparisons
# a. Satellite data
# using processed data after running read_adjust_combine_rasters.R
# add in 11.3 m^2 (extract and summarize) for some variables ...
#USE ^^^ to see how far off everything is

#extract data from each source and do the correlations
#start with FinalCaplesMonitoringPlots2022

vwf_11.37 <- fcmp22data  %>% sf::st_buffer(11.37,endCapStyle = "ROUND") 

canopy_fuel_nbr_dem_RAVG_LIDAR %>%
  exactextractr::exact_extract(., vwf_11.37, c("mean", "median", "min", "max", "count")) -> canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m

canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m %>% dim()

dim(fcmp22data)[[1]] == dim(canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m)[[1]]

#[ ] need a validation step here to make sure it is lined up (too bad the extractr doesn't give rownames for objects!)

cbind(fcmp22data,canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m) -> veg_sat_comp

veg_sat_comp$geometry

veg_sat_comp$UTM_N <- sf::st_coordinates(veg_sat_comp$geometry)[,2]
veg_sat_comp$UTM_E <- sf::st_coordinates(veg_sat_comp$geometry)[,1]

veg_sat_comp$Caples_Severity_Class_dRAVG = cut(veg_sat_comp$mean.ca3872412014620191010_20181118_20191118_rdnbr_cbi4, breaks = c(-1, 0.5, 2, 3, 4), labels=c(0,1,2,3))

#adjust levels for severity class
# dput(unique(veg_sat_comp$Caples_Severity_Class))
veg_sat_comp$Caples_Severity_Class <- ordered(veg_sat_comp$Caples_Severity_Class, levels = c("High", "Mod", "Low", "Unburned", NA))

# #analysis of tree height
# #histograms to see if they make sense
# veg_sat_comp$MaxTreeHeight %>% hist()
# veg_sat_comp$mean.CaplesCanopyHeight2018 %>% hist()

#linear model of one versus the other 
reg<-lm(formula = max.CaplesCanopyHeight2018 ~ MaxTreeHeight + 0,
        data=veg_sat_comp)                      

#get intercept and slope value
coeff<-coefficients(reg)          
slope<-coeff[1]
# slope<- coeff[2]
# intercept<-coeff[1]

reg2<-lm(formula = max.CaplesCanopyHeight2018 ~ MaxTreeHeight,
         data=veg_sat_comp)                      

#get intercept and slope value
coeff<-coefficients(reg2)          
intercept2<-coeff[1]
slope2<-coeff[2]

veg_sat_comp %>%
  filter(!is.na(MaxTreeHeight)) %>%
  filter(Avian_Poin > 0) %>%
  # filter()
  ggplot(.) +
  geom_point(mapping = aes(x = MaxTreeHeight, y = max.CaplesCanopyHeight2018)) +
  geom_abline(mapping = aes(intercept = 0, slope = 1)) + 
  geom_abline(slope = slope, color="red", 
              linetype="dashed", size=1.5) + 
  geom_text(aes(x = MaxTreeHeight, y = jitter(max.CaplesCanopyHeight2018), label = PlotID, hjust = -.25, vjust = -.35)) +
  stat_cor(aes(x = MaxTreeHeight, y = max.CaplesCanopyHeight2018, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 0, label.y = 35) +
  #facet_wrap(~Caples_Severity_Class) +
  #geom_abline(intercept = intercept2, slope = slope2, color="blue", 
  #          linetype="dotted", size=1.5) + 
  #  facet_wrap(~Caples_Severity_Class) +
  labs(y="max Canopy Height before fire (2018 satellite)", x = "Max Tree Height before fire (ground)", title="Correlation between on-the-ground and satellite measurements\n height for vegetation plot (~ 400 meters squared) bird plots") 


summary(reg)  
summary(reg2)  

anova(reg,reg2)







#if we look at the raw corellation, its not so good ... 
cor <- cor.test(veg_sat_comp$Total_Tree_Ht, veg_sat_comp$mean.CaplesCanopyHeight2018)

cor$estimate^2 * 100
# cor 
# 25.90276 

#but in reality we are really asking is how human measurements compare to satellite estimates
#what about cover

#analysis of tree height
#histograms to see if they make sense
veg_sat_comp$Total_Tree_Cover %>% hist()
veg_sat_comp$mean.CaplesCanopyCover2018 %>% hist()

#linear model of one versus the other 
reg<-lm(formula = mean.CaplesCanopyCover2018 ~ Total_Tree_Cover + 0,
        data=veg_sat_comp)                      

#get intercept and slope value
coeff<-coefficients(reg)          
slope<-coeff[1]


reg2<-lm(formula = mean.CaplesCanopyCover2018 ~ Total_Tree_Cover,
        data=veg_sat_comp)                      

#get intercept and slope value
coeff<-coefficients(reg2)          
intercept2<-coeff[1]
slope2<-coeff[2]

veg_sat_comp %>%
  filter(!is.na(Caples_Severity_Class)) %>%
  ggplot(.) +
  geom_point(mapping = aes(x = Total_Tree_Cover, y = mean.CaplesCanopyCover2018)) +
  geom_abline(mapping = aes(intercept = 0, slope = 1)) + 
  geom_abline(slope = slope, color="red", 
              linetype="dashed", size=1.5) + 
  geom_abline(intercept = intercept2, slope = slope2, color="blue", 
              linetype="dotted", size=1.5) + 
  facet_wrap(~Caples_Severity_Class) +
  labs(y="mean Canopy Cover before fire (2018 satellite)", x = "Total Tree Cover before fire (ground)", title="Corelation between on-the-ground and satellite measurements\n cover") 


summary(reg)  
summary(reg2)  

#if we look at the raw corellation, its not so good ... 
cor <- cor.test(veg_sat_comp$Total_Tree_Cover, veg_sat_comp$mean.CaplesCanopyCover2018)

cor$estimate^2 * 100

anova(reg,reg2)

## between years (Before fire 2018 vs 2019)
veg_sat_comp %>%
  filter(!is.na(Caples_Severity_Class)) %>%
  ggplot(.) +
  geom_point(mapping = aes(x = mean.CaplesCanopyHeight2018, y = mean.CaplesCanopyHeight2019)) +
  geom_abline(mapping = aes(intercept = 0, slope = 1)) +
  facet_wrap(~Caples_Severity_Class) 
  

veg_sat_comp %>%
  filter(!is.na(Caples_Severity_Class)) %>%
  ggplot(.) +
  geom_point(mapping = aes(x = mean.CaplesCanopyCover2018, y = mean.CaplesCanopyCover2019)) +
  geom_abline(mapping = aes(intercept = 0, slope = 1)) +
  facet_wrap(~Caples_Severity_Class) 

## after years (2020 vs 2018 or 2019)
veg_sat_comp %>%
  filter(!is.na(Caples_Severity_Class)) %>%
    ggplot(.) +
  geom_point(mapping = aes(x = mean.CaplesCanopyHeight2018, y = mean.CaplesCanopyHeight2020)) +
  geom_abline(mapping = aes(intercept = 0, slope = 1)) +
  facet_wrap(~Caples_Severity_Class) 

veg_sat_comp %>%
  filter(!is.na(Caples_Severity_Class)) %>%
  ggplot(.) +
  geom_point(mapping = aes(x = mean.CaplesCanopyCover2018, y = mean.CaplesCanopyCover2020)) +
  geom_abline(mapping = aes(intercept = 0, slope = 1)) +
  facet_wrap(~Caples_Severity_Class) 

#this suggests looking at RANDOM POINTS might be useful...but at the scale of the point, it is reapeatable!
cor <- cor.test(veg_sat_comp$mean.CaplesCanopyHeight2019, veg_sat_comp$mean.CaplesCanopyHeight2018)

cor$estimate^2 * 100

cor <- cor.test(veg_sat_comp$mean.CaplesCanopyCover2019, veg_sat_comp$mean.CaplesCanopyCover2018)

cor$estimate^2 * 100

###

#Create a data frame from the raster values
df <- cbind.data.frame(values(canopy_fuel_nbr_dem_RAVG_LIDAR$CaplesCanopyHeight2018), values(canopy_fuel_nbr_dem_RAVG_LIDAR$CaplesCanopyHeight2019))
names(df) <- c("CaplesCanopyHeight2018", "CaplesCanopyHeight2019")

#Plotting
library(ggpubr)

ggscatter(data = df, x = "CaplesCanopyHeight2018", y = "CaplesCanopyHeight2019", 
          add = "reg.line") + # Add regression line
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 0, label.y = 1.1)

df1820 <- cbind.data.frame(values(canopy_fuel_nbr_dem_RAVG_LIDAR$CaplesCanopyHeight2018), values(canopy_fuel_nbr_dem_RAVG_LIDAR$CaplesCanopyHeight2020))
names(df1820) <- c("CaplesCanopyHeight2018", "CaplesCanopyHeight2020")

df1820 %>%
  slice_sample(., n = 10^5) %>%
  ggplot(.) +
  geom_point(mapping = aes(x = CaplesCanopyHeight2018, y = CaplesCanopyHeight2020)) +
  geom_abline(mapping = aes(intercept = 0, slope = 1)) +
  stat_cor(aes(x = CaplesCanopyHeight2018, y = CaplesCanopyHeight2020, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 0, label.y = 41.1)

#just for our points at a single 'point' value of 10m diameter (5m radius)
vwf_5 <- veg_plot_points_for_comparison  %>% sf::st_buffer(5,endCapStyle = "SQUARE") 


# 2. Satellite data
# using processed data after running read_adjust_combine_rasters.R
# add in 11.3 m^2 (extract and summarize) for some variables ...
canopy_fuel_nbr_dem_RAVG_LIDAR %>%
  exactextractr::exact_extract(., vwf_5, c("mean", "median", "min", "max", "count")) -> canopy_fuel_nbr_dem_RAVG_LIDAR_extract_vwf_5m 

# add back the plotID
veg_sat_comp_vwf_5m <- cbind(vwf_5,canopy_fuel_nbr_dem_RAVG_LIDAR_extract_vwf_5m)
veg_sat_comp_vwf_5m$Caples_Severity_Class <- ordered(veg_sat_comp_vwf_5m$Caples_Severity_Class, levels = c("High", "Mod", "Low", "Unburned", NA))
canopy_fuel_nbr_dem_RAVG_LIDAR_extract_vwf_5m %>% dim()

## after years (2018 vs 2019)
veg_sat_comp_vwf_5m %>%
  filter(!is.na(Caples_Severity_Class)) %>%
  ggplot(.) +
  geom_point(mapping = aes(x = mean.CaplesCanopyHeight2018, y = mean.CaplesCanopyHeight2019)) +
  geom_abline(mapping = aes(intercept = 0, slope = 1)) +
  stat_cor(aes(x = mean.CaplesCanopyHeight2018, y = mean.CaplesCanopyHeight2019, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 0, label.y = 28) +
  facet_wrap(~Caples_Severity_Class) +
  labs(y="mean Canopy Height before fire (2019 satellite)", 
       x = "mean Canopy Height before fire (2018 satellite)", 
       title="Correlation between years satellite measurements\n canopy height at each point (~1 pixel) pre fire") -> cor1819ht
  
#calc 1 and 2 param models
#linear model of one versus the other 
reg1819ht <- lm(formula = mean.CaplesCanopyHeight2019 ~ mean.CaplesCanopyHeight2018 + 0,
                data=veg_sat_comp_vwf_5m)                      

#get intercept and slope value
coeff<-coefficients(reg1819ht)          
slope<-coeff[1]

#get intercept and slope value
reg1819ht2 <- lm(formula = mean.CaplesCanopyHeight2019 ~ mean.CaplesCanopyHeight2018,
                 data=veg_sat_comp_vwf_5m)                      

#get intercept and slope value
coeff<-coefficients(reg1819ht2)          
intercept2<-coeff[1]
slope2<-coeff[2]

cor1819ht +
  geom_abline(slope = slope, color="red", 
              linetype="dashed", size=1.5) + 
  geom_abline(intercept = intercept2, slope = slope2, color="blue", 
              linetype="dotted", size=1.5) 

anova(reg1819ht,reg1819ht2)

#cover
veg_sat_comp_vwf_5m %>%
  filter(!is.na(Caples_Severity_Class)) %>%
  ggplot(.) +
  geom_point(mapping = aes(x = mean.CaplesCanopyCover2018, y = mean.CaplesCanopyCover2019)) +
  geom_abline(mapping = aes(intercept = 0, slope = 1)) +
  stat_cor(aes(x = mean.CaplesCanopyCover2018, y = mean.CaplesCanopyCover2019, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 0, label.y = 75) +
  facet_wrap(~Caples_Severity_Class) +
  labs(y="mean Canopy Cover before fire (2019 satellite)", 
       x = "mean Canopy Cover before fire (2018 satellite)", 
       title="Correlation between years satellite measurements\n canopy cover at each point (~1 pixel) pre fire") -> cor1819cv
  
#get intercept and slope value
reg1819cv <- lm(formula =  mean.CaplesCanopyCover2019 ~ mean.CaplesCanopyCover2018 + 0,
                  data=veg_sat_comp_vwf_5m)                      
  
#get intercept and slope value
coeff<-coefficients(reg1819cv)          
slope<-coeff[1]
  
#get intercept and slope value
reg1819cv2 <- lm(formula = mean.CaplesCanopyCover2019 ~ mean.CaplesCanopyCover2018,
                   data=veg_sat_comp_vwf_5m)                      
  
#get intercept and slope value
coeff<-coefficients(reg1819cv2)          
intercept2<-coeff[1]
slope2<-coeff[2]
  
cor1819cv +
    geom_abline(slope = slope, color="red", 
                linetype="dashed", size=1.5) + 
    geom_abline(intercept = intercept2, slope = slope2, color="blue", 
                linetype="dotted", size=1.5) 

anova(reg1819cv,reg1819cv2)

#vs 2020
## after years (2020 vs 2018 or 2019)
veg_sat_comp_vwf_5m %>%
  filter(!is.na(Caples_Severity_Class)) %>%
  ggplot(.) +
  geom_point(mapping = aes(x = mean.CaplesCanopyHeight2018, y = mean.CaplesCanopyHeight2020)) +
  geom_abline(mapping = aes(intercept = 0, slope = 1)) +
  stat_cor(aes(x = mean.CaplesCanopyHeight2018, y = mean.CaplesCanopyHeight2020, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 0, label.y = 28) +
  facet_wrap(~Caples_Severity_Class) +
  labs(y="mean Canopy Height after fire (2020 satellite)", 
       x = "mean Canopy Height before fire (2018 satellite)", 
       title="Correlation between years satellite measurements\n canopy height at each point (~1 pixel) before and after fire") -> cor1820ht

#calc 1 and 2 param models
#linear model of one versus the other 
reg1820ht <- lm(formula = mean.CaplesCanopyHeight2020 ~ mean.CaplesCanopyHeight2018 + 0,
                data=veg_sat_comp_vwf_5m)                      

#get intercept and slope value
coeff<-coefficients(reg1820ht)          
slope<-coeff[1]

#get intercept and slope value
reg1820ht2 <- lm(formula = mean.CaplesCanopyHeight2020 ~ mean.CaplesCanopyHeight2018,
                 data=veg_sat_comp_vwf_5m)                      

#get intercept and slope value
coeff<-coefficients(reg1820ht2)          
intercept2<-coeff[1]
slope2<-coeff[2]

cor1820ht +
  geom_abline(slope = slope, color="red", 
              linetype="dashed", size=1.5) + 
  geom_abline(intercept = intercept2, slope = slope2, color="blue", 
              linetype="dotted", size=1.5) 

anova(reg1820ht,reg1820ht2)

#cover 2020 vs 2018
veg_sat_comp_vwf_5m %>%
  filter(!is.na(Caples_Severity_Class)) %>%
  ggplot(.) +
  geom_point(mapping = aes(x = mean.CaplesCanopyCover2018, y = mean.CaplesCanopyCover2020)) +
  geom_abline(mapping = aes(intercept = 0, slope = 1)) +
  stat_cor(aes(x = mean.CaplesCanopyCover2018, y = mean.CaplesCanopyCover2020, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 0, label.y = 75) +
  facet_wrap(~Caples_Severity_Class) +
  labs(y="mean Canopy Cover after fire (2020 satellite)", 
       x = "mean Canopy Cover before fire (2018 satellite)", 
       title="Correlation between years satellite measurements\n canopy Cover at each point (~1 pixel) before and after fire") -> cor1820cv

#calc 1 and 2 param models
#linear model of one versus the other 
reg1820cv <- lm(formula = mean.CaplesCanopyCover2020 ~ mean.CaplesCanopyCover2018 + 0,
                data=veg_sat_comp_vwf_5m)                      

#get intercept and slope value
coeff<-coefficients(reg1820cv)          
slope<-coeff[1]

#get intercept and slope value
reg1820cv2 <- lm(formula = mean.CaplesCanopyCover2020 ~ mean.CaplesCanopyCover2018,
                 data=veg_sat_comp_vwf_5m)                      

#get intercept and slope value
coeff<-coefficients(reg1820cv2)          
intercept2<-coeff[1]
slope2<-coeff[2]

cor1820cv +
  geom_abline(slope = slope, color="red", 
              linetype="dashed", size=1.5) + 
  geom_abline(intercept = intercept2, slope = slope2, color="blue", 
              linetype="dotted", size=1.5) 

anova(reg1820cv,reg1820cv2)

veg_sat_comp_vwf_5m %>%
  filter(!is.na(Caples_Severity_Class)) %>%
  ggplot(.) +
  geom_point(mapping = aes(x = mean.CaplesCanopyHeight2018, y = mean.CaplesCanopyHeight2020)) +
  geom_text(aes(x = mean.CaplesCanopyHeight2018, y = mean.CaplesCanopyHeight2020, label = PlotID, hjust = 1, vjust = -1)) +
  geom_abline(mapping = aes(intercept = 0, slope = 1)) +
  stat_cor(aes(x = mean.CaplesCanopyHeight2018, y = mean.CaplesCanopyHeight2020, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 0, label.y = 30) +
  #facet_wrap(~Caples_Severity_Class) +
  labs(y="mean Canopy Height after fire (2020 satellite)", 
       x = "mean Canopy Height before fire (2018 satellite)", 
       title="Correlation between years satellite measurements\n canopy height at each point (~1 pixel) before and after fire") # -> cor1820cv

