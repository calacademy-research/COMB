#test_of_original_wild_points.R
#TO DO

# 1. search files, emails & GPS/phones for vegetation and wildlife points
# a. Spreadsheets translated to Google sheets
# # CaplesWesternpts_0621_2017_sampled_no.gpx
# CaplesWesternpts_0621_2017_sampled_yes.gpx
# CAPLES_BIRDS_systematic_pts_with_veg_and_utms (2).xlsx
# CAPLES_BIRDS_systematic_pts_with_veg_and_utms.xlsx
# CAPLES_BIRDS_systematic_pts_with_veg_and_utms (2).xlsx
# CaplesWesternpts_0621_2017_sampled_no.kml
# Caples_Avian_PtList_UTM_SizeDensity_20201117.csv
# caples_pts_2018_MKC.xlsx (limited)
# caples_pts_2018_MKC.xlsx (more up to date)
# caples_pts_2018_MKC.xlsx (limited copy 2)
# Wildlife sampling points.xlsx
# 1070_Eldo.gpx & all similar single point GPX exports
# Wildlife sampling points.xlsx
# caples_pts_2018_MKC.xlsx
# ARU_points_20180615.xlsx
# NOTE ALL SAME (on ..46.3 .. 41 GRID) AS 'wild_points.shp' with addition of 1072 see below


# b. Shape files from veg_plot_salo_comp.R
# read in wildlife points
#[ ] update wild_points to include plant plots, see e-mail 2022-06-28 Sarah Jacobes


FinalCaplesMonitoringPlots2022 <- sf::read_sf(here("spatial", "input", "shapefiles", "Monitoring2022", "FinalCaplesMonitoringPlots2022.shp")) %>% # crs not included
  st_set_crs(26910) %>% #were collected using NAD83 coordinate (26910) coordinate reference system
  st_transform(crs(study_area)) #, aoi = study_area$geometry) #defined in other script
#[ ] NOTE PROBLEM NO TRANSFORMATION IS HAPPENING BETWEEN NAD83 AND WGS84 NEED TO FIGURE OUT WHAT'S UP


#run veg_plot_salo_comp.R lines through 70

# 2. name each file with suffix to trace history
# 3. write-down each coordinate reference system CRS
# 4. import to this script
# 5. transform all to same WGS 84 crs (or all to NAD83)
# 6. map each layer by name
# 7. FIGURE Out WTF

# 8. Then correct the veg_plot location
# 9. Redo the QA/QC (accuracy of vegetation poin):
# 9a. @ all vegetation plots x = ground y = satellite: relationship, #params, r^2
# 9b. @ all bird points, view satellite by values ... graphic? 
# 9c. @ all bird points [ ] -> modify (interpret) bird points satellite measures using 9a model if necessary
# 10. Redo Precision 
# 10a. precision of cover between 2018 & 2019 -- all points on study area (bird/veg)
# 10b. precision of height between 2018 & 2019 -- all points on study area (bird/veg)
# 10c. model cover / height change 2018/19 vs 2020 
# 10d. make cover / height decision (average?, if not aggregation function, what to output for Mary)


library(leaflet)
library(readr)
library(tidyverse)

# (from read_points_output_data.R)
# [ ] READ in as close to original as possible

new_wild_points %>% 
  dplyr::select(plotID_av = point_d, plotID_veg = plotID_veg, geometry) -> nwp2

crs(nwp2)

vwf_11.37 <- nwp2  %>% sf::st_buffer(11.37,endCapStyle = "ROUND") 

# 
# 1. Vegetation plots
# 
# Summarize key variables (Sarah)
# using output 'Caples_PlotData_20220225'
# so must source("Caples_Tidy_20220221.R" to make it work)

Caples_PlotData_20220225 <- st_as_sf(Caples_PlotData_20220225, coords = c("Easting", "Northing"), crs = NAD83crs)

#change crs
Caples_PlotData_20220225 <- st_transform(Caples_PlotData_20220225, crs(study_area)) # transformed crs

#select out points for comparisons
Caples_PlotData_20220225 %>%
  filter(SamplingTime == 0) %>%
  dplyr::select(PlotID, geometry, Caples_Severity_Class) %>% 
  distinct() -> veg_plot_points_for_comparison

veg_plot_points_for_comparison %>%
  st_join(nwp2) %>% View()  #way off, no join

View(veg_plot_points_for_comparison)

View(nwp2)

# vwf_11.37 <- veg_plot_points_for_comparison  %>% sf::st_buffer(11.37,endCapStyle = "ROUND") 
#[ ] to do ... redo the chain of custody analysis on the points as data frames and then plot the various options for 
# UTMS on leaflet maps.

# 2. Satellite data
# using processed data after running read_adjust_combine_rasters.R
# add in 11.3 m^2 (extract and summarize) for some variables ...
canopy_fuel_nbr_dem_RAVG_LIDAR %>%
  exactextractr::exact_extract(., vwf_11.37, c("mean", "median", "min", "max", "count")) -> canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m

canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m %>% dim()

#2.5 mapping

point_data <- veg_sat_comp %>% dplyr::select(Avian_Poin,CSE_ID,geometry)


mapa <- leaflet(point_data) %>%
addProviderTiles("Esri.WorldImagery") %>% #addTiles(mapa)
# mapa <-   addPolygons(mapa, color = "#444444", weight = 1, smoothFactor = 0.5,
#                       opacity = 1.0, fillOpacity = 0.5,
#                       fillColor = ~colorQuantile("YlOrRd", Shape_Area)(Shape_Area),
#                       highlightOptions = highlightOptions(color = "red", weight = 2,
#                                                           bringToFront = TRUE))  
addMarkers(lng=sf::st_coordinates(point_data$geometry)[,1], lat=sf::st_coordinates(point_data$geometry)[,2], popup=paste(point_data$plotID_av, point_data$plotID_veg))

mapa

leaflet(point_data) %>% 
  fitBounds(lng1 = min(point_data$lon) - 0.11, 
            lat1 = min(point_data$lat) - 0.11,
            lng2 = max(point_data$lon) + 0.11, 
            lat2 = max(point_data$lat) + 0.11) %>% 
  addCircleMarkers(~lon, ~lat,
                   radius = 3,
                   opacity = 100,
                   color = "white", 
                   label = as.character(point_data$Site), 
                   labelOptions = labelOptions(noHide = TRUE, 
                                               textOnly = TRUE, 
                                               style = list(color = "white"),
                                               direction = "auto",
                                               offset = c(0, -10))) %>% 
  addScaleBar(options = list(imperial = FALSE)) %>%
  addProviderTiles(providers$Esri.WorldImagery)


# 3. Comparisons
Caples_PlotData_20220225 %>%
  filter(SamplingTime == 0) -> CPD_time0

dim(CPD_time0)[[1]] == dim(canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m)[[1]]

#[ ] need a validation step here to make sure it is lined up (too bad the extractr doesn't give rownames for objects!)

cbind(CPD_time0,canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m) -> veg_sat_comp

#adjust levels for severity class
# dput(unique(veg_sat_comp$Caples_Severity_Class))
veg_sat_comp$Caples_Severity_Class <- ordered(veg_sat_comp$Caples_Severity_Class, levels = c("High", "Mod", "Low", "Unburned", NA))

#analysis of tree height
#histograms to see if they make sense
veg_sat_comp$Total_Tree_Ht %>% hist()
veg_sat_comp$mean.CaplesCanopyHeight2018 %>% hist()

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
  ggplot(.) +
  geom_point(mapping = aes(x = MaxTreeHeight, y = max.CaplesCanopyHeight2018)) +
  geom_abline(mapping = aes(intercept = 0, slope = 1)) + 
  geom_abline(slope = slope, color="red", 
              linetype="dashed", size=1.5) + 
  geom_text(aes(x = MaxTreeHeight, y = max.CaplesCanopyHeight2018, label = PlotID, hjust = -.25, vjust = -.35)) +
  stat_cor(aes(x = MaxTreeHeight, y = max.CaplesCanopyHeight2018, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 0, label.y = 35) +
  #facet_wrap(~Caples_Severity_Class) +
  #geom_abline(intercept = intercept2, slope = slope2, color="blue", 
  #          linetype="dotted", size=1.5) + 
  #  facet_wrap(~Caples_Severity_Class) +
  labs(y="max Canopy Height before fire (2018 satellite)", x = "Max Tree Height before fire (ground)", title="Correlation between on-the-ground and satellite measurements\n height for vegetation plot (~ 400 meters squared, avian points)") 


summary(reg)  
summary(reg2)  

anova(reg,reg2)
