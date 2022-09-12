library(leaflet)
library(readr)
library(tidyverse)

# (from read_points_output_data.R)

new_wild_points %>% 
  dplyr::select(plotID_av = Avian_Poin, plotID_veg = CSE_ID, geometry = geometry.x) -> nwp2

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

veg_plot_points_for_comparison %>% as_tibble() 

nwp2 %>% as_tibble() 

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

point_data <- st_transform(nwp2, 4326)


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
  labs(y="max Canopy Height before fire (2018 satellite)", x = "Max Tree Height before fire (ground)", title="Correlation between on-the-ground and satellite measurements\n height for vegetation plot (~ 400 meters squared)") 


summary(reg)  
summary(reg2)  

anova(reg,reg2)
