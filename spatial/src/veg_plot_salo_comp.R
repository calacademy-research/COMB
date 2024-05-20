#read points summarize Veg plot and SALO data for 11.3 m^2 comparison
# veg_plot_salo_comp.R
# 0. Getting started
# 
# [] build new script to fix shape file for ALL points
# FIX POINT LOCATIONS TO BE BECKY'S LATEST vegetation_points
# install.packages("rgdal")
library(rgdal)
library(sf)
library(here)
library(tidyverse)
library(raster)
library(ggpubr)

#1. Read in wildlife points
#[ ] update wild_points to include plant plots, see e-mail 2022-06-28 Sarah Jacobs
#Import
wild_points <- sf::read_sf(here("spatial", "input", "shapefiles", "WildlifePoints.shp")) %>% # crs not included
  st_set_crs(26910) %>% #were collected using NAD83 coordinate (26910) coordinate reference system
  # #save NAD83crs for later 
  # st_crs(wild_points) -> NAD83crs
  st_transform(crs(study_area)) #, aoi = study_area$geometry) #defined in other script
#[X] NOTE PROBLEM NO TRANSFORMATION IS HAPPENING BETWEEN NAD83 AND WGS84 NEED TO FIGURE OUT WHAT'S UP (working now 2023-11-15)

#FIX NAMES
names(wild_points) <- c(
  "point_d", "Cpls_Wt", "VEG_CSE", "AVIAN_S", "WHR_TSD", "SZ_DNS2",
  "Treatment", "geometry"
)


# are points in the fire boundary area?
wild_points$inside_fire_boundary <- as.vector(st_intersects(fire_boundary, wild_points, sparse = FALSE))

#fix with "plotID_UTM.csv" that is now in the input directory (so running "read_points_output_data.R" creates it)
plotID_UTM <- read_csv(here("spatial", "input", "shapefiles", "plotID_UTM.csv"))
  
# plotID_UTM$plotID_av are the bird_points
# plotID_UTM$plotID_veg are the veg_points

#building a bigger table
new_wild_points <- left_join(plotID_UTM, wild_points, by = c("plotID_av" = "point_d"), keep = TRUE)

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
  st_set_crs(crs(study_area)) -> new_wild_points # transformed crs

#could do for all vegetation points (hold off until edited [ ])

crs(new_wild_points)

#write out fixed shapefile
sf::st_write(new_wild_points, dsn = here("spatial", "input", "shapefiles", "WildlifePointsFixed.shp"), append = FALSE)

#1.5 also do the same for FinalCaplesMonitoringPlots2022
FinalCaplesMonitoringPlots2022 <- sf::read_sf(here("spatial", "input", "shapefiles", "Monitoring2022", "FinalCaplesMonitoringPlots2022.shp")) %>% # crs not included
  st_set_crs(26910) %>% #were collected using NAD83 coordinate (26910) coordinate reference system
  st_transform(crs(study_area)) #, aoi = study_area$geometry) #defined in other script

# 2. Vegetation Data
# 
# Summarize key variables (Sarah)
# using output 'Caples_PlotData_20220225'
# so must 
# source(here("vegetation","src","Caples_Tidy_20220221.R"))
#to make it work)

Caples_PlotData_20220225 <- st_as_sf(Caples_PlotData_20220225, coords = c("Easting", "Northing"), crs = NAD83crs)

#change crs
Caples_PlotData_20220225 <- st_transform(Caples_PlotData_20220225, crs(study_area)) # transformed crs

#
Caples_PlotData_20220225$geometry

#select out points for comparisons
Caples_PlotData_20220225 %>%
  filter(SamplingTime == 0) %>%
  distinct() -> veg_plot_points_for_comparison

# 3. Merge new_wild_points and vegetation data

#drop geometry from both:
FinalCaplesMonitoringPlots2022_df <- cbind(st_drop_geometry(FinalCaplesMonitoringPlots2022), st_coordinates(FinalCaplesMonitoringPlots2022))
new_wild_points_df <- cbind(st_drop_geometry(vegetation_points), st_coordinates(vegetation_points))
veg_plot_points_for_comparison_df <- cbind(st_drop_geometry(veg_plot_points_for_comparison), st_coordinates(veg_plot_points_for_comparison))
#[ ]broken here 2022-10-14 but not important // [ ] ACTUALLY is important 2023-11-15!!!
#HERE
# see mismatch between UTM_N/UTM_E, X.x/Y.x,X.y/Y.y. the latter point pairs are off, figure it out!
#merge data based on shared key(s)
FinalCaplesMonitoringPlots2022_df %>%
  left_join(new_wild_points_df, c("CSE_ID" = "PLOTID"), keep = TRUE) %>% View() #used to join by "plotID_veg"
  mutate(CSE_RF_ID = ifelse(is.na(RedFir_ID.x), CSE_ID.x, RedFir_ID.x)) %>% 
  left_join(veg_plot_points_for_comparison_df, by = c("CSE_RF_ID" = "PlotID"),  keep = TRUE) %>% 
  dplyr::select(UTM_N_FCMP22 = UTM_N, UTM_E_FCMP22 = UTM_E, 
                FCMP22_Y.y = Y.x, FCMP22_X.x = X.x, 
                # UTM_N_nwp_Northing = Northing, UTM_N_nwp_Easting = Easting,
                nwp_Y.y = Y.y, nwp_X.y = X.y,
                vppfc_Y = Y, vppfc_X = X, CSE_ID.x:RedFir_ID.x, Cpls_Wt:inside_fire_boundary, PV_ID:WYMO) %>% #View() #yikes, they are all over the place!
  st_as_sf(., coords = c("UTM_E_FCMP22", "UTM_N_FCMP22")) -> fcmp22data #make first two into the geometry

st_crs(fcmp22data) <- 26910
fcmp22data <- st_transform(fcmp22data, crs(study_area))

#extract data from ^^
# all bunk above

#add new locs to fcmp22data as a hack [HERE 2022-10-14] [ ] problematic fcmp22data ... what are the right points?
fcmp22data %>%
  left_join(as.data.frame(.), as.data.frame(vegetation_points), by=c("CSE_ID"="CSE_ID")) %>% View()
  select()

vwf_11.37 <- fcmp22data  %>% sf::st_buffer(11.37,endCapStyle = "ROUND") 

# vwf_11.37 <- vegetation_points  %>% sf::st_buffer(11.37,endCapStyle = "ROUND") 

# 4. Comparisons
# a. Satellite data
# using processed data after running read_adjust_combine_rasters.R
# add in 11.3 m^2 (extract and summarize) for some variables ...
canopy_fuel_nbr_dem_RAVG_LIDAR %>%
  exactextractr::exact_extract(., vwf_11.37, c("mean", "median", "min", "max", "count")) -> canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m

canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m %>% dim()

dim(fcmp22data)[[1]] == dim(canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m)[[1]]

#[ ] need a validation step here to make sure it is lined up (too bad the extractr doesn't give rownames for objects!)

cbind(fcmp22data,canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m) -> veg_sat_comp

veg_sat_comp$geometry

veg_sat_comp$UTM_N <- sf::st_coordinates(veg_sat_comp$geometry)[,2]
veg_sat_comp$UTM_E <- sf::st_coordinates(veg_sat_comp$geometry)[,1]

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
  filter(Avian_Poin > 0) %>%
  # filter()
  ggplot(.) +
    geom_point(mapping = aes(x = MaxTreeHeight, y = max.CaplesCanopyHeight2018)) +
    geom_abline(mapping = aes(intercept = 0, slope = 1)) + #1:1 line
    geom_abline(slope = slope, color="red", #
              linetype="dashed", size=1.5) +
  geom_abline(slope = slope2, intercept = intercept2, color="blue", #
              linetype="dashed", size=1) +
  geom_text(aes(x = MaxTreeHeight, y = jitter(max.CaplesCanopyHeight2018), label = PlotID, hjust = -.25, vjust = -.35)) +
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
  #facet_wrap(~Caples_Severity_Class) +
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
  facet_wrap(~Caples_Severity_Class) +
  labs(y="mean Canopy Height after fire (2020 satellite)", 
       x = "mean Canopy Height before fire (2018 satellite)", 
       title="Correlation between years satellite measurements\n canopy height at each point (~1 pixel) before and after fire") -> cor1820cv

cor1820cv
