#lidar_SALO.R

#canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m depends on 

plot(canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m$max.Dominant_Height,canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m$max.Mean_Height)
plot(canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m$max.Dominant_Height,canopy_fuel_nbr_dem_RAVG_LIDAR_extract_11.37m$max.CaplesCanopyHeight2018)

#linear model of one versus the other 
reg<-lm(formula = max.Dominant_Height ~ MaxTreeHeight + 0,
        data=veg_sat_comp)                      

#get intercept and slope value
coeff<-coefficients(reg)          
slope<-coeff[1]
# slope<- coeff[2]
# intercept<-coeff[1]

reg2<-lm(formula = max.Dominant_Height ~ MaxTreeHeight,
         data=veg_sat_comp)                      

#get intercept and slope value
coeff<-coefficients(reg2)          
intercept2<-coeff[1]
slope2<-coeff[2]


veg_sat_comp %>%
  filter(!is.na(MaxTreeHeight)) %>%
  filter(CSE_ID != 31) %>%
  ggplot(.) +
  geom_point(mapping = aes(x = MaxTreeHeight, y = max.Dominant_Height)) +
  xlim(0, 55) +
  ylim(0,55) +
  geom_abline(mapping = aes(intercept = 0, slope = 1)) + 
  #geom_abline(mapping = aes(intercept = 0, slope = slope, color="red", 
  #            linetype="dashed", size=1.5)) + 
  geom_text(aes(x = MaxTreeHeight, y = max.Dominant_Height, label = Avian_Poin, hjust = -.25, vjust = -.35, cex = 1.0)) +
  stat_cor(aes(x = MaxTreeHeight, y = max.Dominant_Height, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 0, label.y = 50) +
  # facet_wrap(~cut(mean.ca3858612053820210815_20201011_20211016_rdnbr_cbi4,breaks = c(-1, 0.5, 2, 3, 4), labels=c(0,1,2,3))) +
  facet_wrap(~cut(Total_Tree_Cover, 3, labels=c("low","medium","high"))) +
  #geom_abline(intercept = intercept2, slope = slope2, color="blue", 
  #          linetype="dotted", size=1.5) + 
  #  facet_wrap(~Caples_Severity_Class) +
  coord_fixed() +
  labs(y="max Canopy Height before fire (2018(?) LIDAR)", x = "Max Tree Height before fire (ground)", title="Correlation between on-the-ground and satellite measurements\n height for vegetation plot (~ 400 meters squared, avian points)") 
