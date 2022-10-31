#### Sample raster data around the points ####
# eda_forest_data_points.R
# Created by: Durrell D. Kapan
# Created on: 18 April 2022
# modified from Caples_Spatial_EDA.md

#### Load Libraries ####

# load libraries
library(here)
library(tidyverse)
# library(purrr)
library(raster)
#library(googledrive) # connect to google drive

# load functions
source(here("comb_functions.R"))

#------------------------------------------------------------
# Run previous data scripts for rasters and points
#------------------------------------------------------------
# #check if file exists in /spatial/output/rasters && ../tables/
# run read_adjust_combine_rasters.R
# source(here("/spatial/src/read_adjust_combine_rasters.R"))
# -AND-
# run read_points_output_data.R
# source(here("/spatial/src/read_points_output_data.R"))
# **Import raster brick**

#------------------------------------------------------------
# Plots to check sanity of data
#------------------------------------------------------------

# base R histogram

# save pars
par() -> parold

par(mfrow = c(3, 2))
par(oma = c(2, 2, 3, 0))
par(mar = c(3, 4, 1, 1))
hist(wide4havars$Perc_LTg28mHt_2018_4ha, xlim = c(0, .5), main = "", xlab = "2018", ylab = "LT > 28m")
hist(wide4havars$Perc_LTg28mHt_2020_4ha, xlim = c(0, .5), main = "", xlab = "2020", ylab = "")
hist(wide4havars$Perc_LTg25mHt_2018_4ha, xlim = c(0, .5), main = "", xlab = "2018", ylab = "LT > 25m")
hist(wide4havars$Perc_LTg25mHt_2020_4ha, xlim = c(0, .5), main = "", xlab = "2020", ylab = "")
hist(wide4havars$Perc_LTg22mHt_2018_4ha, xlim = c(0, .5), main = "", xlab = "2018", ylab = "LT > 22m")
hist(wide4havars$Perc_LTg22mHt_2020_4ha, xlim = c(0, .5), main = "", xlab = "2020", ylab = "")
title(main = c("Distribution of fraction of plot with tree height >= cutoffs \n (22, 25 & 28m) for 2018 vs 2020 (across 82 4ha avian plots)"), outer = T)

# par() <- parold

# ggplot faceted histogram
tall_forest_variables_4ha %>%
  filter(grepl("LT(.*)Ht", var), Year != 2019, scale == "4ha") %>% #View()
  ggplot(., aes(value, group = var)) +
  geom_histogram(aes(y = stat(density)), bins = 20) +
  # scale_y_continuous(labels = var) +
  facet_wrap(var ~ Year, ncol = 2)

# ggplot faceted histogram
tall_forest_variables_4ha %>%
  filter(grepl("LT(.*)Ht|RAVGcbi4", var), Year != 2019, scale == "4ha") %>%
  pivot_wider(
    names_from = sum_fn:Year,
    # names_glue ="{var}_{.value}", #printf('{%s}_{%s}_{%s', sum_fn, scale, Year
    values_from = value
  ) %>%
  pivot_longer(
    cols = contains("Perc"),
    names_to = c("sum_fn", "var", "Year"),
    names_pattern = "(.*)_(.*)_(.*)",
    values_to = "value",
    values_drop_na = TRUE
  ) %>%
  filter(grepl("LTg22mHt", var)) %>%
  mutate(RAVG = cut(mean_RAVGcbi4_20182019, breaks = c(-1, 0, 1, 2, 3, 4), labels = c("0", "(0,1]", "(1,2]", "(2,3]", "(3,4]"))) -> tfv_input

# plot the input
ggplot(tfv_input, aes(value, group = var)) +
  geom_histogram(aes(y = stat(density)), bins = 15) +
  # scale_y_continuous(labels = var) +
  facet_grid(RAVG ~ Year) +
  ggtitle(c("Large Trees (height >= 22m) versus burn severity (RAVG) before and after fire")) +
  xlab("Fraction of 4ha with tree height >= 22m") +
  ylab("Number of plots in each RAVG level")

#plot each year against each other
tall_forest_variables_4ha %>%
  filter(grepl("LT(.*)Ht|RAVGcbi4", var), Year != 2019, scale == "4ha") %>%  # View()
  pivot_wider(
    names_from = sum_fn:Year,
    # names_glue ="{var}_{.value}", #printf('{%s}_{%s}_{%s', sum_fn, scale, Year)
    values_from = value
  ) %>% 
  mutate(RAVG = cut(mean_RAVGcbi4_20182019, breaks = c(-1, 0, 1, 3, 4), labels = c("unburned 0", "low (0,1]", "moderate (1,3]", "high (3,4]"))) -> ltplotinput #breaks = c(-1, 0, 1, 2, 3, 4), labels = c("0", "(0,1]", "(1,2]", "(2,3]", "(3,4]"))

ltplotinput$RAVG <- ordered(ltplotinput$RAVG, levels = c("high", "moderate", "low", "unburned"))

# plot comparing 2018 to 2020
ggplot(ltplotinput, aes(Perc_LTg22mHt_2018, Perc_LTg22mHt_2020, color = RAVG)) +
  geom_point() +
  xlim(0, .55) +
  ylim(0, .55) +
  coord_equal() +
  theme(aspect.ratio=1) +
  geom_abline(intercept = 0, slope = 1) +
  # scale_y_continuous(labels = var) +
  facet_grid(~RAVG) +
  ggtitle(c("Large Trees (height >= 22m) versus burn severity (RAVG) before and after fire")) +
  xlab("Fraction of 4ha with tree height >= 22m in 2018") +
  ylab("Fraction of 4ha with tree height >= 22m in 2020") 
  

rm(tfvinput)

#these delta's must be calculated on a pixel by pixel basis then look at the central tendency [ ]
wide4havars$deltacov <- (wide4havars$mean_CanopyCover_2020_4ha-wide4havars$mean_CanopyCover_2018_4ha)
wide4havars$deltaht <- (wide4havars$mean_CanopyHeight_2020_4ha-wide4havars$mean_CanopyHeight_2018_4ha)
wide4havars$deltabsht <- (wide4havars$mean_CanopyBaseHeight_2020_4ha-wide4havars$mean_CanopyBaseHeight_2018_4ha)
wide4havars$deltabd <- (wide4havars$mean_CanopyBulkDensity_2020_4ha-wide4havars$mean_CanopyBulkDensity_2018_4ha)
wide4havars$deltalc <- (wide4havars$mean_CanopyLayerCount_2020_4ha-wide4havars$mean_CanopyLayerCount_2018_4ha)
wide4havars$deltalf <- (wide4havars$mean_LadderFuelDensity_2020_4ha-wide4havars$mean_LadderFuelDensity_2018_4ha)
wide4havars$deltasf <- (wide4havars$mean_SurfaceFuels_2020_4ha-wide4havars$mean_SurfaceFuels_2018_4ha)

wide4havars %>%
  mutate(RAVG =cut(mean_RAVGcbi4_20182019, breaks = c(-1, 0, 1, 3, 4), labels = c("unburned 0", "low (0,1]", "moderate (1,3]", "high (3,4]"))) -> wide4havars

wide4havars$RAVG <- ordered(wide4havars$RAVG, levels = c("high", "moderate", "low", "unburned"))

wide4havars %>%
  ggplot(aes(x = avian_point, y = deltaht, color = RAVG))+
  geom_point() +
  geom_abline(intercept = 0, slope = 0) +
  scale_color_manual(values = c("high" = "black",
                                "moderate" = "red",
                                "low" = "orange",
                                "unburned" = "green")) +
  ggtitle("change in mean canopy height 2020-2018") +
  ylab("change in mean canopy height 4ha") +
  facet_wrap(~RAVG)

wide4havars %>%
  ggplot(aes(x = avian_point, y = deltacov, color = RAVG))+
  geom_point() +
  geom_abline(intercept = 0, slope = 0) +
  scale_color_manual(values = c("high" = "black",
                                "moderate" = "red",
                                "low" = "orange",
                                "unburned" = "green")) +
  ggtitle("change in mean canopy cover 2020-2018") +
  ylab("change in mean canopy cover 4ha") +
  facet_wrap(~RAVG)

wide4havars %>%
  ggplot(aes(x = avian_point, y = deltabsht, color = RAVG)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 0) +
  scale_color_manual(values = c("high" = "black",
                                "moderate" = "red",
                                "low" = "orange",
                                "unburned" = "green")) +
  ggtitle("change in mean base height 2020-2018") +
  ylab("change in mean base height 4ha") +
  facet_wrap(~RAVG)

wide4havars %>%
  ggplot(aes(x = avian_point, y = deltalf, color = RAVG)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 0) +
  scale_color_manual(values = c("high" = "black",
                                "moderate" = "red",
                                "low" = "orange",
                                "unburned" = "green")) +
  ggtitle("change in ladder fuels 2020-2018") +
  ylab("change in mean ladder fuels 4ha") +
  facet_wrap(~RAVG)

wide4havars %>%
  ggplot(aes(x = mean_CanopyHeight_2018_4ha, y = mean_CanopyHeight_2020_4ha, color = RAVG)) +
  geom_point() +
  stat_summary(fun.data=mean_cl_normal) +
  geom_smooth(method='lm', formula= y~x) +
  geom_abline(intercept = 0, slope = 1) +
  coord_fixed() +
  theme(aspect.ratio=1) +
  scale_color_manual(values = c("high" = "black",
                                "moderate" = "red",
                                "low" = "orange",
                                "unburned" = "green")) +
  ggtitle(c("Relationship between canopy variables vs \nburn severity (RAVG) before and after fire")) +
  xlab("Mean canopy height 4ha (2018)") +
  ylab("Mean canopy height 4ha (2020)") +
  facet_wrap(~RAVG) 

wide4havars %>%
  ggplot(aes(x = mean_CanopyCover_2018_4ha, y = mean_CanopyCover_2020_4ha, color = RAVG)) +
  geom_point() +
  stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x) +
  geom_abline(intercept = 0, slope = 1) +
  coord_fixed() +
  theme(aspect.ratio=1) +
  scale_color_manual(values = c("high" = "black",
                                "moderate" = "red",
                                "low" = "orange",
                                "unburned" = "green")) +
  ggtitle(c("Relationship between canopy variables vs \nburn severity (RAVG) before and after fire")) +
  xlab("Mean canopy cover 4ha (2018)") +
  ylab("Mean canopy cover 4ha (2020)") +
  facet_wrap(~RAVG) 

wide4havars %>%
  ggplot(aes(x = deltaht, y = deltacov, color = RAVG)) +
  geom_point() +
  stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x) +
  geom_abline(intercept = 0, slope = 1) +
  coord_fixed() +
  theme(aspect.ratio=1) +
  scale_color_manual(values = c("high" = "black",
                                "moderate" = "red",
                                "low" = "orange",
                                "unburned" = "green")) +
  ggtitle(c("Relationship between canopy variables vs \nburn severity (RAVG) before and after fire")) +
  xlab("Change in height (2020-2018)") +
  ylab("Change in cover (2020-2018)") +
  facet_wrap(~RAVG) 

wide4havars %>%
  ggplot(aes(x = deltasf, y = deltabsht, color = RAVG)) +
  geom_point() +
  stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x) +
  scale_color_manual(values = c("high" = "black",
                                "moderate" = "red",
                                "low" = "orange",
                                "unburned" = "green")) +
  ggtitle(c("Relationship between fuel and canopy base height \n vs burn severity (RAVG) before and after fire")) +
  xlab("Change in surface fuel (2020-2018)") +
  ylab("Change in canopy base height (2020-2018)") +
  facet_wrap(~RAVG)

wide4havars %>%
  ggplot(aes(x = deltalf, y = deltabsht, color = RAVG)) +
  geom_point() +
  stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x) +
  # geom_abline(intercept = 0, slope = 1) +
  coord_fixed() +
  theme(aspect.ratio=1) +
  scale_color_manual(values = c("high" = "black",
                                "moderate" = "red",
                                "low" = "orange",
                                "unburned" = "green")) +
  ggtitle(c("Relationship between fuel and canopy base height \n vs burn severity (RAVG) before and after fire")) +
  xlab("Change in ladder fuel (2020-2018)") +
  ylab("Change in canopy base height (2020-2018)") +
  facet_wrap(~RAVG) 

wide4havars %>%
  ggplot(aes(x = deltasf, y = deltalf, color = RAVG)) +
  geom_point() +
  # geom_abline(intercept = 0, slope = 1) +
  stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x) +
  coord_fixed() +
  theme(aspect.ratio=1) +
  scale_color_manual(values = c("high" = "black",
                                "moderate" = "red",
                                "low" = "orange",
                                "unburned" = "green")) +
  ggtitle(c("Relationship between fuel variables vs \nburn severity (RAVG) before and after fire")) +
  ylab("Change in ladder fuel (2020-2018)") +
  xlab("Change in surface fuel (2020-2018)") +
  facet_wrap(~RAVG) 

