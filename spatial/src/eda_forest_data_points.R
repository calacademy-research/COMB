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
library(googledrive) # connect to google drive

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
hist(wide1havars$Perc_LTg28mHt_2018_1ha, xlim = c(0, .5), main = "", xlab = "2018", ylab = "LT > 28m")
hist(wide1havars$Perc_LTg28mHt_2020_1ha, xlim = c(0, .5), main = "", xlab = "2020", ylab = "")
hist(wide1havars$Perc_LTg25mHt_2018_1ha, xlim = c(0, .5), main = "", xlab = "2018", ylab = "LT > 25m")
hist(wide1havars$Perc_LTg25mHt_2020_1ha, xlim = c(0, .5), main = "", xlab = "2020", ylab = "")
hist(wide1havars$Perc_LTg22mHt_2018_1ha, xlim = c(0, .5), main = "", xlab = "2018", ylab = "LT > 22m")
hist(wide1havars$Perc_LTg22mHt_2020_1ha, xlim = c(0, .5), main = "", xlab = "2020", ylab = "")
title(main = c("Distribution of fraction of plot with tree height >= cutoffs \n (22, 25 & 28m) for 2018 vs 2020 (across 82 1ha avian plots)"), outer = T)

par() <- parold

# ggplot faceted histogram
tall_forest_variables %>%
  filter(grepl("LT(.*)Ht", var), Year != 2019, scale == "4ha") %>%
  ggplot(., aes(value, group = var)) +
  geom_histogram(aes(y = stat(density)), bins = 20) +
  # scale_y_continuous(labels = var) +
  facet_wrap(var ~ Year, ncol = 2)

# ggplot faceted histogram
tall_forest_variables %>%
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

rm(tfvinput)
