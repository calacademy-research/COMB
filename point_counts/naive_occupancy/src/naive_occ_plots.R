#### Naive Occupancy Plots for Species ####
# naive_occ_plots-best_DK_update.R
# Created by: Natalie Beckman-Smith 
# Based on earlier work by Mary Clapp & Durrell Kapan
# Created on: 26 July 2021
# Last updated on: 26 July 2021
# modifed by Durrell Kapan to be worked into the main
# codebase 2021-09-30 & 2022-03-20 (added function)

#### Load Libraries ####
library(tidyverse)
library(here)
library(forcats)

#### Load Data ####
# Must first run naive_occupancy_17-21.R to get diff df
#source(here("point_counts","naive_occupancy","src","naive_occupancy_17-21.R")) #if necessary to create the diff file
# [ ] fix "diff" to be read from output if doesn't exist & you don't want to rerun the ^ script to make iit
diff_filename <- "diff.csv"
#Read the google data
PointC <- fread(here("point_counts","naive_occupancy","output", diff_filename))

#define function to plot
diffplot <- function(datadiff, dw, ww, ttl, ttl_size, txt_size, txt_angle){
  ggplot(data = datadiff, 
       aes(x = reorder(birdName, -naive_occ), y = naive_occ, fill = year,
           ymin = as.numeric(min1719), 
           ymax = as.numeric(max1719))) +
  geom_col(position = position_dodge(dw), 
           width = ww) + 
  geom_errorbar(data = datadiff, stat = "identity", width = .5, 
                colour = "black", alpha = 0.8, size = .5, 
                position = position_dodge(dw)) +
  #ylim(0, max(max1719, na.rm=TRUE)) +
  scale_fill_brewer(palette = "Dark2", name = "", 
                    labels = c("Pre-fire\n 2017-19", "Post-fire\n 2020", "Post-fire\n 2021")) +
  labs(title = ttl, x = "Common Name", y = "# of survey points detected") +
  theme_bw() +
  theme(plot.title = element_text(size = ttl_size), axis.title = element_text(size = txt_size), 
        axis.text.x = element_text(angle = txt_angle, hjust = 1, size = txt_size), 
        plot.margin = margin(2, 2, 2, 50, unit = "pt"), 
        legend.position = c(.775, .875), legend.direction = "horizontal")
}

# if desired set the filter as a percentage 
# of the total 82 points (e.g. >8 forces at 
# least 10% of the points to be occupied
# in one of average pre vs post-fire)
# no_occ_filter <- 8

#diff %>% 
#  filter(pre >= no_occ_filter | post >= no_occ_filter) %>% ...
##more common species either before or after

#### Now calculate different scenarios
#### Species that increased in 2020 ####
diff_temp_up20 <- diff %>%
  select("birdName","birdCode_fk","1719mean"="mean1719","2020","2021","min1719","max1719","diff") %>%
  filter(`2020` > max1719) %>% 
  pivot_longer(cols = "1719mean":"2021", 
               names_to = "year",
               values_to = "naive_occ")

#add back NAs
diff_temp_up20$min1719[diff_temp_up20$year == "2020"] <- NA
diff_temp_up20$min1719[diff_temp_up20$year == "2021"] <- NA
diff_temp_up20$max1719[diff_temp_up20$year == "2020"] <- NA
diff_temp_up20$max1719[diff_temp_up20$year == "2021"] <- NA 

ttl <- paste0("Species that increased post-fire")

up20 <- 
  diffplot(datadiff = diff_temp_up20, dw = .8, ww = .8, ttl = ttl, ttl_size = 12, txt_size = 10, txt_angle = 42)

#### Species that decreased in 2020 ####
diff_temp_down20 <- diff %>%
  select("birdName","birdCode_fk","1719mean"="mean1719","2020","2021","min1719","max1719","diff") %>%
  filter(`2020` < min1719) %>%   
  pivot_longer(cols = "1719mean":"2021", 
               names_to = "year",
               values_to = "naive_occ")

#add back NAs
diff_temp_down20$min1719[diff_temp_down20$year == "2020"] <- NA
diff_temp_down20$min1719[diff_temp_down20$year == "2021"] <- NA
diff_temp_down20$max1719[diff_temp_down20$year == "2020"] <- NA
diff_temp_down20$max1719[diff_temp_down20$year == "2021"] <- NA    

ttl <- paste0("Species that decreased post-fire")

down20 <- 
  diffplot(datadiff = diff_temp_down20, dw = .8, ww = .8, ttl = ttl, ttl_size = 12, txt_size = 10, txt_angle = 42)

#### Species that stayed the same in 2020 ####
diff_temp_same20 <- diff %>%
  select("birdName","birdCode_fk","1719mean"="mean1719","2020","2021","min1719","max1719","diff") %>%
  filter(`2020` <= max1719,
         `2020` >= min1719) %>% 
  pivot_longer(cols = "1719mean":"2021", 
               names_to = "year",
               values_to = "naive_occ")

#add back NAs
diff_temp_same20$min1719[diff_temp_same20$year == "2020"] <- NA
diff_temp_same20$min1719[diff_temp_same20$year == "2021"] <- NA
diff_temp_same20$max1719[diff_temp_same20$year == "2020"] <- NA
diff_temp_same20$max1719[diff_temp_same20$year == "2021"] <- NA  

ttl <- paste0("Species that remained the same post-fire")

same20 <- 
  diffplot(datadiff = diff_temp_same20, dw = .8, ww = .8, ttl = ttl, ttl_size = 12, txt_size = 10, txt_angle = 42)
#[ ]^ not a good plot, some are 0=0!

#### Species that decreased in 2021 ####
diff_temp_down21 <- diff %>%
  select("birdName","birdCode_fk","1719mean"="mean1719","2020","2021","min1719","max1719","diff") %>%
  filter(`2021` < min1719) %>%   
  pivot_longer(cols = "1719mean":"2021", 
               names_to = "year",
               values_to = "naive_occ")

#add back NAs
diff_temp_down21$min1719[diff_temp_down21$year == "2020"] <- NA
diff_temp_down21$min1719[diff_temp_down21$year == "2021"] <- NA
diff_temp_down21$max1719[diff_temp_down21$year == "2020"] <- NA
diff_temp_down21$max1719[diff_temp_down21$year == "2021"] <- NA    

ttl <- paste0("Species with distributions that decreased in 2021")

down21 <- 
  diffplot(datadiff = diff_temp_down21, dw = .8, ww = .8, ttl = ttl, ttl_size = 12, txt_size = 10, txt_angle = 42)

#### Species that increased in both 2020 and 2021 ####
diff_temp_up20_up21 <- diff %>%
  select("birdName","birdCode_fk","1719mean"="mean1719","2020","2021","min1719","max1719","diff") %>%
  filter(`2020` > max1719,
         `2021` > max1719) %>% 
  pivot_longer(cols = "1719mean":"2021", 
               names_to = "year",
               values_to = "naive_occ")

#add back NAs
diff_temp_up20_up21$min1719[diff_temp_up20_up21$year == "2020"] <- NA
diff_temp_up20_up21$min1719[diff_temp_up20_up21$year == "2021"] <- NA
diff_temp_up20_up21$max1719[diff_temp_up20_up21$year == "2020"] <- NA
diff_temp_up20_up21$max1719[diff_temp_up20_up21$year == "2021"] <- NA  

ttl <- paste0("Post-fire species occupancy > pre-fire levels in 2020 & 2021")

up20up21 <- 
  diffplot(datadiff = diff_temp_up20_up21, dw = .8, ww = .8, ttl = ttl, ttl_size = 12, txt_size = 10, txt_angle = 42)

#### Species that increased in 2020 and remained the same in 2021 ####
diff_temp_up20_same21 <- diff %>%
  select("birdName","birdCode_fk","1719mean"="mean1719","2020","2021","min1719","max1719","diff") %>%
  filter(`2020` > max1719,
         `2021` <= max1719,
         `2021` >= min1719) %>% 
  pivot_longer(cols = "1719mean":"2021", 
               names_to = "year",
               values_to = "naive_occ")

#add back NAs
diff_temp_up20_same21$min1719[diff_temp_up20_same21$year == "2020"] <- NA
diff_temp_up20_same21$min1719[diff_temp_up20_same21$year == "2021"] <- NA
diff_temp_up20_same21$max1719[diff_temp_up20_same21$year == "2020"] <- NA
diff_temp_up20_same21$max1719[diff_temp_up20_same21$year == "2021"] <- NA  

ttl <- paste0("Post-fire species occupancy > pre-fire levels (2020) \n returning to pre-fire levels in 2021.")

up20same21 <- 
  diffplot(datadiff = diff_temp_up20_same21, dw = .8, ww = .8, ttl = ttl, ttl_size = 12, txt_size = 10, txt_angle = 42)

#### Species that remained the same in 2020 but increased in 2021 ####
diff_temp_same20_up21 <- diff %>%
  select("birdName","birdCode_fk","1719mean"="mean1719","2020","2021","min1719","max1719","diff") %>%
  filter(`2020` <= max1719,
         `2020` >= min1719,
         `2021` > max1719,) %>% 
  pivot_longer(cols = "1719mean":"2021", 
               names_to = "year",
               values_to = "naive_occ")

#add back NAs
diff_temp_same20_up21$min1719[diff_temp_same20_up21$year == "2020"] <- NA
diff_temp_same20_up21$min1719[diff_temp_same20_up21$year == "2021"] <- NA
diff_temp_same20_up21$max1719[diff_temp_same20_up21$year == "2020"] <- NA
diff_temp_same20_up21$max1719[diff_temp_same20_up21$year == "2021"] <- NA  

ttl <- paste0("Post-fire species occupancy within pre-fire levels (2020) \n increasing above pre-fire levels in 2021.")

same20up21 <- 
  diffplot(datadiff = diff_temp_same20_up21, dw = .8, ww = .8, ttl = ttl, ttl_size = 12, txt_size = 10, txt_angle = 42)

#### Species that remained the same in both 2020 and 2021 ####
diff_temp_same20_same21 <- diff %>%
  select("birdName","birdCode_fk","1719mean"="mean1719","2020","2021","min1719","max1719","diff") %>%
  filter(`2020` <= max1719,
         `2020` >= min1719,
         `2021` <= max1719,
         `2021` >= min1719,) %>% 
  pivot_longer(cols = "1719mean":"2021", 
               names_to = "year",
               values_to = "naive_occ")

#add back NAs
diff_temp_same20_same21$min1719[diff_temp_same20_same21$year == "2020"] <- NA
diff_temp_same20_same21$min1719[diff_temp_same20_same21$year == "2021"] <- NA
diff_temp_same20_same21$max1719[diff_temp_same20_same21$year == "2020"] <- NA
diff_temp_same20_same21$max1719[diff_temp_same20_same21$year == "2021"] <- NA  

ttl <- paste0("Post-fire species occupancy within pre-fire range \nin both 2020 & 2021")

same20same21 <- 
  diffplot(datadiff = diff_temp_same20_same21, dw = .8, ww = .8, ttl = ttl, ttl_size = 12, txt_size = 10, txt_angle = 42)

#### Species that decreased in 2020 but returned to pre-fire levels in 2021 or decreaed ####
diff_temp_down20_same21 <- diff %>%
  select("birdName","birdCode_fk","1719mean"="mean1719","2020","2021","min1719","max1719","diff") %>%
  filter(`2020` < min1719,
         `2021` <= max1719,
         `2021` >= min1719,) %>% 
  pivot_longer(cols = "1719mean":"2021", 
               names_to = "year",
               values_to = "naive_occ")

#add back NAs
diff_temp_down20_same21$min1719[diff_temp_down20_same21$year == "2020"] <- NA
diff_temp_down20_same21$min1719[diff_temp_down20_same21$year == "2021"] <- NA
diff_temp_down20_same21$max1719[diff_temp_down20_same21$year == "2020"] <- NA
diff_temp_down20_same21$max1719[diff_temp_down20_same21$year == "2021"] <- NA  

ttl <- paste0("Species that decreased post-fire in 2020 and returned to pre-fire range in 2021")

down20same21 <- 
  diffplot(datadiff = diff_temp_down20_same21, dw = .8, ww = .8, ttl = ttl, ttl_size = 12, txt_size = 10, txt_angle = 42)
  
#### Species that decreased in 2020 and/or 2021 ####
  diff_temp_down20_or_down21 <- diff %>%
    select("birdName","birdCode_fk","1719mean"="mean1719","2020","2021","min1719","max1719","diff") %>%
    filter(`2020` < min1719 | `2021` < min1719) %>% 
    pivot_longer(cols = "1719mean":"2021", 
                 names_to = "year",
                 values_to = "naive_occ")
  
#add back NAs
diff_temp_down20_or_down21$min1719[diff_temp_down20_or_down21$year == "2020"] <- NA
diff_temp_down20_or_down21$min1719[diff_temp_down20_or_down21$year == "2021"] <- NA
diff_temp_down20_or_down21$max1719[diff_temp_down20_or_down21$year == "2020"] <- NA
diff_temp_down20_or_down21$max1719[diff_temp_down20_or_down21$year == "2021"] <- NA  

ttl <- c("Post-fire species occupancy lower than pre-fire range \nin  2020 and / or 2021")
  
down20ordown21 <- 
  diffplot(datadiff = diff_temp_down20_or_down21, dw = .8, ww = .8, ttl = ttl, ttl_size = 12, txt_size = 10, txt_angle = 42)

#### Printing Plots ###
print(up20)
print(same20)
print(down20)

print(up20up21)
print(up20same21)

print(same20up21)
print(same20same21)

print(down20same21)
print(down21)

#decide to save those that are final

#if we removed extra busy labels, could be combined via
#library(patchwork)

#up20 / same20 / down20
#up20up21 + up20same21
#same20up21 + same20same21
#down20same21 + down21

#scale_fill_manual(name = "",
#labels = c("Pre-fire (2017-2019)","2020","2021"),
#values = c("#57C4AD","#EDA247","#EDA247")) +
