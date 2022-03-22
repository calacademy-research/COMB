# naive_occupancy_17-21.R based on previous work:
# firesymp_figs_mean_range_mkc_ddk.R
# code to generate figures for exploratory data analysis pre- 
# and post-fire species effects
# written by M Clapp on 4/30/2021
# last updated 4/30/2021 by M Clapp
# modifications by DDK 5/1/2021 5/2/2021
# mods by MKC 5/3/2021 - added table of burn severity, 
# prelim geographical distributions by spp. [UNDER CONSTRUCTION]
# mods by DDK 5/4/2021 - noting that we can use the min and max to categorize the data
# relies on read_PC.R to generate object "PointC" (start there)
# working copy 3/20/2022
#head(PCTClong)
#what was the distribution of effort at points across the years and 
#are there any data entry errors

# [X] update this step to create naive occ object
# [ ] link to output file of readPC.R 
# [X] load key libraries here 
library(tidyverse)
library(here)
library(lubridate)
library(data.table)
library(fs)

#source(here("point_counts","data_ingest","src","readPC.R")) #if necessary to create the PC output file
#/point_counts/data_ingest/output/PointC_2022-03-19.csv or more up-to-date version

# Creating the folder within inputs that contains the symlinks that point to the output directory from "/data_ingest/src/readPC.R"
if (dir.exists(here("point_counts/naive_occupancy/input/")) == F) {
  dir.create(here("point_counts/naive_occupancy/input/"))
}

#Creating a symlink from the output of /point_count/data_ingest/ to the input of /point_count/naive_occupancy/
PointC_filename <- "PointC_2022-03-21.csv"
link_create(here("point_counts","data_ingest","output", PointC_filename), here("point_counts","naive_occupancy","input", PointC_filename))
#Read the google data
PointC <- fread(here("point_counts","naive_occupancy","input", PointC_filename))

#make a point_count summary dataframe
point_counts <- PointC %>% 
  ungroup() %>%
  mutate(Yr = year(DateTime)) %>%
  select(point_ID_fk, Yr, pointCount_ID_fk) %>% 
  distinct() %>% 
  arrange(point_ID_fk, Yr, pointCount_ID_fk) %>% 
  group_by(point_ID_fk, Yr) %>% 
  summarise(ncount= length(pointCount_ID_fk))

#looks like data entry errors, removed times with ":" in them in FMPRO [ ] need to adjust DB input type to string only ()
PointC %>% ungroup() %>% mutate(Yr = year(DateTime)) %>% filter(is.na(Yr)) %>% select(pointCount_ID_fk, point_ID_fk, DateTime, observer_fk) %>% distinct() -> tofix

#[ ] tofix has three issues:

#first issue
#sum(PointC$pointCount_ID_fk==0)
#[1] 351
PointC %>%
  filter(pointCount_ID_fk > 0) -> PointC
#[x] ^removed one empty record

PointC %>%
  filter(pointCount_ID_fk == 164) %>% View()
#PointC %>%
#  filter(pointCount_ID_fk == 164) %>% View()
#[ ] must fix this for NIKR (duplicate PC ID)
#pointCount_ID_fk, point_ID_fx, DateTime, observer_fk
#164 408 NA NIKR

PointC %>%
  filter(pointCount_ID_fk == 44) %>% View()
#PointC %>%
#  filter(pointCount_ID_fk == 44) %>% View()
#[ ] must fix this for JPD (duplicate PC ID)
#pointCount_ID_fk, point_ID_fx, DateTime, observer_fk
#44 582 NA JPD

PC <- PointC %>% 
  mutate(Yr = year(DateTime)) %>% 
  #filter by <100m, not flying over, not NA -- adding back flyover for overview
  filter(Yr>=2017, dist_detect != "100m+", birdCode_fk != "UNKN", birdCode_fk != "", observer_fk != "MASC", observer_fk != "IASH") %>% #dist_detect != "flyover",
  group_by(point_ID_fk, DateTime, Yr, observer_fk, birdCode_fk) %>%
  summarise(abun=n())

PCwide <- PC %>% 
  spread(key=birdCode_fk, value = abun, fill = 0)

PClong <- PCwide %>%
  pivot_longer(cols = AMDI:WIWA, names_to = "birdCode_fk", values_to = "count") %>% #WTSW
  group_by(Yr, point_ID_fk, birdCode_fk) %>%
  summarise(max.ct = as.integer(max(count))) %>%
  mutate(naive.occ = 0)

PClong$naive.occ[PClong$max.ct>0]<-1

#this is a list of species detected by point counts only
PCTClong %>% 
  select(birdName, birdCode_fk, eBird_6_code) %>% 
  distinct() -> Spp4LetterList

#update diff to include range for 2017-2019, and median or mean value ...
diff <- PClong %>% left_join(Spp4LetterList) %>%
  group_by(Yr, birdName, birdCode_fk) %>% 
  summarise(naive.occ = sum(naive.occ)) %>% 
  pivot_wider(names_from = Yr, values_from = naive.occ) %>% #view()
  mutate(mean1719 = mean(c(`2017`,`2018`,`2019`))) %>% 
  mutate(min1719 = min(c(`2017`,`2018`,`2019`))) %>%  
  mutate(max1719 = max(c(`2017`,`2018`,`2019`))) %>% 
  mutate(mean2021 = mean(c(`2020`,`2021`))) %>% 
  mutate(min2021 = min(c(`2020`,`2021`))) %>%  
  mutate(max2021 = max(c(`2020`,`2021`))) %>% 
  mutate(diff17192021 =`mean2021`-`mean1719`) %>% #View()
  mutate(absdiff17192021 =abs(`mean2021`-`mean1719`)) %>% #View()
  mutate(post_g_max17192021 = `max2021` > max1719) %>%
  mutate(post_l_min17192021 = `min2021` < min1719) %>%
  mutate(post_g_l_n = case_when(`max2021` > max1719 ~ "greater",
                                `max2021` <= max1719 && `min2021` >= min1719  ~ "in_range",
                                `min2021` < min1719 ~ "less")) %>%
  # mutate(pc1721 = mean1719/(mean1719+mean2021)) %>%
  mutate(pc1721 = (mean2021-mean1719)/(mean1719+mean2021)) %>%
  mutate(diff1719 = (`2020`-`mean1719`)) %>% 
  # mutate(mean1819 = mean(c(`2018`,`2019`))) %>% #optional 
  # mutate(min1819 = min(c(`2018`,`2019`))) %>%   #statistics 
  # mutate(max1819 = max(c(`2018`,`2019`))) %>%   #for 
  # mutate(diff1819 =`2020`-`mean1819`) %>%       #different 
  # mutate(diff2019 =`2020`-`2019`) %>%           #comparisons
  arrange(desc(abs(diff17192021)))

#set difference to look at
diff$diff <- diff$diff1719 #change if you compute optional statistic
diff$change <- diff$diff
diff$change <- ifelse(diff$diff > 0, "pos", "neg")
diff$change[diff$diff==0] <- "none"

#[ ] check source of these errors NAT found looking through naive occ best.R:
diff <- diff[-1*c(57),] #DOWO was a duplicate but not recognized by distinct(), unique() or duplicated()
diff <- diff %>% filter(birdCode_fk != "GCCR", #Removed non-species
                        birdCode_fk != "SPHY",
                        birdCode_fk != "0")

library(forcats)

#relevel (reorder) factor levels
lvls <- c("less", "in_range", "greater")
diff %>% 
  ungroup() %>%
  mutate(post_g_l_n = post_g_l_n %>%
           fct_relevel(lvls))

#these are the summary counts for the different cutoffs employed
#see post_g_l_n above
table(diff$post_g_l_n, rev(cut(diff$diff,c(-10,-3,3,8, 20))))
#[ ] make a clearly labeled table for kable

# change in NAIVE occupancy (# points detected) of *most* species pre and post fire. 
#For ease of reading, I filtered results to species detected 5+ times in either year, but we could remove that filter to look at all of them

diff %>% 
  filter(mean1719>=5 |`2020`>=5 |`2021`>=5) %>% 
ggplot() +
  geom_point(aes(x=reorder(birdCode_fk, -diff), y=diff, color=change)) + 
  geom_hline(aes(yintercept=0), size=1, color="blue", alpha=0.5) +
  geom_segment(aes(x=reorder(birdCode_fk, diff), y=0, xend=reorder(birdCode_fk, diff), yend=diff, color=change, size=0.1)) + #[ ] , size = 1
  ylim(-11,30) + 
  labs(x="4-letter avian species code", y="change in # of points detected (avg) pre- and post-fire") +
  scale_color_manual(values=c("#D95F02","#7570B3","#1B9E77")) +
  theme_bw() +
  theme(text = element_text(size=16), 
        axis.text.x = element_text(angle=65, hjust=1),
        legend.position = "none")
#these aren't filtered by whether or not they 
#are unusual given the range of variation 
#see suggested graphics suggested below

#set the 'pre' = previous comparison, choose from commented out possibilities
diff$pre <- diff$mean1719

#build the long dataset
diff %>% 
  mutate(post = `2020`) %>% # 
  pivot_longer(cols = c(pre, post), names_to = "pre_post", values_to = "naive.occ") %>%
  mutate(prop.pts = naive.occ/82, abs.diff = abs(diff)) %>% 
  mutate(min1719 = ifelse(pre_post == "post", `2020`, min1719)) %>% 
  mutate(max1719 = ifelse(pre_post == "post", `2020`, max1719)) -> difflong #View() #uncomment to inspect data table

#what is the frequency distribution of the presence on points
difflong %>%
  mutate(pre_post = fct_relevel(pre_post, rev)) %>%
  ggplot(aes(x=naive.occ))+
  geom_histogram(binwidth = 8) +
  facet_wrap(facets = vars(pre_post))

# [X] generate output file 
# Write diff to output folder
if (dir.exists(here("point_counts/naive_occupancy/output/")) == F) {
  dir.create(here("point_counts/naive_occupancy/output/"))
}

diff_filename <- paste0("diff", ".csv")
fwrite(diff, here("point_counts","naive_occupancy","output", diff_filename))

# remaining plots now superseded by 
# "naive_occ_plots-best_DK_update.R" 

# geographical patterns ---------------------------------------------------

## UNDER CONSTRUCTION MKC 5/3/2021

# [ ] TODO: work on "greying out" the previous years' points on the 2020 map; making it pretty, etc. may have to add separate layers/objects per year to do this
# delete tics ' ' to uncomment the below code which needs work
'

env <- read_csv("data/caples_env_1ha_20210315.csv") #work with Caples_spatial
utms <- read_csv("data/Wildlife_Sampling_Points.csv") #work with Caples_spatial
head(env)

prefire <- env %>% dplyr::select(point_d, mean_Cover_2018_1ha, mean_Height_2018_1ha, mean_NBR_2018_1ha, Perc_NtoE_2020_1ha, Binary_NtoE_2020_1ha, mean_Elevation_2020_1ha) %>% left_join(utms, by=c("point_d"="point_id")) 

geo <- PClong %>% dplyr::select(-naive.occ) %>% 
  filter(point_ID_fk !=397) %>% 
  pivot_wider(names_from = Yr, values_from = max.ct) %>%
  left_join(utms, by = c("point_ID_fk"="point_id")) 

osfl <- geo %>% filter(birdCode_fk=="OSFL") %>%
  mutate(everd = ifelse(`2017` == 0 & `2018`==0 & `2019`==0 & `2020`==0, "no", "yes")) %>% 
  filter(everd=="yes") %>% 
  pivot_longer(cols = 3:6, names_to = "Yr", values_to = "max.ct") %>% 
  filter(max.ct != 0)

osfl$max.ct <- as.integer(osfl$max.ct)

ggplot() +
  geom_point(data=geo, aes(x=X_U10N83, y=Y_U10N83)) +
  geom_point(data=osfl, aes(x=X_U10N83, y=Y_U10N83, color=Yr), size=3) +
  labs(title=paste0("Olive-sided Flycatcher detections")) +
  scale_color_viridis_d() +
  facet_wrap(~Yr)

bbwo <- geo %>% filter(birdCode_fk=="BBWO") %>%
  mutate(everd = ifelse(`2017` == 0 & `2018`==0 & `2019`==0 & `2020`==0, "no", "yes")) %>% 
  filter(everd=="yes") %>% 
  pivot_longer(cols = 3:6, names_to = "Yr", values_to = "max.ct") %>% 
  filter(max.ct != 0)

bbwo$max.ct <- as.integer(bbwo$max.ct)

ggplot() +
  geom_point(data=geo, aes(x=X_U10N83, y=Y_U10N83)) +
  geom_point(data=bbwo, aes(x=X_U10N83, y=Y_U10N83, color=Yr), size=3) +
  labs(title=paste0("Black-backed Woodpecker detections")) +
  scale_color_viridis_d() +
  facet_wrap(~Yr)

amro <- geo %>% filter(birdCode_fk=="AMRO") %>%
  mutate(everd = ifelse(`2017` == 0 & `2018`==0 & `2019`==0 & `2020`==0, "no", "yes")) %>% 
  filter(everd=="yes") %>% 
  pivot_longer(cols = 3:6, names_to = "Yr", values_to = "max.ct") %>% 
  filter(max.ct != 0)

amro$max.ct <- as.integer(amro$max.ct)

ggplot() +
  geom_point(data=geo, aes(x=X_U10N83, y=Y_U10N83)) +
  geom_point(data=amro, aes(x=X_U10N83, y=Y_U10N83, color=Yr, size=max.ct)) +
  labs(title=paste0("American Robin detections")) +
  scale_color_viridis_d() +
  facet_wrap(~Yr)



# burn severity across points ---------------------------------------------

yesCBI<-read_csv(paste0(here(),"/largeIO/CaplesWesternpts_0621_2017_sampled_yes_CBI.csv"))
yesCBI4<-read_csv(paste0(here(),"/largeIO/CaplesWesternpts_0621_2017_sampled_yes_CBI4.csv"))

burnSev<-yesCBI %>% full_join(yesCBI4, by="Name") %>% select(Name,RASTERVALU.x,RASTERVALU.y)

summary(burnSev,RASTERVALU.y)

veg_bird_burnSev <- burnSev %>% left_join(utms, by = c("Name"="point_id")) %>% 
  separate(SZ_DNS2, c("Size","Density"), sep=1) %>%
  separate(WHR_TSD, c("Forest.Type", "TSD"), sep=3) %>% 
  mutate(severity_class = as.factor(RASTERVALU.y))

veg_bird_burnSev %>% group_by(severity_class) %>% summarise(npts = n_distinct(Name))

'

