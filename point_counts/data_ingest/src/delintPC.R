# script to delint point count data from Caples Creek project
# reads in object 'PointC' from readPC.R
# creates data frame 'dfc' whose rows refer to bird counts of each species x point x visit x year.

# TODO: run readPC.R directly from this script

# for now, run readPC.R if you haven't already before starting here
PointC <- fread(here("point_counts/data_ingest/output/PointC_2022-03-23.csv"))

PC <- PointC %>% 
  filter(dist_detect != "flyover", # comment these on/off to restrict distances
         dist_detect !="100m+", 
         dist_detect !="50to100m", 
         point_ID_fk!="1072", # missing veg data for this point
         birdCode_fk != "UNKN", birdCode_fk != "0", !is.na(birdCode_fk),
         birdCode_fk != "XXHU", birdCode_fk != "XXWO", # comment on/off to include XX__ entries
         year(DateTime) != 2017,
         observer_fk != "IASH", observer_fk!="MASC") %>%
  mutate(year = year(DateTime)) %>%
  group_by(pointCount_ID_fk, point_ID_fk, DateTime, year, observer_fk, birdCode_fk) %>%
  summarise(abun=n()) %>% 
  spread(birdCode_fk, value = abun, fill = 0) %>%
  pivot_longer(AMDI:WISA, names_to = "birdCode_fk", values_to = "abun")

visits <- PointC %>% 
  mutate(year = year(DateTime)) %>%
  filter(point_ID_fk!="1072", 
         observer_fk != "IASH", observer_fk!="MASC",
         birdCode_fk != "UNKN", birdCode_fk != "XXHU", birdCode_fk != "XXWO", !is.na(birdCode_fk), 
         year != 2017) %>%
  group_by(point_ID_fk, year, DateTime) %>% 
  summarise(rich = n_distinct(fbirdName)) %>%
  arrange(DateTime, .by_group=TRUE) %>%
  mutate(visit=seq_along(DateTime))

extraVisits <- visits %>% filter(visit>3)

dfc <- full_join(PC, visits) %>% 
  ungroup() %>% 
  dplyr::select(abun, observer_fk, DateTime, point_ID_fk, year, birdCode_fk, visit) %>% 
  filter(visit < 5) 

dups <- dfc[which(duplicated(dfc[,4:7])==TRUE),] %>% group_by(point_ID_fk, year, visit) %>% summarise(nspec=n_distinct(birdCode_fk))

dfc[which(duplicated(dfc[,4:7])==TRUE),] %>% group_by(point_ID_fk, year, visit, DateTime) %>% arrange(DateTime) %>% summarise(nspec=n_distinct(birdCode_fk)) %>% arrange(DateTime) # all of these are double-observer/training counts

dfc <- dfc[which(duplicated(dfc[,4:7])==FALSE),]

# check this is 0
dfc[which(duplicated(dfc[,4:7])==TRUE),] %>% group_by(point_ID_fk, year, visit) %>% summarise(nspec=n_distinct(birdCode_fk)) # okay

dfc <- dfc %>% filter(!is.na(birdCode_fk))

dfc$observer_fk[dfc$observer_fk=="MC"] <- "MKC"
unique(dfc$observer_fk)