library(chron)
library(data.table)
library(forcats)
library(googledrive)
library(here)
library(lubridate)
library(reshape2)
library(tidyverse)

# script to delint point count data from Caples Creek project
# creates data frame 'dfc' whose rows refer to bird counts of each species x point x visit x year.

# for now, run readPC.R if you haven't already before starting here
PointC <- fread(here("point_counts/data_ingest/output/PointC.csv"))

PointC$observer_fk[PointC$observer_fk == "MC"] <- "MKC"

PointC$DateTime <- with_tz(PointC$DateTime, tzone = "America/Los_Angeles")

PC <- PointC %>%
  filter(
    dist_detect != "flyover", # comment these on/off to restrict distances
    dist_detect != "100m+",
    # dist_detect !="50to100m",
    # point_ID_fk != "1072", # missing veg data for this point
    birdCode_fk != "UNKN", birdCode_fk != "0", !is.na(birdCode_fk),
    birdCode_fk != "XXHU", birdCode_fk != "XXWO", birdCode_fk != "XXWA",
    birdCode_fk != "SPHY", birdCode_fk != "GCCR",# comment on/off to include XX__ entries
    year(DateTime) != 2017,
    observer_fk != "IASH", observer_fk != "MASC"
  ) %>%
  mutate(year = year(DateTime)) %>%
  group_by(pointCount_ID_fk, point_ID_fk, DateTime, year, observer_fk, birdCode_fk) %>%
  summarise(abun = n()) %>%
  spread(birdCode_fk, value = abun, fill = 0) %>%
  pivot_longer(AMDI:WIWA, names_to = "birdCode_fk", values_to = "abun")

visits <- PointC %>%
  mutate(year = year(DateTime)) %>%
  filter(
    # point_ID_fk != "1072",
    observer_fk != "IASH", observer_fk != "MASC",
   # birdCode_fk != "UNKN", birdCode_fk != "XXHU", birdCode_fk != "XXWO", !is.na(birdCode_fk),
    year != 2017
  ) %>%
  group_by(point_ID_fk, year, DateTime) %>%
  summarise(rich = n_distinct(fbirdName)) %>%
  arrange(DateTime, .by_group = TRUE) %>%
  mutate(visit = seq_along(DateTime))

extraVisits <- visits %>% filter(visit > 3)

dfc <- full_join(PC, visits) %>%
  ungroup() %>%
  dplyr::select(abun, observer_fk, DateTime, point_ID_fk, year, birdCode_fk, visit) %>%
  filter(visit < 4)

dups <- dfc[which(duplicated(dfc[, 4:7]) == TRUE), ] %>%
  group_by(point_ID_fk, year, visit) %>%
  summarise(nspec = n_distinct(birdCode_fk))

dfc[which(duplicated(dfc[, 4:7]) == TRUE), ] %>%
  group_by(point_ID_fk, year, visit, DateTime) %>%
  arrange(DateTime) %>%
  summarise(nspec = n_distinct(birdCode_fk)) %>%
  arrange(DateTime) # all of these are double-observer/training counts

dfc <- dfc[which(duplicated(dfc[, 4:7]) == FALSE), ]

# check this is 0
dfc[which(duplicated(dfc[, 4:7]) == TRUE), ] %>%
  group_by(point_ID_fk, year, visit) %>%
  summarise(nspec = n_distinct(birdCode_fk)) # okay

dfc <- dfc %>% filter(!is.na(birdCode_fk))

# check 4-letter codes against ML's bird name index:
if (file.exists(here("acoustic/data_ingest/input/birds.txt")) == F) {
  drive_download(as_id("1--EpSKrcM1GxY20MkUM8IL8avmXH7pW37s6l4EK2JW4"), here("acoustic/data_ingest/input/birds.txt"))
}
sixBirdcols <- as.data.frame(fread(here("acoustic/data_ingest/input/birds.txt"), stringsAsFactors = FALSE, header = FALSE))
colnames(sixBirdcols)

setdiff(dfc$birdCode_fk,sixBirdcols$V3) # in point count set but not ML set
setdiff(sixBirdcols$V3, dfc$birdCode_fk) # in ML set but not in point count set

dfc$birdCode_fk[which(dfc$birdCode_fk=="AUWA")] <- "YRWA"
dfc$birdCode_fk[which(dfc$birdCode_fk=="ORJU")] <- "DEJU"
dfc$birdCode_fk[which(dfc$birdCode_fk=="NOPY")] <- "NOPO"

write_csv(dfc, "point_counts/data_ingest/output/PC_delinted.csv")
write_csv(visits, "point_counts/data_ingest/output/PC_visit_metadata.csv")

# #2021-11-29 QA/QC
# visits %>%
#   select(point_ID_fk, year, visit) %>% 
#   group_by(point_ID_fk) %>%
#   rowwise(.) %>%
#   mutate(maxvisit = max(visit)) %>% View()
#   pivot_wider(., names_from = year, values_from = visit, ) %>% View()
#   
# table(visits$,visits$year) < 3 %>% 
#   as_data_frame(.) %>% View() 
#   as_tibble(.) %>% 
#   mutate(sum = rowSums(across(c(1:4)))) %>% View()
  