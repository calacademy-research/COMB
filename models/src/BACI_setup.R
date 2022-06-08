# Libraries ---------------------------------------------------------------

library(tidyverse)
library(chron)
library(lubridate)
library(googledrive)
library(here)
library(data.table)
library(reshape2)

source(here("comb_functions.R"))

# Read in data ------------------------------------------------------------

if (dir.exists(here("models/jags/")) == FALSE) {
  dir.create(here("models/jags"))
}
if (dir.exists(here("models/output/")) == FALSE) {
  dir.create(here("models/output"))
}

# TODO []: link to Google Drive &/or a link to `read_PC.R`. The Google Drive code in readPC.R does not work when using non-CalAcademy login. For now, doing a hard download.

# TODO []: why does the total # of entries differ from that in the Caples git?

PC <- fread(here("point_counts/data_ingest/output/PC_delinted.csv"))
visits <- read_csv("point_counts/data_ingest/output/PC_visit_metadata.csv")

birdTots <- PC %>% group_by(birdCode_fk) %>% summarise(tot = sum(abun)) %>% filter(tot>3)
dfc <- PC %>% filter(birdCode_fk %in% birdTots$birdCode_fk)

y <- melt(dfc, id.var = c("birdCode_fk", "point_ID_fk", "visit", "year"), measure.var = "abun") %>% acast(point_ID_fk ~ visit ~ year ~ birdCode_fk)

# COVARIATES

dates1 <- visits %>%
  filter(visit < 4) %>%
  mutate(JDay = as.numeric(format(DateTime, "%j"))) %>%
  melt(id.var = c("point_ID_fk", "year", "visit"), measure.var = "JDay") %>%
  acast(point_ID_fk ~ visit ~ year)

# standardize dates
mdate <- mean(dates1, na.rm = TRUE)
sddate <- sqrt(var(dates1[1:length(dates1)], na.rm = TRUE))
date1 <- (dates1 - mdate) / sddate

date1[is.na(date1)] <- 0

# TIME matrix [under construction - issue with time zone and converting to list of matrices]

times <- visits %>%
  filter(visit < 4) %>%
  mutate(Time = as.numeric(times(format(DateTime, "%H:%M:%S")))) %>%
  melt(id.var = c("point_ID_fk", "year", "visit"), measure.var = "Time") %>%
  acast(point_ID_fk ~ visit ~ year)

# standardize times
mtime <- mean(times, na.rm = TRUE)
sdtime <- sqrt(var(times[1:length(times)], na.rm = TRUE))
time1 <- (times - mtime) / sdtime

time1[is.na(time1)] <- 0

# [] TODO: update with google drive links &/or symlinks

veg <- fread(here("spatial/output/rasters/wide4havars.csv")) 
siteData4 <- veg %>%
  dplyr::select(point_d, Perc_LTg22mHt_2018_4ha, Perc_LTg22mHt_2020_4ha, mean_CanopyCover_2018_4ha, mean_CanopyCover_2019_4ha, mean_CanopyCover_2020_4ha,  max_Elevation_NA_4ha, mean_CaplesdNBR_Nov18Nov19_4ha, mean_RAVGcbi4_20182019_4ha) %>%
  arrange(point_d)

#utms <- read_csv("data/Wildlife_Sampling_Points.csv")

siteData4$RAVG = cut(siteData4$mean_RAVGcbi4_20182019_4ha, breaks = c(-1, 0.5, 2, 3, 4), labels=c(0,1,2,3))

siteData4 %>% pivot_longer(cols = mean_CanopyCover_2018_4ha:mean_CanopyCover_2020_4ha, names_to= "year", values_to = "canopycov") %>% 
  ggplot() +
  geom_violin(aes(x=RAVG, y=canopycov, color=year)) 

siteData4 %>% pivot_longer(cols = Perc_LTg22mHt_2018_4ha:Perc_LTg22mHt_2020_4ha, names_to= "year", values_to = "lgtree") %>% 
  ggplot() +
  geom_violin(aes(x=RAVG, y=lgtree, color=year)) 

# TODO []: deal with 1072 missing data and remove filters below once found and incorporated into RAVG data

nopc <- setdiff(siteData4$point_d, siteList) # these three aren't point count locations
setdiff(siteList, siteData4$point_d) # 1072 is missing in the metadata for some reason

dfc <- dfc %>% filter(point_ID_fk != 1072)
siteData4 <- siteData4 %>% filter(point_d != 397, point_d != 587, point_d != 845)

siteList <- sort(unique(dfc$point_ID_fk))

nsite <- length(unique(dfc$point_ID_fk))
nvisit <- length(unique(dfc$visit))
nyear <- length(unique(dfc$year))
nspec <- length(unique(dfc$birdCode_fk))

y <- melt(dfc, id.var = c("birdCode_fk", "point_ID_fk", "visit", "year"), measure.var = "abun") %>% acast(point_ID_fk ~ visit ~ year ~ birdCode_fk)

siteData4 %>% ggplot() +
  geom_boxplot(aes(x=RAVG, y=mean_CaplesdNBR_Nov18Nov19_4ha))

Cover <- scale(siteData4$mean_CanopyCover_2018_4ha)
Ht <- scale(siteData4$Perc_LTg22mHt_2018_4ha)
El <- scale(siteData4$max_Elevation_NA_4ha)
Sev <- as.numeric(as.character(siteData4$RAVG))

BA <- matrix(NA, nsite, nyear)
#colnames(BA) <- 2018:2021
#rownames(BA) <- siteList
BA[,1:2] <- 0
BA[,3:4] <- 1

BACI <- BA*Sev
BACI[BACI>0] <- 1

spp <- unique(dfc$birdCode_fk)

occ <- 1
if (occ == 1) {
  y[y > 0] <- 1
} else {
  y <- y
}

data <- list(
  nspec = nspec, # number of species (integer)
  nsite = nsite, # number of sites (integer)
  nsurvey = 3, # number of visits to each site
  nyear = nyear, # number of years (integer)
  y = y, # 4D array of detection/nondetection data
  Time = time1, # 3D array of standardized times [nsite x nsurvey x nyear]
  Date = date1,
  Ht = as.vector(Ht),
  Cov = as.vector(Cover),
  BACI = BACI
  #dNBR = as.vector(dNBR)
  # and any other variables we'd want in either part of the model (aspect, etc... but just starting with these for now)
)

save(data, Ht, Cover, BACI, Sev, siteList, spp, file = "models/output/BACI_input.RData") # need to git ignore this
