# generate visit-level covariates for combined occupancy model
library(lubridate)
library(reshape2)
library(tidyverse)
library(chron)

pc <- read_csv("point_counts/data_ingest/output/PC_delinted.csv")
head(pc)

visits <- pc %>% group_by(year, point_ID_fk, visit, DateTime, observer_fk) %>%
  filter(abun !=0) %>%
  summarise(rich=n_distinct(birdCode_fk))

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

# TIME matrix

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
