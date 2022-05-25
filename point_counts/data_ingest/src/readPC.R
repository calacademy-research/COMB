## ReadPC.R
## Script for reading in point count data
## Durrell D. Kapan, Mary Clapp, Jack Dumbacher & Ian Shriner
## dkapan@calacademy.org - Sept 2019, May 2020 for eBird 6 codes
## July 2020 debugging and cleanup
## June 2021 new dataset

## December 2021 - Reorganized Data Structure
## Mostly copied from old read_PC.R
## Mark Schulist
rm(list = ls())
# Load Libraries and Caples Scripts -----------
library(tidyverse)
library(here)
library(googledrive)
library(lubridate)
library(auk)
library(forcats)
library(data.table)

here <- here()

source(here("comb_functions.R"))


# Syncing file with google drive --------------

drive_sync(here("point_counts/data_ingest/input/"), "https://drive.google.com/drive/u/0/folders/13bXkNCZzFwGC8H4k-CJ4Cf3tYJ3t18zW")

# Reading inputs -------------------

if (dir.exists(here("point_counts/data_ingest/input/")) == FALSE) {
  dir.create(here("point_counts/data_ingest/input/"))
}

PCTC <- fread(here("point_counts/data_ingest/input/2021-06-17_2017-2021_Caples_data.out.csv"))

# Do some QA/QC -------------------
PCTC %>%
  dplyr::select(point_ID_fk, Date) %>%
  mutate(Year = lubridate::year(mdy(Date))) %>%
  distinct() %>%
  group_by(Year, point_ID_fk) %>%
  tally() %>%
  spread(Year, n) %>%
  View()

# Filter out 'bad' or 'nonconforming data' --------
# [ ]

# Making data into "long" data -------------
PCTClong <- PCTC %>%
  group_by(birdCode_fk) %>%
  pivot_longer(c(
    -Date, -point_ID_fk, -timeStart, -observer_fk, -interval_5to10_valid, -Sampling_Unit_fk,
    -`Caples Watershed`, -TRTMT, -CONTROL, -birdCode_fk, -birdName, -TOTAL_individuals_0to10min,
    -taxonCount_ID_pk, -pointCount_ID_fk, -notes
  ), names_to = "orig_str", values_to = "counts")


# PCTClong <- PCTC %>%
#   group_by(birdCode_fk) %>%
#   gather(key = "orig_str", value = "counts", -Date, -point_ID_fk, -timeStart, -observer_fk,
#          -interval_5to10_valid, -Sampling_Unit_fk, -`Caples Watershed`, -TRTMT, -CONTROL, -birdCode_fk, -birdName, -TOTAL_individuals_0to10min, -taxonCount_ID_pk, -pointCount_ID_fk, -notes)
# Re-coding the key-value pairs into different columns -----------
# time_period define
PCTClong$time_period <- PCTClong$orig_str
# recode
PCTClong$time_period[str_detect(PCTClong$orig_str, "0to5min")] <- "0to5min"
PCTClong$time_period[str_detect(PCTClong$orig_str, "5to10min")] <- "5to10min"

# distance define
PCTClong$dist_detect <- PCTClong$orig_str
# recode
PCTClong$dist_detect[str_detect(PCTClong$orig_str, "0to25m")] <- "0to25m"
PCTClong$dist_detect[str_detect(PCTClong$orig_str, "25to50m")] <- "25to50m"
PCTClong$dist_detect[str_detect(PCTClong$orig_str, "50to100m")] <- "50to100m"
PCTClong$dist_detect[str_detect(PCTClong$orig_str, "moreThan100m")] <- "100m+"
PCTClong$dist_detect[str_detect(PCTClong$orig_str, "flyover")] <- "flyover"

# fix levels
PCTClong$dist_detect <- ordered(PCTClong$dist_detect, levels = c("0to25m", "25to50m", "50to100m", "100m+", "flyover"))

# new or continuing
PCTClong$NEWorCONT <- PCTClong$orig_str
# recode
PCTClong$NEWorCONT[str_detect(PCTClong$orig_str, "NEW")] <- "NEW"
PCTClong$NEWorCONT[str_detect(PCTClong$orig_str, "CONT")] <- "CONT"
PCTClong$NEWorCONT[str_detect(PCTClong$orig_str, "NEW|CONT", negate = TRUE)] <- "NA"

# detection type define
PCTClong$detect_type <- PCTClong$orig_str

# recode
# KEY: HO = heard only; seen_call = seen and heard call; seen_song = seen and heard song; SO = seen only; unspec = unspecified
PCTClong$detect_type[str_detect(PCTClong$orig_str, "ALL")] <- "ALL"
PCTClong$detect_type[str_detect(PCTClong$orig_str, "HO_call")] <- "HO_call"
PCTClong$detect_type[str_detect(PCTClong$orig_str, "HO_other")] <- "HO_other"
PCTClong$detect_type[str_detect(PCTClong$orig_str, "HO_song")] <- "HO_song"
PCTClong$detect_type[str_detect(PCTClong$orig_str, "seen_call")] <- "seen_call"
PCTClong$detect_type[str_detect(PCTClong$orig_str, "seen_song")] <- "seen_song"
PCTClong$detect_type[str_detect(PCTClong$orig_str, "SO")] <- "SO"
PCTClong$detect_type[str_detect(PCTClong$orig_str, "unspec")] <- "unspec"

# Making Dates into datetime objects ----------
PCTClong$Date <- as.Date(PCTClong$Date, "%m/%d/%Y")

# Getting the left 2 digits and the right 2 digits to get the complete time regardless of user input (12:34 vs 1234)
PCTClong$timeStart <- PCTClong$timeStart %>%
  as.character() %>%
  str_pad(4, pad = "0")

PCTClong$timeStart <- paste0(left(PCTClong$timeStart, 2), ":", right(PCTClong$timeStart, 2))

PCTClong$DateTime <- as.POSIXct(strptime(paste(PCTClong$Date, PCTClong$timeStart), "%Y-%m-%d %H:%M"))

# Checking for 4 letter banding code mistakes ------------------
PCTClong %>%
  ungroup() %>%
  mutate(birdCode_fk = ifelse(grepl("TRSW", birdCode_fk), "TRES", birdCode_fk)) -> PCTClong

PCTClong %>%
  ungroup() %>%
  mutate(birdCode_fk = ifelse(grepl("RSFL", birdCode_fk), "NOFL", birdCode_fk)) -> PCTClong

PCTClong %>%
  ungroup() %>%
  mutate(counts = ifelse(grepl("1`", counts), "1", counts)) %>%
  mutate(counts = as.numeric(counts)) -> PCTClong

# Fixing the bad species names due to user input error -----------------

# Defining species that need to be fixed
to_fix_spp <- c("Audubon's (Yellow-rumped) Warbler", "Nashville Warbler (Western)", "Oregon (Dark-eyed) Junco", "Pileated Woodpecker (Western)", "Red-shafted (Northern) Flicker", "Warbling Vireo (Western)", "White-breasted Nuthatch (Great Basin)", "Williamson's / Red-breasted Sapsucker", "Wilson's Warbler (Western)", "Unknown / Uncertain", "")
fixed_spp <- c("Yellow-rumped Warbler", "Nashville Warbler", "Dark-eyed Junco", "Pileated Woodpecker", "Northern Flicker", "Warbling Vireo", "White-breasted Nuthatch", "Sapsucker sp.", "Wilson's Warbler", "Uk", NA)

# Adding a duplicate column to do the error correction on
PCTClongf <- PCTClong %>% mutate(fbirdName = birdName)

# quick  loop over these ^ and fix the overall dataframe
for (i in 1:length(fixed_spp)) {
  PCTClongf <- PCTClongf %>% mutate(fbirdName = case_when(birdName %in% to_fix_spp[c(i)] ~ fixed_spp[c(i)], TRUE ~ fbirdName))
}

# Put the column back into PCTClong
PCTClongf %>%
  mutate(birdName = fbirdName) -> PCTClong

# Clean up
rm(PCTClongf, fixed_spp, i, to_fix_spp)

# Add 6 letter codes using eBird ----------------
PCTClong %>%
  mutate(eBird_6_code = ebird_species(birdName, "code")) -> PCTClong

# Write file to output --------------------
# Make sure the output folder exists
if (dir.exists(here("point_counts/data_ingest/output/")) == FALSE) {
  dir.create(here("point_counts/data_ingest/output/"))
}

# Finally writing
# Getting input file name and adding ".long" for output
long_input_file_name <- list.files(here("point_counts/data_ingest/input/")) %>%
  str_sub(end = -5) %>%
  paste0(".long.csv")

# Writing
fwrite(PCTClong, here("point_counts/data_ingest/output/", long_input_file_name))

# PointC ------------------------------------------------------------------


# Below, PointC is an object that summarizes PCTClong into a matrix where every row is an observation of a species within a certain distance band and detected using a certain detection type. Most entries are 1 except where multiple birds of the same species are detected in the same distance bin and using the same type. This data frame can be further summarised to recover a species count per point visit (as is done in ./models/occupancy_setup.R), filtered to restrict distance, select certain detection types, etc.).

# IMPORTANT NOTE: The files saved here contain ALL point count data, including training data and counts made on days of inclement weather, etc. Further post-processing for downstream analyses occur in other scripts.

PointC <- PCTClong %>%
  mutate_if(is.factor, as.character) %>% # changes everything to characters also speedier
  filter(NEWorCONT != "CONT", detect_type != "ALL") %>%
  group_by(pointCount_ID_fk, point_ID_fk, DateTime, observer_fk, time_period, dist_detect, detect_type, NEWorCONT, birdCode_fk, fbirdName) %>%
  summarise(tot = sum(counts)) %>%
  filter(tot > 0) %>%
  arrange(DateTime, point_ID_fk) # puts it in order for the entire dataset

# Write PointC to output folder
PointC_filename <- "PointC.csv"
fwrite(PointC, here("point_counts/data_ingest/output/", PointC_filename))


# find duplicates of pointCount_ID_fk (the record number from Filemaker Pro): it doesn't matter if there are duplicates there as long as the data are summarized using unique combinations of year, site, and date-time.
PointC %>%
  group_by(pointCount_ID_fk, DateTime, point_ID_fk) %>%
  summarise(rich = n_distinct(birdCode_fk)) %>%
  filter(point_ID_fk != 0) %>%
  View()

# old notes ---------------------------------------------------------------


# NOTE check against the final database ...
# some checks
# function to find missing columns
# NotIn <- function(x,y) !(x %in% y)
# #columns in failing dplyr call
# grplst<-c("not", "pointCount_ID_fk", "point_ID_fk", "DateTime", "observer_fk", "time_period", "dist_detect", "detect_type", "NEWorCONT", "birdCode_fk", "birdName")
#
# #quick test (including a true positive 'not')
# NotIn(grplst[1:10],colnames(PCTClong))
