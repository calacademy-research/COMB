#############################################################
## readML.R --
## read in Tom's machine learning second results set March 2021!
## Script to get in huge dataset of Google machine learning of birds
## Clean it up and assign points to it...
## Team Caples (Jack Dumbacher, Mary Clapp, Durrell Kapan, Mark Schulist)
## built on readML.R from CAPLES/scripts/readML.R
## March 2021
## November 2021 - Mark Schulist
## New Version of readML.R without the sudo probs
## Output goes to compareML.R for comparison with selection tables
## February 2022 - Cleaned up script (a little) and made to work with COMB git - Mark Schulist

# Load libraries ------------------------------------------------------------

# libraries required for this script
library(tidyverse)
library(data.table)
library(here)
library(googledrive)
library(lubridate)
library(fs)
library(furrr)

# Load some Caples functions

source(here("comb_functions.R"))

############# OPTIONS FOR SCRIPT ###################
# If you want to run the script below, make sure that you have sufficient memory (< 64 GB) and time
# Otherwise, just download the outputs using the following command
# This will download multiple data sets, each with varying amounts of data truncated
# Name Structure:
#   dataML_m[minimum_logit]_prop[proportion_of_remaining_data]
drive_sync(here("acoustic/data_ingest/output/"), "https://drive.google.com/drive/folders/1eOrXsDmiIW9YqJWrlUWR9-Cgc7hHKD_5")

# Download/organize files from google drive if not already downloaded --------------------------------

# Creating the folder within inputs that contains the raw_files
if (dir.exists(here("acoustic/data_ingest/input/")) == F) {
  dir.create(here("acoustic/data_ingest/input/"))
}

if (dir.exists(here("acoustic/data_ingest/input/raw_files/")) == F) {
  dir.create(here("acoustic/data_ingest/input/raw_files/"))
}

# Getting list of files that need to be downloaded
google_raw_file_names <- drive_ls(as_dribble("https://drive.google.com/drive/folders/10TRgKc4NemrDeYMnBOtFmf-PS6HI5F6F"), pattern = "*tar.gz")
local_tar_raw_file_names <- list.files(here("acoustic/data_ingest/input/raw_files/"), pattern = "*tar.gz")

# Filtering google files that are not already downloaded
needed_google_raw_file_names <- google_raw_file_names %>% filter(!(google_raw_file_names$name %in% local_tar_raw_file_names))

# Downloading raw_data csv files from google drive if not already downloaded
map2(
  needed_google_raw_file_names$id, needed_google_raw_file_names$name,
  ~ drive_download(as_id(.x), path = here("acoustic", "data_ingest", "input", "raw_files", .y))
)

# Unzipping files if not already unzipped (technically untarred?)
local_tar_raw_file_names <- list.files(here("acoustic/data_ingest/input/raw_files/"), pattern = "*tar.gz") %>%
  str_sub(end = -8) # Getting rid of extension as untarring gets rid of extensions
local_untarred_raw_file_names <- list.dirs(here("acoustic/data_ingest/input/raw_files/"), recursive = F, full.names = F)

# Filtering files that have not already been untarred
untar_needed_file_names <- local_tar_raw_file_names[(!(local_tar_raw_file_names %in% local_untarred_raw_file_names))]


# Untarring all of the files which are not already untarred
map(
  paste0(untar_needed_file_names, ".tar.gz"),
  ~ untar(here("acoustic", "data_ingest", "input", "raw_files", .), exdir = here("acoustic/data_ingest/input/raw_files/"))
)

# Renaming the untarred files to match the tar.gz
# untar() removes part of the name for some reason
map2(
  .x = here("acoustic/data_ingest/input/raw_files", str_sub(untar_needed_file_names)), .y = here("acoustic/data_ingest/input/raw_files", str_sub(untar_needed_file_names, start = 13)),
  ~ file.rename(.y, .x)
)

# Cleaning up all as everything is downloaded and not needing in memory
rm(list = ls())
source(here("comb_functions.R"))

# get the names in order from Tom's "info.txt" and download from drive if not already downloaded
if (file.exists(here("acoustic/data_ingest/input/birds.txt")) == F) {
  drive_download(as_id("1--EpSKrcM1GxY20MkUM8IL8avmXH7pW37s6l4EK2JW4"), here("acoustic/data_ingest/input/birds.txt"))
}
sixBirdcols <- fread(here("acoustic/data_ingest/input/birds.txt"), stringsAsFactors = FALSE, header = FALSE)

# Ingest the many small csv files and filter with awk for memory efficiency --------------------------------------------

# Defining all of the small files that we want to ingest
files <- list.files(here("acoustic/data_ingest/input/raw_files/"), recursive = T, full.names = T, pattern = "*.csv")

# Making a string which contains all of the column numbers that awk will filter
# Change the part after the ">" in order to set the minimum logit for each row
awk <- unlist(map(
  5:95,
  ~ paste0("($", .x, ">-3)||")
)) %>%
  paste0(collapse = "") %>%
  str_sub(end = -3)

file_filter <- ""

# If filtering is needed, run the following command:
# file_filter <- "WHITE-2-CAPL"

awk_command <- paste0("awk -F, '{if(", "($1 ~ /", file_filter, "/)&&(", awk, ")){print $0}}'", " ")

# Mapping the awk reading command over all of the small dfs and making one large df with map_df()
# Using library(furrr) in order to do this over multiple cores
plan(multisession, workers = 32)
dataML <- future_map_dfr(
  files,
  ~ fread(cmd = paste0(awk_command, .), col.names = c("File_ID", "Start_Time", "End_Time", trimws(sixBirdcols$V2)))
)

# Add the survey point number where each recording took place ------------------------------------------------------------

# Create symlink from output of file_to_point to input of readML.R
if (file.exists(here("acoustic/data_ingest/input/aru2point.csv")) == F) {
  drive_download(as_id("1e2DV2eykUskpffVqIZBOjmB08x3Z_cbt5NK93B3V6JE"), here("acoustic/data_ingest/input/aru2point.csv"))
}

# READ IN THE ARU to actual point key value pair file from the file: 2020_ARU_data.soundfilelist.csv
aru2point <- fread(here("acoustic/data_ingest/input/aru2point.csv"), header = TRUE, sep = ",") %>%
  separate(Filename, into = c("File_ID", "filetype"), sep = "\\.", remove = FALSE) %>%
  filter(!(is.na(Point))) # to get rid of rows when we don't know where the ARU is located

# add point ID to df
dataML <- right_join(aru2point, dataML, by = "File_ID")

# Change the dates and times to posix date_time codes for analysis -----------------------------------------------------------


dataML <- dataML %>%
  separate(col = File_ID, into = c("ARU_ID", "Date", "Time"), sep = "_", remove = FALSE) %>% # pulling out date and time info
  mutate(Date = ymd(Date)) %>% # reformating date
  unite("Date_Time", Date:Time) %>% # making it into a format readable by the next line
  mutate(Date_Time = ymd_hms(Date_Time))

# Save the new table and clean up ------------------------------------------------------------

if (dir.exists(here("acoustic/data_ingest/output/")) == F) {
  dir.create(here("acoustic/data_ingest/output/"))
}

# Command for making output smaller by only keeping a proportion of the dataset grouped by date and time
dataML <- dataML %>% 
  group_by(Date_Time) %>% 
  slice_sample(prop = 1)

fwrite(dataML, here("acoustic/data_ingest/output/dataML.csv"))