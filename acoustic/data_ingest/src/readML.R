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
#------------------------------------------------------------

#------------------------------------------------------------
# Load libraries
#------------------------------------------------------------

# libraries required for this script
library(tidyverse)
library(data.table)
library(here)
library(googledrive)
library(lubridate)
library(fs)
library(furrr)

#------------------------------------------------------------
# Load some Caples functions
#------------------------------------------------------------

source(here("machine_learning/file_to_point/src/file_to_point.R"))

rm(list = ls())

source(here("caples_functions.R"))

#------------------------------------------------------------
# Download/organize files from google drive if not already downloaded
#------------------------------------------------------------

# Creating the folder within inputs that contains the raw_files
if (dir.exists(here("machine_learning/data_ingest/input/")) == F) {
  dir.create(here("machine_learning/data_ingest/input/"))
}

if (dir.exists(here("machine_learning/data_ingest/input/raw_files/")) == F) {
  dir.create(here("machine_learning/data_ingest/input/raw_files/"))
}

# Getting list of files that need to be downloaded
google_raw_file_names <- drive_ls(as_dribble("https://drive.google.com/drive/u/0/folders/141rV6UByPrTO15UQd9I9OBNgjR9EhW2c"), pattern = "*tar.gz")
local_tar_raw_file_names <- list.files(here("machine_learning/data_ingest/input/raw_files/"), pattern = "*tar.gz")

# Filtering google files that are not already downloaded
needed_google_raw_file_names <- google_raw_file_names %>% filter(!(google_raw_file_names$name %in% local_tar_raw_file_names))

# Downloading raw_data csv files from google drive if not already downloaded
map2(
  needed_google_raw_file_names$id, needed_google_raw_file_names$name,
  ~ drive_download(as_id(.x), path = here("machine_learning", "data_ingest", "input", "raw_files", .y))
)

# Unzipping files if not already unzipped (technically untarred?)
local_tar_raw_file_names <- list.files(here("machine_learning/data_ingest/input/raw_files/"), pattern = "*tar.gz") %>%
  str_replace(".*?_", "") %>% # Getting rid of first part of file name as untar gets rid of it
  str_sub(end = -8) # Getting rid of extension as untarring gets rid of extensions
local_untarred_raw_file_names <- list.dirs(here("machine_learning/data_ingest/input/raw_files/"), recursive = F, full.names = F)

# Filtering files that have not already been untarred
# And adding extensions and prefixes back
untar_needed_file_names <- paste0("predictions_", local_tar_raw_file_names[(!(local_tar_raw_file_names %in% local_untarred_raw_file_names))], ".tar.gz")


# Untarring all of the files which are not already untarred
# Hacky solution to needing to paste the extensions back in previous step
if (untar_needed_file_names != "predictions_.tar.gz") {
  map(
    untar_needed_file_names,
    ~ untar(here("machine_learning", "data_ingest", "input", "raw_files", .), exdir = here("machine_learning/data_ingest/input/raw_files/"))
  )
}

# Cleaning up all as everything is downloaded and not needing in memory
rm(list = ls())

# get the names in order from Tom's "info.txt" and download from drive if not already downloaded
if (file.exists(here("machine_learning/data_ingest/input/birds.txt")) == F) {
  drive_download(as_id("1zPMO_fW4nK7BwDIrp5_eDW8lgSlkRANFz8eH944-ycU"), here("machine_learning/data_ingest/input/birds.txt"))
}
sixBirdcols <- fread(here("machine_learning/data_ingest/input/birds.txt"), stringsAsFactors = FALSE, header = FALSE)

#------------------------------------------------------------
# Ingest the many small csv files and filter with awk for memory efficiency
#------------------------------------------------------------

# Defining all of the small files that we want to ingest
files <- list.files(here("machine_learning/data_ingest/input/raw_files/"), recursive = T, full.names = T, pattern = "*.csv")

# Making a string which contains all of the column numbers that awk will filter
awk <- unlist(map(
  5:95,
  ~ paste0("($", .x, ">-2)||")
)) %>%
  paste0(collapse = "") %>%
  str_sub(end = -3)

file_filter <- ""

# If filtering is needed, run the following command:
# file_filter <- "WHITE-2-CAPL"

awk_command <- paste0("awk -F, '{if(", "($1 ~ /", file_filter, "/)&&(", awk, ")){print $0}}'", " ")

# Mapping the awk reading command over all of the small dfs and making one large df with map_df()
# Using library(furrr) in order to do this over multiple cores
plan(multisession)
dataML <- future_map_dfr(
  files,
  ~ fread(cmd = paste0(awk_command, .), col.names = c("File_ID", "Start_Time", "End_Time", trimws(sixBirdcols$V2)))
) #%>%
  slice_sample(prop = .6)


#------------------------------------------------------------
# add headers and variable names
#------------------------------------------------------------

# now can rename in one fell swoop; e.g.: rename_at(vars(oldnames), ~ newnames) since dimnames and tibble don't mix (ok with dataframes though)
colnames(dataML) <- c("File_ID", "Start_Time", "End_Time", trimws(sixBirdcols$V2))

dataML %>%
  arrange(File_ID, Start_Time) -> dataML

#------------------------------------------------------------
# add the survey point number where each recording took place
#------------------------------------------------------------

# Create symlink from output of file_to_point to input of readML.R
link_create(here("machine_learning/file_to_point/output/aru2point.csv"), here("machine_learning/data_ingest/input/aru2point.csv"))

# READ IN THE ARU to actual point key value pair file from the file: 2020_ARU_data.soundfilelist.csv
aru2point2020 <- fread(here("machine_learning/data_ingest/input/aru2point.csv"), header = TRUE, sep = ",") %>%
  select(filename, point) %>%
  separate(filename, into = c("File_ID", "filetype"), sep = "\\.", remove = FALSE) %>%
  select(point, File_ID) %>%
  filter(!(is.na(point)))

# add point ID to df
dataML <- right_join(aru2point2020, dataML, by = "File_ID")

#-----------------------------------------------------------------------
# change the dates and times to posix date_time codes for analysis
#-----------------------------------------------------------------------


dataML <- dataML %>%
  separate(col = File_ID, into = c("ARU_ID", "Date", "Time"), sep = "_", remove = FALSE) %>% # pulling out date and time info
  mutate(Date = ymd(Date)) %>% # reformating date
  unite("Date_Time", Date:Time) %>% # making it into a format readable by the next line
  mutate(Date_Time = ymd_hms(Date_Time))

#------------------------------------------------------------
# save the new table and clean up
#------------------------------------------------------------

if (dir.exists(here("machine_learning/data_ingest/output/")) == F) {
  dir.create(here("machine_learning/data_ingest/output/"))
}

fwrite(dataML, here("machine_learning/data_ingest/output/dataML.csv"))


rm(list = c("ppp", "hh", "pp", "paths_to_directory_ML", "filesML", "a_file_PC", "paths_to_directory_PC", "sixBirdcols", "aru2point2018f", "aru2point2019f", "aru2point2018d", "aru2point2019d", "aru2point201819d", "aru2point201819df", "CountML", "CountMLdupes"))
