#############################################################
## readML.subset.R --
## read in Tom's machine learning second results set March 2021!
## Script to get in PARTIAL dataset of Google machine learning of birds
## from outputs that have already been read in by readML.R
##  and subsequently filtered.  In this particular case, it loads in
## the dataML_m1.5.csv file, but that can easily be changed by making
## a substitution in the code.  The idea is to get a working copy of the
## ML data in so that other work can be done...
## Team Caples (Jack Dumbacher, Mary Clapp, Durrell Kapan, Mark Schulist)
## built on readML.R from CAPLES/scripts/readML.R
## ############################################################
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
here <- here()

dataML <- fread(paste0(here, "/acoustic/data_ingest/output/dataML_m1.5.csv"))


# Add the survey point number where each recording took place ------------------------------------------------------------

# Create symlink from output of file_to_point to input of readML.R
if (file.exists(here("acoustic/data_ingest/input/aru2point.csv")) == F) {
  drive_download(as_id("1e2DV2eykUskpffVqIZBOjmB08x3Z_cbt5NK93B3V6JE"), here("acoustic/data_ingest/input/aru2point.csv"))
}

# READ IN THE ARU to actual point key value pair file from the file: 2020_ARU_data.soundfilelist.csv
aru2point <- fread(here("acoustic/data_ingest/input/aru2point.csv"), header = TRUE, sep = ",") %>%
  separate(Filename, into = c("File_ID", "filetype"), sep = "\\.", remove = FALSE) %>%
  filter(!is.na(as.numeric(point))) %>% # throws out all that are NA or have text and are not on a point
  filter(as.numeric(point) > 0)

distinct(aru2point[, 4])
distinct(dataML[, 7])
dataML %>% distinct(Point)
distinct(dataML, Point)
distinct(aru2point, Point)

dataML %>%
  arrange(distinct(Point))

cols <- colnames(dataML)[12:100]

filter_values <- paste(cols, ">-2", "| ") %>%
  paste(collapse = "") %>%
  str_sub(end = -4)

dataML_m2.0 <- dataML %>% filter(rlang::eval_tidy(rlang::parse_expr(filter_values)))

cols <- colnames(dataML)[1:100]

#########################################
# to create a subset of data with a given
# minimum threshold across all species
##########################################
cols <- colnames(dataML)[12:100]

# For this example, the cutoff is max logit > -2

filter_values <- paste(cols, ">-0.5", "| ") %>%
  paste(collapse = "") %>%
  str_sub(end = -4)

dataML_m0.5 <- dataML %>% filter(rlang::eval_tidy(rlang::parse_expr(filter_values)))

fwrite(dataML_m0.5, here("acoustic/data_ingest/output/dataML_m0.5.csv"))



###########################################################
# to convert the logit to something more like a probability,
# we can use the formula
#
#       p = exp(logit)/(exp(logit)+1)
#
###########################################################

# define function logit_to_p

logit_to_p <- function(logit) {
  p <- exp(logit) / (exp(logit) + 1)
  return(p)
}

dataML_rebnut <- dataML %>%
  filter(logit_to_p(rebnut) > 0.1)

p <- ggplot(dataML_rebnut, aes(x = rebnut, y = (logit_to_p(rebnut)))) +
  geom_violin()

p

# Rotate the violin plot
p + coord_flip()


##############################################################
# this section of code subsets the dataML dataframe
# to provide raw data for Y[ARUi,j]
# where i is the ith point
# and j is the jth recording made on point i,
# and sum is the total number of detections exceeding
# a particular threshold (in the case below, logit > 0.5)
#
# we are additionally using only recordings taken from
# 06:00 to 10:00am, complete recordings, recordings
# made on a point.
##############################################################
dataML_model <- dataML %>%
  filter(!is.na(as.numeric(point))) %>% # throws out all that are NA or have text and are not on a point
  filter(as.numeric(point) > 0) %>% # removes recordings on point "0"
  # select(filename, File_ID, ARU_ID, Date_Time, point, Start_Time, rebnut) %>% # choose only the focal species
  # mutate(P = logit_to_p(rebnut)) %>%    # if you want to add prob
  filter(year(Date_Time) == 2021) %>% # to filter only 2021
  filter(hour(Date_Time) >= 6, hour(Date_Time) < 10) %>% # to filter morning obs
  filter(minute(Date_Time) == 00 | minute(Date_Time) == 30) # to remove partial recordings

fwrite(dataML_model, here("acoustic/data_ingest/output/dataML_model.csv"))

dataML_sample <- dataML_model %>%
  group_by(File_ID, point, Date_Time) %>%
  summarize(sum = sum(rebnut > 0)) # this restricts data to single species (rebnut) and counts only logits greater than a certain cutoff (logit>0)


dataML_point <- dataML_sample %>%
  group_by(point) %>%
  summarize(n = sum(sum > 0))

p <- ggplot(data = dataML_point, aes(x = n)) +
  geom_histogram(binwidth = 6)
p
