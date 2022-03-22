#############################################################
## caples_functions.R
## Script for functions for general use across project
## Durrell D. Kapan & ...
## dkapan@calacademy.org - Sept 2019
## Mark Schulist - January 2022
#------------------------------------------------------------

#------------------------------------------------------------
# 0 - Function definitions
#------------------------------------------------------------

# excel like functions to utilize substr()
# [ ] maybe use library(stringr) in other cases
# defining some basic functions
left <- function(text, num_char) {
  substr(text, 1, num_char)
}

mid <- function(text, start_num, num_char) {
  substr(text, start_num, start_num + num_char - 1)
}

right <- function(text, num_char) {
  substr(text, nchar(text) - (num_char - 1), nchar(text))
}

# updated read for fread()
# define function to read in data w/filename
read_plus <- function(flnm, skp = 0) {
  fread(flnm, sep = "\t", skip = skp, fill = TRUE, colClasses = c("character")) %>% # did use 'colClasses=c("character")'
    mutate(filename = flnm)
}

# validate column names
# not used but will put in periods . instead of blanks if you want
validate.names <- function(df) {
  rtn <- df
  valid_column_names <- make.names(names = names(df), unique = TRUE, allow_ = TRUE)
  names(rtn) <- valid_column_names
  rtn
}

# easy negation function
NotIn <- function(x, y) !(x %in% y)

# list all objects by type (google it for location)
list.objects <- function(env = .GlobalEnv) {
  if (!is.environment(env)) {
    env <- deparse(substitute(env))
    stop(sprintf('"%s" must be an environment', env))
  }
  obj.type <- function(x) class(get(x, envir = env))
  foo <- sapply(ls(envir = env), obj.type)
  object.name <- names(foo)
  names(foo) <- seq(length(foo))
  dd <- data.frame(
    CLASS = foo, OBJECT = object.name,
    stringsAsFactors = FALSE
  )
  dd[order(dd$CLASS), ]
}

# do not remove (google it for location)
do.not.rm <- function(vars, envir = .GlobalEnv) {
  vars <- c(vars, "dnrm")
  keep <- match(x = vars, table = ls(envir = envir))
  if (any(is.na(keep))) {
    stop(paste(
      "Some of the variables were not found in",
      environmentName(envir)
    ))
  }
  rm(list = ls(envir = envir)[-keep], envir = envir)
  cat("Removed all but", length(keep), "objects from",
    environmentName(envir),
    fill = TRUE
  )
}


# superseded by rm.all.but()

# here is also conflicting with lubridate -- not sure why I had this investigate
# snippet hh
# here::here("${1}")

# readVariableWidthFile
readVariableWidthFile <- function(filePath, head = TRUE, sep = "\t", colnms = FALSE, colClasses = "character", skip = 0) {
  con <- file(filePath)
  lines <- readLines(con)
  close(con)
  slines <- strsplit(lines, sep)
  colCount <- max(unlist(lapply(slines, length)))

  FileContent <- read.table(filePath,
    header = head,
    sep = sep,
    # if(colnms = FALSE) col.names = paste0("V", seq_len(colCount)),
    fill = TRUE,
    colClasses = "character"
  )
  return(FileContent)
}

# example:
# readVariableWidthFile(paste0(paths_to_selection_tables,"/",files_selection_tables[1]))


drive_sync <- function(local_dir, drive_folder, pattern = NULL) {
  # First creating the local_dir if it does not exist
  if (dir.exists(local_dir) == FALSE) {
    dir.create(local_dir)
  }
  
  # Getting info about the directories and the files in them
  if (is.null(pattern) == TRUE) {
    google_files <- drive_ls(as_dribble(drive_folder))
  } else {
    google_files <- drive_ls(as_dribble(drive_folder)) %>%
      filter(str_detect(name, pattern = pattern))
  }

  if (is.null(pattern) == TRUE) {
    local_files <- basename(system(paste0("find ", local_dir, " -mindepth 1 -maxdepth 1 ! -type l"), intern = TRUE))
  } else {
    local_files <- basename(system(paste0("find ", local_dir, " -mindepth 1 -maxdepth 1 ! -type l"), intern = TRUE)) %>%
      str_subset(pattern = pattern)
  }

  # Comparing the two directories
  only_local <- local_files[!(local_files %in% google_files$name)]
  only_google <- google_files %>% filter(!(google_files$name %in% local_files))

  # Uploading the only_local and downloading the only_google
  map(
    only_local,
    ~ drive_upload(paste0(local_dir, .x), path = as_dribble(drive_folder))
  )

  map2(
    only_google$id, only_google$name,
    ~ drive_download(.x, path = paste0(local_dir, .y))
  )
}