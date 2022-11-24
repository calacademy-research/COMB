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
  
  # fixing part where it tries to download folders (which throws an error)
  google_files <- google_files %>% 
    filter(
      unlist(map(1:length(google_files$drive_resource), 
                 ~ google_files$drive_resource[[.x]][["mimeType"]] != "application/vnd.google-apps.folder"))
    )
  
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
    ~ drive_upload(paste0(local_dir, "/", .x), path = as_dribble(drive_folder))
  )

  map2(
    only_google$id, only_google$name,
    ~ drive_download(.x, path = paste0(local_dir, "/", .y))
  )
}


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

# FUNCTIONS FROM JAGSUI ------------
#Functions for manipulating and extracting info from mcmc.list-class objects
#from package rjags/coda

# This is a subset of the functions in mcmc_tools in devel version 1.5.1.9024

###------------------------------------------------------------------------------
#Remove brackets and indices from parameter names in mcmc.list
strip_params <- function(params_raw, unique=FALSE){
  params_strip <- sapply(strsplit(params_raw,'[', fixed=T),'[',1)
  if(unique) return( unique(params_strip) )
  params_strip
}
#------------------------------------------------------------------------------

###------------------------------------------------------------------------------
#Identify which columns in mcmc.list object correspond to a given
#parameter name (useful for non-scalar parameters)
which_params <- function(param, params_raw){
  params_strip <- strip_params(params_raw)
  if( ! param %in% params_strip ){
    return(NULL)
  }
  which(params_strip == param)
}
#------------------------------------------------------------------------------

###------------------------------------------------------------------------------
#Get names of parameters from an mcmc.list
#If simplify=T, also drop brackets/indices
param_names <- function(mcmc_list, simplify=FALSE){
  raw <- colnames(mcmc_list[[1]])
  if(!simplify) return(raw)
  strip_params(raw, unique=T)
}
#------------------------------------------------------------------------------

###------------------------------------------------------------------------------
#Match parameter name to scalar or array versions of parameter name
match_params <- function(params, params_raw){
  unlist(lapply(params, function(x){
    if(x %in% params_raw) return(x)
    if(!x %in% strip_params(params_raw)) return(NULL)
    params_raw[which_params(x, params_raw)]
  }))
}
#------------------------------------------------------------------------------

###------------------------------------------------------------------------------
#Subset cols of mcmc.list (simple version of [.mcmc.list method)
select_cols <- function(mcmc_list, col_inds){
  out <- lapply(1:length(mcmc_list), FUN=function(x){
    mcmc_element <- mcmc_list[[x]][,col_inds,drop=FALSE]
    attr(mcmc_element,'mcpar') <- attr(mcmc_list[[x]], 'mcpar')
    class(mcmc_element) <- 'mcmc'
    mcmc_element
  })
  class(out) <- 'mcmc.list'
  out
}
#------------------------------------------------------------------------------

###------------------------------------------------------------------------------
#Convert one parameter in mcmc.list to matrix, n_iter * n_chains
mcmc_to_mat <- function(samples, param){
  psamples <- select_cols(samples, param)
  n_chain <- length(samples)
  n_iter <- nrow(samples[[1]])
  matrix(unlist(psamples), nrow=n_iter, ncol=n_chain)
}
#------------------------------------------------------------------------------
