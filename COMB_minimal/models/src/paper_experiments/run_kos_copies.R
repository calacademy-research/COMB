library(here)
library(jagsUI)
library(lubridate)
library(tidyverse)
library(chron)

source(here("comb_functions.R"))
source(here("COMB_minimal/models/src/model_read_lib_agg.R"))

# MCMC settings
na <- 100
ni <- 8000
nt <- 1
nb <- 1000
nc <- 6

site_covars <- read_csv("COMB_minimal/models/input/wide4havars.csv") %>%
  dplyr::mutate(Point = avian_point) %>%
  dplyr::select(Point,
                mean_CanopyCover_2020_4ha,
                mean_RAVGcbi_20182019_4ha,
                Perc_LTg22mHt_2020_4ha) %>%
  arrange(Point)

#' Extract time of day as decimal hours
#'
#' Converts an array of time strings (HH:MM:SS) into an array of numeric
#' values representing the hour of the day from midnight. This version is
#' vectorized for efficiency.
#'
#' @param time_mat A character array of time values. The function handles
#'   a special character "O" by converting it to `NA`.
#'
#' @return A numeric array with the same dimensions as `time_mat`, where each
#'   element is the time of day in hours.
get_hour_mat <- function(time_mat) {
  # Store original dimensions to reshape at the end
  original_dims <- dim(time_mat)
  
  time_mat[time_mat == "O"] <- NA
  decimal_hour <- hour(hms(time_mat)) + minute(hms(time_mat)) / 60
  
  # Reshape the resulting vector back to the original array dimensions.
  dim(decimal_hour) <- original_dims
  return(decimal_hour)
}

#' Calculate summary statistics of scores by site
#'
#' This function takes the data object prepared for JAGS and calculates
#' summary statistics (mean, median, max, sd) for the machine learning
#' scores for each site.
#'
#' @param jagsdata A list object returned by `readCombined()`. It must contain
#'   `siteid` (a vector of site indices for each score), `score` (a vector of
#'   scores), and `indices$point` (a data frame mapping site indices to site
#'   IDs).
#'
#' @return A data frame with summary statistics for scores at each site.
#'   The `Point` column from the input is used as row names.
get_score_stats_by_site <- function(jagsdata) {
  tibble::tibble(site_index = jagsdata$siteid,
                 score = jagsdata$score) %>%
    dplyr::left_join(jagsdata$indices$point,
                     by = c("site_index" = "Point_Index")) %>%
    dplyr::group_by(Point) %>%
    dplyr::summarize(mean_score = mean(score, na.rm = TRUE),
                     median_score = median(score, na.rm = TRUE),
                     max_score = max(score, na.rm = TRUE),
                     sd_score = sd(score, na.rm = TRUE),
                     .groups = "drop") %>%
    tibble::column_to_rownames(var = "Point")
}

#' Get naive occupancy and score statistics by site
#'
#' This function calculates naive occupancy estimates from both ARU and point
#' count data, along with summary statistics for machine learning scores. It
#' combines these into a single data frame, with one row per site.
#'
#' @param jagsdata A list object, typically the output of `readCombined()`. It
#'   must contain detection summaries (`y_aru_sum`, `y_pc_sum`), raw detection
#'   data (`y.aru`, `y.pc` with site IDs as row names), and score data used by
#'   `get_score_stats_by_site()`.
#'
#' @return A data frame with sites as rows and columns for naive occupancy,
#'   visit counts, and score statistics. The data frame will have site IDs
#'   (`Point`) as row names.
get_naive_occ_by_site_df <- function(jagsdata) {
  score_df <- get_score_stats_by_site(jagsdata)
  
  tibble::tibble(
    Point = rownames(jagsdata$y.aru),
    y_aru_sum = jagsdata$y_aru_sum,
    naive_occ_aru = as.integer(jagsdata$y_aru_sum > 0),
    y_pc_sum = jagsdata$y_pc_sum,
    naive_occ_pc = as.integer(jagsdata$y_pc_sum > 0),
    naive_occ_both = as.integer(jagsdata$y_aru_sum > 0 | jagsdata$y_pc_sum > 0),
    n_aru_visits = rowSums(!is.na(jagsdata$y.aru)),
    n_pc_visits = rowSums(!is.na(jagsdata$y.pc))
  ) %>%
    dplyr::left_join(tibble::rownames_to_column(score_df, "Point"), 
                     by = "Point") %>%
    dplyr::arrange(as.numeric(Point)) %>%
    tibble::column_to_rownames(var = "Point")
}


#' Run a knockout experiment for a single species, year, and threshold
#'
#' This function serves as a pipeline to run a suite of "knockout" models for a
#' given species, year, and detection threshold. It handles data preparation,
#' runs multiple JAGS models corresponding to different data-source
#' combinations (e.g., ARU only, Point Count only, integrated), saves the full
#' JAGS output for each model, and returns a summary data frame.
#'
#' @details
#' This function has significant side effects:
#' 1.  It saves one `.RData` file for each JAGS model run into `saveDirJags`.
#' 2.  It saves one `.RData` file containing the final occupancy summary data
#'     frame into `saveDirOcc`.
#'
#' The function assumes that the destination directories exist.
#'
#' @param th Numeric. The detection threshold for machine learning scores.
#' @param sp Character. The species code to be analyzed.
#' @param year Numeric. The year of the study.
#' @param begin_hour,end_hour Numeric. The time-of-day window (in hours) to
#'   filter the data.
#' @param daysLimit Numeric. The maximum number of distinct days to include for
#'   ARU visits.
#' @param samplesperdayLimit Numeric. The maximum number of samples per day to
#'   include for ARU visits.
#' @param modDir Character. Path to the directory containing JAGS model files.
#' @param saveDirJags Character. Path to the directory where full JAGS results
#'   will be saved.
#' @param saveDirOcc Character. Path to the directory where the summary
#'   occupancy data frame will be saved.
#' @param include_occ_covar Logical. If `TRUE`, occupancy covariates (burn,
#'   cover) are added to the JAGS data.
#' @param include_time_det_covars Logical. If `TRUE`, time-based detection
#'   covariates (hour of day) are added.
#' @param site_covars_df Data frame. A data frame containing site-level
#'   covariates, with a 'Point' column for merging.
#' @param aru_only_models Logical. If `TRUE`, only a subset of models that use
#'   ARU data are run.
#' @param monitored_params Character vector. A vector of parameter names to
#'   monitor and save from the JAGS output.
#'
#' @return A list containing two elements:
#'   \item{jagsResultlist}{A list of all `jagsUI` output objects, named by
#'     model type.}
#'   \item{occ_df}{A data frame with naive occupancy and model-estimated
#'     occupancy probabilities (`z`) for each site.}
run_ko <- function(th, sp, year, begin_hour, end_hour, daysLimit,
                   samplesperdayLimit,
                   modDir = "COMB_minimal/models/src/paper_experiments/jags_files/no_covars/",
                   saveDirJags = "COMB_minimal/models/src/paper_experiments/jagsResults/jagsResults_no_covars/",
                   saveDirOcc = "COMB_minimal/models/src/paper_experiments/occResults/no_covars/",
                   include_occ_covar = FALSE,
                   include_time_det_covars = FALSE,
                   site_covars_df = site_covars,
                   aru_only_models = FALSE,
                   monitored_params = c(
                     "mean_psi", "beta0", "beta1", "p11", "p_aru11", "p_aru01",
                     "mu", "sigma", "psi", "psi.pred.burn", "NOcc", "PropOcc",
                     "T_pc_obs", "T_pc_sim", "T_aru_obs", "T_aru_sim",
                     "TvegObs", "TvegSim", "Dobs", "Dsim", "z"
                   )) {
  # 1. Read and prepare data
  data <- readCombined(
    species = sp,
    years = c(year),
    beginTime = begin_hour,
    endTime = end_hour,
    visitAggregation = "file",
    daysLimit = daysLimit,
    samplesperdayLimit = samplesperdayLimit,
    thresholdOptions = list(value = th, is.quantile = FALSE),
    squeeze = TRUE,
    logit_col = "max_logit"
  )
  
  naive_occ_df <- get_naive_occ_by_site_df(data)
  
  # 2. Clean data object and add covariates
  if (include_occ_covar) {
    covs <- site_covars_df %>%
      dplyr::filter(Point %in% data$indices$point$Point) %>%
      dplyr::arrange(Point)
    veg <- dplyr::left_join(as.data.frame(data$indices$point), 
                            covs,
                            by = "Point")
    data$cover <- as.numeric(scale(veg$mean_CanopyCover_2020_4ha))
    data$burn  <- veg$mean_RAVGcbi_20182019_4ha
    data$Xburn <- seq(0, 2.5, length.out = 100)
  }
  
  if (include_time_det_covars) {
    data$aru_Hour <- get_hour_mat(data$aru_Time)
    data$pc_Hour <- get_hour_mat(data$pc_Time)
    cols_to_remove <- c("indices", "pc_DateTime", "pc_Time",
                        "aru_DateTime", "aru_Time", "scores_datetime")
  } else {
    cols_to_remove <- c("indices", "pc_DateTime", "pc_YDay", "pc_Time",
                        "aru_DateTime", "aru_YDay", "aru_Time", "scores_datetime")
  }
  if (all(is.na(data$y.pc))){
    cols_to_remove <- c(cols_to_remove,"y.pc")
  }
  jagsData <- within(data, rm(list = cols_to_remove))
  
  # 3. Define models and initial values
  kos <- if (aru_only_models) {
    kos <- c("model_A", "model_AS", "model_HA", "model_HAS")
  } else {
    kos <- c("model_A", "model_AS", "model_H", "model_HA",
             "model_HAS", "model_HS", "model_S")
  }
  
  inits <- function() {
    list(
      mu = c(-2, 0),
      sigma = c(1, 1),
      z = rep(1, jagsData$nsites),
      beta0 = 0,
      beta1 = 0 # Not all models use this, but JAGS handles unused inits
    )
  }
  
  # 4. Run models and save results
  results <- list()
  for (model_type in kos) {
    model_file <- file.path(modDir, paste0(model_type, ".txt"))
    
    jagsResult <- jagsUI::jags(
      data = jagsData,
      inits = inits,
      parameters.to.save = monitored_params,
      model.file = model_file,
      n.adapt = na,
      n.chains = nc,
      n.thin = nt,
      n.iter = ni,
      n.burnin = nb,
      parallel = TRUE
    )
    
    # Append model estimates to the occupancy summary data frame
    naive_occ_df[[paste("mean_z", model_type, sep = "_")]] <- jagsResult$mean$z
    naive_occ_df[[paste("median_z", model_type, sep = "_")]] <- jagsResult$q50$z
    naive_occ_df[[paste("rhat_z", model_type, sep = "_")]] <- jagsResult$Rhat$z
    
    # Save the full JAGS output object
    result_name <- paste("jagsResult", th, sp, year, model_type, Sys.Date(), sep = "_")
    save_path_jags <- file.path(saveDirJags, paste0(result_name, ".RData"))
    save(jagsResult, file = save_path_jags)
    
    results[[model_type]] <- jagsResult
  }
  
  # 5. Save the final summary data frame
  mod_str <- ifelse(aru_only_models, "aru_only_models", "all_models")
  occ_df_name <- paste("occ_df", th, sp, year, mod_str, Sys.Date(), sep = "_")
  save_path_occ <- file.path(saveDirOcc, paste0(occ_df_name, ".RData"))
  save(naive_occ_df, file = save_path_occ)
  
  return(list(jagsResultlist = results, occ_df = naive_occ_df))
}



# Single run as a test

set.seed(123)
test_results_BBWO <- run_ko(th = 0, sp = "BBWO", year=2021, begin_hour=4, end_hour=10, daysLimit=5, 
                            samplesperdayLimit=12,
                            modDir = "COMB_minimal/models/src/paper_experiments/jags_files/no_covars/",
                            saveDirJags = "COMB_minimal/models/src/paper_experiments/jagsResults/jagsResults_no_covars/",
                            saveDirOcc = "COMB_minimal/models/src/paper_experiments/occResults/no_covars/",
                            include_occ_covar = FALSE,
                            include_time_det_covars = FALSE,
                            site_covars_df = site_covars,
                            aru_only_models = FALSE)


########### Loop over species and thresholds######################
# Define the settings
monitored_HAS <- c(
  "mean_psi",
  "beta0",
  "beta1",
  "p11",
  "p_aru11",
  "p_aru01",
  "mu",
  "sigma",
  "psi",
  "psi.pred.burn",
  "NOcc", "PropOcc",
  "T_pc_obs",
  "T_pc_sim",
  "T_aru_obs",
  "T_aru_sim",
  "TvegObs",
  "TvegSim",
  "Dobs",
  "Dsim",
  "z"
)

begin_hour = 4
end_hour = 10
speciesCode <- sort(c("AMRO", "BBWO", "CONI"))
year <- 2021
threshold <- c(-2,-1,0,1,2)
aruVisitLimit <- 60 
daysLimit=5 
samplesperdayLimit=12
modDir = "COMB_minimal/models/src/paper_experiments/jags_files/no_covars/"
saveDirJags = "COMB_minimal/models/src/paper_experiments/jagsResults/jagsResults_no_covars/"
saveDirOcc = "COMB_minimal/models/src/paper_experiments/occResults/no_covars/"
if (!dir.exists(saveDirJags)) {
  dir.create(saveDirJags, recursive = TRUE) 
} 
if (!dir.exists(saveDirOcc)) {
  dir.create(saveDirOcc, recursive = TRUE) 
} 

ko_occ_results <- list()
ko_jags_results <- list()
## Iterate over thresholds and species codes
for (th in threshold){
  for (sp in speciesCode){
    ko_results <- run_ko(th = th, sp = sp, year=year, 
                         begin_hour=begin_hour, 
                         end_hour=end_hour, 
                         daysLimit=daysLimit, 
                         samplesperdayLimit=samplesperdayLimit,
                         modDir = modDir,
                         saveDirJags = saveDirJags,
                         saveDirOcc = saveDirOcc,
                         include_occ_covar = FALSE,
                         include_time_det_covars = FALSE,
                         site_covars_df = site_covars,
                         monitored_params = monitored_HAS,
                         aru_only_models = FALSE)

    ko_occ_results <- list.append(ko_occ_results, 
                                  list(species = sp, 
                                       year = year, 
                                       threshold=th, 
                                       occ_df = ko_results$occ_df))
    ko_jags_results <- list.append(ko_jags_results, ko_results$jagsResultlist)
  }
}



#############Loop over all species at threshold  = 0 ############################
rawaru <- read_csv(dataMlPath)
species_set <- unique(rawaru$species)[1:(length(unique(rawaru$species)) - 2)]
rm(rawaru)
rawpc <- read_csv(pointCountsPath)
species_pc <- unique(rawpc$birdCode_fk[rawpc$abun > 0])

setdiff(species_set, species_pc)
setdiff(species_pc, species_set)

all_species_set <- sort(intersect(species_pc, species_set))



begin_hour = 4
end_hour = 10
year <- 2021
threshold <- c(-2,0)
aruVisitLimit <- 60 
daysLimit=5 
samplesperdayLimit=12
# modDir = "COMB_minimal/models/src/paper_experiments/jags_files/single_covar/"
# saveDirJags = "COMB_minimal/models/src/paper_experiments/jagsResults/jagsResults_single_covar_all_sp/"
# saveDirOcc = "COMB_minimal/models/src/paper_experiments/occResults/single_covar_all_sp/"

modDir = "COMB_minimal/models/src/paper_experiments/jags_files/no_covars/"
saveDirJags = "COMB_minimal/models/src/paper_experiments/jagsResults/jagsResults_no_covar_all_sp_20250919/"
saveDirOcc = "COMB_minimal/models/src/paper_experiments/occResults/no_covar_all_sp_20250919/"
if (!dir.exists(saveDirJags)) {
  dir.create(saveDirJags, recursive = TRUE) 
} 
if (!dir.exists(saveDirOcc)) {
  dir.create(saveDirOcc, recursive = TRUE) 
} 

monitored_HAS <- c(
  "mean_psi",
  "p11",
  "p_aru11",
  "p_aru01",
  "mu",
  "sigma",
  "psi",
  "psi.pred.burn",
  "NOcc", "PropOcc",
  "T_pc_obs",
  "T_pc_sim",
  "T_aru_obs",
  "T_aru_sim",
  "TvegObs",
  "TvegSim",
  "Dobs",
  "Dsim",
  "z",
  "mu_diff_2_1"
)


ko_occ_results = list()
ko_jags_results = list()
threshold_set <- c(-2, 0)
# species_set <- c("AMRO", "BBWO", "CONI", "DUFL" ,"HETH", "HEWA", "MOUQ", "OSFL", "RBSA", "SOGR", "WEWP")
# species_set <- c("HEWA", "MOUQ", "OSFL", "RBSA", "SOGR", "WEWP")
for (sp in all_species_set){
  for (th in threshold_set){
    ko_results <- run_ko(th = th, sp = sp, year=year, 
                         begin_hour=begin_hour, 
                         end_hour=end_hour, 
                         daysLimit=daysLimit, 
                         samplesperdayLimit=samplesperdayLimit,
                         modDir = modDir,
                         saveDirJags = saveDirJags,
                         saveDirOcc = saveDirOcc,
                         include_occ_covar = FALSE,
                         include_time_det_covars = FALSE,
                         site_covars_df = site_covars,
                         monitored_params = monitored_HAS,
                         aru_only_models = FALSE)
  
    # ko_occ_results <- append(ko_occ_results, 
    #                               list(species = sp, 
    #                                    year = year, 
    #                                    threshold=th, 
    #                                    occ_df = ko_results$occ_df))
    # ko_jags_results <- append(ko_jags_results, ko_results$jagsResultlist)
  }
}

data_test <- data <- readCombined(
  species = species_set[2],
  years = c(year),
  beginTime = begin_hour,
  endTime = end_hour,
  visitAggregation = "file",
  daysLimit = daysLimit,
  samplesperdayLimit = samplesperdayLimit,
  thresholdOptions = list(value = th, is.quantile = FALSE),
  squeeze = TRUE,
  logit_col = "max_logit"
)
