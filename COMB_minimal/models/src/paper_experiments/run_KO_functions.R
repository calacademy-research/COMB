library(here)
library(jagsUI)
library(lubridate)
library(tidyverse)
library(chron)

source(here("comb_functions.R"))
source(here("COMB_minimal/models/src/model_read_lib_agg.R")) # Functions were modified to read in ARU data that was pre-aggregated# run "knockout" for each species on the list

set.seed(123)

# MCMC settings
na <- 100
ni <- 8000
nt <- 1
nb <- 1000
nc <- 6

site_covars <- read_csv("COMB_minimal/models/input/wide4havars.csv") %>%
  mutate(Point = avian_point) %>% dplyr::select(Point, mean_CanopyCover_2020_4ha, mean_RAVGcbi_20182019_4ha, Perc_LTg22mHt_2020_4ha) %>% arrange(Point)

## Function to extract the time of day in hours from midnight given a matrix of time values
get_hour_mat <- function(time_mat){
  n_dim = length(dim(time_mat))
  hour_matrix <- apply(time_mat, MARGIN = 1:n_dim, FUN = function(x) hour(hms(ifelse(x == "O", NA, x))))
  minute_matrix <- apply(time_mat, MARGIN = 1:n_dim, FUN = function(x) minute(hms(ifelse(x == "O", NA, x))))
  hour_mat <- hour_matrix + minute_matrix/60
  return(hour_mat)
}

# Function to map scores to siteid (note: requires indices created by readCombined to be present) and get aggregated statistics by site
get_score_stats_by_site <- function(jagsdata){
  score_df <- as.data.frame(cbind(jagsdata$siteid, jagsdata$score))
  names(score_df) <- c("site_index", "score")
  score_df <- merge(score_df, jagsdata$indices$point, by.x="site_index", by.y="Point_Index", all.x = TRUE)
  score_by_site <- score_df %>% 
    group_by(Point) %>%
    summarize(mean_score = mean(score), median_score = median(score), max_score = max(score),
              sd_score = sd(score)) %>%
    column_to_rownames(var = "Point")
  return(score_by_site)
}

# Function to extract the 
get_naive_occ_by_site_df <- function(df){
  occ_df <- as.data.frame(cbind(df$y_aru_sum, 
                                ifelse(df$y_aru_sum > 0, 1, df$y_aru_sum),
                                df$y_pc_sum, 
                                ifelse(df$y_pc_sum > 0, 1, 0),
                                ifelse(df$y_aru_sum + df$y_pc_sum > 0, 1, 0)
  ))
  names(occ_df) <- c("y_aru_sum", "naive_occ_aru", "y_pc_sum", "naive_occ_pc", "naive_occ_both")
  occ_df$n_aru_visits <- rowSums(!is.na(df$y.aru))
  occ_df$n_pc_visits <- rowSums(!is.na(df$y.pc))
  score_df <- get_score_stats_by_site(df)
  merged_occ_df <- merge(occ_df, score_df, by.x = 0, by.y = 0, all.x = TRUE) %>%
    arrange(as.numeric(Row.names)) %>%
    column_to_rownames(var = "Row.names")
  return(merged_occ_df)
}


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
                   )){
  #Run knockout experiment for single species, threshold, and year
  
  # First read the data
  data <- readCombined(
    species = sp,
    years = c(year),
    beginTime = begin_hour,
    endTime = end_hour,
    visitAggregation = "file",
    daysLimit = daysLimit,
    samplesperdayLimit = samplesperdayLimit, # determines how many samples per day to sample
    thresholdOptions = list(value = th,
                            is.quantile = F),
    squeeze = T,
    logit_col = "max_logit" # This is specifying we want the max_logit column from aggregated data
    # scale_datetime = T
  )
  
  naive_occ_df <- get_naive_occ_by_site_df(data)
  
  # Clean data object and add covariates
  if(include_occ_covar){
    covs <- site_covars_df %>% filter(Point %in% data$indices$point$Point) %>% arrange(Point)
    Veg <- left_join(as.data.frame(data$indices$point), covs, by = "Point")
    data$cover <- as.numeric(scale(Veg$mean_CanopyCover_2020_4ha))
    data$burn  <- Veg$mean_RAVGcbi_20182019_4ha
    data$Xburn <- seq(0, 2.5, length.out = 100)
  }
  
  if(include_time_det_covars){
    data$aru_Hour <- get_hour_mat(data$aru_Time)
    data$pc_Hour <- get_hour_mat(data$pc_Time)
    jagsData <- within(data, rm('indices','pc_DateTime', 'pc_Time',
                                "aru_DateTime", "aru_Time","scores_datetime"))
  } else{
    jagsData <- within(data, rm('indices','pc_DateTime','pc_YDay','pc_Time',
                                "aru_DateTime","aru_YDay", "aru_Time","scores_datetime"))
  }
  
  assign(paste("data", th, sp, year, sep="_"), jagsData)
  
  #Establish models
  
  if(aru_only_models){
    kos <- c("model_A", "model_AS", "model_HA", "model_HAS")
  } else{
    kos <- c("model_A", "model_AS", "model_H", "model_HA", "model_HAS", "model_HS", "model_S")
  }
  zst <- rep(1, jagsData$nsites)
  psit <- runif(jagsData$nsamples)
  
  gst <- sample(1:2, jagsData$nsamples, replace = TRUE)
  gst[jagsData$score > 0.0] <- 1
  gst[jagsData$score <= 0.0] <- 2
  
  inits <- function() {
    list(
      mu = c(-2, 0),
      sigma = c(1, 1),
      z = zst,
      g = gst,
      beta0 = 0,
      beta1 = 0,
      reg_parm = 1
    )
  }
  results <- list() 
  # Run models iterating over ko conditions
  for (k in seq_along(kos)) {
    model <- switch(kos[k],
                    "model_A"   = paste0(modDir, "model_A.txt"),
                    "model_AS"  = paste0(modDir, "model_AS.txt"),
                    "model_H"   = paste0(modDir, "model_H.txt"),
                    "model_HA"  = paste0(modDir, "model_HA.txt"),
                    "model_HAS" = paste0(modDir, "model_HAS.txt"),
                    "model_HS"  = paste0(modDir, "model_HS.txt"),
                    paste0(modDir, "model_S.txt"))  # default
    
    jagsResult <- jags(
      data = jagsData,
      inits = inits,
      parameters.to.save = monitored_params,
      model.file = model,
      n.adapt = na,
      n.chains = nc,
      n.thin = nt,
      n.iter = ni,
      n.burnin = nb,
      parallel = TRUE
    )
    naive_occ_df[[paste("mean_z", kos[k], sep = "_")]] <- jagsResult$mean$z 
    naive_occ_df[[paste("median_z", kos[k], sep = "_")]] <- jagsResult$q50$z 
    naive_occ_df[[paste("rhat_z", kos[k], sep = "_")]] <- jagsResult$Rhat$z
    
    result_name <- paste("jagsResult", th, sp, paste(year, collapse = "-"),
                         str_extract(model, "(?<=/)[^/]+(?=\\.)"), Sys.Date(), sep = "_")
    
    assign(result_name, jagsResult)
    results[[result_name]] <- jagsResult
    save(list = result_name, file = paste0(saveDirJags, result_name, ".RData"))
  }
  
  mod_str <- ifelse(aru_only_models, "aru_only_models", "all_models")
  occ_df_name <- paste("occ_df", th, sp, paste(year, collapse = "-"),mod_str, Sys.Date(), sep = "_")
  save(list = naive_occ_df, file = paste0(saveDirOcc, occ_df_name, ".RData"))
  return(list(jagsResultlist = results, occ_df = naive_occ_df))
}



# Single run as a test
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

## Iterate over thresholds and species codes
for (th in threshold){
  for (sp in speciesCode){
    ko_results <- run_ko(th = th, sp = sp, year=year, 
                         begin_hour=begin_hour, 
                         end_hour=end_hour, 
                         daysLimit=daysLimit, 
                         samplesperdayLimit=samplesperdayLimit,
                         include_occ_covar = FALSE,
                         include_time_det_covars = FALSE,
                         site_covars_df = site_covars,
                         monitored_params = monitored_HAS,
                         aru_only_models = FALSE)
    occ_df_name <- paste("occ_df", th, sp, paste(year, collapse = "-"),
                         str_extract(model, "(?<=/)[^/]+(?=\\.)"), Sys.Date(), sep = "_")
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

begin_hour = 4
end_hour = 10
year <- 2021
threshold <- 0
aruVisitLimit <- 60 
daysLimit=5 
samplesperdayLimit=12
modDir = "COMB_minimal/models/src/paper_experiments/jags_files/no_covars/"
saveDirJags = "COMB_minimal/models/src/paper_experiments/jagsResults/jagsResults_no_covars_all_sp/"
saveDirOcc = "COMB_minimal/models/src/paper_experiments/occResults/no_covars_all_sp/"
if (!dir.exists(saveDirJags)) {
  dir.create(saveDirJags, recursive = TRUE) 
} 
if (!dir.exists(saveDirOcc)) {
  dir.create(saveDirOcc, recursive = TRUE) 
} 


for (sp in speciesCode){
  ko_results <- run_ko(th = th, sp = sp, year=year, 
                       begin_hour=begin_hour, 
                       end_hour=end_hour, 
                       daysLimit=daysLimit, 
                       samplesperdayLimit=samplesperdayLimit,
                       include_occ_covar = FALSE,
                       include_time_det_covars = FALSE,
                       site_covars_df = site_covars,
                       monitored_params = monitored_HAS,
                       aru_only_models = FALSE)
  occ_df_name <- paste("occ_df", th, sp, paste(year, collapse = "-"),
                       str_extract(model, "(?<=/)[^/]+(?=\\.)"), Sys.Date(), sep = "_")
  ko_occ_results <- list.append(ko_occ_results, 
                                list(species = sp, 
                                     year = year, 
                                     threshold=th, 
                                     occ_df = ko_results$occ_df))
  ko_jags_results <- list.append(ko_jags_results, ko_results$jagsResultlist)
}


