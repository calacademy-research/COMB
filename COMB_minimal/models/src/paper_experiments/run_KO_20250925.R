# Code to iterate 7 different "knockout" model structures over a range of species
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


speciesCode <- sort(c("AMRO", "BBWO", "CONI"))

year <- 2021
threshold <- c(-2,-1,0,1,2)
aruVisitLimit <- 60 # only consider this many ARU visits per site (ordered). 32=4 full mornings


# read in/prepare field data ----------------------------------------------

# SITE covars
site_covars <- read_csv("COMB_minimal/models/input/wide4havars.csv") %>%
  mutate(Point = avian_point) %>% dplyr::select(Point, mean_CanopyCover_2020_4ha, mean_RAVGcbi_20182019_4ha, Perc_LTg22mHt_2020_4ha) %>% arrange(Point)

# JAGS structuring --------------------------------------------------------

for (sp in speciesCode) {
  begin_hour = 4
  end_hour = 10
  data <- readCombined(
    species = c(sp),
    years = c(year),
    beginTime = begin_hour,
    endTime = end_hour,
    visitLimit = NA,
    visitAggregation = "file",
    daysLimit = 4,
    samplesperdayLimit = 10,
    thresholdOptions = list(value = threshold,
                            is.quantile = F),
    squeeze = T,
    logit_col = "max_logit", # This is specifying we want the max_logit column from aggregated data
    scale_datetime = T
  )
  jagsData <- within(data, rm('indices','pc_DateTime','pc_YDay','pc_Time',"aru_DateTime","aru_YDay", "aru_Time","scores_datetime")) #TODO: change over the date computations to pull directly from 'data' object
  
  jagsData$cover <- as.numeric(scale(Veg$mean_CanopyCover_2020_4ha)) # pull covs of interest from preloaded Veg data
  jagsData$burn <- Veg$mean_RAVGcbi_20182019_4ha
  jagsData$Xburn = seq(0, 2.5, length.out=100) # create a vector to predict on
  
  assign(paste0("data_", sp), jagsData)
}

get_hour_mat <- function(time_mat){
  n_dim = length(dim(time_mat))
  hour_matrix <- apply(time_mat, MARGIN = 1:n_dim, FUN = function(x) hour(hms(ifelse(x == "O", NA, x))))
  minute_matrix <- apply(time_mat, MARGIN = 1:n_dim, FUN = function(x) minute(hms(ifelse(x == "O", NA, x))))
  hour_mat <- hour_matrix + minute_matrix/60
  return(hour_mat)
}

begin_hour = 4
end_hour = 9
sp = "BBWO"
th = -1

data <- readCombined(
  species = c(sp),
  years = c(year),
  beginTime = begin_hour,
  endTime = end_hour,
  visitLimit = NA,
  visitAggregation = "file",
  daysLimit = 5,
  samplesperdayLimit = 10,
  thresholdOptions = list(value = th,
                          is.quantile = F),
  squeeze = T,
  logit_col = "max_logit", # This is specifying we want the max_logit column from aggregated data
  scale_datetime = F
)
  
covs <- site_covars %>% filter(Point %in% data$indices$point$Point) %>% arrange(Point)
Veg <- left_join(as.data.frame(data$indices$point), covs, by = "Point")
data$pc_hour <- get_hour_mat(data$pc_Time)
data$aru_hour <- get_hour_mat(data$aru_Time)
# jagsData <- within(data, rm('indices','pc_DateTime','pc_YDay','pc_Time',"aru_DateTime","aru_YDay", "aru_Time","scores_datetime")) #TODO: change over the date computations to pull directly from 'data' object
jagsData <- within(data, rm('indices','pc_Time','pc_DateTime',"aru_DateTime", "aru_Time","scores_datetime"))
# minimum_pc_yday = min(debug_read_data$pc_YDay[debug_read_data$pc_YDay >0])

jagsData$cover <- as.numeric(scale(Veg$mean_CanopyCover_2020_4ha)) # pull covs of interest from preloaded Veg data
jagsData$burn <- Veg$mean_RAVGcbi_20182019_4ha
jagsData$Xburn = seq(0, 2.5, length.out=100) # create a vector to predict on

assign(paste("data", th, sp, year, sep="_"), jagsData)


datasets <- setNames(
  lapply(speciesCode, function(sp) get(assign(paste("data", th, sp, year, sep="_"), jagsData))),
  speciesCode
)




# JAGS structuring --------------------------------------------------------



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


run_ko <- function(th, sp, year, begin_hour, end_hour, daysLimit, 
                   samplesperdayLimit,
                   modDir = "COMB_minimal/models/src/paper_experiments/jags_files/no_covars/",
                   include_occ_covar = FALSE,
                   include_time_det_covars = FALSE,
                   site_covars_df = site_covars,
                   monitored_params = monitored_HAS,
                   aru_only_models = FALSE){
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
  zst <- rep(1, data$nsites)
  psit <- runif(data$nsamples)
  
  gst <- sample(1:2, data$nsamples, replace = TRUE)
  gst[data$score > 0.0] <- 1
  gst[data$score <= 0.0] <- 2
  
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
    result_name <- paste("jagsResult", th, sp, paste(year, collapse = "-"),
                         str_extract(model, "(?<=/)[^/]+(?=\\.)"), Sys.Date(), sep = "_")
    
    assign(result_name, jagsResult)
    results[[result_name]] <- jagsResult
    save(list = result_name, file = paste0("COMB_minimal/models/src/paper_experiments/jagsResults/jagsResults_single_covar/", result_name, ".RData"))
}

datasets <- setNames(
  lapply(speciesCode, function(sp) get(paste("data", th, sp, sep="_"))),
  speciesCode
)


# monitored elements ------------------------------------------------------
# FULL (HUMAN + ARU + SCORES)
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

#kos <- c("model_A", "model_AS", "model_H", "model_HA", "model_HAS", "model_HS", "model_S")
kos <- c("model_A", "model_AS", "model_HA", "model_HAS") # only ones we need for thresh expt
#kos <- "model_HAS"
#modDir <- "COMB_minimal/models/jags_files/"
modDir <- "COMB_minimal/models/jags_files/no_covars/"

#modDir <- "models/jags_files/"
#model <- paste0(modDir,"model_HAS_multiyear_nocovar.txt")
#sp <- "AMRO"
#data <- data_AMRO

# Inits lists -------------------------------------------------------------

results <- list()  # optional: store all results here

for (sp in names(datasets)) {
  data <- datasets[[sp]]
  
  # Create inits components specific to this data
  zst <- rep(1, data$nsites)
  psit <- runif(data$nsamples)
  
  gst <- sample(1:2, data$nsamples, replace = TRUE)
  gst[data$score > 0.0] <- 1
  gst[data$score <= 0.0] <- 2
  
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
      data = data,
      inits = inits,
      parameters.to.save = monitored_HAS,
      model.file = model,
      n.adapt = na,
      n.chains = nc,
      n.thin = nt,
      n.iter = ni,
      n.burnin = nb,
      parallel = TRUE
    )
    
    # for each result, give it a name and put it in the results list
    result_name <- paste("jagsResult", th, sp, paste(year, collapse = "-"),
                         str_extract(model, "(?<=/)[^/]+(?=\\.)"), Sys.Date(), sep = "_")
    
    assign(result_name, jagsResult)
    results[[result_name]] <- jagsResult
  }
  #}
} # end of species dataset loop
rm(jagsResult)

jags_objects <- ls(pattern = "^jagsResult")

# Loop through the list and save each object as an RData file
for (obj in jags_objects) {
  save(list = obj, file = paste0("COMB_minimal/results/jagsResults/jagsResults_no_covars/", obj, ".RData"))
}

}





for(th in threshold) {

  for (sp in speciesCode) {
    begin_hour = 4
    end_hour = 9
    data <- readCombined(
      species = "AMRO",
      years = c(year),
      beginTime = dhours(begin_hour),
      endTime = dhours(end_hour),
      visitLimit = aruVisitLimit, # total number of files
      visitAggregation = "file",
    #  daysLimit = 4,
    #  samplesperdayLimit = 10, # determines how many samples per day to sample
      thresholdOptions = list(value = th,
                              is.quantile = F),
      squeeze = T,
      logit_col = "max_logit" # This is specifying we want the max_logit column from aggregated data
     # scale_datetime = T
    )
  
    covs <- site_covars %>% filter(Point %in% data$indices$point$Point) %>% arrange(Point)
    
    Veg <- left_join(as.data.frame(data$indices$point), covs, by = "Point")
    
    # Clean data object and add covariates
    jagsData <- within(data, rm('indices','pc_DateTime','pc_YDay','pc_Time',
                                "aru_DateTime","aru_YDay", "aru_Time","scores_datetime"))
    
    jagsData$cover <- as.numeric(scale(Veg$mean_CanopyCover_2020_4ha))
    jagsData$burn  <- Veg$mean_RAVGcbi_20182019_4ha
    jagsData$Xburn <- seq(0, 2.5, length.out = 100)
    
    assign(paste("data", th, sp, sep="_"), jagsData)
  }
  
  datasets <- setNames(
    lapply(speciesCode, function(sp) get(paste("data", th, sp, sep="_"))),
    speciesCode
  )
  
  
  # monitored elements ------------------------------------------------------
  # FULL (HUMAN + ARU + SCORES)
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
  
  #kos <- c("model_A", "model_AS", "model_H", "model_HA", "model_HAS", "model_HS", "model_S")
  kos <- c("model_A", "model_AS", "model_HA", "model_HAS") # only ones we need for thresh expt
  #kos <- "model_HAS"
  #modDir <- "COMB_minimal/models/jags_files/"
  modDir <- "COMB_minimal/models/jags_files/no_covars/"
  
  #modDir <- "models/jags_files/"
  #model <- paste0(modDir,"model_HAS_multiyear_nocovar.txt")
  #sp <- "AMRO"
  #data <- data_AMRO
  
  # Inits lists -------------------------------------------------------------
  
  for (sp in names(datasets)) {
    data <- datasets[[sp]]
    
    # Create inits components specific to this data
    zst <- rep(1, data$nsites)
    psit <- runif(data$nsamples)
    
    gst <- sample(1:2, data$nsamples, replace = TRUE)
    gst[data$score > 0.0] <- 1
    gst[data$score <= 0.0] <- 2
    
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
        data = data,
        inits = inits,
        parameters.to.save = monitored_HAS,
        model.file = model,
        n.adapt = na,
        n.chains = nc,
        n.thin = nt,
        n.iter = ni,
        n.burnin = nb,
        parallel = TRUE
      )
      
      # for each result, give it a name and put it in the results list
      result_name <- paste("jagsResult", th, sp, paste(year, collapse = "-"),
                           str_extract(model, "(?<=/)[^/]+(?=\\.)"), Sys.Date(), sep = "_")
      
      assign(result_name, jagsResult)
      results[[result_name]] <- jagsResult
    }
  #}
  } # end of species dataset loop
  rm(jagsResult)
  
  jags_objects <- ls(pattern = "^jagsResult")
  
    # Loop through the list and save each object as an RData file
    for (obj in jags_objects) {
      save(list = obj, file = paste0("COMB_minimal/results/jagsResults/jagsResults_no_covars/", obj, ".RData"))
    }
}


str_split(jags_objects[1], pattern = "_")