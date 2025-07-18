# Code to iterate 7 different "knockout" model structures over a range of species

library(googledrive)
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


# read in/prepare field data ----------------------------------------------

# parameters --------------------------------------------------------------
speciesCode <- sort(c("BBWO", "WHWO", "NOFL", "HAWO", "PIWO", "RBSA")) # must match prefiltering of dataML_model.csv
year <- 2021
threshold <- 0
aruVisitLimit <- 24 # only consider this many ARU visits per site (ordered). 32=4 mornings

# JAGS structuring --------------------------------------------------------

for (sp in speciesCode) {
  data <- readCombined(
    species = c(sp),
    years = c(year),
    beginTime = dhours(5),
    endTime = dhours(9),
    visitLimit = aruVisitLimit,
    visitAggregation = "file",
    samplesperdayLimit = 24,
    thresholdOptions = list(value = threshold,
                            is.quantile = F),
    squeeze = T,
    logit_col = "max_logit", # This is specifying we want the max_logit column from aggregated data
    scale_datetime = T
  )
  jagsData <- within(data, rm('indices','pc_DateTime','pc_YDay','pc_Time',"aru_DateTime","aru_YDay", "aru_Time","scores_datetime")) #TODO: change over the date computations to pull directly from 'data' object

#  jagsData$cover <- as.numeric(scale(Veg$mean_CanopyCover_2020_4ha)) # pull covs of interest from preloaded Veg data
#  jagsData$burn <- Veg$mean_RAVGcbi_20182019_4ha
#  jagsData$Xburn = seq(0, 2.5, length.out=100) # create a vector to predict on

  assign(paste0("data_", sp), jagsData)
}

datasets <- list(BBWO=data_BBWO, HAWO=data_HAWO, NOFL=data_NOFL, PIWO=data_PIWO, RBSA=data_RBSA, WHWO=data_WHWO)


# Inits lists -------------------------------------------------------------

# FULL
# these initial values
zst <- rep(1, data_HAWO$nsites)
psit <- runif(data_HAWO$nsamples)

# any model containing scores needs:
gst <- sample(1:2, data_HAWO$nsamples, replace = TRUE)
gst[data_HAWO$score > 0.0] <- 1
gst[data_HAWO$score <= 0.0] <- 2

inits <- function() {
  list(
    mu = c(-2, 0),
    sigma = c(1, 1),
    z = zst,
    psi = psit,
 #   p11 = runif(1, 0.2, 0.8),
    #  lambda = runif(1, 1, 20), omega = runif(1, 0, 15),
    g = gst,
#    beta0 = 0,
#    beta1 = 0,
    #   beta2=0,
    reg_parm = 1
  )
}


#PC ONLY

inits.pc <- function() {
  list(
    #   mu = c(-2, 0),
    #    sigma = c(1, 1),
    z = zst,
    psi = psit,
#    p11 = runif(1, 0.2, 0.8),
    #  lambda = runif(1, 1, 20), omega = runif(1, 0, 15),
    #   g = gst,
 #   beta0 = 0,
 #   beta1 = 0,
    #   beta2=0,
    reg_parm = 1
  )
}


# monitored elements ------------------------------------------------------
# FULL (HUMAN + ARU + SCORES)
monitored_HAS <- c(
  "mean_psi",
#  "beta0",
#  "beta1",
  "p11",
  "p_aru11",
  "p_aru01",
  "mu",
  "sigma",
  "psi",
#  "psi.pred.burn",
  "NOcc", "PropOcc",
  "T_pc_obs",
  "T_pc_sim",
  "T_aru_obs",
  "T_aru_sim",
#  "TvegObs",
#  "TvegSim",
  "Dobs",
  "Dsim", "z"
)

# HUMAN ONLY 
monitored_H <- c(
  "mean_psi",
#  "beta0",
#  "beta1",
  "p11",
#  "p_aru11",
#  "p_aru01",
#  "mu",
#  "sigma",
  "psi",
#  "psi.pred.burn",
"NOcc", "PropOcc",
    "T_pc_obs",
    "T_pc_sim", "z"
  #  "T_aru_obs",
  #  "T_aru_sim",
 #   "TvegObs",
 #   "TvegSim"#,
  #  "Dobs",
  #  "Dsim"
)

# ARU + SCORE
monitored_AS <- c(
  "mean_psi",
#  "beta0",
 # "beta1",
  "p11",
  "p_aru11",
  "p_aru01",
  "mu",
  "sigma",
  "psi",
 # "psi.pred.burn",
  "NOcc", "PropOcc",
  #  "T_pc_obs",
  #  "T_pc_sim",
    "T_aru_obs",
    "T_aru_sim",
  #  "TvegObs",
  #  "TvegSim",
    "Dobs",
    "Dsim", "z"
)

# HUMAN + ARU
monitored_HA <- c(
  "mean_psi",
 # "beta0",
 # "beta1",
   "p11",
  "p_aru11",
  "p_aru01",
 # "mu",
 # "sigma",
  "psi",
#  "psi.pred.burn",
 "NOcc", "PropOcc",
    "T_pc_obs",
    "T_pc_sim",
    "T_aru_obs",
    "T_aru_sim", "z"
  #  "TvegObs",
 #   "TvegSim"#,
  #  "Dobs",
  #  "Dsim"
)

# HUMAN + SCORES
monitored_HS <- c(
  "mean_psi",
 # "beta0",
 # "beta1",
  "p11",
 # "p_aru11",
 # "p_aru01",
  "mu",
  "sigma",
  "psi",
#  "psi.pred.burn",
 "NOcc", "PropOcc",
    "T_pc_obs",
    "T_pc_sim",
  #  "T_aru_obs",
  #  "T_aru_sim",
#    "TvegObs",
#    "TvegSim",
    "Dobs",
    "Dsim", "z"
)

# ARU ONLY
monitored_A <- c(
  "mean_psi",
 # "beta0",
 # "beta1",
  #"p11",
   "p_aru11",
   "p_aru01",
  #"mu",
  #"sigma",
  "psi",
 # "psi.pred.burn",
  "NOcc", "PropOcc", "z"#,
  #"T_pc_obs",
  #"T_pc_sim",
  #  "T_aru_obs",
  #  "T_aru_sim",
#  "TvegObs",
#  "TvegSim"
  #"Dobs",
  #"Dsim"
)

# SCORES ONLY
monitored_S <- c(
  "mean_psi",
 # "beta0",
 # "beta1",
  # "p11",
  # "p_aru11",
  # "p_aru01",
  "mu",
  "sigma",
  "psi",
#  "psi.pred.burn",
  "NOcc", "PropOcc",
  #"T_pc_obs",
  #"T_pc_sim",
  #  "T_aru_obs",
  #  "T_aru_sim",
 # "TvegObs",
 # "TvegSim",
  "Dobs",
  "Dsim", "z"
)

monitorList <- list(monitored_A, monitored_AS, monitored_H, monitored_HA, monitored_HAS, monitored_HS, monitored_S)

# Run JAGS, Save Outputs --------------------------------------------------
kos <- c("model_A", "model_AS", "model_H", "model_HA", "model_HAS", "model_HS", "model_S")

# make a list of jags objects to call each model
# ifelse?
modDir <- "COMB_minimal/models/jags_files/no_covars/"

kos

# currently the species dataset loop doesn't work
# need to manually specify species
sp <- "HAWO"
data <- data_HAWO

  # any model containing scores needs:
  gst <- sample(1:2, data$nsamples, replace = TRUE)
  gst[data$score > 0.0] <- 1
  gst[data$score <= 0.0] <- 2

  inits <- function() {
    list(
      mu = c(-2, 0),
      sigma = c(1, 1),
      z = zst,
      #   psi = psit,
      #   p11 = runif(1, 0.2, 0.8),
      #  lambda = runif(1, 1, 20), omega = runif(1, 0, 15),
      g = gst,
  #    beta0 = 0,
  #    beta1 = 0,
      #   beta2=0,
      reg_parm = 1
    )
  }

for (k in 1:length(kos)) {
  model <- ifelse(kos[k]=="model_A", paste0(modDir, "model_A.txt"),
                  ifelse(kos[k]=="model_AS", paste0(modDir, "model_AS.txt"),
                         ifelse(kos[k]=="model_H", paste0(modDir,"model_H.txt"),
                                ifelse(kos[k]=="model_HA", paste0(modDir,"model_HA.txt"), 
                                       ifelse(kos[k]=="model_HAS", paste0(modDir,"model_HAS.txt"), 
                                              ifelse(kos[k]=="model_HS", paste0(modDir, "model_HS.txt"), paste0(modDir,"model_S.txt")))))) )
  
monitored <- "z"
  
jagsResult <- jags(
  data=data,
  inits = inits,
  parameters.to.save = monitored,
  model.file = model,
  n.adapt = na,
  n.chains = nc,
  n.thin = nt,
  n.iter = ni,
  n.burnin = nb,
  parallel = TRUE,
)
assign(paste("jagsResult_no_covars", monitored, sp, year, str_extract(model, "(?<=/)[^/]+(?=\\.)"), Sys.Date(), sep="_"), jagsResult)
} # end of KO loop
  
  
  
#} # end of species dataset loop
rm(jagsResult)
  jags_objects <- ls(pattern = "^jagsResult")
  
  # Loop through the list and save each object as an RData file
  for (obj in jags_objects) {
    save(list = obj, file = paste0(obj, ".RData"))
  }

  