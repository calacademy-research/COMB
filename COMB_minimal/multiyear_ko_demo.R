library(googledrive)
library(here)
library(jagsUI)
library(lubridate)
library(tidyverse)
library(chron)

source(here("comb_functions.R"))
source(here("models/src/model_read_lib_agg.R")) # Functions were modified to read in ARU data that was pre-aggregated# run "knockout" for each species on the list

set.seed(123)


speciesCode <- sort(c("BBWO", "CONI")) # must match prefiltering of dataML_model.csv
# BBWO should start at 5 AM and end at 9 AM, start = 5, end = 10, samples per day = 10, 4 days max
# CONI most likely to be vocal at 4 AM, 6 PM to 11 PM
year <- c(2020,2021)
threshold <- 0
aruVisitLimit <- 32 # only consider this many ARU visits per site (ordered). 32=4 mornings

# JAGS structuring --------------------------------------------------------

for (sp in speciesCode) {
  begin_hour = ifelse(sp == "CONI", 4, 5)
  end_hour = ifelse(sp == "CONI",9,10)
  data <- readCombined(
    species = c(sp),
    years = c(year),
    beginTime = dhours(begin_hour),
    endTime = dhours(end_hour),
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

datasets <- list(BBWO=data_BBWO, CONI=data_CONI)

zst <- matrix(1, data$nyears, data$nsites)
#psit = matrix(runif(data$nyears*data$nsamples), data$nyears, data$nsites)

# any model containing scores needs:
gst <- sample(1:2, data_CONI$nsamples, replace = TRUE)
gst[data_CONI$score > 0.0] <- 1
gst[data_CONI$score <= 0.0] <- 2

inits <- function() {
  list(
    mu = c(-2, 0),
    sigma = c(1, 1),
    z = zst,
    #   psi = psit,
    #   p11 = runif(1, 0.2, 0.8),
    #  lambda = runif(1, 1, 20), omega = runif(1, 0, 15),
    g = gst,
    beta0 = 0,
    beta1 = 0,
    beta2 = 0,
    #   beta2=0,
    reg_parm = 1
  )
}

monitored_HAS <- c(
  "mean_psi",
  "beta0",
  "beta1",
  "beta2",
  "alpha0",
  "alpha1",
  "alpha2",
  "alpha3",
  "alpha4",
  "alpha5",
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
  "Dsim"
)

modDir <- "models/jags_files/"
model <- paste0(modDir,"model_HAS_multiyear.txt")
sp <- "CONI"
data <- data_CONI

jagsResult_my <- jags(
  data=data,
  inits = inits,
  parameters.to.save = monitored_HAS,
  model.file = model,
  n.adapt = na,
  n.chains = nc,
  n.thin = nt,
  n.iter = ni,
  n.burnin = nb,
  parallel = TRUE,
)