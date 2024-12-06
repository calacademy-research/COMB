# libraries ---------------------------------------------------------------
library(googledrive)
library(here)
library(jagsUI)
library(lubridate)
library(tidyverse)
library(chron)

source(here("comb_functions.R"))
source(here("models/src/model_read_lib_agg.R")) # Functions were modified to read in ARU data that was pre-aggregated


# parameters --------------------------------------------------------------
speciesCode <- "BBWO" # must match prefiltering of dataML_model.csv
year <- 2020
threshold <- 0
aruVisitLimit <- 24 # only consider this many ARU visits per site (ordered). 32=4 mornings


# JAGS structuring --------------------------------------------------------
data <- readCombined(
  species = c(speciesCode),
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
  scale_datetime = F
)


data$y.pc # performs a quick check to see if pointcount data were read in properly
data$y.aru # aru detections
data$score

# names(data)
# data$aru_DateTime
# data$indices$visit.aru

# JAGS specification ------------------------------------------------------
modelFile <- tempfile()
cat(
  file = modelFile,
  "
model {

  # Priors
  #p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1) # commenting out p11 bc we are providing parameters for p below 
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)
  
  alpha0 ~ dunif(-5, 5)
  alpha1 ~ dunif(-5, 5)
  alpha2 ~ dunif(-5,5)
  alpha3 ~ dunif(-5, 5)
  alpha4 ~ dunif(-5,5)

  beta0 ~ dunif(-5, 5) # Occupancy intercept on prob. scale
  beta1 ~ dunif(-5, 5)
#  beta2 ~ dunif(-5, 5)

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, .1)
  mu[2] ~ dnorm(-2, .1)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1 & 2: PC and ARU detections
  
  for (i in 1:nsites) { # Loop over sites
    logit(psi[i]) <- beta0 + beta1*burn[i] 
    z[i] ~ dbern(psi[i]) # Latent occupancy states
    
    #GOF - Regression: Simulations
    psiSim[i] <- exp(beta0 + beta1*burn[i] )/(1+exp(beta0 + beta1*burn[i] ))
    zSim[i] ~ dbern(psiSim[i])

    # Point count detection process
   
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p11[i,j]*z[i]) # Observed occ. data (if available)
      logit(p11[i,j]) <- alpha0 + alpha1*Time[i,j] + alpha2*Time[i,j]*Time[i,j] + alpha3*Date[i,j] + alpha4*Date[i,j]*Date[i,j]
    
    # GOF Point Count - Tukey-Freeman Discrepancy
    expected_pc[i,j] = p11[i,j]*z[i] # in j loop because y.ind is parameterized
    y_pc_Sim[i,j] ~ dbern(p11[i,j] * z[i])
    
    }
    expected_sum_pc[i] = sum(expected_pc[i, ])
    T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(expected_sum_pc[i]))^2 
    y_pc_Sim_sum[i] = sum(y_pc_Sim[i,])
    T_pc_sim0[i] <- (sqrt(y_pc_Sim_sum[i]) - sqrt(expected_sum_pc[i]))^2  # ...and for simulated data
    
 # ARU - binomial
    p_aru[i] <- z[i]*p_aru11 + p_aru01 # Detection probability with false positives
    for(j in 1:nsurveys.aru) {
      y.aru[i,j] ~ dbern(p_aru[i])  # n_surveys.aru = Total files processed
      }

    # GOF ARU Count - Tukey-Freeman Discrepancy
    T_aru_obs0[i] <- (sqrt(y_aru_sum[i]) - sqrt(p_aru[i]*n_s[i]))^2  # FT discrepancy for observed data
    y_aru_Sim[i] ~ dbin(p_aru[i], n_s[i])
    T_aru_sim0[i] <- (sqrt(y_aru_Sim[i]) - sqrt(p_aru[i]*n_s[i]))^2  # ...and for simulated data
  }

  # Likelihood part 2: feature score data
  for(k in 1:nsamples) {
    score[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
    #GOF for ARU scores
      LLobs_score[k] <- logdensity.norm(score[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
      score_Sim[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
      LLsim_score[k] <- logdensity.norm(score_Sim[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
  }

  # GOF assessment
  T_pc_obs <- sum(T_pc_obs0)
  T_pc_sim <- sum(T_pc_sim0)
  T_aru_obs <- sum(T_aru_obs0)
  T_aru_sim <- sum(T_aru_sim0)

  #GOF assessment for scores
  Dobs <- -2 * sum(LLobs_score)
  Dsim <- -2 * sum(LLsim_score)

  #GOF - Regression: Difference in mean vegetation
  cz1 <- sum(z)
  cz0 <- sum(1 - z)
  TvegObs <- sum(burn*z) / ifelse(cz1>0,cz1, 1)  - sum(burn*(1-z)) / ifelse(cz0>0,cz0, 1)
  czSim1 <- sum(zSim)
  czSim0 <- sum(1 - zSim)
  TvegSim <- sum(burn*zSim) / ifelse(czSim1>0,czSim1, 1) - sum(burn*(1-zSim))/ifelse(czSim0>0,czSim0, 1)

  mean_psi <- mean(psi)
  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  PropOcc <- NOcc/nsites

  # simulate psi over a range of data
  for(k in 1:100) {
    logit(psi.pred.burn[k]) <- beta0 + beta1 * Xburn[k] # psi predictions for burn over the range of burn
  }

}
")


# inits for FULL MODEL ----------------------------------------------------


zst <- rep(1, data$nsites)
psit = runif(data$nsamples)
gst <- sample(1:2, data$nsamples, replace = TRUE)
gst[data$score > 0.0] <- 1
gst[data$score <= 0.0] <- 2
inits <- function() {
  list(
    mu = c(-2, 0),
    sigma = c(1, 1),
    z = zst,
    # psi = psit,
    #p11 = runif(1, 0.2, 0.8),
    #  lambda = runif(1, 1, 20), omega = runif(1, 0, 15),
    g = gst,
    beta0 = 0,
    beta1 = 0,
    #   beta2=0,
    reg_parm = 1
  )
}

monitored <- c(
  "mean_psi",
  "beta0",
  "beta1",
  #  "beta2",
  "alpha0",
  "alpha1",
  "alpha2",
  "alpha3",
  "alpha4",
  # "p11",
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


jagsData <- within(data, rm('indices','pc_DateTime','pc_YDay','pc_Time',"aru_DateTime","aru_YDay", "aru_Time","scores_datetime"))
#TODO: change over the date computations to pull directly from 'data' object

# Add site and visit covariates -------------------------------------------

# SITE covars
site_covars <- read_csv("./models/input/wide4havars.csv") %>%
  mutate(Point = avian_point) %>%
  filter(Point %in% data$indices$point$Point)

colnames(site_covars) # list of options for site-level covariates-- I chose measures of canopy cover and burn severity at the 4ha scale

covs <- site_covars %>% dplyr::select(Point, mean_CanopyCover_2020_4ha, mean_RAVGcbi_20182019_4ha, Perc_LTg22mHt_2020_4ha) %>% arrange(Point)

Veg <-left_join(as.data.frame(data2$indices$point), covs, by = "Point") 

jagsData$cover <- as.numeric(scale(Veg$mean_CanopyCover_2020_4ha))
jagsData$burn <- Veg$mean_RAVGcbi_20182019_4ha

summary(jagsData$burn)
jagsData$Xburn = seq(0, 2.5, length.out=100)


PC <- read_csv("./models/input/PC_delinted.csv")
PC$DateTime <- with_tz(PC$DateTime, tzone = "America/Los_Angeles")
PC$JDay <- as.numeric(format(PC$DateTime, "%j"))
PC$Time <- times(format(PC$DateTime, "%H:%M:%S"))
PC$JDay.s <- scale(PC$JDay)
PC$Time.s <- scale(PC$Time)

visit_pc <- PC %>% filter(year==2021, birdCode_fk==speciesCode) %>% dplyr::select(year, point_ID_fk, visit, JDay, Time, JDay.s, Time.s)
pcTime <- as.matrix(visit_pc %>% group_by(point_ID_fk, visit) %>% dplyr::select(point_ID_fk, visit, Time.s) %>% pivot_wider(names_from = visit, values_from = Time.s))

pcDate <- as.matrix(visit_pc %>% group_by(point_ID_fk, visit) %>% dplyr::select(point_ID_fk, visit, JDay.s) %>% pivot_wider(names_from = visit, values_from = JDay.s))

jagsData$Time <- pcTime[,2:4]
jagsData$Time[is.na(jagsData$Time)] <- 0

jagsData$Date <- pcDate[,2:4]
jagsData$Date[is.na(jagsData$Date)] <- 0

set.seed(123)

# MCMC settings
na <- 10
ni <- 80
nt <- 1
nb <- 10
nc <- 6


jagsResult_CONI20 <- jags(
  jagsData,
  inits,
  monitored,
  modelFile,
  n.adapt = na,
  n.chains = nc,
  n.thin = nt,
  n.iter = ni,
  n.burnin = nb,
  parallel = TRUE,
)
