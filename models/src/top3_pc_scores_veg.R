# combined.R: Script to fit a (point count)+(ARU) occupancy model
#
# Usage:
#   * Run interactively
#   * or ```   source('models/src/combined.R')   ```


# libraries ---------------------------------------------------------------
library(googledrive)
library(here)
library(jagsUI)
library(lubridate)
library(tidyverse)

source(here("comb_functions.R"))
source(here("models/src/model_read_lib_top3.R"))


# parameters --------------------------------------------------------------
speciesCode <- "BBWO" # must match prefiltering of dataML_model.csv
year <- 2021
threshold <- -100
aruVisitLimit <- 24 # only consider this many ARU visits per site (ordered)


# JAGS structuring --------------------------------------------------------
data <- readCombined(
  species = c(speciesCode),
  years = c(year),
  beginTime = dhours(6),
  endTime = dhours(10),
  visitLimit = aruVisitLimit,
  visitAggregation = "file",
  thresholdOptions = list(
    value = threshold,
    is.quantile = F
  ),
  squeeze = T
)

data$y.pc[,1,,] # performs a quick check to see if pointcount data were read in properly
data$y.aru # aru detections
data$score

# JAGS specification ------------------------------------------------------
modelFile <- tempfile()
cat(file = modelFile, "
model {

  # Priors
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  beta0 ~ dnorm(0, 10) # Intercept for occupancy logistic regression
  beta1 ~ dnorm(0, 10) 
  beta2 ~ dnorm(0, 10) 

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, 3)
  mu[2] ~ dnorm(-2, 3)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1: detection data and ARU counts
  for (i in 1:nsites) { # Loop over sites
    logit(psi[i]) <- beta0 + beta1*cover[i] + beta2*burn[i]
    z[i] ~ dbern(psi[i]) # Latent occupancy states

    # Point count
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }
    # GOF Point Count - Tukey-Freeman Discrepancy
    Tobs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
    ySim[i] ~ dbin(p11 * z[i], n_v[i])
    Tsim0[i] <- (sqrt(ySim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data

    #GOF - Regression: Simulations
    psiSim[i] <- exp(beta0 + beta1*cover[i] + beta2*burn[i])/(1+exp(beta0 + beta1*cover[i] + beta2*burn[i]))
    zSim[i] ~ dbern(psiSim[i])

  }

  # Gaussian mixture model for ML ARU scores
  for(k in 1:nsamples) {
    score[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
    #GOF for ARU scores
      LLobs_score[k] <- logdensity.norm(score[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
      score_Sim[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
      LLsim_score[k] <- logdensity.norm(score_Sim[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
  }


  # GOF assessment for Point Count data
  T_pc_obs <- sum(Tobs0)
  T_pc_sim <- sum(Tsim0)

  #GOF assessment for scores
  Dobs <- -2 * sum(LLobs_score)
  Dsim <- -2 * sum(LLsim_score)

  #GOF - Regression: Difference in mean vegetation
  cz1 <- sum(z)
  cz0 <- sum(1 - z)
  TcoverObs <- sum(cover*z) / ifelse(cz1>0,cz1, 1)  - sum(cover*(1-z)) / ifelse(cz0>0,cz0, 1)
  czSim1 <- sum(zSim)
  czSim0 <- sum(1 - zSim)
  TcoverSim <- sum(cover*zSim) / ifelse(czSim1>0,czSim1, 1) - sum(cover*(1-zSim))/ifelse(czSim0>0,czSim0, 1)

  mean_psi = mean(psi) # Mean occupancy across sites
}
")

# initialization
zst <- rep(1, data$nsites)
psit = runif(data$nsamples)
gst <- sample(1:2, data$nsamples, replace = TRUE)
gst[data$score > 0.0] <- 1
gst[data$score <= 0.0] <- 2
inits <- function() {
  list(
    mu = c(1, -1), sigma = c(1, 1), z = zst,
   # psi = psit,
    p11 = runif(1, 0.2, 0.8),
    lambda = runif(1, 1, 20), omega = runif(1, 0, 15),
    g = gst, beta0 = 0, beta1 = 0,
    reg_parm = 1
  )
}


# JAGS execution ----------------------------------------------------------

monitored <- c("beta0", "beta1", "beta2", "p11", 
               "T_pc_obs", "T_pc_sim", "Dobs", "Dsim","mu", "sigma",
               "TcoverObs", "TcoverSim" , "mean_psi")

# MCMC settings
na <- 1000
ni <- 6000
nt <- 1
nb <- 1000
nc <- 6


n_v_per_site <- rowSums(!is.na(data$y.ind[, 1:3])) # number of point-count visits
y_pc_sum <- rowSums(data$y.ind[, 1:3], na.rm = TRUE)
n_samp_per_site <- rowSums(!is.na(data$y.aru[, 1:24])) # number of samples
y_aru_sum <- rowSums(data$y.aru[, 1:24], na.rm = TRUE)


data2 <- append(data, list(n_v = n_v_per_site, y_pc_sum = y_pc_sum, n_s = n_samp_per_site, y_aru_sum = y_aru_sum))

jagsData <- within(data2, rm(indices))

site_covars <- read_csv("./models/input/wide4havars.csv") %>%
  mutate(Point = avian_point) %>%
  filter(Point %in% data$indices$point$Point)

colnames(site_covars) # list of options for site-level covariates-- I chose measures of canopy cover and burn severity at the 4ha scale
  
covs <- site_covars %>% dplyr::select(Point, mean_CanopyCover_2020_4ha, mean_RAVGcbi4_20202021_4ha)

Veg <-left_join(as.data.frame(data$indices$point), covs, by = "Point")
jagsData$cover <- as.numeric(scale(Veg$mean_CanopyCover_2020_4ha))
jagsData$burn <- as.numeric(scale(Veg$mean_RAVGcbi4_20202021_4ha)) # check out a histogram of this... sort of bimodal 

set.seed(123)

jagsResult <- jags(jagsData, inits, monitored, modelFile,
  n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE,
  verbose=T
)

pp.check(jagsResult, "TcoverObs", "TcoverSim", main=
           paste("Posterior Predictive Check\nLogistic Regression: Difference in Occupied Means\n","Species = ", speciesCode,"; Year = ", year, sep = " ")
)

pp.check(jagsResult, "T_pc_obs", "T_pc_sim", main=
           paste("Posterior Predictive Check\nFreeman-Tukey discrepancy - Point Count\n","Species = ", speciesCode,"; Year = ", year, sep = " ")
)

pp.check(jagsResult, "Dobs", "Dsim", main=
           paste("Posterior Predictive Check\nDeviance - ARU\n","Species = ", speciesCode,"; Year = ", year, sep = " ")
)
