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
source(here("models/src/model_read_lib_slice.R"))


ModelTrial <- function(params) {
  # parameters --------------------------------------------------------------
  speciesCode <- params$speciesCode
  year <- params$year
  threshold <- -100
  aruVisitLimit <- params$nARU # only consider this many ARU visits per site
  PCVisitLimit <- params$nPC
  
  # data --------------------------------------------------------------------
  #drive_auth(email = TRUE) # do not prompt when only one email has token
  #drive_sync(
  #  here("acoustic/data_ingest/output/"),
  #  "https://drive.google.com/drive/folders/1eOrXsDmiIW9YqJWrlUWR9-Cgc7hHKD_5"
  #)
  
  
  # JAGS structuring --------------------------------------------------------
  data <- readCombined(
    species = c(speciesCode),
    years = c(year),
    beginTime = dhours(6),
    endTime = dhours(10),
    visitLimit = aruVisitLimit,
    PCvisitlimit = PCVisitLimit,
    visitAggregation = "file",
    thresholdOptions = list(
      value = threshold,
      is.quantile = F
    ),
    squeeze = T
  )
  
  # JAGS specification ------------------------------------------------------
  modelFile <- tempfile()
  cat(file = modelFile, "
  model {
  
    # Priors
    p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
    beta0 ~ dnorm(0, 10) # Intercept for occupancy logistic regression
    beta1 ~ ddexp(0, sqrt(2.0)) # Slope for occupancy logistic regression with a Laplace prior for L1 regularization
  
    # Parameters of the observation model for the scores
    mu[1] ~ dnorm(-2, 3)
    mu[2] ~ dnorm(-2, 3)
    sigma[1] ~ dunif(0.1, 5)
    tau[1] <- 1 / (sigma[1] * sigma[1])
    sigma[2] ~ dunif(0.1, 5)
    tau[2] <- 1 / (sigma[2] * sigma[2])
  
    # Likelihood part 1: detection data and ARU counts
    for (i in 1:nsites) { # Loop over sites
      logit(psi[i]) <- beta0 + beta1*veg[i]
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
      psiSim[i] <- exp(beta0 + beta1*veg[i])/(1+exp(beta0 + beta1*veg[i]))
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
    TvegObs <- sum(veg*z) / ifelse(cz1>0,cz1, 1)  - sum(veg*(1-z)) / ifelse(cz0>0,cz0, 1)
    czSim1 <- sum(zSim)
    czSim0 <- sum(1 - zSim)
    TvegSim <- sum(veg*zSim) / ifelse(czSim1>0,czSim1, 1) - sum(veg*(1-zSim))/ifelse(czSim0>0,czSim0, 1)
  
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
  
  monitored <- c("beta0", "beta1", "p11", "T_pc_obs", "T_pc_sim", "Dobs", "Dsim", 
                 "mu", "sigma", "TvegObs", "TvegSim" , "mean_psi")
  
  # MCMC settings
  na <- 1000
  ni <- 6000
  nt <- 1
  nb <- 1000
  nc <- 6
  
  
  n_v_per_site <- rowSums(!is.na(data$y.ind[, 1:PCVisitLimit, drop=F])) # number of visits
  y_pc_sum <- rowSums(data$y.ind[, 1:PCVisitLimit, drop=F], na.rm = TRUE)
  n_samp_per_site <- rowSums(!is.na(data$y.aru[, 1:aruVisitLimit, drop=F])) # number of samples
  y_aru_sum <- rowSums(data$y.aru[, 1:aruVisitLimit, drop=F], na.rm = TRUE)
  
  
  data2 <- append(data, list(n_v = n_v_per_site, y_pc_sum = y_pc_sum, n_s = n_samp_per_site, y_aru_sum = y_aru_sum))
  
  jagsData <- within(data2, rm(indices))
  
  covs <- read_csv("./models/input/wide4havars.csv") %>%
    mutate(Point = avian_point) %>%
    dplyr::select(Point, mean_CanopyCover_2020_4ha)
  
  Veg <- left_join(as.data.frame(data$indices$point), covs, by = "Point")
  jagsData$veg <- as.numeric(scale(Veg$mean_CanopyCover_2020_4ha))
  
  set.seed(123)
  
  jags(jagsData, inits, monitored, modelFile,
                     n.adapt = na,
                     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE,
  )
}
