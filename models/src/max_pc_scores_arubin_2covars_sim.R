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
source("models/src/simulation_with_covars.R")
# source(here("models/src/model_read_lib_agg.R")) # Functions were modified to read in ARU data that was pre-aggregated


ModelTrial <- function(params) {
  # parameters --------------------------------------------------------------
  aruVisitLimit <- params$nARU # only consider this many ARU visits per site (ordered)
  PCVisitLimit <- params$nPC
  p11 = params$p11
  p_aru11 = params$p_aru11
  p_aru01 = params$p_aru01
  beta0 <- params$beta0
  beta1 <- params$beta1
  mu <- c(-2.5, -1.5)
  sigma <- c(0.3, 1)
  n_points <- params$n_points
  

  data <- simulation(
    p11 = p11,
    p_aru11 = p_aru11,
    p_aru01 = p_aru01,
    mu = mu,
    sigma = sigma,
    n_visits = PCVisitLimit,
    n_recordings = aruVisitLimit,
    n_points = n_points, 
    beta0 = beta0, 
    beta1 = beta1
  )

  # JAGS specification ------------------------------------------------------
  modelFile <- tempfile()
  cat(
    file = modelFile,
    "
    model {
    
      # Priors
      p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
      p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
      p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)
      beta0 ~ dnorm(0, 0.10) # Intercept for occupancy logistic regression
      beta1 ~ dnorm(0, 0.10)
    
      # Parameters of the observation model for the scores
      mu[1] ~ dnorm(-2, 0.2)
      mu[2] ~ dnorm(-2, 0.2)
      sigma[1] ~ dunif(0.1, 5)
      tau[1] <- 1 / (sigma[1] * sigma[1])
      sigma[2] ~ dunif(0.1, 5)
      tau[2] <- 1 / (sigma[2] * sigma[2])
    
      # Likelihood part 1 & 2: PC and ARU detections
      for (i in 1:nsites) { # Loop over sites
        logit(psi[i]) <- beta0 + beta1*burn[i]
        z[i] ~ dbern(psi[i]) # Latent occupancy states
    
        # Point count
        p[i] <- z[i]*p11 # Detection probability
        for(j in 1:nsurveys.pc) {
          y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
        }
        # GOF Point Count - Tukey-Freeman Discrepancy
        T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
        y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
        T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data
    
        #GOF - Regression: Simulations
        psiSim[i] <- exp(beta0 + beta1*burn[i])/(1+exp(beta0 + beta1*burn[i]))
        zSim[i] ~ dbern(psiSim[i])
    
        # ARU - binomial
        p_aru[i] <- z[i]*p_aru11 + p_aru01 # Detection probability with false positives
        for(j in 1:nsurveys.aru) {
          y.aru[i,j] ~ dbern(p_aru[i])  # n_surveys.aru = Total files processed
        }
        # GOF ARU Count - Tukey-Freeman Discrepancy
        T_aru_obs0[i] <- (sqrt(y_aru_sum[i]) - sqrt(p_aru[i]*n_s[i]))^2  # FT discrepancy for observed data
        y_aru_Sim[i] ~ dbin(p_aru[i], n_v[i])
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
      czSim1 <- sum(zSim)
      czSim0 <- sum(1 - zSim)

      mean_psi <- mean(psi)
      NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
      PropOcc <- NOcc/nsites
    
    }
    "
  )

  # initialization
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
      # psi = psit,
      p11 = runif(1, 0.2, 0.8),
      #  lambda = runif(1, 1, 20), omega = runif(1, 0, 15),
      g = gst,
      beta0 = 0,
      beta1 = 0,
      reg_parm = 1
    )
  }


  # JAGS execution ----------------------------------------------------------

  monitored <- c(
    "beta0",
    "beta1",
    "p11",
    "p_aru11",
    "p_aru01",
    "mu",
    "sigma",
    "mean_psi",
    "NOcc", "PropOcc",
    "T_pc_obs",
    "T_pc_sim",
    "T_aru_obs",
    "T_aru_sim",
    "TvegObs",
    "TvegSim",
    "Dobs",
    "Dsim", 
    "psi", 
    "z"
  )

  # MCMC settings
  na <- 1000
  ni <- 8000
  nt <- 1
  nb <- 1000
  nc <- 6


  n_v_per_site <-
    rowSums(!is.na(data$y.ind[, 1:PCVisitLimit, drop = F])) # number of visits
  y_pc_sum <- rowSums(data$y.ind[, 1:PCVisitLimit, drop = F], na.rm = TRUE)
  n_samp_per_site <-
    rowSums(!is.na(data$y.aru[, 1:aruVisitLimit, drop = F])) # number of samples
  y_aru_sum <- rowSums(data$y.aru[, 1:aruVisitLimit, drop = F], na.rm = TRUE)


  jagsData <-
    append(
      data,
      list(
        n_v = n_v_per_site,
        y_pc_sum = y_pc_sum,
        n_s = n_samp_per_site,
        y_aru_sum = y_aru_sum
      )
    )
  
  # set.seed(123)

  jagsResult <- jags(
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

  jagsResult
}
