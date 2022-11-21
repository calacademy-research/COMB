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
source(here("models/src/model_read_lib.R"))


singleSpeciesCombined <- function(params) {
  # parameters --------------------------------------------------------------
  speciesCode <- params$speciesCode
  year <- params$year
  threshold <- 0.5
  aruVisitLimit <- 24 # only consider this many ARU visits per site (ordered)
  PCVisitLimit <- 3 # only consider this many PC visits per site (ordered)
  
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
    psi ~ dunif(0, 1) # psi = Pr(Occupancy)
    p10 ~ dunif(0, 1) # p10 = Pr(y = 1 | z = 0)
    p11 ~ dunif(0, 1) # p11 = Pr(y = 1 | z = 1)
    lam ~ dunif(0, 1000) # lambda: rate of target-species calls detected
    ome ~ dunif(0, 1000) # omega: rate of non-target detections
    # Parameters of the observation model for the scores
    mu[1] ~ dnorm(0, 0.01)
    mu[2] ~ dnorm(0, 0.01)
    sigma ~ dunif(0, 10)
    tau <- 1 / (sigma * sigma)
    # Likelihood part 1: detection data and ARU counts
    for (i in 1:nsites) { # Loop over sites
      z[i] ~ dbern(psi) # Latent occupancy states

      p[i] <- z[i]*p11 + (1-z[i])*p10 # Detection probability including over-detections

      for(j in 1:nsurveys.pc) { # Loop over occasions
        y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
      }

      site.prob[i] <- lam*z[i]/(lam*z[i]+ome) # Pr(sample is target species) ... probability that a logit score came from a true detection

      for(j in 1:nsurveys.aru) { # Loop over occasions
        y.aru[i,j] ~ dpois(lam*z[i] + ome)  # Total samples processed
      }
    }
    # Likelihood part 2: feature score data
    for(k in 1:nsamples) {
      # Sample specific covariate
      score[k] ~ dnorm(mu[g[k]], tau) # parameters are group specific
      probs[k,1] <- site.prob[siteid[k]]
      probs[k,2] <- 1 - site.prob[siteid[k]] # the prior class probabilities
      g[k] ~ dcat(probs[k,])
    }
  }
  ")
  
  # initialization
  
  zst <- rep(1, data$nsites)
  gst <- sample(1:2, data$nsamples, replace = TRUE)
  gst[data$score > threshold] <- 1
  gst[data$score <= threshold] <- 2
  inits <- function() {
    list(
      mu = c(1, -1), sigma = 0.2, z = zst,
      psi = runif(1), p10 = runif(1, 0, 0.05), p11 = runif(1, 0.5, 0.8),
      lam = runif(1, 1, 2), ome = runif(1, 0, 0.4), g = gst
    )
  }
  
  
  # JAGS execution ----------------------------------------------------------
  
  monitored <- c("psi", "p10", "p11", "lam", "ome", "mu", "sigma")
  
  # MCMC settings
  na <- 1000
  ni <- 2000
  nt <- 1
  nb <- 1000
  nc <- 2
  
  jagsData <- within(data, rm(indices))
  jags(jagsData, inits, monitored, modelFile,
       n.adapt = na,
       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE,
       seed = 123
  )
}