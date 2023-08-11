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
speciesCode <- "HAWO" # must match prefiltering of dataML_model.csv
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


# JAGS specification ------------------------------------------------------
modelFile <- tempfile()
cat(file = modelFile, "
model {

  # Priors
  
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  lam ~ dunif(0,20) # lambda: rate of target-species calls detected per file
  ome ~ dunif(0,10) # omega: rate of non-target detections per file

  # Parameters of the observation model for the scores
  mu[1] ~ dunif(-1, 1.5)
  mu[2] ~ dunif(-3, 0)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1: detection data and ARU counts
  for (i in 1:nsites) { # Loop over sites
    psi[i] ~ dunif(0, 1) # psi = Pr(Occupancy)
    z[i] ~ dbern(psi[i]) # Latent occupancy states
    
    # Point count
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }

    # ARU
    for(j in 1:nsurveys.aru) {
      y.aru[i,j] ~ dpois(lam*z[i] + ome)  # Total samples processed
    }
    # y.aru.sum[i] <- sum(y.aru[i,:])
    
    # GOF assessment - Point Count
    Tobs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
    ySim[i] ~ dbin(p11 * z[i], n_v[i])
    Tsim0[i] <- (sqrt(ySim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data
    
    # GOF assessment - ARU Deviance
    LLobs0[i] <- logdensity.pois(y_aru_sum[i], (lam*z[i] + ome)*n_s[i])
    y_aru_Sim[i] ~ dpois((lam*z[i] + ome)*n_s[i])
    LLsim0[i] <- logdensity.pois(y_aru_Sim[i], (lam*z[i] + ome)*n_s[i])
     # T_aru_obs0[i] <- (sqrt(y_aru_sum[i]) - sqrt((lam*z[i] + ome)*n_s[i]))^2  # FT discrepancy for observed data
    #  y_aru_Sim[i] ~ dpois((lam*z[i] + ome)*n_s[i])
    #  T_aru_sim0[i] <- (sqrt(y_aru_Sim[i]) - sqrt((lam*z[i] + ome)*n_s[i]))^2  # ...and for simulated data
  }

  # Likelihood part 2: feature score data
 # for (i in 1:nsites) {
#   site.prob[i] <- lam*z[i]/(lam*z[i]+ome) # Pr(sample is target species)
 # }
 # for(k in 1:nsamples) {
    # Sample specific covariate
 #   score[k] ~ dnorm(mu[g[k]], tau[g[k]]) # parameters are group specific
 #   probs[k,1] <- site.prob[siteid[k]]
 #   probs[k,2] <- 1 - site.prob[siteid[k]] # the prior class probabilities
 #   g[k] ~ dcat(probs[k,])
 #   N1[k] <- ifelse(g[k]==1, 1, 0)
 # }

  # GOF assessment
  T_pc_obs <- sum(Tobs0)
  T_pc_sim <- sum(Tsim0)
  Dobs <- -2 * sum(LLobs0)
  Dsim <- -2 * sum(LLsim0)
 #  T_aru_obs <- sum(T_aru_obs0)
 #  T_aru_sim <- sum(T_aru_sim0)
  # Derived quantities
 # Npos <- sum(N1[])
}
")

# initialization
zst <- rep(1, data$nsites)
psit = runif(data$nsites)
gst <- sample(1:2, data$nsamples, replace = TRUE)
gst[data$score > 0.0] <- 1
gst[data$score <= 0.0] <- 2
inits <- function() {
  list(
    mu = c(1, -1), sigma = c(1, 1), z = zst,
    psi = psit, p11 = runif(1, 0.2, 0.8),
    lam = runif(1, 1, 10), ome = runif(1, 0, 10), 
    g = gst
  )
}


# JAGS execution ----------------------------------------------------------

monitored <- c("psi", "p11", "lam", "ome", "mu", "sigma", "T_pc_obs", "T_pc_sim", "T_aru_obs", "Dobs", "Dsim")

# MCMC settings
na <- 1000
ni <- 4000
nt <- 1
nb <- 1000
nc <- 6


n_v_per_site <- rowSums(!is.na(data$y.ind[, 1:3])) # number of visits
y_pc_sum <- rowSums(data$y.ind[, 1:3], na.rm = TRUE)
n_samp_per_site <- rowSums(!is.na(data$y.aru[, 1:24])) # number of samples
y_aru_sum <- rowSums(data$y.aru[, 1:24], na.rm = TRUE)


data2 <- append(data, list(n_v = n_v_per_site, y_pc_sum = y_pc_sum, n_s = n_samp_per_site, y_aru_sum = y_aru_sum))
# TODO(matth79): JAGS does not like indices in the list, since it's non-numeric.
# Discuss team preferences on whether to omit it, nest the return value of
# readCombined, or some other alternative.
jagsData <- within(data2, rm(indices))

set.seed(123)

jagsResult <- jags(jagsData, inits, monitored, modelFile,
  n.adapt = na,
  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE,
)

pp.check(jagsResult, "T_pc_obs", "T_pc_sim", main="Posterior Predictive Check\nFreeman-Tukey discrepancy - Point Count")

pp.check(jagsResult, "Dobs", "Dsim", main="Posterior Predictive Check\nDeviance - ARU")
# mean(jagsResult$sims.list$Tsim >jagsResult$sims.list$Tobs)


MASS::eqscplot(
  y_aru_sum,
  jagsResult$sims.list$y_aru_Sim,
  xlim = range(0, 500),
  ylim = range(0, 500),
  xlab = "Observed data",
  ylab = "Simulated data"
)

mean(jagsResult$sims.list$Dobs >jagsResult$sims.list$Dsim, na.rm = TRUE)