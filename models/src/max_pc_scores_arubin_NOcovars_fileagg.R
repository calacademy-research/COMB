# combined.R: Script to fit a (point count)+(ARU) occupancy model
#
# Usage:
#   * Run interactively
#   * or ```   source('models/src/combined.R')   ```

source('models/src/max_pc_scores_arubin_NOcovars_fileagg.R')

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
aruVisitLimit <- 32 # only consider this many ARU visits per site (ordered). 32=4 mornings


# JAGS structuring --------------------------------------------------------
data <- readCombined(
  species = c(speciesCode),
  years = c(year),
  beginTime = dhours(5),
  endTime = dhours(9),
  visitLimit = aruVisitLimit,
  visitAggregation = "file",
  thresholdOptions = list(value = threshold,
                          is.quantile = F),
  squeeze = T,
  logit_col = "max_logit" # This is specifying we want the max_logit column from aggregated data
)
data$y.pc[,1,,] # performs a quick check to see if pointcount data were read in properly
data$y.aru # aru detections
data$y.aru[data$y.aru>1] <- 1
data$score

dim(data$indices$visit.aru)
length(data$score)

# JAGS specification ------------------------------------------------------
modelFile <- tempfile()
cat(
  file = modelFile,
  "
model {

  # Priors
  psi ~ dunif(0, 1)
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, .1)
  mu[2] ~ dnorm(-2, .1)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1 & 2: PC and ARU detections
  
  for (i in 1:nsites) { # Loop over sites
  #  logit(psi[i]) <- beta0 + beta1*burn[i] 
    z[i] ~ dbern(psi) # Latent occupancy states

    # Point count detection process
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }

    # GOF Point Count - Tukey-Freeman Discrepancy
      T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
      y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
      T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data
  
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

  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  PropOcc <- NOcc/nsites

}
")


# inits for FULL MODEL ----------------------------------------------------


zst <- rep(1, data$nsites)
#psit = runif(data$nsamples)
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
 #   beta2=0,
    reg_parm = 1
  )
}

monitored <- c(
  "psi",
  "p11",
  "p_aru11",
  "p_aru01",
  "mu",
  "sigma",
  "NOcc", "PropOcc",
  "T_pc_obs",
  "T_pc_sim",
  "T_aru_obs",
  "T_aru_sim",
  "Dobs",
  "Dsim"
)

#monitored <- c("z")

# MCMC settings
na <- 1000
ni <- 8000
nt <- 1
nb <- 1000
nc <- 6


n_v_per_site <-
  rowSums(!is.na(data$y.ind[, 1:3])) # number of visits
y_pc_sum <- rowSums(data$y.ind[, 1:3], na.rm = TRUE)
n_samp_per_site <-
  rowSums(!is.na(data$y.aru[, 1:aruVisitLimit])) # number of samples
y_aru_sum <- rowSums(data$y.aru[, 1:aruVisitLimit], na.rm = TRUE)

data2 <-
  append(
    data,
    list(
      n_v = n_v_per_site,
      y_pc_sum = y_pc_sum,
      n_s = n_samp_per_site,
      y_aru_sum = y_aru_sum
    )
  )
# TODO(matth79): JAGS does not like indices in the list, since it's non-numeric.
# Discuss team preferences on whether to omit it, nest the return value of
# readCombined, or some other alternative.
jagsData <- within(data2, rm(indices))


# Add site and visit covariates -------------------------------------------
# SEE 'max_pc_scores_arubin_2covars_detcovars.R' for section on adding site and visit covars

set.seed(123)

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


# PPCs ------------------------------------------------------

pp.check(
  jagsResult,
  "T_pc_obs",
  "T_pc_sim",
  main =
    paste(
      "Posterior Predictive Check\nFreeman-Tukey discrepancy - Point Count\n",
      "Species = ",
      speciesCode,
      "; Year = ",
      year,
      sep = " "
    )
)

pp.check(
  jagsResult,
  "T_aru_obs",
  "T_aru_sim",
  main =
    paste(
      "Posterior Predictive Check\nFreeman-Tukey discrepancy - ARU\n",
      "Species = ",
      speciesCode,
      "; Year = ",
      year,
      sep = " "
    )
)
#FAIL

pp.check(
  jagsResult,
  "Dobs",
  "Dsim",
  main =
    paste(
      "Posterior Predictive Check\nDeviance - ARU\n",
      "Species = ",
      speciesCode,
      "; Year = ",
      year,
      sep = " "
    )
)

null0 <- jagsResult

# {PC ONLY MODEL} -----------------------------------------------------------
modelFile <- tempfile()
cat(
  file = modelFile,
  "
model {

  # Priors
  psi ~ dunif(0, 1)
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
#  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
#  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)

  # Parameters of the observation model for the scores
#  mu[1] ~ dnorm(-2, .1)
#  mu[2] ~ dnorm(-2, .1)
#  sigma[1] ~ dunif(0.1, 5)
#  tau[1] <- 1 / (sigma[1] * sigma[1])
#  sigma[2] ~ dunif(0.1, 5)
#  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1 & 2: PC and ARU detections
  
  for (i in 1:nsites) { # Loop over sites
  #  logit(psi[i]) <- beta0 + beta1*burn[i] 
    z[i] ~ dbern(psi) # Latent occupancy states

    # Point count detection process
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }

    # GOF Point Count - Tukey-Freeman Discrepancy
      T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
      y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
      T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data

  }

  # GOF assessment
  T_pc_obs <- sum(T_pc_obs0)
  T_pc_sim <- sum(T_pc_sim0)
  
    NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  PropOcc <- NOcc/nsites

}
")


# inits for PC ONLY model -------------------------------------------------

#zst <- rep(1, data$nsites)
#psit = runif(data$nsamples)
#gst <- sample(1:2, data$nsamples, replace = TRUE)
#gst[data$score > 0.0] <- 1
#gst[data$score <= 0.0] <- 2
inits <- function() {
  list(
  #  mu = c(-2, 0),
  #  sigma = c(1, 1),
    z = zst,
 #    psi = psit,
   p11 = runif(1, 0.2, 0.8),
    #  lambda = runif(1, 1, 20), omega = runif(1, 0, 15),
 #   g = gst,
  #  beta0 = 0,
 #   beta1 = 0,
    #   beta2=0,
    reg_parm = 1
  )
}

monitored <- c(
  "psi",
  "p11",
#  "p_aru11",
#  "p_aru01",
#  "mu",
#  "sigma",
   "NOcc", "PropOcc",
  "T_pc_obs",
  "T_pc_sim",
#  "T_aru_obs",
#  "T_aru_sim",
#  "Dobs",
#  "Dsim"
)

n_v_per_site <-
  rowSums(!is.na(data$y.ind[, 1:3])) # number of visits
y_pc_sum <- rowSums(data$y.ind[, 1:3], na.rm = TRUE)
#n_samp_per_site <-
#  rowSums(!is.na(data$y.aru[, 1:aruVisitLimit])) # number of samples
#y_aru_sum <- rowSums(data$y.aru[, 1:aruVisitLimit], na.rm = TRUE)

data2 <-
  append(
    data,
    list(
      n_v = n_v_per_site,
      y_pc_sum = y_pc_sum,
      n_s = n_samp_per_site,
      y_aru_sum = y_aru_sum
    )
  )
# TODO(matth79): JAGS does not like indices in the list, since it's non-numeric.
# Discuss team preferences on whether to omit it, nest the return value of
# readCombined, or some other alternative.
jagsData <- within(data2, rm(indices))


# Add site and visit covariates -------------------------------------------
# SEE 'max_pc_scores_arubin_2covars_detcovars.R' for section on adding site and visit covars

set.seed(123)

jagsResult.pc <- jags(
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

pp.check(
  jagsResult.pc,
  "T_pc_obs",
  "T_pc_sim",
  main =
    paste(
      "Posterior Predictive Check\nFreeman-Tukey discrepancy - Point Count\n",
      "Species = ",
      speciesCode,
      "; Year = ",
      year,
      sep = " "
    )
)

nullpc <- jagsResult.pc

# {ARU ONLY MODEL} -----------------------------------------------------------
modelFile <- tempfile()
cat(
  file = modelFile,
  "
model {

  # Priors
  psi ~ dunif(0, 1)
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)

  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, .1)
  mu[2] ~ dnorm(-2, .1)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1 & 2: PC and ARU detections
  
  for (i in 1:nsites) { # Loop over sites
  #  logit(psi[i]) <- beta0 + beta1*burn[i] 
    z[i] ~ dbern(psi) # Latent occupancy states

    # Point count detection process
  #  p[i] <- z[i]*p11 # Detection probability
  #  for(j in 1:nsurveys.pc) {
  #    y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
  #  }

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
#  T_pc_obs <- sum(T_pc_obs0)
#  T_pc_sim <- sum(T_pc_sim0)
  T_aru_obs <- sum(T_aru_obs0)
  T_aru_sim <- sum(T_aru_sim0)

  #GOF assessment for scores
  Dobs <- -2 * sum(LLobs_score)
  Dsim <- -2 * sum(LLsim_score)

 NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  PropOcc <- NOcc/nsites

}
")

n_v_per_site <-
  rowSums(!is.na(data$y.ind[, 1:3])) # number of visits
y_pc_sum <- rowSums(data$y.ind[, 1:3], na.rm = TRUE)
n_samp_per_site <-
  rowSums(!is.na(data$y.aru[, 1:aruVisitLimit])) # number of samples
y_aru_sum <- rowSums(data$y.aru[, 1:aruVisitLimit], na.rm = TRUE)

data2 <-
  append(
    data,
    list(
      n_v = n_v_per_site,
      y_pc_sum = y_pc_sum,
      n_s = n_samp_per_site,
      y_aru_sum = y_aru_sum
    )
  )
# TODO(matth79): JAGS does not like indices in the list, since it's non-numeric.
# Discuss team preferences on whether to omit it, nest the return value of
# readCombined, or some other alternative.
jagsData <- within(data2, rm(indices))

zst <- rep(1, data$nsites)
#psit = runif(data$nsamples)
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
  #  beta0 = 0,
 #   beta1 = 0,
    #   beta2=0,
    reg_parm = 1
  )
}

monitored <- c(
  "psi",
  #  "p11",
  "p_aru11",
  "p_aru01",
  "mu",
  "sigma",
   "NOcc", "PropOcc",
  #  "T_pc_obs",
  #  "T_pc_sim",
  "T_aru_obs",
  "T_aru_sim",
  "Dobs",
  "Dsim"
)

jagsResult.aru <- jags(
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


pp.check(
  jagsResult.aru,
  "T_aru_obs",
  "T_aru_sim",
  main =
    paste(
      "Posterior Predictive Check\nFreeman-Tukey discrepancy - ARU\n",
      "Species = ",
      speciesCode,
      "; Year = ",
      year,
      sep = " "
    )
)

pp.check(
  jagsResult.aru,
  "Dobs",
  "Dsim",
  main =
    paste(
      "Posterior Predictive Check\nDeviance - ARU\n",
      "Species = ",
      speciesCode,
      "; Year = ",
      year,
      sep = " "
    )
)

nullaru <- jagsResult.aru

# {ARU model - no scores} ---------------------------------------------------

modelFile <- tempfile()
cat(
  file = modelFile,
  "
model {

  # Priors
  psi ~ dunif(0, 1)
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1) 
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)

  # Likelihood part 1 & 2: PC and ARU detections

  for (i in 1:nsites) { # Loop over sites
  #  logit(psi[i]) <- beta0 + beta1*burn[i] 
    z[i] ~ dbern(psi) # Latent occupancy states

    # Point count detection process
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }
    

    # GOF Point Count - Tukey-Freeman Discrepancy
      T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
      y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
      T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data


    # ARU - binomial
      p_aru[i] <- z[i]*p_aru11 + p_aru01 # Detection probability with false positives
      for(j in 1:nsurveys.aru) {
        y.aru[i,j] ~ dbern(p_aru[i]) 
    }

    # GOF ARU Count - Tukey-Freeman Discrepancy
    T_aru_obs0[i] <- (sqrt(y_aru_sum[i]) - sqrt(p_aru[i]*n_s[i]))^2  # FT discrepancy for observed data
    y_aru_Sim[i] ~ dbin(p_aru[i], n_v[i])
    T_aru_sim0[i] <- (sqrt(y_aru_Sim[i]) - sqrt(p_aru[i]*n_s[i]))^2  # ...and for simulated data
  }

  # Likelihood part 2: feature score data
#  for(k in 1:nsamples) {
#    score[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
    #GOF for ARU scores
#      LLobs_score[k] <- logdensity.norm(score[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
#      score_Sim[k] ~ dnorm(mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
#      LLsim_score[k] <- logdensity.norm(score_Sim[k], mu[z[siteid[k]] + 1], tau[z[siteid[k]] + 1])
#  }

  # GOF assessment
  T_pc_obs <- sum(T_pc_obs0)
  T_pc_sim <- sum(T_pc_sim0)
 T_aru_obs <- sum(T_aru_obs0)
 T_aru_sim <- sum(T_aru_sim0)

  #GOF assessment for scores
 # Dobs <- -2 * sum(LLobs_score)
#  Dsim <- -2 * sum(LLsim_score)

  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  PropOcc <- NOcc/nsites

}
")


# inits for ARU ONLY model -------------------------------------------------

# initialization
zst <- rep(1, data$nsites)
#psit = runif(data$nsamples)
#gst <- sample(1:2, data$nsamples, replace = TRUE)
#gst[data$score > 0.0] <- 1
#gst[data$score <= 0.0] <- 2
inits <- function() {
  list(
#    mu = c(-2, 0),
#    sigma = c(1, 1),
    z = zst,
    #  psi = psit,
     p11 = runif(1, 0.2, 0.8),
#    lambda = runif(1, 1, 20), omega = runif(1, 0, 15),
  #  g = gst,
 #   beta0 = 0,
  #  beta1 = 0,
    #   beta2=0,
    reg_parm = 1
  )
}

monitored <- c(
  "psi",
  "p11",
  "p_aru11",
  "p_aru01",
  "NOcc", "PropOcc",
    "T_pc_obs",
    "T_pc_sim",
    "T_aru_obs",
    "T_aru_sim"
#  "Dobs",
#  "Dsim"
)

jagsResult.noscore <- jags(
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

pp.check(
  jagsResult.noscore,
  "T_pc_obs",
  "T_pc_sim",
  main =
    paste(
      "Posterior Predictive Check\nFreeman-Tukey discrepancy - Point Count\n",
      "Species = ",
      speciesCode,
      "; Year = ",
      year,
      sep = " "
    )
)

pp.check(
  jagsResult.noscore,
  "T_aru_obs",
  "T_aru_sim",
  main =
    paste(
      "Posterior Predictive Check\nFreeman-Tukey discrepancy - ARU\n",
      "Species = ",
      speciesCode,
      "; Year = ",
      year,
      sep = " "
    )
)

nullnoscore <- jagsResult.noscore

# {ARU model - scores, no detections }------------------------------------
modelFile <- tempfile()
cat(
  file = modelFile,
  "
model {

  # Priors
  psi ~ dunif(0,1)
  p11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru11 ~ dbeta(2, 2) # p11 = Pr(y = 1 | z = 1)
  p_aru01 ~ dbeta(1, 3)I(0, 1 - p_aru11) # p11 = Pr(y = 1 | z = 0)
  
  # Parameters of the observation model for the scores
  mu[1] ~ dnorm(-2, 5)
  mu[2] ~ dnorm(-2, 5)
  sigma[1] ~ dunif(0.1, 5)
  tau[1] <- 1 / (sigma[1] * sigma[1])
  sigma[2] ~ dunif(0.1, 5)
  tau[2] <- 1 / (sigma[2] * sigma[2])

  # Likelihood part 1 & 2: PC and ARU detections
  
  for (i in 1:nsites) { # Loop over sites
  #  logit(psi[i]) <- beta0 + beta1*burn[i] 
    z[i] ~ dbern(psi) # Latent occupancy states

    # Point count detection process
    p[i] <- z[i]*p11 # Detection probability
    for(j in 1:nsurveys.pc) {
      y.ind[i,j] ~ dbern(p[i]) # Observed occ. data (if available)
    }

    # GOF Point Count - Tukey-Freeman Discrepancy
      T_pc_obs0[i] <- (sqrt(y_pc_sum[i]) - sqrt(p11*z[i]*n_v[i]))^2  # FT discrepancy for observed data
      y_pc_Sim[i] ~ dbin(p11 * z[i], n_v[i])
      T_pc_sim0[i] <- (sqrt(y_pc_Sim[i]) - sqrt(p11*z[i]*n_v[i]))^2  # ...and for simulated data
  
 # ARU - binomial
    p_aru[i] <- z[i]*p_aru11 + p_aru01 # Detection probability with false positives
#    for(j in 1:nsurveys.aru) {
 #     y.aru[i,j] ~ dbern(p_aru[i])  # n_surveys.aru = Total files processed
  #    }

    # GOF ARU Count - Tukey-Freeman Discrepancy
#    T_aru_obs0[i] <- (sqrt(y_aru_sum[i]) - sqrt(p_aru[i]*n_s[i]))^2  # FT discrepancy for observed data
 #   y_aru_Sim[i] ~ dbin(p_aru[i], n_s[i])
 #   T_aru_sim0[i] <- (sqrt(y_aru_Sim[i]) - sqrt(p_aru[i]*n_s[i]))^2  # ...and for simulated data
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
#  T_aru_obs <- sum(T_aru_obs0)
#  T_aru_sim <- sum(T_aru_sim0)

  #GOF assessment for scores
  Dobs <- -2 * sum(LLobs_score)
  Dsim <- -2 * sum(LLsim_score)

  NOcc <- sum(z[]) # derived quantity for # of sites occupied (to compare with 'naive' sum(y.aru[]) and sum(y.pc[])
  PropOcc <- NOcc/nsites

}
")


# inits for ARU SCORE ONLY model --------------------------------------------

# initialization
zst <- rep(1, data$nsites)
#psit = runif(data$nsamples)
#gst <- sample(1:2, data$nsamples, replace = TRUE)
#gst[data$score > 0.0] <- 1
#gst[data$score <= 0.0] <- 2
inits <- function() {
  list(
        mu = c(-2, 0),
        sigma = c(1, 1),
    z = zst,
    #  psi = psit,
    p11 = runif(1, 0.2, 0.8),
    #    lambda = runif(1, 1, 20), omega = runif(1, 0, 15),
       g = gst,
 #   beta0 = 0,
  #  beta1 = 0,
    #   beta2=0,
    reg_parm = 1
  )
}
monitored <- c(
  "psi",
  "p11",
 # "p_aru11",
 # "p_aru01",
  "mu",
  "sigma",
  "NOcc", "PropOcc",
  "T_pc_obs",
  "T_pc_sim",
 # "T_aru_obs",
 # "T_aru_sim",
  "Dobs",
  "Dsim"
)

jagsResult.scoreonly <- jags(
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




pp.check(
  jagsResult.scoreonly,
  "T_pc_obs",
  "T_pc_sim",
  main =
    paste(
      "Posterior Predictive Check\nFreeman-Tukey discrepancy - Point Count\n",
      "Species = ",
      speciesCode,
      "; Year = ",
      year,
      sep = " "
    )
)


pp.check(
  jagsResult.scoreonly,
  "Dobs",
  "Dsim",
  main =
    paste(
      "Posterior Predictive Check\nDeviance - ARU\n",
      "Species = ",
      speciesCode,
      "; Year = ",
      year,
      sep = " "
    )
)

nullscore <- jagsResult.scoreonly



save(jagsData, null0, nullaru, nullpc, nullnoscore, nullscore, file = "BBWO2020_NULLmodel_knockout.RData")



# explore -----------------------------------------------------------------

latlong <- read.csv("models/input/latlong.csv", header = F, col.names = c("siteid", "lat", "long"))

fullsumm <- data.frame(null0$summary)
arusumm <- data.frame(nullaru$summary)
pcsumm <- data.frame(nullpc$summary)
noscoresumm <- data.frame(nullnoscore$summary)
scoreonlysumm <- data.frame(nullscore$summary)
#dayaggsumm <- data.frame(dayagg$summary)
#dayaggsumm2 <- data.frame(dayagg2$summary)

fullsumm$model <- "full"
arusumm$model <- "aruonly"
pcsumm$model <- "pconly"
noscoresumm$model <- "noscores"
scoreonlysumm$model <- "scoreonly"
#dayaggsumm$model <- "dayagg"
#dayaggsumm2$model <- "dayagg3"

fullsumm$measure <- rownames(fullsumm)
arusumm$measure <- rownames(arusumm)
pcsumm$measure <- rownames(pcsumm)
noscoresumm$measure <- rownames(noscoresumm)
scoreonlysumm$measure <- rownames(scoreonlysumm)
#dayaggsumm$measure <- rownames(dayaggsumm)
#dayaggsumm2$measure <- rownames(dayaggsumm2)

allmods <- rbind(fullsumm, arusumm, pcsumm, noscoresumm, scoreonlysumm) #, dayaggsumm, dayaggsumm2)
allmods$model <- factor(allmods$model, levels=c("full", "pconly", "aruonly", "noscores", "scoreonly"))

# 

# MODEL COMPARISONS -------------------------------------------------------


# overall psi -------------------------------------------------------------

psi.all <- allmods %>% filter(measure=="psi" | measure=="PropOcc")
psi.plot <- ggplot(psi.all) +
  geom_point(aes(x=model, y=mean, color=model)) +
  geom_errorbar(aes(x=model, ymin=X2.5., ymax=X97.5.), linewidth=0.5) +
  theme_bw() +
  labs(y="posterior p(OCCUPANCY)") +
  geom_hline(yintercept = 0.383, linetype="dashed") + # naive detections by either source
  geom_hline(yintercept = 0.123, linetype="twodash", color = "#CD9600") + # naive dets by PC only
  geom_hline(yintercept=0.346, linetype="twodash", color = "#00BE67") + # naive dets by ARU only
  facet_wrap(~measure, scales = "free_y", ncol=1) +
  scale_y_continuous(breaks=seq(0.1, 0.7, 0.1))

ggsave(plot = psi.plot, filename = "models/figures/psi_NULLmods_BBWO20.png",width = 6, height = 5, units = "in", dpi = 300)

print(jagsResult)
print(jagsData)
#(site_pos <- data_frame(naive.pc = apply(data$y.pc[, 1, , ], 1, max, na.rm=T), naive.aru = apply(data$y.aru[,1:aruVisitLimit], 1, max, na.rm=TRUE)))
site_pos <- data_frame(naive.pc.n = rowSums(data$y.pc[,1,,], na.rm=TRUE), naive.aru.n = rowSums(data$y.aru[,], na.rm=TRUE))
site_pos$naive.pc.pa <- ifelse(site_pos$naive.pc.n==0, 0, 1)
site_pos$naive.aru.pa <- ifelse(site_pos$naive.aru.n==0,0,1)
site_pos$both <- site_pos$naive.pc.pa + site_pos$naive.aru.pa
mean(apply(site_pos, 1, max), na.rm = TRUE) # average number of sites with positive detections from either source

site_pos$either <- ifelse(site_pos$both > 0, 1, 0)
jagsResult$mean$mean_psi # Mean of the posterior of the mean occupancy rate psi across sites from model
site_pos$psi <- jagsResult$mean$psi
site_pos$siteid <- as.numeric(rownames(jagsData$y.ind[,]))
site_pos$species <- speciesCode
site_pos$year <- year
tab <- left_join(latlong, site_pos)
#write.csv(tab, file = "models/output/bbwo20.csv")

jagsResult$mean$PropOcc

ggplot(site_pos) +
  geom_point(aes(x=naive.aru.pa, y=psi)) +
  coord_flip()


det_tbl <- data.frame(mean = colMeans(site_pos[3:6], na.rm=T), 
model=colnames(site_pos[3:6]), 
measure="naivePsi") # how close the pc and aru estimates are to the combined model estimates likely differs by species
# do we have a simulation that gives different z matrices for pc, aru, and show that we return a reasonable psi?
det_tbl
psi_tbl <- psi.all %>% filter(measure=="PropOcc") %>% dplyr::select(mean, model, measure)%>% rbind(det_tbl)
write.csv(psi_tbl, file = "models/output/nullmods_vs_naive_psi.csv")



  
detcovars <- allmods %>% filter(measure=="p_aru11" | measure=="p_aru01")

ggplot(detcovars) +
  geom_point(aes(x=model, y=mean, color=measure)) +
  geom_errorbar(aes(x=model, ymin=X2.5., ymax=X97.5., color=measure)) +
  theme_bw() +
  labs(y="true (p_aru11) and false (p_aru01) detection rates") +
  scale_color_manual(values=c("darkred", "darkgreen"))
